# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

"""
Delta-fit truncation for strongly forward-peaked phase matrices
===============================================================

Background
----------
The doubling-adding solver can represent an expansion of at most ~2*nmug terms.
For a sharp forward peak (large ice crystals, big droplets) the true expansion
needs hundreds-to-thousands of terms, so it must be truncated to a small order
M (=> small, cheap nmug). The bundled delta-M truncation
(:func:`pymiedap.tmatrix.delta_m_truncate`) sets the truncated coefficients by
*moment projection* (orthogonality integrals over the full angular range). When
the peak is sharp, that short projected series rings (Gibbs) and goes negative
in the backscatter hemisphere -- which makes the doubling diverge (NaN) and,
where it does run, converges only slowly in the polarized signal.

Delta-fit (Hu, Wielicki, Lin et al. 2000) replaces the projection with a
weighted least-squares **fit** of the short expansion to the tabulated phase
matrix over the *backscatter hemisphere only* (scattering angle >= theta_cut),
deliberately ignoring the forward spike. The removed forward energy is folded
into a delta-scaling of the optical thickness and single-scattering albedo,
exactly as in delta-M. Because the fit is not forced to chase the spike, the
truncated phase function stays smooth and non-negative, and is accurate exactly
where reflected-light radiative transfer needs it.

This module implements the scalar (F11) fit first (Phase 1); the vector
extension to all six elements is added in :func:`deltafit_matrix` (Phase 2).

Conventions
-----------
PyMieDAP stores F11(mu) = sum_l a_l P_l(mu) with P_l the Legendre polynomials
and a_0 = 1 (the phase function is normalised so (1/2) int F11 dmu = 1, and
a_1 = 3 g). The fitted coefficients a_l are exactly these (a_l = alpha1_l).
"""

import numpy as np

__all__ = ["legendre_design", "deltafit_scalar", "deltafit_matrix",
           "fill_aerosol_deltafit", "deltafit_from_baum"]

NCOEFS_MAX = 4001   # match the compiled Fortran / pymiedap.baum


def legendre_design(mu, M):
    """Return the (len(mu), M+1) matrix of Legendre polynomials P_l(mu),
    l = 0..M, via the stable three-term recurrence."""
    mu = np.asarray(mu, float)
    A = np.zeros((mu.size, M + 1))
    A[:, 0] = 1.0
    if M >= 1:
        A[:, 1] = mu
    for l in range(1, M):
        A[:, l + 1] = ((2 * l + 1) * mu * A[:, l] - l * A[:, l - 1]) / (l + 1)
    return A


def deltafit_scalar(theta_deg, F11, M, theta_cut=5.0, weight="iso",
                    enforce_nonneg=True):
    """Delta-fit a forward-peaked phase function F11 to order M.

    Parameters
    ----------
    theta_deg : array
        Scattering angles [deg] (the tabulated grid, e.g. the 498 Baum angles).
    F11 : array
        Phase function, normalised so (1/2) int F11 dmu = 1 (the Baum files are).
    M : int
        Truncation order (number of Legendre terms - 1). Keep M <= 2*nmug.
    theta_cut : float
        Forward-peak cutoff [deg]; only angles >= theta_cut are fitted.
    weight : "iso" | "sin" | array
        Fit weights over the fitted angles. "iso" = uniform in mu; "sin" =
        proportional to sin(theta) (solid-angle weighting); or an explicit array.
    enforce_nonneg : bool
        If True, iteratively reweight/clip so the reconstructed F11 stays >= 0
        on the full grid (a few passes; cheap).

    Returns
    -------
    a : ndarray (M+1,)
        Fitted, renormalised Legendre coefficients (a_0 = 1); these are the
        PyMieDAP alpha1_l for the truncated phase function.
    f : float
        Truncated forward-energy fraction (= 1 - a0_raw), for delta-scaling
        tau' = (1 - f*omega) tau and omega' = (1-f) omega / (1 - f*omega).
    diag : dict
        Diagnostics: 'min_full' (min of reconstructed F11 over all angles),
        'rms_back' (RMS fractional residual over the fitted region), 'g'.
    """
    theta_deg = np.asarray(theta_deg, float)
    F11 = np.asarray(F11, float)
    mu = np.cos(np.radians(theta_deg))
    fitmask = theta_deg >= theta_cut

    if isinstance(weight, str):
        if weight == "sin":
            w = np.sin(np.radians(theta_deg))
        else:                      # "iso": uniform in mu
            w = np.ones_like(theta_deg)
    else:
        w = np.asarray(weight, float)
    w = np.where(fitmask, w, 0.0)

    Afull = legendre_design(mu, M)          # (nang, M+1) over all angles
    A = Afull[fitmask]
    y = F11[fitmask]
    ww = w[fitmask]

    def solve(weights):
        sw = np.sqrt(np.maximum(weights, 0.0))
        coef, *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
        return coef

    a_raw = solve(ww)
    # iterative non-negativity: up-weight angles where the reconstruction
    # dips below zero and refit (smooth, avoids hard constraints).
    if enforce_nonneg:
        wcur = ww.copy()
        for _ in range(8):
            recon_full = Afull @ a_raw
            if recon_full.min() >= 0:
                break
            # penalise the (fitted) angles closest to the negative dips
            neg = np.clip(-recon_full[fitmask], 0, None)
            wcur = wcur + 50.0 * neg / (np.abs(y).mean() + 1e-30)
            a_raw = solve(wcur)

    a0 = a_raw[0]
    f = 1.0 - a0                              # forward energy removed
    a = a_raw / a0                            # renormalise so alpha1_0 = 1

    recon = Afull @ a_raw
    rms_back = float(np.sqrt(np.mean(((recon[fitmask] - y) /
                                      (y + 1e-30)) ** 2)))
    diag = dict(min_full=float(recon.min()), rms_back=rms_back,
                g=float(a[1] / 3.0) if M >= 1 else 0.0, f=float(f), a0=float(a0))
    return a, f, diag


def deltafit_matrix(theta_deg, F11, F12, F22, F33, F34, F44, M,
                    theta_cut=5.0, weight="iso"):
    """Vector delta-fit: least-squares-fit all six independent scattering-matrix
    elements to order M over the backscatter hemisphere.

    Each element is expanded in the generalised spherical function appropriate
    to it (the same basis the doubling-adding solver uses): F11, F44 in P^l_00;
    F12, F34 in P^l_02; (F22 +/- F33) in P^l_{2,+/-2}. F11 is fit with the
    non-negativity-preserving scalar routine; all elements are normalised by the
    same factor so alpha1_0 = 1 and the matrix stays mutually consistent.

    Parameters
    ----------
    theta_deg : array
        Scattering angles [deg].
    F11..F44 : arrays
        Absolute matrix elements (F11 normalised so (1/2) int F11 dmu = 1;
        F34 = -P43*P11 sign already applied, as read_baum_netcdf returns).
    M, theta_cut, weight : see :func:`deltafit_scalar`.

    Returns
    -------
    coefs : ndarray (4, 4, M+1)
        PyMieDAP expansion coefficients of the truncated matrix (alpha1_0 = 1):
        [0,0]=alpha1 [1,1]=alpha2 [2,2]=alpha3 [3,3]=alpha4
        [0,1]=[1,0]=beta1  [2,3]=-[3,2]=beta2.
    f : float
        Truncated forward-energy fraction (delta-scaling, from the F11 fit).
    diag : dict
        Diagnostics from the F11 fit plus per-element backscatter RMS.
    """
    from .baum import _gsf      # P00,P02,P22,P2m2 GSF recurrence (validated)

    theta_deg = np.asarray(theta_deg, float)
    mu = np.cos(np.radians(theta_deg))
    fitmask = theta_deg >= theta_cut
    if isinstance(weight, str):
        w = np.sin(np.radians(theta_deg)) if weight == "sin" else np.ones_like(theta_deg)
    else:
        w = np.asarray(weight, float)
    P = _gsf(mu, M)                       # (4, M+1, n): 0=P00 1=P02 2=P22 3=P2m2
    sw = np.sqrt(np.maximum(np.where(fitmask, w, 0.0), 0.0))[fitmask]

    def fit(basis, y):
        A = P[basis].T[fitmask]           # (nfit, M+1)
        c, *_ = np.linalg.lstsq(A * sw[:, None], np.asarray(y, float)[fitmask] * sw,
                                rcond=None)
        return c

    # F11: non-negativity-preserving scalar fit -> alpha1 (normalised), f, a0.
    a1, f, diag = deltafit_scalar(theta_deg, F11, M, theta_cut, weight)
    a0 = diag["a0"]
    # Other elements: same backscatter fit, divided by a0 to share F11's norm.
    a4 = fit(0, F44) / a0
    b1 = fit(1, F12) / a0
    b2 = fit(1, F34) / a0
    sp = fit(2, np.asarray(F22) + np.asarray(F33)) / a0
    sm = fit(3, np.asarray(F22) - np.asarray(F33)) / a0
    a2 = 0.5 * (sp + sm)
    a3 = 0.5 * (sp - sm)

    coefs = np.zeros((4, 4, M + 1))
    coefs[0, 0] = a1
    coefs[1, 1] = a2
    coefs[2, 2] = a3
    coefs[3, 3] = a4
    coefs[0, 1] = b1; coefs[1, 0] = b1
    coefs[2, 3] = b2; coefs[3, 2] = -b2
    return coefs, f, diag


def fill_aerosol_deltafit(aero, theta_deg, F11, F12, F22, F33, F34, F44,
                          ssalb, M, theta_cut=5.0, weight="iso"):
    """Delta-fit a tabulated phase matrix into an Aerosols object.

    Drop-in alternative to delta-M: fits the matrix at order M (stable, non-
    negative), sets ``aero.coefs``/``ncoefs``/``ssalb``/``ssca``, and returns
    the per-wavelength optical-thickness scale ``(1 - f*omega)`` to multiply the
    layer tau with (the single-scattering albedo is delta-scaled to
    ``omega' = (1-f) omega / (1 - f*omega)``).

    F-elements are (nwvl, nang) or (nang,); ``ssalb`` is (nwvl,) or scalar.
    Inputs need not be pre-normalised -- F11 is normalised internally and the
    other elements scaled by the same factor.
    """
    F11 = np.atleast_2d(F11); F12 = np.atleast_2d(F12); F22 = np.atleast_2d(F22)
    F33 = np.atleast_2d(F33); F34 = np.atleast_2d(F34); F44 = np.atleast_2d(F44)
    ssalb = np.atleast_1d(ssalb).astype(float)
    nwvl = F11.shape[0]
    mu = np.cos(np.radians(np.asarray(theta_deg, float)))
    order = np.argsort(mu)

    coefs_all = np.zeros((nwvl, 4, 4, NCOEFS_MAX), order="F")
    ncoefs = np.zeros(nwvl, order="F")
    ssalb_new = np.zeros(nwvl)
    tau_scale = np.ones(nwvl)
    for z in range(nwvl):
        sc = 0.5 * np.trapz(F11[z][order], mu[order])     # normalise phase fn
        c, f, _ = deltafit_matrix(theta_deg, F11[z] / sc, F12[z] / sc,
                                  F22[z] / sc, F33[z] / sc, F34[z] / sc,
                                  F44[z] / sc, M, theta_cut, weight)
        coefs_all[z, :, :, :M + 1] = c
        ncoefs[z] = M
        w = ssalb[z]
        ssalb_new[z] = (1.0 - f) * w / (1.0 - f * w)
        tau_scale[z] = 1.0 - f * w
    aero.coefs = coefs_all
    aero.ncoefs = ncoefs
    aero.ssalb = ssalb_new
    aero.sext = np.ones(nwvl)
    aero.ssca = aero.ssalb * aero.sext
    aero.f = 1.0
    return tau_scale


def deltafit_from_baum(aero, baum_nc, wavelengths_um, deff_um, M,
                       theta_cut=5.0, weight="iso", wavelength_range_um=(0.2, 2.0)):
    """Read a Baum full-phase-matrix NetCDF and delta-fit its ice into ``aero``.

    Convenience wrapper: reads the tabulated matrix (handles .nc / .nc.gz) at
    the requested effective diameter, selects the nearest cached wavelength for
    each requested one, and calls :func:`fill_aerosol_deltafit`. Returns the
    per-wavelength tau-scale factor.
    """
    from .baum import read_baum_netcdf
    d = read_baum_netcdf(baum_nc, wavelength_range_um=wavelength_range_um,
                         deff_um=deff_um)
    wl = d["wavelength_um"]
    idx = [int(np.argmin(np.abs(wl - w)))
           for w in np.atleast_1d(wavelengths_um)]
    th = d["theta_deg"]
    pick = lambda k: d[k][idx]
    return fill_aerosol_deltafit(aero, th, pick("F11"), pick("F12"),
                                 pick("F22"), pick("F33"), pick("F34"),
                                 pick("F44"), d["ssa"][idx], M, theta_cut, weight)
