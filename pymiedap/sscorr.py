# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

"""
Single-scattering (TMS / Nakajima-Tanaka) correction for forward-peaked layers
==============================================================================

Strongly forward-peaked particles (large ice crystals, big cloud droplets) need
hundreds of expansion terms, and the doubling-adding solver's cost grows as
~nmug^4 because the Gauss quadrature must satisfy ``nmug >~ M/2``. Brute force is
therefore infeasible (D_eff~60 um ice -> nmug~300-500 -> days/wavelength).

The Nakajima & Tanaka (1988) "TMS" correction breaks the nmug<->M coupling:
run the doubling-adding with an *aggressive* delta-M truncation at small M (so
small, cheap nmug) for the multiple-scattering field, then replace the
truncated single-scattering term with the exact full-phase-matrix single
scattering, computed analytically per geometry::

    R_corrected = R_DAP(delta-M, small M)
                  - R_ss(delta-M-scaled tau', omega', truncated phase matrix)
                  + R_ss(true tau, omega, full phase matrix)

The single-scattering reflection of a homogeneous layer, for unpolarised
incident sunlight [1,0,0,0], expressed in PyMieDAP's *local meridian frame*
(i.e. before the disk-integration beta rotation that geos/read_dap applies):

    [I,Q,U,V]_ss = G * [F11(Theta), F12(Theta) cos2i2, F12(Theta) sin2i2, 0]
    G            = omega / (4 (mu + mu0)) * (1 - exp[-tau (1/mu + 1/mu0)])
    cosTheta     = -mu mu0 + sqrt((1-mu^2)(1-mu0^2)) cos(dphi)
    cos i2       = (mu0 + mu cosTheta) / (sinTheta sqrt(1-mu^2))

This local-frame single-scattering reflection has been validated element by
element (I, Q and U) against a thin-cloud doubling-adding run, confirming the
geometry, the i2 rotation and all sign conventions match PyMieDAP exactly.

The truncated single-scattering matrix is reconstructed from the delta-M
coefficients with ``module_readmie.expand`` (valid for M <= 1000, which an
aggressive truncation easily satisfies); the full single-scattering matrix is
taken from the tabulated phase matrix (e.g. the 498-angle Baum data),
interpolated to Theta, so it is exact regardless of how many terms the full
expansion would need.

Validation status
-----------------
* ``single_scattering_local`` (the building block) is VALIDATED: it reproduces
  a thin-cloud doubling-adding run element-by-element (I, Q, U), confirming the
  geometry, i2 rotation and sign conventions.
* ``tms_correct_local`` (the full assembly) is implemented but NOT yet
  validated against a *converged* reference -- that needs a high-nmug solve
  (nmug ~ ncoef/2 ~ 100+ for a peaked cloud), which only a fast machine can do.
  Use ``examples/validate_tms.py`` on such a machine to confirm it recovers the
  converged result before trusting it for production. Until then treat the TMS
  output as provisional.
"""

import numpy as np

__all__ = [
    "scattering_geometry",
    "single_scattering_local",
    "tms_correct_local",
]


def scattering_geometry(sza_deg, emi_deg, dphi_deg):
    """Return (Theta_deg, cos2i2, sin2i2) for reflection geometry.

    Parameters are solar zenith, viewing (emergent) zenith and azimuth
    difference, all in degrees; arrays are accepted element-wise. The rotation
    angle i2 maps the scattering plane onto the local meridian plane of the
    emergent beam (PyMieDAP convention, validated against the DAP).
    """
    sza = np.radians(np.asarray(sza_deg, float))
    emi = np.radians(np.asarray(emi_deg, float))
    dp = np.radians(np.asarray(dphi_deg, float))
    mu0 = np.cos(sza)
    mu = np.cos(emi)
    s0 = np.sqrt(np.clip(1 - mu0 ** 2, 0, 1))
    sm = np.sqrt(np.clip(1 - mu ** 2, 0, 1))
    cosT = np.clip(-mu * mu0 + sm * s0 * np.cos(dp), -1.0, 1.0)
    Theta = np.degrees(np.arccos(cosT))
    sinT = np.sqrt(np.clip(1 - cosT ** 2, 0, 1))
    denom = sinT * sm
    with np.errstate(divide="ignore", invalid="ignore"):
        cosi2 = np.where(denom > 1e-9,
                         np.clip((mu0 + mu * cosT) / denom, -1, 1), 1.0)
    i2 = np.arccos(cosi2)
    return Theta, np.cos(2 * i2), np.sin(2 * i2)


def single_scattering_local(mu, mu0, F11, F12, tau, ssalb, cos2i2, sin2i2):
    """Local-frame single-scattering reflected Stokes [I,Q,U,V].

    All inputs are scalars or matching arrays. ``F11``/``F12`` are the
    scattering-matrix elements at the scattering angle of each geometry.
    """
    mu = np.asarray(mu, float); mu0 = np.asarray(mu0, float)
    G = ssalb / (4.0 * (mu + mu0)) * (1.0 - np.exp(-tau * (1.0 / mu + 1.0 / mu0)))
    I = G * F11
    Q = G * F12 * cos2i2
    U = G * F12 * sin2i2
    V = np.zeros_like(np.asarray(I, float))
    return I, Q, U, V


def _expand_F(coefs_4x4, ncoef, theta_deg):
    """F11, F12 at theta_deg from PyMieDAP coefficients via module_readmie.expand
    (requires the native module; coefs sliced/padded to (4,4,1001))."""
    import module_readmie as _rd
    c = np.zeros((4, 4, 1001), order="F")
    n = min(int(ncoef), 1000)
    c[:, :, :n + 1] = coefs_4x4[:, :, :n + 1]
    F11 = np.empty_like(np.asarray(theta_deg, float))
    F12 = np.empty_like(F11)
    th = np.atleast_1d(theta_deg)
    f = np.zeros(6)
    for i, t in enumerate(th):
        _rd.expand(n, c, float(t), f)
        F11.flat[i] = f[0]; F12.flat[i] = f[4]
    return F11, F12


def tms_correct_local(I_dap, Q_dap, U_dap, sza_deg, emi_deg, dphi_deg,
                      coefs_trunc, ncoef_trunc, tau_dm, ssalb_dm,
                      theta_full_deg, F11_full, F12_full, tau, ssalb):
    """Apply the TMS single-scattering correction at given local geometries.

    Parameters
    ----------
    I_dap, Q_dap, U_dap : arrays
        Local-frame Stokes from the delta-M doubling-adding run (beta=0, i.e.
        before the disk rotation), one per geometry.
    sza_deg, emi_deg, dphi_deg : arrays
        Geometry of each sample (solar zenith, emergent zenith, azimuth diff).
    coefs_trunc, ncoef_trunc : (4,4,>=1001) array, int
        delta-M-truncated expansion coefficients and order used in the DAP run.
    tau_dm, ssalb_dm : float
        delta-M-scaled optical thickness and single-scattering albedo used in
        the DAP run.
    theta_full_deg, F11_full, F12_full : arrays
        Tabulated *full* (untruncated) phase-matrix elements vs scattering
        angle (e.g. the 498-angle Baum data).
    tau, ssalb : float
        True (unscaled) optical thickness and single-scattering albedo.

    Returns
    -------
    I, Q, U : arrays
        TMS-corrected local-frame Stokes (apply the usual beta rotation /
        disk weighting afterwards, exactly as for the raw DAP output).
    """
    sza = np.asarray(sza_deg, float); emi = np.asarray(emi_deg, float)
    mu0 = np.cos(np.radians(sza)); mu = np.cos(np.radians(emi))
    Theta, c2, s2 = scattering_geometry(sza, emi, dphi_deg)

    # truncated single scattering (matches what the delta-M DAP run contains)
    F11t, F12t = _expand_F(np.asarray(coefs_trunc, float), ncoef_trunc, Theta)
    It, Qt, Ut, _ = single_scattering_local(mu, mu0, F11t, F12t,
                                            tau_dm, ssalb_dm, c2, s2)
    # full single scattering (exact, from the tabulated phase matrix)
    order = np.argsort(theta_full_deg)
    thf = np.asarray(theta_full_deg, float)[order]
    F11f = np.interp(Theta, thf, np.asarray(F11_full, float)[order])
    F12f = np.interp(Theta, thf, np.asarray(F12_full, float)[order])
    If, Qf, Uf, _ = single_scattering_local(mu, mu0, F11f, F12f,
                                            tau, ssalb, c2, s2)

    return (np.asarray(I_dap) - It + If,
            np.asarray(Q_dap) - Qt + Qf,
            np.asarray(U_dap) - Ut + Uf)
