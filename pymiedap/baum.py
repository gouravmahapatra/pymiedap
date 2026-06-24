# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

"""
Baum / Yang / Heymsfield ice-cloud bulk scattering models for PyMieDAP
======================================================================

The SSEC (University of Wisconsin) ice cloud bulk scattering models give the
single-scattering properties of randomly-oriented, severely-roughened ice
crystal habit mixtures (general habit mixture, solid columns, aggregate of
solid columns) as a function of effective diameter, at 445 wavelengths from
0.2 to 100 um. See Baum et al. (2011), J. Appl. Meteor. Climatol.

The "full phase matrix" files provide six independent matrix elements at 498
scattering angles:

    P11,  P21/P11,  P22/P11,  P33/P11,  P43/P11,  P44/P11

For a randomly-oriented ensemble with a plane of symmetry the scattering matrix
is block-diagonal,

        | F11  F12   0    0  |
    F = | F12  F22   0    0  |
        |  0    0   F33  F34 |
        |  0    0  -F34  F44 |

so the SSEC elements map to PyMieDAP's elements as

    F11 = P11
    F12 = (P21/P11) * P11
    F22 = (P22/P11) * P11
    F33 = (P33/P11) * P11
    F34 = -(P43/P11) * P11        # NOTE the sign: the file gives P43 = -F34
    F44 = (P44/P11) * P11

Unlike spheres, F22 != F11 and F44 != F33, so the spherical-particle expansion
routine bundled with PyMieDAP (``module_mie.devel``) cannot be used directly —
it assumes F22 = F11 and F44 = F33. This module provides the full six-element
generalized-spherical-function expansion (:func:`expand_scattering_matrix`),
which has been validated to reproduce ``module_mie.devel`` to machine precision
on the spherical case while additionally handling the non-spherical elements.

Workflow
--------
    coefs, ncoef = phase_matrix_to_coeffs(theta_deg, F11, F12, F22, F33, F34, F44)
    # or, straight into an Aerosols object ready for a Layer:
    load_baum_into_aerosol(aero, theta_deg, F11, F12, F22, F33, F34, F44, ssalb)

and for large crystals feed the result through
:func:`pymiedap.tmatrix.delta_m_truncate` to keep the doubling-adding solver
within its angular resolution.
"""

import os

import numpy as np

__all__ = [
    "expand_scattering_matrix",
    "phase_matrix_to_coeffs",
    "load_baum_into_aerosol",
    "read_baum_netcdf",
    "convert_baum_netcdf",
    "load_baum_coeffs",
    "fill_aerosol_from_cache",
]

# Match the compiled Fortran array bound (nmatMAX * ... ); PyMieDAP stores
# expansion coefficients in a (4, 4, ncoefsMAX) array with ncoefsMAX = 4001.
NCOEFS_MAX = 4001


# ---------------------------------------------------------------------------
# Generalized spherical functions  P^l_{00}, P^l_{02}, P^l_{22}, P^l_{2,-2}
# ---------------------------------------------------------------------------
# These recurrences are transcribed verbatim from mie_source/devel.f /
# readmie_source/read_mie_output.f (de Haan, Bosma & Hovenier 1987, "the adding
# paper", Eqs. 77-82), so the convention/normalisation matches PyMieDAP exactly.
def _gsf(u, L):
    """Return arrays P[k, l, i] for k in {P00, P02, P22, P2m2}, l=0..L, at u[i]."""
    u = np.asarray(u, dtype=float)
    n = u.size
    qroot6 = -0.25 * np.sqrt(6.0)
    P = np.zeros((4, L + 1, n))
    P00 = np.zeros((n, 2)); P02 = np.zeros((n, 2))
    P22 = np.zeros((n, 2)); P2m2 = np.zeros((n, 2))
    lnew, lold = 0, 1
    sql41 = 0.0
    for l in range(0, L + 1):
        if l == 0:
            P00[:, lold] = 1.0; P00[:, lnew] = 0.0
            P02[:, lold] = 0.0; P22[:, lold] = 0.0; P2m2[:, lold] = 0.0
            P02[:, lnew] = 0.0; P22[:, lnew] = 0.0; P2m2[:, lnew] = 0.0
        else:
            f1 = (2.0 * l - 1.0) / l
            f2 = (l - 1.0) / l
            P00[:, lold] = f1 * u * P00[:, lnew] - f2 * P00[:, lold]
        if l == 2:
            P02[:, lold] = qroot6 * (1.0 - u * u)
            P22[:, lold] = 0.25 * (1.0 + u) ** 2
            P2m2[:, lold] = 0.25 * (1.0 - u) ** 2
            P02[:, lnew] = 0.0; P22[:, lnew] = 0.0; P2m2[:, lnew] = 0.0
            sql41 = 0.0
        elif l > 2:
            sql4 = sql41
            sql41 = np.sqrt(l * l - 4.0)
            twol1 = 2.0 * l - 1.0
            t1 = twol1 / sql41
            t2 = sql4 / sql41
            den = (l - 1.0) * (l * l - 4.0)
            f1 = twol1 * (l - 1.0) * l / den
            f2 = 4.0 * twol1 / den
            f3 = l * ((l - 1.0) ** 2 - 4.0) / den
            P02[:, lold] = t1 * u * P02[:, lnew] - t2 * P02[:, lold]
            P22[:, lold] = (f1 * u - f2) * P22[:, lnew] - f3 * P22[:, lold]
            P2m2[:, lold] = (f1 * u + f2) * P2m2[:, lnew] - f3 * P2m2[:, lold]
        lnew, lold = lold, lnew
        P[0, l] = P00[:, lnew]
        P[1, l] = P02[:, lnew]
        P[2, l] = P22[:, lnew]
        P[3, l] = P2m2[:, lnew]
    return P


def expand_scattering_matrix(u, w, F11, F22, F33, F44, F12, F34, ncoef):
    """Six-element generalized-spherical-function expansion.

    Parameters
    ----------
    u, w : arrays
        Gauss-Legendre nodes (= cos of scattering angle) and weights on [-1, 1].
    F11, F22, F33, F44, F12, F34 : arrays
        Scattering-matrix elements sampled at ``u``.
    ncoef : int
        Highest expansion order to compute (l = 0 .. ncoef).

    Returns
    -------
    coefs : ndarray, shape (4, 4, ncoef+1)
        Expansion coefficients in PyMieDAP layout:
        ``coefs[0,0]=alpha1``, ``[1,1]=alpha2``, ``[2,2]=alpha3``,
        ``[3,3]=alpha4``, ``[0,1]=[1,0]=beta1``, ``[2,3]=-[3,2]=beta2``.
    """
    # On-the-fly GSF recurrence (ping-pong like devel.f) so memory stays O(n)
    # rather than O(n * ncoef) -- essential for the large ncoef (~4000) that
    # the forward diffraction peak of large ice crystals requires.
    u = np.asarray(u, float); w = np.asarray(w, float)
    qroot6 = -0.25 * np.sqrt(6.0)
    wF11 = w * F11; wF44 = w * F44
    wFp = w * (F22 + F33); wFm = w * (F22 - F33)
    wF12 = w * F12; wF34 = w * F34
    n = u.size
    P00 = np.zeros((n, 2)); P02 = np.zeros((n, 2))
    P22 = np.zeros((n, 2)); P2m2 = np.zeros((n, 2))
    a1 = np.zeros(ncoef + 1); a2 = np.zeros(ncoef + 1); a3 = np.zeros(ncoef + 1)
    a4 = np.zeros(ncoef + 1); b1 = np.zeros(ncoef + 1); b2 = np.zeros(ncoef + 1)
    lnew, lold = 0, 1
    sql41 = 0.0
    for l in range(0, ncoef + 1):
        if l == 0:
            P00[:, lold] = 1.0; P00[:, lnew] = 0.0
            P02[:, lold] = 0.0; P22[:, lold] = 0.0; P2m2[:, lold] = 0.0
            P02[:, lnew] = 0.0; P22[:, lnew] = 0.0; P2m2[:, lnew] = 0.0
        else:
            f1 = (2.0 * l - 1.0) / l; f2 = (l - 1.0) / l
            P00[:, lold] = f1 * u * P00[:, lnew] - f2 * P00[:, lold]
        if l == 2:
            P02[:, lold] = qroot6 * (1.0 - u * u)
            P22[:, lold] = 0.25 * (1.0 + u) ** 2
            P2m2[:, lold] = 0.25 * (1.0 - u) ** 2
            P02[:, lnew] = 0.0; P22[:, lnew] = 0.0; P2m2[:, lnew] = 0.0
            sql41 = 0.0
        elif l > 2:
            sql4 = sql41; sql41 = np.sqrt(l * l - 4.0); twol1 = 2.0 * l - 1.0
            t1 = twol1 / sql41; t2 = sql4 / sql41
            den = (l - 1.0) * (l * l - 4.0)
            f1 = twol1 * (l - 1.0) * l / den
            f2 = 4.0 * twol1 / den
            f3 = l * ((l - 1.0) ** 2 - 4.0) / den
            P02[:, lold] = t1 * u * P02[:, lnew] - t2 * P02[:, lold]
            P22[:, lold] = (f1 * u - f2) * P22[:, lnew] - f3 * P22[:, lold]
            P2m2[:, lold] = (f1 * u + f2) * P2m2[:, lnew] - f3 * P2m2[:, lold]
        lnew, lold = lold, lnew
        fl = l + 0.5
        a1[l] = fl * np.dot(P00[:, lnew], wF11)
        a4[l] = fl * np.dot(P00[:, lnew], wF44)
        ap = np.dot(P22[:, lnew], wFp)
        am = np.dot(P2m2[:, lnew], wFm)
        a2[l] = fl * 0.5 * (ap + am)
        a3[l] = fl * 0.5 * (ap - am)
        b1[l] = fl * np.dot(P02[:, lnew], wF12)
        b2[l] = fl * np.dot(P02[:, lnew], wF34)

    coefs = np.zeros((4, 4, ncoef + 1))
    coefs[0, 0] = a1
    coefs[1, 1] = a2
    coefs[2, 2] = a3
    coefs[3, 3] = a4
    coefs[0, 1] = b1; coefs[1, 0] = b1
    coefs[2, 3] = b2; coefs[3, 2] = -b2
    return coefs


def phase_matrix_to_coeffs(theta_deg, F11, F12, F22, F33, F34, F44,
                           ncoef=None, ngauss=None, normalize=True):
    """Convert a tabulated scattering matrix to PyMieDAP expansion coefficients.

    The six elements are given at arbitrary scattering angles ``theta_deg``
    (e.g. the 498 Baum angles). They are interpolated onto Gauss-Legendre nodes
    and expanded with :func:`expand_scattering_matrix`.

    Parameters
    ----------
    theta_deg : array
        Scattering angles in degrees (need not be uniform; 0 and 180 included).
    F11, F12, F22, F33, F34, F44 : arrays
        Absolute scattering-matrix elements (NOT the P_ij/P11 ratios — convert
        first, remembering F34 = -(P43/P11)*P11).
    ncoef : int, optional
        Number of expansion orders. Default: ``2*ngauss - 1`` (capped at
        ``NCOEFS_MAX - 1``).
    ngauss : int, optional
        Number of Gauss-Legendre nodes. Default 2000. Large ice crystals have
        an extremely sharp forward diffraction peak (size parameter ~10^2-10^3),
        so a high node count is needed for the expansion to reproduce the true
        asymmetry parameter; 2000 nodes match the SSEC files' tabulated ``g`` to
        ~10^-3 across 0.2-2 um. The resulting expansion has ncoef ~ 4000 and is
        meant to be fed through :func:`pymiedap.tmatrix.delta_m_truncate` before
        a doubling-adding run.
    normalize : bool
        If True (default), scale all elements so that ``alpha1[0] = 1`` (the
        phase function is normalised), which is PyMieDAP's convention.

    Returns
    -------
    coefs : ndarray (4, 4, ncoef+1)
    ncoef : int
        The actual highest order returned (also the value to store as ncoefs).
    """
    theta_deg = np.asarray(theta_deg, float)
    if ngauss is None:
        ngauss = 2000
    if ncoef is None:
        ncoef = min(2 * ngauss - 1, NCOEFS_MAX - 1)

    # Gauss-Legendre nodes/weights on [-1, 1]; sort everything by mu.
    nodes, weights = np.polynomial.legendre.leggauss(ngauss)
    mu_tab = np.cos(np.radians(theta_deg))
    order = np.argsort(mu_tab)
    mu_tab = mu_tab[order]

    def onto_nodes(Fel):
        Fel = np.asarray(Fel, float)[order]
        # cubic interpolation if scipy is available, else linear
        try:
            from scipy.interpolate import CubicSpline
            return CubicSpline(mu_tab, Fel)(nodes)
        except Exception:
            return np.interp(nodes, mu_tab, Fel)

    f11 = onto_nodes(F11); f12 = onto_nodes(F12); f22 = onto_nodes(F22)
    f33 = onto_nodes(F33); f34 = onto_nodes(F34); f44 = onto_nodes(F44)

    if normalize:
        norm = 0.5 * np.sum(weights * f11)   # = alpha1[0] before scaling
        if norm <= 0:
            raise ValueError("Non-positive phase-function norm; check inputs.")
        f11 /= norm; f12 /= norm; f22 /= norm
        f33 /= norm; f34 /= norm; f44 /= norm

    coefs = expand_scattering_matrix(nodes, weights, f11, f22, f33, f44,
                                     f12, f34, ncoef)
    return coefs, ncoef


def load_baum_into_aerosol(aero, theta_deg, F11, F12, F22, F33, F34, F44,
                           ssalb, sext=None, ncoef=None, ngauss=None):
    """Fill an :class:`~pymiedap.pymiedap.Aerosols` object from Baum data.

    ``theta_deg`` is shared across wavelengths; the six element arrays and
    ``ssalb`` are either 1-D (single wavelength) or 2-D with shape
    ``(nwvl, nangle)`` / ``(nwvl,)``. Sets ``aero.coefs`` (nwvl, 4, 4,
    NCOEFS_MAX), ``aero.ncoefs``, ``aero.ssalb`` and ``aero.ssca`` so that
    ``mix_aerosols`` reproduces the single-scattering albedo.
    """
    F11 = np.atleast_2d(F11); F12 = np.atleast_2d(F12); F22 = np.atleast_2d(F22)
    F33 = np.atleast_2d(F33); F34 = np.atleast_2d(F34); F44 = np.atleast_2d(F44)
    ssalb = np.atleast_1d(ssalb).astype(float)
    nwvl = F11.shape[0]

    coefs_all = np.zeros((nwvl, 4, 4, NCOEFS_MAX), order='F')
    ncoefs_all = np.zeros(nwvl, order='F')
    for z in range(nwvl):
        c, nc = phase_matrix_to_coeffs(theta_deg, F11[z], F12[z], F22[z],
                                       F33[z], F34[z], F44[z],
                                       ncoef=ncoef, ngauss=ngauss)
        coefs_all[z, :, :, :nc + 1] = c
        ncoefs_all[z] = nc

    aero.coefs = coefs_all
    aero.ncoefs = ncoefs_all
    aero.ssalb = ssalb
    aero.sext = np.ones(nwvl) if sext is None else np.asarray(sext, float)
    aero.ssca = aero.ssalb * aero.sext
    aero.f = 1.0
    return aero


# ---------------------------------------------------------------------------
# NetCDF reader
# ---------------------------------------------------------------------------
# Variable names differ slightly between SSEC file releases, so we look up each
# quantity from a list of candidates and fail with a helpful message listing the
# file's actual variables. Finalised/tested against the real files supplied by
# B. Baum (full-phase-matrix, GHM / solid-column / aggregate models).
# Variable names as found in the SSEC full-phase-matrix files (B. Baum, 2013).
# Lists allow tolerance across releases; the first match wins.
_VAR_CANDIDATES = {
    "wavelength":   ["wavelengths", "wavelength", "wnum", "lambda"],
    "deff":         ["effective_diameter", "Deff", "De", "diameter"],
    "angle":        ["phase_angles", "scattering_angle", "angle", "theta"],
    "ssa":          ["single_scattering_albedo", "ssa", "w0"],
    "asym":         ["asymmetry_parameter", "asym", "g"],
    "ext_over_iwc": ["extinction_coefficient_over_iwc", "ext_over_iwc"],
    "P11":          ["p11_phase_function", "p11", "P11"],
    "P21":          ["p21_phase_function", "p21", "P21"],
    "P22":          ["p22_phase_function", "p22", "P22"],
    "P33":          ["p33_phase_function", "p33", "P33"],
    "P43":          ["p43_phase_function", "p43", "P43"],
    "P44":          ["p44_phase_function", "p44", "P44"],
}


def _find_var(ds, key):
    for name in _VAR_CANDIDATES[key]:
        if name in ds.variables:
            return name
    raise KeyError(
        "Could not find a variable for '{}' (tried {}). "
        "Variables in file: {}".format(
            key, _VAR_CANDIDATES[key], list(ds.variables)))


def read_baum_netcdf(path, wavelength_range_um=(0.2, 2.0), deff_um=None):
    """Read a Baum/SSEC full-phase-matrix NetCDF file.

    Parameters
    ----------
    path : str
        Path to the gunzipped full-phase-matrix NetCDF file.
    wavelength_range_um : (lo, hi)
        Keep only wavelengths within this range (default 0.2-2.0 um).
    deff_um : float or None
        If given, return the single nearest effective diameter; otherwise all.

    Returns
    -------
    dict with keys:
        'wavelength_um' (nW,), 'deff_um' (scalar or (nD,)), 'theta_deg' (nA,),
        'ssa', 'asym', 'ext_over_iwc' (each (nW, nD) or (nW,) for a single deff),
        and absolute matrix elements 'F11','F12','F22','F33','F34','F44'
        each (nW, nD, nA) -- or (nW, nA) if a single deff was requested.
        F34 already carries the corrected sign (= -P43 * P11).

    Requires ``netCDF4`` (``pip install netCDF4``).
    """
    try:
        from netCDF4 import Dataset
    except Exception as e:  # pragma: no cover
        raise ImportError("read_baum_netcdf needs netCDF4: pip install netCDF4") from e

    # Transparently handle gzip-compressed files (the SSEC full-phase-matrix
    # files ship as *.nc.gz). netCDF4 reads from an in-memory buffer.
    if str(path).endswith(".gz"):
        import gzip
        with gzip.open(path, "rb") as fh:
            buf = fh.read()
        ds_ctx = Dataset("inmem.nc", mode="r", memory=buf)
    else:
        ds_ctx = Dataset(path)

    with ds_ctx as ds:
        wname = _find_var(ds, "wavelength")
        dname = _find_var(ds, "deff")
        aname = _find_var(ds, "angle")
        wl = np.asarray(ds.variables[wname][:], float)
        de = np.asarray(ds.variables[dname][:], float)
        th = np.asarray(ds.variables[aname][:], float)
        nW, nD, nA = len(wl), len(de), len(th)
        if len({nW, nD, nA}) != 3:
            raise ValueError("wavelength/deff/angle axis sizes are not all "
                             "distinct ({}, {}, {}); cannot orient by size."
                             .format(nW, nD, nA))

        def _axis_order(shape, sizes):
            """Return the permutation that puts axes in the order given by
            `sizes` (matching by length; lengths are all distinct)."""
            shape = list(shape)
            return [shape.index(s) for s in sizes]

        def read2d(key):
            """Read a (deff, wavelength)-type field and return it as (nW, nD)."""
            a = np.asarray(ds.variables[_find_var(ds, key)][:], float)
            return np.transpose(a, _axis_order(a.shape, [nW, nD]))

        def read3d(key):
            """Read a phase-matrix field and return it as (nW, nD, nA)."""
            a = np.asarray(ds.variables[_find_var(ds, key)][:], float)
            return np.transpose(a, _axis_order(a.shape, [nW, nD, nA]))

        ssa = read2d("ssa"); asym = read2d("asym")
        ext_over_iwc = read2d("ext_over_iwc")
        P11 = read3d("P11"); r21 = read3d("P21"); r22 = read3d("P22")
        r33 = read3d("P33"); r43 = read3d("P43"); r44 = read3d("P44")

    sel = np.where((wl >= wavelength_range_um[0]) &
                   (wl <= wavelength_range_um[1]))[0]
    P11 = P11[sel]; r21 = r21[sel]; r22 = r22[sel]
    r33 = r33[sel]; r43 = r43[sel]; r44 = r44[sel]
    F11 = P11
    F12 = r21 * P11
    F22 = r22 * P11
    F33 = r33 * P11
    F34 = -r43 * P11           # sign correction: file stores P43 = -F34
    F44 = r44 * P11
    out = dict(wavelength_um=wl[sel], deff_um=de, theta_deg=th,
               ssa=ssa[sel], asym=asym[sel], ext_over_iwc=ext_over_iwc[sel],
               F11=F11, F12=F12, F22=F22, F33=F33, F34=F34, F44=F44)

    if deff_um is not None:
        j = int(np.argmin(np.abs(de - deff_um)))
        out["deff_um"] = float(de[j])
        for k in ("F11", "F12", "F22", "F33", "F34", "F44"):
            out[k] = out[k][:, j, :]
        for k in ("ssa", "asym", "ext_over_iwc"):
            out[k] = out[k][:, j]
    return out


# ---------------------------------------------------------------------------
# One-shot converter: NetCDF (.nc / .nc.gz)  ->  PyMieDAP coefficient cache
# ---------------------------------------------------------------------------
def convert_baum_netcdf(nc_path, out_npz=None, deff_um=60.0,
                        wavelength_range_um=(0.2, 2.0), ngauss=2000,
                        ncoef=None, verbose=True):
    """Read a Baum/SSEC full-phase-matrix file and expand every wavelength
    (at one effective diameter) into PyMieDAP expansion coefficients, saved as
    a compressed ``.npz`` cache that :func:`load_baum_coeffs` /
    :func:`fill_aerosol_from_cache` can read back instantly.

    Parameters
    ----------
    nc_path : str
        Path to the ``*FullPhaseMatrix.nc`` or ``.nc.gz`` file.
    out_npz : str, optional
        Output cache path. Default: ``<ncbase>_Deff<deff>.npz`` next to input.
    deff_um : float
        Effective diameter to extract (nearest grid value is used).
    wavelength_range_um : (lo, hi)
        Wavelength window to convert (default 0.2-2.0 um).
    ngauss, ncoef : int
        Expansion resolution (see :func:`phase_matrix_to_coeffs`).

    Returns
    -------
    out_npz : str
        Path to the written cache. The archive holds: ``wavelength_um`` (nW,),
        ``deff_um`` (scalar), ``ncoefs`` (nW,), ``coefs`` (nW, 4, 4, Lmax+1),
        ``ssalb`` (nW,), ``asym`` (nW,), ``ext_over_iwc`` (nW,), ``habit`` (str).
    """
    d = read_baum_netcdf(nc_path, wavelength_range_um=wavelength_range_um,
                         deff_um=deff_um)
    wl = d["wavelength_um"]
    th = d["theta_deg"]
    nW = wl.size
    habit = os.path.basename(str(nc_path)).split("_")[0]

    coefs_list = []
    ncoefs = np.zeros(nW, dtype=int)
    g_expand = np.zeros(nW)
    Lmax = 0
    for z in range(nW):
        c, nc = phase_matrix_to_coeffs(th, d["F11"][z], d["F12"][z], d["F22"][z],
                                       d["F33"][z], d["F34"][z], d["F44"][z],
                                       ncoef=ncoef, ngauss=ngauss)
        coefs_list.append(c)
        ncoefs[z] = nc
        g_expand[z] = c[0, 0, 1] / 3.0
        Lmax = max(Lmax, nc)
        if verbose and (z % 25 == 0 or z == nW - 1):
            print("  [{:3d}/{}] lambda={:.3f} um  g_expand={:.4f} g_file={:.4f}"
                  .format(z + 1, nW, wl[z], g_expand[z], d["asym"][z]))

    coefs = np.zeros((nW, 4, 4, Lmax + 1))
    for z, c in enumerate(coefs_list):
        coefs[z, :, :, :c.shape[2]] = c

    if out_npz is None:
        base = os.path.basename(str(nc_path))
        for suff in (".nc.gz", ".nc"):
            if base.endswith(suff):
                base = base[: -len(suff)]
                break
        out_npz = os.path.join(os.path.dirname(os.path.abspath(str(nc_path))),
                               "{}_Deff{:.0f}.npz".format(base, d["deff_um"]))

    np.savez_compressed(out_npz, wavelength_um=wl, deff_um=d["deff_um"],
                        ncoefs=ncoefs, coefs=coefs, ssalb=d["ssa"],
                        asym=d["asym"], asym_expand=g_expand,
                        ext_over_iwc=d["ext_over_iwc"], habit=habit)
    if verbose:
        print("Wrote {}  ({} wavelengths, Deff={:.0f} um, max g-error={:.4f})"
              .format(out_npz, nW, d["deff_um"],
                      np.max(np.abs(g_expand - d["asym"]))))
    return out_npz


def load_baum_coeffs(npz_path, wavelengths_um=None):
    """Load a cache written by :func:`convert_baum_netcdf`.

    If ``wavelengths_um`` is given, the nearest cached wavelength is selected
    for each requested value and the arrays are reordered/subset accordingly.

    Returns a dict with ``wavelength_um``, ``deff_um``, ``ncoefs``, ``coefs``
    (nW, 4, 4, Lmax+1), ``ssalb``, ``asym``.
    """
    z = np.load(npz_path, allow_pickle=True)
    out = {k: z[k] for k in z.files}
    if wavelengths_um is not None:
        wl = out["wavelength_um"]
        idx = [int(np.argmin(np.abs(wl - w))) for w in np.atleast_1d(wavelengths_um)]
        for k in ("wavelength_um", "ncoefs", "coefs", "ssalb", "asym",
                  "asym_expand", "ext_over_iwc"):
            if k in out and np.ndim(out[k]) >= 1 and out[k].shape[0] == wl.size:
                out[k] = out[k][idx]
    return out


def fill_aerosol_from_cache(aero, cache, wavelengths_um=None):
    """Fill an :class:`~pymiedap.pymiedap.Aerosols` object from a Baum cache
    (path or dict). Sets ``coefs`` (nW, 4, 4, NCOEFS_MAX), ``ncoefs``,
    ``ssalb`` and ``ssca`` ready for a :class:`Layer`. Apply
    :func:`pymiedap.tmatrix.delta_m_truncate` afterwards before a DAP run.
    """
    if isinstance(cache, str):
        cache = load_baum_coeffs(cache, wavelengths_um=wavelengths_um)
    elif wavelengths_um is not None:
        wl = cache["wavelength_um"]
        idx = [int(np.argmin(np.abs(wl - w))) for w in np.atleast_1d(wavelengths_um)]
        cache = {k: (v[idx] if (np.ndim(v) >= 1 and getattr(v, "shape", [0])[0] == wl.size) else v)
                 for k, v in cache.items()}
    c = np.asarray(cache["coefs"], float)
    nW, _, _, L = c.shape
    coefs_all = np.zeros((nW, 4, 4, NCOEFS_MAX), order='F')
    coefs_all[:, :, :, :L] = c
    aero.coefs = coefs_all
    aero.ncoefs = np.asarray(cache["ncoefs"], float)
    aero.ssalb = np.asarray(cache["ssalb"], float)
    aero.sext = np.ones(nW)
    aero.ssca = aero.ssalb * aero.sext
    aero.f = 1.0
    return aero
