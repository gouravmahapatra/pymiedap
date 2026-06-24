"""
Core correlated k-distribution (CKD) logic.

This module provides the pure-Python CKD machinery that replaces the original
monolithic ``ckdis.py`` script.  All functions are stateless (no global
variables) and accept parameters explicitly.

Workflow
--------
1. :func:`slitfunction` — build a wavenumber grid for a wavelength window.
2. :func:`compute_cross_section_grid` — call HAPI via :mod:`.hitran`.
3. :func:`ksort` — sort the cross-section into a k-distribution.
4. :func:`gauss_legendre_points` — get integration nodes and weights.
5. :func:`kspec_layer` — run steps 1–4 for every wavelength, returning the
   k-values at the Gauss–Legendre nodes for a single atmospheric layer.
6. :func:`compute_bmabs` — integrate over all layers to produce the
   band-mean absorption optical depth ``bmabs`` used by PyMieDAP.

Unit conventions
----------------
- Wavelength  : µm
- Wavenumber  : cm⁻¹
- Cross-section : cm² molecule⁻¹
- Pressure    : bar (converted to atm internally before HAPI calls)
- Temperature : K
- bmabs       : dimensionless (molecules m⁻² × cm² molecule⁻¹ × 1e-4)
"""

from __future__ import annotations

import logging
from typing import Callable, Optional, Tuple

import numpy as np
from scipy import interpolate

from .hitran import compute_cross_section, fetch_lines
from .gases import get_gas

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Unit conversion helpers
# ---------------------------------------------------------------------------

def wvl2wvn(wavelength_um: float | np.ndarray) -> float | np.ndarray:
    """Convert wavelength [µm] → wavenumber [cm⁻¹]."""
    return 1.0e4 / np.asarray(wavelength_um, dtype=float)


def wvn2wvl(wavenumber_cm: float | np.ndarray) -> float | np.ndarray:
    """Convert wavenumber [cm⁻¹] → wavelength [µm]."""
    return 1.0e4 / np.asarray(wavenumber_cm, dtype=float)


# ---------------------------------------------------------------------------
# Slit function (wavelength window → wavenumber grid)
# ---------------------------------------------------------------------------

def slitfunction(
    wavel0: float,
    sigma: float,
    nw: int,
    truncw: float,
) -> Tuple[int, float, np.ndarray, np.ndarray, float, float, float, float]:
    """Build a wavenumber grid and integrated slit-function weights.

    Parameters
    ----------
    wavel0 : float
        Central wavelength [µm].
    sigma : float
        Band full-width [µm].
    nw : int
        Number of wavenumber samples per µm (spectral resolution control).
    truncw : float
        Extra wavelength margin beyond the band edges [µm] used to capture
        line-shape wings in the cross-section computation.

    Returns
    -------
    nnv : int
        Number of wavenumber grid points inside [vmin, vmax].
    truncv : float
        Truncation half-width in wavenumbers [cm⁻¹].
    specv : ndarray, shape (nnv,)
        Wavenumber grid [cm⁻¹] for the slit interval.
    speci : ndarray, shape (nnv,)
        Cumulative normalised wavelength-width weights (used for IRF-weighted
        k-distribution sorting).
    vmin, vmax : float
        Min / max wavenumber of the slit interval [cm⁻¹].
    tmin, tmax : float
        Min / max wavenumber including the truncation wing [cm⁻¹].
    """
    wavelmin = wavel0 - 0.5 * sigma
    wavelmax = wavel0 + 0.5 * sigma

    if wavelmin <= 0.0:
        raise ValueError(
            f"wavelmin = {wavelmin:.4f} µm is out of bounds. "
            f"Reduce sigma ({sigma}) or increase wavel0 ({wavel0})."
        )

    nnv = int((wavelmax - wavelmin) * nw)

    # Wavelength → wavenumber (note: max wavel → min wavenumber)
    vmin = float(wvl2wvn(wavelmax))
    vmax = float(wvl2wvn(wavelmin))

    # Truncation in wavenumbers
    truncv = abs(1.0e4 * (1.0 / (wavelmax + truncw) - 1.0 / wavelmax))
    if truncv <= 0.0:
        truncv = 1.0e4 / wavelmax

    tmin = vmin - truncv
    tmax = vmax + truncv

    # Build the wavenumber grid
    specv_list = []
    i = 0
    specvi = 0.0
    while specvi < vmax:
        specvi = vmin + (vmax - vmin) * i / max(nnv, 1)
        specv_list.append(specvi)
        i += 1
    nnv = i
    specv = np.array(specv_list, dtype=float)
    specw = wvn2wvl(specv)

    # Build cumulative normalised wavelength-width array
    tot = 0.0
    speci = np.zeros(nnv)
    for i in range(nnv - 1, 0, -1):
        w1 = specw[i - 1]
        w2 = specw[i]
        plus = w1 - w2
        tot += plus
        speci[i] = tot
    if tot > 0.0:
        speci /= tot

    return nnv, truncv, specv, speci, vmin, vmax, tmin, tmax


# ---------------------------------------------------------------------------
# Gauss–Legendre quadrature on [0, 1]
# ---------------------------------------------------------------------------

def gauss_legendre_points(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return *n*-point Gauss–Legendre nodes and weights on [0, 1].

    The standard ``np.polynomial.legendre.leggauss`` returns nodes on [-1, 1].
    This function rescales to [0, 1] and returns the weights unchanged (the
    caller must scale by (b-a)/2 = 0.5 when integrating).

    Parameters
    ----------
    n : int
        Number of quadrature points.

    Returns
    -------
    gp : ndarray, shape (n,)
        Gauss–Legendre nodes on [0, 1].
    gw : ndarray, shape (n,)
        Gauss–Legendre weights (sum = 2; caller scales by 0.5).
    """
    x, w = np.polynomial.legendre.leggauss(n)
    gp = 0.5 * (x + 1.0)   # map [-1, 1] → [0, 1]
    return gp, w


# ---------------------------------------------------------------------------
# k-sort (absorption cross-section → CK distribution)
# ---------------------------------------------------------------------------

def ksort(
    spec: np.ndarray,
    specv: np.ndarray,
    irf_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sort an absorption cross-section spectrum into a k-distribution.

    Parameters
    ----------
    spec : ndarray, shape (N,)
        Absorption cross-section values [cm² molecule⁻¹] on the wavenumber
        grid *specv*.
    specv : ndarray, shape (N,)
        Wavenumber grid [cm⁻¹].
    irf_fn : callable or None
        Instrument response function.  If provided, it must accept *specv*
        and return a non-negative weight array of the same length.  When None,
        the slit is treated as a flat top-hat (uniform weights).

    Returns
    -------
    xscale : ndarray, shape (N+1,)
        Normalised g-space x-axis (uniform spacing, range [0, 1]).
    xscale_irf : ndarray, shape (N+1,)
        Normalised g-space x-axis weighted by the IRF.
    abs_spec_sorted : ndarray, shape (N+1,)
        Absorption cross-sections sorted by increasing magnitude.
    """
    sort_idx = np.argsort(spec)
    N = len(spec)

    # Sorted cross-sections (prepend a zero for the lower bound)
    abs_spec_sorted = np.empty(N + 1)
    abs_spec_sorted[0] = 0.0
    abs_spec_sorted[1:] = spec[sort_idx]

    # Uniform x-scale (flat slit)
    xwidths = np.ones(N) * abs(specv[0] - specv[1]) if N > 1 else np.ones(N)
    xscale = np.zeros(N + 1)
    running = 0.0
    for i, dw in enumerate(xwidths):
        running += dw
        xscale[i + 1] = running
    total = xscale[-1]
    if total > 0.0:
        xscale /= total

    # IRF-weighted x-scale
    if irf_fn is not None:
        irf_vals = irf_fn(specv)
        irf_vals = np.abs(irf_vals)
        irf_max = irf_vals.max()
        if irf_max > 0.0:
            irf_vals /= irf_max
        xwidths_irf = xwidths * irf_vals
        xwidths_irf_sorted = xwidths_irf[sort_idx]
    else:
        xwidths_irf_sorted = xwidths[sort_idx]

    xscale_irf = np.zeros(N + 1)
    running = 0.0
    for i, dw in enumerate(xwidths_irf_sorted):
        running += dw
        xscale_irf[i + 1] = running
    total_irf = xscale_irf[-1]
    if total_irf > 0.0:
        xscale_irf /= total_irf

    return xscale, xscale_irf, abs_spec_sorted


# ---------------------------------------------------------------------------
# k-spectrum for a single atmospheric layer
# ---------------------------------------------------------------------------

def kspec_layer(
    wav: np.ndarray,
    pres_bar: float,
    temp_K: float,
    gauss_points: np.ndarray,
    molecule: str,
    cache_dir: str,
    sigma_um: float = 0.1,
    truncw_um: float = 0.2,
    nw: int = 50000,
    irf_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    use_irf: bool = False,
    chi_fn: Optional[Callable] = None,
    wing_cm: float = 50.0,
) -> np.ndarray:
    """Compute the CK-distributed absorption cross-sections at Gauss nodes.

    For each central wavelength in *wav*, this function:

    1. Builds the wavenumber grid via :func:`slitfunction`.
    2. Computes σ(ν) via HAPI (:func:`.hitran.compute_cross_section`) or the
       sub-Lorentzian routine (:func:`.hitran.compute_cross_section_subL`) when
       *chi_fn* is set.
    3. Sorts into a k-distribution via :func:`ksort`.
    4. Interpolates onto *gauss_points*.

    Parameters
    ----------
    wav : ndarray, shape (W,)
        Central wavelengths [µm].
    pres_bar : float
        Layer pressure [bar].
    temp_K : float
        Layer temperature [K].
    gauss_points : ndarray, shape (G,)
        Gauss–Legendre nodes on [0, 1] from :func:`gauss_legendre_points`.
    molecule : str
        Absorbing gas (e.g. ``'CO2'``).
    cache_dir : str
        HAPI line-list cache directory.
    sigma_um : float
        Band full-width for each spectral window [µm].
    truncw_um : float
        Truncation margin [µm].
    nw : int
        Wavenumber resolution [samples per µm].
    irf_fn : callable or None
        Instrument response function.  Signature: ``irf_fn(specv) → weights``.
    use_irf : bool
        If True, interpolate using the IRF-weighted x-scale; otherwise use the
        flat x-scale.
    chi_fn : callable or None, optional
        Sub-Lorentzian χ-factor function — ``chi_fn(delta_nu) → ndarray``.
        When set, :func:`.hitran.compute_cross_section_subL` is used instead of
        the standard HAPI Voigt routine.  Recommended for layers with
        P > 1 bar (Venus lower atmosphere).
        See :data:`.lineprofiles.CHI_FACTORS` for available models.
    wing_cm : float, optional
        Far-wing truncation [cm⁻¹] for the sub-Lorentzian sum; ignored when
        *chi_fn* is ``None``.  Default 50.

    Returns
    -------
    layer_kdis : ndarray, shape (W, G)
        k-distributed cross-sections [cm² molecule⁻¹] at each Gauss node for
        each wavelength.
    """
    from .hitran import compute_cross_section_subL

    layer_kdis = np.zeros((len(wav), len(gauss_points)))

    for iw, wavel0 in enumerate(wav):
        # 1. Build the wavenumber grid for this window
        nnv, truncv, specv, speci, vmin, vmax, tmin, tmax = slitfunction(
            wavel0, sigma_um, nw, truncw_um
        )

        # 2. Compute wvn_step.
        wvn_range = vmax - vmin
        wvn_step  = max(wvn_range / 500.0, 0.005)  # never finer than 0.005 cm⁻¹

        # 3. Compute absorption cross-section
        try:
            if chi_fn is not None:
                # Sub-Lorentzian LBL sum with χ-factor correction
                nu, sigma = compute_cross_section_subL(
                    molecule, temp_K, pres_bar,
                    wvn_min=tmin, wvn_max=tmax,
                    wvn_step=wvn_step,
                    cache_dir=cache_dir,
                    chi_fn=chi_fn,
                    wing_cm=wing_cm,
                )
            else:
                # Standard HAPI Voigt (faster)
                nu, sigma = compute_cross_section(
                    molecule, temp_K, pres_bar,
                    wvn_min=tmin, wvn_max=tmax,
                    wvn_step=wvn_step,
                    cache_dir=cache_dir,
                )
            # Slice to the slit interval [vmin, vmax]
            mask = (nu >= vmin) & (nu <= vmax)
            sigma_slit = sigma[mask]
            specv_slit = nu[mask]
        except Exception as exc:
            log.warning(
                "kspec_layer: cross-section failed at λ=%.5f µm, P=%.4f bar, "
                "T=%.1f K: %s.  Using zeros.",
                wavel0, pres_bar, temp_K, exc
            )
            sigma_slit = np.zeros(nnv)
            specv_slit = specv

        if len(sigma_slit) == 0:
            log.warning(
                "kspec_layer: empty spectrum at λ=%.5f µm. Using zeros.", wavel0
            )
            sigma_slit = np.zeros(nnv)
            specv_slit = specv

        # 4. Sort into k-distribution
        xscale, xscale_irf, abs_sorted = ksort(sigma_slit, specv_slit, irf_fn)

        # 5. Choose the appropriate x-scale for interpolation
        x_interp = xscale_irf if use_irf else xscale

        # 6. Interpolate onto Gauss nodes
        gp_clipped = np.clip(gauss_points, x_interp[0], x_interp[-1])
        k_fn = interpolate.interp1d(x_interp, abs_sorted, kind='linear',
                                    bounds_error=False, fill_value=(0.0, abs_sorted[-1]))
        layer_kdis[iw, :] = k_fn(gp_clipped)

    return layer_kdis


# ---------------------------------------------------------------------------
# Atmosphere-level bmabs computation
# ---------------------------------------------------------------------------

def compute_bmabs(
    atmosphere,
    molecule: str,
    wav: np.ndarray,
    gauss_points: np.ndarray,
    cache_dir: str,
    sigma_um: float = 0.1,
    truncw_um: float = 0.2,
    nw: int = 50000,
    irf_fn: Optional[Callable] = None,
    use_irf: bool = False,
    verbose: bool = True,
    chi_fn: Optional[Callable] = None,
    wing_cm: float = 50.0,
) -> np.ndarray:
    """Compute band-mean absorption optical depth (bmabs) for all layers.

    This is the main entry point for the CKD pipeline.  It loops over all
    atmospheric levels, calls :func:`kspec_layer` for each, averages adjacent
    levels into layers, and multiplies by the column density to obtain bmabs.

    Parameters
    ----------
    atmosphere : :class:`.atmosphere.Atmosphere`
        Atmospheric profile.  The ``gas_vmr`` dict may contain VMR profiles;
        if *molecule* is absent from ``gas_vmr``, a VMR of 1.0 is assumed.
    molecule : str
        Absorbing gas, e.g. ``'CO2'``.
    wav : ndarray, shape (W,)
        Central wavelengths [µm].
    gauss_points : ndarray, shape (G,)
        Gauss–Legendre nodes on [0, 1].
    cache_dir : str
        HAPI line-list cache directory.
    sigma_um : float
        Spectral window full-width [µm].
    truncw_um : float
        Truncation margin [µm].
    nw : int
        Wavenumber resolution [samples per µm].
    irf_fn : callable or None
        Instrument response function; passed to :func:`kspec_layer`.
    use_irf : bool
        Whether to use the IRF-weighted k-scale.
    verbose : bool
        Print progress to stdout.
    chi_fn : callable or None, optional
        Sub-Lorentzian χ-factor — passed through to :func:`kspec_layer`.
        When set, :func:`.hitran.compute_cross_section_subL` is used in place
        of the standard HAPI Voigt routine.

        Recommended usage for Venus modelling::

            from pymiedap.ckdistribution.lineprofiles import chi_tonkov96
            bmabs = compute_bmabs(..., chi_fn=chi_tonkov96)

        Use ``chi_fn=None`` (default) for the standard Voigt profile.
    wing_cm : float, optional
        Far-wing truncation [cm⁻¹] for the sub-Lorentzian sum; ignored when
        *chi_fn* is ``None``.  Default 50.

    Returns
    -------
    bmabs_kdis : ndarray, shape (nlayer, W, G)
        Band-mean absorption optical depth at each Gauss node, for each
        wavelength, for each atmospheric layer.  ``nlayer = nlev - 1``.

        * ``bmabs_kdis[:, :, 0]``  — lowest Gauss point (weakest absorption)
        * ``bmabs_kdis[:, :, -1]`` — highest Gauss point (strongest absorption,
          used to reproduce Figure 3 of Mahapatra et al. 2024)
    """
    gas_info = get_gas(molecule)
    mol_mass_g = gas_info['molar_mass']

    # Ensure line list is fetched
    wvn_max = float(wvl2wvn(min(wav) - 0.5 * sigma_um))
    wvn_min = float(wvl2wvn(max(wav) + 0.5 * sigma_um))
    fetch_lines(molecule, wvn_min, wvn_max, cache_dir)

    nlev = atmosphere.nlev
    nlay = nlev - 1
    nwav = len(wav)
    ng   = len(gauss_points)

    mode_str = f"sub-Lorentzian χ={chi_fn.__name__}" if chi_fn is not None else "Voigt"
    if verbose:
        print(f"  Profile mode: {mode_str}")

    # Cross-sections at every level: shape (nlev, nwav, ng)
    sigma_levels = np.zeros((nlev, nwav, ng))

    for i in range(nlev):
        pres = float(atmosphere.pressure[i])
        temp = float(atmosphere.temperature[i])
        if verbose:
            print(f"  Level {i+1:3d}/{nlev}  P={pres:.4f} bar  T={temp:.1f} K")
        sigma_levels[i] = kspec_layer(
            wav, pres, temp, gauss_points, molecule, cache_dir,
            sigma_um=sigma_um, truncw_um=truncw_um, nw=nw,
            irf_fn=irf_fn, use_irf=use_irf,
            chi_fn=chi_fn, wing_cm=wing_cm,
        )

    # Build bmabs for each layer (average of top + bottom cross-sections × Nd)
    bmabs_kdis = np.zeros((nlay, nwav, ng))

    for i in range(nlay):
        # Layer-averaged cross-section [cm² molecule⁻¹]
        sigma_avg = 0.5 * (sigma_levels[i] + sigma_levels[i + 1])

        # Column density [molecules m⁻²]
        Nd = atmosphere.column_density(molecule, i, mol_mass_g)

        # bmabs [dimensionless]: σ [cm²/mol] × Nd [mol/m²] → convert cm²→m²: ×1e-4
        bmabs_kdis[i] = sigma_avg * Nd * 1.0e-4

    return bmabs_kdis
