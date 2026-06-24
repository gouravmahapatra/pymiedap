"""
HAPI-based HITRAN absorption cross-section interface.

This module replaces the original Fortran ``salpha.f`` + ``spectrum.f``
pipeline with pure-Python HAPI calls.  HAPI handles all molecules in the
HITRAN database and computes correct partition functions for temperature
scaling via its built-in PYTIPS look-up tables.

HAPI reference:
    Kochanov et al., J. Quant. Spectrosc. Radiat. Transfer 177, 15–30 (2016).
    https://hitran.org/hapi

Usage example::

    from pymiedap.ckdistribution.hitran import fetch_lines, compute_cross_section
    fetch_lines('CO2', 6000, 7000, cache_dir='/tmp/hitran_cache')
    nu, sigma = compute_cross_section('CO2', T=250., P_bar=0.5,
                                      wvn_min=6200., wvn_max=6400.,
                                      wvn_step=0.02, cache_dir='/tmp/hitran_cache')
"""

from __future__ import annotations

import os
import logging
from typing import Callable, Optional, Tuple

import numpy as np

from .gases import get_gas
from .constants import BAR_TO_ATM

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# HAPI import helper
# ---------------------------------------------------------------------------

def _get_hapi():
    """Import HAPI, raising a clear error if it is not installed."""
    try:
        import hapi as _hapi
        return _hapi
    except ImportError:
        raise ImportError(
            "HAPI (HITRAN Application Programming Interface) is required for "
            "absorption cross-section computations but is not installed.\n"
            "Install it with:  pip install hitran-api\n"
            "or download hapi.py from https://hitran.org/hapi and place it on "
            "your Python path."
        )


# ---------------------------------------------------------------------------
# Line-list caching / fetching
# ---------------------------------------------------------------------------

def fetch_lines(
    molecule: str,
    wvn_min: float,
    wvn_max: float,
    cache_dir: str,
    *,
    margin: float = 50.0,
    force_download: bool = False,
    iso_ids: Optional[list] = None,
) -> None:
    """Download HITRAN line parameters for *molecule* if not already cached.

    The line list is stored as ``{cache_dir}/{molecule}.data`` and
    ``{cache_dir}/{molecule}.header``.  Subsequent calls with the same
    molecule and cache directory are no-ops (unless *force_download* is True).

    Parameters
    ----------
    molecule : str
        Molecule name, e.g. ``'CO2'``, ``'CH4'``.  Must be in the gases
        registry (see :mod:`.gases`).
    wvn_min, wvn_max : float
        Wavenumber range of interest [cm⁻¹].  A margin of *margin* cm⁻¹ is
        added on each side to capture pressure-broadened wings.
    cache_dir : str
        Directory where HAPI stores its ``.data`` / ``.header`` files.
    margin : float, optional
        Extra wavenumber range to fetch on each side [cm⁻¹].  Default 50.
    force_download : bool, optional
        Re-download even if the file already exists.
    iso_ids : list of int, optional
        HITRAN *global* isotopologue IDs to include.  Defaults to all
        isotopologues listed in the gases registry.
    """
    hapi = _get_hapi()
    os.makedirs(cache_dir, exist_ok=True)

    # Activate the HAPI working directory
    hapi.db_begin(cache_dir)

    table = molecule.upper()
    data_file = os.path.join(cache_dir, table + '.data')

    if not force_download and os.path.isfile(data_file):
        # Load the table into HAPI's in-memory store if not already there
        try:
            known = hapi.tableList()
        except Exception:
            known = []
        if table not in known:
            log.debug("Loading cached HITRAN table '%s' from %s", table, cache_dir)
            hapi.db_begin(cache_dir)
        else:
            log.debug("Table '%s' already loaded.", table)
        return

    # Need to download
    gas_info = get_gas(molecule)
    ids = iso_ids if iso_ids is not None else gas_info['global_iso_ids']

    log.info(
        "Downloading HITRAN lines for %s (iso_ids=%s) "
        "in wavenumber range [%.1f, %.1f] cm⁻¹ ...",
        molecule, ids, wvn_min - margin, wvn_max + margin
    )
    hapi.fetch_by_ids(table, ids, wvn_min - margin, wvn_max + margin)
    log.info("Download complete → %s", data_file)


# ---------------------------------------------------------------------------
# Cross-section computation
# ---------------------------------------------------------------------------

def compute_cross_section(
    molecule: str,
    T: float,
    P_bar: float,
    wvn_min: float,
    wvn_max: float,
    wvn_step: float,
    cache_dir: str,
    *,
    wing_hw: float = 25.0,
    gamma_l: str = 'gamma_air',
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Voigt absorption cross-section spectrum using HAPI.

    Parameters
    ----------
    molecule : str
        Molecule name (e.g. ``'CO2'``).
    T : float
        Temperature [K].
    P_bar : float
        Total pressure [bar].
    wvn_min, wvn_max : float
        Output wavenumber range [cm⁻¹].
    wvn_step : float
        Output wavenumber step [cm⁻¹].
    cache_dir : str
        Directory containing the downloaded ``.data`` / ``.header`` files.
        Must have been populated via :func:`fetch_lines` first.
    wing_hw : float, optional
        Half-width of the Voigt wing included beyond the output range [cm⁻¹].
        Default 25.
    gamma_l : str, optional
        Lorentzian broadening parameter to use — ``'gamma_air'`` (default) or
        ``'gamma_self'``.

    Returns
    -------
    nu : ndarray, shape (N,)
        Wavenumber grid [cm⁻¹].
    sigma : ndarray, shape (N,)
        Absorption cross-section [cm² molecule⁻¹].

    Raises
    ------
    RuntimeError
        If HAPI returns an empty or all-zero spectrum.
    """
    hapi = _get_hapi()
    hapi.db_begin(cache_dir)

    table = molecule.upper()
    P_atm = float(P_bar) * BAR_TO_ATM

    nu, coef = hapi.absorptionCoefficient_Voigt(
        SourceTables=table,
        Environment={'T': float(T), 'p': float(P_atm)},
        OmegaRange=[float(wvn_min), float(wvn_max)],
        OmegaStep=float(wvn_step),
        OmegaWingHW=float(wing_hw),
        GammaL=gamma_l,
        HITRAN_units=True,   # output in cm²/molecule
    )

    nu    = np.asarray(nu, dtype=float)
    sigma = np.asarray(coef, dtype=float)

    if sigma.size == 0:
        raise RuntimeError(
            f"HAPI returned an empty spectrum for {molecule} at "
            f"T={T} K, P={P_bar} bar, wvn=[{wvn_min}, {wvn_max}] cm⁻¹.\n"
            f"Have you called fetch_lines() first?"
        )

    return nu, sigma


# ---------------------------------------------------------------------------
# Convenience: cross-section averaged over a wavenumber band
# ---------------------------------------------------------------------------

def mean_cross_section(
    molecule: str,
    T: float,
    P_bar: float,
    wvl_um: float,
    bandwidth_um: float,
    cache_dir: str,
    *,
    nw: int = 50000,
    truncw_um: float = 0.2,
    gamma_l: str = 'gamma_air',
) -> float:
    """Band-averaged absorption cross-section [cm² molecule⁻¹].

    Parameters
    ----------
    molecule : str
    T : float  [K]
    P_bar : float  [bar]
    wvl_um : float
        Central wavelength [µm].
    bandwidth_um : float
        Band full-width [µm].
    cache_dir : str
    nw : int
        Number of wavenumber points per µm for the internal grid.
    truncw_um : float
        Truncation margin added to each side of the band [µm].
    gamma_l : str

    Returns
    -------
    float
        Mean cross-section [cm² molecule⁻¹].
    """
    wvl_min = wvl_um - 0.5 * bandwidth_um
    wvl_max = wvl_um + 0.5 * bandwidth_um
    if wvl_min <= 0:
        raise ValueError("wvl_um too small or bandwidth_um too large.")

    wvn_max = 1.0e4 / wvl_min
    wvn_min = 1.0e4 / wvl_max
    wvn_step = 1.0 / nw  # ~0.02 cm⁻¹ at SWIR if nw≈50000/µm

    nu, sigma = compute_cross_section(
        molecule, T, P_bar, wvn_min, wvn_max, wvn_step,
        cache_dir, gamma_l=gamma_l
    )
    return float(np.mean(sigma))


# ---------------------------------------------------------------------------
# Sub-Lorentzian cross-section (line-by-line + χ-factor)
# ---------------------------------------------------------------------------

#: Second radiation constant c₂ = hc/k [cm K]
_C2 = 1.4387769

#: HITRAN standard reference temperature [K]
_T_REF = 296.0


def compute_cross_section_subL(
    molecule: str,
    T: float,
    P_bar: float,
    wvn_min: float,
    wvn_max: float,
    wvn_step: float,
    cache_dir: str,
    *,
    chi_fn: Optional[Callable] = None,
    wing_cm: float = 50.0,
    gamma_l: str = 'gamma_air',
) -> Tuple[np.ndarray, np.ndarray]:
    """Absorption cross-section with optional sub-Lorentzian χ-factor correction.

    Performs a direct line-by-line sum over HITRAN lines, applying the χ-factor
    to suppress Lorentzian far wings.  This is the recommended routine for
    Venus atmosphere modelling where pressures > 1 bar cause the standard
    Voigt/Lorentz profile to overestimate far-wing absorption.

    The cross-section at wavenumber ν is:

    .. code-block:: none

        σ(ν) = Σ_k  S_k(T) · V(ν − ν₀ₖ; γ_D,k, γ_L,k) · χ(|ν − ν₀ₖ|)

    where:

    * ``S_k(T)`` — temperature-scaled HITRAN intensity [cm/molecule]
    * ``V``       — normalised Voigt profile [1/cm⁻¹] via Faddeeva function
    * ``χ``       — sub-Lorentzian correction (``chi_fn``); 1 if ``chi_fn=None``

    Temperature scaling uses the standard HITRAN formula:

    .. code-block:: none

        S(T) = S(T_ref) · Q(T_ref)/Q(T)
                         · exp(−c₂ E" (1/T − 1/T_ref))
                         · [1 − exp(−c₂ ν₀/T)] / [1 − exp(−c₂ ν₀/T_ref)]

    with partition sums from HAPI's built-in PYTIPS tables.

    Parameters
    ----------
    molecule : str
        Absorbing gas (e.g. ``'CO2'``).
    T : float [K]
        Layer temperature.
    P_bar : float [bar]
        Layer total pressure.
    wvn_min, wvn_max : float [cm⁻¹]
        Output wavenumber range.
    wvn_step : float [cm⁻¹]
        Output wavenumber spacing.
    cache_dir : str
        HAPI line-list cache directory (populated via :func:`fetch_lines`).
    chi_fn : callable or None, optional
        χ-factor function — signature ``chi_fn(delta_nu: ndarray) → ndarray``.
        Use :data:`.lineprofiles.chi_tonkov96` for Venus 1.4 µm modelling.
        Pass ``None`` (default) to recover the standard Voigt result.
    wing_cm : float, optional
        Lines within ±*wing_cm* [cm⁻¹] of the output range are included.
        Beyond ~50 cm⁻¹ the Tonkov96 χ suppresses contributions to < 10⁻¹⁵;
        default 50.0 is safe for all built-in χ-factor models.
    gamma_l : str, optional
        Broadening parameter: ``'gamma_air'`` (default) or ``'gamma_self'``.
        For Venus (CO₂ VMR ≈ 0.965) ``'gamma_self'`` is physically appropriate
        but requires self-broadening data in the HITRAN table.

    Returns
    -------
    nu : ndarray [cm⁻¹]
    sigma : ndarray [cm²/molecule]

    Notes
    -----
    This function is ~10–50× slower than :func:`compute_cross_section` because
    it iterates over individual spectral lines in Python.  For typical CKD
    windows (25 cm⁻¹, ~500 grid points, ~500 lines) a single call takes 0.5–2 s
    on a modern CPU.  Numba or Cython could accelerate this by 10–100×.

    See Also
    --------
    :func:`compute_cross_section` : Standard HAPI Voigt (faster, no χ-factor)
    :data:`.lineprofiles.CHI_FACTORS` : Registry of available χ-factor models
    :func:`.lineprofiles.chi_tonkov96` : Tonkov et al. (1996) — default for CO₂
    """
    from .lineprofiles import voigt_profile, doppler_hwhm

    hapi = _get_hapi()
    hapi.db_begin(cache_dir)

    table  = molecule.upper()
    P_atm  = float(P_bar) * BAR_TO_ATM
    T      = float(T)

    # ── Read line parameters from HAPI in-memory cache ────────────────────────
    try:
        cache = hapi.LOCAL_TABLE_CACHE[table]['data']
    except KeyError:
        raise RuntimeError(
            f"HAPI table '{table}' not found in cache '{cache_dir}'. "
            "Call fetch_lines() first."
        )

    nu0      = np.asarray(cache['nu'],       dtype=float)
    S_ref    = np.asarray(cache['sw'],       dtype=float)
    e_lower  = np.asarray(cache['elower'],   dtype=float)
    n_t      = np.asarray(cache['n_air'],    dtype=float)
    molec_id = np.asarray(cache['molec_id'], dtype=int)
    iso_id   = np.asarray(cache['local_iso_id'], dtype=int)

    # Broadening parameter
    if gamma_l == 'gamma_self' and 'gamma_self' in cache:
        gL_ref = np.asarray(cache['gamma_self'], dtype=float)
    else:
        gL_ref = np.asarray(cache['gamma_air'], dtype=float)

    # ── Filter to lines within the relevant wavenumber window ─────────────────
    in_range = (nu0 >= wvn_min - wing_cm) & (nu0 <= wvn_max + wing_cm)
    if not in_range.any():
        nu_grid = np.arange(wvn_min, wvn_max + 0.5 * wvn_step, wvn_step)
        return nu_grid, np.zeros(len(nu_grid))

    nu0      = nu0[in_range]
    S_ref    = S_ref[in_range]
    e_lower  = e_lower[in_range]
    n_t      = n_t[in_range]
    gL_ref   = gL_ref[in_range]
    molec_id = molec_id[in_range]
    iso_id   = iso_id[in_range]

    # ── Temperature-scale line intensities ────────────────────────────────────
    # Boltzmann factor: exp(−c₂ E" (1/T − 1/T_ref))
    boltz = np.exp(-_C2 * e_lower * (1.0 / T - 1.0 / _T_REF))

    # Stimulated-emission factor: [1 − exp(−c₂ ν₀/T)] / [1 − exp(−c₂ ν₀/T_ref)]
    stim = ((1.0 - np.exp(-_C2 * nu0 / T))
            / (1.0 - np.exp(-_C2 * nu0 / _T_REF)))

    # Partition function ratio Q(T_ref) / Q(T)
    # hapi.partitionSum(M, I, T) uses molec_id and local_iso_id
    try:
        Q_ratio = np.array([
            hapi.partitionSum(int(m), int(i), _T_REF) /
            hapi.partitionSum(int(m), int(i), T)
            for m, i in zip(molec_id, iso_id)
        ])
    except Exception as exc:
        log.warning(
            "Partition function lookup failed (%s); falling back to "
            "approximate T^1 scaling.", exc
        )
        Q_ratio = np.full(len(nu0), _T_REF / T)

    S_T = S_ref * Q_ratio * boltz * stim

    # ── Line-shape parameters at (T, P) ───────────────────────────────────────
    # Pressure-broadened Lorentzian HWHM [cm⁻¹]
    gamma_L = P_atm * gL_ref * (_T_REF / T) ** n_t

    # Doppler HWHM [cm⁻¹]
    gas_info  = get_gas(molecule)
    gamma_D   = doppler_hwhm(nu0, T, gas_info['molar_mass'])

    # ── Build output grid and accumulate cross-section ────────────────────────
    nu_grid = np.arange(wvn_min, wvn_max + 0.5 * wvn_step, wvn_step)
    sigma   = np.zeros(len(nu_grid))

    nlines = len(nu0)
    log.debug(
        "compute_cross_section_subL: %d lines, %d grid points, T=%.1f K, P=%.4f bar",
        nlines, len(nu_grid), T, P_bar,
    )

    for k in range(nlines):
        delta   = nu_grid - nu0[k]          # signed detuning [cm⁻¹]
        abs_d   = np.abs(delta)

        # Restrict to lines within wing_cm (some may already be filtered above)
        within = abs_d <= wing_cm
        if not within.any():
            continue

        # Voigt profile on the relevant grid points
        V = voigt_profile(delta[within], gamma_D[k], gamma_L[k])

        # Apply χ-factor if requested
        if chi_fn is not None:
            V = V * chi_fn(abs_d[within])

        sigma[within] += S_T[k] * V

    return nu_grid, sigma
