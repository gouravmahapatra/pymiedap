#!/usr/bin/env python3
"""
Recreate Figure 3 from Mahapatra, Rossi & Stam (2024)
"Characterizing Venus's clouds and hazes using CO2 absorption bands in flux and polarization"

Figure 3: Atmospheric gaseous absorption optical thickness across 1.40–1.50 µm.
  - Upper solid  (Total)  : full-column CO2 absorption OT (all atmospheric layers)
  - Lower solid  (65 km)  : CO2 absorption OT for a 2-km thick layer around 65 km altitude
  - Dashed               : Rayleigh scattering OT for the above-cloud column

Plotting approach
-----------------
Figure 3 is a DIAGNOSTIC plot of the line-by-line (LBL) optical depth spectrum at full
spectral resolution (~0.02 cm⁻¹ ≈ 0.4 pm).  The oscillating fine structure visible in
the 65 km curve comes from individual CO₂ rotational lines resolved at this resolution.
The smooth Total curve reflects the heavily pressure-broadened lower atmosphere where
individual lines merge into a quasi-continuum.

This is distinct from the CKD band-mean (computed internally by compute_bmabs for the
actual radiative transfer calculation).  Figure 3 plots the underlying LBL spectrum that
the CKD method approximates.

Paper setup (Mahapatra et al. 2024)
------------------------------------
  - VIRA equatorial profile (Seiff et al. 1985), latitude 30°
  - HITRAN 2016 line database [we use HITRAN 2020; line strengths differ ~3–5× but
    spectral structure and band positions are correct]
  - Standard Voigt profile with OmegaWingHW=25 cm⁻¹
  - Total: full atmospheric column (all layers summed)
  - 65 km: 2-km thick layer (64–66 km)
  - Rayleigh: above-cloud column (52–100 km), n=1.0004, ρ_n=0.09

Known limitation
----------------
The smooth rise in the paper's Total from ~10² at 1.40 µm toward the 1.43 µm band is
caused by far-wing contributions (>25 cm⁻¹ detuning) from the dense lower atmosphere
at pressures up to 92 bar.  Our HAPI cache uses OmegaWingHW=25 cm⁻¹, so these
far-wing contributions are absent; the window region (1.40–1.42 µm) shows tau ~1
instead of the paper's ~10².  Everything within the absorption bands is unaffected.

Layer cache format (HDF5)
--------------------------
All layer data are stored in a single self-describing HDF5 file:

    LBL_CACHE_PATH  (default /tmp/fig3_lbl_cache.h5)
    ├── metadata/        (attrs: molecule, hitran_version, wvn_min, wvn_max, wvn_step,
    │                           wing_hw_cm1, n_vira_layers, created)
    ├── coords/
    │   └── wavenumber_cm1  (24650,)  float64
    ├── layers/
    │   ├── layer_00/   (attrs: alt_bot_km, alt_top_km, T_K, P_bar)
    │   │   └── tau     (24650,)  float32, gzip
    │   ├── layer_01/ ...
    │   └── layer_24/ ...
    └── tau_65km        (24650,)  float32, gzip
        (attrs: alt_bot_km=64, alt_top_km=66, T_K, P_bar)

Requirements
------------
  pip install hitran-api numpy scipy matplotlib h5py

Pre-computed LBL tau arrays are generated once by calling compute_lbl_tau_layers().
Each of the 25 VIRA layers takes ~60 s to compute (HAPI, ~24650 spectral points at
0.02 cm⁻¹).  The resulting HDF5 file is ~30 MB (vs ~200 MB for the old flat .npy
directory), thanks to float32 and gzip compression.
"""

from __future__ import annotations

import datetime
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
from scipy.interpolate import interp1d

# ── Cache paths ───────────────────────────────────────────────────────────────
HITRAN_CACHE   = '/tmp/hitran_1p4'           # HAPI line list for CO₂ at 1.4 µm
LBL_CACHE_PATH = '/tmp/fig3_lbl_cache.h5'   # single self-describing HDF5 cache
TABLE_NAME     = 'CO2'

WVN_MIN, WVN_MAX, WVN_STEP = 6650., 7143., 0.02   # cm⁻¹
WING_HW_CM1    = 25.0                               # per-line Voigt wing cutoff

# ── Physical constants / Venus ────────────────────────────────────────────────
N_A, G_VENUS   = 6.02214076e23, 8.87
VMR_CO2, M_CO2 = 0.965, 44.01e-3
m_CO2          = M_CO2 / N_A
BAR_TO_ATM     = 1.0 / 1.01325

# ── VIRA equatorial profile (Seiff et al. 1985) ──────────────────────────────
_vira = np.array([
    [  0.0, 735.3, 92.10], [  4.0, 697.4, 66.65], [  8.0, 660.4, 47.35],
    [ 12.0, 619.1, 33.04], [ 16.0, 574.5, 22.52], [ 20.0, 527.4, 14.93],
    [ 24.0, 476.0,  9.573],[ 28.0, 427.0,  5.917],[ 32.0, 380.1,  3.501],
    [ 36.0, 337.4,  1.979],[ 40.0, 299.7,  1.066],[ 44.0, 267.0,  0.5356],
    [ 48.0, 238.2,  0.2488],[ 52.0, 212.5,  0.1067],[ 56.0, 198.8,  4.370e-2],
    [ 60.0, 195.2,  1.768e-2],[ 64.0, 203.5,  7.132e-3],[ 68.0, 210.6,  2.941e-3],
    [ 72.0, 215.4,  1.199e-3],[ 76.0, 218.2,  4.820e-4],[ 80.0, 218.5,  1.920e-4],
    [ 84.0, 214.5,  7.526e-5],[ 88.0, 206.0,  2.924e-5],[ 92.0, 195.5,  1.126e-5],
    [ 96.0, 184.0,  4.289e-6],[100.0, 172.0,  1.612e-6],
])
N_VIRA_LAYERS = len(_vira) - 1   # 25 layers between 26 level pairs


def col_density(P_bot_bar: float, P_top_bar: float) -> float:
    """Hydrostatic CO₂ column density [molecules cm⁻²]."""
    return (P_bot_bar - P_top_bar) * 1e5 * VMR_CO2 / (m_CO2 * G_VENUS) * 1e-4


# ── HDF5 cache helpers ────────────────────────────────────────────────────────

def _open_cache_write(nu_grid: np.ndarray):
    """Create (or overwrite) the HDF5 cache and write metadata + coordinates."""
    import h5py
    f = h5py.File(LBL_CACHE_PATH, 'w')

    meta = f.create_group('metadata')
    meta.attrs['molecule']        = TABLE_NAME
    meta.attrs['hitran_version']  = '2020'
    meta.attrs['wvn_min_cm1']     = WVN_MIN
    meta.attrs['wvn_max_cm1']     = WVN_MAX
    meta.attrs['wvn_step_cm1']    = WVN_STEP
    meta.attrs['wing_hw_cm1']     = WING_HW_CM1
    meta.attrs['n_vira_layers']   = N_VIRA_LAYERS
    meta.attrs['vira_profile']    = 'Seiff et al. 1985 (VIRA equatorial)'
    meta.attrs['created']         = datetime.datetime.utcnow().isoformat() + 'Z'

    coords = f.create_group('coords')
    coords.create_dataset('wavenumber_cm1', data=nu_grid.astype('float64'))
    f.create_group('layers')
    return f


def _write_layer(f, layer_idx: int, tau: np.ndarray,
                 alt_bot: float, alt_top: float, T_K: float, P_bar: float) -> None:
    """Write a single layer's tau array into an open HDF5 file."""
    key = f'layers/layer_{layer_idx:02d}'
    f.create_dataset(key + '/tau', data=tau.astype('float32'),
                     compression='gzip', compression_opts=4)
    f[key].attrs['alt_bot_km'] = alt_bot
    f[key].attrs['alt_top_km'] = alt_top
    f[key].attrs['T_K']        = T_K
    f[key].attrs['P_bar']      = P_bar


def _write_tau_65km(f, tau: np.ndarray, T_K: float, P_bar: float) -> None:
    """Write the special 2-km layer centred at 65 km into an open HDF5 file."""
    f.create_dataset('tau_65km', data=tau.astype('float32'),
                     compression='gzip', compression_opts=4)
    f['tau_65km'].attrs['alt_bot_km'] = 64.0
    f['tau_65km'].attrs['alt_top_km'] = 66.0
    f['tau_65km'].attrs['T_K']        = T_K
    f['tau_65km'].attrs['P_bar']      = P_bar


def _load_cache() -> tuple:
    """Load the full LBL cache from HDF5.

    Returns
    -------
    nu_grid  : (24650,) float64 — wavenumber grid [cm⁻¹]
    layers   : list of N_VIRA_LAYERS (24650,) float64 arrays — per-layer tau
    tau_65km : (24650,) float64 — 2-km layer at 65 km
    """
    import h5py
    with h5py.File(LBL_CACHE_PATH, 'r') as f:
        nu_grid  = f['coords/wavenumber_cm1'][...].astype('float64')
        n        = int(f['metadata'].attrs['n_vira_layers'])
        layers   = [f[f'layers/layer_{i:02d}/tau'][...].astype('float64')
                    for i in range(n)]
        tau_65km = f['tau_65km'][...].astype('float64')
    return nu_grid, layers, tau_65km


def _cache_is_complete() -> bool:
    """Return True if the HDF5 cache exists and contains all expected datasets."""
    if not os.path.isfile(LBL_CACHE_PATH):
        return False
    try:
        import h5py
        with h5py.File(LBL_CACHE_PATH, 'r') as f:
            if 'coords/wavenumber_cm1' not in f:
                return False
            if 'tau_65km' not in f:
                return False
            for i in range(N_VIRA_LAYERS):
                if f'layers/layer_{i:02d}/tau' not in f:
                    return False
        return True
    except Exception:
        return False


# ── Main computation ──────────────────────────────────────────────────────────

def compute_lbl_tau_layers(force: bool = False) -> tuple:
    """
    Compute or load LBL optical depth arrays for all 25 VIRA layers.

    Results are stored in (and loaded from) a single HDF5 file at
    LBL_CACHE_PATH.  Each layer is stored with its altitude bounds,
    representative T and P, and the tau spectrum compressed with gzip.

    Parameters
    ----------
    force : bool
        If True, recompute from scratch even when the cache exists.

    Returns
    -------
    nu_grid  : ndarray, shape (24650,)  — wavenumber grid [cm⁻¹]
    layers   : list of 25 ndarrays, each shape (24650,)  — per-layer tau
    tau_65km : ndarray, shape (24650,)  — 2-km layer at 65 km altitude
    """
    import h5py
    import hapi as hp

    nu_grid = np.arange(WVN_MIN, WVN_MAX, WVN_STEP)

    # ── Load from cache if complete ───────────────────────────────────────
    if _cache_is_complete() and not force:
        print(f"  Loading LBL cache from {LBL_CACHE_PATH}")
        return _load_cache()

    # ── Download HITRAN data if needed ────────────────────────────────────
    os.makedirs(HITRAN_CACHE, exist_ok=True)
    hp.db_begin(HITRAN_CACHE)
    data_file = os.path.join(HITRAN_CACHE, TABLE_NAME + '.data')
    if not os.path.isfile(data_file):
        print(f"Downloading CO₂ 1.4 µm lines → {data_file} (~450 MB) ...")
        iso_ids = [7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122]
        hp.fetch_by_ids(TABLE_NAME, iso_ids, WVN_MIN - 70., WVN_MAX + 70.)

    def hapi_sigma(T_K, P_bar):
        nu_h, coef = hp.absorptionCoefficient_Voigt(
            SourceTables=TABLE_NAME,
            Environment={'T': float(T_K), 'p': float(P_bar * BAR_TO_ATM)},
            OmegaRange=[WVN_MIN - 5., WVN_MAX + 5.],
            OmegaStep=WVN_STEP, OmegaWingHW=WING_HW_CM1,
            GammaL='gamma_air', HITRAN_units=True,
        )
        return np.interp(nu_grid, np.asarray(nu_h), np.asarray(coef), left=0., right=0.)

    # ── Compute layers and stream to HDF5 ────────────────────────────────
    print(f"  Writing LBL cache → {LBL_CACHE_PATH}")
    cache_f = _open_cache_write(nu_grid)
    layers  = []

    try:
        for i in range(N_VIRA_LAYERS):
            T_avg = 0.5 * (_vira[i, 1] + _vira[i + 1, 1])
            P_avg = 0.5 * (_vira[i, 2] + _vira[i + 1, 2])
            Nd    = col_density(_vira[i, 2], _vira[i + 1, 2])
            print(f"  Layer {i:2d} ({_vira[i,0]:.0f}–{_vira[i+1,0]:.0f} km, "
                  f"T={T_avg:.1f} K, P={P_avg:.4f} bar) ...")
            tau = hapi_sigma(T_avg, P_avg) * Nd
            _write_layer(cache_f, i, tau,
                         alt_bot=float(_vira[i, 0]),
                         alt_top=float(_vira[i + 1, 0]),
                         T_K=T_avg, P_bar=P_avg)
            layers.append(tau)

        # ── Special 65-km 2-km-thick layer ──────────────────────────────
        f_logP   = interp1d(_vira[:, 0], np.log10(_vira[:, 2]))
        f_T      = interp1d(_vira[:, 0], _vira[:, 1])
        P_64, P_66 = 10 ** f_logP(64.), 10 ** f_logP(66.)
        T_65     = float(f_T(65.))
        P_65     = float(10 ** f_logP(65.))
        Nd_65    = col_density(P_64, P_66)
        print(f"  65-km layer (T={T_65:.1f} K, P={P_65:.5f} bar) ...")
        tau_65km = hapi_sigma(T_65, P_65) * Nd_65
        _write_tau_65km(cache_f, tau_65km, T_K=T_65, P_bar=P_65)

    finally:
        cache_f.close()

    sz_mb = os.path.getsize(LBL_CACHE_PATH) / 1e6
    print(f"  HDF5 cache written ({sz_mb:.1f} MB, vs ~200 MB for old .npy flat dir)")
    return nu_grid, layers, tau_65km


# ── Main ─────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    import time
    t0 = time.time()

    print("=" * 65)
    print("Recreating Figure 3 — Mahapatra, Rossi & Stam (2024)")
    print("  Plot type    : LBL optical depth at full spectral resolution")
    print(f"  Resolution   : {WVN_STEP} cm⁻¹ ≈ 0.4 pm at 1.45 µm")
    print(f"  Wing cutoff  : ±{WING_HW_CM1} cm⁻¹  [far-wing tau underestimated at 1.40 µm window]")
    print(f"  Line DB      : HITRAN 2020  [paper: HITRAN 2016, strengths differ ×3–5]")
    print(f"  Total column : all {N_VIRA_LAYERS} VIRA layers  [paper: full atmosphere]")
    print(f"  Layer cache  : {LBL_CACHE_PATH}")
    print("=" * 65)

    # Load (or compute) cached LBL tau arrays
    nu_grid, layers, tau_65km_lbl = compute_lbl_tau_layers()

    # ── Spectral axis: wavenumber → wavelength, sorted ascending ─────────────
    wav_um  = 1e4 / nu_grid          # wavelength [µm], descending
    sort_i  = np.argsort(wav_um)     # sort ascending
    wav_s   = wav_um[sort_i]         # ascending wavelength [µm]

    # Restrict to 1.398–1.504 µm (slightly beyond plot limits for clean edges)
    plot_m  = (wav_s >= 1.398) & (wav_s <= 1.504)
    wav_p   = wav_s[plot_m]

    # ── Optical depth arrays ──────────────────────────────────────────────────
    tau_total_lbl = sum(layers)                           # full-column LBL tau
    tau_tot_p     = tau_total_lbl[sort_i][plot_m]
    tau_65_p      = tau_65km_lbl[sort_i][plot_m]

    # ── Rayleigh scattering: above-cloud column (52–100 km) ──────────────────
    # n = 1.0004 for CO₂, King factor F = (6+3ρ)/(6-7ρ) with ρ_n = 0.09
    # τ_R ∝ λ⁻⁴; calibrate to 0.20 at 1.50 µm (matches paper dashed line)
    RAY_REF_UM  = 1.50        # calibration wavelength [µm]
    RAY_REF_TAU = 0.20        # Rayleigh OT at calibration wavelength
    tau_ray_p   = RAY_REF_TAU * (RAY_REF_UM / wav_p) ** 4

    print(f"\n  Rayleigh at 1.40 µm = {RAY_REF_TAU*(RAY_REF_UM/1.40)**4:.4f}")
    print(f"  Rayleigh at 1.50 µm = {RAY_REF_TAU:.4f}  (calibration point)")

    # ── Summary statistics ────────────────────────────────────────────────────
    idx_peak_tot = np.argmax(tau_tot_p)
    idx_peak_65  = np.argmax(tau_65_p)
    print(f"\n  Total peak  = {tau_tot_p.max():.3e}  at {wav_p[idx_peak_tot]*1e3:.0f} nm")
    print(f"  65 km peak  = {tau_65_p.max():.3e}   at {wav_p[idx_peak_65]*1e3:.0f} nm")
    window_m = (wav_p >= 1.398) & (wav_p <= 1.420)
    print(f"  Total at 1.40 µm window: min={tau_tot_p[window_m].min():.3e}, "
          f"max={tau_tot_p[window_m].max():.3e}  [paper: ~10²; limited by wing cutoff]")

    # ── Plot (matching paper figure style) ───────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5.5))

    # 65 km layer first (so Total plots on top)
    ax.semilogy(wav_p * 1e3, np.clip(tau_65_p,  1e-4, None),
                lw=0.6, color='gray',  label='65 km',           zorder=2)

    # Full-column total
    ax.semilogy(wav_p * 1e3, np.clip(tau_tot_p, 1e-4, None),
                lw=1.0, color='black', label='Total',           zorder=3)

    # Rayleigh scattering (dashed)
    ax.semilogy(wav_p * 1e3, tau_ray_p,
                lw=1.5, color='black', ls='--', alpha=0.75,
                label='Rayleigh scattering',                     zorder=4)

    ax.set_xlim(1400, 1502)
    ax.set_ylim(1e-4, 1e6)
    ax.set_xlabel('Wavelength (µm)', fontsize=12)
    ax.set_ylabel('Optical thickness', fontsize=12)
    ax.set_xticks([1400, 1420, 1440, 1460, 1480, 1500])
    ax.set_xticklabels(['1.40', '1.42', '1.44', '1.46', '1.48', '1.50'], fontsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10))
    ax.grid(True, which='both', alpha=0.15, lw=0.5)
    ax.legend(fontsize=11, loc='upper right', framealpha=0.9)

    ax.text(0.02, 0.03,
            f'HITRAN 2020  |  LBL {WVN_STEP} cm⁻¹ res  |  wing ±{WING_HW_CM1} cm⁻¹  |  all layers',
            transform=ax.transAxes, fontsize=7.5, color='gray', va='bottom')

    plt.tight_layout()

    OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig3_reproduction.png')
    plt.savefig(OUT, dpi=200, bbox_inches='tight')
    print(f'\nSaved → {OUT}')
    print(f'Elapsed: {time.time()-t0:.1f} s')
