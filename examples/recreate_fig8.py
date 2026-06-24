#!/usr/bin/env python3
"""Reproduce Figure 8 (top row) of Mahapatra, Rossi & Stam (2024).

Reflects sunlight from a Venus-like atmosphere in the 1.40–1.50 µm CO₂
absorption band.  Two models are computed:

  1. "Only cloud"    – cloud aerosol slab (ba=30, top at 65 km)
  2. "Cloud + haze"  – same cloud + haze slab (ba=0.1, top at 67 km)

Geometry : nadir viewing (θ=0°), SZA θ₀=30°, φ−φ₀=0°
Spectral  : 1.40–1.50 µm, box 1 nm bins, 10 Gauss–Legendre points
Gas       : CO₂ via HITRAN 2016 correlated-k distribution
Particles : mode-2 cloud (r_g=1.05 µm, σ_g=1.21, n=1.40)
            mode-1 haze  (r_g=0.15 µm, σ_g=1.91, n=1.40)
Surface   : black (asurf=0)

Speed note: the CKD step calls HAPI ONCE per layer for the full
1.40–1.50 µm range, then sorts into bins (numpy), rather than
one HAPI call per (layer × bin).  This reduces 3600 HAPI calls
to 8 (one per coarsened atmospheric layer).
"""

from __future__ import annotations

import os
import sys
import warnings
import time
warnings.filterwarnings('ignore')

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np
from scipy import interpolate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Physical / model constants  (matching the paper exactly)
# ---------------------------------------------------------------------------

WAV_MIN_UM   = 1.400
WAV_MAX_UM   = 1.500
BIN_WIDTH_UM = 0.001      # 1 nm bins
N_GAUSS      = 10         # Gauss–Legendre points per bin

SZA_DEG   = 30.0
EMISS_DEG = 0.0
AZ_DEG    = 0.0
# Phase angle = SZA for nadir viewing (the angle at the planet between sun and observer).
# NOT 180-SZA: that would be the scattering angle, which is a different quantity.
# Confirmed by venus_co2_hires_spectropol.py: PHASE = SZA = 45° for nadir+SZA=45°.
PHASE_DEG = SZA_DEG   # = 30°

N_R_CLOUD = 1.40
N_I_CLOUD = 0.0
N_R_HAZE  = 1.40
N_I_HAZE  = 0.0

# Mode-2 cloud: log-normal r_g=1.05 µm, σ_g=1.21 → effective gamma params
_sg_c = np.log(1.21)
CLOUD_R_EFF = 1.05 * np.exp(2.5 * _sg_c**2)   # ≈ 1.15 µm
CLOUD_V_EFF = np.exp(_sg_c**2) - 1              # ≈ 0.037
CLOUD_TAU   = 30.0
CLOUD_TOP_KM = 65.0
CLOUD_BOT_KM = 47.0

# Mode-1 haze: log-normal r_g=0.15 µm, σ_g=1.91 → effective gamma params
_sg_h = np.log(1.91)
HAZE_R_EFF = 0.15 * np.exp(2.5 * _sg_h**2)    # ≈ 0.43 µm
HAZE_V_EFF = np.exp(_sg_h**2) - 1               # ≈ 0.52
HAZE_TAU    = 0.1
HAZE_TOP_KM = 67.0    # 2 km above cloud top
HAZE_BOT_KM = CLOUD_TOP_KM   # haze sits directly on top of cloud

VMR_CO2  = 0.965
G_VENUS  = 8.87       # m s⁻²
M_CO2_KG = 44.01e-3 / 6.022e23   # kg/molecule
N_A      = 6.022e23
N_CO2_N  = 1.0004    # refractive index CO₂ at 1.4–1.5 µm
RHO_CO2  = 0.09      # depolarization factor

NMU_MIE  = 16
NMU_DAP  = 20
NSUBR    = 20

HITRAN_CACHE = '/tmp/hitran_fig8'          # OK in /tmp: re-downloadable
DAP_DIR      = Path('/tmp/fig8_dap')       # OK in /tmp: working files only
OUT_PNG      = REPO_ROOT / 'examples' / 'mahapatra2024_fig8_toprow.png'

# ---------------------------------------------------------------------------
# Persistent cache — one HDF5 file per unique parameter set
# ---------------------------------------------------------------------------
# Stored in .cache/fig8/ inside the repo so it survives workspace resets.
# The filename embeds an 8-char SHA256 fingerprint of every model parameter
# that affects the result, so changing SZA / particle sizes / etc. always
# produces a fresh computation rather than silently reusing stale data.

CACHE_DIR = REPO_ROOT / '.cache' / 'fig8'


def _params_fingerprint() -> str:
    """Return an 8-char hex digest that uniquely identifies these model params."""
    import hashlib, json
    params = dict(
        wav_min=WAV_MIN_UM, wav_max=WAV_MAX_UM, bin_width=BIN_WIDTH_UM,
        n_gauss=N_GAUSS, sza=SZA_DEG, emiss=EMISS_DEG, az=AZ_DEG, phase=PHASE_DEG,
        n_r_cloud=N_R_CLOUD,  n_i_cloud=N_I_CLOUD,
        cloud_r_eff=round(CLOUD_R_EFF, 6), cloud_v_eff=round(CLOUD_V_EFF, 6),
        cloud_tau=CLOUD_TAU,
        n_r_haze=N_R_HAZE, n_i_haze=N_I_HAZE,
        haze_r_eff=round(HAZE_R_EFF, 6), haze_v_eff=round(HAZE_V_EFF, 6),
        haze_tau=HAZE_TAU,
        n_layers=N_LAYERS, layer_bounds=list(LAYER_BOUNDS),
        nmu_mie=NMU_MIE, nmu_dap=NMU_DAP, nsubr=NSUBR,
    )
    digest = hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()
    return digest[:8]


def _cache_path() -> Path:
    """Path to the single HDF5 file holding all intermediate results."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / f'run_{_params_fingerprint()}.h5'


def _h5_write_atomic(path: Path, writer_fn):
    """Write an HDF5 file atomically: write to .tmp then rename."""
    import h5py
    tmp = path.with_suffix('.tmp')
    with h5py.File(tmp, 'w') as f:
        writer_fn(f)
    tmp.rename(path)   # atomic on POSIX filesystems


def _validate_bmabs(bmabs: np.ndarray, expected_shape: tuple) -> None:
    """Raise if the loaded bmabs looks wrong (wrong shape or all zeros)."""
    if bmabs.shape != expected_shape:
        raise ValueError(
            f'bmabs shape mismatch: got {bmabs.shape}, expected {expected_shape}. '
            'Delete the cache file and rerun.'
        )
    if bmabs.max() == 0.0:
        raise ValueError(
            'bmabs is all zeros — likely cached from a run where HITRAN data was '
            'missing. Delete the cache file and rerun (or use --force).'
        )

# ---------------------------------------------------------------------------
# VIRA profile interpolators
# ---------------------------------------------------------------------------

_VIRA = np.array([
    [  0.0,735.3,92.10],[  4.0,697.4,66.65],[  8.0,660.4,47.35],
    [ 12.0,619.1,33.04],[ 16.0,574.5,22.52],[ 20.0,527.4,14.93],
    [ 24.0,476.0, 9.573],[ 28.0,427.0, 5.917],[ 32.0,380.1, 3.501],
    [ 36.0,337.4, 1.979],[ 40.0,299.7, 1.066],[ 44.0,267.0, 0.5356],
    [ 48.0,238.2, 0.2488],[ 52.0,212.5, 0.1067],[ 56.0,198.8, 4.370e-2],
    [ 60.0,195.2, 1.768e-2],[ 64.0,203.5, 7.132e-3],[ 68.0,210.6, 2.941e-3],
    [ 72.0,215.4, 1.199e-3],[ 76.0,218.2, 4.820e-4],[ 80.0,218.5, 1.920e-4],
    [ 84.0,214.5, 7.526e-5],[ 88.0,206.0, 2.924e-5],[ 92.0,195.5, 1.126e-5],
    [ 96.0,184.0, 4.289e-6],[100.0,172.0, 1.612e-6],
])
_T_fn = interpolate.interp1d(_VIRA[:,0], _VIRA[:,1], kind='linear', fill_value='extrapolate')
_P_fn = interpolate.interp1d(_VIRA[:,0], np.log(_VIRA[:,2]), kind='linear', fill_value='extrapolate')

def vira_T(z): return float(_T_fn(z))
def vira_P(z): return float(np.exp(_P_fn(z)))

# ---------------------------------------------------------------------------
# Dynamic layer construction
#
# LAYER_BOUNDS is derived from the cloud/haze geometry rather than hardcoded.
# Boundaries are automatically placed at the cloud bottom, cloud top (= haze
# bottom), and haze top.  GAS_LAYER_TOPS_KM defines the additional break
# points above the haze; add or remove entries to change vertical resolution.
# ---------------------------------------------------------------------------

# Break points for the gas layers above the haze (in km).
# Edit this list to change vertical resolution above the haze layer.
GAS_LAYER_TOPS_KM = [75., 82., 88., 94., 100.]


def make_layer_bounds() -> np.ndarray:
    """Derive layer boundaries from model geometry.

    Automatically includes the surface (0 km), cloud bottom, cloud top /
    haze bottom, haze top, and all GAS_LAYER_TOPS_KM entries above the haze.
    Deduplication and sorting are handled automatically, so you cannot create
    an invalid boundary set by editing the constants.
    """
    key_alts = sorted(set(
        [0.0, CLOUD_BOT_KM, CLOUD_TOP_KM, HAZE_TOP_KM]
        + [z for z in GAS_LAYER_TOPS_KM if z > HAZE_TOP_KM]
    ))
    return np.array(key_alts)


LAYER_BOUNDS = make_layer_bounds()
N_LAYERS     = len(LAYER_BOUNDS) - 1


def layer_midpoint(z_bot: float, z_top: float):
    """Return (T [K], P [bar]) at the midpoint of layer [z_bot, z_top] km."""
    T = 0.5 * (vira_T(z_bot) + vira_T(z_top))
    P = 0.5 * (vira_P(z_bot) + vira_P(z_top))
    return T, P


def _layer_has_cloud(z_bot: float, z_top: float) -> bool:
    """True when the layer [z_bot, z_top] overlaps with the cloud altitude slab."""
    return z_bot < CLOUD_TOP_KM and z_top > CLOUD_BOT_KM


def _layer_has_haze(z_bot: float, z_top: float) -> bool:
    """True when the layer [z_bot, z_top] overlaps with the haze altitude slab."""
    return z_bot < HAZE_TOP_KM and z_top > HAZE_BOT_KM


def _layer_name(z_bot: float, z_top: float, has_cloud: bool, has_haze: bool) -> str:
    """Auto-generate a descriptive layer name from altitude range and content."""
    if has_cloud and has_haze:
        tag = 'cloud_haze'
    elif has_cloud:
        tag = 'cloud'
    elif has_haze:
        tag = 'haze'
    else:
        tag = 'gas'
    return f'{tag}_{z_bot:.0f}_{z_top:.0f}'


def co2_column_cm2(z_bot_km, z_top_km):
    """CO₂ number column density [molecules cm⁻²] via hydrostatic eq."""
    P_bot = vira_P(z_bot_km) * 1.0e5   # bar → Pa
    P_top = vira_P(z_top_km) * 1.0e5
    N_col = (P_bot - P_top) * VMR_CO2 / (M_CO2_KG * G_VENUS)   # molec m⁻²
    return max(N_col * 1.0e-4, 0.0)   # m⁻² → cm⁻²


def rayleigh_tau(z_bot_km, z_top_km, wav_um):
    """CO₂ Rayleigh optical depth for a layer (formula from Eq. 4 of paper)."""
    n   = N_CO2_N
    lam = wav_um * 1.0e-4   # µm → cm
    rho = RHO_CO2
    fac = (24.0 * np.pi**3 * (n**2 - 1)**2 * (6.0 + 3.0*rho)
           / (lam**4 * (n**2 + 2.0)**2 * (6.0 - 7.0*rho)))
    N_L = 2.6868e19   # Loschmidt number [cm⁻³]
    sigma_ray = fac / N_L**2  # cm²/molecule
    N_col = co2_column_cm2(z_bot_km, z_top_km)
    return float(sigma_ray * N_col)


# ---------------------------------------------------------------------------
# Gauss–Legendre points on [0, 1]
# ---------------------------------------------------------------------------

def gauss_legendre_01(n):
    """Gauss–Legendre nodes and weights on [0, 1]."""
    nodes, weights = np.polynomial.legendre.leggauss(n)
    nodes_01   = 0.5 * (nodes + 1.0)
    weights_01 = 0.5 * weights
    return nodes_01, weights_01


# ---------------------------------------------------------------------------
# Fast CKD: one HAPI call per layer for the full band, then sort per bin
# ---------------------------------------------------------------------------

def _ensure_hitran(cache_dir, wvn_lo, wvn_hi):
    """Download HITRAN CO2 lines if not already cached.  Raises on failure."""
    import hapi as hp, os
    data_file = os.path.join(cache_dir, 'CO2.data')
    os.makedirs(cache_dir, exist_ok=True)
    hp.db_begin(cache_dir)
    if not os.path.exists(data_file) or os.path.getsize(data_file) < 1000:
        print(f'  HITRAN CO2 data missing — downloading from hitran.org …', flush=True)
        hp.fetch('CO2', 2, 1, wvn_lo - 10, wvn_hi + 10)   # main isotopologue
        print(f'  Download complete.', flush=True)
    else:
        hp.db_begin(cache_dir)   # just register the table


def _hapi_lbl_full(molecule, T_K, P_bar, wvn_min, wvn_max, wvn_step, cache_dir):
    """Compute full LBL absorption cross-section via HAPI Voigt."""
    import hapi as hp
    hp.db_begin(cache_dir)
    table = molecule.upper()
    P_atm = P_bar / 1.01325
    nu, sigma = hp.absorptionCoefficient_Voigt(
        SourceTables=table,
        Environment={'T': float(T_K), 'p': float(P_atm)},
        OmegaRange=[float(wvn_min), float(wvn_max)],
        OmegaStep=float(wvn_step),
        OmegaWingHW=25.0,
        GammaL='gamma_air',
        HITRAN_units=True,
    )
    nu_arr = np.asarray(nu, dtype='float64')
    sig_arr = np.asarray(sigma, dtype='float64')
    if sig_arr.max() == 0.0:
        raise RuntimeError(
            f'HAPI returned all-zero cross-sections at T={T_K:.0f}K P={P_bar:.3e}bar. '
            'HITRAN table may be missing or empty.'
        )
    return nu_arr, sig_arr


def compute_bmabs_fast(
    cache_path: Path,
    wav_centers,    # [µm], shape (nwav,)
    gauss_points,   # [0,1], shape (ngauss,)
    bin_width_um = BIN_WIDTH_UM,
    force        = False,
):
    """Compute bmabs[nlayer, nwav, ngauss] via one LBL call per layer.

    Results are stored in the /bmabs group of cache_path (a shared HDF5 file).
    The group attribute 'done_layers' tracks resumability; when it reaches
    N_LAYERS the dataset is considered complete.

    On load, shape and non-zero content are validated — a silent all-zero
    result (e.g. from a missing HITRAN table) raises immediately rather than
    being silently used.
    """
    import h5py

    nwav   = len(wav_centers)
    ngauss = len(gauss_points)
    expected_shape = (N_LAYERS, nwav, ngauss)

    # ── Try loading from cache ────────────────────────────────────────────
    if cache_path.exists() and not force:
        with h5py.File(cache_path, 'r') as f:
            if 'bmabs' in f:
                done_layers = int(f['bmabs'].attrs.get('done_layers', 0))
                if done_layers >= N_LAYERS:
                    print(f'CKD: loading from cache ({cache_path.name})')
                    bmabs = f['bmabs'][...].astype('float64')
                    _validate_bmabs(bmabs, expected_shape)   # raises on zeros/wrong shape
                    return bmabs
                else:
                    # Partial: resume from where we stopped
                    bmabs = f['bmabs'][...].astype('float64')
                    print(f'CKD: resuming from layer {done_layers}/{N_LAYERS}')
            else:
                bmabs = np.zeros(expected_shape, dtype='float64')
                done_layers = 0
    else:
        bmabs = np.zeros(expected_shape, dtype='float64')
        done_layers = 0

    # ── Compute remaining layers ──────────────────────────────────────────
    MARGIN_CM = 5.0
    WVN_HI = 1.0e4 / WAV_MIN_UM + MARGIN_CM
    WVN_LO = 1.0e4 / WAV_MAX_UM - MARGIN_CM
    WVN_STEP = 0.10

    _ensure_hitran(HITRAN_CACHE, WVN_LO, WVN_HI)
    print(f'CKD: LBL [{WVN_LO:.0f}–{WVN_HI:.0f} cm⁻¹  step={WVN_STEP}]')

    for il in range(done_layers, N_LAYERS):
        z_bot = float(LAYER_BOUNDS[il]); z_top = float(LAYER_BOUNDS[il + 1])
        T, P  = layer_midpoint(z_bot, z_top)
        N_col = co2_column_cm2(z_bot, z_top)
        print(f'  Layer {il}: {z_bot:.0f}–{z_top:.0f} km  T={T:.1f}K  '
              f'P={P:.3e}bar  N={N_col:.2e} cm⁻²', flush=True)

        t0 = time.time()
        nu, sigma = _hapi_lbl_full('CO2', T, P, WVN_LO, WVN_HI, WVN_STEP, HITRAN_CACHE)
        print(f'    HAPI: {time.time()-t0:.1f}s  pts={len(nu):,}  '
              f'σmax={sigma.max():.2e}', flush=True)

        tau_lbl = sigma * N_col
        for iw, wav_c in enumerate(wav_centers):
            wvn_c  = 1.0e4 / wav_c
            dnu    = 0.5e4 * bin_width_um / wav_c**2
            mask   = (nu >= wvn_c - dnu) & (nu <= wvn_c + dnu)
            tau_bin = tau_lbl[mask]
            if len(tau_bin) == 0:
                continue
            tau_sorted = np.sort(tau_bin)
            x_cdf = np.linspace(0.0, 1.0, len(tau_sorted))
            for ig, gp in enumerate(gauss_points):
                bmabs[il, iw, ig] = np.interp(gp, x_cdf, tau_sorted)

        # ── Atomic partial save after each layer ─────────────────────────
        tmp = cache_path.with_suffix('.tmp')
        with h5py.File(tmp, 'a') as f:
            if 'bmabs' in f:
                del f['bmabs']
            ds = f.create_dataset('bmabs', data=bmabs.astype('float32'),
                                   compression='gzip')
            ds.attrs['done_layers'] = il + 1
            ds.attrs['shape'] = list(expected_shape)
        tmp.rename(cache_path)
        print(f'    [saved layers 0–{il}]', flush=True)

    print(f'CKD complete: bmabs max={bmabs.max():.3e}')
    return bmabs


# ---------------------------------------------------------------------------
# PyMieDAP model builder
# ---------------------------------------------------------------------------

def _dummy_mixed(nwav):
    """Transparent-layer mixed_aerosols placeholder (gas-only layers)."""
    import pymiedap.pymiedap as pmd
    ma = pmd.Aerosols()
    ma.typ     = 'G'
    ma.coefs   = np.zeros((nwav, 4, 4, 1), dtype='float64')
    ma.ncoefs  = np.ones(nwav, dtype='float64')
    ma.ssalb   = np.zeros(nwav)
    ma.sext    = np.zeros(nwav)
    ma.ssca    = np.zeros(nwav)
    ma.col_dens = 0.0
    return ma


def build_model(wav_list, with_haze=False, layer_bounds=None):
    """Build a PyMieDAP Model for one wavelength (or a small list).

    Which layers receive cloud/haze aerosols is determined purely by altitude
    overlap with the CLOUD_BOT_KM–CLOUD_TOP_KM and HAZE_BOT_KM–HAZE_TOP_KM
    slabs — no hardcoded layer indices.  Layer names are auto-generated from
    altitude ranges so they remain meaningful after any resolution change.

    Parameters
    ----------
    layer_bounds : array-like, optional
        Override LAYER_BOUNDS for this call.  Defaults to the module-level
        LAYER_BOUNDS derived from make_layer_bounds().
    """
    import pymiedap.pymiedap as pmd

    if layer_bounds is None:
        layer_bounds = LAYER_BOUNDS

    nwav     = len(wav_list)
    n_layers = len(layer_bounds) - 1

    if n_layers > 50:
        raise ValueError(
            f'PyMieDAP Fortran limit is 50 layers; got {n_layers}. '
            'Coarsen GAS_LAYER_TOPS_KM.'
        )

    model = pmd.Model()
    model.wvl_list = list(wav_list)
    model.asurf    = 0.0

    # Remove default Layers() attributes
    for nm in list(vars(model.layers).keys()):
        delattr(model.layers, nm)

    for il in range(n_layers):
        z_bot = float(layer_bounds[il])
        z_top = float(layer_bounds[il + 1])
        _, P_mid = layer_midpoint(z_bot, z_top)

        tau_ray_w = [rayleigh_tau(z_bot, z_top, w) for w in wav_list]

        # Determine aerosol content from altitude overlap — no index magic
        is_cloud = _layer_has_cloud(z_bot, z_top)
        is_haze  = with_haze and _layer_has_haze(z_bot, z_top)

        aero_tau = (CLOUD_TAU if is_cloud else 0.0) + (HAZE_TAU if is_haze else 0.0)

        lyr = pmd.Layer(
            tau    = [aero_tau] * nwav,
            tau_g  = [0.0]      * nwav,
            tau_ray= tau_ray_w,
            rayscat= False,
            press  = P_mid,
        )

        if is_cloud:
            lyr.aerosols = pmd.Aerosols(
                nr=[N_R_CLOUD] * nwav, ni=[N_I_CLOUD] * nwav,
                r_eff=CLOUD_R_EFF, v_eff=CLOUD_V_EFF, psd='2', typ='C',
            )
        if is_haze:
            # psd='5' (log-normal r_eff/v_eff): psd='2' (gamma) crashes when
            # v_eff > 1/3 because alpha = 1/v_eff - 3 < 0 → gammln failure.
            lyr.aerosols_haze = pmd.Aerosols(
                nr=[N_R_HAZE] * nwav, ni=[N_I_HAZE] * nwav,
                r_eff=HAZE_R_EFF, v_eff=HAZE_V_EFF, psd='5', typ='H',
            )
        if not is_cloud and not is_haze:
            lyr.mixed_aerosols = _dummy_mixed(nwav)

        setattr(model.layers, _layer_name(z_bot, z_top, is_cloud, is_haze), lyr)

    return model


def prepare_mie(model, wav_list):
    """Run mie_code + mix_aerosols for cloud/haze layers only."""
    import pymiedap.pymiedap as pmd
    for lyr_name, lyr in vars(model.layers).items():
        if hasattr(lyr, 'mixed_aerosols'):
            continue   # already a gas layer with dummy mixed_aerosols
        for aero_name, aero in vars(lyr).items():
            if isinstance(aero, pmd.Aerosols):
                pmd.mie_code(aero, wav_list, ngaur=NMU_MIE, nsubr=NSUBR)
        lyr.mix_aerosols()


# ---------------------------------------------------------------------------
# RT spectrum loop
# ---------------------------------------------------------------------------

def run_spectrum(cache_path: Path, wav_centers, gauss_points, gauss_weights,
                 bmabs, with_haze: bool, tag: str):
    """Compute F and P spectra, resumable per wavelength.

    Progress is stored in the /{tag}/I_sum and /{tag}/Q_sum datasets of
    cache_path, with a 'done_wav' attribute tracking how many bins are done.
    Returns (F_spec, P_spec) when complete, (None, None) when the wall-time
    budget is reached mid-computation (caller should rerun).
    """
    import h5py
    import pymiedap.pymiedap as pmd

    nwav   = len(wav_centers)
    ngauss = len(gauss_points)

    # ── Load any partial progress from cache ─────────────────────────────
    I_sum = np.zeros(nwav); Q_sum = np.zeros(nwav); done_wav = 0
    if cache_path.exists():
        with h5py.File(cache_path, 'r') as f:
            if tag in f and 'I_sum' in f[tag]:
                done_wav = int(f[tag]['I_sum'].attrs.get('done_wav', 0))
                if done_wav > 0:
                    I_sum = f[tag]['I_sum'][...].astype('float64')
                    Q_sum = f[tag]['Q_sum'][...].astype('float64')
                    if done_wav < nwav:
                        print(f'  [{tag}] Resuming from λ {done_wav}/{nwav}')

    if done_wav >= nwav:
        print(f'  [{tag}] All wavelengths already in cache.')
        return np.pi * I_sum, -Q_sum / np.maximum(I_sum, 1e-30)

    # ── RT loop ──────────────────────────────────────────────────────────
    DAP_DIR.mkdir(parents=True, exist_ok=True)
    phase_arr = np.array([PHASE_DEG]); sza_arr = np.array([SZA_DEG])
    emiss_arr = np.array([EMISS_DEG]); az_arr  = np.array([AZ_DEG])
    # beta=0: nadir+az=0 → scattering plane = meridian plane → no rotation.
    # get_cosbeta() wrongly returns 90° for emission=0 (sin(0)=0 in denominator),
    # which flips the sign of Q.  Provide the correct value explicitly.
    beta_arr = np.zeros(1)

    t_start = time.time()
    MAX_WALL = 38.0

    for iw in range(done_wav, nwav):
        wav = wav_centers[iw]
        print(f'  [{tag}] λ={wav:.4f} µm  ({iw+1}/{nwav})', end='  ', flush=True)
        t0 = time.time()

        model = build_model([wav], with_haze=with_haze)
        prepare_mie(model, [wav])
        model.geom = pmd.Geom(sza=[SZA_DEG], emission=[EMISS_DEG], azimuth=[AZ_DEG])
        model.geom.phase = [PHASE_DEG]

        I_band = 0.0; Q_band = 0.0
        for ig in range(ngauss):
            for il, (_, lyr) in enumerate(vars(model.layers).items()):
                lyr.tau_g = np.array([float(bmabs[il, iw, ig])])

            pmd.dap_code(model, output_name=f'f8_{tag}_w{iw:03d}_g{ig:02d}',
                         path_output=str(DAP_DIR) + '/', nmug=NMU_DAP)
            fou_file = str(DAP_DIR / f'fou_{wav:.7f}.dat')
            try:
                I0, Q0, _, _ = pmd.read_dap_output(
                    phase_arr, sza_arr, emiss_arr, fou_file,
                    phi=az_arr, beta=beta_arr)
                I_band += float(gauss_weights[ig]) * float(I0[0])
                Q_band += float(gauss_weights[ig]) * float(Q0[0])
            except Exception as exc:
                print(f'\n    read_dap_output failed: {exc}')

        I_sum[iw] = I_band;  Q_sum[iw] = Q_band
        print(f'done ({time.time()-t0:.1f}s)', flush=True)

        # ── Atomic partial save after each wavelength ─────────────────────
        tmp = cache_path.with_suffix('.tmp')
        import shutil
        if cache_path.exists():
            shutil.copy2(cache_path, tmp)   # preserve other groups
        with h5py.File(tmp, 'a') as f:
            grp = f.require_group(tag)
            for key, arr in [('I_sum', I_sum), ('Q_sum', Q_sum)]:
                if key in grp:
                    del grp[key]
                ds = grp.create_dataset(key, data=arr.astype('float64'),
                                        compression='gzip')
                ds.attrs['done_wav'] = iw + 1
        tmp.rename(cache_path)

        if time.time() - t_start > MAX_WALL:
            print(f'  [{tag}] Wall-time budget reached ({iw+1}/{nwav}). '
                  'Run again to continue.', flush=True)
            return None, None

    return np.pi * I_sum, -Q_sum / np.maximum(I_sum, 1e-30)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def make_plot(wav_centers, F_cloud, P_cloud, F_haze, P_haze, out_path):
    wav_nm = wav_centers * 1000.0

    fig, (ax_F, ax_P) = plt.subplots(1, 2, figsize=(10, 5))

    ls_cloud = dict(color='black', lw=1.6, label='Only cloud')
    ls_haze  = dict(color='black', lw=1.6, ls='--', label='Cloud + haze')

    ax_F.plot(wav_nm, F_cloud, **ls_cloud)
    ax_F.plot(wav_nm, F_haze,  **ls_haze)
    ax_F.set_xlabel('Wavelength (nm)', fontsize=12)
    ax_F.set_ylabel('F', fontsize=12)
    ax_F.set_xlim(1400, 1500)
    ax_F.set_ylim(bottom=0)
    ax_F.legend(frameon=False, fontsize=10)
    ax_F.tick_params(direction='in', top=True, right=True)
    ax_F.set_title('(a) Flux F', fontsize=11)

    ax_P.plot(wav_nm, 100.0*P_cloud, **ls_cloud)
    ax_P.plot(wav_nm, 100.0*P_haze,  **ls_haze)
    ax_P.axhline(0, color='grey', lw=0.7, ls=':')
    ax_P.set_xlabel('Wavelength (nm)', fontsize=12)
    ax_P.set_ylabel('P (%)', fontsize=12)
    ax_P.set_xlim(1400, 1500)
    ax_P.legend(frameon=False, fontsize=10)
    ax_P.tick_params(direction='in', top=True, right=True)
    ax_P.set_title('(b) Degree of polarization P', fontsize=11)

    fig.suptitle(
        r'Venus CO$_2$ 1.4–1.5 µm band  ·  SZA=30°  ·  Nadir  ·  1 nm  ·  '
        'Mahapatra, Rossi & Stam (2024) Fig. 8 top row',
        fontsize=10,
    )
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved → {out_path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(force: bool = False) -> int:
    """
    Main entry point.  Pass --force (or call with force=True) to ignore all
    cached results and recompute everything from scratch.

    Cache design
    ────────────
    All intermediate results live in a *single* HDF5 file:
        .cache/fig8/run_{fingerprint}.h5
    whose name encodes every model parameter that affects the result.  If you
    change SZA, particle sizes, or anything else, the fingerprint changes and
    the old cache is simply never found — no manual cleanup needed.

    File layout
    ───────────
        /bmabs           float32 (8, 100, 10)   attrs: done_layers
        /cloud/I_sum     float64 (100,)          attrs: done_wav
        /cloud/Q_sum     float64 (100,)          attrs: done_wav
        /haze/I_sum      float64 (100,)          attrs: done_wav
        /haze/Q_sum      float64 (100,)          attrs: done_wav

    Validation
    ──────────
    After loading bmabs the code asserts shape and max>0, so a silent
    all-zero result (e.g. from a missing HITRAN table) raises loudly.
    """
    import h5py

    print('=== Mahapatra et al. (2024) – Figure 8 top-row ===')
    print(f'Cloud: r_eff={CLOUD_R_EFF:.3f} µm  v_eff={CLOUD_V_EFF:.4f}  τ={CLOUD_TAU}')
    print(f'Haze:  r_eff={HAZE_R_EFF:.3f} µm  v_eff={HAZE_V_EFF:.4f}  τ={HAZE_TAU}')

    cache = _cache_path()
    fp    = _params_fingerprint()
    print(f'Cache: {cache.name}  (fingerprint={fp})')
    if force and cache.exists():
        cache.unlink()
        print('  --force: deleted old cache.')
    print()

    n_bins = int(round((WAV_MAX_UM - WAV_MIN_UM) / BIN_WIDTH_UM))
    wav    = WAV_MIN_UM + BIN_WIDTH_UM * (np.arange(n_bins) + 0.5)
    gp, gw = gauss_legendre_01(N_GAUSS)

    # ── Step 1: CKD ──────────────────────────────────────────────────────
    print('[1/3] Band-mean absorption (CKD)')
    t0    = time.time()
    bmabs = compute_bmabs_fast(cache, wav, gp, force=force)
    # CKD incomplete if the last layer hasn't been written yet
    if cache.exists():
        with h5py.File(cache, 'r') as f:
            done_layers = int(f['bmabs'].attrs.get('done_layers', 0)) if 'bmabs' in f else 0
    else:
        done_layers = 0
    if done_layers < N_LAYERS:
        print(f'  CKD incomplete ({done_layers}/{N_LAYERS} layers). Re-run to continue.')
        return 1
    print(f'  bmabs {bmabs.shape}  max={bmabs.max():.3e}  ({time.time()-t0:.0f}s)')

    # ── Step 2: Cloud RT ─────────────────────────────────────────────────
    print('\n[2/3] RT spectrum: Only cloud')
    F_cloud, P_cloud = run_spectrum(cache, wav, gp, gw, bmabs,
                                     with_haze=False, tag='cloud')
    if F_cloud is None:
        print('  Cloud RT incomplete. Re-run to continue.')
        return 1

    # ── Step 3: Haze RT ──────────────────────────────────────────────────
    print('\n[3/3] RT spectrum: Cloud + haze')
    F_haze, P_haze = run_spectrum(cache, wav, gp, gw, bmabs,
                                   with_haze=True, tag='haze')
    if F_haze is None:
        print('  Haze RT incomplete. Re-run to continue.')
        return 1

    # ── Plot ─────────────────────────────────────────────────────────────
    make_plot(wav, F_cloud, P_cloud, F_haze, P_haze, OUT_PNG)
    print('=== DONE ===')
    return 0


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='Reproduce Fig. 8 of Mahapatra+2024')
    ap.add_argument('--force', action='store_true',
                    help='Delete cached results and recompute from scratch')
    args = ap.parse_args()
    raise SystemExit(main(force=args.force))
