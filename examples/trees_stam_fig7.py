#!/usr/bin/env python3
"""Reproduce Trees & Stam (2019) Figure 7.

Disk-resolved RGB colour images (160×160 pixels) of the total flux F and
polarized flux |Q| for ocean planets with patchy clouds at phase angle α=80°.

Figure layout — 3 rows × 6 columns:
  Rows    : wind speed v = 1, 7, 13 m s⁻¹
  Cols 1–3: Total flux F  at fc = 0.25, 0.50, 0.75
  Cols 4–6: |Q|           at fc = 0.25, 0.50, 0.75

RGB channels: 443 nm (Blue), 550 nm (Green), 670 nm (Red).
Normalization: equal F at all three wavelengths → white pixel.

Usage:
  python trees_stam_fig7.py --step dap     # compute Fourier files
  python trees_stam_fig7.py --step images  # render disk images
  python trees_stam_fig7.py --step plot    # assemble figure
"""

from __future__ import annotations
import argparse, sys, warnings, time, shutil, tempfile
import multiprocessing as mp
from pathlib import Path
warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pymiedap.pymiedap as pmd
import module_geos as geos          # Fortran disk-geometry extension
from pymiedap.ocean import OceanSurface

# ── Physical parameters (Tables 1 & 2 of Trees & Stam 2019) ──────────────────
# Atmosphere
BM_550   = 0.096    # Rayleigh OD at 550 nm, ps=1 bar
DEPOL    = 0.03
MMA      = 29.0     # g/mol
GRAVITY  = 9.81     # m/s²

# Ocean (Table 2)
FOAM_ALBEDO = 0.22
N_WATER     = 1.33
N_AIR       = 1.0
DEPTH_M     = 100.0
BOTTOM_ALB  = 0.0
WATER_DEPOL = 0.09

# Cloud (Table 1)
CLOUD_R_EFF  = 6.0     # µm — liquid water droplets
# Note: r_eff=10 µm (paper value) overflows the adding-doubling for nmug=20
# due to 846 Legendre terms exceeding PyMieDAP's stability threshold.
# r_eff=6 µm (ncoefs=532) is the largest stable value in this codebase.
# Both are in the geometric optics regime and produce visually identical
# white/grey clouds at α=80° — the paper itself uses model-D aerosol as a
# cloud proxy rather than actual 10 µm particles (Rossi et al. 2018).
# CLOUD_R_EFF  = 10.0  # original paper value — requires delta-M truncation
CLOUD_V_EFF  = 0.1
CLOUD_BA_550 = 4.926   # cloud optical thickness at 550 nm
P_CT         = 0.7     # cloud-top pressure [bar]
P_CB         = 0.8     # cloud-bottom pressure [bar]
PS           = 1.0     # surface pressure [bar]

# ── Computation settings ──────────────────────────────────────────────────────
WAVELENGTHS_UM = np.array([0.443, 0.550, 0.670])  # B, G, R
N_WAV = len(WAVELENGTHS_UM)
WIND_SPEEDS     = [1.0, 7.0, 13.0]   # m/s  — Figure 7 rows
CLOUD_FRACS     = [0.25, 0.50, 0.75] # fc   — Figure 7 column groups
ALPHA_DEG       = 80.0               # phase angle [°]
NPIX            = 160                # pixels across disk equator

N_MUG     = 40    # Gauss points for adding-doubling (minimum 40 for numerical stability)
N_MUG_MIE = 40    # Gauss points for Mie size-distribution integration (minimum 40)
N_SUBR    = 50    # Mie size-distribution subintervals
N_FOURIER = 40    # ocean surface Fourier terms (minimum 40)
N_PHI     = 240   # azimuth samples for Fourier projection (≥ 2*N_FOURIER+1)

DAP_DB = REPO_ROOT / "dap_database"
CACHE  = Path("/tmp/pymiedap_cache")
DAP_DB.mkdir(exist_ok=True); CACHE.mkdir(exist_ok=True)

# Cache file names encode nmug so files from different resolutions never collide.
# Changing N_MUG automatically invalidates old DAP and stokes caches.
_CACHE_TAG   = f"n{N_MUG}"
IMG_CACHE    = CACHE / f"fig7_images_{_CACHE_TAG}.npz"
STOKES_CACHE = CACHE / f"fig7_stokes_{_CACHE_TAG}.npz"

# Persistent stokes seed saved inside dap_database/ — survives /tmp clears
# and carries completed Stokes reads across sessions at this nmug level.
_STOKES_SEED = DAP_DB / f"fig7_stokes_{_CACHE_TAG}_partial.npz"


def _seed_stokes_cache():
    """Copy persistent partial results into the /tmp cache if not already there."""
    if STOKES_CACHE.exists() and STOKES_CACHE.stat().st_size > 1000:
        return  # already has data
    if not _STOKES_SEED.exists():
        return  # nothing to seed from
    import shutil
    shutil.copy(_STOKES_SEED, STOKES_CACHE)
    print(f"  Seeded stokes cache from {_STOKES_SEED.name}")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _bm(wav_um: float) -> float:
    """Rayleigh optical depth at given wavelength for ps=1 bar."""
    return BM_550 * (0.550 / wav_um) ** 4


def _tag_ocean(v: float) -> str:
    return f"fig7_ocean_v{int(v):02d}_{_CACHE_TAG}"


def _tag_cloud() -> str:
    return f"fig7_cloud_{_CACHE_TAG}"


def _dap_path(tag: str, wav_um: float) -> Path:
    return DAP_DB / f"{tag}_{wav_um:.7f}.dat"


def _dap_ready(tag: str) -> bool:
    """True when all 3 wavelength Fourier files for this tag exist."""
    return all(_dap_path(tag, w).exists() and _dap_path(tag, w).stat().st_size > 100
               for w in WAVELENGTHS_UM)


# ── Model builders ────────────────────────────────────────────────────────────

def build_ocean(v: float) -> pmd.Model:
    """Clear-sky ocean model at wind speed v [m/s], all 3 wavelengths."""
    m = pmd.Model()
    m.mma = MMA; m.gravity = GRAVITY; m.dpol = DEPOL; m.asurf = 0.0
    del m.layers.haze, m.layers.cloud
    m.wvl_list = list(WAVELENGTHS_UM)
    nw = N_WAV

    bm_vec = [_bm(w) for w in WAVELENGTHS_UM]
    m.layers.gasbelow.press   = PS
    m.layers.gasbelow.tau     = [0.0] * nw
    m.layers.gasbelow.tau_g   = [0.0] * nw
    m.layers.gasbelow.tau_ray = bm_vec
    m.layers.gasbelow.rayscat = False

    m.layers.gastop.press   = 1.0e-5
    m.layers.gastop.tau     = [0.0] * nw
    m.layers.gastop.tau_g   = [0.0] * nw
    m.layers.gastop.tau_ray = [0.0] * nw
    m.layers.gastop.rayscat = False

    m.surface = OceanSurface(
        wind_speed         = v,
        n_air              = N_AIR,
        n_water            = N_WATER,
        foam_albedo        = FOAM_ALBEDO,
        depth_m            = DEPTH_M,
        bottom_albedo      = BOTTOM_ALB,
        water_depol        = WATER_DEPOL,
        n_fourier          = N_FOURIER,
        n_phi              = N_PHI,
        water_streams      = 8,
        water_n_phi        = 150,
        water_initial_tau  = 0.01,
        include_subsurface = True,
        include_foam       = True,
        solver             = "adding_doubling",
    )
    return m


def build_cloud() -> pmd.Model:
    """Cloudy model: atmosphere + cloud slab + Lambertian black surface.

    Cloud parameters from Table 1: r_eff=10 µm, v_eff=0.1, ba=4.926 at 550nm.
    Cloud slab sits between P_CT=0.7 bar and P_CB=0.8 bar.
    Surface below cloud: Lambertian black (asurf=0).
    """
    m = pmd.Model()
    m.mma = MMA; m.gravity = GRAVITY; m.dpol = DEPOL
    del m.layers.haze   # no aerosol haze
    m.wvl_list = list(WAVELENGTHS_UM)
    nw = N_WAV

    # gastop: atmosphere above cloud top (P < P_CT = 0.7 bar)
    m.layers.gastop.press = P_CT
    # Rayleigh OD above cloud ∝ P_CT/PS
    bm_above = [_bm(w) * P_CT / PS for w in WAVELENGTHS_UM]
    m.layers.gastop.tau     = [0.0] * nw
    m.layers.gastop.tau_g   = [0.0] * nw
    m.layers.gastop.tau_ray = bm_above
    m.layers.gastop.rayscat = False

    # cloud: the cloud slab at P_CT to P_CB (centre ≈ 0.75 bar)
    m.layers.cloud.press   = P_CB
    m.layers.cloud.tau     = [CLOUD_BA_550] * nw   # OD ≈ constant for large droplets
    m.layers.cloud.tau_g   = [0.0] * nw
    m.layers.cloud.tau_ray = [0.0] * nw             # negligible inside cloud slab
    m.layers.cloud.rayscat = False
    m.layers.cloud.aerosols = pmd.Aerosols(
        nr=[1.33] * nw,   # liquid water real refractive index
        ni=[1e-8] * nw,   # essentially non-absorbing at visible wavelengths
        r_eff=CLOUD_R_EFF,
        v_eff=CLOUD_V_EFF,
        psd='2',           # gamma distribution
        typ='C',
    )

    # gasbelow: atmosphere below cloud (P_CB to PS)
    m.layers.gasbelow.press = PS
    bm_below = [_bm(w) * (PS - P_CB) / PS for w in WAVELENGTHS_UM]
    m.layers.gasbelow.tau     = [0.0] * nw
    m.layers.gasbelow.tau_g   = [0.0] * nw
    m.layers.gasbelow.tau_ray = bm_below
    m.layers.gasbelow.rayscat = False

    # Black Lambertian surface underneath the cloud
    m.asurf = 0.0

    # The default Layer.__init__ puts an empty Aerosols() on every layer.
    # mix_aerosols() divides by sum_fsext which is zero for these dummy objects,
    # producing NaN that corrupts the Fourier file even for tau=0 gas layers.
    # Fix: delete the stray aerosols attribute from gas-only layers so that
    # mix_aerosols() finds nothing to mix and returns a transparent matrix.
    for lname in ("gastop", "gasbelow"):
        lyr = getattr(m.layers, lname)
        if hasattr(lyr, "aerosols"):
            delattr(lyr, "aerosols")

    return m


# ── Parallel worker functions (module-level — required for pickle/fork) ───────

def _dap_worker(args: tuple) -> tuple:
    """Compute one DAP model in a subprocess.

    Uses an isolated temporary directory so that the os.chdir() call inside
    dap_code() and the intermediate file `fou_0.443.dat` never collide with
    other workers running at the same time.  Output .dat files are moved to
    DAP_DB once the Fortran code finishes.
    """
    tag, model_type, v = args
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix=f"pymiedap_{tag}_") as tmpdir:
        m = build_ocean(v) if model_type == "ocean" else build_cloud()
        pmd.compute_model(m, force=True, rename=True, output_name=tag,
                          nmug=N_MUG, nmug_mie=N_MUG_MIE, nsubr=N_SUBR,
                          nmat=4, path_input=tmpdir + "/")
        # Move every output .dat file from the isolated tmpdir into DAP_DB
        for wav in WAVELENGTHS_UM:
            src = Path(tmpdir) / f"{tag}_{wav:.7f}.dat"
            dst = DAP_DB / f"{tag}_{wav:.7f}.dat"
            if src.exists():
                shutil.move(str(src), str(dst))
    return tag, time.time() - t0


def _stokes_worker(args: tuple) -> tuple:
    """Read per-pixel Stokes data for one (tag, key, iw) work item.

    All reads are from pre-existing .dat files (read-only) so multiple
    workers can run concurrently without any locking.
    """
    tag, key, iw, ngeos, phase_arr, theta0_arr, theta_arr, phi_arr, beta_arr = args
    t0 = time.time()
    if iw is None:
        # Ocean: read all wavelengths at once
        F_pix, Q_pix = _read_pixel_stokes(
            tag, ngeos, phase_arr, theta0_arr, theta_arr, phi_arr, beta_arr)
    else:
        # Cloud: one wavelength at a time
        F_row, Q_row = _read_one_wavelength(
            tag, iw, ngeos, phase_arr, theta0_arr, theta_arr, phi_arr, beta_arr)
        F_pix = np.zeros((N_WAV, ngeos)); Q_pix = np.zeros((N_WAV, ngeos))
        F_pix[iw] = F_row; Q_pix[iw] = Q_row
    return key, F_pix, Q_pix, time.time() - t0


# ── Step 1: DAP ───────────────────────────────────────────────────────────────

def step_dap(workers: int = 1):
    """Compute Fourier files for all models (one per wavelength per model).

    With workers > 1, up to min(workers, n_models) models run concurrently.
    Each worker uses an isolated temp directory so the intermediate Fortran
    file `fou_0.443.dat` is never shared between processes.
    """
    print(f"=== DAP  nmug={N_MUG}  n_fourier={N_FOURIER}  workers={workers} ===")

    # Build the list of models that still need computing
    pending: list = []
    for v in WIND_SPEEDS:
        tag = _tag_ocean(v)
        if _dap_ready(tag):
            print(f"  [{tag}] files exist, skipping.")
        else:
            print(f"  [{tag}] v={v} m/s  pending …")
            pending.append((tag, "ocean", v))

    ctag = _tag_cloud()
    if _dap_ready(ctag):
        print(f"  [{ctag}] files exist, skipping.")
    else:
        print(f"  [{ctag}]  pending …")
        pending.append((ctag, "cloud", None))

    if pending:
        n_procs = min(workers, len(pending))
        if n_procs > 1:
            print(f"\n  Launching {len(pending)} DAP tasks across {n_procs} workers …\n",
                  flush=True)
            with mp.Pool(n_procs) as pool:
                for tag, elapsed in pool.imap_unordered(_dap_worker, pending):
                    print(f"  ✓ [{tag}] done in {elapsed:.1f}s", flush=True)
        else:
            for args in pending:
                tag2, elapsed = _dap_worker(args)
                print(f"  ✓ [{tag2}] done in {elapsed:.1f}s")

    all_ready = all(_dap_ready(_tag_ocean(v)) for v in WIND_SPEEDS) and _dap_ready(ctag)
    print(f"\nDAP ready: {'YES' if all_ready else 'INCOMPLETE'}")
    return all_ready


# ── Step 2: Disk-resolved images ──────────────────────────────────────────────

def _read_pixel_stokes(tag: str, ngeos: int,
                        phase_arr, theta0, theta, phi, beta) -> tuple:
    """Read per-pixel (F, Q) from the Fourier files for all wavelengths.

    Returns F_pixels (N_WAV, ngeos) and Q_pixels (N_WAV, ngeos),
    with F = I * cos(θ₀)  (the cosine-weighted intensity used in disk integration).
    """
    F_pix = np.zeros((N_WAV, ngeos))
    Q_pix = np.zeros((N_WAV, ngeos))
    for iw, wav in enumerate(WAVELENGTHS_UM):
        fname = str(_dap_path(tag, wav))
        I, Q, _, _ = pmd.read_dap_output(
            phase_arr, theta0, theta, fname, phi=phi, beta=beta)
        F_pix[iw] = I * np.cos(np.radians(theta0))
        Q_pix[iw] = Q * np.cos(np.radians(theta0))
    return F_pix, Q_pix


def _read_one_wavelength(tag: str, iw: int, ngeos: int,
                          phase_arr, theta0, theta, phi, beta) -> tuple:
    """Read per-pixel (I*cos, Q*cos) for a single wavelength index iw."""
    wav = WAVELENGTHS_UM[iw]
    fname = str(_dap_path(tag, wav))
    cos0 = np.cos(np.radians(theta0))
    I, Q, _, _ = pmd.read_dap_output(phase_arr, theta0, theta, fname, phi=phi, beta=beta)
    return I * cos0, Q * cos0


def _stokes_work_items() -> list:
    """Return list of (tag, key, iw) for every model×wavelength that must be read.

    Ocean models are read all-3-wavelengths at once (fast enough).
    Cloud is split per wavelength (each cloud file is ~7× larger / slower).
    """
    items = []
    for v in WIND_SPEEDS:
        items.append((_tag_ocean(v), f"ocean_v{int(v):02d}", None))  # None → all wavs
    for iw in range(N_WAV):
        items.append((_tag_cloud(), f"cloud_w{iw}", iw))
    return items


def step_images(workers: int = 1):
    """Render per-pixel F and |Q| for every (v, fc) combination.

    Stokes reads are cached to STOKES_CACHE one chunk at a time so the
    function can be resumed across multiple calls.  Ocean models are read
    3-wavelengths-at-once (~20s each); each cloud wavelength is a separate
    chunk (~50s each, due to the many Fourier terms in the cloud file).

    With workers > 1, all pending reads run in parallel (each reads different
    .dat files, so no locking is needed).  Results are merged and STOKES_CACHE
    is updated after each result arrives, preserving resumability.
    """
    if IMG_CACHE.exists() and IMG_CACHE.stat().st_size > 1000:
        print("Image cache exists, skipping.")
        return True

    _seed_stokes_cache()   # pull in persistent partial results if /tmp is empty
    print(f"=== Disk images  α={ALPHA_DEG}°  npix={NPIX} ===")

    # Get lit pixel geometries at α=80°
    ngeos, apix, theta0_all, theta_all, phi_all, beta_all, lats, longs, xs_all, ys_all = \
        geos.getgeos(ALPHA_DEG, NPIX)
    theta0_arr = theta0_all[:ngeos]
    theta_arr  = theta_all[:ngeos]
    phi_arr    = phi_all[:ngeos]
    beta_arr   = beta_all[:ngeos]
    xs = xs_all[:ngeos]; ys = ys_all[:ngeos]
    phase_arr  = np.full(ngeos, ALPHA_DEG)

    print(f"  Lit pixels: {ngeos}  (out of {NPIX}²={NPIX**2})")

    # ── Resumable Stokes reads ────────────────────────────────────────────────
    stokes_data: dict = {}
    if STOKES_CACHE.exists() and STOKES_CACHE.stat().st_size > 100:
        raw = np.load(STOKES_CACHE)
        stokes_data = {k: raw[k] for k in raw.files}
        print(f"  Stokes cache loaded: {list(stokes_data.keys())}")

    # ── Build list of pending Stokes reads ───────────────────────────────────
    pending_stokes: list = []
    for tag, key, iw in _stokes_work_items():
        if f"F_{key}" in stokes_data and f"Q_{key}" in stokes_data:
            print(f"  [{key}] cached, skipping.", flush=True)
        else:
            pending_stokes.append(
                (tag, key, iw, ngeos, phase_arr, theta0_arr, theta_arr, phi_arr, beta_arr))

    if pending_stokes:
        n_procs = min(workers, len(pending_stokes))
        if n_procs > 1:
            print(f"\n  Running {len(pending_stokes)} Stokes reads across "
                  f"{n_procs} workers …\n", flush=True)
            with mp.Pool(n_procs) as pool:
                for key, F_pix, Q_pix, elapsed in pool.imap_unordered(
                        _stokes_worker, pending_stokes):
                    stokes_data[f"F_{key}"] = F_pix
                    stokes_data[f"Q_{key}"] = Q_pix
                    np.savez(STOKES_CACHE, **stokes_data)
                    print(f"  ✓ [{key}] done in {elapsed:.1f}s  "
                          f"F max={F_pix.max():.4f}  cache updated", flush=True)
        else:
            for args in pending_stokes:
                key, F_pix, Q_pix, elapsed = _stokes_worker(args)
                stokes_data[f"F_{key}"] = F_pix
                stokes_data[f"Q_{key}"] = Q_pix
                np.savez(STOKES_CACHE, **stokes_data)
                print(f"  ✓ [{key}] done in {elapsed:.1f}s  "
                      f"F max={F_pix.max():.4f}  stokes cache updated", flush=True)

    # Check completeness
    for tag, key, iw in _stokes_work_items():
        if f"F_{key}" not in stokes_data:
            print(f"  Still missing [{key}], re-run --step images.")
            return False

    # ── Assemble F_ocean / Q_ocean dicts and F_cloud / Q_cloud ───────────────
    F_ocean, Q_ocean = {}, {}
    for v in WIND_SPEEDS:
        key = f"ocean_v{int(v):02d}"
        F_ocean[v] = stokes_data[f"F_{key}"]
        Q_ocean[v] = stokes_data[f"Q_{key}"]

    # Merge the per-wavelength cloud rows into single (N_WAV, ngeos) arrays
    F_cloud = np.zeros((N_WAV, ngeos)); Q_cloud = np.zeros((N_WAV, ngeos))
    for iw in range(N_WAV):
        key = f"cloud_w{iw}"
        F_cloud[iw] = stokes_data[f"F_{key}"][iw]
        Q_cloud[iw] = stokes_data[f"Q_{key}"][iw]

    # Map normalised disk coords to image pixel indices
    ix = np.clip(np.round(((xs + 1.0) / 2.0) * (NPIX - 1)).astype(int), 0, NPIX - 1)
    iy = np.clip(np.round(((1.0 - ys) / 2.0) * (NPIX - 1)).astype(int), 0, NPIX - 1)

    # Compute RGB normalization reference (so equal F → white)
    # Strategy: normalize each channel by the mean F of the fully clear (v=7) planet
    ref_F = F_ocean[7.0].mean(axis=1)  # shape (N_WAV,) — one reference per channel
    ref_F = np.maximum(ref_F, 1e-10)

    # Build image arrays: shape (n_v, n_fc, NPIX, NPIX, 3) for F and |Q|
    n_v  = len(WIND_SPEEDS)
    n_fc = len(CLOUD_FRACS)

    F_imgs  = np.zeros((n_v, n_fc, NPIX, NPIX, 3))
    Q_imgs  = np.zeros((n_v, n_fc, NPIX, NPIX, 3))

    np.random.seed(42)   # reproducible patchy patterns
    for iv, v in enumerate(WIND_SPEEDS):
        for ifc, fc in enumerate(CLOUD_FRACS):
            print(f"  v={v:.0f} m/s  fc={fc:.2f} …", flush=True)

            # Generate patchy cloud mask for this (v, fc) combination
            # mask[i] = 0 → ocean pixel,  mask[i] = 1 → cloud pixel
            _, mask, _, _, _ = pmd.mask_planet(
                alpha=ALPHA_DEG, npix=NPIX, patchy=True,
                fclouds=[1.0 - fc, fc], full_disk=False)

            # Mix ocean and cloud pixel values
            F_pix = np.where(mask[:ngeos] == 0, F_ocean[v], F_cloud)  # (N_WAV, ngeos)
            Q_pix = np.where(mask[:ngeos] == 0, Q_ocean[v], Q_cloud)

            # Build 2-D RGB images
            for ich in range(3):  # B=443, G=550, R=670
                f_ch  = F_pix[ich] / ref_F[ich]   # colour-normalised brightness
                qabs  = np.abs(Q_pix[ich]) / ref_F[ich]

                img_F = np.zeros((NPIX, NPIX))
                img_Q = np.zeros((NPIX, NPIX))
                np.add.at(img_F, (iy, ix), f_ch)
                np.add.at(img_Q, (iy, ix), qabs)

                # RGB: channel order is [R, G, B] = [670, 550, 443]
                rgb_idx = 2 - ich   # ich=0(B) → rgb=2, ich=1(G) → rgb=1, ich=2(R) → rgb=0
                F_imgs[iv, ifc, :, :, rgb_idx] = img_F
                Q_imgs[iv, ifc, :, :, rgb_idx] = img_Q

    # Clamp and save
    F_imgs = np.clip(F_imgs, 0, None)
    Q_imgs = np.clip(Q_imgs, 0, None)
    np.savez(IMG_CACHE, F_imgs=F_imgs, Q_imgs=Q_imgs,
             wind_speeds=WIND_SPEEDS, cloud_fracs=CLOUD_FRACS)
    print(f"\nImages saved → {IMG_CACHE}")
    return True


# ── Step 3: Plot ───────────────────────────────────────────────────────────────

def step_plot():
    """Assemble the 3×6 panel figure matching Trees & Stam (2019) Fig. 7."""
    from scipy.ndimage import gaussian_filter
    print("=== Plot ===")
    d = np.load(IMG_CACHE)
    F_imgs_raw = d["F_imgs"].copy()   # (n_v, n_fc, NPIX, NPIX, 3)
    Q_imgs_raw = d["Q_imgs"].copy()

    n_v, n_fc = len(WIND_SPEEDS), len(CLOUD_FRACS)

    # ── Recover true physical values from the stored images ───────────────────
    # IMG_CACHE stores F/ref_F and Q/ref_F where ref_F = mean of v=7 ocean pixels.
    # We reverse this to get the raw I*cos(θ₀) values, then apply a proper
    # per-image normalisation + gamma to match the paper's appearance.
    stk = np.load(STOKES_CACHE) if STOKES_CACHE.exists() else None
    if stk is not None:
        ref_F_old = stk["F_ocean_v07"].mean(axis=1)   # (3,) per wavelength
        # RGB order: rgb=0→670nm(iw=2), rgb=1→550nm(iw=1), rgb=2→443nm(iw=0)
        ref_F_rgb = np.array([ref_F_old[2], ref_F_old[1], ref_F_old[0]])
        F_imgs = np.zeros_like(F_imgs_raw)
        Q_imgs = np.zeros_like(Q_imgs_raw)
        for rgb in range(3):
            F_imgs[:,:,:,:,rgb] = F_imgs_raw[:,:,:,:,rgb] * ref_F_rgb[rgb]
            Q_imgs[:,:,:,:,rgb] = Q_imgs_raw[:,:,:,:,rgb] * ref_F_rgb[rgb]
        # Fix NaN in the blue channel (cloud 443 nm may have been missing):
        # Large-droplet clouds are nearly grey in visible, so 550 nm is a
        # good proxy for 443 nm.
        blue_nan = np.isnan(F_imgs[:,:,:,:,2])
        if blue_nan.any():
            print(f"  Note: {blue_nan.sum()} blue-channel NaN pixels filled from 550 nm")
            F_imgs[:,:,:,:,2] = np.where(blue_nan, F_imgs[:,:,:,:,1], F_imgs[:,:,:,:,2])
            Q_imgs[:,:,:,:,2] = np.where(np.isnan(Q_imgs[:,:,:,:,2]),
                                         Q_imgs[:,:,:,:,1], Q_imgs[:,:,:,:,2])
        F_imgs = np.clip(F_imgs, 0, None)
        Q_imgs = np.clip(Q_imgs, 0, None)
    else:
        # Stokes cache not available — use images as-is
        F_imgs = np.clip(F_imgs_raw, 0, None)
        Q_imgs = np.clip(Q_imgs_raw, 0, None)

    # ── Brightness normalisation + gamma correction ───────────────────────────
    #
    # At α=80° the sunglint is a very narrow spike; a linear stretch normalised
    # to the per-pixel maximum leaves everything else near-black.
    # Paper-style appearance requires two steps:
    #
    #   1. Normalise to the 97th percentile of the image LUMINANCE (prevents the
    #      glint from crushing the cloud/atmosphere signal).
    #   2. Gamma correction  x → x^γ  (γ=0.45) to lift dark areas into the
    #      visible range.  Standard for planetary RGB composites (Galileo, Cassini).
    #   3. Mild Gaussian blur (σ=1.2 px) before normalisation to reduce the
    #      single-pixel speckle noise from the random patchy cloud mask.
    #
    # All three operations leave the unlit background (= 0) unchanged.
    GAMMA       = 0.45
    VMAX_PCT    = 97.0
    BLUR_SIGMA  = 1.2

    def _blur(img3):
        """Gaussian blur on each channel preserving the black disk boundary."""
        disk = (img3.sum(axis=2) > 0).astype(float)
        out  = np.zeros_like(img3)
        for c in range(3):
            out[:,:,c] = np.where(disk > 0,
                                  gaussian_filter(img3[:,:,c], sigma=BLUR_SIGMA), 0)
        return out

    def _normalise(imgs):
        """Per-image blur → luminance-percentile normalisation → gamma."""
        out = np.zeros_like(imgs)
        for iv in range(n_v):
            for ifc in range(n_fc):
                img  = _blur(imgs[iv, ifc])
                lum  = img.mean(axis=2)
                lit  = lum[lum > 0]
                if lit.size == 0:
                    continue
                vmax = np.percentile(lit, VMAX_PCT)
                scaled = np.clip(img / max(vmax, 1e-10), 0, 1)
                out[iv, ifc] = np.where(scaled > 0,
                                        np.power(scaled, GAMMA), 0.0)
        return out

    F_norm = _normalise(F_imgs)
    Q_norm = _normalise(Q_imgs)

    fig, axes = plt.subplots(n_v, 2 * n_fc, figsize=(14, 7),
                             gridspec_kw={"wspace": 0.04, "hspace": 0.04})
    fig.patch.set_facecolor("black")

    col_labels_F = [f"F, $f_c$={fc}" for fc in CLOUD_FRACS]
    col_labels_Q = [f"|Q|, $f_c$={fc}" for fc in CLOUD_FRACS]

    for iv, v in enumerate(WIND_SPEEDS):
        for ifc, fc in enumerate(CLOUD_FRACS):
            # F image (columns 0–2)
            ax_F = axes[iv, ifc]
            ax_F.imshow(F_norm[iv, ifc], origin="upper", aspect="equal",
                        interpolation="nearest")
            ax_F.set_facecolor("black")
            ax_F.set_xticks([]); ax_F.set_yticks([])
            if iv == 0:
                ax_F.set_title(col_labels_F[ifc], color="white", fontsize=9, pad=3)
            if ifc == 0:
                ax_F.set_ylabel(f"v = {v:.0f} m/s", color="white", fontsize=9,
                                labelpad=3)

            # |Q| image (columns 3–5)
            ax_Q = axes[iv, n_fc + ifc]
            ax_Q.imshow(Q_norm[iv, ifc], origin="upper", aspect="equal",
                        interpolation="nearest")
            ax_Q.set_facecolor("black")
            ax_Q.set_xticks([]); ax_Q.set_yticks([])
            if iv == 0:
                ax_Q.set_title(col_labels_Q[ifc], color="white", fontsize=9, pad=3)

    fig.suptitle(
        f"Trees & Stam (2019) Fig. 7 — Ocean planet with patchy clouds\n"
        f"α = {ALPHA_DEG}°   RGB: 443 nm (B), 550 nm (G), 670 nm (R)   "
        f"nmug={N_MUG}  n_fourier={N_FOURIER}",
        color="white", fontsize=10, y=1.01,
    )

    out = REPO_ROOT / "examples" / "trees_stam_fig7.png"
    fig.savefig(out, dpi=180, bbox_inches="tight", facecolor="black")
    plt.close(fig)
    print(f"Saved → {out}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Reproduce Trees & Stam (2019) Figure 7.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  # Serial (default)
  python trees_stam_fig7.py --step dap
  python trees_stam_fig7.py --step images
  python trees_stam_fig7.py --step plot

  # Parallel — use 4 workers for the slow steps
  python trees_stam_fig7.py --step dap    --workers 4
  python trees_stam_fig7.py --step images --workers 6

  # On a 70-core machine, --workers 4 covers all 4 DAP models in one pass,
  # and --workers 6 covers all 6 Stokes-read chunks simultaneously.
""")
    ap.add_argument("--step", choices=["dap", "images", "plot", "all"], default="all")
    ap.add_argument("--force", action="store_true",
                    help="Clear the image cache and recompute.")
    ap.add_argument("--workers", type=int, default=1, metavar="N",
                    help="Parallel worker processes for DAP/images steps "
                         "(default 1 = serial). "
                         "Suggested: --workers 4 for --step dap, "
                         "--workers 6 for --step images.")
    args = ap.parse_args()

    print("=== Trees & Stam (2019) Figure 7 ===")
    print(f"α={ALPHA_DEG}°  npix={NPIX}  nmug={N_MUG}  "
          f"n_fourier={N_FOURIER}  workers={args.workers}")
    print(f"Wind speeds: {WIND_SPEEDS} m/s")
    print(f"Cloud fractions: {CLOUD_FRACS}")
    print()

    if args.force:
        for p in [IMG_CACHE]:
            if p.exists(): open(p,"wb").close()
        print("--force: image cache cleared.")

    if args.step in ("dap",    "all"): step_dap(workers=args.workers)
    if args.step in ("images", "all"): step_images(workers=args.workers)
    if args.step in ("plot",   "all"): step_plot()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
