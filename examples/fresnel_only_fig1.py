#!/usr/bin/env python3
"""Reproduce the Fresnel-only dashed curve from Trees & Stam (2019) Fig. 1.

Planet: rough Cox-Munk Fresnel surface (n_air=1.0, n_water=1.33, v=7 m/s)
        with a black surface below.  No sub-surface ocean, no atmosphere above.

This is the simplest case and a clean validation of the polarized ocean-surface
physics: the reflected signal is purely from the wind-roughened air-water interface.

Run sequentially:
  python fresnel_only_fig1.py --step dap     # compute Fourier files (once)
  python fresnel_only_fig1.py --step integ   # disk integration (batches of 10 phase angles)
  python fresnel_only_fig1.py --step plot    # assemble figure

Or all at once (slow):
  python fresnel_only_fig1.py --step all
"""

from __future__ import annotations
import argparse, sys, warnings, time
from pathlib import Path
warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pymiedap.pymiedap as pmd
from pymiedap.ocean import OceanSurface

# ── Parameters ─────────────────────────────────────────────────────────────────
# The Fresnel surface uses a fixed n_water=1.33 with no wavelength dependence,
# so all wavelengths give identical F, Q, P.  We compute ONE representative
# wavelength and tile the result to all five for the final plot — this cuts
# the runtime by 5× (from ~20 min to ~4 min at nmug=50, npix=40).
WAVELENGTHS_UM = np.array([0.350, 0.443, 0.550, 0.670, 0.865])
WAVELENGTHS_NM = WAVELENGTHS_UM * 1000
N_WAV = len(WAVELENGTHS_UM)

# Single representative wavelength used for computation
COMPUTE_WAV_UM = np.array([0.550])   # one wavelength — result applies to all

# Trees & Stam Table 2 — ocean surface
WIND_SPEED  = 7.0     # m/s
N_WATER     = 1.33
N_AIR       = 1.0
FOAM_ALBEDO = 0.22

# Computation settings
N_MUG      = 120     # Gauss points for adding-doubling (atmosphere + surface)

# ── Paths ──────────────────────────────────────────────────────────────────────
DAP_DB   = REPO_ROOT / "dap_database"
CACHE    = Path("/tmp/pymiedap_cache")
OUT_PNG  = REPO_ROOT / "examples" / f"fresnel_only_fig1_nmug{N_MUG}.png"
OUT_NPZ  = CACHE / f"fresnel_only_nmug{N_MUG}.npz"
PART_NPZ = CACHE / f"fresnel_only_nmug{N_MUG}.partial.npz"
DAP_DB.mkdir(exist_ok=True); CACHE.mkdir(exist_ok=True)
N_MUG_MIE  = 16      # Gauss points for Mie (not used here — no aerosols)
N_SUBR     = 50      # size-distribution subdivisions (unused for gas-only)
# n_fourier=120 causes numerical cancellation at small phase angles where F≈0
# (Gibbs-like noise from 121 terms summing to near-zero). The glint is fully
# converged at n_fourier=80 for v=7 m/s; higher terms add only noise.
N_FOURIER  = 80      # ocean surface Fourier terms (converged; 120 gives noise)
N_PHI      = 360     # azimuth samples  (≥ 2*N_FOURIER+1 = 161)

ALPHA_DEG  = np.arange(0., 181., 3.)   # 61 phase angles 0-180°
N_PIX      = 40      # pixels across the equator (paper: ≥ 40)

WAV_COLORS = {
    0.350: "#BF5FBF",   # violet
    0.443: "#2980B9",   # blue
    0.550: "#27AE60",   # green
    0.670: "#E74C3C",   # red
    0.865: "#7F8C8D",   # grey
}
WAV_LABELS = {w: f"{int(w*1000)} nm" for w in WAVELENGTHS_UM}


# ── Model builder ──────────────────────────────────────────────────────────────

def build_model(wav_list=None) -> pmd.Model:
    """Fresnel surface + black bottom, negligible atmosphere."""
    if wav_list is None:
        wav_list = list(WAVELENGTHS_UM)

    m = pmd.Model()
    m.mma = 29; m.gravity = 9.81; m.dpol = 0.03; m.asurf = 0.0
    del m.layers.haze, m.layers.cloud
    m.wvl_list = list(wav_list)
    nw = len(wav_list)

    # Negligible atmosphere: set Rayleigh OD = 0, very small pressure
    for lname in ("gasbelow", "gastop"):
        lyr = getattr(m.layers, lname)
        lyr.press   = 1e-8 if lname == "gasbelow" else 1e-10
        lyr.tau     = [0.0] * nw
        lyr.tau_g   = [0.0] * nw
        lyr.tau_ray = [0.0] * nw
        lyr.rayscat = False

    m.surface = OceanSurface(
        wind_speed         = WIND_SPEED,
        n_air              = N_AIR,
        n_water            = N_WATER,
        foam_albedo        = FOAM_ALBEDO,
        n_fourier          = N_FOURIER,
        n_phi              = N_PHI,
        include_subsurface = False,   # pure Fresnel reflection only
        include_foam       = True,
        solver             = "adding_doubling",
    )
    return m


# ── Step 1: pre-compute one DAP Fourier file per wavelength ───────────────────

def _dap_name(iw: int, wav: float) -> str:
    return f"fonly{N_MUG}_w{iw}"

def _dap_file(iw: int, wav: float) -> Path:
    return DAP_DB / f"{_dap_name(iw, wav)}_{wav:.7f}.dat"


def step_dap():
    """Compute the Fourier coefficient file for the single representative wavelength."""
    iw, wav = 0, float(COMPUTE_WAV_UM[0])
    fpath = _dap_file(iw, wav)
    if fpath.exists() and fpath.stat().st_size > 100:
        print(f"  λ={wav:.3f} µm: Fourier file exists, skipping.")
    else:
        print(f"  λ={wav:.3f} µm: computing (nmug={N_MUG}, n_fourier={N_FOURIER}) …",
              flush=True)
        t0 = time.time()
        m = build_model([wav])
        pmd.compute_model(
            m,
            force        = True,
            rename       = True,
            output_name  = _dap_name(iw, wav),
            nmug         = N_MUG,
            nmug_mie     = N_MUG_MIE,
            nsubr        = N_SUBR,
            nmat         = 4,
            path_input   = str(DAP_DB) + "/",
        )
        print(f"    done in {time.time()-t0:.1f}s  ->  {m.name[0]}")

    ready = fpath.exists() and fpath.stat().st_size > 100
    print(f"\nDAP file ready: {int(ready)}/1")
    return ready


# ── Step 2: disk integration over all phase angles (resumable) ────────────────

def step_integ():
    """Disk-integrate F and Q over all phase angles using a single wavelength.

    Since the Fresnel surface has no wavelength dependence, we compute at
    COMPUTE_WAV_UM only (1 wavelength instead of 5) then broadcast to all
    five wavelengths for the final plot — 5× speedup.
    """
    iw0, wav0 = 0, float(COMPUTE_WAV_UM[0])
    if not (_dap_file(iw0, wav0).exists() and _dap_file(iw0, wav0).stat().st_size > 100):
        print(f"Missing Fourier file for λ={wav0:.3f} µm. Run --step dap first.")
        return False

    n_alpha = len(ALPHA_DEG)
    # Store as (1, n_alpha) — single wavelength
    F1 = np.full((1, n_alpha), np.nan)
    Q1 = np.full((1, n_alpha), np.nan)
    start_idx = 0

    if PART_NPZ.exists() and PART_NPZ.stat().st_size > 100:
        dp = np.load(PART_NPZ)
        F1 = dp["F"]; Q1 = dp["Q"]
        done = ~np.isnan(F1[0])
        start_idx = int(np.sum(done))
        print(f"Resuming from angle index {start_idx}/{n_alpha}")

    if start_idx >= n_alpha:
        # Broadcast single wavelength to all five for the final save
        F = np.tile(F1, (N_WAV, 1)); Q = np.tile(Q1, (N_WAV, 1))
        np.savez(OUT_NPZ, F=F, Q=Q, alpha_deg=ALPHA_DEG, wavelengths_um=WAVELENGTHS_UM)
        PART_NPZ.unlink(missing_ok=True)
        print(f"All angles done. Saved -> {OUT_NPZ}")
        return True

    # Process one angle at a time and save after each — safe for 44s bash timeout
    ANGLE_BUDGET = 3   # max angles per call (~11s each at nmug=80, 3 = ~33s)
    t_call = time.time()
    angles_done_this_call = 0

    for idx in range(start_idx, n_alpha):
        alph = ALPHA_DEG[idx]
        print(f"  α={alph:.0f}°  [{idx+1}/{n_alpha}] …", flush=True)

        m = build_model([wav0])
        m.name = [str(_dap_file(iw0, wav0))]

        t0 = time.time()
        pmd.planet_integrated(
            [m],
            alpha        = [alph],
            npix         = N_PIX,
            output_names = ["fonly50"],
            nmug         = N_MUG,
            nmug_mie     = N_MUG_MIE,
            nsubr        = N_SUBR,
            nmat         = 4,
            force        = False,   # skip compute_model — file already exists
        )
        dt = time.time() - t0

        # m.I / m.Q have shape (1, 1) — single wavelength, single angle
        F1[0, idx] = float(np.array(m.I).ravel()[0])
        Q1[0, idx] = float(np.array(m.Q).ravel()[0])

        np.savez(PART_NPZ, F=F1, Q=Q1, alpha_deg=ALPHA_DEG, wavelengths_um=COMPUTE_WAV_UM)
        print(f"    done in {dt:.1f}s  F={F1[0,idx]:.5f}  Q={Q1[0,idx]:.5f}  ({idx+1}/{n_alpha} saved)", flush=True)

        angles_done_this_call += 1
        if angles_done_this_call >= ANGLE_BUDGET:
            remaining = n_alpha - (idx + 1)
            if remaining > 0:
                print(f"  {remaining} angles remaining — re-run to continue.")
                return False
            break

    # All done — broadcast single wavelength to all five
    F_all = np.tile(F1, (N_WAV, 1)); Q_all = np.tile(Q1, (N_WAV, 1))
    np.savez(OUT_NPZ, F=F_all, Q=Q_all, alpha_deg=ALPHA_DEG, wavelengths_um=WAVELENGTHS_UM)
    PART_NPZ.unlink(missing_ok=True)
    print(f"All angles done! (single λ broadcast to {N_WAV}) Saved -> {OUT_NPZ}")
    return True


# ── Step 3: plot ──────────────────────────────────────────────────────────────

def step_plot():
    if not OUT_NPZ.exists():
        print("No data yet. Run --step dap then --step integ first.")
        return

    d = np.load(OUT_NPZ)
    F  = d["F"];  Q = d["Q"];  alpha = d["alpha_deg"]

    # P = -Q/F, masked where F is negligible
    P = np.full_like(F, np.nan)
    for i in range(N_WAV):
        Fmax = np.nanmax(np.abs(F[i]))
        good = np.abs(F[i]) > 5e-4 * max(Fmax, 1e-10)
        P[i, good] = -Q[i, good] / F[i, good]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    fig.suptitle(
        f"Trees & Stam (2019) Fig. 1 dashed curve — Fresnel-only\n"
        f"(v={WIND_SPEED} m/s, foam, n_water={N_WATER}, nmug={N_MUG}, "
        f"n_fourier={N_FOURIER}, npix={N_PIX})",
        fontsize=10,
    )

    titles = ["Total flux $F / (\\pi F_0)$",
              "Polarized flux $Q / (\\pi F_0)$",
              "Degree of polarization $P$"]
    ylims  = [(0, 0.30), (-0.06, 0.02), (-0.05, 1.0)]

    for col, (ax, title, ylim) in enumerate(zip(axes, titles, ylims)):
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0, 180); ax.set_xticks([0,30,60,90,120,150,180])
        ax.set_ylim(*ylim)
        ax.tick_params(direction="in", top=True, right=True, labelsize=8)
        ax.set_xlabel("Phase angle α (°)", fontsize=9)
        if col == 1:
            ax.axhline(0, color="grey", lw=0.6, ls=":")
        if col == 2:
            ax.axhline(0, color="grey", lw=0.6, ls=":")

        for iw, wl in enumerate(WAVELENGTHS_UM):
            kw = dict(color=WAV_COLORS[wl], lw=1.6, label=WAV_LABELS[wl])
            if col == 0:
                ax.plot(alpha, F[iw], **kw)
            elif col == 1:
                ax.plot(alpha, Q[iw], **kw)
            else:
                ax.plot(alpha, P[iw], **kw)

    axes[0].legend(fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {OUT_PNG}")

    # Print key values for comparison with paper
    print("\n=== Key values vs paper's dashed curve ===")
    idx90 = np.argmin(np.abs(alpha - 90))
    idx60 = np.argmin(np.abs(alpha - 60))
    print(f"{'wl':>6}  F(0°)    F(90°)   Q(60°)   P(90°)")
    for i, wl in enumerate(WAVELENGTHS_NM):
        print(f"{int(wl):>6}  {F[i,0]:.4f}  {F[i,idx90]:.4f}  "
              f"{Q[i,idx60]:.5f}  {P[i,idx90]:.3f}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", choices=["dap","integ","plot","all"], default="all")
    args = ap.parse_args()

    print(f"=== Fresnel-only (Trees & Stam dashed curve) ===")
    print(f"nmug={N_MUG}  n_fourier={N_FOURIER}  n_phi={N_PHI}  npix={N_PIX}")

    if args.step in ("dap", "all"):
        done = step_dap()
        if args.step == "all" and not done:
            return 1

    if args.step in ("integ", "all"):
        done = step_integ()
        if not done:
            return 1   # more angles to go; re-run

    if args.step in ("plot", "all") or OUT_NPZ.exists():
        if OUT_NPZ.exists():
            step_plot()

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
