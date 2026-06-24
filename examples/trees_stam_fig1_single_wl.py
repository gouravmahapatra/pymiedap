#!/usr/bin/env python3
"""Reproduce Trees & Stam (2019) Figure 1 at a single wavelength.

Three models are computed for the chosen wavelength:
  1. Ocean planet  — Rayleigh atmosphere + rough Fresnel + polarized water body
  2. Black surface — same atmosphere, asurf = 0 (Lambertian, no ocean)
  3. Fresnel-only  — rough Fresnel surface, no atmosphere, no water body (dashed)

Usage:
  python trees_stam_fig1_single_wl.py --wl 350  --step dap
  python trees_stam_fig1_single_wl.py --wl 350  --step integ
  python trees_stam_fig1_single_wl.py --wl 350  --step plot
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
# Atmosphere (Tables 1 & 2 of the paper)
BM_550    = 0.096      # Rayleigh OD at 550 nm
DEPOL     = 0.03
MMA       = 29.0       # g/mol
GRAVITY   = 9.81       # m/s²

# Ocean (Table 2)
WIND_SPEED  = 7.0
FOAM_ALBEDO = 0.22
N_WATER     = 1.33
N_AIR       = 1.0
DEPTH_M     = 100.0
BOTTOM_ALB  = 0.0
WATER_DEPOL = 0.09

# Computation
N_MUG      = 80
N_MUG_MIE  = 16
N_SUBR     = 50
N_FOURIER  = 80
N_PHI      = 360
ALPHA_DEG  = np.arange(0., 181., 3.)   # 61 phase angles
N_PIX      = 40

DAP_DB = REPO_ROOT / "dap_database"
CACHE  = Path("/tmp/pymiedap_cache")
DAP_DB.mkdir(exist_ok=True); CACHE.mkdir(exist_ok=True)

MODELS = ("ocean", "black", "fonly")


def _bm(wav_um: float) -> float:
    """Rayleigh optical depth at wavelength wav_um [µm]."""
    return BM_550 * (0.550 / wav_um) ** 4


# ── Model builders ─────────────────────────────────────────────────────────────

def _base_atm(wav_um: float) -> pmd.Model:
    """Cloud-free Earth-like atmosphere at the given wavelength."""
    m = pmd.Model()
    m.mma = MMA; m.gravity = GRAVITY; m.dpol = DEPOL; m.asurf = 0.0
    del m.layers.haze, m.layers.cloud
    m.wvl_list = [wav_um]

    bm = _bm(wav_um)
    m.layers.gasbelow.press   = 1.0
    m.layers.gasbelow.tau     = [0.0]
    m.layers.gasbelow.tau_g   = [0.0]
    m.layers.gasbelow.tau_ray = [bm]
    m.layers.gasbelow.rayscat = False

    m.layers.gastop.press   = 1.0e-5
    m.layers.gastop.tau     = [0.0]
    m.layers.gastop.tau_g   = [0.0]
    m.layers.gastop.tau_ray = [0.0]
    m.layers.gastop.rayscat = False
    return m


def build_ocean(wav_um: float) -> pmd.Model:
    """Ocean planet: atmosphere + rough Fresnel + polarized pure-water body."""
    m = _base_atm(wav_um)
    m.surface = OceanSurface(
        wind_speed         = WIND_SPEED,
        wavelength_um      = wav_um,
        n_air              = N_AIR,
        n_water            = N_WATER,
        foam_albedo        = FOAM_ALBEDO,
        depth_m            = DEPTH_M,
        bottom_albedo      = BOTTOM_ALB,
        water_depol        = WATER_DEPOL,
        n_fourier          = N_FOURIER,
        n_phi              = N_PHI,
        water_streams      = 8,
        water_n_phi        = 180,
        water_initial_tau  = 0.01,
        include_subsurface = True,
        include_foam       = True,
        solver             = "adding_doubling",
    )
    return m


def build_black(wav_um: float) -> pmd.Model:
    """Black-surface planet: same atmosphere, asurf=0."""
    m = _base_atm(wav_um)
    m.asurf = 0.0
    return m


def build_fonly(wav_um: float) -> pmd.Model:
    """Fresnel-only: rough surface, no atmosphere, no water body."""
    m = pmd.Model()
    m.mma = MMA; m.gravity = GRAVITY; m.dpol = DEPOL; m.asurf = 0.0
    del m.layers.haze, m.layers.cloud
    m.wvl_list = [wav_um]
    for lname in ("gasbelow", "gastop"):
        lyr = getattr(m.layers, lname)
        lyr.press = 1e-8 if lname == "gasbelow" else 1e-10
        lyr.tau = [0.0]; lyr.tau_g = [0.0]; lyr.tau_ray = [0.0]
        lyr.rayscat = False
    m.surface = OceanSurface(
        wind_speed         = WIND_SPEED,
        n_air              = N_AIR,
        n_water            = N_WATER,
        foam_albedo        = FOAM_ALBEDO,
        n_fourier          = N_FOURIER,
        n_phi              = N_PHI,
        include_subsurface = False,
        include_foam       = True,
        solver             = "adding_doubling",
    )
    return m


BUILDERS = {"ocean": build_ocean, "black": build_black, "fonly": build_fonly}


# ── File helpers ──────────────────────────────────────────────────────────────

def _dap_path(tag: str, wav_nm: int) -> Path:
    return DAP_DB / f"ts19_{tag}_{wav_nm}nm_{N_MUG}mug_{N_FOURIER}fou_{wav_nm/1000:.7f}.dat"

def _npz_path(tag: str, wav_nm: int) -> Path:
    return CACHE / f"ts19_{tag}_{wav_nm}nm_{N_MUG}mug_{N_FOURIER}fou.npz"

def _part_path(tag: str, wav_nm: int) -> Path:
    return CACHE / f"ts19_{tag}_{wav_nm}nm_{N_MUG}mug_{N_FOURIER}fou.partial.npz"


# ── Step 1: DAP ───────────────────────────────────────────────────────────────

def step_dap(wav_um: float):
    wav_nm = int(round(wav_um * 1000))
    print(f"=== DAP step  λ={wav_nm} nm  nmug={N_MUG}  n_fourier={N_FOURIER} ===")
    for tag in MODELS:
        fpath = _dap_path(tag, wav_nm)
        if fpath.exists() and fpath.stat().st_size > 100:
            print(f"  [{tag}] Fourier file exists, skipping.")
            continue
        print(f"  [{tag}] computing …", flush=True)
        t0 = time.time()
        m = BUILDERS[tag](wav_um)
        oname = f"ts19_{tag}_{wav_nm}nm_{N_MUG}mug_{N_FOURIER}fou"
        pmd.compute_model(m, force=True, rename=True, output_name=oname,
                          nmug=N_MUG, nmug_mie=N_MUG_MIE, nsubr=N_SUBR,
                          nmat=4, path_input=str(DAP_DB) + "/")
        print(f"    done in {time.time()-t0:.1f}s  ->  {m.name[0]}")
    ready = all(_dap_path(t, wav_nm).exists() and _dap_path(t, wav_nm).stat().st_size > 100
                for t in MODELS)
    print(f"\nDAP files ready: {sum(_dap_path(t,wav_nm).exists() and _dap_path(t,wav_nm).stat().st_size>100 for t in MODELS)}/{len(MODELS)}")
    return ready


# ── Step 2: Integration ───────────────────────────────────────────────────────

def step_integ(wav_um: float):
    wav_nm = int(round(wav_um * 1000))
    print(f"=== Integration  λ={wav_nm} nm ===")
    n_alpha = len(ALPHA_DEG)
    remaining = []

    for tag in MODELS:
        if _npz_path(tag, wav_nm).exists() and _npz_path(tag, wav_nm).stat().st_size > 100:
            print(f"  [{tag}] already complete, skipping.")
            continue

        F1 = np.full((1, n_alpha), np.nan); Q1 = np.full((1, n_alpha), np.nan)
        start_idx = 0
        pp = _part_path(tag, wav_nm)
        if pp.exists() and pp.stat().st_size > 100:
            dp = np.load(pp); F1 = dp["F"]; Q1 = dp["Q"]
            done = ~np.isnan(F1[0]); start_idx = int(done.sum())
            print(f"  [{tag}] Resuming from angle {start_idx}/{n_alpha}")

        if start_idx >= n_alpha:
            np.savez(_npz_path(tag, wav_nm), F=F1, Q=Q1, alpha_deg=ALPHA_DEG)
            pp.unlink(missing_ok=True)
            print(f"  [{tag}] All done, saved.")
            continue

        # Batch: ~2 angles per 44s call (each angle ~20s at 350 nm)
        BATCH = 2
        batch_alpha = list(ALPHA_DEG[start_idx:start_idx + BATCH])
        end_idx = start_idx + len(batch_alpha)

        print(f"  [{tag}] α={batch_alpha[0]:.0f}°–{batch_alpha[-1]:.0f}°  "
              f"({end_idx}/{n_alpha}) …", flush=True)

        m = BUILDERS[tag](wav_um)
        m.name = [str(_dap_path(tag, wav_nm))]
        t0 = time.time()
        pmd.planet_integrated([m], alpha=batch_alpha, npix=N_PIX,
                              output_names=[f"ts19_{tag}_{wav_nm}nm"],
                              nmug=N_MUG, nmug_mie=N_MUG_MIE, nsubr=N_SUBR,
                              nmat=4, force=False)
        F1[:, start_idx:end_idx] = np.array(m.I)
        Q1[:, start_idx:end_idx] = np.array(m.Q)
        np.savez(pp, F=F1, Q=Q1, alpha_deg=ALPHA_DEG)
        print(f"    batch done in {time.time()-t0:.0f}s")

        if end_idx >= n_alpha:
            np.savez(_npz_path(tag, wav_nm), F=F1, Q=Q1, alpha_deg=ALPHA_DEG)
            pp.unlink(missing_ok=True)
            print(f"  [{tag}] All angles done!  Saved -> {_npz_path(tag, wav_nm)}")
        else:
            remaining.append(tag)

    if remaining:
        print(f"\nAngles remaining for: {remaining} — re-run to continue.")
        return False
    all_done = all(_npz_path(t, wav_nm).exists() and _npz_path(t,wav_nm).stat().st_size > 100
                   for t in MODELS)
    return all_done


# ── Step 3: Plot ──────────────────────────────────────────────────────────────

def step_plot(wav_um: float):
    wav_nm = int(round(wav_um * 1000))
    print(f"=== Plot  λ={wav_nm} nm ===")

    data = {}
    for tag in MODELS:
        fp = _npz_path(tag, wav_nm)
        if not (fp.exists() and fp.stat().st_size > 100):
            print(f"  Missing data for [{tag}]. Run --step integ first.")
            return
        d = np.load(fp); data[tag] = (d["F"][0], d["Q"][0])

    alpha = np.load(_npz_path("ocean", wav_nm))["alpha_deg"]
    bm = _bm(wav_um)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    fig.suptitle(
        f"Trees & Stam (2019) Fig. 1 — λ = {wav_nm} nm  "
        f"[bm = {bm:.3f}, nmug={N_MUG}, n_fourier={N_FOURIER}]\n"
        f"v=7 m/s · depth=100 m · ps=1 bar · cloud-free", fontsize=10)

    styles = {
        "ocean": dict(color="#2980B9", lw=1.8, label=f"Ocean planet ({wav_nm} nm)"),
        "black": dict(color="#e74c3c",  lw=1.8, label=f"Black surface ({wav_nm} nm)"),
        "fonly": dict(color="black",    lw=1.4, ls="--",
                      label="Fresnel-only (no atm, no ocean)"),
    }

    for col, (ytitle, ylim, key) in enumerate([
        ("Total flux $F/(\\pi F_0)$",     (0, 0.30),     "F"),
        ("Polarized flux $Q/(\\pi F_0)$", (-0.06, 0.02), "Q"),
        ("Degree of polarization $P$",    (-0.05, 1.00), "P"),
    ]):
        ax = axes[col]
        ax.set_title(ytitle, fontsize=10)
        ax.set_xlim(0, 180); ax.set_xticks([0,30,60,90,120,150,180])
        ax.set_ylim(*ylim)
        ax.tick_params(direction="in", top=True, right=True, labelsize=8)
        ax.set_xlabel("Phase angle α (°)", fontsize=9)
        if col in (1, 2):
            ax.axhline(0, color="grey", lw=0.6, ls=":")

        for tag in MODELS:
            F, Q = data[tag]
            if key == "F":
                ax.plot(alpha, F, **styles[tag])
            elif key == "Q":
                ax.plot(alpha, Q, **styles[tag])
            else:
                good = np.abs(F) > 5e-4 * np.nanmax(np.abs(F))
                P = np.full_like(F, np.nan)
                P[good] = -Q[good] / F[good]
                ax.plot(alpha, P, **styles[tag])

    axes[0].legend(fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()

    out_png = REPO_ROOT / "examples" / f"trees_stam_fig1_{wav_nm}nm.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out_png}")

    # Print key values
    print(f"\nKey values at λ={wav_nm} nm:")
    idx90  = np.argmin(np.abs(alpha - 90))
    idx150 = np.argmin(np.abs(alpha - 150))
    print(f"  {'Model':<10}  {'F(α=0°)':>9}  {'F(α=90°)':>9}  {'Q(α=90°)':>10}  {'P(α=90°)':>9}")
    for tag in MODELS:
        F, Q = data[tag]
        F0  = F[0]; F90 = F[idx90]; Q90 = Q[idx90]
        P90 = -Q90/F90 if abs(F90) > 1e-6 else float("nan")
        print(f"  {tag:<10}  {F0:>9.5f}  {F90:>9.5f}  {Q90:>10.6f}  {P90:>+9.3f}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wl",   type=float, default=350.0,
                    help="Wavelength in nm (default: 350)")
    ap.add_argument("--step", choices=["dap","integ","plot","all"], default="all")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    wav_um = args.wl / 1000.0
    wav_nm = int(round(args.wl))

    print(f"=== Trees & Stam (2019) Fig. 1 — λ={wav_nm} nm ===")
    print(f"Rayleigh bm = {_bm(wav_um):.4f}  (vs 0.096 at 550 nm)")

    if args.force:
        for tag in MODELS:
            for p in [_dap_path(tag, wav_nm), _npz_path(tag, wav_nm),
                      _part_path(tag, wav_nm)]:
                if p.exists(): open(p,"wb").close()
        print("  --force: caches cleared.")

    if args.step in ("dap",  "all"): step_dap(wav_um)
    if args.step in ("integ","all"):
        done = step_integ(wav_um)
        if not done: return 1
    if args.step in ("plot", "all"):
        step_plot(wav_um)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
