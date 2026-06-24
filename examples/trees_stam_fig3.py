#!/usr/bin/env python3
"""Reproduce Trees & Stam (2019) Figure 3.

F, Q, P vs phase angle for the OCEAN PLANET at three surface pressures:
  top row    : ps = 0.5 bar
  middle row : ps = 5.0 bar
  bottom row : ps = 10.0 bar

Five wavelengths (350, 443, 550, 670, 865 nm) are shown per panel.
No black surface, no Fresnel-only dashed line.

Usage:
  python trees_stam_fig3.py --step dap      # DAP for all 3 pressures (3 calls)
  python trees_stam_fig3.py --step integ    # disk integration (resumable)
  python trees_stam_fig3.py --step plot     # assemble 3×3 figure
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

# ── Model parameters (Tables 1 & 2) ───────────────────────────────────────────
WAVELENGTHS_UM = np.array([0.350, 0.443, 0.550, 0.670, 0.865])
WAVELENGTHS_NM = (WAVELENGTHS_UM * 1000).astype(int)
N_WAV = len(WAVELENGTHS_UM)

SURFACE_PRESSURES = [0.5, 5.0, 10.0]   # bar — Figure 3 rows

BM_550   = 0.096   # Rayleigh OD at 550 nm, ps=1 bar
DEPOL    = 0.03
MMA      = 29.0
GRAVITY  = 9.81

WIND_SPEED  = 7.0
FOAM_ALBEDO = 0.22
N_WATER     = 1.33
N_AIR       = 1.0
DEPTH_M     = 100.0
BOTTOM_ALB  = 0.0
WATER_DEPOL = 0.09

# Computation settings (user-specified)
N_MUG      = 50
N_MUG_MIE  = 16
N_SUBR     = 50
N_FOURIER  = 50
N_PHI      = 300   # min = 2*50+1 = 101; 300 gives good accuracy
ALPHA_DEG  = np.arange(0., 181., 3.)
N_PIX      = 40

DAP_DB = REPO_ROOT / "dap_database"
CACHE  = Path("/tmp/pymiedap_cache")
DAP_DB.mkdir(exist_ok=True); CACHE.mkdir(exist_ok=True)

# Wavelength colours matching Trees & Stam (magenta→blue→green→red→brown)
WAV_COLORS = {350: "#C0392B", 443: "#2980B9", 550: "#27AE60",
              670: "#E67E22", 865: "#7D3C98"}
WAV_LABELS = {350: "350 nm", 443: "443 nm", 550: "550 nm",
              670: "670 nm", 865: "865 nm"}


# ── Helpers ────────────────────────────────────────────────────────────────────

def _bm(wav_um: float, ps: float) -> float:
    """Rayleigh optical depth at wavelength wav_um [µm] and surface pressure ps [bar]."""
    return BM_550 * ps * (0.550 / wav_um) ** 4


def _ps_tag(ps: float) -> str:
    return f"ps{ps:.1f}".replace(".", "p")


def _dap_names_file(ps: float) -> Path:
    """Text file storing the 5 DAP Fourier-file paths created by compute_model."""
    return CACHE / f"fig3_ocean_{_ps_tag(ps)}_{N_MUG}mug_{N_FOURIER}fou.names"


def _dap_ready(ps: float) -> bool:
    """True when all 5 wavelength Fourier files exist and are non-empty."""
    nf = _dap_names_file(ps)
    if not (nf.exists() and nf.stat().st_size > 0):
        return False
    for p in nf.read_text().strip().split("\n"):
        if not (Path(p).exists() and Path(p).stat().st_size > 100):
            return False
    return True


def _load_dap_names(ps: float):
    return _dap_names_file(ps).read_text().strip().split("\n")


def _npz_path(ps: float) -> Path:
    return CACHE / f"fig3_ocean_{_ps_tag(ps)}_{N_MUG}mug_{N_FOURIER}fou.npz"


def _part_path(ps: float) -> Path:
    return CACHE / f"fig3_ocean_{_ps_tag(ps)}_{N_MUG}mug_{N_FOURIER}fou.partial.npz"


def _dap_name(ps: float) -> str:
    return f"fig3_ocean_{_ps_tag(ps)}_{N_MUG}mug_{N_FOURIER}fou"


# ── Model builder ──────────────────────────────────────────────────────────────

def build_ocean(ps: float) -> pmd.Model:
    """Ocean planet at surface pressure ps [bar], all 5 wavelengths."""
    m = pmd.Model()
    m.mma = MMA; m.gravity = GRAVITY; m.dpol = DEPOL; m.asurf = 0.0
    del m.layers.haze, m.layers.cloud
    m.wvl_list = list(WAVELENGTHS_UM)

    # Rayleigh OD scaled by surface pressure
    bm_vec = [_bm(w, ps) for w in WAVELENGTHS_UM]
    m.layers.gasbelow.press   = float(ps)
    m.layers.gasbelow.tau     = [0.0] * N_WAV
    m.layers.gasbelow.tau_g   = [0.0] * N_WAV
    m.layers.gasbelow.tau_ray = bm_vec
    m.layers.gasbelow.rayscat = False

    m.layers.gastop.press   = 1.0e-5
    m.layers.gastop.tau     = [0.0] * N_WAV
    m.layers.gastop.tau_g   = [0.0] * N_WAV
    m.layers.gastop.tau_ray = [0.0] * N_WAV
    m.layers.gastop.rayscat = False

    m.surface = OceanSurface(
        wind_speed         = WIND_SPEED,
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


# ── DAP helpers (per-wavelength resumable) ────────────────────────────────────

def _dap_file_for_wav(ps: float, wav_um: float) -> Path:
    """Expected path of the Fourier .dat file for one (ps, wavelength) pair."""
    return DAP_DB / f"{_dap_name(ps)}_{wav_um:.7f}.dat"


def _dap_wav_done(ps: float, wav_um: float) -> bool:
    p = _dap_file_for_wav(ps, wav_um)
    return p.exists() and p.stat().st_size > 100


def _dap_ready(ps: float) -> bool:
    """True when all 5 wavelength Fourier files exist, are non-empty, and names file written."""
    if not all(_dap_wav_done(ps, w) for w in WAVELENGTHS_UM):
        return False
    nf = _dap_names_file(ps)
    if not (nf.exists() and nf.stat().st_size > 0):
        # Regenerate names file from known paths
        paths = [str(_dap_file_for_wav(ps, w)) for w in WAVELENGTHS_UM]
        nf.write_text("\n".join(paths))
    return True


def build_ocean_1wav(ps: float, wav_um: float) -> pmd.Model:
    """Ocean planet at surface pressure ps [bar], single wavelength."""
    m = pmd.Model()
    m.mma = MMA; m.gravity = GRAVITY; m.dpol = DEPOL; m.asurf = 0.0
    del m.layers.haze, m.layers.cloud
    m.wvl_list = [float(wav_um)]

    bm_vec = [_bm(wav_um, ps)]
    m.layers.gasbelow.press   = float(ps)
    m.layers.gasbelow.tau     = [0.0]
    m.layers.gasbelow.tau_g   = [0.0]
    m.layers.gasbelow.tau_ray = bm_vec
    m.layers.gasbelow.rayscat = False

    m.layers.gastop.press   = 1.0e-5
    m.layers.gastop.tau     = [0.0]
    m.layers.gastop.tau_g   = [0.0]
    m.layers.gastop.tau_ray = [0.0]
    m.layers.gastop.rayscat = False

    m.surface = OceanSurface(
        wind_speed         = WIND_SPEED,
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


# ── Step 1: DAP ────────────────────────────────────────────────────────────────

def step_dap():
    """Compute DAP Fourier files one wavelength at a time (resumable).

    Each call computes the next missing wavelength for the first incomplete
    pressure, then returns.  Re-run until "DAP files ready: 3/3".
    """
    print(f"=== DAP  nmug={N_MUG}  n_fourier={N_FOURIER} ===")

    did_work = False
    for ps in SURFACE_PRESSURES:
        if _dap_ready(ps):
            print(f"  ps={ps:.1f} bar: all wavelengths done, skipping.")
            continue

        # Find the first missing wavelength
        for wav_um in WAVELENGTHS_UM:
            if _dap_wav_done(ps, wav_um):
                continue
            wav_nm = int(round(wav_um * 1000))
            print(f"  ps={ps:.1f} bar  λ={wav_nm}nm: computing …", flush=True)
            t0 = time.time()
            m = build_ocean_1wav(ps, wav_um)
            pmd.compute_model(m, force=True, rename=True,
                              output_name=_dap_name(ps),
                              nmug=N_MUG, nmug_mie=N_MUG_MIE, nsubr=N_SUBR,
                              nmat=4, path_input=str(DAP_DB) + "/")
            print(f"    done in {time.time()-t0:.1f}s  → {m.name[0]}")
            did_work = True
            break   # one wavelength per call; re-run for the next

        if did_work:
            break   # one wavelength total per call

    # Rebuild names files for any now-complete pressures
    for ps in SURFACE_PRESSURES:
        if all(_dap_wav_done(ps, w) for w in WAVELENGTHS_UM):
            paths = [str(_dap_file_for_wav(ps, w)) for w in WAVELENGTHS_UM]
            _dap_names_file(ps).write_text("\n".join(paths))

    ready = sum(_dap_ready(ps) for ps in SURFACE_PRESSURES)
    print(f"\nDAP files ready: {ready}/{len(SURFACE_PRESSURES)}")
    return ready == len(SURFACE_PRESSURES)


# ── Step 2: Integration ────────────────────────────────────────────────────────

def step_integ():
    """Disk-integrate F and Q for all pressures (resumable batches)."""
    print(f"=== Integration  nmug={N_MUG}  n_fourier={N_FOURIER} ===")
    n_alpha = len(ALPHA_DEG)
    remaining = []

    for ps in SURFACE_PRESSURES:
        # Check complete
        if _npz_path(ps).exists() and _npz_path(ps).stat().st_size > 100:
            print(f"  ps={ps:.1f} bar: already complete, skipping.")
            continue

        # Check DAP files
        if not _dap_ready(ps):
            print(f"  ps={ps:.1f} bar: DAP files missing — run --step dap first.")
            remaining.append(ps)
            continue

        # Load partial progress
        F = np.full((N_WAV, n_alpha), np.nan)
        Q = np.full((N_WAV, n_alpha), np.nan)
        start_idx = 0
        pp = _part_path(ps)
        if pp.exists() and pp.stat().st_size > 100:
            dp = np.load(pp); F = dp["F"]; Q = dp["Q"]
            start_idx = int((~np.isnan(F[0])).sum())
            print(f"  ps={ps:.1f} bar: Resuming from angle {start_idx}/{n_alpha}")

        if start_idx >= n_alpha:
            np.savez(_npz_path(ps), F=F, Q=Q,
                     alpha_deg=ALPHA_DEG, wavelengths_um=WAVELENGTHS_UM, ps=ps)
            pp.unlink(missing_ok=True)
            print(f"  ps={ps:.1f} bar: all done, saved.")
            continue

        # Process one batch
        BATCH = 2   # ~20s/angle at nmug=50 → 2 angles (~39s) fits in 44s
        batch = list(ALPHA_DEG[start_idx:start_idx + BATCH])
        end_idx = start_idx + len(batch)
        print(f"  ps={ps:.1f} bar: α={batch[0]:.0f}°–{batch[-1]:.0f}°  "
              f"({end_idx}/{n_alpha}) …", flush=True)

        m = build_ocean(ps)
        m.name = _load_dap_names(ps)   # 5 per-wavelength Fourier files

        t0 = time.time()
        pmd.planet_integrated([m], alpha=batch, npix=N_PIX,
                              output_names=[_dap_name(ps)],
                              nmug=N_MUG, nmug_mie=N_MUG_MIE, nsubr=N_SUBR,
                              nmat=4, force=False)
        F[:, start_idx:end_idx] = np.array(m.I)
        Q[:, start_idx:end_idx] = np.array(m.Q)
        np.savez(pp, F=F, Q=Q, alpha_deg=ALPHA_DEG,
                 wavelengths_um=WAVELENGTHS_UM, ps=ps)
        print(f"    batch done in {time.time()-t0:.0f}s")

        if end_idx >= n_alpha:
            np.savez(_npz_path(ps), F=F, Q=Q,
                     alpha_deg=ALPHA_DEG, wavelengths_um=WAVELENGTHS_UM, ps=ps)
            pp.unlink(missing_ok=True)
            print(f"  ps={ps:.1f} bar: ALL ANGLES DONE  →  {_npz_path(ps)}")
        else:
            remaining.append(ps)

    if remaining:
        print(f"\nStill running: ps={remaining} — re-run to continue.")
        return False
    all_done = all(_npz_path(ps).exists() and _npz_path(ps).stat().st_size > 100
                   for ps in SURFACE_PRESSURES)
    return all_done


# ── Step 3: Plot ───────────────────────────────────────────────────────────────

def step_plot():
    """Produce the 3×3 Figure 3 panel."""
    print("=== Plot ===")
    data = {}
    for ps in SURFACE_PRESSURES:
        fp = _npz_path(ps)
        if not (fp.exists() and fp.stat().st_size > 100):
            print(f"  Missing data for ps={ps:.1f} bar. Run --step integ first.")
            return
        d = np.load(fp)
        data[ps] = (d["F"], d["Q"])

    alpha = np.load(_npz_path(SURFACE_PRESSURES[0]))["alpha_deg"]

    fig, axes = plt.subplots(3, 3, figsize=(13, 11), sharex=True)
    fig.suptitle(
        f"Trees & Stam (2019) Fig. 3 — Ocean planet, varying surface pressure\n"
        f"nmug={N_MUG}  n_fourier={N_FOURIER}  v=7 m/s  depth=100 m  cloud-free",
        fontsize=11,
    )

    row_labels  = ["$p_s = 0.5$ bar", "$p_s = 5$ bar", "$p_s = 10$ bar"]
    col_titles  = ["Total flux $F/(\\pi F_0)$",
                   "Polarized flux $Q/(\\pi F_0)$",
                   "Degree of polarization $P$"]
    ylims = [(0, 0.30), (-0.06, 0.02), (0.0, 1.0)]

    for row, ps in enumerate(SURFACE_PRESSURES):
        F, Q = data[ps]

        for col in range(3):
            ax = axes[row, col]
            if row == 0:
                ax.set_title(col_titles[col], fontsize=10)
            if col == 0:
                ax.set_ylabel(row_labels[row] + "\n", fontsize=9)
            ax.set_xlim(0, 180)
            ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
            ax.set_ylim(*ylims[col])
            ax.tick_params(direction="in", top=True, right=True, labelsize=8)
            if col in (1, 2):
                ax.axhline(0, color="grey", lw=0.6, ls=":")
            if row == 2:
                ax.set_xlabel("Phase angle α (°)", fontsize=9)

            for iw, (wl_um, wl_nm) in enumerate(zip(WAVELENGTHS_UM, WAVELENGTHS_NM)):
                kw = dict(color=WAV_COLORS[wl_nm], lw=1.5,
                          label=WAV_LABELS[wl_nm] if row == 0 else "_")
                if col == 0:
                    ax.plot(alpha, F[iw], **kw)
                elif col == 1:
                    ax.plot(alpha, Q[iw], **kw)
                else:
                    good = np.abs(F[iw]) > 5e-4 * np.nanmax(np.abs(F[iw]))
                    P = np.full(len(alpha), np.nan)
                    P[good] = -Q[iw, good] / F[iw, good]
                    ax.plot(alpha, P, **kw)

    # Single legend at top
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5,
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.01))

    fig.tight_layout(rect=[0, 0.04, 1, 1])
    out = REPO_ROOT / "examples" / f"trees_stam_fig3_nmug{N_MUG}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")

    # Key values
    print("\nKey values (P at α=90°):")
    alpha90 = np.argmin(np.abs(alpha - 90))
    print(f"  {'ps (bar)':<10}  " + "  ".join(f"{wl:>6}nm" for wl in WAVELENGTHS_NM))
    for ps in SURFACE_PRESSURES:
        F, Q = data[ps]
        Ps = [-Q[iw, alpha90] / F[iw, alpha90] for iw in range(N_WAV)]
        print(f"  {ps:<10.1f}  " + "  ".join(f"{p:>7.3f}" for p in Ps))


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", choices=["dap", "integ", "plot", "all"], default="all")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    print(f"=== Trees & Stam (2019) Figure 3 — nmug={N_MUG}  n_fourier={N_FOURIER} ===")
    print(f"Surface pressures: {SURFACE_PRESSURES} bar")
    print(f"Wavelengths: {list(WAVELENGTHS_NM)} nm")
    print()

    if args.force:
        for ps in SURFACE_PRESSURES:
            for p in [_dap_names_file(ps), _npz_path(ps), _part_path(ps)]:
                if p.exists(): open(p, "wb").close()
            # Zero individual Fourier files if names file exists
            nf = _dap_names_file(ps)
            if nf.exists():
                for fp in nf.read_text().strip().split("\n"):
                    if Path(fp).exists(): open(fp, "wb").close()
        print("  --force: caches cleared.")

    if args.step in ("dap",   "all"): step_dap()
    if args.step in ("integ", "all"):
        done = step_integ()
        if not done: return 1
    if args.step in ("plot",  "all"): step_plot()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
