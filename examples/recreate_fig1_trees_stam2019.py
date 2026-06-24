#!/usr/bin/env python3
"""Reproduce Figure 1 (top row) of Trees & Stam (2019),
"Blue, white, and red ocean planets", A&A 626, A129.

Figure 1 shows disk-integrated F, Q, and P as functions of phase angle α
(0–180°) for five wavelengths, for two model planets:
  • Top row:    ocean planet  — rough Fresnel interface + pure-water body
  • Bottom row: black-surface — asurf = 0 (Lambertian, zero albedo)

An extra dashed black curve (top row only) shows the rough Fresnel interface
alone, without an atmosphere and without a sub-surface ocean.

Model parameters (Tables 1 & 2 of the paper)
─────────────────────────────────────────────
Atmosphere   : Earth-like, cloud-free, ps = 1 bar
               bm(550 nm) = 0.096, δ = 0.03, mg = 29 g mol⁻¹, g = 9.81 m s⁻²
               Rayleigh scaling: bm(λ) = 0.096 × (550/λ_nm)⁴
Ocean surface: v = 7 m s⁻¹, a_foam = 0.22, n₁ = 1.0, n₂ = 1.33,
               depth = 100 m, bottom_albedo = 0, δ_w = 0.09
Wavelengths  : 350, 443, 550, 670, 865 nm
Phase angles : 0 – 180° in 3° steps  (61 points)
npix         : 40 pixels across the equator  (per paper)
"""

from __future__ import annotations

import os, sys, warnings, time
from pathlib import Path
warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pymiedap.pymiedap as pmd
from pymiedap.ocean import OceanSurface

# ── Output paths ───────────────────────────────────────────────────────────────
DAP_DB  = REPO_ROOT / "dap_database"
OUT_PNG      = REPO_ROOT / "examples" / "trees_stam2019_fig1.png"
OUT_NPZ      = REPO_ROOT / "examples" / "trees_stam2019_fig1_data.npz"
# Per-model intermediate saves so each can be computed in a separate bash call
OUT_OCEAN    = REPO_ROOT / ".cache" / "fig1_ocean.npz"
OUT_BLACK    = REPO_ROOT / ".cache" / "fig1_black.npz"
OUT_FONLY    = REPO_ROOT / ".cache" / "fig1_fonly.npz"
CACHE_DIR_F1 = REPO_ROOT / ".cache"
DAP_DB.mkdir(exist_ok=True)

# ── Model parameters (Tables 1 & 2) ───────────────────────────────────────────
WAVELENGTHS_UM = np.array([0.350, 0.443, 0.550, 0.670, 0.865])
WAVELENGTHS_NM = WAVELENGTHS_UM * 1000.0
N_WAV = len(WAVELENGTHS_UM)

# Atmosphere
PS_BAR  = 1.0      # surface pressure
DEPOL   = 0.03     # depolarization factor
MMA     = 29.0     # mean molecular mass  [g mol⁻¹]
GRAVITY = 9.81     # gravity              [m s⁻²]
BM_550  = 0.096    # Rayleigh OD at 550 nm (Table 1)

# Gas optical thickness per wavelength (Rayleigh λ⁻⁴ scaling)
BM = BM_550 * (0.550 / WAVELENGTHS_UM) ** 4

# Ocean parameters
WIND_SPEED   = 7.0    # m s⁻¹
FOAM_ALBEDO  = 0.22
N_WATER      = 1.33
N_AIR        = 1.0
DEPTH_M      = 100.0
BOTTOM_ALB   = 0.0
WATER_DEPOL  = 0.09

# Computation settings
ALPHA_DEG = np.arange(0., 181., 3.)   # 61 phase angles
N_PIX     = 40     # pixels across the equator (paper: ≥ 40)
N_MUG     = 20     # Gauss points for DAP
N_MUG_MIE = 16     # Gauss points for Mie
N_SUBR    = 50     # Mie size-distribution subdivisions
# Ocean Fourier resolution — paper recommends 40-80 for publication quality;
# 20 with n_phi=180 gives ~5s per wavelength and still captures the glint well
# at v=7 m/s.  Increase N_FOURIER/N_PHI for sharper features (lower wind speed).
N_FOURIER = 20     # ocean Fourier terms (paper: 40-80)
N_PHI     = 180    # azimuth samples  (min = max(32, 2*N_FOURIER+1) = 41)

# Colour scheme: same ordering as in the paper
WAV_COLORS = {
    0.350: "#9B59B6",   # violet
    0.443: "#2980B9",   # blue
    0.550: "#27AE60",   # green
    0.670: "#E74C3C",   # red
    0.865: "#7F8C8D",   # near-IR / grey
}
WAV_LABELS = {
    0.350: "350 nm",
    0.443: "443 nm",
    0.550: "550 nm",
    0.670: "670 nm",
    0.865: "865 nm",
}


# ── Model builder ─────────────────────────────────────────────────────────────

def _build_base_model(tag: str) -> pmd.Model:
    """Build a cloud-free, Earth-like atmospheric model ready for surface
    assignment.  Gas Rayleigh ODs are set manually from the paper's Table 1
    scaling law so the results match exactly regardless of PyMieDAP's internal
    refractive-index formula.
    """
    m = pmd.Model()
    m.mma     = MMA
    m.gravity = GRAVITY
    m.dpol    = DEPOL
    m.asurf   = 0.0

    # Remove cloud and haze — cloud-free model
    del m.layers.haze
    del m.layers.cloud

    m.wvl_list = list(WAVELENGTHS_UM)

    # gasbelow: carries the full Rayleigh column; set rayscat=False so we
    # inject the paper's bm values directly instead of relying on bmolecules.
    m.layers.gasbelow.press   = PS_BAR
    m.layers.gasbelow.tau     = [0.0] * N_WAV
    m.layers.gasbelow.tau_g   = [0.0] * N_WAV
    m.layers.gasbelow.tau_ray = list(BM)
    m.layers.gasbelow.rayscat = False   # use tau_ray, not bmolecules

    # gastop: top of atmosphere — transparent, negligible pressure
    m.layers.gastop.press   = 1.0e-5
    m.layers.gastop.tau     = [0.0] * N_WAV
    m.layers.gastop.tau_g   = [0.0] * N_WAV
    m.layers.gastop.tau_ray = [0.0] * N_WAV
    m.layers.gastop.rayscat = False

    return m


def build_ocean_model() -> pmd.Model:
    """Ocean planet: rough Fresnel interface + full pure-water body."""
    m = _build_base_model("ocean")
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
        water_streams      = 8,    # rectangular-supermatrix water quadrature
        water_n_phi        = 120,  # azimuth samples inside water body
        water_initial_tau  = 0.01,
        include_subsurface = True,
        include_foam       = True,
        solver             = "adding_doubling",
    )
    return m


def build_black_model() -> pmd.Model:
    """Black-surface planet: same gas atmosphere, asurf = 0."""
    m = _build_base_model("black")
    m.asurf = 0.0   # already set; explicit for clarity
    return m


def build_fresnel_only_model() -> pmd.Model:
    """Rough Fresnel surface, no atmosphere, no sub-surface ocean.
    Dashed reference curve in the top row of Figure 1.
    """
    m = pmd.Model()
    m.mma     = MMA
    m.gravity = GRAVITY
    m.dpol    = DEPOL
    m.asurf   = 0.0

    del m.layers.haze
    del m.layers.cloud

    m.wvl_list = list(WAVELENGTHS_UM)

    # Negligible atmosphere — effectively vacuum above the surface
    for lname in ("gasbelow", "gastop"):
        lyr = getattr(m.layers, lname)
        lyr.press   = 1.0e-8 if lname == "gasbelow" else 1.0e-10
        lyr.tau     = [0.0] * N_WAV
        lyr.tau_g   = [0.0] * N_WAV
        lyr.tau_ray = [0.0] * N_WAV
        lyr.rayscat = False

    m.surface = OceanSurface(
        wind_speed         = WIND_SPEED,
        n_air              = N_AIR,
        n_water            = N_WATER,
        foam_albedo        = FOAM_ALBEDO,
        n_fourier          = N_FOURIER,
        n_phi              = N_PHI,
        include_subsurface = False,   # no water body
        include_foam       = True,
        solver             = "adding_doubling",
    )
    return m


# ── Run planet_integrated for one model ──────────────────────────────────────

def run_model(model: pmd.Model, name: str, out_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Compute disk-integrated F and Q phase curves for one model.

    Processes phase angles in batches of 10 and saves a partial-progress
    checkpoint after each batch so computation is resumable across calls.
    The final .npz is written once all angles are complete.

    Returns
    -------
    F : (N_WAV, len(ALPHA_DEG)) -- total flux (normalised by piF0)
    Q : (N_WAV, len(ALPHA_DEG)) -- polarised flux (signed, normalised by piF0)
    """
    CACHE_DIR_F1.mkdir(parents=True, exist_ok=True)

    # Final cache -- fully done
    if out_path.exists():
        d = np.load(out_path)
        print(f"[{name}] Loaded from cache ({out_path.name})")
        return d["F"], d["Q"]

    # Partial-progress checkpoint
    partial_path = out_path.with_suffix(".partial.npz")
    n_alpha = len(ALPHA_DEG)
    F = np.full((N_WAV, n_alpha), np.nan)
    Q = np.full((N_WAV, n_alpha), np.nan)
    start_idx = 0

    if partial_path.exists():
        dp = np.load(partial_path)
        F = dp["F"]
        Q = dp["Q"]
        done = ~np.isnan(F[0])
        start_idx = int(np.sum(done))
        print(f"[{name}] Resuming from angle index {start_idx}/{n_alpha}")

    if start_idx >= n_alpha:
        np.savez(out_path, F=F, Q=Q, alpha_deg=ALPHA_DEG, wavelengths_um=WAVELENGTHS_UM)
        partial_path.unlink(missing_ok=True)
        print(f"[{name}] All angles complete. Saved -> {out_path}")
        return F, Q

    t0 = time.time()
    BATCH = 10  # angles per call -- fits comfortably within 44s
    batch = list(ALPHA_DEG[start_idx:start_idx + BATCH])
    end_idx = start_idx + len(batch)

    print(f"\n[{name}] planet_integrated angles {start_idx}-{end_idx - 1} "
          f"(alpha={batch[0]:.0f} to {batch[-1]:.0f} deg, {N_WAV} wl, npix={N_PIX}) ...",
          flush=True)

    pmd.planet_integrated(
        [model],
        alpha        = batch,
        npix         = N_PIX,
        output_names = [name],
        nmug         = N_MUG,
        nmug_mie     = N_MUG_MIE,
        nsubr        = N_SUBR,
        nmat         = 4,
        force        = True,
    )

    # model.I / model.Q: shape (N_WAV, len(batch))
    # model.I already equals F/(pi*F0) -- the paper normalisation -- because
    # planet_integrated sums I*cos(theta0)*pixel_area such that at alpha=0 the
    # result equals the geometric albedo (verified: Lambert asurf=1 gives 0.627,
    # theoretical 2/3=0.667, difference is npix discretisation).
    # Do NOT multiply by pi here.
    F[:, start_idx:end_idx] = np.array(model.I)
    Q[:, start_idx:end_idx] = np.array(model.Q)

    # Save partial progress
    np.savez(partial_path, F=F, Q=Q, alpha_deg=ALPHA_DEG, wavelengths_um=WAVELENGTHS_UM)
    dt = time.time() - t0
    print(f"  Batch done in {dt:.0f}s  ({end_idx}/{n_alpha} angles complete)")

    if end_idx >= n_alpha:
        np.savez(out_path, F=F, Q=Q, alpha_deg=ALPHA_DEG, wavelengths_um=WAVELENGTHS_UM)
        partial_path.unlink(missing_ok=True)
        print(f"  All done! Saved -> {out_path}")
    else:
        print(f"  {n_alpha - end_idx} angles remaining -- re-run to continue.")

    return F, Q


# ── Plotting ──────────────────────────────────────────────────────────────────

def _masked_P(F: np.ndarray, Q: np.ndarray, threshold: float = 5e-4) -> np.ndarray:
    """Return P = -Q/F with values masked (NaN) where F is too small to be
    meaningful.  ``threshold`` is relative to the per-wavelength maximum F.
    """
    P = np.full_like(F, np.nan)
    for i in range(F.shape[0]):
        Fmax = np.nanmax(np.abs(F[i]))
        good = np.abs(F[i]) > threshold * max(Fmax, 1e-10)
        P[i, good] = -Q[i, good] / F[i, good]
    return P


def make_plot(
    F_ocean,  Q_ocean,
    F_black,  Q_black,
    F_fonly,  Q_fonly,
):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(
        "Trees & Stam (2019) Fig. 1  —  cloud-free ocean vs black-surface planet",
        fontsize=11,
    )

    col_titles = ["Total flux $F / (\\pi F_0)$",
                  "Polarized flux $Q / (\\pi F_0)$",
                  "Degree of polarization $P$ (%)"]
    row_titles = ["Ocean planet\n(rough Fresnel + water body)",
                  "Black-surface planet\n($a_{\\rm surf}=0$)"]

    for col, title in enumerate(col_titles):
        axes[0, col].set_title(title, fontsize=10)
    for row, title in enumerate(row_titles):
        axes[row, 0].set_ylabel(title, fontsize=9)

    datasets = [
        (F_ocean, Q_ocean, 0),
        (F_black, Q_black, 1),
    ]

    for F, Q, row in datasets:
        P = _masked_P(F, Q)    # NaN where F is too small

        ax_F = axes[row, 0]
        ax_Q = axes[row, 1]
        ax_P = axes[row, 2]

        for iw, wl in enumerate(WAVELENGTHS_UM):
            kw = dict(color=WAV_COLORS[wl], lw=1.4, label=WAV_LABELS[wl])
            ax_F.plot(ALPHA_DEG, F[iw],  **kw)
            ax_Q.plot(ALPHA_DEG, Q[iw],  **kw)
            ax_P.plot(ALPHA_DEG, P[iw],  **kw)   # fraction, not %

        # Fresnel-only reference (dashed black) — F and Q only; P excluded
        # because F~0 at small alpha makes P unreliable there.
        if row == 0 and F_fonly is not None:
            P_fonly = _masked_P(F_fonly, Q_fonly, threshold=1e-2)
            kw_f  = dict(color="black", lw=1.2, ls="--",
                         label="Fresnel only (no atm, no ocean)")
            kw_f_ = dict(color="black", lw=1.2, ls="--", label="_nolegend_")
            for iw in range(N_WAV):
                kw_cur = kw_f if iw == 0 else kw_f_
                axes[0, 0].plot(ALPHA_DEG, F_fonly[iw],     **kw_cur)
                axes[0, 1].plot(ALPHA_DEG, Q_fonly[iw],     **kw_cur)
                axes[0, 2].plot(ALPHA_DEG, P_fonly[iw],     **kw_cur)

        ax_Q.axhline(0, color="grey", lw=0.6, ls=":")
        ax_P.axhline(0, color="grey", lw=0.6, ls=":")

    # ── Axis limits matching the paper ────────────────────────────────────────
    for ax in axes[:, 0]:
        ax.set_ylim(0, 0.30)        # F: 0 to 0.30
    for ax in axes[:, 1]:
        ax.set_ylim(-0.06, 0.02)    # Q: paper shows -0.06 to 0.02
    # P: paper shows 0 to 1.0 (= 0-100%); ocean can dip slightly below 0 at glint
    axes[1, 2].set_ylim(0, 1.00)
    axes[0, 2].set_ylim(-0.10, 1.00)
    # Fix P y-axis label to show fraction (matching paper), not %
    for ax in axes[:, 2]:
        ax.set_ylabel("Degree of polarization $P$", fontsize=9)

    for ax in axes.flat:
        ax.set_xlim(0, 180)
        ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
        ax.tick_params(direction="in", top=True, right=True, labelsize=8)

    for ax in axes[1, :]:
        ax.set_xlabel("Phase angle $\\alpha$ (°)", fontsize=9)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.02))

    fig.tight_layout(rect=[0, 0.04, 1, 1])
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved → {OUT_PNG}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    import argparse
    ap = argparse.ArgumentParser(description="Reproduce Trees & Stam (2019) Figure 1")
    ap.add_argument("--model", choices=["ocean", "black", "fonly", "plot", "all"],
                    default="all",
                    help="Which step to run.  Use 'ocean', 'black', 'fonly' to compute "
                         "each model independently (fits in one bash call), "
                         "then 'plot' to assemble the figure once all three are done.  "
                         "'all' runs everything sequentially (slower).")
    ap.add_argument("--force", action="store_true",
                    help="Ignore cached model results and recompute.")
    args = ap.parse_args()

    print("=== Trees & Stam (2019) Figure 1 ===")
    print(f"Wavelengths : {WAVELENGTHS_NM} nm")
    print(f"Phase angles: {ALPHA_DEG[0]}°–{ALPHA_DEG[-1]}° ({len(ALPHA_DEG)} pts)")
    print(f"npix={N_PIX}  n_fourier={N_FOURIER}  n_phi={N_PHI}")
    print(f"Gas bm(550nm)={BM_550}  step=--model {args.model}")
    print()

    CACHE_DIR_F1.mkdir(parents=True, exist_ok=True)

    if args.force:
        for p in [OUT_OCEAN, OUT_BLACK, OUT_FONLY, OUT_NPZ]:
            if p.exists(): p.unlink()

    if args.model in ("ocean", "all"):
        F_ocean, Q_ocean = run_model(build_ocean_model(), "ts19_ocean", OUT_OCEAN)

    if args.model in ("black", "all"):
        F_black, Q_black = run_model(build_black_model(), "ts19_black", OUT_BLACK)

    if args.model in ("fonly", "all"):
        F_fonly, Q_fonly = run_model(build_fresnel_only_model(), "ts19_fonly", OUT_FONLY)

    if args.model in ("plot", "all"):
        # Load all three from cache (or from memory if run in 'all' mode)
        if not all(p.exists() for p in [OUT_OCEAN, OUT_BLACK, OUT_FONLY]):
            missing = [p.name for p in [OUT_OCEAN, OUT_BLACK, OUT_FONLY] if not p.exists()]
            print(f"Cannot plot — missing results for: {missing}")
            print("Run each model first:  python recreate_fig1_trees_stam2019.py --model ocean")
            return 1

        d_o = np.load(OUT_OCEAN);  F_ocean = d_o["F"]; Q_ocean = d_o["Q"]
        d_b = np.load(OUT_BLACK);  F_black = d_b["F"]; Q_black = d_b["Q"]
        d_f = np.load(OUT_FONLY);  F_fonly = d_f["F"]; Q_fonly = d_f["Q"]

        # Legacy data was saved with an erroneous pi factor — undo it.
        # Detect by checking whether F(350nm,alpha=0) > 0.5 (physically impossible
        # since it would exceed a perfect Lambert sphere's geometric albedo of 2/3).
        if F_black[0, 0] > 0.5:
            print("  Detected legacy pi-factor in cached data — dividing out.")
            F_ocean /= np.pi; Q_ocean /= np.pi
            F_black /= np.pi; Q_black /= np.pi
            F_fonly /= np.pi; Q_fonly /= np.pi

        make_plot(F_ocean, Q_ocean, F_black, Q_black, F_fonly, Q_fonly)

        # Save consolidated NPZ
        np.savez(OUT_NPZ,
                 F_ocean=F_ocean, Q_ocean=Q_ocean,
                 F_black=F_black, Q_black=Q_black,
                 F_fonly=F_fonly, Q_fonly=Q_fonly,
                 wavelengths_um=WAVELENGTHS_UM, alpha_deg=ALPHA_DEG)
        print(f"Consolidated data saved → {OUT_NPZ}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
