#!/usr/bin/env python3
"""Approximate Hansen & Hovenier (1974) Figure 4 with PyMieDAP.

This example recreates the paper's disk-integrated Venus polarization curves
at 0.55 um for five effective radii. The optical inputs follow the figure
caption: n_r = 1.44, b = 0.07, and a Rayleigh contribution f_r = 0.045.

The original paper uses a thick, homogeneous Venus atmosphere and its own
phase-matrix construction. In this repo, the closest self-contained proxy is a
single optically thick aerosol layer over a black surface, with user-specified
Rayleigh scattering scaled to 0.55 um. That means this script is intended as a
qualitative reproduction of the Figure 4 curve family, not a point-by-point
reconstruction of the published calculations.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

if TYPE_CHECKING:
    import pymiedap.pymiedap as pmd


WAVELENGTH_UM = 0.55
REFRACTIVE_INDEX_REAL = 1.44
EFFECTIVE_VARIANCE = 0.07
RAYLEIGH_FR_AT_0365 = 0.045
REFERENCE_WAVELENGTH_UM = 0.365
EFFECTIVE_RADII_UM = (0.6, 0.9, 1.05, 1.2, 1.5)
DEFAULT_PHASE_COUNT = 61
DEFAULT_CLOUD_TAU = 30.0
DEFAULT_NPIX = 18
DEFAULT_NMUG = 16
DEFAULT_NMUG_MIE = 16
DEFAULT_NSUBR = 20
FIGURE_PATH = REPO_ROOT / "examples" / "hansen_hovenier_1974_fig4.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Approximate Figure 4 of Hansen & Hovenier (1974) with "
            "disk-integrated PyMieDAP runs."
        )
    )
    parser.add_argument("--phase-count", type=int, default=DEFAULT_PHASE_COUNT)
    parser.add_argument("--npix", type=int, default=DEFAULT_NPIX)
    parser.add_argument("--nmug", type=int, default=DEFAULT_NMUG)
    parser.add_argument("--nmug-mie", type=int, default=DEFAULT_NMUG_MIE)
    parser.add_argument("--nsubr", type=int, default=DEFAULT_NSUBR)
    parser.add_argument("--cloud-tau", type=float, default=DEFAULT_CLOUD_TAU)
    parser.add_argument("--curve-only", action="store_true")
    parser.add_argument("--radius", type=float, default=None)
    parser.add_argument("--curve-output", type=Path, default=None)
    parser.add_argument(
        "--output",
        type=Path,
        default=FIGURE_PATH,
        help="Path of the PNG figure to write.",
    )
    return parser.parse_args()


def rayleigh_ratio_at_wavelength(
    wavelength_um: float,
    fr_at_reference: float = RAYLEIGH_FR_AT_0365,
    reference_wavelength_um: float = REFERENCE_WAVELENGTH_UM,
) -> float:
    return fr_at_reference * (reference_wavelength_um / wavelength_um) ** 4


def build_model(
    effective_radius_um: float,
    cloud_tau: float,
    rayleigh_ratio: float,
) -> "pmd.Model":
    import pymiedap.pymiedap as pmd

    model = pmd.Model()
    model.wvl_list = [WAVELENGTH_UM]
    del model.layers.gastop
    del model.layers.haze
    del model.layers.cloud

    # The paper uses a homogeneous, optically thick atmosphere. A single,
    # thick scattering layer over a black surface is a close PyMieDAP proxy.
    layer = model.layers.gasbelow
    layer.tau = [cloud_tau]
    layer.rayscat = False
    layer.tau_ray = [cloud_tau * rayleigh_ratio]
    layer.aerosols = pmd.Aerosols(
        nr=[REFRACTIVE_INDEX_REAL],
        ni=[0.0],
        r_eff=effective_radius_um,
        v_eff=EFFECTIVE_VARIANCE,
        psd="2",
    )
    model.asurf = 0.0
    return model


def compute_curve(
    effective_radius_um: float,
    phases_deg: np.ndarray,
    cloud_tau: float,
    npix: int,
    nmug: int,
    nmug_mie: int,
    nsubr: int,
    run_id: str,
) -> np.ndarray:
    import pymiedap.pymiedap as pmd

    rayleigh_ratio = rayleigh_ratio_at_wavelength(WAVELENGTH_UM)
    model = build_model(effective_radius_um, cloud_tau, rayleigh_ratio)
    tag = f"hh74f4{run_id}_{str(effective_radius_um).replace('.', '')}"
    pmd.planet_integrated(
        [model],
        alpha=phases_deg,
        npix=npix,
        output_names=[tag],
        nmug=nmug,
        nmug_mie=nmug_mie,
        nsubr=nsubr,
    )
    return 100.0 * model.P[0]


def style_axes(ax: plt.Axes) -> None:
    ax.set_xlim(0, 180)
    ax.set_ylim(-5.5, 7.0)
    ax.set_xticks(np.arange(0, 181, 20))
    ax.set_yticks(np.arange(-4, 7, 2))
    ax.set_xlabel("Phase Angle")
    ax.set_ylabel("% Polarization")
    ax.tick_params(direction="in", top=True, right=True)


def plot_curves(phases_deg: np.ndarray, curves: dict[float, np.ndarray], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    line_styles = {
        1.5: {"linestyle": (0, (9, 4)), "linewidth": 1.2},
        1.2: {"linestyle": (0, (1, 2)), "linewidth": 1.5},
        1.05: {"linestyle": "-", "linewidth": 1.4},
        0.9: {"linestyle": (0, (6, 3)), "linewidth": 1.2},
        0.6: {"linestyle": "-", "linewidth": 0.8, "alpha": 0.65},
    }

    fig, ax = plt.subplots(figsize=(6.2, 6.8), dpi=160)
    for radius in sorted(curves.keys(), reverse=True):
        ax.plot(
            phases_deg,
            curves[radius],
            color="black",
            label=f"a = {radius} um",
            **line_styles[radius],
        )

    style_axes(ax)
    ax.legend(loc="upper left", frameon=False)
    ax.text(
        0.97,
        0.94,
        "\n".join(
            (
                "lambda = 0.55 um",
                "n_r = 1.44",
                "b = 0.07",
            )
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox={"facecolor": "white", "edgecolor": "black", "boxstyle": "square,pad=0.35"},
    )
    ax.set_title("Hansen & Hovenier (1974) Figure 4 Approximation", pad=12)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def run_curve_subprocess(
    radius: float,
    args: argparse.Namespace,
) -> np.ndarray:
    with tempfile.TemporaryDirectory(prefix="hh74f4_") as tmpdir:
        curve_output = Path(tmpdir) / f"curve_{str(radius).replace('.', '')}.npz"
        cmd = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--curve-only",
            "--radius",
            str(radius),
            "--curve-output",
            str(curve_output),
            "--phase-count",
            str(args.phase_count),
            "--npix",
            str(args.npix),
            "--nmug",
            str(args.nmug),
            "--nmug-mie",
            str(args.nmug_mie),
            "--nsubr",
            str(args.nsubr),
            "--cloud-tau",
            str(args.cloud_tau),
        ]
        subprocess.run(cmd, check=True)
        data = np.load(curve_output)
        return data["polarization_percent"]


def main() -> int:
    args = parse_args()
    phases_deg = np.linspace(0.0, 180.0, args.phase_count)
    curves: dict[float, np.ndarray] = {}
    run_id = f"{os.getpid() % 1000:03d}{int(time.time()) % 10000:04d}"

    if args.curve_only:
        if args.radius is None or args.curve_output is None:
            raise ValueError("--curve-only requires --radius and --curve-output.")
        curve = compute_curve(
            effective_radius_um=args.radius,
            phases_deg=phases_deg,
            cloud_tau=args.cloud_tau,
            npix=args.npix,
            nmug=args.nmug,
            nmug_mie=args.nmug_mie,
            nsubr=args.nsubr,
            run_id=run_id,
        )
        args.curve_output.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            args.curve_output,
            phases_deg=phases_deg,
            polarization_percent=curve,
        )
        finite = np.isfinite(curve)
        if np.any(finite):
            print(
                f"a_eff = {args.radius:.2f} um "
                f"range = {np.nanmin(curve):.3f} to {np.nanmax(curve):.3f} %"
            )
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(0)

    print("Using PyMieDAP proxy for Hansen & Hovenier (1974) Figure 4")
    print(f"wavelength = {WAVELENGTH_UM:.3f} um")
    print(f"n_r = {REFRACTIVE_INDEX_REAL:.2f}")
    print(f"b = {EFFECTIVE_VARIANCE:.2f}")
    print(
        "effective Rayleigh ratio at 0.55 um = "
        f"{rayleigh_ratio_at_wavelength(WAVELENGTH_UM):.6f}"
    )

    for radius in EFFECTIVE_RADII_UM:
        print(f"Computing a_eff = {radius:.2f} um")
        curves[radius] = run_curve_subprocess(radius, args)
        finite = np.isfinite(curves[radius])
        if np.any(finite):
            print(
                f"  range = {np.nanmin(curves[radius]):.3f} to "
                f"{np.nanmax(curves[radius]):.3f} %"
            )

    plot_curves(phases_deg, curves, args.output)
    print(f"Saved figure to {args.output}")
    print(
        "Note: this reproduces the paper's thick homogeneous atmosphere "
        "qualitatively using a black-surface, optically thick PyMieDAP proxy."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
