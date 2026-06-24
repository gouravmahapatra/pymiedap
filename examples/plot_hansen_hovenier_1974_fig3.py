#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pymiedap.pymiedap as pmd


WAVELENGTH_UM = 1.0
REFRACTIVE_INDEX_REAL = 1.44
REFRACTIVE_INDEX_IMAG = 0.0
EFFECTIVE_VARIANCES = (0.05, 0.07, 0.15)
CONTOUR_LEVELS = np.array(
    [-10, -5, -2, -1, 0, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80],
    dtype=float,
)


@contextmanager
def working_directory(path: Path):
    original = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original)


def compute_single_scattering_polarization(
    xeff_grid: np.ndarray,
    phase_grid_deg: np.ndarray,
    veff: float,
    nsubr: int,
    ngaur: int,
) -> np.ndarray:
    polarization_grid = np.zeros((xeff_grid.size, phase_grid_deg.size), dtype=float)

    with tempfile.TemporaryDirectory(prefix="hh74_fig4_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        with working_directory(tmpdir_path):
            for row, xeff in enumerate(xeff_grid):
                reff = xeff * WAVELENGTH_UM / (2.0 * np.pi)
                aerosol = pmd.Aerosols(
                    nr=[REFRACTIVE_INDEX_REAL],
                    ni=[REFRACTIVE_INDEX_IMAG],
                    r_eff=reff,
                    v_eff=veff,
                    psd="2",
                    typ="HH74",
                )

                pmd.mie_code(
                    aerosol,
                    [WAVELENGTH_UM],
                    output=True,
                    nsubr=nsubr,
                    ngaur=ngaur,
                )

                scattering_file = tmpdir_path / f"{aerosol.typ}.sc.{WAVELENGTH_UM:06.3f}"
                theta_deg, pl = pmd.read_mie_output(
                    os.fspath(scattering_file),
                    full_output=False,
                    nameout=os.fspath(tmpdir_path / "mie_output.dat"),
                )

                phase_deg = 180.0 - theta_deg
                sort_idx = np.argsort(phase_deg)
                phase_sorted = phase_deg[sort_idx]
                pl_sorted = 100.0 * pl[sort_idx]
                polarization_grid[row, :] = np.interp(phase_grid_deg, phase_sorted, pl_sorted)

    return polarization_grid


def make_plot(
    xeff_grid: np.ndarray,
    phase_grid_deg: np.ndarray,
    panels: list[np.ndarray],
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(len(EFFECTIVE_VARIANCES), 1, figsize=(9, 12), sharex=True)

    if len(EFFECTIVE_VARIANCES) == 1:
        axes = [axes]

    radius_grid = xeff_grid * WAVELENGTH_UM / (2.0 * np.pi)
    phase_mesh, xeff_mesh = np.meshgrid(phase_grid_deg, xeff_grid)

    for idx, (ax, veff, polarization) in enumerate(zip(axes, EFFECTIVE_VARIANCES, panels)):
        contours = ax.contour(
            phase_mesh,
            xeff_mesh,
            polarization,
            levels=CONTOUR_LEVELS,
            colors="black",
            linewidths=0.9,
        )
        ax.clabel(contours, inline=True, fontsize=7, fmt="%g")
        ax.set_ylabel(r"$x_{\mathrm{eff}} = 2\pi r_{\mathrm{eff}} / \lambda$")
        ax.set_title(
            f"({chr(ord('a') + idx)}) n = {REFRACTIVE_INDEX_REAL:.2f} + {REFRACTIVE_INDEX_IMAG:.0f}i, "
            f"v_eff = {veff:.2f}, lambda = {WAVELENGTH_UM:.1f} um",
            loc="left",
        )
        ax.grid(True, alpha=0.18)

        secx = ax.secondary_xaxis("top", functions=(lambda p: 180.0 - p, lambda s: 180.0 - s))
        secx.set_xlabel("Scattering Angle [deg]")

        secy = ax.secondary_yaxis(
            "right",
            functions=(
                lambda xeff: xeff * WAVELENGTH_UM / (2.0 * np.pi),
                lambda reff: reff * (2.0 * np.pi) / WAVELENGTH_UM,
            ),
        )
        secy.set_ylabel(r"$r_{\mathrm{eff}}$ [um] at $\lambda=1$ um")
        secy.set_yticks(np.arange(0.0, max(radius_grid) + 0.5, 0.5))

    axes[-1].set_xlabel("Phase Angle [deg]")
    axes[-1].set_xlim(0.0, 180.0)
    axes[-1].set_ylim(xeff_grid.min(), xeff_grid.max())
    fig.suptitle("Hansen & Hovenier (1974) Figure 4 Reproduction", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce the Figure 4 single-scattering polarization contours from Hansen & Hovenier (1974)."
    )
    parser.add_argument("--nradii", type=int, default=28, help="Number of effective-size samples (default: 28).")
    parser.add_argument("--ngaur", type=int, default=120, help="Number of Gaussian angles in Mie runs (default: 120).")
    parser.add_argument("--nsubr", type=int, default=40, help="Number of radius subintervals in Mie runs (default: 40).")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    phase_grid_deg = np.linspace(0.0, 180.0, 181)
    xeff_grid = np.linspace(0.5, 35.0, args.nradii)

    panels = []
    for veff in EFFECTIVE_VARIANCES:
        print(f"Computing v_eff={veff:.2f} ...")
        panel = compute_single_scattering_polarization(
            xeff_grid=xeff_grid,
            phase_grid_deg=phase_grid_deg,
            veff=veff,
            nsubr=args.nsubr,
            ngaur=args.ngaur,
        )
        panels.append(panel)

    output_path = Path(__file__).with_name("hansen_hovenier_1974_fig4.png")
    make_plot(xeff_grid=xeff_grid, phase_grid_deg=phase_grid_deg, panels=panels, output_path=output_path)

    print(f"Saved figure to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
