#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pymiedap.pymiedap as pmd


def main() -> int:
    model = pmd.Model()
    model.wvl_list = [0.7]
    del model.layers.gastop
    del model.layers.haze
    del model.layers.cloud
    model.layers.gasbelow.rayscat = False
    model.layers.gasbelow.tau = [0.0]
    model.surface[0, 0] = 1.0

    alphas = np.linspace(0.0, np.pi, 80)
    alphas_deg = np.degrees(alphas)
    theta = np.pi - alphas
    lambert_phase = 2 * (np.sin(theta) - theta * np.cos(theta)) / (3.0 * np.pi)

    pmd.planet_integrated([model], npix=60, alpha=alphas_deg, output_names=["plotLambert"])

    intensity = model.I[0, :]
    residual = intensity - lambert_phase

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    axes[0].plot(alphas_deg, intensity, label="PyMieDAP", linewidth=2)
    axes[0].plot(alphas_deg, lambert_phase, label="Analytical Lambert", linewidth=2, linestyle="--")
    axes[0].set_ylabel("Disk-Integrated I")
    axes[0].set_title("Lambertian Phase Curve Benchmark")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(alphas_deg, residual, color="black", linewidth=1.5)
    axes[1].axhline(0.0, color="tab:red", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Phase Angle [deg]")
    axes[1].set_ylabel("Residual")
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()

    output_path = Path(__file__).with_name("lambert_phase_curve.png")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

    print(f"Saved figure to {output_path}")
    print(f"Max absolute residual: {np.max(np.abs(residual)):.8e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
