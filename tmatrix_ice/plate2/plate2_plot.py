#!/usr/bin/env python3
"""
Assemble the Plate 2 figure (Mishchenko, Travis & Mackowski, JQSRT 1998).

Degree of linear polarization for unpolarized incident light, -F21/F11 (%),
versus scattering angle (0-180 deg) and size parameter (0-30), at m=1.53+0.008i:

  (a) monodisperse spheres                        <- spher_p2a.f  -> panelA.dat
  (b) oblate spheroids eps=1.7, axis || beam      <- amplq_plate2 -> panelB.dat
  (c) oblate spheroids eps=1.7, axis perp. beam   <- amplq_plate2 -> panelC.dat
  (d) oblate spheroids eps=1.7, random orientation<- tmq_p2d      -> panelD.dat

Each panel*.dat file has three columns:  size_parameter  angle_deg  (-F21/F11)
(the third column is a fraction; multiplied by 100 below for percent).

Run from the folder holding the panel*.dat files:
    python3 plate2_plot.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(fn):
    d = np.loadtxt(fn)
    sizes = np.unique(d[:, 0])
    angles = np.unique(d[:, 1])
    Z = np.full((len(sizes), len(angles)), np.nan)
    si = {s: i for i, s in enumerate(sizes)}
    ai = {a: i for i, a in enumerate(angles)}
    for s, a, v in d:
        Z[si[s], ai[a]] = v * 100.0       # fraction -> percent
    return angles, sizes, Z


def main():
    panels = [("panelA.dat", "(a) SPHERES"),
              ("panelB.dat", "(b) SPHEROIDS, axis PARALLEL to beam"),
              ("panelC.dat", "(c) SPHEROIDS, axis PERPENDICULAR to beam"),
              ("panelD.dat", "(d) SPHEROIDS, RANDOM ORIENTATION")]
    fig, axs = plt.subplots(2, 2, figsize=(13, 10), sharex=True, sharey=True)
    fig.suptitle("Plate 2 replication - degree of linear polarization  "
                 r"$-F_{21}/F_{11}$ (%)" + "\n"
                 "m = 1.53 + 0.008i;  spheroids: oblate, aspect ratio 1.7",
                 fontsize=13, fontweight="bold")
    levels = np.arange(-100, 100.01, 10)
    im = None
    for ax, (fn, title) in zip(axs.ravel(), panels):
        ang, siz, Z = load(fn)
        im = ax.contourf(ang, siz, Z, levels=levels, cmap="RdBu_r", extend="both")
        ax.set_title(title, fontsize=11)
        ax.set_xlim(0, 180); ax.set_xticks(range(0, 181, 30)); ax.set_ylim(0, 30)
    for ax in axs[1, :]:
        ax.set_xlabel("Scattering angle (deg)")
    for ax in axs[:, 0]:
        ax.set_ylabel("Size parameter")
    fig.tight_layout(rect=[0, 0, 0.9, 0.95])
    cax = fig.add_axes([0.92, 0.12, 0.02, 0.76])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Degree of linear polarization (%)")
    fig.savefig("plate2_replication.png", dpi=150)
    print("wrote plate2_replication.png")


if __name__ == "__main__":
    main()
