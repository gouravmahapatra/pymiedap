#!/usr/bin/env python3
"""
Assemble the Plate 3 and Plate 4 figures
(Mishchenko, Travis & Mackowski, JQSRT 60, 309-324, 1998).

Both are polydisperse diagrams over scattering angle x EFFECTIVE size parameter,
power-law size distribution with v_eff = 0.1, m = 1.53 + 0.008i.

Plate 3: degree of linear polarization -F21/F11 (%) for randomly oriented
         polydisperse spheroids, six shapes (prolate/oblate x eps=1.4,1.7,2).
         Data: plate3_<tag>.dat  (cols: x_eff  angle  -F21/F11)
Plate 4: scattering-matrix elements for polydisperse spheres.
         Data: plate4.dat (cols: x_eff angle F11 F33/F11 F12/F11 F34/F11)

Run from the folder holding the .dat files:
    python3 plate34_plot.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def grid3(fn):
    d = np.loadtxt(fn)
    s = np.unique(d[:, 0]); a = np.unique(d[:, 1])
    Z = np.full((len(s), len(a)), np.nan)
    si = {v: i for i, v in enumerate(s)}; ai = {v: i for i, v in enumerate(a)}
    for sv, av, vv in d:
        Z[si[sv], ai[av]] = vv * 100
    return a, s, Z


def plate3():
    shapes = [('pro14', r'PROLATE, $\epsilon$=1.4'), ('obl14', r'OBLATE, $\epsilon$=1.4'),
              ('pro17', r'PROLATE, $\epsilon$=1.7'), ('obl17', r'OBLATE, $\epsilon$=1.7'),
              ('pro20', r'PROLATE, $\epsilon$=2'),   ('obl20', r'OBLATE, $\epsilon$=2')]
    fig, axs = plt.subplots(3, 2, figsize=(11, 12), sharex=True)
    fig.suptitle('Plate 3 replication - degree of linear polarization '
                 r'$-F_{21}/F_{11}$ (%)' + '\n'
                 'randomly oriented polydisperse spheroids, m=1.53+0.008i, '
                 r'$v_{eff}$=0.1', fontsize=12, fontweight='bold')
    lev = np.arange(-100, 100.01, 10)
    im = None
    for ax, (tag, title) in zip(axs.ravel(), shapes):
        a, s, Z = grid3(f'plate3_{tag}.dat')
        im = ax.contourf(a, s, Z, levels=lev, cmap='RdBu_r', extend='both')
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0, 180); ax.set_xticks(range(0, 181, 30)); ax.set_ylim(0, 30)
    for ax in axs[2, :]:
        ax.set_xlabel('Scattering angle (deg)')
    for ax in axs[:, 0]:
        ax.set_ylabel('Effective size parameter')
    fig.tight_layout(rect=[0, 0, 0.9, 0.95])
    cax = fig.add_axes([0.92, 0.12, 0.02, 0.76])
    plt.colorbar(im, cax=cax).set_label('DLP (%)')
    fig.savefig('plate3_replication.png', dpi=150)
    print('wrote plate3_replication.png')


def plate4():
    d = np.loadtxt('plate4.dat')
    s = np.unique(d[:, 0]); a = np.unique(d[:, 1])
    si = {v: i for i, v in enumerate(s)}; ai = {v: i for i, v in enumerate(a)}
    F11 = np.full((len(s), len(a)), np.nan); F33 = F11.copy(); F12 = F11.copy(); F34 = F11.copy()
    for r in d:
        i, j = si[r[0]], ai[r[1]]
        F11[i, j], F33[i, j], F12[i, j], F34[i, j] = r[2], r[3] * 100, r[4] * 100, r[5] * 100
    fig, axs = plt.subplots(2, 3, figsize=(15, 9), sharex=True)
    fig.suptitle('Plate 4 replication - scattering-matrix elements, polydisperse '
                 'spheres, m=1.53+0.008i, ' + r'$v_{eff}$=0.1', fontsize=12, fontweight='bold')
    lev = np.arange(-100, 100.01, 10)

    def panel(ax, Z, title, log=False):
        if log:
            im = ax.contourf(a, s, np.log10(np.clip(Z, 1e-3, None)),
                             levels=np.linspace(-1, 3, 17), cmap='viridis', extend='both')
        else:
            im = ax.contourf(a, s, Z, levels=lev, cmap='RdBu_r', extend='both')
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0, 180); ax.set_xticks(range(0, 181, 30))
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    panel(axs[0, 0], F11, r'$\log_{10}(F_{11})$', log=True)
    panel(axs[0, 1], np.full_like(F11, 100), r'$F_{22}/F_{11}$ (%)')
    panel(axs[0, 2], F33, r'$F_{33}/F_{11}$ (%)')
    panel(axs[1, 0], F33, r'$F_{44}/F_{11}$ (%)')   # spheres: F44=F33
    panel(axs[1, 1], F12, r'$F_{12}/F_{11}$ (%)')
    panel(axs[1, 2], F34, r'$F_{34}/F_{11}$ (%)')
    for ax in axs[1, :]:
        ax.set_xlabel('Scattering angle (deg)')
    for ax in axs[:, 0]:
        ax.set_ylabel('Effective size parameter')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig('plate4_replication.png', dpi=150)
    print('wrote plate4_replication.png')


if __name__ == '__main__':
    plate3()
    plate4()
