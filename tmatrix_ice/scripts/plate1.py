#!/usr/bin/env python3
"""
Reproduce the Plate 1 figures from Mishchenko, Travis & Mackowski (JQSRT 1996/1998).

Inputs (produced by scripts/build_and_run.sh, expected in the current directory):
    spher.print               - Mie single sphere, x=5, m=1.5+0.005i
    bisphere_d2r.print        - bispheres, d = 2r (touching)
    bisphere_d2p5r.print      - bispheres, d = 2.5r
    bisphere_d4r.print        - bispheres, d = 4r
    bisphere_d8r.print        - bispheres, d = 8r

Outputs:
    plate1_mie_reference.png  - the single-sphere "black curve" only (6 panels)
    plate1_mie_reference.csv  - per-angle values for the sphere
    plate1_full_replication.png - full Plate 1: 4 bispheres + sphere reference

Optional verification (skipped automatically if miepython is not installed):
    cross-checks spher.f against miepython for x=5, m=1.5+0.005i.

Requirements: numpy, matplotlib  (miepython optional, for the cross-check).
Run from the folder holding the *.print files:
    python3 scripts/plate1.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# parsers for the Mishchenko-format scattering-matrix tables
# ----------------------------------------------------------------------
def _find_table(lines, needed):
    for i, l in enumerate(lines):
        if all(tok in l for tok in needed):
            return i + 1
    raise RuntimeError("scattering-matrix header not found")


def parse_spher(fn):
    """spher.f MATR prints: angle  F11  F33  F12  F34  (spheres: F22=F11, F44=F33)."""
    lines = open(fn).read().splitlines()
    s = _find_table(lines, ("F11", "F33", "F12", "F34"))
    rows = []
    for l in lines[s:]:
        p = l.split()
        if len(p) >= 5:
            try:
                rows.append(list(map(float, p[:5])))
            except ValueError:
                break
        else:
            break
    a = np.array(rows)  # ang F11 F33 F12 F34
    ang, F11, F33, F12, F34 = a.T
    return dict(ang=ang, F11=F11, F22=F11, F33=F33, F44=F33, F12=F12, F34=F34)


def parse_bisphere(fn):
    """bisphere.f MATR prints: angle  F11  F22  F33  F44  F12  F34."""
    lines = open(fn).read().splitlines()
    s = _find_table(lines, ("F11", "F22", "F44"))
    rows = []
    for l in lines[s:]:
        p = l.split()
        if len(p) >= 7:
            try:
                rows.append(list(map(float, p[:7])))
            except ValueError:
                break
        else:
            break
    a = np.array(rows)  # ang F11 F22 F33 F44 F12 F34
    ang, F11, F22, F33, F44, F12, F34 = a.T
    return dict(ang=ang, F11=F11, F22=F22, F33=F33, F44=F44, F12=F12, F34=F34)


def ratios(d):
    """Return the six panel quantities: F11 and the five % ratios to F11."""
    F11 = d["F11"]
    return (d["ang"], F11,
            100 * d["F22"] / F11, 100 * d["F33"] / F11,
            100 * d["F44"] / F11, 100 * d["F12"] / F11, 100 * d["F34"] / F11)


# ----------------------------------------------------------------------
# optional independent cross-check of spher.f
# ----------------------------------------------------------------------
def verify_with_miepython(mie):
    try:
        import miepython
    except ImportError:
        print("[verify] miepython not installed -- skipping cross-check")
        return
    x, m = 5.0, 1.5 - 0.005j     # miepython uses n - k i
    qext, qsca, qback, g = miepython.efficiencies_mx(m, x)
    print(f"[verify] miepython:  Qext={qext:.5f}  Qsca={qsca:.5f}  "
          f"albedo={qsca/qext:.5f}  g={g:.5f}")
    mu = np.cos(np.radians(mie["ang"]))
    S1, S2 = miepython.S1_S2(m, x, mu, norm="albedo")
    S1, S2 = np.asarray(S1), np.asarray(S2)
    F11 = 0.5 * (np.abs(S1) ** 2 + np.abs(S2) ** 2)
    F12 = 0.5 * (np.abs(S2) ** 2 - np.abs(S1) ** 2)
    F33 = (S2 * np.conj(S1)).real
    F34 = (S1 * np.conj(S2)).imag
    for name, sp, mp in [("F33/F11", mie["F33"] / mie["F11"], F33 / F11),
                         ("F12/F11", mie["F12"] / mie["F11"], F12 / F11),
                         ("F34/F11", mie["F34"] / mie["F11"], F34 / F11)]:
        print(f"[verify] max|d {name}| (spher.f vs miepython) = "
              f"{np.nanmax(np.abs(sp - mp)):.2e}")


# ----------------------------------------------------------------------
# plotting
# ----------------------------------------------------------------------
PANELS = [("F11 (phase function)", "F11", "log"),
          ("F22/F11 (%)", "F22r", None), ("F33/F11 (%)", "F33r", None),
          ("F44/F11 (%)", "F44r", None), ("F12/F11 (%)", "F12r", None),
          ("F34/F11 (%)", "F34r", None)]


def _unpack(d):
    ang, F11, F22r, F33r, F44r, F12r, F34r = ratios(d)
    return dict(ang=ang, F11=F11, F22r=F22r, F33r=F33r,
                F44r=F44r, F12r=F12r, F34r=F34r)


def plot_reference(mie, out_png, out_csv):
    p = _unpack(mie)
    fig, ax = plt.subplots(2, 3, figsize=(13, 7.5))
    fig.suptitle("Plate 1 black reference curve — single Mie sphere, "
                 "x = 5, m = 1.5 + 0.005i", fontsize=13, fontweight="bold")
    for a, (title, key, scale) in zip(ax.ravel(), PANELS):
        (a.semilogy if scale == "log" else a.plot)(p["ang"], p[key], "k-", lw=1.8)
        a.set_title(title, fontsize=11)
        a.set_xlim(0, 180); a.set_xticks(range(0, 181, 30))
        a.grid(alpha=0.3); a.set_xlabel("Scattering angle (deg)")
        if scale != "log":
            a.axhline(0, color="gray", lw=0.6)
    ax[0, 1].set_ylim(99, 101)
    for c in (2, 3, 4, 5):
        ax.ravel()[c].set_ylim(-100, 100)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=150)
    print("wrote", out_png)

    out = np.column_stack([p["ang"], p["F11"], p["F22r"], p["F33r"],
                           p["F44r"], p["F12r"], p["F34r"]])
    np.savetxt(out_csv, out, delimiter=",", fmt="%.6g",
               header="angle_deg,F11,F22_over_F11_pct,F33_over_F11_pct,"
                      "F44_over_F11_pct,F12_over_F11_pct,F34_over_F11_pct",
               comments="")
    print("wrote", out_csv)


def plot_full(bispheres, mie, out_png):
    styles = {"d2r": ("d = 2r (touching)", "#d62728"),
              "d2p5r": ("d = 2.5r", "#ff7f0e"),
              "d4r": ("d = 4r", "#2ca02c"),
              "d8r": ("d = 8r", "#1f77b4")}
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle("Plate 1 (Mishchenko et al. 1996, JQSRT) — replicated\n"
                 "Randomly oriented bispheres, component x = 5, m = 1.5 + 0.005i; "
                 "black = single Mie sphere", fontsize=12, fontweight="bold")
    pm = _unpack(mie)
    for a, (title, key, scale) in zip(ax.ravel(), PANELS):
        for tag, (lbl, c) in styles.items():
            p = _unpack(bispheres[tag])
            (a.semilogy if scale == "log" else a.plot)(
                p["ang"], p[key], color=c, lw=1.6, label=lbl)
        (a.semilogy if scale == "log" else a.plot)(
            pm["ang"], pm[key], "k-", lw=2.0, label="single sphere")
        a.set_title(title, fontsize=11)
        a.set_xlim(0, 180); a.set_xticks(range(0, 181, 30))
        a.grid(alpha=0.3); a.set_xlabel("Scattering angle (deg)")
        if scale != "log":
            a.axhline(0, color="gray", lw=0.6)
    ax[0, 1].set_ylim(60, 102)
    for c in (2, 3, 4, 5):
        ax.ravel()[c].set_ylim(-100, 100)
    ax.ravel()[0].legend(fontsize=8, loc="upper right")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_png, dpi=150)
    print("wrote", out_png)


def main():
    mie = parse_spher("spher.print")
    verify_with_miepython(mie)
    plot_reference(mie, "plate1_mie_reference.png", "plate1_mie_reference.csv")

    tags = {"d2r": "bisphere_d2r.print", "d2p5r": "bisphere_d2p5r.print",
            "d4r": "bisphere_d4r.print", "d8r": "bisphere_d8r.print"}
    if all(os.path.exists(f) for f in tags.values()):
        bispheres = {k: parse_bisphere(v) for k, v in tags.items()}
        plot_full(bispheres, mie, "plate1_full_replication.png")
    else:
        print("[full] bisphere_*.print not all present -- "
              "skipping full Plate 1 figure")


if __name__ == "__main__":
    main()
