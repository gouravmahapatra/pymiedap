#!/usr/bin/env python3
# This file is part of PyMieDAP, released under GNU General Public License.
"""
bench_nmug.py -- measure doubling-adding cost vs nmug and extrapolate.

Runs a single-wavelength, forward-peaked cloud (Mie r_eff=10 um + delta-M, a
stand-in for the Baum ice workload) at a series of `nmug` values, timing the
DAP step and recording peak memory. It then fits a power law and extrapolates
the per-wavelength wall-time to a target `nmug`, and prints the read_dap_output
`rfou` memory footprint for a target (nmuMAX, nfouMAX) build.

Use this to decide whether a target nmug (e.g. 400-500) is affordable on your
machine and how many wavelengths you can run in parallel.

    python examples/bench_nmug.py                         # default sweep
    python examples/bench_nmug.py --nmugs 50,100,150,200 --target 500
    python examples/bench_nmug.py --target 500 --nmu-max 512 --nfou-max 1024

Note: `nmug` cannot exceed the compiled `nmuMAX` (201 by default; raise it with
rebuild_highres_nmug.py). The extrapolation works from whatever points fit.
"""

import argparse
import os
import resource
import sys
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import pymiedap.pymiedap as pmd                       # noqa: E402
from pymiedap.tmatrix import delta_m_truncate         # noqa: E402


def _transparent():
    g = pmd.Aerosols(); g.typ = 'G'
    g.coefs = np.zeros((1, 4, 4, 1)); g.ncoefs = np.ones(1)
    g.ssalb = np.zeros(1); g.sext = np.zeros(1); g.ssca = np.zeros(1)
    g.col_dens = 0.0
    return g


def make_base_cloud(reff=8.0, wav=0.55, ngaur=240):
    """Compute the (expensive, nmug-independent) Mie cloud once; return its
    full expansion so each timed run only does delta-M + DAP."""
    a = pmd.Aerosols(nr=[1.33], ni=[0.0], r_eff=reff, v_eff=0.1, psd='2', typ='C')
    pmd.mie_code(a, [wav], ngaur=ngaur, nsubr=50)
    return dict(coefs=a.coefs.copy(), ncoefs=a.ncoefs.copy(),
                ssalb=a.ssalb.copy(), sext=a.sext.copy(), wav=wav)


def one_run(nmug, base, mfac=2.0):
    """Time one DAP solve at the given nmug, reusing the precomputed Mie cloud."""
    wav = base["wav"]
    a = pmd.Aerosols(typ='C')
    a.coefs = base["coefs"].copy(); a.ncoefs = base["ncoefs"].copy()
    a.ssalb = base["ssalb"].copy(); a.sext = base["sext"].copy()
    a.ssca = a.ssalb * a.sext; a.f = 1.0
    M = min(int(mfac * nmug) - 1, int(a.ncoefs[0]))
    ts = delta_m_truncate(a, M)
    m = pmd.Model(wvl_list=np.array([wav]))
    del m.layers.gasbelow, m.layers.haze
    g = m.layers.gastop
    g.rayscat = False; g.tau = [0.0]; g.tau_g = [0.0]; g.tau_ray = [0.0]
    m.layers.cloud.rayscat = False
    m.layers.cloud.tau = [5.0 * ts[0]]
    m.layers.cloud.tau_g = [0.0]; m.layers.cloud.tau_ray = [0.0]
    m.surface[0, 0] = 0.05
    m.layers.cloud.aerosols = a
    g.mixed_aerosols = _transparent()
    m.layers.cloud.mix_aerosols(); m.name = [""]
    t0 = time.time()
    pmd.dap_code(m, rename=True, output_name='bench_nmug', nmug=nmug, nmat=4)
    return time.time() - t0, M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmugs", default="40,80,120,160",
                    help="Comma-separated nmug values to time (<= compiled nmuMAX).")
    ap.add_argument("--target", type=int, default=500,
                    help="nmug to extrapolate the per-wavelength time to.")
    ap.add_argument("--nmu-max", type=int, default=512,
                    help="Target compiled nmuMAX for the rfou memory estimate.")
    ap.add_argument("--nfou-max", type=int, default=1024,
                    help="Target compiled nfouMAX for the rfou memory estimate.")
    args = ap.parse_args()

    nmugs = [int(x) for x in args.nmugs.split(",")]
    print("Computing base Mie cloud once (r_eff=8 um) ...")
    base = make_base_cloud()
    print("Timing DAP (delta-M M~2*nmug), one wavelength:\n")
    print("  {:>6}  {:>6}  {:>10}".format("nmug", "M", "DAP[s]"))
    ts, good = [], []
    for n in nmugs:
        try:
            dt, M = one_run(n, base)
            print("  {:6d}  {:6d}  {:10.1f}".format(n, M, dt))
            ts.append(dt); good.append(n)
        except Exception as e:
            print("  {:6d}   --     FAILED ({})".format(n, type(e).__name__))

    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # KB->GB on Linux
    print("\nPeak process memory so far: {:.2f} GB".format(peak_gb))

    if len(good) >= 2:
        # power-law fit  t = c * nmug^p   (log-log linear fit)
        p, logc = np.polyfit(np.log(good), np.log(ts), 1)
        c = np.exp(logc)
        t_target = c * args.target ** p
        print("\nPower-law fit: DAP time ~ nmug^{:.2f}".format(p))
        print("Extrapolated per-wavelength DAP time at nmug={}: {:.0f} s "
              "(~{:.1f} min)".format(args.target, t_target, t_target / 60))
    else:
        print("\nNeed >=2 successful points to extrapolate.")

    # rfou memory for the target build (read_dap_output allocation)
    rfou_gb = 4 * args.nmu_max * args.nmu_max * (args.nfou_max + 1) * 8 / 1e9
    print("\nread_dap_output rfou footprint at nmuMAX={}, nfouMAX={}: {:.1f} GB "
          "per process".format(args.nmu_max, args.nfou_max, rfou_gb))
    print("(Each parallel wavelength task needs about this much RAM, plus the "
          "DAP supermatrices ~ {:.2f} GB.)".format(
              (4 * args.nmu_max) ** 2 * 8 * 8 / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
