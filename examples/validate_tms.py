#!/usr/bin/env python3
# This file is part of PyMieDAP, released under GNU General Public License.
"""
validate_tms.py -- validate the TMS single-scattering correction.

Checks that an AGGRESSIVE delta-M solve (small M, small nmug) PLUS the TMS
single-scattering correction (pymiedap.sscorr.tms_correct_local) recovers the
result of a fully CONVERGED doubling-adding solve -- and does so far more
cheaply.

This needs a converged reference, which for a forward-peaked cloud requires
nmug >~ ncoef/2 (~100+). That is minutes on a fast workstation/cluster but
exceeds a constrained sandbox, which is why this validation lives here rather
than being run inline.

    python examples/validate_tms.py --reff 4 --tau 1.0 \
        --nmug-ref 150 --m-cheap 24 --nmug-cheap 16

Interpretation: TMS error (vs the converged reference) should be a few percent
and clearly smaller than the bare delta-M error, across geometries and in both
I and Q. If so, the correction is working and you can run the Baum ice at a
small, affordable nmug with TMS instead of brute-forcing nmug=500.
"""

import argparse
import os
import sys
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import module_readmie as rd                    # noqa: E402
import pymiedap.pymiedap as pmd                # noqa: E402
from pymiedap.tmatrix import delta_m_truncate  # noqa: E402
from pymiedap.sscorr import tms_correct_local  # noqa: E402


def _transparent():
    t = pmd.Aerosols(); t.typ = 'G'
    t.coefs = np.zeros((1, 4, 4, 1)); t.ncoefs = np.ones(1)
    t.ssalb = np.zeros(1); t.sext = np.zeros(1); t.ssca = np.zeros(1)
    t.col_dens = 0.0
    return t


def _build(aero, tauv, nmug, name):
    m = pmd.Model(wvl_list=np.array([0.55]))
    del m.layers.gasbelow, m.layers.haze
    g = m.layers.gastop
    g.rayscat = False; g.tau = [0.]; g.tau_g = [0.]; g.tau_ray = [0.]
    m.layers.cloud.rayscat = False
    m.layers.cloud.tau = [tauv]; m.layers.cloud.tau_g = [0.]; m.layers.cloud.tau_ray = [0.]
    m.surface[0, 0] = 0.0
    m.layers.cloud.aerosols = aero
    m.layers.cloud.mix_aerosols()
    g.mixed_aerosols = _transparent()
    m.name = [""]
    pmd.dap_code(m, rename=True, output_name=name, nmug=nmug, nmat=4)
    return m.name[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reff", type=float, default=4.0)
    ap.add_argument("--wav", type=float, default=0.55)
    ap.add_argument("--tau", type=float, default=1.0)
    ap.add_argument("--ngaur", type=int, default=200, help="Mie angle count.")
    ap.add_argument("--nmug-ref", type=int, default=150,
                    help="nmug for the converged reference (>= ncoef/2).")
    ap.add_argument("--m-cheap", type=int, default=24, help="delta-M order.")
    ap.add_argument("--nmug-cheap", type=int, default=16, help="cheap nmug.")
    args = ap.parse_args()

    a = pmd.Aerosols(nr=[1.33], ni=[0.0], r_eff=args.reff, v_eff=0.1,
                     psd='2', typ='C')
    pmd.mie_code(a, [args.wav], ngaur=args.ngaur, nsubr=50)
    full = a.coefs.copy(); ncf = int(a.ncoefs[0]); ssalb = float(a.ssalb[0])
    g = full[0, 0, 0, 1] / 3
    print("reff=%.1f  ncoef=%d  g=%.3f  ssalb=%.5f" % (args.reff, ncf, g, ssalb))
    if args.nmug_ref < ncf / 2:
        print("WARNING: nmug-ref=%d < ncoef/2=%d -> reference NOT converged; "
              "raise --nmug-ref." % (args.nmug_ref, ncf // 2))

    t0 = time.time()
    fn_ref = _build(a, args.tau, args.nmug_ref, 'tms_ref')
    t_ref = time.time() - t0

    ac = pmd.Aerosols(typ='C')
    ac.coefs = full.copy(); ac.ncoefs = np.array([ncf], float)
    ac.ssalb = a.ssalb.copy(); ac.sext = np.ones(1); ac.ssca = np.ones(1); ac.f = 1.0
    ts = delta_m_truncate(ac, args.m_cheap)
    tau_dm = args.tau * ts[0]; ssalb_dm = float(ac.ssalb[0])
    trc = np.asarray(ac.coefs[0])
    t0 = time.time()
    fn_cheap = _build(ac, tau_dm, args.nmug_cheap, 'tms_cheap')
    t_cheap = time.time() - t0
    print("converged ref: nmug=%d (%.1fs)   cheap delta-M: M=%d nmug=%d (%.1fs)  "
          "speedup x%.0f" % (args.nmug_ref, t_ref, args.m_cheap, args.nmug_cheap,
                             t_cheap, t_ref / max(t_cheap, 1e-6)))

    # full tabulated phase matrix for the exact single-scatter term
    thg = np.linspace(0, 180, 721)
    F11f = np.zeros_like(thg); F12f = np.zeros_like(thg)
    cf = np.zeros((4, 4, 1001), order='F'); nn = min(ncf, 1000)
    cf[:, :, :nn + 1] = full[0, :, :, :nn + 1]; f = np.zeros(6)
    for i, t in enumerate(thg):
        rd.expand(nn, cf, float(t), f); F11f[i] = f[0]; F12f[i] = f[4]

    print("\n sza emi dphi |  errI bare  errI tms  |  errQ bare  errQ tms")
    geos = [(40, 30, 40), (55, 35, 80), (35, 55, 120), (50, 50, 160),
            (60, 25, 100), (45, 45, 60)]
    aggE = {"Ib": [], "It": [], "Qb": [], "Qt": []}
    for sza, emi, dphi in geos:
        Ir, Qr, Ur, _ = pmd.read_dap_output(np.array([0.]), np.array([sza]),
                                            np.array([emi]), fn_ref,
                                            phi=np.array([dphi]), beta=np.array([0.]))
        Ib, Qb, Ub, _ = pmd.read_dap_output(np.array([0.]), np.array([sza]),
                                            np.array([emi]), fn_cheap,
                                            phi=np.array([dphi]), beta=np.array([0.]))
        Ic, Qc, Uc = tms_correct_local(Ib[0], Qb[0], Ub[0], sza, emi, dphi,
                                       trc, args.m_cheap, tau_dm, ssalb_dm,
                                       thg, F11f, F12f, args.tau, ssalb)
        eIb = abs(Ib[0] - Ir[0]) / abs(Ir[0]) * 100
        eIt = abs(Ic - Ir[0]) / abs(Ir[0]) * 100
        dQ = max(abs(Qr[0]), 1e-6)
        eQb = abs(Qb[0] - Qr[0]) / dQ * 100
        eQt = abs(Qc - Qr[0]) / dQ * 100
        aggE["Ib"].append(eIb); aggE["It"].append(eIt)
        aggE["Qb"].append(eQb); aggE["Qt"].append(eQt)
        print(" %3d %3d %4d | %7.1f%% %7.1f%% | %7.1f%% %7.1f%%"
              % (sza, emi, dphi, eIb, eIt, eQb, eQt))

    print("\nmean |error|:  I bare %.1f%% -> tms %.1f%% ;  Q bare %.1f%% -> tms %.1f%%"
          % (np.mean(aggE["Ib"]), np.mean(aggE["It"]),
             np.mean(aggE["Qb"]), np.mean(aggE["Qt"])))
    ok = np.mean(aggE["It"]) < np.mean(aggE["Ib"]) and np.mean(aggE["It"]) < 5
    print("VERDICT:", "TMS improves on bare delta-M and is within a few %% -- "
          "correction validated." if ok else
          "TMS did NOT clearly beat bare delta-M -- inspect (is the reference "
          "converged? is M_cheap too small?).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
