#!/usr/bin/env python3
# This file is part of PyMieDAP, released under GNU General Public License.
"""
ice_thin_sensitivity.py
=============================================================================
Does a realistic, optically-thin Baum/Yang ice cloud actually change the
disk-integrated Earth-like (ocean + liquid water cloud) spectrum -- and does
the ACCURACY of the ice forward-peak treatment matter?

This is the decision experiment after concluding that TMS is not the production
route for D_eff~60 um ice (the doubling-adding NaNs at the small truncation
order TMS would need, and where it is stable bare delta-M already converges).

Cases (liquid cloud + surface identical in all three; only the ice changes):

    A : liquid water cloud only, NO ice
    B : liquid + ice, COARSE stable delta-M truncation of the ice (M_ICE_B)
    C : liquid + ice, FINER  stable delta-M truncation of the ice (M_ICE_C)

Comparison, per wavelength and phase angle:
    F = I  (disk-integrated reflectance)   Q   P = -Q/I   Pl = sqrt(Q^2+U^2)/I

Decision:
    If  B vs C  differ by  < ~1% in F, < ~2-3% in Q, < ~0.5 pp in P
        -> a cheap stable ice treatment is good enough; use it.
    If larger (esp. near cloudbow / glint / high phase)
        -> the ice forward peak matters; develop delta-fit truncation.
    If B and/or C come back NaN (the ice is too forward-peaked to truncate
        stably at that M) the script reports the minimum stable M it found --
        if even that is expensive, delta-fit is required.

Brute-force nmug=200-500 is deliberately NOT used.

Prerequisite: a Baum coefficient cache for the chosen D_eff, built with
    python examples/convert_baum_to_pymiedap.py --deff 60 --wlmin 0.2 --wlmax 2.0

Usage:
    python examples/ice_thin_sensitivity.py            # defaults below
    python examples/ice_thin_sensitivity.py --deff 60 --m-ice-b 24 --m-ice-c 60 \
        --m-liq 80 --nmug 50
"""

import argparse
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import pymiedap.pymiedap as pmd                       # noqa: E402
from pymiedap.tmatrix import delta_m_truncate          # noqa: E402
from pymiedap.baum import fill_aerosol_from_cache      # noqa: E402

# --- Earth-like configuration (Roccetti et al. 2025-ish, ocean scenario) ----
LIQ_REFF = 9.0;  LIQ_VEFF = 0.10;  LIQ_TAU = 5.44
LIQ_NR = 1.335;  LIQ_NI = 1.0e-8
ICE_TAU = 0.6
P0 = 1.013;  SCALE_H_KM = 8.0
P_GASTOP = 1e-3
P_ICE = P0 * np.exp(-3.88 / SCALE_H_KM)
P_LIQ = P0 * np.exp(-1.55 / SCALE_H_KM)
P_SURF = P0
DAP_DIR = os.path.join(REPO, "dap_database")


def _transparent(nwvl):
    g = pmd.Aerosols(); g.typ = 'G'
    g.coefs = np.zeros((1, 4, 4, 1)); g.ncoefs = np.ones(1)
    g.ssalb = np.zeros(1); g.sext = np.zeros(1); g.ssca = np.zeros(1)
    g.col_dens = 0.0
    return g


def build_and_compute(wvls, case, baum_cache, m_ice, m_liq, nmug, nmug_mie,
                      asurf, name):
    """Build the cloudy column for a case and run the DAP. Returns (model, ok)
    where ok=False if the disk reflectance comes back non-finite."""
    wvls = np.asarray(wvls, float)
    m = pmd.Model(wvl_list=wvls); m.asurf = asurf
    L = m.layers; nw = len(wvls)
    L.gastop.press = P_GASTOP; L.gastop.tau = [0.0]
    L.gasbelow.press = P_SURF; L.gasbelow.tau = [0.0]
    L.haze.press = P_ICE
    L.cloud.press = P_LIQ; L.cloud.tau = [LIQ_TAU]
    a = L.cloud.aerosols
    a.typ = 'C'; a.psd = '2'; a.r_eff = LIQ_REFF; a.v_eff = LIQ_VEFF
    a.nr = [LIQ_NR] * max(2, nw); a.ni = [LIQ_NI] * max(2, nw); a.layered = False
    L.haze.tau = [0.0 if case == 'A' else ICE_TAU]
    m.wvl_list = m.wvl_list

    for lay, layer in vars(m.layers).items():
        if hasattr(layer, 'mixed_aerosols'):
            del layer.mixed_aerosols
        if lay == 'cloud':                                   # liquid: Mie+delta-M
            pmd.mie_code(layer.aerosols, wvls, ngaur=nmug_mie, nsubr=50)
            ts = delta_m_truncate(layer.aerosols, m_liq)
            t = np.atleast_1d(layer.tau).astype(float)
            if t.size != nw:
                t = t[0] * np.ones(nw)
            layer.tau = t * ts
            layer.mix_aerosols()
        elif lay == 'haze':                                  # ice
            if case == 'A':
                layer.mixed_aerosols = _transparent(nw)      # no ice
            else:
                fill_aerosol_from_cache(layer.aerosols, baum_cache,
                                        wavelengths_um=wvls)
                layer.aerosols.f = 1.0
                ts = delta_m_truncate(layer.aerosols, m_ice)
                t = np.atleast_1d(layer.tau).astype(float)
                if t.size != nw:
                    t = t[0] * np.ones(nw)
                layer.tau = t * ts
                layer.mix_aerosols()
        else:                                                # gas
            for an, aero in vars(layer).items():
                if isinstance(aero, pmd.Aerosols):
                    pmd.mie_code(aero, wvls, ngaur=20, nsubr=50)
            layer.mix_aerosols()

    pmd.dap_code(m, rename=True, output_name=name, nmug=nmug, nmat=4,
                 path_output=DAP_DIR)
    return m


def build_clear(wvls, nmug, asurf, name):
    """Cloud-free companion column (gas + surface), computed once and reused.
    planet_integrated needs a [cloudy, clear] pair; single-model full_disk
    returns zero in this build."""
    wvls = np.asarray(wvls, float)
    m = pmd.Model(wvl_list=wvls); m.asurf = asurf
    L = m.layers
    L.gastop.press = P_GASTOP; L.gastop.tau = [0.0]
    L.haze.press = P_ICE; L.haze.tau = [0.0]
    L.cloud.press = P_LIQ; L.cloud.tau = [0.0]
    L.gasbelow.press = P_SURF; L.gasbelow.tau = [0.0]
    m.wvl_list = m.wvl_list
    pmd.compute_model(m, force=True, output_name=name, path_input=DAP_DIR,
                      nmug=nmug)
    return m


def disk(cloudy, clear, alphas, npix, cloud_cover, nmug):
    """Disk-integrate the cloudy column (mixed with the clear companion at the
    given cloud cover) at the requested phase angles."""
    cc = min(max(cloud_cover, 1e-3), 1.0 - 1e-3)
    pmd.planet_integrated([cloudy, clear], alpha=list(alphas), force=False,
                          rename=True, nmug=nmug, npix=npix, patchy=False,
                          fclouds=[cc, 1.0 - cc])
    return cloudy.I, cloudy.Q, cloudy.U   # shape (nwvl, nphase)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deff", type=float, default=60.0)
    ap.add_argument("--m-ice-b", type=int, default=24)
    ap.add_argument("--m-ice-c", type=int, default=60)
    ap.add_argument("--m-liq", type=int, default=80)
    ap.add_argument("--nmug", type=int, default=50)
    ap.add_argument("--nmug-mie", type=int, default=120)
    ap.add_argument("--npix", type=int, default=15)
    ap.add_argument("--wvls", default="0.40,0.55,0.70,0.85,1.00")
    ap.add_argument("--phases", default="60,90,120,150")
    ap.add_argument("--cloud-cover", type=float, default=1.0,
                    help="Fraction of the disk covered by the cloudy column "
                         "(default 1.0 = fully cloudy, the most stringent ice "
                         "test). Mixed against a clear companion column.")
    ap.add_argument("--asurf", type=float, default=0.05,
                    help="Lambertian albedo (dark-ocean proxy). The B-vs-C "
                         "difference is surface-independent; use a real ocean "
                         "BPDF only if you also want glint realism.")
    ap.add_argument("--cache", default=None,
                    help="Baum .npz cache (default: examples/baum_cache for --deff).")
    args = ap.parse_args()

    wvls = np.array([float(x) for x in args.wvls.split(",")])
    phases = np.array([float(x) for x in args.phases.split(",")])
    cache = args.cache or os.path.join(
        HERE, "baum_cache",
        "GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix"
        "_Deff{:.0f}.npz".format(args.deff))
    if not os.path.exists(cache):
        print("Baum cache not found:", cache)
        print("Build it first:  python examples/convert_baum_to_pymiedap.py "
              "--deff {:.0f} --wlmin {:.2f} --wlmax {:.2f}".format(
                  args.deff, wvls.min(), wvls.max()))
        return 1

    print("Earth-like ice sensitivity | D_eff=%.0f um | nmug=%d  M_liq=%d  "
          "M_ice: B=%d C=%d" % (args.deff, args.nmug, args.m_liq,
                                args.m_ice_b, args.m_ice_c))
    clear = build_clear(wvls, args.nmug, args.asurf, "ice_sens_clear")
    cases = {}
    specs = [("A", None), ("B", args.m_ice_b), ("C", args.m_ice_c)]
    for tag, mice in specs:
        m = build_and_compute(wvls, tag, cache, mice, args.m_liq, args.nmug,
                              args.nmug_mie, args.asurf, "ice_sens_%s" % tag)
        I, Q, U = disk(m, clear, phases, args.npix, args.cloud_cover, args.nmug)
        finite = np.isfinite(I).all() and np.isfinite(Q).all()
        cases[tag] = dict(I=np.array(I), Q=np.array(Q), U=np.array(U), ok=finite)
        print("  case %s computed; finite=%s%s" % (
            tag, finite, "" if finite else "  <-- UNSTABLE (raise its M)"))

    A, B, C = cases["A"], cases["B"], cases["C"]
    if not (B["ok"] and C["ok"]):
        print("\nB and/or C are UNSTABLE at these M for D_eff=%.0f um ice." % args.deff)
        print("=> there is no cheap *stable* delta-M order for this ice; either")
        print("   raise --m-ice-* (and --nmug>=M/2) until finite, or develop the")
        print("   delta-fit positivity-preserving truncation.")
        return 0

    def P(c):  # signed degree of pol, and total
        return -c["Q"] / c["I"], np.sqrt(c["Q"]**2 + c["U"]**2) / c["I"]
    PA, _ = P(A); PB, _ = P(B); PC, _ = P(C)

    print("\n wl[um] alpha |   F_A     F_B     F_C   | dF(B-C) | P_A    P_B    P_C  | dP(B-C)pp | dF(A-B)")
    maxdF_bc = maxdQ_bc = maxdP_bc = 0.0
    for i, w in enumerate(wvls):
        for j, al in enumerate(phases):
            FA, FB, FC = A["I"][i, j], B["I"][i, j], C["I"][i, j]
            dF_bc = abs(FB - FC) / abs(FC) * 100
            dQ_bc = abs(B["Q"][i, j] - C["Q"][i, j]) / max(abs(C["Q"][i, j]), 1e-9) * 100
            dP_bc = abs(PB[i, j] - PC[i, j]) * 100   # percentage points
            dF_ab = abs(FA - FB) / abs(FB) * 100
            maxdF_bc = max(maxdF_bc, dF_bc); maxdQ_bc = max(maxdQ_bc, dQ_bc)
            maxdP_bc = max(maxdP_bc, dP_bc)
            print(" %5.2f %5.0f | %.5f %.5f %.5f | %6.2f%% | %+.3f %+.3f %+.3f | %7.3f | %6.2f%%"
                  % (w, al, FA, FB, FC, dF_bc, PA[i, j], PB[i, j], PC[i, j], dP_bc, dF_ab))

    print("\nMax B-vs-C differences:  dF=%.2f%%  dQ=%.2f%%  dP=%.3f pp" %
          (maxdF_bc, maxdQ_bc, maxdP_bc))
    ok = maxdF_bc < 1.0 and maxdQ_bc < 3.0 and maxdP_bc < 0.5
    print("\nDECISION:", (
        "PASS -- cheap stable ice (M~%d) is good enough; the ice forward-peak "
        "accuracy does NOT materially change the Earth-like spectrum." % args.m_ice_b
        if ok else
        "FAIL -- the ice forward-peak treatment DOES matter (B != C beyond "
        "tolerance); proceed to develop delta-fit / positivity-preserving "
        "truncation."))
    print("(A-vs-B above shows how much the ice layer changes the spectrum at all.)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
