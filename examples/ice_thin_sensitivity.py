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


FLOOR = 1e-4   # reflectance below this == non-converged (collapsed) DAP solve


def compute_case(wvls, phases, tag, m_ice, args, clear, cache, name):
    """Build+disk one case. nmug auto-raised so nmug >= M_ice/2 for ice cases."""
    nmug = args.nmug if m_ice is None else max(args.nmug, (m_ice + 1) // 2 + 5)
    m = build_and_compute(wvls, tag, cache, m_ice, args.m_liq, nmug,
                          args.nmug_mie, args.asurf, name)
    I, Q, U = disk(m, clear, phases, args.npix, args.cloud_cover, nmug)
    I = np.array(I); Q = np.array(Q); U = np.array(U)
    ok = np.isfinite(I).all() and bool((np.abs(I) >= FLOOR).all())
    badwl = [float(wvls[k]) for k in range(len(wvls))
             if np.any(np.abs(I[k]) < FLOOR)]
    return dict(I=I, Q=Q, U=U, ok=ok, m_ice=m_ice, nmug=nmug, badwl=badwl)


def _pol(c):
    return -c["Q"] / c["I"]


def _compare(x, y):
    """Max reflectance diff [%] and max polarization diff [pp] between cases."""
    dF = float(np.max(np.abs(x["I"] - y["I"]) / np.abs(y["I"])) * 100)
    dP = float(np.max(np.abs(_pol(x) - _pol(y))) * 100)
    return dF, dP


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deff", type=float, default=60.0)
    ap.add_argument("--m-ice-b", type=int, default=24)
    ap.add_argument("--m-ice-c", type=int, default=60)
    ap.add_argument("--liq-reff", type=float, default=4.0,
                    help="Liquid cloud effective radius [um]. Held identical in "
                         "A/B/C so it cancels in B-vs-C; keep it small enough to "
                         "stay stable at the bluest wavelength (r_eff=9 um is "
                         "unstable below ~0.7 um at affordable M_liq/nmug).")
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
    ap.add_argument("--auto", action="store_true",
                    help="Auto-bracket: raise M_ice on a ladder until the "
                         "polarization converges (dP<tol) or --max-m is hit.")
    ap.add_argument("--m-start", type=int, default=60, help="auto: first M_ice.")
    ap.add_argument("--m-step", type=int, default=30, help="auto: M_ice increment.")
    ap.add_argument("--max-m", type=int, default=150, help="auto: largest M_ice.")
    ap.add_argument("--tol-df", type=float, default=1.0, help="F convergence [%].")
    ap.add_argument("--tol-dp", type=float, default=0.5, help="P convergence [pp].")
    args = ap.parse_args()

    global LIQ_REFF
    LIQ_REFF = args.liq_reff
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

    print("Earth-like ice sensitivity | D_eff=%.0f um | liq r_eff=%.1f um "
          "tau_ice=%.2f | M_liq=%d nmug=%d" % (args.deff, LIQ_REFF, ICE_TAU,
                                               args.m_liq, args.nmug))
    clear = build_clear(wvls, args.nmug, args.asurf, "ice_sens_clear")

    # Case A: liquid only (no ice). Must be valid or nothing else is meaningful.
    A = compute_case(wvls, phases, "A", None, args, clear, cache, "ice_sens_A")
    print("  Case A (liquid only): valid=%s%s" % (
        A["ok"], "" if A["ok"] else "  bad lambda=" +
        ",".join("%.2f" % w for w in A["badwl"])))
    if not A["ok"]:
        print("\nINCONCLUSIVE -- liquid-only Case A collapsed (numerics, not an "
              "ice result). The liquid is identical in all cases so it cancels "
              "in the comparison: lower --liq-reff or raise --m-liq/--nmug until "
              "Case A is valid at every wavelength, then re-run.")
        return 0

    # -------------------------------------------------------------------------
    if args.auto:
        print("\nAUTO convergence ladder: M_ice = %d, %d, ... <= %d "
              "(stop when dP<%.2f pp and dF<%.2f%%)" %
              (args.m_start, args.m_start + args.m_step, args.max_m,
               args.tol_dp, args.tol_df))
        print("  M_ice nmug | ice-vs-noice: dF%  dP(pp) | step-to-step: dF%  dP(pp)")
        prev = None; converged = None; last = None
        M = args.m_start
        while M <= args.max_m:
            c = compute_case(wvls, phases, "I", M, args, clear, cache,
                             "ice_sens_M%d" % M)
            if not c["ok"]:
                print("  %5d %4d | UNSTABLE (collapsed) -- skipping" % (M, c["nmug"]))
                M += args.m_step; continue
            dFA, dPA = _compare(c, A)            # how much ice changes the spectrum
            if prev is None:
                print("  %5d %4d | %5.1f %6.2f |   (baseline)" %
                      (M, c["nmug"], dFA, dPA))
            else:
                dF, dP = _compare(c, prev)       # step-to-step convergence
                print("  %5d %4d | %5.1f %6.2f | %5.2f %6.3f" %
                      (M, c["nmug"], dFA, dPA, dF, dP))
                if dF < args.tol_df and dP < args.tol_dp:
                    converged = (prev["m_ice"], M); last = c; break
            prev = c; last = c; M += args.m_step

        print()
        if converged:
            print("CONVERGED: polarization stable between M_ice=%d and %d "
                  "(dP<%.2f pp). Production ice setting: M_ice=%d, nmug=%d. "
                  "delta-fit NOT needed." %
                  (converged[0], converged[1], args.tol_dp, converged[1],
                   last["nmug"]))
        else:
            print("NOT CONVERGED up to M_ice=%d: the polarization keeps drifting "
                  "more than %.2f pp per step. Brute-force refinement is too slow "
                  "=> develop delta-fit / positivity-preserving truncation." %
                  (args.max_m, args.tol_dp))
        return 0

    # -------------------------------------------------------------------------
    # Explicit two-order comparison (B vs C), with the dP(pp) verdict.
    B = compute_case(wvls, phases, "B", args.m_ice_b, args, clear, cache, "ice_sens_B")
    C = compute_case(wvls, phases, "C", args.m_ice_c, args, clear, cache, "ice_sens_C")
    for tag, c in (("B", B), ("C", C)):
        print("  Case %s (ice M=%d, nmug=%d): valid=%s%s" % (
            tag, c["m_ice"], c["nmug"], c["ok"], "" if c["ok"] else
            "  bad lambda=" + ",".join("%.2f" % w for w in c["badwl"])))
    if not (B["ok"] and C["ok"]):
        stable = [c["m_ice"] for c in (B, C) if c["ok"]]
        print("\nUnstable ice order(s) present. Stable: %s. Re-run with two "
              "stable orders, or use --auto to find the converged M "
              "automatically." % (stable or "none"))
        return 0

    PA, PB, PC = _pol(A), _pol(B), _pol(C)
    print("\n wl[um] alpha |   F_A     F_B     F_C   | dF(B-C) | P_A    P_B    P_C  | dP(B-C)pp | dF(A-B)")
    maxdF = maxdP = 0.0
    for i, w in enumerate(wvls):
        for j, al in enumerate(phases):
            dF = abs(B["I"][i, j] - C["I"][i, j]) / abs(C["I"][i, j]) * 100
            dP = abs(PB[i, j] - PC[i, j]) * 100
            dFab = abs(A["I"][i, j] - B["I"][i, j]) / abs(B["I"][i, j]) * 100
            maxdF = max(maxdF, dF); maxdP = max(maxdP, dP)
            print(" %5.2f %5.0f | %.5f %.5f %.5f | %6.2f%% | %+.3f %+.3f %+.3f | %7.3f | %6.2f%%"
                  % (w, al, A["I"][i, j], B["I"][i, j], C["I"][i, j], dF,
                     PA[i, j], PB[i, j], PC[i, j], dP, dFab))
    print("\nB-vs-C convergence:  max dF=%.2f%% (tol %.1f%%)   max dP=%.3f pp "
          "(tol %.2f pp)" % (maxdF, args.tol_df, maxdP, args.tol_dp))
    if maxdF < args.tol_df and maxdP < args.tol_dp:
        print("DECISION: CONVERGED -- M_ice=%d is good enough (reflectance AND "
              "polarization stable). Use it for production." % args.m_ice_b)
    else:
        print("DECISION: NOT converged at M_ice=%d (polarization still drifting). "
              "Push higher M (or run --auto); if it won't converge by ~M=150, "
              "develop delta-fit." % args.m_ice_b)
    print("(dF(A-B) shows how much the ice layer changes the spectrum at all.)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
