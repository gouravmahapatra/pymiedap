#!/usr/bin/env python3
# This file is part of PyMieDAP, released under GNU General Public License.
"""
validate_deltafit.py  (Phase 4)
=============================================================================
End-to-end test of the delta-fit ice truncation against the expensive delta-M
reference, in the full Earth-like disk-integrated geometry.

Two ice treatments, identical everything else (liquid cloud, surface):
    REF : Baum ice via delta-M at a high, near-converged order (M=180, nmug~95)
          -- this is the slow "truth" from the convergence ladder.
    DF  : Baum ice via DELTA-FIT at a small order (M=60, nmug~35) -- cheap.

If DF reproduces REF to dF<1% and dP<0.5 pp across phase/wavelength, then
delta-fit gives the converged polarized spectrum at a fraction of the cost
(~(35/95)^4 ~ 2%). That validates delta-fit as the production ice route.

Prereqs on the cluster (in the pymiedap conda env, from ~/pymiedap):
  * the Baum NetCDF (or .nc.gz) present, and a delta-M cache built for REF:
      python examples/convert_baum_to_pymiedap.py --deff 60 --wlmin 0.4 --wlmax 1.0
  * pass the NetCDF path with --baum-nc (for the delta-fit branch).

Usage (run backgrounded; REF at nmug~95 is the slow part):
  nohup python examples/validate_deltafit.py --deff 60 --liq-reff 4 \
      --baum-nc GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc.gz \
      --wvls 0.55,0.85 --phases 60,90,120,150 > df.log 2>&1 &
"""

import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

# Reuse the validated Earth-like machinery.
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "its", os.path.join(HERE, "ice_thin_sensitivity.py"))
its = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(its)

import pymiedap.pymiedap as pmd                       # noqa: E402
from pymiedap.tmatrix import delta_m_truncate          # noqa: E402
from pymiedap.baum import fill_aerosol_from_cache      # noqa: E402
from pymiedap.deltafit import deltafit_from_baum       # noqa: E402


def build_ice_case(wvls, method, args, name):
    """Build+compute the liquid+ice column with the ice done by `method`
    ('deltam' from the cache at M=180, or 'deltafit' from the NetCDF at M=60)."""
    its.LIQ_REFF = args.liq_reff
    wvls = np.asarray(wvls, float)
    nmug = args.nmug_ref if method == "deltam" else args.nmug_df
    m = pmd.Model(wvl_list=wvls); m.asurf = args.asurf
    L = m.layers; nw = len(wvls)
    L.gastop.press = its.P_GASTOP; L.gastop.tau = [0.0]
    L.gasbelow.press = its.P_SURF; L.gasbelow.tau = [0.0]
    L.haze.press = its.P_ICE; L.haze.tau = [its.ICE_TAU]
    L.cloud.press = its.P_LIQ; L.cloud.tau = [its.LIQ_TAU]
    a = L.cloud.aerosols
    a.typ = 'C'; a.psd = '2'; a.r_eff = args.liq_reff; a.v_eff = its.LIQ_VEFF
    a.nr = [its.LIQ_NR] * max(2, nw); a.ni = [its.LIQ_NI] * max(2, nw)
    m.wvl_list = m.wvl_list

    for lay, layer in vars(m.layers).items():
        if hasattr(layer, 'mixed_aerosols'):
            del layer.mixed_aerosols
        if lay == 'cloud':
            pmd.mie_code(layer.aerosols, wvls, ngaur=args.nmug_mie, nsubr=50)
            ts = delta_m_truncate(layer.aerosols, args.m_liq)
            layer.tau = np.atleast_1d(layer.tau).astype(float)
            if layer.tau.size != nw:
                layer.tau = layer.tau[0] * np.ones(nw)
            layer.tau = layer.tau * ts
            layer.mix_aerosols()
        elif lay == 'haze':
            if method == "deltam":
                fill_aerosol_from_cache(layer.aerosols, args.cache, wavelengths_um=wvls)
                layer.aerosols.f = 1.0
                ts = delta_m_truncate(layer.aerosols, args.m_ref)
            else:  # deltafit
                ts = deltafit_from_baum(layer.aerosols, args.baum_nc, wvls,
                                        args.deff, args.m_df, theta_cut=args.theta_cut)
            layer.tau = np.atleast_1d(layer.tau).astype(float)
            if layer.tau.size != nw:
                layer.tau = layer.tau[0] * np.ones(nw)
            layer.tau = layer.tau * ts
            layer.mix_aerosols()
        else:
            for an, aero in vars(layer).items():
                if isinstance(aero, pmd.Aerosols):
                    pmd.mie_code(aero, wvls, ngaur=20, nsubr=50)
            layer.mix_aerosols()

    pmd.dap_code(m, rename=True, output_name=name, nmug=nmug, nmat=4,
                 path_output=its.DAP_DIR)
    return m, nmug


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deff", type=float, default=60.0)
    ap.add_argument("--liq-reff", type=float, default=4.0)
    ap.add_argument("--m-liq", type=int, default=80)
    ap.add_argument("--m-ref", type=int, default=180, help="delta-M reference order.")
    ap.add_argument("--nmug-ref", type=int, default=95)
    ap.add_argument("--m-df", type=int, default=60, help="delta-fit order.")
    ap.add_argument("--nmug-df", type=int, default=35)
    ap.add_argument("--theta-cut", type=float, default=5.0)
    ap.add_argument("--nmug-mie", type=int, default=120)
    ap.add_argument("--npix", type=int, default=15)
    ap.add_argument("--cloud-cover", type=float, default=1.0)
    ap.add_argument("--asurf", type=float, default=0.05)
    ap.add_argument("--wvls", default="0.55,0.85")
    ap.add_argument("--phases", default="60,90,120,150")
    ap.add_argument("--cache", default=None, help="delta-M cache for REF.")
    ap.add_argument("--baum-nc", required=True, help="Baum NetCDF (.nc/.nc.gz) for delta-fit.")
    args = ap.parse_args()

    wvls = np.array([float(x) for x in args.wvls.split(",")])
    phases = np.array([float(x) for x in args.phases.split(",")])
    args.cache = args.cache or os.path.join(
        HERE, "baum_cache",
        "GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix"
        "_Deff{:.0f}.npz".format(args.deff))

    print("Phase 4: delta-fit (M=%d, nmug=%d) vs delta-M reference (M=%d, nmug=%d)"
          % (args.m_df, args.nmug_df, args.m_ref, args.nmug_ref))
    clear = its.build_clear(wvls, max(args.nmug_ref, args.nmug_df), args.asurf,
                            "df_clear")

    mref, nref = build_ice_case(wvls, "deltam", args, "df_ref")
    Iref, Qref, Uref = its.disk(mref, clear, phases, args.npix, args.cloud_cover, nref)
    Iref, Qref = np.array(Iref), np.array(Qref)
    mdf, ndf = build_ice_case(wvls, "deltafit", args, "df_fit")
    Idf, Qdf, Udf = its.disk(mdf, clear, phases, args.npix, args.cloud_cover, ndf)
    Idf, Qdf = np.array(Idf), np.array(Qdf)

    Pref, Pdf = -Qref / Iref, -Qdf / Idf
    print("\n wl[um] alpha |  F_ref    F_df   dF%  | P_ref   P_df   dP(pp)")
    maxdF = maxdP = 0.0
    for i, w in enumerate(wvls):
        for j, al in enumerate(phases):
            dF = abs(Idf[i, j] - Iref[i, j]) / abs(Iref[i, j]) * 100
            dP = abs(Pdf[i, j] - Pref[i, j]) * 100
            maxdF = max(maxdF, dF); maxdP = max(maxdP, dP)
            print(" %5.2f %5.0f | %.5f %.5f %5.2f | %+.3f %+.3f %6.3f"
                  % (w, al, Iref[i, j], Idf[i, j], dF, Pref[i, j], Pdf[i, j], dP))
    print("\nmax dF=%.2f%%  max dP=%.3f pp  | cost ratio df/ref ~ %.1f%%"
          % (maxdF, maxdP, (args.nmug_df / args.nmug_ref) ** 4 * 100))
    if maxdF < 1.0 and maxdP < 0.5:
        print("DECISION: delta-fit VALIDATED -- reproduces the M=%d reference to "
              "<0.5 pp at a fraction of the cost. Use delta-fit (M=%d, nmug=%d) "
              "as the production ice route." % (args.m_ref, args.m_df, args.nmug_df))
    else:
        print("DECISION: delta-fit differs from the reference beyond tolerance; "
              "raise --m-df / lower --theta-cut and retest.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
