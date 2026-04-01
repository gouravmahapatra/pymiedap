#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pymiedap.pymiedap as pmd


NOTEBOOK_TOL = 5e-5
LAMBERT_TOL = 2e-4


def compare(name: str, actual: np.ndarray, expected: np.ndarray, tol: float = NOTEBOOK_TOL) -> dict[str, float | bool]:
    diff = np.abs(np.asarray(actual) - np.asarray(expected))
    max_abs = float(np.max(diff))
    return {
        "name": name,
        "max_abs_diff": max_abs,
        "passed": bool(max_abs <= tol),
    }


def build_aerosols() -> tuple[pmd.Aerosols, pmd.Aerosols, pmd.Aerosols, pmd.Aerosols]:
    aer_a = pmd.Aerosols(nr=[1.45], ni=[0], r_eff=0.23, v_eff=0.18)
    aer_b = pmd.Aerosols(nr=[1.44], ni=[0], r_eff=1.05, v_eff=0.07)
    aer_c = pmd.Aerosols(nr=[1.33], ni=[0], v_eff=0.07, r_eff=2, par3=0.5, psd="7")
    aer_d = pmd.Aerosols(nr=[1.33], ni=[0], r_eff=2.2, v_eff=0.07)

    pmd.mie_code(aer_a, [0.55], nsubr=40, ngaur=500)
    pmd.mie_code(aer_b, [0.55], nsubr=40, ngaur=500)
    pmd.mie_code(aer_c, [0.70], nsubr=40, ngaur=500)
    pmd.mie_code(aer_d, [0.70], nsubr=40, ngaur=500)
    return aer_a, aer_b, aer_c, aer_d


def build_notebook_models(aer_c: pmd.Aerosols) -> tuple[pmd.Model, pmd.Model]:
    model_b = pmd.Model()
    model_b.wvl_list = [0.7]
    del model_b.layers.gasbelow
    del model_b.layers.haze
    model_b.layers.gastop.rayscat = False
    model_b.layers.gastop.tau_ray = [0.1]
    model_b.layers.cloud.rayscat = False
    model_b.layers.cloud.tau = [0.4]
    model_b.layers.cloud.tau_ray = [0.1]
    model_b.layers.cloud.aerosols = aer_c
    model_b.dpol = 0.0279
    model_b.surface[0, 0] = 0.1

    model_a = pmd.Model()
    model_a.wvl_list = [0.7]
    del model_a.layers.gastop
    del model_a.layers.haze
    del model_a.layers.cloud
    model_a.layers.gasbelow.rayscat = False
    model_a.layers.gasbelow.tau = [1.0]
    model_a.layers.gasbelow.aerosols = aer_c
    model_a.surface[0, 0] = 0.0

    return model_a, model_b


def geometry_arrays() -> tuple[np.ndarray, ...]:
    mus = np.array([0.1, 0.5, 1.0, 0.1, 0.5, 1.0, 0.1, 0.5, 1.0, 0.1, 0.5, 1.0])
    emissions = np.degrees(np.arccos(mus))
    mu0s = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    szas = np.degrees(np.arccos(mu0s))
    dphi = np.radians([0, 0, 0, 30, 30, 30, 0, 0, 0, 30, 30, 30])
    alphas = mus * mu0s + np.sqrt(1 - mus**2) * np.sqrt(1 - mu0s**2) * np.cos(dphi)
    phases = np.degrees(np.arccos(alphas))
    azimuth = pmd.calc_azimuth(phases, szas, emissions, deg=True)
    return mus, emissions, mu0s, szas, dphi, phases, azimuth


def run() -> int:
    coordsx = [0, 1, 2, 3, 0, 2]
    coordsy = [0, 1, 2, 3, 1, 3]

    aer_a, aer_b, aer_c, aer_d = build_aerosols()
    model_a, model_b = build_notebook_models(aer_c)
    mus, emissions, mu0s, szas, dphi, phases, azimuth = geometry_arrays()

    pmd.compute_model(model_a, output_name="benchA", rename=True, nmat=4, nmug=40)
    pmd.compute_model(model_b, output_name="benchB", rename=True, nmat=4, nmug=40)

    ia, qa, ua, va = pmd.read_dap_output(
        phases,
        szas,
        emissions,
        model_a.name[0],
        beta=np.zeros_like(dphi),
        phi=np.degrees(dphi),
    )
    ib, qb, ub, vb = pmd.read_dap_output(
        phases,
        szas,
        emissions,
        model_b.name[0],
        beta=np.zeros_like(dphi),
        phi=np.degrees(dphi),
    )

    notebook_checks = [
        compare("aerA asym", aer_a.asym, np.array([0.72100229])),
        compare("aerB asym", aer_b.asym, np.array([0.7179972])),
        compare("aerC asym", aer_c.asym, np.array([0.804201])),
        compare("aerD asym", aer_d.asym, np.array([0.80188248])),
        compare(
            "aerA coef l0",
            aer_a.coefs[0, coordsx, coordsy, 0],
            np.array([1.0, 0.0, 0.0, 0.92901809, 0.0, 0.0]),
        ),
        compare(
            "aerA coef l10",
            aer_a.coefs[0, coordsx, coordsy, 10],
            np.array([0.08802913, 0.11121631, 0.10151652, 0.0832367, -0.00358954, 0.02614555]),
        ),
        compare(
            "aerA coef l15",
            aer_a.coefs[0, coordsx, coordsy, 15],
            np.array([0.00228455, 0.00218155, 0.00224082, 0.00255047, 0.00033707, 0.00046371]),
        ),
        compare(
            "aerA coef l24",
            aer_a.coefs[0, coordsx, coordsy, 24],
            np.array([5.48015333e-06, 6.77440344e-06, 2.37713944e-06, 1.97521301e-06, -7.37991466e-07, 3.53992457e-06]),
        ),
        compare(
            "aerB coef l0",
            aer_b.coefs[0, coordsx, coordsy, 0],
            np.array([1.0, 0.0, 0.0, 0.86462372, 0.0, 0.0]),
        ),
        compare(
            "aerB coef l10",
            aer_b.coefs[0, coordsx, coordsy, 10],
            np.array([5.06996479, 5.33408747, 5.31811435, 5.12294652, 0.03998552, 0.02898462]),
        ),
        compare(
            "aerB coef l15",
            aer_b.coefs[0, coordsx, coordsy, 15],
            np.array([4.29325014, 4.22823484, 4.16063179, 4.35190093, -0.04425792, 0.31830627]),
        ),
        compare(
            "aerB coef l24",
            aer_b.coefs[0, coordsx, coordsy, 24],
            np.array([1.48138437, 1.66704814, 1.50932983, 1.42450347, -0.04210442, 0.24849267]),
        ),
        compare(
            "aerC coef l0",
            aer_c.coefs[0, coordsx, coordsy, 0],
            np.array([1.0, 0.0, 0.0, 0.95461472, 0.0, 0.0]),
        ),
        compare(
            "aerC coef l10",
            aer_c.coefs[0, coordsx, coordsy, 10],
            np.array([1.17713563, 1.2847585, 1.2659546, 1.17411335, -0.02230815, 0.09498393]),
        ),
        compare(
            "aerC coef l15",
            aer_c.coefs[0, coordsx, coordsy, 15],
            np.array([0.42561199, 0.42569832, 0.4231388, 0.4298012, -0.0055369, 0.04263078]),
        ),
        compare(
            "aerC coef l80",
            aer_c.coefs[0, coordsx, coordsy, 80],
            np.array([1.28782521e-06, 1.40241107e-06, 1.04194570e-06, 9.98265353e-07, 8.49592384e-08, 4.62722533e-07]),
        ),
        compare(
            "aerD coef l0",
            aer_d.coefs[0, coordsx, coordsy, 0],
            np.array([1.0, 0.0, 0.0, 0.91852089, 0.0, 0.0]),
        ),
        compare(
            "aerD coef l10",
            aer_d.coefs[0, coordsx, coordsy, 10],
            np.array([7.20704159, 7.35295547, 7.25693617, 7.13468861, -0.029961, 0.09669141]),
        ),
        compare(
            "aerD coef l15",
            aer_d.coefs[0, coordsx, coordsy, 15],
            np.array([7.98539793, 8.00874554, 7.98884822, 7.99970188, 0.02916147, 0.16475289]),
        ),
        compare(
            "aerD coef l80",
            aer_d.coefs[0, coordsx, coordsy, 80],
            np.array([3.04667428e-03, 3.18324854e-03, 2.83399090e-03, 2.77620381e-03, 2.74463843e-05, 5.59438200e-04]),
        ),
        compare(
            "geometry phases",
            phases,
            np.array([24.26082952, 8.53773646e-07, 60.0, 37.2274111, 25.9050793, 60.0, 0.0, 24.26082952, 84.26082952, 29.8461183, 37.2274111, 84.26082952]),
        ),
        compare(
            "geometry emissions",
            emissions,
            np.array([84.26082952, 60.0, 0.0, 84.26082952, 60.0, 0.0, 84.26082952, 60.0, 0.0, 84.26082952, 60.0, 0.0]),
        ),
        compare(
            "geometry sza",
            szas,
            np.array([60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 84.26082952, 84.26082952, 84.26082952, 84.26082952, 84.26082952, 84.26082952]),
        ),
        compare(
            "geometry azimuth",
            azimuth,
            np.array([8.53773646e-07, 8.53773646e-07, 0.0, 30.0, 30.0, 0.0, 0.0, 8.53773646e-07, 0.0, 30.0, 30.0, 0.0]),
        ),
        compare("modelA bmsca", model_a.bmsca, np.zeros(50)),
        compare("modelA basca", model_a.basca, np.array([1.0] + [0.0] * 49)),
        compare("modelB bmsca", model_b.bmsca, np.array([0.1, 0.1] + [0.0] * 48)),
        compare("modelB basca", model_b.basca, np.array([0.4] + [0.0] * 49)),
        compare(
            "modelA I",
            ia * mu0s * np.pi,
            np.array([1.10269, 0.31943, 0.033033, 0.66414, 0.25209, 0.033033, 2.93214, 0.22054, 0.009287, 0.76910, 0.132828, 0.009287]),
        ),
        compare(
            "modelA Q",
            qa * mu0s * np.pi,
            np.array([0.004604, -0.002881, -0.002979, 0.000303, -0.001444, -0.001489, 0.009900, 0.000976, -0.000815, -0.003758, 0.000220, -0.000408]),
        ),
        compare(
            "modelA U",
            ua * mu0s * np.pi,
            np.array([0.0, 0.0, 0.0, -0.002770, -0.004141, -0.002580, 0.0, 0.0, 0.0, 0.003124, -0.000525, -0.000706]),
        ),
        compare(
            "modelA V",
            va * mu0s * np.pi,
            np.array([0.0, 0.0, 0.0, 0.000038, 0.000017, 0.0, 0.0, 0.0, 0.0, 0.000012, 0.000007, 0.0]),
        ),
        compare(
            "modelB I",
            ib * mu0s * np.pi,
            np.array([0.53295, 0.20843, 0.093680, 0.41814, 0.18497, 0.093680, 0.52277, 0.106590, 0.026009, 0.27630, 0.083628, 0.026009]),
        ),
        compare(
            "modelB Q",
            qb * mu0s * np.pi,
            np.array([-0.028340, -0.036299, -0.024156, -0.000058, -0.019649, -0.012078, 0.011506, -0.005186, -0.014984, 0.034368, 0.003839, -0.007492]),
        ),
        compare(
            "modelB U",
            ub * mu0s * np.pi,
            np.array([0.0, 0.0, 0.0, -0.073105, -0.041401, -0.020920, 0.0, 0.0, 0.0, -0.016042, -0.014492, -0.012976]),
        ),
        compare(
            "modelB V",
            vb * mu0s * np.pi,
            np.array([0.0, 0.0, 0.0, 0.000106, 0.000040, 0.0, 0.0, 0.0, 0.0, 0.000027, 0.000017, 0.0]),
        ),
    ]

    lambert_model = pmd.Model()
    lambert_model.wvl_list = [0.7]
    del lambert_model.layers.gastop
    del lambert_model.layers.haze
    del lambert_model.layers.cloud
    lambert_model.layers.gasbelow.rayscat = False
    lambert_model.layers.gasbelow.tau = [0.0]
    lambert_model.surface[0, 0] = 1.0

    alphas = np.linspace(0, np.pi, 80)
    alphas_deg = np.degrees(alphas)
    theta = np.pi - alphas
    lambert_expected = 2 * (np.sin(theta) - theta * np.cos(theta)) / (3.0 * np.pi)

    pmd.planet_integrated([lambert_model], npix=60, alpha=alphas_deg, output_names=["benchLambert"])
    lambert_diff = lambert_model.I[0, :] - lambert_expected
    finite_mask = np.isfinite(lambert_diff)
    lambert_max_abs = float(np.max(np.abs(lambert_diff[finite_mask])))
    lambert_nonfinite = np.where(~finite_mask)[0]

    print("Notebook-value checks")
    failed = False
    for result in notebook_checks:
        status = "PASS" if result["passed"] else "FAIL"
        print(f"{status:4} {result['name']:<20} max_abs_diff={result['max_abs_diff']:.8e}")
        failed = failed or (not result["passed"])

    print("")
    print("Lambertian check")
    print(f"finite max_abs_diff={lambert_max_abs:.8e}")
    print(f"nonfinite tail indices={lambert_nonfinite.tolist()}")
    if lambert_max_abs > LAMBERT_TOL:
        print(f"FAIL finite Lambertian residual exceeds tolerance {LAMBERT_TOL:.1e}")
        failed = True
    else:
        print(f"PASS finite Lambertian residual within tolerance {LAMBERT_TOL:.1e}")

    if len(lambert_nonfinite) > 0:
        print("FAIL Lambertian phase curve contains non-finite values.")
        failed = True

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(run())
