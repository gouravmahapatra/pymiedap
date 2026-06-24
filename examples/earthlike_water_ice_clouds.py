#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
earthlike_water_ice_clouds.py
=============================================================================
Set up a PyMieDAP two-cloud-layer (liquid water + ice) Earth-like atmosphere
and compute the disk-integrated reflected-light SPECTRUM and PHASE CURVE in
intensity and polarization.

The cloud parameters reproduce the disk-averaged Earth-like values of

    Roccetti, Emde, Sterzik, Manev, Seidel & Bagnulo (2025),
    "Planet Earth in reflected and polarized light. I.", A&A,
    arXiv:2504.02048  --  Table 2 (last row, "average").

Table 2 (Earth-like planet, alpha = 90 deg) disk-averaged clouds
--------------------------------------------------------------------
                       cloud cover   height H   r_eff      tau
    liquid water clouds   38.0 %      1.55 km    8.6 um     5.44
    ice    water clouds   51.5 %      3.88 km   44.3 um     0.62
    all clouds            57.0 %        --        --         --

PyMieDAP layers used here (top -> bottom):
    gastop   : optically thin gas above the clouds
    cloud_ice: ice cloud, tau ~ 0.62, scattering matrix from T-matrix coeffs
    cloud_liq: liquid water cloud, tau ~ 5.44, Mie (r_eff = 8.6 um)
    gasbelow : gas column below the clouds, down to the surface

=============================================================================
IMPORTANT MODELLING CAVEATS  (read before trusting any number)
=============================================================================
1. ICE OPTICS / WAVELENGTH COVERAGE.
   The only T-matrix ice scattering file shipped with the repo is
   tmatrix_ice/ice_oblate_0.5um.coeffs : lambda = 0.5 um, r_eff ~ 1 um,
   oblate spheroid (a/b = 2).  The T-matrix (EBCM) method is numerically
   UNSTABLE for the paper's ~44 um ice crystals (size parameter ~550), which
   is precisely why the paper uses geometric-optics / the six-habit HEY
   parameterization for ice.  So:
     * The PHASE CURVE is computed at 0.5 um, where the T-matrix data is
       physically valid, and it carries the non-sphericity signature
       (F22 != F11) that Mie spheres cannot reproduce.
     * The SPECTRUM holds the ice single-scattering properties fixed at the
       0.5 um T-matrix values across the VIS-NIR band.  This is a deliberate
       "grey ice" approximation, justified because ice is essentially
       non-absorbing in 0.4-0.9 um (single-scattering albedo W ~ 0.99998),
       but it does NOT capture the size/wavelength dependence of real cirrus.
   To do better, generate one .coeffs file per wavelength with the Fortran
   tool in tmatrix_ice/ (see its README: change MRR/MRI for the ice
   refractive index, AXMAX/RAT for size) and pass them via ICE_COEFF_FILES.

2. PLANE-PARALLEL COLUMNS.  PyMieDAP integrates independent 1-D columns over
   the disk.  It cannot reproduce the 3-D sub-grid cloud variability that the
   paper found dominant; this script emulates the paper's MEAN cloud field.

3. CLOUD COVER.  The paper gives different cover for liquid (38%) and ice
   (51.5%).  Here both cloud layers sit in one "cloudy" column and the disk
   integration mixes that column with a clear column at the combined
   all-cloud cover (57%).  Refine with separate columns if you need the
   liquid/ice cover split.
=============================================================================

Run:
    python examples/earthlike_water_ice_clouds.py

Outputs (written next to this script):
    earthlike_spectrum.csv / .png        reflectance + P vs wavelength
    earthlike_phasecurve_0.5um.csv / .png reflectance + P vs phase angle
"""

import os
import sys
import numpy as np

# --------------------------------------------------------------------------
# Repo layout: this file lives in <repo>/examples/. PyMieDAP and the compiled
# native modules (module_mie, module_dap, ...) live in <repo>/.
# --------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import pymiedap.pymiedap as pmd          # noqa: E402

# =========================================================================
# 1.  PAPER PARAMETERS  (Roccetti et al. 2025, Table 2, Earth-like average)
# =========================================================================
LIQ_TAU    = 5.44      # liquid water cloud optical thickness (at 550 nm)
LIQ_REFF   = 8.6       # um, effective radius
LIQ_VEFF   = 0.10      # effective variance (typical for water clouds)
LIQ_H_KM   = 1.55      # km, mean cloud height
LIQ_NR     = 1.335     # water refractive index, real part (~VIS)
LIQ_NI     = 1.0e-8    # water, imaginary part (negligible absorption in VIS)

ICE_TAU    = 0.62      # ice water cloud optical thickness (at 550 nm)
ICE_H_KM   = 3.88      # km, mean cloud height

CC_ALL     = 0.57      # combined all-cloud cover fraction

# Pressure (bar) at the bottom of each layer. Rough hydrostatic mapping of the
# paper's cloud heights (US standard atmosphere): p ~ p0 * exp(-z/H), H~8 km.
P0          = 1.013                       # surface pressure [bar]
SCALE_H_KM  = 8.0
def p_at_km(z_km):
    return P0 * np.exp(-z_km / SCALE_H_KM)

P_GASTOP    = 1e-3                        # top of atmosphere
P_ICE_BASE  = p_at_km(ICE_H_KM)          # ~0.62 bar
P_LIQ_BASE  = p_at_km(LIQ_H_KM)          # ~0.83 bar
P_SURFACE   = P0                          # bottom of gasbelow

SURF_ALBEDO = 0.05                        # dark Lambertian land/ocean proxy

# Path to the T-matrix ice expansion-coefficient file (single wavelength).
ICE_TMATRIX_COEFFS = os.path.join(REPO, "tmatrix_ice", "ice_oblate_0.5um.coeffs")

# =========================================================================
# 1b.  NUMERICS  --  the single most important practical knob
# =========================================================================
# Large cloud droplets give a strongly forward-peaked phase function whose
# Mie expansion needs hundreds of terms (ncoefs ~ 440-900 for r_eff=8.6 um).
# The doubling-adding solver can only represent ~2*nmug terms, so we apply
# delta-M truncation (Wiscombe 1977; vector form below) to order M_TRUNC.
# For numerical stability the truncation fraction f must stay < ~0.2, which
# requires M_TRUNC >= ~90 for r_eff=8.6 um, hence nmug >= ~45 (since the
# quadrature must satisfy M_TRUNC <= 2*nmug). That is accurate but heavy
# (~1-2 min per wavelength) -- the OFFLINE_* defaults below.
#
# The DEMO_* values use a smaller r_eff that converges in ~10 s/wavelength,
# for quick end-to-end runs / CI / sandboxes.  Switch with MODE.
MODE = os.environ.get("EWIC_MODE", "demo")    # "demo" or "offline"

if MODE == "offline":
    LIQ_REFF = 8.6     # the paper's exact value
    NMUG     = 60      # DAP Gauss points (>= M_TRUNC/2)
    NMUG_MIE = 120     # Mie Gauss points (resolve the large-particle scattering)
    M_TRUNC  = 120     # delta-M truncation order -> f ~ 0.13 (stable)
else:                  # "demo": fast, converged, slightly smaller droplets
    LIQ_REFF = 3.0     # tractable in a 45 s sandbox; cloudbow still present
    NMUG     = 30
    NMUG_MIE = 100
    M_TRUNC  = 60      # f ~ 0.05 for r_eff=3 um (negligible truncation)


# =========================================================================
# 2.  T-MATRIX SUPPORT  (now part of the PyMieDAP package)
# =========================================================================
# The T-matrix .coeffs -> PyMieDAP converter and the vector delta-M truncation
# used below are first-class functions in pymiedap.tmatrix.  See that module's
# docstring for the two file formats and the delta-M derivation.
from pymiedap.tmatrix import (              # noqa: E402
    tmatrix_to_pymiedap_coeffs,
    delta_m_truncate,
)


# =========================================================================
# 3.  BUILD THE EARTH-LIKE CLOUDY MODEL
# =========================================================================
def build_cloudy_model(wvl_list):
    """Two cloud layers (ice over liquid) embedded in a gas column."""
    m = pmd.Model(wvl_list=np.asarray(wvl_list, dtype=float))
    m.asurf = SURF_ALBEDO

    L = m.layers
    nwvl = len(m.wvl_list)

    # ---- gas layers (Rayleigh only) ----
    L.gastop.press   = P_GASTOP
    L.gastop.tau     = [0.0]
    L.gasbelow.press = P_SURFACE
    L.gasbelow.tau   = [0.0]

    # The default Layers() has a single 'cloud' layer + a 'haze' layer.
    # Repurpose: haze -> ICE cloud (loaded coeffs), cloud -> LIQUID cloud (Mie).
    # -------------------- ICE cloud (top) --------------------
    ice = L.haze
    ice.press = P_ICE_BASE
    ice.tau   = [ICE_TAU]
    ice.aerosols.typ = 'I'
    ice.aerosols.layered = False          # value irrelevant; Mie is bypassed
    # -------------------- LIQUID cloud (bottom) --------------
    liq = L.cloud
    liq.press = P_LIQ_BASE
    liq.tau   = [LIQ_TAU]
    a = liq.aerosols
    a.typ   = 'C'
    a.psd   = '2'                         # two-parameter gamma (reff, veff)
    a.r_eff = LIQ_REFF
    a.v_eff = LIQ_VEFF
    a.nr    = [LIQ_NR] * max(2, nwvl)
    a.ni    = [LIQ_NI] * max(2, nwvl)
    a.layered = False

    # propagate wvl-dependent arrays
    m.wvl_list = m.wvl_list
    return m


def build_clear_model(wvl_list):
    """Cloud-free companion column (gas + surface only)."""
    m = pmd.Model(wvl_list=np.asarray(wvl_list, dtype=float))
    m.asurf = SURF_ALBEDO
    L = m.layers
    L.gastop.press   = P_GASTOP;   L.gastop.tau   = [0.0]
    L.haze.press     = P_ICE_BASE; L.haze.tau     = [0.0]
    L.cloud.press    = P_LIQ_BASE; L.cloud.tau    = [0.0]
    L.gasbelow.press = P_SURFACE;  L.gasbelow.tau = [0.0]
    m.wvl_list = m.wvl_list
    return m


# =========================================================================
# 4.  CUSTOM COMPUTE  (Mie+delta-M for liquid, loaded T-matrix coeffs for ice)
# =========================================================================
# We cannot use pmd.compute_model directly: it runs Mie on EVERY aerosol and
# would overwrite the ice coefficients we load.  We reproduce its steps but
# inject the ice scattering matrix by hand and apply delta-M to the liquid
# cloud so the doubling-adding stays stable for large droplets.
# --------------------------------------------------------------------------
def compute_cloudy(model, ice_coeff_file, ice_albedo, output_name,
                   nmug=NMUG, nmug_mie=NMUG_MIE, nsubr=50, nmat=4,
                   m_trunc=M_TRUNC,
                   path_input=os.path.join(REPO, "dap_database")):
    wvl = model.wvl_list
    nwvl = len(wvl)

    for lay_name, layer in vars(model.layers).items():
        if hasattr(layer, 'mixed_aerosols'):
            del layer.mixed_aerosols

        if lay_name == 'haze':                       # ICE: load coeffs
            ice = layer.aerosols
            # one .coeffs file per wavelength; here the same 0.5um file is
            # reused for every wavelength (grey-ice approximation, see header).
            ice.load_coefs([ice_coeff_file] * nwvl)
            # mix_aerosols needs cross sections; only the RATIO (=albedo)
            # matters once tau is user-set, so use unit extinction.
            ice.sext  = np.ones(nwvl)
            ice.ssca  = np.ones(nwvl) * ice_albedo
            ice.ssalb = np.ones(nwvl) * ice_albedo
            ice.f     = 1.0
            layer.mix_aerosols()
        elif lay_name == 'cloud':                    # LIQUID: Mie + delta-M
            a = layer.aerosols
            pmd.mie_code(a, wvl, ngaur=nmug_mie, nsubr=nsubr)
            tau_scale = delta_m_truncate(a, m_trunc)
            layer.tau = np.atleast_1d(layer.tau).astype(float)
            if layer.tau.size != nwvl:
                layer.tau = layer.tau[0] * np.ones(nwvl)
            layer.tau = layer.tau * tau_scale         # conserve scaled energy
            layer.mix_aerosols()
        else:                                        # transparent gas layers
            for aero_name, aero in vars(layer).items():
                if isinstance(aero, pmd.Aerosols):
                    pmd.mie_code(aero, wvl, ngaur=20, nsubr=nsubr)
            layer.mix_aerosols()

    pmd.dap_code(model, rename=True, output_name=output_name,
                 nmug=nmug, nmat=nmat, path_output=path_input)


# =========================================================================
# 5.  DRIVERS  (spectrum + phase curve)
# =========================================================================
def run_phase_curve(ice_coeff_file, ice_albedo,
                    lam_um=0.5, alphas=None):
    """Disk-integrated reflectance + polarization vs phase angle at one lambda."""
    if alphas is None:
        alphas = np.arange(0.0, 170.0, 10.0)

    cloudy = build_cloudy_model([lam_um])
    clear  = build_clear_model([lam_um])

    # Pre-compute the cloudy column with the injected ice coefficients and
    # set cloudy.name, so planet_integrated(force=False) will NOT recompute it
    # (which would overwrite the ice scattering matrix via Mie). The clear
    # column is left fresh for planet_integrated to compute itself.
    compute_cloudy(cloudy, ice_coeff_file, ice_albedo, output_name='earth_cloudy')

    # planet_integrated stores results on the FIRST model and returns None.
    # patchy=False -> deterministic uniform fractional cloud cover, so the
    # result is reproducible (no random cloud realisation) and consistent
    # across wavelengths/phase angles.  Use patchy=True only to sample the
    # cloud-position variability (the paper's 1-sigma bands).
    pmd.planet_integrated([cloudy, clear], alpha=list(alphas),
                          force=False, rename=True, nmug=NMUG,
                          output_names=['earth_cloudy', 'earth_clear'],
                          fclouds=[CC_ALL, 1.0 - CC_ALL], patchy=False)
    I, Q, U, P = cloudy.I[0], cloudy.Q[0], cloudy.U[0], cloudy.P[0]
    return np.asarray(alphas), I, Q, U, P


def run_spectrum(ice_coeff_file, ice_albedo,
                 wvls=None, alpha=90.0):
    """Disk-integrated reflectance + polarization vs wavelength at one phase."""
    if wvls is None:
        wvls = np.arange(0.40, 0.901, 0.05)     # um, VIS-NIR

    cloudy = build_cloudy_model(wvls)
    clear  = build_clear_model(wvls)

    compute_cloudy(cloudy, ice_coeff_file, ice_albedo, output_name='earth_cloudy_spec')

    # patchy=False -> deterministic uniform fractional cloud cover (see
    # run_phase_curve note), giving a smooth, reproducible spectrum.
    pmd.planet_integrated([cloudy, clear], alpha=[alpha],
                          force=False, rename=True, nmug=NMUG,
                          output_names=['earth_cloudy_spec', 'earth_clear_spec'],
                          fclouds=[CC_ALL, 1.0 - CC_ALL], patchy=False)
    I, Q, U, P = cloudy.I[:, 0], cloudy.Q[:, 0], cloudy.U[:, 0], cloudy.P[:, 0]
    return np.asarray(wvls), I, Q, U, P


# =========================================================================
# 6.  MAIN
# =========================================================================
def main():
    # ---- prepare the ice coefficient file in PyMieDAP format ----
    converted = os.path.join(HERE, "ice_oblate_0.5um_pymiedap.coeffs")
    albedo, nrows, lmax = tmatrix_to_pymiedap_coeffs(
        ICE_TMATRIX_COEFFS, converted,
        lam_um=0.5, nr=1.3117, ni=1.0e-8, reff=1.0, veff=0.1)
    print("Converted ice coeffs: albedo={:.6f}, {} orders (Lmax decl. {})"
          .format(albedo, nrows, lmax))

    # ---- phase curve at 0.5 um (T-matrix ice physically valid here) ----
    print("\n=== Phase curve @ 0.5 um ===")
    a, Ia, Qa, Ua, Pa = run_phase_curve(converted, albedo, lam_um=0.5)
    np.savetxt(os.path.join(HERE, "earthlike_phasecurve_0.5um.csv"),
               np.column_stack([a, Ia, Qa, Ua, Pa]),
               header="alpha_deg, I(reflectance), Q, U, P=-Q/I", delimiter=",")

    # ---- spectrum at alpha = 90 deg (grey-ice approximation, see header) ----
    print("\n=== Spectrum @ alpha = 90 deg ===")
    w, Iw, Qw, Uw, Pw = run_spectrum(converted, albedo, alpha=90.0)
    np.savetxt(os.path.join(HERE, "earthlike_spectrum.csv"),
               np.column_stack([w, Iw, Qw, Uw, Pw]),
               header="wavelength_um, I(reflectance), Q, U, P=-Q/I", delimiter=",")

    # ---- optional plots ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(a, Ia, 'o-'); ax[0].set_xlabel("phase angle [deg]")
        ax[0].set_ylabel("reflectance I"); ax[0].set_title("Phase curve @ 0.5 um")
        ax[1].plot(a, 100*Pa, 's-', color='C3'); ax[1].set_xlabel("phase angle [deg]")
        ax[1].set_ylabel("P = -Q/I [%]"); ax[1].set_title("Polarization")
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "earthlike_phasecurve_0.5um.png"), dpi=130)

        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(w*1000, Iw, 'o-'); ax[0].set_xlabel("wavelength [nm]")
        ax[0].set_ylabel("reflectance I"); ax[0].set_title("Spectrum @ a=90 deg")
        ax[1].plot(w*1000, 100*Pw, 's-', color='C3'); ax[1].set_xlabel("wavelength [nm]")
        ax[1].set_ylabel("P = -Q/I [%]"); ax[1].set_title("Polarization")
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "earthlike_spectrum.png"), dpi=130)
        print("\nWrote CSVs and PNGs to", HERE)
    except Exception as e:
        print("Plotting skipped:", e)


if __name__ == "__main__":
    main()
