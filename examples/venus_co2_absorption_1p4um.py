#!/usr/bin/env python3
"""
Venus CO2 absorption spectrum — 1.40–1.50 µm
=============================================

Simulates the apparent reflected-light spectrum of Venus as seen by a
nadir-viewing satellite using the Hansen-Hovenier viewing geometry:

    SZA = 45°,  emission angle = 0°,  Venus cloud albedo ≈ 0.75

The CO2 absorption is computed line-by-line (LBL) using HITRAN 2020 line
parameters retrieved via HAPI.  The vertical integration follows the
Ignatiev et al. Venus mesosphere profile (57–88 km, above cloud top).

Physics
-------
The two-stream reflectance approximation for a Lambert cloud deck beneath
an absorbing gas column:

    I(λ) / F = A_cloud × exp[-τ_CO₂(λ) × (1/μ₀ + 1/μ)]

where τ_CO₂(λ) is the CO₂ column optical depth summed from 57 to 88 km,
μ₀ = cos(SZA) = 0.707, and μ = cos(0°) = 1.

Output
------
Saves ``venus_co2_absorption_1p4um.png`` in the same directory.

Dependencies
------------
  pip install hitran-api numpy scipy matplotlib
"""

from __future__ import annotations

import os
import sys
import warnings
warnings.filterwarnings('ignore')

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.ndimage import gaussian_filter1d

import hapi as hp
from pymiedap.ckdistribution.constants import N_A, G_VENUS, BAR_TO_ATM

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CACHE_DIR  = os.path.join(os.path.expanduser('~'), '.pymiedap', 'hitran')
TABLE_NAME = 'CO2_14um'           # HAPI table name for the 1.4 µm region
VMR_CO2    = 0.965                # Venus CO2 VMR
M_CO2      = 44.01e-3             # kg/mol
m_CO2      = M_CO2 / N_A         # kg/molecule

# Viewing geometry (Hansen-Hovenier style)
SZA_DEG  = 45.0
EMI_DEG  =  0.0
mu0      = np.cos(np.radians(SZA_DEG))   # cos solar zenith angle
mu       = np.cos(np.radians(EMI_DEG))   # cos emission angle
A_CLOUD  = 0.75                          # Venus cloud SWIR reflectance

# Spectral grid
WVN_MIN  = 6650.0    # cm-1
WVN_MAX  = 7151.0    # cm-1
WVN_STEP = 0.02      # cm-1 (fine enough to resolve individual lines)

# Instrument resolution for convolved spectrum
IRF_FWHM_NM = 5.0   # nm FWHM (typical SWIR spectrometer)

# Venus atmospheric levels above cloud deck (Ignatiev et al. profile)
#   columns:  altitude [km],  pressure [bar],  temperature [K]
LEVELS = np.array([
    [57.0,   0.0500,  232.0],
    [61.0,   0.0200,  214.0],
    [65.0,   0.0085,  205.0],
    [69.0,   0.0035,  208.0],
    [73.0,   0.0014,  213.0],
    [77.0,   0.0006,  218.0],
    [82.0,   0.00022, 212.0],
    [88.0,   0.00008, 204.0],
])

# ---------------------------------------------------------------------------
# Step 1: Fetch HITRAN line list (downloads once, then cached)
# ---------------------------------------------------------------------------
os.makedirs(CACHE_DIR, exist_ok=True)
hp.db_begin(CACHE_DIR)

data_file = os.path.join(CACHE_DIR, TABLE_NAME + '.data')
if not os.path.isfile(data_file):
    print(f"Downloading CO2 lines (6580–7220 cm⁻¹) → {data_file} ...")
    # All major CO2 isotopologues (global HITRAN IDs)
    iso_ids = [7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122]
    hp.fetch_by_ids(TABLE_NAME, iso_ids, WVN_MIN - 70., WVN_MAX + 70.)
else:
    print(f"Using cached line list: {data_file}")
    known = hp.tableList() if 'tableList' in dir(hp) else []
    if TABLE_NAME not in known:
        hp.db_begin(CACHE_DIR)

# ---------------------------------------------------------------------------
# Step 2: LBL cross-sections and optical depth
# ---------------------------------------------------------------------------
nu_grid   = np.arange(WVN_MIN, WVN_MAX, WVN_STEP)
wvl_nm    = 1.0e7 / nu_grid           # nm (decreasing)

tau_total    = np.zeros_like(nu_grid)
sigma_list   = []

print(f"\nComputing LBL CO₂ optical depth ({len(LEVELS)-1} layers):")
for i in range(len(LEVELS) - 1):
    alt_bot, P_bot, T_bot = LEVELS[i]
    alt_top, P_top, T_top = LEVELS[i + 1]

    P_avg = 0.5 * (P_bot + P_top)
    T_avg = 0.5 * (T_bot + T_top)
    p_atm = P_avg * BAR_TO_ATM

    print(f"  Layer {i+1}  alt={alt_bot:.0f}–{alt_top:.0f} km  "
          f"P={P_avg:.4f} bar  T={T_avg:.0f} K", end=' ', flush=True)

    nu_h, coef = hp.absorptionCoefficient_Voigt(
        SourceTables=TABLE_NAME,
        Environment={'T': float(T_avg), 'p': float(p_atm)},
        OmegaRange=[float(WVN_MIN - 5.), float(WVN_MAX + 5.)],
        OmegaStep=float(WVN_STEP),
        OmegaWingHW=20.0,
        GammaL='gamma_air',
        HITRAN_units=True,          # cm² molecule⁻¹
    )
    sigma_interp = np.interp(nu_grid, np.asarray(nu_h), np.asarray(coef), left=0., right=0.)
    sigma_list.append(sigma_interp.copy())

    # Hydrostatic column density [molecules m⁻²]
    dP_Pa = (P_bot - P_top) * 1.0e5
    N_d   = dP_Pa * VMR_CO2 / (m_CO2 * G_VENUS)

    tau_layer  = sigma_interp * N_d * 1.0e-4   # σ[cm²] × N[m⁻²] × 1e-4
    tau_total += tau_layer
    print(f"  τ_peak = {tau_layer.max():.2f}")

sigma_arr = np.array(sigma_list)   # shape (nlayer, nnu)

# ---------------------------------------------------------------------------
# Step 3: Apparent reflectance (two-stream approximation)
# ---------------------------------------------------------------------------
exponent    = tau_total * (1.0 / mu0 + 1.0 / mu)
I_native    = A_CLOUD * np.exp(-exponent)

# Convolve with a Gaussian IRF (simulate finite instrument resolution)
fwhm_cm_at_centre = 1.0e4 * (IRF_FWHM_NM * 1.0e-3) / (1.45 ** 2)
sigma_px    = (fwhm_cm_at_centre / WVN_STEP) / (2.0 * np.sqrt(2.0 * np.log(2.0)))
I_convolved = gaussian_filter1d(I_native, sigma=sigma_px)
I_clear     = np.full_like(nu_grid, A_CLOUD)

# ---------------------------------------------------------------------------
# Step 4: Three-panel figure
# ---------------------------------------------------------------------------
GREY  = '#2d2d2d'
BLUE  = '#1f5fa6'
RED   = '#c0392b'
ORNG  = '#e67e22'
GRN   = '#27ae60'

fig, axes = plt.subplots(3, 1, figsize=(11, 11),
                          gridspec_kw={'hspace': 0.48})

# --- Ascending wavelength for all panels ---
order    = np.argsort(wvl_nm)
wvl_plot = wvl_nm[order]
nu_plot  = nu_grid[order]

# Panel 1: Cross-section at cloud top
ax1 = axes[0]
ax1.semilogy(wvl_plot, sigma_arr[0][order], lw=0.5, color=BLUE, alpha=0.85)
ax1.set_ylabel(r'$\sigma_\mathrm{CO_2}$ [cm$^2$ molecule$^{-1}$]', fontsize=11)
ax1.set_title(
    r'CO$_2$ absorption cross-section  ($P = 0.05$ bar, $T = 232$ K, cloud top $\approx 57$ km)',
    fontsize=10.5
)
ax1.set_xlim(1400, 1505)
ax1.grid(True, which='both', alpha=0.2)
ax1.tick_params(labelbottom=False)

# Secondary x-axis in wavenumber
ax1t = ax1.twiny()
ax1t.set_xlim(1.0e7 / 1505., 1.0e7 / 1400.)
ax1t.set_xlabel(r'Wavenumber [cm$^{-1}$]', fontsize=9, labelpad=3)
ax1t.tick_params(labelsize=8)

# Panel 2: Column optical depth
ax2 = axes[1]
ax2.semilogy(wvl_plot, np.clip(tau_total[order], 1e-5, None),
             lw=0.7, color=GREY, alpha=0.9)
ax2.axhline(1.0, ls='--', lw=1.3, color=RED,  alpha=0.75, label=r'$\tau = 1$')
ax2.axhline(0.1, ls=':',  lw=1.3, color=ORNG, alpha=0.75, label=r'$\tau = 0.1$')
ax2.set_ylabel(r'CO$_2$ column optical depth $\tau(\lambda)$', fontsize=11)
ax2.set_title(
    r'Column optical depth — Venus mesosphere, 57–88 km  (CO$_2$ VMR = 96.5 %)',
    fontsize=10.5
)
ax2.set_xlim(1400, 1505)
ax2.set_ylim(1e-4, 3e3)
ax2.legend(fontsize=9, loc='upper left')
ax2.grid(True, which='both', alpha=0.2)
ax2.tick_params(labelbottom=False)

# Mark band positions
band_labels = {1435.: r'$3\nu_1\!+\!\nu_3$', 1469.: r'$2\nu_1\!+\!2\nu_2\!+\!\nu_3$'}
for wl, lbl in band_labels.items():
    if 1400 < wl < 1505:
        ax2.axvline(wl, ls=':', lw=1.0, color='purple', alpha=0.55)
        ax2.text(wl + 0.8, 500., lbl, color='purple', fontsize=7.5, va='top', rotation=90)
        axes[2].axvline(wl, ls=':', lw=1.0, color='purple', alpha=0.35)

# Panel 3: Apparent reflectance
ax3 = axes[2]
ax3.fill_between(wvl_plot, I_native[order], I_clear[order],
                 alpha=0.15, color=BLUE, label=r'CO$_2$ absorption')
ax3.plot(wvl_plot, I_clear[order],    ls='--', lw=1.5, color=GRN,  alpha=0.8,
         label='Clear sky (no gas)')
ax3.plot(wvl_plot, I_native[order],   lw=0.5,  color=GREY, alpha=0.5,
         label='LBL (native res.)')
ax3.plot(wvl_plot, I_convolved[order], lw=1.8,  color=BLUE,
         label=rf'Convolved ($\Delta\lambda \approx {IRF_FWHM_NM:.0f}$ nm)')
ax3.set_xlabel(r'Wavelength [nm]', fontsize=12)
ax3.set_ylabel(r'Normalised reflectance $I/F$', fontsize=11)
ax3.set_title(
    r'Apparent reflectance — Venus satellite view'
    r'  ($\theta_0 = 45°$, nadir,  $A_\mathrm{cloud} = 0.75$)',
    fontsize=10.5
)
ax3.set_xlim(1400, 1505)
ax3.set_ylim(0, 0.82)
ax3.legend(fontsize=9, loc='lower left', ncol=2)
ax3.grid(True, alpha=0.25)

fig.text(
    0.5, 0.01,
    r'Hansen-Hovenier geometry  ·  HITRAN 2020 CO$_2$ lines via HAPI  ·  pymiedap.ckdistribution',
    ha='center', fontsize=8, color='grey', style='italic'
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'venus_co2_absorption_1p4um.png')
plt.savefig(OUT, dpi=180, bbox_inches='tight')
print(f'\nSaved → {OUT}')
