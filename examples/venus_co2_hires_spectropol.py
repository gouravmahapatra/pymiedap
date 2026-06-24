"""
Venus CO2 absorption spectrum at 1 nm spectropolarimeter resolution
1420–1460 nm  |  Hansen-Hovenier viewing geometry  |  SZA = 45°, nadir
"""
import sys, os, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/sessions/eloquent-admiring-ptolemy/mnt/pymiedap')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from scipy.ndimage import gaussian_filter1d
import hapi as hp
from pymiedap.ckdistribution.constants import N_A, G_VENUS, BAR_TO_ATM

# ─── Setup ──────────────────────────────────────────────────────────────────
CACHE   = '/tmp/hitran_test'
TABLE   = 'CO2_14'          # already cached from previous run

VMR_CO2 = 0.965
M_CO2   = 44.01e-3
m_CO2   = M_CO2 / N_A

# Geometry: satellite nadir, SZA = 45°
MU0     = np.cos(np.radians(45.))
MU      = 1.0
PHASE   = 45.                           # degrees
A_CLOUD = 0.75

# Spectral range  (1420–1460 nm  →  6849–7042 cm-1)
WVL_LO_NM  = 1420.
WVL_HI_NM  = 1460.
WVN_LO     = 1e7 / WVL_HI_NM           # 6849.3 cm-1
WVN_HI     = 1e7 / WVL_LO_NM           # 7042.3 cm-1
WVN_STEP   = 0.005                      # cm-1 (ultra-fine LBL grid)

# Instrument: 1 nm FWHM
IRF_NM     = 1.0

# Venus levels above cloud top (57 – 88 km, Ignatiev profile)
LEVELS = np.array([
    [57.,  0.0500,  232.],
    [61.,  0.0200,  214.],
    [65.,  0.0085,  205.],
    [69.,  0.0035,  208.],
    [73.,  0.0014,  213.],
    [77.,  0.0006,  218.],
    [82.,  0.00022, 212.],
    [88.,  0.00008, 204.],
])

# CO2 vibrational band assignments in 1420–1460 nm region
# (combination / overtone bands; labels follow HITRAN notation)
BANDS = {
    1427.6:  r'$3\nu_1\!+\!\nu_3$  (main)',
    1432.5:  r'hot band',
    1438.4:  r'$3\nu_1\!+\!\nu_3$  (²⁰C¹⁸O)',
    1443.7:  r'$2\nu_1\!+\!2\nu_2\!+\!\nu_3$',
    1453.2:  r'hot band',
}

hp.db_begin(CACHE)

# ─── LBL cross-section at each layer ────────────────────────────────────────
nu_grid   = np.arange(WVN_LO - 2., WVN_HI + 2., WVN_STEP)
wvl_all   = 1.e7 / nu_grid                 # nm, decreasing

tau_total  = np.zeros_like(nu_grid)
sig_levels = []

print(f"Grid: {len(nu_grid):,} pts  Δν = {WVN_STEP} cm⁻¹  "
      f"({WVL_LO_NM:.0f}–{WVL_HI_NM:.0f} nm)\n")

for i in range(len(LEVELS) - 1):
    alt_b, P_b, T_b = LEVELS[i]
    alt_t, P_t, T_t = LEVELS[i+1]
    P_avg = 0.5 * (P_b + P_t);  T_avg = 0.5 * (T_b + T_t)
    p_atm = P_avg * BAR_TO_ATM

    print(f"  Layer {i+1}:  {alt_b:.0f}–{alt_t:.0f} km  "
          f"P={P_avg:.4f} bar  T={T_avg:.0f} K", end='  ', flush=True)

    nu_h, coef = hp.absorptionCoefficient_Voigt(
        SourceTables=TABLE,
        Environment={'T': float(T_avg), 'p': float(p_atm)},
        OmegaRange=[float(WVN_LO - 3.), float(WVN_HI + 3.)],
        OmegaStep=float(WVN_STEP),
        OmegaWingHW=25.0,
        GammaL='gamma_air',
        HITRAN_units=True,
    )
    sigma = np.interp(nu_grid, np.asarray(nu_h), np.asarray(coef), left=0., right=0.)
    sig_levels.append(sigma.copy())

    dP_Pa = (P_b - P_t) * 1.e5
    Nd    = dP_Pa * VMR_CO2 / (m_CO2 * G_VENUS)
    tau_l = sigma * Nd * 1.e-4
    tau_total += tau_l
    print(f"τ_peak = {tau_l.max():.2f}")

sig_levels = np.array(sig_levels)

# ─── Clip to exact wavelength window ────────────────────────────────────────
mask   = (wvl_all >= WVL_LO_NM) & (wvl_all <= WVL_HI_NM)
nu_w   = nu_grid[mask]
wvl_w  = wvl_all[mask]                     # nm, decreasing
tau_w  = tau_total[mask]
sig_w  = sig_levels[:, mask]

# sort ascending in wavelength for plotting
idx    = np.argsort(wvl_w)
wvl    = wvl_w[idx]                        # ascending nm
nu_s   = nu_w[idx]                         # cm-1
tau    = tau_w[idx]
sigma0 = sig_w[0][idx]                     # cloud-top level

# ─── Reflectance ────────────────────────────────────────────────────────────
expo       = tau * (1./MU0 + 1./MU)
I_native   = A_CLOUD * np.exp(-expo)
I_clear    = np.full_like(wvl, A_CLOUD)

# 1 nm Gaussian IRF
fwhm_cm    = 1.e4 * (IRF_NM * 1.e-3) / (1.44**2)     # ~4.8 cm-1
sig_px     = (fwhm_cm / WVN_STEP) / (2.*np.sqrt(2.*np.log(2.)))
I_conv     = gaussian_filter1d(I_native, sigma=sig_px)

# ─── Stokes Q/I  (spectropolarimeter channel) ────────────────────────────────
# Physical model: two contributions to the observed polarized light:
#   1. Cloud deck (Mie): polarization P_cloud at this phase angle.
#      H&H 1974 Fig 3 gives P ≈ +4–6 % for Venus-like particles (reff~1 µm,
#      veff=0.07, n=1.44) at phase φ = 45°.  We use P_cloud = +0.05.
#   2. CO2 Rayleigh above cloud: depolarization factor ρ_CO2 = 0.0802
#      (Sneep & Ubachs 2005).  Single-scatter polarisation:
#      P_Ray(φ) = -(1-ρ)sin²φ / [(1-ρ) + (1+ρ)cos²φ]
#      At φ=45°: P_Ray = -0.920×0.5 / (0.920 + 1.080×0.5) = -0.315
#      Rayleigh single-scatter optical depth above cloud (P<0.05 bar):
#      τ_Ray(λ) = τ_Ray_ref × (550 nm / λ)^4  ;  τ_Ray_ref(550 nm) ≈ 0.04
#      (for CO2 column with P < 0.05 bar)
#   Composite Q/I  (linear combination weighted by emerging intensities):
#      Q = [P_cloud × I_cloud + P_Ray × I_Ray] / (I_cloud + I_Ray)
P_CLOUD = 0.050
RHO_CO2 = 0.0802
P_RAY   = -(1-RHO_CO2)*np.sin(np.radians(PHASE))**2 / \
           ((1-RHO_CO2) + (1+RHO_CO2)*np.cos(np.radians(PHASE))**2)
print(f"\nP_Rayleigh at {PHASE}°: {P_RAY:.4f}  ({100*P_RAY:.1f} %)")

TAU_RAY_550 = 0.038     # CO2 Rayleigh OD above cloud at 550 nm
tau_ray     = TAU_RAY_550 * (0.550 / (wvl * 1.e-3))**4   # at each λ

# Cloud contribution (transmission through CO2 column)
I_cloud_frac = np.exp(-(tau / MU0 + tau / MU))          # double-pass
# Rayleigh single-scatter above cloud (upward only, sun-side)
I_ray_frac   = 0.5 * (1. - np.exp(-tau_ray / MU0))      # weak

# Weighted Q/I
I_total_model = I_cloud_frac + I_ray_frac
Q_total_model = P_CLOUD * I_cloud_frac + P_RAY * I_ray_frac
QoI_native    = Q_total_model / np.maximum(I_total_model, 1.e-30)
QoI_conv      = gaussian_filter1d(QoI_native, sigma=sig_px)

# Also express as percentage
QoI_pct       = 100. * QoI_conv

# ─── Figure ─────────────────────────────────────────────────────────────────
BLUE  = '#1f5fa6'
DKBL  = '#0a3060'
RED   = '#c0392b'
ORNG  = '#e06c00'
GRN   = '#1a7a4a'
PURP  = '#6c2b8f'
GREY  = '#2d2d2d'
LG    = '#bbbbbb'

fig   = plt.figure(figsize=(12, 13))
gs    = fig.add_gridspec(4, 1, hspace=0.52, top=0.95, bottom=0.06)
axes  = [fig.add_subplot(g) for g in gs]

# ── Panel 1: LBL cross-section (cloud-top level) ──────────────────────────
ax = axes[0]
ax.semilogy(wvl, sigma0, lw=0.35, color=BLUE, alpha=0.85)
ax.set_ylabel(r'$\sigma_\mathrm{CO_2}$ [cm$^2$/mol]', fontsize=10.5)
ax.set_title(
    r'CO$_2$ LBL cross-section  ($P\!=\!0.05$ bar, $T\!=\!232$ K, cloud top $\!\approx\!57$ km)',
    fontsize=10)
ax.set_xlim(WVL_LO_NM, WVL_HI_NM)
ax.set_ylim(1e-30, 1e-22)
ax.grid(True, which='both', alpha=0.18, lw=0.5)
ax.tick_params(labelbottom=False)
# secondary wavenumber axis
at = ax.twiny()
at.set_xlim(1e7/WVL_HI_NM, 1e7/WVL_LO_NM)
at.set_xlabel(r'Wavenumber [cm$^{-1}$]', fontsize=8.5, labelpad=2)
at.tick_params(labelsize=8)

# ── Panel 2: column optical depth with band annotations ───────────────────
ax = axes[1]
tau_clipped = np.clip(tau, 1e-5, 5e3)
ax.fill_between(wvl, tau_clipped, 1e-5, alpha=0.10, color=BLUE)
ax.semilogy(wvl, tau_clipped, lw=0.55, color=GREY, alpha=0.9)
ax.axhline(1.0, ls='--', lw=1.3, color=RED,  alpha=0.8, label=r'$\tau=1$')
ax.axhline(0.1, ls=':',  lw=1.3, color=ORNG, alpha=0.8, label=r'$\tau=0.1$')
ax.set_ylabel(r'$\tau_\mathrm{CO_2}(\lambda)$   [column, 57–88 km]', fontsize=10.5)
ax.set_title(r'CO$_2$ column optical depth — band assignments', fontsize=10)
ax.set_xlim(WVL_LO_NM, WVL_HI_NM)
ax.set_ylim(1e-4, 5e3)
ax.legend(fontsize=9, loc='upper right')
ax.grid(True, which='both', alpha=0.18, lw=0.5)
ax.tick_params(labelbottom=False)

band_items = sorted(BANDS.items())
for k, (wl0, lbl) in enumerate(band_items):
    if WVL_LO_NM <= wl0 <= WVL_HI_NM:
        ax.axvline(wl0, ls=':', lw=0.9, color=PURP, alpha=0.7)
        ypos = [800., 200., 50., 20., 5.][k % 5]
        ax.text(wl0 + 0.25, ypos, lbl, color=PURP, fontsize=6.8,
                va='center', rotation=90, alpha=0.9)

# ── Panel 3: I/F  (intensity channel) ────────────────────────────────────
ax = axes[2]
ax.fill_between(wvl, I_native, I_clear, alpha=0.12, color=BLUE,
                label=r'CO$_2$ absorption')
ax.plot(wvl, I_clear,   ls='--', lw=1.4, color=GRN,  alpha=0.8,
        label='Clear sky  (no CO$_2$)')
ax.plot(wvl, I_native,  lw=0.4, color=LG, alpha=0.7,
        label='LBL  (Δν = 0.005 cm⁻¹)')
ax.plot(wvl, I_conv,    lw=2.0, color=BLUE,
        label=r'Convolved  ($\Delta\lambda = 1$ nm)')
# Shade the narrow spectral windows where τ < 1
win_mask = tau < 1.0
if win_mask.any():
    ax.fill_between(wvl, 0, I_conv, where=win_mask,
                    alpha=0.25, color=GRN, label='Window  ($\\tau < 1$)')
ax.set_ylabel(r'Normalised reflectance  $I/F$', fontsize=10.5)
ax.set_title(
    rf'Intensity channel  —  SZA $= {int(np.degrees(np.arccos(MU0)))}°$, '
    rf'nadir,  $A_\mathrm{{cloud}} = {A_CLOUD}$',
    fontsize=10)
ax.set_xlim(WVL_LO_NM, WVL_HI_NM)
ax.set_ylim(-0.01, 0.80)
ax.legend(fontsize=8.5, loc='upper left', ncol=2)
ax.grid(True, alpha=0.22, lw=0.5)
ax.tick_params(labelbottom=False)
for wl0, _ in band_items:
    if WVL_LO_NM <= wl0 <= WVL_HI_NM:
        ax.axvline(wl0, ls=':', lw=0.7, color=PURP, alpha=0.35)

# ── Panel 4: Q/I  (polarisation channel) ─────────────────────────────────
ax = axes[3]
ax.axhline(100.*P_CLOUD, ls='--', lw=1.2, color=GRN, alpha=0.75,
           label=rf'Cloud continuum  ($P_\mathrm{{cloud}} = {100*P_CLOUD:.0f}$%)')
ax.axhline(100.*P_RAY,   ls=':',  lw=1.2, color=RED, alpha=0.75,
           label=rf'Rayleigh limit  ($P_\mathrm{{Ray}} = {100*P_RAY:.1f}$%)')
ax.plot(wvl, QoI_pct, lw=1.8, color=DKBL,
        label=r'$Q/I$ estimate  (convolved, 1 nm)')
ax.fill_between(wvl, QoI_pct, 100.*P_CLOUD, alpha=0.12, color=DKBL)
ax.set_xlabel(r'Wavelength [nm]', fontsize=12)
ax.set_ylabel(r'$Q/I$  [%]', fontsize=10.5)
ax.set_title(
    r'Polarisation channel  —  H&H cloud ($n\!=\!1.44$, $r_\mathrm{eff}\!=\!1\,\mu$m) '
    r'+ CO$_2$ Rayleigh',
    fontsize=10)
ax.set_xlim(WVL_LO_NM, WVL_HI_NM)
ax.set_ylim(-35., 12.)
ax.legend(fontsize=8.5, loc='lower left', ncol=2)
ax.grid(True, alpha=0.22, lw=0.5)
for wl0, _ in band_items:
    if WVL_LO_NM <= wl0 <= WVL_HI_NM:
        ax.axvline(wl0, ls=':', lw=0.7, color=PURP, alpha=0.35)

# ── Shared band label strip between panels 3 and 4 ───────────────────────
for wl0, lbl in band_items:
    if WVL_LO_NM <= wl0 <= WVL_HI_NM:
        axes[2].annotate('', xy=(wl0, -0.01), xytext=(wl0, 0.03),
                         xycoords='data', textcoords='data',
                         arrowprops=dict(arrowstyle='->', color=PURP, lw=0.8))

# ── Super-title and footer ────────────────────────────────────────────────
fig.suptitle(
    r'Venus CO$_2$ absorption spectrum  ·  1420–1460 nm  ·  1 nm spectropolarimeter',
    fontsize=13, fontweight='bold', y=0.985)
fig.text(
    0.5, 0.015,
    r'Hansen-Hovenier geometry  ·  HITRAN 2020 via HAPI  ·  Ignatiev et al. (2009) P–T profile  ·  pymiedap.ckdistribution',
    ha='center', fontsize=7.5, color='grey', style='italic')

OUT = '/sessions/eloquent-admiring-ptolemy/mnt/pymiedap/examples/venus_co2_hires_spectropol.png'
plt.savefig(OUT, dpi=200, bbox_inches='tight')
print(f'\nSaved → {OUT}')
