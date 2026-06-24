#!/usr/bin/env python3
"""
Earth SWIR absorption spectrum — CO₂ and CH₄ bands
====================================================

Demonstrates the pymiedap.ckdistribution subpackage without running the full
PyMieDAP radiative-transfer code.  It shows:

  1. How to set up an Earth standard atmosphere.
  2. How to fetch HITRAN line parameters via HAPI.
  3. How to compute absorption cross-sections and k-distributions.
  4. How to plot the resulting band-mean bmabs spectrum.

Running this script requires:
  - hapi (``pip install hitran-api``)
  - numpy, scipy, matplotlib

The HITRAN line files are downloaded on the first run (~50–100 MB) and cached
locally at ``~/.pymiedap/hitran/``.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# ---- Import the subpackage ------------------------------------------------
from pymiedap.ckdistribution import (
    earth_standard,
    gauss_legendre_points,
    fetch_lines,
    compute_cross_section,
    slitfunction,
    kspec_layer,
    compute_bmabs,
    wvl2wvn,
    wvn2wvl,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CACHE_DIR = os.path.join(os.path.expanduser('~'), '.pymiedap', 'hitran')
os.makedirs(CACHE_DIR, exist_ok=True)

# Wavelength grid: CO₂ 1.6 µm band + CH₄ 1.67 µm band
WAV_UM = np.arange(1.56, 1.70, 0.005)   # µm

SIGMA_UM   = 0.005   # spectral window full-width [µm]
TRUNCW_UM  = 0.10    # truncation margin [µm]
NW         = 20000   # wavenumber resolution [pts per µm]
N_GAUSS    = 10      # Gauss–Legendre quadrature points

# Atmosphere
NLEV = 20
atm  = earth_standard(nlev=NLEV)

# ---------------------------------------------------------------------------
# Step 1: fetch HITRAN lines for CO₂ and CH₄
# ---------------------------------------------------------------------------
print("=" * 60)
print("Fetching HITRAN line lists …")
print("=" * 60)

# Wavenumber range from the wavelength grid + truncation
wvn_max = float(wvl2wvn(WAV_UM.min() - 0.5 * SIGMA_UM))
wvn_min = float(wvl2wvn(WAV_UM.max() + 0.5 * SIGMA_UM))

fetch_lines('CO2', wvn_min, wvn_max, CACHE_DIR)
fetch_lines('CH4', wvn_min, wvn_max, CACHE_DIR)
print("Line lists ready.\n")

# ---------------------------------------------------------------------------
# Step 2: plot absorption cross-sections at a single level (P=0.5 bar, T=250 K)
# ---------------------------------------------------------------------------
print("Computing cross-sections at P=0.5 bar, T=250 K …")
wvl_test = 1.60   # µm
nnv, truncv, specv, speci, vmin, vmax, tmin, tmax = slitfunction(
    wvl_test, SIGMA_UM, NW, TRUNCW_UM
)
wvn_step = (vmax - vmin) / max(nnv - 1, 1)

nu_co2, sig_co2 = compute_cross_section(
    'CO2', T=250., P_bar=0.5,
    wvn_min=tmin, wvn_max=tmax, wvn_step=wvn_step,
    cache_dir=CACHE_DIR
)
nu_ch4, sig_ch4 = compute_cross_section(
    'CH4', T=250., P_bar=0.5,
    wvn_min=tmin, wvn_max=tmax, wvn_step=wvn_step,
    cache_dir=CACHE_DIR
)

fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
axes[0].semilogy(wvn2wvl(nu_co2) * 1e3, sig_co2, lw=0.5, color='C0')
axes[0].set_ylabel(r'CO$_2$  $\sigma_{abs}$ [cm$^2$ molecule$^{-1}$]')
axes[0].set_title(r'Absorption cross-sections at $P=0.5$ bar, $T=250$ K')
axes[0].grid(True, alpha=0.3)

axes[1].semilogy(wvn2wvl(nu_ch4) * 1e3, sig_ch4, lw=0.5, color='C1')
axes[1].set_ylabel(r'CH$_4$  $\sigma_{abs}$ [cm$^2$ molecule$^{-1}$]')
axes[1].set_xlabel('Wavelength [nm]')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('swir_cross_sections.pdf', dpi=150)
print("Saved: swir_cross_sections.pdf\n")

# ---------------------------------------------------------------------------
# Step 3: compute bmabs for a single level and plot the CK distribution
# ---------------------------------------------------------------------------
print("Computing k-distribution for CO₂ at P=0.5 bar, T=250 K …")
gp, gw = gauss_legendre_points(N_GAUSS)

# kspec_layer for a single (P, T) level
layer_kdis = kspec_layer(
    WAV_UM, pres_bar=0.5, temp_K=250., gauss_points=gp,
    molecule='CO2', cache_dir=CACHE_DIR,
    sigma_um=SIGMA_UM, truncw_um=TRUNCW_UM, nw=NW
)
# layer_kdis.shape == (nwav, ngauss)

fig, ax = plt.subplots(figsize=(10, 4))
for ig in range(N_GAUSS):
    ax.semilogy(WAV_UM * 1e3, layer_kdis[:, ig],
                lw=1, alpha=0.7, label=f'g={gp[ig]:.2f}')
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel(r'$k_g$ [cm$^2$ molecule$^{-1}$]')
ax.set_title(r'CO$_2$ k-distribution at Gauss points — $P=0.5$ bar, $T=250$ K')
ax.legend(ncol=5, fontsize=7, loc='upper right')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('swir_kdistribution.pdf', dpi=150)
print("Saved: swir_kdistribution.pdf\n")

# ---------------------------------------------------------------------------
# Step 4: full bmabs computation over all levels
# ---------------------------------------------------------------------------
print("Computing bmabs for the full Earth atmosphere (CO₂) …")
print(f"  nlev = {atm.nlev}, nwav = {len(WAV_UM)}, ngauss = {N_GAUSS}")

bmabs_co2 = compute_bmabs(
    atmosphere=atm,
    molecule='CO2',
    wav=WAV_UM,
    gauss_points=gp,
    cache_dir=CACHE_DIR,
    sigma_um=SIGMA_UM,
    truncw_um=TRUNCW_UM,
    nw=NW,
    verbose=True,
)
# bmabs_co2.shape == (nlayer, nwav, ngauss)
print(f"\nbmabs_co2.shape = {bmabs_co2.shape}")

# ---------------------------------------------------------------------------
# Step 5: plot the band-mean bmabs spectrum
# ---------------------------------------------------------------------------
# Gauss-weight-averaged bmabs per wavelength (sum over Gauss points × weight × 0.5)
bmabs_mean = 0.5 * np.sum(bmabs_co2 * gw[np.newaxis, np.newaxis, :], axis=2)
# Sum over layers to get the column bmabs
bmabs_col  = np.sum(bmabs_mean, axis=0)   # shape (nwav,)

fig, ax = plt.subplots(figsize=(10, 4))
ax.semilogy(WAV_UM * 1e3, bmabs_col, lw=1.5, color='C0',
            label=r'CO$_2$ (Earth standard, 20 levels)')
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel(r'Column bmabs (Gauss-averaged)')
ax.set_title(r'CO$_2$ band-mean absorption optical depth — Earth SWIR')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('swir_bmabs_spectrum.pdf', dpi=150)
print("Saved: swir_bmabs_spectrum.pdf")
print("\nDone.")
