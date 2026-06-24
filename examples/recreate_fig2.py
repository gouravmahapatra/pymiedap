#!/usr/bin/env python3
"""
Recreate Figure 2 from Mahapatra, Rossi & Stam (2024)
"Characterizing Venus's clouds and hazes using CO2 absorption bands in flux and polarization"

Figure 2: Pressure and temperature profiles of the Venus model atmosphere (VIRA / Seiff et al. 1985)
with horizontal layer boundary lines (4 km below 40 km, 1 km above 40 km).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatter
from scipy.interpolate import interp1d

# ----------------------------------------------------------------
# VIRA equatorial profile anchor points (from venus_ignatiev in atmosphere.py)
# alt [km], T [K], P [bar]
# ----------------------------------------------------------------
_ig = np.array([
    [  0.0,  735.3, 92.10],
    [  4.0,  697.4, 66.65],
    [  8.0,  660.4, 47.35],
    [ 12.0,  619.1, 33.04],
    [ 16.0,  574.5, 22.52],
    [ 20.0,  527.4, 14.93],
    [ 24.0,  476.0,  9.573],
    [ 28.0,  427.0,  5.917],
    [ 32.0,  380.1,  3.501],
    [ 36.0,  337.4,  1.979],
    [ 40.0,  299.7,  1.066],
    [ 44.0,  267.0,  0.5356],
    [ 48.0,  238.2,  0.2488],
    [ 52.0,  212.5,  0.1067],
    [ 56.0,  198.8,  4.370e-2],
    [ 60.0,  195.2,  1.768e-2],
    [ 64.0,  203.5,  7.132e-3],
    [ 68.0,  210.6,  2.941e-3],
    [ 72.0,  215.4,  1.199e-3],
    [ 76.0,  218.2,  4.820e-4],
    [ 80.0,  218.5,  1.920e-4],
    [ 84.0,  214.5,  7.526e-5],
    [ 88.0,  206.0,  2.924e-5],
    [ 92.0,  195.5,  1.126e-5],
    [ 96.0,  184.0,  4.289e-6],
    [100.0,  172.0,  1.612e-6],
])

alt_anchor = _ig[:, 0]   # km
T_anchor   = _ig[:, 1]   # K
P_bar      = _ig[:, 2]   # bar → convert to atm
P_atm      = P_bar / 1.01325

# Interpolate to the model layer grid:
# Levels at 0, 4, 8, ... 40 km (4 km spacing), then 41, 42, ... 100 km (1 km spacing)
alt_below40 = np.arange(0,  41, 4)   # 0..40 km, 4 km steps → layer boundaries
alt_above40 = np.arange(40, 101, 1)  # 40..100 km, 1 km steps → layer boundaries
alt_all     = np.unique(np.concatenate([alt_below40, alt_above40]))

# Smooth interpolation for plotting the continuous profile
f_T = interp1d(alt_anchor, T_anchor, kind='cubic')
f_P = interp1d(alt_anchor, np.log10(P_atm), kind='cubic')  # interpolate in log space

alt_fine = np.linspace(0, 100, 2000)
T_fine   = f_T(alt_fine)
P_fine   = 10 ** f_P(alt_fine)   # back to linear atm

# ----------------------------------------------------------------
# Build the figure — single panel with dual x-axis
# ----------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(5.5, 8))

# --- Bottom x-axis: Pressure (log scale) ---
# Pressure decreases with altitude, so we plot P_fine on ax1
p_line, = ax1.semilogx(P_fine, alt_fine, color='black', lw=1.8, label='Pressure')
ax1.set_xlabel('Pressure (atm)', fontsize=12)
ax1.set_xlim(1e2, 1e-6)   # high pressure at left (surface), low at right (top)
ax1.set_ylim(0, 100)
ax1.set_ylabel('Altitude (km)', fontsize=12)

# Set log ticks explicitly
ax1.set_xticks([1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
ax1.xaxis.set_major_formatter(
    matplotlib.ticker.LogFormatterSciNotation(base=10)
)
ax1.tick_params(axis='x', which='both', labelsize=9)
ax1.tick_params(axis='y', labelsize=10)

# --- Top x-axis: Temperature (linear scale) ---
ax2 = ax1.twiny()
t_line, = ax2.plot(T_fine, alt_fine, color='black', lw=1.8, ls='--', label='Temperature')
ax2.set_xlabel('Temperature (K)', fontsize=12)
ax2.set_xlim(150, 800)   # wide enough to show 200–700 K
ax2.set_xticks([200, 300, 400, 500, 600, 700])
ax2.tick_params(axis='x', labelsize=10)

# --- Horizontal layer boundary lines ---
# Below 40 km: every 4 km; above 40 km: every 1 km
# Draw thin grey horizontal lines at each level
for alt_lev in alt_all:
    ax1.axhline(alt_lev, color='grey', lw=0.4, alpha=0.6, zorder=0)

# Also add y-axis ticks every 10 km
ax1.set_yticks(np.arange(0, 101, 10))
ax1.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5))

# --- Legend ---
lines = [p_line, t_line]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper right', fontsize=10, framealpha=0.8)

ax1.grid(False)
ax2.grid(False)

plt.tight_layout()
out_path = '/sessions/eloquent-admiring-ptolemy/mnt/outputs/fig2_reproduction.png'
plt.savefig(out_path, dpi=200, bbox_inches='tight')
print(f"Saved: {out_path}")
plt.close()
