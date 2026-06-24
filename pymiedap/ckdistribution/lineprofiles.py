"""
Sub-Lorentzian line shape profiles and χ-factor corrections for CO₂.

Background
----------
The standard Voigt (or Lorentzian) profile significantly overestimates
absorption in the far wings of CO₂ lines when the bath gas pressure exceeds
~1 bar.  Venus's lower atmosphere reaches 92 bar at the surface, making this
correction critical for accurate radiative transfer in near-IR windows.

The correction is expressed through a dimensionless χ-factor (chi-factor):

    f_subL(Δν) = f_Voigt(Δν) × χ(|Δν|)

where Δν = ν − ν₀ is the detuning from the line centre.  By definition
χ(0) = 1; it falls smoothly below 1 for large |Δν|, suppressing the
unphysically strong Lorentzian far wings.

Implemented χ-factor models
----------------------------
``chi_tonkov96``
    Tonkov et al. (1996), derived from room-temperature CO₂ self-broadened
    measurements at 1.4 and 1.6 µm.  Piecewise-exponential in three regions.
    **Recommended for Venus 1.4 µm window modelling.**

``chi_perrin89``
    Perrin & Hartmann (1989), derived for CO₂–N₂ mixtures at 4.3 µm.
    Stronger far-wing suppression; use as a sensitivity check.

References
----------
Tonkov, M.V. et al. (1996). Measurements and empirical modeling of pure CO₂
    absorption in the 2.3- and 1.6-µm region at room temperature: far wings,
    allowed and collision-induced bands.
    J. Quant. Spectrosc. Radiat. Transfer **56**, 783–794.
    https://doi.org/10.1016/0022-4073(96)00082-3

Perrin, M.Y. & Hartmann, J.-M. (1989). Temperature-dependent measurements
    and modeling of absorption by CO₂–N₂ mixtures in the far line-wings of the
    4.3-µm CO₂ band.
    J. Quant. Spectrosc. Radiat. Transfer **42**, 311–318.
    https://doi.org/10.1016/0022-4073(89)90077-0

Meadows, V.S. & Crisp, D. (1996). Ground-based near-infrared observations of
    the Venus nightside: The thermal structure and water abundance near the
    surface.
    J. Geophys. Res. Planets **101**, 4595–4622.
    https://doi.org/10.1029/95JE03567
"""

from __future__ import annotations

import logging
from typing import Callable, Optional

import numpy as np
from scipy.special import wofz

log = logging.getLogger(__name__)

# ── Numerical constants ───────────────────────────────────────────────────────
_LN2     = np.log(2.0)
_SQRT2   = np.sqrt(2.0)
_SQRTPI  = np.sqrt(np.pi)
_SQRT2LN2 = np.sqrt(2.0 * _LN2)   # = sqrt(2 ln 2) ≈ 1.1774


# ═══════════════════════════════════════════════════════════════════════════════
# χ-factor implementations
# ═══════════════════════════════════════════════════════════════════════════════

def chi_tonkov96(delta_nu: np.ndarray) -> np.ndarray:
    r"""Tonkov et al. (1996) sub-Lorentzian χ-factor for CO₂.

    Piecewise-exponential parameterisation derived from room-temperature
    self-broadened CO₂ measurements at 1.4 and 1.6 µm (see Table 1 and
    Eq. 3 of the reference):

    .. math::

        \chi(\Delta\nu) = \begin{cases}
            1                                                     & |\Delta\nu| \le 3 \text{ cm}^{-1} \\
            \exp\!\bigl(-0.0214\,(|\Delta\nu| - 3)\bigr)         & 3 < |\Delta\nu| \le 30 \text{ cm}^{-1} \\
            \chi(30)\cdot\exp\!\bigl(-0.393\,(|\Delta\nu| - 30)\bigr) & |\Delta\nu| > 30 \text{ cm}^{-1}
        \end{cases}

    where :math:`\chi(30) = \exp(-0.0214 \times 27) \approx 0.561`.

    Parameters
    ----------
    delta_nu : array-like
        Absolute detuning from line centre :math:`|\nu - \nu_0|` [cm⁻¹].
        May be a scalar or any-shape ndarray.

    Returns
    -------
    chi : ndarray  (same shape as *delta_nu*)
        Dimensionless correction factor, :math:`0 < \chi \le 1`.

    Notes
    -----
    The parameters (0.0214 cm, 0.393 cm, breakpoint at 30 cm⁻¹) are specific
    to CO₂ self-broadening.  For Venus (VMR_CO2 ≈ 0.965) this is the dominant
    broadening mechanism and these values are directly applicable.
    """
    d = np.asarray(delta_nu, dtype=float)
    chi = np.ones_like(d)

    # Region 2: intermediate far wing (3 – 30 cm⁻¹)
    m2 = (d > 3.0) & (d <= 30.0)
    if m2.any():
        chi[m2] = np.exp(-0.0214 * (d[m2] - 3.0))

    # Region 3: deep far wing (> 30 cm⁻¹)
    m3 = d > 30.0
    if m3.any():
        chi_at_30 = np.exp(-0.0214 * 27.0)          # χ(30) ≈ 0.561
        chi[m3] = chi_at_30 * np.exp(-0.393 * (d[m3] - 30.0))

    return chi


def chi_perrin89(delta_nu: np.ndarray) -> np.ndarray:
    r"""Perrin & Hartmann (1989) sub-Lorentzian χ-factor for CO₂.

    Derived for CO₂–N₂ mixtures at 4.3 µm but widely used as a sensitivity
    benchmark.  Applies stronger suppression beyond 30 cm⁻¹ than Tonkov96:

    .. math::

        \chi(\Delta\nu) = \begin{cases}
            1                                                       & |\Delta\nu| \le 3 \text{ cm}^{-1} \\
            \exp\!\bigl(-0.0302\,(|\Delta\nu| - 3)\bigr)           & 3 < |\Delta\nu| \le 30 \text{ cm}^{-1} \\
            \chi(30)\cdot\exp\!\bigl(-0.600\,(|\Delta\nu| - 30)\bigr) & |\Delta\nu| > 30 \text{ cm}^{-1}
        \end{cases}

    Parameters
    ----------
    delta_nu : array-like
        :math:`|\nu - \nu_0|` [cm⁻¹].

    Returns
    -------
    chi : ndarray
        Correction factor, :math:`0 < \chi \le 1`.
    """
    d = np.asarray(delta_nu, dtype=float)
    chi = np.ones_like(d)

    a1, dnu1, a2 = 0.0302, 30.0, 0.600

    m2 = (d > 3.0) & (d <= dnu1)
    if m2.any():
        chi[m2] = np.exp(-a1 * (d[m2] - 3.0))

    m3 = d > dnu1
    if m3.any():
        chi_at_dnu1 = np.exp(-a1 * (dnu1 - 3.0))
        chi[m3] = chi_at_dnu1 * np.exp(-a2 * (d[m3] - dnu1))

    return chi


# ── Registry of available χ-factor models ────────────────────────────────────

#: Dictionary mapping model name → callable (or ``None`` for no correction).
#:
#: Use with :func:`~pymiedap.ckdistribution.hitran.compute_cross_section_subL`:
#:
#: .. code-block:: python
#:
#:     from pymiedap.ckdistribution.lineprofiles import CHI_FACTORS
#:     nu, sig = compute_cross_section_subL(..., chi_fn=CHI_FACTORS['tonkov96'])
CHI_FACTORS: dict[str, Optional[Callable]] = {
    'tonkov96': chi_tonkov96,
    'perrin89': chi_perrin89,
    'none':     None,          # standard Voigt, no correction
}


# ═══════════════════════════════════════════════════════════════════════════════
# Voigt profile utility
# ═══════════════════════════════════════════════════════════════════════════════

def voigt_profile(
    delta_nu: np.ndarray,
    gamma_D:  float,
    gamma_L:  float,
) -> np.ndarray:
    r"""Normalised Voigt line profile [cm] via the Faddeeva function.

    The Voigt profile is the convolution of a Gaussian (Doppler broadening)
    and a Lorentzian (pressure broadening):

    .. math::

        V(\Delta\nu;\,\gamma_D,\,\gamma_L)
        = \frac{\operatorname{Re}\bigl[w(z)\bigr]}{\sigma_D\,\sqrt{2\pi}}

    where

    .. math::

        z = \frac{\Delta\nu + i\,\gamma_L}{\sigma_D\,\sqrt{2}},
        \quad \sigma_D = \frac{\gamma_D}{\sqrt{2\ln 2}}

    and :math:`w(z) = e^{-z^2}\operatorname{erfc}(-iz)` is the Faddeeva
    (complex error) function, computed via :func:`scipy.special.wofz`.

    The profile is normalised so that :math:`\int_{-\infty}^{+\infty} V(\Delta\nu)\,
    d(\Delta\nu) = 1` when integrated in cm⁻¹ units.

    Parameters
    ----------
    delta_nu : ndarray [cm⁻¹]
        Signed detuning :math:`\nu - \nu_0` from the line centre.
    gamma_D : float [cm⁻¹]
        Doppler HWHM.
    gamma_L : float [cm⁻¹]
        Lorentzian (pressure-broadened) HWHM.

    Returns
    -------
    V : ndarray [cm]  (= 1/cm⁻¹)
        Normalised line profile. Multiply by line intensity S(T) [cm/molecule]
        to obtain the absorption cross-section σ(ν) [cm²/molecule].

    Notes
    -----
    For pure Doppler (γ_L → 0): V(Δν) → Gaussian with HWHM γ_D.
    For pure Lorentzian (γ_D → 0): V(Δν) → γ_L / [π(Δν² + γ_L²)].
    Both limits integrate to 1 over ℝ.

    Accuracy is limited by :func:`scipy.special.wofz`, which returns values
    correct to 14 significant digits.
    """
    # Guard against degenerate widths
    if gamma_D <= 0.0:
        gamma_D = 1e-30
    if gamma_L < 0.0:
        gamma_L = 0.0

    sigma_D = gamma_D / _SQRT2LN2          # Gaussian std dev [cm⁻¹]
    z = (np.asarray(delta_nu, dtype=complex) + 1j * gamma_L) / (sigma_D * _SQRT2)
    return np.real(wofz(z)) / (sigma_D * _SQRT2 * _SQRTPI)


def doppler_hwhm(nu0: np.ndarray, T: float, molar_mass_g: float) -> np.ndarray:
    r"""Doppler (thermal) HWHM for spectral lines [cm⁻¹].

    .. math::

        \gamma_D = \nu_0 \sqrt{\frac{8\ln 2\, R\, T}{M\, c^2}}

    Parameters
    ----------
    nu0 : ndarray [cm⁻¹]
        Line centre wavenumbers.
    T : float [K]
        Temperature.
    molar_mass_g : float [g/mol]
        Molar mass of the absorbing molecule.

    Returns
    -------
    gamma_D : ndarray [cm⁻¹]
    """
    # R [J/(mol K)] = 8.31446  ;  c [cm/s] = 2.99792458e10
    # γ_D/ν₀ = √(8 ln2 · R · T / (M_kg · c²))
    #         = √(8 ln2 · 8.31446 / (M_kg · (2.99792458e10)²)) · √T
    # Factor at 296 K for CO2 (M=44.01 g/mol): ~4.19e-5 × √(T)
    # Use c in m/s (consistent with R in J/(mol K) = kg m²/(s² mol K)).
    # HWHM formula: γ_D = ν₀ √(2 ln2 · R · T / (M_kg · c_m²))
    # Using 8 ln2 would give FWHM; mixing c in cm/s would give a ~50× error.
    M_kg = molar_mass_g * 1e-3
    R    = 8.31446        # J/(mol K)
    c_m  = 2.99792458e8   # m/s  ← must be SI to match R
    return np.asarray(nu0) * np.sqrt(2.0 * _LN2 * R * T / (M_kg * c_m**2))
