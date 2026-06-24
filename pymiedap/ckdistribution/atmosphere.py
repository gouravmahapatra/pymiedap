"""
Atmosphere class and planet preset profiles.

The ``Atmosphere`` object holds the 1-D vertical profile (altitude, pressure,
temperature) together with volume-mixing-ratio (VMR) profiles for one or more
absorbing gases.  All atmospheric computations in ckdistribution work against
this abstraction so no planet-specific constants are hardcoded elsewhere.
"""

from __future__ import annotations

import numpy as np
from typing import Dict, Optional
from scipy import interpolate

from .constants import G_EARTH, G_VENUS, G_MARS, N_A


class Atmosphere:
    """Vertical atmospheric profile.

    Parameters
    ----------
    altitude : array_like
        Layer boundary altitudes [km], monotonically increasing, length N.
    pressure : array_like
        Pressure at each altitude level [bar], length N.
    temperature : array_like
        Temperature at each altitude level [K], length N.
    gas_vmr : dict, optional
        Mapping ``{molecule_name: vmr_array}`` where *vmr_array* has length N.
        Volume mixing ratios are dimensionless fractions (0 → 1).
    gravity : float, optional
        Gravitational acceleration [m s⁻²].  Defaults to Earth standard.
    name : str, optional
        Human-readable label for this profile.

    Attributes
    ----------
    altitude, pressure, temperature : ndarray
        As passed, converted to numpy arrays.
    gas_vmr : dict
        VMR profiles (same length as altitude).
    gravity : float
    name : str
    nlev : int
        Number of altitude levels.
    """

    def __init__(
        self,
        altitude: np.ndarray,
        pressure: np.ndarray,
        temperature: np.ndarray,
        gas_vmr: Optional[Dict[str, np.ndarray]] = None,
        gravity: float = G_EARTH,
        name: str = "atmosphere",
    ):
        self.altitude = np.asarray(altitude, dtype=float)
        self.pressure = np.asarray(pressure, dtype=float)
        self.temperature = np.asarray(temperature, dtype=float)
        self.gas_vmr = {k: np.asarray(v, dtype=float)
                        for k, v in (gas_vmr or {}).items()}
        self.gravity = float(gravity)
        self.name = name

        if not (self.altitude.shape == self.pressure.shape == self.temperature.shape):
            raise ValueError("altitude, pressure, and temperature must have the same length.")

    @property
    def nlev(self) -> int:
        """Number of vertical levels."""
        return len(self.altitude)

    # ------------------------------------------------------------------
    # Column-density helper
    # ------------------------------------------------------------------

    def column_density(self, molecule: str, layer_idx: int,
                       mol_mass_g_per_mol: float) -> float:
        """Number column density [molecules m⁻²] for a given layer.

        Computed from the pressure difference between layer boundaries using
        the hydrostatic equation:

            N = (Δp / (m_molecule * g)) × VMR

        Parameters
        ----------
        molecule : str
            Gas name as in ``gas_vmr``.
        layer_idx : int
            Index *i* for the layer between levels *i* and *i+1* (0-based).
        mol_mass_g_per_mol : float
            Molar mass of the molecule [g/mol].

        Returns
        -------
        float
            Column density in molecules per m².
        """
        p_bot = self.pressure[layer_idx] * 1.0e5       # bar → Pa
        p_top = self.pressure[layer_idx + 1] * 1.0e5   # bar → Pa
        if p_bot < p_top:
            raise ValueError(
                f"Bottom pressure ({p_bot:.3f} Pa) < top pressure ({p_top:.3f} Pa) "
                f"at layer {layer_idx}.  Pressure array must be decreasing with altitude."
            )
        m_mol = mol_mass_g_per_mol * 1.0e-3  # g/mol → kg/mol
        m_molecule = m_mol / N_A              # kg/molecule
        vmr = self.gas_vmr.get(molecule, np.ones(self.nlev))
        vmr_avg = 0.5 * (vmr[layer_idx] + vmr[layer_idx + 1])
        Nd = (p_bot - p_top) / (m_molecule * self.gravity)  # molecules/m²
        return Nd * vmr_avg

    # ------------------------------------------------------------------
    # Layer geometry helpers
    # ------------------------------------------------------------------

    def layer_bounds(self, layer_idx: int) -> dict:
        """Return the physical bounds and average state for a given layer.

        A layer is the atmospheric slab between altitude levels *layer_idx*
        (bottom) and *layer_idx + 1* (top).

        Parameters
        ----------
        layer_idx : int
            Zero-based layer index.  Must satisfy ``0 <= layer_idx < nlev - 1``.

        Returns
        -------
        dict with keys:
            * ``alt_bot_km``  – altitude of layer bottom [km]
            * ``alt_top_km``  – altitude of layer top [km]
            * ``T_avg_K``     – layer-mean temperature [K]
            * ``P_avg_bar``   – layer-mean pressure [bar]
            * ``P_bot_bar``   – pressure at layer bottom [bar]
            * ``P_top_bar``   – pressure at layer top [bar]

        Raises
        ------
        IndexError
            If *layer_idx* is out of range.
        """
        nlayer = self.nlev - 1
        if not (0 <= layer_idx < nlayer):
            raise IndexError(
                f"layer_idx {layer_idx} is out of range for an atmosphere with "
                f"{self.nlev} levels ({nlayer} layers)."
            )
        i = layer_idx
        return {
            'alt_bot_km': float(self.altitude[i]),
            'alt_top_km': float(self.altitude[i + 1]),
            'T_avg_K':    float(0.5 * (self.temperature[i] + self.temperature[i + 1])),
            'P_avg_bar':  float(0.5 * (self.pressure[i] + self.pressure[i + 1])),
            'P_bot_bar':  float(self.pressure[i]),
            'P_top_bar':  float(self.pressure[i + 1]),
        }

    @property
    def nlayer(self) -> int:
        """Number of atmospheric layers (= nlev - 1)."""
        return max(0, self.nlev - 1)

    # ------------------------------------------------------------------
    # Interpolation
    # ------------------------------------------------------------------

    def interpolate_to(self, new_altitude: np.ndarray) -> "Atmosphere":
        """Return a new ``Atmosphere`` interpolated onto *new_altitude* [km]."""
        alt = self.altitude
        T_fn = interpolate.interp1d(alt, self.temperature, kind='linear',
                                    fill_value='extrapolate')
        P_fn = interpolate.interp1d(alt, np.log(self.pressure), kind='linear',
                                    fill_value='extrapolate')
        new_T = T_fn(new_altitude)
        new_P = np.exp(P_fn(new_altitude))
        new_vmr = {}
        for mol, vmr in self.gas_vmr.items():
            fn = interpolate.interp1d(alt, vmr, kind='linear',
                                      fill_value='extrapolate')
            new_vmr[mol] = np.clip(fn(new_altitude), 0.0, 1.0)
        return Atmosphere(new_altitude, new_P, new_T, gas_vmr=new_vmr,
                          gravity=self.gravity, name=self.name)

    def __repr__(self) -> str:
        gases = list(self.gas_vmr.keys()) or ['none']
        return (f"Atmosphere(name={self.name!r}, nlev={self.nlev}, "
                f"p=[{self.pressure.min():.3g}..{self.pressure.max():.3g}] bar, "
                f"T=[{self.temperature.min():.1f}..{self.temperature.max():.1f}] K, "
                f"gases={gases})")


# ---------------------------------------------------------------------------
# Planet presets
# ---------------------------------------------------------------------------

def earth_standard(nlev: int = 50) -> Atmosphere:
    """U.S. Standard Atmosphere 1976 — simplified 7-region model.

    Parameters
    ----------
    nlev : int
        Number of altitude levels to sample (linearly spaced 0–86 km).

    Returns
    -------
    Atmosphere
        Earth standard atmosphere with CO2 (400 ppm), CH4 (1.8 ppm),
        H2O (column-mean ~0.5%), N2O (320 ppb), CO (0.1 ppm).
    """
    # Anchor points: (altitude km, T K, P bar)
    _anchors = np.array([
        [  0.0,  288.15, 1.01325],
        [ 11.0,  216.65, 2.26320e-1],
        [ 20.0,  216.65, 5.47489e-2],
        [ 32.0,  228.65, 8.68019e-3],
        [ 47.0,  270.65, 1.10906e-3],
        [ 51.0,  270.65, 6.69388e-4],
        [ 71.0,  214.65, 3.95642e-5],
        [ 86.0,  186.87, 3.73383e-6],
    ])
    alt_anchor = _anchors[:, 0]
    T_anchor   = _anchors[:, 1]
    P_anchor   = _anchors[:, 2]

    alt = np.linspace(0.0, 86.0, nlev)
    T_fn = interpolate.interp1d(alt_anchor, T_anchor, kind='linear')
    P_fn = interpolate.interp1d(alt_anchor, np.log(P_anchor), kind='linear')
    T = T_fn(alt)
    P = np.exp(P_fn(alt))

    # Simple VMR profiles (constant with altitude — coarse but usable)
    vmr = {
        'CO2': np.full(nlev, 400.0e-6),    # 400 ppm
        'CH4': np.full(nlev, 1.8e-6),      # 1.8 ppm
        'H2O': np.full(nlev, 5.0e-3),      # ~0.5% (troposphere mean)
        'N2O': np.full(nlev, 320.0e-9),    # 320 ppb
        'CO':  np.full(nlev, 0.1e-6),      # 0.1 ppm
        'O3':  np.full(nlev, 50.0e-9),     # 50 ppb (stratospheric peak ignored)
        'O2':  np.full(nlev, 0.2095),      # 20.95%
    }

    return Atmosphere(alt, P, T, gas_vmr=vmr, gravity=G_EARTH,
                      name='Earth Standard 1976')


def venus_ignatiev(nlev: int = 71) -> Atmosphere:
    """Venus atmosphere based on the Ignatiev et al. 38-layer profile.

    The profile is parameterised analytically from the tabulated pressures
    and temperatures.  CO2 is the dominant gas (~96.5%), with SO2, N2, CO
    as minor species.  This matches the profile used in the original
    ck-distribution code (``VenusIg.dat`` / ``vertical_profile_ignatiev_38lays.htp``).

    Parameters
    ----------
    nlev : int
        Number of output levels (default 71 as in the original 71-layer model).

    Returns
    -------
    Atmosphere
        Venus atmosphere with CO2 VMR = 0.965.
    """
    # Ignatiev et al. anchor points (altitude km, T K, P bar)
    # Derived from Table 1 of Ignatiev et al. (2009, J. Geophys. Res.)
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
    alt_anchor = _ig[:, 0]
    T_anchor   = _ig[:, 1]
    P_anchor   = _ig[:, 2]

    # Low: 0–40 km at 4 km resolution; high: 41–100 km at 1 km resolution
    alt_low  = np.array([0., 4., 8., 12., 16., 20., 24., 28., 32., 36., 40.])
    alt_high = np.arange(41, 101, 1, dtype=float)
    alt = np.concatenate([alt_low, alt_high])
    alt = alt[:nlev]

    T_fn = interpolate.interp1d(alt_anchor, T_anchor, kind='linear',
                                fill_value='extrapolate')
    P_fn = interpolate.interp1d(alt_anchor, np.log(P_anchor), kind='linear',
                                fill_value='extrapolate')
    T = T_fn(alt)
    P = np.exp(P_fn(alt))

    vmr = {
        'CO2': np.full(len(alt), 0.965),   # 96.5 %
        'SO2': np.full(len(alt), 1.0e-4),  # ~100 ppm
        'N2':  np.full(len(alt), 0.035),   # ~3.5 %
        'CO':  np.full(len(alt), 2.0e-5),  # ~20 ppm
    }

    return Atmosphere(alt, P, T, gas_vmr=vmr, gravity=G_VENUS,
                      name='Venus Ignatiev')


def mars_mgs(nlev: int = 36) -> Atmosphere:
    """Simple Mars profile based on MGS radio-occultation climatology.

    CO2 dominates (~95.3 %).  The profile covers 0–120 km.

    Parameters
    ----------
    nlev : int
        Number of output levels.

    Returns
    -------
    Atmosphere
        Mars atmosphere with CO2 VMR = 0.953.
    """
    _mars = np.array([
        [   0.,  210., 6.36e-3],
        [  10.,  195., 2.84e-3],
        [  20.,  185., 1.24e-3],
        [  30.,  175., 5.21e-4],
        [  40.,  165., 2.12e-4],
        [  50.,  155., 8.32e-5],
        [  60.,  148., 3.12e-5],
        [  70.,  143., 1.13e-5],
        [  80.,  140., 3.94e-6],
        [  90.,  139., 1.32e-6],
        [ 100.,  143., 4.28e-7],
        [ 110.,  151., 1.35e-7],
        [ 120.,  165., 4.14e-8],
    ])
    alt_anchor = _mars[:, 0]
    T_anchor   = _mars[:, 1]
    P_anchor   = _mars[:, 2]

    alt = np.linspace(0., 120., nlev)
    T_fn = interpolate.interp1d(alt_anchor, T_anchor, kind='linear',
                                fill_value='extrapolate')
    P_fn = interpolate.interp1d(alt_anchor, np.log(P_anchor), kind='linear',
                                fill_value='extrapolate')
    T = T_fn(alt)
    P = np.exp(P_fn(alt))

    vmr = {
        'CO2': np.full(nlev, 0.953),
        'N2':  np.full(nlev, 0.027),
        'Ar':  np.full(nlev, 0.016),
        'CO':  np.full(nlev, 7.0e-4),
    }

    return Atmosphere(alt, P, T, gas_vmr=vmr, gravity=G_MARS,
                      name='Mars MGS')
