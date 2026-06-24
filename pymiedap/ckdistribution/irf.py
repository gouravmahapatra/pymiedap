"""
Instrument Response Functions (IRF) for spectral convolution.

Each function returns a *callable* that maps a wavenumber array [cm⁻¹] to a
non-negative weight array of the same length.  Pass the returned callable as
the *irf_fn* argument to :func:`.kdistribution.kspec_layer` or
:func:`.kdistribution.compute_bmabs`.

Available IRFs
--------------
- :func:`gaussian_irf` — symmetric Gaussian with a given FWHM.
- :func:`box_irf` — uniform (top-hat) over a given half-width.
- :func:`spicav_ir_irf` — SPICAV-IR instrument PSF loaded from a data file.
- :func:`no_irf` — returns flat weights (equivalent to using no IRF).
"""

from __future__ import annotations

import numpy as np
from typing import Callable, Optional
from scipy import interpolate


# ---------------------------------------------------------------------------
# Gaussian IRF
# ---------------------------------------------------------------------------

def gaussian_irf(
    center_wvn: float,
    fwhm_nm: float,
    wvl_um: float,
) -> Callable[[np.ndarray], np.ndarray]:
    """Return a Gaussian IRF centred at *center_wvn*.

    Parameters
    ----------
    center_wvn : float
        Centre wavenumber [cm⁻¹] of the Gaussian.
    fwhm_nm : float
        Full-width at half-maximum of the instrument slit [nm].
    wvl_um : float
        Reference wavelength used to convert nm → cm⁻¹ [µm].

    Returns
    -------
    Callable[[ndarray], ndarray]
        Function that takes a wavenumber array and returns Gaussian weights.
    """
    # Convert FWHM in nm to FWHM in cm⁻¹ at the reference wavelength
    wvl_cm = wvl_um * 1.0e-4  # µm → cm
    fwhm_cm = fwhm_nm * 1.0e-7  # nm → cm
    fwhm_wvn = fwhm_cm / (wvl_cm ** 2)  # Δν = Δλ / λ²

    sigma_wvn = fwhm_wvn / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    def _irf(specv: np.ndarray) -> np.ndarray:
        return np.exp(-0.5 * ((specv - center_wvn) / sigma_wvn) ** 2)

    return _irf


# ---------------------------------------------------------------------------
# Box (top-hat) IRF
# ---------------------------------------------------------------------------

def box_irf(
    center_wvn: float,
    half_width_wvn: float,
) -> Callable[[np.ndarray], np.ndarray]:
    """Return a flat top-hat IRF of total width 2 × *half_width_wvn*.

    Parameters
    ----------
    center_wvn : float
        Centre wavenumber [cm⁻¹].
    half_width_wvn : float
        Half-width of the box [cm⁻¹].

    Returns
    -------
    Callable[[ndarray], ndarray]
        Function that returns 1 inside the box and 0 outside.
    """
    def _irf(specv: np.ndarray) -> np.ndarray:
        weights = np.zeros_like(specv, dtype=float)
        weights[np.abs(specv - center_wvn) <= half_width_wvn] = 1.0
        return weights

    return _irf


# ---------------------------------------------------------------------------
# No IRF (flat weights)
# ---------------------------------------------------------------------------

def no_irf() -> Callable[[np.ndarray], np.ndarray]:
    """Return a flat-weight function (equivalent to no IRF).

    Returns
    -------
    Callable[[ndarray], ndarray]
        Function that returns an array of ones.
    """
    def _irf(specv: np.ndarray) -> np.ndarray:
        return np.ones_like(specv, dtype=float)

    return _irf


# ---------------------------------------------------------------------------
# SPICAV-IR IRF loaded from a PSF file
# ---------------------------------------------------------------------------

def spicav_ir_irf(
    psf_file: str,
    perch: float = 0.0,
) -> Callable[[np.ndarray], np.ndarray]:
    """Load the SPICAV-IR instrument PSF from a text file.

    The PSF file must contain three whitespace-separated columns::

        dnu[cm⁻¹]   freq[kHz]   normalized_response

    This matches the format of SPICAV-IR calibration files produced by the
    LATMOS / IKI team (e.g. ``psf_lw_all_O2_1270_desc.txt``).

    Parameters
    ----------
    psf_file : str
        Path to the PSF data file.
    perch : float, optional
        Percentage perturbation applied to the central peak of the PSF
        (used for sensitivity studies).  Default 0 (no perturbation).

    Returns
    -------
    Callable[[ndarray], ndarray]
        Function that takes a relative-wavenumber array ``Δν = ν - ν₀`` [cm⁻¹]
        and returns the (possibly perturbed) normalised SPICAV-IR response.

    Notes
    -----
    The returned callable takes an absolute wavenumber array *specv* as input.
    You must supply the central wavenumber *center_wvn* at call time by using
    :func:`spicav_ir_irf_at_center` which wraps this function.
    """
    import numpy as _np
    from scipy import interpolate as _interp

    irf_data = _np.loadtxt(psf_file, skiprows=1)
    dnu  = irf_data[:, 0]   # cm⁻¹ relative wavenumber
    ch0  = irf_data[:, 2]   # normalised response

    # Apply optional perturbation to the central peak region
    if perch != 0.0:
        ch0_perturbed = ch0.copy()
        central_mask = (dnu > -4.5) & (dnu < 6)
        ch0_perturbed[central_mask] = 0.0
        ch0 = ch0 + perch * ch0_perturbed / 100.0

    # Extend the domain to ±100 cm⁻¹ to avoid edge artefacts
    dnu_ext  = _np.arange(-100, 100, abs(dnu[0] - dnu[1]))
    irf_base = _interp.Akima1DInterpolator(dnu, ch0)
    ch0_ext  = irf_base(dnu_ext)
    ch0_ext  = _np.nan_to_num(ch0_ext)

    irf_interp = _interp.Akima1DInterpolator(dnu_ext, ch0_ext)

    def _make_irf(center_wvn: float) -> Callable[[_np.ndarray], _np.ndarray]:
        """Create an IRF callable centred at *center_wvn* [cm⁻¹]."""
        def _irf(specv: _np.ndarray) -> _np.ndarray:
            dnu_arr = specv - center_wvn
            response = irf_interp(dnu_arr)
            response = _np.nan_to_num(response)
            response = _np.clip(response, 0.0, None)
            return response

        return _irf

    return _make_irf


def spicav_ir_irf_at_center(
    psf_file: str,
    center_wvn: float,
    perch: float = 0.0,
) -> Callable[[np.ndarray], np.ndarray]:
    """Convenience wrapper: return a SPICAV-IR IRF fixed at *center_wvn*.

    Parameters
    ----------
    psf_file : str
        Path to the SPICAV-IR PSF data file.
    center_wvn : float
        Wavenumber of the spectral channel centre [cm⁻¹].
    perch : float, optional
        Perturbation percentage for sensitivity studies.

    Returns
    -------
    Callable[[ndarray], ndarray]
        IRF function taking an absolute wavenumber array [cm⁻¹].
    """
    make_irf = spicav_ir_irf(psf_file, perch=perch)
    return make_irf(center_wvn)


# ---------------------------------------------------------------------------
# Convolution utility (for post-processing spectra)
# ---------------------------------------------------------------------------

def convolve_spectrum(
    nu: np.ndarray,
    flux: np.ndarray,
    irf_fn: Callable[[np.ndarray], np.ndarray],
    center_wvn: Optional[float] = None,
    output_nu: Optional[np.ndarray] = None,
) -> tuple:
    """Convolve a spectrum with an IRF in wavenumber space.

    Parameters
    ----------
    nu : ndarray
        Input wavenumber grid [cm⁻¹].
    flux : ndarray
        Input spectrum (any units) on *nu*.
    irf_fn : callable
        IRF function ``irf_fn(nu_array) → weights``.  If *center_wvn* is
        provided, the function is evaluated at ``nu - center_wvn``; otherwise
        the IRF is evaluated at *nu* directly (absolute wavenumbers).
    center_wvn : float or None
        Centre wavenumber for a relative IRF [cm⁻¹].
    output_nu : ndarray or None
        Wavenumber grid for the output.  If None, the output is on *nu*.

    Returns
    -------
    out_nu : ndarray
        Output wavenumber grid.
    out_flux : ndarray
        Convolved spectrum.
    """
    weights = irf_fn(nu)
    if np.sum(weights) == 0:
        weights = np.ones_like(nu)

    # Convolve via weighted integration
    if output_nu is None:
        output_nu = nu.copy()

    out_flux = np.zeros_like(output_nu, dtype=float)
    for i, nu_out in enumerate(output_nu):
        dnu = np.abs(nu - nu_out)
        w = irf_fn(nu - nu_out if center_wvn is None else nu - nu_out)
        denom = np.trapz(w, nu)
        if denom > 0:
            out_flux[i] = np.trapz(w * flux, nu) / denom
        else:
            # Nearest-neighbour fallback
            out_flux[i] = flux[np.argmin(dnu)]

    return output_nu, out_flux
