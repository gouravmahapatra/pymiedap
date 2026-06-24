"""
PyMieDAP integration: reflected-light spectrum with gaseous absorption.

This module ties the CK-distribution machinery to PyMieDAP's radiative
transfer engine.  The key idea is:

1. Compute ``bmabs[layer, wavelength, gauss_point]`` via the CKD pipeline.
2. For each wavelength *w* and Gauss point *g*, set ``layer.tau_g = bmabs[:, w, g]``
   on the PyMieDAP model and call ``model.dap_code()`` to compute the
   Stokes vector.
3. Accumulate the results weighted by the Gauss–Legendre weights.

The Gauss–Legendre integration approximates the band-averaged Stokes vector:

    I_band(λ) ≈ (1/2) Σ_g  w_g  I(k_g(λ))

where the factor of 1/2 comes from rescaling [−1, 1] → [0, 1].

Usage example::

    import numpy as np
    import pymiedap.pymiedap as pmd
    from pymiedap.ckdistribution import (
        earth_standard, compute_bmabs, gauss_legendre_points, compute_reflected_spectrum
    )
    from pathlib import Path

    atm = earth_standard(nlev=20)
    wav = np.arange(1.56, 1.65, 0.005)   # µm, CO2 / CH4 band
    gp, gw = gauss_legendre_points(10)

    I_band, Q_band, U_band, V_band = compute_reflected_spectrum(
        model_factory=my_model_factory,
        bmabs=my_bmabs,
        wav=wav,
        gauss_points=gp,
        gauss_weights=gw,
        phase=np.array([0.]),
        sza=np.array([30.]),
        emission=np.array([0.]),
        output_dir=Path('/tmp/dap_output'),
    )
"""

from __future__ import annotations

import datetime
import logging
import os
from pathlib import Path
from typing import Callable, Dict, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Core integration loop
# ---------------------------------------------------------------------------

def compute_reflected_spectrum(
    model_factory: Callable[[], object],
    bmabs: np.ndarray,
    wav: np.ndarray,
    gauss_points: np.ndarray,
    gauss_weights: np.ndarray,
    phase: np.ndarray,
    sza: np.ndarray,
    emission: np.ndarray,
    output_dir: Path,
    *,
    azimuth: Optional[np.ndarray] = None,
    beta: Optional[np.ndarray] = None,
    rename: bool = False,
    verbose: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute the CK-averaged reflected Stokes spectrum.

    Parameters
    ----------
    model_factory : callable
        Zero-argument callable that returns a freshly configured
        :class:`pymiedap.pymiedap.Model` (or sub-class) instance.  The model
        must have its cloud / aerosol / surface layers already set up; only
        ``tau_g`` (gaseous absorption) will be modified per Gauss point.
    bmabs : ndarray, shape (nlayer, nwav, ngauss)
        Band-mean absorption optical depth array from
        :func:`.kdistribution.compute_bmabs`.  Index order must match:
        layer index (surface → TOA), wavelength index, Gauss-point index.
    wav : ndarray, shape (nwav,)
        Central wavelengths [µm].
    gauss_points : ndarray, shape (ngauss,)
        Gauss–Legendre nodes on [0, 1].
    gauss_weights : ndarray, shape (ngauss,)
        Gauss–Legendre weights (sum = 2 on [−1, 1], so multiply by 0.5).
    phase, sza, emission : ndarray
        Geometry arrays [degrees], same length (ngeos).
    output_dir : Path
        Directory for temporary DAP Fourier output files.
    azimuth, beta : ndarray or None
        Additional geometry arrays.  If None, zeros are used.
    rename : bool
        Passed through to ``model.dap_code()``; controls output filename scheme.
    verbose : bool
        Print progress.

    Returns
    -------
    I_band, Q_band, U_band, V_band : ndarray, shape (nwav, ngeos)
        Gauss–Legendre-averaged Stokes parameters for each wavelength and
        geometry.

    Notes
    -----
    The factor of 0.5 in the weights arises from the change of variables
    from the standard Gauss–Legendre interval [−1, 1] to the k-distribution
    domain [0, 1]:

        ∫₀¹ f(g) dg ≈ (1/2) Σ_j w_j f(g_j)
    """
    import pymiedap.pymiedap as pmd

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    nwav   = len(wav)
    ngauss = len(gauss_points)
    ngeos  = len(phase)

    if azimuth is None:
        azimuth = np.zeros(ngeos)
    if beta is None:
        beta = np.zeros(ngeos)

    I_band = np.zeros((nwav, ngeos))
    Q_band = np.zeros((nwav, ngeos))
    U_band = np.zeros((nwav, ngeos))
    V_band = np.zeros((nwav, ngeos))

    for iw, wavel in enumerate(wav):
        if verbose:
            print(f"  Wavelength {iw+1:3d}/{nwav}  λ = {wavel:.5f} µm")

        for ig, (gp, gw) in enumerate(zip(gauss_points, gauss_weights)):
            # Build a fresh model for this Gauss point
            model = model_factory()

            # Inject gaseous absorption optical depth into each layer.
            # tau_g must remain a 1-element ndarray so that dap_code can index
            # it as layer.tau_g[z] with z=0 (single-wavelength call).
            for ilayer, layer in enumerate(model.layers):
                if ilayer < bmabs.shape[0]:
                    layer.tau_g = np.array([float(bmabs[ilayer, iw, ig])])
                else:
                    layer.tau_g = np.array([0.0])

            # Run DAP
            try:
                model.dap_code(
                    wav=[wavel],
                    output_dir=output_dir,
                    rename=rename,
                )
            except Exception as exc:
                log.warning(
                    "DAP failed at λ=%.5f µm, Gauss point %d/%d: %s. "
                    "Contribution set to zero.",
                    wavel, ig + 1, ngauss, exc
                )
                continue

            # Read DAP output
            try:
                fname = model.name[0]
                I0, Q0, U0, V0 = pmd.read_dap_output(
                    phase, sza, emission, fname,
                    beta=beta, phi=azimuth
                )
            except Exception as exc:
                log.warning(
                    "read_dap_output failed at λ=%.5f µm, Gauss %d/%d: %s.",
                    wavel, ig + 1, ngauss, exc
                )
                continue

            # Accumulate with Gauss–Legendre weight (×0.5 for [0,1] rescaling)
            weight = 0.5 * gw
            I_band[iw] += weight * I0
            Q_band[iw] += weight * Q0
            U_band[iw] += weight * U0
            V_band[iw] += weight * V0

    return I_band, Q_band, U_band, V_band


# ---------------------------------------------------------------------------
# High-level convenience wrapper
# ---------------------------------------------------------------------------

def run_ckd_spectrum(
    molecule: str,
    atmosphere,
    wav: np.ndarray,
    model_factory: Callable,
    phase: np.ndarray,
    sza: np.ndarray,
    emission: np.ndarray,
    output_dir: Path,
    cache_dir: Optional[str] = None,
    *,
    n_gauss: int = 10,
    sigma_um: float = 0.005,
    truncw_um: float = 0.2,
    nw: int = 50000,
    irf_fn=None,
    use_irf: bool = False,
    azimuth: Optional[np.ndarray] = None,
    beta: Optional[np.ndarray] = None,
    verbose: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """End-to-end CKD reflected-spectrum computation.

    Combines :func:`.kdistribution.compute_bmabs` and
    :func:`compute_reflected_spectrum` into a single call.

    Parameters
    ----------
    molecule : str
        Absorbing gas (e.g. ``'CO2'``).
    atmosphere : :class:`.atmosphere.Atmosphere`
        Atmospheric profile.
    wav : ndarray, shape (W,)
        Central wavelengths [µm].
    model_factory : callable
        Returns a configured :class:`pymiedap.pymiedap.Model`.
    phase, sza, emission : ndarray
        Observation geometry [degrees].
    output_dir : Path
        Directory for DAP output files.
    cache_dir : str or None
        HAPI line-list cache directory.  Defaults to ``~/.pymiedap/hitran``.
    n_gauss : int
        Number of Gauss–Legendre quadrature points.
    sigma_um : float
        Spectral window full-width [µm].
    truncw_um : float
        Truncation margin [µm].
    nw : int
        Wavenumber resolution [samples per µm].
    irf_fn : callable or None
        Instrument response function.
    use_irf : bool
        Use IRF-weighted k-scale.
    azimuth, beta : ndarray or None
        Additional geometry.
    verbose : bool
        Print progress.

    Returns
    -------
    wav : ndarray, shape (W,)
        Wavelengths [µm].
    I_band, Q_band, U_band, V_band : ndarray, shape (W, ngeos)
        Band-averaged Stokes parameters.
    """
    from .kdistribution import compute_bmabs, gauss_legendre_points

    if cache_dir is None:
        cache_dir = os.path.join(os.path.expanduser('~'), '.pymiedap', 'hitran')

    if verbose:
        print(f"[ckdistribution] Computing CKD spectrum for {molecule}")
        print(f"  Atmosphere:  {atmosphere}")
        print(f"  Wavelengths: {wav[0]:.4f} – {wav[-1]:.4f} µm  ({len(wav)} points)")
        print(f"  Gauss pts:   {n_gauss}")
        print(f"  Cache dir:   {cache_dir}")

    gp, gw = gauss_legendre_points(n_gauss)

    if verbose:
        print("\n[1/2] Computing bmabs ...")
    bmabs = compute_bmabs(
        atmosphere, molecule, wav, gp, cache_dir,
        sigma_um=sigma_um, truncw_um=truncw_um, nw=nw,
        irf_fn=irf_fn, use_irf=use_irf, verbose=verbose
    )

    if verbose:
        print("\n[2/2] Running PyMieDAP radiative transfer ...")
    I, Q, U, V = compute_reflected_spectrum(
        model_factory, bmabs, wav, gp, gw,
        phase, sza, emission, output_dir,
        azimuth=azimuth, beta=beta, verbose=verbose
    )

    return wav, I, Q, U, V


# ---------------------------------------------------------------------------
# bmabs I/O utilities  (HDF5-backed, self-describing)
# ---------------------------------------------------------------------------

def save_bmabs(
    bmabs: np.ndarray,
    wav: np.ndarray,
    gauss_weights: np.ndarray,
    filepath: str,
    *,
    molecule: str = "",
    atmosphere_name: str = "",
    gauss_order: int = 0,
    extra_meta: Optional[Dict[str, str]] = None,
) -> None:
    """Save a bmabs array and provenance metadata to an HDF5 file.

    The file layout is::

        /metadata/
            molecule          (str attr)
            atmosphere_name   (str attr)
            gauss_order       (int attr)
            wav_min_um        (float attr)
            wav_max_um        (float attr)
            n_wav             (int attr)
            n_layer           (int attr)
            n_gauss           (int attr)
            created           (ISO-8601 timestamp attr)
            <extra_meta keys> (str attrs)
        /coords/
            wav_um            (nwav,)   float64
            gauss_weights     (ngauss,) float64
        /bmabs                (nlayer, nwav, ngauss) float32, gzip compressed

    Parameters
    ----------
    bmabs : ndarray, shape (nlayer, nwav, ngauss)
    wav : ndarray, shape (nwav,)  [µm]
    gauss_weights : ndarray, shape (ngauss,)
    filepath : str
        Output path.  A ``.h5`` suffix is appended when absent.
    molecule : str
        Absorbing gas name (e.g. ``'CO2'``), stored in metadata.
    atmosphere_name : str
        Human-readable label from ``Atmosphere.name``, stored in metadata.
    gauss_order : int
        Number of Gauss–Legendre quadrature points, stored in metadata.
    extra_meta : dict or None
        Any additional key/value strings to store under ``/metadata/``.
    """
    try:
        import h5py
    except ImportError as exc:
        raise ImportError(
            "h5py is required for HDF5 bmabs storage.  "
            "Install it with:  pip install h5py"
        ) from exc

    path = Path(filepath)
    if path.suffix not in {'.h5', '.hdf5'}:
        path = path.with_suffix(path.suffix + '.h5')

    nlayer, nwav, ngauss = bmabs.shape
    with h5py.File(path, 'w') as f:
        # ── Metadata group ────────────────────────────────────────────────
        meta = f.create_group('metadata')
        meta.attrs['molecule']        = molecule
        meta.attrs['atmosphere_name'] = atmosphere_name
        meta.attrs['gauss_order']     = int(gauss_order) or ngauss
        meta.attrs['wav_min_um']      = float(wav.min())
        meta.attrs['wav_max_um']      = float(wav.max())
        meta.attrs['n_wav']           = nwav
        meta.attrs['n_layer']         = nlayer
        meta.attrs['n_gauss']         = ngauss
        meta.attrs['created']         = datetime.datetime.utcnow().isoformat() + 'Z'
        for k, v in (extra_meta or {}).items():
            meta.attrs[str(k)] = str(v)

        # ── Coordinate datasets ───────────────────────────────────────────
        coords = f.create_group('coords')
        coords.create_dataset('wav_um',       data=wav.astype('float64'))
        coords.create_dataset('gauss_weights', data=gauss_weights.astype('float64'))

        # ── Main bmabs dataset (float32 + gzip) ───────────────────────────
        f.create_dataset(
            'bmabs',
            data=bmabs.astype('float32'),
            compression='gzip',
            compression_opts=4,
            chunks=(min(nlayer, 4), min(nwav, 64), ngauss),
        )

    log.info("Saved bmabs → %s  (shape %s, molecule=%r, atm=%r)",
             path, bmabs.shape, molecule, atmosphere_name)


def load_bmabs(
    filepath: str,
    *,
    expected_wav: Optional[np.ndarray] = None,
    rtol: float = 1e-5,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load a bmabs array from an HDF5 file written by :func:`save_bmabs`.

    Parameters
    ----------
    filepath : str
        Path to ``.h5`` file (or without extension — ``.h5`` is appended).
    expected_wav : ndarray or None
        If given, the stored wavelength grid is compared against this array.
        A :class:`ValueError` is raised when they differ beyond *rtol*,
        preventing silent cache mismatches.
    rtol : float
        Relative tolerance for wavelength comparison (default 1 × 10⁻⁵).

    Returns
    -------
    bmabs : ndarray, shape (nlayer, nwav, ngauss)  float64
    wav : ndarray, shape (nwav,)  [µm]
    gauss_weights : ndarray, shape (ngauss,)
    """
    try:
        import h5py
    except ImportError as exc:
        raise ImportError(
            "h5py is required for HDF5 bmabs storage.  "
            "Install it with:  pip install h5py"
        ) from exc

    path = Path(filepath)
    if not path.exists() and path.suffix not in {'.h5', '.hdf5'}:
        path = path.with_suffix(path.suffix + '.h5')

    with h5py.File(path, 'r') as f:
        bmabs        = f['bmabs'][...].astype('float64')
        wav          = f['coords/wav_um'][...].astype('float64')
        gauss_weights = f['coords/gauss_weights'][...].astype('float64')

        # Log provenance for diagnostics
        meta = f['metadata']
        log.info(
            "Loaded bmabs from %s  (shape %s, molecule=%r, atm=%r, created=%r)",
            path, bmabs.shape,
            meta.attrs.get('molecule', ''),
            meta.attrs.get('atmosphere_name', ''),
            meta.attrs.get('created', ''),
        )

    # ── Consistency check ─────────────────────────────────────────────────
    if expected_wav is not None:
        expected_wav = np.asarray(expected_wav, dtype=float)
        if wav.shape != expected_wav.shape:
            raise ValueError(
                f"Cached bmabs has {wav.shape[0]} wavelengths but expected "
                f"{expected_wav.shape[0]}.  Delete {path} and recompute."
            )
        max_dev = np.max(np.abs(wav - expected_wav) / (np.abs(expected_wav) + 1e-30))
        if max_dev > rtol:
            raise ValueError(
                f"Cached wavelength grid differs from expected by {max_dev:.2e} "
                f"(rtol={rtol}).  Delete {path} and recompute."
            )

    return bmabs, wav, gauss_weights


def write_bmabs_text(
    bmabs: np.ndarray,
    wav: np.ndarray,
    gauss_weights: np.ndarray,
    filepath: str,
) -> None:
    """Write bmabs to a human-readable text file (original ckdis.py format).

    .. deprecated::
        The text format carries no provenance metadata and produces very large
        files.  Prefer :func:`save_bmabs` (HDF5) for new work.

    Parameters
    ----------
    bmabs : ndarray, shape (nlayer, nwav, ngauss)
    wav : ndarray, shape (nwav,)  [µm]
    gauss_weights : ndarray, shape (ngauss,)
    filepath : str
    """
    nlay, nwav, ngauss = bmabs.shape
    with open(filepath, 'w') as f:
        f.write(f'Number of wavelength values:\n{nwav}\n')
        f.write(f'Number of layers:\n{nlay}\n')
        f.write(f'Number of gauss points:\n{ngauss}\n')
        f.write('lamb_no    gausspoint    layer    bmabs_kd   lambda (microns)\n')
        for iw in range(nwav):
            for ig in range(ngauss):
                for il in range(nlay):
                    f.write(
                        f'   {iw+1:2.0f}    {ig+1:2.0f}    {il+1:2.0f}'
                        f'    {bmabs[il, iw, ig]:12.10E}'
                        f'   {wav[iw]:8.4f}\n'
                    )
    log.info("Wrote bmabs text file to %s", filepath)
