# PyMieDAP changelog

## Unreleased

### New features
- Added `pymiedap.baum`, integrating the SSEC / Baum-Yang-Heymsfield ice cloud
  bulk scattering models (severely-roughened, randomly-oriented habit mixtures;
  Baum et al. 2011):
  - `read_baum_netcdf()` reads the full-phase-matrix `.nc`/`.nc.gz` files
    (transparent gzip), applies the `F34 = -P43*P11` sign convention, and
    subsets wavelengths.
  - `expand_scattering_matrix()` / `phase_matrix_to_coeffs()` perform the full
    six-element generalized-spherical-function expansion (non-spherical
    `F22 != F11`, `F44 != F33`), validated to reproduce `module_mie.devel` to
    machine precision on the spherical case and the files' own asymmetry
    parameter to <0.5% across 0.2-2 um.
  - `convert_baum_netcdf()` + `examples/convert_baum_to_pymiedap.py` build
    `.npz` coefficient caches; `load_baum_coeffs()` / `fill_aerosol_from_cache()`
    load them into an `Aerosols` object.
  - `earthlike_water_ice_clouds.py` gained an `EWIC_ICE=tmatrix|baum` switch to
    use Baum GHM ice. Note: D_eff~60 um crystals are so forward-peaked that a
    stable delta-M run needs nmug beyond the compiled `nmuMAX=201`; rebuild with
    a larger `nmuMAX` to run the Baum ice offline.
- Added `rebuild_highres_nmug.py`: a helper that raises `nmuMAX` (default 512,
  for nmug up to 500) consistently across `dap_source/max_incl`,
  `geos_source/max_incl` and the matching `pymiedap.py` constants
  (`_RFOU_DIM*`, default `nmuMAX=` arguments), then recompiles the native
  modules. Supports `--nmug-max`, `--nfou-max`, `--patch-only`, `--restore`.
  Intended for the multi-core Linux cluster build.
- Added `pymiedap.tmatrix`, a module integrating the bundled T-matrix ice code
  (`tmatrix_ice/`) into the package:
  - `tmatrix_to_pymiedap_coeffs()` converts a T-matrix `.coeffs` file to the
    Meerhoff-Mie expansion-coefficient format read by `module_readmie`.
  - `load_tmatrix_into_aerosol()` loads converted coefficients (one file per
    wavelength) straight into an `Aerosols` object.
  - `delta_m_truncate()` applies vector delta-M truncation so strongly
    forward-peaked phase functions (large droplets/crystals) stay within the
    angular resolution of the doubling-adding solver.
  - `run_tmatrix()` runs a compiled `tmatrix_ice` binary.
- Added `examples/earthlike_water_ice_clouds.py`: a worked example that builds
  a two-cloud-layer (liquid water via Mie + ice via T-matrix) Earth-like
  atmosphere from the disk-averaged cloud parameters of Roccetti et al. (2025,
  A&A, arXiv:2504.02048, Table 2) and computes disk-integrated reflected- and
  polarized-light spectra and phase curves. Includes a digitized overlay
  against that paper's Fig. 13.

## V 1.5.0  2019-04-15

This version has a few bug fixes:
### Bug fixes
- There was an error in the calculation of the scattering and absorption
  optical thickness. This was no a problem when the single-scattering albedo
    was 1, but would yield incorrect result with absorbing particles
- There was an issue that could appear for users running Python 2.7, creating
  erroneous results when using computations with multiple wavelentgths
- The use of an odd number of pixels was giving some unexpected behaviour for a
  calculation of asymmetry, now fixed.
- The calculation of the asymmetry of the cloud pattern was incorrect.

### Improvements
- Some optimization of the patchy cloud generation algorithm. Improves cloud
  coverage on the edges of the disk.

## V 1.0.0 2018-01-15
