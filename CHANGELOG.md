# PyMieDAP changelog

## Unreleased

### New features
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
