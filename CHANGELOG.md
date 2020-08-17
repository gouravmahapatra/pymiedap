# PyMieDAP changelog

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
