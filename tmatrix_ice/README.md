# T-matrix ice scattering — oblate spheroids, visible (0.5 µm)

Mishchenko & Travis T-matrix code (extended/quad precision) configured for ice
particles and made to compile and run under modern **gfortran**.

## The ice case computed

| Parameter | Value | Notes |
|---|---|---|
| Wavelength `LAM` | 0.500 µm | visible |
| Refractive index `MRR + MRI·i` | 1.3117 + 1.0e-8 i | ice at ~0.5 µm (Warren & Brandt 2008); MRI set to a tiny effectively-non-absorbing floor |
| Shape | oblate spheroid, `NP=-1`, `EPS=2.0` | aspect ratio a/b = 2 (flattened — proxy for ice plates) |
| Size distribution | power law (`NDISTR=3`), Hansen–Travis | |
| Effective radius `reff` | 1.0 µm (equal-surface-area-sphere, `RAT=0.5`) | `VEFF=0.1` |
| Orientation | random | analytically averaged |
| Accuracy `DDELT` | 1e-3 | recommended value |

These are set in the `INPUT DATA` block of `tmq.lp.f` (search for `RAT=0.5D0`).
Change `MRR`/`MRI` for other wavelengths, `EPS` for aspect ratio (`<1` = prolate
columns), `NP=-2` for finite cylinders, `AXMAX` for size.

## Results (see `ice_oblate_0.5um.out`)

- **van der Mee & Hovenier physical-correctness test: SATISFIED**
- Single-scattering albedo `W = 0.99998` (near 1 — ice is essentially
  non-absorbing in the visible)
- Asymmetry parameter `<cos> = 0.809`
- Cross sections `CEXT = CSCA = 6.32 µm²`
- 59 generalized-spherical-function expansion coefficients (ALPHA1–4, BETA1–2)
- Full normalized scattering matrix F11, F22, F33, F44, F12, F34 vs angle.
  Note **F22 ≠ F11** — the nonsphericity signature that Mie spheres cannot
  reproduce (for spheres F22 ≡ F11 identically).

## Output files

- `ice_oblate_0.5um.out` — human-readable diagnostics + expansion coefficients
  + scattering matrix (Fortran unit 6, originally `test`).
- `ice_oblate_0.5um.coeffs` — machine-readable: first line is `albedo  Lmax`,
  followed by the six expansion coefficients per order (Fortran unit 10,
  originally `tmatr.write`). This is the format consumed downstream.

## Using these coefficients in PyMieDAP

The `.coeffs` files are wired into the package via `pymiedap.tmatrix`:

```python
import pymiedap.pymiedap as pmd
from pymiedap.tmatrix import load_tmatrix_into_aerosol, delta_m_truncate

ice = pmd.Aerosols(typ='I')
load_tmatrix_into_aerosol(ice, ['tmatrix_ice/ice_oblate_0.5um.coeffs'])
# -> ice.coefs / ice.ncoefs / ice.ssalb are now set; drop `ice` into a Layer.
```

`pymiedap.tmatrix.tmatrix_to_pymiedap_coeffs()` does the format conversion on
its own, `run_tmatrix()` runs a compiled binary, and `delta_m_truncate()`
keeps strongly forward-peaked expansions within the doubling-adding solver's
angular resolution. A complete worked example (liquid Mie + ice T-matrix
clouds, Earth-as-exoplanet spectra) is in
`examples/earthlike_water_ice_clouds.py`.

## Changes made to compile under gfortran

The uploaded `tmq.lp.f` carried a block of sample console output prepended as
plain text (not Fortran), and both files used vendor-specific quad-precision
spellings from the original (Compaq/DEC) compiler. To build with gfortran:

1. Removed the prepended sample-output header (159 lines) from `tmq.lp.f`.
2. Replaced quad intrinsic names with gfortran's generic equivalents
   (case-insensitive), in both files:
   `QCONJG→CONJG`, `QFLOAT→DFLOAT`, `QSQRT→SQRT`, `QCOS→COS`, `QSIN→SIN`,
   `QCOSH→COSH`, `QSINH→SINH`, `QARCOS→ACOS`, `QATAN→ATAN`, `QABS→ABS`,
   `QMAX1→MAX`, `QIMAG→AIMAG`.
3. `QCMPLX(re,im)` → `CMPLX(re,im,16)` (one occurrence) to preserve quad kind.
4. `Q10.3` quad edit descriptor in a FORMAT → `E10.3`.

Identifiers that merely start with Q (QEXT, QSCA, QGAUSS, QR, QI, …) were left
untouched. `REAL*16`/`COMPLEX*32` and `…q0` literals are accepted by gfortran
as-is.

## How to rebuild and run

```sh
gfortran -std=legacy -O2 -w -c lpq.f    -o lpq.o
gfortran -std=legacy -O2 -w -c tmq.lp.f -o tmq.lp.o   # needs tmq.par.f in cwd (INCLUDE)
gfortran -std=legacy -O2 tmq.lp.o lpq.o -o tmatrix_ice
./tmatrix_ice          # writes ./test and ./tmatr.write
```

(Binary not included — it is Linux/aarch64-specific; rebuild on your machine.)

## Caveats

- The T-matrix (extended boundary condition) method becomes unstable for very
  large size parameters or very large aspect ratios; quad precision raises the
  ceiling but does not remove it. This compact a/b=2, reff≈1 µm case converges
  easily.
- Spheroids/cylinders are smooth rotationally-symmetric proxies — they cannot
  represent true hexagonal ice facets or halos. For large faceted crystals use
  geometric-optics/IGOM or DDA databases (Yang & Liou, Baran).
- The `time = … min` line in the output is a cosmetic artifact of the legacy
  `mclock` timer on this platform; ignore it.
- To feed the `refl.f → interp.f` reflection pipeline, the `.coeffs` file matches
  the generalized-spherical-function expansion of `spher.write`, but check the
  exact line format `refl.f` expects before chaining.
