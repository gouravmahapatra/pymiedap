# PyMieDAP

PyMieDAP (Python Mie Doubling Adding Program) is a package to make light
scattering computations with Mie scattering and radiative transfer computations
with full orders of scattering and taking into account the polarization of the
light scattered.

Full planet modeling at any phase angle is possible.  With the subpackage
exopy, it is also possible to simulate systems with a star, a planet and a
possible moon.

## Development status

**PyMieDAP** is not longer maintained. I may reply to questions regarding the code, but no updates will be made until further notice.

## Getting Started

The most reliable way to run this repository today is:

1. install the system dependencies with Homebrew
2. create a project-local virtual environment
3. install the Python packages into that virtual environment
4. build the Fortran extension modules in place
5. run the examples, tests, or notebooks from that same virtual environment

This repository has been tested on Unix-like systems. macOS and Linux are the
recommended targets.

## Recommended Setup

The commands below assume:

* macOS
* Homebrew is installed
* the repository is checked out locally

The recommended Python version for this repository is **Python 3.11**.

## Quick Start

### 1. Install System Dependencies

Install the native build tools and runtime libraries with Homebrew:

```bash
xcode-select --install
brew update
brew install python@3.11 gcc pkgconf openblas
```

What these provide:

* `python@3.11`: the recommended Python interpreter for this repository
* `gcc`: provides `gfortran`, which is required to build the Fortran modules
* `pkgconf`: compiler/linker metadata helper
* `openblas`: optimized BLAS/LAPACK library

### 2. Create and Activate a Virtual Environment

Run all PyMieDAP commands from inside a virtual environment. This keeps the
installation isolated from other Python projects.

```bash
cd /path/to/PyMieDAP
$(brew --prefix)/bin/python3.11 -m venv .venv
source .venv/bin/activate
```

After activation, your shell prompt should show `(.venv)`.

Whenever you come back to the project later, reactivate it with:

```bash
cd /path/to/PyMieDAP
source .venv/bin/activate
```

### 3. Install Python Dependencies Into the Virtual Environment

With the virtual environment activated:

```bash
python -m pip install --upgrade pip wheel "setuptools<60"
python -m pip install numpy scipy matplotlib pillow ipykernel jupyterlab
```

### 4. Build the Native PyMieDAP Modules

PyMieDAP depends on compiled native modules such as `module_mie` and
`module_dap`. Build them from the repository root:

```bash
python setup.py build_ext --inplace
```

If the build succeeds, you should see compiled files such as:

```bash
module_mie.cpython-311-darwin.so
module_mieshell.cpython-311-darwin.so
module_readmie.cpython-311-darwin.so
module_dap.cpython-311-darwin.so
module_geos.cpython-311-darwin.so
```

### 5. Verify the Installation

Check that the native modules import correctly:

```bash
python -c "import module_mie, module_mieshell, module_readmie, module_dap, module_geos; print('native modules ok')"
```

Then verify the package imports:

```bash
python -c "import pymiedap.pymiedap as pmd; import pymiedap.exopy as exopy; print('pymiedap ok')"
```

## Running the Code

### Run the Test Suite

From the repository root, with the virtual environment activated:

```bash
python -m unittest discover -s tests -v
```

### Run the Benchmark Script

```bash
python examples/run_pymiedap_benchmark.py
```

### Generate the Lambertian Phase-Curve Figure

```bash
python examples/plot_lambert_phase_curve.py
```

This writes a figure to:

```text
examples/lambert_phase_curve.png
```

### Run a Notebook

Register the virtual environment as a Jupyter kernel:

```bash
python -m ipykernel install --user --name pymiedap-venv --display-name "PyMieDAP (.venv)"
```

Start Jupyter from the repository root:

```bash
jupyter notebook examples/pymiedap_benchmark_updated.ipynb
```

or:

```bash
jupyter lab examples/pymiedap_benchmark_updated.ipynb
```

In Jupyter, select the kernel:

```text
PyMieDAP (.venv)
```

If you use the wrong kernel, the notebook will usually fail with errors such
as:

```text
ModuleNotFoundError: No module named 'module_mie'
```

That means the notebook is not running inside the same virtual environment that
was used to build the native modules.

## Typical Workflow

For day-to-day use, the minimal workflow is:

```bash
cd /path/to/PyMieDAP
source .venv/bin/activate
python setup.py build_ext --inplace
python -m unittest discover -s tests -v
python examples/run_pymiedap_benchmark.py
```

## Non-spherical particles: T-matrix support

PyMieDAP computes single-scattering properties of **spherical** particles
internally with Mie theory. For **non-spherical** particles (e.g. ice
crystals), the scattering matrix has to come from elsewhere. The `tmatrix_ice/`
directory ships the Mishchenko & Travis T-matrix (extended boundary condition)
Fortran code, configured for randomly-oriented ice spheroids, which writes a
machine-readable expansion-coefficient file (`*.coeffs`). The `pymiedap.tmatrix`
module integrates those coefficients into the normal PyMieDAP workflow:

| Function | Purpose |
|---|---|
| `tmatrix_to_pymiedap_coeffs()` | Convert a T-matrix `.coeffs` file to the Meerhoff-Mie expansion-coefficient format read by `module_readmie`. |
| `load_tmatrix_into_aerosol()` | Load converted coefficients (one file per wavelength) straight into an `Aerosols` object. |
| `delta_m_truncate()` | Apply vector delta-M truncation so strongly forward-peaked phase functions (large droplets/crystals) stay within the doubling-adding solver's angular resolution. |
| `run_tmatrix()` | Run a compiled `tmatrix_ice` binary and return the produced `.coeffs` path. |

Minimal use — drop a T-matrix ice scattering matrix into a layer:

```python
import pymiedap.pymiedap as pmd
from pymiedap.tmatrix import load_tmatrix_into_aerosol

ice = pmd.Aerosols(typ='I')
load_tmatrix_into_aerosol(ice, ['tmatrix_ice/ice_oblate_0.5um.coeffs'])
# -> ice.coefs / ice.ncoefs / ice.ssalb are now set; assign `ice` to a Layer.
```

To generate coefficients for other wavelengths/sizes/shapes, edit the
`INPUT DATA` block of the Fortran source and rebuild as described in
`tmatrix_ice/README.md`.

### Earth-as-an-exoplanet example

`examples/earthlike_water_ice_clouds.py` is a complete worked example that
builds a two-cloud-layer atmosphere — a **liquid water cloud** (Mie) below an
**ice cloud** (T-matrix) — from the disk-averaged cloud parameters of Roccetti
et al. (2025, A&A, [arXiv:2504.02048](https://arxiv.org/abs/2504.02048),
Table 2), then computes disk-integrated reflected- and polarized-light spectra
and phase curves and overlays them on that paper's digitized Fig. 13.

```bash
source .venv/bin/activate
EWIC_MODE=demo    python examples/earthlike_water_ice_clouds.py   # fast, r_eff=3 um
EWIC_MODE=offline python examples/earthlike_water_ice_clouds.py   # paper's r_eff=8.6 um
```

`demo` mode runs in seconds with slightly smaller droplets; `offline` mode uses
the paper's exact 8.6 um droplets (heavier: large `nmug`/`M_TRUNC` for stability,
~1-2 min per wavelength).

## Troubleshooting

### `ModuleNotFoundError: No module named 'module_mie'`

Cause:

* the native extensions were not built, or
* the notebook/script is running under the wrong Python interpreter

Fix:

```bash
source .venv/bin/activate
python setup.py build_ext --inplace
python -c "import module_mie; print('module_mie ok')"
```

If this happens inside Jupyter, switch the kernel to `PyMieDAP (.venv)` and
restart the notebook kernel.

### `gfortran` not found

Install the compiler with:

```bash
brew install gcc
```

### Jupyter is not available

Install it into the virtual environment:

```bash
source .venv/bin/activate
python -m pip install jupyterlab ipykernel
```

## Tutorials and Examples

Useful entry points in this repository:

* `pymiedap_tutorial.ipynb`
* `pymiedap_benchmark.ipynb`
* `examples/pymiedap_benchmark_updated.ipynb`
* `examples/run_pymiedap_benchmark.py`
* `examples/plot_lambert_phase_curve.py`
* `examples/earthlike_water_ice_clouds.py` — liquid (Mie) + ice (T-matrix)
  Earth-like clouds; disk-integrated spectra/phase curves vs Roccetti et al.
  (2025) Fig. 13 (see the T-matrix section above)

## Authors

* **Loïc Rossi** - TU Delft/LATMOS - *Initial work, Python and Fortran interface* -
    [Gitlab](https://gitlab.com/loic.cg.rossi), [Website](http://loic.cg.rossi.gitlab.io)
* **Daphne Stam** - TU Delft - *Initial work, Fortran code* -
* **Javier Bersoza** - TU Delft - *Initial work, Exopy* - 

## License

This project is licensed under the GNU GPL and CeCILL-B License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

If you want to use this code in a scientific publication, please cite:
* Rossi, L., Berzosa-Molina, J., Stam, D. M., _PyMieDAP: a Python-Fortran tool for computing fluxes and polarization signals
  of (exo)planets_. Astronomy & Astrophysics, Vol. 616, A147.

  Link is [here](https://doi.org/10.1051/0004-6361/201832859), and
  [also on arXiv](https://arxiv.org/abs/1804.08357)

If you use Exopy, please refer to:
* Berzosa-Molina, J., Rossi, L. and Stam, D. M.; _Traces of exomoons in
  computed flux and polarization phase curves of starlight reflected by
exoplanets_. Astronomy & Astrophysics, in press.
        [Link here](https://doi.org/10.1051/0004-6361/201833320).



## References

The method used for Mie and Doubling-Adding calculations can be found in the
following references.
For the Mie scattering:
* de Rooij, W. A. & van der Stap, C. C. A. H. _Expansion of Mie scattering
matrices in generalized spherical functions_, A&A, 1984, 131, 237-248
For the Doubling-Adding:
* de Haan, J. F.; Bosma, P. B. & Hovenier, J. W. _The adding method for multiple
 scattering calculations of polarized light_, A&A, 1987, 183, 371-391

Some examples of use of PyMieDAP can be found in the following papers:
* Fauchez, T.; Rossi, L. & Stam, D. M. _The O2 A-Band in the Fluxes and
Polarization of Starlight Reflected by Earth-Like Exoplanets_, The Astrophysical
Journal, 2017, 842, 41
* Rossi, L. and Stam, D. M. _Using polarimetry to retrieve cloud coverage of
  Earth-like exoplanets_, Astronomy and Astrophysics, 2017, 607, A57
