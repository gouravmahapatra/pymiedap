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

These instructions will get you a copy of the project up and running on your
local machine for development and testing purposes. See deployment for notes on
how to deploy the project on a live system.

### Prerequisites

PyMieDAP requires 
* Python (2.7 or 3.x)
* numpy
* matplotlib
* scipy
* a FORTRAN compiler ('gfortran' is known to work fine)

but if you are a scientist, you might already have those installed.

Be careful: this has only been tested for Unices, Windows systems might face
issues with compilers.

### Installing

To install, you can use the setup.py script

```
python setup.py install
```

this will install PyMieDAP where python itself is installed. If you want to
install it in your home folder, you can use

```
python setup.py install --home=~

```
If you install it locally, you`ll need to add ~/lib/python (or the path you
chose) to your PYTHONPATH.
Add this to your .bashrc

```
export PYTHONPATH=$PYTHONPATH:~/lib/python
```
then run the install.
Be careful that on Mac, you might need to use another file than .bashrc. Also,
for **Spyder** users, you can set it with Spyder's PYTHONPATH manager.

Once installed and your path set correctly, you can import PyMieDAP and Exopy in a
(i)Python terminal or a script using:

```
import pymiedap.pymiedap as pmd
import pymiedap.exopy as exopy
```

To know more about how to use PyMieDAP, you can refer to the notebook
`pymiedap_tutorial.ipynb`.

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
