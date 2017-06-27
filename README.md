# PyMieDAP

PyMieDAP (Python Mie Doubling Adding Program) is a package to make light
scattering computations with Mie scattering and radiative transfer computations
with full orders of scattering and taking into account the polarization of the
light scattered.

Full planet modeling at any phase angle is possible.  With the subpackage
exopy, it is also possible to simulate systems with a star, a planet and a
possible moon.

## Getting Started

These instructions will get you a copy of the project up and running on your
local machine for development and testing purposes. See deployment for notes on
how to deploy the project on a live system.

### Prerequisites

PyMieDAP requires 
* numpy
* matplotlib
* scipy
* a FORTRAN compiler ('gfortran' is known to work fine)

but if you are a scientist, you might already have those installed.

Be careful: this has only been tested for Unices and with Python2.7.

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

Once installed and your path set correctly, you can import PyMieDAP and Exopy in a
(i)Python terminal or a script using:

```
import pymiedap.pymiedap as pmd
import pymiedap.exopy as exopy
```

To know more about how to use PyMieDAP, you can refer to the notebook or to the
script 'pymiedap_demo.py'.

## Contributing

Please read
[CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for
details on our code of conduct, and the process for submitting pull requests to
us.

## Authors

* **Loïc Rossi** - TU Delft - *Initial work, Python and Fortran interface* -
* [Gitlab](https://gitlab.com/loic.cg.rossi), [TU Delft](http://homepage.tudelft.nl/k6v71/)
* **Daphne Stam** - TU Delft - *Initial work, Fortran code* -
* **Javier Bersoza** - TU Delft - *Initial work, Exopy* - 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

If you want to use this code in a scientific publication, it would be
appreciated if you cite us:
* 

## References

The method used for Mie and Doubling-Adding calculations can be found in the
following references.
For the Mie scattering:
* de Rooij, W. A. & van der Stap, C. C. A. H. _Expansion of Mie scattering
matrices in generalized spherical functions_, A&A, 1984, 131, 237-248
For the Doubling-Adding:
* de Haan, J. F.; Bosma, P. B. & Hovenier, J. W. _The adding method for multiple
 scattering calculations of polarized light_, A&A, 1987, 183, 371-391
The geometry code is refered to in the following papers:
* Fauchez, T.; Rossi, L. & Stam, D. M. _The O2 A-Band in the Fluxes and
Polarization of Starlight Reflected by Earth-Like Exoplanets_, The Astrophysical
Journal, 2017, 842, 41
