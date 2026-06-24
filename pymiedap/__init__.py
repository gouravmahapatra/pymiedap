# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

"""PyMieDAP package."""

# T-matrix support (pure-Python; safe to import without the native modules).
from . import tmatrix
from .tmatrix import (
    tmatrix_to_pymiedap_coeffs,
    load_tmatrix_into_aerosol,
    delta_m_truncate,
    run_tmatrix,
)

__all__ = [
    "__version__",
    "tmatrix",
    "tmatrix_to_pymiedap_coeffs",
    "load_tmatrix_into_aerosol",
    "delta_m_truncate",
    "run_tmatrix",
]
__version__ = "0.1.1"
