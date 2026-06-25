# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

"""PyMieDAP package."""

# T-matrix and Baum ice-cloud support (pure-Python; safe to import without the
# native modules).
from . import tmatrix
from .tmatrix import (
    tmatrix_to_pymiedap_coeffs,
    load_tmatrix_into_aerosol,
    delta_m_truncate,
    run_tmatrix,
)
from . import baum
from .baum import (
    expand_scattering_matrix,
    phase_matrix_to_coeffs,
    load_baum_into_aerosol,
    read_baum_netcdf,
    convert_baum_netcdf,
    load_baum_coeffs,
    fill_aerosol_from_cache,
)
from . import sscorr
from .sscorr import (
    scattering_geometry,
    single_scattering_local,
    tms_correct_local,
)

__all__ = [
    "__version__",
    "tmatrix",
    "tmatrix_to_pymiedap_coeffs",
    "load_tmatrix_into_aerosol",
    "delta_m_truncate",
    "run_tmatrix",
    "baum",
    "expand_scattering_matrix",
    "phase_matrix_to_coeffs",
    "load_baum_into_aerosol",
    "read_baum_netcdf",
    "convert_baum_netcdf",
    "load_baum_coeffs",
    "fill_aerosol_from_cache",
]
__version__ = "0.1.1"
