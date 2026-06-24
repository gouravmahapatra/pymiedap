# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

"""
T-matrix support for PyMieDAP
=============================

PyMieDAP computes the single-scattering properties of *spherical* particles
internally with Mie theory (:func:`pymiedap.pymiedap.mie_code`).  For
*non-spherical* particles -- e.g. ice crystals -- the scattering matrix has to
come from elsewhere.  The ``tmatrix_ice/`` directory of this repository ships
the Mishchenko & Travis T-matrix (extended boundary condition) Fortran code,
configured for randomly-oriented ice spheroids, which writes a machine-readable
expansion-coefficient file (``*.coeffs``).

This module makes those coefficients first-class PyMieDAP citizens:

* :func:`tmatrix_to_pymiedap_coeffs` converts a T-matrix ``.coeffs`` file into
  the "Meerhoff Mie" expansion-coefficient text format that
  :func:`pymiedap.pymiedap.read_mie_output` / ``module_readmie.file2coefs``
  expect.
* :func:`load_tmatrix_into_aerosol` fills an :class:`~pymiedap.pymiedap.Aerosols`
  object with the converted coefficients and a single-scattering albedo, ready
  to drop into a :class:`~pymiedap.pymiedap.Layer`.
* :func:`delta_m_truncate` applies vector delta-M truncation to any aerosol's
  expansion so that strongly forward-peaked phase functions (large droplets or
  crystals) stay within the angular resolution of the doubling-adding solver.
* :func:`run_tmatrix` is a thin convenience wrapper that runs a compiled
  ``tmatrix_ice`` binary (see ``tmatrix_ice/README.md`` for how to build it).

File formats
------------
**T-matrix ``.coeffs``** (Mishchenko-Travis; see ``tmatrix_ice/README.md``)::

    <single-scattering-albedo>D+00      <Lmax>
    alpha1  alpha2  alpha3  alpha4  beta1  beta2     # order l = 0
    alpha1  alpha2  alpha3  alpha4  beta1  beta2     # order l = 1
    ...

The six columns are the generalised-spherical-function expansion coefficients
in the same order PyMieDAP uses, and they carry the ``(2l+1)`` factor
(``alpha1_0 = 1``).  Adjacent columns can be glued together when a value is
negative (``...D+01-0.67...D-01``), so parsing is done with a regex.

**PyMieDAP ("Meerhoff Mie") format** read by ``readsc``::

    line 1               : title, first 23 chars == " EXPANSION COEFFICIENTS"
    lines 2-7  (6 lines) : free-text header (skipped)
    line 8               : single-scattering albedo, FORMAT(32X, E25.14)
    lines 9-14 (6 lines) : free-text header (skipped; last = column header)
    line 15+             : FORMAT(i4, 6f19.14) -> l a1 a2 a3 a4 b1 b2
                           terminated by a line with l = -1

and the matrix is filled as::

    alfbet(1,1,l)=a1  alfbet(2,2,l)=a2  alfbet(3,3,l)=a3  alfbet(4,4,l)=a4
    alfbet(1,2,l)=alfbet(2,1,l)=b1      alfbet(3,4,l)=-alfbet(4,3,l)=b2
"""

import os
import re
import subprocess

import numpy as np

__all__ = [
    "tmatrix_to_pymiedap_coeffs",
    "load_tmatrix_into_aerosol",
    "delta_m_truncate",
    "run_tmatrix",
]

# Matches Fortran "D"/"E" floats, including ones glued to the previous column
# by a leading minus sign.
_NUM_RE = re.compile(r"[+-]?\d*\.\d+[DdEe][+-]?\d+")


def _to_floats(line):
    return [float(t.replace("D", "E").replace("d", "E"))
            for t in _NUM_RE.findall(line)]


def tmatrix_to_pymiedap_coeffs(src, dst, lam_um=0.5, nr=1.3117, ni=1.0e-8,
                               reff=1.0, veff=0.1):
    """Convert a T-matrix ``.coeffs`` file to a PyMieDAP coefficient file.

    Parameters
    ----------
    src : str
        Path to the T-matrix ``.coeffs`` file.
    dst : str
        Path to write the PyMieDAP-format expansion-coefficient file.
    lam_um, nr, ni, reff, veff : float
        Metadata written into the (otherwise ignored) header lines, for
        provenance only.  They do **not** affect the coefficients.

    Returns
    -------
    albedo : float
        Single-scattering albedo read from the T-matrix file.
    nrows : int
        Number of expansion orders written (l = 0 .. nrows-1).
    lmax_declared : int
        The Lmax value declared on the first line of the source file.
    """
    with open(src) as f:
        first = f.readline().split()
        albedo = float(first[0].replace("D", "E").replace("d", "E"))
        lmax_declared = int(first[1])
        rows = []
        for line in f:
            vals = _to_floats(line)
            if len(vals) >= 6:
                rows.append(vals[:6])
    nrows = len(rows)

    # The albedo is read with FORMAT(32X, E25.14): columns 1-32 are the label,
    # the number occupies columns 33-57.
    ssa_prefix = (" single scattering albedo  a  =").ljust(32)
    ssa_line = ssa_prefix + " {:.14E}".format(albedo)

    with open(dst, "w") as f:
        f.write(" EXPANSION COEFFICIENTS SCATTERING MATRIX\n")
        f.write(" CONVERTED FROM T-MATRIX (Mishchenko-Travis)\n")
        f.write(" lambda= {:11.7f} Re(m) = {:11.7f} Im(m) = {:11.7f}\n"
                .format(lam_um, nr, ni))
        f.write(" Size distribution : (see source)\n")
        f.write(" reff  = {:11.7f} veff  = {:11.7f}\n".format(reff, veff))
        f.write(" rmin  =   0.0000000 rmax  =   0.0000000\n")
        f.write(" (header line - ignored by reader)\n")
        # ^ exactly 6 filler lines between title and the albedo line, so that
        #   readsc's "DO 10 i=1,6" lands the albedo read on the 8th line.
        f.write(ssa_line + "\n")
        f.write(" asymmetry parameter <cos th> =  0.0\n")
        f.write(" extinction efficiency  Qext  =  0.0\n")
        f.write(" scattering efficiency  Qsca  =  0.0\n")
        f.write(" geometrical cross section G  =  0.0\n")
        f.write(" average volume            V  =  0.0\n")
        f.write("   l       alpha1             alpha2             alpha3"
                "             alpha4             beta1              beta2\n")
        for l, row in enumerate(rows):
            f.write("{:4d}".format(l)
                    + "".join("{:19.14f}".format(v) for v in row) + "\n")
        f.write("{:4d}".format(-1)
                + "".join("{:19.14f}".format(0.0) for _ in range(6)) + "\n")
    return albedo, nrows, lmax_declared


def load_tmatrix_into_aerosol(aero, coeff_files, albedos=None,
                              workdir=None, **convert_kwargs):
    """Load T-matrix ``.coeffs`` files into an :class:`Aerosols` object.

    One coefficient file is required per wavelength of the model.  The files
    are converted to PyMieDAP format (in ``workdir``, default: alongside the
    sources) and loaded via ``Aerosols.load_coefs``.  The single-scattering
    albedo of each file is stored on the aerosol so that ``mix_aerosols``
    reproduces it.

    Parameters
    ----------
    aero : pymiedap.pymiedap.Aerosols
        Aerosol object to fill (its ``coefs``/``ncoefs``/``ssalb``/``ssca``
        are overwritten).
    coeff_files : list of str
        T-matrix ``.coeffs`` files, one per wavelength.
    albedos : list of float, optional
        Single-scattering albedos.  If omitted, they are read from the files.
    workdir : str, optional
        Directory for the converted files.  Defaults to the directory of each
        source file.
    **convert_kwargs :
        Passed through to :func:`tmatrix_to_pymiedap_coeffs` (lam_um, nr, ...).

    Returns
    -------
    aero : pymiedap.pymiedap.Aerosols
        The same object, for chaining.
    """
    nwvl = len(coeff_files)
    converted = []
    read_albedos = []
    for cf in coeff_files:
        base = os.path.splitext(os.path.basename(cf))[0]
        out_dir = workdir or os.path.dirname(os.path.abspath(cf))
        dst = os.path.join(out_dir, base + "_pymiedap.coeffs")
        alb, _, _ = tmatrix_to_pymiedap_coeffs(cf, dst, **convert_kwargs)
        converted.append(dst)
        read_albedos.append(alb)

    aero.load_coefs(converted)
    if albedos is None:
        albedos = read_albedos
    albedos = np.asarray(albedos, dtype=float)
    aero.sext = np.ones(nwvl)
    aero.ssca = albedos.copy()
    aero.ssalb = albedos.copy()
    aero.f = 1.0
    return aero


def delta_m_truncate(aero, M):
    """Apply vector delta-M truncation (Wiscombe 1977) to ``aero.coefs``.

    A forward-peaked phase function needs O(size parameter) Legendre terms; the
    doubling-adding solver represents at most ~2*nmug.  Delta-M subtracts an
    (unpolarised) forward delta-peak of fractional weight ``f = chi_M`` and
    rescales the optical thickness and single-scattering albedo to conserve the
    truncated energy.  Only the (1,1) phase-function moments get the
    ``(2l+1)*f`` shift; all other Greek-matrix coefficients are scaled by
    ``1/(1-f)``.  The expansion coefficients carry the ``(2l+1)`` factor, so
    ``f = alpha1_M / (2M+1)``.

    The aerosol's ``coefs``, ``ncoefs``, ``ssalb`` and ``ssca`` are modified
    in place.  Returns the per-wavelength optical-thickness scaling factor
    ``(1 - f*omega)`` to multiply the layer ``tau`` with.

    Note
    ----
    For numerical stability of the doubling the truncation fraction ``f``
    should stay below ~0.2 (a negative truncated phase function otherwise makes
    the doubling diverge), which sets a lower bound on ``M`` -- and hence on
    ``nmug`` (>= M/2) -- for large particles.
    """
    nwvl = aero.coefs.shape[0]
    tau_scale = np.ones(nwvl)
    for z in range(nwvl):
        nc = int(aero.ncoefs[z])
        Mu = min(M, nc - 1)
        l = np.arange(Mu)
        f = aero.coefs[z, 0, 0, Mu] / (2 * Mu + 1)        # f = chi_M
        c = aero.coefs[z].copy()
        c[:, :, :Mu] /= (1.0 - f)                          # scale all elements
        c[0, 0, :Mu] = (aero.coefs[z, 0, 0, :Mu]
                        - (2 * l + 1) * f) / (1.0 - f)     # shift (1,1) only
        c[:, :, Mu:] = 0.0
        aero.coefs[z] = c
        w = aero.ssalb[z]
        w_new = (1.0 - f) * w / (1.0 - f * w)
        aero.ssalb[z] = w_new
        aero.ssca[z] = w_new * aero.sext[z]   # keep ssca consistent for mixing
        aero.ncoefs[z] = Mu
        tau_scale[z] = (1.0 - f * w)
    return tau_scale


def run_tmatrix(binary, workdir=None, output_basename=None):
    """Run a compiled ``tmatrix_ice`` binary and return its output paths.

    The Fortran T-matrix code is not an f2py extension; it is a standalone
    program whose input is hard-coded in its ``INPUT DATA`` block (wavelength,
    refractive index, shape, size).  Build it as described in
    ``tmatrix_ice/README.md``::

        gfortran -std=legacy -O2 -w -c lpq.f    -o lpq.o
        gfortran -std=legacy -O2 -w -c tmq.lp.f -o tmq.lp.o
        gfortran -std=legacy -O2 tmq.lp.o lpq.o -o tmatrix_ice

    It writes two files in the working directory: ``test`` (human-readable
    diagnostics) and ``tmatr.write`` (the machine-readable coefficients, the
    ``.coeffs`` format consumed by :func:`tmatrix_to_pymiedap_coeffs`).

    Parameters
    ----------
    binary : str
        Path to the compiled ``tmatrix_ice`` executable.
    workdir : str, optional
        Directory to run in (defaults to the binary's directory).
    output_basename : str, optional
        If given, ``tmatr.write`` is renamed to ``<output_basename>.coeffs``
        and its path returned; otherwise the raw ``tmatr.write`` path.

    Returns
    -------
    coeffs_path : str
        Path to the produced ``.coeffs`` (or ``tmatr.write``) file.
    """
    binary = os.path.abspath(binary)
    workdir = workdir or os.path.dirname(binary)
    subprocess.run([binary], cwd=workdir, check=True)
    raw = os.path.join(workdir, "tmatr.write")
    if not os.path.exists(raw):
        raise FileNotFoundError(
            "tmatrix run did not produce 'tmatr.write' in {}".format(workdir))
    if output_basename:
        dst = os.path.join(workdir, output_basename + ".coeffs")
        os.replace(raw, dst)
        return dst
    return raw
