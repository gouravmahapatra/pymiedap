#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_baum_to_pymiedap.py
=============================================================================
Drop-in converter: place the SSEC / Baum ice-cloud bulk-scattering files in a
folder, run this script, and get PyMieDAP-ready scattering-coefficient caches.

The SSEC "full phase matrix" files
(`*_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc[.gz]`, habit models:
GeneralHabitMixture, SolidColumns, AggregateSolidColumns) give the six
independent scattering-matrix elements at 498 scattering angles, for 445
wavelengths (0.2-100 um) and 23 effective diameters (10-120 um).

For each file found, this tool reads it (gzip is handled transparently),
expands the phase matrix into generalized-spherical-function expansion
coefficients in PyMieDAP's exact convention (validated to reproduce the file's
own asymmetry parameter to <0.5%), and writes a compact ``.npz`` cache holding
``coefs`` / ``ncoefs`` / ``ssalb`` / ``asym`` per wavelength at the requested
effective diameter.

Use the cache in a model with::

    import pymiedap.pymiedap as pmd
    from pymiedap.baum import fill_aerosol_from_cache
    from pymiedap.tmatrix import delta_m_truncate

    ice = pmd.Aerosols(typ='I')
    fill_aerosol_from_cache(ice, 'GeneralHabitMixture_..._Deff60.npz',
                            wavelengths_um=my_wvls)
    delta_m_truncate(ice, M=...)        # tame the forward peak for doubling-adding

Usage
-----
    python examples/convert_baum_to_pymiedap.py [--dir DIR] [--deff 60]
        [--wlmin 0.2] [--wlmax 2.0] [--ngauss 2000] [--out OUTDIR]

With no arguments it scans the repository root for the three habit files and
writes caches into ``examples/baum_cache/`` at D_eff = 60 um, 0.2-2.0 um.
"""

import argparse
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

from pymiedap.baum import convert_baum_netcdf   # noqa: E402

PATTERN = "*FullPhaseMatrix.nc*"   # matches .nc and .nc.gz


def main():
    ap = argparse.ArgumentParser(description="Convert SSEC/Baum ice NetCDF files "
                                             "to PyMieDAP coefficient caches.")
    ap.add_argument("--dir", default=REPO,
                    help="Folder to scan for Baum *FullPhaseMatrix.nc[.gz] files "
                         "(default: repository root).")
    ap.add_argument("--out", default=os.path.join(HERE, "baum_cache"),
                    help="Output directory for the .npz caches.")
    ap.add_argument("--deff", type=float, default=60.0,
                    help="Effective diameter in um (nearest grid value used).")
    ap.add_argument("--wlmin", type=float, default=0.2, help="Min wavelength [um].")
    ap.add_argument("--wlmax", type=float, default=2.0, help="Max wavelength [um].")
    ap.add_argument("--ngauss", type=int, default=2000,
                    help="Gauss-Legendre nodes for the expansion (forward peak).")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.dir, PATTERN)))
    if not files:
        print("No Baum files matching '{}' found in {}".format(PATTERN, args.dir))
        print("Place the GHM / SolidColumns / AggregateSolidColumns "
              "*FullPhaseMatrix.nc[.gz] files there and re-run.")
        return 1

    os.makedirs(args.out, exist_ok=True)
    print("Found {} file(s); converting at D_eff={:.0f} um, {:.2f}-{:.2f} um:"
          .format(len(files), args.deff, args.wlmin, args.wlmax))
    for f in files:
        print("\n== {} ==".format(os.path.basename(f)))
        base = os.path.basename(f)
        for suff in (".nc.gz", ".nc"):
            if base.endswith(suff):
                base = base[: -len(suff)]
                break
        out_npz = os.path.join(args.out, "{}_Deff{:.0f}.npz".format(base, args.deff))
        convert_baum_netcdf(f, out_npz=out_npz, deff_um=args.deff,
                            wavelength_range_um=(args.wlmin, args.wlmax),
                            ngauss=args.ngauss)
    print("\nDone. Caches written to {}".format(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
