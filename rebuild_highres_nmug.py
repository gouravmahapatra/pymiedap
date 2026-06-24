#!/usr/bin/env python3
# This file is part of PyMieDAP, released under GNU General Public License.
"""
rebuild_highres_nmug.py
=============================================================================
Raise PyMieDAP's doubling-adding angular-resolution limit so that `nmug` (the
number of Gaussian division points) can go up to a chosen maximum, then
recompile the native Fortran modules.

WHY: strongly forward-peaked particles -- in particular the Baum/SSEC ice
habit mixtures at D_eff ~ 60 um (asymmetry g ~ 0.8) -- need a high `nmug`
(several hundred) for a stable delta-M run. The shipped build is compiled with
`nmuMAX = 201`, which caps `nmug` at 201. This tool bumps `nmuMAX` (default to
512, enough for `nmug = 500`) consistently across the Fortran include files AND
the Python constants that must match the compiled array sizes, then rebuilds.

WHAT IT PATCHES (atomically, with .bak backups):
  * dap_source/max_incl    : nmuMAX
  * geos_source/max_incl   : nmuMAX, nfouMAX
  * pymiedap/pymiedap.py    : every `nmuMAX=201` default, and the hardcoded
                              _RFOU_DIM1/2/3 that must equal the compiled
                              (nmatMAX*nmuMAX, nmuMAX, nfouMAX+1) rfou shape.

MEMORY: read_dap_output allocates rfou of shape
    (4*nmuMAX, nmuMAX, nfouMAX+1)  doubles.
  nmuMAX=512, nfouMAX=1024  ->  ~8.6 GB per process.
  nmuMAX=512, nfouMAX=500   ->  ~4.2 GB per process.
Pick nfouMAX to match the largest delta-M order M you will use (Fourier terms
~ M); lower it to save memory. On a multi-core cluster you run one process per
wavelength, so size for (per-process RAM) x (concurrent processes).

USAGE (run on the target machine, where gfortran/Python are available):
    python rebuild_highres_nmug.py                 # nmuMAX=512, nfouMAX=1024
    python rebuild_highres_nmug.py --nmug-max 512 --nfou-max 1024
    python rebuild_highres_nmug.py --patch-only    # edit files, don't compile
    python rebuild_highres_nmug.py --restore       # revert to .bak and rebuild

This MUST be paired with a recompile -- patching the Python constants without
rebuilding module_geos would corrupt read_dap_output. The default flow does
both; --patch-only is for inspection.
=============================================================================
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

REPO = os.path.dirname(os.path.abspath(__file__))
DAP_INCL = os.path.join(REPO, "dap_source", "max_incl")
GEOS_INCL = os.path.join(REPO, "geos_source", "max_incl")
PYMIEDAP = os.path.join(REPO, "pymiedap", "pymiedap.py")
TARGETS = [DAP_INCL, GEOS_INCL, PYMIEDAP]


def _backup(path):
    bak = path + ".bak"
    if not os.path.exists(bak):
        shutil.copy2(path, bak)


def _restore(path):
    bak = path + ".bak"
    if os.path.exists(bak):
        shutil.copy2(bak, path)
        print("  restored", os.path.relpath(path, REPO))


def _sub_in_file(path, replacements):
    with open(path) as f:
        txt = f.read()
    for pat, repl, n_expected in replacements:
        txt, n = re.subn(pat, repl, txt)
        tag = os.path.relpath(path, REPO)
        if n_expected is not None and n != n_expected:
            print("  WARNING [{}]: pattern {!r} matched {} times "
                  "(expected {})".format(tag, pat, n, n_expected))
        else:
            print("  [{}] {} substitution(s): {!r}".format(tag, n, pat))
    with open(path, "w") as f:
        f.write(txt)


def patch(nmug_max, nfou_max):
    nsup = 4 * nmug_max          # nmatMAX(=4) * nmuMAX
    for p in TARGETS:
        _backup(p)

    # --- Fortran include files ---
    _sub_in_file(DAP_INCL, [
        (r"nmuMAX=\d+", "nmuMAX={}".format(nmug_max), 1),
    ])
    _sub_in_file(GEOS_INCL, [
        (r"nmuMAX=\d+", "nmuMAX={}".format(nmug_max), 1),
        (r"nfouMAX=\d+", "nfouMAX={}".format(nfou_max), 1),
    ])

    # --- Python constants that MUST match the recompiled module_geos ---
    _sub_in_file(PYMIEDAP, [
        (r"nmuMAX=201\b", "nmuMAX={}".format(nmug_max), None),
        (r"_RFOU_DIM1 = \d+", "_RFOU_DIM1 = {}".format(nsup), 1),
        (r"_RFOU_DIM2 = \d+", "_RFOU_DIM2 = {}".format(nmug_max), 1),
        (r"_RFOU_DIM3 = \d+", "_RFOU_DIM3 = {}".format(nfou_max + 1), 1),
    ])

    gb = nsup * nmug_max * (nfou_max + 1) * 8 / 1e9
    print("\nPatched for nmug up to {} (nmuMAX={}, nfouMAX={}).".format(
        nmug_max, nmug_max, nfou_max))
    print("read_dap_output rfou array will be ~{:.1f} GB per process.".format(gb))


def rebuild():
    print("\nRebuilding native modules (python setup.py build_ext --inplace) ...")
    # force a clean rebuild so the patched max_incl is picked up
    for so in os.listdir(REPO):
        if so.endswith(".so"):
            os.remove(os.path.join(REPO, so))
    build_dir = os.path.join(REPO, "build")
    if os.path.isdir(build_dir):
        shutil.rmtree(build_dir)
    subprocess.run([sys.executable, "setup.py", "build_ext", "--inplace"],
                   cwd=REPO, check=True)


def smoke_test():
    print("\nSmoke test ...")
    code = (
        "import module_mie, module_dap, module_geos, module_readmie, "
        "module_mieshell;"
        "import pymiedap.pymiedap as pmd;"
        "print('  native modules OK');"
        "print('  pmd.Model() layers:', list(vars(pmd.Model().layers).keys()))"
    )
    subprocess.run([sys.executable, "-c", code], cwd=REPO, check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nmug-max", type=int, default=512,
                    help="New nmuMAX (max allowed nmug). Default 512 (>=500).")
    ap.add_argument("--nfou-max", type=int, default=1024,
                    help="New geos nfouMAX (Fourier terms ~ delta-M order M). "
                         "Default 1024. Lower it to save memory.")
    ap.add_argument("--patch-only", action="store_true",
                    help="Edit the files but do not recompile.")
    ap.add_argument("--restore", action="store_true",
                    help="Restore the original files from .bak and rebuild.")
    args = ap.parse_args()

    if args.restore:
        print("Restoring original sources from .bak ...")
        for p in TARGETS:
            _restore(p)
        rebuild()
        smoke_test()
        print("\nReverted to the original nmuMAX build.")
        return 0

    if args.nmug_max < 1 or args.nfou_max < 1:
        print("nmug-max and nfou-max must be positive.")
        return 1

    patch(args.nmug_max, args.nfou_max)
    if args.patch_only:
        print("\n--patch-only: skipped recompile. Run "
              "'python setup.py build_ext --inplace' yourself, or re-run "
              "without --patch-only.")
        return 0
    rebuild()
    smoke_test()
    print("\nDone. nmug can now go up to {}.".format(args.nmug_max))
    print("Revert any time with:  python rebuild_highres_nmug.py --restore")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
