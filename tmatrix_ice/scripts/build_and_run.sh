#!/usr/bin/env bash
# Build the Mie + bisphere codes and run all cases needed for the Plate 1 figures.
#
# Requirements: gfortran (any recent version).
# Run from the directory that contains spher.f and bisphere.f, e.g.:
#     bash scripts/build_and_run.sh
#
# Source files expected (already configured for the Plate 1 cases):
#   spher.f     - Lorenz-Mie code, set to monodisperse sphere x=5, m=1.5+0.005i
#   bisphere.f  - Mishchenko & Mackowski two-sphere superposition T-matrix
#
# The 159-line (spher.f) / 89-line (bisphere.f) sample-output headers that ship
# at the top of the original files have ALREADY been stripped in the copies in
# this folder. If you start from a pristine NASA download, strip them first:
#   sed -i '1,71d'  spher.f      # remove prepended test output
#   sed -i '1,89d'  bisphere.f
set -euo pipefail
FFLAGS="-std=legacy -O2 -w"

echo ">> Building Mie reference (spher.f)"
gfortran $FFLAGS spher.f -o spher_mie
./spher_mie                      # writes spher.print (scattering matrix table)

echo ">> Building bisphere.f"
gfortran $FFLAGS bisphere.f -o bisphere

echo ">> Sweeping centre separation d = 2r, 2.5r, 4r, 8r"
# component radius = 5 (LAM=2*pi in the source, so size parameter == radius).
# d/r = R12 / R(1). R12 is set on the line  '      R12=...D0'  in bisphere.f.
for spec in "10D0:d2r" "12.5D0:d2p5r" "20D0:d4r" "40D0:d8r"; do
  val="${spec%%:*}"; tag="${spec##*:}"
  sed -E "s/^      R12=[0-9.]+D0/      R12=${val}/" bisphere.f > bs_${tag}.f
  gfortran $FFLAGS bs_${tag}.f -o bs_${tag}
  ./bs_${tag} >/dev/null 2>&1
  cp bisphere.print bisphere_${tag}.print
  echo "   ${tag}: R12=${val}  ->  bisphere_${tag}.print"
done

echo ">> Done. Now run:  python3 scripts/plate1.py"
