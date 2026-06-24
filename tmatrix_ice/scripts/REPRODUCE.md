# Reproducing the Plate 1 figures

Two steps, run from the `tmatrix_ice/` folder (which holds `spher.f` and
`bisphere.f`, already patched to compile under gfortran and pre-configured for
the Plate 1 cases).

```sh
bash scripts/build_and_run.sh     # compiles Fortran, runs Mie + 4 bisphere cases
python3 scripts/plate1.py         # parses outputs, verifies, makes the figures
```

## What each step does

`build_and_run.sh`
- compiles `spher.f` -> `spher_mie`, runs it -> `spher.print`
  (monodisperse sphere, size parameter x = 5, m = 1.5 + 0.005i)
- compiles `bisphere.f` and runs it four times, editing the centre separation
  `R12` for d = 2r, 2.5r, 4r, 8r -> `bisphere_d2r/d2p5r/d4r/d8r.print`

`plate1.py`
- parses the Mishchenko-format scattering-matrix tables
- (optional) cross-checks `spher.f` against `miepython` if it is installed
- writes `plate1_mie_reference.png` + `.csv` (the black single-sphere curve)
- writes `plate1_full_replication.png` (4 bispheres + sphere reference)

## Requirements

- `gfortran` (any recent version)
- Python: `numpy`, `matplotlib`; `miepython` optional (only for the cross-check)

```sh
pip install numpy matplotlib miepython
```

## Verification built in

- All bisphere runs print "TEST OF VAN DER MEE & HOVENIER IS SATISFIED".
- The d = 4r run reproduces the reference table shipped at the top of the
  original `bisphere.f` to ~1e-5 (integral quantities match exactly).
- `plate1.py` reports max|Δ| between `spher.f` and `miepython` for the
  normalized ratios (~6e-5 in this case).

## Starting from a pristine NASA download

The copies here have had their prepended sample-output blocks removed (they are
plain text, not Fortran). For a fresh download, strip them first:

```sh
sed -i '1,71d' spher.f       # remove the 71-line sample output
sed -i '1,89d' bisphere.f    # remove the 89-line sample output
```

`spher.f` additionally needed its vendor quad-intrinsic spellings translated for
gfortran — but note the Mie `spher.f` used here is double precision and compiles
as-is; the quad-precision translation applies to the separate single-particle
T-matrix code `tmq.lp.f` (see ../README.md).
