# Reproducing Plate 2

Plate 2 (Mishchenko, Travis & Mackowski, *JQSRT* 60, 309-324, 1998) is four
colour diagrams of the **degree of linear polarization for unpolarized incident
light**, −F21/F11 (%), over **scattering angle (0-180 deg) × size parameter
(0-30)**, all at refractive index **m = 1.53 + 0.008i**:

| panel | particle | code used |
|---|---|---|
| (a) | monodisperse spheres | `spher_p2a.f` (Lorenz-Mie) |
| (b) | oblate spheroid ε=1.7, axis **parallel** to incident beam (fixed) | `amplq_plate2.lp.f` |
| (c) | oblate spheroid ε=1.7, axis **perpendicular** to incidence (fixed) | `amplq_plate2.lp.f` |
| (d) | oblate spheroid ε=1.7, **random** orientation | `tmq_p2d.f` |

Panels (b)/(c) need the fixed-orientation amplitude-matrix code (`amplq.lp.f`),
panel (d) the orientation-averaged code (`tmq.lp.f`). Both are quad precision and
were sanitized to build under gfortran (see ../README.md). The single-particle
codes here read the **size parameter from stdin** and were configured for the
Plate 2 physical parameters.

## Build (run from this folder; sources copied from ../ as needed)

```sh
gfortran -std=legacy -O2 -w -c ../lpq_amplq.f -o lpqa.o     # quad LAPACK (amplq companion)
gfortran -std=legacy -O2 -w -c ../lpq.f        -o lpq.o     # quad LAPACK (tmq companion)
cp ../amplq.par.f ../tmq.par.f .                            # INCLUDE files

gfortran -std=legacy -O2 -w spher_p2a.f                 -o spher_p2a   # panel a
gfortran -std=legacy -O2 -w -c amplq_plate2.lp.f -o ap2.o
gfortran -std=legacy -O2 ap2.o lpqa.o                   -o amplq_p2    # panels b,c
gfortran -std=legacy -O2 -w -c tmq_p2d.f         -o td.o
gfortran -std=legacy -O2 td.o lpq.o                     -o tmq_p2d     # panel d
```

amplq_plate2 reads `AXI IPANEL` (IPANEL=2 -> panel b, =3 -> panel c) and writes
`ampl.out` (angle, Z11, Z21, -Z21/Z11). spher_p2a and tmq_p2d read `AXI` and
write their scattering matrix to `spher.print` / `test`.

## Run the sweeps

The driver scripts loop the size parameter 0.5-30 (step 0.5) and append
`size angle dlp` to the panel data files. The T-matrix panels are chunked to
keep each invocation short; `run_chunk*.sh` process a `todo*.txt` list of sizes
until ~35 s elapse, then resume on the next call (call repeatedly until "DONE").

```sh
# panel a (fast, one shot)
: > panelA.dat
for x in $(seq 0.5 0.5 30); do echo $x | ./spher_p2a >/dev/null 2>&1
  awk -v x=$x 'f{print x,$1,(-$4/$2)} /<.*F11.*F33.*F12.*F34/{f=1}' spher.print >> panelA.dat; done

# panels b, c  (uses run_amplq.sh for the cheap end, run_chunk.sh for x>15)
: > panelB.dat; bash run_amplq.sh 2 panelB.dat 0.5 15 0.5
seq 15.5 0.5 30 > todoB.txt; until bash run_chunk.sh 2 panelB.dat todoB.txt | grep -q DONE; do :; done
: > panelC.dat; bash run_amplq.sh 3 panelC.dat 0.5 15 0.5
seq 15.5 0.5 30 > todoC.txt; until bash run_chunk.sh 3 panelC.dat todoC.txt | grep -q DONE; do :; done

# panel d
: > panelD.dat; seq 0.5 0.5 30 > todoD.txt
until bash run_chunk_d.sh | grep -q DONE; do :; done
```

## Plot

```sh
python3 plate2_plot.py     # -> plate2_replication.png
```

## Definitions / conventions

- Degree of linear polarization for unpolarized incidence = −Z21/Z11, where Z is
  the 4×4 phase (Stokes scattering) matrix. For spheres and random orientation
  Z21 = F12, matching the plate's "−F21/F11".
- Spheroid size axis = equal-surface-area-sphere size parameter (RAT computed by
  SAREA; wavelength set to 2π so that size parameter equals the input radius).
- Fixed-orientation geometry: incidence along the lab z-axis; particle symmetry
  axis at Euler angle β = 0 (panel b, axis ∥ beam) or β = 90° (panel c, axis ⟂
  beam); scattering traced in the φ = 0 meridian, so scattering angle = θ.

## Verification

- Panel (a) cross-checked against `miepython` (x=10, m=1.53+0.008i): the −F21/F11
  curve matches to ~0.02 %.
- The fixed-orientation build (`amplq.lp.f`) reproduces the amplitude + phase
  matrix shipped in the original file's test header exactly.
- Qualitative check vs the paper: (a) shows the sphere interference fringes;
  (b)/(c) are more intricate and mutually very different; (d) is markedly
  smoother — exactly the message of Plate 2.

## Resolution note

These diagrams use a size-parameter step of 0.5 (60 rows) and a 1° angle grid.
Halve the step for a closer match to the published smoothness (panel d and b/c
runtime scale roughly with the number of sizes; x≈30 is ~15 s per size).
