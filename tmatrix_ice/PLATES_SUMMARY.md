# Reproducing Plates 1-4 of Mishchenko, Travis & Mackowski (JQSRT 1996/1998)

A record of what was reproduced from the review paper "T-matrix computations of
light scattering by nonspherical particles," using the original NASA GISS
Fortran codes (sanitized to build under gfortran) and verified independently
where possible.

## The codes used

| code | role | panels it produced |
|---|---|---|
| `spher.f` | Lorenz-Mie (homogeneous spheres) | Plate 1 black curve, Plate 2a, Plate 4 |
| `bisphere.f` | two-sphere superposition T-matrix | Plate 1 (bispheres) |
| `tmq.lp.f` | orientation-averaged T-matrix (single axisymmetric particle) | Plate 2d, Plate 3 |
| `amplq.lp.f` | fixed-orientation amplitude-matrix T-matrix | Plate 2b, Plate 2c |

All are quad/double-precision Fortran 77. Sanitizing meant stripping prepended
sample-output blocks and translating vendor quad intrinsics (`QSQRT`, `QCONJG`,
`QCMPLX`, `QACOS`, ...) to gfortran's generic names. Every build was checked
against the reference output shipped inside each source file.

## Plate 1 - scattering matrix of randomly oriented bispheres

Component spheres x=5, m=1.5+0.005i, centre separations d = 2r, 2.5r, 4r, 8r,
with the single Mie sphere (x=5) as the black reference.

Result: full reproduction. **Key finding** - F22/F11 is identically 100% for the
single sphere but drops below it for the bispheres, deepest for the touching
pair (d=2r, ~60% near 120 deg) and recovering toward 100% as d grows; by d=8r the
two spheres scatter almost independently. The single-sphere interference
structure dominates every element except F22/F11. This is the signature of
cooperative (dependent) scattering.

Verification: all cases pass the van der Mee & Hovenier test; the d=4r run
reproduces the reference table inside `bisphere.f` to 1e-5.

## Plate 2 - degree of linear polarization, monodisperse particles

-F21/F11 (%) over scattering angle x size parameter (0-30), m=1.53+0.008i.
(a) spheres; (b) oblate spheroid eps=1.7, axis parallel to beam; (c) same,
axis perpendicular; (d) oblate spheroid eps=1.7, random orientation.

Result: full reproduction of all four panels. **Key finding** - the monodisperse
diagrams are dominated by a dense field of sharp interference maxima/minima.
For the fixed-orientation spheroids (b, c) the pattern is even more complex and
the two orientations look completely different. Random orientation (d) is
markedly smoother - orientation averaging already removes much of the
interference structure. This is the motivation for using polydisperse, randomly
oriented ensembles rather than monodisperse fixed particles.

Verification: panel (a) matches `miepython` to ~0.02%; the `amplq.lp.f` build
reproduces its shipped amplitude+phase-matrix test exactly.

## Plate 3 - degree of linear polarization, polydisperse spheroids

-F21/F11 (%) over scattering angle x EFFECTIVE size parameter, randomly oriented
polydisperse spheroids, six shapes (prolate/oblate x eps=1.4, 1.7, 2),
power-law size distribution with v_eff = 0.1 (so x1 = 0.8903 x_eff,
x2 = 1.5654 x_eff), m=1.53+0.008i.

Result: reproduced over the computationally reachable range. **Key finding** -
size averaging on top of orientation averaging removes essentially all residual
interference structure: the diagrams are smooth, with a weak, broadly negative
(blue) polarization through the side-scattering angles and shape-dependent
positive features near backscatter (strongest for oblate eps=2). The
prolate/oblate distinction and the eps trend are visible. Contrasting Plates 2,
3 and 4 confirms the paper's thesis: averaging over sizes and orientations
yields smooth diagrams that enable meaningful comparison of particle types.

Scope limit (honest): this was run on a single machine with a per-call time
ceiling. The T-matrix cost rises steeply with size and asphericity, so the
effective size parameter was capped where a single polydisperse point still
fit in time: oblate (all eps) and prolate eps=1.4 to x_eff=15, prolate eps=1.7
to x_eff=13, prolate eps=2 to x_eff=12. The published plate goes to x_eff=30;
the upper region needs substantially more compute (Mishchenko used the improved
code and significant CPU). Size integration used a modest number of Gaussian
points (NKMAX=6); increasing it and the x_eff range would sharpen the diagrams.

## Plate 4 - scattering-matrix elements, polydisperse spheres

Six elements (log10 F11, F22/F11, F33/F11, F44/F11, F12/F11, F34/F11) over
scattering angle x effective size parameter, polydisperse spheres, same
distribution (v_eff=0.1) and refractive index.

Result: full reproduction. **Key finding** - because of size averaging the
diagrams are smooth (no interference fringes). F11 is strongly forward-peaked;
F22/F11 is identically 100% (the defining property of spheres - contrast with
Plate 1 bispheres and Plate 3 spheroids where it departs from 100%);
F44/F11 = F33/F11 for spheres; F12, F33, F34 show smooth, size-dependent
structure. This is the spherical-particle baseline against which Plates 3, 5-9
measure the effect of nonsphericity.

## Overall

The four plates together tell one story: **monodisperse, fixed-orientation
scattering (Plate 2 a-c) is buried under interference structure that makes
particle comparison hopeless; averaging over orientation (2d), and then over
size (Plates 3, 4), progressively smooths it away.** Nonsphericity then shows up
cleanly - most diagnostically in F22/F11, which is exactly 100% only for
spheres (Plate 4) and departs from it for aggregates (Plate 1) and spheroids
(Plate 3).

Folders: `plate1_full/`, `plate2/`, `plate34/` hold the figures, per-point data,
the configured Fortran variants, the sweep drivers, and the plotting scripts.
Each reproduction was regenerated from saved data in a clean directory as a
check.
