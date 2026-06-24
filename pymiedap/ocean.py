# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

"""Ocean surface support for PyMieDAP.

Physical overview
-----------------
A wind-ruffled ocean surface is a two-component optical system:

1. **Rough air-water interface** — specular reflection and refraction at
   wind-tilted microfacets whose slope distribution follows Cox & Munk (1954).
   Each facet contributes a Fresnel reflection/transmission weighted by the
   slope probability density function (PDF) and a Smith/Sancer shadowing
   factor.  The result is an angle-dependent, polarizing operator.

2. **Subinterface water body** — the water column below the interface scatters
   and absorbs photons.  Pure seawater scattering is strongly forward-peaked
   with a Rayleigh-like polarisation signature (anisotropic density fluctuation
   model).  The water body is solved here with a polarized adding-doubling
   scheme on its own discrete-ordinate angular grid.

These two components are coupled following the formalism of
**Trees & Stam (2019, A&A 626, A129)**, whose Appendix A gives the clean-ocean
reflection operator:

    R_CO = R_I + T_I* · R_W · (I − R_I* · R_W)^{-1} · T_I       (TS19, A.1)

where:
  - R_I  = rough air→water reflection operator (upwelling, viewed from above)
  - T_I  = rough air→water downward transmission operator
  - R_I* = rough water→air reflection operator (total-internal-reflection
           effects included; this is the "starred" conjugate path)
  - T_I* = rough water→air upward transmission operator
  - R_W  = polarized reflection operator of the water body (maps downwelling
           water streams to upwelling water streams)
  - I    = identity on the water-stream space

The factor (I − R_I* · R_W)^{-1} accounts for the infinite series of
multiple reflections between the interface (from below) and the water body.

The final ocean operator mixes the clean-ocean operator with a Lambertian,
non-polarizing whitecap (foam) term using the Monahan & O'Muircheartaigh
(1980) wind-speed-dependent whitecap fraction *f_w*:

    R_ocean = f_w · R_foam + (1 − f_w) · R_CO

Polarisation sign convention
----------------------------
PyMieDAP uses the IAU/Stokes [I, Q, U, V] convention where Q > 0 means
electric-field oscillations preferentially parallel to the scattering/
reference plane and Q < 0 means oscillations preferentially perpendicular.
This matches the classical DAP Rayleigh convention (Hovenier & van der Mee
1983, A&A 128, 1; de Haan et al. 1987).

For Fresnel reflection, s-polarisation (perpendicular to the plane of
incidence) is always reflected more strongly than p-polarisation, so the
reflected beam acquires Q < 0.  Throughout this module the Fresnel off-
diagonal element is defined as:

    b = (R_p − R_s) / 2    →  b < 0 because R_p < R_s

so that the Mueller element M[1,0] = b < 0 maps unpolarised I to negative Q,
consistent with the Rayleigh matrix convention used in the rest of PyMieDAP.

Reference-plane (meridian-plane) rotations
-------------------------------------------
All Stokes vectors in PyMieDAP are defined relative to the local meridian
plane: the half-plane that contains the local vertical (z-axis) and the ray
direction.  Whenever a scattering or reflection event changes the reference
plane (i.e. the pre- and post-event meridian planes differ), the Stokes vector
must be rotated before and after applying the interaction Mueller matrix.

The rotation operator is L(ξ), the standard Mueller rotation matrix:

    L(ξ) = diag(1) ⊕ [[cos2ξ, sin2ξ], [−sin2ξ, cos2ξ]] ⊕ diag(1)

Two rotation angles are needed for each interaction:
  - ξ₁ : from the incoming meridian plane to the scattering/facet plane
  - ξ₂ : from the scattering/facet plane to the outgoing meridian plane

The full Mueller operator is  L(−ξ₂) · M_interaction · L(−ξ₁).

This module implements a Trees & Stam (2019)-style ocean bottom boundary
for the doubling-adding solver.  The user-facing object is
:class:`OceanSurface`, which can be assigned to ``model.surface`` before
calling ``pmd.compute``.

Compared with the first rough-ocean branch, this version replaces the
lightweight diffuse water-colour closure with an internal polarized water-body
solver and an interface/water coupling step.  The implemented clean-ocean
operator follows Appendix A of Trees & Stam (2019):

    R_CO = R_I + T_I^* R_W (I - R_I^* R_W)^-1 T_I

where R_I and T_I are rough air-water interface reflection/transmission
operators, R_I^* and T_I^* are the corresponding water-to-air operators, and
R_W is the polarized reflection operator of the subinterface water body.  The
final ocean operator adds wind-generated foam as a Lambertian, non-polarizing
whitecap term.

The water body is solved as a homogeneous, pure-seawater slab over a black or
Lambertian bottom by a discrete-ordinate adding-doubling calculation.  Pure
seawater scattering is approximated with a depolarized Rayleigh Mueller matrix.
The implementation uses rectangular supermatrices internally: the atmosphere
uses the DAP angular streams, while the ocean body uses its own, usually denser,
water-stream quadrature.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from typing import Callable, Literal, Tuple

import numpy as np


_EPS = 1.0e-14  # Small floor value used throughout to avoid division by zero


# ---------------------------------------------------------------------------
# Generic Stokes/vector geometry helpers
# ---------------------------------------------------------------------------


def _rotation_matrix(angle: float) -> np.ndarray:
    """Mueller rotation matrix L(angle) for the [I, Q, U, V] Stokes convention.

    Physical meaning
    ----------------
    When a Stokes vector is defined relative to one reference plane (e.g. the
    meridian plane of the incoming ray) and we wish to express it relative to a
    rotated reference plane (e.g. the scattering plane), we must apply L(ξ)
    where ξ is the signed rotation angle from the old plane to the new plane,
    measured around the ray direction.

    The rotation acts only on the linear-polarisation components Q and U; it
    leaves total intensity I and circular polarisation V unchanged:

        L(ξ) = | 1    0       0      0 |
                | 0  cos2ξ  sin2ξ    0 |
                | 0 -sin2ξ  cos2ξ    0 |
                | 0    0       0      1 |

    Note the factor of 2 in the trigonometric arguments: this arises because
    linear polarisation is invariant under a 180° rotation of the reference
    plane (it is a spin-2 quantity), so a physical rotation by ξ rotates the
    Stokes (Q, U) pair by 2ξ.

    Parameters
    ----------
    angle : float
        Rotation angle ξ in radians (positive = counterclockwise when viewed
        along the direction of light propagation).

    Returns
    -------
    L : np.ndarray, shape (4, 4)
        Mueller rotation matrix.

    References
    ----------
    Hovenier & van der Mee (1983), A&A 128, 1, eq. (2.6).
    """

    c = np.cos(2.0 * angle)   # cos(2ξ) — factor of 2 for spin-2 quantity
    s = np.sin(2.0 * angle)   # sin(2ξ)
    return np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, c, s, 0.0],
            [0.0, -s, c, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )


def _safe_normalize(vec: np.ndarray) -> np.ndarray:
    """Return a unit vector parallel to ``vec``, or zero if ``vec`` is tiny.

    This guard prevents NaN propagation when degenerate geometries arise (e.g.
    a ray directed exactly along the z-axis, where the cross product with z
    vanishes, making the meridian plane undefined).

    Parameters
    ----------
    vec : np.ndarray, shape (3,)
        Input 3-D vector.

    Returns
    -------
    np.ndarray, shape (3,)
        Normalised vector, or zero vector if ``|vec| < _EPS``.
    """
    norm = np.linalg.norm(vec)
    if norm < _EPS:
        return np.zeros_like(vec)
    return vec / norm


def _signed_angle_between_plane_normals(
    old_normal: np.ndarray, new_normal: np.ndarray, ray_direction: np.ndarray
) -> float:
    """Signed Stokes-reference rotation angle ξ around ``ray_direction``.

    Physical meaning
    ----------------
    Each of the two planes is specified by its normal vector (a vector
    perpendicular to the plane and perpendicular to the ray direction).
    ``old_normal`` is the normal to the current reference plane; ``new_normal``
    is the normal to the target reference plane.  The function returns the
    signed angle ξ such that rotating ``old_normal`` by ξ around
    ``ray_direction`` gives ``new_normal``.

    The sign is determined using the triple product
        y = ray · (old_normal × new_normal)
    and the usual ``arctan2`` convention: positive ξ means a
    counterclockwise rotation when viewed along the direction of propagation.

    This angle is then passed to ``_rotation_matrix(ξ)`` to produce the
    appropriate Mueller rotation matrix L(ξ).

    Parameters
    ----------
    old_normal : np.ndarray, shape (3,)
        Normal to the current (old) reference plane.
    new_normal : np.ndarray, shape (3,)
        Normal to the desired (new) reference plane.
    ray_direction : np.ndarray, shape (3,)
        Direction of light propagation (defines the rotation axis).

    Returns
    -------
    float
        Signed angle ξ in radians in the range (−π, π].
    """

    old_normal = _safe_normalize(old_normal)
    new_normal = _safe_normalize(new_normal)
    ray_direction = _safe_normalize(ray_direction)
    # After normalisation the zero-vector guard fires if the input was degenerate.
    if np.linalg.norm(old_normal) < _EPS or np.linalg.norm(new_normal) < _EPS:
        return 0.0
    # x = cos ξ: projection of new_normal onto old_normal
    x = float(np.clip(np.dot(old_normal, new_normal), -1.0, 1.0))
    # y = sin ξ: signed component of the cross product along the ray axis
    y = float(np.dot(ray_direction, np.cross(old_normal, new_normal)))
    return float(np.arctan2(y, x))


def _meridian_plane_normal(ray_direction: np.ndarray) -> np.ndarray:
    """Return the normal to the meridian plane of ``ray_direction``.

    The meridian plane is the vertical half-plane containing both the local
    zenith direction (z-axis) and the ray.  Its outward normal is:

        n_mer = (z × k) / |z × k|

    where k is the ray direction and z = [0, 0, 1].

    This normal serves as the reference-plane normal for Stokes vectors: Q > 0
    means the E-field oscillates in the direction of ``n_mer`` (within the
    meridian plane), and Q < 0 means oscillations are perpendicular to it.

    Edge case: when the ray points exactly along ±z (zenith/nadir), the cross
    product vanishes and the meridian plane is undefined.  By convention the
    y-z plane normal [0, 1, 0] is returned, giving a stable, reproducible
    reference for the vertical directions that are always used for the
    extra μ = 1 stream in PyMieDAP.

    Parameters
    ----------
    ray_direction : np.ndarray, shape (3,)
        Unit vector giving the propagation direction of the ray.

    Returns
    -------
    np.ndarray, shape (3,)
        Unit normal to the meridian plane.
    """
    z = np.array([0.0, 0.0, 1.0])
    normal = np.cross(z, ray_direction)
    if np.linalg.norm(normal) < _EPS:
        # At zenith/nadir the meridian plane is undefined.  The x-z plane is a
        # stable convention and avoids numerical singularities.
        normal = np.array([0.0, 1.0, 0.0])
    return _safe_normalize(normal)


def _direction(mu: float, azimuth: float, hemi: Literal["up", "down"]) -> np.ndarray:
    """Construct a unit direction vector from cosine-zenith and azimuth angle.

    Uses the standard spherical-coordinate parametrisation where μ = cos θ and
    φ = azimuth.  The z-component sign encodes the hemisphere:
      - "up"   → z > 0 (upwelling, toward zenith)
      - "down" → z < 0 (downwelling, toward nadir)

    Parameters
    ----------
    mu : float
        |cos θ|, i.e. the absolute cosine of the polar angle (0 = horizontal,
        1 = vertical).  Clipped to [0, 1].
    azimuth : float
        Azimuth angle φ in radians.  The incoming reference direction used
        throughout this module has azimuth = 0.
    hemi : {"up", "down"}
        Which hemisphere the direction belongs to.

    Returns
    -------
    np.ndarray, shape (3,)
        Unit direction vector [sin θ cos φ, sin θ sin φ, ±cos θ].
    """
    mu = float(np.clip(mu, 0.0, 1.0))
    s = np.sqrt(max(0.0, 1.0 - mu * mu))   # sin θ = sqrt(1 − μ²)
    z = mu if hemi == "up" else -mu         # sign encodes hemisphere
    return np.array([s * np.cos(azimuth), s * np.sin(azimuth), z], dtype=float)


def _stokes_block(mat: np.ndarray, nmat: int) -> np.ndarray:
    """Extract the leading (nmat × nmat) sub-block of a (4 × 4) Mueller matrix.

    PyMieDAP can run in reduced-Stokes modes: nmat=1 (intensity only, I),
    nmat=3 (linear polarisation, [I, Q, U]), or nmat=4 (full, [I, Q, U, V]).
    This helper extracts the appropriate block so that the Fourier projection
    loops do not carry unnecessary zeros for V when nmat < 4.

    Parameters
    ----------
    mat : np.ndarray, shape (4, 4)
        Full 4×4 Mueller matrix.
    nmat : int
        Number of Stokes components to retain (1, 3, or 4).

    Returns
    -------
    np.ndarray, shape (nmat, nmat)
        Upper-left sub-block of ``mat``.
    """
    if nmat == 4:
        return mat
    if nmat == 3:
        return mat[:3, :3]
    if nmat == 1:
        return mat[:1, :1]
    raise ValueError("nmat must be 1, 3, or 4")


# ---------------------------------------------------------------------------
# Fresnel, Cox-Munk, shadowing and whitecaps
# ---------------------------------------------------------------------------


def _fresnel_amplitude_coefficients(mu_i: float, n1: float, n2: float) -> Tuple[float, float, float, bool]:
    """Compute Fresnel amplitude reflection coefficients for s and p polarisations.

    Physical background
    -------------------
    At a planar interface between two dielectric media with real refractive
    indices n1 (incident side) and n2 (transmitted side), the Fresnel
    equations give the complex reflection amplitude coefficients r_s and r_p
    for s-polarisation (electric field perpendicular to the plane of incidence)
    and p-polarisation (electric field parallel to the plane of incidence):

        r_s = (n1 cos θ_i − n2 cos θ_t) / (n1 cos θ_i + n2 cos θ_t)
        r_p = (n2 cos θ_i − n1 cos θ_t) / (n2 cos θ_i + n1 cos θ_t)

    where θ_i and θ_t are the angles of incidence and transmission, related by
    Snell's law: n1 sin θ_i = n2 sin θ_t.

    For real refractive indices and incidence below the critical angle (no TIR),
    the amplitudes are real.  For TIR (sin²θ_t ≥ 1), both amplitudes have unit
    modulus and the function flags the condition rather than returning phases.

    Parameters
    ----------
    mu_i : float
        cos θ_i — cosine of the angle of incidence (must be positive; its
        absolute value is taken internally).
    n1 : float
        Refractive index of the incident medium (typically 1.0 for air).
    n2 : float
        Refractive index of the transmitted medium (typically ~1.33 for water).

    Returns
    -------
    rs_amp : float
        Fresnel amplitude reflection coefficient for s-polarisation.
    rp_amp : float
        Fresnel amplitude reflection coefficient for p-polarisation.
    mu_t : float
        cos θ_t — cosine of the transmission angle (0.0 under TIR).
    tir : bool
        True if total internal reflection occurs (n1 > n2 and θ_i > θ_c).
    """
    mu_i = float(np.clip(abs(mu_i), 0.0, 1.0))
    sin_i2 = max(0.0, 1.0 - mu_i * mu_i)             # sin²θ_i = 1 − cos²θ_i
    sin_t2 = (n1 / n2) ** 2 * sin_i2                  # Snell: sin²θ_t = (n1/n2)² sin²θ_i
    if sin_t2 >= 1.0:
        # Total internal reflection: no transmitted ray; return canonical TIR
        # values (rs_amp = -1, rp_amp = +1 by convention at grazing incidence).
        return -1.0, 1.0, 0.0, True
    mu_t = np.sqrt(max(0.0, 1.0 - sin_t2))            # cos θ_t
    rs_amp = (n1 * mu_i - n2 * mu_t) / (n1 * mu_i + n2 * mu_t + _EPS)
    rp_amp = (n2 * mu_i - n1 * mu_t) / (n2 * mu_i + n1 * mu_t + _EPS)
    return float(rs_amp), float(rp_amp), float(mu_t), False


def fresnel_reflection_mueller(mu_i: float, n1: float = 1.0, n2: float = 1.33) -> np.ndarray:
    """Mueller matrix for Fresnel reflection at a flat dielectric interface.

    Physical background
    -------------------
    For a planar interface the Fresnel Mueller matrix for reflection takes the
    block-diagonal form (in the plane-of-incidence reference frame):

        M_refl = | a  b  0  0 |
                 | b  a  0  0 |
                 | 0  0  c  0 |
                 | 0  0  0  c |

    where:
        R_s = r_s²,  R_p = r_p²         (reflectances for s and p polarisation)
        a = (R_s + R_p) / 2              (average reflectance → maps I to I)
        b = (R_p − R_s) / 2             (differential reflectance)
        c = sqrt(R_s · R_p)             (geometric mean → cross-polarisation term)

    Sign convention for b (critical for Q-sign consistency)
    --------------------------------------------------------
    Physical Fresnel optics at an air-water interface gives R_s > R_p (the
    s-polarised component is more strongly reflected), which means b < 0.
    This is intentional: the Mueller element M[1,0] = b generates Q < 0 from
    unpolarised incident light (I), meaning the reflected beam is preferentially
    s-polarised (perpendicular to the plane of incidence).

    PyMieDAP adopts the convention where Q < 0 means "preferentially
    perpendicular" (matching the Hovenier & van der Mee / DAP Rayleigh matrix
    convention where M[1,0] < 0 at 90° scattering).  Using b = (R_p − R_s)/2
    preserves this sign throughout the code.

    Parameters
    ----------
    mu_i : float
        cos θ_i — cosine of the angle of incidence on the facet normal.
        Must be non-negative (the magnitude is taken internally).
    n1 : float
        Refractive index of the incident medium (default: 1.0, air).
    n2 : float
        Refractive index of the transmitted medium (default: 1.33, water).

    Returns
    -------
    M : np.ndarray, shape (4, 4)
        4×4 Fresnel reflection Mueller matrix in the plane-of-incidence frame.
        The matrix must be rotated to/from meridian planes before use in a
        full scattering calculation (see ``_rotation_matrix``).

    Notes
    -----
    Under total internal reflection (n1 > n2, θ_i > θ_c), all incident light
    is reflected.  The matrix reduces to diag(1, 1, 1, −1) in the amplitude
    representation; here R_s = R_p = 1 so a=1, b=0, c=1.
    """
    rs_amp, rp_amp, _mu_t, _tir = _fresnel_amplitude_coefficients(mu_i, n1, n2)
    rs = rs_amp * rs_amp    # R_s = |r_s|² — reflectance for s-polarisation
    rp = rp_amp * rp_amp    # R_p = |r_p|² — reflectance for p-polarisation
    a = 0.5 * (rs + rp)     # average reflectance; maps I → a·I
    # b = (R_p − R_s)/2.  At air-water interface R_s > R_p, so b < 0.
    # This correctly produces Q < 0 (perpendicular polarisation) for reflected
    # unpolarised light, consistent with the PyMieDAP/DAP Rayleigh convention.
    b = 0.5 * (rp - rs)
    c = np.sqrt(max(rs * rp, 0.0))  # geometric mean; governs U and V coupling
    # The sign of U/V reflection depends on convention.  PyMieDAP's local
    # convention follows the historical DAP sign choice used here.
    return np.array(
        [
            [a, b, 0.0, 0.0],
            [b, a, 0.0, 0.0],
            [0.0, 0.0, c, 0.0],
            [0.0, 0.0, 0.0, c],
        ],
        dtype=float,
    )


def fresnel_transmission_mueller(mu_i: float, n1: float = 1.0, n2: float = 1.33) -> np.ndarray:
    """Power-normalised Mueller matrix for Fresnel transmission at a flat interface.

    Physical background
    -------------------
    By energy conservation at a lossless interface, the power transmittance for
    each polarisation component satisfies:

        T_s = 1 − R_s,   T_p = 1 − R_p

    Note: this is the *power* transmittance (fraction of incident irradiance
    transmitted), not the amplitude transmittance.  It already incorporates the
    refractive-index contrast and the obliquity (cos θ_t / cos θ_i) factor that
    arises when converting between Stokes parameters (which scale as irradiance)
    across an interface.

    The rough-interface transmission operator ``_transmission_local`` additionally
    applies the solid-angle Jacobian  n2² μ_t μ_i / (n2 μ_t − n1 μ_i)²  that
    accounts for the refraction of the solid-angle element when light crosses the
    interface (Zhai et al. 2010, Trees & Stam 2019, Appendix A).

    The Mueller matrix has the same block structure as the reflection matrix:

        M_trans = | a  b  0  0 |   with  a = (T_s + T_p)/2
                  | b  a  0  0 |         b = (T_p − T_s)/2
                  | 0  0  c  0 |         c = sqrt(T_s · T_p)
                  | 0  0  0  c |

    Sign convention for b
    ----------------------
    At an air-water interface T_s < T_p (since R_s > R_p, less s-polarised light
    is transmitted), so b = (T_p − T_s)/2 > 0.  This means the transmitted beam
    acquires Q > 0 (parallel polarisation); the reflected beam acquires Q < 0.
    Both signs are consistent with the PyMieDAP/DAP Rayleigh convention.

    Parameters
    ----------
    mu_i : float
        cos θ_i — cosine of the angle of incidence.
    n1 : float
        Refractive index of the incident medium (default: 1.0, air).
    n2 : float
        Refractive index of the transmitted medium (default: 1.33, water).

    Returns
    -------
    M : np.ndarray, shape (4, 4)
        4×4 power-normalised Fresnel transmission Mueller matrix in the plane-
        of-incidence frame.  Returns a zero matrix under TIR.
    """
    rs_amp, rp_amp, _mu_t, tir = _fresnel_amplitude_coefficients(mu_i, n1, n2)
    if tir:
        # Total internal reflection: no transmitted power in any polarisation.
        return np.zeros((4, 4), dtype=float)
    rs = rs_amp * rs_amp    # R_s
    rp = rp_amp * rp_amp    # R_p
    ts = max(0.0, 1.0 - rs) # T_s = 1 − R_s (energy conservation)
    tp = max(0.0, 1.0 - rp) # T_p = 1 − R_p
    a = 0.5 * (ts + tp)     # average transmittance
    # b = (T_p − T_s)/2.  T_p > T_s at air-water interface → b > 0,
    # giving Q > 0 in the transmitted beam (parallel polarisation preferred).
    # Same (R_p−R_s)/2 sign convention used in fresnel_reflection_mueller.
    b = 0.5 * (tp - ts)
    c = np.sqrt(max(ts * tp, 0.0))  # geometric mean; governs U and V coupling
    return np.array(
        [
            [a, b, 0.0, 0.0],
            [b, a, 0.0, 0.0],
            [0.0, 0.0, c, 0.0],
            [0.0, 0.0, 0.0, c],
        ],
        dtype=float,
    )


def fresnel_unpolarized_reflectance(mu_i: float, n1: float = 1.0, n2: float = 1.33) -> float:
    """Scalar Fresnel reflectance for unpolarised light (average of R_s and R_p).

    This is the [0, 0] element of the Fresnel reflection Mueller matrix, i.e.
    the fraction of incident intensity that is reflected regardless of the
    polarisation state.  Useful for quick energy-budget estimates.

    Parameters
    ----------
    mu_i : float
        cos θ_i — cosine of the angle of incidence.
    n1, n2 : float
        Refractive indices of the incident and transmitted media.

    Returns
    -------
    float
        Unpolarised reflectance in [0, 1].
    """
    return float(fresnel_reflection_mueller(mu_i, n1=n1, n2=n2)[0, 0])


def cox_munk_slope_variance(wind_speed: float) -> float:
    """Cox & Munk (1954) isotropic mean-square wave-slope variance σ².

    The empirical linear relationship between wind speed and slope variance was
    established by Cox & Munk (1954, J. Mar. Res. 13, 198) from sun-glitter
    photographs taken from aircraft:

        σ² = 0.003 + 0.00512 · U_{10}

    where U_{10} is the wind speed at 10 m height in m/s.  Here the isotropic
    (azimuth-averaged) form is used; the original paper also provides an
    anisotropic decomposition into upwind/crosswind components which is not
    implemented here.

    Parameters
    ----------
    wind_speed : float
        Wind speed U_{10} in m/s at 10 m height above the ocean surface.

    Returns
    -------
    float
        Isotropic slope variance σ² (dimensionless).  The RMS slope is σ = √σ².
    """
    return 0.003 + 0.00512 * float(wind_speed)


def monahan_whitecap_fraction(wind_speed: float) -> float:
    """Monahan & O'Muircheartaigh (1980) whitecap (foam) fractional coverage.

    Whitecaps are wind-generated patches of breaking waves covered by air
    bubbles and foam.  They scatter light roughly Lambertian-like with a high
    broadband albedo (~0.22).  The empirical power law from Monahan &
    O'Muircheartaigh (1980, J. Phys. Oceanogr. 10, 2094) relates coverage to
    wind speed:

        f_w = 2.95 × 10⁻⁶ · U_{10}^{3.52}

    At typical ocean conditions (U ≈ 7 m/s) this gives f_w ≈ 0.5%.  The
    fraction grows rapidly at high wind speeds (e.g. ~3% at 14 m/s).

    Parameters
    ----------
    wind_speed : float
        Wind speed U_{10} in m/s.

    Returns
    -------
    float
        Whitecap fractional coverage f_w in [0, 1] (dimensionless).
    """
    return 2.95e-6 * float(wind_speed) ** 3.52


def cox_munk_pdf(mu_n: float, sigma: float) -> float:
    """Cox & Munk (1954) probability density for the facet-normal zenith angle.

    The isotropic Cox-Munk model assumes that ocean-surface microfacet slopes
    follow a 2-D isotropic Gaussian distribution with variance σ².  For a
    facet whose outward normal has zenith cosine μ_n = cos θ_n (θ_n is the
    tilt angle of the facet away from the horizontal), the probability density
    function for μ_n is:

        P(μ_n) = exp(−tan²θ_n / σ²) / (π σ² μ_n³)

    where tan²θ_n = (1 − μ_n²) / μ_n².

    This PDF describes the fraction of ocean surface area occupied by facets
    tilted by θ_n.  It peaks at μ_n = 1 (flat facet) and falls off for large
    tilts.  The scale factor 1/(π σ² μ_n³) normalises the distribution over
    the hemisphere of facet normals.

    Parameters
    ----------
    mu_n : float
        cos θ_n — cosine of the facet-normal zenith angle.  Must be in (0, 1].
        Returns 0 for μ_n ≤ 0 (facet pointing downward — unphysical).
    sigma : float
        RMS slope σ = √σ² (square root of the Cox-Munk slope variance).
        Must be positive.

    Returns
    -------
    float
        Probability density P(μ_n) in units of sr⁻¹ (per steradian).
    """
    if mu_n <= 0.0 or sigma <= 0.0:
        return 0.0
    mu_n = max(float(mu_n), _EPS)
    # tan²θ_n = (1 − μ_n²) / μ_n² — squared tangent of the facet tilt angle
    return float(
        np.exp(-((1.0 - mu_n * mu_n) / (sigma * sigma * mu_n * mu_n)))
        / (np.pi * sigma * sigma * mu_n**3)
    )


def smith_sancer_shadowing(mu: float, mu0: float, sigma: float) -> float:
    """Smith (1967) / Sancer (1969) bi-directional shadowing-masking factor.

    Physical meaning
    ----------------
    On a rough surface, some facets are shadowed (not illuminated by the
    incoming ray) or masked (hidden from the detector).  The shadowing factor
    S(μ, μ₀, σ) is the fraction of facets that are simultaneously illuminated
    and visible, i.e. not shadowed by neighbouring surface waves from either
    direction.

    The Smith (1967)/Sancer (1969) approximation treats the two shadowing
    contributions as statistically independent:

        S = 1 / (1 + Λ(μ) + Λ(μ₀))

    where the Lambda function is:

        Λ(γ) = 0.5 [ (σ/γ) √((1−γ²)/π) exp(−γ²/(σ²(1−γ²))) − erfc(γ/(σ√(1−γ²))) ]

    with γ = μ (the cosine of the observation or illumination zenith angle).

    For a flat surface (σ → 0) or near-nadir geometry (μ, μ₀ → 1), S → 1
    (no shadowing).  At grazing angles (μ or μ₀ → 0), S → 0 (complete
    shadowing).

    Parameters
    ----------
    mu : float
        cos θ_r — cosine of the reflected/observation zenith angle.
    mu0 : float
        cos θ_i — cosine of the incidence zenith angle.
    sigma : float
        RMS slope σ = √σ² (Cox-Munk slope standard deviation).

    Returns
    -------
    float
        Shadowing factor S in [0, 1].  Multiplied into the reflection/
        transmission kernel to suppress contributions at grazing angles.
    """
    def lam(gamma: float) -> float:
        """Single-direction Lambda function for the Smith/Sancer model."""
        gamma = float(np.clip(gamma, _EPS, 1.0 - _EPS))
        root = np.sqrt(max(1.0 - gamma * gamma, 0.0))  # sin θ
        if root < _EPS:
            return 0.0
        arg = gamma / (sigma * root)  # γ / (σ sin θ) = cot θ / σ
        import math

        return 0.5 * (
            (sigma / gamma) * np.sqrt((1.0 - gamma * gamma) / np.pi) * np.exp(-(arg * arg))
            - math.erfc(arg)
        )

    return 1.0 / (1.0 + lam(mu) + lam(mu0))


# ---------------------------------------------------------------------------
# Pure-water optical coefficients
# ---------------------------------------------------------------------------

# Compact pure-water optical coefficient lookup table.
# Absorption values (aw) are based on Pope & Fry (1997, Appl. Opt. 36, 8710)
# and Smith & Baker (1981, Appl. Opt. 20, 177) at visible/NIR wavelengths.
# Scattering values (bw) follow the λ⁻⁴·³² Rayleigh-like fit of
# Morel (1974) / Shifrin (1988).  Users who need exact wavelength-specific
# values should pass aw and bw directly to WaterBody or OceanSurface.
_WATER_WAVELENGTH_NM = np.array(
    [350, 380, 400, 425, 443, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 800, 865],
    dtype=float,
)
_WATER_ABS_M1 = np.array(
    [0.010, 0.0046, 0.0066, 0.011, 0.0145, 0.017, 0.033, 0.045, 0.060, 0.090, 0.220, 0.275, 0.350, 0.430, 0.650, 1.15, 2.60, 4.60, 12.0],
    dtype=float,
)
_WATER_SCAT_M1 = np.array(
    [0.0100, 0.0075, 0.0063, 0.0049, 0.0040, 0.0031, 0.0025, 0.0020, 0.0017, 0.0014, 0.0012, 0.0010, 0.0009, 0.0008, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003],
    dtype=float,
)


def pure_water_coefficients(wavelength_um: float) -> Tuple[float, float]:
    """Return approximate pure-water absorption and scattering coefficients.

    Interpolates the built-in lookup table (Pope & Fry / Morel) at the
    requested wavelength.  For wavelengths outside the table range
    [350, 865 nm], numpy's ``np.interp`` clamps to the nearest endpoint value.

    Parameters
    ----------
    wavelength_um : float
        Wavelength in microns (e.g. 0.550 for green light at 550 nm).

    Returns
    -------
    aw : float
        Absorption coefficient a_w in m⁻¹.  Rises steeply in the NIR (>600 nm)
        and has a strong minimum around 420 nm (≈ 0.004 m⁻¹).
    bw : float
        Scattering coefficient b_w in m⁻¹.  Decreases monotonically with
        wavelength following an approximate λ⁻⁴·³² power law; at 500 nm
        b_w ≈ 0.0025 m⁻¹.
    """
    wl_nm = 1000.0 * float(wavelength_um)   # convert μm → nm for table lookup
    aw = float(np.interp(wl_nm, _WATER_WAVELENGTH_NM, _WATER_ABS_M1))
    bw = float(np.interp(wl_nm, _WATER_WAVELENGTH_NM, _WATER_SCAT_M1))
    return aw, bw


def pure_water_diffuse_albedo(
    wavelength_um: float,
    depth_m: float = 100.0,
    bottom_albedo: float = 0.0,
    aw: float | None = None,
    bw: float | None = None,
) -> float:
    """Legacy scalar diffuse water-column albedo (unpolarised, hemispheric closure).

    This function is a fast, closed-form approximation retained for backward
    compatibility and quick diagnostics.  It is used by the ``diffuse_closure``
    solver path.  The proper polarized adding-doubling solver (``WaterBody``)
    should be used for publication-quality calculations.

    The water-leaving component uses the Gordons (1973) / Morel & Prieur (1977)
    single-scattering albedo approximation:

        A_water = 0.20 · ω₀ · (1 − exp(−τ)) / (1 + √(3(1−ω₀)))

    where ω₀ = b_w / (a_w + b_w) is the single-scattering albedo and
    τ = (a_w + b_w) · depth is the total optical depth.  The factor 0.20 is an
    empirical diffuse-geometry constant.

    The bottom contribution is attenuated by round-trip propagation:

        A_bottom = A_bot · exp(−2τ)

    Parameters
    ----------
    wavelength_um : float
        Wavelength in microns.
    depth_m : float
        Geometric depth of the water column in metres.
    bottom_albedo : float
        Lambertian bottom albedo (0 = black, 1 = perfectly white sand).
    aw, bw : float | None
        Optional override for absorption and scattering coefficients in m⁻¹.
        If None, the built-in pure-water table is used.

    Returns
    -------
    float
        Effective diffuse water-column albedo in [0, 1].
    """
    if aw is None or bw is None:
        aw0, bw0 = pure_water_coefficients(wavelength_um)
        aw = aw0 if aw is None else aw
        bw = bw0 if bw is None else bw

    c = max(float(aw) + float(bw), _EPS)            # total attenuation coefficient c = a + b
    omega = np.clip(float(bw) / c, 0.0, 1.0)        # single-scattering albedo ω₀ = b/c
    tau = c * max(float(depth_m), 0.0)               # total vertical optical depth τ = c·z
    # Water-leaving irradiance fraction (approximate diffuse closure)
    water_leaving = 0.20 * omega * (1.0 - np.exp(-tau)) / (1.0 + np.sqrt(max(0.0, 3.0 * (1.0 - omega))))
    # Bottom contribution attenuated by two-way propagation through the water column
    bottom = float(bottom_albedo) * np.exp(-2.0 * tau)
    return float(np.clip(water_leaving + bottom, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Fourier projection and Rayleigh scattering
# ---------------------------------------------------------------------------


def _quadrature(n: int, include_extra_mu_one: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gauss-Legendre quadrature nodes and weights on the interval (0, 1].

    The Doubling-Adding-Polarization (DAP) solver performs hemisphere-integrated
    radiance calculations using Gauss-Legendre quadrature over the polar cosine
    μ = cos θ ∈ (0, 1].  The standard Gauss-Legendre nodes lie on [−1, 1]; this
    function rescales them to [0, 1] to span the upward (or downward) hemisphere:

        x ∈ [−1, 1]  →  μ = (x + 1) / 2 ∈ [0, 1]
        w_{GL}       →  w = w_{GL} / 2

    The "smf" (scaling/metric factor) pre-multiplies μ-space Stokes vectors to
    account for the solid-angle element dΩ = 2π μ dμ that appears in the
    radiative transfer equation integral operators:

        smf_i = √(2 μ_i w_i)

    Supermatrices formed from outer products ``smf_i * smf_j`` embed this
    quadrature weight into the matrix elements, so the adding-doubling
    operations on supermatrices are equivalent to quadrature-weighted integrals.

    Parameters
    ----------
    n : int
        Number of Gauss-Legendre streams per hemisphere.
    include_extra_mu_one : bool
        When True, append an extra stream at μ = 1 (nadir/zenith direction).
        This is needed in DAP for the direct-solar-beam direction.
        The extra weight is set to 0.5 (following DAP's setmu.f convention).

    Returns
    -------
    mus : np.ndarray, shape (n,) or (n+1,)
        Quadrature nodes μ_i ∈ (0, 1].
    weights : np.ndarray, shape (n,) or (n+1,)
        Quadrature weights w_i (positive, sum to ~0.5).
    smf : np.ndarray, shape (n,) or (n+1,)
        Scaling factors smf_i = √(2 μ_i w_i) for supermatrix construction.
    """
    if n < 1:
        raise ValueError("number of streams must be positive")
    x, w = np.polynomial.legendre.leggauss(n)     # nodes/weights on [−1, 1]
    mus = 0.5 * (x + 1.0)                          # rescale to (0, 1]
    weights = 0.5 * w                              # rescale weights accordingly
    if include_extra_mu_one:
        mus = np.concatenate([mus, np.array([1.0])])
        # DAP setmu.f uses w=0.5 for the extra vertical direction.
        weights = np.concatenate([weights, np.array([0.5])])
    smf = np.sqrt(2.0 * mus * weights)             # supermatrix metric factor
    return mus, weights, smf


def _use_sine_block(row: int, col: int) -> bool:
    """Determine whether a Fourier block uses sine (True) or cosine (False) expansion.

    In the azimuthal Fourier decomposition of the reflection/transmission
    operators, different Stokes combinations transform differently under
    azimuth reversal φ → −φ:

    - I, Q (row/col 0, 1) are even in φ → cosine expansion
    - U, V (row/col 2, 3) have mixed parity in φ → sine expansion for
      cross-terms (I↔U, I↔V, Q↔U, Q↔V), cosine for self-terms

    The parity rule is: use sine if exactly one of {row, col} is in {2, 3}
    (the XOR of the two indicator functions).  This is the same convention
    used in the original DAP code and in the earlier ocean branch.

    Parameters
    ----------
    row : int
        Row Stokes index (0=I, 1=Q, 2=U, 3=V).
    col : int
        Column Stokes index.

    Returns
    -------
    bool
        True if the (row, col) block expands in sine terms; False for cosine.
    """
    # Same parity rule as the earlier ocean branch.  It matches the common DAP
    # convention where U/V terms are sine-like under azimuth reversal.
    return ((row in (2, 3)) ^ (col in (2, 3)))


def _signed_angle_vec(
    old_normals: np.ndarray,
    new_normals: np.ndarray,
    ray_directions: np.ndarray,
) -> np.ndarray:
    """Vectorised signed Stokes-reference rotation angle for an array of geometries.

    This is the batch (all-phi-at-once) counterpart of
    ``_signed_angle_between_plane_normals``.  It computes the signed angle ξ_i
    between ``old_normals[i]`` and ``new_normals[i]`` around
    ``ray_directions[i]`` for all i simultaneously, enabling the vectorised
    reflection/transmission paths (_reflection_vec, _transmission_vec) to
    avoid a Python loop over azimuth angles.

    The calculation is identical to the scalar version:
        cos ξ_i = old_n_i · new_n_i      (dot product after normalisation)
        sin ξ_i = ray_i · (old_n_i × new_n_i)   (triple product for sign)
        ξ_i     = arctan2(sin ξ_i, cos ξ_i)

    Parameters
    ----------
    old_normals : np.ndarray, shape (n, 3)
        Array of old reference-plane normals (one per azimuth sample).
    new_normals : np.ndarray, shape (n, 3)
        Array of new reference-plane normals.
    ray_directions : np.ndarray, shape (n, 3)
        Array of ray-propagation direction vectors (rotation axes).

    Returns
    -------
    np.ndarray, shape (n,)
        Signed rotation angles ξ_i in radians, each in (−π, π].
    """
    mag_o = np.linalg.norm(old_normals, axis=1, keepdims=True)
    mag_n = np.linalg.norm(new_normals, axis=1, keepdims=True)
    # Safe normalisation: set to zero vector if the magnitude is negligible
    on = np.where(mag_o > _EPS, old_normals / (mag_o + _EPS), 0.0)
    nn = np.where(mag_n > _EPS, new_normals / (mag_n + _EPS), 0.0)
    rd = ray_directions / (np.linalg.norm(ray_directions, axis=1, keepdims=True) + _EPS)
    # cos ξ: dot product of normalised plane normals, clipped to [-1, 1]
    x = np.clip(np.einsum("ni,ni->n", on, nn), -1.0, 1.0)
    # sin ξ: signed component of the cross product projected onto the ray axis
    cross = np.cross(on, nn)   # (n, 3)
    y = np.einsum("ni,ni->n", rd, cross)
    return np.arctan2(y, x)


def _project_operator_fourier(
    dst_mus: np.ndarray,
    dst_smf: np.ndarray,
    src_mus: np.ndarray,
    src_smf: np.ndarray,
    nmat: int,
    n_fourier: int,
    n_phi: int,
    local_matrix: Callable[[float, float, float], np.ndarray],
    local_matrix_vec=None,
) -> np.ndarray:
    """Project an azimuth-dependent Mueller-kernel operator into DAP supermatrices.

    What this function computes
    ---------------------------
    In the Doubling-Adding-Polarization (DAP) method, reflection and transmission
    operators are represented as *supermatrices* whose (i·nmat : (i+1)·nmat,
    j·nmat : (j+1)·nmat) block is the nmat×nmat Mueller matrix coupling the j-th
    source stream (at μ_src[j]) to the i-th destination stream (at μ_dst[i]).

    Because the ocean surface is azimuthally symmetric about the vertical (z)
    axis, the operators can be expanded in a Fourier cosine/sine series in the
    relative azimuth Δφ = φ_dst − φ_src:

        M(μ_dst, μ_src, Δφ) = Σ_{m=0}^{N} M^(m)(μ_dst, μ_src) cos(m Δφ)
                                                       [or sin for U/V blocks]

    This function evaluates M at n_phi equidistant azimuth samples, computes
    the Fourier coefficients M^(m) by numerical integration (trapezoidal rule
    with equal spacing = discrete Fourier transform), and assembles the
    result into the output supermatrix array of shape (nd·nmat, ns·nmat, N+1).

    The supermatrix elements include the quadrature weighting:

        out[i·nmat:, j·nmat:, m] = smf_dst[i] · M^(m)(μ_dst[i], μ_src[j]) · smf_src[j]

    so that contracting with another smf-scaled supermatrix implements a
    hemisphere-integrated polarised radiative transfer operator.

    Vectorised vs scalar evaluation path
    --------------------------------------
    Two evaluation paths exist for the local Mueller matrix M(μ_dst, μ_src, φ):

    **Scalar path** (``local_matrix_vec is None``): calls ``local_matrix`` once
    per (μ_dst, μ_src, φ) triple in a Python loop over phi.  Correct but slow
    for large n_phi.

    **Vectorised path** (``local_matrix_vec`` provided): calls
    ``local_matrix_vec(μ_dst, μ_src, phis_array)`` once, returning the full
    (n_phi, 4, 4) array in a single NumPy call.  ~50× faster.  Used by the
    ``_reflection_vec`` / ``_transmission_vec`` methods of ``RoughInterface``.

    Einsum details
    --------------
    After collecting all azimuth samples into ``mats`` of shape
    (n_phi, nmat, nmat), the Fourier coefficients are computed as:

        mats_T  = mats.transpose(1, 2, 0)       # → (nmat, nmat, n_phi)
        coeffs_cos[m, r, c] = Σ_k mats_T[r, c, k] · cos(m φ_k) / n_phi
        coeffs_sin[m, r, c] = Σ_k mats_T[r, c, k] · sin(m φ_k) / n_phi

    using ``np.einsum("rck,mk->mrc", mats_T, trig_terms)``.  The axis labels
    are: r = Stokes row, c = Stokes column, k = phi-sample index,
    m = Fourier order.  The output has shape (n_fourier+1, nmat, nmat).

    Finally, a boolean mask ``use_sine_mask[r, c]`` selects between the cosine
    and sine coefficient arrays per Stokes block (U/V cross-terms use sine,
    all others use cosine).

    Parameters
    ----------
    dst_mus : np.ndarray, shape (nd,)
        Cosines of the destination (outgoing) polar angles.
    dst_smf : np.ndarray, shape (nd,)
        Supermatrix scaling factors for the destination streams.
    src_mus : np.ndarray, shape (ns,)
        Cosines of the source (incoming) polar angles.
    src_smf : np.ndarray, shape (ns,)
        Supermatrix scaling factors for the source streams.
    nmat : int
        Number of Stokes components (1, 3, or 4).
    n_fourier : int
        Highest Fourier order N retained.
    n_phi : int
        Number of equidistant azimuth quadrature samples.  Should satisfy
        n_phi ≥ max(32, 2·n_fourier + 1) to avoid aliasing.
    local_matrix : Callable[[float, float, float], np.ndarray]
        Function that returns the local 4×4 Mueller matrix for scalar
        (μ_dst, μ_src, phi) input.  Used when ``local_matrix_vec`` is None.
    local_matrix_vec : Callable or None
        Optional vectorised version of ``local_matrix`` that accepts
        (μ_dst, μ_src, phis_array) and returns (n_phi, 4, 4).

    Returns
    -------
    out : np.ndarray, shape (nd·nmat, ns·nmat, n_fourier+1), Fortran order
        Fourier-decomposed, smf-scaled supermatrix.  The third axis indexes
        Fourier order m = 0, 1, …, n_fourier.
    """
    nd = len(dst_mus)
    ns = len(src_mus)
    nsup_d = nd * nmat
    nsup_s = ns * nmat
    out  = np.zeros((nsup_d, nsup_s, n_fourier + 1), dtype=np.float64, order="F")
    # Equidistant azimuth samples on [0, 2π); endpoint=False avoids double-counting φ=0=2π
    phis = np.linspace(0.0, 2.0 * np.pi, int(n_phi), endpoint=False)
    # Pre-compute trig basis for all Fourier orders: shape (n_fourier+1, n_phi)
    cos_terms = np.array([np.cos(m * phis) for m in range(n_fourier + 1)])
    sin_terms = np.array([np.sin(m * phis) for m in range(n_fourier + 1)])

    # Sine/cosine selector per (row, col): True → sin, False → cos.
    # This 2-D boolean mask has shape (nmat, nmat).
    use_sine_mask = np.array(
        [[_use_sine_block(r, c) for c in range(nmat)] for r in range(nmat)],
        dtype=bool,
    )  # (nmat, nmat)

    for i, mu_dst in enumerate(dst_mus):
        for j, mu_src in enumerate(src_mus):
            # --- Collect Mueller matrices at all azimuth angles ---
            if local_matrix_vec is not None:
                # Vectorised path: one call returns (n_phi, 4, 4)
                raw  = local_matrix_vec(float(mu_dst), float(mu_src), phis)
                mats = raw[:, :nmat, :nmat]           # (n_phi, nmat, nmat)
            else:
                # Scalar path: n_phi separate calls (slow but always available)
                mats = np.empty((len(phis), nmat, nmat), dtype=float)
                for ip, phi in enumerate(phis):
                    mats[ip] = _stokes_block(
                        local_matrix(float(mu_dst), float(mu_src), float(phi)), nmat
                    )

            # ── Vectorised Fourier projection ──────────────────────────────
            # Instead of 3 nested Python loops (n_fourier × nmat × nmat),
            # use einsum to compute all Fourier coefficients at once:
            #   coeffs_cos[m, r, c] = (1/n_phi) Σ_k mats[k,r,c] · cos(m·φ_k)
            #   coeffs_sin[m, r, c] = (1/n_phi) Σ_k mats[k,r,c] · sin(m·φ_k)
            # Axis key: r=Stokes row, c=Stokes col, k=phi-sample, m=Fourier order.
            mats_T = np.ascontiguousarray(mats.transpose(1, 2, 0))   # → (nmat, nmat, n_phi)
            coeffs_cos = np.einsum("rck,mk->mrc", mats_T, cos_terms) / n_phi
            coeffs_sin = np.einsum("rck,mk->mrc", mats_T, sin_terms) / n_phi

            # Select sine or cosine coefficient per Stokes block using the
            # (n_fourier+1, nmat, nmat) boolean broadcast of use_sine_mask.
            coeffs = np.where(use_sine_mask[np.newaxis], coeffs_sin, coeffs_cos)

            # Write into the output supermatrix, scaled by smf factors.
            # coeffs has shape (n_fourier+1, nmat, nmat); transpose to
            # (nmat, nmat, n_fourier+1) to match the output layout.
            out[i * nmat:(i + 1) * nmat,
                j * nmat:(j + 1) * nmat, :] = (
                dst_smf[i] * coeffs.transpose(1, 2, 0) * src_smf[j]
            )
    return out


def rayleigh_mueller(cos_scatter: float, depolarization: float = 0.0) -> np.ndarray:
    """Depolarized Rayleigh scattering Mueller matrix in the scattering plane.

    Physical background
    -------------------
    Pure seawater scattering arises from random density fluctuations (Einstein-
    Smoluchowski theory).  The angular scattering pattern is very close to that
    of classical Rayleigh scattering but is slightly depolarized due to the
    anisotropy of the water molecule's polarisability tensor.

    The standard Rayleigh phase matrix (in the scattering plane frame) is:

        M_Rayleigh(θ) = (3/4) ·
            | 1+cos²θ   −sin²θ     0        0       |
            | −sin²θ    1+cos²θ    0        0       |
            |    0          0    2cos θ      0       |
            |    0          0      0      2cos θ     |

    (normalised so that ∫ P11 dΩ = 4π).

    Depolarization correction (King factor)
    ----------------------------------------
    With a King-factor depolarization δ ∈ [0, 0.5], the polarizing off-diagonal
    elements are reduced by a factor:

        pol = (1 − δ) / (1 + δ/2)

    This is the standard depolarization correction from King (1923) /
    Chandrasekhar (1960).  For pure seawater δ ≈ 0.09 (Morel 1974).  Setting
    δ = 0 recovers the ideal Rayleigh matrix.

    Parameters
    ----------
    cos_scatter : float
        Cosine of the scattering angle θ (−1 to 1).  cos θ = k̂_in · k̂_out.
    depolarization : float
        King depolarization factor δ in [0, 0.5].  Default 0 = no depolarization.

    Returns
    -------
    M : np.ndarray, shape (4, 4)
        Rayleigh scattering Mueller matrix expressed in the scattering plane
        reference frame.  Must be sandwiched between rotation matrices before
        use in a full meridian-plane calculation.

    Notes
    -----
    P12 < 0 for all scattering angles (the scattered light tends to have Q < 0
    = perpendicular polarisation), consistent with the Rayleigh sign convention
    used throughout PyMieDAP.
    """
    c = float(np.clip(cos_scatter, -1.0, 1.0))
    s2 = max(0.0, 1.0 - c * c)   # sin²θ
    # King-like depolarization reduction factor for off-diagonal polarising elements.
    # For delta=0 this reduces to the conventional Rayleigh matrix with integral P11 = 4π.
    delta = float(np.clip(depolarization, 0.0, 0.5))
    pol = (1.0 - delta) / (1.0 + 0.5 * delta)   # reduction factor ∈ [0, 1]
    p11 = 0.75 * (1.0 + c * c)    # (3/4)(1 + cos²θ) — governs I→I
    p12 = -0.75 * pol * s2        # −(3/4)·pol·sin²θ — always ≤ 0 → Q < 0 for unpolarised I
    p22 = p11                     # for Rayleigh, P22 = P11 exactly
    p33 = 1.5 * pol * c           # (3/2)·pol·cos θ — governs U→U and V→V coupling
    p44 = p33                     # P44 = P33 for Rayleigh
    return np.array(
        [
            [p11, p12, 0.0, 0.0],
            [p12, p22, 0.0, 0.0],
            [0.0, 0.0, p33, 0.0],
            [0.0, 0.0, 0.0, p44],
        ],
        dtype=float,
    )


def _scattering_matrix_local(k_in: np.ndarray, k_out: np.ndarray, depolarization: float) -> np.ndarray:
    """Rayleigh Mueller matrix rotated to the meridian-plane reference frame.

    This is the full single-scattering phase matrix for one photon direction
    change from ``k_in`` to ``k_out`` inside the water body, expressed in the
    meridian-plane (local Stokes) reference frame of each ray.

    The procedure is:
    1. Find the scattering angle θ from k̂_in · k̂_out = cos θ.
    2. Identify the scattering plane normal n_scat = (k_in × k_out) / |...|.
    3. Compute ξ₁: signed angle from the incoming meridian plane to the
       scattering plane, around the incoming ray direction k_in.
    4. Compute ξ₂: signed angle from the scattering plane to the outgoing
       meridian plane, around the outgoing ray direction k_out.
    5. Return  L(−ξ₂) · M_Rayleigh(θ) · L(−ξ₁).

    This sandwich of rotation matrices ensures the Stokes vector undergoes the
    physically correct reference-plane transformation at the scattering event.

    Parameters
    ----------
    k_in : np.ndarray, shape (3,)
        Incoming ray direction unit vector (need not be normalised; normalised
        internally).
    k_out : np.ndarray, shape (3,)
        Outgoing (scattered) ray direction unit vector.
    depolarization : float
        King depolarization factor δ for the Rayleigh matrix.

    Returns
    -------
    M : np.ndarray, shape (4, 4)
        Full single-scattering Mueller matrix in the meridian-plane frame.
    """
    k_in = _safe_normalize(k_in)
    k_out = _safe_normalize(k_out)
    # Cosine of the scattering angle
    cos_theta = float(np.clip(np.dot(k_in, k_out), -1.0, 1.0))
    # Normal to the scattering plane (perpendicular to both incoming and outgoing rays)
    scatter_plane = _safe_normalize(np.cross(k_in, k_out))
    if np.linalg.norm(scatter_plane) < _EPS:
        # Degenerate case: forward or backward scattering (k_in ∥ k_out).
        # Use the outgoing meridian plane as the scattering plane.
        scatter_plane = _meridian_plane_normal(k_out)
    mer_in  = _meridian_plane_normal(k_in)
    mer_out = _meridian_plane_normal(k_out)
    # ξ₁ = angle from incoming meridian → scattering plane (measured around k_in)
    xi1 = _signed_angle_between_plane_normals(mer_in, scatter_plane, k_in)
    # ξ₂ = angle from scattering plane → outgoing meridian (measured around k_out)
    xi2 = _signed_angle_between_plane_normals(scatter_plane, mer_out, k_out)
    # Full Mueller matrix: rotate incoming, apply Rayleigh, rotate outgoing
    return _rotation_matrix(-xi2) @ rayleigh_mueller(cos_theta, depolarization) @ _rotation_matrix(-xi1)


# ---------------------------------------------------------------------------
# Water-body adding-doubling solver
# ---------------------------------------------------------------------------


@dataclass
class WaterBody:
    """Polarized subinterface water-column solver using adding-doubling.

    Physical model
    --------------
    The ocean water column below the rough air-water interface is modelled as a
    *plane-parallel, homogeneous* slab of pure seawater.  The slab has:

    - **Absorption** a_w and **scattering** b_w coefficients from the built-in
      pure-water table (Pope & Fry / Morel) or user-supplied values.
    - **Total attenuation** c = a_w + b_w and **optical depth** τ = c · depth.
    - **Single-scattering albedo** ω₀ = b_w / c (fraction of attenuated photons
      that are scattered rather than absorbed).
    - A **depolarized Rayleigh phase matrix** to describe the angular shape of
      scattering.  Pure seawater scatters almost like an anisotropic Rayleigh
      scatterer (depolarization factor δ ≈ 0.09).

    The lower boundary is either a black bottom (bottom_albedo=0, default) or a
    Lambertian reflector.

    Adding-doubling algorithm
    -------------------------
    The water body is solved for its polarized *reflection supermatrix* R_W
    using the **adding-doubling method** (van de Hulst 1963, de Haan et al. 1987):

    1. **Initialise a thin slab**: start with a layer of optical depth
       τ₀ = τ / 2^N where N = ceil(log₂(τ / τ_init)) is chosen so that
       τ₀ ≤ ``initial_tau`` (typically 0.01).  For such a thin slab the
       single-scattering approximation is very accurate.

       For a thin layer of optical depth τ₀, the m-th Fourier order of the
       reflection and transmission supermatrices are:

           R^(m)_{ij} = (ω₀/4) · h_R(i,j) · P_ref^(m)(μ_i, μ_j)
           T^(m)_{ij} = (ω₀/4) · h_T(i,j) · P_tra^(m)(μ_i, μ_j)  +  δ_{ij} exp(−τ₀/μ_i)

       where h_R(i,j) = (1 − e^{−τ₀/μ_i} · e^{−τ₀/μ_j}) / (μ_i + μ_j)
       and   h_T(i,j) = (e^{−τ₀/μ_j} − e^{−τ₀/μ_i}) / (μ_i − μ_j)
       (L'Hôpital limit when μ_i = μ_j).

       The diagonal term in T is the **direct-beam transmittance** exp(−τ₀/μ_i).

    2. **Double N times**: apply the doubling formula to combine two identical
       layers, each of optical depth τ/2^k, into one layer of depth τ/2^{k−1}:

           R₂ = R + T · (I − R²)^{−1} · R · T
           T₂ = T · (I − R²)^{−1} · T

       After N doublings the slab has optical depth τ.

    3. **Add the bottom boundary** (if non-black): combine the water slab with
       a Lambertian lower boundary using the adding equations.

    The result is the upward-directed reflection supermatrix R_W mapping
    downwelling water streams to upwelling water streams, needed by the
    Trees & Stam coupling formula.

    Rectangular supermatrices
    -------------------------
    The water body uses its own, user-configurable Gaussian quadrature
    (``n_streams``) which can be denser than the atmosphere's DAP streams.
    This "rectangular supermatrix" approach (separate angular grids for the
    interface and the water body) is described in Trees & Stam (2019),
    Appendix A.  The interface transmission operators T_I and T_I* couple the
    two grids.

    Parameters
    ----------
    wavelength_um : float
        Wavelength in microns.  Used to look up a_w and b_w if not provided.
    depth_m : float
        Geometric water depth in metres.  Trees & Stam (2019) use 100 m.
    bottom_albedo : float
        Lambertian lower-boundary albedo in [0, 1].  Default 0 = black bottom.
    depolarization : float
        King depolarization factor δ for the Rayleigh phase matrix.
        Default 0.09 matches Morel (1974) for pure seawater.
    n_streams : int
        Number of Gaussian-Legendre streams per hemisphere inside the water.
        More streams → more accurate but slower.  16 is adequate for most uses.
    initial_tau : float
        Maximum optical thickness of the starting thin-slab approximation.
        Smaller values improve accuracy at the cost of more doubling iterations.
    n_phi : int
        Number of azimuth quadrature samples for the Fourier projection of the
        Rayleigh phase matrix.  240 captures up to Fourier order ~120 without
        aliasing.
    aw, bw : float | None
        Optional user-supplied absorption/scattering coefficients in m⁻¹.
        When None, the built-in interpolated table values are used.
    """

    wavelength_um: float
    depth_m: float = 100.0
    bottom_albedo: float = 0.0
    depolarization: float = 0.09
    n_streams: int = 16
    initial_tau: float = 0.01
    n_phi: int = 240
    aw: float | None = None
    bw: float | None = None

    def optical_properties(self) -> Tuple[float, float, float, float]:
        """Compute fundamental optical properties of the water column.

        Looks up (or uses the provided) absorption and scattering coefficients,
        then derives the total attenuation, optical depth, and single-scattering
        albedo.

        Returns
        -------
        aw : float
            Absorption coefficient in m⁻¹.
        bw : float
            Scattering coefficient in m⁻¹.
        tau : float
            Total vertical optical depth τ = (a_w + b_w) · depth_m.
            A larger τ means more strongly attenuating water.
        omega : float
            Single-scattering albedo ω₀ = b_w / (a_w + b_w) ∈ [0, 1].
            ω₀ = 1 means purely scattering (no absorption); ω₀ = 0 means
            purely absorbing.  For pure seawater at 500 nm, ω₀ ≈ 0.07.
        """
        aw, bw = pure_water_coefficients(self.wavelength_um)
        # Override with user-supplied values if provided
        if self.aw is not None:
            aw = float(self.aw)
        if self.bw is not None:
            bw = float(self.bw)
        c = max(aw + bw, _EPS)                           # total attenuation c = a + b
        tau = c * max(float(self.depth_m), 0.0)          # optical depth τ = c · z
        omega = float(np.clip(bw / c, 0.0, 1.0))        # single-scattering albedo ω₀
        return aw, bw, tau, omega

    def _phase_fourier(self, mus: np.ndarray, smf: np.ndarray, nmat: int, n_fourier: int, hemi_out: str) -> np.ndarray:
        """Compute the Fourier-projected Rayleigh phase supermatrix.

        Evaluates the depolarized Rayleigh scattering matrix for every pair of
        stream directions (μ_in downward, μ_out in hemisphere ``hemi_out``),
        then projects over azimuth using ``_project_operator_fourier``.

        This is used twice in ``_thin_layer``:
          - ``hemi_out="up"``   → reflection phase matrix P_ref
          - ``hemi_out="down"`` → transmission phase matrix P_tra

        Parameters
        ----------
        mus : np.ndarray
            Quadrature cosines for the water streams.
        smf : np.ndarray
            Supermatrix scaling factors for the water streams.
        nmat : int
            Number of Stokes components.
        n_fourier : int
            Highest Fourier order.
        hemi_out : str
            "up" (reflection) or "down" (transmission).

        Returns
        -------
        np.ndarray, shape (n_str·nmat, n_str·nmat, n_fourier+1)
            Fourier-projected, smf-scaled phase supermatrix.
        """
        def local(mu_out: float, mu_in: float, dphi: float) -> np.ndarray:
            # Incident light is always downwelling inside the water column.
            # For reflection the scattered direction is upwelling;
            # for transmission through the layer it remains downwelling.
            k_in  = _direction(mu_in, 0.0, "down")
            k_out = _direction(mu_out, dphi, hemi_out)  # type: ignore[arg-type]
            return _scattering_matrix_local(k_in, k_out, self.depolarization)

        return _project_operator_fourier(mus, smf, mus, smf, nmat, n_fourier, self.n_phi, local)

    def _thin_layer(self, mus: np.ndarray, smf: np.ndarray, nmat: int, n_fourier: int, tau0: float, omega: float) -> Tuple[list[np.ndarray], list[np.ndarray]]:
        """Compute reflection and transmission supermatrices for a thin slab.

        This is the initialisation step of the adding-doubling algorithm.  For
        a slab of optical depth τ₀ ≪ 1, the single-scattering approximation
        gives analytically integrable path integrals:

        **Reflection kernel** (upward scattered by a downwelling photon):
            h_R(i,j) = [1 − exp(−τ₀/μ_i) · exp(−τ₀/μ_j)] / (μ_i + μ_j)

        This integral arises from integrating the source function over the slab
        depth z ∈ [0, τ₀], accounting for attenuation of both the incoming beam
        (factor exp(−z/μ_j)) and the outgoing scattered beam (factor exp(−(τ₀−z)/μ_i)
        for upward direction → integrates to the expression above).

        **Transmission kernel** (downward scattered, exiting bottom):
            h_T(i,j) = [exp(−τ₀/μ_j) − exp(−τ₀/μ_i)] / (μ_i − μ_j)   for μ_i ≠ μ_j
                     = (τ₀/μ_i²) exp(−τ₀/μ_i)                            for μ_i = μ_j  (L'Hôpital)

        The L'Hôpital form for equal angles avoids 0/0 on the diagonal.

        Both kernels are then multiplied by (ω₀/4) · phase_matrix^(m) to give
        the multiple-scattering-free reflection/transmission operators.

        For the m=0 Fourier term, the direct-transmission diagonal is added:
            T[i·nmat+k, i·nmat+k] += exp(−τ₀/μ_i)

        This represents photons that pass straight through without being
        scattered.

        Parameters
        ----------
        mus : np.ndarray, shape (n_str,)
            Stream cosines for the water-body quadrature.
        smf : np.ndarray, shape (n_str,)
            Supermatrix scaling factors.
        nmat : int
            Number of Stokes components.
        n_fourier : int
            Highest Fourier order.
        tau0 : float
            Optical depth of the thin starting slab (≪ 1 for accuracy).
        omega : float
            Single-scattering albedo ω₀ of the water.

        Returns
        -------
        R_list : list of np.ndarray, length (n_fourier+1)
            R_list[m] = m-th Fourier order of the reflection supermatrix.
        T_list : list of np.ndarray, length (n_fourier+1)
            T_list[m] = m-th Fourier order of the transmission supermatrix.
        """
        nsup = len(mus) * nmat
        # Phase matrices for reflection (upward) and transmission (downward)
        phase_ref = self._phase_fourier(mus, smf, nmat, n_fourier, "up")
        phase_tra = self._phase_fourier(mus, smf, nmat, n_fourier, "down")

        # ── Vectorised hR / hT over all stream pairs at once ────────────
        # e_mu[i] = exp(−τ₀/μ_i) = direct transmittance for stream i
        e_mu = np.exp(-tau0 / np.clip(mus, _EPS, None))    # shape (n_str,)

        # hR[i,j] = (1 - e_i · e_j) / (μ_i + μ_j)
        # Numerator: 1 − exp(−τ₀/μ_i − τ₀/μ_j)  (two-path attenuation)
        # Denominator: μ_i + μ_j  (geometrical factor for reflection path)
        hR = (1.0 - np.outer(e_mu, e_mu)) / (mus[:, None] + mus[None, :] + _EPS)

        # hT[i,j] = (e_j − e_i) / (μ_i − μ_j)  with L'Hôpital when μ_i = μ_j.
        # Note: the ordering (e_i[:, None] − e_j[None, :]) uses i=out, j=in.
        dmu = mus[:, None] - mus[None, :]                   # (n_str, n_str)
        off_diag = np.abs(dmu) > 1.0e-10                    # True where μ_i ≠ μ_j
        hT = np.where(
            off_diag,
            # Off-diagonal: standard formula
            (e_mu[:, None] - e_mu[None, :]) / np.where(off_diag, dmu, 1.0),
            # Diagonal (μ_i = μ_j): L'Hôpital limit = τ₀/μ_i² · exp(−τ₀/μ_i)
            tau0 * e_mu[:, None] / (mus[:, None] ** 2 + _EPS),
        )

        # Expand scalar hR/hT matrices from (n_str, n_str) to (nsup, nsup) by
        # tiling each element into an (nmat, nmat) block via Kronecker product.
        # hR_block[i·nmat:(i+1)·nmat, j·nmat:(j+1)·nmat] = hR[i,j] · ones(nmat,nmat)
        hR_block = np.kron(hR, np.ones((nmat, nmat)))   # (nsup, nsup)
        hT_block = np.kron(hT, np.ones((nmat, nmat)))

        R_list: list[np.ndarray] = []
        T_list: list[np.ndarray] = []

        for m in range(n_fourier + 1):
            # Single-scattering reflection/transmission: (ω₀/4) · h · P^(m)
            # The factor 1/4 = 1/(4π) normalisation constant for the phase function
            # in the discrete-ordinate quadrature scheme.
            R = 0.25 * omega * hR_block * phase_ref[:, :, m]
            T = 0.25 * omega * hT_block * phase_tra[:, :, m]
            if m == 0:
                # The m=0 Fourier term includes the azimuth-independent
                # direct-beam transmittance on the diagonal: exp(−τ₀/μ_i)
                # for each Stokes component k of stream i.
                idx = np.arange(len(mus)) * nmat
                for k in range(nmat):
                    T[idx + k, idx + k] += e_mu
            R_list.append(R)
            T_list.append(T)
        return R_list, T_list

    @staticmethod
    def _double_layer(R: np.ndarray, T: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Apply one doubling step: combine two identical layers into one.

        The doubling formulas (van de Hulst 1963, de Haan et al. 1987) combine
        two identical layers (each with reflection R and transmission T) into a
        single layer of twice the optical depth:

            R₂ = R + T · (I − R²)^{-1} · R · T
            T₂ = T · (I − R²)^{-1} · T

        The factor (I − R²)^{-1} = (I − R·R)^{-1} accounts for the infinite
        series of multiple reflections between the two sub-layers:
            (I − R·R)^{-1} = I + R² + R⁴ + R⁶ + ...

        Each term in the series represents one additional round-trip bounce
        between the top surface of the lower sub-layer and the bottom surface
        of the upper sub-layer.

        The matrix solve ``np.linalg.solve(I − R², I)`` is used instead of
        explicit matrix inversion for numerical stability (avoids forming the
        inverse explicitly).

        Parameters
        ----------
        R : np.ndarray, shape (nsup, nsup)
            Reflection supermatrix of a single sub-layer.
        T : np.ndarray, shape (nsup, nsup)
            Transmission supermatrix of a single sub-layer.

        Returns
        -------
        R2 : np.ndarray, shape (nsup, nsup)
            Reflection supermatrix of the doubled layer.
        T2 : np.ndarray, shape (nsup, nsup)
            Transmission supermatrix of the doubled layer.
        """
        I = np.eye(R.shape[0])
        # Solve (I − R²) X = I  →  X = (I − R²)^{-1}; more stable than np.linalg.inv
        den = np.linalg.solve(I - R @ R, np.eye(R.shape[0]))
        R2 = R + T @ den @ R @ T   # combined reflection
        T2 = T @ den @ T           # combined transmission
        return R2, T2

    def _add_bottom(self, R: np.ndarray, T: np.ndarray, mus: np.ndarray, smf: np.ndarray, nmat: int, m: int) -> np.ndarray:
        """Add a Lambertian bottom boundary to the water-column reflection operator.

        Uses the adding formulas to combine the water-column operators (R, T)
        with a Lambertian bottom boundary of albedo A_bot:

            R_total = R + T · (I − R_bot · R)^{-1} · R_bot · T

        The Lambertian bottom has a non-polarizing, isotropic (m=0 only)
        reflection operator in the supermatrix basis:

            R_bot[i·nmat, j·nmat] = A_bot · smf[i] · smf[j]   (intensity only)

        All higher Fourier orders (m ≠ 0) are zero for an isotropic bottom, so
        the function returns ``R`` unchanged for m > 0.

        Parameters
        ----------
        R : np.ndarray, shape (nsup, nsup)
            Reflection supermatrix of the water column alone.
        T : np.ndarray, shape (nsup, nsup)
            Transmission supermatrix of the water column.
        mus : np.ndarray
            Water stream cosines (used for sizing).
        smf : np.ndarray
            Supermatrix scaling factors for the water streams.
        nmat : int
            Number of Stokes components.
        m : int
            Current Fourier order.  Only m=0 is modified.

        Returns
        -------
        np.ndarray, shape (nsup, nsup)
            Updated reflection supermatrix including the bottom contribution.
        """
        if self.bottom_albedo <= 0.0 or m != 0:
            # No bottom contribution, or Fourier order > 0 (Lambertian is isotropic).
            return R
        nsup = len(mus) * nmat
        # Build Lambertian bottom reflection supermatrix (non-polarizing, m=0 only).
        # Only the I→I block (indices divisible by nmat) is non-zero.
        Rb = np.zeros((nsup, nsup), dtype=float)
        idx = np.arange(len(mus)) * nmat   # indices of the I-component of each stream
        Rb[np.ix_(idx, idx)] = float(self.bottom_albedo) * np.outer(smf, smf)
        I = np.eye(nsup)
        # Adding formula: R_total = R_col + T · (I − R_bot · R_col)^{-1} · R_bot · T
        den = np.linalg.solve(I - Rb @ R, np.eye(nsup))
        return R + T @ den @ Rb @ T

    def reflection_supermatrices(self, nmat: int = 4, n_fourier: int = 2) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute the full polarized reflection supermatrix R_W for the water body.

        This is the main entry point for ``WaterBody``.  It executes the
        complete adding-doubling calculation and returns the result in a format
        ready for the Trees & Stam coupling step in ``OceanSurface``.

        Algorithm summary
        -----------------
        1. Compute quadrature nodes, weights, and smf factors for n_streams.
        2. Retrieve optical properties (a_w, b_w, τ, ω₀).
        3. Determine the number of doubling steps N = ceil(log₂(τ / τ_init)).
        4. Compute the thin-slab operators R^(m), T^(m) for τ₀ = τ / 2^N.
        5. Double N times: apply ``_double_layer`` N times per Fourier term.
        6. Add the bottom boundary using ``_add_bottom`` (for m=0 only).
        7. Assemble the result into the output array R_W.

        Parameters
        ----------
        nmat : int
            Number of Stokes components (1, 3, or 4).  Default 4 = full Stokes.
        n_fourier : int
            Number of Fourier terms to compute.  For pure Rayleigh scattering
            only m=0, 1, 2 are non-zero, but higher terms may be requested if
            needed for the atmosphere coupling.  Extra terms are numerically zero.

        Returns
        -------
        Rw : np.ndarray, shape (n_streams·nmat, n_streams·nmat, n_fourier+1), Fortran order
            Fourier-decomposed, smf-scaled reflection supermatrix of the water body.
            Maps downwelling water streams to upwelling water streams.
        mus : np.ndarray, shape (n_streams,)
            Quadrature cosines used for the water streams.
        weights : np.ndarray, shape (n_streams,)
            Quadrature weights.
        smf : np.ndarray, shape (n_streams,)
            Supermatrix scaling factors.
        """
        mus, weights, smf = _quadrature(self.n_streams, include_extra_mu_one=False)
        aw, bw, tau, omega = self.optical_properties()
        nsup = len(mus) * nmat
        Rw = np.zeros((nsup, nsup, n_fourier + 1), dtype=np.float64, order="F")

        if tau <= 0.0 or omega <= 0.0:
            # Optically trivial case: no optical depth or purely absorbing water.
            # The water body reflects nothing; return a zero operator.
            return Rw, mus, weights, smf

        # Determine number of doubling iterations: N = ceil(log2(τ / τ₀_init))
        # so that the starting slab thickness τ₀ = τ / 2^N ≤ initial_tau.
        ndouble = max(0, int(np.ceil(np.log2(max(tau / max(self.initial_tau, _EPS), 1.0)))))
        tau0 = tau / (2**ndouble)   # optical depth of the thin starting slab

        # Step 1: initialise the thin-slab operators for each Fourier term
        R_list, T_list = self._thin_layer(mus, smf, nmat, n_fourier, tau0, omega)

        for m in range(n_fourier + 1):
            R = R_list[m]
            T = T_list[m]
            # Steps 2–N: double the layer thickness N times
            for _ in range(ndouble):
                R, T = self._double_layer(R, T)
            # Step N+1: add the bottom Lambertian boundary (m=0 only)
            R = self._add_bottom(R, T, mus, smf, nmat, m)
            Rw[:, :, m] = R
        return Rw, mus, weights, smf


# ---------------------------------------------------------------------------
# Rough refractive interface and clean-ocean coupling
# ---------------------------------------------------------------------------


@dataclass
class RoughInterface:
    """Rough Cox-Munk air-water interface polarized reflection and transmission operators.

    Physical model
    --------------
    The ocean surface is modelled as an ensemble of randomly oriented planar
    microfacets whose slope distribution follows the isotropic Cox & Munk (1954)
    Gaussian model.  Each facet reflects and refracts light according to the
    Fresnel equations for a flat dielectric interface between air (n=n_air) and
    water (n=n_water).

    The bidirectional reflection distribution function (BRDF) for a Cox-Munk
    surface is (Beckmann & Spizzichino 1963, Cox & Munk 1954):

        f_r(μ_out, φ_out; μ_in) = S(μ_out, μ_in; σ) · π · P(μ_n; σ) ·
                                    M_Fresnel(μ_i) / (4 μ_out μ_in μ_n)

    where:
      - μ_n = cos θ_n is the cosine of the facet-normal zenith angle,
        determined by the geometry: n̂ = (k̂_out − k̂_in)/|k̂_out − k̂_in|
        for reflection.
      - P(μ_n; σ) is the Cox-Munk slope PDF.
      - S(μ_out, μ_in; σ) is the Smith/Sancer shadowing factor.
      - M_Fresnel(μ_i) is the 4×4 Fresnel Mueller matrix evaluated at the
        local incidence angle μ_i = |k̂_in · n̂|.
      - The prefactor (4 μ_out μ_in μ_n)^{-1} converts from the slope PDF
        (defined per unit solid angle of facet normals) to the BRDF
        (defined per unit projected area × unit solid angle of reflection).

    For transmission (refraction), the same geometrical framework applies but
    the facet normal is determined by Snell's law (vector form):

        n̂ = (n1 k̂_in − n2 k̂_out) / |n1 k̂_in − n2 k̂_out|

    and an additional solid-angle Jacobian

        J = n2² μ_t μ_i / (n2 μ_t − n1 μ_i)²

    accounts for the compression/expansion of the solid-angle element as light
    crosses the refractive interface.

    Polarisation-accurate Mueller matrix
    -------------------------------------
    Because the facet normal is oriented relative to the ray pair (not the
    global vertical), the Fresnel Mueller matrix must be sandwiched between two
    Mueller rotation matrices to transform to/from the meridian-plane
    (global Stokes) reference frame:

        M_local = L(−ξ₂) · M_Fresnel · L(−ξ₁)

    where ξ₁ = angle from incoming meridian plane to facet plane (around k_in),
    and ξ₂ = angle from facet plane to outgoing meridian plane (around k_out).

    Vectorised vs scalar paths
    --------------------------
    Each method has two implementations:
      - **Scalar** ``_reflection_local`` / ``_transmission_local``: evaluates
        the matrix for a single (μ_out, μ_in, φ) triple.  Used by the public
        wrappers and when building supermatrices with the slow path.
      - **Vectorised** ``_reflection_vec`` / ``_transmission_vec``: evaluates
        the matrix for a single (μ_out, μ_in) pair but all azimuth angles φ
        simultaneously, returning an (n_phi, 4, 4) array in one NumPy call.
        ~50× faster; used as the ``local_matrix_vec`` argument to
        ``_project_operator_fourier``.

    Parameters
    ----------
    wind_speed : float
        Wind speed U_{10} in m/s.  Controls the Cox-Munk slope variance σ².
    n_air : float
        Refractive index of the air-side medium (usually 1.0).
    n_water : float
        Refractive index of the water-side medium (usually ~1.33 for seawater).

    Related classes
    ---------------
    ``WaterBody`` — provides R_W, the water-column reflection operator.
    ``OceanSurface`` — couples ``RoughInterface`` with ``WaterBody`` using the
    Trees & Stam (2019) adding formula to produce the final ocean operator.
    """

    wind_speed: float = 7.0
    n_air: float = 1.0
    n_water: float = 1.33

    @property
    def sigma(self) -> float:
        """RMS surface slope σ = √σ² derived from the Cox-Munk wind-speed relation.

        Returns
        -------
        float
            Root-mean-square wave slope (dimensionless).
        """
        return float(np.sqrt(cox_munk_slope_variance(self.wind_speed)))

    def _reflection_local(self, mu_out: float, mu_in: float, dphi: float, from_below: bool = False) -> np.ndarray:
        """Compute the local Cox-Munk reflection Mueller matrix for one (μ_out, μ_in, φ) triple.

        This is the scalar (single-azimuth) reflection kernel.  For a batch
        over all azimuth angles use ``_reflection_vec``.

        Geometry
        --------
        For air-side reflection (from_below=False):
          - k_in  = downwelling ray at azimuth φ=0 (reference direction)
          - k_out = upwelling reflected ray at azimuth φ=dphi
          - facet normal n̂ = normalise(k_out − k_in)  (bisector of incident + reflected)

        For water-side reflection (from_below=True, TIR path):
          - k_in  = upwelling ray (in water) at azimuth φ=0
          - k_out = downwelling reflected ray at azimuth φ=dphi
          - facet normal n̂ = normalise(k_in − k_out)  (upward-pointing bisector)

        The scale factor is:
            scale = S(μ_out, μ_in; σ) · π · P(μ_n; σ) / (4 μ_out μ_in μ_n)

        Parameters
        ----------
        mu_out : float
            |cos θ_out| — cosine of the outgoing (reflected) polar angle.
        mu_in : float
            |cos θ_in| — cosine of the incoming polar angle.
        dphi : float
            Relative azimuth angle φ_out − φ_in in radians.
        from_below : bool
            If True, computes the water→air (conjugate) reflection R_I*.
            If False, computes the air→water reflection R_I.

        Returns
        -------
        np.ndarray, shape (4, 4)
            4×4 Mueller matrix (in the meridian-plane frame) for this geometry.
            Returns zeros if the geometry is degenerate (μ_n ≤ 0).
        """
        mu_out = float(np.clip(mu_out, _EPS, 1.0))
        mu_in  = float(np.clip(mu_in,  _EPS, 1.0))
        if from_below:
            # Water-to-air path: incoming ray upwells, outgoing ray downwells.
            k_in  = _direction(mu_in,  0.0,  "up")
            k_out = _direction(mu_out, dphi, "down")
            # Facet normal: bisector of in-going and reversed outgoing (upward-pointing)
            normal = _safe_normalize(k_in - k_out)
            n1, n2 = self.n_water, self.n_air
        else:
            # Air-to-water path: incoming ray downwells, outgoing ray upwells.
            k_in  = _direction(mu_in,  0.0,  "down")
            k_out = _direction(mu_out, dphi, "up")
            # Facet normal: bisector of outgoing and reversed incoming (upward-pointing)
            normal = _safe_normalize(k_out - k_in)
            n1, n2 = self.n_air, self.n_water

        mu_n = float(normal[2])          # cos θ_n = z-component of the facet normal
        if mu_n <= 0.0:
            # Facet pointing downward or horizontal — unphysical; return zero.
            return np.zeros((4, 4), dtype=float)
        mu_i = float(np.clip(abs(np.dot(k_in, normal)), 0.0, 1.0))  # local incidence angle on facet
        if mu_i <= 0.0:
            return np.zeros((4, 4), dtype=float)

        # BRDF prefactor = shadowing × Cox-Munk PDF / geometric factor
        pdf    = cox_munk_pdf(mu_n, self.sigma)
        shadow = smith_sancer_shadowing(mu_out, mu_in, self.sigma)
        scale  = shadow * np.pi * pdf / (4.0 * mu_out * mu_in * mu_n + _EPS)

        # Rotation angles to/from the facet reference plane
        mer_in   = _meridian_plane_normal(k_in)
        mer_out  = _meridian_plane_normal(k_out)
        # facet_plane normal = n̂ × k̂_out (perpendicular to both normal and k_out)
        facet_plane = _safe_normalize(np.cross(normal, k_out))
        if np.linalg.norm(facet_plane) < _EPS:
            facet_plane = mer_out  # degenerate: use outgoing meridian as fallback
        # ξ₁ = rotation from incoming meridian → facet plane (around k_in)
        xi1 = _signed_angle_between_plane_normals(mer_in, facet_plane, k_in)
        # ξ₂ = rotation from facet plane → outgoing meridian (around k_out)
        xi2 = _signed_angle_between_plane_normals(facet_plane, mer_out, k_out)
        rf = fresnel_reflection_mueller(mu_i, n1, n2)
        # Full Mueller matrix: rotate to facet frame, apply Fresnel, rotate back
        return scale * (_rotation_matrix(-xi2) @ rf @ _rotation_matrix(-xi1))

    def _transmission_local(self, mu_out: float, mu_in: float, dphi: float, upward: bool = False) -> np.ndarray:
        """Compute the local Cox-Munk transmission Mueller matrix for one (μ_out, μ_in, φ) triple.

        This is the scalar (single-azimuth) transmission kernel.  For a batch
        over all azimuth angles use ``_transmission_vec``.

        ``upward=False`` maps air-downwelling → water-downwelling (T_I, downward refraction).
        ``upward=True``  maps water-upwelling → air-upwelling  (T_I*, upward refraction).

        Facet normal via vector Snell's law
        ------------------------------------
        For a refracting interface the facet normal is not simply the bisector
        of k_in and k_out; instead it is determined by Snell's law in vector form:

            n1 k̂_in − n2 k̂_out ∥ n̂

        (this can be derived from the boundary condition that the tangential
        components of the wave vector are continuous across the interface).
        The sign is chosen so that n̂ points upward (positive z component).

        Solid-angle Jacobian
        --------------------
        When light is refracted, the solid angle element dΩ changes according
        to the Jacobian:

            J = n2² μ_t μ_i / (n2 μ_t − n1 μ_i)²

        This factor arises from the change in the mapping between the incoming
        and outgoing angular coordinates and must multiply the scale factor.
        See Zhai et al. (2010, Appl. Opt. 49, 2441), eq. (5).

        Parameters
        ----------
        mu_out : float
            |cos θ_out| — cosine of the outgoing (transmitted) polar angle.
        mu_in : float
            |cos θ_in| — cosine of the incoming polar angle.
        dphi : float
            Relative azimuth angle in radians.
        upward : bool
            Direction of transmission:
              - False (default): air → water (T_I, downwelling refraction)
              - True: water → air (T_I*, upwelling refraction)

        Returns
        -------
        np.ndarray, shape (4, 4)
            4×4 Mueller matrix (meridian-plane frame) for this geometry.
            Returns zeros if TIR applies or the geometry is degenerate.
        """
        mu_out = float(np.clip(mu_out, _EPS, 1.0))
        mu_in  = float(np.clip(mu_in,  _EPS, 1.0))
        if upward:
            # Water → air: both incoming and outgoing rays point upward
            k_in  = _direction(mu_in,  0.0,  "up")
            k_out = _direction(mu_out, dphi, "up")
            n1, n2 = self.n_water, self.n_air
        else:
            # Air → water: both incoming and outgoing rays point downward
            k_in  = _direction(mu_in,  0.0,  "down")
            k_out = _direction(mu_out, dphi, "down")
            n1, n2 = self.n_air, self.n_water

        # Vector Snell: facet normal is parallel to (n1 k̂_in − n2 k̂_out)
        normal = _safe_normalize(n1 * k_in - n2 * k_out)
        if normal[2] < 0.0:
            normal = -normal    # enforce upward-pointing normal convention
        mu_n = float(normal[2])
        if mu_n <= 0.0:
            return np.zeros((4, 4), dtype=float)
        mu_i = float(np.clip(abs(np.dot(k_in, normal)), 0.0, 1.0))  # local incidence on facet
        rs_amp, rp_amp, mu_t, tir = _fresnel_amplitude_coefficients(mu_i, n1, n2)
        if tir or mu_t <= 0.0:
            # Total internal reflection: no transmitted ray possible.
            return np.zeros((4, 4), dtype=float)

        pdf    = cox_munk_pdf(mu_n, self.sigma)
        shadow = smith_sancer_shadowing(mu_out, mu_in, self.sigma)
        # Solid-angle Jacobian for refraction (Zhai et al. 2010, eq. 5)
        jac   = (n2 * n2 * mu_t * mu_i) / ((n2 * mu_t - n1 * mu_i) ** 2 + _EPS)
        # Transmission scale factor: shadowing × π × PDF × Jacobian / (geometry)
        scale = shadow * np.pi * pdf * jac / (mu_in * mu_out * mu_n + _EPS)

        # Rotation angles to/from the facet reference frame
        mer_in   = _meridian_plane_normal(k_in)
        mer_out  = _meridian_plane_normal(k_out)
        facet_plane = _safe_normalize(np.cross(normal, k_out))
        if np.linalg.norm(facet_plane) < _EPS:
            facet_plane = mer_out
        xi1 = _signed_angle_between_plane_normals(mer_in, facet_plane, k_in)
        xi2 = _signed_angle_between_plane_normals(facet_plane, mer_out, k_out)
        tf = fresnel_transmission_mueller(mu_i, n1, n2)
        # Full Mueller matrix: rotate to facet frame, apply Fresnel transmission, rotate back
        return scale * (_rotation_matrix(-xi2) @ tf @ _rotation_matrix(-xi1))

    def reflection_air(self, mu_out: float, mu_in: float, dphi: float) -> np.ndarray:
        """Air-side (R_I): Cox-Munk reflection from air-to-water, viewed from above.

        Computes the reflection operator for a downwelling incident ray (in air)
        scattering into an upwelling reflected ray (in air).  This is the
        standard sun-glint reflection operator.

        Parameters
        ----------
        mu_out : float
            |cos θ_out| of the upwelling reflected ray.
        mu_in : float
            |cos θ_in| of the downwelling incident ray.
        dphi : float
            Relative azimuth angle φ_out − φ_in in radians.

        Returns
        -------
        np.ndarray, shape (4, 4)
            Mueller matrix for specular reflection off the rough ocean surface.
        """
        return self._reflection_local(mu_out, mu_in, dphi, from_below=False)

    def reflection_water(self, mu_out: float, mu_in: float, dphi: float) -> np.ndarray:
        """Water-side (R_I*): conjugate Cox-Munk reflection from water-to-air viewed from below.

        Computes the reflection operator for an upwelling incident ray (in water)
        that is totally internally reflected (or partially reflected) back downward
        (in water).  This is the R_I* operator in the Trees & Stam formula, which
        represents the interface's opacity to upwelling water-body radiation.

        Parameters
        ----------
        mu_out : float
            |cos θ_out| of the downwelling reflected ray (in water).
        mu_in : float
            |cos θ_in| of the upwelling incident ray (in water).
        dphi : float
            Relative azimuth angle in radians.

        Returns
        -------
        np.ndarray, shape (4, 4)
            Mueller matrix for reflection at the interface from below.
        """
        return self._reflection_local(mu_out, mu_in, dphi, from_below=True)

    def transmission_down(self, mu_out_water: float, mu_in_air: float, dphi: float) -> np.ndarray:
        """Downward transmission T_I: maps air-downwelling to water-downwelling.

        Computes the refraction operator for a downwelling ray in air entering
        the water (Snell refraction, air → water).  This is the T_I operator in
        the Trees & Stam formula, driving sunlight into the water body.

        Parameters
        ----------
        mu_out_water : float
            |cos θ_out| of the downwelling refracted ray inside water.
        mu_in_air : float
            |cos θ_in| of the downwelling incident ray in air.
        dphi : float
            Relative azimuth angle in radians.

        Returns
        -------
        np.ndarray, shape (4, 4)
            Mueller matrix for downward refraction at the rough interface.
        """
        return self._transmission_local(mu_out_water, mu_in_air, dphi, upward=False)

    def transmission_up(self, mu_out_air: float, mu_in_water: float, dphi: float) -> np.ndarray:
        """Upward transmission T_I*: maps water-upwelling to air-upwelling.

        Computes the refraction operator for an upwelling ray in water exiting
        into air (Snell refraction, water → air).  This is the T_I* operator in
        the Trees & Stam formula, transmitting water-leaving radiance to the
        observer.

        Parameters
        ----------
        mu_out_air : float
            |cos θ_out| of the upwelling refracted ray in air.
        mu_in_water : float
            |cos θ_in| of the upwelling incident ray inside water.
        dphi : float
            Relative azimuth angle in radians.

        Returns
        -------
        np.ndarray, shape (4, 4)
            Mueller matrix for upward refraction at the rough interface.
        """
        return self._transmission_local(mu_out_air, mu_in_water, dphi, upward=True)

    # ------------------------------------------------------------------
    # Vectorized-over-phi methods (same physics, ~50× faster than scalar
    # loop when called for all azimuth angles at once).
    # ------------------------------------------------------------------

    def _reflection_vec(
        self, mu_out: float, mu_in: float, phis: np.ndarray, from_below: bool = False
    ) -> np.ndarray:
        """Reflection Mueller matrix vectorised over all azimuth angles simultaneously.

        This method computes the same physics as ``_reflection_local`` but for
        the full array of n_phi azimuth angles in a single NumPy call, avoiding
        a Python loop over phi.  It is ~50× faster for large n_phi and is
        supplied to ``_project_operator_fourier`` as ``local_matrix_vec``.

        Implementation strategy
        -----------------------
        1. Build k_out as an (n_phi, 3) array by stacking the azimuth-varying
           x/y components with a constant ±μ_out z-component.
        2. Compute the facet normals as (n_phi, 3) arrays.
        3. Compute the scale factors (PDF, shadowing) as (n_phi,) arrays.
        4. Vectorise the Fresnel coefficients over the (n_phi,) array of local
           incidence angles μ_i.
        5. Compute the rotation angles ξ₁, ξ₂ using ``_signed_angle_vec``.
        6. Build the (n_phi, 4, 4) product L(−ξ₂) · M_Fresnel · L(−ξ₁) without
           a Python loop using pre-expanded matrix algebra.

        Vectorised matrix multiply detail (step 6)
        -------------------------------------------
        Rather than calling ``np.einsum`` or a batch matrix multiply, the
        product L(−ξ₂) · R_F · L(−ξ₁) is assembled analytically by exploiting
        the sparse structure of L(ξ) and R_F:

        First compute RL = R_F · L(−ξ₁):
          - L(−ξ₁) acts on COLUMNS 1 and 2 (the Q and U columns):
              col_1_new = c1 · col_1_old − s1 · col_2_old  (note: but L mixes differently)
              col_2_new = s1 · col_1_old + c1 · col_2_old
          - Since R_F has the block form [[a,b],[b,a]] in (0,1):
              RL[0,0]=a, RL[0,1]=b·c1, RL[0,2]=b·s1
              RL[1,0]=b, RL[1,1]=a·c1, RL[1,2]=a·s1
              RL[2,1]=−c_f·s1, RL[2,2]=c_f·c1, RL[3,3]=c_f

        Then compute L(−ξ₂) · RL:
          - L(−ξ₂) acts on ROWS 1 and 2:
              row_1_new =  c2 · old_row_1 + s2 · old_row_2
              row_2_new = −s2 · old_row_1 + c2 · old_row_2

        This is done in-place on the ``result`` copy.

        Parameters
        ----------
        mu_out : float
            |cos θ_out| of the outgoing ray.
        mu_in : float
            |cos θ_in| of the incoming ray.
        phis : np.ndarray, shape (n_phi,)
            Array of relative azimuth angles φ in radians.
        from_below : bool
            False = air-side reflection R_I; True = water-side reflection R_I*.

        Returns
        -------
        np.ndarray, shape (n_phi, 4, 4)
            Array of 4×4 Mueller matrices, one per azimuth angle.
            Zero matrices are returned for degenerate geometries.
        """
        n_phi = len(phis)
        mu_out = float(np.clip(mu_out, _EPS, 1.0))
        mu_in  = float(np.clip(mu_in,  _EPS, 1.0))
        sin_in  = float(np.sqrt(max(0.0, 1.0 - mu_in  * mu_in)))   # sin θ_in
        sin_out = float(np.sqrt(max(0.0, 1.0 - mu_out * mu_out)))   # sin θ_out

        if from_below:
            # Water→air path: k_in upwelling (z>0), k_out downwelling (z<0)
            k_in = np.array([sin_in, 0.0,  mu_in],  dtype=float)
            k_out = np.column_stack([sin_out * np.cos(phis),
                                     sin_out * np.sin(phis),
                                     np.full(n_phi, -mu_out)])
            n1, n2 = self.n_water, self.n_air
            # Facet normal: bisector pointing upward for TIR geometry
            diff = k_in[np.newaxis, :] - k_out          # normal = norm(k_in - k_out)
        else:
            # Air→water path: k_in downwelling (z<0), k_out upwelling (z>0)
            k_in  = np.array([sin_in, 0.0, -mu_in], dtype=float)
            k_out = np.column_stack([sin_out * np.cos(phis),
                                     sin_out * np.sin(phis),
                                     np.full(n_phi, mu_out)])
            n1, n2 = self.n_air, self.n_water
            diff = k_out - k_in[np.newaxis, :]           # normal = norm(k_out - k_in)

        # --- Facet normal geometry ---
        mag = np.linalg.norm(diff, axis=1, keepdims=True)           # (n_phi, 1)
        valid = (mag.squeeze() > _EPS)                               # boolean mask
        normal = np.where(mag > _EPS, diff / (mag + _EPS), 0.0)     # (n_phi, 3)

        mu_n = normal[:, 2]          # z-component of facet normal = cos θ_n  (n_phi,)
        valid &= (mu_n > 0.0)        # require upward-pointing facet

        # Local angle of incidence on the facet: μ_i = |k̂_in · n̂|
        mu_i_arr = np.abs(np.einsum('i,pi->p', k_in, normal))       # (n_phi,)
        valid &= (mu_i_arr > 0.0)

        # --- Scale factors ---
        # Shadowing: scalar (independent of azimuth for isotropic Cox-Munk)
        shadow = smith_sancer_shadowing(mu_out, mu_in, self.sigma)
        sig = self.sigma

        # Cox-Munk PDF (varies with μ_n across azimuths)
        mu_n_s = np.where(valid, np.clip(mu_n, _EPS, None), 1.0)
        tan2 = (1.0 - mu_n_s * mu_n_s) / (mu_n_s * mu_n_s * sig * sig + _EPS)  # tan²θ_n / σ²
        pdf  = np.where(valid, np.exp(-tan2) / (np.pi * sig * sig * mu_n_s ** 3 + _EPS), 0.0)
        scale = np.where(valid, shadow * np.pi * pdf / (4.0 * mu_out * mu_in * mu_n_s + _EPS), 0.0)

        # --- Fresnel reflection coefficients (vectorised over all azimuths) ---
        sin_i2 = np.maximum(0.0, 1.0 - mu_i_arr * mu_i_arr)       # sin²θ_i
        sin_t2 = (n1 / n2) ** 2 * sin_i2                           # sin²θ_t (Snell)
        tir    = sin_t2 >= 1.0                                      # TIR mask
        mu_t   = np.sqrt(np.maximum(0.0, 1.0 - sin_t2))            # cos θ_t
        # R_s = r_s²,  R_p = r_p² (set to 1 under TIR)
        rs = np.where(tir, 1.0, ((n1 * mu_i_arr - n2 * mu_t) / (n1 * mu_i_arr + n2 * mu_t + _EPS)) ** 2)
        rp = np.where(tir, 1.0, ((n2 * mu_i_arr - n1 * mu_t) / (n2 * mu_i_arr + n1 * mu_t + _EPS)) ** 2)

        # --- Meridian planes and facet plane (all as (n_phi, 3) arrays) ---
        z = np.array([0.0, 0.0, 1.0])
        # Incoming meridian plane normal: constant over phi (k_in is fixed at phi=0)
        mer_in_v   = np.cross(z, k_in)
        mer_in_mag = np.linalg.norm(mer_in_v)
        mer_in_v   = mer_in_v / mer_in_mag if mer_in_mag > _EPS else np.array([0.0, 1.0, 0.0])
        mer_in_rep = np.tile(mer_in_v, (n_phi, 1))                  # (n_phi, 3): broadcast

        # Outgoing meridian plane normal: varies with azimuth
        mer_out_c = np.cross(z[np.newaxis, :], k_out)               # (n_phi, 3)
        mer_out_m = np.linalg.norm(mer_out_c, axis=1, keepdims=True)
        mer_out   = np.where(mer_out_m > _EPS, mer_out_c / (mer_out_m + _EPS),
                             np.array([0.0, 1.0, 0.0]))

        # Facet plane normal = n̂ × k̂_out
        facet_c = np.cross(normal, k_out)                            # (n_phi, 3)
        facet_m = np.linalg.norm(facet_c, axis=1, keepdims=True)
        facet_p = np.where(facet_m > _EPS, facet_c / (facet_m + _EPS), mer_out)

        k_in_rep = np.tile(k_in, (n_phi, 1))                        # (n_phi, 3)
        xi1 = _signed_angle_vec(mer_in_rep, facet_p, k_in_rep)      # (n_phi,)
        xi2 = _signed_angle_vec(facet_p,    mer_out, k_out)         # (n_phi,)

        # --- Vectorised matrix product: scale · L(−ξ₂) · R_F · L(−ξ₁) ---
        c1, s1 = np.cos(-2.0 * xi1), np.sin(-2.0 * xi1)   # trig for L(−ξ₁)
        c2, s2 = np.cos(-2.0 * xi2), np.sin(-2.0 * xi2)   # trig for L(−ξ₂)
        a_f = 0.5 * (rs + rp)           # Fresnel average reflectance
        b_f = 0.5 * (rp - rs)           # (R_p − R_s)/2 < 0 → correct Q < 0 sign
        c_f = np.sqrt(np.maximum(rs * rp, 0.0))

        # Step 1: RL = R_F · L(−ξ₁)
        # L(−ξ₁) rotates columns 1 and 2 of R_F by mixing with trig factors c1/s1.
        # R_F is block-diagonal: [[a,b],[b,a]] for I/Q and [[c,0],[0,c]] for U/V.
        RL = np.zeros((n_phi, 4, 4))
        RL[:, 0, 0] = a_f;       RL[:, 0, 1] = b_f * c1;   RL[:, 0, 2] = b_f * s1
        RL[:, 1, 0] = b_f;       RL[:, 1, 1] = a_f * c1;   RL[:, 1, 2] = a_f * s1
        RL[:, 2, 1] = -c_f * s1; RL[:, 2, 2] = c_f * c1;   RL[:, 3, 3] = c_f
        # (rows 0 and 3 of RL are correct above; only rows 1 and 2 mix in step 2)

        # Step 2: result = L(−ξ₂) · RL
        # L(−ξ₂) acts on ROWS 1 and 2 of RL (rows 0 and 3 are unchanged).
        # New row 1 =  c2 · old_row_1 + s2 · old_row_2
        # New row 2 = −s2 · old_row_1 + c2 · old_row_2
        result = RL.copy()
        result[:, 1, :] =  c2[:, np.newaxis] * RL[:, 1, :] + s2[:, np.newaxis] * RL[:, 2, :]
        result[:, 2, :] = -s2[:, np.newaxis] * RL[:, 1, :] + c2[:, np.newaxis] * RL[:, 2, :]
        result *= scale[:, np.newaxis, np.newaxis]   # apply BRDF prefactor
        result[~valid] = 0.0                         # zero out degenerate geometries
        return result

    def _transmission_vec(
        self, mu_out: float, mu_in: float, phis: np.ndarray, upward: bool = False
    ) -> np.ndarray:
        """Transmission Mueller matrix vectorised over all azimuth angles simultaneously.

        This is the batch counterpart of ``_transmission_local``.  It computes
        the same refraction physics for a fixed (μ_out, μ_in) pair but all
        azimuth angles φ at once, returning an (n_phi, 4, 4) array.

        Key differences from ``_reflection_vec``
        -----------------------------------------
        1. **Facet normal via vector Snell**: diff = n1 k̂_in − n2 k̂_out; the
           normal points in the direction of this vector (with sign chosen so
           z > 0).  Note the ``flip`` operation that ensures the upward convention.

        2. **TIR guard**: entries where sin²θ_t ≥ 1 are flagged invalid and
           zeroed in the output.

        3. **Solid-angle Jacobian**: the transmission scale additionally includes
           ``jac = n2² μ_t μ_i / (n2 μ_t − n1 μ_i)²``, accounting for the
           solid-angle compression/expansion on refraction.

        4. **Fresnel transmission amplitudes**: the power transmittances T_s and
           T_p are computed from the Fresnel amplitude formulas:
               t_s = 1 + r_s,   t_p = (1 + r_p) · n1/n2
           giving T_s = t_s² and T_p = t_p².  The sign convention for
           b = (T_p − T_s)/2 follows that of ``fresnel_transmission_mueller``.

        Parameters
        ----------
        mu_out : float
            |cos θ_out| of the outgoing (refracted) ray.
        mu_in : float
            |cos θ_in| of the incoming ray.
        phis : np.ndarray, shape (n_phi,)
            Array of relative azimuth angles φ in radians.
        upward : bool
            False = air → water (downward refraction, T_I).
            True  = water → air (upward refraction, T_I*).

        Returns
        -------
        np.ndarray, shape (n_phi, 4, 4)
            Array of 4×4 Mueller matrices, one per azimuth angle.
            Zero for TIR and degenerate geometries.
        """
        n_phi = len(phis)
        mu_out = float(np.clip(mu_out, _EPS, 1.0))
        mu_in  = float(np.clip(mu_in,  _EPS, 1.0))
        sin_in  = float(np.sqrt(max(0.0, 1.0 - mu_in  * mu_in)))
        sin_out = float(np.sqrt(max(0.0, 1.0 - mu_out * mu_out)))

        if upward:
            # Water → air: both rays point upward (z > 0)
            k_in  = np.array([sin_in, 0.0, mu_in], dtype=float)
            k_out = np.column_stack([sin_out * np.cos(phis), sin_out * np.sin(phis),
                                     np.full(n_phi, mu_out)])
            n1, n2 = self.n_water, self.n_air
        else:
            # Air → water: both rays point downward (z < 0)
            k_in  = np.array([sin_in, 0.0, -mu_in], dtype=float)
            k_out = np.column_stack([sin_out * np.cos(phis), sin_out * np.sin(phis),
                                     np.full(n_phi, -mu_out)])
            n1, n2 = self.n_air, self.n_water

        # --- Facet normal via vector Snell: n̂ ∝ n1 k̂_in − n2 k̂_out ---
        diff = n1 * k_in[np.newaxis, :] - n2 * k_out                # (n_phi, 3)
        # Enforce upward-pointing normal (z > 0 convention)
        flip = diff[:, 2] < 0.0
        diff[flip] *= -1.0
        mag   = np.linalg.norm(diff, axis=1, keepdims=True)
        valid = (mag.squeeze() > _EPS)
        normal = np.where(mag > _EPS, diff / (mag + _EPS), 0.0)     # (n_phi, 3)

        mu_n     = normal[:, 2]                                      # cos θ_n
        valid   &= (mu_n > 0.0)
        mu_i_arr = np.abs(np.einsum('i,pi->p', k_in, normal))       # local incidence angle
        valid   &= (mu_i_arr > 0.0)

        # --- Fresnel transmission coefficients ---
        sin_i2 = np.maximum(0.0, 1.0 - mu_i_arr * mu_i_arr)
        sin_t2 = (n1 / n2) ** 2 * sin_i2
        tir    = sin_t2 >= 1.0    # total internal reflection mask
        valid &= ~tir             # TIR entries are invalid (no transmitted ray)
        mu_t   = np.sqrt(np.maximum(0.0, 1.0 - sin_t2))
        valid &= (mu_t > 0.0)
        # Fresnel amplitude reflection coefficients (needed to get transmission amplitudes)
        rs_amp = (n1 * mu_i_arr - n2 * mu_t) / (n1 * mu_i_arr + n2 * mu_t + _EPS)
        rp_amp = (n2 * mu_i_arr - n1 * mu_t) / (n2 * mu_i_arr + n1 * mu_t + _EPS)
        # Power transmission amplitudes: t_s = 1 + r_s; t_p = (1 + r_p) · n1/n2
        ts_amp = 1.0 + rs_amp
        tp_amp = (1.0 + rp_amp) * n1 / (n2 + _EPS)
        ts = ts_amp * ts_amp   # T_s = |t_s|²
        tp = tp_amp * tp_amp   # T_p = |t_p|²

        # --- Scale factors ---
        shadow = smith_sancer_shadowing(mu_out, mu_in, self.sigma)
        sig    = self.sigma
        mu_n_s = np.where(valid, np.clip(mu_n, _EPS, None), 1.0)
        tan2   = (1.0 - mu_n_s * mu_n_s) / (mu_n_s * mu_n_s * sig * sig + _EPS)
        pdf    = np.where(valid, np.exp(-tan2) / (np.pi * sig * sig * mu_n_s ** 3 + _EPS), 0.0)
        # Solid-angle Jacobian for refraction (Zhai et al. 2010, eq. 5)
        jac    = np.where(valid, (n2 * n2 * mu_t * mu_i_arr) / ((n2 * mu_t - n1 * mu_i_arr) ** 2 + _EPS), 0.0)
        scale  = np.where(valid, shadow * np.pi * pdf * jac / (mu_in * mu_out * mu_n_s + _EPS), 0.0)

        # --- Meridian planes and facet plane ---
        z = np.array([0.0, 0.0, 1.0])
        mer_in_v = np.cross(z, k_in)
        m = np.linalg.norm(mer_in_v)
        mer_in_v   = mer_in_v / m if m > _EPS else np.array([0.0, 1.0, 0.0])
        mer_in_rep = np.tile(mer_in_v, (n_phi, 1))

        mer_out_c = np.cross(z[np.newaxis, :], k_out)
        mer_out_m = np.linalg.norm(mer_out_c, axis=1, keepdims=True)
        mer_out   = np.where(mer_out_m > _EPS, mer_out_c / (mer_out_m + _EPS),
                             np.array([0.0, 1.0, 0.0]))
        facet_c = np.cross(normal, k_out)
        facet_m = np.linalg.norm(facet_c, axis=1, keepdims=True)
        facet_p = np.where(facet_m > _EPS, facet_c / (facet_m + _EPS), mer_out)

        k_in_rep = np.tile(k_in, (n_phi, 1))
        xi1 = _signed_angle_vec(mer_in_rep, facet_p, k_in_rep)     # (n_phi,)
        xi2 = _signed_angle_vec(facet_p, mer_out, k_out)            # (n_phi,)

        # --- Vectorised matrix product: scale · L(−ξ₂) · M_trans · L(−ξ₁) ---
        # Transmission Mueller matrix elements
        a_f = 0.5 * (ts + tp)           # average transmittance
        b_f = 0.5 * (tp - ts)           # differential; (T_p−T_s)/2; for air→water b_f > 0
        c_f = np.sqrt(np.maximum(ts * tp, 0.0))

        c1, s1 = np.cos(-2.0 * xi1), np.sin(-2.0 * xi1)
        c2, s2 = np.cos(-2.0 * xi2), np.sin(-2.0 * xi2)

        # Step 1: RL = M_trans · L(−ξ₁) (same sparse structure as reflection)
        RL = np.zeros((n_phi, 4, 4))
        RL[:, 0, 0] = a_f;       RL[:, 0, 1] = b_f * c1;    RL[:, 0, 2] = b_f * s1
        RL[:, 1, 0] = b_f;       RL[:, 1, 1] = a_f * c1;    RL[:, 1, 2] = a_f * s1
        RL[:, 2, 1] = -c_f * s1; RL[:, 2, 2] = c_f * c1;    RL[:, 3, 3] = c_f

        # Step 2: result = L(−ξ₂) · RL — mixes rows 1 and 2 only
        result = RL.copy()
        result[:, 1, :] =  c2[:, np.newaxis] * RL[:, 1, :] + s2[:, np.newaxis] * RL[:, 2, :]
        result[:, 2, :] = -s2[:, np.newaxis] * RL[:, 1, :] + c2[:, np.newaxis] * RL[:, 2, :]
        result *= scale[:, np.newaxis, np.newaxis]   # apply scale prefactor
        result[~valid] = 0.0                         # zero TIR and degenerate entries
        return result


# ---------------------------------------------------------------------------
# Public OceanSurface
# ---------------------------------------------------------------------------


@dataclass
class OceanSurface:
    """Trees & Stam (2019) style complete ocean surface operator for PyMieDAP.

    Physical model
    --------------
    ``OceanSurface`` assembles the full ocean reflectance operator by combining
    three components:

    1. **Rough interface** (``RoughInterface``): specular Fresnel reflection and
       refraction at a Cox-Munk wind-roughened surface.  Produces the air-side
       reflection R_I, downward transmission T_I, water-side reflection R_I*, and
       upward transmission T_I*.

    2. **Water body** (``WaterBody``): polarized multiple-scattering reflection
       operator R_W of the pure-seawater column below the interface, computed by
       adding-doubling.

    3. **Whitecap foam**: a Lambertian, non-polarizing contribution representing
       wind-broken wave foam, weighted by the Monahan fraction f_w.

    These are coupled using the Trees & Stam (2019, A&A 626, A129) formula
    (Appendix A, eq. A.1):

        R_CO = R_I + T_I* · R_W · (I − R_I* · R_W)^{-1} · T_I

    Physical interpretation of each term:
      - R_I           : light reflected at the surface without entering the water
      - T_I* · [...] · T_I : light that enters the water (T_I), bounces around
                            inside and emerges from the water body (R_W term),
                            passes back through the interface (T_I*); the factor
                            (I − R_I* · R_W)^{-1} sums the infinite series of
                            multiple back-reflections between the interface and
                            the water body.

    The final ocean operator is:
        R_ocean = f_w · R_foam + (1 − f_w) · R_CO

    Solver options
    --------------
    ``solver="adding_doubling"`` (default): proper polarized solver implementing
    the full Trees & Stam coupling formula.

    ``solver="diffuse_closure"`` (fast): legacy scalar closure that adds a
    diffuse water-colour albedo on top of the rough-interface reflection, without
    full Stokes coupling.  Retained for backward compatibility and speed tests.

    Rectangular supermatrices
    -------------------------
    The atmosphere and the ocean water body are discretised on independent
    angular grids:
      - The atmosphere uses DAP's ``nmug`` Gauss-Legendre streams.
      - The water body uses ``water_streams`` streams (typically denser).

    The coupling operators T_I (mapping air streams → water streams) and T_I*
    (mapping water streams → air streams) are therefore *rectangular*
    supermatrices, not square.  This allows independent optimisation of the
    angular resolution in each medium.

    Parameters
    ----------
    wind_speed : float
        Isotropic wind speed U_{10} in m/s.  Default 7 m/s matches the baseline
        in Trees & Stam (2019).  Affects σ² and the whitecap fraction.
    wavelength_um : float | None
        Wavelength in microns.  Must be set (or the model will set it via
        ``with_wavelength``) before calling ``fourier_supermatrix``.
    n_air, n_water : float
        Real refractive indices of the air and water sides of the interface.
    foam_albedo : float
        Lambertian whitecap single-scattering albedo.  Default 0.22 is typical
        for fresh foam (Kokhanovsky 2004).
    depth_m : float
        Geometrical water depth in metres (passed to ``WaterBody``).
    bottom_albedo : float
        Lambertian lower boundary of the water column (0 = black, default).
    water_depol : float
        King depolarization factor for the pure-water Rayleigh phase matrix.
    water_streams : int
        Number of Gaussian streams in the water body (independent of DAP).
    water_initial_tau : float
        Starting optical thickness for the water-body doubling initialisation.
    n_fourier : int
        Number of Fourier terms in the ocean surface operators.  Must be large
        enough to resolve the sun-glint specular spike (typically 80–200).
    n_phi : int
        Number of azimuth quadrature samples for the interface Fourier projection.
        Must satisfy n_phi ≥ max(32, 2·n_fourier + 1) to avoid aliasing.
    water_n_phi : int
        Number of azimuth quadrature samples for the water-body phase-matrix
        Fourier projection (separate from n_phi; 240 is typically sufficient).
    include_subsurface : bool
        If False, omit the water body entirely and return only R_I (+ foam).
        Useful for quick surface-glint-only calculations.
    include_foam : bool
        If False, the whitecap fraction is forced to zero (no foam mixing).
    aw, bw : float | None
        User-supplied absorption/scattering coefficients in m⁻¹ for the water
        body.  None = use built-in table.
    solver : {"adding_doubling", "diffuse_closure"}
        Solver algorithm for the water-body / interface coupling.
    """

    wind_speed: float = 7.0
    wavelength_um: float | None = None
    n_air: float = 1.0
    n_water: float = 1.33
    foam_albedo: float = 0.22
    depth_m: float = 100.0
    bottom_albedo: float = 0.0
    water_depol: float = 0.09
    water_streams: int = 16
    water_initial_tau: float = 0.01
    n_fourier: int = 80
    n_phi: int = 720
    water_n_phi: int = 240
    include_subsurface: bool = True
    include_foam: bool = True
    aw: float | None = None
    bw: float | None = None
    solver: Literal["adding_doubling", "diffuse_closure"] = "adding_doubling"

    @property
    def whitecap_fraction(self) -> float:
        """Monahan whitecap fractional area coverage f_w.

        Returns 0 if ``include_foam`` is False, otherwise calls
        ``monahan_whitecap_fraction`` with the current wind speed.

        Returns
        -------
        float
            Whitecap fraction in [0, 1].
        """
        if not self.include_foam:
            return 0.0
        return monahan_whitecap_fraction(self.wind_speed)

    @property
    def slope_variance(self) -> float:
        """Cox-Munk isotropic wave-slope variance σ² for the current wind speed.

        Returns
        -------
        float
            Slope variance σ² (dimensionless).
        """
        return cox_munk_slope_variance(self.wind_speed)

    @property
    def subsurface_albedo(self) -> float:
        """Legacy scalar water-column albedo (unpolarised diffuse closure).

        Used only by the ``diffuse_closure`` solver path and the diagnostic
        ``reflection_matrix`` method.  Returns 0 when ``include_subsurface``
        is False or wavelength is not yet set.

        Returns
        -------
        float
            Approximate hemispheric water-leaving + bottom albedo in [0, 1].
        """
        if not self.include_subsurface or self.wavelength_um is None:
            return 0.0
        return pure_water_diffuse_albedo(
            self.wavelength_um,
            depth_m=self.depth_m,
            bottom_albedo=self.bottom_albedo,
            aw=self.aw,
            bw=self.bw,
        )

    def with_wavelength(self, wavelength_um: float) -> "OceanSurface":
        """Return a copy of this ``OceanSurface`` with a concrete wavelength assigned.

        Because ``OceanSurface`` is a frozen dataclass, this creates a new instance
        via ``dataclasses.replace``.  Useful for building wavelength-specific surface
        operators inside a wavelength loop without mutating the original object.

        Parameters
        ----------
        wavelength_um : float
            Wavelength in microns.

        Returns
        -------
        OceanSurface
            New instance identical to ``self`` except for ``wavelength_um``.
        """
        return replace(self, wavelength_um=float(wavelength_um))

    def _interface(self) -> RoughInterface:
        """Construct a ``RoughInterface`` object from the current surface parameters.

        Returns
        -------
        RoughInterface
            Configured with this surface's wind speed and refractive indices.
        """
        return RoughInterface(wind_speed=self.wind_speed, n_air=self.n_air, n_water=self.n_water)

    def _water_body(self) -> WaterBody:
        """Construct a ``WaterBody`` object from the current surface parameters.

        Raises
        ------
        ValueError
            If ``wavelength_um`` is None (must be set before calling this).

        Returns
        -------
        WaterBody
            Configured with this surface's depth, optical, and quadrature parameters.
        """
        if self.wavelength_um is None:
            raise ValueError("OceanSurface.wavelength_um must be set before building the water body")
        return WaterBody(
            wavelength_um=float(self.wavelength_um),
            depth_m=self.depth_m,
            bottom_albedo=self.bottom_albedo,
            depolarization=self.water_depol,
            n_streams=self.water_streams,
            initial_tau=self.water_initial_tau,
            n_phi=self.water_n_phi,
            aw=self.aw,
            bw=self.bw,
        )

    def reflection_matrix(self, mu: float, mu0: float, dphi: float) -> np.ndarray:
        """Diagnostic local ocean reflection Mueller matrix for a single angle triple.

        Computes the ocean reflectance for viewing direction (μ, φ₀+Δφ) and
        solar direction (μ₀, φ₀=0) using:
          - The rough-interface reflection R_I (always included).
          - An optional scalar diffuse water-colour contribution (unpolarised)
            using the legacy ``pure_water_diffuse_albedo`` closure.
          - Whitecap mixing using the Monahan fraction.

        This is a diagnostic tool only.  The proper polarized water coupling
        (Trees & Stam formula) is implemented in ``fourier_supermatrix`` via
        ``_adding_doubling_supermatrix``.

        Parameters
        ----------
        mu : float
            |cos θ| of the viewing direction.
        mu0 : float
            |cos θ₀| of the solar direction.
        dphi : float
            Relative azimuth Δφ = φ_view − φ_sun in radians.

        Returns
        -------
        np.ndarray, shape (4, 4)
            Local ocean reflection Mueller matrix (meridian-plane frame).
            Note: only [0,0] includes the water-colour contribution; the full
            polarized water-body coupling is only in ``fourier_supermatrix``.
        """
        interface = self._interface()
        clean = interface.reflection_air(mu, mu0, dphi)

        if self.include_subsurface and self.wavelength_um is not None:
            water = np.zeros((4, 4), dtype=float)
            a_sub = self.subsurface_albedo
            if a_sub > 0.0:
                # Approximate water contribution: intensity only, two-way Fresnel transmittance
                t_down = 1.0 - fresnel_unpolarized_reflectance(mu0, self.n_air, self.n_water)
                t_up   = 1.0 - fresnel_unpolarized_reflectance(mu,  self.n_air, self.n_water)
                water[0, 0] = a_sub * t_down * t_up
            clean = clean + water

        # Whitecap foam: Lambertian, non-polarizing (only I component)
        foam = np.zeros((4, 4), dtype=float)
        foam[0, 0] = self.foam_albedo
        q = float(np.clip(self.whitecap_fraction, 0.0, 1.0))
        # Linear mixture: f_w · foam + (1 − f_w) · clean_ocean
        return q * foam + (1.0 - q) * clean

    def _foam_supermatrix(self, mus: np.ndarray, smf: np.ndarray, nmat: int, n_fourier: int) -> np.ndarray:
        """Build the Lambertian whitecap-foam Fourier supermatrix.

        Whitecap foam is modelled as an isotropic (Lambertian) non-polarizing
        reflector with constant albedo ``foam_albedo``.  Its Mueller matrix is
        simply ``foam_albedo · diag(1, 0, 0, 0)`` (only I component; no Q, U, V).

        In the Fourier supermatrix representation, a Lambertian surface
        contributes only to the m=0 Fourier term (azimuth-independent), and only
        to the I–I block of the supermatrix:

            foam[i·nmat, j·nmat, m=0] = foam_albedo · smf[i] · smf[j]

        All other entries are zero.

        Parameters
        ----------
        mus : np.ndarray
            Quadrature cosines for the air streams.
        smf : np.ndarray
            Supermatrix scaling factors.
        nmat : int
            Number of Stokes components.
        n_fourier : int
            Highest Fourier order (foam contributes only to m=0).

        Returns
        -------
        np.ndarray, shape (nsup, nsup, n_fourier+1), Fortran order
            Foam reflection supermatrix.
        """
        nsup = len(mus) * nmat
        foam = np.zeros((nsup, nsup, n_fourier + 1), dtype=np.float64, order="F")
        idx = np.arange(len(mus)) * nmat   # indices of the I-component of each stream
        # Index the m=0 slice (a 2-D view), then fill the I–I block
        foam[:, :, 0][np.ix_(idx, idx)] = float(self.foam_albedo) * np.outer(smf, smf)
        return foam

    def _diffuse_closure_supermatrix(self, mus: np.ndarray, smf: np.ndarray, nmat: int, n_fourier: int) -> np.ndarray:
        """Build the rough-interface + scalar-water-colour supermatrix (legacy solver).

        This is the legacy ``diffuse_closure`` solver path.  It computes:
          1. The polarized rough-interface reflection R_I as a Fourier supermatrix.
          2. Optionally adds a scalar diffuse water-colour contribution to the
             m=0 I–I block (using ``pure_water_diffuse_albedo``).

        The water contribution is intensity-only and uses unpolarised Fresnel
        transmittances.  No polarization coupling between the interface and the
        water is modelled.

        Parameters
        ----------
        mus : np.ndarray
            Quadrature cosines for the air streams.
        smf : np.ndarray
            Supermatrix scaling factors.
        nmat : int
            Number of Stokes components.
        n_fourier : int
            Highest Fourier order.

        Returns
        -------
        np.ndarray, shape (nsup, nsup, n_fourier+1), Fortran order
            Clean-ocean supermatrix (interface only + optional diffuse water).
        """
        interface = self._interface()
        # Compute the full Fourier-decomposed rough-interface reflection supermatrix R_I
        rough = _project_operator_fourier(
            mus, smf, mus, smf, nmat, n_fourier, self.n_phi,
            interface.reflection_air,
            local_matrix_vec=lambda mo, mi, ps: interface._reflection_vec(mo, mi, ps, from_below=False),
        )
        if self.include_subsurface and self.wavelength_um is not None:
            a_sub = self.subsurface_albedo   # scalar diffuse water-column albedo
            if a_sub > 0.0:
                # For each stream, compute the Fresnel transmittance into the water.
                # This is 1 − R_Fresnel(μ_i) for unpolarised light.
                t_vec = np.array([
                    1.0 - fresnel_unpolarized_reflectance(mu, self.n_air, self.n_water)
                    for mu in mus
                ])                                           # (n_str,)
                # Add water-colour contribution to the m=0 I–I block:
                # ΔR[i·nmat, j·nmat, 0] += a_sub · T_i · smf[i] · smf[j] · T_j
                idx = np.arange(len(mus)) * nmat
                rough[:, :, 0][np.ix_(idx, idx)] += (
                    a_sub * np.outer(smf * t_vec, smf * t_vec)
                )
        return rough

    def _adding_doubling_supermatrix(self, air_mus: np.ndarray, air_smf: np.ndarray, nmat: int, n_fourier: int) -> np.ndarray:
        """Build the full polarized clean-ocean supermatrix via the Trees & Stam coupling formula.

        This implements eq. (A.1) of Trees & Stam (2019):

            R_CO^(m) = R_I^(m) + T_I*^(m) · R_W^(m) · (I − R_I*^(m) · R_W^(m))^{-1} · T_I^(m)

        independently for each Fourier order m = 0, 1, …, n_fourier.

        Step-by-step
        ------------
        1. Compute R_I (air–air reflection supermatrix) using
           ``_project_operator_fourier`` with the vectorised reflection kernel.

        2. If ``include_subsurface`` is False, return R_I immediately.

        3. Compute R_W (water-body reflection supermatrix) via
           ``WaterBody.reflection_supermatrices``.

        4. Compute the four rectangular interface operators:
           - T_I    : (n_water×nmat, n_air×nmat)  — air streams → water streams
           - R_I*   : (n_water×nmat, n_water×nmat) — water streams → water (TIR)
           - T_I*   : (n_air×nmat,   n_water×nmat) — water streams → air streams

        5. For each Fourier order m, solve the coupling equation:
               D = (I − R_I*^(m) · R_W^(m))^{-1} · T_I^(m)     [np.linalg.solve for stability]
               U = R_W^(m) · D
               R_water_leaving = T_I*^(m) · U
               R_CO^(m) = R_I^(m) + R_water_leaving

        Note on solve vs inv
        --------------------
        ``np.linalg.solve(A, B)`` solves AX = B without forming A^{-1} explicitly.
        This is more numerically stable and faster than ``np.linalg.inv(A) @ B``
        for the typically well-conditioned (I − R_I* R_W) systems encountered here.

        Parameters
        ----------
        air_mus : np.ndarray
            Quadrature cosines of the DAP atmosphere streams.
        air_smf : np.ndarray
            Supermatrix scaling factors for the atmosphere streams.
        nmat : int
            Number of Stokes components.
        n_fourier : int
            Highest Fourier order.

        Returns
        -------
        np.ndarray, shape (nsup_air, nsup_air, n_fourier+1), Fortran order
            Complete clean-ocean reflection supermatrix R_CO.
        """
        interface = self._interface()
        # Step 1: Rough-interface air-side reflection supermatrix R_I
        # Shape: (n_air·nmat, n_air·nmat, n_fourier+1)
        R_air = _project_operator_fourier(
            air_mus, air_smf, air_mus, air_smf, nmat, n_fourier, self.n_phi,
            interface.reflection_air,
            local_matrix_vec=lambda mo, mi, ps: interface._reflection_vec(mo, mi, ps, from_below=False),
        )

        if not self.include_subsurface:
            # Surface-only mode: return just the rough-interface reflection.
            return R_air

        # Step 2: Water-body reflection supermatrix R_W on the water stream grid
        # Shape: (n_water·nmat, n_water·nmat, n_fourier+1)
        water_body = self._water_body()
        Rw, w_mus, _w_weights, w_smf = water_body.reflection_supermatrices(nmat=nmat, n_fourier=n_fourier)

        # Step 3: Rectangular interface coupling operators
        # T_I = downward transmission: maps air streams → water streams
        # Shape: (n_water·nmat, n_air·nmat, n_fourier+1)
        T_down = _project_operator_fourier(
            w_mus, w_smf, air_mus, air_smf, nmat, n_fourier, self.n_phi,
            interface.transmission_down,
            local_matrix_vec=lambda mo, mi, ps: interface._transmission_vec(mo, mi, ps, upward=False),
        )
        # R_I* = water-side interface reflection: maps water streams → water streams
        # Shape: (n_water·nmat, n_water·nmat, n_fourier+1)
        R_int_star = _project_operator_fourier(
            w_mus, w_smf, w_mus, w_smf, nmat, n_fourier, self.n_phi,
            interface.reflection_water,
            local_matrix_vec=lambda mo, mi, ps: interface._reflection_vec(mo, mi, ps, from_below=True),
        )
        # T_I* = upward transmission: maps water streams → air streams
        # Shape: (n_air·nmat, n_water·nmat, n_fourier+1)
        T_up = _project_operator_fourier(
            air_mus, air_smf, w_mus, w_smf, nmat, n_fourier, self.n_phi,
            interface.transmission_up,
            local_matrix_vec=lambda mo, mi, ps: interface._transmission_vec(mo, mi, ps, upward=True),
        )

        # Step 4: Trees & Stam coupling for each Fourier order
        nsup_air = len(air_mus) * nmat
        clean = np.zeros((nsup_air, nsup_air, n_fourier + 1), dtype=np.float64, order="F")
        for m in range(n_fourier + 1):
            if np.max(np.abs(Rw[:, :, m])) <= 0.0:
                # Water body has no reflectance at this Fourier order
                # (or at all, for a zero-ω₀ water body); skip the coupling.
                clean[:, :, m] = R_air[:, :, m]
                continue
            # Q1 = R_I*^(m) · R_W^(m)
            Q1 = R_int_star[:, :, m] @ Rw[:, :, m]
            Iw = np.eye(Q1.shape[0])
            # D = (I − R_I*^(m) · R_W^(m))^{-1} · T_I^(m)
            # Using solve instead of inv for numerical stability.
            D = np.linalg.solve(Iw - Q1, T_down[:, :, m])
            # U = R_W^(m) · D  (upwelling water radiance after round-trip)
            U = Rw[:, :, m] @ D
            # Water-leaving contribution: T_I*^(m) · U
            R_water_leaving = T_up[:, :, m] @ U
            # Total clean-ocean: surface glint + water-leaving
            clean[:, :, m] = R_air[:, :, m] + R_water_leaving
        return clean

    def fourier_supermatrix(
        self,
        nmug: int,
        nmat: int = 4,
        nmuMAX: int = 201,
        nmatMAX: int = 4,
    ) -> Tuple[np.ndarray, int, int, np.ndarray]:
        """Build and return DAP-ready Fourier supermatrices for the ocean surface.

        This is the main entry point called by the PyMieDAP model setup code
        (``pmd.compute`` via the surface-check hook).  It performs all the
        computationally expensive steps: Fourier projection of the rough
        interface, water-body adding-doubling, Trees & Stam coupling, and
        whitecap mixing.

        The returned ``surfmat`` array can be passed directly to the patched
        Fortran DAP surface routine as a bottom boundary condition.

        Parameters
        ----------
        nmug : int
            Number of Gauss-Legendre atmospheric streams per hemisphere (the
            DAP stream count).  The quadrature includes an extra μ=1 stream,
            so the actual grid has nmug+1 nodes.
        nmat : int
            Number of Stokes components (1, 3, or 4).
        nmuMAX : int
            Maximum allowed number of angular streams (Fortran array bound).
            Must satisfy nmug+1 ≤ nmuMAX.
        nmatMAX : int
            Maximum allowed number of Stokes components (Fortran array bound).
            Must satisfy nmat ≤ nmatMAX.

        Returns
        -------
        surfmat : np.ndarray, shape (nsup, nsup, nsurfou+1), Fortran order
            Full ocean surface Fourier supermatrix including glint, water body,
            and foam contributions.
        nsup : int
            Size of the supermatrix: (nmug+1) · nmat.
        nsurfou : int
            Number of Fourier orders = ``self.n_fourier``.
        mus : np.ndarray, shape (nmug+1,)
            Quadrature cosines of the DAP stream grid (with the extra μ=1 node).

        Raises
        ------
        ValueError
            If any of the input size parameters fail their bounds checks, or if
            ``n_phi < max(32, 2·n_fourier + 1)`` (would cause Fourier aliasing).
        """
        if nmat not in (1, 3, 4):
            raise ValueError("nmat must be 1, 3 or 4")
        if nmug < 1:
            raise ValueError("nmug must be positive")
        if self.n_fourier < 0:
            raise ValueError("n_fourier must be non-negative")
        if self.n_phi < max(32, 2 * self.n_fourier + 1):
            raise ValueError("n_phi should be at least max(32, 2*n_fourier + 1)")
        if nmug + 1 > nmuMAX:
            raise ValueError("nmug+1 exceeds nmuMAX")
        if nmat > nmatMAX:
            raise ValueError("nmat exceeds nmatMAX")

        # Build the DAP quadrature grid (nmug Gauss-Legendre nodes + μ=1 extra)
        air_mus, _air_weights, air_smf = _quadrature(nmug, include_extra_mu_one=True)
        nsup    = len(air_mus) * nmat   # total supermatrix dimension
        nsurfou = int(self.n_fourier)   # number of Fourier orders

        # Choose solver: full polarized coupling or legacy diffuse closure
        if self.solver == "diffuse_closure":
            clean = self._diffuse_closure_supermatrix(air_mus, air_smf, nmat, nsurfou)
        elif self.solver == "adding_doubling":
            clean = self._adding_doubling_supermatrix(air_mus, air_smf, nmat, nsurfou)
        else:
            raise ValueError("solver must be 'adding_doubling' or 'diffuse_closure'")

        # Mix clean ocean with whitecap foam using the Monahan fraction f_w
        q    = float(np.clip(self.whitecap_fraction, 0.0, 1.0))
        foam = self._foam_supermatrix(air_mus, air_smf, nmat, nsurfou)
        # R_ocean = f_w · R_foam + (1 − f_w) · R_CO
        surfmat = q * foam + (1.0 - q) * clean
        return np.asfortranarray(surfmat), nsup, nsurfou, air_mus
