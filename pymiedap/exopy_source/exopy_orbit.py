# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

# -*- coding: utf-8 -*-
#"""
#==================================================================
#EXOPY module: exopy_orbit.py
#Delft University of Technology
#------------------------------------------------------------------
#Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
#Date: 2016-2017
#------------------------------------------------------------------
#
#Dependences:
#
#DESCRIPTION
#------------------------------------------------------------------
#Script containing the functions required for the computation of
#the orbits of the bodies under the assumption of a nested two
#body problem.
#
#LIST OF FUNCTIONS
#------------------------------------------------------------------
# - kepler: Calculates the eccentric anomaly given a value of mean
#	   anomaly and eccentricity.
# - nested2bp: Computes the orbit of each body under the assumption
#	      of a nested two body problem.
#
#REFERENCES
#------------------------------------------------------------------
#  [1] Kipping, D. M. (2011). luna: an algorithm for generating dynamic
#      planetary moon transits. Monthly Notices of the Royal Astronomical
#      Society, 416(1), 689-709.
#  [2] Wakker, K. F. (2015). Fundamentals of astrodynamics. TU Delft Library.
#
#"""

# Required modules
import math, numpy as np
import exopy_functions as fun

def kepler(M, e, method = 'Newton-Raphson'):
    """
    ==================================================================
    EXOPY function: kepler()
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    Dependences:

    DESCRIPTION
    ------------------------------------------------------------------
    Calculates the Eccentric anomaly given a value of Mean anomaly and
    eccentricity, based on [1].

    INPUTS
    ------------------------------------------------------------------
    - M: Mean anomaly [rad=0:2pi] (float)
    - e: Eccentricity [-] (float)
    - method: 'Newton-Raphson' or 'Merkley'

    OUTPUTS
    ------------------------------------------------------------------
    - E: Eccentric anomaly [rad] (float)

    REFERENCES
    ------------------------------------------------------------------
    [1] Wakker, K. F. (2015). Fundamentals of astrodynamics. TU Delft Library.

    """
    if method == 'Merkley':

        alpha = 3*math.pi**2 +1.6*math.pi*(math.pi-abs(M))/(1+e)/(math.pi**2-6)
        d = 3*(1-e) + alpha*e
        r = 3*alpha*d*(d-1+e)*M+M**3
        q = 2*alpha*d*(1-e)-M**2
        w = (abs(r)+(q**3+r**2)**0.5)**(2/3)

        E = (2*r*w/(w**2+w*q+q**2)+M)/d

        f  = E - e*math.sin(E) - M
        f1 = 1 - e*math.cos(E)
        f2 = e*math.sin(E)
        f3 = 1 - f1
        f4 = -f2

        d3 = -f/(f1-0.5*f*f2/f1)
        d4 = -f/(f1+0.5*d3*f2+d3*d3*f3/6)
        d5 = -f/(f1+0.5*d4*f2+f3*d4*d4/6+d4*d4*d4*f4/24)

        E = E+d5

    elif method == 'Newton-Raphson':
        E = M
        error = 1
        while error > 1e-10:

            E_prev = E
            E = E - (E-e*math.sin(E)-M)/(1-e*math.cos(E))
            error = abs(E-E_prev)

    return E




def nested2bp(moon, planet, star, dt, tf):
    """
    ==================================================================
    EXOPY function: nested2bp()
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    Dependences:

    DESCRIPTION
    ------------------------------------------------------------------
    Computes the orbit of each body under the assumption of a nested
    two body problem, based on [1].

    INPUTS
    ------------------------------------------------------------------
    - moon: body object of type 'moon' [-] ('body' object)
    - planet: body object of type 'planet' [-] ('body' object)
    - star: body object of type 'star' [-] ('body' object)
    - dt: time step [s] (float)
    - tf: final time to be computed [s] (float)

    OUTPUTS
    ------------------------------------------------------------------
    - moon: updated body object of type 'moon' [-] ('body' object)
    - planet: updated body object of type 'planet' [-] ('body' object)
    - star: updated body object of type 'star' [-] ('body' object)

    REFERENCES
    ------------------------------------------------------------------
    [1] Kipping, D. M. (2011). luna: an algorithm for generating
        dynamic planetary moon transits. Monthly Notices of the Royal
        Astronomical Society, 416(1), 689-709.


    """

    #==============================================================================
    #
    #  ZERO STEP: Load the problem parameters from the input data
    #
    #==============================================================================
    G = 6.674e-11                              # Universal gravitational constant
                                               # [N⋅m2/kg2]
    m_m      = moon.properties.m               # Mass of moon   [kg]
    m_p      = planet.properties.m             # Mass of planet [kg]
    m_s      = star.properties.m               # Mass of star   [kg]

    a_mb     = moon.orbital_elements.a         # Moon-baryc. semimajor axis [m]
    e_mb     = moon.orbital_elements.e         # Moon orbit eccentricity [-]
    i_mb     = moon.orbital_elements.i         # Moon orbit inclination [deg]
    Omega_mb = moon.orbital_elements.Omega     # Moon RAAN [deg]
    omega_mb = moon.orbital_elements.omega     # Argument of moon periapsis [deg]
    t0_m     = moon.orbital_elements.t0	    # Lunar time of last per. passage [s]

    a_bs     = planet.orbital_elements.a_b     # Barycenter-star semim.axis [m]
    e_bs     = planet.orbital_elements.e_b     # Barycenter eccentricity [-]
    i_bs     = planet.orbital_elements.i_b     # Barycenter inclination [deg]
    Omega_bs = planet.orbital_elements.Omega_b # Barycenter RAAN [deg]
    omega_bs = planet.orbital_elements.omega_b # Argument of bar. periapsis [deg]
    t0_b     = planet.orbital_elements.t0_b    # Bar. time of last per. passage [s]

    #==============================================================================
    #
    #  FIRST STEP: Create time vector and true anomaly vector
    #
    #==============================================================================

    mu_m   = G*(m_p + m_m)       # Moon-planet standard gravitational parameter
    mu_b   = G*(m_m + m_p + m_s) # Moon-planet-barycenter stnd. grav. parameter
                                 # [N⋅m2/kg]
    time   = np.arange(0, tf, dt, dtype=np.float)  # Time vector [s]

    n_mb   = (mu_m/(a_mb*(m_p+m_m)/m_p)**3)**0.5 # Moon-bar. mean motion [1/s]
    n_bs   = (mu_b/a_bs**3)**0.5   # Barycenter-star mean motion [1/s]

    if t0_b is None and t0_m is None:

        f0_m     = moon.orbital_elements.f0_m
        f0_b     = planet.orbital_elements.f0_b

        E0m  = 2* np.atan( np.sqrt((1-e_mb)/(1+e_mb)) * np.tan(f0_m/2) )
        E0b  = 2* np.atan( np.sqrt((1-e_bs)/(1+e_bs)) * np.tan(f0_b/2) )

        M0m  = E0m - e_mb*np.sin(E0m)
        M0b  = E0b - e_bs*np.sin(E0b)

        t0_m = -M0m/n_mb
        t0_b = -M0b/n_bs

    M_mb   = n_mb * (time - t0_m)  # Moon-barycenter mean anomaly [rad]
    M_bs   = n_bs * (time - t0_b)  # Barycenter-star mean anomaly [rad]

    # Eccentric anomaly and true anomaly are initialized
    E_mb   = np.zeros(len(time))
    E_bs   = np.zeros(len(time))
    f_mb   = np.zeros(len(time))
    f_bs   = np.zeros(len(time))

    for index in range(len(E_mb)):

        # E is calculated through the Kepler equation. Then, f can be computed
        # The fun.symp function is used in order to reduce the angular range
        # from 0 to 2pi
        E_mb[index] = kepler(fun.symp(M_mb[index]), e_mb, method = 'Newton-Raphson')
        E_bs[index] = kepler(fun.symp(M_bs[index]), e_bs, method = 'Newton-Raphson')
        f_mb[index]   = math.atan2((math.sin(E_mb[index])*(1-e_mb**2)**0.5)/
                        (1-e_mb*math.cos(E_mb[index])),((math.cos(E_mb[index])
                        -e_mb)/(1-e_mb*math.cos(E_mb[index]))))
        f_bs[index]   = math.atan2((math.sin(E_bs[index])*(1-e_bs**2)**0.5)/
                        (1-e_bs*math.cos(E_bs[index])),((math.cos(E_bs[index])
                        -e_bs)/(1-e_bs*math.cos(E_bs[index]))))
        f_mb[index] = fun.symp(f_mb[index])
        f_bs[index] = fun.symp(f_bs[index])

    #==============================================================================
    #
    #  SECOND STEP: Keplerian orbit of moon (s) and planet (p) around
    #               the moon-planet system barycenter (b).
    #
    #==============================================================================

    # Conversion of angles from degrees to radians

    i_mb     = i_mb     * math.pi / 180 # [rad]
    Omega_mb = Omega_mb * math.pi / 180 # [rad]
    omega_mb = omega_mb * math.pi / 180 # [rad]

    # Orbit equation for the moon around the planet

    r_mb = a_mb*(1-math.pow(e_mb,2))/(1+e_mb*np.cos(f_mb)) # [meters]

    # Two-dimensional cartesian coordinates of the moon motion around the
    # barycenter <position2D_mb = (x_mb, y_mb, z_mb=0)>

    x_mb    = r_mb * np.cos(f_mb)    # [meters]
    y_mb    = r_mb * np.sin(f_mb)    # [meters]

    position2D_mb    = np.zeros((3,len(x_mb)))  # [meters]
    position2D_mb[0] = x_mb                     # [meters]
    position2D_mb[1] = y_mb                     # [meters]

    # Two-dimensional cartesian coordinates of the planet motion around the
    # barycenter, and planet and moon motion around each other.
    # <position2D_mp = (x_mp, y_mp, z_mp=0)>
    # <position2D_pm = (x_pm, y_pm, z_pm=0)>
    # <position2D_pb = (x_pb, y_pb, z_pb=0)>
    #
    # Being the origin of the reference frame the center of mass of the system,
    # the following relation is met:
    # <m_s·position2D_mb + m_p·position2D_pb = 0>
    # while <position2D_mp = position2D_mb - position2D_pb>

    position2D_pb = - position2D_mb * m_m / m_p          # [meters]
    r_pb = np.linalg.norm(position2D_pb, axis = 0)       # [meters]

    position2D_mp =   position2D_mb - position2D_pb      # [meters]
    position2D_pm =  -position2D_mb + position2D_pb      # [meters]

    r_mp = np.linalg.norm(position2D_mp, axis = 0)       # [meters]
    r_pm = np.linalg.norm(position2D_pm, axis = 0)       # [meters]

    # Three-dimensional cartesian coordinates of the moon motion around the
    # planet <position3D_mb = (X_mb, Y_mb, Z_mb)>

    # The effect of orbit inclination, longitude of the ascending node and the
    # argument of periapsis are introduced by mean of three rotation matrices.

    Mz_Omega_mb = np.array([[ np.cos(Omega_mb), np.sin(Omega_mb),   0   ],
                            [-np.sin(Omega_mb), np.cos(Omega_mb),   0   ],
                            [        0        ,        0        ,   1   ]]).T

    Mz_omega_mb = np.array([[ np.cos(omega_mb), np.sin(omega_mb),   0   ],
                            [-np.sin(omega_mb), np.cos(omega_mb),   0   ],
                            [        0        ,        0        ,   1   ]]).T

    Mx_i_mb     = np.array([[      1      ,        0        ,      0      ],
                            [      0      ,  np.cos(i_mb)   , np.sin(i_mb)],
                            [      0      , -np.sin(i_mb)   , np.cos(i_mb)]]).T

    position3D_mb =Mz_Omega_mb.dot(Mx_i_mb).dot(Mz_omega_mb).dot(position2D_mb)

    # Three-dimensional cartesian coordinates of the planet motion around the
    # barycenter, and planet and moon motion around each other.
    # <position3D_mp = (X_mp, Y_mp, Z_mp)>
    # <position3D_pm = (X_pm, Y_pm, Z_pm)>
    # <position3D_pb = (X_pb, Y_pb, Z_pb)>
    #
    # Being the origin of the reference frame the center of mass of the system,
    # the following relation is met:
    # <m_s·position3D_mb + m_p·position3D_pb = 0>
    # while <position3D_mp = position3D_mb - position3D_pb>

    position3D_pb =  -position3D_mb * m_m / m_p          # [meters]
    position3D_mp =   position3D_mb - position3D_pb      # [meters]
    position3D_pm =  -position3D_mb + position3D_pb      # [meters]


    #==============================================================================
    #
    #  THIRD STEP: Keplerian orbit of the moon-planet system barycenter (b) around
    #              the hosting star (s).
    #
    #==============================================================================

    # Conversion of angles from degrees to radians

    i_bs     = i_bs     * math.pi / 180 # [rad]
    Omega_bs = Omega_bs * math.pi / 180 # [rad]
    omega_bs = omega_bs * math.pi / 180 # [rad]

    # Orbit equation for the barycenter

    r_bs = a_bs*(1-math.pow(e_bs,2))/(1+e_bs*np.cos(f_bs)) # [meters]

    # 2-dimensional cartesian coordinates of the moon-planet barycenter motion
    # around the hosting star <position2D_bs = (x_bs, y_bs, z_bs=0)>

    x_bs    = r_bs * np.cos(f_bs)    # [meters]
    y_bs    = r_bs * np.sin(f_bs)    # [meters]

    position2D_bs    = np.zeros((3,len(x_bs)))  # [meters]
    position2D_bs[0] = x_bs                     # [meters]
    position2D_bs[1] = y_bs                     # [meters]


    # 3-dimensional cartesian coordinates of the moon-planet barycenter around
    # the star <position3D_bs = (X_bs, Y_bs, Z_bs)>
    #
    # The effect of orbit inclination, longitude of the ascending node and the
    # argument of periapsis are introduced by mean of three rotation matrices.

    Mz_Omega_bs = np.array([[ np.cos(Omega_bs), np.sin(Omega_bs),   0   ],
                            [-np.sin(Omega_bs), np.cos(Omega_bs),   0   ],
                            [        0        ,        0        ,   1   ]]).T

    Mz_omega_bs = np.array([[ np.cos(omega_bs), np.sin(omega_bs),   0   ],
                            [-np.sin(omega_bs), np.cos(omega_bs),   0   ],
                            [        0        ,        0        ,   1   ]]).T

    Mx_i_bs     = np.array([[     1     ,        0        ,      0      ],
                            [     0     ,  np.cos(i_bs)   , np.sin(i_bs)],
                            [     0     , -np.sin(i_bs)   , np.cos(i_bs)]]).T

    position3D_bs =Mz_Omega_bs.dot(Mx_i_bs).dot(Mz_omega_bs).dot(position2D_bs)

    # Three-dimensional motion of moon and planet with respect to the hosting
    # star

    position3D_ms = Mz_Omega_bs.dot(Mx_i_bs).dot(Mz_omega_bs).dot(
                                                   position2D_bs+position3D_mb)
    position3D_ps = Mz_Omega_bs.dot(Mx_i_bs).dot(Mz_omega_bs).dot(
                                                   position2D_bs+position3D_pb)

    r_ms = np.linalg.norm(position3D_ms, axis = 0)     # [meters]
    r_ps = np.linalg.norm(position3D_ps, axis = 0)     # [meters]


    # The angle between planet and moon with vertex at the star is calculated
    # to have an indication of the proximityof eclipse events.

    b2b_angle = np.rad2deg(np.arccos(np.round(np.diagonal(np.dot(
                position3D_ps.T,position3D_ms))/np.linalg.norm(
                position3D_ps, axis = 0)/np.linalg.norm(
                position3D_ms,axis=0),6)))


    # The angle between the star and the third body with vertex at the body
    # under observation is calculated to have an indication of the proximity
    # of eclipse events.

    aux1 = position3D_ps
    aux2 = position3D_ms - position3D_ps

    beta_p = np.rad2deg(np.arccos((np.diagonal(np.dot(
                -aux1.T,aux2))/np.linalg.norm(
                aux2, axis = 0)/np.linalg.norm(
                aux1,axis=0))))

    aux1 = position3D_ms
    aux2 = position3D_ps - position3D_ms

    beta_m = np.rad2deg(np.arccos((np.diagonal(np.dot(
                -aux1.T,aux2))/np.linalg.norm(
                aux1, axis = 0)/np.linalg.norm(
                aux2,axis=0))))

    # Export results to body objects

    moon.ephemeris.time           = time
    moon.orbital_elements.nu_mb   = f_mb
    moon.orbital_elements.M_mb    = M_mb
    moon.orbital_elements.E_mb    = E_mb
    moon.ephemeris.r_b            = r_mb
    moon.ephemeris.r_p            = r_mp
    moon.ephemeris.r_s            = r_ms
    moon.ephemeris.position2D_mp  = position2D_mp
    moon.ephemeris.position3D_mp  = position3D_mp
    moon.ephemeris.position2D_mb  = position2D_mb
    moon.ephemeris.position3D_mb  = position3D_mb
    moon.ephemeris.position3D_s   = position3D_ms
    moon.ephemeris.period         = 2*math.pi/(n_mb)
    moon.geometry.b2b_angle       = b2b_angle
    moon.geometry.beta            = beta_m

    planet.ephemeris.time         = time
    planet.orbital_elements.nu_pm = f_mb
    planet.orbital_elements.nu_bm = f_mb
    planet.orbital_elements.M_pb  = M_mb
    planet.orbital_elements.E_pb  = E_mb
    planet.ephemeris.r_b          = r_pb
    planet.ephemeris.r_s          = r_ps
    planet.ephemeris.r_m          = r_pm
    planet.ephemeris.position2D_pb= position2D_pb
    planet.ephemeris.position2D_pm= position2D_pm
    planet.ephemeris.position3D_pb= position3D_pb
    planet.ephemeris.position3D_s = position3D_ps
    planet.ephemeris.position3D_pm= position3D_pm
    planet.ephemeris.period       = 2*math.pi/(n_mb)
    planet.ephemeris.period_bs    = 2*math.pi/(n_bs)
    planet.geometry.b2b_angle     = b2b_angle
    planet.geometry.beta          = beta_p

    star.ephemeris.time           = time
    star.orbital_elements.nu_bs   = f_bs
    star.orbital_elements.M_bs    = M_bs
    star.orbital_elements.E_bs    = E_bs
    star.ephemeris.r_b            = r_bs
    star.ephemeris.r_m            = r_ms
    star.ephemeris.r_p            = r_ps
    star.ephemeris.position2D_bs  = position2D_bs
    star.ephemeris.position3D_bs  = position3D_bs
    star.ephemeris.period_bs      = 2*math.pi/(n_bs)
    star.ephemeris.position3D_s   = np.zeros([3,len(time)])

    print('    ✓ The trajectories of ' + moon.name + ', ' + planet.name + ', and ' + star.name + ' have been calculated.\n')

    return moon, planet, star


def kepler_orbit(planet, star, dt, tf):
    """
    ==================================================================
        EXOPY function: kepler_orbit()
        Delft University of Technology
        ------------------------------------------------------------------
        Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
        Date: 2016-2017
        ------------------------------------------------------------------
        Dependences:
        DESCRIPTION
        ------------------------------------------------------------------
        Computes the orbit of the planet under the assumption of a Kepler
        problem.
        INPUTS
        ------------------------------------------------------------------
        - planet: body object of type 'planet' [-] ('body' object)
        - star: body object of type 'star' [-] ('body' object)
        - dt: time step [s] (float)
        - tf: final time to be computed [s] (float)
        OUTPUTS
        ------------------------------------------------------------------
        - planet: updated body object of type 'planet' [-] ('body' object)
        - star: updated body object of type 'star' [-] ('body' object)
        REFERENCES
        ------------------------------------------------------------------

    """
    #===========================================================
    #
    #  ZERO STEP: Load the problem parameters from the input data
    #
    #===========================================================

    G = 6.674e-11                             # Universal gravitational constant # [N⋅m2/kg2]
    m_p      = planet.properties.m            # Mass of planet [kg]
    m_s      = star.properties.m              # Mass of star   [kg]
    a_s     = planet.orbital_elements.a_b     # Planet-star semim.axis [m]
    e_s     = planet.orbital_elements.e_b     # Planet eccentricity [-]
    i_s     = planet.orbital_elements.i_b     # Planet inclination [deg]
    Omega_s = planet.orbital_elements.Omega_b # Planet RAAN [deg]
    omega_s = planet.orbital_elements.omega_b # Argument of planet periapsis [deg]
    t0      = planet.orbital_elements.t0_b    # Planet time of last per.  passage [s]

    #===========================================================
    #
    #  FIRST STEP: Create time vector and true anomaly vector
    #
    #===========================================================
    mu_s   = G*(m_p + m_s)       # planet-star stnd. grav. parameter [N⋅m2/kg]

    time   = np.arange(0, tf, dt, dtype=np.float)  # Time vector [s]
    n_s   = (mu_s/a_s**3)**0.5   # Barycenter-star mean motion [1/s]

    if t0 is None and t0 is None:
         f0     = planet.orbital_elements.f0_b
         E0s  = 2* np.atan( np.sqrt((1-e_s)/(1+e_s)) * np.tan(f0/2) )
         M0s  = E0s - e_s*np.sin(E0s)
         t0 = -M0s/n_s
    M_s = n_s * (time - t0)  # Barycenter-star mean anomaly [rad]
    # Eccentric anomaly and true anomaly are initialized
    E_s = np.zeros(len(time))
    f_s = np.zeros(len(time))

    for index in range(len(E_s)):
        # E is calculated through the Kepler equation. Then, f can be computed
        # The fun.symp function is used in order to reduce the angular range
        # from 0 to 2pi

        E_s[index] = kepler(fun.symp(M_s[index]), e_s, method = 'Newton-Raphson')
        f_s[index] = math.atan2((math.sin(E_s[index])*(1-e_s**2)**0.5)/
                                (1-e_s*math.cos(E_s[index])),( (math.cos(E_s[index])
                                     -e_s)/(1-e_s*math.cos(E_s[index]))))
        f_s[index] = fun.symp(f_s[index])
    #===========================================================
    #
    #  SECOND STEP: Keplerian orbit of planet (p) around the hosting star (s).
    #
    #===========================================================
    # Conversion of angles from degrees to radians
    i_s     = i_s     * math.pi / 180 # [rad]
    Omega_s = Omega_s * math.pi / 180 # [rad]
    omega_s = omega_s * math.pi / 180 # [rad]

    # Orbit equation for the barycenter
    r_s = a_s*(1-math.pow(e_s,2))/(1+e_s*np.cos(f_s)) # [meters]

    # 2-dimensional cartesian coordinates of the planet position
    # around the hosting star <position2D_s = (x_s, y_s, z_s=0)>
    x_s    = r_s * np.cos(f_s)    # [meters]
    y_s    = r_s * np.sin(f_s)    # [meters]
    position2D_s    = np.zeros((3,len(x_s)))  # [meters]
    position2D_s[0] = x_s                     # [meters]
    position2D_s[1] = y_s                     # [meters]

    # 3-dimensional cartesian coordinates of the planet position around
    # the star <position3D_s = (X_s, Y_s, Z_s)>
    #
    # The effect of orbit inclination, longitude of the ascending node and the
    # argument of periapsis are introduced by mean of three rotation matrices.

    Mz_Omega_s = np.array([[ np.cos(Omega_s), np.sin(Omega_s),   0   ],
                            [-np.sin(Omega_s), np.cos(Omega_s),   0   ],
                            [        0        ,        0        ,   1 ]]).T
    Mz_omega_s = np.array([[ np.cos(omega_s), np.sin(omega_s),   0   ],
                            [-np.sin(omega_s), np.cos(omega_s),   0   ],
                            [        0        ,        0        ,   1 ]]).T
    Mx_i_s     = np.array([[     1     ,        0        ,      0      ],
                            [     0     ,  np.cos(i_s)   , np.sin(i_s)],
                            [     0     , -np.sin(i_s)   , np.cos(i_s)]]).T

    # Three-dimensional motion of the planet with respect to the hosting
    # star

    position3D_s =Mz_Omega_s.dot(Mx_i_s).dot(Mz_omega_s).dot(position2D_s)
    r_ps = np.linalg.norm(position3D_s, axis = 0)     # [meters]

    # Export results to body objects
    planet.ephemeris.time         = time
    planet.orbital_elements.nu_ps = f_s
    planet.orbital_elements.M_ps  = M_s
    planet.orbital_elements.E_ps  = E_s
    planet.ephemeris.r_s          = r_ps
    planet.ephemeris.position2D_s = position2D_s
    planet.ephemeris.position3D_s = position3D_s
    planet.ephemeris.period_ps    = 2*math.pi/(n_s)
    star.ephemeris.time           = time
    star.orbital_elements.nu_ps   = f_s
    star.orbital_elements.M_ps    = M_s
    star.orbital_elements.E_ps    = E_s
    star.ephemeris.r_p            = r_ps
    star.ephemeris.position2D_ps  = position2D_s
    star.ephemeris.position3D_ps  = position3D_s
    star.ephemeris.period_ps      = 2*math.pi/(n_s)
    star.ephemeris.position3D_s   = np.zeros([3,len(time)])
    print('    ✓ The trajectories of ' + planet.name + ' and ' + star.name + ' have been calculated.\n')
    return planet, star


###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
###-------------------------- End of script --------------------------------'''
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
