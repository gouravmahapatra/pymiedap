# -*- coding: utf-8 -*-

# ==================================================================
# EXOPY module: exopy_compute.py
# Delft University of Technology
# ------------------------------------------------------------------
# Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
# Date: 2016-2017
# ------------------------------------------------------------------
#
# Dependences:
#
# DESCRIPTION
# ------------------------------------------------------------------
# The 'expoy_compute' script contain a series of functions which
# serve as interface between the user and the exopy functions.
#
# In essence, these organize the information, identify potential
# errors on the input data and takes the output back to the user.
#
# LIST OF FUNCTIONS
# ------------------------------------------------------------------
#  - orbit: Function computing the orbits of the bodies under the
# 	  assumption of a nested two body problem.
#  - geometry: Function computing the geometric parameters involved
# 	     in the motion of the planetary system.
#  - phases: Function computing the darkenned pixels according to the
# 	   phase angle.
#  - transits: Function computing the shadowed pixels due to transits
#  - eclipses: Function computing the shadowed pixels due to eclipses
#  - int_radiance: Function integrating the reflected radiance of
# 		 each body.
#  - combine: Function combining the signal of the different bodies
# 	    into a single one.
#
#


import exopy_config                    as _cfg
from exopy_geometry import geometry    as _geom
from exopy_phase    import phase       as _phase
from exopy_eclipse  import eclipse     as _eclipse
from exopy_transit  import transits    as _transit
from exopy_radiance import integration as _integ
from exopy_radiance import combine     as _comb
from exopy_orbit    import nested2bp   as _2bp



def orbit(moon, planet, star, delta_t, final_t):
    '''
    ==================================================================
    EXOPY function: exopy.compute.orbit
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    INPUTS
    ------------------------------------------------------------------
    - moon: moon object [-] ('body' object)
    - planet: planet object [-] ('body' object)
    - star: star object [-] ('body' object)
    - delta_t: time interval between two consecutive computed epochs [s] (float)
    - final_f: final computation time [s] (float)

    OUTPUTS
    ------------------------------------------------------------------
    - moon: updated moon object [-] ('body' object)
    - planet: updated planet object [-] ('body' object)
    - star: updated star object [-] ('body' object)

    DESCRIPTION
    ------------------------------------------------------------------
    Serves as interface between the user and the function computing
    the orbits of the extrasolar planetary system.

    '''

    moon, planet, star = _2bp(moon, planet, star, delta_t, final_t)

    return moon, planet, star


def geometry(bodies):
    '''
    ==================================================================
    EXOPY function: exopy.compute.geometry
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    INPUTS
    ------------------------------------------------------------------
    - bodies: list comprising a planet moon & star object [-] (list)

    OUTPUTS
    ------------------------------------------------------------------
    - bodies: updated list comprising a planet moon & star object [-]
        (list)

    DESCRIPTION
    ------------------------------------------------------------------
    Serves as interface between the user and the function computing
    the geometry parameters of the extrasolar planetary system.


    '''

    import sys

    if _cfg.ref_line == None:

        if type(bodies) != list:
            bodies = _geom(bodies)
        else:
            for i in range(len(bodies)):
                bodies[i] = _geom(bodies[i])

    else:

        if type(bodies) != list:
            if _cfg.ref_line != bodies.name:
                sys.exit('Error: value for ref_line does not match the name of any input body')

            bodies = _geom(bodies, ref_line_angle = 'fix')
        else:
            names = []
            for body in bodies: names.append(body.name)

            if _cfg.ref_line not in names:
                sys.exit('Error: value for ref_line does not match the name of any input body')
            else: index = names.index(_cfg.ref_line)
            bodies[index] = _geom(bodies[index], ref_line_angle = 'fix')
            for i in range(len(bodies)):
                if i == index: continue
                bodies[i] = _geom(bodies[i], ref_line_angle = bodies[index].geometry.ref_line_angle)

    return bodies


def phases(bodies, star):
    '''
    ==================================================================
    EXOPY function: exopy.compute.phases
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    INPUTS
    ------------------------------------------------------------------
    - bodies: list comprising a planet & moon object [-] (list)
    - star: star type of body object [-] (body object)

    OUTPUTS
    ------------------------------------------------------------------
    - bodies: updated list comprising a planet & moon object [-]
        (list)

    DESCRIPTION
    ------------------------------------------------------------------
    Serves as interface between the user and the function computing
    the phase shadowing of the extrasolar planetary system.


    '''

    import sys

    if type(bodies) != list:
        bodies = _phase(bodies, star)
    else:
        for i in range(len(bodies)):
            bodies[i] = _phase(bodies[i], star)
    #
    #    else:
    #
    #        if type(bodies) != list:
    #            if _cfg.ref_line != bodies.name:
    #                sys.exit('Error: value for ref_line does not match the name of any input body')
    #
    #            bodies = geom.phase(bodies, star, ref_line_angle = 'fix')
    #        else:
    #            names = []
    #            for body in bodies: names.append(body.name)
    #
    #            if _cfg.ref_line not in names:
    #                sys.exit('Error: value for ref_line does not match the name of any input body')
    #            else: index = names.index(_cfg.ref_line)
    #            bodies[index] = geom.phase(bodies[index], star)
    #            for i in range(len(bodies)):
    #                if i == index: continue
    #                bodies[i] = geom.phase(bodies[i], star)

    return bodies


def transits(bodies):
    '''
    ==================================================================
    EXOPY function: exopy.compute.transits
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    INPUTS
    ------------------------------------------------------------------
    - bodies: list comprising a series of body objects [-] (list)

    OUTPUTS
    ------------------------------------------------------------------
    - bodies: updated list of body objects [-] (list)

    DESCRIPTION
    ------------------------------------------------------------------
    Serves as interface between the user and the function computing
    the transits shadowing of the extrasolar planetary system.


    '''

    import numpy as np
    import sys

    names = []
    for body in bodies:
        if not hasattr(body.grid, 'shadow'):
           body.grid.shadow = np.ones([len(body.ephemeris.time), len(body.grid.nodes)])

        names.append(body.name)

    if _cfg.ref_line == None:
        ref_line_angle = np.zeros_like(bodies[0].ephemeris.time)
    else:

        if _cfg.ref_line not in names:
            sys.exit('Error: value for ref_line does not match the name of any input body')
        else: index = names.index(_cfg.ref_line)

        ref_line_angle = bodies[index].geometry.ref_line_angle

    bodies = _transit(bodies, ref_line_angle)

    return bodies


def eclipses(bodies, star):
    '''
    ==================================================================
    EXOPY function: exopy.compute.transits
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    INPUTS
    ------------------------------------------------------------------
    - bodies: list comprising a series of planet and moon objects [-]
        (list)
    - star: star type of body  object [-] (body object)

    OUTPUTS
    ------------------------------------------------------------------
    - bodies: updated list of planet and moon objects [-] (list)

    DESCRIPTION
    ------------------------------------------------------------------
    Serves as interface between the user and the function computing
    the eclipses shadowing of the extrasolar planetary system.


    '''

    bodies = _eclipse(bodies, star)

    return bodies


def int_radiance(bodies, path_input = './dap_database/', nmug = 20, nmug_mie = 20, nmat=4, nsubr=50):
    '''
    ==================================================================
    EXOPY function: exopy.compute.int_radiance
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    INPUTS
    ------------------------------------------------------------------
    - bodies: list comprising a series of planet and moon objects [-]
        (list)
    - path_input: Path to the Fourier files storage folder [-]
            ('body' object)
    - nmug: number of Gauss points for DAP calculations [-] (int)
    - nmug_mie: number of Gauss points for MIE calculations [-] (int)
    - nmat: number of Stokes elements to compute [-] (int)
    - nsubr: number of subintervals for the distribution [-] (int)

    OUTPUTS
    ------------------------------------------------------------------
    - bodies: updated list comprising a series of planet and moon
        objects [-] (list)

    DESCRIPTION
    ------------------------------------------------------------------
    Serves as interface between the user and the function integrating
    the reflected radiance of each body.


    '''

    for body in bodies:
        body = _integ(body, path_input = path_input, nmug = nmug, nmug_mie = nmug_mie, nmat=nmat, nsubr=nsubr)

    return bodies


def combine(bodies):
    '''
    ==================================================================
    EXOPY function: exopy.compute.combine
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    INPUTS
    ------------------------------------------------------------------
    - bodies: list comprising a series of planet and moon objects [-]
        (list)

    OUTPUTS
    ------------------------------------------------------------------
    - I: First stokes vector: flux [normalized] (numpy array)
    - Q: Second stokes vector: linear polarization [normalized] (numpy
        array)
    - U: Third stokes vector: linear polarization [normalized] (numpy
        array)
    - V: Fourth stokes vector: circular polarization [normalized]
        (numpy array)

    DESCRIPTION
    ------------------------------------------------------------------
    Serves as interface between the user and the function combining
    the reflected radiance by different bodies.


    '''

    for body in bodies:
        if body.name == _cfg.ref_body:
            I,Q,U,V = _comb(bodies,body)
            break

    return I,Q,U,V
