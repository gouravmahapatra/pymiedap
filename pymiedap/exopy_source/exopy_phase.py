# -*- coding: utf-8 -*-
"""
==================================================================
EXOPY module: exopy_phase.py
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
------------------------------------------------------------------

Dependences:

DESCRIPTION
------------------------------------------------------------------
Script containing the functions required for calculating the phase
shadow of a body.

LIST OF FUNCTIONS
------------------------------------------------------------------
 - phase: Function computing the darkenned pixels according to the
	  phase angle.


"""

import numpy as np
import exopy_config as _cfg
#from exopy_functions import PolyArea



def phase(body, star, conf):
    """
==================================================================
EXOPY function: phase()
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
------------------------------------------------------------------

Dependences:

DESCRIPTION
------------------------------------------------------------------
Computes the darkenned pixels according to the phase angle for
each time epoch.

INPUTS
------------------------------------------------------------------
 - body: Planet or moon type of body object [-] ('body' object)
 - star: Star type of body object [-] ('body' object)

OUTPUTS
------------------------------------------------------------------
 - body: Planet or moon type of body object [-] ('body' object)


    """

    print('\n    ... phase of ' + body.type + ' ' + body.name + '\n')

    approach = conf.approach
    N        = conf.N

    time = body.ephemeris.time
    T = len(time)
    t_count = np.array(range(len(time)))

    Rb  = body.properties.R
    Rs  = star.properties.R
    r_s = body.ephemeris.r_s

    alpha = body.geometry.phase_angle#alpha
    phi_s = body.geometry.solar_azimuth_angle

    if approach is 'parallel':
        xp = np.zeros_like(time)
        Rp = np.ones_like(time)
        angle = np.ones(T)*np.radians(90)
    elif approach is 'conical':
        xp = (Rs/Rb-1)*Rb/r_s
        Rp = (1-(xp)**2)**0.5
        angle = np.arctan(Rp/xp)

    Rob = body.ephemeris.position3D_s_ob

    body.flag.phase_d = True
    angle_nodes = np.arccos(np.einsum('Ni,it->tN', body.grid.nodes_xyz_rot, -Rob)/(0.5*r_s[:,np.newaxis]))
    body.grid.illuminated_nodes[angle_nodes>angle[:,np.newaxis]] = 0
    body.grid.shadow = body.grid.illuminated_nodes.astype(float)
    body.geometry.phase_area_d = np.degrees(np.sum(body.grid.area.reshape([1,len(body.grid.area)]).repeat(len(time),0)* ( ~body.grid.illuminated_nodes ), 1))*4

    return body
