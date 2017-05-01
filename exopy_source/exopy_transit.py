# -*- coding: utf-8 -*-
"""
==================================================================
EXOPY module: exopy_transit.py
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
------------------------------------------------------------------

Dependences:

DESCRIPTION
------------------------------------------------------------------
The 'expoy_transit' script contain the functions required for the
computation of the pixel darkening of the planetary system bodies 
due to shadowing through transits.

LIST OF FUNCTIONS
------------------------------------------------------------------
 - transits: Function computing the transits shadowing of the extra-
	     solar planetary system.


"""


import exopy_config as _cfg
import numpy as np
#from exopy_functions import PolyArea


def transits(bodies, ref_line_angle = None):   
    """
==================================================================
EXOPY function: eclipse()
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
------------------------------------------------------------------

Dependences:

DESCRIPTION
------------------------------------------------------------------
Computes the degree of pixel darkening due to eclipses shadowing 
among the different bodies.

INPUTS
------------------------------------------------------------------
 - bodies: list of body type of objects [-] (list)
 - ref_line_angle: shift angle towards reference line for the com-
		   putation of reflected radiance [rad] (numpy 
		   array)

OUTPUTS
------------------------------------------------------------------
 - bodies: updated body objects list [-] (list)


    """

#    '''
#    bodies = [Earth,Moon,Sun]
#    compute = 'cd'
#    N = 200
#    el = 0
#    az = 0
#    ref_line_angle = None
#    
#    for body in bodies:
#        if not hasattr(body.grid, 'shadow'):
#           body.grid.shadow = np.zeros([len(body.ephemeris.time), len(body.grid.nodes)])  
#    '''
    
    import itertools
    
    # Time vector
    time = bodies[0].ephemeris.time
    T = len(time)
    t_count = np.array(range(len(time)))
    
    N       = _cfg.N
    
    if ref_line_angle is None: ref_line_angle = np.zeros_like(time)
    
    # Algorithm control flags
    #merge_eclipses  = False

    # Create list for radii, re-initialize flags, and initialize points lists
    R = [bodies[i].properties.R for i in range(len(bodies))]
    for i in range(len(bodies)): bodies[i].flag.transit = []

    for i in range(len(bodies)): 
        bodies[i].flag.transit_d = True
                           
    # Bodies are arranged according to their radii, so that that R1 <= R2 <= R3
    bodies = [i for (r,i) in sorted(zip(R,bodies))]
    R.sort()
    
    # (n)
    # (2)  posibilities

    for i, j in itertools.combinations(range(len(bodies)),2):   

        bodies[i].flag.transit.append( [ np.zeros_like(time, dtype=bool), bodies[j].name ] )
        bodies[j].flag.transit.append( [ np.zeros_like(time, dtype=bool), bodies[i].name ] )
           
        print('\n    ... transits between ' + bodies[i].type + ' ' + bodies[i].name + ' and ' + bodies[j].type + ' ' + bodies[j].name + '\n')
                
        # Minimum distance for existance of transit
        L  = R[j] + R[i]
        M  = R[j] - R[i]

        # Ratio of radii
        r = R[i]/R[j]
            
        # A reference frame centered at bodyj is considered. The vector 
        # from body2 to bodyi is denominated 'd'.
        
        d = np.array([ bodies[i].ephemeris.position3D_s_ob[0,:] - 
                                       bodies[j].ephemeris.position3D_s_ob[0,:],
                       bodies[i].ephemeris.position3D_s_ob[1,:] - 
                                       bodies[j].ephemeris.position3D_s_ob[1,:],
                       bodies[i].ephemeris.position3D_s_ob[2,:] - 
                                       bodies[j].ephemeris.position3D_s_ob[2,:] ])    
        
        # Is bodyi in front of bodyj?
        bodyi_infront = np.array([False]).repeat(len(time))
        bodyi_infront[d[2,:]>0] = True
        
        # D: Distance from body2 to body1 in 2D plane
        D = np.linalg.norm([d[0,:], d[1,:]],axis = 0)     
        
        # Flags are assigned
        bodies[i].flag.transit[-1][0][np.invert(bodyi_infront)*(D<=L)] = True
        bodies[j].flag.transit[-1][0][(bodyi_infront)*(D<=L)] = True

        b = (np.linalg.norm( -bodies[j].grid.nodes.T[:,np.newaxis,:]*2*R[j] - d[0:2,bodyi_infront & (D<=L),np.newaxis]         , axis=0)<R[i])
        a = (np.linalg.norm( -bodies[i].grid.nodes.T[:,np.newaxis,:]*2*R[i] + d[0:2,np.invert(bodyi_infront)&(D<=L),np.newaxis], axis=0)<R[j])
	
	bodies[j].grid.shadow[bodyi_infront * (D<=L), :] = bodies[j].grid.shadow[bodyi_infront * (D<=L), :] * np.invert(b).astype(int)
	bodies[i].grid.shadow[np.invert(bodyi_infront) * (D<=L), :] = bodies[i].grid.shadow[np.invert(bodyi_infront) * (D<=L), :] * np.invert(a).astype(int)
                                        
    return bodies
