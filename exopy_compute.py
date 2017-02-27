# -*- coding: utf-8 -*-

"""
Created on Thu Dec  1 01:36:41 2016

@author: javier
"""


import exopy_config                    as _cfg
from exopy_geometry import geometry    as _geom
from exopy_phase    import phase       as _phase
from exopy_eclipse  import eclipse     as _eclipse
from exopy_transit  import transits    as _transit
from exopy_radiance import integration as _integ
from exopy_radiance import combine     as _comb
from exopy_orbit    import nested2bp   as _2bp

def orbit(moon, planet, star, delta_t, final_t, t0_moon = 0, t0_planet = 0):
    import module_orbit as orbit
    
    moon, planet, star = _2bp(moon, planet, star, t0_moon, t0_planet, delta_t, final_t)
    
    return moon, planet, star


def geometry(bodies):
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
    
    bodies = _eclipse(bodies, star)  
    
    return bodies

def int_radiance(bodies):
    
    for body in bodies:
	body = _integ(body)

    return bodies


def combine(bodies):

    for body in bodies:
	if body.name == _cfg.ref_body:
	    I,Q,U,V = _comb(bodies,body)
	    break

    return I,Q,U,V
