# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 01:33:54 2016

@author: javier
"""

def detail_radiance(bodies,I,Q,U,V, save = False):
	import module_plot as plot
	plot.detail_radiance(bodies[0],bodies[1],I,Q,U,V, save)

def radiance(body, phase = False, save = False):
	import module_plot as plot
	plot.radiance(body, phase, save)

def grid(body):
    import module_plot as plot

    plot.plot_grid(body)
    
def shadow_d(body, t=0, save = False, dots = False):
    import module_plot as plot
    import numpy as np
    if not hasattr(body.grid, 'shadow'):
       body.grid.shadow = np.ones([len(body.ephemeris.time), len(body.grid.nodes)])  
    
    plot.plot_shadow_d(body, t, save, dots) 
  

def shadow_dd(body, t=[0,0,0], save = False, dots = False):
    import module_plot as plot
    import numpy as np
    if not hasattr(body.grid, 'shadow'):
       body.grid.shadow = np.ones([len(body.ephemeris.time), len(body.grid.nodes)])  
    
    plot.plot_shadow_dd(body, t, save, dots) 

def geometry_d(body, t=0, save = False, dots = False):
    import module_plot as plot
    plot.plot_geometry_d(body, t, save, dots) 

def radiance_d(body, t=0, save = False, dots = False):
    import module_plot as plot
    plot.plot_radiance_d(body, t, save, dots) 

def I_d(body, t=0, dots = False):
    import module_plot as plot
    plot.plot_I_d(body, t, dots) 

def U_d(body, t=0, dots = False):
    import module_plot as plot
    plot.plot_U_d(body, t, dots)

def Q_d(body, t=0, dots = False):
    import module_plot as plot
    plot.plot_Q_d(body, t, dots)

def V_d(body, t=0, dots = False):
    import module_plot as plot
    plot.plot_V_d(body, t, dots)

def shadow_cd(body, t=0, dots = False):
    import module_plot as plot
    import numpy as np
    if not hasattr(body.grid, 'shadow'):
       body.grid.shadow = np.ones([len(body.ephemeris.time), len(body.grid.nodes)])  
    
    plot.plot_shadow_cd(body, t, dots)    

def shadow_c(body, t=0):
    import module_plot as plot

    plot.plot_shadow_c(body, t)
    
def shadow_cvsd(body, t=0, dots = False):
    import module_plot as plot
    import numpy as np
    if not hasattr(body.grid, 'shadow'):
       body.grid.shadow = np.ones([len(body.ephemeris.time), len(body.grid.nodes)])  
    
    plot.plot_shadow_cvsd(body, t, dots)

def anim_shadow_d(bodies, dots=False, time = 'all'):
    import module_plot as plot
    import sys
        
    if type(bodies) != list:
                
        plot = plot.plot_anim_shadow_d(bodies, dots, time)
            
    elif len(bodies) == 2:
               
        plot = plot.plot_anim_shadow_d2(bodies[0],bodies[1], dots, time)

    else:
        sys.exit('Error: Maximum number of bodies admitted in plot function is 2')

    return plot


def anim_shadow_c(bodies, time='all'):
    import module_plot as plot
    import sys
    
    if type(bodies) != list:
        
        plot = plot.plot_anim_shadow_c(bodies, time)

    elif len(bodies) == 2:
        
        plot = plot.plot_anim_shadow_c2(bodies[0],bodies[1],time)

    else:
        sys.exit('Error: Maximum number of bodies admitted in plot function is 2')

    return plot
    
def anim_shadow_cd(bodies, dots = False, time ='all'):
    import module_plot as plot
    import sys
    
    if type(bodies) != list:
        
        plot = plot.plot_anim_shadow_cd(bodies, dots, time)

    else:
        sys.exit('Error: Maximum number of bodies admitted in plot function is 2')

    return plot
    
def anim_shadow_cvsd(bodies, dots = False, time ='all'):
    import module_plot as plot
    import sys
    
    if type(bodies) != list:
        
        plot = plot.plot_anim_shadow_cvsd(bodies, dots, time)

    else:
        sys.exit('Error: Maximum number of bodies admitted in plot function is 2')

    return plot

def phasearea_cd(body):
    import module_plot as plot
    import sys
    
    if not body.flag.phase_d:
        sys.exit('Error in plotanim_shadow_cd: ' + body.type + ' ' + body.name + ' has flag phase_d = 0 --> no information to be plotted')
    if not body.flag.phase_c:
        sys.exit('Error in plotanim_shadow_cd: ' + body.type + ' ' + body.name + ' has flag phase_c = 0 --> no information to be plotted')        
    else:        
        plot.plotvs_phasearea_cd(body)
        
        
def anim_orbit(bodies):
    import module_plot as plot

    plot.plot_XYZorbitanimation(bodies[0].ephemeris.time,bodies[0].ephemeris.position3D_s,bodies[1].ephemeris.position3D_s)

def anim_planetmoon(time,bodies):
    import module_plot as plot

    plot.plot_planetmoon(time,bodies[0],bodies[1])



