# -*- coding: utf-8 -*-

"""
Exopy: Exoplanets and exomoons model
-------------------------------------
Author: J. Berzosa Molina

Date: 2016-2017
-------------------------------------

Model created in the framework of my master thesis work on exoplanets and exomoons.

Recipe: 

  1) Create a three body system.
  2) ...

  J. Berzosa

"""

# ==============
# IMPORT MODULES
# ==============

import exopy_plot as plot
import exopy_config as cfg
import exopy_compute as compute
from matplotlib import rc as _rc
_rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
_rc('text', usetex=True)
import matplotlib.pyplot as _plt
_plt.rc('text', usetex=True)
_plt.rc('font', family='serif')

# ======================================
# ======================================

if __name__ == "__main__":
    print(__doc__)  


def test_integration(body, scene='clear', plot = 'no'):

	import time as t
	import pymiedap.pymiedap as pmd
	import numpy as np
	import matplotlib.pyplot as plt


	file1 = {'clear':'exopy/files/table_clear.dat','cloudy':'exopy/files/table_cloudy.dat','gas':'exopy/files/table_gas.dat'  }
	file2 = {'clear':'exopy/files/clear_0.500.dat','cloudy':'exopy/files/cloudy_0.500.dat','gas':'exopy/files/fou_jup_100.dat'}
	
	t11 = t.time()
	
	time = body.ephemeris.time

	Ip = np.zeros_like(time)
	Qp = np.zeros_like(time)
	Up = np.zeros_like(time)

	output = np.zeros([4,len(body.ephemeris.time),body.grid.N_points])
	
	area  = np.repeat(body.grid.area[:,np.newaxis],4,1).T
	phase = np.degrees(body.geometry.phase_angle)
	sza   = np.degrees(body.grid.solar_zenith_angle)
	emission = np.degrees(body.grid.observer_zenith_angle)
	beta  = np.degrees(body.grid.beta)
    	phi   = np.degrees(body.grid.azimuth)
	
	for i,j in enumerate(time):
	    t1 = t.time()
	    print i+1, ' out of ', len(time)
	    A = body.grid.shadow[i,:]>10E-10

	    I,Q,U,V = pmd.read_dap_output(np.repeat(phase[i], sum(A)), sza[i,A], emission[A], file2[scene], beta=beta[i,A], phi=phi[i,A])

	    output[0,i,A] = 4*np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[0,A]*I/np.pi
	    output[1,i,A] = 4*np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[1,A]*Q/np.pi
	    output[2,i,A] = 4*np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[2,A]*U/np.pi
	    print t.time()-t1

	Ip = np.nansum(output[0,:],1)
	Qp = np.nansum(output[1,:],1)
	Up = np.nansum(output[2,:],1)

#==========================================================================
#==========================================================================
#==========================================================================
	if plot != 'no':

        	data = np.loadtxt(file1[scene])

        	phase_angle = data[:,0]
        	I           = data[:,1]
        	Q           = data[:,2]
        	U           = data[:,3]

     	  	del data,

     		f, (ax1, ax2) = plt.subplots(1,2)
        	ax1.plot(phase_angle,I,'b-')
        	ax1.plot(phase_angle,Q,'g-')
        	ax1.plot(phase_angle,U,'r-')
        	ax1.set_xlabel('Phase angle [deg]')

		ax1.plot(np.rad2deg(body.geometry.phase_angle), Ip,'o-b')
		ax1.plot(np.rad2deg(body.geometry.phase_angle), Qp,'og-')
		ax1.plot(np.rad2deg(body.geometry.phase_angle), Up,'o-r')

		ax1.plot(phase,np.interp(phase,phase_angle,I) ,'xb')
		ax1.plot(phase,np.interp(phase,phase_angle,U) ,'xr')
		ax1.plot(phase,np.interp(phase,phase_angle,Q) ,'xg')
		ax1.legend(['Reference I','Reference Q','Reference U','Computed I','Computed Q','Computed U'])

		ax2.plot(phase,Ip/np.interp(phase,phase_angle,I),'o-b')
		ax2.plot(phase,Qp/np.interp(phase,phase_angle,Q),'o-g')
		ax2.plot(phase,Up/np.interp(phase,phase_angle,U),'o-r')
		ax2.set_xlabel('Phase angle [deg]')
		ax2.set_ylabel('Difference')
		ax2.legend(['$\Delta I$','$\Delta Q$','$\Delta U$'])

		plt.tight_layout()

		t22 = t.time()

		print t22-t11

	
	return [Ip, Qp, Up]








def new_body(names, types):
    '''
===============================================================================
  Function: model.new_body        November 2016, Javier B.M., TU Delft
-------------------------------------------------------------------------------

The function creates a set of new body objects

Parameters
----------

names : str/list
    name of the body(ies) to be created
types : str/list
    type of body(ies) to be created: 'planet', 'moon', 'star'   

Returns
-------

bodies : model.body object/list
    Output generated body(ies)
   
Body object description
----------------------
    
A body object serves as binder element for all the information related to \
a celestial body. The following subclasses are associated to a body-type \
object:
    
* Properties: contains the main physical properties of the body

* Orbital elements: contains the orbital elements describing the orbit of the \
body  
    
* Ephemeris: contains the position of each body for each epoch

* Geometry: contains the angles defining the geometry for each epoch, as well \
as the shadowing conditions

* Flags: contains a series of flag indicators as a function of time

* Grid: contains the grid properties for discretizing the body

===============================================================================    
    '''
    import exopy_functions as fun
    import numpy as np
    
    if type(names)==str:
        bodies = fun.body(names, types)
    else:        
        bodies = list(np.zeros(len(names)))
        for i in range(len(bodies)):
            bodies[i] = fun.body(names[i], types[i])

    return bodies


def new_system(identifier):
    '''
===============================================================================
  Function: model.new_system        November 2016, Javier B.M., TU Delft
-------------------------------------------------------------------------------

The function creates a set of new body objects, based on a predefined system \
stored in the 'model.py' file.

Parameters
----------

identifier : str
    identification string of the system of bodies to be loaded   

Returns
-------

bodies : model.body object/list
    Output generated body(ies)

Body object description
----------------------

A body object serves as binder element for all the information related to \
a celestial body. The following subclasses are associated to a body-type \
object:

* Properties: contains the main physical properties of the body

* Orbital elements: contains the orbital elements describing the orbit of the \
body

* Ephemeris: contains the position of each body for each epoch

* Geometry: contains the angles defining the geometry for each epoch, as well \
as the shadowing conditions

* Flags: contains a series of flag indicators as a function of time

* Grid: contains the grid properties for discretizing the body

===============================================================================    
    '''    
    import exopy_functions as fun
    
    if identifier == 'transit':
    
        Planet = fun.body('Planet', 'planet')
        Moon   = fun.body('Moon', 'moon')
        Star   = fun.body('Star', 'star')
        
        Planet .properties.m = 5.972e24  # [kg]
        Moon.properties.m  = 7.342e22  # [kg]
        Star.properties.m   = 2e30      # [kg]
        
        Planet .properties.R = 6.371008e6  # [m]
        Moon.properties.R  = 1.7374e6    # [m]
        Star.properties.R   = 6.957e8     # [m]
        
        Moon.orbital_elements.a      = 3.844e8*5.9723e24/(5.972e24+7.342e22) # [m]
        Moon.orbital_elements.e      = 0#0.0549  # 0.077 # [-]
        Moon.orbital_elements.i      = 90   # [deg]
        Moon.orbital_elements.Omega  = 0     # [deg]
        Moon.orbital_elements.omega  = 90     # [deg]
        
        Star.orbital_elements.a_b     = 1.5e11  # [m]
        Star.orbital_elements.e_b     = 0#0.0167  # [-]
        Star.orbital_elements.i_b     = 0#23.5    # [deg]
        Star.orbital_elements.Omega_b = 0.0     # [deg]
        Star.orbital_elements.omega_b = 0.0     # [deg]
            
        bodies = [Moon, Planet , Star]


    if identifier == 'edge-on':

        Planet = fun.body('Planet', 'planet')
        Moon  = fun.body('Moon', 'moon')
        Star= fun.body('Star', 'star')

        Planet.properties.m = 5.972e24  # [kg]
        Moon.properties.m  = 7.342e22  # [kg]
        Star.properties.m   = 2e30      # [kg]

        Planet.properties.R = 6.371008e6  # [m]
        Moon.properties.R  = 1.7374e6    # [m]
        Star.properties.R   = 6.957e8     # [m]

        Moon.orbital_elements.a      = 3.844e8*5.9723e24/(5.972e24+7.342e22) # [m]
        Moon.orbital_elements.e      = 0#0.0549  # 0.077 # [-]
        Moon.orbital_elements.i      = 0   # [deg]
        Moon.orbital_elements.Omega  = 0     # [deg]
        Moon.orbital_elements.omega  = 90     # [deg]

        Star.orbital_elements.a_b     = 1.5e11  # [m]
        Star.orbital_elements.e_b     = 0#0.0167  # [-]
        Star.orbital_elements.i_b     = 90#23.5    # [deg]
        Star.orbital_elements.Omega_b = 0.0     # [deg]
        Star.orbital_elements.omega_b = 0.0     # [deg]

        bodies = [Moon, Planet, Star]


    if identifier == 'face-on':

        Planet= fun.body('Planet', 'planet')
        Moon  = fun.body('Moon', 'moon')
        Star = fun.body('Star', 'star')

        Planet.properties.m = 5.972e24  # [kg]
        Moon.properties.m  = 7.342e22  # [kg]
        Star.properties.m   = 2e30      # [kg]

        Planet.properties.R = 6.371008e6  # [m]
        Moon.properties.R  = 1.7374e6    # [m]
        Star.properties.R   = 6.957e8     # [m]

        Moon.orbital_elements.a      = 3.844e8*5.9723e24/(5.972e24+7.342e22) # [m]
        Moon.orbital_elements.e      = 0#0.0549  # 0.077 # [-]
        Moon.orbital_elements.i      = 0   # [deg]
        Moon.orbital_elements.Omega  = 0     # [deg]
        Moon.orbital_elements.omega  = 90     # [deg]

        Star.orbital_elements.a_b     = 1.5e11  # [m]
        Star.orbital_elements.e_b     = 0#0.0167  # [-]
        Star.orbital_elements.i_b     = 0#23.5    # [deg]
        Star.orbital_elements.Omega_b = 0.0     # [deg]
        Star.orbital_elements.omega_b = 0.0     # [deg]

        bodies = [Moon, Planet, Star]


        
    if identifier == 's_system':
    
        Earth = fun.body('Earth', 'planet')
        Moon  = fun.body('Moon', 'moon')
        Sun   = fun.body('Sun', 'star')
        
        Earth.properties.m = 5.972e24  # [kg]
        Moon.properties.m  = 7.342e22  # [kg]
        Sun.properties.m   = 2e30      # [kg]
        
        Earth.properties.R = 6.371008e6  # [m]
        Moon.properties.R  = 1.7374e6    # [m]
        Sun.properties.R   = 6.957e8     # [m]
        
        Moon.orbital_elements.a      = 3.844e8*5.9723e24/(5.972e24+7.342e22) # [m]
        Moon.orbital_elements.e      = 0.0549  # 0.077 # [-]
        Moon.orbital_elements.i      = 5.145   # [deg]
        Moon.orbital_elements.Omega  = 0.0     # [deg]
        Moon.orbital_elements.omega  = 0.0     # [deg]
        
        Sun.orbital_elements.a_b     = 1.5e11  # [m]
        Sun.orbital_elements.e_b     = 0.0167  # [-]
        Sun.orbital_elements.i_b     = 23.5    # [deg]
        Sun.orbital_elements.Omega_b = 0.0     # [deg]
        Sun.orbital_elements.omega_b = 0.0     # [deg]
            
        bodies = [Moon, Earth, Sun]
    
    if identifier == 'pool_ball':
        
        Earth = fun.body('Earth', 'planet')
        Moon  = fun.body('Moon', 'moon')
        Sun   = fun.body('Sun', 'star')
        
        Earth.properties.m = 6e24  # [kg]
        Moon.properties.m  = 0     # [kg]
        Sun.properties.m   = 2e30  # [kg]
        
        Earth.properties.R = 6.371e4  # [m]
        Moon.properties.R  = 1.737e6  # [m]
        Sun.properties.R   = 6.957e8  # [m]
        
        Moon.orbital_elements.a      = 3.8e8  # [m]
        Moon.orbital_elements.e      = 0.055  # [-]
        Moon.orbital_elements.i      = 5.145  # [deg]
        Moon.orbital_elements.Omega  = 0.0    # [deg]
        Moon.orbital_elements.omega  = 0.0    # [deg]
        
        Sun.orbital_elements.a_b     = 1e5    # [m]
        Sun.orbital_elements.e_b     = 0.017  # [-]
        Sun.orbital_elements.i_b     = 60     # [deg]
        Sun.orbital_elements.Omega_b = 270    # [deg]
        Sun.orbital_elements.omega_b = 0.0    # [deg]
        #  60*60*24*1e-9,60*60*24*3e-7)
        bodies = [Moon, Earth, Sun]
    
    if identifier == 'Moon_Earth_Sun_mod_eclipses_s_phase_meeting':
    
        Earth = fun.body('Earth', 'planet')
        Moon  = fun.body('Moon', 'moon')
        Sun   = fun.body('Sun', 'star')
        
        Earth.properties.m = 5.972e24  # [kg]
        Moon.properties.m  = 7.342e22  # [kg]
        Sun.properties.m   = 2e30      # [kg]
        
        Earth.properties.R = 6.371e6  # [m]
        Moon.properties.R  = 1.737e6  # [m]
        Sun.properties.R   = 6.957e8  # [m]
        
        Moon.orbital_elements.a      = 3.84748e8* 6e24/(6e24+7e22) # [m]
        Moon.orbital_elements.e      = 0 # [-]
        Moon.orbital_elements.i      = 90  # [deg]
        Moon.orbital_elements.Omega  = 0    # [deg]
        Moon.orbital_elements.omega  = 90    # [deg]
        
        Sun.orbital_elements.a_b     = 1.5e11 # [m]
        Sun.orbital_elements.e_b     = 0  # [-]
        Sun.orbital_elements.i_b     = 0     # [deg]
        Sun.orbital_elements.Omega_b = 0   # [deg]
        Sun.orbital_elements.omega_b = 0    # [deg]      
    
        bodies = [Moon, Earth, Sun]
    
    if identifier == 'Moon_Earth_Sun_mod_eclipses_s_meeting':
    
        Earth = fun.body('Earth', 'planet')
        Moon  = fun.body('Moon', 'moon')
        Sun   = fun.body('Sun', 'star')
        
        Earth.properties.m = 5.972e24  # [kg]
        Moon.properties.m  = 7.342e22  # [kg]
        Sun.properties.m   = 2e30      # [kg]
        
        Earth.properties.R = 6.371e6  # [m]
        Moon.properties.R  = 1.737e6  # [m]
        Sun.properties.R   = 6.957e8  # [m]
        
        Moon.orbital_elements.a      = 3.84748e8* 6e24/(6e24+7e22) # [m]
        Moon.orbital_elements.e      = 0.5 # [-]
        Moon.orbital_elements.i      = 0  # [deg]
        Moon.orbital_elements.Omega  = 0    # [deg]
        Moon.orbital_elements.omega  = 180    # [deg]
        
        Sun.orbital_elements.a_b     = 1.5e11 # [m]
        Sun.orbital_elements.e_b     = 0  # [-]
        Sun.orbital_elements.i_b     = 50     # [deg]
        Sun.orbital_elements.Omega_b = 0   # [deg]
        Sun.orbital_elements.omega_b = 300    # [deg]  
        
        bodies = [Moon, Earth, Sun]
        
    if identifier == 'Moon_Earth_Sun_mod_eclipses_ob_meeting':
        
        Earth = fun.body('Earth', 'planet')
        Moon  = fun.body('Moon', 'moon')
        Sun   = fun.body('Sun', 'star')
        
        Earth.properties.m = 5.972e24  # [kg]
        Moon.properties.m  = 7.342e22  # [kg]
        Sun.properties.m   = 2e30      # [kg]
        
        Earth.properties.R = 6.371e6  # [m]
        Moon.properties.R  = 1.737e6  # [m]
        Sun.properties.R   = 6.957e8  # [m]
        
        Moon.orbital_elements.a      = 3.84748e8* 6e24/(6e24+7e22) # [m]
        Moon.orbital_elements.e      = 0 # [-]
        Moon.orbital_elements.i      = 0  # [deg]
        Moon.orbital_elements.Omega  = 0    # [deg]
        Moon.orbital_elements.omega  = -1.3    # [deg]
        
        Sun.orbital_elements.a_b     = 1.5e11 # [m]
        Sun.orbital_elements.e_b     = 0  # [-]
        Sun.orbital_elements.i_b     = 90     # [deg]
        Sun.orbital_elements.Omega_b = 0   # [deg]
        Sun.orbital_elements.omega_b = 89.73    # [deg]          
        
        bodies = [Moon, Earth, Sun]
        
        
    return bodies
    

def run_simulation(body1, body2 ,star,dt, tf, t0_moon=0, t0_planet=0):
    if body1.type == 'moon':
        compute.orbit(body1, body2, star, dt, tf,t0_moon,t0_planet);
    elif body2.type == 'moon':
        compute.orbit(body2, body1, star, dt, tf,t0_moon,t0_planet);
    compute.geometry([body1, body2, star])
    compute.phases([body1,body2], star)
    compute.transits([body1, body2, star])
    compute.eclipses([body1,body2], star)

    body1.atmosphere.name = 'Hola'
    body2.atmosphere.name = 'Adios'
    compute.int_radiance([body1, body2])    
    I,Q,U,V = compute.combine([body1, body2])

    return I,Q,U,V



