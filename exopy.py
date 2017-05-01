# -*- coding: utf-8 -*-

"""
==================================================================
EXOplanets & EXOmoons PYthon (EXOPY)
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
--------------------------------------- .  '  *     _..._       .  . '
                _______                   .  *    .' .::::.  * -+-  
       ___o |==|__(o(__D              .    * .   :  ::::::::    '  *
     /\  \/       /|\                     * .    :  ::::::::  ' .  .
    / /          / | \                 *   *  .  `. '::::::'     . 
    ` `         /  |  \                  '   *     `-.::''

Dependences: numpy, 
exopy_plot, exopy_config, exopy_compute

DESCRIPTION
------------------------------------------------------------------
This code is dedicated to defining an extrasolar planetary system
comprised of one planet and one moon hosted by a central star and 
computing theflux and polarization state of the reflected starlight
with the help of the PYMIEDAP code [reference].

In particular, the following actions can be conducted:

 * Computation of the three-dimensional orbits of the panet and
   moon around the hosting star.
 * Determination of the planet and moon shadowed region as a func-
   tion of phase angle and observer position.
 * Determination of shadows due to transits, i.e. one body putting 
   before a targeted body and the observer. 
 * Determination of shadows due to eclipses, i.e. one body putting
   before a targeted body and the star.
 * Computation of the individual reflected signals of the bodies 
   and combination of the results into a standalone output signal.

The different tools provided by exopy are given by the instances:
 
 * new_body: Used to create objects of type 'body', i.e. a planet,
   moon, or star.
 * new_system: Used to load a pre-defined planetary system.
 * run_simulation: Used to run a set of pre-defined simulation ins-
   tructions.
 * compute: Used to run any of the different computation modules
   available in EXOPY.
 * cfg: Used to specify a series of generic simulation settings.

Detailed instructions on how to use the EXOPY code and any of its
modules as well as example scripts are provided in the documenta-
tion of the various functions and in:
  ../exopy_source/EXOPY_how_to_use.html

"""

# ==============
# IMPORT MODULES
# ==============

import exopy_source.exopy_plot as plot
import exopy_source.exopy_config as cfg
import exopy_source.exopy_compute as compute
import exopy_source.exopy_atm_models as atm_models
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




def new_body(names, types):
    '''
==================================================================
EXOPY function: exopy.new_body
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
------------------------------------------------------------------

INPUTS
------------------------------------------------------------------
 - names: names of the bodies to be created [-] (str/list)
 - types: types of the bodies to be created [-] (str/list)
	  'planet', 'moon', 'star'

OUTPUTS
------------------------------------------------------------------
 - bodies: created bodies [-] ('body' object/list)

DESCRIPTION
------------------------------------------------------------------
Function of the EXOPY tool, creating a set of new 'body' objects
which serve as binding element for all information related to a 
celestial body.

    ''' 

    import exopy_source.exopy_functions as fun
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
==================================================================
EXOPY function: exopy.new_system
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
------------------------------------------------------------------

INPUTS
------------------------------------------------------------------
 - identifier: names of the system to be loaded [-] (str)

OUTPUTS
------------------------------------------------------------------
 - bodies: List of created bodies [-] (list)

DESCRIPTION
------------------------------------------------------------------
Function of the EXOPY tool, creating a set of new 'body' objects
which serve as binding element for all information related to a 
celestial body.

    ''' 
    import exopy_source.exopy_functions as fun

    if identifier == 'transit':

        Planet = fun.body('Planet', 'planet')
        Moon   = fun.body('Moon', 'moon')
        Star   = fun.body('Star', 'star')

        Planet .properties.m          = 5.972e24  # [kg]
        Moon.properties.m             = 7.342e22  # [kg]
        Star.properties.m             = 2e30      # [kg]

        Planet .properties.R          = 6.371008e6  # [m]
        Moon.properties.R             = 1.7374e6    # [m]
        Star.properties.R             = 6.957e8     # [m]

        Moon.orbital_elements.a       = 3.844e8*5.9723e24/(5.972e24+7.342e22) # [m]
        Moon.orbital_elements.e       = 0#0.0549  # 0.077 # [-]
        Moon.orbital_elements.i       = 90   # [deg]
        Moon.orbital_elements.Omega   = 0    # [deg]
        Moon.orbital_elements.omega   = 90   # [deg]
        Moon.orbital_elements.t0      = 0    # [s]

        Planet.orbital_elements.a_b     = 1.5e11    # [m]
        Planet.orbital_elements.e_b     = 0#0.0167  # [-]
        Planet.orbital_elements.i_b     = 0#23.5    # [deg]
        Planet.orbital_elements.Omega_b = 0.0       # [deg]
        Planet.orbital_elements.omega_b = 0.0       # [deg]
        Planet.orbital_elements.t0_b    = 0         # [s]

        bodies = [Moon, Planet , Star]


    if identifier == 'edge-on':

        Planet = fun.body('Planet', 'planet')
        Moon  = fun.body('Moon', 'moon')
        Star= fun.body('Star', 'star')

        Planet.properties.m = 5.972e24  # [kg]
        Moon.properties.m   = 7.342e22  # [kg]
        Star.properties.m   = 2e30      # [kg]

        Planet.properties.R = 6.371008e6  # [m]
        Moon.properties.R   = 1.7374e6    # [m]
        Star.properties.R   = 6.957e8     # [m]

        Moon.orbital_elements.a      = 3.844e8*5.9723e24/(5.972e24+7.342e22) # [m]
        Moon.orbital_elements.e      = 0#0.0549  # 0.077 # [-]
        Moon.orbital_elements.i      = 0         # [deg]
        Moon.orbital_elements.Omega  = 0         # [deg]
        Moon.orbital_elements.omega  = 90        # [deg]
        Moon.orbital_elements.t0     = 0         # [s]

        Planet.orbital_elements.a_b     = 1.5e11   # [m]
        Planet.orbital_elements.e_b     = 0#0.0167 # [-]
        Planet.orbital_elements.i_b     = 90#23.5  # [deg]
        Planet.orbital_elements.Omega_b = 0.0      # [deg]
        Planet.orbital_elements.omega_b = 0.0      # [deg]
        Planet.orbital_elements.t0_b    = 0        # [s]


        bodies = [Moon, Planet, Star]


    if identifier == 'face-on':

        Planet= fun.body('Planet', 'planet')
        Moon  = fun.body('Moon', 'moon')
        Star = fun.body('Star', 'star')

        Planet.properties.m = 5.972e24  # [kg]
        Moon.properties.m  = 7.342e22   # [kg]
        Star.properties.m   = 2e30      # [kg]

        Planet.properties.R = 6.371008e6  # [m]
        Moon.properties.R  = 1.7374e6     # [m]
        Star.properties.R   = 6.957e8     # [m]

        Moon.orbital_elements.a      = 3.844e8*5.9723e24/(5.972e24+7.342e22) # [m]
        Moon.orbital_elements.e      = 0#0.0549  # 0.077 # [-]
        Moon.orbital_elements.i      = 0      # [deg]
        Moon.orbital_elements.Omega  = 0      # [deg]
        Moon.orbital_elements.omega  = 90     # [deg]
        Moon.orbital_elements.t0     = 0      # [s]

        Planet.orbital_elements.a_b     = 1.5e11  # [m]
        Planet.orbital_elements.e_b     = 0#0.0167  # [-]
        Planet.orbital_elements.i_b     = 0#23.5    # [deg]
        Planet.orbital_elements.Omega_b = 0.0     # [deg]
        Planet.orbital_elements.omega_b = 0.0     # [deg]
        Planet.orbital_elements.t0_b    = 0       # [s]


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
        Moon.orbital_elements.omega  = 180.0   # [deg]
        Moon.orbital_elements.t0     = 0       # [s]

        Earth.orbital_elements.a_b     = 1.5e11  # [m]
        Earth.orbital_elements.e_b     = 0.0167  # [-]
        Earth.orbital_elements.i_b     = 23.5    # [deg]
        Earth.orbital_elements.Omega_b = 0.0     # [deg]
        Earth.orbital_elements.omega_b = 0.0     # [deg]
        Earth.orbital_elements.t0_b    = 0       # [s]


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
        Moon.orbital_elements.t0     = 0      # [s]

        Planet.orbital_elements.a_b     = 1e5    # [m]
        Planet.orbital_elements.e_b     = 0.017  # [-]
        Planet.orbital_elements.i_b     = 60     # [deg]
        Planet.orbital_elements.Omega_b = 270    # [deg]
        Planet.orbital_elements.omega_b = 0.0    # [deg]
        Planet.orbital_elements.t0_b    = 0      # [s]
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
        Moon.orbital_elements.e      = 0      # [-]
        Moon.orbital_elements.i      = 90     # [deg]
        Moon.orbital_elements.Omega  = 0      # [deg]
        Moon.orbital_elements.omega  = 90     # [deg]
        Moon.orbital_elements.t0     = 0      # [s]

        Planet.orbital_elements.a_b     = 1.5e11    # [m]
        Planet.orbital_elements.e_b     = 0         # [-]
        Planet.orbital_elements.i_b     = 0         # [deg]
        Planet.orbital_elements.Omega_b = 0         # [deg]
        Planet.orbital_elements.omega_b = 0         # [deg]
        Planet.orbital_elements.t0_b    = 0         # [s]


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
        Moon.orbital_elements.e      = 0.5     # [-]
        Moon.orbital_elements.i      = 0       # [deg]
        Moon.orbital_elements.Omega  = 0       # [deg]
        Moon.orbital_elements.omega  = 180     # [deg]
        Moon.orbital_elements.t0     = 0       # [s]

        Planet.orbital_elements.a_b     = 1.5e11    # [m]
        Planet.orbital_elements.e_b     = 0         # [-]
        Planet.orbital_elements.i_b     = 50        # [deg]
        Planet.orbital_elements.Omega_b = 0         # [deg]
        Planet.orbital_elements.omega_b = 300       # [deg]
        Planet.orbital_elements.t0_b    = 0         # [s]


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
        Moon.orbital_elements.e      = 0        # [-]
        Moon.orbital_elements.i      = 0        # [deg]
        Moon.orbital_elements.Omega  = 0        # [deg]
        Moon.orbital_elements.omega  = -1.3     # [deg]
        Moon.orbital_elements.t0     = 0        # [s]

        Planet.orbital_elements.a_b     = 1.5e11   # [m]
        Planet.orbital_elements.e_b     = 0        # [-]
        Planet.orbital_elements.i_b     = 90       # [deg]
        Planet.orbital_elements.Omega_b = 0        # [deg]
        Planet.orbital_elements.omega_b = 89.73    # [deg]
        Planet.orbital_elements.t0_b    = 0        # [s]


        bodies = [Moon, Earth, Sun]


    return bodies


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
        	print(i+1, ' out of ', len(time))
        	A = body.grid.shadow[i,:]>10E-10

        	I,Q,U,V = pmd.read_dap_output(np.repeat(phase[i], sum(A)), sza[i,A], emission[A], file2[scene], beta=beta[i,A], phi=phi[i,A])

        	output[0,i,A] = 4*np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[0,A]*I/np.pi
        	output[1,i,A] = 4*np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[1,A]*Q/np.pi
        	output[2,i,A] = 4*np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[2,A]*U/np.pi
        	print(t.time()-t1)

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

        	print(t22-t11)


	return [Ip, Qp, Up]



def run_simulation(body1, body2 ,star,dt, tf, flag_transits=True, flag_eclipses=True, flag_radiance=True, path_input = '../dap_database/', nmug=20, nmug_mie=20, nmat=4, nsubr=50):
    '''
==================================================================
EXOPY function: exopy.run_simulation
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
------------------------------------------------------------------

INPUTS
------------------------------------------------------------------
 - body1: planet or moon object [-] ('body' object)
 - body2: planet or moon object [-] ('body' object)
 - star: star object [-] ('body' object)
 - dt: time interval between two consecutive computed epochs [s] (float)
 - tf: final computation time [s] (float)
 - flag_transits: turns ON/OFF the computation of transits [True/False] (bool)
 - flag_eclipses: turns ON/OFF the computation of eclipses [True/False] (bool)
 - flag_radiance: turns ON/OFF the computation of reflected radiance [True/False] (bool)
 - path_input: location of the Fourier files database [-] (str)
 - nmug: number of Gauss points for DAP calculations [-] (int)
 - nmug_mie: number of Gauss points for MIE calculations [-] (int)
 - nmat: number of Stokes elements to compute [-] (int)
 - nsubr: number of subintervals for the distribution [-] (int)

OUTPUTS
------------------------------------------------------------------
 - body1: planet or moon object [-] ('body' object)
 - body2: planet or moon object [-] ('body' object)
 - star: star object [-] ('body' object)
 - I: First stokes vector: flux [normalized] (numpy array)
 - Q: Second stokes vector: linear polarization [normalized] (numpy
      array)
 - U: Third stokes vector: linear polarization [normalized] (numpy 
      array)
 - V: Fourth stokes vector: circular polarization [normalized] 
      (numpy array)

DESCRIPTION
------------------------------------------------------------------
Function of the EXOPY tool, conducting the following operations:

 1. Computes the orbits of the extrasolar planetary system through
    the compute.orbit subfunction.
 2. Computes the geometries involved in the motion of the bodies 
    at each time epoch.
 3. Computes the shadowed region of the bodies as a function of the
    phase angle and the observer's position
 4. Computes the shadowed region of the bodies due to transits.
 5. Computes the shadowed region of the bodies due to eclipses.
 6. Computes the integrated reflected radiance of each body.
 7. Combines the signal of the different bodies into a single one.


    ''' 
    if body1.type == 'moon':
        compute.orbit(body1, body2, star, dt, tf);
    elif body2.type == 'moon':
        compute.orbit(body2, body1, star, dt, tf);
    compute.geometry([body1, body2, star])
    compute.phases([body1,body2], star)
    if flag_transits:    compute.transits([body1, body2, star])
    if flag_eclipses:    compute.eclipses([body1,body2], star)
    if flag_radiance:
        compute.int_radiance([body1, body2], path_input = path_input, nmug=nmug, nmug_mie = nmug_mie, nmat=nmat, nsubr=nsubr)
        I,Q,U,V = compute.combine([body1, body2])

    return I,Q,U,V



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
#-------------------------- End of script --------------------------------'''
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
