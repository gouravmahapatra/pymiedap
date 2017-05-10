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


def save_pickle(bodies, tf, dt, I, Q, U, V, exopy, directory = None,
                name = None, description = None):
    ''' 
    ==================================================================
    EXOPY function: exopy.save_pickle
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------
    
    INPUTS
    ------------------------------------------------------------------
     - bodies: list of body type of objects [-] (list)
     - tf: final time to be computed [s] (float)  
     - dt: time step [s] (float)
     - I: First stokes vector: flux [normalized] (numpy array)
     - Q: Second stokes vector: linear polarization [normalized] 
       (numpy array)
     - U: Third stokes vector: linear polarization [normalized] 
       (numpy array)
     - V: Fourth stokes vector: circular polarization [normalized]
       (numpy array)
     - exopy: names of the system to be loaded [-] (str)
     - directoy: path towards ouput folder [-] (str)
     - name: output file name [-] (str)
     - description: notes on stored data [-] (str)
    
    OUTPUTS
    ------------------------------------------------------------------
     - pickle file: output .pickle file [-] (.pickle file)
    
    DESCRIPTION
    ------------------------------------------------------------------
    Function creating and storing a .pickle file containing all 
    information on the current simulation.

    ''' 
    
    import os
    import pickle
    import os.path

    
    if directory is None or name is None:
        print(' ')
        print('Current directory: '+os.getcwd()+'/')
        dirs = [d for d in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), d))]
        for i in dirs:
            print i
        print(' ')
        directory = raw_input('Save directory: ')
        print(' ')
        print('List of existing files:')
        files = os.listdir(os.getcwd()+'/'+directory)
        for i in files:
            print i
        print('')
        name = raw_input('Name of the file: ')
    
    if description is None:
        print('')
        description = raw_input('Simulation description: ')
    
    save_data = [
        bodies,
        tf,
        dt,
        I,
        Q,
        U,
        V,
        exopy.cfg.az,
        exopy.cfg.el,
        exopy.cfg.approach,
        exopy.cfg.case,
        exopy.cfg.ref_body,
        exopy.cfg.ref_line,
        exopy.cfg.plot_color,
        exopy.cfg.N,
        description ]
    
    with open(os.getcwd()+'/'+directory+'/'+name+'.pickle', 'wb') as f:
         pickle.dump(save_data, f)
    
    f = open(os.getcwd()+'/'+directory+'/'+name+'.txt', 'w+')
    f.write('Description of the simulation: \n')
    f.write(description)
    f.write('\n')
    f.close
    
    print('\nData has been saved.')
    
#    del save_data, files, dirs, directory, name


def load_pickle(directory = None, name = None):
    ''' 
    ==================================================================
    EXOPY function: exopy.load_pickle
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------
    
    INPUTS
    ------------------------------------------------------------------
     - directoy: path towards ouput folder [-] (str)
     - name: output file name [-] (str) 
    
    OUTPUTS
    ------------------------------------------------------------------
     - bodies: list of body type of objects [-] (list)
     - tf: final time to be computed [s] (float)  
     - dt: time step [s] (float)
     - I: first stokes vector: flux [normalized] (numpy array)
     - Q: second stokes vector: linear polarization [normalized] 
       (numpy array)
     - U: third stokes vector: linear polarization [normalized] 
       (numpy array)
     - V: fourth stokes vector: circular polarization [normalized]
       (numpy array)
     - cfg: EXOPy configuration object [-] (exopy.cfg object)
     - description: notes on stored data [-] (str)
    
    DESCRIPTION
    ------------------------------------------------------------------
    Function loading a .pickle file containing all information on a 
    previous simulation.

    ''' 
    
    import os
    import pickle
    import os.path
    
    
    if directory is None:
        print(' ')
        print('Current directory: '+os.getcwd()+'/')
        dirs = [d for d in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), d))]
        for i in dirs:
            print i
        print(' ')
        directory = raw_input('Load directory: ')
        
    if name is None:
        print(' ')
        print('List of existing files:')
        files = os.listdir(os.getcwd()+'/'+directory)
        for i in files:
            print i
        print('')
        name = raw_input('Name of the file: ')
            
    with open(os.getcwd()+'/'+directory+'/'+name+'.pickle', 'rb') as f:
         load_data = pickle.load(f)
    
    conf = cfg
    
    bodies          = load_data[0]
    tf              = load_data[1]
    dt              = load_data[2]
    I               = load_data[3]
    Q               = load_data[4]
    U               = load_data[5]
    V               = load_data[6]
    conf.az         = load_data[7]
    conf.el         = load_data[8]
    conf.approach   = load_data[9]
    conf.case       = load_data[10]
    conf.ref_body   = load_data[11]
    conf.ref_line   = load_data[12]
    conf.plot_color = load_data[13]
    conf.N          = load_data[14]
    description     = load_data[15]

    print('\nData has been loaded.')

    return bodies, tf, dt, I, Q, U, V, cfg, description, 
    
    
    
def print_txt(bodies, tf, dt, I, Q, U, V, exopy, directory = None,
              name = None, description = None):
    ''' 
    ==================================================================
    EXOPY function: exopy.print_txt
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------
    
    INPUTS
    ------------------------------------------------------------------
     - bodies: list of body type of objects [-] (list)
     - tf: final time to be computed [s] (float)
     - dt: time step [s] (float)
     - I: First stokes vector: flux [normalized] (numpy array)
     - Q: Second stokes vector: linear polarization [normalized] 
       (numpy array)
     - U: Third stokes vector: linear polarization [normalized] 
       (numpy array)
     - V: Fourth stokes vector: circular polarization [normalized]
       (numpy array)
     - exopy: names of the system to be loaded [-] (str)
     - directoy: path towards ouput folder [-] (str)
     - name: output file name [-] (str)
     - description: notes on stored data [-] (str)
    
    OUTPUTS
    ------------------------------------------------------------------
     - text file: output .txt file [-] (.txt file)
    
    DESCRIPTION
    ------------------------------------------------------------------
    Function creating and storing a .txt file containing all 
    information on the current simulation.

    ''' 
    
    import os
    import os.path
    import numpy as np
    import time
    
    if directory is None or name is None:
        print(' ')
        print('Current directory: '+os.getcwd()+'/')
        dirs = [d for d in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), d))]
        for i in dirs:
        	print i
        print(' ')
        directory = raw_input('Save directory: ')
        print(' ')
        print('List of existing files:')
        files = os.listdir(os.getcwd()+'/'+directory)
        for i in files:
        	print i
        print('')
        name = raw_input('Name of the file: ')
    
    if description is None:
        print('')
        description = raw_input('Simulation description: ')
    
    file = open(os.getcwd()+'/'+directory+'/'+name+'.txt', 'w+')
    
    
    file.write('-- EXOPY OUTPUT DATA -- ' + time.strftime('%H:%M:%S') + ' -- ' + time.strftime('%d\%m\%Y') + ' --\n\n')
    
    file.write('Description:     \t' + description + '\n')
    file.write('Flag eclipses:   \t' + str(bodies[0].flag.eclipse_d) + '\nFlag transits:      \t' + str(bodies[0].flag.transit_d) + '\n')
    file.write('Reference line: \t' + exopy.cfg.ref_line + '\nReference body: \t' + exopy.cfg.ref_body + '\n')
    file.write('Computed epochs: \t' + str(len(bodies[0].ephemeris.time)) + '\n')
    file.write('Number of wavelengths: \t' + str(len(bodies[0].atmosphere.wvl_list)) + '\n\n')
    
    file.write('-- ORBITAL ELEMENTS --\n\n')
    
    file.write('           \t a[m]                  \t e[-]     \t i[deg] \t omega[deg] \t Omega[deg] \t t0[s]         \t Period[s]\n')
    file.write('Barycentre \t %-.10e \t %-.5e \t %-06.3f \t %-07.3f \t %-6.3f \t %-10.1f \t %-10.1f\n'   % (bodies[0].orbital_elements.a_b, bodies[0].orbital_elements.e_b, bodies[0].orbital_elements.i_b, bodies[0].orbital_elements.omega_b, bodies[0].orbital_elements.Omega_b, bodies[0].orbital_elements.t0_b, bodies[0].ephemeris.period_bs))
    file.write('Moon       \t %-.10e \t %-.5e \t %-06.3f \t %-07.3f \t %-6.3f \t %-10.1f \t %-10.1f\n\n' % (bodies[1].orbital_elements.a  , bodies[1].orbital_elements.e  , bodies[1].orbital_elements.i  , bodies[1].orbital_elements.omega  , bodies[1].orbital_elements.Omega  , bodies[1].orbital_elements.t0  , bodies[1].ephemeris.period   ))
    
    file.write('-- PROPERTIES OF BODIES --\n\n')
    
    file.write('           \t Radius[m]      \t Mass[kg]         \t Neq  \t Npoints  \t Fourier coeff.\n')
    file.write('Planet     \t %-.10e \t %-.10e \t %-3d \t %-6d \t %s \n'   % (bodies[0].properties.R, bodies[0].properties.m, bodies[0].grid.Nsq, bodies[0].grid.N_points, bodies[0].atmosphere.name[0]))
    file.write('Moon       \t %-.10e \t %-.10e \t %-3d \t %-6d \t %s \n'   % (bodies[1].properties.R, bodies[1].properties.m, bodies[1].grid.Nsq, bodies[1].grid.N_points, bodies[1].atmosphere.name[0]))
    file.write('Star       \t %-.10e \t %-.10e \t - \t - \t\t - \n\n' % (bodies[2].properties.R, bodies[2].properties.m))
    
    
    file.write('-- ATMOSPHERES --\n\n')
    
    file.write('               == Planet ==\n\n')
    file.write('mma [AMU] \t dpol [-] \t albedo [-] \t no. layers\n')
    file.write('%-8d \t %f   \t %f \t %d\n\n'  % (bodies[0].atmosphere.mma, bodies[0].atmosphere.dpol, bodies[0].atmosphere.surface[0,0], len(vars(bodies[0].atmosphere.layers).items())))
    file.write('Layer \t\t Type  \t Level \t P [bar] \t rayscat \t wvl [nm] \t tau [-] \t tau_gas [-]  \n')
    
    n_wvl = np.size(bodies[0].atmosphere.wvl_list)
    for layer_name, layer in vars(bodies[0].atmosphere.layers).items():
        if hasattr(layer,'mixed_aerosols') == True:
             file.write('%-10s \t %s \t %d \t %-0.2e \t %-7s \t %-06.4f \t %-06.3f \t %-06.3f \n' % (str(layer_name), layer.mixed_aerosols.typ, layer.level, layer.press, str(layer.rayscat), bodies[0].atmosphere.wvl_list[0], layer.tau[0], layer.tau_g[0]))
    	else:
         file.write('%-10s \t %s \t %d \t %0.2e \t %-7s \t %6.4f \t %06.3f \t %06.3f \n' % (str(layer_name), layer.aerosols.typ, layer.level, layer.press, str(layer.rayscat), bodies[0].atmosphere.wvl_list[0], layer.tau[0], layer.tau_g[0]))
    
    	for i in range(n_wvl-1):
         file.write('\t\t\t\t\t\t\t\t %6.4f \t %06.3f \t %06.3f \n' % (bodies[0].atmosphere.wvl_list[i+1], layer.tau[i+1], layer.tau_g[i+1]))
    
    file.write('               == Moon ==\n\n')
    file.write('mma [AMU] \t dpol [-] \t albedo [-] \t no. layers\n')
    file.write('%-8d \t %f   \t %f \t %d\n\n'  % (bodies[1].atmosphere.mma, bodies[1].atmosphere.dpol, bodies[1].atmosphere.surface[0,0], len(vars(bodies[1].atmosphere.layers).items())))
    file.write('Layer  \t\t Type  \t Level \t P [bar] \t rayscat \t wvl [nm] \t tau [-] \t tau_gas [-] \n')
    
    n_wvl = np.size(bodies[1].atmosphere.wvl_list)
    for layer_name, layer in vars(bodies[1].atmosphere.layers).items():
        if hasattr(layer,'mixed_aerosols') == True:
            file.write('%-10s \t %s \t %d \t %0.2e \t %-7s \t %6.4f \t %06.3f \t %06.3f \n' % (str(layer_name), layer.mixed_aerosols.typ, layer.level, layer.press, str(layer.rayscat), bodies[1].atmosphere.wvl_list[0], layer.tau[0], layer.tau_g[0]))
        else:
            file.write('%-10s \t %s \t %d \t %0.2e \t %-7s \t %6.4f \t %06.3f \t %06.3f \n' % (str(layer_name), layer.aerosols.typ, layer.level, layer.press, str(layer.rayscat), bodies[1].atmosphere.wvl_list[0], layer.tau[0], layer.tau_g[0]))
    
        for i in range(n_wvl-1):
            file.write('\t\t\t\t\t\t\t\t %6.4f \t %06.3f \t %06.3f \n' % (bodies[1].atmosphere.wvl_list[i+1], layer.tau[i+1], layer.tau_g[i+1]))
    
    
    file.write('\n-- POLARIZED SIGNAL --\n\n')
    
    file.write('Time [s] \t F_system [normalized] \t Q_system [normalized] \t U_system [normalized] \t V_system [normalized] \t P[%]        \t Chi[deg] \n\n')
    for i,t in enumerate(bodies[0].ephemeris.time):
        P   = np.sqrt(Q[0,i]**2 + U[0,i]**2)/I[0,i]*100
        Chi = 0.5*np.rad2deg(np.arctan(U[0,i]/Q[0,i])) 
        file.write('%10.1f \t %.10e \t %.10e \t %.10e \t %.10e \t %9.6f \t %.5e \n' % (bodies[0].ephemeris.time[i], I[0,i], Q[0,i], U[0,i], V[0,i], P, Chi ))
    file.write('\n\n')
    
    file.write('Time [s]\tAlpha [deg]\tP-S dist. [m]\tF_planet [-]\tQ_planet [-]\tU_planet [-]\tV_planet [-]\n\n')
    for i,t in enumerate(bodies[0].ephemeris.time):
        file.write('%10.0f\t%6.3f\t\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\n' % (bodies[0].ephemeris.time[i], np.rad2deg(bodies[0].geometry.phase_angle[i]), bodies[0].ephemeris.r_s[i], bodies[0].radiance.I[0,i], bodies[0].radiance.Q[0,i], bodies[0].radiance.U[0,i], bodies[0].radiance.V[0,i]))
    file.write('\n\n')
    
    file.write('Time [s]\tAlpha [deg]\tM-S dist. [m]\tF_moon [norm]\tQ_moon [-]\tU_moon [-]\tV_moon [-]\n\n')
    for i,t in enumerate(bodies[1].ephemeris.time):
        file.write('%10.0f\t%6.3f\t\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\n' % (bodies[1].ephemeris.time[i], np.rad2deg(bodies[1].geometry.phase_angle[i]), bodies[1].ephemeris.r_s[i], bodies[1].radiance.I[0,i], bodies[1].radiance.Q[0,i], bodies[1].radiance.U[0,i], bodies[1].radiance.V[0,i]))
    
    
    file.write('\n\n -- END OF FILE --')
    
    file.close()
    
#    del name, dirs, directory, files, i, t, description, n_wvl, layer, layer_name, P, Chi, 
    
    print('\nData has been printed.')
    
    
def read_txt(directory = None, name = None):
    ''' 
    ==================================================================
    EXOPY function: exopy.read_txt
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------
    
    INPUTS
    ------------------------------------------------------------------
     - directoy: path towards ouput folder [-] (str)
     - name: output file name [-] (str) 
    
    OUTPUTS
    ------------------------------------------------------------------
     - data: object containing information on loaded simulation [-] 
       (exopy.read_txt.reading_class class)
    
    DESCRIPTION
    ------------------------------------------------------------------
    Function loading a .txt file containing all information on a 
    previous simulation.

    '''     
    
    import os
    import os.path
    import numpy as np
    
    ref1 = 8 
    ref2 = ref1 + 5 
    ref3 = ref2 + 6 
    ref4 = ref3 + 8 + 14 + 1
    ref5 = ref4 + 4 
    ref6 = ref5 + 4 
    
    if directory is None:
        print(' ')
        print('Current directory: '+os.getcwd()+'/')
        dirs = [d for d in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), d))]
        for i in dirs:
                print i
        print(' ')
        directory = raw_input('Load directory: ')
        
    if name is None:
        print(' ')
        print('List of existing files:')
        files = os.listdir(os.getcwd()+'/'+directory)
        for i in files:
                print i
        print('')
        name = raw_input('Name of the file: ')
        
    file = open(os.getcwd() + '/' + directory + '/' + name +'.txt').read()
    lines = file.split('\n')
    totalline = len(lines)
    
    Nepochs = int(lines[ref1-1].split('\t')[1])
    Nwvl    = int(lines[ref1].split('\t')[1])
    
    class reading_class:
        def __init__(self,Nepochs):
            self.Rp         = 0
            self.Rm         = 0
            self.Rs         = 0
            self.Mp         = 0
            self.Mm         = 0
            self.Ms         = 0
            self.Neq_p      = 0
            self.Neq_m      = 0
            self.Npoint_p   = 0
            self.Npoint_m   = 0
            self.atm_name_p = 'None'
            self.atm_name_m = 'None'
            self.a_b        = 0
            self.a_m        = 0
            self.e_b        = 0
            self.e_m        = 0
            self.i_b        = 0
            self.i_m        = 0
            self.omega_b    = 0
            self.omega_m    = 0
            self.Omega_b    = 0
            self.Omega_m    = 0
            self.t0_p       = 0
            self.t0_m       = 0
            self.Period_p   = 0
            self.Period_m   = 0
            self.time       = np.zeros(Nepochs)
            self.F          = np.zeros(Nepochs)
            self.Q          = np.zeros(Nepochs)
            self.U          = np.zeros(Nepochs)
            self.V          = np.zeros(Nepochs)
            self.P          = np.zeros(Nepochs)
            self.Chi        = np.zeros(Nepochs)
            self.alpha_p    = np.zeros(Nepochs)
            self.alpha_m    = np.zeros(Nepochs)
            self.d_ps       = np.zeros(Nepochs)
            self.d_ms       = np.zeros(Nepochs)
            self.F_p        = np.zeros(Nepochs)
            self.Q_p        = np.zeros(Nepochs)
            self.U_p        = np.zeros(Nepochs)
            self.V_p        = np.zeros(Nepochs)
            self.F_m        = np.zeros(Nepochs)
            self.Q_m        = np.zeros(Nepochs)
            self.U_m        = np.zeros(Nepochs)
            self.V_m        = np.zeros(Nepochs)
    
    data = reading_class(Nepochs)
    
    aux             = lines[ref2].split('\t')
    data.a_b        = float(aux[1])
    data.e_b        = float(aux[2])
    data.i_b        = float(aux[3])
    data.omega_b    = float(aux[4])
    data.Omega_b    = float(aux[5])
    data.t0_b       = float(aux[6])
    data.Period_b   = float(aux[7])
    
    aux             = lines[ref2+1].split('\t')
    data.a_m        = float(aux[1])
    data.e_m        = float(aux[2])
    data.i_m        = float(aux[3])
    data.omega_m    = float(aux[4])
    data.Omega_m    = float(aux[5])
    data.t0_m       = float(aux[6])
    data.Period_m   = float(aux[7])
    
    aux             = lines[ref3].split('\t')
    data.Rp         = float(aux[1])
    data.Mp         = float(aux[2])
    data.Neq_p      = int(aux[3])
    data.Npoints_p  = int(aux[4])
    data.atm_name_p = aux[5]
    
    aux             = lines[ref3+1].split('\t')
    data.Rm         = float(aux[1])
    data.Mm         = float(aux[2])
    data.Neq_m      = int(aux[3])
    data.Npoints_m  = int(aux[4])
    data.atm_name_m = aux[5]
    
    aux     = lines[ref3+2].split('\t')
    data.Rs = float(aux[1])
    data.Ms = float(aux[2])
    
    layers_planet = int(lines[ref3+9].split('\t')[3])
    layers_moon   = int(lines[ref3+16-1+layers_planet*Nwvl].split('\t')[3])
    extra         = layers_planet*Nwvl + layers_moon*Nwvl
    
    for i in range(Nepochs):
        aux             = lines[ref4 + i + extra].split('\t')
        data.time[i]    = float(aux[0])
        data.F[i]       = float(aux[1])
        data.Q[i]       = float(aux[2])
        data.U[i]       = float(aux[3])
        data.V[i]       = float(aux[4])
        data.P[i]       = float(aux[5])
        data.Chi[i]     = float(aux[6])
    
    for i in range(Nepochs):
        aux             = lines[Nepochs + ref5 + i + extra].split('\t')
        data.alpha_p[i] = float(aux[1])
        data.d_ps[i]    = float(aux[3])
        data.F_p[i]     = float(aux[4])
        data.Q_p[i]     = float(aux[5])
        data.U_p[i]     = float(aux[6])
        data.V_p[i]     = float(aux[7])
    
    for i in range(Nepochs):
        aux             = lines[2*Nepochs + ref6 + i + extra].split('\t')
        data.alpha_m[i] = float(aux[1])
        data.d_ms[i]    = float(aux[3])
        data.F_m[i]     = float(aux[4])
        data.Q_m[i]     = float(aux[5])
        data.U_m[i]     = float(aux[6])
        data.V_m[i]     = float(aux[7])
    

    return data
    
    

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
        Moon.orbital_elements.omega  = 0         # [deg]
        Moon.orbital_elements.t0     = 0         # [s]

        Planet.orbital_elements.a_b     = 1.5e11   # [m]
        Planet.orbital_elements.e_b     = 0#0.0167 # [-]
        Planet.orbital_elements.i_b     = 90#23.5  # [deg]
        Planet.orbital_elements.Omega_b = 0.0      # [deg]
        Planet.orbital_elements.omega_b = 270.0    # [deg]
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

        Earth = fun.body('Planet', 'planet')
        Moon  = fun.body('Moon', 'moon')
        Sun   = fun.body('Sun', 'star')

        Earth.properties.m = 5.972e24  # [kg]
        Moon.properties.m  = 7.342e22  # [kg]
        Sun.properties.m   = 1.989e30      # [kg]

        Earth.properties.R = 6.371008e6  # [m]
        Moon.properties.R  = 1.7374e6    # [m]
        Sun.properties.R   = 6.957e8     # [m]

        Moon.orbital_elements.a      = 3.844e8*5.9723e24/(5.972e24+7.342e22) # [m]
        Moon.orbital_elements.e      = 0.0549  # 0.077 # [-]
        Moon.orbital_elements.i      = 5.145   # [deg]
        Moon.orbital_elements.Omega  = 0.0     # [deg]
        Moon.orbital_elements.omega  = 0.0   # [deg]
        Moon.orbital_elements.t0     = 0.0       # [s]

        Earth.orbital_elements.a_b     = 1.5e11  # [m]
        Earth.orbital_elements.e_b     = 0.0167  # [-]
        Earth.orbital_elements.i_b     = 90.0 #+23.5  # [deg]
        Earth.orbital_elements.Omega_b = 0.0     # [deg]
        Earth.orbital_elements.omega_b = 270.0   # [deg]
        Earth.orbital_elements.t0_b    = 0.0       # [s]


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



def run_simulation(body1, body2 ,star,dt, tf, flag_transits=True, flag_eclipses=True, flag_radiance=True, path_input = '../dap_database/'):
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
        compute.int_radiance([body1, body2], path_input = path_input)
        I,Q,U,V,P,Chi = compute.combine([body1, body2])

    return I,Q,U,V,P,Chi



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
#-------------------------- End of script --------------------------------'''
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
