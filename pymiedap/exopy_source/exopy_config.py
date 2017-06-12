# -*- coding: utf-8 -*-

"""
==================================================================
Exopy Configuration File: exopy_config.py
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
------------------------------------------------------------------

Basic configuration parameters for tuning some performance aspects
of the EXOPY tool. These are loaded when importing the EXOPY tool
and are accesible through exopy.doc.

LIST OF PARAMETERS
------------------------------------------------------------------
 - az: Azimuthal angle of the observer's position wrt the exosystem
       reference frame [deg?] (float)
 - el: Elevation angle of the observer's position wrt the exosystem 
       reference frame [deg?] (float)
 - approach: 'conical' or 'parallel' rays approach [-] (str)
 - case: ¿?
 - N:    ¿?
 - ref_body: reference body for the combination of the various ra-
	     diance signals [-] (str)
 - ref_line: reference body line for the combination of the various 
	     radiance signals [-] (str)
 - plot_color: plot color identifier [-] (int)
 - plot_faces: sets ON/OFF the plot of pixel borders [True/False] (bool)


"""

az = 0
el = 0
approach = 'conical'
case     = 'd'
N        = 200
ref_body = None
ref_line = None
plot_color = 0
plot_faces = True
