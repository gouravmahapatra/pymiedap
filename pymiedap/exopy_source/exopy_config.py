# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

# -*- coding: utf-8 -*-

#"""
#==================================================================
#Exopy Configuration File: exopy_config.py
#Delft University of Technology
#------------------------------------------------------------------
#Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
#Date: 2016-2017
#------------------------------------------------------------------
#
#Basic configuration parameters for tuning some performance aspects
#of the EXOPY tool. These are loaded when importing the EXOPY tool
#and are accesible through exopy.doc.
#
#LIST OF PARAMETERS
#------------------------------------------------------------------
# - az: Azimuthal angle of the observer's position wrt the exosystem
#       reference frame [deg?] (float)
# - el: Elevation angle of the observer's position wrt the exosystem
#       reference frame [deg?] (float)
# - approach: 'conical' or 'parallel' rays approach [-] (str)
# - case: 
# - N:    
# - ref_body: reference body for the combination of the various ra-
#	     diance signals [-] (str)
# - ref_line: reference body line for the combination of the various
#	     radiance signals [-] (str)
# - plot_color: plot color identifier [-] (int)
# - plot_faces: sets ON/OFF the plot of pixel borders [True/False] (bool)
#
#
#"""

import numpy as _np

#az = 0
#el = 0
#approach = 'conical'
#case     = 'd'
#N        = 200
#ref_body = None
#ref_line = None
#plot_color = 0
#plot_faces = True

class Settings():
    ''' A Class to define computation settings'''

    def __init__(self, az=0, el=0, approach='conical', case='d', N=200,
                 ref_body=None, ref_line=None, plot_color=0, plot_faces=True):

        self.az=az
        self.el=el
        self.approach=approach
        self.case=case
        self.N=N
        self.ref_body=ref_body
        self.ref_line=ref_line
        self.plot_color=plot_color
        self.plot_faces=plot_faces

    def _plot_color(self, string):

        color = _np.array([[ 'faces'      , '#aaadab', 'b'  ],
                        [ 'nodes'      , '#aaadab', 'r'  ],
                        [ 'circle1'    , '#D3D3D3', 'k'  ],
                        [ 'circle'     , '#444544', 'r'  ],
                        [ 'sun'        ,  'y'     , 'y'  ],
                        [ 'umbra_c'    ,  'k'     , 'k'  ],
                        [ 'antumbra_c' ,  '0.4'   , '0.4'],
                        [ 'phase_c'    ,  'g'     , 'g'  ],
                        [ 'penumbra_c' ,  '0.7'   , '0.7'],
                        [ 'transit_c'  ,  'y'     , 'y'  ],
                        [ 'umbra'      ,  'k'     , 'k'  ],
                        [ 'antumbra'   ,  '0.4'   , '0.4'],
                        [ 'phase'      ,  'k'     , 'k'  ],
                        [ 'penumbra'   ,  '0.7'   , '0.7'],
                        [ 'transit'    ,  'k'     , 'k'  ],
                        [ 'background' ,  'w'     , 'w'  ],
                        [ 'border'     ,  'r'     , 'm'  ]])

        return color[_np.where(color==string)[0][0],1+self.plot_color]

