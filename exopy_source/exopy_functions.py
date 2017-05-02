
# -*- coding: utf-8 -*-

"""
==================================================================
EXOPY module: exopy_functions.py
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
------------------------------------------------------------------

Dependences:

Module containing a set of common functions which are accessed by
the rest of the modules.

LIST OF CLASSES
------------------------------------------------------------------
 - body: class for the creation of planet, moon and star objects
 - properties: contains information regarding bodies' properties
 - ephemeris: contains information regarding bodies' ephemeris
 - geometry: contains information regarding bodies' geometry wrt
	     the extrasolar planetary system
 - radiance: contains information regarding the bodies' reflected
 	     starlight
 - orbital_elements: contains information regarding the orbital ele-
	  	     ments of the bodies' orbits
 - flag: contains flag indicators of various aspects related to a
	 body
 - grid: contains information regarding the bodies' grid
 - ProgressBar: class for the creation of a progress bar
 - Arrow3D: class for the creation of 3D arrow objects


LIST OF FUNCTIONS
------------------------------------------------------------------
 - anim_to_html: Converts animation to HTML for compatibility with
                 jupyter notebook.
 - display_animation: Displays animation in jupyter notebook
 - plot_xy: Predefined function for 2D plots.
 - quesry_yes_no: Asks the user a yes/no question.
 - grid_area: -- TBC --

REFERENCES
------------------------------------------------------------------

  [1] jakevdp.github.io/blog/2013/05/12/embedding-matplotlib-animations/
  [2] en.wikipedia.org/wiki/Shoelace_formula#
  [3] Recipe 577058 from activestate.com. Adapted to Python 3.
  [4] stackoverflow.com/questions/29188612/arrows-in-matplotlib-using-mplot3d


"""

# Required modules
from __future__ import print_function
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import sys, numpy as np
import math as m
import matplotlib.pyplot as plt
from IPython.display import HTML
from tempfile import NamedTemporaryFile
import base64
import re
import exopy_grid as grd
from .. import pymiedap as pmd

plt.ion()

#from shapely.geometry import Polygon as Poly
#from descartes import PolygonPatch
#from matplotlib.collections import PatchCollection

VIDEO_TAG = """<video controls>
 <source src="data:video/x-m4v;base64,{0}" type="video/mp4">
 Your browser does not support the video tag.
</video>"""


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
#--------------------------    CLASSES    --------------------------------'''
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''


class body():
    '''
    ==================================================================
    EXOPY class: body
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    The class 'body' provides a framework for the storage of all the
    information regaring the planet, moon, or star type of bodies

    METHODS
    ------------------------------------------------------------------
    - atmosphere: PYMIEDAP class for the definition and characteriza-
            tion of the body's surface and atmosphere.
    - ephemeris: class storing ephemeris information on the body's or-
            bit.
    - flag: class storing flag indicators.
    - geometry: class for the storage of geometry features of a body
            along its orbit which do not depend on the pixel dis-
            cretization.
    - grid: class for the characterization of the pixel discratization
        of a body and the storage of geometry information which is
        pixel dependant.
    - name: string name associated with the body.
    - note: string note providing more information if neded.
    - orbital_elements: class for the characterization of the body's
                orbit through its Keplerian elements.
    - properties: class storing the natural properties of a body.
    - radiance: class storing the output Stokes elements.
    - reset: !!??
    - show: method displaying a list of all created bodies.
    - type: string type associated with the body.

    '''
    __track       = [['#', 'Variable name', 'Type of body']]

    def __repr__(self):

        import numpy as np
        text = ['\n']
        text.append('Body object:\n' )
        text.append('   Name of body:  '   + self.name )
        text.append('   Type of body:  '   + self.type )
        text.append('   Note:          '   + self.note  + '\n' )
        text.append('   Radius:  '         + str(self.properties.R) + ' m' )
        text.append('   Mass:    '         + str(self.properties.m) + ' kg' )
        if self.properties.m != 0 and self.properties.R != 0:
                text.append('   Density: '         + str(self.properties.m/(4./3.*np.pi*self.properties.R**3)) + ' kg/m3'  )
        text.append('   Surface albedo: '  + str(self.atmosphere.surface[0,0]) + '\n')

        if self.type != 'star':
                text.append('   Wavelength [nm]\tFourier file')
                for i,w in enumerate(self.atmosphere.wvl_list):
                        text.append('       '+str(self.atmosphere.wvl_list[i])+'\t\t'+self.atmosphere.name[i])

                text.append(' ')
                text[-1] = str(self.flag)

                text.append(' ')
                if self.flag.orbit == False: aux=0
                else: aux = (len(str(np.size(self.ephemeris.time)))+1)*3
                if self.type == 'planet': aux1=12
                else: aux1 = 0
                text[-1] = str(self.orbital_elements)[0:-112-aux-aux1]

                #text.append('\n   Period: ' + str(self.ephemeris.period) + ' s  (' + str(self.ephemeris.period/31557600) + ' years)')
                #text.append('   Semi-major axis (a): '    + str(self.orbital_elements.a) + ' m')
                #text.append('   Eccentricity (e):    '    + str(self.orbital_elements.e))
                #text.append('   Inclination (i):     '    + str(self.orbital_elements.i) + ' deg')
                #text.append('   Argument of periapsis (omega): '  + str(self.orbital_elements.omega) + ' deg')
                #text.append('   Right ascension of the ascending node (Omega): '  + str(self.orbital_elements.Omega) + ' deg')
                #text.append('   Time from periapsis passage (t0): '  + str(self.orbital_elements.t0) + 's')

                text.append('\n   Grid type:\t'     + str(self.grid.type))
                text.append('   Grid points:\t'     + str(self.grid.N_points))
                text.append('   Grid eq. points: '  + str(self.grid.Nsq))

                text.append('\nRelation of properties at time index \'i\' available through body(i).')


        return '\n'.join(text)

    def __call__(self, index = 0):
        import numpy as np
        text = ['\n']
        text.append('Body object:\n' )
        text.append('   Name of body:  '   + self.name )
        text.append('   Type of body:  '   + self.type )
        text.append('   Note:          '   + self.note  + '\n' )
        text.append('   Radius:  '         + str(self.properties.R) + ' m' )
        text.append('   Mass:    '         + str(self.properties.m) + ' kg' )
        if self.properties.m != 0 and self.properties.R != 0:
                text.append('   Density: '         + str(self.properties.m/(4./3.*np.pi*self.properties.R**3)) + ' kg/m3'  )
                text.append('   Surface albedo: '  + str(self.atmosphere.surface[0,0]) + '\n')

        if self.type != 'star':
                text.append('   Wavelength [nm]\tFourier file')
                for i,w in enumerate(self.atmosphere.wvl_list):
                        text.append('       '+str(self.atmosphere.wvl_list[i])+'\t\t'+self.atmosphere.name[i])

        print('\n'.join(text))

        self.flag(index)
        self.orbital_elements(index=index)
        self.ephemeris(index=index)
        self.geometry (index=index,  unit='deg')



    def __init__(self, name, Type, note = ''):

        body.__updatelist(self,name)

        '''The new body is characterized by the setup function'''
        body.__setup(self,name,Type,note)
        print('    ✓ New',Type, name,'created!')

    def __updatelist(self,name):

        # If a body already exists in the bodies database, the older entries
        # are deleted

        while next((i for i, j in enumerate(body.__track)
                                              if name in j), None) is not None:
            i = next((i for i, j in enumerate(body.__track) if name in j),None)
            del(body.__track[i])

            ' The list is re-numbered '
            for index in range(len(body.__track)-i):
                body.__track[index+i][0] = i+index

            print('    ! The object has been overwritten in the body list...')

    def __setup(self,name,Type,note):
        self.note = note
        self.name = name
        self.type = Type
        body.__track.append([len(body.__track), name, Type])

        self.properties       = properties()
        self.ephemeris        = ephemeris(Type)
        self.orbital_elements = orbital_elements(Type)
        self.geometry         = geometry(Type)
        self.flag             = flag(Type)
        self.grid             = grid(Type)
        self.radiance         = radiance(Type)
        self.atmosphere       = pmd.Model()

    def show(self):
        print('\nList of bodies:\n')

        for index in range(len(body.__track)):
            print('%8s  %12s  %12s'
                % (body.__track[index][0], body.__track[index][1],
                   body.__track[index][2]))
        if len(body.__track)==1:
            print('%25s' %('EMPTY'))
        print('\n')

    def reset():
        print('\n')
        if query_yes_no('Note: Reseting the database does not delete the \
        objects from the memory.\n\ Do you want to continue?', default="no"):
            del(body.__track)
            body.__track = [['#', 'Variable name', 'Type of body']]
            print('\nDatabase reseted!\n')
            body.show()
        else:
            print('\nReset aborted.\n')


class properties():
    '''
    ==================================================================
    EXOPY class: body.properties
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    The class 'properties' stores all relevant information regarding
    the natural properties of a body.

    METHODS
    ------------------------------------------------------------------
    - m: mass of the body [kg] (float)
    - R: radius of the body [m] (float)


    '''

    def __repr__(self):
        text = ['\n']
        text.append('Relation of body properties:\n' )
        text.append('   Mass: m = ' + str(self.m) + ' kg' )
        text.append('   Radius: R = ' + str(self.R) + ' m')
        return '\n'.join(text)

    def __init__(self):

        self.m   = 0
        self.R   = 0
        self.fourier_scene = 'clear' # !!!

    def __call__(self):
        text = ['\n']
        text.append('Relation of body properties:\n' )
        text.append('   Mass: m = ' + str(self.m) + ' kg' )
        text.append('   Radius: R = ' + str(self.R) + ' m')

        print('\n'.join(text))


class ephemeris():
    '''
    ==================================================================
    EXOPY class: body.ephemeris
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    The class 'ephemeris' stores all relevant information regarding
    the ephemeris of a body.

    METHODS
    ------------------------------------------------------------------
    - period: period of the body's orbit [s] (float)
    - position2D_ij: 2D position vector of body i wrt body j [m]
                    (numpy array) *
    - position3D_ij: 3D position vector of body i wrt body j [m]
                    (numpy array) *
    - position3D_s: 3D position vector wrt star [m] (numpy array)
    - position3D_s_ob: 3D position vector wrt star in observer's ref.
                frame [m] (numpy array)
    - r_b: distance to planet-moon system barycentre [m] (numpy array)
    - r_m: distance to moon [m] (numpy array)
    - r_s: distance to planet [m] (numpy array)
    - time: Array of computed time epochs [s] (numpy array)

    * i,j = p, m, b, s
    with p=planet, m=moon, b=planet-moon barycentre, s=star

    A relation of ephemeris values at time epoch 'i' can be retrieved
    by body.ephemeris(i).

    '''

    def __repr__(self):
        import numpy as np
        text = ['\n']
        text.append('Relation of body ephemeris:\n' )
        if self.__type == 'planet':
                text.append('   period = ' + str(self.period)  + ' s  (' + str(self.period/31557600) + ' years)')
                text.append('   position2D_pb   = ' + str(np.shape(self.position2D_pb))  + ' array [m]' )
                text.append('   position2D_pm   = ' + str(np.shape(self.position2D_pm))  + ' array [m]' )
                text.append('   position3D_pb   = ' + str(np.shape(self.position3D_pb))  + ' array [m]' )
                text.append('   position3D_pm   = ' + str(np.shape(self.position3D_pm))  + ' array [m]' )
                text.append('   position3D_s    = ' + str(np.shape(self.position3D_s))   + ' array [m]' )
                text.append('   position3D_s_ob = ' + str(np.shape(self.position3D_s_ob))+ ' array [m]' )
                text.append('   r_b  = '            + str(np.shape(self.r_b))            + ' array [m]' )
                text.append('   r_m  = '            + str(np.shape(self.r_m))            + ' array [m]' )
                text.append('   r_s  = '            + str(np.shape(self.r_s))            + ' array [m]' )

        if self.__type == 'moon':
                text.append('   period = ' + str(self.period)  + ' s  (' + str(self.period/31557600) + ' years)')
                text.append('   position2D_mb   = ' + str(np.shape(self.position2D_mb))  + ' array [m]' )
                text.append('   position2D_mp   = ' + str(np.shape(self.position2D_mp))  + ' array [m]' )
                text.append('   position3D_mb   = ' + str(np.shape(self.position3D_mb))  + ' array [m]' )
                text.append('   position3D_mp   = ' + str(np.shape(self.position3D_mp))  + ' array [m]' )
                text.append('   position3D_s    = ' + str(np.shape(self.position3D_s))   + ' array [m]' )
                text.append('   position3D_s_ob = ' + str(np.shape(self.position3D_s_ob))+ ' array [m]' )
                text.append('   r_b  = '            + str(np.shape(self.r_b))            + ' array [m]' )
                text.append('   r_p  = '            + str(np.shape(self.r_p))            + ' array [m]' )
                text.append('   r_s  = '            + str(np.shape(self.r_s))            + ' array [m]' )

        if self.__type == 'star':
                text.append('   period_bs = ' + str(self.period_bs)  + ' s  (' + str(self.period_bs/31557600) + ' years)')
                text.append('   position2D_bs   = '+ str(np.shape(self.position2D_bs))  + ' array [m]' )
                text.append('   position3D_bs   = '+ str(np.shape(self.position3D_bs))  + ' array [m]' )
                text.append('   position3D_s    = '+ str(np.shape(self.position3D_s))   + ' array [m]' )
                text.append('   position3D_s_ob = '+ str(np.shape(self.position3D_s_ob))+ ' array [m]' )
                text.append('   r_b  = '           + str(np.shape(self.r_b))            + ' array [m]' )
                text.append('   r_m  = '           + str(np.shape(self.r_m))            + ' array [m]' )
                text.append('   r_p  = '           + str(np.shape(self.r_p))            + ' array [m]' )
        text.append('   time = ' + str(np.shape(self.time))  + ' array [s]' )
        return '\n'.join(text)

    def __init__(self,Type):
        self.time   = None
        self.__type  = Type
        if Type == 'planet':
                self.period 		= 0.
                self.position2D_pb   = None
                self.position2D_pm   = None
                self.position3D_pb   = None
                self.position3D_pm   = None
                self.position3D_s    = None
                self.position3D_s_ob = None
                self.r_b             = None
                self.r_m 		= None
                self.r_s 		= None

        if Type == 'moon':
                self.period 		= 0.
                self.position2D_mb   = None
                self.position2D_mp   = None
                self.position3D_mb   = None
                self.position3D_mp   = None
                self.position3D_s    = None
                self.position3D_s_ob = None
                self.r_b             = None
                self.r_p 		= None
                self.r_s 		= None

        if Type == 'star':
                self.period_bs	= 0.
                self.position2D_bs   = None
                self.position3D_bs   = None
                self.position3D_s    = None
                self.position3D_s_ob = None
                self.r_b             = None
                self.r_m 		= None
                self.r_p 		= None

    def __call__(self, index = 0):

        print('\nRelation of body ephemeris variables at time index ' + str(index) + ':\n')

        if self.__type == 'planet':
                print('   period = ' + str(self.period)  + ' s  (' + str(self.period/31557600) + ' years)')
                print('   position2D_pb   = ' + str(self.position2D_pb[:,index])  + ' m' )
                print('   position2D_pm   = ' + str(self.position2D_pm[:,index])  + ' m' )
                print('   position3D_pb   = ' + str(self.position3D_pb[:,index])  + ' m' )
                print('   position3D_pm   = ' + str(self.position3D_pm[:,index])  + ' m' )
                print('   position3D_s    = ' + str(self.position3D_s[:,index])   + ' m' )
                print('   position3D_s_ob = ' + str(self.position3D_s_ob[:,index])+ ' m' )
                print('   r_b  = '             + str(self.r_b[index])              + ' m' )
                print('   r_m  = '             + str(self.r_m[index])              + ' m' )
                print('   r_s  = '             + str(self.r_s[index])              + ' m' )

        if self.__type == 'moon':
                print('   period = ' + str(self.period)  + ' s  (' + str(self.period/31557600) + ' years)')
                print('   position2D_mb   = ' + str(self.position2D_mb[:,index])  + ' m' )
                print('   position2D_mp   = ' + str(self.position2D_mp[:,index])  + ' m' )
                print('   position3D_mb   = ' + str(self.position3D_mb[:,index])  + ' m' )
                print('   position3D_mp   = ' + str(self.position3D_mp[:,index])  + ' m' )
                print('   position3D_s    = ' + str(self.position3D_s[:,index])   + ' m' )
                print('   position3D_s_ob = ' + str(self.position3D_s_ob[:,index])+ ' m' )
                print('   r_b  = '             + str(self.r_b[index])              + ' m' )
                print('   r_p  = '             + str(self.r_p[index])              + ' m' )
                print('   r_s  = '             + str(self.r_s[index])              + ' m' )

        if self.__type == 'star':
                print('   period_bs = ' + str(self.period_bs)  + ' s  (' + str(self.period_bs/31557600) + ' years)')
                print('   position2D_bs   = ' + str(self.position2D_bs[:,index])  + ' m' )
                print('   position3D_bs   = ' + str(self.position3D_bs[:,index])  + ' m' )
                print('   position3D_s    = ' + str(self.position3D_s[:,index])   + ' m' )
                print('   position3D_s_ob = ' + str(self.position3D_s_ob[:,index])+ ' m' )
                print('   r_b  = '             + str(self.r_b[index])              + ' m' )
                print('   r_m  = '             + str(self.r_m[index])              + ' m' )
                print('   r_p  = '             + str(self.r_p[index])              + ' m' )

        print('   time = ' + str(self.time[index])  + ' s' )




class geometry():
    '''
    ==================================================================
    EXOPY class: body.geometry
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    The class 'geometry' stores all relevant information regarding the
    geometry of a body in the extrasolar planetary system which do not
    depend on the pixel discretization of the disk.

    METHODS
    ------------------------------------------------------------------
    - phase_angle: angle between the star, planet, and observer [rad]
                (numpy array)
    - ref_line_angle: shift angle to centre of reference body [rad]
            (numpy array)
    - ref_plane_angle: shift angle to reference angle for radiance
                calc. [rad] (numpy array)
    - ref_plane_to_ref_line_angle: difference between the ref. line
                    and plane angles (numpy array)
    - solar_azimuth_angle: azimuth angle of star as seen from the body
                in observer reference frame [rad]
                (numpy array)
    - alpha: supplementary angle to phase angle [rad] (numpy array)

    A relation of ephemeris values at time epoch 'i' can be retrieved
    by body.geometry(i).

    '''

    def __repr__(self):
        import numpy as np
        text = ['\n']
        text.append('Relation of body geometry variables:\n')
        if self.__type != 'star':
                text.append('   phase_angle      = ' + str(np.shape(self.phase_angle))     + ' array [rad]' )
                text.append('   ref_line_angle   = ' + str(np.shape(self.ref_line_angle))  + ' array [rad]' )
                text.append('   ref_plane_angle  = ' + str(np.shape(self.ref_plane_angle)) + ' array [rad]' )
                text.append('   ref_plane_to_ref_line_angle  = ' + str(np.shape(self.ref_plane_to_ref_line_angle))  + ' array [rad]' )
                text.append('   solar_azimuth_angle  = ' + str(np.shape(self.solar_azimuth_angle))  + ' array [rad]' )
                text.append('   alpha = ' + str(np.shape(self.alpha)) + ' array [rad]' )
        else:
                text.append('   --No data available--')

        return '\n'.join(text)

    def __init__(self,Type):
        if Type!= 'star':
                self.phase_angle    			= None
                self.ref_line_angle 		  	= None
                self.ref_plane_angle  		= None
                self.ref_plane_to_ref_line_angle  	= None
                self.solar_azimuth_angle          	= None
                self.alpha			  	= None
        self.__type  = Type

    def __call__(self, index = 0, unit = 'rad'):
        import numpy as np

        scale = 1
        if unit == 'deg':
                scale = 180./np.pi
        elif unit != 'rad':
                raise ValueError('Input value for input \'unit\' must be \'rad\' (default) or \'deg\'.')

        text = ['\n']

        if self.__type != 'star':

                if np.ndim(self.phase_angle)== 1:

                        if index > np.size(self.phase_angle,0)-1:
                                text.append('   ! Index is larger than arrays\' dimension')
                        else:
                                text.append('Relation of body geometry variables at time index ' + str(index) + ':\n')

                                text.append('   phase_angle      = ' + str(self.phase_angle[index]*scale)     + ' ' + unit)
                                text.append('   ref_line_angle   = ' + str(self.ref_line_angle[index]*scale)  + ' ' + unit)
                                text.append('   ref_plane_angle  = ' + str(self.ref_plane_angle[index]*scale) + ' ' + unit)
                                text.append('   ref_plane_to_ref_line_angle  = ' + str(self.ref_plane_to_ref_line_angle[index]*scale)  + ' ' + unit)
                                text.append('   solar_azimuth_angle  = ' + str(self.solar_azimuth_angle[index]*scale)  + ' ' + unit)
                                text.append('   alpha = ' + str(self.alpha[index]*scale) + ' ' + unit)
                                if unit == 'rad': print('\n   Note: show angles in deg via input parameter unit = \'deg\'.')

                else:

                        text.append('Relation of body geometry variables:\n')

                        text.append('   phase_angle      = ' + str(self.phase_angle)     + ' ' + unit)
                        text.append('   ref_line_angle   = ' + str(self.ref_line_angle)  + ' ' + unit)
                        text.append('   ref_plane_angle  = ' + str(self.ref_plane_angle) + ' ' + unit)
                        text.append('   ref_plane_to_ref_line_angle  = ' + str(self.ref_plane_to_ref_line_angle)  + ' ' + unit)
                        text.append('   solar_azimuth_angle  = ' + str(self.solar_azimuth_angle)  + ' ' + unit)
                        text.append('   alpha = ' + str(self.alpha) + ' ' + unit)
                        if unit == 'rad': print('\n   Note: show angles in deg via input parameter unit = \'deg\'.')

        else:
                text.append('   --No data available--')

        print('\n'.join(text))



class radiance():
    '''
    ==================================================================
    EXOPY class: body.radiance
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    The class 'radiance' stores all relevant information regarding the
    body's reflected starlight as produced by the PYMIEDAP code. Pixel
    by pixel contributions to the total flux are stored under the
    'grid' category.

    METHODS
    ------------------------------------------------------------------
    - I: First Stokes element: reflected starlight flux [normalized]
        (numpy array)
    - U: Second Stokes element: linear polarization of reflected
        starlight [normalized] (numpy array)
    - V: Third Stokes element: circular polarization of reflected
        starlight [normalized] (numpy array)
    - Q: Fourth Stokes element: circular polarization of reflected
        starlight [normalized] (numpy array)
    - I_ref: First Stokes element wrt reference plane: reflected
        starlight flux [normalized] (numpy array)
    - U_ref: Second Stokes element wrt reference plane: linear pola-
        rization of reflected
            starlight [normalized] (numpy array)
    - V_ref: Third Stokes element wrt reference plane: circular pola-
        rization of reflected starlight [normalized]
        (numpy array)
    - Q_ref: Fourth Stokes element wrt reference plane: circular po-
        larization of reflected starlight [normalized]
        (numpy array)

    A relation of radiance values at time epoch 'i' can be retrieved
    by body.radiance(i).

    '''

    def __init__(self, Type):

        if Type != 'star':
            self.I 		= None
            self.U 		= None
            self.V 		= None
            self.Q 		= None
            self.I_ref 	= None
            self.U_ref 	= None
            self.V_ref 	= None
            self.Q_ref 	= None
        self.__type = Type


    def __repr__(self):
        import numpy as np
        text = ['\n']
        text.append('Relation of body radiance variables:\n')
        if self.__type != 'star':
            text.append('   I = ' + str(np.shape(self.I)) + ' array [normalized]' )
            text.append('   U = ' + str(np.shape(self.U)) + ' array [normalized]' )
            text.append('   V = ' + str(np.shape(self.V)) + ' array [normalized]' )
            text.append('   Q = ' + str(np.shape(self.Q)) + ' array [normalized]' )
            text.append('   I_ref = ' + str(np.shape(self.I_ref)) + ' array [normalized]' )
            text.append('   U_ref = ' + str(np.shape(self.U_ref)) + ' array [normalized]' )
            text.append('   V_ref = ' + str(np.shape(self.V_ref)) + ' array [normalized]' )
            text.append('   Q_ref = ' + str(np.shape(self.Q_ref)) + ' array [normalized]' )
        else:
            text.append('   --No data available--')
        return '\n'.join(text)


    def __call__(self, index = 0, wvl = 0):
        import numpy as np
        text = ['\n']

        if self.__type != 'star':

            if np.ndim(self.I)== 2:

                if index > np.size(self.I,1)-1:
                    text.append('   ! Index is larger than arrays\' dimension')
                else:
                    text.append('Relation of body radiance variables at time index ' + str(index) + ' and wavelength index ' + str(wvl) + ':\n')

                    text.append('   I = ' + str(self.I[wvl,index]) + ' [-]' )
                    text.append('   U = ' + str(self.U[wvl,index]) + ' [-]' )
                    text.append('   V = ' + str(self.V[wvl,index]) + ' [-]' )
                    text.append('   Q = ' + str(self.Q[wvl,index]) + ' [-]' )
                    text.append('   I_ref = ' + str(self.I_ref[wvl,index]) + ' [-]' )
                    text.append('   U_ref = ' + str(self.U_ref[wvl,index]) + ' [-]' )
                    text.append('   V_ref = ' + str(self.V_ref[wvl,index]) + ' [-]' )
                    text.append('   Q_ref = ' + str(self.Q_ref[wvl,index]) + ' [-]' )

            else:

                text.append('Relation of body radiance variables:\n')
                text.append('   I = ' + str(self.I) + ' [-]' )
                text.append('   U = ' + str(self.U) + ' [-]' )
                text.append('   V = ' + str(self.V) + ' [-]' )
                text.append('   Q = ' + str(self.Q) + ' [-]' )
                text.append('   I_ref = ' + str(self.I_ref) + ' [-]' )
                text.append('   U_ref = ' + str(self.U_ref) + ' [-]' )
                text.append('   V_ref = ' + str(self.V_ref) + ' [-]' )
                text.append('   Q_ref = ' + str(self.Q_ref) + ' [-]' )

        else:
           text.append('Relation of body radiance variables:\n')
           text.append('   --No data available--')

        print('\n'.join(text))


class orbital_elements():
    '''
    ==================================================================
    EXOPY class: body.orbital_elements
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    The class 'orbital elements' stores all input/output orbital ele-
    ments involved in the computation of the bodies' orbits.

    METHODS
    ------------------------------------------------------------------
    - a_b: semi-major axis of the planet-moon system barycentre's or-
            bit around the star [m] (float)
    - e_b: eccentricity of the planet-moon system barycentre's orbit
            around the star [-] (float)
    - i_b: inclination of the planet-moon system barycentre's orbit
            around the star [deg] (float)
    - Omega_b: Right Ascension of the Ascending Node of the planet-
                moon system barycentre's orbit around the star [deg]
                (float)
    - omega_b: Argument of periapsis of the planet-moon system bary-
                centre's orbit around the star [deg] (float)
    - t0_b: time since last periapsis passage of the planet-moon sys-
            tem barycentre's orbit around the star [s] (float)
    - E_bs: Eccentric anomaly of the planet-moon system barycentre's
            orbit around the star [rad] (numpy array)
    - M_bs: Mean anomaly of the planet-moon system barycentre's
            orbit around the star [rad] (numpy array)
    - nu_bs: True anomaly of the planet-moon system barycentre's
            orbit around the star [rad] (numpy array)
    - a: semi-major axis of the moon's orbit around the planet-moon
        system barycentre [m] (float)
    - e: eccentricity of the moon's orbit around the planet-moon
        system barycentre [-] (float)
    - i: inclination of the moon's orbit around the planet-moon
        system barycentre [deg] (float)
    - Omega: Right Ascension of the Ascending Node of the moon's orbit
            around the planet-moon system barycentre [deg] (float)
    - omega: Argument of periapsis of the moon's orbit around the planet
            -moon system barycentre [deg] (float)
    - t0: time since last periapsis passage of the moon's orbit around
        the planet-moon system barycentre [s] (float)
    - E_mb: Eccentric anomaly of the moon's orbit around the planet
            -moon system barycentre [rad] (numpy array)
    - M_mb: Mean anomaly of the moon's orbit around the planet
            -moon system barycentre [rad] (numpy array)
    - nu_mb: True anomaly of the moon's orbit around the planet
            -moon system barycentre [rad] (numpy array)


    A relation of orbital elements values at time epoch 'i' can be re-
    trieved by body.orbital_elements(i).

    '''

    def __init__(self,Type):

        self.__type = Type

        if Type == 'moon':
                self.a     = None
                self.e     = None
                self.i     = None
                self.Omega = None
                self.omega = None
                self.t0    = None
                self.E_mb  = None
                self.M_mb  = None
                self.nu_mb = None
                self.nu_mp = None


        elif Type == 'planet':
                self.a_b     = None
                self.e_b     = None
                self.i_b     = None
                self.Omega_b = None
                self.omega_b = None
                self.t0_b    = None
                self.E_bs  = None
                self.M_bs  = None
                self.nu_bs = None


    def __call__(self, index = 0):
        import numpy as np
        text = ['\n']

        if self.__type == 'planet':
                if np.ndim(self.E_bs)== 1:

                        if index > np.size(self.E_bs,0)-1:
                                text.append('   ! Index is larger than arrays\' dimension')
                        else:
                                text.append('Relation of orbital elements at time index ' + str(index) + ':\n')

                                text.append('   Bar. semi-major axis (a_b): '    + str(self.a_b) + ' m'  )
                                text.append('   Bar. eccentricity (e_b):    '    + str(self.e_b) + '[-]' )
                                text.append('   Bar. inclination (i_b):     '    + str(self.i_b) + ' deg')
                                text.append('   Bar. argument of periapsis (omega_b):    '  + str(self.omega_b) + ' deg')
                                text.append('   Bar. right ascension of the ascending node (Omega_b): '  + str(self.Omega_b) + ' deg')
                                text.append('   Bar. time from periapsis passage (t0_b): '  + str(self.t0_b) + ' s')
                                text.append('   Bar. eccentric anomaly (E_bs) = ' + str(self.E_bs[index]) + ' rad' )
                                text.append('   Bar. mean anomaly (M_bs)  = ' + str(self.M_bs[index])  + ' rad' )
                                text.append('   Bar. true anomaly (nu_bs) = ' + str(self.nu_bs[index]) + ' rad' )

                else:

                        text.append('Relation of orbital elements:\n')

                        text.append('   Bar. semi-major axis (a_b): '    + str(self.a_b) + ' m'  )
                        text.append('   Bar. eccentricity (e_b):    '    + str(self.e_b) + '[-]' )
                        text.append('   Bar. inclination (i_b):     '  + str(self.i_b)   + ' deg')
                        text.append('   Bar. argument of periapsis (omega_b):    '  + str(self.omega_b) + ' deg')
                        text.append('   Bar. right ascension of the ascending node (Omega_b): '  + str(self.Omega_b) + ' deg')
                        text.append('   Bar. time from periapsis passage (t0_b): '  + str(self.t0_b) + ' s')
                        text.append('   Bar. eccentric anomaly (E_bs) = ' + str(np.shape(self.E_bs)) + ' array' )
                        text.append('   Bar. mean anomaly (M_bs)  = ' + str(np.shape(self.M_bs)) + ' array' )
                        text.append('   Bar. true anomaly (nu_bs) = ' + str(np.shape(self.nu_bs))+ ' array' )

        elif self.__type == 'moon':

                if np.ndim(self.E_mb)== 1:

                        if index > np.size(self.E_mb,0)-1:
                                text.append('   ! Index is larger than arrays\' dimension')
                        else:
                                text.append('Relation of orbital elements at time index ' + str(index) + ':\n')

                                text.append('   Semi-major axis (a): '  + str(self.a) + ' m'  )
                                text.append('   Eccentricity (e):    '  + str(self.e) + '[-]' )
                                text.append('   Inclination (i):     '  + str(self.i) + ' deg')
                                text.append('   Argument of periapsis (omega):    '  + str(self.omega) + ' deg')
                                text.append('   Right ascension of the ascending node (Omega): '  + str(self.Omega) + ' deg')
                                text.append('   Time from periapsis passage (t0): '  + str(self.t0) + ' s')
                                text.append('   Eccentric anomaly (E_mb) = ' + str(self.E_mb[index]) + ' rad' )
                                text.append('   Mean anomaly (M_mb)  = ' + str(self.M_mb[index])  + ' rad' )
                                text.append('   True anomaly (nu_mb) = ' + str(self.nu_mb[index]) + ' rad' )

                else:
                        text.append('Relation of orbital elements:\n')

                        text.append('   Semi-major axis (a): '  + str(self.a) + ' m'  )
                        text.append('   Eccentricity (e):    '  + str(self.e) + '[-]' )
                        text.append('   Inclination (i):     '  + str(self.i) + ' deg')
                        text.append('   Argument of periapsis (omega):    '  + str(self.omega) + ' deg')
                        text.append('   Right ascension of the ascending node (Omega): '  + str(self.Omega) + ' deg')
                        text.append('   Time from periapsis passage (t0): '  + str(self.t0) + ' s')
                        text.append('   Eccentric anomaly (E_mb) = ' + str(np.shape(self.E_mb)) + ' rad' )
                        text.append('   Mean anomaly (M_mb)  = ' + str(np.shape(self.M_mb))  + ' rad' )
                        text.append('   True anomaly (nu_mb) = ' + str(np.shape(self.nu_mb)) + ' rad' )


        else:
                text.append('Relation of body radiance variables:\n')
                text.append('   --No data available--')

        print('\n'.join(text))



    def __repr__(self):
        import numpy as np
        text = ['\n']

        if self.__type == 'planet':

                text.append('Relation of orbital elements:\n')

                text.append('   Bar. semi-major axis (a_b): '    + str(self.a_b) + ' m'  )
                text.append('   Bar. eccentricity (e_b):    '    + str(self.e_b) + '[-]' )
                text.append('   Bar. inclination (i_b):     '  + str(self.i_b)   + ' deg')
                text.append('   Bar. argument of periapsis (omega_b):    '  + str(self.omega_b) + ' deg')
                text.append('   Bar. right ascension of the ascending node (Omega_b): '  + str(self.Omega_b) + ' deg')
                text.append('   Bar. time from periapsis passage (t0_b): '  + str(self.t0_b) + ' s')
                text.append('   Bar. eccentric anomaly (E_bs) = ' + str(np.shape(self.E_bs)) + ' array' )
                text.append('   Bar. mean anomaly (M_bs)  = ' + str(np.shape(self.M_bs)) + ' array' )
                text.append('   Bar. true anomaly (nu_bs) = ' + str(np.shape(self.nu_bs))+ ' array' )

        elif self.__type == 'moon':

                text.append('Relation of orbital elements:\n')

                text.append('   Semi-major axis (a): '  + str(self.a) + ' m'  )
                text.append('   Eccentricity (e):    '  + str(self.e) + '[-]' )
                text.append('   Inclination (i):     '  + str(self.i) + ' deg')
                text.append('   Argument of periapsis (omega):    '  + str(self.omega) + ' deg')
                text.append('   Right ascension of the ascending node (Omega): '  + str(self.Omega) + ' deg')
                text.append('   Time from periapsis passage (t0): '  + str(self.t0)  + ' s')
                text.append('   Eccentric anomaly (E_mb) = ' + str(np.shape(self.E_mb)) + ' rad' )
                text.append('   Mean anomaly (M_mb)  = ' + str(np.shape(self.M_mb))  + ' rad' )
                text.append('   True anomaly (nu_mb) = ' + str(np.shape(self.nu_mb)) + ' rad' )
        else:

                text.append('Relation of body radiance variables:\n')
                text.append('   --No data available--')


        return '\n'.join(text)



class flag():
    '''
    ==================================================================
    EXOPY class: body.flag
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    The class 'flag' stores all relevant indicators on the computation
    status of the body object.

    METHODS
    ------------------------------------------------------------------
    - eclipse_d: flag indicating the computation status of eclipse
            shadow [True/False] (bool)
    - transit_d: flag indicating the computation status of transit
            shadow [True/False] (bool)
    - phase_d: flag indicating the computation status of phase shadow
            [True/False] (bool)
    - radiance: flag indicating the computation status of reflected
            starlight [True/False] (bool)
    - combine: flag indicating the computation status of radiance sig-
            nals combination [True/False] (bool)
    - orbit: flag indicating the computation status of the orbital
        geometry of the body [True/False] (bool)
    - transit: flag indicating the occurrence of transit shadowing [-]
            (list)
    - antumbra: flag indicating the occurrence of antumbra shadowing
            [-] (list)
    - umbra: flag indicating the occurrence of umbra shadowing [-]
        (bool)
    - penumbra: flag indicating the occurrence of penumbra shadowing
            [-] (list)


    A relation of flag values at time epoch 'i' can be retrieved
    by body.flag(i).

    '''

    def __init__(self, Type):

        if  Type != 'star':
                self.transit          = []
                self.antumbra         = []
                self.umbra            = []
                self.penumbra         = []
                self.transit_d        = False
                self.phase_d          = False
                self.eclipse_d        = False
                self.radiance	      = False
                self.combine	      = False
                self.orbit	      = False
                self.__type           = Type

    def __call__(self, index = 0):
        import numpy as np
        text = ['\n']
        text.append('Relation of body flag indicators at time index ' + str(index) + ':\n')
        if self.__type != 'star':

                text.append('   Orbit geometry calculated:\t'  + str(self.orbit))
                text.append('   Phase shadow calculated:\t'    + str(self.phase_d))
                text.append('   Transits shadow calculated:\t' + str(self.transit_d))
                text.append('   Eclipses shadow calculated:\t' + str(self.eclipse_d))
                text.append('   Reflected light calculated:\t' + str(self.radiance))
                text.append('   Reflected light combined:\t'   + str(self.combine) + '\n')

                for i in range(np.size(self.transit,0)):
                        existance = np.bool(self.transit[i][0][index])
                        text.append('   Transit events with ' + self.transit[i][1]        + ': ' + str(existance))

                for i in range(np.size(self.umbra,0)):
                        existance = np.bool(self.umbra[i][0][index])
                        text.append('   Umbral shadowing due to ' + self.umbra[i][1]    + ': ' + str(existance))

                for i in range(np.size(self.antumbra,0)):
                        existance = np.bool(self.antumbra[i][0][index])
                        text.append('   Antumbral shadowing due to ' + self.antumbra[i][1] + ': ' + str(existance))

                for i in range(np.size(self.penumbra,0)):
                        existance = np.bool(self.penumbra[i][0][index])
                        text.append('   Penumbral shadowing due to ' + self.penumbra[i][1] + ': ' + str(existance))

        else:
                text.append('   --No data available--')
        print( '\n'.join(text))



    def __repr__(self):
        import numpy as np
        text = ['\n']
        text.append('Relation of body flag indicators:\n')
        if self.__type != 'star':

                text.append('   Orbit geometry calculated:\t'  + str(self.orbit))
                text.append('   Phase shadow calculated:\t'    + str(self.phase_d))
                text.append('   Transits shadow calculated:\t' + str(self.transit_d))
                text.append('   Eclipses shadow calculated:\t' + str(self.eclipse_d))
                text.append('   Reflected light calculated:\t' + str(self.radiance))
                text.append('   Reflected light combined:\t'   + str(self.combine) + '\n')

                for i in range(np.size(self.transit,0)):
                        existance = np.bool(np.sum(self.transit[i][0]))
                        text.append('   Transit events with ' + self.transit[i][1]        + ': ' + str(existance))

                for i in range(np.size(self.umbra,0)):
                        existance = np.bool(np.sum(self.umbra[i][0]))
                        text.append('   Umbral shadowing due to ' + self.umbra[i][1]    + ': ' + str(existance))

                for i in range(np.size(self.antumbra,0)):
                        existance = np.bool(np.sum(self.antumbra[i][0]))
                        text.append('   Antumbral shadowing due to ' + self.antumbra[i][1] + ': ' + str(existance))

                for i in range(np.size(self.penumbra,0)):
                        existance = np.bool(np.sum(self.penumbra[i][0]))
                        text.append('   Penumbral shadowing due to ' + self.penumbra[i][1] + ': ' + str(existance))

        else:
                text.append('   --No data available--')
        return '\n'.join(text)


class grid():
    '''
    ==================================================================
    EXOPY cl 	ass: body.grid
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    The class 'grid' allows to define the disk discretization and sto-
    res all relevant information on the body grid, as well as radiance
    and geometry information wich varies pixel-to-pixel.

    METHODS
    ------------------------------------------------------------------
    - N_points: Total number of pixels in grid [-] (int)
    - Nsq: Number of pixels along the equator [-] (int)
    - area: Area of each pixel along the gris [-] (numpy array)
    - type: Type of discretization [-] (str)
    - nodes: 2D Cartesian coordinates of each pixel in grid reference
        frame [-] (numpy array)
    - nodes_xyz: 3D Cartesian coordinates of each pixel in grid refe-
            rence frame [-] (numpy array)
    - nodes_xyz_rot: Rotated 3D Cartesian coordinates of each pixe
            in grid reference frame [-] (numpy array)
    - I: First Stokes element: reflected starlight flux per pixel
        [normalized] (numpy array)
    - U: Second Stokes element: linear polarization of reflected
        starlight per pixel [normalized] (numpy array)
    - V: Third Stokes element: circular polarization of reflected
        starlight per pixel [normalized] (numpy array)
    - Q: Fourth Stokes element: circular polarization of reflected
        starlight per pixel [normalized] (numpy array)
    - faces: Cartesian points delimiting the vertices of each pixel
        [-] (numpy array)
    - azimuth: Azimuth angle of each pixel at each time epoch [rad]
            (numpy array)
    - beta: Beta angle of each pixel at each time epoch [rad] (numpy
        array)
    - distance_nodes_ob: Distance from star to each pixel [m] (numpy
                array)
    - illuminated_nodes: Array indicating the illumination status of
                each pixel at each time epoch [True/False]
                (numpy array)
    - observer_zenith_angle: Angle between the observer position,
                the centre of the pixel and the zenith
                direction [rad] (numpy array)
    - solar_zenith_angle: Angle between the star and the zenith di-
                rection, with centre at the pixel centre
                for each time epoch [rad] (numpy array)
    - position_nodes_ob: 3D Cartesian coordinates of each pixel in
                observer reference frame [m] (numpy array)
    - set_grid: Function for the definition of the grid.
    - shadow: Shadowing status of each pixel at each time expressed
        from 0 to 1.
    - show_grid: Function for plotting an overview of the grid used.

    A relation of grid values at time epoch 'i' can be retrieved
    by body.grid(i).

    '''

    def __init__(self,Type):

        self.Nsq      = 100
        self.type     = 'Square'
        self.I	      = None
        self.Q	      = None
        self.U	      = None
        self.V	      = None
        self.azimuth  = None
        self.beta     = None
        self.distance_nodes_ob     = None
        self.illuminated_nodes     = None
        self.nodes_xyz_rot         = None
        self.observer_zenith_angle = None
        self.position_nodes_ob     = None
        self.solar_zenith_angle	   = None
        self.__type 		   = Type
        self.nodes, self.faces, self.area, self.N_points, self.nodes_xyz = grd.square(self.Nsq)

    def __call__(self,index=0, wvl=0):
        import numpy as np
        text = ['\n']

        if self.__type != 'star':

                if np.ndim(self.illuminated_nodes) == 2:

                        if index > np.size(self.illuminated_nodes,0)-1:
                                text.append('   ! Index is larger than arrays\' dimension')
                        else:
                                text.append('Relation of body grid variables at time index ' + str(index) + ' and wavelength index ' + str(wvl) + ':\n')

                                text.append('   Grid type:\t'     + str(self.type))
                                text.append('   Grid points:\t'     + str(self.N_points))
                                text.append('   Grid eq. points: '  + str(self.Nsq) + '\n')
                                text.append('   I = ' + str(np.shape(self.I[wvl,index])) + ' subarray [normalized]' )
                                text.append('   Q = ' + str(np.shape(self.Q[wvl,index])) + ' subarray [normalized]' )
                                text.append('   U = ' + str(np.shape(self.U[wvl,index])) + ' subarray [normalized]' )
                                text.append('   V = ' + str(np.shape(self.V[wvl,index])) + ' subarray [normalized]' )
                                text.append('   azimuth = ' + str(np.shape(self.azimuth )) + ' array [rad]' )
                                text.append('   beta    = ' + str(np.shape(self.beta    )) + ' array [rad]' )
                                text.append('   distance_nodes_ob = ' + str(np.shape(self.distance_nodes_ob[index])) + ' array [m]' )
                                text.append('   illuminated_nodes = ' + str(np.shape(self.illuminated_nodes[index])) + ' array [bool]' )
                                text.append('   nodes_xyz_rot = ' + str(np.shape(self.nodes_xyz_rot)) + ' array [normalized]' )
                                text.append('   observer_zenith_angle = ' + str(np.shape(self.observer_zenith_angle)) + ' array [rad]' )
                                text.append('   nodes = ' + str(np.shape(self.nodes)) + ' array [normalized]' )
                                text.append('   faces = ' + str(np.shape(self.faces)) + ' array [normalized]' )
                                text.append('   area = ' + str(np.shape(self.area)) + ' array [normalized]' )
                                text.append('   nodes_xyz = ' + str(np.shape(self.nodes_xyz)) + ' array [normalized]\n' )

                                text.append('   show_grid = Grid plot function')
                                text.append('   set_grid = Grid definition function')

                else:

                        text.append('Relation of body grid variables:\n')

                        text.append('   Grid type:\t'     + str(self.type))
                        text.append('   Grid points:\t'     + str(self.N_points))
                        text.append('   Grid eq. points: '  + str(self.Nsq) + '\n')
                        text.append('   I = ' + str(np.shape(self.I)) + ' subarray [normalized]' )
                        text.append('   Q = ' + str(np.shape(self.Q)) + ' subarray [normalized]' )
                        text.append('   U = ' + str(np.shape(self.U)) + ' subarray [normalized]' )
                        text.append('   V = ' + str(np.shape(self.V)) + ' subarray [normalized]' )
                        text.append('   azimuth = ' + str(np.shape(self.azimuth )) + ' array [rad]' )
                        text.append('   beta    = ' + str(np.shape(self.beta    )) + ' array [rad]' )
                        text.append('   distance_nodes_ob = ' + str(np.shape(self.distance_nodes_ob)) + ' array [m]' )
                        text.append('   illuminated_nodes = ' + str(np.shape(self.illuminated_nodes)) + ' array [bool]' )
                        text.append('   nodes_xyz_rot = ' + str(np.shape(self.nodes_xyz_rot)) + ' array [normalized]' )
                        text.append('   observer_zenith_angle = ' + str(np.shape(self.observer_zenith_angle)) + ' array [rad]' )
                        text.append('   nodes = ' + str(np.shape(self.nodes)) + ' array [normalized]' )
                        text.append('   faces = ' + str(np.shape(self.faces)) + ' array [normalized]' )
                        text.append('   area = ' + str(np.shape(self.area)) + ' array [normalized]' )
                        text.append('   nodes_xyz = ' + str(np.shape(self.nodes_xyz)) + ' array [normalized]\n' )

                        text.append('   show_grid = Grid plot function')
                        text.append('   set_grid = Grid definition function')

        else:
                text.append('Relation of body grid variables:\n')
                text.append('   --No data available--')

        print('\n'.join(text))


    def __repr__(self):
        import numpy as np
        text = ['\n']
        text.append('Relation of body grid variables:\n')
        if self.__type != 'star':

            text.append('\n   Grid type:\t'     + str(self.type))
            text.append('   Grid points:\t'     + str(self.N_points))
            text.append('   Grid eq. points: '  + str(self.Nsq) + '\n')
            text.append('   I = ' + str(np.shape(self.I)) + ' array [normalized]' )
            text.append('   Q = ' + str(np.shape(self.Q)) + ' array [normalized]' )
            text.append('   U = ' + str(np.shape(self.U)) + ' array [normalized]' )
            text.append('   V = ' + str(np.shape(self.V)) + ' array [normalized]' )
            text.append('   azimuth = ' + str(np.shape(self.azimuth )) + ' array [rad]' )
            text.append('   beta    = ' + str(np.shape(self.beta    )) + ' array [rad]' )
            text.append('   distance_nodes_ob = ' + str(np.shape(self.distance_nodes_ob)) + ' array [m]' )
            text.append('   illuminated_nodes = ' + str(np.shape(self.illuminated_nodes)) + ' array [bool]' )
            text.append('   nodes_xyz_rot = ' + str(np.shape(self.nodes_xyz_rot)) + ' array [normalized]' )
            text.append('   observer_zenith_angle = ' + str(np.shape(self.observer_zenith_angle)) + ' array [rad]' )
            text.append('   nodes = ' + str(np.shape(self.nodes)) + ' array [normalized]' )
            text.append('   faces = ' + str(np.shape(self.faces)) + ' array [normalized]' )
            text.append('   area = ' + str(np.shape(self.area)) + ' array [normalized]' )
            text.append('   nodes_xyz = ' + str(np.shape(self.nodes_xyz)) + ' array [normalized]\n' )

            text.append('   show_grid = Grid plot function')
            text.append('   set_grid = Grid definition function')

        else:
            text.append('   --No data available--')
        return '\n'.join(text)


    def show_grid(self):
        '''
        ==================================================================
        EXOPY class: body.grid.show_grid()
        Delft University of Technology
        ------------------------------------------------------------------
        Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
        Date: 2016-2017
        ------------------------------------------------------------------

        DESCRIPTION
        ------------------------------------------------------------------
        The function 'show_grid' allows the user to easily plot an over-
        view of the grid defined on a certain body.

        INPUT
        ------------------------------------------------------------------
        - No inputs -

        '''

        fig , ax = plt.subplots()
        #plt.plot(self.nodes[:,0], self.nodes[:,1], 'or')
        for i in range(self.N_points):
            plt.plot(self.faces[i,0,:], self.faces[i,1,:], 'b-', linewidth = '1')
            circle1 = plt.Circle((0, 0), 0.5, color = 'k', fill=False, zorder=1)
            ax.add_artist(circle1)
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim([-0.55, 0.55])
            ax.set_ylim([-0.55, 0.6])
            ax.text(-0.45,  0.525,'Points: '+ str(self.N_points), fontsize= 15)

    #def show_shadow(self,t):

    #    polygon = list(np.zeros(self.N_points))
    #    p   = list(np.zeros(self.N_points))

    #    fig , ax = plt.subplots()
    #    plt.plot(self.nodes[:,0], self.nodes[:,1], 'or')
    #    for i in range(self.N_points):
    #        plt.plot(self.faces[i,0,:], self.faces[i,1,:], 'b-', linewidth = '2',zorder=2)
    #        polygon[i] = PolygonPatch(Poly(self.faces[i,:,:].T))
    #        p[i] = PatchCollection([polygon[i]], alpha=1, color=str(int(self.shadow[t,i])),edgecolor='none', zorder=1)
    #        ax.add_collection(p[i])

    #    circle1 = plt.Circle((0, 0), 0.5, color = 'r', fill=False, zorder=3)
    #    ax.add_artist(circle1)
    #    ax.set_aspect('equal', adjustable='box')
    #    ax.set_xlim([-0.55, 0.55])
    #    ax.set_ylim([-0.55, 0.6])
    #    ax.text(-0.4,  0.525,'Points: '+ str(self.N_points), fontsize= 15)
    #    ax.text(0.14,  0.525,'$\^{A}$: '+ str(np.round(np.sum(self.area)/(m.pi*0.5**2),3)), fontsize= 15)



    def set_grid(self, grid_type='square', Nlon=7, Nlat=7, cos1 = 80, cos2 = 60, diff = 0.009, Nsq = 7, Nang=15, Nrad=3, max_area=0.03, min_angle=30, circle_edges = 20, centre = False):
    	'''
        ==================================================================
        EXOPY class: body.grid.set_grid()
        Delft University of Technology
        ------------------------------------------------------------------
        Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
        Date: 2016-2017
        ------------------------------------------------------------------

        DESCRIPTION
        ------------------------------------------------------------------
        The function 'set_grid' allows the user to discretize the disk of
        the body.

        INPUT
        ------------------------------------------------------------------
        - grid_type: Type of grid ('square', 'radial', 'triangular',
                'sphere')
        - Nlon:
        - Nlat:
        - cos1:
        - cos2:
        - diff:
        - Nsq: ['square'] Number of pixels along the equator [-] (int)
        - Nang:
        - Nrad:
        - max_area:
        - min_angle:
        - circle_edges:
        - centre:

            TO BE COMPLETED

    	'''

        if grid_type == 'triangular':
            self.nodes, self.faces, self.area, self.N_points, self.nodes_xyz = grd.triangl(max_area, min_angle, circle_edges, centre)
            print('  Triangular grid created with Meshpy')
        elif grid_type == 'square':
            self.nodes, self.faces, self.area, self.N_points, self.nodes_xyz = grd.square(Nsq)
            self.Nsq  =  Nsq
            print('  Square grid created!')
        elif grid_type == 'radial':
            self.nodes, self.faces, self.area, self.N_points = grd.radial(Nang, Nrad)
            print('  Radial grid created!')
        elif grid_type == 'sphere':
            self.nodes, self.faces, self.area, self.N_points = grd.sphere(Nlon, Nlat, cos1, cos2, diff)
            print('  Spherical grid created!')
        else:
            sys.exit('ERROR: Grid type not recognized!')


class ProgressBar(object):
    '''
    ==================================================================
    EXOPY class: ProgressBar
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    Provides a progress bar to be shown during computations.

    '''

    DEFAULT = '    Progress: %(bar)s %(percent)3d%%'
    FULL = '    %(bar)s %(current)d/%(total)d (%(percent)3d%%) %(remaining)d to go'

    def __init__(self, total, width=20, fmt=DEFAULT, symbol='█',
                 output=sys.stderr):
        assert len(symbol) == 1

        self.total = total
        self.width = width
        self.symbol = symbol
        self.output = output
        self.fmt = re.sub(r'(?P<name>%\(.+?\))d',
            r'\g<name>%dd' % len(str(total)), fmt)

        self.current = 0

    def __call__(self):
        percent = self.current / float(self.total)
        size = int(self.width * percent)
        remaining = self.total - self.current
        bar = '|' + self.symbol * size + ' ' * (self.width - size) + '|'

        args = {
            'total': self.total,
            'bar': bar,
            'current': self.current,
            'percent': percent * 100,
            'remaining': remaining
        }
        print('\r' + self.fmt % args, file=self.output, end='')

    def done(self):
        self.current = self.total
        self()
        print('', file=self.output)


class Arrow3D(FancyArrowPatch):
    '''
    ==================================================================
    EXOPY class: Arrow3D
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    Class for the creation of three dimensional arrow objects, based
    on [4]

    [4] stackoverflow.com/questions/29188612/arrows-in-matplotlib-using-mplot3d


        '''

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
#--------------------------   FUNCTIONS   --------------------------------'''
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''




def anim_to_html(anim):
    '''
    ==================================================================
    EXOPY function: anim_to_html
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    Converts a python animation object to a HTML format compatible with
    Jupyter Notebook. Based on [1]

    [1] jakevdp.github.io/blog/2013/05/12/embedding-matplotlib-animations/

    INPUTS
    ------------------------------------------------------------------
    - anim: Animation object.

    OUTPUTS
    ------------------------------------------------------------------
    - Encoded video animation.


    '''

    if not hasattr(anim, '_encoded_video'):
        with NamedTemporaryFile(suffix='.mp4') as f:
            anim.save(f.name, fps=20, extra_args=['-vcodec', 'libx264'])
            video = open(f.name, "rb").read()
            anim._encoded_video = base64.b64encode(video).decode('utf-8')

    return VIDEO_TAG.format(anim._encoded_video)


def display_animation(anim):
    '''
    ==================================================================
    EXOPY function: anim_to_html
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    Displays an animation in Jupyter Notebook. Based on [1].

    [1] jakevdp.github.io/blog/2013/05/12/embedding-matplotlib-animations/

    INPUTS
    ------------------------------------------------------------------
    - anim: Animation object.

    OUTPUTS
    ------------------------------------------------------------------
    - HTML video.


    '''

    plt.close(anim._fig)
    #return HTML(anim_to_html(anim))
    return HTML(anim.to_html5_video())



def plot_xy(x, y1, y2 = None, info = ['x [-]', 'y [-]', ' ', None]):
    '''
    ==================================================================
    EXOPY function: plot_xy
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    Predefined function for 2D plots, e.g. phase vs time.

    INPUTS
    ------------------------------------------------------------------
    - x: Horizontal coordinate.
    - y1: Vertical coordinate.
    - y2: Second vertical coordinate.
    - info :Information on lables
        0. x-lable (str)
        1. y-lable (str)
        2. title   (str)
        3. legend  [(str),(str)]

    OUTPUTS
    ------------------------------------------------------------------
    - Figure.

    '''

    plt.figure()
    plt.plot(x,y1,linewidth=2, linestyle="-", label="Curve 1")
    if y2 is not None: # If a second curve is required...
        plt.plot(x,y2,linewidth=2, linestyle="-", label="Curve 2")
    plt.xlabel(info[0])
    plt.ylabel(info[1])
    plt.xlim([min(x),max(x)])
    plt.title(info[2])
    if y2 is not None:
        if info[3] is not None:
            plt.legend([info[3][0],info[3][1]], loc='upper left',frameon=False)
        else:
            plt.legend(loc='upper left', frameon=False)
    else:
        if info[3] is not None:
            plt.legend([info[3]], loc='upper left', frameon=False)
    plt.grid()
    plt.show()


def query_yes_no(question, default="no"):
    '''
    ==================================================================
    EXOPY function: plot_xy
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    Asks a yes/no question via input() and return their answer, based
    on [3]

    [3] Recipe 577058 from activestate.com. Adapted to Python 3.

    INPUTS
    ------------------------------------------------------------------
    - question: String to be presented to the user [-] (str)
    - default: Presumed answer if user presses <Enter>. Values admi-
            tted: 'yes', 'no' (default), None (meaning an answer
            is required from the side of the user) [-] (str)

    OUTPUTS
    ------------------------------------------------------------------
    - True or False for positive and negative answers accordingly.

        '''

    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        aux = question+prompt
        sys.stdout.write(aux)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")



def grid_area(grid,t):
    '''
    ==================================================================
    EXOPY function: grid_area
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    -- TBC --

    INPUTS
    ------------------------------------------------------------------
    - grid: Grid type of object [-] (grid class)
    - t: time epoch index [-] (int)

    OUTPUTS
    ------------------------------------------------------------------
    - area: -- TBC --

    '''


    if t== 'all':
        area = np.degrees(4*np.einsum('i,ji->j',grid.area,grid.shadow!=1))
    else:
        area = np.degrees(4*np.dot(grid.area,grid.shadow[t,:]!=1))

    return area


def symp(x, deg = False):
    '''
    ==================================================================
    EXOPY function: symp
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    DESCRIPTION
    ------------------------------------------------------------------
    Converts an input angle to the range 0, 2pi.

    INPUTS
    ------------------------------------------------------------------
    - x: Input angle [rad/deg] (float)
    - deg: True if x given in degrees, False if x given in radians
        [True/False] (bool)

    OUTPUTS
    ------------------------------------------------------------------
    - x: Converted angle

    '''


    if deg is False:
        x = x * 180 / m.pi
    if x > 360:
        n = x // 360
        x = x - 360 * n

    elif x < 0:
        n = abs(x) // 360
        x = x + 360 * (n+1)
    if deg is False:
        x = x * m.pi / 180
    return x




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
#-------------------------- End of script --------------------------------'''
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''

'''

def symp(x, deg = False):
    """
    ==============================================================================
       October 2016, Javier B.M., TU Delft
    ------------------------------------------------------------------------------
       Convert an input angle to the range 0, 2pi.

       Inputs:
           - x: Input angle.
           - deg: True,  if x is given in degrees.
                  False, if x is given in radians. (Default)
           - y2: Second vertical coordinate.

       Outputs:
           - x: Output angle.
    ==============================================================================
    """

    if deg is False:
        x = x * 180 / m.pi
    if x > 360:
        n = x // 360
        x = x - 360 * n

    elif x < 0:
        n = abs(x) // 360
        x = x + 360 * (n+1)
    if deg is False:
        x = x * m.pi / 180
    return x


def symp1(x, deg = False):
    """
    ==============================================================================
       October 2016, Javier B.M., TU Delft
    ------------------------------------------------------------------------------
       Convert an input angle to the range -pi, pi.

       Inputs:
           - x: Input angle.
           - deg: True,  if x is given in degrees.
                  False, if x is given in radians. (Default)
           - y2: Second vertical coordinate.

       Outputs:
           - x: Output angle.
    ==============================================================================
    """

    if deg is False:
        x = x * 180 / m.pi

    while x < -180:
        x += 360
    while x > 180:
        x -= 360

    if deg is False:
        x = x * m.pi / 180
    return x






def PolyArea(x,y): # Shoelace formula
    """
    ==============================================================================
       October 2016, Javier B.M., TU Delft
    ------------------------------------------------------------------------------
       Calculates the area inside a simple non-intersecting polygon, based on the
        formula [2].

       Inputs:
           - anim: Animation object.
       Outputs:
           - HTML video.

       [2] https://en.wikipedia.org/wiki/Shoelace_formula
    ==============================================================================
    """
    return 0.5*np.abs(np.einsum('tN,tN->t', x, np.roll(y,1,axis=1))-np.einsum('tN,tN->t', y, np.roll(x,1,axis=1)))

'''
