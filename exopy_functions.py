
# -*- coding: utf-8 -*-

#==============================================================================
#                               MT_FUNCTIONS.PY
#==============================================================================
#   October 2016
#   Javier Berzosa Molina
#   Delf University of Technology
#   Astrodynamics & Space Missions and Planetary Exploration
#------------------------------------------------------------------------------
#   Module containing the common functions which are accessed by the rest of
#   the modules.
#
#   List of functions:
#       - anim_to_html: Converts animation to HTML for compatibility with
#                       jupyter notebook.
#       - display_animation: Displays animation in jupyter notebooko
#       - PolyArea: Calculates the area inside a non-intersecting polygon using
#                   the Shoelace formula.
#       - plot_xy: Predefined function for 2D plots, e.g. phase vs time.
#       - symp: Convert an input angle to the range 0, 2pi.
#       - quesry_yes_no: Asks the user a yes/no question.
#
#   List of classes:
#
#       - body: Class for the creation of planet, moon and star objects.
#       - Arrow3D: Class for the creation of three dimensional arrow objects.
#
#   References:
#   [1] jakevdp.github.io/blog/2013/05/12/embedding-matplotlib-animations/
#   [2] en.wikipedia.org/wiki/Shoelace_formula#
#   [3] Recipe 577058 from activestate.com. Adapted to Python 3.
#   [4] stackoverflow.com/questions/29188612/arrows-in-matplotlib-using-mplot3d
#
#==============================================================================
#                                                      ___
#                                                   ,o88888
#                                                ,o8888888'
#                          ,:o:o:oooo.        ,8O88Pd8888"
#                      ,.::.::o:ooooOoOoO. ,oO8O8Pd888'"
#                    ,.:.::o:ooOoOoOO8O8OOo.8OOPd8O8O"
#                   , ..:.::o:ooOoOOOO8OOOOo.FdO8O8"
#                  , ..:.::o:ooOoOO8O888O8O,COCOO"
#                 , . ..:.::o:ooOoOOOO8OOOOCOCO"
#                  . ..:.::o:ooOoOoOO8O8OCCCC"o
#                     . ..:.::o:ooooOoCoCCC"o:o
#                     . ..:.::o:o:,cooooCo"oo:o:
#                  `   . . ..:.:cocoooo"'o:o:::'
#                  .`   . ..::ccccoc"'o:o:o:::'
#                 :.:.    ,c:cccc"':.:.:.:.:.'
#               ..:.:"'`::::c:"'..:.:.:.:.:.'
#             ...:.'.:.::::"'    . . . . .'
#            .. . ....:."' `   .  . . ''
#          . . . ...."'
#          .. . ."'
#         .
#==============================================================================

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
import pymiedap as pmd

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
    ==============================================================================
    October 2016, Javier B.M., TU Delft
    ------------------------------------------------------------------------------
    Class for the creation of 'planet', 'moon', and 'star' objects.

    Internal variables:
    - body.__track: Numerated list of created bodies.
    Subclasses:
        - properties: contains mass and radius attributes
        - ephemeris: contains the position of each body for each epoch
        - geometry: contains the angles defining the geometry for each epoch
        - orbital_elements: contains the orbital elements describing each orbit
    Functions:
        - body.documentation(): Shows documentation.\n\
        - body.show(): Shows list of bodies in the database.\n\
        - body.reset(): Clears the bodies database.\n\n')

    More info available at the class documentation: >> body.documentation()
    ==============================================================================
    '''
    __track       = [['#', 'Variable name', 'Type of body']]

    def __init__(self, name, Type):

        body.__updatelist(self,name)

        '''The new body is characterized by the setup function'''
        body.__setup(self,name,Type)
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

    def __setup(self,name,Type):
        self.name = name
        self.type = Type
        body.__track.append([len(body.__track), name, Type])

        self.properties = properties()
        self.ephemeris = ephemeris()
        self.orbital_elements = orbital_elements()
        self.geometry = geometry()
        self.flag = flag()
        self.grid = grid()
        self.radiance = radiance()
        self.atmosphere = pmd.Model()

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


class settings():

    def __init__(self):

        self.trial = None


class properties():

    def __init__(self):

        self.m   = None
        self.R   = None
	self.fourier_scene = 'clear'

    def __call__(self):

        print('\nRelation of body properties:')

        print('   Mass (m) = ', self.m, ' kg')
        print('   Radius (R) = ', self.R, 'm')


class ephemeris():

    def __init__(self):

        self.time   = [None]

    def __call__(self, index = 0):

        print('\nRelation of body ephemeris:')

        print('   Time (time) = ', self.time[index], 's')


class geometry():

    def __init__(self):

        self.phase   = [None]

    def __call__(self, index = 0):

        print('\nRelation of body geometry variables:')

        print('   Phase (phase) = ', self.phase[index], ' deg')


class radiance():

    def __init__(self):

	self.I = [None]
	self.U = [None]
	self.V = [None]
	self.Q = [None]


class orbital_elements():

    def __init__(self):

        self.a     = None
        self.e     = None
        self.i     = None
        self.Omega = None
        self.omega = None

    def __call__(self):

        print('\nRelation of orbital elements:')

        print('   Semi-major axis (a) = ', self.a, 'm')
        print('   Eccentricity (e) = ', self.e)
        print('   Inclination (i) = ', self.i, ' deg')
        print('   Right ascension of the ascending node (Omega) = ', self.Omega, ' deg')
        print('   Argument of periapsis (omega) = ', self.omega, ' deg')


class flag():

    def __init__(self):

        self.eclipse          = []
        self.transit          = []
        self.antumbra         = []
        self.umbra            = []
        self.penumbra         = []
        self.transit_d        = False
        self.transit_c        = False
        self.phase_d          = False
        self.phase_c          = False
        self.eclipse_d        = False
        self.eclipse_c        = False



    def __call__(self, index = 0):

        print('\nRelation of body flags:')

        print('   Eclipse exists (eclipse) = ', self.eclipse[index])
        print('   Transit 1 exists (transit1) = ', self.transit1[index])
        print('   Transit 2 exists (transit2) = ', self.transit2[index])


class grid():

    def __init__(self):

        self.nodes = None
        self.faces = None
        self.area  = None
        self.N_points = None

        self.nodes, self.faces, self.area, self.N_points, self.nodes_xyz = grd.square(6)

    def __call__(self):

        self.show_grid()


    def show_grid(self):

        fig , ax = plt.subplots()
        plt.plot(self.nodes[:,0], self.nodes[:,1], 'or')
        for i in range(self.N_points):
            plt.plot(self.faces[i,0,:], self.faces[i,1,:], 'b-', linewidth = '2')
        circle1 = plt.Circle((0, 0), 0.5, color = 'k', fill=False, zorder=1)
        ax.add_artist(circle1)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim([-0.55, 0.55])
        ax.set_ylim([-0.55, 0.6])
        ax.text(-0.4,  0.525,'Points: '+ str(self.N_points), fontsize= 15)
        ax.text(0.14,  0.525,'$\^{A}$: '+ str(np.round(np.sum(self.area)/(m.pi*0.5**2),3)), fontsize= 15)

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

        if grid_type == 'triangular':
            self.nodes, self.faces, self.area, self.N_points, self.nodes_xyz = grd.triangl(max_area, min_angle, circle_edges, centre)
            print('  Triangular grid created with Meshpy')
        elif grid_type == 'square':
            self.nodes, self.faces, self.area, self.N_points, self.nodes_xyz = grd.square(Nsq)
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
    """
    #==============================================================================
    #   October 2016, Javier B.M., TU Delft
    #------------------------------------------------------------------------------
    #   Class for the creation of three dimensional arrow objects, based on [4]
    #
    #   [4] stackoverflow.com/questions/29188612/arrows-in-matplotlib-using-mplot3d
    #==============================================================================
    """

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
    """
    ==============================================================================
       October 2016, Javier B.M., TU Delft
    ------------------------------------------------------------------------------
       Converts a python animation object to a HTML format compatible with
       Jupyter Notebook. Based on [1].

       Inputs:
           - anim: Animation object.
       Outputs:
           - Encoded video animation.

       [1] jakevdp.github.io/blog/2013/05/12/embedding-matplotlib-animations/
    ==============================================================================
    """

    if not hasattr(anim, '_encoded_video'):
        with NamedTemporaryFile(suffix='.mp4') as f:
            anim.save(f.name, fps=20, extra_args=['-vcodec', 'libx264'])
            video = open(f.name, "rb").read()
            anim._encoded_video = base64.b64encode(video).decode('utf-8')

    return VIDEO_TAG.format(anim._encoded_video)


def display_animation(anim):
    """
    ==============================================================================
       October 2016, Javier B.M., TU Delft
    ------------------------------------------------------------------------------
       Displays an animation in Jupyter Notebook. Based on [1].

       Inputs:
           - anim: Animation object.
       Outputs:
           - HTML video.

       [1] jakevdp.github.io/blog/2013/05/12/embedding-matplotlib-animations/
    ==============================================================================
    """

    plt.close(anim._fig)
    #return HTML(anim_to_html(anim))
    return HTML(anim.to_html5_video())


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


def plot_xy(x, y1, y2 = None, info = ['x [-]', 'y [-]', ' ', None]):
    """
        # X axis label, Y axis label, Title, Legend
    ==============================================================================
       October 2016, Javier B.M., TU Delft
    ------------------------------------------------------------------------------
       Predefined function for 2D plots, e.g. phase vs time.

       Inputs:
           - x: Horizontal coordinate.
           - y1: Vertical coordinate.
           - y2: Second vertical coordinate.
           - info = (0) x-label
                    (1) y-label
                    (2) title
                    (3) legend
       Outputs:
           - Plot.
    ==============================================================================
    """

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


def query_yes_no(question, default="no"):
    """
    ==============================================================================
    October 2016, Javier B.M., TU Delft
    ------------------------------------------------------------------------------
    Ask a yes/no question via input() and return their answer, based on [3]

    Inputs:
        - question: String to be presented to the user.
        - default: Presumed answer if user presses <Enter>. Values
                    admitted: 'yes', 'no' (default), None (meaning an answer is
                    required of the user).
    Outputs:
        - True or False for positive and negative answers accordingly.

    [3] Recipe 577058 from activestate.com. Adapted to Python 3.
    ==============================================================================
    """

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

    if t== 'all':
        area = np.degrees(4*np.einsum('i,ji->j',grid.area,grid.shadow!=1))
    else:
        area = np.degrees(4*np.dot(grid.area,grid.shadow[t,:]!=1))

    return area


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
#-------------------------- End of script --------------------------------'''
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'''
