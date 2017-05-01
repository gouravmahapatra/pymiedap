# -*- coding: utf-8 -*-

"""
==================================================================
EXOPY module: exopy_plot.py
Delft University of Technology
------------------------------------------------------------------
Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
Date: 2016-2017
------------------------------------------------------------------

Dependences:

DESCRIPTION
------------------------------------------------------------------
Script containing the functions pre-defined to display the results
from EXOPY.

LIST OF FUNCTIONS
------------------------------------------------------------------
 - -- TBC --


"""

import matplotlib as _matplotlib
import matplotlib.pyplot as _plt
from matplotlib.patches import Polygon as _Polygon
from matplotlib.collections import PatchCollection as _PatchCollection
from matplotlib import animation as _animation
#from shapely.geometry import Polygon as _Poly
import numpy as _np
from exopy_functions import grid_area as _grid_area
import time as _time
import exopy_config as _cfg


def IQ2(bodies, t = 0, wvl=0, phase = False, save = False):
    from matplotlib.patches import Rectangle
    # create a figure with subplots
    fig = _plt.figure(figsize=(15,7))
    ax1 = _plt.subplot2grid((2,4), (0,0), colspan=2)
    ax2 = _plt.subplot2grid((2,4), (1,0), sharex=ax1, colspan=2)
    ax3 = _plt.subplot2grid((2,4), (0,2))
    ax5 = _plt.subplot2grid((2,4), (0,3))
    ax4 = _plt.subplot2grid((2,4), (1,2))
    ax6 = _plt.subplot2grid((2,4), (1,3))

    if phase==False:
        xlabelstr = 'Time [Earth days]'
    else:
        xlabelstr = 'Phase angle [deg]'
	
    I0ref = max(abs(bodies[0].radiance.I_ref[wvl,:]))
    I1ref = max(abs(bodies[1].radiance.I_ref[wvl,:]))
    Q0ref = max(abs(bodies[0].radiance.Q_ref[wvl,:]))
    Q1ref = max(abs(bodies[1].radiance.Q_ref[wvl,:]))

    ax1.grid()
    ax1.set_title('Flux')
    ax1.set_xlabel(xlabelstr)
    ax1.set_ylabel('$\hat I~[-]$')
    I0 = bodies[0].radiance.I[wvl,:]#/I0ref
    I0[_np.isnan(I0)] = 0
    I1 = bodies[1].radiance.I[wvl,:]#/I1ref
    I1[_np.isnan(I1)] = 0
    if phase is False:
        ax1.plot(bodies[0].ephemeris.time/3600/24, I0,'-k', linewidth = 1.5, label = bodies[0].name)
        ax1.plot(bodies[1].ephemeris.time/3600/24, I1,'k' , dashes=[10, 5, 10, 5], linewidth = 1.5, label = bodies[1].name)
    else:
        ax1.plot(_np.rad2deg(bodies[0].geometry.phase_angle), I0,'-k', linewidth = 1.5, label = bodies[0].name)
        ax1.plot(_np.rad2deg(bodies[1].geometry.phase_angle), I1,'k' , dashes=[10, 10, 10, 10], linewidth = 1.5, label = bodies[1].name)

    #ax1.legend(loc=1)

    ax2.grid()
    ax2.set_title('Q parameter')
    ax2.set_xlabel(xlabelstr)
    ax2.set_ylabel('$\hat Q~[-]$')
    Q0 = bodies[0].radiance.Q[wvl,:]#/Q0ref
    Q0[_np.isnan(Q0)] = 0
    Q1 = bodies[1].radiance.Q[wvl,:]#/Q1ref
    Q1[_np.isnan(Q1)] = 0    
    if phase is False:
        ax2.plot(bodies[0].ephemeris.time/3600/24, Q0,'-k', linewidth = 1.5, label = bodies[0].name)
        ax2.plot(bodies[1].ephemeris.time/3600/24, Q1,'k' , dashes=[10, 5, 10, 5], linewidth = 1.5, label = bodies[1].name)
    else:
        ax2.plot(_np.rad2deg(bodies[0].geometry.phase_angle), Q0,'-k', linewidth = 1.5, label = bodies[0].name)
        ax2.plot(_np.rad2deg(bodies[1].geometry.phase_angle), Q1,'k' , dashes=[10, 5, 10, 5], linewidth = 1.5, label = bodies[1].name)

    ax2.legend(loc=4)



    circle1 = _plt.Circle((0, 0), 1, facecolor = _plot_color('circle1'), edgecolor = _plot_color('circle'), linewidth = 2.5)
    circle2 = _plt.Circle((0, 0), 1, edgecolor = _plot_color('circle'), facecolor = _plot_color('circle1'), fill=True , linewidth = 2.5)
    circle3 = _plt.Circle((0, 0), 1, edgecolor = _plot_color('circle'), facecolor = _plot_color('circle1'), fill=True , linewidth = 2.5)
    circle4 = _plt.Circle((0, 0), 1, edgecolor = _plot_color('circle'), facecolor = _plot_color('circle1'), fill=True , linewidth = 2.5)
    ax3.add_artist(circle1)
    ax4.add_artist(circle2)
    ax5.add_artist(circle3)
    ax6.add_artist(circle4)

    #faces = _plot_color('faces')
    faces = "none"

    N 	= _np.sum(bodies[0].grid.shadow[t,:]!=0)
    faces00	= bodies[0].grid.faces[bodies[0].grid.shadow[t,:]!=0,0,0]
    faces10	= bodies[0].grid.faces[bodies[0].grid.shadow[t,:]!=0,1,0]
    faces01	= bodies[0].grid.faces[bodies[0].grid.shadow[t,:]!=0,0,1]
    I0	= bodies[0].grid.I[wvl,t,bodies[0].grid.shadow[t,:]!=0]#/I0ref
    Q0	= bodies[0].grid.Q[wvl,t,bodies[0].grid.shadow[t,:]!=0]#/Q0ref
    
    patches = []
    for i in range(N):
        square = Rectangle( (faces00[i]*2,faces10[i]*2),-2*faces00[i]+2*faces01[i],-2*faces00[i]+2*faces01[i] )
        patches.append(square)

    p1 = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)
    p2 = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)

    p1.set_array(I0)
    p2.set_array(Q0)
    aux3 = ax3.add_collection(p1)
    fig.colorbar(aux3,ax=ax3,fraction=0.046, pad=0.04)
    aux4 = ax4.add_collection(p2)
    fig.colorbar(aux4,ax=ax4,fraction=0.046, pad=0.04)

    ax3.set_title('Planet I parameter')
    ax4.set_title('Planet Q parameter')

    _plot_config_r(ax3)
    _plot_config_r(ax4)

    _plt.tight_layout()

    circle5 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    circle6 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax3.add_artist(circle5)
    ax4.add_artist(circle6)
	



    N 	= _np.sum(bodies[1].grid.shadow[t,:]!=0)
    faces00	= bodies[1].grid.faces[bodies[1].grid.shadow[t,:]!=0,0,0]
    faces10	= bodies[1].grid.faces[bodies[1].grid.shadow[t,:]!=0,1,0]
    faces01	= bodies[1].grid.faces[bodies[1].grid.shadow[t,:]!=0,0,1]
    I1	= bodies[1].grid.I[wvl,t,bodies[1].grid.shadow[t,:]!=0]#/I1ref
    Q1	= bodies[1].grid.Q[wvl,t,bodies[1].grid.shadow[t,:]!=0]#/Q1ref
    
    patches = []
    for i in range(N):
        square = Rectangle( (faces00[i]*2,faces10[i]*2),-2*faces00[i]+2*faces01[i],-2*faces00[i]+2*faces01[i] )
        patches.append(square)

    p3 = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)
    p4 = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)

    p3.set_array(I1)
    p4.set_array(Q1)
    aux5 = ax5.add_collection(p3)
    fig.colorbar(aux5,ax=ax5,fraction=0.046, pad=0.04)
    aux6 = ax6.add_collection(p4)
    fig.colorbar(aux6,ax=ax6,fraction=0.046, pad=0.04)

    ax5.set_title('Moon I parameter')
    ax6.set_title('Moon Q parameter')

    _plot_config_r(ax5)
    _plot_config_r(ax6)

    circle7 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    circle8 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax5.add_artist(circle7)
    ax6.add_artist(circle8)

    _plt.tight_layout()


    if save:
        filename = 'IQ_' + '_'  + _time.strftime("%d-%m-%Y") + '_' + _time.strftime("%X")

        fig.savefig('Images/'+filename + '.eps')
        fig.savefig('Images/'+filename + '.png')



def geometry_d(body, t = 0, save = False, dots = False):

    print('\n    ⇒ Plotting geometry parameters of ' + body.name+' at t = '+str(body.ephemeris.time[t])+' seconds')
    from matplotlib.patches import Rectangle

    # create a figure with subplots
    fig = _plt.figure(figsize=(10,8))
    ax1 = _plt.subplot2grid((2,2), (0,0))
    ax2 = _plt.subplot2grid((2,2), (0,1))
    ax3 = _plt.subplot2grid((2,2), (1,0))
    ax4 = _plt.subplot2grid((2,2), (1,1))


    if _cfg.plot_faces == True:
    	faces = _plot_color('faces')
    else:
        faces = "none"
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)

    patches = []
    for i in range(body.grid.N_points):
        square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
        patches.append(square)

    p1 = _PatchCollection(patches,cmap=_matplotlib.cm.YlOrRd, alpha = 1, edgecolor = faces)
    p2 = _PatchCollection(patches,cmap=_matplotlib.cm.YlOrRd, alpha = 1, edgecolor = faces)
    p3 = _PatchCollection(patches,cmap=_matplotlib.cm.coolwarm, alpha = 1, edgecolor = faces)
    p4 = _PatchCollection(patches,cmap=_matplotlib.cm.coolwarm, alpha = 1, edgecolor = faces)


    p1.set_array(_np.degrees(body.grid.solar_zenith_angle[t,:]))
    p2.set_array(_np.degrees(body.grid.observer_zenith_angle))
    p3.set_array(_np.degrees(body.grid.beta[t,:]))
    p4.set_array(_np.degrees(body.grid.azimuth[t,:]))
    aux1 = ax1.add_collection(p1)
    fig.colorbar(aux1,ax=ax1,fraction=0.046, pad=0.04)
    aux2 = ax2.add_collection(p2)
    fig.colorbar(aux2,ax=ax2,fraction=0.046, pad=0.04)
    aux3 = ax3.add_collection(p3)
    fig.colorbar(aux3,ax=ax3,fraction=0.046, pad=0.04)
    aux4 = ax4.add_collection(p4)
    fig.colorbar(aux4,ax=ax4,fraction=0.046, pad=0.04)

    circle1 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    circle2 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    circle3 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    circle4 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax1.add_artist(circle1)
    ax2.add_artist(circle2)
    ax3.add_artist(circle3)
    ax4.add_artist(circle4)

    ax1.set_title('SZA')
    ax2.set_title('EMISSION')
    ax3.set_title('BETA')
    ax4.set_title('PHI')

    _plot_config_r(ax1)
    _plot_config_r(ax2)
    _plot_config_r(ax3)
    _plot_config_r(ax4)
    _plt.tight_layout()

    if save:

        filename = 'geometry-pixels_' + planet1.properties.fourier_scene + '_alpha-' + str(int(_np.round(_np.degrees(body.geometry.phase_angle[t])))) + '_'  + _time.strftime("%d-%m-%Y") + '_' + _time.strftime("%X")

        fig.savefig('Images/'+filename + '.eps')
        fig.savefig('Images/'+filename + '.png')



def radiance(body, wvl=0, phase = False, save = False):
    # create a figure with subplots
    fig = _plt.figure(figsize=(9,12))
    ax1 = _plt.subplot2grid((4,1), (0,0))
    ax2 = _plt.subplot2grid((4,1), (1,0),sharex=ax1)
    ax3 = _plt.subplot2grid((4,1), (2,0),sharex=ax1)
    ax4 = _plt.subplot2grid((4,1), (3,0),sharex=ax1)

    if phase==False:
        xlabelstr = 'Time [Earth days]'
    else:
        xlabelstr = 'Phase angle'

    ax1.grid()
    ax1.set_title(body.name + ' I parameter')
    ax1.set_xlabel(xlabelstr)
    if phase is False:
        ax1.plot(body.ephemeris.time/3600/24, body.radiance.I[wvl,:],'-b')
    else:
        ax1.plot(_np.rad2deg(body.geometry.phase_angle), body.radiance.I[wvl,:],'-b')

    ax2.grid()
    ax2.set_title(body.name + ' Q parameter')
    ax2.set_xlabel(xlabelstr)
    if phase is False:
        ax2.plot(body.ephemeris.time/3600/24, body.radiance.Q[wvl,:],'-g')
    else:
        ax2.plot(_np.rad2deg(body.geometry.phase_angle), body.radiance.Q[wvl,:],'-g')
        ax3.grid()
    ax3.set_title(body.name + ' U parameter')
    ax3.set_xlabel(xlabelstr)
    if phase is False:
        ax3.plot(body.ephemeris.time/3600/24, body.radiance.U[wvl,:],'-r')
    else:
        ax3.plot(_np.rad2deg(body.geometry.phase_angle), body.radiance.U[wvl,:],'-r')

    ax4.grid()
    ax4.set_title(body.name + ' V parameter')
    ax4.set_xlabel(xlabelstr)
    if phase is False:
        ax4.plot(body.ephemeris.time/3600/24, body.radiance.V[wvl,:],'-k')
    else:
        ax4.plot(_np.rad2deg(body.geometry.phase_angle), body.radiance.V[wvl,:],'-k')

    _plt.tight_layout()

    if save:
        filename = 'radiance_' + body.properties.fourier_scene + '_alpha-' + str(int(_np.round(_np.degrees(body.geometry.phase_angle[t])))) + '_'  + _time.strftime("%d-%m-%Y") + '_' + _time.strftime("%X")

        fig.savefig('Images/'+filename + '.eps')
        fig.savefig('Images/'+filename + '.png')

    _plt.show()


def detail_radiance(bodies,I,Q,U,V, wvl = 0, save = False):

        body1 = bodies[0]
        body2 = bodies[1]

        # create a figure with subplots
        ax3 = _plt.subplot2grid((2,2), (0,0))
        ax4 = _plt.subplot2grid((2,2), (0,1),sharey  = ax3)
        ax5 = _plt.subplot2grid((2,2), (1,0),colspan = 2  )

        ax3.grid()
        ax3.set_title(body1.name + ' reflected light')
        ax3.set_xlabel('Time [Earth days]')
        ax3.plot(body1.ephemeris.time/3600/24, body1.radiance.I_ref[wvl,:],'-')
        #ax3.plot(body2.ephemeris.time, body2.radiance.I)
        ax3.plot(body1.ephemeris.time/3600/24, body1.radiance.Q_ref[wvl,:],'-')
        #ax3.plot(body2.ephemeris.time, body2.radiance.Q)
        ax3.plot(body1.ephemeris.time/3600/24, body1.radiance.U_ref[wvl,:],'-')
        #ax3.plot(body2.ephemeris.time, body2.radiance.U)

        ax4.grid()
        ax4.set_title(body2.name + ' reflected light')
        ax4.set_xlabel('Time [Earth days]')
        _plt.setp(ax4.get_yticklabels(), visible=False)
        #ax4.plot(body1.ephemeris.time, body1.radiance.I_ref)
        ax4.plot(body2.ephemeris.time/3600/24, body2.radiance.I_ref[wvl,:],'-')
        #ax4.plot(body1.ephemeris.time, body1.radiance.Q_ref)
        ax4.plot(body2.ephemeris.time/3600/24, body2.radiance.Q_ref[wvl,:],'-')
        #ax4.plot(body1.ephemeris.time, body1.radiance.U_ref)
        ax4.plot(body2.ephemeris.time/3600/24, body2.radiance.U_ref[wvl,:],'-')

        ax5.grid()
        ax5.set_title('Total reflected light')
        ax5.set_xlabel('Time [Earth days]')
        ax5.plot(body2.ephemeris.time/3600/24, I[wvl,:],'-')
        ax5.plot(body2.ephemeris.time/3600/24, Q[wvl,:],'-')
        ax5.plot(body2.ephemeris.time/3600/24, U[wvl,:],'-')
        ax5.legend(['I','Q','U'])
        _plt.tight_layout()

        if save:

                filename = 'detail-radiance_' + body1.properties.fourier_scene +'-'+ body2.properties.fourier_scene + '_alpha-' + str(int(_np.round(_np.degrees(body.geometry.phase_angle[t])))) + '_'  + _time.strftime("%d-%m-%Y") + '_' + _time.strftime("%X")

                fig.savefig('Images/'+filename + '.eps')
                fig.savefig('Images/'+filename + '.png')

        _plt.show()

def _move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = _matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        #f.canvas.manager.window.move(x, y)
        None
    _plt.show()

def _plot_color(string):

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

    return color[_np.where(color==string)[0][0],1+_cfg.plot_color]

def _initialize_indicators(ax, m=2):
    circles = [0,0,0,0,0]
    # Phase
    circles[0] = _plt.Circle((0.54*m, -0.54*m), 0.02*m, color = _plot_color('phase_c')    ,ec='k'   , fill=False, zorder=1)
    # Umbra
    circles[1] = _plt.Circle((0.54*m, -0.44*m), 0.01*m, color = _plot_color('umbra_c')    ,ec='none', fill=False, zorder=2)
    # Antumbra
    circles[2] = _plt.Circle((0.54*m, -0.44*m), 0.01*m, color = _plot_color('antumbra_c') ,ec='none', fill=False, zorder=2)
    # Penumbra
    circles[3] = _plt.Circle((0.54*m, -0.44*m), 0.02*m, color = _plot_color('penumbra_c') ,ec='k'   , fill=False, zorder=1)
    # Transit
    circles[4] = _plt.Circle((0.54*m, -0.49*m), 0.02*m, color = _plot_color('transit_c')  ,ec='k'   , fill=False, zorder=1)

    ax.add_artist(circles[0])
    ax.add_artist(circles[1])
    ax.add_artist(circles[2])
    ax.add_artist(circles[3])
    ax.add_artist(circles[4])

    return circles


def _update_indicators(ax,flags,t,circles):

    if flags.phase[t]:     circles[0].fill  = True
    else:                  circles[0].fill  = False
    if flags.umbra[t]:     circles[1].fill  = True
    else:                  circles[1].fill  = False
    if flags.antumbra[t]:  circles[2].fill  = True
    else:                  circles[2].fill  = False
    if flags.penumbra[t]:  circles[3].fill  = True
    else:                  circles[3].fill  = False
    if flags.transit[t]:   circles[4].fill  = True
    else:                  circles[4].fill  = False


def _load_flags(body,Type,time):
    flags = plot_flags(time)

    if Type == 'd':
        flags.PHASE   = body.flag.phase_d
        flags.TRANSIT = body.flag.transit_d
        flags.ECLIPSE = body.flag.eclipse_d

        if not flags.PHASE:
            print('    ... Caution: discrete phase data not available for '+body.type+' '+body.name+'!')
        if not flags.TRANSIT:
            print('    ... Caution: discrete transit data not available for '+body.type+' '+body.name+'!')
        if not flags.ECLIPSE:
            print('    ... Caution: discrete eclipse data not available for '+body.type+' '+body.name+'!')

    elif Type == 'c':
        flags.PHASE   = body.flag.phase_c
        flags.TRANSIT = body.flag.transit_c
        flags.ECLIPSE = body.flag.eclipse_c

        if not flags.PHASE:
            print('    ... Caution: continuous phase data not available for '+body.type+' '+body.name+'!')
        if not flags.TRANSIT:
            print('    ... Caution: continuous transit data not available for '+body.type+' '+body.name+'!')
        if not flags.ECLIPSE:
            print('    ... Caution: continuous eclipse data not available for '+body.type+' '+body.name+'!')

    aux = _np.zeros_like(body.ephemeris.time, dtype = bool)
    if flags.PHASE:                  flags.phase    = _np.ones(len(range(time[0],time[1])), dtype=bool)
    if len(body.flag.umbra)    != 0:
        aux = aux * 0
        for i in range(len(body.flag.umbra)):
            aux = aux + body.flag.umbra[i][0]
        flags.umbra  = (aux[time[0]:time[1]]*flags.ECLIPSE)

    if len(body.flag.antumbra) != 0:
        aux = aux * 0
        for i in range(len(body.flag.antumbra)):
            aux = aux + body.flag.antumbra[i][0]
        flags.antumbra  = (aux[time[0]:time[1]]*flags.ECLIPSE)

    if len(body.flag.penumbra) != 0:
        aux = aux * 0
        for i in range(len(body.flag.penumbra)):
            aux = aux + body.flag.penumbra[i][0]
        flags.penumbra  = (aux[time[0]:time[1]]*flags.ECLIPSE)

    if len(body.flag.transit)  != 0:
        aux = aux * 0
        for i in range(len(body.flag.transit)):
            aux = aux + body.flag.transit[i][0]
        flags.transit  = (aux[time[0]:time[1]]*flags.TRANSIT)

    return flags


def _plot_config(ax):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_axis_bgcolor((_plot_color('background')))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([-1.14, 1.14])
    ax.set_ylim([-1.14, 1.20])

def _plot_config_r(ax):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_axis_bgcolor((_plot_color('background')))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([-1.01, 1.01])
    ax.set_ylim([-1.01, 1.01])
    #ax.patch.set_visible(False)
    ax.axis('off')

def _initialize_text(ax, Type):

    if Type == 'r':
        text_r    = list(range(6))
        #text_r[3] = ax.text(-1.01, -1.02, '', fontsize = 12)
        #text_r[4] = ax.text(-1.01, -0.96, '', fontsize = 12)
        return text_r


    if Type == 'd':
        text_d    = list(range(6))
        text_d[0] = ax.text( 1.12,  1.12, '', fontsize = 12, horizontalalignment='right')
        text_d[1] = ax.text( 1.12,  1.03, '', fontsize = 12, horizontalalignment='right')
        text_d[2] = ax.text(-1.12,  1.12, '', fontsize = 12)
        text_d[3] = ax.text(-1.11, -1.12, '', fontsize = 12)
        text_d[4] = ax.text(-1.11, -1.04, '', fontsize = 12)
        text_d[5] = ax.text(-1.12, -0.96, '', fontsize = 12)
        return text_d

    if Type == 'c':
        text_c    = list(range(3))
#        text_c[0] = ax.text(-1.12,  1.12, '', fontsize = 12)
        text_c[0] = ax.text(-1.105,-1.12, '', fontsize = 12)
        text_c[1] = ax.text(-1.105,-1.04, '', fontsize = 12)
        text_c[2] = ax.text(-1.12, -0.96, '', fontsize = 12)
        return text_c

def _update_text(text, body, t, Type):

    #if Type == 'r':
        #text[3].set_text('Time   = %4d days, %02d hours (%d)'%
        #             (int(body.ephemeris.time[t]/60/60/24), (body.ephemeris.time[t]/24/60/60-
        #              int(body.ephemeris.time[t]/60/60/24))*24,t))
        #text[4].set_text('Phase = %04.2f$^{\circ}$' % (_np.degrees(body.geometry.phase_angle[t]),))

    if Type == 'd':
        text[0].set_text('Points: '+ str(body.grid.N_points))
        #text[1].set_text('$\^{A}$= '+ str(_np.round(_np.sum(body.grid.area)/(m.pi*0.5**2),3)))
        text[2].set_text('$A_{shadow}$= %4.2f [-]' % (_grid_area(body.grid,t)))
        text[3].set_text('t   = %4d days, %02d hours (%d)'%
                     (int(body.ephemeris.time[t]/60/60/24), (body.ephemeris.time[t]/24/60/60-
                      int(body.ephemeris.time[t]/60/60/24))*24,t))
        text[4].set_text('$Phase$  = %04.2f$^{\circ}$' % (_np.degrees(body.geometry.phase_angle[t]),))
        text[5].set_text('$\\varphi_S$ = %04.2f$^{\circ}$' % (_np.degrees(body.geometry.solar_azimuth_angle[t]),))
#
#        text[0].set_text('$Points: '+ str(body.grid.N_points)+'$')
#        text[1].set_text('$\^{A}= '+ str(_np.round(_np.sum(body.grid.area)/(m.pi*0.5**2),3))+'$')
#        text[2].set_text('$A_{shadow}= \, %4.2f \, [-]$' % (_grid_area(body.grid,t)))
#        text[3].set_text('$t \,\,\,\,\, = \, %4d \, days, \, %02d \, hours \, (%d)$'%
#                     (int(body.ephemeris.time[t]/60/60/24), (body.ephemeris.time[t]/24/60/60-
#                      int(body.ephemeris.time[t]/60/60/24))*24,t))
#        text[4].set_text('$\\alpha \,\,\, = \, %04.2f^{\circ}$' % (_np.degrees(body.geometry.phase_angle[t]),))
#        text[5].set_text('$\\varphi_S = \, %04.2f^{\circ}$' % (_np.degrees(body.geometry.solar_azimuth_angle[t]),))

    if Type == 'c':
        text[0].set_text('t   = %4d days, %02d hours (%d)'%
                     (int(body.ephemeris.time[t]/60/60/24), (body.ephemeris.time[t]/24/60/60-
                      int(body.ephemeris.time[t]/60/60/24))*24,t))
        text[1].set_text('$Phase$  = %04.2f$^{\circ}$' % (_np.degrees(body.geometry.phase_angle[t]),))
        text[2].set_text('$\\varphi_S$ = %04.2f$^{\circ}$' % (_np.degrees(body.geometry.solar_azimuth_angle[t]),))
#
#        text[0].set_text('$t \,\,\,\,\, = \, %4d \, days, \, %02d \, hours \, (%d)$'%
#                     (int(body.ephemeris.time[t]/60/60/24), (body.ephemeris.time[t]/24/60/60-
#                      int(body.ephemeris.time[t]/60/60/24))*24,t))
#        text[1].set_text('$Phase \,\,\, = \, %04.2f^{\circ}$' % (_np.degrees(body.geometry.phase_angle[t]),))
#        text[2].set_text('$\\varphi_S = \, %04.2f^{\circ}$' % (_np.degrees(body.geometry.solar_azimuth_angle[t]),))


class plot_flags():

    def __init__(self,time):

        self.PHASE   = False
        self.TRANSIT = False
        self.ECLIPSE = False
        self.phase      = _np.zeros(len(range(time[0],time[1])), dtype=bool)
        self.umbra      = _np.zeros(len(range(time[0],time[1])), dtype=bool)
        self.antumbra   = _np.zeros(len(range(time[0],time[1])), dtype=bool)
        self.penumbra   = _np.zeros(len(range(time[0],time[1])), dtype=bool)
        self.transit    = _np.zeros(len(range(time[0],time[1])), dtype=bool)


def grid(body):

    print('\n    ⇒ Plotting grid of ' + body.type + ' ' + body.name)

    fig , ax = _plt.subplots()
    _plt.plot(body.grid.nodes[:,0], body.grid.nodes[:,1], 'o', color = _plot_color('nodes'), markersize=3)
    for i in range(body.grid.N_points):
        _plt.plot(body.grid.faces[i,0,:], body.grid.faces[i,1,:], _plot_color('faces'), linewidth = 2)
    circle1 = _plt.Circle((0, 0), 0.5, color = _plot_color('circle'), fill=False, zorder=1)
    ax.add_artist(circle1)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([-0.55, 0.55])
    ax.set_ylim([-0.55, 0.6])
    ax.text(-0.4,  0.525,'Points: '+ str(body.grid.N_points), fontsize= 15)
    #ax.text(0.14,  0.525,'$\^{A}$: '+ str(_np.round(_np.sum(body.grid.area)/(m.pi*0.5**2),3)), fontsize= 15)


def shadow_d(body, t = 0, save = False, dots = False):

    if not hasattr(body.grid, 'shadow'):
    	body.grid.shadow = _np.ones([len(body.ephemeris.time), len(body.grid.nodes)])

    print('\n    ⇒ Plotting shadow of ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')

    flags = _load_flags(body,'d',[t,t+1])

    # Drawing
    fig = _plt.figure()
    #_move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)
    cell        = []#list(_np.zeros(body.grid.N_points))
    patch_cells = []#list(_np.zeros(body.grid.N_points))

    if _cfg.plot_faces == True:
    	ax.plot(body.grid.faces[:,0,:].T*2, body.grid.faces[:,1,:].T*2, color= _plot_color('faces'), linewidth = 0.5,zorder=2)
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)
    for i in range(body.grid.N_points):
        cell.append( [_Polygon(2*body.grid.faces[i,:,:].T)] )
        patch_cells.append( _PatchCollection(cell[i], alpha=1, color=str(body.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax.add_collection(patch_cells[i])


    # Body circle
    circle = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle)

    # Indicators
    circles = _initialize_indicators(ax)
    _update_indicators(ax,flags,0, circles)

    # Sun point
    P = [ 1.08*_np.cos(body.geometry.solar_azimuth_angle[t]) , 1.08*_np.sin(body.geometry.solar_azimuth_angle[t]+_np.pi)]
    circleP = _plt.Circle((P[0], P[1]), 0.03, color = _plot_color('sun'), fill=True, zorder=1)
    ax.add_artist(circleP)

    # Text
    text_d = _initialize_text(ax, 'd')
    _update_text(text_d, body, t, 'd')

    ax.set_title('Discretized '+body.type+' '+body.name)

    _plot_config(ax)
#    _plt.tight_layout()

    if save:

    	filename = 'shadow_' + body.properties.fourier_scene + '_alpha-' + str(int(_np.round(_np.degrees(body.geometry.phase_angle[t])))) + '_'  + _time.strftime("%d-%m-%Y") + '_' + _time.strftime("%X")

    	fig.savefig('Images/'+filename + '.eps')
    	fig.savefig('Images/'+filename + '.png')

    _plt.show()

def shadow_dd(body, t = [0,0,0], save = False, dots = False):

        if not hasattr(body.grid, 'shadow'):
                body.grid.shadow = _np.ones([len(body.ephemeris.time), len(body.grid.nodes)])

        print('\n    ⇒ Plotting shadow of ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t[0]])+', t ='+str(body.ephemeris.time[t[1]])+', and t ='+str(body.ephemeris.time[t[2]])+' seconds')
        from matplotlib.patches import Rectangle

        # create a figure with subplots
        fig = _plt.figure(figsize=(11,4))
        ax1 = _plt.subplot2grid((1,3), (0,0))
        ax2 = _plt.subplot2grid((1,3), (0,1))
        ax3 = _plt.subplot2grid((1,3), (0,2))

        if _cfg.plot_faces == True:
    	        faces = _plot_color('faces')
        else:
                faces = "none"
        if dots:
                ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)

        patches = []
        for i in range(body.grid.N_points):
                square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
                patches.append(square)

        p1 = _PatchCollection(patches, alpha = 1, cmap = 'Greys_r', edgecolor = faces)
        p2 = _PatchCollection(patches, alpha = 1, cmap = 'Greys_r', edgecolor = faces)
        p3 = _PatchCollection(patches, alpha = 1, cmap = 'Greys_r', edgecolor = faces)

        p1.set_array(body.grid.shadow[t[0],:])
        p2.set_array(body.grid.shadow[t[1],:])
        p3.set_array(body.grid.shadow[t[2],:])
        aux1 = ax1.add_collection(p1)
        aux2 = ax2.add_collection(p2)
        aux3 = ax3.add_collection(p3)

        circle1 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)        
        circle2 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
        circle3 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
        ax1.add_artist(circle1)
        ax2.add_artist(circle2)
        ax3.add_artist(circle3)

        ax1.set_title('t = %1d days %1d hours'%
                     (int(body.ephemeris.time[t[0]]/60/60/24), (body.ephemeris.time[t[0]]/24/60/60-
                      int(body.ephemeris.time[t[0]]/60/60/24))*24))
        ax2.set_title('t = %1d days %1d hours'%
                     (int(body.ephemeris.time[t[1]]/60/60/24), (body.ephemeris.time[t[1]]/24/60/60-
                      int(body.ephemeris.time[t[1]]/60/60/24))*24))
        ax3.set_title('t = %1d days %1d hours'%
                     (int(body.ephemeris.time[t[2]]/60/60/24), (body.ephemeris.time[t[2]]/24/60/60-
                      int(body.ephemeris.time[t[2]]/60/60/24))*24))

        _plot_config_r(ax1)
        _plot_config_r(ax2)
        _plot_config_r(ax3)
        _plt.tight_layout()

        if save:

                filename = 'shadow3_' + body.properties.fourier_scene + '_alpha-' + str(int(_np.round(_np.degrees(body.geometry.phase_angle[t[0]])))) + '-'+str(int(_np.round(_np.degrees(body.geometry.phase_angle[t[1]])))) + '-' + str(int(_np.round(_np.degrees(body.geometry.phase_angle[t[2]])))) + '_'  + _time.strftime("%d-%m-%Y") + '_' + _time.strftime("%X")

                fig.savefig('Images/'+filename + '.eps')
                fig.savefig('Images/'+filename + '.png')

        _plt.show()

def I_d(body, t = 0, wvl=0, dots = False):

    print('\n    ⇒ Plotting I parameter ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    from matplotlib.patches import Rectangle

    flags = _load_flags(body,'d',[t,t+1])

    # Drawing
    #fig = _plt.figure(figsize=(12,10))
    fig = _plt.figure(figsize=(10,8))
    _move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)

    if _cfg.plot_faces == True:
    	faces = _plot_color('faces')
    else:
        faces = "none"
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)

    patches = []
    for i in range(body.grid.N_points):
        square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
        patches.append(square)

    p = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces)
    p.set_array(body.grid.I[wvl, t,:])
    ax.add_collection(p)
    _plt.colorbar(p)

    # Body circle
    circle = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle)

    ax.set_title('I parameter')

    _plot_config_r(ax)
    _plt.tight_layout()
    #fig.patch.set_visible(False)


def Q_d(body, t = 0, wvl=0, dots = False):

    print('\n    ⇒ Plotting Q parameter ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    from matplotlib.patches import Rectangle

    flags = _load_flags(body,'d',[t,t+1])

    # Drawing
    #fig = _plt.figure(figsize=(12,10))
    fig = _plt.figure(figsize=(10,8))
    _move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)

    if _cfg.plot_faces == True:
        faces = _plot_color('faces')
    else:
        faces = "none"
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)

    patches = []
    for i in range(body.grid.N_points):
        square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
        patches.append(square)

    p = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces)
    p.set_array(body.grid.Q[wvl,t,:])
    ax.add_collection(p)
    _plt.colorbar(p)

    # Body circle
    circle = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle)

    ax.set_title('Q parameter')

    _plot_config_r(ax)
    _plt.tight_layout()


def U_d(body, t = 0, wvl=0, dots = False):

    print('\n    ⇒ Plotting U parameter ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    from matplotlib.patches import Rectangle

    flags = _load_flags(body,'d',[t,t+1])

    # Drawing
    #fig = _plt.figure(figsize=(12,10))
    fig = _plt.figure(figsize=(10,8))
    _move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)

    if _cfg.plot_faces == True:
        faces = _plot_color('faces')
    else:
        faces = "none"
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)

    patches = []
    for i in range(body.grid.N_points):
        square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
        patches.append(square)

    p = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces)
    p.set_array(body.grid.U[wvl,t,:])
    ax.add_collection(p)
    _plt.colorbar(p)

    # Body circle
    circle = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle)

    ax.set_title('U parameter')

    _plot_config_r(ax)
    _plt.tight_layout()


def V_d(body, t = 0, wvl=0, dots = False):

    print('\n    ⇒ Plotting V parameter ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    from matplotlib.patches import Rectangle

    flags = _load_flags(body,'d',[t,t+1])

    # Drawing
    #fig = _plt.figure(figsize=(12,10))
    fig = _plt.figure(figsize=(10,8))
    _move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)

    if _cfg.plot_faces == True:
        faces = _plot_color('faces')
    else:
        faces = "none"
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)

    patches = []
    for i in range(body.grid.N_points):
        square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
        patches.append(square)

    p = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces)
    p.set_array(body.grid.V[wvl,t,:])
    ax.add_collection(p)
    _plt.colorbar(p)

    # Body circle
    circle = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle)

    ax.set_title('V parameter')

    _plot_config_r(ax)
    _plt.tight_layout()


def radiance_d(body, t = 0, wvl=0, save = False, dots = False):

    print('\n    ⇒ Plotting stokes parameters of ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    from matplotlib.patches import Rectangle

    # create a figure with subplots
    fig = _plt.figure(figsize=(10,8))
    ax1 = _plt.subplot2grid((2,2), (0,0))
    ax2 = _plt.subplot2grid((2,2), (0,1))
    ax3 = _plt.subplot2grid((2,2), (1,0))
    ax4 = _plt.subplot2grid((2,2), (1,1))

    circle1 = _plt.Circle((0, 0), 1, facecolor = _plot_color('circle1'), edgecolor = _plot_color('circle'), linewidth = 2.5)
    circle2 = _plt.Circle((0, 0), 1, edgecolor = _plot_color('circle'), facecolor = _plot_color('circle1'), fill=True , linewidth = 2.5)
    circle3 = _plt.Circle((0, 0), 1, edgecolor = _plot_color('circle'), facecolor = _plot_color('circle1'), fill=True , linewidth = 2.5)
    circle4 = _plt.Circle((0, 0), 1, edgecolor = _plot_color('circle'), facecolor = _plot_color('circle1'), fill=True , linewidth = 2.5)
    ax1.add_artist(circle1)
    ax2.add_artist(circle2)
    ax3.add_artist(circle3)
    ax4.add_artist(circle4)

    if _cfg.plot_faces == True:
        faces = _plot_color('faces')
    else:
        faces = "none"
        if dots:
            ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)

    N 	= _np.sum(body.grid.shadow[t,:]!=0)
    faces00	= body.grid.faces[body.grid.shadow[t,:]!=0,0,0]
    faces10	= body.grid.faces[body.grid.shadow[t,:]!=0,1,0]
    faces01	= body.grid.faces[body.grid.shadow[t,:]!=0,0,1]
    I	= body.grid.I[wvl,t,body.grid.shadow[t,:]!=0]
    Q	= body.grid.Q[wvl,t,body.grid.shadow[t,:]!=0]
    U	= body.grid.U[wvl,t,body.grid.shadow[t,:]!=0]
    V	= body.grid.V[wvl,t,body.grid.shadow[t,:]!=0]

    patches = []
    for i in range(N):
        square = Rectangle( (faces00[i]*2,faces10[i]*2),-2*faces00[i]+2*faces01[i],-2*faces00[i]+2*faces01[i] )
        patches.append(square)

    p1 = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)
    p2 = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)
    p3 = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)
    p4 = _PatchCollection(patches,cmap=_matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)

    Q_aux = -Q/I
    U_aux = U/I
    V_aux = V/I

    Q_aux[_np.isnan(Q_aux)] = 0
    U_aux[_np.isnan(U_aux)] = 0
    V_aux[_np.isnan(V_aux)] = 0

    p1.set_array(I)
    p2.set_array(Q_aux)
    p3.set_array(U_aux)
    p4.set_array(V_aux)
    aux1 = ax1.add_collection(p1)
    fig.colorbar(aux1,ax=ax1,fraction=0.046, pad=0.04)
    aux2 = ax2.add_collection(p2)
    fig.colorbar(aux2,ax=ax2,fraction=0.046, pad=0.04)
    aux3 = ax3.add_collection(p3)
    fig.colorbar(aux3,ax=ax3,fraction=0.046, pad=0.04)
    aux4 = ax4.add_collection(p4)
    fig.colorbar(aux4,ax=ax4,fraction=0.046, pad=0.04)

    ax1.set_title('I parameter')
    ax2.set_title('Q parameter')
    ax3.set_title('U parameter')
    ax4.set_title('V parameter')

    _plot_config_r(ax1)
    _plot_config_r(ax2)
    _plot_config_r(ax3)
    _plot_config_r(ax4)
    _plt.tight_layout()

    circle5 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    circle6 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    circle7 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    circle8 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax1.add_artist(circle5)
    ax2.add_artist(circle6)
    ax3.add_artist(circle7)
    ax4.add_artist(circle8)

    if save:
        filename = 'radiance-pixels_' + body.properties.fourier_scene + '_alpha-' + str(int(_np.round(_np.degrees(body.geometry.phase_angle[t])))) + '_'  + _time.strftime("%d-%m-%Y") + '_' + _time.strftime("%X")
        fig.savefig('Images/'+filename + '.eps')
        fig.savefig('Images/'+filename + '.png')


def anim_shadow_1(body, dots = False, time  = 'all', info = [15, False]):

#==============================================================================
#   October 2016, Javier B.M., TU Delft
#------------------------------------------------------------------------------
#   Predefined
#
#   Inputs:
#       - body1: body object for plotting phase
#       - body2: second body object for plotting phase
#       - info = (0) x-label
#                (1) y-label
#                (2) interval
#                (3) store? True or False
#   Outputs:
#       - Plot
#==============================================================================
    global patch_cells
    global t_vector

    if time == 'all': time = [0,len(body.ephemeris.time)]

    print('\n    ⇒ Plotting shadow animation of ' + body.type + ' ' + body.name+' from t ='+str(body.ephemeris.time[time[0]])+' seconds to t = '+str(body.ephemeris.time[time[1]-1])+' seconds')

    flags = _load_flags(body,'d', time)

    # Drawing t=0
    fig = _plt.figure(figsize=(13,6))
    _move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)
    cell        = []#list(_np.zeros(body.grid.N_points))
    patch_cells = []#list(_np.zeros(body.grid.N_points))

    t_i = 0
    t = time[t_i]
    t_vector = _np.arange(time[0],time[1])

    ax.plot(body.grid.faces[:,0,:].T*2, body.grid.faces[:,1,:].T*2, _plot_color('faces'), linewidth = 0.5,zorder=2)
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)
    for i in range(body.grid.N_points):
        cell.append( [_Polygon(2*body.grid.faces[i,:,:].T)] )
        patch_cells.append( _PatchCollection(cell[i], alpha=1, color=str(body.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax.add_collection(patch_cells[i])

    # Body circle
    circle = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle)

    # Indicators
    circles = _initialize_indicators(ax)
    _update_indicators(ax,flags,t, circles)

    # Sun point
    P = [ 1.08*_np.cos(body.geometry.solar_azimuth_angle[time[0]:time[1]]+_np.pi) , 1.08*_np.sin(body.geometry.solar_azimuth_angle[time[0]:time[1]]+_np.pi)]

    circleP = _plt.Circle((P[0][t], P[1][t]), 0.03, color = _plot_color('sun'), fill=True, zorder=1)
    ax.add_artist(circleP)

    # Text
    text_d = _initialize_text(ax, 'd')
    _update_text(text_d, body, t, 'd')

    ax.set_title('Discretized '+body.type+' '+body.name)

    _plot_config(ax)

    def animate(t_i):
        global t_vector
        global patch_cells

        t_i %= len(_np.arange(time[0],time[1]))
        t    = t_vector[t_i]

        circleP.center = (P[0][t_i], P[1][t_i])
        _update_indicators(ax,flags,t_i, circles)
        _update_text(text_d, body, t, 'd')

        for i in range(body.grid.N_points):
            patch_cells[i].set_color(str(body.grid.shadow[t,i])) # set new color colors

    anim = _animation.FuncAnimation(fig, animate, interval=info[0], blit=False)

    if info[1] is True:
        _plt.rcParams['_animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
        mywriter = _animation.FFMpegWriter()
        anim.save('phase_anim.mp4', writer=mywriter)

    _plt.show()
    return anim

def anim_shadow_2(body1, body2, dots = False, time  = 'all', info = [ 15, False]):

#==============================================================================
#   October 2016, Javier B.M., TU Delft
#------------------------------------------------------------------------------
#   Predefined
#
#   Inputs:
#       - body1: body object for plotting phase
#       - body2: second body object for plotting phase
#       - info = (0) x-label
#                (1) y-label
#                (2) interval
#                (3) store? True or False
#   Outputs:
#       - Plot
#==============================================================================
    global t_vector
    global patch_cells1
    global patch_cells2

    if time == 'all': time = [0,len(body1.ephemeris.time)]

    print('\n    ⇒ Plotting shadow animation of ' + body1.type + ' ' + body1.name+' and '+ body2.type + ' ' + body2.name +' from t ='+str(body1.ephemeris.time[time[0]])+' seconds to t = '+str(body2.ephemeris.time[time[1]-1])+' seconds')

    flags1 = _load_flags(body1,'d', time)
    flags2 = _load_flags(body2,'d', time)

    # Drawing t=0
    fig = _plt.figure(figsize=(13,6))
    _move_figure(fig, 155, 110)
    ax1  = fig.add_subplot(1, 2, 1)
    ax2  = fig.add_subplot(1, 2, 2)
    cell1        = []#list(_np.zeros(body1.grid.N_points))
    patch_cells1 = []#list(_np.zeros(body1.grid.N_points))
    cell2        = []#list(_np.zeros(body1.grid.N_points))
    patch_cells2 = []#list(_np.zeros(body1.grid.N_points))

    t_i = 0
    t = time[t_i]
    t_vector = _np.arange(time[0],time[1])

    ax1.plot(body1.grid.faces[:,0,:].T*2, body1.grid.faces[:,1,:].T*2, _plot_color('faces'), linewidth = 0.5,zorder=2)
    if dots:
        ax1.plot(body1.grid.nodes[:,0]*2, body1.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)
    for i in range(body1.grid.N_points):
        cell1.append( [_Polygon(2*body1.grid.faces[i,:,:].T)] )
        patch_cells1.append( _PatchCollection(cell1[i], alpha=1, color=str(body1.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax1.add_collection(patch_cells1[i])
    ax2.plot(body2.grid.faces[:,0,:].T*2, body2.grid.faces[:,1,:].T*2, _plot_color('faces'), linewidth = 0.5,zorder=2)
    if dots:
        ax2.plot(body2.grid.nodes[:,0]*2, body2.grid.nodes[:,1]*2, 'o', color = _plot_color('nodes'), markersize=3)
    for i in range(body2.grid.N_points):
        cell2.append( [_Polygon(2*body2.grid.faces[i,:,:].T)] )
        patch_cells2.append( _PatchCollection(cell2[i], alpha=1, color=str(body2.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax2.add_collection(patch_cells2[i])

    # Body circle
    circle1 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=10,linewidth = 2.5)
    ax1.add_artist(circle1)
    circle2 = _plt.Circle((0, 0), 1, color = _plot_color('circle'), fill=False, zorder=10,linewidth = 2.5)
    ax2.add_artist(circle2)

    # Indicators
    circles1 = _initialize_indicators(ax1)
    _update_indicators(ax1,flags1,t_i, circles1)
    circles2 = _initialize_indicators(ax2)
    _update_indicators(ax2,flags2,t_i, circles2)

    # Sun point
    P1 = [ 1.08*_np.cos(body1.geometry.solar_azimuth_angle[time[0]:time[1]]+_np.pi) , 1.08*_np.sin(body1.geometry.solar_azimuth_angle[time[0]:time[1]]+_np.pi)]
    circleP1 = _plt.Circle((P1[0][t_i], P1[1][t_i]), 0.03, color = _plot_color('sun'), fill=True, zorder=1)
    ax1.add_artist(circleP1)
    P2 = [ 1.08*_np.cos(body2.geometry.solar_azimuth_angle[time[0]:time[1]]+_np.pi) , 1.08*_np.sin(body2.geometry.solar_azimuth_angle[time[0]:time[1]]+_np.pi)]
    circleP2 = _plt.Circle((P2[0][t_i], P2[1][t_i]), 0.03, color = _plot_color('sun'), fill=True, zorder=1)
    ax2.add_artist(circleP2)

    # Text
    text_d1 = _initialize_text(ax1, 'd')
    _update_text(text_d1, body1, t, 'd')
    text_d2 = _initialize_text(ax2, 'd')
    _update_text(text_d2, body2, t, 'd')

    ax1.set_title('Discretized '+body1.type+' '+body1.name)
    ax2.set_title('Discretized '+body2.type+' '+body2.name)

    _plot_config(ax1)
    _plot_config(ax2)


    def animate(t_i):
        global t_vector
        global patch_cells1
        global patch_cells2

        t_i %= len(_np.arange(time[0],time[1]))
        t    = t_vector[t_i]

        circleP1.center = (P1[0][t_i], P1[1][t_i])
        _update_indicators(ax1,flags1,t_i, circles1)
        _update_text(text_d1, body1, t, 'd')
        circleP2.center = (P2[0][t_i], P2[1][t_i])
        _update_indicators(ax2,flags2,t_i, circles2)
        _update_text(text_d2, body2, t, 'd')

        for i in range(body1.grid.N_points):
            patch_cells1[i].set_color(str(body1.grid.shadow[t,i])) # set new color colors
        for i in range(body2.grid.N_points):
            patch_cells2[i].set_color(str(body2.grid.shadow[t,i])) # set new color colors

    anim = _animation.FuncAnimation(fig, animate, interval=info[0], blit=False)

    if info[1] is True:
        _plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
        mywriter = _animation.FFMpegWriter()
        anim.save('phase_anim.mp4', writer=mywriter)

    return anim





def ecl(t, body1, body2, star):

    import sys
    sys.path.append('/home/javier/anaconda3/lib/python3.5/site-packages')
    import matplotlib.pyplot as _plt
    #from shapely.geometry import Polygon as Poly
    from matplotlib.collections import _PatchCollection
    import math as m
    #from descartes import PolygonPatch
    import numpy as _np

    R1 = _np.zeros(3)
    R2 = _np.zeros(3)

    R1[0] =  body1.ephemeris.position3D_s_ob[0][t]
    R1[1] =  body1.ephemeris.position3D_s_ob[1][t]
    R1[2] =  body1.ephemeris.position3D_s_ob[2][t]

    R2[0] =  body2.ephemeris.position3D_s_ob[0][t]
    R2[1] =  body2.ephemeris.position3D_s_ob[1][t]
    R2[2] =  body2.ephemeris.position3D_s_ob[2][t]

    R1abs = body1.ephemeris.r_s[t]
    R2abs = body2.ephemeris.r_s[t]

    ds = R1.dot(R2)/R2abs
    d2 = ds -body2.ephemeris.r_s[t]
    d1 = (R1abs**2-(R2abs+d2)**2)**0.5

    rs = star.properties.R
    r2 = body2.properties.R
    r1 = body1.properties.R

    d2s = R2abs

    psi = m.asin((rs-r2)/d2s)
    print('PSI', _np.degrees(psi))
    omega = m.asin((r2+rs)/d2s)
    print('OMEGA', _np.degrees(omega))

    O2A6 = (r1+r2)/m.sin(omega)
    print('O2A6: ',O2A6)
    print('d2: ',d2)
    print('d1: ',d1)
    rho = m.atan(d1/(d2+O2A6))
    print('RHO', _np.degrees(rho))
    if rho>omega:
        print('rho larger than omega')
    else:
        print('omega larger than rho')


    fig = _plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    x_aux1 = d2-r2/m.sin(psi)
    x_aux2 = d2+r2/m.sin(omega)
    r = m.tan(psi)/m.tan(omega)
    x_inter = (x_aux2+r*x_aux1)/(r+1)

    penumbra = []
    #penumbra1 = Poly([(d2-r2/m.sin(psi),0),(x_inter,m.tan(psi)*(x_inter-d2+r2/m.sin(psi))),(d2+r2/m.sin(omega)-2*ds*m.cos(omega),2*ds*m.sin(omega)),(d2-r2/m.sin(psi)-2*ds*m.cos(psi),2*ds*m.sin(psi))])
    #penumbra2 = Poly([(d2-r2/m.sin(psi),0),(x_inter,-m.tan(psi)*(x_inter-d2+r2/m.sin(psi))),(d2+r2/m.sin(omega)-2*ds*m.cos(omega),-2*ds*m.sin(omega)),(d2-r2/m.sin(psi)-2*ds*m.cos(psi),-2*ds*m.sin(psi))])
    #penumbra.append(PolygonPatch(penumbra1))
    #penumbra.append(PolygonPatch(penumbra2))
    #p = _PatchCollection(penumbra, alpha=1, color='0.5', edgecolor=None, zorder = 1)
    #ax.add_collection(p)

    antumbra = []
    #antumbra1 = Poly([(d2-r2/m.sin(psi),0),(d2-r2/m.sin(psi)-2*ds*m.cos(psi),2*ds*m.sin(psi)),(d2-r2/m.sin(psi)-2*ds*m.cos(psi),-2*ds*m.sin(psi))])
    #antumbra.append(PolygonPatch(antumbra1))
    #p = _PatchCollection(antumbra, alpha=1, color='0.75', edgecolor=None, zorder = 1)
    #ax.add_collection(p)

    umbra = []
    #umbra1 = Poly([(d2-r2/m.sin(psi),0) , ( d2-r2*m.sin(psi),r2*m.cos(psi) ) , (d2,0) , ( d2-r2*m.sin(psi),-r2*m.cos(psi) )])
    #umbra.append(PolygonPatch(umbra1))
    #p = _PatchCollection(umbra, alpha=1, color='0.1', edgecolor=None, zorder = 1)
    #ax.add_collection(p)

    circ=_plt.Circle((ds,0), radius=rs, color='b', ec= 'k', fill=True)
    ax.add_patch(circ)
    circ=_plt.Circle((d2,0), radius=r2, color='r', ec = 'k', fill=True)
    ax.add_patch(circ)
    circ=_plt.Circle((0,d1), radius=r1, color='g', ec = 'k', fill=True)
    ax.add_patch(circ)

    ax.arrow(d2-r2/m.sin(psi), 0, 2*ds*m.cos(psi),2*ds*m.sin(psi), head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.arrow(d2-r2/m.sin(psi), 0, 2*ds*m.cos(psi),-2*ds*m.sin(psi), head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.arrow(d2-r2/m.sin(psi), 0, -2*ds*m.cos(psi),2*ds*m.sin(psi), head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.arrow(d2-r2/m.sin(psi), 0, -2*ds*m.cos(psi),-2*ds*m.sin(psi), head_width=0.05, head_length=0.1, fc='k', ec='k')

    ax.arrow(d2+r2/m.sin(omega), 0, 2*ds*m.cos(omega),2*ds*m.sin(omega), head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.arrow(d2+r2/m.sin(omega), 0, 2*ds*m.cos(omega),-2*ds*m.sin(omega), head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.arrow(d2+r2/m.sin(omega), 0, -2*ds*m.cos(omega),-2*ds*m.sin(omega), head_width=0.05, head_length=0.1, fc='k', ec='k')
    ax.arrow(d2+r2/m.sin(omega), 0, -2*ds*m.cos(omega),2*ds*m.sin(omega), head_width=0.05, head_length=0.1, fc='k', ec='k')

    #patches = []
    #patches.append(Polygon(, True))
    #p = _PatchCollection(patches, color='red', alpha=0.4)
    #ax.add_collection(p)

    ax.set_xlim([-r1, d2])
    ax.set_ylim([-r2*20,  r2*20])
    #ax.set_aspect('equal', adjustable='box')
    ax.set_axis_bgcolor('0.9')

    _plt.show()







def xyorbit(position2D_1, position2D_2 = None,  info =
        ['x axis [m]', 'y axis [m]','2-D orbit', None, 'O']):
        # X axis label, Y axis label, Title, Legend, Origin name

#==============================================================================
#   October 2016, Javier B.M., TU Delft
#------------------------------------------------------------------------------
#   Predefined function for 2D orbit plots (up to two orbits)
#
#   Inputs:
#       - position2D_1: Vector position of first body (x,y)
#       - position2D_2: Vector position of second body (x,y)
#       - info = (0) x-label
#                (1) y-label
#                (2) title
#                (3) legend
#                (4) origin name
#   Outputs:
#       - Plot
#==============================================================================

    fig = _plt.figure()

    # The orbits of the two bodies are plotted, if available
    _plt.plot(position2D_1[0], position2D_1[1],'r', linewidth=1,
                                                 linestyle="-", label="Body 1")
    if position2D_2 is not None:
        _plt.plot(position2D_2[0], position2D_2[1],'b', linewidth=1,
                                                 linestyle="-", label="Body 2")

    # Origin point and name
    _plt.plot(0,0, 'ko')
    ax = fig.gca()
    ax.annotate(info[4], xy=(0, 0), xytext=(0,0),
                arrowprops = dict(facecolor='black', shrink=0.05),
                horizontalalignment='right', verticalalignment='top',
                color='black', size= 'large',)

    # Axes labels and title are established
    _plt.xlabel(info[0])
    _plt.ylabel(info[1])
    _plt.title(info[2])

    # A same scale is applied to both axes
    _plt.gca().set_aspect('equal', adjustable='box')
    _plt.grid()

    # Axes limits and legend depend on the number of curves plotted
    if position2D_2 is not None:
        aux1 = max([position2D_1.max(),position2D_2.max()])
        _plt.xlim([min([position2D_1[0].min(),position2D_2[0].min()]),
                           max([position2D_1[0].max(),position2D_2[0].max()])])
        _plt.ylim([min([position2D_1[1].min(),position2D_2[1].min()]),
                           max([position2D_1[1].max(),position2D_2[1].max()])])

        if info[3] is not None:
            _plt.legend([info[3][0],info[3][1]],loc='upper left', frameon=False)
        else:
            _plt.legend(loc='upper left', frameon=False)
    else:
        aux1 = position2D_1.max()
        _plt.xlim([position2D_1[0].min(),position2D_1[0].max()])
        _plt.ylim([position2D_1[1].min(),position2D_1[1].max()])

        if info[3] is not None:
            _plt.legend([info[3]],loc='upper left', frameon=False)

    # Reference frame arrows are added
    ax.arrow( 0, 0, aux1/3, 0, fc="r", ec="r",
                                     head_width=aux1/20, head_length=aux1/20 )
    ax.arrow( 0, 0, 0, aux1/3, fc="g", ec="g",
                                     head_width=aux1/20, head_length=aux1/20 )

    _plt.show()




def XYZorbit(position3D_1, position3D_2 = None, info =
        ['x [m]',' y [m]','z [m]','3-D orbit',None,'O','general']):
   # X axis label, Y axis label, Z axis label, Title, Legend, Origin name, view

#==============================================================================
#   October 2016, Javier B.M., TU Delft
#------------------------------------------------------------------------------
#   Predefined function for 3D orbit plots (up to two orbits)
#
#   Inputs:
#       - position2D_1: Vector position of first body (x,y)
#       - position2D_2: Vector position of second body (x,y)
#       - info = (0) x-label
#                (1) y-label
#                (2) z-label
#                (3) title
#                (4) legend
#                (5) origin name
#                (6) view: close or general
#   Outputs:
#       - Plot
#==============================================================================

    import statistics as stat
    import exopy_functions as fun

    fig = _plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # The orbits of the two bodies are plotted, if available
    ax.plot(position3D_1[0], position3D_1[1],
                                           position3D_1[2],'r', label='Body 1')
    if position3D_2 is not None:
        ax.plot(position3D_2[0], position3D_2[1], position3D_2[2],'b', label='Body 2')

    # Axes labels and title are established
    ax = fig.gca()
    ax.set_xlabel(info[0])
    ax.set_ylabel(info[1])
    ax.set_zlabel(info[2])
    _plt.title(info[3])

    fig.gca().set_aspect('equal', adjustable='box')
    _plt.grid()

    # Axes limits and legend depend on the number of curves plotted
    if position3D_2 is not None:
        if info[6] is 'close':

            aux1 = max([position3D_1[0].max(),position3D_2[0].max()])
            aux2 = min([position3D_1[0].min(),position3D_2[0].min()])

            aux3 = max([position3D_1[1].max(),position3D_2[1].max()])
            aux4 = min([position3D_1[1].min(),position3D_2[1].min()])

            aux5 = max([position3D_1[2].max(),position3D_2[2].max()])
            aux6 = min([position3D_1[2].min(),position3D_2[2].min()])

            aux  = max([aux1 - aux2, aux3 - aux4, aux5 - aux6])

            ax.set_xlim3d(stat.mean([aux1, aux2]) - aux/2,
                                               stat.mean([aux1, aux2]) + aux/2)
            ax.set_ylim3d(stat.mean([aux3, aux4]) - aux/2,
                                               stat.mean([aux3, aux4]) + aux/2)
            ax.set_zlim3d(stat.mean([aux5, aux6]) - aux/2,
                                               stat.mean([aux5, aux6]) + aux/2)

        if info[6] is 'general':

            aux1 = max([position3D_1[0].max(),position3D_2[0].max(), 0])
            aux2 = min([position3D_1[0].min(),position3D_2[0].min(), 0])

            aux3 = max([position3D_1[1].max(),position3D_2[1].max(), 0])
            aux4 = min([position3D_1[1].min(),position3D_2[1].min(), 0])

            aux5 = max([position3D_1[2].max(),position3D_2[2].max(), 0])
            aux6 = min([position3D_1[2].min(),position3D_2[2].min(), 0])

            aux  = max([aux1 - aux2, aux3 - aux4, aux5 - aux6])

            ax.set_xlim3d(stat.mean([aux1, aux2]) - aux/2,
                                               stat.mean([aux1, aux2]) + aux/2)
            ax.set_ylim3d(stat.mean([aux3, aux4]) - aux/2,
                                               stat.mean([aux3, aux4]) + aux/2)
            ax.set_zlim3d(stat.mean([aux5, aux6]) - aux/2,
                                               stat.mean([aux5, aux6]) + aux/2)

        if info[4] is not None:
            _plt.legend([info[4][0],info[4][1]],loc='center left',frameon=False)
        else:
            _plt.legend(loc='center left', frameon=False)
    else:
        if info[6] is 'close':
            aux1 = max([position3D_1[0].max()])
            aux2 = min([position3D_1[0].min()])

            aux3 = max([position3D_1[1].max()])
            aux4 = min([position3D_1[1].min()])

            aux5 = max([position3D_1[2].max()])
            aux6 = min([position3D_1[2].min()])

            aux  = max([aux1 - aux2, aux3 - aux4, aux5 - aux6])

            ax.set_xlim3d(stat.mean([aux1, aux2]) - aux/2,
                                               stat.mean([aux1, aux2]) + aux/2)
            ax.set_ylim3d(stat.mean([aux3, aux4]) - aux/2,
                                               stat.mean([aux3, aux4]) + aux/2)
            ax.set_zlim3d(stat.mean([aux5, aux6]) - aux/2,
                                               stat.mean([aux5, aux6]) + aux/2)

        if info[6] is 'general':
            aux1 = max([position3D_1[0].max(), 0])
            aux2 = min([position3D_1[0].min(), 0])

            aux3 = max([position3D_1[1].max(), 0])
            aux4 = min([position3D_1[1].min(), 0])

            aux5 = max([position3D_1[2].max(), 0])
            aux6 = min([position3D_1[2].min(), 0])

            aux  = max([aux1 - aux2, aux3 - aux4, aux5 - aux6])

            ax.set_xlim3d(stat.mean([aux1, aux2]) - aux/2,
                                               stat.mean([aux1, aux2]) + aux/2)
            ax.set_ylim3d(stat.mean([aux3, aux4]) - aux/2,
                                               stat.mean([aux3, aux4]) + aux/2)
            ax.set_zlim3d(stat.mean([aux5, aux6]) - aux/2,
                                               stat.mean([aux5, aux6]) + aux/2)

        if info[4] is not None:
            _plt.legend([info[4]],loc='center left', frameon=False)

    # Reference frame arrows are added
    Xvector = fun.Arrow3D([0,aux/6],[0,  0  ],[0,  0  ],
                          mutation_scale=20, lw=1, arrowstyle="-|>", color="r")
    Yvector = fun.Arrow3D([0,  0  ],[0,aux/6],[0,  0  ],
                          mutation_scale=20, lw=1, arrowstyle="-|>", color="g")
    Zvector = fun.Arrow3D([0,  0  ],[0,  0  ],[0,aux/6],
                          mutation_scale=20, lw=1, arrowstyle="-|>", color="y")
    ax.add_artist(Xvector)
    ax.add_artist(Yvector)
    ax.add_artist(Zvector)

    # Origin point and text
    ax.scatter([0],[0],[0], color="k",s=40)
    ax.text(aux/100, aux/100, aux/100, info[5], [0,1,0])

    _plt.show()






def anim_orbit(time, position3D_1, position3D_2 = None, info =
 ['x axis [m]','y axis [m]', 'z axis [m]','3-D orbit in cartesian coordinates',
  None,'O','close',False, 60]):
   # X label,Y label,Z label,Title,Legend,Origin name, view,store?,interval

#==============================================================================
#   October 2016, Javier B.M., TU Delft
#------------------------------------------------------------------------------
#   Predefined function for 3D orbit plots (up to two orbits)
#
#   Inputs:
#       - position2D_1: Vector position of first body (x,y)
#       - position2D_2: Vector position of second body (x,y)
#       - info = (0) x-label
#                (1) y-label
#                (2) z-label
#                (3) title
#                (4) legend
#                (5) origin name
#                (6) view: close or general
#                (7) store animation?: True or False
#                (8) interval
#   Outputs:
#       - Plot
#==============================================================================

    import matplotlib.pyplot as __plt
    import mpl_toolkits.mplot3d.axes3d as p3
    import matplotlib.animation as animation
    import statistics as stat
    import exopy_functions as fun

    # Update animation function
    def update_lines(num, dataLines, lines, pts) :
        for line, pt, data in zip(lines, pts, dataLines) :
            line.set_data(data[0:2,:num])
            line.set_3d_properties(data[2,:num])
            pt.set_data(data[0:2,num-1:num])
            pt.set_3d_properties(data[2,num-1:num])
            time_text.set_text('Time: %4.2f Solar days' % (time[num]/60/60/24))

        return lines

    # Attach 3D axis to the figure
    fig = __plt.figure()
    ax = p3.Axes3D(fig)

    # Axes labels and title are established
    ax = fig.gca()
    ax.set_xlabel(info[0])
    ax.set_ylabel(info[1])
    ax.set_zlabel(info[2])
    ax.set_title(info[3])

    fig.gca().set_aspect('equal', adjustable='box')
    __plt.grid()

    leg = [0, 0]

    # Axes limits and legend depend on the number of curves plotted
    if position3D_2 is not None:
        if info[6] is 'close':

            aux1 = max([position3D_1[0].max(),position3D_2[0].max()])
            aux2 = min([position3D_1[0].min(),position3D_2[0].min()])

            aux3 = max([position3D_1[1].max(),position3D_2[1].max()])
            aux4 = min([position3D_1[1].min(),position3D_2[1].min()])

            aux5 = max([position3D_1[2].max(),position3D_2[2].max()])
            aux6 = min([position3D_1[2].min(),position3D_2[2].min()])

            aux  = max([aux1 - aux2, aux3 - aux4, aux5 - aux6])

            ax.set_xlim3d(stat.mean([aux1, aux2]) - aux/2,
                                               stat.mean([aux1, aux2]) + aux/2)
            ax.set_ylim3d(stat.mean([aux3, aux4]) - aux/2,
                                               stat.mean([aux3, aux4]) + aux/2)
            ax.set_zlim3d(stat.mean([aux5, aux6]) - aux/2,
                                               stat.mean([aux5, aux6]) + aux/2)

        elif info[6] is 'general':

            aux1 = max([position3D_1[0].max(),position3D_2[0].max(), 0])
            aux2 = min([position3D_1[0].min(),position3D_2[0].min(), 0])

            aux3 = max([position3D_1[1].max(),position3D_2[1].max(), 0])
            aux4 = min([position3D_1[1].min(),position3D_2[1].min(), 0])

            aux5 = max([position3D_1[2].max(),position3D_2[2].max(), 0])
            aux6 = min([position3D_1[2].min(),position3D_2[2].min(), 0])

            aux  = max([aux1 - aux2, aux3 - aux4, aux5 - aux6])

            ax.set_xlim3d(stat.mean([aux1, aux2]) - aux/2,
                                               stat.mean([aux1, aux2]) + aux/2)
            ax.set_ylim3d(stat.mean([aux3, aux4]) - aux/2,
                                               stat.mean([aux3, aux4]) + aux/2)
            ax.set_zlim3d(stat.mean([aux5, aux6]) - aux/2,
                                               stat.mean([aux5, aux6]) + aux/2)

        if info[4] is None:
            leg[0] = 'Body 1'
            leg[1] = 'Body 2'
        else:
            leg[0] = info[4][0]
            leg[1] = info[4][1]
    else:
        if info[6] is 'close':
            print('a')
            aux1 = max([position3D_1[0].max()])
            aux2 = min([position3D_1[0].min()])

            aux3 = max([position3D_1[1].max()])
            aux4 = min([position3D_1[1].min()])

            aux5 = max([position3D_1[2].max()])
            aux6 = min([position3D_1[2].min()])

            aux  = max([aux1 - aux2, aux3 - aux4, aux5 - aux6])

            ax.set_xlim3d(stat.mean([aux1, aux2]) - aux/2,
                                               stat.mean([aux1, aux2]) + aux/2)
            ax.set_ylim3d(stat.mean([aux3, aux4]) - aux/2,
                                               stat.mean([aux3, aux4]) + aux/2)
            ax.set_zlim3d(stat.mean([aux5, aux6]) - aux/2,
                                               stat.mean([aux5, aux6]) + aux/2)

        elif info[6] is 'general':
            aux1 = max([position3D_1[0].max(), 0])
            aux2 = min([position3D_1[0].min(), 0])

            aux3 = max([position3D_1[1].max(), 0])
            aux4 = min([position3D_1[1].min(), 0])

            aux5 = max([position3D_1[2].max(), 0])
            aux6 = min([position3D_1[2].min(), 0])

            aux  = max([aux1 - aux2, aux3 - aux4, aux5 - aux6])

            ax.set_xlim3d(stat.mean([aux1, aux2]) - aux/2,
                                               stat.mean([aux1, aux2]) + aux/2)
            ax.set_ylim3d(stat.mean([aux3, aux4]) - aux/2,
                                               stat.mean([aux3, aux4]) + aux/2)
            ax.set_zlim3d(stat.mean([aux5, aux6]) - aux/2,
                                               stat.mean([aux5, aux6]) + aux/2)

        if info[4] is None:
            leg[0] = None
            leg[1] = None
        else:
            leg[0] = info[4][0]
            leg[1] = None


    # Reference frame arrows are added
    Xvector = fun.Arrow3D([0,aux/6],[0,  0  ],[0,  0  ],
                          mutation_scale=20, lw=1, arrowstyle="-|>", color="r")
    Yvector = fun.Arrow3D([0,  0  ],[0,aux/6],[0,  0  ],
                          mutation_scale=20, lw=1, arrowstyle="-|>", color="g")
    Zvector = fun.Arrow3D([0,  0  ],[0,  0  ],[0,aux/6],
                          mutation_scale=20, lw=1, arrowstyle="-|>", color="y")
    ax.add_artist(Xvector)
    ax.add_artist(Yvector)
    ax.add_artist(Zvector)

    # Origin point and text
    ax.scatter([0],[0],[0], color="k",s=40)
    ax.text(aux/100, aux/100, aux/100, info[5], [0,1,0])

    # Set up animation data
    n = len(position3D_1[0])

    if position3D_2 is not None:
        data = [_np.vstack((position3D_1[0], position3D_1[1], position3D_1[2])),
                _np.vstack((position3D_2[0], position3D_2[1], position3D_2[2])),
                _np.vstack((position3D_1[0], position3D_1[1], position3D_1[2])),
                _np.vstack((position3D_2[0], position3D_2[1], position3D_2[2]))]
        lines = [ax.plot(data[0][0,0:1], data[0][1,0:1], data[0][2,0:1],
                                                      'r' , label = leg[0])[0],
                 ax.plot(data[0][0,0:1], data[0][1,0:1], data[0][2,0:1],
                                                      'b' , label = leg[1])[0]]
        pts   = [ax.plot(data[0][0,0:1],data[0][1,0:1],data[0][2,0:1],'ro')[0],
                 ax.plot(data[0][0,0:1],data[0][1,0:1],data[0][2,0:1],'bo')[0]]

    else:

        data = [_np.vstack((position3D_1[0], position3D_1[1], position3D_1[2])),
                _np.vstack((0, 0, 0)),
                _np.vstack((position3D_1[0], position3D_1[1], position3D_1[2])),
                _np.vstack((0, 0, 0))]
        lines = [ax.plot(data[0][0,0:1], data[0][1,0:1], data[0][2,0:1],
                                                       'r', label = leg[0])[0],
                 ax.plot(data[0][0,0:1], data[0][1,0:1], data[0][2,0:1],
                                                       'b', label = leg[1])[0]]
        pts   = [ax.plot(data[0][0,0:1],data[0][1,0:1],data[0][2,0:1],'ro')[0],
                 ax.plot(data[0][0,0:1],data[0][1,0:1],data[0][2,0:1],'bo')[0]]

    time_text = ax.text2D(0.2, 0.2, "2D Text", transform=ax.transAxes)
    __plt.legend(loc="center left")

    # Creating the Animation object
    ani = animation.FuncAnimation( fig, update_lines,n,fargs=(data,lines,pts),
                                                  interval=info[8], blit=False)

    if info[7] is 'store_yes':
        __plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
        mywriter = animation.FFMpegWriter()
        ani.save('mymovie.mp4', writer=mywriter)

    return ani









def plot_planetmoon(time, planet, moon):

    fig = _plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")

    #draw sphere
    u, v = _np.mgrid[0:2*_np.pi:20j, 0:_np.pi:10j]
    x=_np.cos(u)*_np.sin(v)
    y=_np.sin(u)*_np.sin(v)
    z=_np.cos(v)

    pos = moon.ephemeris.position3D_s[:,time]-planet.ephemeris.position3D_s[:,time]
    x0 = pos[0]
    y0 = pos[1]
    z0 = pos[2]

#    planet.properties.R = 1
#    moon.properties.R = 1

    ax.plot_wireframe(planet.properties.R*x, planet.properties.R*y, planet.properties.R*z, color="r")
    ax.plot_wireframe(x0+moon.properties.R*x,y0+moon.properties.R*y,z0+moon.properties.R*z, color="b")

    m = max(x0+moon.properties.R*x.max(),y0+moon.properties.R*y.max(),z0+moon.properties.R*z.max())

    ax.set_xlim3d(0, m)
    ax.set_ylim3d(0, m)
    ax.set_zlim3d(0, m)



    return



