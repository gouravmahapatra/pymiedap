# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 22:01:58 2016

@author: javier
"""
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon    
from matplotlib.collections import PatchCollection
from matplotlib import animation
#from shapely.geometry import Polygon as Poly
import numpy as np
import math as m
from module_functions import grid_area
import exopy_cfg as _cfg
import time


def plot_geometry_d(body, t = 0, save = False, dots = False):

    	print('\n    ⇒ Plotting geometry parameters of ' + body.name+' at t = '+str(body.ephemeris.time[t])+' seconds')
    	from matplotlib.patches import Rectangle

        # create a figure with subplots
    	fig = plt.figure(figsize=(10,8))
        ax1 = plt.subplot2grid((2,2), (0,0))
        ax2 = plt.subplot2grid((2,2), (0,1))
	ax3 = plt.subplot2grid((2,2), (1,0))
	ax4 = plt.subplot2grid((2,2), (1,1))
    
    	if _cfg.plot_faces == True:
    		faces = plot_color('faces')
    	else:
		faces = "none"  
    	if dots:
        	ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)

	

    	patches = []
    	for i in range(body.grid.N_points):
		square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
		patches.append(square)

   	p1 = PatchCollection(patches,cmap=matplotlib.cm.YlOrRd, alpha = 1, edgecolor = faces)
   	p2 = PatchCollection(patches,cmap=matplotlib.cm.YlOrRd, alpha = 1, edgecolor = faces)
   	p3 = PatchCollection(patches,cmap=matplotlib.cm.coolwarm, alpha = 1, edgecolor = faces)
   	p4 = PatchCollection(patches,cmap=matplotlib.cm.coolwarm, alpha = 1, edgecolor = faces)


	p1.set_array(np.degrees(body.grid.solar_zenith_angle[t,:]))
	p2.set_array(np.degrees(body.grid.observer_zenith_angle))
	p3.set_array(np.degrees(body.grid.beta[t,:]))
	p4.set_array(np.degrees(body.grid.azimuth[t,:]))
    	aux1 = ax1.add_collection(p1)
    	fig.colorbar(aux1,ax=ax1,fraction=0.046, pad=0.04) 
    	aux2 = ax2.add_collection(p2)
    	fig.colorbar(aux2,ax=ax2,fraction=0.046, pad=0.04) 
    	aux3 = ax3.add_collection(p3)
    	fig.colorbar(aux3,ax=ax3,fraction=0.046, pad=0.04) 
    	aux4 = ax4.add_collection(p4)
    	fig.colorbar(aux4,ax=ax4,fraction=0.046, pad=0.04) 

	circle1 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
	circle2 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
	circle3 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
	circle4 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    	ax1.add_artist(circle1)   
    	ax2.add_artist(circle2)   
    	ax3.add_artist(circle3)   
    	ax4.add_artist(circle4)   

        ax1.set_title('SZA')
        ax2.set_title('EMISSION')
        ax3.set_title('BETA')
        ax4.set_title('PHI')
     
    	plot_config_r(ax1)
    	plot_config_r(ax2)
    	plot_config_r(ax3)
    	plot_config_r(ax4)
    	plt.tight_layout()

	if save:

		filename = 'geometry-pixels_' + planet1.properties.fourier_scene + '_alpha-' + str(int(np.round(np.degrees(body.geometry.phase_angle[t])))) + '_'  + time.strftime("%d-%m-%Y") + '_' + time.strftime("%X")

		fig.savefig('Images/'+filename + '.eps')
		fig.savefig('Images/'+filename + '.png')



def radiance(body, phase = False, save = False):

        # create a figure with subplots
    	fig = plt.figure(figsize=(9,12))
        ax1 = plt.subplot2grid((4,1), (0,0))
        ax2 = plt.subplot2grid((4,1), (1,0),sharex=ax1)
	ax3 = plt.subplot2grid((4,1), (2,0),sharex=ax1)
	ax4 = plt.subplot2grid((4,1), (3,0),sharex=ax1)

	ax1.grid()
	ax1.set_title(body.name + ' I parameter')
	ax1.set_xlabel('Time [Earth days]')
	if phase is False:
		ax1.plot(body.ephemeris.time/3600/24, body.radiance.I_ref,'-b')
	else:
		ax1.plot(np.rad2deg(body.geometry.phase_angle), body.radiance.I_ref,'-b')

	ax2.grid()
	ax2.set_title(body.name + ' Q parameter')
	ax2.set_xlabel('Time [Earth days]')
	if phase is False:
		ax2.plot(body.ephemeris.time/3600/24, body.radiance.Q_ref,'-g')
	else:
		ax2.plot(np.rad2deg(body.geometry.phase_angle), body.radiance.Q_ref,'-g')
	ax3.grid()
	ax3.set_title(body.name + ' U parameter')
        ax3.set_xlabel('Time [Earth days]')
	if phase is False:
		ax3.plot(body.ephemeris.time/3600/24, body.radiance.U_ref,'-r')
	else:
		ax3.plot(np.rad2deg(body.geometry.phase_angle), body.radiance.U_ref,'-r')
	
	ax4.grid()
	ax4.set_title(body.name + ' V parameter')
        ax4.set_xlabel('Time [Earth days]')
	if phase is False:
		ax4.plot(body.ephemeris.time/3600/24, body.radiance.V_ref,'-k')
	else:
		ax4.plot(np.rad2deg(body.geometry.phase_angle), body.radiance.V_ref,'-k')

	plt.tight_layout()

	if save:

		filename = 'radiance_' + body.properties.fourier_scene + '_alpha-' + str(int(np.round(np.degrees(body.geometry.phase_angle[t])))) + '_'  + time.strftime("%d-%m-%Y") + '_' + time.strftime("%X")

		fig.savefig('Images/'+filename + '.eps')
		fig.savefig('Images/'+filename + '.png')



def detail_radiance(body1,body2,I,Q,U,V, save = False):

        # create a figure with subplots
        ax3 = plt.subplot2grid((2,2), (0,0))
        ax4 = plt.subplot2grid((2,2), (0,1),sharey=ax3)
	ax5 = plt.subplot2grid((2,2), (1,0),colspan = 2)


	ax3.grid()
	ax3.set_title(body1.name + ' reflected light')
	ax3.set_xlabel('Time [Earth days]')
	ax3.plot(body1.ephemeris.time/3600/24, body1.radiance.I_ref,'-')
        #ax3.plot(body2.ephemeris.time, body2.radiance.I)
        ax3.plot(body1.ephemeris.time/3600/24, body1.radiance.Q_ref,'-')
        #ax3.plot(body2.ephemeris.time, body2.radiance.Q)
        ax3.plot(body1.ephemeris.time/3600/24, body1.radiance.U_ref,'-')
        #ax3.plot(body2.ephemeris.time, body2.radiance.U)
	
	ax4.grid()
	ax4.set_title(body2.name + ' reflected light')
	ax4.set_xlabel('Time [Earth days]')
	plt.setp(ax4.get_yticklabels(), visible=False)
	#ax4.plot(body1.ephemeris.time, body1.radiance.I_ref)
        ax4.plot(body2.ephemeris.time/3600/24, body2.radiance.I_ref,'-')
        #ax4.plot(body1.ephemeris.time, body1.radiance.Q_ref)
        ax4.plot(body2.ephemeris.time/3600/24, body2.radiance.Q_ref,'-')
        #ax4.plot(body1.ephemeris.time, body1.radiance.U_ref)
        ax4.plot(body2.ephemeris.time/3600/24, body2.radiance.U_ref,'-')

	ax5.grid()
        ax5.set_title('Total reflected light')
        ax5.set_xlabel('Time [Earth days]')
	ax5.plot(body2.ephemeris.time/3600/24, I,'-')
        ax5.plot(body2.ephemeris.time/3600/24, Q,'-')
        ax5.plot(body2.ephemeris.time/3600/24, U,'-')
	ax5.legend(['I','Q','U'])
	plt.tight_layout()

	if save:

		filename = 'detail-radiance_' + body1.properties.fourier_scene +'-'+ body2.properties.fourier_scene + '_alpha-' + str(int(np.round(np.degrees(body.geometry.phase_angle[t])))) + '_'  + time.strftime("%d-%m-%Y") + '_' + time.strftime("%X")

		fig.savefig('Images/'+filename + '.eps')
		fig.savefig('Images/'+filename + '.png')


def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        #f.canvas.manager.window.move(x, y)
        None
    plt.show()
    
def plot_color(string):
    
    color = np.array([[ 'faces'      , '#aaadab', 'b'  ],
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
    
    return color[np.where(color==string)[0][0],1+_cfg.plot_color]
    
def initialize_indicators(ax, m=2):
    circles = [0,0,0,0,0]
    # Phase
    circles[0] = plt.Circle((0.54*m, -0.54*m), 0.02*m, color = plot_color('phase_c')    ,ec='k'   , fill=False, zorder=1)
    # Umbra
    circles[1] = plt.Circle((0.54*m, -0.44*m), 0.01*m, color = plot_color('umbra_c')    ,ec='none', fill=False, zorder=2)
    # Antumbra
    circles[2] = plt.Circle((0.54*m, -0.44*m), 0.01*m, color = plot_color('antumbra_c') ,ec='none', fill=False, zorder=2)
    # Penumbra
    circles[3] = plt.Circle((0.54*m, -0.44*m), 0.02*m, color = plot_color('penumbra_c') ,ec='k'   , fill=False, zorder=1)
    # Transit
    circles[4] = plt.Circle((0.54*m, -0.49*m), 0.02*m, color = plot_color('transit_c')  ,ec='k'   , fill=False, zorder=1)
    
    ax.add_artist(circles[0])
    ax.add_artist(circles[1])
    ax.add_artist(circles[2])
    ax.add_artist(circles[3])
    ax.add_artist(circles[4])    
    
    return circles
    
    
def update_indicators(ax,flags,t,circles):
    
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
    
    
def load_flags(body,Type,time):
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
            
    aux = np.zeros_like(body.ephemeris.time, dtype = bool)
    if flags.PHASE:                  flags.phase    = np.ones(len(range(time[0],time[1])), dtype=bool)
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
    
    
def plot_config(ax):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False) 
    ax.set_axis_bgcolor((plot_color('background')))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([-1.14, 1.14])
    ax.set_ylim([-1.14, 1.20])    
    
def plot_config_r(ax):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False) 
    ax.set_axis_bgcolor((plot_color('background')))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([-1.01, 1.01])
    ax.set_ylim([-1.01, 1.01])  
    #ax.patch.set_visible(False)
    ax.axis('off')

def initialize_text(ax, Type):
    
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
        
def update_text(text, body, t, Type):
    
    #if Type == 'r':
        #text[3].set_text('Time   = %4d days, %02d hours (%d)'%
        #             (int(body.ephemeris.time[t]/60/60/24), (body.ephemeris.time[t]/24/60/60-
        #              int(body.ephemeris.time[t]/60/60/24))*24,t))
        #text[4].set_text('Phase = %04.2f$^{\circ}$' % (np.degrees(body.geometry.phase_angle[t]),))

    if Type == 'd':
        text[0].set_text('Points: '+ str(body.grid.N_points))
        #text[1].set_text('$\^{A}$= '+ str(np.round(np.sum(body.grid.area)/(m.pi*0.5**2),3)))
        text[2].set_text('$A_{shadow}$= %4.2f [-]' % (grid_area(body.grid,t)))
        text[3].set_text('t   = %4d days, %02d hours (%d)'%
                     (int(body.ephemeris.time[t]/60/60/24), (body.ephemeris.time[t]/24/60/60-
                      int(body.ephemeris.time[t]/60/60/24))*24,t))
        text[4].set_text('$Phase$  = %04.2f$^{\circ}$' % (np.degrees(body.geometry.phase_angle[t]),))
        text[5].set_text('$\\varphi_S$ = %04.2f$^{\circ}$' % (np.degrees(body.geometry.solar_azimuth_angle[t]),))   
#        
#        text[0].set_text('$Points: '+ str(body.grid.N_points)+'$')
#        text[1].set_text('$\^{A}= '+ str(np.round(np.sum(body.grid.area)/(m.pi*0.5**2),3))+'$')
#        text[2].set_text('$A_{shadow}= \, %4.2f \, [-]$' % (grid_area(body.grid,t)))
#        text[3].set_text('$t \,\,\,\,\, = \, %4d \, days, \, %02d \, hours \, (%d)$'%
#                     (int(body.ephemeris.time[t]/60/60/24), (body.ephemeris.time[t]/24/60/60-
#                      int(body.ephemeris.time[t]/60/60/24))*24,t))
#        text[4].set_text('$\\alpha \,\,\, = \, %04.2f^{\circ}$' % (np.degrees(body.geometry.phase_angle[t]),))
#        text[5].set_text('$\\varphi_S = \, %04.2f^{\circ}$' % (np.degrees(body.geometry.solar_azimuth_angle[t]),))   
        
    if Type == 'c':
        text[0].set_text('t   = %4d days, %02d hours (%d)'%
                     (int(body.ephemeris.time[t]/60/60/24), (body.ephemeris.time[t]/24/60/60-
                      int(body.ephemeris.time[t]/60/60/24))*24,t))
        text[1].set_text('$Phase$  = %04.2f$^{\circ}$' % (np.degrees(body.geometry.phase_angle[t]),))
        text[2].set_text('$\\varphi_S$ = %04.2f$^{\circ}$' % (np.degrees(body.geometry.solar_azimuth_angle[t]),))   
#        
#        text[0].set_text('$t \,\,\,\,\, = \, %4d \, days, \, %02d \, hours \, (%d)$'%
#                     (int(body.ephemeris.time[t]/60/60/24), (body.ephemeris.time[t]/24/60/60-
#                      int(body.ephemeris.time[t]/60/60/24))*24,t))
#        text[1].set_text('$Phase \,\,\, = \, %04.2f^{\circ}$' % (np.degrees(body.geometry.phase_angle[t]),))
#        text[2].set_text('$\\varphi_S = \, %04.2f^{\circ}$' % (np.degrees(body.geometry.solar_azimuth_angle[t]),))   


class plot_flags():
    
    def __init__(self,time):
    
        self.PHASE   = False
        self.TRANSIT = False
        self.ECLIPSE = False
        self.phase      = np.zeros(len(range(time[0],time[1])), dtype=bool)
        self.umbra      = np.zeros(len(range(time[0],time[1])), dtype=bool)
        self.antumbra   = np.zeros(len(range(time[0],time[1])), dtype=bool)
        self.penumbra   = np.zeros(len(range(time[0],time[1])), dtype=bool)
        self.transit    = np.zeros(len(range(time[0],time[1])), dtype=bool)


def plot_grid(body):
    
    print('\n    ⇒ Plotting grid of ' + body.type + ' ' + body.name)
    
    fig , ax = plt.subplots()
    plt.plot(body.grid.nodes[:,0], body.grid.nodes[:,1], 'o', color = plot_color('nodes'), markersize=3)
    for i in range(body.grid.N_points):
        plt.plot(body.grid.faces[i,0,:], body.grid.faces[i,1,:], plot_color('faces'), linewidth = 2)
    circle1 = plt.Circle((0, 0), 0.5, color = plot_color('circle'), fill=False, zorder=1)
    ax.add_artist(circle1)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([-0.55, 0.55])
    ax.set_ylim([-0.55, 0.6])
    ax.text(-0.4,  0.525,'Points: '+ str(body.grid.N_points), fontsize= 15)
    #ax.text(0.14,  0.525,'$\^{A}$: '+ str(np.round(np.sum(body.grid.area)/(m.pi*0.5**2),3)), fontsize= 15)
    

def plot_shadow_d(body, t = 0, save = False, dots = False):
          
    print('\n    ⇒ Plotting shadow of ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
          
    flags = load_flags(body,'d',[t,t+1])
    
    # Drawing
    fig = plt.figure(figsize=(13,6))
    move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)
    cell        = []#list(np.zeros(body.grid.N_points))
    patch_cells = []#list(np.zeros(body.grid.N_points))

    if _cfg.plot_faces == True:
    	ax.plot(body.grid.faces[:,0,:].T*2, body.grid.faces[:,1,:].T*2, color= plot_color('faces'), linewidth = 0.5,zorder=2)    
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)
    for i in range(body.grid.N_points):
        cell.append( [Polygon(2*body.grid.faces[i,:,:].T)] )
        patch_cells.append( PatchCollection(cell[i], alpha=1, color=str(body.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax.add_collection(patch_cells[i])     
        
        
    # Body circle 
    circle = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle) 
    
    # Indicators
    circles = initialize_indicators(ax)
    update_indicators(ax,flags,0, circles)
    
    # Sun point
    P = [ 1.08*np.cos(body.geometry.solar_azimuth_angle[t]) , 1.08*np.sin(body.geometry.solar_azimuth_angle[t]+np.pi)]
    circleP = plt.Circle((P[0], P[1]), 0.03, color = plot_color('sun'), fill=True, zorder=1)
    ax.add_artist(circleP)    
    
    # Text
    text_d = initialize_text(ax, 'd')
    update_text(text_d, body, t, 'd')

    ax.set_title('Discretized '+body.type+' '+body.name)
    
    plot_config(ax)
    plt.tight_layout()

    if save:

    	filename = 'shadow_' + body.properties.fourier_scene + '_alpha-' + str(int(np.round(np.degrees(body.geometry.phase_angle[t])))) + '_'  + time.strftime("%d-%m-%Y") + '_' + time.strftime("%X")

    	fig.savefig('Images/'+filename + '.eps')
    	fig.savefig('Images/'+filename + '.png')



def plot_shadow_dd(body, t = [0,0,0], save = False, dots = False):

    	print('\n    ⇒ Plotting shadow of ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t[0]])+', t ='+str(body.ephemeris.time[t[1]])+', and t ='+str(body.ephemeris.time[t[2]])+' seconds')
    	from matplotlib.patches import Rectangle

        # create a figure with subplots
    	fig = plt.figure(figsize=(11,4))
        ax1 = plt.subplot2grid((1,3), (0,0))
        ax2 = plt.subplot2grid((1,3), (0,1))
	ax3 = plt.subplot2grid((1,3), (0,2))
   
    	if _cfg.plot_faces == True:
    		faces = plot_color('faces')
    	else:
		faces = "none"  
    	if dots:
        	ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)

    	patches = []
    	for i in range(body.grid.N_points):
		square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
		patches.append(square)

   	p1 = PatchCollection(patches, alpha = 1, cmap = 'Greys_r', edgecolor = faces)
   	p2 = PatchCollection(patches, alpha = 1, cmap = 'Greys_r', edgecolor = faces)
   	p3 = PatchCollection(patches, alpha = 1, cmap = 'Greys_r', edgecolor = faces)

	p1.set_array(body.grid.shadow[t[0],:])
	p2.set_array(body.grid.shadow[t[1],:])
	p3.set_array(body.grid.shadow[t[2],:])
    	aux1 = ax1.add_collection(p1)
    	aux2 = ax2.add_collection(p2)
    	aux3 = ax3.add_collection(p3)

	circle1 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
	circle2 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
	circle3 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
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
     
    	plot_config_r(ax1)
    	plot_config_r(ax2)
    	plot_config_r(ax3)
    	plt.tight_layout()

        if save:

    		filename = 'shadow3_' + body.properties.fourier_scene + '_alpha-' + str(int(np.round(np.degrees(body.geometry.phase_angle[t[0]])))) + '-'+str(int(np.round(np.degrees(body.geometry.phase_angle[t[1]])))) + '-' + str(int(np.round(np.degrees(body.geometry.phase_angle[t[2]])))) + '_'  + time.strftime("%d-%m-%Y") + '_' + time.strftime("%X")
		
    		fig.savefig('Images/'+filename + '.eps')
    		fig.savefig('Images/'+filename + '.png')



def plot_I_d(body, t = 0, dots = False):
          
    print('\n    ⇒ Plotting I parameter ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    from matplotlib.patches import Rectangle

    flags = load_flags(body,'d',[t,t+1])
    
    # Drawing
    #fig = plt.figure(figsize=(12,10))
    fig = plt.figure(figsize=(10,8))
    move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)

    if _cfg.plot_faces == True:
    	faces = plot_color('faces')
    else:
	faces = "none"  
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)

    patches = []
    for i in range(body.grid.N_points):
	square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
	patches.append(square)

    p = PatchCollection(patches,cmap=matplotlib.cm.jet, alpha = 1, edgecolor = faces)
    p.set_array(body.grid.I[t,:])
    ax.add_collection(p)
    plt.colorbar(p) 

    # Body circle 
    circle = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle)   
    
    ax.set_title('I parameter')
    
    plot_config_r(ax)
    plt.tight_layout()
    #fig.patch.set_visible(False)




def plot_Q_d(body, t = 0, dots = False):
          
    print('\n    ⇒ Plotting Q parameter ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    from matplotlib.patches import Rectangle

    flags = load_flags(body,'d',[t,t+1])
    
    # Drawing
    #fig = plt.figure(figsize=(12,10))
    fig = plt.figure(figsize=(10,8))
    move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)

    if _cfg.plot_faces == True:
    	faces = plot_color('faces')
    else:
	faces = "none"  
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)

    patches = []
    for i in range(body.grid.N_points):
	square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
	patches.append(square)

    p = PatchCollection(patches,cmap=matplotlib.cm.jet, alpha = 1, edgecolor = faces)
    p.set_array(body.grid.Q[t,:])
    ax.add_collection(p)
    plt.colorbar(p) 

    # Body circle 
    circle = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle)   

    ax.set_title('Q parameter')
    
    plot_config_r(ax)
    plt.tight_layout()




def plot_U_d(body, t = 0, dots = False):
          
    print('\n    ⇒ Plotting U parameter ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    from matplotlib.patches import Rectangle

    flags = load_flags(body,'d',[t,t+1])
    
    # Drawing
    #fig = plt.figure(figsize=(12,10))
    fig = plt.figure(figsize=(10,8))
    move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)

    if _cfg.plot_faces == True:
    	faces = plot_color('faces')
    else:
	faces = "none"  
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)

    patches = []
    for i in range(body.grid.N_points):
	square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
	patches.append(square)

    p = PatchCollection(patches,cmap=matplotlib.cm.jet, alpha = 1, edgecolor = faces)
    p.set_array(body.grid.U[t,:])
    ax.add_collection(p)
    plt.colorbar(p) 

    # Body circle 
    circle = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle)   

    ax.set_title('U parameter')
    
    plot_config_r(ax)
    plt.tight_layout()



    

def plot_V_d(body, t = 0, dots = False):
          
    print('\n    ⇒ Plotting V parameter ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    from matplotlib.patches import Rectangle

    flags = load_flags(body,'d',[t,t+1])
    
    # Drawing
    #fig = plt.figure(figsize=(12,10))
    fig = plt.figure(figsize=(10,8))
    move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)

    if _cfg.plot_faces == True:
    	faces = plot_color('faces')
    else:
	faces = "none"  
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)

    patches = []
    for i in range(body.grid.N_points):
	square = Rectangle( (body.grid.faces[i,0,0]*2,body.grid.faces[i,1,0]*2),-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1],-2*body.grid.faces[i,0,0]+2*body.grid.faces[i,0,1] )
	patches.append(square)

    p = PatchCollection(patches,cmap=matplotlib.cm.jet, alpha = 1, edgecolor = faces)
    p.set_array(body.grid.V[t,:])
    ax.add_collection(p)
    plt.colorbar(p) 

    # Body circle 
    circle = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle)   

    ax.set_title('V parameter')
    
    plot_config_r(ax)
    plt.tight_layout()





def plot_radiance_d(body, t = 0, save = False, dots = False):

    	print('\n    ⇒ Plotting stokes parameters of ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    	from matplotlib.patches import Rectangle

        # create a figure with subplots
    	fig = plt.figure(figsize=(10,8))
        ax1 = plt.subplot2grid((2,2), (0,0))
        ax2 = plt.subplot2grid((2,2), (0,1))
	ax3 = plt.subplot2grid((2,2), (1,0))
	ax4 = plt.subplot2grid((2,2), (1,1))
    
	circle1 = plt.Circle((0, 0), 1, facecolor = plot_color('circle1'), edgecolor = plot_color('circle'), linewidth = 2.5)
	circle2 = plt.Circle((0, 0), 1, edgecolor = plot_color('circle'), facecolor = plot_color('circle1'), fill=True , linewidth = 2.5)
	circle3 = plt.Circle((0, 0), 1, edgecolor = plot_color('circle'), facecolor = plot_color('circle1'), fill=True , linewidth = 2.5)
	circle4 = plt.Circle((0, 0), 1, edgecolor = plot_color('circle'), facecolor = plot_color('circle1'), fill=True , linewidth = 2.5)
    	ax1.add_artist(circle1)   
    	ax2.add_artist(circle2)   
    	ax3.add_artist(circle3)   
    	ax4.add_artist(circle4)   

    	if _cfg.plot_faces == True:
    		faces = plot_color('faces')
    	else:
		faces = "none"  
    	if dots:
        	ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)

	N 	= np.sum(body.grid.shadow[t,:]!=0)
	faces00	= body.grid.faces[body.grid.shadow[t,:]!=0,0,0] 
	faces10	= body.grid.faces[body.grid.shadow[t,:]!=0,1,0] 
	faces01	= body.grid.faces[body.grid.shadow[t,:]!=0,0,1] 
	I	= body.grid.I[t,body.grid.shadow[t,:]!=0]
	Q	= body.grid.Q[t,body.grid.shadow[t,:]!=0]
	U	= body.grid.U[t,body.grid.shadow[t,:]!=0]
	V	= body.grid.V[t,body.grid.shadow[t,:]!=0]

    	patches = []
    	for i in range(N):
		square = Rectangle( (faces00[i]*2,faces10[i]*2),-2*faces00[i]+2*faces01[i],-2*faces00[i]+2*faces01[i] )
		patches.append(square)

   	p1 = PatchCollection(patches,cmap=matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)
   	p2 = PatchCollection(patches,cmap=matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)
   	p3 = PatchCollection(patches,cmap=matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)
   	p4 = PatchCollection(patches,cmap=matplotlib.cm.jet, alpha = 1, edgecolor = faces, zorder = 2)

	Q_aux = -Q/I
	U_aux = U/I
	V_aux = V/I

	Q_aux[np.isnan(Q_aux)] = 0
	U_aux[np.isnan(U_aux)] = 0
	V_aux[np.isnan(V_aux)] = 0

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
     
    	plot_config_r(ax1)
    	plot_config_r(ax2)
    	plot_config_r(ax3)
    	plot_config_r(ax4)
    	plt.tight_layout()

	circle5 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
	circle6 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
	circle7 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
	circle8 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    	ax1.add_artist(circle5)   
    	ax2.add_artist(circle6)   
    	ax3.add_artist(circle7)   
    	ax4.add_artist(circle8) 

	if save:
	
		filename = 'radiance-pixels_' + body.properties.fourier_scene + '_alpha-' + str(int(np.round(np.degrees(body.geometry.phase_angle[t])))) + '_'  + time.strftime("%d-%m-%Y") + '_' + time.strftime("%X")

		fig.savefig('Images/'+filename + '.eps')
		fig.savefig('Images/'+filename + '.png')








	

def plot_shadow_cd(body, t, dots=False):
    
    print('\n    ⇒ Plotting shadow of ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    
    flags_d = load_flags(body,'d',[t,t+1])
    flags_c = load_flags(body,'c',[t,t+1])
    
    fig = plt.figure(figsize=(13,6))
    move_figure(fig, 155, 110)
    ax_d = fig.add_subplot(1, 2, 1)
    ax_c = fig.add_subplot(1, 2, 2)

    plot_config(ax_d)
    plot_config(ax_c)
    
    # Drawing
    cell        = []#list(np.zeros(body.grid.N_points))
    patch_cells = []#list(np.zeros(body.grid.N_points))

#    ax_d.plot(body.grid.faces[:,0,:].T*2, body.grid.faces[:,1,:].T*2, '-',color=plot_color('faces'), linewidth = 0.5,zorder=2)    
#    if dots:
#        ax_d.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color= plot_color('nodes'), markersize=3)
#    for i, x in enumerate(np.arange(body.grid.N_points)[body.grid.shadow[t,:]==0]):
#        cell.append( Polygon(2*body.grid.faces[x,:,:].T) )
#    patch_cells.append( PatchCollection(cell, alpha=1, color=str(0),edgecolor='none', zorder=1) )
#    ax_d.add_collection(patch_cells[0])       
#    for i, x in enumerate(np.arange(body.grid.N_points)[(body.grid.shadow[t,:]!=1)&(body.grid.shadow[t,:]!=0)]):
#        cell.append( Polygon(2*body.grid.faces[x,:,:].T) )
#        patch_cells.append( PatchCollection(cell[x+1], alpha=1, color=str(body.grid.shadow[t,x]),edgecolor='none', zorder=1) )
#        ax_d.add_collection(patch_cells[x+1])  
        
        
    ax_d.plot(body.grid.faces[:,0,:].T*2, body.grid.faces[:,1,:].T*2, plot_color('faces'), linewidth = 0.5,zorder=2)    
    if dots:
        ax_d.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)
    for i in range(body.grid.N_points):
        cell.append( [Polygon(2*body.grid.faces[i,:,:].T)] )
        patch_cells.append( PatchCollection(cell[i], alpha=1, color=str(body.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax_d.add_collection(patch_cells[i])    
            
    if flags_c.PHASE:
        phase=PatchCollection([Polygon(body.geometry.phase_points[0:2,t,:].T)],
                               alpha=1,color=plot_color('phase'),edgecolor='none', zorder=4) 
        ax_c.add_collection(phase)    

    if flags_c.TRANSIT:
        aux_polygon = []
        for i in range(len(body.geometry.transit_points)):
            aux_polygon.append(Polygon(body.geometry.transit_points[i][0][0:2,t,:].T))
        transit=PatchCollection(aux_polygon, alpha=1,color=plot_color('transit'),edgecolor='none', zorder=4) 
        ax_c.add_collection(transit)     

    if flags_c.ECLIPSE:
        aux_antumbra = []
        aux_umbra    = []
        aux_penumbra = []
        for i in range(len(body.geometry.umbra_points)):        
            aux_umbra.append( Polygon(body.geometry.umbra_points[i][0][0:2,t,:].T) )
            aux_antumbra.append( Polygon(body.geometry.antumbra_points[i][0][0:2,t,:].T) )
            aux_penumbra.append( Polygon(body.geometry.penumbra_points[i][0][0:2,t,:].T) )                                      
        umbra    = PatchCollection(aux_umbra   , alpha=1,color=plot_color('umbra')   ,edgecolor='none', zorder=3)
        antumbra = PatchCollection(aux_antumbra, alpha=1,color=plot_color('antumbra'),edgecolor='none', zorder=3)
        penumbra = PatchCollection(aux_penumbra, alpha=1,color=plot_color('penumbra'),edgecolor='none', zorder=2)
        ax_c.add_collection(umbra)           
        ax_c.add_collection(antumbra)           
        ax_c.add_collection(penumbra)    
               
                   
    # Body circle 
    circle_d = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=5,linewidth = 2.5)
    ax_d.add_artist(circle_d) 
    circle_c = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=5,linewidth = 2.5)
    ax_c.add_artist(circle_c) 
    
    # Indicators
    circles_d = initialize_indicators(ax_d)
    update_indicators(ax_d,flags_d,0, circles_d)
    circles_c = initialize_indicators(ax_c)
    update_indicators(ax_c,flags_c,0, circles_c)
    
    # Sun point
    P = [ 1.08*np.cos(body.geometry.solar_azimuth_angle[t]+np.pi) , 1.08*np.sin(body.geometry.solar_azimuth_angle[t]+np.pi)]
    circleP_d = plt.Circle((P[0], P[1]), 0.03, color = plot_color('sun'), fill=True, zorder=5)
    circleP_c = plt.Circle((P[0], P[1]), 0.03, color = plot_color('sun'), fill=True, zorder=5)
    ax_d.add_artist(circleP_d)              
    ax_c.add_artist(circleP_c)      
           
    # Text
    text_d = initialize_text(ax_d, 'd')
    text_c = initialize_text(ax_c, 'c')
    update_text(text_d, body, t, 'd')
    update_text(text_c, body, t, 'c')
   
    # Title
    ax_d.set_title('Discretized '+body.type+' '+body.name)           
    ax_c.set_title('Continuous ' +body.type+' '+body.name)           
           

def plot_shadow_c(body, t):
    
    print('\n    ⇒ Plotting shadow of ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
    
    flags_c = load_flags(body,'c',[t,t+1])
    
    fig = plt.figure(figsize=(13,6))
    move_figure(fig, 155, 110)
    ax_c = fig.add_subplot(1, 1, 1)

    plot_config(ax_c)
    
    # Drawing
            
    if flags_c.PHASE:
        phase=PatchCollection([Polygon(body.geometry.phase_points[0:2,t,:].T)],
                               alpha=1,color=plot_color('phase'),edgecolor='none', zorder=4) 
        ax_c.add_collection(phase)    

    if flags_c.TRANSIT:
        aux_polygon = []
        for i in range(len(body.geometry.transit_points)):
            aux_polygon.append(Polygon(body.geometry.transit_points[i][0][0:2,t,:].T))
        transit=PatchCollection(aux_polygon, alpha=1,color=plot_color('transit'),edgecolor='none', zorder=4) 
        ax_c.add_collection(transit)     

    if flags_c.ECLIPSE:
        aux_antumbra = []
        aux_umbra    = []
        aux_penumbra = []
        for i in range(len(body.geometry.umbra_points)):        
            aux_umbra.append( Polygon(body.geometry.umbra_points[i][0][0:2,t,:].T) )
            aux_antumbra.append( Polygon(body.geometry.antumbra_points[i][0][0:2,t,:].T) )
            aux_penumbra.append( Polygon(body.geometry.penumbra_points[i][0][0:2,t,:].T) )                                      
        umbra    = PatchCollection(aux_umbra   , alpha=1,color=plot_color('umbra')   ,edgecolor='none', zorder=3)
        antumbra = PatchCollection(aux_antumbra, alpha=1,color=plot_color('antumbra'),edgecolor='none', zorder=3)
        penumbra = PatchCollection(aux_penumbra, alpha=1,color=plot_color('penumbra'),edgecolor='none', zorder=2)
        ax_c.add_collection(umbra)           
        ax_c.add_collection(antumbra)           
        ax_c.add_collection(penumbra)    
                  
                   
    # Body circle 
    circle_c = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=5,linewidth = 2.5)
    ax_c.add_artist(circle_c) 
    
    # Indicators
    circles_c = initialize_indicators(ax_c)
    update_indicators(ax_c,flags_c,0, circles_c)
    
    # Sun point
    P = [ 1.08*np.cos(body.geometry.solar_azimuth_angle[t]+np.pi) , 1.08*np.sin(body.geometry.solar_azimuth_angle[t]+np.pi)]
    circleP_c = plt.Circle((P[0], P[1]), 0.03, color = plot_color('sun'), fill=True, zorder=5)
    ax_c.add_artist(circleP_c)      
           
    # Text
    text_c = initialize_text(ax_c, 'c')
    update_text(text_c, body, t, 'c')
   
    # Title
    ax_c.set_title('Continuous ' +body.type+' '+body.name)           
           
           

def plot_shadow_cvsd(body, t, dots=False):
    
    print('\n    ⇒ Plotting shadow of ' + body.type + ' ' + body.name+' at t ='+str(body.ephemeris.time[t])+' seconds')
       
    flags_d = load_flags(body,'d',[t,t+1])
    flags_c = load_flags(body,'c',[t,t+1])
    
    fig = plt.figure(figsize=(13,6))
    move_figure(fig, 155, 110)
    ax_d = fig.add_subplot(1, 1, 1)

    plot_config(ax_d)
    
    # Drawing
    cell        = []#list(np.zeros(body.grid.N_points))
    patch_cells = []#list(np.zeros(body.grid.N_points))

#    ax_d.plot(body.grid.faces[:,0,:].T*2, body.grid.faces[:,1,:].T*2, plot_color('faces'), linewidth = 0.5,zorder=2)    
#    if dots:
#        ax_d.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)
#    for i, x in enumerate(np.arange(body.grid.N_points)[body.grid.shadow[t,:]==0]):
#        cell.append( Polygon(2*body.grid.faces[x,:,:].T) )
#    patch_cells.append( PatchCollection(cell, alpha=1, color=str(0),edgecolor='none', zorder=1) )
#    ax_d.add_collection(patch_cells[0])       
#    for i, x in enumerate(np.arange(body.grid.N_points)[(body.grid.shadow[t,:]!=1)&(body.grid.shadow[t,:]!=0)]):
#        cell.append( Polygon(2*body.grid.faces[x,:,:].T) )
#        patch_cells.append( PatchCollection([cell[x+1]], alpha=1, color=str(body.grid.shadow[t,x]),edgecolor='none', zorder=1) )
#        ax_d.add_collection(patch_cells[1+i])    
        
        
    ax_d.plot(body.grid.faces[:,0,:].T*2, body.grid.faces[:,1,:].T*2, plot_color('faces'), linewidth = 0.5,zorder=2)    
    if dots:
        ax_d.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)
    for i in range(body.grid.N_points):
        cell.append( [Polygon(2*body.grid.faces[i,:,:].T)] )
        patch_cells.append( PatchCollection(cell[i], alpha=1, color=str(body.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax_d.add_collection(patch_cells[i])            
            
    if flags_c.PHASE:
        ax_d.plot(body.geometry.phase_points[0,t,:], body.geometry.phase_points[1,t,:], color = plot_color('border'), linewidth = 2) 

    if flags_c.TRANSIT:
        for i in range(len(body.geometry.transit_points)):
            ax_d.plot(body.geometry.transit_points[i][0][0,t,:], body.geometry.transit_points[i][0][1,t,:], color = plot_color('border'), linewidth = 2)  

    if flags_c.ECLIPSE:
        for i in range(len(body.geometry.umbra_points)):        
            ax_d.plot(body.geometry.antumbra_points[i][0][0,t,:], body.geometry.antumbra_points[i][0][1,t,:], color = plot_color('border'), linewidth = 2) 
            ax_d.plot(body.geometry.umbra_points[i][0][0,t,:]   , body.geometry.umbra_points[i][0][1,t,:]   , color = plot_color('border'), linewidth = 2) 
            ax_d.plot(body.geometry.penumbra_points[i][0][0,t,:], body.geometry.penumbra_points[i][0][1,t,:], color = plot_color('border'), linewidth = 2)
       
    # Body circle 
    circle_d = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=5,linewidth = '2.5')
    ax_d.add_artist(circle_d) 
    
    # Indicators
    circles_d = initialize_indicators(ax_d)
    update_indicators(ax_d,flags_d,0, circles_d)
    
    # Sun point
    P = [ 1.08*np.cos(body.geometry.solar_azimuth_angle[t]+np.pi) , 1.08*np.sin(body.geometry.solar_azimuth_angle[t]+np.pi)]
    circleP_d = plt.Circle((P[0], P[1]), 0.03, color = plot_color('sun'), fill=True, zorder=5)
    ax_d.add_artist(circleP_d)              
           
    # Text
    text_d = initialize_text(ax_d, 'd')
    update_text(text_d, body, t, 'd')
   
    # Title
    ax_d.set_title('Discretized vs continuous '+body.type+' '+body.name)           

    
    
    
def plot_anim_shadow_d(body, dots = False, time  = 'all', info = [ 15, False]):            

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

    flags = load_flags(body,'d', time)
    
    # Drawing t=0
    fig = plt.figure(figsize=(13,6))
    move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)
    cell        = []#list(np.zeros(body.grid.N_points))
    patch_cells = []#list(np.zeros(body.grid.N_points))

    t_i = 0    
    t = time[t_i]
    t_vector = np.arange(time[0],time[1])
    
    ax.plot(body.grid.faces[:,0,:].T*2, body.grid.faces[:,1,:].T*2, plot_color('faces'), linewidth = 0.5,zorder=2)    
    if dots:
        ax.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)
    for i in range(body.grid.N_points):
        cell.append( [Polygon(2*body.grid.faces[i,:,:].T)] )
        patch_cells.append( PatchCollection(cell[i], alpha=1, color=str(body.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax.add_collection(patch_cells[i])       
        
    # Body circle 
    circle = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax.add_artist(circle) 
    
    # Indicators
    circles = initialize_indicators(ax)
    update_indicators(ax,flags,t, circles)
    
    # Sun point
    P = [ 1.08*np.cos(body.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi) , 1.08*np.sin(body.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi)]

    circleP = plt.Circle((P[0][t], P[1][t]), 0.03, color = plot_color('sun'), fill=True, zorder=1)
    ax.add_artist(circleP)    
    
    # Text
    text_d = initialize_text(ax, 'd')
    update_text(text_d, body, t, 'd')

    ax.set_title('Discretized '+body.type+' '+body.name)
    
    plot_config(ax)

    def animate(t_i):
        global t_vector
        global patch_cells
        
        t_i %= len(np.arange(time[0],time[1]))
        t    = t_vector[t_i]

        circleP.center = (P[0][t_i], P[1][t_i])
        update_indicators(ax,flags,t_i, circles)
        update_text(text_d, body, t, 'd')
        
        for i in range(body.grid.N_points):
            patch_cells[i].set_color(str(body.grid.shadow[t,i])) # set new color colors
         
    anim = animation.FuncAnimation(fig, animate, interval=info[0], blit=False)
    
    if info[1] is True:
        plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
        mywriter = animation.FFMpegWriter()
        anim.save('phase_anim.mp4', writer=mywriter)
    
    plt.show()
    return anim       

def plot_anim_shadow_d2(body1, body2, dots = False, time  = 'all', info = [ 15, False]):            

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

    flags1 = load_flags(body1,'d', time)
    flags2 = load_flags(body2,'d', time)
    
    # Drawing t=0
    fig = plt.figure(figsize=(13,6))
    move_figure(fig, 155, 110)
    ax1  = fig.add_subplot(1, 2, 1)
    ax2  = fig.add_subplot(1, 2, 2)
    cell1        = []#list(np.zeros(body1.grid.N_points))
    patch_cells1 = []#list(np.zeros(body1.grid.N_points))
    cell2        = []#list(np.zeros(body1.grid.N_points))
    patch_cells2 = []#list(np.zeros(body1.grid.N_points))
    
    t_i = 0    
    t = time[t_i]
    t_vector = np.arange(time[0],time[1])

    ax1.plot(body1.grid.faces[:,0,:].T*2, body1.grid.faces[:,1,:].T*2, plot_color('faces'), linewidth = 0.5,zorder=2)    
    if dots:
        ax1.plot(body1.grid.nodes[:,0]*2, body1.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)
    for i in range(body1.grid.N_points):
        cell1.append( [Polygon(2*body1.grid.faces[i,:,:].T)] )
        patch_cells1.append( PatchCollection(cell1[i], alpha=1, color=str(body1.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax1.add_collection(patch_cells1[i])       
    ax2.plot(body2.grid.faces[:,0,:].T*2, body2.grid.faces[:,1,:].T*2, plot_color('faces'), linewidth = 0.5,zorder=2)    
    if dots:
        ax2.plot(body2.grid.nodes[:,0]*2, body2.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)
    for i in range(body2.grid.N_points):
        cell2.append( [Polygon(2*body2.grid.faces[i,:,:].T)] )
        patch_cells2.append( PatchCollection(cell2[i], alpha=1, color=str(body2.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax2.add_collection(patch_cells2[i])   
    
    # Body circle 
    circle1 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=10,linewidth = 2.5)
    ax1.add_artist(circle1) 
    circle2 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=10,linewidth = 2.5)
    ax2.add_artist(circle2) 

    # Indicators
    circles1 = initialize_indicators(ax1)
    update_indicators(ax1,flags1,t_i, circles1)
    circles2 = initialize_indicators(ax2)
    update_indicators(ax2,flags2,t_i, circles2)

    # Sun point
    P1 = [ 1.08*np.cos(body1.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi) , 1.08*np.sin(body1.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi)]
    circleP1 = plt.Circle((P1[0][t_i], P1[1][t_i]), 0.03, color = plot_color('sun'), fill=True, zorder=1)
    ax1.add_artist(circleP1)
    P2 = [ 1.08*np.cos(body2.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi) , 1.08*np.sin(body2.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi)]
    circleP2 = plt.Circle((P2[0][t_i], P2[1][t_i]), 0.03, color = plot_color('sun'), fill=True, zorder=1)
    ax2.add_artist(circleP2)    
    
    # Text
    text_d1 = initialize_text(ax1, 'd')
    update_text(text_d1, body1, t, 'd')
    text_d2 = initialize_text(ax2, 'd')
    update_text(text_d2, body2, t, 'd')

    ax1.set_title('Discretized '+body1.type+' '+body1.name)
    ax2.set_title('Discretized '+body2.type+' '+body2.name)
    
    plot_config(ax1)
    plot_config(ax2)


    def animate(t_i):
        global t_vector
        global patch_cells1
        global patch_cells2
        
        t_i %= len(np.arange(time[0],time[1]))
        t    = t_vector[t_i]
        
        circleP1.center = (P1[0][t_i], P1[1][t_i])
        update_indicators(ax1,flags1,t_i, circles1)
        update_text(text_d1, body1, t, 'd')
        circleP2.center = (P2[0][t_i], P2[1][t_i])
        update_indicators(ax2,flags2,t_i, circles2)
        update_text(text_d2, body2, t, 'd')
        
        for i in range(body1.grid.N_points):
            patch_cells1[i].set_color(str(body1.grid.shadow[t,i])) # set new color colors
        for i in range(body2.grid.N_points):
            patch_cells2[i].set_color(str(body2.grid.shadow[t,i])) # set new color colors
            
    anim = animation.FuncAnimation(fig, animate, interval=info[0], blit=False)
    
    if info[1] is True:
        plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
        mywriter = animation.FFMpegWriter()
        anim.save('phase_anim.mp4', writer=mywriter)
        
    return anim       

    
    
    
    
def plot_anim_shadow_c(body, dots = False, time  = 'all', info = [ 15, False]):            

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
    global phase, transit, umbra, antumbra, penumbra
    global t_vector

    if time == 'all': time = [0,len(body.ephemeris.time)]

    print('\n    ⇒ Plotting shadow animation of ' + body.type + ' ' + body.name+' from t ='+str(body.ephemeris.time[time[0]])+' seconds to t = '+str(body.ephemeris.time[time[1]-1])+' seconds')

    flags = load_flags(body,'c', time)
    
    # Drawing t=0
    fig = plt.figure(figsize=(13,6))
    move_figure(fig, 155, 110)
    ax  = fig.add_subplot(1, 1, 1)

    t_i = 0    
    t = time[t_i]
    t_vector = np.arange(time[0],time[1])
    
    if flags.PHASE:
        phase=PatchCollection([Polygon(body.geometry.phase_points[0:2,t,:].T)],
                               alpha=1,color=plot_color('phase'),edgecolor='none', zorder=4) 
        ax.add_collection(phase)    

    if flags.TRANSIT:
        aux_polygon = []
        for i in range(len(body.geometry.transit_points)):
            aux_polygon.append(Polygon(body.geometry.transit_points[i][0][0:2,t,:].T))
        transit=PatchCollection(aux_polygon, alpha=1,color=plot_color('transit'),edgecolor='none', zorder=4) 
        ax.add_collection(transit)     

    if flags.ECLIPSE:
        aux_antumbra = []
        aux_umbra    = []
        aux_penumbra = []
        for i in range(len(body.geometry.umbra_points)):        
            aux_umbra.append( Polygon(body.geometry.umbra_points[i][0][0:2,t,:].T) )
            aux_antumbra.append( Polygon(body.geometry.antumbra_points[i][0][0:2,t,:].T) )
            aux_penumbra.append( Polygon(body.geometry.penumbra_points[i][0][0:2,t,:].T) )                                      
        umbra    = PatchCollection(aux_umbra   , alpha=1,color=plot_color('umbra')   ,edgecolor='none', zorder=3)
        antumbra = PatchCollection(aux_antumbra, alpha=1,color=plot_color('antumbra'),edgecolor='none', zorder=3)
        penumbra = PatchCollection(aux_penumbra, alpha=1,color=plot_color('penumbra'),edgecolor='none', zorder=2)
        ax.add_collection(umbra)           
        ax.add_collection(antumbra)           
        ax.add_collection(penumbra)    
         
        
    # Body circle 
    circle = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=10,linewidth = 2.5)
    ax.add_artist(circle) 
    
    # Indicators
    circles = initialize_indicators(ax)
    update_indicators(ax,flags,t, circles)
    
    # Sun point
    P = [ 1.08*np.cos(body.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi) , 1.08*np.sin(body.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi)]
    circleP = plt.Circle((P[0][t], P[1][t]), 0.03, color = plot_color('sun'), fill=True, zorder=1)
    ax.add_artist(circleP)    
    
    # Text
    text_d = initialize_text(ax, 'c')
    update_text(text_d, body, t, 'c')

    ax.set_title('Continuous '+body.type+' '+body.name)
    
    plot_config(ax)


    def animate(t_i):
        global t_vector
        global phase, transit, umbra, antumbra, penumbra
        
        t_i %= len(np.arange(time[0],time[1]))
        t    = t_vector[t_i]
        
        circleP.center = (P[0][t_i], P[1][t_i])
        update_indicators(ax,flags,t_i, circles)
        update_text(text_d, body, t, 'c')

        if flags.PHASE:
            ax.collections.remove(phase)
            phase=PatchCollection([Polygon(body.geometry.phase_points[0:2,t,:].T)],
                                   alpha=1,color=plot_color('phase'),edgecolor='none', zorder=4) 
            ax.add_collection(phase)    
    
        if flags.TRANSIT:
            ax.collections.remove(transit)
            aux_polygon = []
            for i in range(len(body.geometry.transit_points)):
                aux_polygon.append(Polygon(body.geometry.transit_points[i][0][0:2,t,:].T))
            transit=PatchCollection(aux_polygon, alpha=1,color=plot_color('transit'),edgecolor='none', zorder=4) 
            ax.add_collection(transit)     
    
        if flags.ECLIPSE:
            ax.collections.remove(umbra)
            ax.collections.remove(antumbra)
            ax.collections.remove(penumbra)
            aux_antumbra = []
            aux_umbra    = []
            aux_penumbra = []
            for i in range(len(body.geometry.umbra_points)):        
                aux_umbra.append( Polygon(body.geometry.umbra_points[i][0][0:2,t,:].T) )
                aux_antumbra.append( Polygon(body.geometry.antumbra_points[i][0][0:2,t,:].T) )
                aux_penumbra.append( Polygon(body.geometry.penumbra_points[i][0][0:2,t,:].T) )                                      
            umbra    = PatchCollection(aux_umbra   , alpha=1,color=plot_color('umbra')   ,edgecolor='none', zorder=3)
            antumbra = PatchCollection(aux_antumbra, alpha=1,color=plot_color('antumbra'),edgecolor='none', zorder=3)
            penumbra = PatchCollection(aux_penumbra, alpha=1,color=plot_color('penumbra'),edgecolor='none', zorder=2)
            ax.add_collection(umbra)           
            ax.add_collection(antumbra)           
            ax.add_collection(penumbra)    
                     
         
    anim = animation.FuncAnimation(fig, animate, interval=info[0], blit=False)
    
    if info[1] is True:
        plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
        mywriter = animation.FFMpegWriter()
        anim.save('phase_anim.mp4', writer=mywriter)
        
    return anim       

def plot_anim_shadow_c2(body1, body2, dots = False, time  = 'all', info = [ 15, False]):            

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
    global phase1, transit1, umbra1, antumbra1, penumbra1
    global phase2, transit2, umbra2, antumbra2, penumbra2
    global t_vector

    if time == 'all': time = [0,len(body1.ephemeris.time)]
        
    print('\n    ⇒ Plotting shadow animation of ' + body1.type + ' ' + body1.name+' and '+ body2.type + ' ' + body2.name +' from t ='+str(body1.ephemeris.time[time[0]])+' seconds to t = '+str(body2.ephemeris.time[time[1]-1])+' seconds')

    flags1 = load_flags(body1,'c', time)
    flags2 = load_flags(body2,'c', time)
    
    # Drawing t=0
    fig = plt.figure(figsize=(13,6))
    move_figure(fig, 155, 110)
    ax1  = fig.add_subplot(1, 2, 1)
    ax2  = fig.add_subplot(1, 2, 2)
    
    t_i = 0    
    t = time[t_i]
    t_vector = np.arange(time[0],time[1])

    if flags1.PHASE:
        phase1=PatchCollection([Polygon(body1.geometry.phase_points[0:2,t,:].T)],
                               alpha=1,color=plot_color('phase'),edgecolor='none', zorder=4) 
        ax1.add_collection(phase1)    

    if flags1.TRANSIT:
        aux_polygon = []
        for i in range(len(body1.geometry.transit_points)):
            aux_polygon.append(Polygon(body1.geometry.transit_points[i][0][0:2,t,:].T))
        transit1=PatchCollection(aux_polygon, alpha=1,color=plot_color('transit'),edgecolor='none', zorder=4) 
        ax1.add_collection(transit1)     

    if flags1.ECLIPSE:
        aux_antumbra = []
        aux_umbra    = []
        aux_penumbra = []
        for i in range(len(body1.geometry.umbra_points)):        
            aux_umbra.append( Polygon(body1.geometry.umbra_points[i][0][0:2,t,:].T) )
            aux_antumbra.append( Polygon(body1.geometry.antumbra_points[i][0][0:2,t,:].T) )
            aux_penumbra.append( Polygon(body1.geometry.penumbra_points[i][0][0:2,t,:].T) )                                      
        umbra1    = PatchCollection(aux_umbra   , alpha=1,color=plot_color('umbra')   ,edgecolor='none', zorder=3)
        antumbra1 = PatchCollection(aux_antumbra, alpha=1,color=plot_color('antumbra'),edgecolor='none', zorder=3)
        penumbra1 = PatchCollection(aux_penumbra, alpha=1,color=plot_color('penumbra'),edgecolor='none', zorder=2)
        ax1.add_collection(umbra1)           
        ax1.add_collection(antumbra1)           
        ax1.add_collection(penumbra1)    
         
        
    if flags2.PHASE:
        phase2=PatchCollection([Polygon(body2.geometry.phase_points[0:2,t,:].T)],
                               alpha=1,color=plot_color('phase'),edgecolor='none', zorder=4) 
        ax2.add_collection(phase2)    

    if flags2.TRANSIT:
        aux_polygon = []
        for i in range(len(body2.geometry.transit_points)):
            aux_polygon.append(Polygon(body2.geometry.transit_points[i][0][0:2,t,:].T))
        transit2=PatchCollection(aux_polygon, alpha=1,color=plot_color('transit'),edgecolor='none', zorder=4) 
        ax2.add_collection(transit2)     

    if flags2.ECLIPSE:
        aux_antumbra = []
        aux_umbra    = []
        aux_penumbra = []
        for i in range(len(body2.geometry.umbra_points)):        
            aux_umbra.append( Polygon(body2.geometry.umbra_points[i][0][0:2,t,:].T) )
            aux_antumbra.append( Polygon(body2.geometry.antumbra_points[i][0][0:2,t,:].T) )
            aux_penumbra.append( Polygon(body2.geometry.penumbra_points[i][0][0:2,t,:].T) )                                      
        umbra2    = PatchCollection(aux_umbra   , alpha=1,color=plot_color('umbra')   ,edgecolor='none', zorder=3)
        antumbra2 = PatchCollection(aux_antumbra, alpha=1,color=plot_color('antumbra'),edgecolor='none', zorder=3)
        penumbra2 = PatchCollection(aux_penumbra, alpha=1,color=plot_color('penumbra'),edgecolor='none', zorder=2)
        ax2.add_collection(umbra2)           
        ax2.add_collection(antumbra2)           
        ax2.add_collection(penumbra2)    
           
        
    
    # Body circle 
    circle1 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax1.add_artist(circle1) 
    circle2 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax2.add_artist(circle2) 

    # Indicators
    circles1 = initialize_indicators(ax1)
    update_indicators(ax1,flags1,t_i, circles1)
    circles2 = initialize_indicators(ax2)
    update_indicators(ax2,flags2,t_i, circles2)

    # Sun point
    P1 = [ 1.08*np.cos(body1.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi) , 1.08*np.sin(body1.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi)]
    circleP1 = plt.Circle((P1[0][t_i], P1[1][t_i]), 0.03, color = plot_color('sun'), fill=True, zorder=1)
    ax1.add_artist(circleP1)
    P2 = [ 1.08*np.cos(body2.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi) , 1.08*np.sin(body2.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi)]
    circleP2 = plt.Circle((P2[0][t_i], P2[1][t_i]), 0.03, color = plot_color('sun'), fill=True, zorder=1)
    ax2.add_artist(circleP2)    
    
    # Text
    text_d1 = initialize_text(ax1, 'c')
    update_text(text_d1, body1, t, 'c')
    text_d2 = initialize_text(ax2, 'c')
    update_text(text_d2, body2, t, 'c')

    ax1.set_title('Continuous '+body1.type+' '+body1.name)
    ax2.set_title('Continuous '+body2.type+' '+body2.name)
    
    plot_config(ax1)
    plot_config(ax2)


    def animate(t_i):
        global phase1, transit1, umbra1, antumbra1, penumbra1
        global phase2, transit2, umbra2, antumbra2, penumbra2
        global t_vector
        
        t_i %= len(np.arange(time[0],time[1]))
        t    = t_vector[t_i]
        
        circleP1.center = (P1[0][t_i], P1[1][t_i])
        update_indicators(ax1,flags1,t_i, circles1)
        update_text(text_d1, body1, t, 'c')
        circleP2.center = (P2[0][t_i], P2[1][t_i])
        update_indicators(ax2,flags2,t_i, circles2)
        update_text(text_d2, body2, t, 'c')
        
        if flags1.PHASE:
            ax1.collections.remove(phase1)
            phase1=PatchCollection([Polygon(body1.geometry.phase_points[0:2,t,:].T)],
                                   alpha=1,color=plot_color('phase'),edgecolor='none', zorder=4) 
            ax1.add_collection(phase1)    
    
        if flags1.TRANSIT:
            ax1.collections.remove(transit1)
            aux_polygon = []
            for i in range(len(body1.geometry.transit_points)):
                aux_polygon.append(Polygon(body1.geometry.transit_points[i][0][0:2,t,:].T))
            transit1=PatchCollection(aux_polygon, alpha=1,color=plot_color('transit'),edgecolor='none', zorder=4) 
            ax1.add_collection(transit1)     
    
        if flags1.ECLIPSE:
            ax1.collections.remove(umbra1)
            ax1.collections.remove(antumbra1)
            ax1.collections.remove(penumbra1)
            aux_antumbra = []
            aux_umbra    = []
            aux_penumbra = []
            for i in range(len(body1.geometry.umbra_points)):        
                aux_umbra.append( Polygon(body1.geometry.umbra_points[i][0][0:2,t,:].T) )
                aux_antumbra.append( Polygon(body1.geometry.antumbra_points[i][0][0:2,t,:].T) )
                aux_penumbra.append( Polygon(body1.geometry.penumbra_points[i][0][0:2,t,:].T) )                                      
            umbra1    = PatchCollection(aux_umbra   , alpha=1,color=plot_color('umbra')   ,edgecolor='none', zorder=3)
            antumbra1 = PatchCollection(aux_antumbra, alpha=1,color=plot_color('antumbra'),edgecolor='none', zorder=3)
            penumbra1 = PatchCollection(aux_penumbra, alpha=1,color=plot_color('penumbra'),edgecolor='none', zorder=2)
            ax1.add_collection(umbra1)           
            ax1.add_collection(antumbra1)           
            ax1.add_collection(penumbra1)    
            
            
        if flags2.PHASE:
            ax2.collections.remove(phase2)
            phase2=PatchCollection([Polygon(body2.geometry.phase_points[0:2,t,:].T)],
                                   alpha=1,color=plot_color('phase'),edgecolor='none', zorder=4) 
            ax2.add_collection(phase2)    
    
        if flags2.TRANSIT:
            ax2.collections.remove(transit2)
            aux_polygon = []
            for i in range(len(body2.geometry.transit_points)):
                aux_polygon.append(Polygon(body2.geometry.transit_points[i][0][0:2,t,:].T))
            transit2=PatchCollection(aux_polygon, alpha=1,color=plot_color('transit'),edgecolor='none', zorder=4) 
            ax2.add_collection(transit2)     
    
        if flags2.ECLIPSE:
            ax2.collections.remove(umbra2)
            ax2.collections.remove(antumbra2)
            ax2.collections.remove(penumbra2)
            aux_antumbra = []
            aux_umbra    = []
            aux_penumbra = []
            for i in range(len(body2.geometry.umbra_points)):        
                aux_umbra.append( Polygon(body2.geometry.umbra_points[i][0][0:2,t,:].T) )
                aux_antumbra.append( Polygon(body2.geometry.antumbra_points[i][0][0:2,t,:].T) )
                aux_penumbra.append( Polygon(body2.geometry.penumbra_points[i][0][0:2,t,:].T) )                                      
            umbra2    = PatchCollection(aux_umbra   , alpha=1,color=plot_color('umbra')   ,edgecolor='none', zorder=3)
            antumbra2 = PatchCollection(aux_antumbra, alpha=1,color=plot_color('antumbra'),edgecolor='none', zorder=3)
            penumbra2 = PatchCollection(aux_penumbra, alpha=1,color=plot_color('penumbra'),edgecolor='none', zorder=2)
            ax2.add_collection(umbra2)           
            ax2.add_collection(antumbra2)           
            ax2.add_collection(penumbra2)    
              
            
    anim = animation.FuncAnimation(fig, animate, interval=info[0], blit=False)
    
    if info[1] is True:
        plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
        mywriter = animation.FFMpegWriter()
        anim.save('phase_anim.mp4', writer=mywriter)
        
    return anim       
    
    
    
    
def plot_anim_shadow_cd(body,dots = False, time  = 'all', info = [ 15, False]):          

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
    global patch_cells1
    global phase2, transit2, umbra2, antumbra2, penumbra2
    global t_vector

    if time == 'all': time = [0,len(body.ephemeris.time)]
        
    print('\n    ⇒ Plotting shadow animation of ' + body.type + ' ' + body.name+' from t ='+str(body.ephemeris.time[time[0]])+' seconds to t = '+str(body.ephemeris.time[time[1]-1])+' seconds')

    flags1 = load_flags(body,'d', time)
    flags2 = load_flags(body,'c', time)
    
    # Drawing t=0
    fig = plt.figure(figsize=(13,6))
    move_figure(fig, 155, 110)
    ax1  = fig.add_subplot(1, 2, 1)
    ax2  = fig.add_subplot(1, 2, 2)
    
    t_i = 0    
    t = time[t_i]
    t_vector = np.arange(time[0],time[1])
    patch_cells1 = []
    cell1 = []
    
    ax1.plot(body.grid.faces[:,0,:].T*2, body.grid.faces[:,1,:].T*2, plot_color('faces'), linewidth = 0.5,zorder=2)    
    if dots:
        ax1.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)
    for i in range(body.grid.N_points):
        cell1.append( [Polygon(2*body.grid.faces[i,:,:].T)] )
        patch_cells1.append( PatchCollection(cell1[i], alpha=1, color=str(body.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax1.add_collection(patch_cells1[i])  
        
    if flags2.PHASE:
        phase2=PatchCollection([Polygon(body.geometry.phase_points[0:2,t,:].T)],
                               alpha=1,color=plot_color('phase'),edgecolor='none', zorder=4) 
        ax2.add_collection(phase2)    

    if flags2.TRANSIT:
        aux_polygon = []
        for i in range(len(body.geometry.transit_points)):
            aux_polygon.append(Polygon(body.geometry.transit_points[i][0][0:2,t,:].T))
        transit2=PatchCollection(aux_polygon, alpha=1,color=plot_color('transit'),edgecolor='none', zorder=4) 
        ax2.add_collection(transit2)     

    if flags2.ECLIPSE:
        aux_antumbra = []
        aux_umbra    = []
        aux_penumbra = []
        for i in range(len(body.geometry.umbra_points)):        
            aux_umbra.append( Polygon(body.geometry.umbra_points[i][0][0:2,t,:].T) )
            aux_antumbra.append( Polygon(body.geometry.antumbra_points[i][0][0:2,t,:].T) )
            aux_penumbra.append( Polygon(body.geometry.penumbra_points[i][0][0:2,t,:].T) )                                      
        umbra2    = PatchCollection(aux_umbra   , alpha=1,color=plot_color('umbra')   ,edgecolor='none', zorder=3)
        antumbra2 = PatchCollection(aux_antumbra, alpha=1,color=plot_color('antumbra'),edgecolor='none', zorder=3)
        penumbra2 = PatchCollection(aux_penumbra, alpha=1,color=plot_color('penumbra'),edgecolor='none', zorder=2)
        ax2.add_collection(umbra2)           
        ax2.add_collection(antumbra2)           
        ax2.add_collection(penumbra2)    
                  
        
    
    # Body circle 
    circle1 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax1.add_artist(circle1) 
    circle2 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax2.add_artist(circle2) 

    # Indicators
    circles1 = initialize_indicators(ax1)
    update_indicators(ax1,flags1,t_i, circles1)
    circles2 = initialize_indicators(ax2)
    update_indicators(ax2,flags2,t_i, circles2)

    # Sun point
    P1 = [ 1.08*np.cos(body.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi) , 1.08*np.sin(body.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi)]
    circleP1 = plt.Circle((P1[0][t_i], P1[1][t_i]), 0.03, color = plot_color('sun'), fill=True, zorder=1)
    ax1.add_artist(circleP1)
    circleP2 = plt.Circle((P1[0][t_i], P1[1][t_i]), 0.03, color = plot_color('sun'), fill=True, zorder=1)
    ax2.add_artist(circleP2)    
    
    # Text
    text_d1 = initialize_text(ax1, 'd')
    update_text(text_d1, body, t, 'd')
    text_d2 = initialize_text(ax2, 'c')
    update_text(text_d2, body, t, 'c')

    ax1.set_title('Discretized '+body.type+' '+body.name)
    ax2.set_title('Continuous '+body.type+' '+body.name)
    
    plot_config(ax1)
    plot_config(ax2)


    def animate(t_i):
        global patch_cells1
        global phase2, transit2, umbra2, antumbra2, penumbra2
        global t_vector
        
        t_i %= len(np.arange(time[0],time[1]))
        t    = t_vector[t_i]
        
        circleP1.center = (P1[0][t_i], P1[1][t_i])
        update_indicators(ax1,flags1,t_i, circles1)
        update_text(text_d1, body, t, 'd')
        circleP2.center = (P1[0][t_i], P1[1][t_i])
        update_indicators(ax2,flags2,t_i, circles2)
        update_text(text_d2, body, t, 'c')
        
        for i in range(body.grid.N_points):
            patch_cells1[i].set_color(str(body.grid.shadow[t,i])) # set new color colors
                
        if flags2.PHASE:
            ax2.collections.remove(phase2)
            phase2=PatchCollection([Polygon(body.geometry.phase_points[0:2,t,:].T)],
                                   alpha=1,color=plot_color('phase'),edgecolor='none', zorder=4) 
            ax2.add_collection(phase2)    
    
        if flags2.TRANSIT:
            ax2.collections.remove(transit2)
            aux_transit = []
            for i in range(len(body.geometry.transit_points)):
                aux_transit.append(Polygon(body.geometry.transit_points[i][0][0:2,t,:].T))
            transit2=PatchCollection(aux_transit, alpha=1,color=plot_color('transit'),edgecolor='none', zorder=4) 
            ax2.add_collection(transit2)     
    
        if flags2.ECLIPSE:
            ax2.collections.remove(umbra2)
            ax2.collections.remove(antumbra2)
            ax2.collections.remove(penumbra2)
            aux_antumbra = []
            aux_umbra    = []
            aux_penumbra = []
            for i in range(len(body.geometry.umbra_points)):        
                aux_umbra.append( Polygon(body.geometry.umbra_points[i][0][0:2,t,:].T) )
                aux_antumbra.append( Polygon(body.geometry.antumbra_points[i][0][0:2,t,:].T) )
                aux_penumbra.append( Polygon(body.geometry.penumbra_points[i][0][0:2,t,:].T) )                                      
            umbra2    = PatchCollection(aux_umbra   , alpha=1,color=plot_color('umbra')   ,edgecolor='none', zorder=3)
            antumbra2 = PatchCollection(aux_antumbra, alpha=1,color=plot_color('antumbra'),edgecolor='none', zorder=3)
            penumbra2 = PatchCollection(aux_penumbra, alpha=1,color=plot_color('penumbra'),edgecolor='none', zorder=2)
            ax2.add_collection(umbra2)           
            ax2.add_collection(antumbra2)           
            ax2.add_collection(penumbra2)    
            
    anim = animation.FuncAnimation(fig, animate, interval=info[0], blit=False)
    
    if info[1] is True:
        plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
        mywriter = animation.FFMpegWriter()
        anim.save('phase_anim.mp4', writer=mywriter)
        
    return anim      




    
def plot_anim_shadow_cvsd(body,dots = False, time  = 'all', info = [ 15, False]):          

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
    global line_phase, line_umbra, line_antumbra, line_penumbra, line_transit
    global patch_cells1
    global t_vector

    if time == 'all': time = [0,len(body.ephemeris.time)]
        
    print('\n    ⇒ Plotting shadow animation of ' + body.type + ' ' + body.name+' from t ='+str(body.ephemeris.time[time[0]])+' seconds to t = '+str(body.ephemeris.time[time[1]-1])+' seconds')

    flags1 = load_flags(body,'d', time)
    flags2 = load_flags(body,'c', time)
    
    # Drawing t=0
    fig = plt.figure(figsize=(13,6))
    move_figure(fig, 155, 110)
    ax1  = fig.add_subplot(1, 1, 1)
    
    t_i = 0    
    t = time[t_i]
    t_vector = np.arange(time[0],time[1])
    patch_cells1 = []
    cell1 = []
    
    ax1.plot(body.grid.faces[:,0,:].T*2, body.grid.faces[:,1,:].T*2, plot_color('faces'), linewidth = 0.5,zorder=2)    
    if dots:
        ax1.plot(body.grid.nodes[:,0]*2, body.grid.nodes[:,1]*2, 'o', color = plot_color('nodes'), markersize=3)
    for i in range(body.grid.N_points):
        cell1.append( [Polygon(2*body.grid.faces[i,:,:].T)] )
        patch_cells1.append( PatchCollection(cell1[i], alpha=1, color=str(body.grid.shadow[t,i]), edgecolor='none', zorder=1) )
        ax1.add_collection(patch_cells1[i])  
        
    if flags2.PHASE:
        line_phase, = ax1.plot(body.geometry.phase_points[0,t,:], body.geometry.phase_points[1,t,:], color = plot_color('border'), linewidth = 2) 

    if flags2.TRANSIT:
        line_transit = [] 
        for i in range(len(body.geometry.transit_points)):
            aux, = ax1.plot(body.geometry.transit_points[i][0][0,t,:], body.geometry.transit_points[i][0][1,t,:], color = plot_color('border'), linewidth = 2)  
            line_transit.append(aux)
            
    if flags2.ECLIPSE:
        line_antumbra = []
        line_umbra    = []
        line_penumbra = []
        for i in range(len(body.geometry.umbra_points)):        
            aux_antumbra, = ax1.plot(body.geometry.antumbra_points[i][0][0,t,:], body.geometry.antumbra_points[i][0][1,t,:], color = plot_color('border'), linewidth = 2) 
            aux_umbra,    = ax1.plot(body.geometry.umbra_points[i][0][0,t,:]   , body.geometry.umbra_points[i][0][1,t,:]   , color = plot_color('border'), linewidth = 2) 
            aux_penumbra, = ax1.plot(body.geometry.penumbra_points[i][0][0,t,:], body.geometry.penumbra_points[i][0][1,t,:], '-', color = plot_color('border'), linewidth = 2) 
            line_antumbra.append(aux_antumbra)
            line_umbra.append(aux_umbra)
            line_penumbra.append(aux_penumbra)
        
    # Body circle 
    circle1 = plt.Circle((0, 0), 1, color = plot_color('circle'), fill=False, zorder=3,linewidth = 2.5)
    ax1.add_artist(circle1) 

    # Indicators
    circles1 = initialize_indicators(ax1)
    update_indicators(ax1,flags1,t_i, circles1)

    # Sun point
    P1 = [ 1.08*np.cos(body.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi) , 1.08*np.sin(body.geometry.solar_azimuth_angle[time[0]:time[1]]+np.pi)]
    circleP1 = plt.Circle((P1[0][t_i], P1[1][t_i]), 0.03, color = plot_color('sun'), fill=True, zorder=1)
    ax1.add_artist(circleP1)   
    
    # Text
    text_d1 = initialize_text(ax1, 'd')
    update_text(text_d1, body, t, 'd')

    ax1.set_title('Discretized and continuous '+body.type+' '+body.name)
    
    plot_config(ax1)

    def animate(t_i):
        global line_phase, line_umbra, line_antumbra, line_penumbra, line_transit
        global patch_cells1
        global t_vector
        
        t_i %= len(np.arange(time[0],time[1]))
        t    = t_vector[t_i]
        
        circleP1.center = (P1[0][t_i], P1[1][t_i])
        update_indicators(ax1,flags1,t_i, circles1)
        update_text(text_d1, body, t, 'd')
        
        for i in range(body.grid.N_points):
            patch_cells1[i].set_color(str(body.grid.shadow[t,i])) # set new color colors
                
        if flags2.PHASE:
            line_phase.set_data(body.geometry.phase_points[0,t,:], body.geometry.phase_points[1,t,:])
            
        if flags2.TRANSIT:
            for i in range(len(body.geometry.transit_points)):
                line_transit[i].set_data(body.geometry.transit_points[i][0][0,t,:],body.geometry.transit_points[i][0][1,t,:])
    
        if flags2.ECLIPSE:
            for i in range(len(body.geometry.umbra_points)):               
                line_antumbra[i].set_data( body.geometry.antumbra_points[i][0][0,t,:], body.geometry.antumbra_points[i][0][1,t,:])
                line_umbra[i].set_data   ( body.geometry.umbra_points[i][0][0,t,:]   , body.geometry.umbra_points[i][0][1,t,:]   )
                line_penumbra[i].set_data( body.geometry.penumbra_points[i][0][0,t,:], body.geometry.penumbra_points[i][0][1,t,:])
            
        return# line_phase, line_transit[1], line_transit[0], line_umbra, line_antumbra, line_penumbra,
        
    anim = animation.FuncAnimation(fig, animate, interval=info[0], blit=False)
    
    if info[1] is True:
        plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
        mywriter = animation.FFMpegWriter()
        anim.save('phase_anim.mp4', writer=mywriter)
        
    return anim




        
    
    
    
    
def ecl(t, body1, body2, star):
    
    import sys
    sys.path.append('/home/javier/anaconda3/lib/python3.5/site-packages')
    import matplotlib.pyplot as plt
    #from shapely.geometry import Polygon as Poly
    from matplotlib.collections import PatchCollection
    import math as m
    #from descartes import PolygonPatch
    import numpy as np

    R1 = np.zeros(3)
    R2 = np.zeros(3)
    
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
    print('PSI', np.degrees(psi))
    omega = m.asin((r2+rs)/d2s)
    print('OMEGA', np.degrees(omega))
    
    O2A6 = (r1+r2)/m.sin(omega) 
    print('O2A6: ',O2A6)
    print('d2: ',d2)
    print('d1: ',d1)
    rho = m.atan(d1/(d2+O2A6))
    print('RHO', np.degrees(rho))
    if rho>omega:
        print('rho larger than omega')
    else:
        print('omega larger than rho')
    
    
    fig = plt.figure()
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
    #p = PatchCollection(penumbra, alpha=1, color='0.5', edgecolor=None, zorder = 1)
    #ax.add_collection(p)  
    
    antumbra = []
    #antumbra1 = Poly([(d2-r2/m.sin(psi),0),(d2-r2/m.sin(psi)-2*ds*m.cos(psi),2*ds*m.sin(psi)),(d2-r2/m.sin(psi)-2*ds*m.cos(psi),-2*ds*m.sin(psi))])
    #antumbra.append(PolygonPatch(antumbra1))
    #p = PatchCollection(antumbra, alpha=1, color='0.75', edgecolor=None, zorder = 1)
    #ax.add_collection(p)  
    
    umbra = []
    #umbra1 = Poly([(d2-r2/m.sin(psi),0) , ( d2-r2*m.sin(psi),r2*m.cos(psi) ) , (d2,0) , ( d2-r2*m.sin(psi),-r2*m.cos(psi) )])
    #umbra.append(PolygonPatch(umbra1))
    #p = PatchCollection(umbra, alpha=1, color='0.1', edgecolor=None, zorder = 1)
    #ax.add_collection(p)  
    
    circ=plt.Circle((ds,0), radius=rs, color='b', ec= 'k', fill=True)
    ax.add_patch(circ)
    circ=plt.Circle((d2,0), radius=r2, color='r', ec = 'k', fill=True)
    ax.add_patch(circ)
    circ=plt.Circle((0,d1), radius=r1, color='g', ec = 'k', fill=True)
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
    #p = PatchCollection(patches, color='red', alpha=0.4)
    #ax.add_collection(p)
    
    ax.set_xlim([-r1, d2])
    ax.set_ylim([-r2*20,  r2*20])        
    #ax.set_aspect('equal', adjustable='box')
    ax.set_axis_bgcolor('0.9')
    
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def plot_xyorbit(position2D_1, position2D_2 = None,  info =
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

    fig = plt.figure()
    
    # The orbits of the two bodies are plotted, if available
    plt.plot(position2D_1[0], position2D_1[1],'r', linewidth=2,
                                                 linestyle="-", label="Body 1")
    if position2D_2 is not None:
        plt.plot(position2D_2[0], position2D_2[1],'b', linewidth=2,
                                                 linestyle="-", label="Body 2")
    
    # Origin point and name
    plt.plot(0,0, 'ko') 
    ax = fig.gca()
    ax.annotate(info[4], xy=(0, 0), xytext=(0,0),   
                arrowprops = dict(facecolor='black', shrink=0.05),
                horizontalalignment='right', verticalalignment='top',
                color='black', size= 'large',)
                
    # Axes labels and title are established
    plt.xlabel(info[0])
    plt.ylabel(info[1])
    plt.title(info[2])
    
    # A same scale is applied to both axes
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid()
    
    # Axes limits and legend depend on the number of curves plotted
    if position2D_2 is not None:
        aux1 = max([position2D_1.max(),position2D_2.max()])
        plt.xlim([min([position2D_1[0].min(),position2D_2[0].min()]),
                           max([position2D_1[0].max(),position2D_2[0].max()])])
        plt.ylim([min([position2D_1[1].min(),position2D_2[1].min()]),
                           max([position2D_1[1].max(),position2D_2[1].max()])])
    
        if info[3] is not None:
            plt.legend([info[3][0],info[3][1]],loc='upper left', frameon=False)
        else:
            plt.legend(loc='upper left', frameon=False)
    else:
        aux1 = position2D_1.max()                
        plt.xlim([position2D_1[0].min(),position2D_1[0].max()])
        plt.ylim([position2D_1[1].min(),position2D_1[1].max()])
    
        if info[3] is not None:
            plt.legend([info[3]],loc='upper left', frameon=False)
    
    # Reference frame arrows are added
    ax.arrow( 0, 0, aux1/3, 0, fc="r", ec="r", 
                                     head_width=aux1/20, head_length=aux1/20 )
    ax.arrow( 0, 0, 0, aux1/3, fc="g", ec="g",
                                     head_width=aux1/20, head_length=aux1/20 )
    
    plt.show()
    
    
    
    
def plot_XYZorbit(position3D_1, position3D_2 = None, info = 
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
    import module_functions as fun

    fig = plt.figure()
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
    plt.title(info[3])
    
    fig.gca().set_aspect('equal', adjustable='box')
    plt.grid()

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
            plt.legend([info[4][0],info[4][1]],loc='center left',frameon=False)
        else:
            plt.legend(loc='center left', frameon=False)
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
            plt.legend([info[4]],loc='center left', frameon=False)
            
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
    ax.scatter([0],[0],[0], color="k",s=100)
    ax.text(aux/100, aux/100, aux/100, info[5], [0,1,0])

    plt.show()
     
     
     
     
     
     
def plot_XYZorbitanimation(time, position3D_1, position3D_2 = None, info = 
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
     
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d.axes3d as p3
    import matplotlib.animation as animation
    import statistics as stat
    import module_functions as fun 
     
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
    fig = plt.figure()
    ax = p3.Axes3D(fig)
         
    # Axes labels and title are established
    ax = fig.gca()
    ax.set_xlabel(info[0])
    ax.set_ylabel(info[1])
    ax.set_zlabel(info[2])
    ax.set_title(info[3])
    
    fig.gca().set_aspect('equal', adjustable='box')
    plt.grid()

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
    ax.scatter([0],[0],[0], color="k",s=100)
    ax.text(aux/100, aux/100, aux/100, info[5], [0,1,0])

    # Set up animation data
    n = len(position3D_1[0])
    
    if position3D_2 is not None:   
        data = [np.vstack((position3D_1[0], position3D_1[1], position3D_1[2])),
                np.vstack((position3D_2[0], position3D_2[1], position3D_2[2])),
                np.vstack((position3D_1[0], position3D_1[1], position3D_1[2])),
                np.vstack((position3D_2[0], position3D_2[1], position3D_2[2]))]
        lines = [ax.plot(data[0][0,0:1], data[0][1,0:1], data[0][2,0:1],
                                                      'r' , label = leg[0])[0],
                 ax.plot(data[0][0,0:1], data[0][1,0:1], data[0][2,0:1],
                                                      'b' , label = leg[1])[0]]
        pts   = [ax.plot(data[0][0,0:1],data[0][1,0:1],data[0][2,0:1],'ro')[0],
                 ax.plot(data[0][0,0:1],data[0][1,0:1],data[0][2,0:1],'bo')[0]]
        
    else:
        
        data = [np.vstack((position3D_1[0], position3D_1[1], position3D_1[2])),
                np.vstack((0, 0, 0)), 
                np.vstack((position3D_1[0], position3D_1[1], position3D_1[2])),
                np.vstack((0, 0, 0))]
        lines = [ax.plot(data[0][0,0:1], data[0][1,0:1], data[0][2,0:1],
                                                       'r', label = leg[0])[0],
                 ax.plot(data[0][0,0:1], data[0][1,0:1], data[0][2,0:1],
                                                       'b', label = leg[1])[0]]
        pts   = [ax.plot(data[0][0,0:1],data[0][1,0:1],data[0][2,0:1],'ro')[0],
                 ax.plot(data[0][0,0:1],data[0][1,0:1],data[0][2,0:1],'bo')[0]]
        
    time_text = ax.text2D(0.2, 0.2, "2D Text", transform=ax.transAxes)
    plt.legend(loc="center left")
    
    # Creating the Animation object
    ani = animation.FuncAnimation( fig, update_lines,n,fargs=(data,lines,pts),
                                                  interval=info[8], blit=False)
    
    if info[7] is 'store_yes':
        plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
        mywriter = animation.FFMpegWriter()
        ani.save('mymovie.mp4', writer=mywriter)
    
    return ani
    
    
    
    
    
    
    
    
    
def plot_planetmoon(time, planet, moon): 
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')    
    ax.set_aspect("equal")
    
    #draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
def plotvs_phasearea_cd(body):
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

    plt.figure()
    area = np.zeros(len(body.ephemeris.time))
    for i in range(len(body.ephemeris.time)):
        area[i] =  np.sum(body.grid.area[body.grid.shadow[i,:]==0])*180/m.pi*4
    
    # The orbits of the two bodies are plotted, if available
    plt.plot(body.ephemeris.time/3600/24, body.geometry.area_phase_c, 'r', linewidth=2,
                                                 linestyle="-", label="Continuous")

    #plt.plot(body.ephemeris.time/3600/24, body.geometry.area_phase_d,'b', linewidth=2, linestyle="-", label="discrete")
    plt.plot(body.ephemeris.time/3600/24, area,'g', linewidth=2, linestyle="-", label="discrete")
                  
    # Axes labels and title are established
    plt.xlabel('Time [days]')
    plt.ylabel('Area [-]')
    plt.title('Area vs Time')
    plt.legend()
    plt.grid()
    
    
