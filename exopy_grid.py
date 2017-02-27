# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:41:28 2016

@author: javier
"""

from __future__ import division
from __future__ import absolute_import
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/javier/anaconda3/lib/python3.5/site-packages')
import math as m
import numpy as np
import numpy.linalg as la
from six.moves import range
#import meshpy.triangle as triangle




def square(Nsq):
    
    if Nsq == 0:
        sys.exit('Error 1 in grid creation: Nsq must be > 0!')
    
    a = 0.5 + np.sqrt(2)/(2*Nsq)
    
    mesh_nodes = np.zeros([Nsq*Nsq,2])
    mesh_faces = np.zeros([Nsq,Nsq,2,5])
    mesh_area  = np.zeros([Nsq,Nsq])
    mesh_nodes_xyz = np.zeros([Nsq,Nsq,3])
    
    x = np.linspace(-0.5+1/(2*Nsq),0.5-1/(2*Nsq),Nsq)
    y = np.linspace(-0.5+1/(2*Nsq),0.5-1/(2*Nsq),Nsq)
    xv, yv = np.meshgrid(x, y)
    mesh_nodes[:,0] = xv.flatten()
    mesh_nodes[:,1] = yv.flatten()
    
    r = np.zeros([Nsq,Nsq])
    r = np.zeros([Nsq,Nsq])
    
    for i in range(Nsq):
        for j in range(Nsq):
            r[i,j] = (x[i]**2 + y[j]**2)**0.5
            mesh_faces[i,j,0,:] = [-1/(2*Nsq), 1/(2*Nsq), 1/(2*Nsq), -1/(2*Nsq), -1/(2*Nsq)]+ x[i] 
            mesh_faces[i,j,1,:] = [-1/(2*Nsq), -1/(2*Nsq), 1/(2*Nsq), 1/(2*Nsq), -1/(2*Nsq)]+ y[j]  
            mesh_area[i,j]      = (x[1]-x[0])*(y[1]-y[0])
            theta = m.atan2(y[j],x[i])
            if r[i,j]<=0.5:
                mesh_nodes_xyz[i,j,:] = np.array([r[i,j]*np.cos(theta), r[i,j]*np.sin(theta), np.sqrt(0.5**2-r[i,j]**2)])
    mesh_nodes = np.array([yv[(r <= 0.5)], xv[(r <= 0.5)]]).T
    mesh_faces = mesh_faces[r <= 0.5, :,:]
    #mesh_area  = 4*mesh_area [r <= 0.5]
    mesh_area  = mesh_area[r <= 0.5]
    mesh_area  = np.pi/np.sum(r<=0.5) * np.ones_like(mesh_area)
    N_points   = len(mesh_nodes)
    Area       = np.sum(mesh_area)
    mesh_nodes_xyz = mesh_nodes_xyz[r<=0.5,:]
    
#    fig , ax = plt.subplots()
#    plt.plot(mesh_nodes[:,0], mesh_nodes[:,1], 'or')
#    for i in range(N_points):
#        plt.plot(mesh_faces[i,0,:], mesh_faces[i,1,:], 'b-', linewidth = '2')
#    circle1 = plt.Circle((0, 0), 0.5, color = 'k', fill=False, zorder=1)
#    ax.add_artist(circle1)   
#    ax.set_aspect('equal', adjustable='box')
#    ax.set_xlim([-0.55, 0.55])
#    ax.set_ylim([-0.55, 0.6])     
#    fig.suptitle('Regular square grid', fontsize=20)
#    ax.text(-0.4,  0.525,'Points: '+ str(N_points), fontsize= 15)
#    ax.text(0.14,  0.525,'$\^{A}$: '+ str(np.round(Area/(m.pi*0.5**2),3)), fontsize= 15)

    return mesh_nodes, mesh_faces, mesh_area, N_points, mesh_nodes_xyz



















def sphere(Nlon=7, Nlat=7, cos1 = 80, cos2 = 60, diff = 0.009):
    
    cos1 = cos1 * m.pi/180
    cos2 = cos2 * m.pi/180
    
    mesh_faces  = np.zeros([Nlon*Nlat,2,5])
    mesh_nodes  = np.zeros([Nlon*Nlat,2])
    mesh_area   = np.zeros([Nlon*Nlat])
    
    a = [-0.5 + 1/(2*Nlon) + i/Nlon for i in range(Nlon)]
    y = [-0.5 + 1/(2*Nlat) + i/Nlat for i in range(Nlat)]
    
    if cos1 != False:
        yy = [(0.5-diff)*m.sin(i*(2*cos1)/Nlat-cos1)/m.sin(cos1) for i in range(Nlat+1)]
    else: 
        yy = [-0.5+diff + (1-2*diff)*i/(Nlat) for i in range(Nlat+1)]
        
    if cos2 != False:
        aa = [0.5*m.sin(i*(2*cos2)/Nlon-cos2)/m.sin(cos2) for i in range(Nlon+1)]
    else: 
        aa = [-0.5 + i/Nlon for i in range(Nlon+1)]
        
    av, yv = np.meshgrid(a, y)
    xv = np.cos( np.arcsin(yv*2) )*av
        
    for i in range(Nlon):
        for j in range(Nlat):    
            mesh_faces[i*Nlat+j,1,:] = [ yy[j],yy[j],yy[j+1],yy[j+1],yy[j] ] 
            mesh_faces[i*Nlat+j,0,:] = [ m.cos(m.asin(2*mesh_faces[i*Nlat+j,1,0]))*aa[i], m.cos(m.asin(2*mesh_faces[i*Nlat+j,1,1]))*aa[i+1], m.cos(m.asin(2*mesh_faces[i*Nlat+j,1,2]))*aa[i+1], m.cos(m.asin(2*mesh_faces[i*Nlat+j,1,3]))*aa[i], m.cos(m.asin(2*mesh_faces[i*Nlat+j,1,0]))*aa[i] ]
            mesh_nodes[i*Nlat+j,:]   = [ np.mean(mesh_faces[i*Nlat+j,0,0:4]), np.mean(mesh_faces[i*Nlat+j,1,0:4]) ]
            mesh_area[i*Nlat+j]      = ( abs(mesh_faces[i*Nlat+j,0,1]-mesh_faces[i*Nlat+j,0,0])+abs(mesh_faces[i*Nlat+j,0,2]-mesh_faces[i*Nlat+j,0,3]) ) * 0.5 * abs( mesh_faces[i*Nlat+j,1,2] - mesh_faces[i*Nlat+j,1,1] )
            
    mesh_nodes_prime = np.array([xv.flatten(), yv.flatten()]).T
        
    N_points = Nlon*Nlat
    Area     = np.sum(mesh_area)
    
#    fig , ax = plt.subplots()
#    plt.plot(mesh_nodes[:,0], mesh_nodes[:,1], 'or')
#    for i in range(N_points):
#        plt.plot(mesh_faces[i,0,:], mesh_faces[i,1,:], 'b-', linewidth = '2')
#    circle1 = plt.Circle((0, 0), 0.5, color = 'k', fill=False, zorder=1)
#    ax.add_artist(circle1)   
#    ax.set_aspect('equal', adjustable='box')
#    ax.set_xlim([-0.55, 0.55])
#    ax.set_ylim([-0.55, 0.6])     
#    fig.suptitle('Spherical grid', fontsize=20)
#    ax.text(-0.4,  0.525,'Points: '+ str(N_points), fontsize= 15)
#    ax.text(0.14,  0.525,'$\^{A}$: '+ str(np.round(Area/(m.pi*0.5**2),3)), fontsize= 15)
    
    return mesh_nodes, mesh_faces, mesh_area, N_points



def radial(Nang=15, Nrad=3):
    
    mesh_faces = np.zeros([Nang*Nrad,2,5])
    mesh_nodes = np.zeros([Nang*Nrad,2])
    mesh_area  = np.zeros([Nang*Nrad])

    angle  = [2*m.pi/Nang * (i+0.5) for i in range(Nang)]
    radius = [1/(2*Nrad) * (0.5 + i) * m.cos(m.pi/Nang) for i in range(Nrad)]

    a = 0.5*m.cos(m.pi/Nang)/Nrad
    b = m.pi/Nang
    for i in range(Nang):
        for j in range(Nrad):
            mesh_nodes[i*Nrad+j,:]   = [radius[j]*m.cos(angle[i]), radius[j]*m.sin(angle[i])]   
            mesh_faces[i*Nrad+j,0,:] = [(1 + j)/(2*Nrad)*m.cos(angle[i]-b), (1 + j)/(2*Nrad)*m.cos(angle[i]+b), (j)/(2*Nrad)*m.cos(angle[i]+b), (j)/(2*Nrad)*m.cos(angle[i]-b), (1 + j)/(2*Nrad)*m.cos(angle[i]-b)] 
            mesh_faces[i*Nrad+j,1,:] = [(1 + j)/(2*Nrad)*m.sin(angle[i]-b), (1 + j)/(2*Nrad)*m.sin(angle[i]+b), (j)/(2*Nrad)*m.sin(angle[i]+b), (j)/(2*Nrad)*m.sin(angle[i]-b), (1 + j)/(2*Nrad)*m.sin(angle[i]-b)] 
            mesh_area[i*Nrad+j]      = m.sin(m.pi/Nang)*m.cos(m.pi/Nang)*( ((1 + j)/(2*Nrad))**2 - ((j)/(2*Nrad))**2 )
        
    N_points = Nang*Nrad
    Area     = np.sum(mesh_area)
    
#    fig , ax = plt.subplots()
#    plt.plot(mesh_nodes[:,0], mesh_nodes[:,1], 'or')
#    for i in range(N_points):
#        plt.plot(mesh_faces[i,0,:], mesh_faces[i,1,:], 'b-', linewidth = '2')
#    circle1 = plt.Circle((0, 0), 0.5, color = 'k', fill=False, zorder=1)
#    ax.add_artist(circle1)   
#    ax.set_aspect('equal', adjustable='box')
#    ax.set_xlim([-0.55, 0.55])
#    ax.set_ylim([-0.55, 0.6])     
#    fig.suptitle('Polar grid', fontsize=20)
#    ax.text(-0.4,  0.525,'Points: '+ str(N_points), fontsize= 15)
#    ax.text(0.14,  0.525,'$\^{A}$: '+ str(np.round(Area/(m.pi*0.5**2),3)), fontsize= 15)

    return mesh_nodes, mesh_faces, mesh_area, N_points


    

def triangl(max_area=0.03, min_angle=30, circle_edges = 20, centre = False):
    
    def round_trip_connect(start, end):
        return [(i, i+1) for i in range(start, end)] + [(end, start)]
        
    #a = 0.2
    if centre == True:
        points = [(0,0)] # [(a, 0), (a, a), (-a, a), (-a, -a), (a, -a), (a, 0)]
        facets = round_trip_connect(0, len(points)-1)
        markers = [2]#    points = []
    else:
        points  = []
        facets  = []
        markers = []
    circ_start = len(points)
    points.extend(
            (0.5 * np.cos(angle), 0.5 * np.sin(angle))
            for angle in np.linspace(0, 2*np.pi, circle_edges, endpoint=False))
    facets.extend(round_trip_connect(circ_start, len(points)-1))

    #markers = [2,2,2,2,2,2]

    markers.extend(list(np.ones(circle_edges, dtype='int')))
    markers = [int(i) for i in markers]

    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(facets, facet_markers=markers)
    #
    info.regions.resize(1)
    # points [x,y] in region, + region number, + regional area constraints
    info.regions[0] = ([0,0] + [1,0.05])

    mesh = triangle.build(info, volume_constraints=True, min_angle = min_angle, max_volume=max_area)

    mesh_points = np.array(mesh.points)
    mesh_tris = np.array(mesh.elements)
    mesh_attr = np.array(mesh.point_markers)
    
    mesh_nodes = np.array([((mesh_points[i[0],0]+mesh_points[i[1],0]+mesh_points[i[2],0])/3, (mesh_points[i[0],1]+mesh_points[i[1],1]+mesh_points[i[2],1])/3) for i in mesh_tris])
    mesh_area  = np.array([ 0.5 * abs( mesh_points[i[0],0] * (mesh_points[i[1],1]-mesh_points[i[2],1]) + mesh_points[i[1],0]*(mesh_points[i[2],1]-mesh_points[i[0],1]) +mesh_points[i[2],0]*(mesh_points[i[0],1]-mesh_points[i[1],1]) ) for i in mesh_tris ])
    mesh_faces = np.array([ [[ mesh_points[i[0],0], mesh_points[i[1],0], mesh_points[i[2],0], mesh_points[i[0],0] ], [ mesh_points[i[0],1], mesh_points[i[1],1], mesh_points[i[2],1], mesh_points[i[0],1] ] ] for i in mesh_tris ])
    
    N_points = len(mesh_tris)
    Area     = np.sum(mesh_area)
    
    mesh_nodes_xyz = np.array([mesh_nodes[:,0],mesh_nodes[:,1],np.sqrt(0.5**2-mesh_nodes[:,0]**2-mesh_nodes[:,1]**2)]).T

    
#    fig , ax = plt.subplots()
#    plt.plot(mesh_nodes[:,0], mesh_nodes[:,1], 'or')
#    plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, 'b-', linewidth = '2')
#    circle1 = plt.Circle((0, 0), 0.5, color = 'k', fill=False, zorder=1)
#    ax.add_artist(circle1)   
#    ax.set_aspect('equal', adjustable='box')
#    ax.set_xlim([-0.55, 0.55])
#    ax.set_ylim([-0.55, 0.6])     
#    fig.suptitle('Triangular grid', fontsize=20)
#    ax.text(-0.4,  0.525,'Points: '+ str(N_points), fontsize= 15)
#    ax.text(0.14,  0.525,'$\^{A}$: '+ str(np.round(Area/(m.pi*0.5**2),3)), fontsize= 15)
         
    return mesh_nodes, mesh_faces, mesh_area, N_points, mesh_nodes_xyz

     
     
#if __name__ == "__main__":
#    sph()
