# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

# -*- coding: utf-8 -*-
# ==================================================================
# EXOPY module: exopy_geometry.py
# Delft University of Technology
# ------------------------------------------------------------------
# Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
# Date: 2016-2017
# ------------------------------------------------------------------
#
# Dependences:
#
# DESCRIPTION
# ------------------------------------------------------------------
# Script containing the functions required for the computation of
# the geometry parameters involved in the motion of the planetary
# system.
#
# LIST OF FUNCTIONS
# ------------------------------------------------------------------
#  - geometry: Calculates the geometry parameters.
#


# Required modules
#import mt_functions as fun, matplotlib.patches as patches
#import sys
#sys.path.append('/home/javier/anaconda3/lib/python3.5/site-packages')
#from exopy_functions import PolyArea
import numpy as np
import exopy_config as _cfg

def geometry(body, conf, ref_line_angle = None):
    '''
    ==================================================================
    EXOPY function: geometry()
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    INPUTS
    ------------------------------------------------------------------
    - body: planet/moon type of body object [-] ('body' object)
    - ref_line_angle: shift angle towards reference line for the com-
            putation of reflected radiance [rad] (numpy
            array)

    OUTPUTS
    ------------------------------------------------------------------
    - body: updated body object [-] ('body' object)


    DESCRIPTION
    ------------------------------------------------------------------
    Function computing the geometry parameters that are driven by the
    orbital motion of the extrasolar planetary system.


        '''


    # The current function accounts for the calculation of relevant variables
    # in the observer reference frame as well as computing the geometry angles

    print('\n    ... geometry of ' + body.type + ' ' + body.name + '\n')

    if body.type == 'star':
        body.ephemeris.position3D_s_ob = body.ephemeris.position3D_s
    else:

        el = conf.el
        az = conf.az

        time = body.ephemeris.time

        # Star-to-body vector: R
        R = body.ephemeris.position3D_s    # [m]


        # Conversion from orbital reference frame to observer reference frame
        # Rotation matrix Mrot:
        Mrot = np.array([[ np.cos(np.radians(el))*np.cos(np.radians(az)), np.cos(np.radians(el))*np.sin(np.radians(az)), -np.sin(np.radians(el)) ],
                         [-np.sin(np.radians(az))                       , np.cos(np.radians(az))                       ,            0            ],
                         [ np.sin(np.radians(el))*np.cos(np.radians(az)), np.sin(np.radians(el))*np.sin(np.radians(az)),  np.cos(np.radians(el)) ]])

        # Star-to-body vector in observer reference frame:
        Rob = np.einsum('ij,jk->ik',Mrot,R)

        body.geometry.ref_plane_angle = np.arctan2(Rob[1,:],Rob[0,:])

        key = 1
        if ref_line_angle is 'fix':
            key = 0
            body.geometry.ref_line_angle = np.arctan2(Rob[1,:],Rob[0,:])
            body.geometry.ref_plane_to_ref_line_angle = np.zeros_like(time)

        elif ref_line_angle is None:
            body.geometry.ref_line_angle = np.zeros_like(time)
            body.geometry.ref_plane_to_ref_line_angle = body.geometry.ref_plane_angle - body.geometry.ref_line_angle

        else:
            body.geometry.ref_line_angle = ref_line_angle
            body.geometry.ref_plane_to_ref_line_angle = body.geometry.ref_plane_angle - body.geometry.ref_line_angle


        Mrot = np.array([[ np.cos(body.geometry.ref_line_angle), np.sin(body.geometry.ref_line_angle), np.zeros_like(time)],
                         [-np.sin(body.geometry.ref_line_angle), np.cos(body.geometry.ref_line_angle), np.zeros_like(time)],
                         [ np.zeros_like(time)                 , np.zeros_like(time)                 , np.ones_like(time) ]])

        body.ephemeris.position3D_s_ob = np.einsum('ijt,jt->it',Mrot,Rob)
        # Conversion from orbital reference frame to observer reference frame
        #    # Rotation matrix Mrot:
        #    Mrot = np.array([[ np.cos(np.radians(el))*np.cos(np.radians(az)), np.cos(np.radians(el))*np.sin(np.radians(az)), -np.sin(np.radians(el)) ],
        #                     [-np.sin(np.radians(az))                       , np.cos(np.radians(az))                       ,            0            ],
        #                     [ np.sin(np.radians(el))*np.cos(np.radians(az)), np.sin(np.radians(el))*np.sin(np.radians(az)),  np.cos(np.radians(el)) ]])

        #    # Star-to-body vector in observer reference frame:
        #    Rob = np.einsum('ij,jk->ik',Mrot,R)

        #    phi_s = np.arctan2(Rob[1,:],Rob[0,:])
        #    body.geometry.ref_line_angle  = phi_s
        #    ref_line_angle = phi_s

        #elif ref_line_angle is None:

        #    ref_line_angle = np.zeros_like(time)

        #Mrot = np.array([[ np.cos(ref_line_angle)*np.cos(np.radians(el))*np.cos(np.radians(az))-np.sin(ref_line_angle)*np.sin(np.radians(az)), np.cos(ref_line_angle)*np.cos(np.radians(el))*np.sin(np.radians(az)) + np.sin(ref_line_angle)*np.cos(np.radians(az)), -np.cos(ref_line_angle)*np.sin(np.radians(el)) ],
        #                 [-np.sin(ref_line_angle)*np.cos(np.radians(el))*np.cos(np.radians(az))-np.cos(ref_line_angle)*np.sin(np.radians(az)),-np.sin(ref_line_angle)*np.cos(np.radians(el))*np.sin(np.radians(az)) + np.cos(ref_line_angle)*np.cos(np.radians(az)),  np.sin(ref_line_angle)*np.sin(np.radians(el)) ],
        #                 [ np.sin(np.radians(el))*np.cos(np.radians(az))*np.ones_like(time)                                                  , np.sin(np.radians(el))*np.sin(np.radians(az))*np.ones_like(time)                                                    ,  np.cos(np.radians(el))*np.ones_like(time)     ]])

        # Star-to-body vector in observer reference frame:
        #Rob = np.einsum('ijk,jk->ik',Mrot,R)
        #body.ephemeris.position3D_s_ob = Rob


        # Phase angle
        body.geometry.alpha = np.arccos(body.ephemeris.position3D_s_ob[2,:]/np.linalg.norm(body.ephemeris.position3D_s_ob, axis=0))
        body.geometry.phase_angle = np.arccos(-body.ephemeris.position3D_s_ob[2,:]/np.linalg.norm(body.ephemeris.position3D_s_ob, axis=0))

        # Azimuth angle
        #body.geometry.solar_azimuth_angle = np.mod(np.arctan2(Rob[1,:],Rob[0,:]),2*np.pi)*key

        body.geometry.solar_azimuth_angle = np.mod(np.arctan2(body.ephemeris.position3D_s_ob[1,:],body.ephemeris.position3D_s_ob[0,:]),2*np.pi)*key

        # Grid geometry
        if conf.case in ['d', 'cd', 'dc']:

            rotation = np.array([[-1, 0, 0], [ 0,-1, 0], [ 0, 0, 1]])

            body.grid.nodes_xyz_rot = np.einsum('Ni,ij->Nj',body.grid.nodes_xyz,rotation)

            body.grid.position_nodes_ob = body.ephemeris.position3D_s_ob.T[:,:,np.newaxis] + 2*body.properties.R*body.grid.nodes_xyz_rot.T[np.newaxis,:,:]
            #body.grid.position_nodes_reference_plane = body.ephemeris.position3D_reference_plane.T[:,:,np.newaxis] + 2*body.properties.R*body.grid.nodes_xyz_rot
            body.grid.distance_nodes_ob = np.linalg.norm(body.grid.position_nodes_ob,axis=1)
            #body.grid.position_nodes_ob = Rob[np.newaxis,:,:] + body.grid.nodes_xyz[:,:,np.newaxis]
            #body.grid.solar_zenith_angle = np.arccos(np.einsum('Ni,tiN->tN', 2*body.grid.nodes_xyz_rot, -body.grid.position_nodes_ob)/(np.linalg.norm(body.grid.position_nodes_ob,axis=1)))
            body.grid.solar_zenith_angle = np.arccos(np.einsum('Ni,tiN->tN', 2*body.grid.nodes_xyz_rot, -body.ephemeris.position3D_s_ob.T[:,:,np.newaxis])/(np.linalg.norm(body.ephemeris.position3D_s_ob.T[:,:,np.newaxis],axis=1)))
            body.grid.observer_zenith_angle = np.arccos(2*body.grid.nodes_xyz_rot[:,2])

            #angle_nodes = np.arccos(np.einsum('Ni,it->tN', body.grid.nodes_xyz_rot, Rob)/(0.5*Rnorm[:,np.newaxis]))
            body.grid.illuminated_nodes = np.ones([np.size(time),body.grid.N_points], dtype=bool)

            #body.grid.phase_angle = np.arccos(-body.grid.position_nodes_ob[:,2]/np.linalg.norm(body.grid.position_nodes_ob,axis=1))

            Mrot = np.array([[ np.cos(body.geometry.ref_plane_to_ref_line_angle), np.sin(body.geometry.ref_plane_to_ref_line_angle)],
                             [-np.sin(body.geometry.ref_plane_to_ref_line_angle), np.cos(body.geometry.ref_plane_to_ref_line_angle)]])

            xy = np.einsum('ijt,gj->itg',Mrot,body.grid.nodes_xyz_rot[:,0:2])
            x =  xy[0,:,:]#body.grid.nodes_xyz[:,0]*np.cos(body.geometry.ref_line_angle)+body.grid.nodes_xyz[:,1]*np.sin(body.geometry.ref_line_angle)
            y =  xy[1,:,:]#-body.grid.nodes_xyz[:,0]*np.sin(body.geometry.ref_line_angle)+body.grid.nodes_xyz[:,1]*np.cos(body.geometry.ref_line_angle)

            #body.grid.position_nodes_reference_plane = body.ephemeris.position3D_reference_plane.T[:,:,np.newaxis] + 2*body.properties.R*body.grid.nodes_xyz.T[np.newaxis,:,:]
            #body.grid.distance_nodes = np.linalg.norm(body.grid.position_nodes_reference_plane,axis=1)

            #body.grid.solar_zenith_angle = np.arccos(np.einsum('Ni,tiN->tN', 2*body.grid.nodes_xyz, -body.ephemeris.position3D_reference_plane.T[:,:,np.newaxis])/(np.linalg.norm(body.ephemeris.position3D_reference_plane.T[:,:,np.newaxis],axis=1)))
            #body.grid.observer_zenith_angle = np.arccos(2*body.grid.nodes_xyz[:,2])

            #body.grid.illuminated_nodes = np.ones([np.size(time),body.grid.N_points], dtype=bool)

            #x = body.grid.nodes_xyz[:,0]
            #y = body.grid.nodes_xyz[:,1]

            #===================
            # phi - phi_o angle
            #===================
            #term1 = np.cos(body.grid.phase_angle) - np.cos(body.grid.observer_zenith_angle[np.newaxis,:])*np.cos(body.grid.solar_zenith_angle)

            term1 = np.cos(body.geometry.phase_angle[:,np.newaxis]) - np.cos(body.grid.observer_zenith_angle[np.newaxis,:])*np.cos(body.grid.solar_zenith_angle)

            term2 = np.sin(body.grid.observer_zenith_angle[np.newaxis,:])*np.sin(body.grid.solar_zenith_angle)

            phi = np.zeros_like(body.grid.solar_zenith_angle)

            epsilon = 10E-9

            term3 = term1[np.abs(term2)>=epsilon]/term2[np.abs(term2)>=epsilon]
            term3[np.isnan(term3)] = 0
            #print(np.sum(term3>1))
            #print(np.sum(term3<-1))
            term3[term3> 1] =  1
            term3[term3<-1] = -1

            phi[np.abs(term2)>=epsilon] = -( np.pi - np.arccos(term3))

            phi[y<0] = -phi[y<0]

            body.grid.azimuth = phi

            #===================
            # beta angle
            #===================

            beta = np.zeros_like(body.grid.solar_zenith_angle)

            beta[(x!=0) & (y!=0)] = (np.arcsin(y[(x!=0) & (y!=0)]/np.sqrt(x[(x!=0) & (y!=0)]**2+y[(x!=0) & (y!=0)]**2)))

            beta[(y< 0) & (x>=0)] = beta[(y< 0) & (x>=0)] + np.pi
            beta[(y< 0) & (x< 0)] =-beta[(y< 0) & (x< 0)]
            beta[(y>=0) & (x< 0)] =-beta[(y>=0) & (x< 0)] + np.pi

            body.grid.beta = beta

    return body

