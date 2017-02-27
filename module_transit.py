# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 17:49:51 2016

@author: javier
"""


import exopy_cfg as _cfg
import numpy as np
from module_functions import PolyArea


def transits(bodies, ref_line_angle = None):   
    '''
    bodies = [Earth,Moon,Sun]
    compute = 'cd'
    N = 200
    el = 0
    az = 0
    ref_line_angle = None
    
    for body in bodies:
        if not hasattr(body.grid, 'shadow'):
           body.grid.shadow = np.zeros([len(body.ephemeris.time), len(body.grid.nodes)])  
    '''    
    import itertools
    
    # Time vector
    time = bodies[0].ephemeris.time
    T = len(time)
    t_count = np.array(range(len(time)))
    
    compute = _cfg.case
    N       = _cfg.N
    
    if ref_line_angle is None: ref_line_angle = np.zeros_like(time)
    
    # Algorithm control flags
    #merge_eclipses  = False

    # Create list for radii, re-initialize flags, and initialize points lists
    R = [bodies[i].properties.R for i in range(len(bodies))]
    for i in range(len(bodies)): bodies[i].flag.transit = []
    if compute in ['c','cd','dc']:
        for i in range(len(bodies)): 
            bodies[i].geometry.transit_points = []
            bodies[i].geometry.transit_area   = []
            bodies[i].flag.transit_c = True
    if compute in ['d','cd','dc']:
        for i in range(len(bodies)): 
            bodies[i].flag.transit_d = True
                           
    # Bodies are arranged according to their radii, so that that R1 <= R2 <= R3
    bodies = [i for (r,i) in sorted(zip(R,bodies))]
    R.sort()
    
#    # Conversion from orbital reference frame to observer reference frame
#    # Rotation matrix Mrot:
#    Mrot1 = np.array([[ np.cos(np.radians(el))*np.cos(np.radians(az)), np.cos(np.radians(el))*np.sin(np.radians(az)), -np.sin(np.radians(el)) ],
#                      [-np.sin(np.radians(az))                       , np.cos(np.radians(az))                       ,            0            ],
#                      [ np.sin(np.radians(el))*np.cos(np.radians(az)), np.sin(np.radians(el))*np.sin(np.radians(az)),  np.cos(np.radians(el)) ]])
#    
#    # HOW TO ACCOUNT FOR REFERENCE LINE? Rotation of ref line angle around z axis!
#    Mrot2 = np.array([[ np.cos(ref_line_angle)       , np.sin(ref_line_angle)       , np.zeros_like(ref_line_angle) ],
#                      [-np.sin(ref_line_angle)       , np.cos(ref_line_angle)       , np.zeros_like(ref_line_angle) ],
#                      [ np.zeros_like(ref_line_angle), np.zeros_like(ref_line_angle), np.ones_like(ref_line_angle)  ]])
#    
    # (n)
    # (2)  posibilities

    for i, j in itertools.combinations(range(len(bodies)),2):   

        bodies[i].flag.transit.append( [ np.zeros_like(time, dtype=bool), bodies[j].name ] )
        bodies[j].flag.transit.append( [ np.zeros_like(time, dtype=bool), bodies[i].name ] )
           
        print('\n    ... transits between ' + bodies[i].type + ' ' + bodies[i].name + ' and ' + bodies[j].type + ' ' + bodies[j].name + '\n')
                
        # Minimum distance for existance of transit
        L  = R[j] + R[i]
        M  = R[j] - R[i]

        # Ratio of radii
        r = R[i]/R[j]
            
        # A reference frame centered at bodyj is considered. The vector 
        # from body2 to bodyi is denominated 'd'.
        
        d = np.array([ bodies[i].ephemeris.position3D_s_ob[0,:] - 
                                       bodies[j].ephemeris.position3D_s_ob[0,:],
                       bodies[i].ephemeris.position3D_s_ob[1,:] - 
                                       bodies[j].ephemeris.position3D_s_ob[1,:],
                       bodies[i].ephemeris.position3D_s_ob[2,:] - 
                                       bodies[j].ephemeris.position3D_s_ob[2,:] ])    
        
        # Is bodyi in front of bodyj?
        bodyi_infront = np.array([False]).repeat(len(time))
        bodyi_infront[d[2,:]>0] = True
        
        # D: Distance from body2 to body1 in 2D plane
        D = np.linalg.norm([d[0,:], d[1,:]],axis = 0)     
        
        # Flags are assigned
        bodies[i].flag.transit[-1][0][np.invert(bodyi_infront)*(D<=L)] = True
        bodies[j].flag.transit[-1][0][(bodyi_infront)*(D<=L)] = True

        
        if compute in ['d', 'cd', 'dc']: 
#            print i,j

#            mesh_ri = np.ones_like(bodies[i].grid.shadow)*R[j]*1.5
#            mesh_rj = np.ones_like(bodies[j].grid.shadow)*R[i]*1.5
#            print(tt.time()-a)
#            mesh_rj[bodyi_infront & (D<=L),:] = np.linalg.norm( bodies[j].grid.nodes.T[:,np.newaxis,:]*2*R[j] - \
#                                d[0:2,bodyi_infront & (D<=L),np.newaxis], axis=0)
#            mesh_ri[np.invert(bodyi_infront)&(D<=L),:] = np.linalg.norm( bodies[i].grid.nodes.T[:,np.newaxis,:]*2*R[i] + \
#                                d[0:2,np.invert(bodyi_infront)&(D<=L),np.newaxis], axis=0)
#            bodies[j].grid.shadow[mesh_rj<R[i]] = 0
#            bodies[i].grid.shadow[mesh_ri<R[j]] = 0
            
#            print((np.linalg.norm( bodies[j].grid.nodes.T[:,np.newaxis,:]*2*R[j] - d[0:2,bodyi_infront & (D<=L),np.newaxis]         , axis=0)<R[i]).shape)            
#            print((np.linalg.norm( bodies[i].grid.nodes.T[:,np.newaxis,:]*2*R[i] + d[0:2,np.invert(bodyi_infront)&(D<=L),np.newaxis], axis=0)<R[j]).shape)
            b = (np.linalg.norm( -bodies[j].grid.nodes.T[:,np.newaxis,:]*2*R[j] - d[0:2,bodyi_infront & (D<=L),np.newaxis]         , axis=0)<R[i])
            a = (np.linalg.norm( -bodies[i].grid.nodes.T[:,np.newaxis,:]*2*R[i] + d[0:2,np.invert(bodyi_infront)&(D<=L),np.newaxis], axis=0)<R[j])
	
	    #print(b)		
	    #print(np.invert(b))
	    #print np.invert(b).astype(int).shape

	    bodies[j].grid.shadow[bodyi_infront * (D<=L), :] = bodies[j].grid.shadow[bodyi_infront * (D<=L), :] * np.invert(b).astype(int)
	    bodies[i].grid.shadow[np.invert(bodyi_infront) * (D<=L), :] = bodies[i].grid.shadow[np.invert(bodyi_infront) * (D<=L), :] * np.invert(a).astype(int)
	    #print b.shape

	    #print d.shape, bodies[j].grid.nodes.shape, bodies[j].grid.nodes.T[:,np.newaxis,:].shape, d[0:2,bodyi_infront * (D<=L),np.newaxis].shape 
	    #print b.shape
            #bodies[j].grid.shadow[np.linalg.norm( bodies[j].grid.nodes.T[:,np.newaxis,:]*2*R[j] - d[0:2,bodyi_infront * (D<=L),np.newaxis]         , axis=0)<R[i]] = 0
            #bodies[i].grid.shadow[np.linalg.norm( bodies[i].grid.nodes.T[:,np.newaxis,:]*2*R[i] + d[0:2,np.invert(bodyi_infront)*(D<=L),np.newaxis], axis=0)<R[j]] = 0
                        
        if compute in ['c','cd','dc']:
            
            bodies[i].geometry.transit_points.append( [ np.zeros([2,T,N]), bodies[j].name ] )
            bodies[i].geometry.transit_area.append(   [ np.zeros(T)      , bodies[j].name ] )            
            bodies[j].geometry.transit_points.append( [ np.zeros([2,T,N]), bodies[i].name ] )
            bodies[j].geometry.transit_area.append(   [ np.zeros(T)      , bodies[i].name ] )           
            
            if any(D<L):
            
                theta = np.linspace(0,360,N,endpoint=True)
                xy_circ1 = np.array([np.cos(np.radians(theta)),np.sin(np.radians(theta))])
                
                if R[j] == R[i]:
                    bodies[i].geometry.transit_points[-1][0][:,(D==0)&(~bodyi_infront),:] = xy_circ1[:,np.newaxis,:]
                    bodies[j].geometry.transit_points[-1][0][:,(D==0)&( bodyi_infront),:] = xy_circ1[:,np.newaxis,:]
                
                bodies[i].geometry.transit_points[-1][0][:,(D<=M)&(~bodyi_infront),:] = xy_circ1[:,np.newaxis,:]
                bodies[j].geometry.transit_points[-1][0][:,(D<=M)&( bodyi_infront),:] = d[0:2,( bodyi_infront)&(D<=M),np.newaxis]/R[j] + r*xy_circ1[:,np.newaxis,:]              
                
                bodies[i].geometry.transit_area[-1][0][(D<=M)&(~bodyi_infront)] = 180
                bodies[j].geometry.transit_area[-1][0][(D<=M)&( bodyi_infront)] = \
                PolyArea(bodies[j].geometry.transit_points[-1][0][0,(D<=M)&( bodyi_infront),:],bodies[j].geometry.transit_points[-1][0][1,(D<=M)&( bodyi_infront),:])
                
                if any(D>M):
                
                    n     = np.zeros([2,T])
                    geom  = np.zeros([3,T]) # Contains a,b,h
                    OMEGA = np.zeros([2,T]) # Contains OMEGA1, OMEGA2
                    PHI   = np.zeros(T)
                    MPHI  = np.zeros([2,2,T])
    
                    int_i  = (D<L)&(D>M)&(bodyi_infront)
                    int_j  = (D<L)&(D>M)&(~bodyi_infront)
                    int_ij = (D<L)&(D>M)
                    
                    geom[0:2,int_ij] = ((D[int_ij]**2+R[j]**2-R[i]**2)/(2*D[int_ij]), (D[int_ij]**2+R[i]**2-R[j]**2)/(2*D[int_ij]))
                    geom[2, int_ij]  = np.sqrt(R[j]**2-geom[0,int_ij]**2)
                    OMEGA[:, int_ij] = (np.mod(np.arctan2(geom[2,int_ij],geom[1,int_ij]),2*np.pi), np.mod(np.arctan2(geom[2,int_ij],geom[0,int_ij]),2*np.pi)) 
                    
                    n[0,int_ij]      = np.round(N*R[i]*OMEGA[0,int_ij]/(R[i]*OMEGA[0,int_ij]+R[j]*OMEGA[1,int_ij]))
                    n[1,int_ij]      = N - n[0,int_ij]
                
                    PHI[int_i]  = np.arctan2(d[1,int_i],d[0,int_i])
                    PHI[int_j]  = np.arctan2(d[1,int_j],d[0,int_j])+np.pi ####
                    MPHI[:,:,int_ij] = np.array([[ np.cos(PHI[int_ij]),  np.sin(PHI[int_ij])],
                                                 [ np.sin(PHI[int_ij]), -np.cos(PHI[int_ij])]])
                                            
                    omega = [0 if ~int_ij[t] else [np.linspace(-OMEGA[1,t],OMEGA[1,t],n[1,t]),np.linspace(-OMEGA[0,t],OMEGA[0,t],n[0,t]) ]     for t in t_count]

                    aux1 = np.array([np.concatenate(( np.einsum('ij,jk->ik',MPHI[:,:,t],np.array([np.cos(omega[t][0]), np.sin(omega[t][0])])), D[t]/R[j]*np.array([np.cos(PHI[t]),np.sin(PHI[t])])[:,np.newaxis]-r*np.einsum('ij,jk->ik',MPHI[:,:,t],np.array([np.cos(omega[t][1]),np.sin(omega[t][1])]))  ),axis=1)
                    for t in t_count[int_i]])
                        
                    aux2 = np.array([np.concatenate((np.einsum('ij,jk->ik',MPHI[:,:,t],np.array([np.cos(omega[t][1]), np.sin(omega[t][1])])), D[t]/R[i]*np.array([np.cos(PHI[t]),np.sin(PHI[t])])[:,np.newaxis]-(1/r)*np.einsum('ij,jk->ik',MPHI[:,:,t],np.array([np.cos(omega[t][0]),np.sin(omega[t][0])]))  ),axis=1)
                    for t in t_count[int_j]])
            
                    if len(aux1)!=0: 
                        bodies[j].geometry.transit_points[-1][0][:,int_i,:] = aux1.swapaxes(0,1)
                        bodies[j].geometry.transit_area[-1][0][int_i] = \
                            PolyArea(bodies[j].geometry.transit_points[-1][0][0,int_i,:],bodies[j].geometry.transit_points[-1][0][1,int_i,:])
                
                    if len(aux2)!=0: 
                        bodies[i].geometry.transit_points[-1][0][:,int_j,:] = aux2.swapaxes(0,1)
                        bodies[j].geometry.transit_area[-1][0][int_j] = \
                            PolyArea(bodies[j].geometry.transit_points[-1][0][0,int_j,:],bodies[j].geometry.transit_points[-1][0][1,int_j,:])
                
    return bodies
