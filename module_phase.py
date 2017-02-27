# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 17:49:09 2016

@author: javier
"""
import numpy as np
import exopy_cfg as _cfg
from module_functions import PolyArea

def phase(body, star):
    '''
    body = Earth
    star = Sun
    '''
    
    print('\n    ... phase of ' + body.type + ' ' + body.name + '\n')
    
    approach = _cfg.approach
    compute  = _cfg.case
    N        = _cfg.N

    time = body.ephemeris.time
    T = len(time)
    t_count = np.array(range(len(time)))
            
    Rb  = body.properties.R
    Rs  = star.properties.R
    r_s = body.ephemeris.r_s
    
    alpha = body.geometry.phase_angle#alpha
    phi_s = body.geometry.solar_azimuth_angle
    
    if approach is 'parallel':
        xp = np.zeros_like(time)
        Rp = np.ones_like(time)
        angle = np.ones(T)*np.radians(90)
    elif approach is 'conical':
        xp = (Rs/Rb-1)*Rb/r_s
        Rp = (1-(xp)**2)**0.5   
        angle = np.arctan(Rp/xp)
            
    if compute in ['d','cd','dc']:
        
        Rob = body.ephemeris.position3D_s_ob  
        
        body.flag.phase_d = True
        angle_nodes = np.arccos(np.einsum('Ni,it->tN', body.grid.nodes_xyz_rot, -Rob)/(0.5*r_s[:,np.newaxis]))
        body.grid.illuminated_nodes[angle_nodes>angle[:,np.newaxis]] = 0       
        body.grid.shadow = body.grid.illuminated_nodes.astype(float)
        body.geometry.phase_area_d = np.degrees(np.sum(body.grid.area.reshape([1,len(body.grid.area)]).repeat(len(time),0)* ( ~body.grid.illuminated_nodes ), 1))*4

    if compute in ['c','cd','dc']:

        body.flag.phase_c = True    
            
        # Azimuthal angle for phase edge definition
        # Starting at theta = 0 allows a proper collocation of elements in the 
        # cartesian coordinate vectors
        theta = np.append(np.linspace(0,2*np.pi,N-1, endpoint = False), 0)
        
        Ralpha = np.array([[np.sin(alpha)       , np.zeros_like(alpha) , -np.cos(alpha)       ],
                           [np.zeros_like(alpha), np.ones_like(alpha)  ,  np.zeros_like(alpha)],
                           [np.cos(alpha)       , np.zeros_like(alpha) ,  np.sin(alpha)       ]])
        
        Rphi_s = np.array([[  np.cos(phi_s)       , -np.sin(phi_s)       , np.zeros_like(phi_s) ],
                           [  np.sin(phi_s)       ,  np.cos(phi_s)       , np.zeros_like(phi_s) ],  
                           [  np.zeros_like(phi_s),  np.zeros_like(phi_s), np.ones_like(phi_s) ]])  
        
        # Calculation of phase edge with reference x axis
        # Phase edge cartesian coordinates in solar directed frame. v [3,time,N]

        v = np.array([np.repeat(xp[:,np.newaxis],N,axis=1), Rp[:,np.newaxis]*np.sin(theta[:,np.newaxis].T), Rp[:,np.newaxis]*np.cos(theta[:,np.newaxis].T)])
        # Phase edge cartesian coordinates in intermediate ref. frame
        #    x, y, z = Ralpha * v
        xyz   = np.einsum('ijk,jkl->ikl',Ralpha,v)
        
        flag  = np.sum(xyz[2]>0,1)
        
        #if approach is 'radial_sphere':
        
        aux1 = np.diff(np.concatenate( ( np.zeros([T,1]) , xyz[2]<0 , np.zeros([T,1]) ), axis = 1 ),axis=1)
        aux2 = [np.where(aux1[i,:] != 0 )[0] for i in t_count]
        n_jumps = np.array([np.size(aux2[i]) for i in t_count])
        
        angles2 = [0]*T
        angles4 = [0]*T
        
        angles2 = [ np.linspace( np.arctan2(xyz[1,i,np.mod(aux2[i][0]-1,N)],xyz[0,i,np.mod(aux2[i][0]-1,N)]), np.arctan2(xyz[1,i,aux2[i][1]], xyz[0,i,aux2[i][1]]), np.diff(aux2[i]) , endpoint = True) if (flag[i]!=0)&((n_jumps[i]==2)|(n_jumps[i]==4)) else 0 for i in t_count]
        angles4 = [ np.linspace( np.arctan2(xyz[1,i,np.mod(aux2[i][2]-1,N)],xyz[0,i,np.mod(aux2[i][2]-1,N)]), np.arctan2(xyz[1,i,aux2[i][3]], xyz[0,i,aux2[i][3]]), np.diff(aux2[i]) , endpoint = True) if n_jumps[i]== 4 else 0 for i in t_count]
        
        x2 = [ np.cos(angles2[i]) if (flag[i]!=0)&((n_jumps[i]==2)|(n_jumps[i]==4)) else 0 for i in t_count]
        y2 = [ np.sin(angles2[i]) if (flag[i]!=0)&((n_jumps[i]==2)|(n_jumps[i]==4)) else 0 for i in t_count]       
        x4 = [ np.cos(angles4[i]) if (n_jumps[i]==4) else 0 for i in t_count]
        y4 = [ np.sin(angles4[i]) if (n_jumps[i]==4) else 0 for i in t_count]
          
#        plt.plot(xyz[0,1000,:],xyz[1,1000,:],'o-b')
#        plt.plot(xyz[0,1000,aux2[1000][0]],xyz[1,1000,aux2[1000][0]],'o-g',zorder = 5)
#        plt.plot(xyz[0,1000,aux2[1000][1]],xyz[1,1000,aux2[1000][1]],'o-m',zorder = 5)
#        print(xyz[2,1000,aux2[1000][1]])
#        print(xyz[2,1000,aux2[1000][0]])
        
        for i in t_count[flag==0]:                               xyz[:,i,:]                     = np.zeros([3,N])
        for i in t_count[(flag!=0)&((n_jumps==2)|(n_jumps==4))]: xyz[:,i,aux2[i][0]:aux2[i][1]] = (x2[i],y2[i],np.zeros_like(x2[i]))
        for i in t_count[n_jumps==4]:                            xyz[:,i,aux2[i][2]:aux2[i][3]] = (x4[i],y4[i],np.zeros_like(x4[i]))                           

        body.geometry.phase_points = np.einsum('ijk,jkl->ikl', Rphi_s, xyz)
        body.geometry.phase_area_c = np.degrees(PolyArea(xyz[0,:,:],xyz[1,:,:]))
        
#        plt.plot(xyz[0,1000,:],xyz[1,1000,:],'o-r')
              
    return body
