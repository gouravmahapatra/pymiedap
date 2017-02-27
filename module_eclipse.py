# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 17:50:52 2016

@author: javier
"""
import numpy as np
import exopy_cfg as _cfg
from module_functions import PolyArea

def eclipse(bodies, star = None):
    '''
    bodies = [Moon,Earth]o
    star = Sunp

    
    for body in bodies:
        if not hasattr(body.grid, 'shadow'):
           body.grid.shadow = np.zeros([len(body.ephemeris.time), len(body.grid.nodes)])  
    '''    
    import time as tt
    a = tt.time()        
    
    import itertools
    
    # Time vector
    time = bodies[0].ephemeris.time
    T = len(time)
    t_count = np.array(range(len(time)))
    
    N = _cfg.N      
    compute = _cfg.case
     
    # Create list for radii, re-initialize flags, and initialize points lists
#    R = [bodies[i].properties.R for i in range(len(bodies))]
    for i in range(len(bodies)): 
        bodies[i].flag.umbra = []
        bodies[i].flag.antumbra = []
        bodies[i].flag.penumbra = []
    if compute in ['c','cd','dc']:
        for i in range(len(bodies)):
            bodies[i].geometry.umbra_points    = []
            bodies[i].geometry.antumbra_points = []
            bodies[i].geometry.penumbra_points = []
            bodies[i].geometry.eclipse_area    = []
            bodies[i].geometry.eclipse_area    = []
            bodies[i].geometry.eclipse_area    = []
            bodies[i].flag.eclipse_c = True
    if compute in ['d','cd','dc']:
        for i in range(len(bodies)):
            bodies[i].flag.eclipse_d = True
                           
#    # Bodies are arranged according to their radii, so that that R1 <= R2 <= R3
#    bodies = [i for (r,i) in sorted(zip(R,bodies))]
#    R.sort()

    # (n)
    # (2)  posibilities

    for i, j in itertools.permutations(range(len(bodies)),2):  

        bodies[i].flag.umbra.append(    [ np.zeros_like(time, dtype=bool), bodies[j].name ] )
        bodies[i].flag.antumbra.append( [ np.zeros_like(time, dtype=bool), bodies[j].name ] )
        bodies[i].flag.penumbra.append( [ np.zeros_like(time, dtype=bool), bodies[j].name ] )
           
        print('\n    ... eclipses between ' + bodies[i].type + ' ' + bodies[i].name + ' and ' + bodies[j].type + ' ' + bodies[j].name + '\n')
            
        # Which body is eclipsing and which one is eclipsed?
        # It is assumed that bodyi is being eclipsed
    
        Rs = star.properties.R
    
        Di = bodies[i].ephemeris.r_s
        Dj = bodies[j].ephemeris.r_s
    
        bodyi_behind = np.ones_like(time, dtype = bool)
        bodyi_behind[Dj>Di] = False
        
        di  = bodies[i].ephemeris.position3D_s_ob
        dj  = bodies[j].ephemeris.position3D_s_ob
        
        Radi = bodies[i].properties.R
        Radj = bodies[j].properties.R
                
        # IT IS IMPORTANT TO SEPARATE BOTH CASES, SINCE THE NUMBER OF NODES IN 
        # EACH MESH CAN BE DIFFERENT

        d_i = np.zeros([2,3,T])
        D_i = np.zeros([2,T])
        dijpp = np.zeros([3,T])

        d_i[:,:, bodyi_behind] = np.array([di[:, bodyi_behind], dj[:, bodyi_behind]])
        D_i[:, bodyi_behind] = np.array([Di[bodyi_behind], Dj[bodyi_behind]])
        
        
        # Changing reference frame
        Rphi_s    = np.zeros([3,3,bodyi_behind.size])        #                   d2
        Ralpha    = np.zeros([3,3,bodyi_behind.size])        #                   d2
        
        alpha = bodies[j].geometry.alpha
        phi_s = bodies[j].geometry.solar_azimuth_angle
        
        Rphi_s[:,:,bodyi_behind] = \
                 np.array([[-np.cos(phi_s[bodyi_behind])       , -np.sin(phi_s[bodyi_behind])       , np.zeros_like(phi_s[bodyi_behind]) ],
                           [ np.sin(phi_s[bodyi_behind])       , -np.cos(phi_s[bodyi_behind])       , np.zeros_like(phi_s[bodyi_behind]) ],
                           [ np.zeros_like(phi_s[bodyi_behind]),  np.zeros_like(phi_s[bodyi_behind]), np.ones_like(phi_s[bodyi_behind])  ]])
                          
        Ralpha[:,:,bodyi_behind] = \
                 np.array([[  np.sin(alpha[bodyi_behind])       , np.zeros_like(alpha[bodyi_behind]), np.cos(alpha[bodyi_behind])       ],
                           [  np.zeros_like(alpha[bodyi_behind]), np.ones_like(alpha[bodyi_behind]) , np.zeros_like(alpha[bodyi_behind])],
                           [ -np.cos(alpha[bodyi_behind])       , np.zeros_like(alpha[bodyi_behind]), np.sin(alpha[bodyi_behind])       ]])

        dijaux  = np.einsum('ijk,jk->ik', Rphi_s[:,:,bodyi_behind], -di[:,bodyi_behind]+dj[:,bodyi_behind])  
        dijp    = np.einsum('jik,jk->ik', Ralpha[:,:,bodyi_behind], dijaux)  
        
        gamma = np.arctan2(dijp[2,:], dijp[1,:])
                           
        Rgamma = np.array([[ np.ones_like(gamma) ,  np.zeros_like(gamma),  np.zeros_like(gamma) ],
                           [ np.zeros_like(gamma), -np.cos(gamma)       , -np.sin(gamma)        ],
                           [ np.zeros_like(gamma),  np.sin(gamma)       , -np.cos(gamma)        ]])  
                           
        dijpp[:,bodyi_behind] = np.einsum('ijk,jk->ik', Rgamma, dijp)          

                
        #########  CHANGE TO APPROACH GIVEN BY THE NEXT THREE LINES, I.E. USED DOT PRODUCT INSTEAD OF ATAN2 AND ROTATION MATRIXES

        d1_i       = np.zeros([T])          #                     Bi
        d2_i       = np.zeros([T])          #                     /\
        sinOMEGA_i = np.zeros([T])          #                    / \ 
        cosOMEGA_i = np.zeros([T])          #                   /  \  d1
        O2A6_i     = np.zeros([T])          #                B2/___\
        R_A6P_i    = np.zeros([2,T])        #                   d2
        sinrho_i   = np.zeros([T])       
        
#        d2_i[bodyi_behind] = np.einsum('iT,iT->T', di[:,bodyi_behind], dj[:,bodyi_behind])/Dj[bodyi_behind] - Dj[bodyi_behind]
#        d1_i[bodyi_behind] = np.sqrt(D_i[0, bodyi_behind]**2-(Dj[bodyi_behind]+d2_i[ bodyi_behind])**2)
        
        d2_i[bodyi_behind] = np.abs(dijpp[0,bodyi_behind])
        d1_i[bodyi_behind] = np.abs(dijpp[1,bodyi_behind])
        
        sinOMEGA_i[bodyi_behind]  =  (Rs + Radj) / Dj[bodyi_behind]
        cosOMEGA_i[bodyi_behind]  =   np.sqrt(1-(sinOMEGA_i[bodyi_behind])**2)
        
        O2A6_i[bodyi_behind]      =  (Radi + Radj) / sinOMEGA_i[bodyi_behind]
        R_A6P_i[:,bodyi_behind]   =   np.array([-d2_i[bodyi_behind] - O2A6_i[bodyi_behind] , - d1_i[bodyi_behind]])
        sinrho_i[bodyi_behind]    = - R_A6P_i[1,bodyi_behind] / ( np.linalg.norm(R_A6P_i[:,bodyi_behind], axis = 0) )
 
        eclipsed = (sinrho_i < sinOMEGA_i)
 
        if compute in ['d', 'cd', 'dc']: 
            
            if any(eclipsed): # Body i eclipsed at some time
                            
                # Update illuminated nodes indicator -> Non eclipse epochs are 
                # discarded
                i_illuminated = (bodies[i].grid.illuminated_nodes)&(eclipsed[:,np.newaxis])
                
                sinPSI = np.ones(T)*np.nan
                cosPSI = np.ones(T)
                tanPSI = np.ones(T)*np.nan

                sinPSI[eclipsed]   = (Rs - Radj) / Dj[eclipsed]
                cosPSI[eclipsed]   =  np.sqrt(Dj[eclipsed]**2-(Rs-Radj)**2)/Dj[eclipsed]
                tanPSI[eclipsed]   =  sinPSI[eclipsed]/cosPSI[eclipsed]
            
                ###############################################################
                ###############################################################
                ######            LOOK FOR FULL UMBRA NODES              ######
                ###############################################################
                ###############################################################

                O2A1   = np.zeros(T)
                R_A1P  = np.ones([2,T])
                coszeta= np.zeros(T)
                
                O2A1[eclipsed]    = -(Radj - Radi) / sinPSI[eclipsed]
                R_A1P[:,eclipsed] =   np.array([- d2_i[eclipsed] - O2A1[eclipsed], - d1_i[eclipsed] ])
                coszeta[eclipsed] =   R_A1P[0,eclipsed] / ( np.linalg.norm(R_A1P[:,eclipsed],axis=0) )

                # if coszeta > cosPSI --> Full umbra!
                full_umbra_times = eclipsed & (coszeta>cosPSI)
                bodies[i].grid.shadow[ full_umbra_times[:,np.newaxis] & (i_illuminated)  ] = 0
                                
                # full umbra flags are set                
                bodies[i].flag.umbra[-1][0][(full_umbra_times)&(eclipsed)] = True
                
                # Update illuminated nodes indicator
                i_illuminated = (i_illuminated)&(~full_umbra_times[:,np.newaxis])

                # Times for continuing the eclipse search:
                remaining = eclipsed&~full_umbra_times&bodyi_behind

                # The rest of the cases will require doing an node-wise analysis
                # So far d1 and d2 accounted for the center of the bodies, now, 
                # d1_nodes and d2_nodes do exactly the same but for each node

                # For the construction of d_nodes, it is taken into account that
                # now, di accounts for the position of the nodes and dj is still
                # the position of the center of the eclipsing body
                
                d_nod       = bodies[i].grid.position_nodes_ob.swapaxes(0,1)
#                d_nod_norm  = bodies[i].grid.distance_nodes_ob
                d2_nod      = np.zeros([T,bodies[i].grid.N_points])
                d1_nod      = np.zeros([T,bodies[i].grid.N_points])
#                d2_nod1      = np.zeros([T,bodies[i].grid.N_points])
#                d1_nod1      = np.zeros([T,bodies[i].grid.N_points])
#                aux         = np.zeros([T,bodies[i].grid.N_points])
                
                print len(remaining)
                print np.shape(Rphi_s)
                dijaux  = np.einsum('ijtl,jtl->itl', Rphi_s[:,:,remaining,np.newaxis], -d_nod[:,remaining,:] + dj[:,remaining,np.newaxis])  
                dijp    = np.einsum('jitl,jtl->itl', Ralpha[:,:,remaining,np.newaxis], dijaux)                  
                
                d2_nod[remaining,:] = np.abs(dijp[0,:,:])
                d1_nod[remaining,:] = np.linalg.norm(dijp[1:,:,:],axis=0)
                
                
#                aux[remaining,:]    = np.einsum('iTN,iT->TN', d_nod[:,remaining,:], dj[:,remaining])#/Dj[:,np.newaxis] - Dj[:,np.newaxis]
#                d2_nod1[remaining,:] = aux[remaining,:]/Dj[remaining,np.newaxis] - Dj[remaining,np.newaxis]
#                del aux
#                d1_nod1[remaining,:] = np.sqrt(d_nod_norm[remaining,:]**2-(Dj[remaining,np.newaxis]+d2_nod[remaining,:])**2)   
##                
                # FROM d1_nod AND d2_nod, ANY ANGLE CAN BE CALCULATED FOR EACH POINT IN THE MESH
                # LET'S CONTINUE LOOKING FOR OTHER ECLIPSE CASES
        
                ###############################################################
                ###############################################################
                ######           LOOK FOR PARTIAL UMBRA NODES            ######
                ###############################################################
                ###############################################################

                O2A2   = np.zeros(T)
                O2A3   = np.zeros(T)
                R_A2P  = np.zeros([2,T])
                R_A3P  = np.zeros([2,T])
                costheta = np.zeros(T)
              
                O2A2[remaining]    = -(Radi + Radj) / sinPSI[remaining]
                O2A3[remaining]    = - Radj / sinPSI[remaining]
                
                R_A2P[:,remaining] =   np.array([ - d2_i[remaining] - O2A2[remaining], d1_i[remaining] ])
                R_A3P[:,remaining] =   np.array([ - d2_i[remaining] - O2A3[remaining], d1_i[remaining] ])
                
                R_A2P_norm = np.linalg.norm( R_A2P, axis=0 )
                R_A3P_norm = np.linalg.norm( R_A3P, axis=0 )
                costheta[remaining]=   R_A2P[0,remaining]/R_A2P_norm[remaining]
                      
                # if costheta > cosPSI and (RA2P>R1/tanPSI or RA3P<R1)  --> Partial umbra!
                
                partial_umbra_times = (remaining) & (costheta>cosPSI) & ( (R_A2P_norm > Radi/tanPSI) | (R_A3P_norm < Radi) )
                partial_umbra_nodes = (partial_umbra_times[:,np.newaxis]) & (i_illuminated)                       
                
                # Now we know which nodes at which times are suitable for being in umbra.
                # A final decision needs calculating the beta angle:
                
                R_A3P_nod  = np.zeros([2,T,bodies[i].grid.N_points])#*np.nan   # I commented the nan product, since it multiplies the time by ~o(2)
                R_A3P_nod_norm = np.zeros([T,bodies[i].grid.N_points])
                O2A3       = O2A3[:,np.newaxis].repeat(bodies[i].grid.N_points, axis= 1)
                sinbeta    = np.zeros([T,bodies[i].grid.N_points])
                
                R_A3P_nod[0,partial_umbra_nodes] = - d2_nod[partial_umbra_nodes] - O2A3[partial_umbra_nodes]  #O2A3[partial_umbra,np.newaxis]
                R_A3P_nod[1,partial_umbra_nodes] =   d1_nod[partial_umbra_nodes]########## he quitado un signo negativo -
                R_A3P_nod_norm[partial_umbra_nodes] = np.linalg.norm(R_A3P_nod[:,partial_umbra_nodes],axis=0) #np.sqrt( d1_nod[partial_umbra_nodes]**2  +  (d2_nod[partial_umbra_nodes] + O2A3[partial_umbra_nodes])**2 )
                
                sinbeta[partial_umbra_nodes]  = R_A3P_nod[1,partial_umbra_nodes] / R_A3P_nod_norm[partial_umbra_nodes]
            
                partial_umbra = (sinbeta < sinPSI[:,np.newaxis]) & (partial_umbra_nodes)
                
                bodies[i].grid.shadow[ partial_umbra ] = 0
                
                # partial umbra flags are set                
                bodies[i].flag.umbra[-1][0][np.sum(partial_umbra,axis=1,dtype=bool)] = True   
                
                # Update illuminated nodes indicator
                i_illuminated = i_illuminated&~partial_umbra
                                            
                ###############################################################
                ###############################################################
                ######             LOOK FOR ANTUMBRA NODES               ######
                ###############################################################
                ###############################################################
                                            
                # First, the center of the body needs to fulfill:
                # 180-zeta  <  PSI
                # coszeta was already calulated at the beginning

                # if coszeta < -cosPSI --> Antumbra!
                antumbra_times = remaining & (coszeta<-cosPSI)
                
                # Which nodes are to be evaluated?
                antumbra_nodes = antumbra_times[:,np.newaxis] & i_illuminated 
                
                # Now we know which nodes at which times are suitable for being in antumbra.
                # A final decision needs calculating the beta angle:
                cosbeta = np.zeros([T,bodies[i].grid.N_points])

                R_A3P_nod[0,antumbra_nodes] = - d2_nod[antumbra_nodes] - O2A3[antumbra_nodes]  #O2A3[partial_umbra,np.newaxis]
                R_A3P_nod[1,antumbra_nodes] = - d1_nod[antumbra_nodes]
                R_A3P_nod_norm[antumbra_nodes] = np.sqrt( d1_nod[antumbra_nodes]**2  +  (d2_nod[antumbra_nodes] + O2A3[antumbra_nodes])**2 )
                
                cosbeta[antumbra_nodes]  = R_A3P_nod[0,antumbra_nodes] / R_A3P_nod_norm[antumbra_nodes] 
            
                antumbra = (cosbeta < -cosPSI[:,np.newaxis]) & (antumbra_nodes)
                
                angle_sun = np.arcsin( Rs   / bodies[i].grid.distance_nodes_ob[antumbra] )
                angle_j   = np.arcsin( Radj / np.sqrt(d1_nod[antumbra]**2 + d2_nod[antumbra]**2) )
                
                bodies[i].grid.shadow[ antumbra ] = np.ones_like(angle_j) - (angle_j/angle_sun)**2
                
                # Penumbra flags are set                
                bodies[i].flag.antumbra[-1][0][np.sum(antumbra,axis=1,dtype=bool)] = True    
                
                # Update illuminated nodes indicator
                i_illuminated = i_illuminated&~antumbra                
                             
                ###############################################################
                ###############################################################
                ######             LOOK FOR PENUMBRA NODES               ######
                ###############################################################
                ###############################################################

                # The remaining points with i_illuminated == True are candidates
                # for the penumbra region
               
                O2A4 = np.zeros([T])
                O2A4[remaining]  =  Radj / sinOMEGA_i[remaining]

                R_A4P_nod  = np.zeros([2,T,bodies[i].grid.N_points])#*np.nan   # I commented the nan product, since it multiplies the time by ~o(2)
                R_A4P_nod_norm = np.zeros([T,bodies[i].grid.N_points])
                O2A4 = O2A4[:,np.newaxis].repeat(bodies[i].grid.N_points, axis= 1)
                cosomega = np.zeros([T,bodies[i].grid.N_points])
                
                penumbra_nodes = i_illuminated # Remaining illuminated nodes --> Non-eclipse nodes were set to False at the beginning               
                
                R_A4P_nod[0,penumbra_nodes] = - d2_nod[penumbra_nodes] - O2A4[penumbra_nodes]  #O2A3[partial_umbra,np.newaxis]
                R_A4P_nod[1,penumbra_nodes] = - d1_nod[penumbra_nodes]
                R_A4P_nod_norm[penumbra_nodes] = np.sqrt( d1_nod[penumbra_nodes]**2  +  (d2_nod[penumbra_nodes] + O2A4[penumbra_nodes])**2 )
                
                cosomega[penumbra_nodes]  = - R_A4P_nod[0,penumbra_nodes] / R_A4P_nod_norm[penumbra_nodes] 
            
                penumbra = (cosomega > cosOMEGA_i[:,np.newaxis]) & (penumbra_nodes)

                angle_sun = np.arcsin( Rs   / bodies[i].grid.distance_nodes_ob[penumbra] )
                angle_j   = np.arcsin( Radj / np.sqrt(d1_nod[penumbra]**2 + d2_nod[penumbra]**2) )

                dj = dj[:,:,np.newaxis].repeat(bodies[i].grid.N_points, axis= 2)

                aux = np.zeros([T,bodies[i].grid.N_points])
                delta_angle = np.zeros([T,bodies[i].grid.N_points])
                theta1 = np.zeros([T,bodies[i].grid.N_points])
                theta2 = np.zeros([T,bodies[i].grid.N_points])
                aux[remaining,:] = np.einsum('itN,itN->tN', -d_nod[:,remaining,:], -d_nod[:,remaining,:] + dj[:,remaining,:])
                aux[penumbra] = aux[penumbra]/bodies[i].grid.distance_nodes_ob[penumbra]/np.linalg.norm(-d_nod[:,penumbra] + dj[:,penumbra],axis=0)
                delta_angle[penumbra] = (np.arccos(aux[penumbra]))

		cord = 1/delta_angle[penumbra]*np.sqrt(4*delta_angle[penumbra]**2*angle_sun**2 - (delta_angle[penumbra]**2+angle_sun**2-angle_j**2)**2)
		xc_j     = (delta_angle[penumbra]**2+angle_j**2-angle_sun**2)/(2*delta_angle[penumbra])
		xc_sun   = (delta_angle[penumbra]**2+angle_sun**2-angle_j**2)/(2*delta_angle[penumbra])

		#bodies[i].cord = cord
		#bodies[i].miralo = (delta_angle[penumbra]**2+angle_j**2-angle_sun**2)
		#bodies[i].angle_j = angle_j
		#bodies[i].angle_sun = angle_sun
		#bodies[i].delta_angle = delta_angle[penumbra]

                theta1 = 2*np.arcsin( cord/2/angle_j  )
		theta1p = np.copy(theta1)
		theta1p[(angle_sun>angle_j)&(xc_sun>delta_angle[penumbra])] = 2*np.pi - theta1[(angle_sun>angle_j)&(xc_sun>delta_angle[penumbra])]
                theta2 = 2*np.arcsin( cord/2/angle_sun)
		theta2p = np.copy(theta2)
		theta2p[(angle_j>angle_sun)&(xc_j  >delta_angle[penumbra])] = 2*np.pi - theta2[(angle_j>angle_sun)&(xc_j  >delta_angle[penumbra])]

		sign_j = np.ones_like(theta1)
		sign_sun = np.ones_like(theta1)
		sign_j[(angle_sun>angle_j)&(xc_sun>delta_angle[penumbra])] = -1
		sign_sun[(angle_j>angle_sun)&(xc_j>delta_angle[penumbra])] = -1

                #theta1 = np.arccos( (angle_j**2/(2*delta_angle[penumbra]*angle_sun))*(-1+(angle_sun/angle_j)**2+(delta_angle[penumbra]/angle_j)**2) )
                #theta2 = np.arcsin(angle_sun/angle_j*np.sin(theta1))

                #bodies[i].residual = delta_angle[penumbra]**2 - angle_j**2 - angle_sun**2

		#bodies[i].xc_j = xc_j 
		#bodies[i].xc_sun = xc_sun
		#bodies[i].sign_j = sign_j
		#bodies[i].sign_sun = sign_sun
	
		#bodies[i].theta1 = theta1
		#bodies[i].theta2 = theta2
		

		A2 = (theta2p-sign_sun*np.sin(theta2))/2/np.pi
		A1 = ((angle_j/angle_sun))**2*(theta1p-sign_j*np.sin(theta1))/2/np.pi
		R = (angle_j/angle_sun)**2

		#bodies[i].A1 = A1
		#bodies[i].A2 = A2
		#bodies[i].R  = R

                aux1 = A1 + A2
                #aux1 = (theta2-np.sin(theta2))/np.pi+R*(theta1-np.sin(theta1))/np.pi
                bodies[i].grid.shadow[penumbra] = (np.ones_like(aux1) - aux1)
                
                # Penumbra flags are s
                bodies[i].flag.penumbra[-1][0][np.sum(penumbra,axis=1,dtype=bool)] = True              
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        if compute in ['c', 'cd', 'dc']: 
            
            if any(eclipsed): # Body i eclipsed at some time
                 
                bodies[i].geometry.umbra_points.append(    [ np.zeros([2,T,N]), bodies[j].name ] )
                bodies[i].geometry.antumbra_points.append( [ np.zeros([2,T,N]), bodies[j].name ] )
                bodies[i].geometry.penumbra_points.append( [ np.zeros([2,T,N]), bodies[j].name ] )
                bodies[i].geometry.eclipse_area.append(    [ np.zeros(T)      , bodies[j].name ] ) 
                bodies[i].geometry.eclipse_area.append(    [ np.zeros(T)      , bodies[j].name ] ) 
                bodies[i].geometry.eclipse_area.append(    [ np.zeros(T)      , bodies[j].name ] )    
                
#                alpha  = np.zeros([T])
#                phi_s  = np.zeros([T])
#                gamma  = np.zeros([T])
#                Rphi_s = np.zeros([3,3,T])
#                Ralpha = np.zeros([3,3,T])
                
#                di  = np.zeros([3,T])
#                dj  = np.zeros([3,T])
#                dij = np.zeros([3,T])
#                dijp= np.zeros([3,T])
                
                
                # Changing reference frame
                
                
##                 Rotation of reference frame matrices
#                alpha[eclipsed] = bodies[j].geometry.alpha[eclipsed]
#                phi_s[eclipsed] = bodies[j].geometry.solar_azimuth_angle[eclipsed]
#                
#                Rphi_s[:,:,eclipsed] = np.array([[-np.cos(phi_s[eclipsed])       , -np.sin(phi_s[eclipsed])       , np.zeros_like(phi_s[eclipsed]) ],
#                                                 [ np.sin(phi_s[eclipsed])       , -np.cos(phi_s[eclipsed])       , np.zeros_like(phi_s[eclipsed]) ],
#                                                 [ np.zeros_like(phi_s[eclipsed]),  np.zeros_like(phi_s[eclipsed]), np.ones_like(phi_s[eclipsed])  ]])
#                                  
#                Ralpha[:,:,eclipsed] = np.array([[  np.sin(alpha[eclipsed])       , np.zeros_like(alpha[eclipsed]), np.cos(alpha[eclipsed])       ],
#                                                 [  np.zeros_like(alpha[eclipsed]), np.ones_like(alpha[eclipsed]) , np.zeros_like(alpha[eclipsed])],
#                                                 [ -np.cos(alpha[eclipsed])       , np.zeros_like(alpha[eclipsed]), np.sin(alpha[eclipsed])       ]])
#        
#                di[:,eclipsed]  = bodies[i].ephemeris.position3D_s_ob[:,eclipsed]
#                dj[:,eclipsed]  = bodies[j].ephemeris.position3D_s_ob[:,eclipsed]
#                        
#                dij[:,eclipsed] = dj[:,eclipsed] - di[:,eclipsed] # From i to j
##                dji = di[eclipsed] - dj[eclipsed]        
#        
#                dijp[:,eclipsed]  = np.einsum('ijk,jk->ik', Rphi_s[:,:,eclipsed], dij[:,eclipsed])  
#                dijp[:,eclipsed]    = np.einsum('jik,jk->ik', Ralpha[:,:,eclipsed], dijp[:,eclipsed])  
#                
#                gamma = np.arctan2(dijp[2,:], dijp[1,:])
                                   
#                Rgamma = np.array([[ np.ones_like(gamma) ,  np.zeros_like(gamma),  np.zeros_like(gamma) ],
#                                   [ np.zeros_like(gamma), -np.cos(gamma)       , -np.sin(gamma)        ],
#                                   [ np.zeros_like(gamma),  np.sin(gamma)       , -np.cos(gamma)        ]])                  
                
#                Rt = np.zeros([3,3,T])
                
#                Rt[:,:,eclipsed] = np.array([[ - np.sin(alpha[eclipsed])*np.cos(phi_s[eclipsed])  , - np.cos(alpha[eclipsed])*np.sin(gamma[eclipsed])*np.cos(phi_s[eclipsed]) - np.cos(gamma[eclipsed])*np.sin(phi_s[eclipsed]) ,   np.cos(alpha[eclipsed])*np.cos(gamma[eclipsed])*np.cos(phi_s[eclipsed]) + np.sin(gamma[eclipsed])*np.sin(phi_s[eclipsed])  ],
#                                             [ - np.sin(alpha[eclipsed])*np.sin(phi_s[eclipsed])  , - np.cos(alpha[eclipsed])*np.sin(gamma[eclipsed])*np.sin(phi_s[eclipsed]) + np.cos(gamma[eclipsed])*np.cos(phi_s[eclipsed]) ,   np.cos(alpha[eclipsed])*np.cos(gamma[eclipsed])*np.sin(phi_s[eclipsed]) - np.sin(gamma[eclipsed])*np.cos(phi_s[eclipsed])  ],
#                                             [   np.cos(alpha[eclipsed])                          , - np.sin(gamma[eclipsed])*np.sin(alpha[eclipsed])                                                                           , - np.sin(alpha[eclipsed])*np.cos(gamma[eclipsed])                                                                            ]])
                R_1 = np.zeros([3,3,T])
                
                R_1[:,:,eclipsed] = np.array([[ np.sin(alpha[eclipsed])       ,   np.cos(alpha[eclipsed])*np.sin(gamma[eclipsed]), - np.cos(alpha[eclipsed])*np.cos(gamma[eclipsed]) ],
                                              [ np.zeros_like(alpha[eclipsed]), - np.cos(gamma[eclipsed])                        ,   np.sin(gamma[eclipsed])                         ],
                                              [ np.cos(alpha[eclipsed])       , - np.sin(gamma[eclipsed])*np.sin(alpha[eclipsed]), - np.sin(alpha[eclipsed])*np.cos(gamma[eclipsed]) ]])
                                
                Rphi_s = np.zeros([2,2,T])
                Rphi_s[:,:,eclipsed] = np.array([[-np.cos(phi_s[eclipsed]),  np.sin(phi_s[eclipsed]) ],
                                                 [-np.sin(phi_s[eclipsed]), -np.cos(phi_s[eclipsed]) ]])
                                
                # Unit circle is created here to be assigned in the steps ahead
                theta = np.append(np.linspace(0,2*np.pi,N-1, endpoint = False), 0)
                circle = np.array([np.cos(theta), np.sin(theta)])
#                theta = np.linspace(0,2*np.pi,N-1, endpoint = False)
                
                
                #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
                ######            LOOK FOR FULL UMBRA EPOCHS             ######
                ###############################################################
                
                if compute == 'c': # Otherwise, we take the values from the discrete calculation
                                    
                    # Update illuminated nodes indicator -> Non eclipse epochs are 
                    # discarded
#                    i_illuminated = (bodies[i].grid.illuminated_nodes)&(eclipsed[:,np.newaxis])
                    
                    sinPSI = np.ones(T)*np.nan
                    cosPSI = np.ones(T)
                    tanPSI = np.ones(T)*np.nan
    
                    sinPSI[eclipsed]   = (Rs - Radj) / Dj[eclipsed]
                    cosPSI[eclipsed]   =  np.sqrt(Dj[eclipsed]**2-(Rs-Radj)**2)/Dj[eclipsed]
                    tanPSI[eclipsed]   =  sinPSI[eclipsed]/cosPSI[eclipsed]
                
                    O2A1   = np.zeros(T)
                    R_A1P  = np.ones([2,T])
                    coszeta= np.zeros(T)
                    
                    O2A1[eclipsed]    = -(Radj - Radi) / sinPSI[eclipsed]
                    R_A1P[:,eclipsed] =   np.array([- d2_i[eclipsed] - O2A1[eclipsed], - d1_i[eclipsed] ])
                    coszeta[eclipsed] =   R_A1P[0,eclipsed] / ( np.linalg.norm(R_A1P[:,eclipsed],axis=0) )
    
                    # if coszeta > cosPSI --> Full umbra!
                    full_umbra_times = eclipsed & (coszeta>cosPSI)                
                
                bodies[i].geometry.umbra_points[-1][0][:,full_umbra_times,:] = circle[:,np.newaxis,:]

                # Times for continuing the eclipse search: The full umbra epochs are, 
                # together with the full antumbra epochs, the only ones that can
                # be exluded from the rest of the analysis, the partial cases cannot, since
                # they can host penumbra cases too! 
                # However, no type of antumbra can occur at the same time than 
                # umbra! This means that we can save time by defining a time 
                # array for the calculation of penumbra and a different one for
                # the upcoming antumbra, where the partial umbra cases are to 
                # be deleted.
                remaining = eclipsed&~full_umbra_times 
                
                
                #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
                ######          LOOK FOR PARTIAL UMBRA EPOCHS            ######
                ###############################################################
                
                if compute == 'c': # Otherwise, we take the values from the discrete calculation
                
                    O2A2   = np.zeros(T)
                    O2A3   = np.zeros(T)
                    R_A2P  = np.zeros([2,T])
                    R_A3P  = np.zeros([2,T])
                    costheta = np.zeros(T)
                  
                    O2A2[remaining]    = -(Radi + Radj) / sinPSI[remaining]
                    O2A3[remaining]    = - Radj / sinPSI[remaining]
                    
                    R_A2P[:,remaining] =   np.array([ - d2_i[remaining] - O2A2[remaining], -d1_i[remaining] ])
                    R_A3P[:,remaining] =   np.array([ - d2_i[remaining] - O2A3[remaining], -d1_i[remaining] ])
                    
                    R_A2P_norm = np.linalg.norm( R_A2P, axis=0 )
                    R_A3P_norm = np.linalg.norm( R_A3P, axis=0 )
                    costheta[remaining]=   R_A2P[0,remaining]/R_A2P_norm[remaining]
                    
                    # if costheta > cosPSI and (RA2P>R1/tanPSI or RA3P<R1)  --> Partial umbra!
                    partial_umbra_times = (remaining) & (costheta>cosPSI) & ( (R_A2P_norm > Radi/tanPSI) | (R_A3P_norm < Radi) )
                    
                # UMBRA SPOT
                umbra_spot = partial_umbra_times &  (coszeta < -cosPSI)
                
                
                
                if np.sum(umbra_spot)>0:
                    ################################
                    a = R_A3P[1,umbra_spot]
                    b = -R_A3P[0,umbra_spot]
                    c = sinPSI[umbra_spot]    
                    d = cosPSI[umbra_spot]
                    cc= tanPSI[umbra_spot]
                    # If we enter this function, it is 100% sure that we have intersection
                    # In order to avoid the singularity at a = 0...
                    #a[a==0] = 1e-30
        
                    A  = 1/(2*a*d**2)
                    x0 = b*c**2
                    y0 = ( (c*b)**2 - (Radi**2 - a**2) ) / (2*a)                
    
                    x = (+2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta)[np.newaxis,:]/A[:,np.newaxis] + ( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta)[np.newaxis,:]/A[:,np.newaxis])**2 - 4*(x0[:,np.newaxis]**2 + y0[:,np.newaxis]/A[:,np.newaxis] + cc[:,np.newaxis]*b[:,np.newaxis]*np.cos(theta)[np.newaxis,:]/A[:,np.newaxis] ) )**0.5) / (2*Radi)
                    y = (cc[:,np.newaxis] * (x*Radi - b[:,np.newaxis]) * np.cos(theta)[np.newaxis,:] - a[:,np.newaxis]) / Radi
                    z = (cc[:,np.newaxis] * (x*Radi - b[:,np.newaxis]) * np.sin(theta)[np.newaxis,:]    ) / Radi
                    
                    xyz   = np.einsum('ijk,jkl->ikl',R_1[:,:,umbra_spot],np.array([x,y,z]))
                    
    #                flag  = np.sum(xyz[2]>0,1)
                    t = np.size(xyz,1)
                    
                    aux1 = np.diff((xyz[2]<0).astype(int),axis=1)
                    aux2 = [np.where(aux1[ii,:] != 0 )[0] for ii in range(t)]
    #                ref_line = np.array([np.arctan2(xyz[1,:,0],xyz[0,:,0]), xyz[0,:,0],xyz[1,:,0]])
                    #exterior= np.array([np.array([aux1[ii,aux2[ii][0]], -aux1[ii,aux2[ii][-1]]]) for ii in range(t)])
    #                n_patches = np.array([2 if (np.size(aux2[ii])==3 or (np.size(aux2[ii])==2 and exterior[ii,0]==-1)  or (np.size(aux2[ii])==4 and exterior[ii,0]==1) ) else \
    #                                      3 if (np.size(aux2[ii])==4 and exterior[ii,0]==-1) else  \
    #                                      1 for ii in range(t)])
    
                    coord = np.zeros([2,t,N])
                                          
                    coord[:,xyz[2,:,:]>0] = xyz[:2,xyz[2,:,:]>0]
                    
    #                coord[0,:,:] = np.array([  np.concatenate((xyz[0,ii,0:aux2[ii][0]+1], np.cos( ref_line[0,ii] + np.linspace(np.mod(np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1])-ref_line[0,ii],2*np.pi),np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]],xyz[0,ii,aux2[ii][1]])-ref_line[0,ii],2*np.pi),aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:])) if (np.size(aux2[ii])==2 and exterior[ii,0]==1) else \
    #                                           np.zeros([N])   for ii in range(t)  ])
                                               
                    coord[0,:,:] = -np.array([  xyz[0,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[0,ii,aux2[ii][3]+1:]))                                                                                                                                      if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],2*np.pi),xyz[0,ii,aux2[ii][3]+1]) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])
                    coord[1,:,:] = np.array([  xyz[1,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi), N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[1,ii,aux2[ii][3]+1:]))                                                                                                                                                   if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],xyz[0,ii,aux2[ii][3]+1]),2*np.pi) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])                                            
                                               
    #                coord[1,:,:] = np.array([  np.concatenate((xyz[1,ii,0:aux2[ii][0]+1], np.sin( ref_line[0,ii] + np.linspace(np.mod(np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1])-ref_line[0,ii],2*np.pi),np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]],xyz[0,ii,aux2[ii][1]])-ref_line[0,ii],2*np.pi),aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:])) if (np.size(aux2[ii])==2 and exterior[ii,0]==1) else \
    #                                           np.zeros([N])   for ii in range(t)  ])                                      
                                          
                    coord   = np.einsum('jik,jkl->ikl',Rphi_s[:,:,umbra_spot],coord)
                                          
                                          
#                    fig = plt.figure()
#                    ax = fig.add_subplot(111, projection='3d')            
#                    plt.plot(x[20,:],y[20,:],z[20,:],'o-')
#                    plt.plot(xyz[0,20,:],xyz[1,20,:],xyz[2,20,:],'o-')
#                    ax.set_xlim([-1, 1])
#                    ax.set_ylim([-1, 1])
#                    ax.set_zlim([-1, 1])   
#                    ax.set_xlabel('x')   
#                    ax.set_ylabel('y')   
#                    ax.set_zlabel('z')                     
#                    circle1 = plt.Circle((0,0), 1, color = 'k', fill=False, zorder=1)
#                    ax.add_artist(circle1)                                          
#                                          
#                                          
                    bodies[i].geometry.umbra_points[-1][0][:,umbra_spot,:] = coord
#                    
#                    fig , ax1 = plt.subplots()
#                    ax1.set_ylim([-1,1])
#                    ax1.set_xlim([-1,1])
#                    ax1.set_aspect('equal', adjustable='box')
#                    circle1 = plt.Circle((0,0), 1, color = 'k', fill=False, zorder=1)
#                    ax1.add_artist(circle1)
#                    plt.plot(coord[0,20,:],coord[1,20,:],'ko-')
#                    plt.plot(xyz[0,20,:],xyz[1,20,:],'bo-')
#                    plt.plot(x[20,:],y[20,:],'ro-')

                

                # UMBRA EDGE
                umbra_edge = partial_umbra_times & ~(coszeta < -cosPSI)
                
                if np.sum(umbra_edge)>0:
                    ################################
                    a = R_A3P[1,umbra_edge]
                    b = -R_A3P[0,umbra_edge]
                    c = sinPSI[umbra_edge]    
                    d = cosPSI[umbra_edge]
                    cc= tanPSI[umbra_edge]
                    # If we enter this function, it is 100% sure that we have intersection
                    # In order to avoid the singularity at a = 0...
                    #a[a==0] = 1e-30
                    
                    t = np.size(cc)
                    
                    A  = 1/(2*a*d**2)
                    x0 = b*c**2
                    y0 = ( (c*b)**2 - (Radi**2 - a**2) ) / (2*a)
                    
                    theta0 = np.arccos( 2*A/c * (b-x0 + ((x0-b)**2 + y0/A )**0.5 ) )
                
#                    theta1 = np.zeros([round(N/2),t])
#                    theta2 = np.zeros([N-round(N/2)+2,t])
                    
                    i1 = np.arange(round(N/2))
                    i2 = np.arange(N-round(N/2)+2)
                    
                    theta1 = -theta0[:,np.newaxis]*np.cos(i1*np.pi/(round(N/2)-1))[np.newaxis,:]
                    theta2 =  theta0[:,np.newaxis]*np.cos(i2*np.pi/(N-round(N/2)+1))[np.newaxis,:]
                    
#                    for i in range(len(theta)):
#                        theta[i] = -theta0 * m.cos(i*m.pi/(round(N/2)-1))
#                        #theta[i] = -theta0 * m.cos(i*m.pi/(N-1))
#                    for i in range(len(theta1)):
#                        theta1[i] = theta0 * m.cos(i*m.pi/(len(theta1)-1))
#                        #theta1[i] = theta0 * m.cos(i*m.pi/(N-1))
                    theta2 = np.delete(theta2,np.size(theta2,1)-1,1)
                    theta2 = np.delete(theta2,0,1)  
                    
                    aux_1 = ( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta1)/A[:,np.newaxis])**2 - 4*(x0[:,np.newaxis]**2 + y0[:,np.newaxis]/A[:,np.newaxis] + cc[:,np.newaxis]*b[:,np.newaxis]*np.cos(theta1)/A[:,np.newaxis] ))
                    aux_1[aux_1<0]=0
                    aux_2 = ( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta2)/A[:,np.newaxis])**2 - 4*(x0[:,np.newaxis]**2 + y0[:,np.newaxis]/A[:,np.newaxis] + cc[:,np.newaxis]*b[:,np.newaxis]*np.cos(theta2)/A[:,np.newaxis] ))
                    aux_2[aux_2<0]=0
                    
                    x = np.concatenate(( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta1)/A[:,np.newaxis] + ( aux_1 )**0.5) / (2*Radi), (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta2)/A[:,np.newaxis] - ( aux_2 )**0.5) / (2*Radi)  ),1)
                    y = np.concatenate(( (cc[:,np.newaxis] * np.abs(x[:,0:round(N/2)]*Radi - b[:,np.newaxis]) * np.cos(theta1) - a[:,np.newaxis]) / Radi , ( cc[:,np.newaxis] * np.abs(x[:,round(N/2):]*Radi - b[:,np.newaxis]) * np.cos(theta2) - a[:,np.newaxis])/Radi ),1)
                    z = np.concatenate(( (cc[:,np.newaxis] * (x[:,0:round(N/2)]*Radi - b[:,np.newaxis]) * np.sin(theta1)       ) / Radi , ( cc[:,np.newaxis] * (x[:,round(N/2):]*Radi - b[:,np.newaxis]) * np.sin(theta2)   )/Radi  ),1)
                    
                    
#                    fig = plt.figure()
#                    ax = fig.add_subplot(111, projection='3d')            
#                    plt.plot(x[0,:],y[0,:],z[0,:],'o-')
#                    ax.set_xlim([-1, 1])
#                    ax.set_ylim([-1, 1])
#                    ax.set_zlim([-1, 1])   
#                    ax.set_xlabel('x')   
#                    ax.set_ylabel('y')   
#                    ax.set_zlabel('z')                     
#                    circle1 = plt.Circle((0,0), 1, color = 'k', fill=False, zorder=1)
#                    ax.add_artist(circle1)   
                    
                    xyz   = np.einsum('ijk,jkl->ikl',R_1[:,:,umbra_edge],np.array([x,y,z]))
                    
    #                flag  = np.sum(xyz[2]>0,1)
                    
                    aux1 = np.diff((xyz[2]<0).astype(int),axis=1)
                    aux2 = [np.where(aux1[ii,:] != 0 )[0] for ii in range(t)]
    #                ref_line = np.array([np.arctan2(xyz[1,:,0],xyz[0,:,0]), xyz[0,:,0],xyz[1,:,0]])
                    #exterior= np.array([np.array([aux1[ii,aux2[ii][0]], -aux1[ii,aux2[ii][-1]]]) for ii in range(t)])
    #                n_patches = np.array([2 if (np.size(aux2[ii])==3 or (np.size(aux2[ii])==2 and exterior[ii,0]==-1)  or (np.size(aux2[ii])==4 and exterior[ii,0]==1) ) else \
    #                                      3 if (np.size(aux2[ii])==4 and exterior[ii,0]==-1) else  \
    #                                      1 for ii in range(t)])
    
                    coord = np.zeros([2,t,N])
                                                              
    #                coord[0,:,:] = np.array([  np.concatenate((xyz[0,ii,0:aux2[ii][0]+1], np.cos( ref_line[0,ii] + np.linspace(np.mod(np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1])-ref_line[0,ii],2*np.pi),np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]],xyz[0,ii,aux2[ii][1]])-ref_line[0,ii],2*np.pi),aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:])) if (np.size(aux2[ii])==2 and exterior[ii,0]==1) else \
    #                                           np.zeros([N])   for ii in range(t)  ])
                                               
                    coord[0,:,:] = -np.array([  xyz[0,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[0,ii,aux2[ii][3]+1:]))                                                                                                                                      if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],2*np.pi),xyz[0,ii,aux2[ii][3]+1]) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])
                    coord[1,:,:] = np.array([  xyz[1,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi), N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[1,ii,aux2[ii][3]+1:]))                                                                                                                                                   if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],xyz[0,ii,aux2[ii][3]+1]),2*np.pi) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])                                            
                                          
                    coord   = np.einsum('jik,jkl->ikl',Rphi_s[:,:,umbra_edge],coord)                            
                               
                    bodies[i].geometry.umbra_points[-1][0][:,umbra_edge,:] = coord                        
                    
            
                # Times for continuing the eclipse search:
                remaining_antumbra = remaining &~ partial_umbra_times
                
                #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
                ######          LOOK FOR FULL ANTUMBRA EPOCHS            ######
                ###############################################################
                
                if compute == 'c': # Otherwise, we take the values from the discrete calculation
                                    
                    O2A2   = np.zeros(T)
                    R_A2P  = np.zeros([2,T])
                    costheta = np.zeros(T)
                  
                    O2A2[remaining_antumbra]    = -(Radi + Radj) / sinPSI[remaining_antumbra]
                    
                    R_A2P[:,remaining_antumbra] =   np.array([ - d2_i[remaining_antumbra] - O2A2[remaining_antumbra], -d1_i[remaining_antumbra] ])
                    R_A2P_norm = np.linalg.norm( R_A2P, axis=0 )
     
                    costheta[remaining_antumbra]=   R_A2P[0,remaining_antumbra]/R_A2P_norm[remaining_antumbra]
    
                # if costheta < - cosPSI --> Full antumbra!
                full_antumbra_times = remaining_antumbra & (costheta<-cosPSI)
                    
                bodies[i].geometry.antumbra_points[-1][0][:,full_antumbra_times,:] = circle[:,np.newaxis,:]

                # Times for continuing the eclipse search:
                remaining = remaining &~ full_antumbra_times
                remaining_antumbra = remaining_antumbra &~ full_antumbra_times
                
                #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
                ######         LOOK FOR PARTIAL ANTUMBRA EPOCHS          ######
                ###############################################################
                
                partial_antumbra_times = remaining_antumbra & (coszeta < -cosPSI)
                
                # ANTUMBRA SPOT
                antumbra_spot = partial_antumbra_times & (costheta>cosPSI)
                
                
                if np.sum(antumbra_spot)>0:
                    ################################
                    a = R_A3P[1,antumbra_spot]
                    b = R_A3P[0,antumbra_spot]
                    c = sinPSI[antumbra_spot]    
                    d = cosPSI[antumbra_spot]
                    cc= tanPSI[antumbra_spot]
                    # If we enter this function, it is 100% sure that we have intersection
                    # In order to avoid the singularity at a = 0...
                    #a[a==0] = 1e-30


                    A  = 1/(2*a*d**2)
                    x0 = b*c**2
                    y0 = ( (c*b)**2 - (Radi**2 - a**2) ) / (2*a)                
    
                    x = (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta)[np.newaxis,:]/A[:,np.newaxis] - ( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta)[np.newaxis,:]/A[:,np.newaxis])**2 - 4*(x0[:,np.newaxis]**2 + y0[:,np.newaxis]/A[:,np.newaxis] + cc[:,np.newaxis]*b[:,np.newaxis]*np.cos(theta)[np.newaxis,:]/A[:,np.newaxis] ) )**0.5) / (2*Radi)
                    y = (cc[:,np.newaxis] * (x*Radi - b[:,np.newaxis]) * np.cos(theta)[np.newaxis,:] - a[:,np.newaxis]) / Radi
                    z = -(cc[:,np.newaxis] * (x*Radi - b[:,np.newaxis]) * np.sin(theta)[np.newaxis,:]    ) / Radi
                    
                    xyz   = np.einsum('ijk,jkl->ikl',R_1[:,:,antumbra_spot],np.array([x,y,z]))
                    
    #                flag  = np.sum(xyz[2]>0,1)
                    t = np.size(xyz,1)
                    
                    aux1 = np.diff((xyz[2]<0).astype(int),axis=1)
                    aux2 = [np.where(aux1[ii,:] != 0 )[0] for ii in range(t)]
    #                ref_line = np.array([np.arctan2(xyz[1,:,0],xyz[0,:,0]), xyz[0,:,0],xyz[1,:,0]])
                    #exterior= np.array([np.array([aux1[ii,aux2[ii][0]], -aux1[ii,aux2[ii][-1]]]) for ii in range(t)])
    #                n_patches = np.array([2 if (np.size(aux2[ii])==3 or (np.size(aux2[ii])==2 and exterior[ii,0]==-1)  or (np.size(aux2[ii])==4 and exterior[ii,0]==1) ) else \
    #                                      3 if (np.size(aux2[ii])==4 and exterior[ii,0]==-1) else  \
    #                                      1 for ii in range(t)])
    
                    coord = np.zeros([2,t,N])
                                          
#                    coord[:,xyz[2,:,:]>0] = xyz[:2,xyz[2,:,:]>0]
                    
    #                coord[0,:,:] = np.array([  np.concatenate((xyz[0,ii,0:aux2[ii][0]+1], np.cos( ref_line[0,ii] + np.linspace(np.mod(np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1])-ref_line[0,ii],2*np.pi),np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]],xyz[0,ii,aux2[ii][1]])-ref_line[0,ii],2*np.pi),aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:])) if (np.size(aux2[ii])==2 and exterior[ii,0]==1) else \
    #                                           np.zeros([N])   for ii in range(t)  ])
                                               
                    coord[0,:,:] = -np.array([  xyz[0,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[0,ii,aux2[ii][3]+1:]))                                                                                                                                      if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],2*np.pi),xyz[0,ii,aux2[ii][3]+1]) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])
                    coord[1,:,:] = np.array([  xyz[1,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi), N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[1,ii,aux2[ii][3]+1:]))                                                                                                                                                   if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],xyz[0,ii,aux2[ii][3]+1]),2*np.pi) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])                                          
                                               
    #                coord[1,:,:] = np.array([  np.concatenate((xyz[1,ii,0:aux2[ii][0]+1], np.sin( ref_line[0,ii] + np.linspace(np.mod(np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1])-ref_line[0,ii],2*np.pi),np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]],xyz[0,ii,aux2[ii][1]])-ref_line[0,ii],2*np.pi),aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:])) if (np.size(aux2[ii])==2 and exterior[ii,0]==1) else \
    #                                           np.zeros([N])   for ii in range(t)  ])
                                          
                    coord   = np.einsum('jik,jkl->ikl',Rphi_s[:,:,antumbra_spot],coord)
                                          
                                          
#                    fig = plt.figure()
#                    ax = fig.add_subplot(111, projection='3d')            
#                    plt.plot(x[20,:],y[20,:],z[20,:],'o-')
#                    plt.plot(xyz[0,20,:],xyz[1,20,:],xyz[2,20,:],'o-')
#                    ax.set_xlim([-1, 1])
#                    ax.set_ylim([-1, 1])
#                    ax.set_zlim([-1, 1])   
#                    ax.set_xlabel('x')   
#                    ax.set_ylabel('y')   
#                    ax.set_zlabel('z')                     
#                    circle1 = plt.Circle((0,0), 1, color = 'k', fill=False, zorder=1)
#                    ax.add_artist(circle1)                                          

                    bodies[i].geometry.antumbra_points[-1][0][:,antumbra_spot,:] = coord
                    
#                    fig , ax1 = plt.subplots()
#                    ax1.set_ylim([-1,1])
#                    ax1.set_xlim([-1,1])
#                    ax1.set_aspect('equal', adjustable='box')
#                    circle1 = plt.Circle((0,0), 1, color = 'k', fill=False, zorder=1)
#                    ax1.add_artist(circle1)
#                    plt.plot(coord[0,20,:],coord[1,20,:],'ko-')
#                    plt.plot(xyz[0,20,:],xyz[1,20,:],'bo-')
#                    plt.plot(x[20,:],y[20,:],'ro-')

                
                # ANTUMBRA EDGE
                antumbra_edge = partial_antumbra_times & ~(costheta>cosPSI)         
                
                if np.sum(antumbra_edge)>0:
                    ################################
                    a = R_A3P[1,antumbra_edge]
                    b = R_A3P[0,antumbra_edge]
                    c = sinPSI[antumbra_edge]    
                    d = cosPSI[antumbra_edge]
                    cc= tanPSI[antumbra_edge]
                    
                    # If we enter this function, it is 100% sure that we have intersection
                    # In order to avoid the singularity at a = 0...
                    #a[a==0] = 1e-30
                    
                    t = np.size(cc)
                    
                    A  = 1/(2*a*d**2)
                    x0 = b*c**2
                    y0 = ( (c*b)**2 - (Radi**2 - a**2) ) / (2*a)
                    
                    theta0 = np.arccos( 2*A/c * (b-x0 + ((x0-b)**2 + y0/A )**0.5 ) )
                
#                    theta1 = np.zeros([round(N/2),t])
#                    theta2 = np.zeros([N-round(N/2)+2,t])
                    
                    i1 = np.arange(round(N/2))
                    i2 = np.arange(N-round(N/2)+2)
                    
                    theta1 = -theta0[:,np.newaxis]*np.cos(i1*np.pi/(round(N/2)-1))[np.newaxis,:]
                    theta2 =  theta0[:,np.newaxis]*np.cos(i2*np.pi/(N-round(N/2)+1))[np.newaxis,:]
                    
#                    for i in range(len(theta)):
#                        theta[i] = -theta0 * m.cos(i*m.pi/(round(N/2)-1))
#                        #theta[i] = -theta0 * m.cos(i*m.pi/(N-1))
#                    for i in range(len(theta1)):
#                        theta1[i] = theta0 * m.cos(i*m.pi/(len(theta1)-1))
#                        #theta1[i] = theta0 * m.cos(i*m.pi/(N-1))
                    theta2 = np.delete(theta2,np.size(theta2,1)-1,1)
                    theta2 = np.delete(theta2,0,1)  
                    
                    aux_1 = ( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta1)/A[:,np.newaxis])**2 - 4*(x0[:,np.newaxis]**2 + y0[:,np.newaxis]/A[:,np.newaxis] + cc[:,np.newaxis]*b[:,np.newaxis]*np.cos(theta1)/A[:,np.newaxis] ))
                    aux_1[aux_1<0]=0
                    aux_2 = ( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta2)/A[:,np.newaxis])**2 - 4*(x0[:,np.newaxis]**2 + y0[:,np.newaxis]/A[:,np.newaxis] + cc[:,np.newaxis]*b[:,np.newaxis]*np.cos(theta2)/A[:,np.newaxis] ))
                    aux_2[aux_2<0]=0
                    
                    x = np.concatenate(( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta1)/A[:,np.newaxis] + ( aux_1 )**0.5) / (2*Radi), (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta2)/A[:,np.newaxis] - ( aux_2 )**0.5) / (2*Radi)  ),1)
                    y = np.concatenate(( (cc[:,np.newaxis] * np.abs(x[:,0:round(N/2)]*Radi - b[:,np.newaxis]) * np.cos(theta1) - a[:,np.newaxis]) / Radi , ( cc[:,np.newaxis] * np.abs(x[:,round(N/2):]*Radi - b[:,np.newaxis]) * np.cos(theta2) - a[:,np.newaxis])/Radi ),1)
                    z = -np.concatenate(( (cc[:,np.newaxis] * (x[:,0:round(N/2)]*Radi - b[:,np.newaxis]) * np.sin(theta1)       ) / Radi , ( cc[:,np.newaxis] * (x[:,round(N/2):]*Radi - b[:,np.newaxis]) * np.sin(theta2)   )/Radi  ),1)
                    
                    
#                    fig = plt.figure()
#                    ax = fig.add_subplot(111, projection='3d')            
#                    plt.plot(x[0,:],y[0,:],z[0,:],'o-')
#                    ax.set_xlim([-1, 1])
#                    ax.set_ylim([-1, 1])
#                    ax.set_zlim([-1, 1])   
#                    ax.set_xlabel('x')   
#                    ax.set_ylabel('y')   
#                    ax.set_zlabel('z')                     
#                    circle1 = plt.Circle((0,0), 1, color = 'k', fill=False, zorder=1)
#                    ax.add_artist(circle1)   
                    
                    xyz   = np.einsum('ijk,jkl->ikl',R_1[:,:,antumbra_edge],np.array([x,y,z]))
                    
    #                flag  = np.sum(xyz[2]>0,1)
                    
                    aux1 = np.diff((xyz[2]<0).astype(int),axis=1)
                    aux2 = [np.where(aux1[ii,:] != 0 )[0] for ii in range(t)]
    #                ref_line = np.array([np.arctan2(xyz[1,:,0],xyz[0,:,0]), xyz[0,:,0],xyz[1,:,0]])
                    #exterior= np.array([np.array([aux1[ii,aux2[ii][0]], -aux1[ii,aux2[ii][-1]]]) for ii in range(t)])
    #                n_patches = np.array([2 if (np.size(aux2[ii])==3 or (np.size(aux2[ii])==2 and exterior[ii,0]==-1)  or (np.size(aux2[ii])==4 and exterior[ii,0]==1) ) else \
    #                                      3 if (np.size(aux2[ii])==4 and exterior[ii,0]==-1) else  \
    #                                      1 for ii in range(t)])
    
                    coord = np.zeros([2,t,N])
                                                              
    #                coord[0,:,:] = np.array([  np.concatenate((xyz[0,ii,0:aux2[ii][0]+1], np.cos( ref_line[0,ii] + np.linspace(np.mod(np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1])-ref_line[0,ii],2*np.pi),np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]],xyz[0,ii,aux2[ii][1]])-ref_line[0,ii],2*np.pi),aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:])) if (np.size(aux2[ii])==2 and exterior[ii,0]==1) else \
    #                                           np.zeros([N])   for ii in range(t)  ])
                                               
                    coord[0,:,:] = -np.array([  xyz[0,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[0,ii,aux2[ii][3]+1:]))                                                                                                                                      if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],2*np.pi),xyz[0,ii,aux2[ii][3]+1]) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])
                    coord[1,:,:] = np.array([  xyz[1,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi), N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[1,ii,aux2[ii][3]+1:]))                                                                                                                                                   if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],xyz[0,ii,aux2[ii][3]+1]),2*np.pi) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])                                             
                                          
                    coord   = np.einsum('jik,jkl->ikl',Rphi_s[:,:,antumbra_edge],coord)                            
                               
                    bodies[i].geometry.antumbra_points[-1][0][:,antumbra_edge,:] = coord     


                #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
                ######          LOOK FOR FULL PENUMBRA EPOCHS            ######
                ###############################################################
                
                                    
                O2A5   = np.zeros(T)
                R_A5P  = np.zeros([2,T])
                cosphi = np.zeros(T)
              
                O2A5[remaining]    = (Radj - Radi) / sinOMEGA_i[remaining]
                R_A5P[:,remaining] =   np.array([ - d2_i[remaining] - O2A5[remaining], -d1_i[remaining] ])
                R_A5P_norm = np.linalg.norm( R_A5P, axis=0 )
                
                cosphi[remaining]=   -R_A5P[0,remaining]/R_A5P_norm[remaining]
    
                # if costheta < - cosPSI --> Full penumbra!
                full_penumbra_times = remaining & (cosphi > cosOMEGA_i)
                    
                bodies[i].geometry.penumbra_points[-1][0][:,full_penumbra_times,:] = circle[:,np.newaxis,:]

                # Times for continuing the eclipse search:
                remaining = remaining &~ full_penumbra_times
           
                #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
                ######         LOOK FOR PARTIAL PENUMBRA EPOCHS          ######
                ###############################################################
                
                partial_penumbra_times = remaining
                
                # PENUMBRA SPOT
                penumbra_spot = partial_penumbra_times & (cosphi < -cosOMEGA_i)
                
                if np.sum(penumbra_spot)>0:
                    
                    O2A4 = np.zeros([T])
                    O2A4[penumbra_spot]  =  Radj / sinOMEGA_i[penumbra_spot]
    
                    R_A4P  = np.zeros([2,T])#*np.nan   # I commented the nan product, since it multiplies the time by ~o(2)
                                        
                    R_A4P[0,penumbra_spot] = - d2_i[penumbra_spot] - O2A4[penumbra_spot]  #O2A3[partial_umbra,np.newaxis]
                    R_A4P[1,penumbra_spot] =   d1_i[penumbra_spot]                    
                    
                    
                    ################################
                    a = R_A4P[1,penumbra_spot]
                    b = R_A4P[0,penumbra_spot]
                    c = sinOMEGA_i[penumbra_spot]    
                    d = cosOMEGA_i[penumbra_spot]
                    cc= sinOMEGA_i[penumbra_spot]/cosOMEGA_i[penumbra_spot]
                    # If we enter this function, it is 100% sure that we have intersection
                    # In order to avoid the singularity at a = 0...
                    #a[a==0] = 1e-30


                    A  = 1/(2*a*d**2)
                    x0 = b*c**2
                    y0 = ( (c*b)**2 - (Radi**2 - a**2) ) / (2*a)                
    
                    x = (+2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta)[np.newaxis,:]/A[:,np.newaxis] - ( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta)[np.newaxis,:]/A[:,np.newaxis])**2 - 4*(x0[:,np.newaxis]**2 + y0[:,np.newaxis]/A[:,np.newaxis] + cc[:,np.newaxis]*b[:,np.newaxis]*np.cos(theta)[np.newaxis,:]/A[:,np.newaxis] ) )**0.5) / (2*Radi)
                    y = (cc[:,np.newaxis] * (x*Radi - b[:,np.newaxis]) * np.cos(theta)[np.newaxis,:] - a[:,np.newaxis]) / Radi
                    z = -(cc[:,np.newaxis] * (x*Radi - b[:,np.newaxis]) * np.sin(theta)[np.newaxis,:]    ) / Radi
                    
                    xyz   = np.einsum('ijk,jkl->ikl',R_1[:,:,penumbra_spot],np.array([x,y,z]))
                    
    #                flag  = np.sum(xyz[2]>0,1)
                    t = np.size(xyz,1)
                    
                    aux1 = np.diff((xyz[2]<0).astype(int),axis=1)
                    aux2 = [np.where(aux1[ii,:] != 0 )[0] for ii in range(t)]
    #                ref_line = np.array([np.arctan2(xyz[1,:,0],xyz[0,:,0]), xyz[0,:,0],xyz[1,:,0]])
                    #exterior= np.array([np.array([aux1[ii,aux2[ii][0]], -aux1[ii,aux2[ii][-1]]]) for ii in range(t)])
    #                n_patches = np.array([2 if (np.size(aux2[ii])==3 or (np.size(aux2[ii])==2 and exterior[ii,0]==-1)  or (np.size(aux2[ii])==4 and exterior[ii,0]==1) ) else \
    #                                      3 if (np.size(aux2[ii])==4 and exterior[ii,0]==-1) else  \
    #                                      1 for ii in range(t)])
    
    
#                    fig = plt.figure()
#                    ax = fig.add_subplot(111, projection='3d')            
#                    plt.plot(x[20,:],y[20,:],z[20,:],'o-')
#                    plt.plot(xyz[0,20,:],xyz[1,20,:],xyz[2,20,:],'o-')
#                    ax.set_xlim([-1, 1])
#                    ax.set_ylim([-1, 1])
#                    ax.set_zlim([-1, 1])   
#                    ax.set_xlabel('x')   
#                    ax.set_ylabel('y')   
#                    ax.set_zlabel('z')                     
#                    circle1 = plt.Circle((0,0), 1, color = 'k', fill=False, zorder=1)
#                    ax.add_artist(circle1)      
    
    
                    coord = np.zeros([2,t,N])
                                          
                    coord[0,xyz[2,:,:]>0] = -xyz[0,xyz[2,:,:]>0]
                    coord[1,xyz[2,:,:]>0] =  xyz[1,xyz[2,:,:]>0]
                    coord[0,:] =-xyz[0,:]
                    coord[1,:] = xyz[1,:]
                    
                                               
                    coord[0,:,:] = -np.array([  xyz[0,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[0,ii,aux2[ii][3]+1:]))                                                                                                                                      if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],2*np.pi),xyz[0,ii,aux2[ii][3]+1]) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])
                    coord[1,:,:] = np.array([  xyz[1,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi), N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[1,ii,aux2[ii][3]+1:]))                                                                                                                                                   if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],xyz[0,ii,aux2[ii][3]+1]),2*np.pi) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])                                            
                                               
    #                coord[1,:,:] = np.array([  np.concatenate((xyz[1,ii,0:aux2[ii][0]+1], np.sin( ref_line[0,ii] + np.linspace(np.mod(np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1])-ref_line[0,ii],2*np.pi),np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]],xyz[0,ii,aux2[ii][1]])-ref_line[0,ii],2*np.pi),aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:])) if (np.size(aux2[ii])==2 and exterior[ii,0]==1) else \
    #                                           np.zeros([N])   for ii in range(t)  ])                                      
                                          
                    coord   = np.einsum('jik,jkl->ikl',Rphi_s[:,:,penumbra_spot],coord)
                    
#                    fig , ax1 = plt.subplots()
#                    ax1.set_ylim([-1,1])
#                    ax1.set_xlim([-1,1])
#                    ax1.set_aspect('equal', adjustable='box')
#                    circle1 = plt.Circle((0,0), 1, color = 'k', fill=False, zorder=1)
#                    ax1.add_artist(circle1)
#                    plt.plot(coord[0,20,:],coord[1,20,:],'ko-')


                    bodies[i].geometry.penumbra_points[-1][0][:,penumbra_spot,:] = coord



                # PENUMBRA EDGE
                penumbra_edge = partial_penumbra_times & ~(cosphi < -cosOMEGA_i)
                
                
                if np.sum(penumbra_edge)>0:
                    ################################
                    O2A4 = np.zeros([T])
                    O2A4[penumbra_edge]  =  Radj / sinOMEGA_i[penumbra_edge]
    
                    R_A4P  = np.zeros([2,T])#*np.nan   # I commented the nan product, since it multiplies the time by ~o(2)
                                        
                    R_A4P[0,penumbra_edge] = - d2_i[penumbra_edge] - O2A4[penumbra_edge]  #O2A3[partial_umbra,np.newaxis]
                    R_A4P[1,penumbra_edge] =   d1_i[penumbra_edge]                    
                    
                    
                    ################################
                    a = R_A4P[1,penumbra_edge]
                    b = R_A4P[0,penumbra_edge]
                    c = sinOMEGA_i[penumbra_edge]    
                    d = cosOMEGA_i[penumbra_edge]
                    cc= sinOMEGA_i[penumbra_edge]/cosOMEGA_i[penumbra_edge]
                    
                    # If we enter this function, it is 100% sure that we have intersection
                    # In order to avoid the singularity at a = 0...
                    #a[a==0] = 1e-30
                    
                    t = np.size(cc)
                    
                    A  = 1/(2*a*d**2)
                    x0 = b*c**2
                    y0 = ( (c*b)**2 - (Radi**2 - a**2) ) / (2*a)
                    
                    theta0 = np.arccos( 2*A/c * (b-x0 + ((x0-b)**2 + y0/A )**0.5 ) )
                
#                    theta1 = np.zeros([round(N/2),t])
#                    theta2 = np.zeros([N-round(N/2)+2,t])
                    
                    i1 = np.arange(round(N/2))
                    i2 = np.arange(N-round(N/2)+2)
                    
                    theta1 = -theta0[:,np.newaxis]*np.cos(i1*np.pi/(round(N/2)-1))[np.newaxis,:]
                    theta2 =  theta0[:,np.newaxis]*np.cos(i2*np.pi/(N-round(N/2)+1))[np.newaxis,:]
                    
#                    for i in range(len(theta)):
#                        theta[i] = -theta0 * m.cos(i*m.pi/(round(N/2)-1))
#                        #theta[i] = -theta0 * m.cos(i*m.pi/(N-1))
#                    for i in range(len(theta1)):
#                        theta1[i] = theta0 * m.cos(i*m.pi/(len(theta1)-1))
#                        #theta1[i] = theta0 * m.cos(i*m.pi/(N-1))
                    theta2 = np.delete(theta2,np.size(theta2,1)-1,1)
                    theta2 = np.delete(theta2,0,1)  
                    
                    aux_1 = ( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta1)/A[:,np.newaxis])**2 - 4*(x0[:,np.newaxis]**2 + y0[:,np.newaxis]/A[:,np.newaxis] + cc[:,np.newaxis]*b[:,np.newaxis]*np.cos(theta1)/A[:,np.newaxis] ))
                    aux_1[aux_1<0]=0
                    aux_2 = ( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta2)/A[:,np.newaxis])**2 - 4*(x0[:,np.newaxis]**2 + y0[:,np.newaxis]/A[:,np.newaxis] + cc[:,np.newaxis]*b[:,np.newaxis]*np.cos(theta2)/A[:,np.newaxis] ))
                    aux_2[aux_2<0]=0
                    
                    x = np.concatenate(( (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta1)/A[:,np.newaxis] + ( aux_1 )**0.5) / (2*Radi), (2*x0[:,np.newaxis] + cc[:,np.newaxis]*np.cos(theta2)/A[:,np.newaxis] - ( aux_2 )**0.5) / (2*Radi)  ),1)
                    y = np.concatenate(( (cc[:,np.newaxis] * np.abs(x[:,0:round(N/2)]*Radi - b[:,np.newaxis]) * np.cos(theta1) - a[:,np.newaxis]) / Radi , ( cc[:,np.newaxis] * np.abs(x[:,round(N/2):]*Radi - b[:,np.newaxis]) * np.cos(theta2) - a[:,np.newaxis])/Radi ),1)
                    z = np.concatenate(( (cc[:,np.newaxis] * (x[:,0:round(N/2)]*Radi - b[:,np.newaxis]) * np.sin(theta1)       ) / Radi , ( cc[:,np.newaxis] * (x[:,round(N/2):]*Radi - b[:,np.newaxis]) * np.sin(theta2)   )/Radi  ),1)
                    
                    
#                    fig = plt.figure()
#                    ax = fig.add_subplot(111, projection='3d')            
#                    plt.plot(x[0,:],y[0,:],z[0,:],'o-')
#                    ax.set_xlim([-1, 1])
#                    ax.set_ylim([-1, 1])
#                    ax.set_zlim([-1, 1])   
#                    ax.set_xlabel('x')   
#                    ax.set_ylabel('y')   
#                    ax.set_zlabel('z')                     
#                    circle1 = plt.Circle((0,0), 1, color = 'k', fill=False, zorder=1)
#                    ax.add_artist(circle1)   
                    
                    xyz   = np.einsum('ijk,jkl->ikl',R_1[:,:,penumbra_edge],np.array([x,y,z]))
                    
    #                flag  = np.sum(xyz[2]>0,1)
                    
                    aux1 = np.diff((xyz[2]<0).astype(int),axis=1)
                    aux2 = [np.where(aux1[ii,:] != 0 )[0] for ii in range(t)]
    #                ref_line = np.array([np.arctan2(xyz[1,:,0],xyz[0,:,0]), xyz[0,:,0],xyz[1,:,0]])
                    #exterior= np.array([np.array([aux1[ii,aux2[ii][0]], -aux1[ii,aux2[ii][-1]]]) for ii in range(t)])
    #                n_patches = np.array([2 if (np.size(aux2[ii])==3 or (np.size(aux2[ii])==2 and exterior[ii,0]==-1)  or (np.size(aux2[ii])==4 and exterior[ii,0]==1) ) else \
    #                                      3 if (np.size(aux2[ii])==4 and exterior[ii,0]==-1) else  \
    #                                      1 for ii in range(t)])
    
                    coord = np.zeros([2,t,N])

                    coord[0,:,:] = -np.array([  xyz[0,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[0,ii,aux2[ii][3]+1:]))                                                                                                                                      if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],2*np.pi),xyz[0,ii,aux2[ii][3]+1]) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[0,ii,0:aux2[ii][0]+1], np.cos( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[0,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.cos( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[0,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.cos( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[0,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])
                    coord[1,:,:] = np.array([  xyz[1,ii,:] if (np.size(aux2[ii])==0) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                        if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod(np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi), N-aux2[ii][1]  ) )))                                                                                                                                                                                                                                                                                   if (np.size(aux2[ii])==2 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][3]], xyz[0,ii,aux2[ii][3]]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,aux2[ii][3]-aux2[ii][2]  )), xyz[1,ii,aux2[ii][3]+1:]))                                                                                                                                                   if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:aux2[ii][3]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][3]+1], xyz[0,ii,aux2[ii][3]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][3]+1],xyz[0,ii,aux2[ii][3]+1]),2*np.pi) ,N-aux2[ii][3]  ) ))) if (np.size(aux2[ii])==4 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,N-aux2[ii][0]  ) )))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:]))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 if (np.size(aux2[ii])==1 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.concatenate(( xyz[1,ii,0:aux2[ii][0]+1], np.sin( np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][1]], xyz[0,ii,aux2[ii][1]]) - np.arctan2(xyz[1,ii,aux2[ii][0]+1],xyz[0,ii,aux2[ii][0]+1]),2*np.pi) ,aux2[ii][1]-aux2[ii][0]  ) )   , xyz[1,ii,aux2[ii][1]+1:aux2[ii][2]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,-1], xyz[0,ii,-1]) - np.arctan2(xyz[1,ii,aux2[ii][2]+1],xyz[0,ii,aux2[ii][2]+1]),2*np.pi) ,N-aux2[ii][2]  ))))                                                                                                                                                                                            if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]== 1) else \
                                               np.concatenate(( np.sin( np.arctan2(xyz[1,ii,0],xyz[0,ii,0]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][0]], xyz[0,ii,aux2[ii][0]]) - np.arctan2(xyz[1,ii,0],xyz[0,ii,0]),2*np.pi) ,aux2[ii][0]  ) ),  xyz[1,ii,aux2[ii][0]+1:aux2[ii][1]+1],  np.sin( np.arctan2(xyz[1,ii,aux2[ii][1]+1], xyz[0,ii,aux2[ii][1]+1]) + np.linspace(0, np.mod( np.arctan2(xyz[1,ii,aux2[ii][2]], xyz[0,ii,aux2[ii][2]]) - np.arctan2(xyz[1,ii,aux2[ii][1]+1],xyz[0,ii,aux2[ii][1]+1]),2*np.pi) ,aux2[ii][2]-aux2[ii][1] )), xyz[1,ii,aux2[ii][2]+1:]))                                                                                                                                                                                                                               if (np.size(aux2[ii])==3 and aux1[ii,aux2[ii][0]]==-1) else \
                                               np.zeros([N])   for ii in range(t)  ])                                    

                    coord   = np.einsum('jik,jkl->ikl',Rphi_s[:,:,penumbra_edge],coord)                            
                                               
                bodies[i].geometry.penumbra_points[-1][0][:,penumbra_edge,:] = coord

                # Times for continuing the eclipse search:
                remaining = eclipsed &~ partial_penumbra_times                              
               
                print(np.sum(remaining))
                
                
                
                
    return bodies
                