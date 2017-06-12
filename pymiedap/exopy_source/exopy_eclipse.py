# -*- coding: utf-8 -*-

# ==================================================================
# EXOPY module: exopy_eclipse.py
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
# The 'expoy_eclipse' script contain the functions required for the
# computation of the pixel darkening of the planetary system bodies
# due to shadowing through eclipses.
#
# LIST OF FUNCTIONS
# ------------------------------------------------------------------
#  - eclipse: Function computing the eclipses shadowing of the extra-
# 	    solar planetary system.
#
#


import numpy as np
import exopy_config as _cfg
#from exopy_functions import PolyArea

def eclipse(bodies, star = None):
    """
    ==================================================================
    EXOPY function: eclipse()
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    Dependences:

    DESCRIPTION
    ------------------------------------------------------------------
    Computes the degree of pixel darkening due to eclipses shadowing
    among the different bodies.

    INPUTS
    ------------------------------------------------------------------
    - bodies: list comprising a planet & moon object [-] (list)
    - star: star type of body object [-] (body object)

    OUTPUTS
    ------------------------------------------------------------------
    - bodies: updated list comprising a planet & moon object [-]
        (list)


    """

    #bodies = [Moon,Earth]o
    #star = Sunp

    #for body in bodies:
    #    if not hasattr(body.grid, 'shadow'):
    #       body.grid.shadow = np.zeros([len(body.ephemeris.time), len(body.grid.nodes)])

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
		######            LOOK FOR FULL UMBRA NODES              ######
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
        #   d_nod_norm  = bodies[i].grid.distance_nodes_ob
		d2_nod      = np.zeros([T,bodies[i].grid.N_points])
		d1_nod      = np.zeros([T,bodies[i].grid.N_points])
        #   d2_nod1     = np.zeros([T,bodies[i].grid.N_points])
        #   d1_nod1     = np.zeros([T,bodies[i].grid.N_points])
        #   aux         = np.zeros([T,bodies[i].grid.N_points])

        #print len(remaining)
        #print np.shape(Rphi_s)
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
		######           LOOK FOR PARTIAL UMBRA NODES            ######
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


    return bodies

