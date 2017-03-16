import numpy as np
import pymiedap as pmd
import module_geos as geos

def combine(bodies, reference):

	print('\n    ... combining radiance results\n')

	I = np.zeros_like(reference.ephemeris.time)
	Q = np.zeros_like(I)
	U = np.zeros_like(I)
	V = np.zeros_like(I)

	ones  = np.ones_like(reference.ephemeris.time)
	zeros = np.zeros_like(reference.ephemeris.time)

	for body in bodies:

		size_scale = body.properties.R/reference.properties.R
		distance_scale = (reference.ephemeris.r_s/body.ephemeris.r_s)**2

		IQUV = [body.radiance.I, body.radiance.Q, body.radiance.U, body.radiance.V]

        	# this rotation could be written in a much much simpler way!
        	L = np.array([[ones , zeros, zeros, zeros],
                      [zeros,
                       np.cos(2*body.geometry.ref_plane_to_ref_line_angle),
                       np.sin(2*body.geometry.ref_plane_to_ref_line_angle),
                       zeros],
                      [zeros,-np.sin(2*body.geometry.ref_plane_to_ref_line_angle),
                       np.cos(2*body.geometry.ref_plane_to_ref_line_angle),
                       zeros],
                      [zeros, zeros, zeros, ones]])

        	IQUV  = np.einsum('ijt,jt->it',L,IQUV)

        	body.radiance.I_ref = distance_scale * size_scale * IQUV[0,:]
        	body.radiance.Q_ref = distance_scale * size_scale * IQUV[1,:]
        	body.radiance.U_ref = distance_scale * size_scale * IQUV[2,:]
        	body.radiance.V_ref = distance_scale * size_scale * IQUV[3,:]

        	I = I + body.radiance.I_ref
       		Q = Q + body.radiance.Q_ref
        	U = U + body.radiance.U_ref
        	V = V + body.radiance.V_ref

	return I, Q, U, V,


def integration(body):

    ngeosMAX=200000

    print('\n    ... integrating radiance on ' + body.type + ' ' + body.name + ' disk \n')

    #files = dict(np.genfromtxt('exopy/scenes.dat', dtype='str'))

    time = body.ephemeris.time
    wvl_list = body.atmosphere.wvl_list
    nwvl = len(wvl_list)

    body.grid.I = np.zeros((nwvl,np.size(body.grid.shadow)))
    body.grid.Q = np.zeros((nwvl,np.size(body.grid.shadow)))
    body.grid.U = np.zeros((nwvl,np.size(body.grid.shadow)))
    body.grid.V = np.zeros((nwvl,np.size(body.grid.shadow)))

    Ip = np.zeros((nwvl,np.size(time)))
    Qp = np.zeros((nwvl,np.size(time)))
    Up = np.zeros((nwvl,np.size(time)))
    Vp = np.zeros((nwvl,np.size(time)))

    area  = np.repeat(body.grid.area[:,np.newaxis],4,1).T
    phase = np.degrees(body.geometry.phase_angle)
    sza   = np.degrees(body.grid.solar_zenith_angle)
    emission = np.degrees(body.grid.observer_zenith_angle)
    beta  = np.degrees(body.grid.beta)
    phi   = np.degrees(body.grid.azimuth)

    scene = body.properties.fourier_scene


    for l,wvl in enumerate(wvl_list):
        for i,j in enumerate(time):

            print i+1, ' out of ', len(time)
            A = body.grid.shadow[i,:]>10E-10

            IQUV = np.zeros([4,body.grid.N_points])

            ######################################################################################
            #######################       PMD.READ_DAP_OUTPUT      ###############################
            ######################################################################################

            I,Q,U,V = pmd.read_dap_output(np.repeat(phase[i], sum(A)), sza[i,A],
                                        emission[A], body.atmosphere.name[l],
                                        beta=beta[i,A],
                                        phi=phi[i,A])

            # renormalizing
            IQUV[0,A] = 4*np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[0,A]*I/np.pi
            IQUV[1,A] = 4*np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[1,A]*Q/np.pi
            IQUV[2,A] = 4*np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[2,A]*U/np.pi
            IQUV[3,A] = 4*np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[3,A]*V/np.pi

            # storing disk resolved case
            body.grid.I[l,i,:] = IQUV[0,:]
            body.grid.Q[l,i,:] = IQUV[1,:]
            body.grid.U[l,i,:] = IQUV[2,:]
            body.grid.V[l,i,:] = IQUV[3,:]

            # storing disk integrated case
            Ip[l,i] = np.nansum(IQUV[0,:])
            Qp[l,i] = np.nansum(IQUV[1,:])
            Up[l,i] = np.nansum(IQUV[2,:])
            Vp[l,i] = np.nansum(IQUV[3,:])

        body.radiance.I = Ip
        body.radiance.Q = Qp
        body.radiance.U = Up
        body.radiance.V = Vp

    return body,

