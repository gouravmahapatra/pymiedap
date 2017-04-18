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

        # Rotation of the Stokes vectors into observer plane
        nQ, nU = pmd.rotate_stokes(body.radiance.Q, body.radiance.U,
                                   body.geometry.ref_plane_to_ref_line_angle)

        body.radiance.Q = nQ
        body.radiance.U = nU

        body.radiance.I_ref = distance_scale * size_scale * body.radiance.I
        body.radiance.Q_ref = distance_scale * size_scale * body.radiance.Q
        body.radiance.U_ref = distance_scale * size_scale * body.radiance.U
        body.radiance.V_ref = distance_scale * size_scale * body.radiance.V

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
    wvl_list = np.array(body.atmosphere.wvl_list)
    nwvl = len(wvl_list)

    ngrids = np.shape(body.grid.shadow)
    ngrids = (nwvl,) + ngrids
    body.grid.I = np.zeros(ngrids)
    body.grid.Q = np.zeros(ngrids)
    body.grid.U = np.zeros(ngrids)
    body.grid.V = np.zeros(ngrids)

    ngrids = np.shape(time)
    ngrids = (nwvl,) + ngrids
    Ip = np.zeros(ngrids)
    Qp = np.zeros(ngrids)
    Up = np.zeros(ngrids)
    Vp = np.zeros(ngrids)

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

            # If atmosphere not calculated yet
            if body.atmosphere.name[l] == '':
                pmd.compute_model(body.atmosphere)

            ######################################################################################
            #######################       PMD.READ_DAP_OUTPUT      ###############################
            ######################################################################################

            I,Q,U,V = pmd.read_dap_output(np.repeat(phase[i], sum(A)), sza[i,A],
                                        emission[A], body.atmosphere.name[l],
                                        beta=beta[i,A],
                                        phi=phi[i,A])

            # renormalizing
            IQUV[0,A] = np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[0,A]*I
            IQUV[1,A] = np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[1,A]*Q
            IQUV[2,A] = np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[2,A]*U
            IQUV[3,A] = np.cos(np.radians(sza[i,A]))*body.grid.shadow[i,A]*area[3,A]*V

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

