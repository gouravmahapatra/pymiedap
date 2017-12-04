# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

# -*- coding: utf-8 -*-
# ==================================================================
# EXOPY module: exopy_radiance.py
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
# Script containing the functions required for computing the indivi-
# dual and combined reflected signals at each time epoch.
#
# LIST OF FUNCTIONS
# ------------------------------------------------------------------
#  - combine: Function merging the radiance outputs of a planet and
# 	    moon bodies.
#  - integration: Function computing the Stokes elements at each pi-
# 		xel and time epoch and integrating them along the
# 		bodies' disks.
#
#

import time as t
import numpy as np
import pymiedap.pymiedap as pmd

def combine(bodies, reference):
    """
    ==================================================================
    EXOPY function: combine()
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    Dependences:

    DESCRIPTION
    ------------------------------------------------------------------
    Function dedicated to merge the Stokes elements for a set of
    planet and moon bodies.

    INPUTS
    ------------------------------------------------------------------
    - bodies: List of two body objects [-] (list)
    - reference: Reference body [-] ('body' object)

    OUTPUTS
    ------------------------------------------------------------------
    - I: First stokes vector: flux [normalized] (numpy array)
    - Q: Second stokes vector: linear polarization [normalized] (numpy
    array)
    - U: Third stokes vector: linear polarization [normalized] (numpy
    array)
    - V: Fourth stokes vector: circular polarization [normalized]
    (numpy array)
    - P: Degree of linear polarization [%] (numpy array)
    - Chi: Angle of polarization [deg] (numpy array)

    """

    print('\n    ... combining radiance results\n')

    I   = np.zeros_like(reference.radiance.I)
    Q   = np.zeros_like(I)
    U   = np.zeros_like(I)
    V   = np.zeros_like(I)
    P   = np.zeros_like(I)
    Chi = np.zeros_like(I)

    #ones  = np.ones_like(reference.ephemeris.time)
    #zeros = np.zeros_like(reference.ephemeris.time)

    for body in bodies:
        body.radiance.I_ref = np.zeros_like(I)
        body.radiance.Q_ref = np.zeros_like(I)
        body.radiance.U_ref = np.zeros_like(I)
        body.radiance.V_ref = np.zeros_like(I)

    for wvl,j in enumerate(reference.atmosphere.wvl_list):
        for body in bodies:

            size_scale = (body.properties.R/reference.properties.R)**2
            distance_scale = (reference.ephemeris.r_s/body.ephemeris.r_s)**2

            #IQUV = [body.radiance.I[wvl,:], body.radiance.Q[wvl,:], body.radiance.U[wvl,:], body.radiance.V[wvl,:]]

            # Rotation of the Stokes vectors into observer plane
            nQ, nU = pmd.rotate_stokes(body.radiance.Q[wvl,:], body.radiance.U[wvl,:],
                                       body.geometry.ref_plane_to_ref_line_angle)

            body.radiance.Q[wvl,:] = nQ
            body.radiance.U[wvl,:] = nU

            body.radiance.I_ref[wvl,:] = distance_scale * size_scale * body.radiance.I[wvl,:]
            body.radiance.Q_ref[wvl,:] = distance_scale * size_scale * body.radiance.Q[wvl,:]
            body.radiance.U_ref[wvl,:] = distance_scale * size_scale * body.radiance.U[wvl,:]
            body.radiance.V_ref[wvl,:] = distance_scale * size_scale * body.radiance.V[wvl,:]

            I[wvl,:] = I[wvl,:] + body.radiance.I_ref[wvl,:]
            Q[wvl,:] = Q[wvl,:] + body.radiance.Q_ref[wvl,:]
            U[wvl,:] = U[wvl,:] + body.radiance.U_ref[wvl,:]
            V[wvl,:] = V[wvl,:] + body.radiance.V_ref[wvl,:]

        P[wvl,:]   = np.sqrt(Q[wvl,:]**2 + U[wvl,:]**2)/I[wvl,:]*100
        Chi[wvl,:] = 0.5*np.rad2deg(np.arctan(U[wvl,:]/Q[wvl,:]))

    return I, Q, U, V, P, Chi


def integration(body, path_input = './dap_database/', nmug = 20, nmug_mie = 20, nmat=4, nsubr=50):
    """
    ==================================================================
    EXOPY function: integration()
    Delft University of Technology
    ------------------------------------------------------------------
    Author: Javier Berzosa Molina, Loic Rossi, Daphne Stam
    Date: 2016-2017
    ------------------------------------------------------------------

    Dependences:

    DESCRIPTION
    ------------------------------------------------------------------
    Function computing the Stokes elements at each pixel and time
    epoch and integrating them along the bodies' disks.

    INPUTS
    ------------------------------------------------------------------
    - body: Planet or moon type of body object [-] ('body' object)
    - path_input: Path to the Fourier files storage folder [-]
            ('body' object)
    - nmug: number of Gauss points for DAP calculations [-] (int)
    - nmug_mie: number of Gauss points for MIE calculations [-] (int)
    - nmat: number of Stokes elements to compute [-] (int)
    - nsubr: number of subintervals for the distribution [-] (int)

    OUTPUTS
    ------------------------------------------------------------------
    - body: Planet or moon type of body object [-] ('body' object)


    """

    ngeosMAX=200000

    print('\n    ... integrating radiance on ' + body.type + ' ' + body.name + ' disk \n')

    #files = dict(np.genfromtxt('exopy/scenes.dat', dtype='str'))

    t0 = t.time()

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

    for l,wvl in enumerate(wvl_list):
        for i,j in enumerate(time):

            print(i+1, ' out of ', len(time))
            A = body.grid.shadow[i,:]>10E-10

            IQUV = np.zeros([4,body.grid.N_points])

            # If atmosphere not calculated yet
            if body.atmosphere.name[l] == '':
                if hasattr(body.atmosphere, 'tag'):
                    if body.atmosphere.tag != '':
                        pmd.compute_model(body.atmosphere, rename=True,
                                          output_name=body.atmosphere.tag,
                                          path_input=path_input, nmug=nmug,
                                          nmug_mie=nmug_mie, nmat=nmat,
                                          nsubr=nsubr)
                else:
                    pmd.compute_model(body.atmosphere,
                                    path_input=path_input, nmug=nmug,
                                    nmug_mie=nmug_mie, nmat=nmat,
                                    nsubr=nsubr)

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

    body.radiance.I   = Ip
    body.radiance.Q   = Qp
    body.radiance.U   = Up
    body.radiance.V   = Vp
    body.radiance.P   = np.sqrt(Qp**2 + Up**2)/Ip*100
    body.radiance.Chi = 0.5*np.rad2deg(np.arctan(Up/Qp))

    body.radiance.t = t.time() - t0

    return body,
