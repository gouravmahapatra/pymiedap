# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

import numpy as np
import numpy.random as npr
import module_mie as mie
import module_mieshell as mieshell
import module_readmie as readmie
import module_dap as dap
import module_geos as geos
import os
import sys
import os.path
import matplotlib.pyplot as mpl

from .classes import Model, Layers, Layer, Geom, Aerosols
from .utils import calc_azimuth, get_cosbeta, sunblackbody

def dap_code(model, rename=False, output_name='modelA',
             path_output='./dap_database/',
             nlaysMAX=50, ncoefsMAX=4001, nfouMAX=4001,
             nmuMAX=201, nmatMAX=4, nmat=4, nmug=20):
    """ This function launches the DAP code to calculate the supermatrices
    produced by the doubling-adding code. Reads input from a model class.
    Requires the module module_dap.

    Parameters
    ----------
    model : Model object
        an input model containing all the modeling parameters
    rename: Bool
        if True, the user set output file is used for the resulting coefficient files
    output_name : string
        base of filename to be used for the Fourier output files
    path_output : string
        Path of folder where to store the output Fourier files. Folder is
        created if not existing.
    nlaysMAX : integer, optional
        Maximum number of layers. Default 50
    ncoefsMAX: int, optional
        max number of coefs for the Mie expansion. Default 4001.
        Warning: changing this might conflict with Fortan modules
    nfouMAX : integer, optional
        max number of Fourier coefficients. Default is 4001
    nmug : integer, optional
        number of Gauss points used in the calculations. Default 20
    nnuMAX: int, optional
        maximal number of Gauss points
    nmatMAX: int, optional
        max number of Stokes elements. Should be 1, 3 or 4.
        Default 4
    nmat: int, optional
        number of Stokes elements to compute. Should be 1, 3 or 4.
        Default 4

    Returns
    -------
        produces an output file containing Fourier coefficients for the given
        model

    """

    # Checking that the output directory exists
    # if not, creating it
    try:
        os.makedirs(os.path.normpath(path_output))
    except OSError:
        if not os.path.isdir(os.path.normpath(path_output)):
            raise

    #-----------------------------------------------------------------------
    #     Some parameter values:
    #
    #     nmat: size of Stokes parameters (1, 3, or 4)
    #     nmug: number of Gauss points
    #     asurf: surface albedo
    #     nlays: the number of atmospheric layers
    #     nmat: size of Stokes parameters (1, 3, or 4)
    #     nmug: number of Gauss points
    #     asurf: surface albedo
    #-----------------------------------------------------------------------
    nsupMAX = nmuMAX * nmatMAX

    # Reading basic properties of the model
    asurf = model.asurf
    surfmat = surface_check(model,nmug, nmat=nmat, nmuMAX=nmuMAX)
    nlays = len(vars(model.layers))  # number of atmospheric layers
    nwvl = len(model.wvl_list)

    ncoefs = np.zeros(nlaysMAX, order='F')
    wavels = model.wvl_list
    nwavels = len(wavels)

    # Creating array to receive the coefficients
    coefin = np.zeros((nmatMAX,nmatMAX,ncoefsMAX,nlaysMAX), order='F')
    coefs = np.zeros((nmatMAX,nmatMAX,ncoefsMAX,nlaysMAX), order='F')
    ncoefin = np.zeros((nlaysMAX), order='F')
    model.name = ['']*len(model.wvl_list)

    # some global parameters
    gravity = model.gravity
    dpol = model.dpol
    mma = model.mma

    print('Beginning of DAP program')

    # ---------------------
    # Loop on wavelengths
    # ---------------------
    for z, wav in enumerate(wavels):
        #get coefs for each wavelength
        taus = np.zeros(nlaysMAX, order='F')  # all values of tau at wvl z
        taus_g = np.zeros(nlaysMAX, order='F')  # all values of tau_g at wvl z
        basca = np.zeros(nlaysMAX, order='F')  # all values of basca
        baabs = np.zeros(nlaysMAX, order='F')  # all values of baabs
        ssalbs = np.zeros(nlaysMAX, order='F')  # all values of s.scat albedoes
        pres = np.zeros(nlaysMAX, order='F')  # pressure levels for each layer

        print('Wavelength {:06.7f} microns'.format(wav))

        # Loop on layers
        l = 0  # layer number
        for layer_name, layer in vars(model.layers).items():

            #  reading information contained in all layers
            # reading the coefficients of each layer for the current wavelength
            taus[l] = layer.tau[z]  #value of tau at wvl z in layer l
            taus_g[l] = layer.tau_g[z]  #value of tau_g at wvl z in layer l
            pres[l] = layer.press

            # if the layer is transparent, then ignore the aerosols
            # coefficients
            if max(layer.tau)!=0:
                coefin[:,:,:,l] = layer.mixed_aerosols.coefs[z,:,:,:]
                ncoefin[l] = layer.mixed_aerosols.ncoefs[z]
                ssalbs[l] = layer.mixed_aerosols.ssalb[z]
                # Calculate the aerosol scattering and absorption optical
                # thicknesses,
                basca[l] = ssalbs[l] * taus[l]
                baabs[l] = (1. - ssalbs[l]) * taus[l]

            print('{}.sc.{:06.7f}'.format(layer.mixed_aerosols.typ, wav))
            l = l + 1

        layorder = np.argsort(pres[:l]) # finding the order of pressures
        layorder = layorder[::-1]  # putting layers in order (highest pres below)

        taus[:l] = taus[layorder]  # putting the taus in right order
        taus_g[:l] = taus_g[layorder]  # putting the taus_g in right order
        pres[:l] = pres[layorder]  # putting the pressures in right order

        # More reordering
        coefin[:,:,:,:l] = coefin[:,:,:,layorder]
        ncoefin[:l] = ncoefin[layorder]
        ssalbs[:l] = ssalbs[layorder]
        basca[:l] = basca[layorder]
        baabs[:l] = baabs[layorder]

        # Calculate the molecular parameters of the atmosphere:
        ri = model.rindex_gas[z]
        bmsca, bmabs, coefsm = dap.bmolecules(wav, nlays, pres, dpol, ri, mma, gravity)

        model.coefsm = coefsm

        # Storing the effective scattering and absorption opacities
        for layer_name, layer in vars(model.layers).items():
            # identify position of the current layer
            lev = (pres==layer.press)

            # force user-define rayleigh opacity
            if layer.rayscat is False:
                bmsca[lev] = layer.tau_ray[z]

            layer.bmsca[z] = bmsca[lev]
            layer.bmabs[z] = bmabs[lev]
            layer.basca[z] = basca[lev]
            layer.baabs[z] = baabs[lev]

        #---------------------------------------------------------------------
        #     Open the Fourier coefficients output file:
        #---------------------------------------------------------------------
        outputname = 'fou_' + '{:4.3f}'.format(wav) + '.dat'

        #---------------------------------------------------------------
        #     Calculate the combined expansion coefficients
        #---------------------------------------------------------------
        for i in np.arange(nlays):
            ncoefs[i] = max(ncoefin[i], 2)
            # multiply all coefs by the associated optical thickness
            # in each layer
            com = bmsca[i,np.newaxis,np.newaxis] * coefsm
            coa = basca[i,np.newaxis,np.newaxis,np.newaxis] * coefin[:,:,:,i]
            # if the opacities are too small, nullify coefs
            if ((bmsca[i]+basca[i]) < 1e-10):
                coefs[:,:,:,i]= 0.
            else:
                # otherwise, scale coefs according relative opacities
                coefs[:,:,:,i]= (com+coa)/(bmsca[i]+basca[i])

        model.coefst = coefs

        #-----------------------------------------------------------------------
        #     Calculate the total optical thicknesses and albedos:
        #-----------------------------------------------------------------------
        bmabs = taus_g
        b = bmsca + basca + baabs + bmabs
        model.bmsca = bmsca
        model.basca = basca
        model.bmabs = bmabs
        model.baabs = baabs
        model.b = b
        a = np.zeros(len(b))
        # If b is zero, then a should also be 0
        #a[b==0] = 0.0
        a[b!=0] = (bmsca[b!=0]+basca[b!=0])/b[b!=0]
        model.a = a

        #-----------------------------------------------------------------------
        #     Call the doubling-adding routine:
        #-----------------------------------------------------------------------
        dap.adding(outputname,a,b,coefs,ncoefs,nlays,nmug,nmat,surfmat)

        # Naming the model with check for Windows paths
        print('fou_{:4.7f}.dat'.format(wav))
        if rename is True:
            output_file = path_output + output_name + '_{:4.7f}.dat'.format(wav)
            output_file = os.path.normpath(output_file)
            model.name[z] = output_file
            os.rename('fou_{:4.3f}.dat'.format(wav),output_file)
        else:
            output_file = path_output + 'fou_{:4.7f}.dat'.format(wav)
            output_file = os.path.normpath(output_file)
            model.name[z] = output_file
            os.rename(outputname,output_file)
        print('End of DAP program')


def read_dap_output(phase, sza, emission, filename, beta=None, phi=None,
                    ngeosMAX=100000, nmuMAX=300, nfouMAX=2000, nmatMAX=4):
    """ Reads the supermatrices coefficients for given geometry

    Parameters
    ----------
    phase : float or array
        phase angle (deg)
    sza : float or array
        solar zenith angle (deg)
    emission : float or array
        emission angle (deg)
    filename : string
        name of the Fourier file to be read
    beta : None, float or array, optional
        angle between the meridian plane and the scattering plane (deg)
    phi : None, float or array, optional
        azimuthal angle (deg)
    nfouMAX : integer, optional
        max number of Fourier coefficients. Default is 4001
    ngeosMAX : integer, optional
        max number of geometries. Default 100000
    nnuMAX: int, optional
        maximal number of Gauss points
    nmatMAX: int, optional
        max number of Stokes elements. Should be 1, 3 or 4.
        Default 4


    Returns
    -------
    I : array (same shape as phase)
        Stokes element I, normalised with input flux unity (not Pi)
    Q : array (same shape as phase)
        Stokes element Q, normalised with input flux unity (not Pi)
    U : array (same shape as phase)
        Stokes element U, normalised with input flux unity (not Pi)
    V : array (same shape as phase)
        Stokes element V, normalised with input flux unity (not Pi)

    All elements are given assuming normal reflection (i.e. multiply by
    cos(theta0) if you want real observed flux)

    """

    ngeos = len(phase)
    betaF = np.zeros(ngeosMAX, order='F')
    azimuthF = np.zeros(ngeosMAX, order='F')

    if phi is None:
        # Getting needed geometry
        azimuth = calc_azimuth(phase, sza, emission, deg=False)  # azimuth angle in rads
        azimuth = np.pi - azimuth  # corr. due to diff definitions
        #azimuth = azimuth  # corr. due to diff definitions
        azimuthF[:ngeos] = np.degrees(azimuth)
    else:
        azimuth = phi[:ngeos]
        azimuthF[:ngeos] = phi[:ngeos]


    if beta is None:
        azimuth = azimuth[:ngeos]
        beta = get_cosbeta(phase, sza, emission, azimuth)  # rotation angle beta
        beta = np.arccos(beta) #warning: still in radians here
        betaF[:ngeos] = np.degrees(beta)
    else:
        betaF[:ngeos] = beta[:ngeos]


    ngeos = len(phase)

    #Preparing vectors for FORTRAN function
    szaF = np.zeros(ngeosMAX, order='F')
    emissionF = np.zeros(ngeosMAX, order='F')

    szaF[:ngeos] = sza
    emissionF[:ngeos] = emission
    # make sure all input angles are in degrees

    # Reading Stoke vector
    rfou = np.zeros((nmatMAX*nmuMAX,nmuMAX,nfouMAX+1), order='F')
    Sv = geos.read_dap(filename, ngeos, szaF, emissionF, azimuthF, betaF, rfou)
    del(rfou)

    # storing output in proper Stokes elements
    I = Sv[0,:ngeos]
    Q = Sv[1,:ngeos]
    U = Sv[2,:ngeos]
    V = Sv[3,:ngeos]

    # scaling unit input flux
    I = I/np.pi
    Q = Q/np.pi
    U = U/np.pi
    V = V/np.pi

    return I,Q,U,V

def surface_check(model,nmug, nmat=4, nmuMAX=201):
    """ A function to check what type of surface is defined by the user, and
    act accordingly for the DAP code"""

    if type(model.surface)==np.ndarray:
        mus, smf, Lfin = fourier_matrix(nmug=nmug, surf_mat=model.surface, nmat=nmat, nmuMAX=nmuMAX)
    elif model.surface==str:
        print('read fourier file for surface')

    return Lfin

def fourier_matrix(nmug=20, surf_mat=np.diag([1,0,0,0]), nmat=4, nmuMAX=201, nmatMAX=4):
    """ A function to Fourier develop a scattering matrix for a surface
    """

    # Getting the Gauss points
    xs, w = dap.gauleg(201,nmug,0.,1.)
    smf = np.sqrt(2*xs*w)
    # adding extra nadir direction
    xs[nmug] = 1.0
    w[nmug] = 1.0
    smf[nmug] = 1.0
    nmu = nmug+1

    #mu = xs
    #mup = xs
    L = np.zeros((nmat,nmat,nmu,nmu), order='F')

    # if nmat<4, we take the corresponding submatrix of the input surface
    # matrix
    L[0:,0:,:,:] = surf_mat[:nmat,:nmat,np.newaxis,np.newaxis]
    #Bmp = np.zeros((nmat,nmat,nphi,nm))
    #Bmm = np.zeros((nmat,nmat,nphi,nm))

    #for im,mm in enumerate(m):
    #    for ip,pp in enumerate(phi):
    #        Bmp[0:,0:,ip,im] = np.diag([np.cos(mm*pp),np.cos(mm*pp),np.sin(mm*pp),np.sin(mm*pp)])
    #        Bmm[0:,0:,ip,im] = np.diag([-np.sin(mm*pp),-np.sin(mm*pp),np.cos(mm*pp),np.cos(mm*pp)])

    #Bm = Bmp + Bmm

    #Lm = np.einsum('jcbm,ijklp->ijklpm',Bm,L)

    #L2 = scpinteg.simps(Lm,phi,axis=4)
    #L2 = L2/(phi[-1]-phi[0])

    Lfin = np.zeros((nmatMAX*nmuMAX,nmatMAX*nmuMAX),order='F')
    for k in np.arange(nmat):
        for l in np.arange(nmat):
            for i in np.arange(nmu):
                for j in np.arange(nmu):
                    Lfin[nmat*i+l,nmat*j+k] = smf[i]*L[l,k,i,j]*smf[j]

    return xs,smf,Lfin


