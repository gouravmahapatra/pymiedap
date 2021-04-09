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

def mie_code(aerosols, wavelengths, output=False, delta=1e-8, cutoff=1e-8, thmin=0, thmax=180,
             nsubr=50, ngaur=60, ncoefsMAX=4001,
             nfouMAX=4001, nmuMAX=201, nmatMAX=4):
    """ Compute Mie expansion coefficients for Aerosol object.
    Requires the module_mie module.

    Parameters
    ----------
    aerosols : Aerosol object
        an input aerosol type model containing all the modeling parameters
    wavelengths : array
        list of wavelengths for which computations should be performed
    delta: float, optional
        truncation of the Mie sum. Default 1e-8
    cutoff: float, optional
        cutoff value for the particle size distribution. Default 1e-8
    thmin: float, optional
        minimal value of phase angle (in degrees). Default 0
    thmax: float, optional
        maximal phase angle. Default 180.
    nsubr: integer, optional
        number of subintervals for the distribution. Default 50
    ngaur: integer, optional
        number of Gauss points used in the calculations. Default 60
    ncoefsMAX: int, optional
        max number of coefs for the Mie expansion. Default 4001.
        Warning: changing this might conflict with Fortan modules
    nmatMAX: int, optional
        number of Stokes elements to compute. Should be 1, 3 or 4.
        Default 4
    output : Bool, optional
        if True, expansion coefficients are written in a file with name in the form
        scfile_name = aerosols.specie + '.sc.' + '{:06.3f}'.format(wav)

    Returns
    -------
    Stores expansion coefficients in aerosols object

    """

    # ---------------------------
    # Mie calculations parameters
    # ---------------------------

    weight2 = 0.3

    nwavels = len(wavelengths)

    # Creating array to receive the coefficients
    supercoefin = np.zeros((nwavels,nmatMAX,nmatMAX,ncoefsMAX), order='F')
    superncoefin = np.zeros(nwavels, order='F')

    qexts = np.zeros(nwavels)
    sexts = np.zeros(nwavels)
    qscas = np.zeros(nwavels)
    sscas = np.zeros(nwavels)
    asyms = np.zeros(nwavels)

    #Getting aerosols parameters
    r_eff = aerosols.r_eff
    v_eff = aerosols.v_eff
    par3 = aerosols.par3
    specie = aerosols.typ
    idis = aerosols.psd  # index of the particle size dist

    # Getting extrema for radii
    rmin, rmax = mie.rminmax(idis, r_eff, v_eff, par3, weight2, cutoff)
    # -----
    # LOOP ON WAVELENGTHS
    # -----

    print('Beginning of Mie program')

    for i, wav in enumerate(wavelengths):
        print('Wavelength {:06.7f}'.format(wav))
        nr = aerosols.nr[i]
        ni = aerosols.ni[i]
        m = nr - 1j * ni

        scfile_name = specie + '.sc.' + '{:06.3f}'.format(wav)

        # calculation of the scattering matrix
        u, wg, F, miec, nangle = mie.scatmat(m, wav, idis, nsubr, ngaur, rmin,
                                             rmax, r_eff, v_eff, par3, weight2,
                                             delta)

        ncoefs = nangle

        # Generalization to 6 matrix elements
        Fshape = np.shape(F)
        F2 = np.zeros((6,Fshape[1]), order='F')
        F2[0,:] = F[0,:]  # F11
        F2[1,:] = F[0,:]  # F22
        F2[2,:] = F[2,:]  # F33
        F2[3,:] = F[2,:]  # F44
        F2[4,:] = F[1,:]  # F12
        F2[5,:] = F[3,:]  # F34

        #expansion of the matrix
        coefs = matrix_expansion(ncoefs, nangle, u, wg, F2)

        # Store the coefficients for each wvl
        supercoefin[i,:,:,:] = coefs
        superncoefin[i] = ncoefs

        qexts[i] = miec[3] # get Qext somewhere
        sexts[i] = miec[1] # get Sext somewhere
        qscas[i] = miec[1] # get Qext somewhere
        sscas[i] = miec[0] # get Sext somewhere

        #-----------------------------------------------------------------------
        #        Calculate the asymmetry parameter:
        #-----------------------------------------------------------------------
        if nangle >= 1:
            cosbar = coefs[0, 0, 1] / 3.
        else:
            cosbar = 0.

        asyms[i] = cosbar

        #------------------------------------------------------------------
        #        Write the coefficients to the output file:
        #------------------------------------------------------------------
        if output is True:
            mie.writsc(scfile_name, idis, nsubr, ngaur, coefs, ncoefs, cosbar,
                        miec, wav, nr, ni, rmin, rmax, r_eff, v_eff, par3)

        # end of wvl loop


    # Storing the coefficients in the aerosols object
    aerosols.coefs = supercoefin
    aerosols.ncoefs = superncoefin
    aerosols.asym = asyms
    aerosols.sext = sexts
    aerosols.qext = qexts
    aerosols.ssca = sscas
    aerosols.qsca = qscas
    aerosols.ssalb = sscas/sexts

    print('End of Mie program')


def mie_shell(aerosols, wavelengths, output=False, delta=1e-8, cutoff=1e-8, thmin=0, thmax=180,
              nsubr=20, ngaur=20, nlaysMAX=50, ncoefsMAX=4001,
              nfouMAX=4001, nmuMAX=201, nmatMAX=4):
    """ Compute Layered spheres Mie expansion coefficients for Aerosol object.
    Requires the module_mieshell module.

    Parameters
    ----------
    aerosols : Aerosol object
        an input aerosol type model containing all the modeling parameters
    wavelengths : array
        list of wavelengths for which computations should be performed
    delta: float, optional
        truncation of the Mie sum. Default 1e-8
    cutoff: float, optional
        cutoff value for the particle size distribution. Default 1e-8
    thmin: float, optional
        minimal value of phase angle (in degrees). Default 0
    thmax: float, optional
        maximal phase angle. Default 180.
    nsubr: integer, optional
        number of subintervals for the distribution. Default 50
    ngaur: integer, optional
        number of Gauss points used in the calculations. Default 60
    ncoefsMAX: int, optional
        max number of coefs for the Mie expansion. Default 4001.
        Warning: changing this might conflict with Fortan modules
    nmatMAX: int, optional
        number of Stokes elements to compute. Should be 1, 3 or 4.
        Default 4
    output : Bool, optional
        if True, expansion coefficients are written in a file with name in the form
        scfile_name = aerosols.specie + '.sc.' + '{:06.3f}'.format(wav)

    Returns
    -------
    Stores expansion coefficients in aerosols object

    """

    # ---------------------------
    # Mie calculations parameters
    # ---------------------------
    par3 = 0.25  # this last parameter is only used for some PSD
    weight2 = 0.0

    nwavels = len(wavelengths)

    # Creating array to receive the coefficients
    supercoefin = np.zeros((nwavels,nmatMAX,nmatMAX,ncoefsMAX), order='F')
    superncoefin = np.zeros(nwavels, order='F')

    qexts = np.zeros(nwavels)
    sexts = np.zeros(nwavels)
    qscas = np.zeros(nwavels)
    sscas = np.zeros(nwavels)
    asyms = np.zeros(nwavels)

    #Getting aerosols parameters
    r_eff = aerosols.r_eff
    v_eff = aerosols.v_eff
    specie = aerosols.typ
    idis = aerosols.psd  # index of the particle size dist

    # Getting extrema for radii
    rmin, rmax = mie.rminmax(idis, r_eff, v_eff, par3, weight2,
                                cutoff)
    # -----
    # LOOP ON WAVELENGTHS
    # -----

    print('Beginning of Mie program')

    for i, wav in enumerate(wavelengths):
        print('Wavelength {:06.3f}'.format(wav))
        nr_mantle = aerosols.nr[i]
        ni_mantle = aerosols.ni[i]
        m1 = (nr_mantle + 1j * ni_mantle)

        nr_core = aerosols.nr_core[i]
        ni_core = aerosols.ni_core[i]
        m2 = (nr_core + 1j * ni_core)
        # Re + i*Im is necessary for BHCOAT

        ratio = aerosols.rcoremant

        scfile_name = specie + '.sc.' + '{:06.3f}'.format(wav)

        # calculation of the scattering matrix
        u, wg, F, miec, nangle = mieshell.scatmat(m1,m2, wav, idis, nsubr,
                                                  ngaur, rmin, rmax, r_eff,
                                                  v_eff, par3, ratio, weight2,
                                                  delta)

        ncoefs = nangle

        # Generalization to 6 matrix elements
        Fshape = np.shape(F)
        F2 = np.zeros((6,Fshape[1]), order='F')
        F2[0,:] = F[0,:]  # F11
        F2[1,:] = F[0,:]  # F22
        F2[2,:] = F[2,:]  # F33
        F2[3,:] = F[2,:]  # F44
        F2[4,:] = F[1,:]  # F12
        F2[5,:] = F[3,:]  # F34

        #expansion of the matrix
        coefs = matrix_expansion(ncoefs, nangle, u, wg, F2)

        # Store the coefficients for each wvl
        supercoefin[i,:,:,:] = coefs
        superncoefin[i] = ncoefs

        qexts[i] = miec[3] # get Qext somewhere
        sexts[i] = miec[1] # get Sext somewhere
        qscas[i] = miec[2] # get Qsca somewhere
        sscas[i] = miec[0] # get Ssca somewhere

        #-----------------------------------------------------------------------
        #        Calculate the asymmetry parameter:
        #-----------------------------------------------------------------------
        if nangle >= 1:
            cosbar = coefs[0, 0, 1] / 3.
        else:
            cosbar = 0.

        asyms[i] = cosbar

        #------------------------------------------------------------------
        #        Write the coefficients to the output file:
        #------------------------------------------------------------------
        if output is True:
            mie.writsc(scfile_name, idis, nsubr, ngaur, coefs, ncoefs, cosbar,
                        miec, wav, nr_mantle, ni_mantle, rmin, rmax, r_eff, v_eff, par3)

        # end of wvl loop


    # Storing the coefficients in the aerosols object
    aerosols.coefs = supercoefin
    aerosols.ncoefs = superncoefin
    aerosols.asym = asyms
    aerosols.sext = sexts
    aerosols.qext = qexts
    aerosols.ssca = sscas
    aerosols.qsca = qscas
    aerosols.ssalb = sscas/sexts

    print('End of Mie program')
    #return u,F, miec


def matrix_expansion(ncoefs, nangle, u, wg, F):
    """ Expansion of a scattering matrix in fourier coefficients.

    Used only as an alias for mie.devel

    Parameters
    ----------
    ncoefs :
        the number of coefficients to be used in the expansion
    nangle:
        number of scattering angles
    u:
        cosine of scattering angles
    wg:
        gaussian weights associated with u
    F:
        scattering matrix

    Returns
    -------
    coefs : array (nmatMAX, nmatMAX, ncoefsMAX)
        the expansion coefficients

    """

    coefs = mie.devel(ncoefs, nangle, u, wg, F)
    return coefs


def read_mie_output(filename, full_output=False, nameout='stuff.dat'):
    """ Read Mie expansion coefficients from file

    Parameters
    ----------
    filename : string
        the name of the file containing the expansion coefficients
    full_output : Bool, optional
        if True, returns phase angles and full scattering matrix. If False,
        returns only phase and -Q/I
        Default is False
    nameout : string
        name of the output file

    Returns
    -------
    theta : array
        scattering angle
    Pl : array
        degree of linear polarization. If full_output is False
    F : array
        a (6,nangles) array with F11, F22, F33, F44, F12, F34 elements of
        the scattering matrix and degree of polarisation
        Only if full_output is True

    Also produces an output file with name nameout

    """

    theta, F = readmie.readmieoutput(filename, nameout)

    Pl = - F[4,:] / F[0,:]

    if full_output is True:
        return theta, F
    else:
        return theta, Pl


