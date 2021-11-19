# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

"""
# PYthon MIE DAP (PYMIEDAP)
# This code is used to make computations of Mie scattering along with radiative
# transfer calculations with polarization.
# Dependencies : numpy, matplotlib and scipy
# module_mie, module_mieshell, module_readmie, module_dap, module_geos
#
# Authors : Loic Rossi, Daphne Stam
# Date : 2013 - 2016
# Licence for the Python elements: GNU/GPL & CeCILL
# http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt
# http://www.gnu.org/copyleft/gpl.html
#
# If you use this code, please refer to
# de Rooij et al. 1984, A&A
# de Haan et al. 1987, A&A
# Stam et al. 2006, A&A
"""

# ==============
# IMPORT MODULES
# ==============
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
from PIL import Image

from .classes import Model, Layers, Layer, Geom, Aerosols
from .utils import calc_azimuth, get_cosbeta, sunblackbody
from .mie import mie_code, mie_shell, matrix_expansion, read_mie_output
from .dap import dap_code, read_dap_output

# --------------
# FUNCTIONS
#--------------

def compute_model(atm_model, force=False,
               path_input='./dap_database/', set_taus=False, rename=False,
               output_name='modelA', nmug_mie=20, nmug=20, nsubr=50, nmat=4):
    """
    Compute the Fourier files associated with a Model object.

    Runs mie_code for all Aerosol objects and dap_code on the model

    Parameters
    ----------
    atm_model : Model object
        an input model containing all the modeling parameters
    force : bool, optional
        if True, existing Fourier files are overwritten
        if False, Fourier files are computed only if atm_model has no
        associated file yet.
    path_output : string
        Path of folder where to store the output Fourier files. Folder is
        created if not existing.
    set_taus : bool
        if True, the opacities in atm_model are computed using the column
        densities, instead of the user-defined tau
    rename: Bool
        if True, the user set output file is used for the resulting coefficient files
    output_name : string
        base of filename to be used for the Fourier output files
    nmug : integer, optional
        number of Gauss points used in the DAP calculations. Default 20
    nmug_mie : integer, optional
        number of Gauss points used in the Mie calculations. Default 20
    nsubr : int, optional
        Number of subintervals for Mie scattering.
        Default 50
    nmat: int, optional
        number of Stokes elements to compute. Should be 1, 3 or 4.
        Default 4

    Returns
    -------
    Computes the Fourier files related to the given model. Also stores
    their names in the model object.

    """

    # Get wvl list
    wvl = atm_model.wvl_list

    # If the model is not yet computed or is forced to
    if atm_model.name[0] == '' or force is True:
        # Execute Mie on all aerosols types on all layers
        for lay, layer in vars(atm_model.layers).items():
            print('In layer '+lay+':')
            #If there is already an aerosol mix, we overwrite it
            if hasattr(layer,'mixed_aerosols') is True:
                del(layer.mixed_aerosols)
            for aero_name, aero in vars(layer).items():
                if isinstance(aero, Aerosols):
                    if aero.layered is False:
                        mie_code(aero, atm_model.wvl_list, ngaur=nmug_mie, nsubr=nsubr)
                    else:
                        mie_shell(aero, atm_model.wvl_list, ngaur=nmug_mie, nsubr=nsubr)

            layer.mix_aerosols() #mix aerosols

        # Set the opacities
        if set_taus is True:
            for lay, layer in vars(atm_model.layers).items():
                layer.tau = layer.col_dens * layer.mixed_aerosols.sext

        #execute DAP
        # making sure that the directory we write to is the same we'll read
        # from
        dap_code(atm_model, rename=rename, output_name=output_name, nmug=nmug,
                 nmat=nmat, path_output=path_input)


def read_model(atm_model,data,step=20, force=False,
               path_input='./dap_database/', set_taus=False, rename=False,
               output_name='modelA', nmug_mie=20, nmug=20, nsubr=50, nmat=4):
    """ to read the files associated with a model.
    INPUTS:
        atm_model : a Model object with all the input parameters
        data : a Data object with the observations
    KEYWORDS:
        step : step btw two points used for the fits
        force : if False, existing Fourier files are not overwritten; if True
            existing files are replaced by newer versions
        path_input : path of the fourier DAP files
        set_taus: if True, will set opacities following scattering cross section and column density
        rename: if true, output_name is used
        output_name: custom name radical for the output files of the DAP code
        nmug: number of Gauss points for Mie and DAP calculations
        nmat: number of Stokes elements to compute
    OUTPUT: returns the reflected Stokes vectors for the exact geometry of the
    observation given by data and for the model considered
    """

    # Get wvl list
    wvl = atm_model.wvl_list

    #Get geom
    atm_model.geom.phase = data.phase[::step]
    atm_model.phase = data.phase[::step]
    atm_model.geom.sza = data.geo.sza[::step]
    atm_model.geom.emission = data.geo.emission[::step]
    atm_model.geom.calc_azimuth()
    atm_model.geom.calc_beta()

    # Create table to store output
    n_pts = len(data.phase[::step])
    It = np.zeros((len(wvl),n_pts))
    Qt = np.zeros((len(wvl),n_pts))
    Ut = np.zeros((len(wvl),n_pts))
    Vt = np.zeros((len(wvl),n_pts))

    # compute Fourier file
    compute_model(atm_model, force=force, path_input=path_input,
                  set_taus=set_taus, rename=rename, output_name=output_name,
                  nmug_mie=nmug_mie, nmug=nmug, nsubr=nsubr, nmat=nmat)

    #read files and store result
    for j,w in enumerate(wvl):
        print('Reading {}'.format(atm_model.name[j]))
        I,Q,U,V = read_dap_output(data.phase[::step],data.geo.sza[::step],data.geo.emission[::step],atm_model.name[j],beta=atm_model.geom.beta, phi=atm_model.geom.azimuth)
        It[j,:] = I*np.cos(np.radians(atm_model.geom.sza))
        Qt[j,:] = Q*np.cos(np.radians(atm_model.geom.sza))
        Ut[j,:] = U*np.cos(np.radians(atm_model.geom.sza))
        Vt[j,:] = V*np.cos(np.radians(atm_model.geom.sza))
    atm_model.I = It

    # Compute adjusted radiance
    B = sunblackbody(np.array(atm_model.wvl_list)*1e-6, Ts=atm_model.Ts)
    omegas = np.pi * (atm_model.Rs/atm_model.Dps)**2
    atm_model.I2 = atm_model.I * B[:,np.newaxis] * omegas

    atm_model.Q = Qt
    atm_model.U = Ut
    atm_model.V = Vt
    atm_model.P = -Qt/It


def gen_clouds():
    model = Model()
    model.wvl_list = [1.101]
    model.layers.haze.aerosols.r_eff = 0.25
    model.layers.haze.col_dens = 0.2
    model.layers.haze.aerosols.typ = 'H'
    model.layers.gastop.tau[:] = 0.0
    model.layers.gastop.col_dens = 0.0
    model.layers.gasbelow.tau[:] = 0.0
    model.layers.gasbelow.col_dens = 0.0
    return model


def rotate_stokes(Q,U,beta):
    """ Compute values of Q and U after rotation from one reference
	plane to another plane separated by an angle beta

    Parameters
    ----------
    Q : float or array
        Stokes Q
    U : float or array
        Stokes U
    beta : float or array
        Angle between two reference planes (radians)

    Returns
    -------
    newQ : float or array
        rotated Q
    newU : float or array
        rotated U

    """

    newQ = np.cos(2*beta)*Q + np.sin(2*beta)*U
    newU = -np.sin(2*beta)*Q + np.cos(2*beta)*U

    return newQ, newU

def calc_radial(model):
    """ Function to compute the radial polarization
    Input:
        a model object with a pixel resolved computation
    Output: 
        adds the radial polarization
    """

    # check if x and y are defined
    if not hasattr(model.geom, 'x'):
        raise Exception('This model object does not have pixels defined')

    # compute phi angle
    phi = np.arctan(model.geom.x/model.geom.y)

    #make the rotation and store the radial values
    Qr, Ur = rotate_stokes(model.Q, model.U, phi)
    model.Qr = Qr
    model.Ur = Ur


def binned_average(x,y,xbins, errmean=True, weighted=True, sigmas=1.):
    """ This function computes binned averages
    Input:
        x,y : data to be binned
        xbins : bins in which to average
        errmean : if True, the std returned is the error of the mean
        and not the dispersion
    Output :
        mean and std
    """

    # values of averaged bins
    xmid = xbins[:-1]+0.5*np.diff(xbins)

    # Inverse of sigmas squared
    sigm2 = 1/sigmas**2

    if weighted is True:
        # this is \sum 1/sigma_i^2
        n = np.histogram(x,bins=xbins,weights=sigm2)[0]
        # this is \sum x_i / sigma_i^2
        sy = np.histogram(x,bins=xbins,weights=y*sigm2)[0]
        moy = sy/n
        # this is T^2 = \sum (x_i-moy)**2/ sigma_i^2
        #sigma_intermed = np.histogram(x,bins=xbins,weights=sigm2*(y-moy)**2)[0]
        #sigma : sqrt(T**2/n)
        #sigma = np.sqrt(sigma_intermed/n)
        sigma = np.sqrt(1./n)
    else:
        n = np.histogram(x,bins=xbins)[0]
        sy = np.histogram(x,bins=xbins,weights=y)[0]
        sy2 = np.histogram(x,bins=xbins,weights=y*y)[0]
        moy = sy/n
        sigma = np.sqrt(sy2/n - moy**2)
        if errmean is True:
            sigma = sigma/np.sqrt(n)

    return xmid,moy,sigma


def planet_pixels(models, alpha=[10], npix=15, force=False, set_taus=False, rename=True,
                  output_names=['modelA','modelB'], fixed_pattern=False,
                  input_pattern=None, cusp=False, thresh_lat=50., patchy=True,
                  xscale=0.1, yscale=0.01,
                  bands=False, bands_lats=[-90,90],
                  fclouds=[0.5,0.5], constant_fcloud=False, sscloud=False,
                  sigma_c=10., delta_c=[0.], nmug_mie=20, nmug=20, nsubr=50,
                  nmat=4, pixscaler=1, adaptive_pixels=False):
    """ Generate disk-resolved images of a planet according to model

    Parameters
    ----------
    models : list of Model objects
        models to use in computations
    alpha : array
        phase angles for which to compute
    npix : int
        number of pixels (total number of pixels will be npix**2)
    force : bool, optional
        if True, will force recalculation of model
    set_taus : bool, optional
        if True, will set opacities following scattering cross section and column density
    rename : bool, optional
        if True, model output files will be renamed
    output_names : list of strings
        list of names of the models, used for the name of the Fourier files
    fixed_pattern : bool, optional
        if True, a cloud pattern is generated at start and then
        reused for all phase angle after.
    input_pattern: array, optional
        an existing pattern that can be used as a starter
        (caution: must have size nphase*npix*npix)
    cusp : bool, optional
        if True, polar cusps are created
    thresh_lat : float, optional
        defines the latitude above which the cusps exist
    patchy : bool, optional
        if True, patchy clouds are generated
    fcloud : array, optional
        fraction of the planet to be covered with clouds
        should have same length as models input list
    constant_fcloud : bool, optional
        if True, the factor fcloud applies not to the whole
        planet but to the lit part of the planet
    sscloud : bool, optional
        if True, a subsolar cloud is created
    sigma_c : float, optional
        extend in degrees of the subsolar cloud with respect to the SSP. Cloud
        exists for SZA<sigma_c
    delta_c : float, optional
        offset in degrees the position of the subsolar cloud with
        respect to subsolar point.
    bands : bool, optional
        if True, defines latitudinal bands
    bands_lats : array, optional
        array listing the borders of the bands.
        Ex: [-90, 45, 90] defines two bands
    nmug : int, optional
        number of Gauss point for DAP code
    nmug_mie : int, optional
        number of Gauss point for Mie code
    nmat : int, optional
        number of Stokes elements to compute
    nsubr : int, optional
        number of divisions for size dist in Mie calculations
    adaptive_pixels : bool, optional
        if True, npix increases with increasing phase angle (in sin**2 of
        alpha/2)
    pixscaler : float, optional
        factor used in combination with adaptive_pixels to set the
        rate of increase in pix number. Default 1
    xscale : float, optional
        for patchy clouds gives the typical size on x-axis, as a function of npix
    yscale : float, optional
        for patchy clouds gives the typical size on y-axis, as a function of npix

    Returns
    -------
    I,Q,U,V : arrays
        Stokes elements. I(alpha=0) being the geometric albedo
    P : array
        P is -Q/I
    Pqmin,Pqmax : arrays
        min and max values of -Q/I
    Plmin,Plmax : arrays
        min and max values of Pl
    Ptmin,Plmax : arrays
        min and max values of total polarization
    Imin,Imax : arrays
        min and max values of intensity

    Those parameters being stored in the first model object given as input.

    """

    # If onyl one model is given, assumes a full cover.
    if len(models)==1:
        full_disk=True
        patchy=False
    else:
        full_disk=False

    atm_model = models[0]
    wvl = atm_model.wvl_list
    nwvl = len(atm_model.wvl_list)
    nalpha = len(alpha)
    mpl.ioff()

    # Computing the model atmosphere
    # ------------------------------

    for M, model in enumerate(models):
        compute_model(model, force=force,
                      set_taus=set_taus, rename=rename, output_name=output_names[M],
                      nmug_mie=nmug_mie, nmug=nmug, nsubr=nsubr, nmat=nmat)

    #At start, no specific cloud cover
    picture_full = None
    atm_model.fcloud = np.zeros(len(alpha))
    atm_model.asym = np.zeros(len(alpha))

    # Preparing arrays
    # ------------------
    If = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    Qf = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    Uf = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    Vf = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))

    phaf = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    szaf = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    emif = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    azif = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    betf = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    latf = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    lonf = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    xf = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))
    yf = np.nan*np.zeros((nwvl,nalpha,(2*npix)**2))

    if len(alpha)!=len(delta_c):
        delta_c = delta_c[0]*np.ones(len(alpha))

    # Loop on wvl
    # -------------
    for j,w in enumerate(wvl):

        # Loop on phase angle
        # -------------------
        for A,alph in enumerate(alpha):

            if fixed_pattern is False:
                picture_full=None

            if input_pattern is not None:
                picture_full=input_pattern[A,:,:]

            if adaptive_pixels is True:
                npix2 = np.ceil(npix * (1 + np.sin(np.radians(alph)/2.)**2))
                npix2 = int(npix2)
                print('npix=',npix2)
            else:
                npix2 = npix

            #Get geom
            ngeos, apix, theta0, theta, phi, beta, lats, longs, xs, ys = geos.getgeos(alph, npix2)

            theta0 = theta0[:ngeos]
            theta = theta[:ngeos]
            phi = phi[:ngeos]
            beta = beta[:ngeos]
            lats = lats[:ngeos]
            longs = longs[:ngeos]
            phase = np.ones(ngeos)*alph
            x = xs[:ngeos]
            y = ys[:ngeos]

            phaf[j,A,:ngeos] = alph*np.ones(ngeos)
            szaf[j,A,:ngeos] = theta0[:ngeos]
            emif[j,A,:ngeos] = theta[:ngeos]
            azif[j,A,:ngeos] = phi[:ngeos]
            betf[j,A,:ngeos] = beta[:ngeos]
            latf[j,A,:ngeos] = lats[:ngeos]
            lonf[j,A,:ngeos] = longs[:ngeos]
            xf[j,A,:ngeos] = xs[:ngeos]
            yf[j,A,:ngeos] = ys[:ngeos]


            # Create table to store output
            It = np.zeros((len(wvl),ngeos))
            Qt = np.zeros((len(wvl),ngeos))
            Ut = np.zeros((len(wvl),ngeos))
            Vt = np.zeros((len(wvl),ngeos))

            # call maskd_planet
            if npix!=npix2:
                picture_full = None

            picture, mask, picture_full, ncloud, asym = mask_planet(alpha=alph, npix=npix2,
                                                                    fixed_cover=picture_full,
                                                                    cusp=cusp,
                                                                    thresh_lat=thresh_lat,
                                                                    bands=bands,
                                                                    bands_lats=bands_lats,
                                                                    patchy=patchy,
                                                                    xscale=xscale,
                                                                    yscale=yscale,
                                                                    fclouds=fclouds,
                                                                    constant_fcloud=constant_fcloud,
                                                                    full_disk=full_disk,
                                                                    sscloud=sscloud,
                                                                    sigma_c=sigma_c,
                                                                    delta_c=delta_c[A])

            for pixtype,model in enumerate(models):
                print('Reading {}'.format(model.name[j]))
                phaseB = phase[mask==pixtype]
                theta0B = theta0[mask==pixtype]
                thetaB = theta[mask==pixtype]
                phiB = phi[mask==pixtype]
                betaB = beta[mask==pixtype]

                IB,QB,UB,VB = read_dap_output(phaseB,theta0B,thetaB,model.name[j],phi=phiB, beta=betaB)

                model.picture = np.copy(picture_full)

                It[j,mask==pixtype] = IB*np.cos(np.radians(theta0B))
                Qt[j,mask==pixtype] = QB*np.cos(np.radians(theta0B))
                Ut[j,mask==pixtype] = UB*np.cos(np.radians(theta0B))
                Vt[j,mask==pixtype] = VB*np.cos(np.radians(theta0B))

            If[j,A,:ngeos] = It[j,:]
            Qf[j,A,:ngeos] = Qt[j,:]
            Uf[j,A,:ngeos] = Ut[j,:]
            Vf[j,A,:ngeos] = Vt[j,:]

        atm_model.geom.phase = phaf
        atm_model.phase = phaf
        atm_model.geom.sza = szaf
        atm_model.geom.emission = emif
        atm_model.geom.azimuth = azif
        atm_model.geom.beta = betf
        atm_model.geom.latitude = latf
        atm_model.geom.longitude = lonf
        atm_model.geom.x = xf
        atm_model.geom.y = yf
        atm_model.npix = npix

        atm_model.I = If

        # Compute adjusted radiance
        B = sunblackbody(np.array(atm_model.wvl_list)*1e-6, Ts=atm_model.Ts)
        Rs = 696342000.
        Dvs = 108208930000.0
        omegas = np.pi * (atm_model.Rs/atm_model.Dps)**2
        atm_model.I2 = atm_model.I * B[:,np.newaxis,np.newaxis] * omegas

        atm_model.Q = Qf
        atm_model.U = Uf
        atm_model.V = Vf
        atm_model.P = -Qf/If
        atm_model.Pt = np.sqrt(Qf**2 + Uf**2 + Vf**2)/If
        atm_model.Pl = np.sqrt(Qf**2 + Uf**2)/If
        atm_model.Pv = Vf/If

        # PLOT
        font_size=14
        # BUG PIXEL SIZE?
        #pixsize = max(np.diff(lats))

        figsize = 850
        dpi = 90

    mpl.ion()


def plot_pixels(model, wvl_idx=0, display='grid', stokes='P', phase_idx=0,
                title='Polarization', cmap='YlOrRd',vmin=None,vmax=None, data_scale=1.,
                font_size=12, figsize=8, dpi=100):
    """ Function to nicely plot a resolved planet based on Model object

    Parameters
    ----------
    model : Model object
        a model object already computed and read with pmd.planet_pixels
    wvl_idx : int
        index of the wvl to be plotted
    phase_idx : int, optional
        index in the phase angle array to be plotted. Default is 0
    display : string, optional
        if 'grid', displays the planet as an orthographic projection of a
        sphere at a given phase angle.
        If 'map' displays results as function of latitude/longitude
        Default is 'grid'
    stokes : string, optional
        which Stokes element to plot. Allowed are the names of variables in the
        model object, likely to be in 'P' (-Q/I), 'I', 'Q', 'U', 'Qr', 'Qu',
        'V', 'Pt' (total polarization), 'Pl' total linear pol, 'Pv' for V/I.
        Default is 'P'
    title : string, optional
        title of the figure
        Default is 'Polarization'
    cmap : string, optional
        a matplotlib colormap name, default is 'YlOrRd'
    vmin, vmax: floats or None, optional
        minimum and maximum range of values to plot, default are None
    data_scale :
        multiplier for the displayed quantity. plotted output is data_scale*data
        For example if you want P_l in percents, use data_scale=100.
    font_size : int, optional
        size of the font for the figure, default is 12
    figsize : float, optional
        size of the figure in inches
    dpi : int, optional
        dots per inch, resolution of the figure

    Returns
    -------
    Returns a figure with axes

    """

    npix=model.npix

    fig = mpl.figure(figsize=(figsize,figsize), dpi=dpi)
    ax = fig.add_subplot(111, aspect=1)


    if display == 'grid':
        X = model.geom.x[wvl_idx,phase_idx,:]
        Y = model.geom.y[wvl_idx,phase_idx,:]
        circ = mpl.Circle((0,0),1,color='gray')
        ax.add_patch(circ)
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        axh = bbox.height
        axw = bbox.width
        scalingx = np.ones_like(X)
        scalingy = np.ones_like(Y)
        area = (axh*dpi*axw*dpi)/(1.5*npix)**2
    if display == 'map':
        X = model.geom.longitude[wvl_idx,phase_idx,:]
        Y = model.geom.latitude[wvl_idx,phase_idx,:]
        ax.set_xlim(-90,90)
        ax.set_ylim(-90,90)
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        axh = bbox.height
        axw = bbox.width
        scalingx = 1./np.cos(np.radians(X))
        scalingy = 1./np.cos(np.radians(Y))
        area = (axh*dpi*axw*dpi)/(2*1.5*npix)**2

    #getting the required Stokes element
    Ztmp = getattr(model,stokes)
    #selecting the wvl and phase
    Z = Ztmp[wvl_idx,phase_idx,:]

    Z = data_scale * Z #scaling factor

    ax.set_title(title)
    sc = ax.scatter(X, Y, c=Z,lw=0, marker='s',
                    s=area*scalingx*scalingy,
                    cmap=cmap, zorder=10,
                    vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=1.2)
    cb = fig.colorbar(sc,pad=0.02, extend='both')
    cb.set_label(stokes,size=font_size)
    ax.set_aspect('equal')

    return fig,ax


def planet_integrated(models, alpha=[10], npix=15, force=False, set_taus=False,
                      rename=True, output_names=['modelA','modelB'], fixed_pattern=False,
                      input_pattern=None, cusp=False, thresh_lat=50., full_disk=False,
                      patchy=True, fclouds=[0.5,0.5], constant_fcloud=False,
                      xscale=0.1, yscale=0.01,
                      bands=False, bands_lats=[-90,90],
                      sscloud=False, sigma_c=10., delta_c=[0.], nmug_mie=20,
                      niter=1, nmug=20, nsubr=50, nmat=4,
                      adaptive_pixels=False):

    """ Function to generate disk-integrated images of a planet according to model

    Parameters
    ----------
    models : list of Model objects
        models to use in computations
    alpha : array
        phase angles for which to compute
    npix : int
        number of pixels (total number of pixels will be npix**2)
    force : bool
        if True, will force recalculation of model
    set_taus : bool
        if True, will set opacities following scattering cross section and column density
    rename : bool
        if True, model output files will be renamed
    output_names : list of strings
        list of names of the models, used for the name of the Fourier files
    fixed_pattern : bool
        if True, a cloud pattern is generated at start and then
        reused for all phase angle after.
    input_pattern: array
        an existing pattern that can be used as a starter
        (caution: must have size nphase*npix*npix)
    cusp : bool
        if True, polar cusps are created
    thresh_lat : float
        defines the latitude above which the cusps exist
    patchy : bool
        if True, patchy clouds are generated
    fcloud : array
        fraction of the planet to be covered with clouds
        should have same length as models input list
    constant_fcloud : bool
        if True, the factor fcloud applies not to the whole
        planet but to the lit part of the planet
    sscloud : bool
        if True, a subsolar cloud is created
    sigma_c : float
        extend in degrees of the subsolar cloud with respect to the SSP. Cloud
        exists for SZA<sigma_c
    delta_c : float
        offset in degrees the position of the subsolar cloud with
        respect to subsolar point.
    bands : bool
        if True, defines latitudinal bands
    bands_lats : array
        array listing the borders of the bands.
        Ex: [-90, 45, 90] defines two bands
    nmug : int
        number of Gauss point for DAP code
    nmug_mie : int
        number of Gauss point for Mie code
    nmat : int
        number of Stokes elements to compute
    nsubr : int
        number of divisions for size dist in Mie calculations
    adaptive_pixels : bool
        if True, npix increases with increasing phase angle (in sin**2 of
        alpha/2)
    pixscaler : float, optional
        factor used in combination with adaptive_pixels to set the
        rate of increase in pix number. Default 1
    xscale : float
        for patchy clouds gives the typical size on x-axis, as a function of npix
    yscale : float
        for patchy clouds gives the typical size on y-axis, as a function of npix

    Returns
    -------
    I,Q,U,V : arrays
        Stokes elements. I(alpha=0) being the geometric albedo
    P : array
        P is -Q/I
    Pqmin,Pqmax : arrays
        min and max values of -Q/I
    Plmin,Plmax : arrays
        min and max values of Pl
    Ptmin,Plmax : arrays
        min and max values of total polarization
    Imin,Imax : arrays
        min and max values of intensity

    Those parameters being stored in the first model object given as input.

    """

    # If onyl one model is given, assumes a full cover.
    if len(models)==1:
        full_disk=True
        patchy=False

    # Reading some parameters
    atm_model = models[0]
    wvl = atm_model.wvl_list
    nalpha = len(alpha)
    phases = np.zeros(nalpha)
    nwvl = len(wvl)
    ntypes = len(models)

    Nsteps = nalpha * niter * ntypes

    # Create table to store final output
    It = np.zeros((len(wvl),nalpha))
    Qt = np.zeros((len(wvl),nalpha))
    Ut = np.zeros((len(wvl),nalpha))
    Vt = np.zeros((len(wvl),nalpha))

    # Create tables to store raw variations
    Iall = np.ones((len(wvl),nalpha,niter))
    Qall = np.ones((len(wvl),nalpha,niter))
    Uall = np.ones((len(wvl),nalpha,niter))
    Vall = np.ones((len(wvl),nalpha,niter))

    # loop on models to compute
    for M, model in enumerate(models):
        compute_model(model, force=force,
                      set_taus=set_taus, rename=rename, output_name=output_names[M],
                      nmug_mie=nmug_mie, nmug=nmug, nsubr=nsubr, nmat=nmat)


    #At start, no specific cloud cover
    picture_full = None
    vecfcloud = np.zeros((len(alpha),len(fclouds)))
    vecasym = np.zeros((nalpha,niter))

    if len(alpha)!=len(delta_c):
        delta_c = delta_c[0]*np.ones(len(alpha))

    # ===================
    # Loop on phase angle
    # ===================

    for a,alph in enumerate(alpha):

        if fixed_pattern is False:
            picture_full=None

        if input_pattern is not None:
            picture_full=input_pattern[a,:,:]

        #Get geom of pixels
        if adaptive_pixels is True:
            npix2 = np.ceil(npix * (1 + np.sin(np.radians(alph)/2.)**2))
            npix2 = int(npix2)
            print('npix=',npix2)
        else:
            npix2 = npix

        ngeos, apix, theta0, theta, phi, beta, lats, longs, xs, ys = geos.getgeos(alph, npix2)
        phases[a] = alph

        theta0 = theta0[:ngeos]
        theta = theta[:ngeos]
        phi = phi[:ngeos]
        beta = beta[:ngeos]
        lats = lats[:ngeos]
        longs = longs[:ngeos]
        phase = np.ones(ngeos)*alph
        x = xs[:ngeos]
        y = ys[:ngeos]
        atm_model.geom.phase = phase
        atm_model.phase = phase
        atm_model.geom.sza = theta0
        atm_model.geom.emission = theta
        atm_model.geom.azimuth = phi
        atm_model.geom.beta = beta
        atm_model.geom.latitude = lats
        atm_model.geom.longitude = longs

        # store pixels values for each type of pixel
        Is = np.zeros((len(models),len(wvl),ngeos))
        Qs = np.zeros((len(models),len(wvl),ngeos))
        Us = np.zeros((len(models),len(wvl),ngeos))
        Vs = np.zeros((len(models),len(wvl),ngeos))

        # First mask
        picture, mask, picture_full, ncloud, asym = mask_planet(alpha=alph, npix=npix2,
                                                                fixed_cover=picture_full,
                                                                cusp=cusp,
                                                                thresh_lat=thresh_lat,
                                                                bands=bands,
                                                                bands_lats=bands_lats,
                                                                patchy=patchy,
                                                                xscale=xscale,
                                                                yscale=yscale,
                                                                fclouds=fclouds,
                                                                constant_fcloud=constant_fcloud,
                                                                full_disk=full_disk,
                                                                sscloud=sscloud,
                                                                sigma_c=sigma_c,
                                                                delta_c=delta_c[a])



        for pixtype,model in enumerate(models): #for each pixel type
            for j,w in enumerate(wvl): # and each wvl

                # if only one pattern read only relevant pixels
                if niter==1:
                    if pixtype in np.unique(mask):
                        phaseB = phase[mask==pixtype]
                        theta0B = theta0[mask==pixtype]
                        thetaB = theta[mask==pixtype]
                        phiB = phi[mask==pixtype]
                        betaB = beta[mask==pixtype]
                        print('Reading {}'.format(model.name[j]))
                        I,Q,U,V = read_dap_output(phaseB,theta0B,thetaB,model.name[j],phi=phiB, beta=betaB)
                        Is[pixtype,j,mask==pixtype] = I*np.cos(np.radians(theta0B))
                        Qs[pixtype,j,mask==pixtype] = Q*np.cos(np.radians(theta0B))
                        Us[pixtype,j,mask==pixtype] = U*np.cos(np.radians(theta0B))
                        Vs[pixtype,j,mask==pixtype] = V*np.cos(np.radians(theta0B)) # store current pixel type output
                    else:
                        I,Q,U,V = (0,0,0,0)
                else:
                    # if multiple patterns read all pixels
                    print('Reading {}'.format(model.name[j]))
                    I,Q,U,V = read_dap_output(phase,theta0,theta,model.name[j],phi=phi, beta=beta)
                    Is[pixtype,j,:] = I*np.cos(np.radians(theta0))
                    Qs[pixtype,j,:] = Q*np.cos(np.radians(theta0))
                    Us[pixtype,j,:] = U*np.cos(np.radians(theta0))
                    Vs[pixtype,j,:] = V*np.cos(np.radians(theta0)) # store current pixel type output

        # Create table to store output of all pixels
        Ix = np.zeros((len(wvl),ngeos))
        Qx = np.zeros((len(wvl),ngeos))
        Ux = np.zeros((len(wvl),ngeos))
        Vx = np.zeros((len(wvl),ngeos))

        # ====================
        # Loop on wavelength
        # ====================

        #for j,w in enumerate(wvl):

            # ==================
            # loop on iterations
            # ==================

        for citer in np.arange(niter):

            if citer>0 or npix!=npix2:
                picture_full = None
                # Generate a new pixel mask
                picture, mask, picture_full, ncloud, asym = mask_planet(alpha=alph, npix=npix2,
                                                                        fixed_cover=picture_full,
                                                                        cusp=cusp,
                                                                        thresh_lat=thresh_lat,
                                                                        bands=bands,
                                                                        bands_lats=bands_lats,
                                                                        patchy=patchy,
                                                                        xscale=xscale,
                                                                        yscale=yscale,
                                                                        fclouds=fclouds,
                                                                        constant_fcloud=constant_fcloud,
                                                                        full_disk=full_disk,
                                                                        sscloud=sscloud,
                                                                        sigma_c=sigma_c,
                                                                        delta_c=delta_c[a])


            #===============
            # loop on pixels types
            #===============

            if len(mask) !=0: # if there is some visible pixels
                for pixtype,model in enumerate(models):
                    IB = Is[pixtype,:,:]
                    IB = IB[:,mask==pixtype]
                    QB = Qs[pixtype,:,:]
                    QB = QB[:,mask==pixtype]
                    UB = Us[pixtype,:,:]
                    UB = UB[:,mask==pixtype]
                    VB = Vs[pixtype,:,:]
                    VB = VB[:,mask==pixtype]
                    theta0B = theta0[mask==pixtype]

                    # save some information in the model
                    vecfcloud[a,:] = fclouds
                    vecasym[a,citer] = asym
                    model.picture = np.copy(picture_full)

                    Ix[:,mask==pixtype] = IB
                    Qx[:,mask==pixtype] = QB
                    Ux[:,mask==pixtype] = UB
                    Vx[:,mask==pixtype] = VB

                #===============
                # end of loop on pixels types
                #===============

                # Integrating over planet
                Iall[:,a,citer] = np.nansum(Ix,axis=1)*apix
                Qall[:,a,citer] = np.nansum(Qx,axis=1)*apix
                Uall[:,a,citer] = np.nansum(Ux,axis=1)*apix
                Vall[:,a,citer] = np.nansum(Vx,axis=1)*apix

            else:
                # save some information in the model
                vecfcloud[a,:] = fclouds
                vecasym[a,citer] = asym
                model.picture = np.copy(picture_full)

                # Integrating over planet
                Iall[:,a,citer] = np.nan
                Qall[:,a,citer] = np.nan
                Uall[:,a,citer] = np.nan
                Vall[:,a,citer] = np.nan



            # ==================
            # end of loop on iterations
            # ==================


        # ====================
        # end of wvl loop
        # ====================

        progress = (a+1)*(citer+1)*(pixtype+1) #for the progress bar
        sys.stdout.write('\r{:2f}% done\n'.format(100.*progress/Nsteps))
        sys.stdout.flush()

    # ===================
    # end of loop on phase angle
    # ===================

    # Compute mean value
    It = Iall.mean(axis=2)
    Qt = Qall.mean(axis=2)
    Ut = Uall.mean(axis=2)
    Vt = Vall.mean(axis=2)

    # store result in first input model
    atm_model.I = It
    atm_model.Q = Qt
    atm_model.U = Ut
    atm_model.V = Vt
    atm_model.P = (-Qt/It)
    atm_model.Pl = np.sqrt( (Qt**2+Ut**2)/It**2 )
    atm_model.Pt = np.sqrt( (Qt**2+Ut**2+Vt**2)/It**2 )
    atm_model.Pu = Ut/It
    atm_model.Pv = Vt/It

    # store global results in first input model
    atm_model.Iall = Iall
    atm_model.Qall = Qall
    atm_model.Uall = Uall
    atm_model.Vall = Vall
    atm_model.Pqall = (-Qall/Iall)
    atm_model.Plall = np.sqrt( (Qall**2+Uall**2)/Iall**2 )
    atm_model.Ptall = np.sqrt( (Qall**2+Uall**2+Vall**2)/Iall**2 )
    atm_model.Puall = Uall/Iall
    atm_model.Pvall = Vall/Iall

    # saving dispersion
    atm_model.Pqstd = np.std(atm_model.Pqall, axis=2)
    atm_model.Pustd = np.std(atm_model.Puall, axis=2)
    atm_model.Pvstd = np.std(atm_model.Pvall, axis=2)
    atm_model.Plstd = np.std(atm_model.Plall, axis=2)
    atm_model.Ptstd = np.std(atm_model.Ptall, axis=2)
    atm_model.Istd = np.std(atm_model.Iall, axis=2)
    atm_model.Pqmin1s = atm_model.P-1*atm_model.Pqstd
    atm_model.Pqmin2s = atm_model.P-2*atm_model.Pqstd
    atm_model.Pqmin3s = atm_model.P-3*atm_model.Pqstd
    atm_model.Pumin1s = atm_model.Pu-1*atm_model.Pustd
    atm_model.Pumin2s = atm_model.Pu-2*atm_model.Pustd
    atm_model.Pumin3s = atm_model.Pu-3*atm_model.Pustd
    atm_model.Plmin1s = atm_model.Pl-1*atm_model.Plstd
    atm_model.Plmin2s = atm_model.Pl-2*atm_model.Plstd
    atm_model.Plmin3s = atm_model.Pl-3*atm_model.Plstd
    atm_model.Ptmin1s = atm_model.Pt-1*atm_model.Ptstd
    atm_model.Ptmin2s = atm_model.Pt-2*atm_model.Ptstd
    atm_model.Ptmin3s = atm_model.Pt-3*atm_model.Ptstd
    atm_model.Imin1s = atm_model.I-1*atm_model.Istd
    atm_model.Imin2s = atm_model.I-2*atm_model.Istd
    atm_model.Imin3s = atm_model.I-3*atm_model.Istd

    atm_model.Pqmax1s = atm_model.P+1*atm_model.Pqstd
    atm_model.Pqmax2s = atm_model.P+2*atm_model.Pqstd
    atm_model.Pqmax3s = atm_model.P+3*atm_model.Pqstd
    atm_model.Pumax1s = atm_model.Pu+1*atm_model.Pustd
    atm_model.Pumax2s = atm_model.Pu+2*atm_model.Pustd
    atm_model.Pumax3s = atm_model.Pu+3*atm_model.Pustd
    atm_model.Plmax1s = atm_model.Pl+1*atm_model.Plstd
    atm_model.Plmax2s = atm_model.Pl+2*atm_model.Plstd
    atm_model.Plmax3s = atm_model.Pl+3*atm_model.Plstd
    atm_model.Ptmax1s = atm_model.Pt+1*atm_model.Ptstd
    atm_model.Ptmax2s = atm_model.Pt+2*atm_model.Ptstd
    atm_model.Ptmax3s = atm_model.Pt+3*atm_model.Ptstd
    atm_model.Imax1s = atm_model.I+1*atm_model.Istd
    atm_model.Imax2s = atm_model.I+2*atm_model.Istd
    atm_model.Imax3s = atm_model.I+3*atm_model.Istd

    atm_model.Pqmin = atm_model.Pqall.min(axis=2)
    atm_model.Pqmax = atm_model.Pqall.max(axis=2)
    atm_model.Pumax = atm_model.Puall.max(axis=2)
    atm_model.Pumin = atm_model.Puall.min(axis=2)
    atm_model.Plmin = atm_model.Plall.min(axis=2)
    atm_model.Plmax = atm_model.Plall.max(axis=2)
    atm_model.Ptmin = atm_model.Ptall.min(axis=2)
    atm_model.Ptmax = atm_model.Ptall.max(axis=2)
    atm_model.Imin = atm_model.Iall.min(axis=2)
    atm_model.Imax = atm_model.Iall.max(axis=2)

    # Compute adjusted radiance
    B = sunblackbody(np.array(atm_model.wvl_list)*1e-6,Ts=atm_model.Ts)
    Rs = 696342000.
    Dvs = 108208930000.0
    omegas = np.pi * (Rs/Dvs)**2
    atm_model.I2 = atm_model.I * B[:,np.newaxis] * omegas

    atm_model.geom.phase = phases
    atm_model.phase = phases
    atm_model.fcloud = vecfcloud
    atm_model.asym = vecasym
    atm_model.npix = npix


def mask_planet(alpha=0, npix=20, cusp=False, thresh_lat=50., patchy=True,
                 full_disk=False, fclouds=[0.5,0.5], fixed_cover=None,
                 constant_fcloud=False, sscloud=False, sigma_c=10.,
                 delta_c=0., xscale=0.1, yscale=0.01,
                bands=False, bands_lats=[-90,90]):
    """ Generates a mask that can be used for inhomogeneous planets

    Parameters
    ----------
    alpha : int,optional
        phase angle at which the calculation is made
        default is 0
    npix : int, optional
        number of pixels on each axis
        default is 20
    cusp : bool, optional
        if True, polar cusps are added and pixels above thresh_lat have value 0
        default is False
    thresh_lat : float, optional
        latitude threshold above which the cusps extend
        default is 50.
    patchy : bool
        if True, generates a random patchy cloud cover. First type of patches
        is given value 0, second is 1, etc.
        Default is True
    bands : bool
        If True, generates bands that are defined by their latitudes. First
        band has mask value 0, second is 1, etc.
        default is False
    bands_lats : array
        array with limits of the bands. Should start at -90d and end at 90d. Must be in
        increasing order. Default is [-90,90]
    fclouds : array
        List of fractions of the planet that should be covered with pixels of
        each model
        Default is [0.5, 0.5]
    fixed_cover : None or array
        if None, a new cloud cover is generated. If a table is
        given, it will be used as the cloud cover. Default is None.
    constant_fcloud : Bool
        if True, the cloud cover fraction is calculated for
        the given phase angle and not for the whole disk
        Default is False
    sscloud : bool
        if True, a subsolar cloud is generated
    sigma_c : float
        extend in degrees of the subsolar cloud. Points between the
        subsolar point and the points with SZA=alpha+sigma are given value 0.
    delta_c : float
        longitudinal offset for the subsolar cloud, in degrees
    xscale : float
        for patchy clouds gives the typical size on x-axis, as a function of npix
    yscale : float
        for patchy clouds gives the typical size on y-axis, as a function of npix

    Returns
    -------
    grid_lit : 2d array
        array corresponding to the points of the generated pattern that are lit
    grid_out: 1d array
        flat array with the pixels that are actually lit and on the planet
    gird_full : 2d array
        mask of the whole disk, including the non-lit part of the planet
    nb_cloud : float
        fraction of the planet covered with clouds
    asym : float
        asymetry parameter: amount of pixels of the mask that don't match their
        image by symmetry through the equatorial axis

    """

    # if no specific pattern, assume full_disk
    if patchy is False and cusp is False and sscloud is False and bands is False:
        full_disk=True

    # read the pixel geometries
    ngeos, apix, theta0, theta, phi, beta, lats, longs, xs, ys = geos.getgeos(alpha, npix)

    # prepare grids
    #grid = np.zeros((npix,npix))
    grid_lit = np.zeros((npix,npix))
    grid_full = np.zeros((npix,npix))
    grid_lit[:] = np.nan

    # get angles
    theta0 = theta0[:ngeos]
    theta = theta[:ngeos]
    phi = phi[:ngeos]
    beta = beta[:ngeos]
    lats = lats[:ngeos]
    longs = longs[:ngeos]
    phase = np.ones(ngeos)*alpha
    xs = xs[:ngeos]
    ys = ys[:ngeos]
    xs = xs.round(2)
    ys = ys.round(2)

    # X and Y axis
    step = 2./npix
    X = -1 + 0.5*step + np.arange(0,npix)*step
    Y = -1 + 0.5*step + np.arange(0,npix)*step
    X = X.round(2)
    Y = Y.round(2)
    xv, yv = np.meshgrid(X,Y)
    # remove outside of disk
    grid_full[xv*xv+yv*yv>1]=np.nan
    #grid_lit[xv*xv+yv*yv>1]=np.nan

    # find where correct pixels are on a square grid
    xidx = [np.where(X==item)[0][0] for i, item in enumerate(xs) if item in X]
    yidx = [np.where(Y==item)[0][0] for i, item in enumerate(ys) if item in Y]

    if full_disk is True:
        # remove outside of disk
        grid_full[xv*xv+yv*yv>1]=np.nan
        # validate pixels lit
        grid_lit[yidx,xidx] = 0. # validate those
        grid_full[yidx,xidx] = 0. # validate those
        # wARNING! arrays have shape (nlines, ncols), hence the grid[y,x]!


    if sscloud is True:
        grid_full[:] = 1.
        # validate pixels lit
        grid_lit[yidx,xidx] = 1. # validate those

        # coordinates of subsolar point + offset
        ssidx = np.where((theta0+delta_c)<sigma_c)[0]
        newxs = xs[ssidx]
        newys = ys[ssidx]

        ssxidx = [np.where(X==item)[0][0] for i, item in enumerate(newxs) if item in X]
        ssyidx = [np.where(Y==item)[0][0] for i, item in enumerate(newys) if item in Y]
        grid_full[ssyidx,ssxidx] = 0. # validate those
        grid_lit[ssyidx,ssxidx] = 0. # validate those
        # wARNING! arrays have shape (nlines, ncols), hence the grid[y,x]!

    # If polar cusps
    if cusp is True:
        grid_full[:] = 1.
        # validate pixels lit
        grid_lit[yidx,xidx] = 1. # validate those

        # coordinates of subsolar point + offset
        cuspidx = np.where(abs(lats)>thresh_lat)[0]
        newxs = xs[cuspidx]
        newys = ys[cuspidx]

        cuspxidx = [np.where(X==item)[0][0] for i, item in enumerate(newxs) if item in X]
        cuspyidx = [np.where(Y==item)[0][0] for i, item in enumerate(newys) if item in Y]
        #grid[cuspxidx,cuspyidx] = 0. # validate those
        grid_lit[cuspyidx,cuspxidx] = 0. # validate those
        grid_full[cuspyidx,cuspxidx] = 0. # validate those
        # wARNING! arrays have shape (nlines, ncols), hence the grid[y,x]!

    if bands is True:
        grid_full[:] = 1.
        # validate pixels lit
        grid_lit[yidx,xidx] = 1. # validate those

        for latidx,lat in enumerate(bands_lats[:-1]):
            lower_bound = bands_lats[latidx]
            upper_bound = bands_lats[latidx+1]

            bandidx = np.where((lats>lower_bound)*(lats<upper_bound))[0]
            newxs = xs[bandidx]
            newys = ys[bandidx]
            bandxidx = [np.where(X==item)[0][0] for i, item in enumerate(newxs) if item in X]
            bandyidx = [np.where(Y==item)[0][0] for i, item in enumerate(newys) if item in Y]
            grid_lit[bandyidx,bandxidx] = latidx # validate those
            grid_full[bandyidx,bandxidx] = latidx # validate those

        # wARNING! arrays have shape (nlines, ncols), hence the grid[y,x]!

    if patchy is True:
        #if no fixed cover is wanted
        # n types
        ntypes = len(fclouds)
        fclouds = np.array(fclouds).astype(float) #avoids issues if integers are given
        fclouds = fclouds/sum(fclouds) #renormalization

        if fixed_cover is None:
            #compute a new one
            #starting from a no-type cover (repr. with -1)

            grid_full[:] = -1.
            # validate pixels lit
            grid_lit[yidx,xidx] = -1. # validate those

            #loop on types
            total_fcloud = 0. #total cloud coverage

            for T in np.arange(ntypes):
                total_fcloud += fclouds[T]
                nb_cloud=0

                #fill the planet with pixels
                while nb_cloud<total_fcloud:
                    # generate several multivariate gaussians on the grid
                    moy = (npr.randint(0,npix),npr.randint(0,npix))
                    cov = np.diag([npix*yscale,npix*xscale])
                    x,y = npr.multivariate_normal(moy,cov,50).T
                    # Warning: here x is N/S axis and y E/W axis
                    x = np.round(x)
                    y = np.round(y)
                    x = x.astype('int')
                    y = y.astype('int')
                    # if they go beyond the grid, wrap them around
                    x[x>=npix] = x[x>=npix]-npix
                    y[y>=npix] = y[y>=npix]-npix
                    x[x<-npix] = x[x<-npix]+npix
                    y[y<-npix] = y[y<-npix]+npix
                    # if a pixel is not already taken, give the value of the current type
                    #grid[x,y] = np.where(grid[x,y]==-1, T, grid[x,y])
                    grid_full[x,y] = np.where(grid_full[x,y]==-1, T, grid_full[x,y])
                    grid_lit[x,y] = np.where(grid_lit[x,y]==-1, T, grid_lit[x,y])
                    # wARNING! arrays have shape (nlines, ncols), hence the
                    # grid[x,y]!

                    # lit part of the planet and remove out-of-planet pixels
                    grid_full[xv*xv+yv*yv>1]=np.nan

                    # get current cloud coverage at given phase angle
                    if constant_fcloud is True:
                        cl = np.where(grid_lit>=0)[0].size
                        lit = np.where(~np.isnan(grid_lit))[0].size
                        nb_cloud = float(cl)/(lit)
                    else:
                        cl = np.where(grid_full>=0)[0].size
                        ondisk = np.where(~np.isnan(grid_full))[0].size
                        nb_cloud = float(cl)/(ondisk)

        else:
            # else take existing one
            grid_lit = np.zeros((npix,npix))
            grid_full = np.zeros((npix,npix))
            grid_full[:,:] = np.nan
            grid_lit[:,:] = np.nan
            grid_lit[yidx,xidx] = fixed_cover[yidx,xidx]
            grid_full[yidx,xidx] = fixed_cover[yidx,xidx]
            grid_lit[xv*xv+yv*yv>1]=np.nan
            grid_full[xv*xv+yv*yv>1]=np.nan

    # get current cloud coverage at given phase angle
    cl = np.where(grid_lit>=0)[0].size
    lit = np.where(~np.isnan(grid_lit))[0].size
    nb_cloud = float(cl)/(lit+1) #+1 to avoid division by 0

    #compute asymmetry factor
    #difference north south hemispheres
    N = int(max(npix/2, (npix+1)/2)) #safety measure to avoid issues with odd npix
    n = int(min(npix/2, (npix+1)/2))
    diffgrid = grid_lit[:N,:] - (grid_lit[n:,:])[::-1,:]
    # where is the difference in cloud cover?
    asympix = np.where(diffgrid[~np.isnan(diffgrid)]!=0)[0]
    # count the asymetry
    if np.size(diffgrid[~np.isnan(diffgrid)])==0:
        asym=0
    else:
        asym = np.size(asympix)/np.float(np.size(diffgrid[~np.isnan(diffgrid)]))

    # Flatten the grid with only lit points
    grid_out = grid_lit[~np.isnan(grid_lit)]
    grid_out = grid_out.flatten()

    return grid_lit, grid_out, grid_full, nb_cloud, asym




def orthographic_projection(center=np.array([0,0]), npix=20, input_img='./earth_contour.png'):
    """ function to compute the orthographic proj given coordinates"""

    img = Image.open(input_img)
    img.thumbnail((4*npix,2*npix))
    maps = np.array(img)
    maps = maps/np.max(maps)
    maps[maps>0.5] = 1.0
    maps[maps<0.5] = 0.0
    maps = maps.astype('float')
    nlat, nlon = np.shape(maps)
    latitudes = np.linspace(-90,90,nlat)
    longitudes = np.linspace(-180,180,nlon)
    maps2 = np.copy(maps)
    lats = np.radians(latitudes)
    lons = np.radians(longitudes)
    center = np.radians(center)
    center[0] = -center[0]
    x = np.cos(lats[:,np.newaxis])*np.sin(lons[np.newaxis,:]-center[1])
    y = (np.cos(center[0])*np.sin(lats[:,np.newaxis]) -
         np.sin(center[0])*np.cos(lats[:,np.newaxis])*np.cos(lons[np.newaxis,:]-center[1]))
    cosc = (np.sin(center[0])*np.sin(lats[:,np.newaxis]) +
            np.cos(center[0])*np.cos(lats[:,np.newaxis])*np.cos(lons[np.newaxis,:]-center[1]))
    x[cosc<0.] = np.nan
    y[cosc<0.] = np.nan
    maps[cosc<0.] = np.nan

    # X and Y axis
    step = 2./npix
    X = -1 + 0.5*step + np.arange(0,npix)*step
    Y = -1 + 0.5*step + np.arange(0,npix)*step
    X = X.round(2)
    Y = Y.round(2)

    # prepare grids
    grid_lit = np.zeros((npix,npix))
    grid_full = np.zeros((npix,npix))
    grid_lit[:] = np.nan
    grid_full[:] = np.nan

    # find where correct pixels are on a square grid
    xidx = [np.argmin(abs(X-item)) for i, item in enumerate(x.flatten()) if ~np.isnan(item)]
    yidx = [np.argmin(abs(Y-item)) for i, item in enumerate(y.flatten()) if ~np.isnan(item)]

    for i,ii in enumerate(x.flatten()):
        if ~np.isnan(ii):
            xtmp = np.argmin(abs(X-ii))
            jj = y.flatten()[i]
            ytmp = np.argmin(abs(Y-jj))
            grid_full[xtmp,ytmp] = maps.flatten()[i]

    return x,y,maps2, xidx, yidx, grid_full


# -----------
# MAIN
# -----------

if __name__ == '__main__':

    print('ok')
