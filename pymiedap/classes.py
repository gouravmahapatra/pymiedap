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

# ---------
# CLASSES DEFINITION
# ---------

class Layer():
    """ This class is intended to describe a layer for the Doubling-Adding
    program. It contains the basic parameters of the model

    Parameters
    ----------
    tau : array
        optical thickness (for each lambda)
    tau_g : array
        optical thickness related to gaseous absorption (for each lambda)
    tau_ray : array
        optical thickness related to rayleigh scattering set by user (for each lambda)
    rayscat: Boolean
        if True, rayleigh scattering is computed. If false, tau_ray is used instead.
    press : float
        pressure at the bottom of the layer [bars]
    aerosols :
        an object containing the properties of a type of aerosols.
        Several of these aerosol objects can coexist in a layer, but they should have
        different names.
    col_dens: float
        particular column density in particles per square micrometers
    """

    def __init__(self, tau=[30], tau_g=[0.], press=30e-3,
                 mix_factor=0., bmsca=[0], bmabs=[0],
                 tau_ray=[0.], rayscat=True, col_dens=0.,
                 basca=[0], baabs=[0]):
        """ Initializes the model object with default values
        aerosols: a subclass to describe the properties of the aerosols
        """
        self.aerosols = Aerosols()
        self.tau = tau
        self.tau_g = tau_g
        self.tau_ray = tau_ray
        self.rayscat = rayscat
        self.press = press
        self.col_dens = col_dens
        self.bmsca = bmsca
        self.bmabs = bmabs
        self.basca = basca
        self.baabs = baabs

    def update_layer(self, nitems):
        """ If the number of working wavelengths is changed, this method
        updates the vectors that depend on wavelength
        """
        self.tau = self.tau[0] * np.ones(nitems)
        self.tau_g = self.tau_g[0] * np.ones(nitems)
        self.tau_ray = self.tau_ray[0] * np.ones(nitems)
        self.bmsca = self.bmsca[0] * np.ones(nitems)
        self.bmabs = self.bmabs[0] * np.ones(nitems)
        self.basca = self.basca[0] * np.ones(nitems)
        self.baabs = self.baabs[0] * np.ones(nitems)

    def mix_aerosols(self):
        """ Mixes the aerosols by combining their scattering matrices"""
        # Checking if there is already a mixed aerosol object
        if hasattr(self,'mixed_aerosols') is True:
            del(self.mixed_aerosols)

        # Preparing some variables
        sum_f = 0
        sum_fsext = np.zeros(len(self.tau)) #sum of sext*f
        sum_fssca = np.zeros(len(self.tau)) #sum of ssca*f
        max_ncoefs = 0
        sum_tau_sca = np.zeros(len(self.tau)) #preparing sums

        sum_coefs = 0
        typ = ''

        # normalise the factors
        for aero_name, aero in vars(self).items():
            if isinstance(aero, Aerosols):
                sum_f += aero.f
                max_ncoefs = np.maximum(max_ncoefs, aero.ncoefs)
                typ = typ + aero.typ

        # retrieve the total column density
        for aero_name, aero in vars(self).items():
            if isinstance(aero, Aerosols):
                aero.f = aero.f / float(sum_f)
                sum_fsext += aero.f * aero.sext
                sum_fssca += aero.f * aero.ssca
        N = self.tau[0] / sum_fsext[0]

        # use it to combine coefs
        for aero_name, aero in vars(self).items():
            if isinstance(aero, Aerosols):
                sum_tau_sca += aero.f * aero.ssca

        for aero_name, aero in vars(self).items():
            if isinstance(aero, Aerosols):
                sum_coefs += (aero.f * aero.coefs)
        mix_coefs = sum_coefs #/ sum_tau_sca[:,np.newaxis, np.newaxis, np.newaxis]

        # filling the mixed object with the result
        self.mixed_aerosols = Aerosols()
        self.mixed_aerosols.coefs = mix_coefs
        self.mixed_aerosols.ncoefs = max_ncoefs
        self.mixed_aerosols.typ = typ
        self.mixed_aerosols.col_dens = N
        self.mixed_aerosols.sext = sum_fsext
        self.mixed_aerosols.ssca = sum_fssca
        self.mixed_aerosols.ssalb = sum_fssca/sum_fsext
        #Warning! Might not be okay!
        #self.mixed_aerosols.nr = aero.nr
        #self.mixed_aerosols.ni = aero.ni
        #self.mixed_aerosols.nr_core= aero.nr_core
        #self.mixed_aerosols.ni_core= aero.ni_core

        print("Aerosols mixed!")


class Layers:
    """ This class contains layer objects describing the atmospheric
    structure desired"""

    def __init__(self):

        self.gastop = Layer(tau=[0.0], press=1e-5)
        self.haze = Layer(press=10e-3, tau=[0.01])
        self.cloud = Layer(press=1.)
        self.gasbelow = Layer(tau=[0.0], press=100)


class Model(object):
    """ This class defines a model planet.

    Parameters
    ----------
    Layers : Layers class object
        contains several layers
    wvl_list : array
        list of wvl to compute. Each time the list is changed (all at
        once) the others vectors within the layers are updated.
    gravity : float
        value of the acceleration of gravity (in m/s^2)
        default 9.81
    asurf : float
        surface albedo for Lambertian surface. If changed, surface[0,0] changes
        too (see below).
        Default is 0
    mma : float
        mean molecular mass of the gas (in atomic mass units)
    dpol : float
        depolarisation factor for the gas medium
        default is 0.09
    rindex_gas : array
        refractive index of the gas (same length as wvl_list). Default is for air.
    Ts : float
        temperature of the star in K
    Rs : float
        Radius of the star in meters
    Dps : float
        distance planet-star in meters
    surface : 4x4 array or string
        an array describing a constant reflection matrix or a filename
        with the Fourier coefficients of a more complicated surface. If surface
        is updated, asurf changes accordingly.

    Returns
    -------
    I, Q, U, V : nd arrays
        Stokes elements
    P : nd array
        degree of linear polarization (-Q/I)
    Pl : nd array
        total degree of linear polarization (sqrt(Q**2+U**2)/I)
    Pt : nd array
        total degree of polarization
    Pu : nd array
        U/I
    Pv : nd array
        V/I
    Pxstd : ndarray
        Standard deviation of quantity x in case of variable pixel masks
    Pxmax1s : ndarray
        Standard deviation of quantity x in case of variable pixel masks
    Pxmax : ndarray
        Max of quantity x in case of variable pixel masks
    Pxmin : ndarray
        Min of quantity x in case of variable pixel masks
    PxminNs : ndarray
        Value of quantity x - N times the standard deviation in case of
        variable pixel masks
    PxmaxNs : ndarray
        Value of quantity x + N times the standard deviation in case of
        variable pixel masks
    fcloud : float
        cloud cover when doing resolved and integrated simulations
    asym : float
        asymmetry of the planet
    picture : 2d array
        image of cloud cover

    The returned attributes depend on the function used on the model object

    """

    def __init__(self, wvl_list=np.array([1.101]), gravity=9.81, dpol=0.09,
                 mma=44, I=[], Q=[], U=[], V=[], P=[], phase=[], asurf=0,
                 Ts=5750, Rs =696342000.,
                 Dps=108208930000.0, rindex_gas=np.array([0]),
                 fcloud=1., asym=0., picture=[], surface=[]):
        """ Generic properties of the model """
        self.gravity = gravity
        self.dpol = dpol
        self.mma = mma
        self.Ts = Ts
        self.Rs = Rs
        self.Dps = Dps

        self.layers = Layers()
        self.geom = Geom()
        self.phase = phase

        self.I = I
        self.Q = Q
        self.U = U
        self.V = V
        self.P = P  # ! warning P is -Q/I
        self._wvl_list = wvl_list
        self.set_rindex_gas()

        self.fcloud = fcloud
        self.asym = asym
        self.picture = picture

        self._surface = np.diag([asurf,0.,0.,0.])
        self._asurf = asurf #surface albedo for lambertian surf

        self.name = ['']

    # auto update of wvl list in layers
    @property
    def wvl_list(self):
        return self._wvl_list

    @wvl_list.setter
    def wvl_list(self, wlist):
        self._wvl_list = wlist
        self.set_rindex_gas()
        for layername,layer in vars(self.layers).items():
            layer.update_layer(len(wlist))
            for aero_name, aero in vars(layer).items():
                if isinstance(aero, Aerosols):
                    aero.update_arrays(len(wlist))

    # mutual update of asurf and surface matrix
    @property
    def asurf(self):
        return self._asurf

    @asurf.setter
    def asurf(self, alb):
        self._asurf = alb
        self._surface[0,0] = self._asurf

    @property
    def surface(self):
        return self._surface

    @surface.setter
    def surface(self, val):
        self._surface = val
        self._asurf = self._surface[0,0]



    def __repr__(self):
        """ Custom display of basic model parameters"""

        strfin = ''
        print('Model:')
        wvl_str = ('**Operating wavelengths:\n'+
                    str(self.wvl_list) + ' microns\n\n')
        planet_str = ("**Planet data:**\n "+
                      "g={:2.2f} m/s^2; ".format(self.gravity) +
                      "surf.alb.={:2.2f}\n \n".format(self.asurf))
        gas_str = ("**Gas data**\n mma={:2.2f} ".format(self.mma) +
                   "dpol={:2.2f}\n".format(self.dpol))
        lays_str = '\n **Layers** \n'
        for layer_name, layer in vars(self.layers).items():
            if hasattr(layer,'mixed_aerosols') is True:
                strout = ('LAYER ' + str(layer_name) +'\n'+
                          ' Type:' + layer.mixed_aerosols.typ +
                          ', P=' + str(layer.press) +
                          ', tau=' + str(layer.tau) +
                          ', tau_gas=' + str(layer.tau_g) + '\n')
                strout += layer.mixed_aerosols.__repr__()
            else:
                for aero_name, aero in vars(layer).items():
                    if isinstance(aero, Aerosols):
                        strout = ('LAYER '+str(layer_name) +'\n'+
                          ' Type:' + aero.typ +
                          ', P=' + str(layer.press) +
                          ', tau=' + str(layer.tau) +
                          ', tau_gas=' + str(layer.tau_g) + '\n')
                        strout += aero.__repr__()
            strfin = strfin + strout
        return wvl_str+planet_str+gas_str+lays_str+strfin


    def model_atm(self,  H=5., z_top=74., k_max=2., r_c = 1.0, v_c = 0.07,
                  n_c=1.42, profile='vertical_profile_ignatiev_38lays.htp',path='./spicavpol/'):
        """New model of Venus'atmosphere with 37 layers for more precise description
        H : scale height
        z_top : height of the top of clouds
        k_max : maximum value of k_ext
        z : height
        dz : thickness of a layer
        z_int : height at the middle of a layer
        z_cut : height below which k_ext is constant
        k_extA : extinction coefficient for z_int >= z_cut
        k_0 : extinction coefficient for z_int <= 47 km
        k_ext : extinction coefficient forall z
        tau_cloud : optical thickness of cloud
        T : temperature
        P : pression
        """

        wvl = self.wvl_list[0]
        nwav = len(self.wvl_list)

        #Reading the file and stocking the data in an array
        tab = np.genfromtxt(path+profile,skip_header=1)

        #Computing the interesting variables
        z = tab[:,0]
        dz = np.diff(z)
        z_int = z[:-1] + dz/2.
        z_cut = z_top - H*np.log(H*k_max)

        k_extA = (1/H) * np.exp(-(z_int-z_top)/H)
        k_0 = np.zeros(len(z_int))
        k_ext = ( (z_int>=z_cut)*k_extA
                 + (z_int<z_cut)*(z_int>47.)*k_max
                 + (z_int<=47.)*k_0 )
						#for z_int >= z_cut        : k_ext = k_extA (mask)
						#for 47 km < z_int < z_cut : k_ext = k_max  (mask)
						#for z_int <= 47 km        : k_ext = 0 km-1 (no cloud below) (mask)
        tau_cloud = k_ext*dz
        T = tab[:,1]
        P = tab[:,2] / 1e3   #pressures in bars!

        #Deleting the old model's layers
        del self.layers.gastop
        del self.layers.haze
        #del self.layers.cloudtop
        del self.layers.cloud
        del self.layers.gasbelow

        #Loading tau_g data and stocking them in a list
        T = []

        for i in np.arange(37):
            t = np.load(path+'dtau_LW/dtau_{:03d}_LW.npz'.format(i))
            T.append(t)

        tab_tau_c = np.zeros((37, nwav))
        tab_tau_g = np.zeros((37, nwav))
        for z,w in enumerate(self.wvl_list):
            # fill table of cloud opacities
            tab_tau_c[:,z] = tau_cloud
            # test: are we in the band?
            if (w<T[0]['lam0'][0]) or (w>T[0]['lam0'][-1]):
                tab_tau_g[:,z] = 0. # if not gaz opacity=0
            else:
                #else give each layer the gaz opacity of the correct wvl
                for l in np.arange(37):
                    index_wvl = np.where(abs(T[l]['lam0']-w) == np.nanmin(abs(T[l]['lam0']-w)))[0]
                    print(index_wvl)
                    tab_tau_g[l,z] = T[l]['tau0'][index_wvl]


        self.layers.layer0 = Layer(tau=tab_tau_c[0,:], tau_g=tab_tau_g[0,:], press=P[0], level=1, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer1 = Layer(tau=tab_tau_c[1,:], tau_g=tab_tau_g[1,:], press=P[1], level=2, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav),typ='C')
        self.layers.layer2 = Layer(tau=tab_tau_c[2,:], tau_g=tab_tau_g[2,:], press=P[2], level=3, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer3 = Layer(tau=tab_tau_c[3,:], tau_g=tab_tau_g[3,:], press=P[3], level=4, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer4 = Layer(tau=tab_tau_c[4,:], tau_g=tab_tau_g[4,:], press=P[4], level=5, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer5 = Layer(tau=tab_tau_c[5,:], tau_g=tab_tau_g[5,:], press=P[5], level=6, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer6 = Layer(tau=tab_tau_c[6,:], tau_g=tab_tau_g[6,:], press=P[6], level=7, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer7 = Layer(tau=tab_tau_c[7,:], tau_g=tab_tau_g[7,:], press=P[7], level=8, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer8 = Layer(tau=tab_tau_c[8,:], tau_g=tab_tau_g[8,:], press=P[8], level=9, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer9 = Layer(tau=tab_tau_c[9,:], tau_g=tab_tau_g[9,:], press=P[9], level=10, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer10 = Layer(tau=tab_tau_c[10,:], tau_g=tab_tau_g[10,:], press=P[10], level=11, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer11 = Layer(tau=tab_tau_c[11,:], tau_g=tab_tau_g[11,:], press=P[11], level=12, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer12 = Layer(tau=tab_tau_c[12,:], tau_g=tab_tau_g[12,:], press=P[12], level=13, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer13 = Layer(tau=tab_tau_c[13,:], tau_g=tab_tau_g[13,:], press=P[13], level=14, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer14 = Layer(tau=tab_tau_c[14,:], tau_g=tab_tau_g[14,:], press=P[14], level=15, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer15 = Layer(tau=tab_tau_c[15,:], tau_g=tab_tau_g[15,:], press=P[15], level=16, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer16 = Layer(tau=tab_tau_c[16,:], tau_g=tab_tau_g[16,:], press=P[16], level=17, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer17 = Layer(tau=tab_tau_c[17,:], tau_g=tab_tau_g[17,:], press=P[17], level=18, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer18 = Layer(tau=tab_tau_c[18,:], tau_g=tab_tau_g[18,:], press=P[18], level=19, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer19 = Layer(tau=tab_tau_c[19,:], tau_g=tab_tau_g[19,:], press=P[19], level=20, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer20 = Layer(tau=tab_tau_c[20,:], tau_g=tab_tau_g[20,:], press=P[20], level=21, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer21 = Layer(tau=tab_tau_c[21,:], tau_g=tab_tau_g[21,:], press=P[21], level=22, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer22 = Layer(tau=tab_tau_c[22,:], tau_g=tab_tau_g[22,:], press=P[22], level=23, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer23 = Layer(tau=tab_tau_c[23,:], tau_g=tab_tau_g[23,:], press=P[23], level=24, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer24 = Layer(tau=tab_tau_c[24,:], tau_g=tab_tau_g[24,:], press=P[24], level=25, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer25 = Layer(tau=tab_tau_c[25,:], tau_g=tab_tau_g[25,:], press=P[25], level=26, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer26 = Layer(tau=tab_tau_c[26,:], tau_g=tab_tau_g[26,:], press=P[26], level=27, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer27 = Layer(tau=tab_tau_c[27,:], tau_g=tab_tau_g[27,:], press=P[27], level=28, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer28 = Layer(tau=tab_tau_c[28,:], tau_g=tab_tau_g[28,:], press=P[28], level=29, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer29 = Layer(tau=tab_tau_c[29,:], tau_g=tab_tau_g[29,:], press=P[29], level=30, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer30 = Layer(tau=tab_tau_c[30,:], tau_g=tab_tau_g[30,:], press=P[30], level=31, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer31 = Layer(tau=tab_tau_c[31,:], tau_g=tab_tau_g[31,:], press=P[31], level=32, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer32 = Layer(tau=tab_tau_c[32,:], tau_g=tab_tau_g[32,:], press=P[32], level=33, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer33 = Layer(tau=tab_tau_c[33,:], tau_g=tab_tau_g[33,:], press=P[33], level=34, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer34 = Layer(tau=tab_tau_c[34,:], tau_g=tab_tau_g[34,:], press=P[34], level=35, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer35 = Layer(tau=tab_tau_c[35,:], tau_g=tab_tau_g[35,:], press=P[35], level=36, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')
        self.layers.layer36 = Layer(tau=tab_tau_c[36,:], tau_g=tab_tau_g[36,:], press=P[36], level=37, r_eff=r_c, v_eff=v_c, nr=n_c*np.ones(nwav), ni=1e-8*np.ones(nwav), typ='C')

        return z_int, k_ext, tau_cloud, tab_tau_g


    def set_taus(self):
        """ modifies the values of the opacity vector according the column
        density and the scattering cross section
        """
        for layername, layer in vars(self.layers).items():
            layer.tau = layer.col_dens * layer.mixed_aerosols.sext
            # should it be sigma_ext or sigma_sca ?

    def set_nrs(self,slope,k):
        """ this function computes the refractive indices as a function of wvl
        given a slope and a constant as in the following equation:
            n = k + slope*lambda
        """

        wvl = np.array(self.wvl_list)
        nrs = k + slope*wvl
        for layername, layer in vars(self.layers).items():
            layer.nr = nrs

    def set_rindex_gas(self, gas='air'):
        """ Computes the refractive index of a gas and sets the refractive
        indices of the model
        INPUTS:
            wvl: wavelength (in microns)
            gas: choose between 'air', 'CO2', 'N2', 'He' and 'H2'
        """
        S = 1./np.array(self.wvl_list)

        if gas=='air':
            # see Ciddor et al. 1996
            rindex = 1 + (0.05792105/(238.0185 - S**2)) + (0.00167917/(57.362 - S**2))
        if gas=='CO2':
            # see Bideau-Mehu el. 1973 (0.1807 -- 1.69 um)
            rindex = 1 + ((6.99100*1e-2/(166.175 - S**2)) + (1.44720*1e-3/(79.609 - S**2)) +
                    (6.42941*1e-5/(56.3064 - S**2)) + (5.21306*1e-5/(46.0196 - S**2))
                    + (1.46847*1e-6/(0.0584738 - S**2)))
        if gas=='H2':
            # Peck and Hung 1977 (0.168 -- 1.6945 um)
            rindex = 1 + (0.0148956/(180.7-S**2)) + (0.0049037/(92 - S**2))
        if gas=='He':
            # Mansfield and Peck 1969
            rindex = 1 + (0.01470091/(423.98 - S**2))
        if gas=='N2':
            # Peck and Khanna 1966 (0.47 -- 2.06 um)
            rindex = 1 + 6.8552*1e-5 + (3.243157*1e-2/(144-S**2))

        self.rindex_gas = rindex


    def summary(self):
        """Writes a summary of model parameters in a .info file with name of
        current model"""

        if self.name==['']:
            print('Model not titled yet!')
        else:
            for w,wvl in enumerate(self.wvl_list):
                filename = self.name[w].split('/')[-1]
                filename = filename.split('.')[:-1]
                filename = '.'.join(filename)
                filename += '.info'
                fich = open(filename,'w')

                strfin = ''
                wvl_str = ('{}:\n'+
                           '**Operating wavelengths:\n'+
                        str(self.wvl_list) + ' microns\n\n')
                planet_str = ("**Planet data:**\n "+
                            "g={:2.2f} m/s^2; ".format(self.gravity) +
                            "surf.alb.={:2.2f}\n \n".format(self.asurf))
                gas_str = ("**Gas data**\n mma={:2.2f} ".format(self.mma) +
                        "dpol={:2.2f}\n".format(self.dpol))
                lays_str = '\n **Layers** \n'
                fich.write(wvl_str)
                fich.write(planet_str)
                fich.write(gas_str)
                fich.write(lays_str)

                for layer_name, layer in vars(self.layers).items():
                    if hasattr(layer,'mixed_aerosols') is True:
                        strout = ('LAYER '+str(layer_name) +'\n'+
                                ' Type:' + layer.mixed_aerosols.typ +
                                ', P=' + str(layer.press) +
                                ', tau=' + str(layer.tau) +
                                ', tau_gas=' + str(layer.tau_g) + '\n')
                        strout += layer.aerosols.__repr__()
                        strfin = strfin + strout
                    else:
                        strout = ('LAYER '+str(layer_name) +'\n'+
                                ' Type:' + layer.aerosols.typ +
                                ', P=' + str(layer.press) +
                                ', tau=' + str(layer.tau) +
                                ', tau_gas=' + str(layer.tau_g) + '\n')
                        strout += layer.aerosols.__repr__()
                        strfin = strfin + strout
                fich.write(strfin)


class Geom():
    """ subclass containing geometrical informations
    Geom : a class containing all the geometrical informations
    sza : solar zenith angle, obtained from data
    emission : emission angle
    azimuth : azimtuh angle
    latitude: latitude
    longitude: longitude
    local_time
    phase : phase angle
    beta : rotation angle btw local meridian and scattering plane (see Hovenier
    1969)
    """

    def __init__(self, sza=[], emission=[],
                 azimuth=[], phase=[], beta=[],
                 latitude=[], longitude=[], local_time=[]):
        self.sza = sza
        self.emission = emission
        self.azimuth = azimuth
        self.phase = phase
        self.beta = beta
        self.latitude = latitude
        self.longitude = longitude
        self.local_time = local_time

    def calc_phase(self):
        """ Computes the phase angle (in degrees) from azimuth, SZA and
        emission angle (also in degrees)"""
        rsza = np.radians(self.sza)
        remission = np.radians(self.emission)
        razimuth = np.radians(self.azimuth)

        cospha = np.cos(rsza)*np.cos(remission) + np.sin(rsza)*np.sin(remission)*np.cos(razimuth)
        rphase = np.arccos(cospha)
        phase = np.degrees(rphase)

        self.phase = phase

    def calc_azimuth(self, deg=True):
        """This method computes the azimuth angle from the geometric
        data. To be used once all geo data have been read and treated.
        INPUTS:
            phase : phase angle
            sza : solar zenith angle
            emission : emission angle
        OUTPUT:
            the azimuthal angle
            If deg==1, the output is given in degrees
            Radians otherwise
        """
        sza = self.sza
        emission = self.emission
        phase = self.phase

        theta = np.radians(sza)
        thetap = np.radians(emission)
        alpha = np.radians(phase)

        t1 = np.cos(alpha) - (np.cos(theta)*np.cos(thetap))
        t2 = np.sin(theta)*np.sin(thetap)
        c_delta_phi = t1/t2
        c_delta_phi[t2<1e-6] = 1. # if denominator is too small, set cos(phi) to 1.
        c_delta_phi[c_delta_phi>1.] = 1.
        c_delta_phi[c_delta_phi<-1.] = -1.

        if deg==1:
            delta_phi = np.degrees(np.arccos(c_delta_phi))
            azimuth = 180. - delta_phi
            azimuth[azimuth<1e-5] = 0.
        else:
            delta_phi = np.arccos(c_delta_phi)
            azimuth = np.pi - delta_phi
            azimuth[azimuth<1e-5] = 0.

        self.azimuth = azimuth

    def calc_beta(self):
        ''' Calculates the rotation angle between the local meridian plane and the
        scattering plane. (see eqs. (7) and (8) in Hovenier 1969)
        Inputs : SZA (deg)
                EMI (deg)
                PHA (deg)
                AZI (deg)
        Returns: BETA (deg)'''
        SZA = self.sza
        EMI = self.emission
        PHA = self.phase
        AZI = self.azimuth

        sgn = np.ones(len(AZI))
        #sgn[AZI<0.] = 1
        #sgn[AZI>0.] = -1

        num = np.cos(np.pi*SZA/180.)-np.cos(np.pi*EMI/180.)*np.cos(np.pi*PHA/180.)
        denom = (sgn*np.sin(np.pi*EMI/180.)*np.sin(np.pi*PHA/180.))
        cb = num/denom
        cb[denom==0] = 0.
        cb[cb>1.] = 1.
        cb[cb<0.] = 0.

        self.beta = np.degrees(np.arccos(cb))


class Aerosols():
    """subclass containing aerosols properties
    nr, ni : real, imaginary parts of refractive index (for each wvl)
    nr_core, ni_core : real, imaginary parts of refractive index for the inner
        core in case of layered spheres (for each wvl)
    rcoremant: ratio between the radius of the outer sphere and the inner core;
        only relevant for layered spheres
    r_eff, v_eff : effective radius and variance (constant with lambda)
    par3: value that can be used for some size distributions
    typ : string indicating type of aerosols ex:'C' for clouds, 'H' for hazes;
        just for user's reference, doesn't change the results in any way
    layered: if True, the calculation is made with Mie scattering for layered
        sphere (see Bohren and Huffmann)
    psd : type of particle size distribution fct 2 stands for modified gamma
    f: is the mix ratio of this type of aerosol in the layer. If f=0.5, half the
    particules are this type.
    * OUTPUT variables:
    qext : extinction coefficients (for each wvl)
    sext : extinction cross-section (for each wvl)
    qsca : scattering coefficients (for each wvl)
    ssca : scattering cross-section (for each wvl)
    coefs: array containing the expansion corefficients from the single
        scattering. For the combined aerosols and for each wvl.

    * SIZE Distributions (value of psd):
        par1 refers to Aerosols.reff
        par2 refers to Aerosols.veff
        par3 refers to Aerosols.par3

        1: TWO PARAMETER GAMMA with alpha (par1) and b (par2) given
        2: TWO PARAMETER GAMMA with reff (par1) and veff (par2) given
        3: Bimodal gamma with equal mode weights
        4: Log normal with rg (par1) and sigma (par2) given
        5: Log normal with reff (par1) and veff (par2) given
        6: Power law with alpha (par1), rmin (par2), rmax(par3)
        7: MODIFIED GAMMA with alpha (par1), rc (par2) and gamma (par3) given
        8: MODIFIED GAMMA with alpha (par1), b (par2) and gamma (par3) given

    """
    def __init__(self, r_eff=1.05, v_eff=0.07, par3=1., rcoremant=0.1, nr=[1.42, 1.41],
                 ni=[1e-8, 1e-8], qext=[0.1, 0.1], sext=[0.1,0.1],
                 ni_core=[1e-8, 1e-8], nr_core=[1.42, 1.41],
                 qsca=[0.1, 0.1], ssca=[0.1,0.1],
                 col_dens=0.6, typ='C', psd='2', layered=False):
        """ Initializes the aerosols object with default values"""
        self.r_eff = r_eff
        self.v_eff = v_eff
        self.par3 = par3
        self.nr = nr
        self.ni = ni
        self.nr_core = nr_core
        self.ni_core = ni_core
        self.rcoremant = rcoremant
        self.typ = typ
        self.psd = psd
        self.qext = qext
        self.sext = sext
        self.qsca = qsca
        self.ssca = ssca
        self.col_dens = col_dens
        self.coefs = []
        self.ncoefs = []
        self.f = 1
        self.layered = layered

    def update_arrays(self, nitems):
        """ adds a new item to each list linked to wvl"""
        self.nr = self.nr[0] * np.ones(nitems)
        self.ni = self.ni[0] * np.ones(nitems)
        self.nr_core = self.nr_core[0] * np.ones(nitems)
        self.ni_core = self.ni_core[0] * np.ones(nitems)
        self.qext = self.qext[0] * np.ones(nitems)
        self.sext = self.sext[0] * np.ones(nitems)
        self.qsca = self.qsca[0] * np.ones(nitems)
        self.ssca = self.ssca[0] * np.ones(nitems)


    def __repr__(self):
        """ Displayed informations when the object is printed """

        if self.layered is True:
            str0 = "Layered spherical particles\n"
            strA = ("nr_core =" + str(self.nr_core) + " )) nr_mantle = "
                    + str(self.nr) + "))\n")
            strB = ("ni_core =" + str(self.ni_core) + " )) ni_mantle = "
                    + str(self.ni) + "))\n")
            strC = ("   R_eff = {:1.3f} \n".format(self.r_eff))
            strD = ("   R core/mantle = {:1.3f} \n".format(self.rcoremant))
            strE = (" Type: " + str(self.typ) + "\n")
            return (str0 + strA + strB + strC + strD + strE)
        if self.layered  is False:
            str0 = "Spherical particles\n"
            strA = ("nr =" + str(self.nr) + "))\n")
            strB = ("ni =" + str(self.ni) + "))\n")
            strC = ("   R_eff = {:1.3f} \n".format(self.r_eff))
            strD = (" Type: " + str(self.typ) + "\n")
            return (str0 + strA + strB + strC + strD)

    def load_coefs(self,filename_list, ncoefsMAX=4001, nmatMAX=4):
        """ Method to load files with expansion coeficients into the Aerosol
        object."""

        nwavels = len(filename_list)
        # Creating array to receive the coefficients
        supercoefin = np.zeros((nwavels,nmatMAX,nmatMAX,ncoefsMAX), order='F')
        superncoefin = np.zeros(nwavels, order='F')

        for i,filename in enumerate(filename_list):
            ncoef, coefs = readmie.file2coefs(filename)

            # Store the coefficients for each wvl
            supercoefin[i,:,:,:] = coefs
            superncoefin[i] = ncoef

        self.coefs = supercoefin
        self.ncoefs = superncoefin



