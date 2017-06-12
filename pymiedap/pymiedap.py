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
import scipy.interpolate as scpint
import scipy.integrate as scpinteg
import scipy.stats as st
import os
import sys
import os.path
import copy
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as mpl
#from lmfit import minimize, Parameters
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

# ---------
# CLASSES DEFINITION
# ---------

class Layer():
    """ This class is intended to describe a layer for the Doubling-Adding
    program. It contains the basic parameters of the model
    tau : optical thickness (for each lambda)
    tau_g : optical thickness related to absorption (for each lambda)
    tau_ray : optical thickness related to rayleigh scattering set by user (for each lambda)
    rayscat: if True, rayleigh scattering is computed. If false, tau_ray is used instead.
    press : pressure at the bottom of the layer
    level : level of the layer with respect to the others from bottom to top (starts at 0)
    aerosols : an object containing the properties of a type of aerosols.
    col_dens: particular column density in particles per square micrometers
    Several of these aerosol objects can coexist in a layer, but they should have
    different names.
    """

    def __init__(self, tau=[30, 30], tau_g=[0.,0.], press=30e-3, psd='2',
                 level=0, mix_factor=0., bmsca=[0, 0], bmabs=[0,0],
                 tau_ray=[0.,0.], rayscat=True,
                 basca=[0,0], baabs=[0,0]):
        """ Initializes the model object with default values
        aerosols: a subclass to describe the properties of the aerosols
        """
        self.aerosols = Aerosols()
        self.tau = tau
        self.tau_g = tau_g
        self.tau_ray = tau_ray
        self.rayscat = rayscat
        self.press = press
        self.level = level
        self.col_dens = 3.2
        self.bmsca = bmsca
        self.bmabs = bmabs
        self.basca = basca
        self.baabs = baabs

    def update_layer(self, nitems):
        """ If the number of working wavelengths is changed, this method
        updates the vectors that depend on wavelength"""
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
        if hasattr(self,'mixed_aerosols')==True:
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
                sum_coefs += (aero.ssca[:,np.newaxis,np.newaxis,np.newaxis] *  aero.f * aero.coefs)
        mix_coefs = sum_coefs / sum_tau_sca[:,np.newaxis, np.newaxis, np.newaxis]

        # filling the mixed object with the result
        self.mixed_aerosols = Aerosols()
        self.mixed_aerosols.coefs = mix_coefs
        self.mixed_aerosols.ncoefs = max_ncoefs
        self.mixed_aerosols.typ = typ
        self.mixed_aerosols.col_dens = N
        self.mixed_aerosols.sext = sum_fsext
        self.mixed_aerosols.ssca = sum_fssca
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

        self.gastop = Layer(tau=[0.0, 0.0], press=1e-5, level=4)
        self.haze = Layer(press=10e-3, tau=[0.01, 0.01], level=3)
        self.cloud = Layer(level=2, press=1., psd='3')
        self.gasbelow = Layer(tau=[0.0, 0.0], press=100, level=1)


class Model(object):
    """ This class defines a model planet.
    Layers : a Layer class object containing several layers
    I, Q, U, V : Stokes vectors
    P : degree of linear polarization (-Q/I)
    wvl_list : list of wvl to compute. Each time the list is changed (all at
        once) the others vectors within the layers are updated.
    gravity: value of the acceleration of gravity (in m/s^2)
    asurf: surface albedo for Lambertian surface
    mma: molecular mass of the gas (in atomic mass units)
    dpol: depolarisation factor for the gas medium
    Ts: temperature of the star in K
    Rs: Radius of the star in meters
    Dps: distance planet-star in meters
    fcloud: cloud cover when doing resolved and integrated simulations
    asym: asymmetry in cloud cover
    picture: image of cloud cover
    surface: an array describing a constant reflection matrix or a filename
        with the Fourier coefficients of a more complicated surface
    """

    def __init__(self, wvl_list=np.array([1.101]), gravity=9.81, dpol=0.09,
                 mma=44, I=[], Q=[], U=[], V=[], P=[], phase=[], asurf=0,
                 Ts=5750, Rs =696342000.,
                 Dps=108208930000.0,
                 fcloud=1., asym=0., picture=[], surface=[]):
        """ Generic properties of the model """
        self.gravity = gravity
        self.dpol = dpol
        self.mma = mma
        self.asurf = asurf #surface albedo for lambertian surf
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

        self.fcloud = fcloud
        self.asym = asym
        self.picture = picture

        self.surface = np.diag([0.,0,0,0])

        self.name = ['']

    @property
    def wvl_list(self):
        return self._wvl_list

    @wvl_list.setter
    def wvl_list(self, wlist):
        self._wvl_list = wlist
        for layername,layer in vars(self.layers).items():
            layer.update_layer(len(wlist))
            for aero_name, aero in vars(layer).items():
                if isinstance(aero, Aerosols):
                    aero.update_arrays(len(wlist))

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
            if hasattr(layer,'mixed_aerosols')==True:
                strout = (str(layer_name) +
                          ' Type:' + layer.mixed_aerosols.typ +
                          ', level=' + str(layer.level) +
                          ', P=' + str(layer.press) +
                          ', tau=' + str(layer.tau) +
                          ', tau_gas=' + str(layer.tau_g) + '\n')
            else:
                strout = (str(layer_name) +
                          ' Type:' + layer.aerosols.typ +
                          ', level=' + str(layer.level) +
                          ', P=' + str(layer.press) +
                          ', tau=' + str(layer.tau) +
                          ', tau_gas=' + str(layer.tau_g) + '\n')
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

    def summary(self):
        """Writes a summary of model parameters in a .info file with name of
        current model"""

        if self.name==['']:
            print('Model not titled yet!')
        else:
            for w,wvl in enumerate(self.wvl_list):
                filename = self.name[w].split('/')[-1]
                filename = filename.split('_')[0]
                filename = '_'.join(filename)
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
                    if hasattr(layer,'mixed_aerosols')==True:
                        strout = (str(layer_name) +
                                ' Type:' + layer.mixed_aerosols.typ +
                                ', level=' + str(layer.level) +
                                ', P=' + str(layer.press) +
                                ', tau=' + str(layer.tau) +
                                ', tau_gas=' + str(layer.tau_g) + '\n')
                    else:
                        strout = (str(layer_name) +
                                ' Type:' + layer.aerosols.typ +
                                ', level=' + str(layer.level) +
                                ', P=' + str(layer.press) +
                                ', tau=' + str(layer.tau) +
                                ', tau_gas=' + str(layer.tau_g) + '\n')
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
        else:
            delta_phi = np.arccos(c_delta_phi)
            azimuth = np.pi - delta_phi

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

        self.beta = np.degrees(np.arccos(cb))


class Aerosols():
    """subclass containing aerosols properties
    nr, ni : real, imaginary parts of refractive index (for each wvl)
    nr_core, ni_core : real, imaginary parts of refractive index for the inner
        core in case of layered spheres (for each wvl)
    rcoremant: ratio between the radius of the outer sphere and the inner core;
        only relevant for layered spheres
    r_eff, v_eff : effective radius and variance (constant with lambda)
    qext : extinction coefficients (for each wvl)
    sext : extinction cross-section (for each wvl)
    qsca : scattering coefficients (for each wvl)
    ssca : scattering cross-section (for each wvl)
    typ : string indicating type of aerosols ex:'C' for clouds, 'H' for hazes;
        just for user's reference, doesn't change the results in any way
    layered: if True, the calculation is made with Mie scattering for layered
        sphere (see Bohren and Huffmann)
    psd : type of particle size distribution fct 2 stands for modified gamma
    f: is the mix ratio of this type of aerosol in the layer. If f=0.5, half the
    particules are this type.
    coefs: array containing the expansion corefficients from the single
        scattering. For the combined aerosols and for each wvl.
    """
    def __init__(self, r_eff=1.05, v_eff=0.07, rcoremant=0.1, nr=[1.42, 1.41],
                 ni=[1e-8, 1e-8], qext=[0.1, 0.1], sext=[0.1,0.1],
                 ni_core=[1e-8, 1e-8], nr_core=[1.42, 1.41],
                 qsca=[0.1, 0.1], ssca=[0.1,0.1],
                 col_dens=0.6, typ='C', psd='2', layered=False):
        """ Initializes the aerosols object with default values"""
        self.r_eff = r_eff
        self.v_eff = v_eff
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

        if self.layered == True:
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



# --------------
# FUNCTIONS
#--------------

# Small utilities
def calc_azimuth(phase, sza, emission, deg=True):
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
    theta = np.radians(sza)
    thetap = np.radians(emission)
    alpha = np.radians(phase)
    t1 = np.cos(alpha) - (np.cos(theta)*np.cos(thetap))
    t2 = np.sin(theta)*np.sin(thetap)

    c_delta_phi = t1/t2
    c_delta_phi[t2<1e-6] = 1. # if denominator is too small, set cos(phi) to 1.
    c_delta_phi[c_delta_phi>1.] = 1.
    c_delta_phi[c_delta_phi<-1.] = -1.

    if deg is True:
        delta_phi = np.degrees(np.arccos(c_delta_phi))
    else:
        delta_phi = np.arccos(c_delta_phi)

    azimuth = delta_phi
    return azimuth


def get_cosbeta(PHA,SZA,EMI,AZI):

    ''' Calculates the rotation angle between the local meridian plane and the
    scattering plane.
    Inputs : SZA (deg)
             EMI (deg)
             PHA (deg)
             AZI (deg)
    Returns: cos BETA (rad) '''

    sgn = np.ones(len(AZI))
    #sgn[AZI<0.] = 1
    #sgn[AZI>0.] = -1
    #sgn = (-1.* AZI>=np.pi) + (1. * AZI<np.pi)
    #sgn = 1.
    num = np.cos(np.pi*SZA/180.)-np.cos(np.pi*EMI/180.)*np.cos(np.pi*PHA/180.)
    denom = (sgn*np.sin(np.pi*EMI/180.)*np.sin(np.pi*PHA/180.))
    cb = num/denom
    cb[denom==0] = 0.

    return cb


def sunblackbody(L, Ts=5750):
    """ Computes the blackbody radiance of the Sun at given
    wavelength L in metres"""

    h = 6.63e-34  # Planck's constant
    c = 2.99e8  # Speed of light
    kb = 1.38e-23  # Boltzmann cst
    #Ts = 5750  # Sun's surf temp.

    A = (2*h*c*c) / (L**5)
    AA = 1. / (np.exp((h*c) / (L*kb*Ts)) - 1)
    B = A * AA  # units : W/m2/sr/m
    B = B * 1e-6  # W/m2/sr/um

    return B


# Main functions
def mie_code(aerosols, wavelengths, output=False, delta=1e-8, cutoff=1e-8, thmin=0, thmax=180,
             step=1, nsubr=50, ngaur=60, nlaysMAX=50, ncoefsMAX=4001,
             nfouMAX=4001, nmuMAX=201, nmatMAX=4):
    """ Takes an input Model object and computes the Mie expansion coefficients
    for the different aerosols types.
    Requires the module_mie module.
    INPUT:
        aerosols : an input aerosol type model containing all the modeling parameters
    OPTIONAL INPUT:
        delta: truncation of the Mie sum
        cutoff: cutoff value for the particle size distribution
        thmin: minimal value of phase angle (in degrees)
        thmax: maximal phase angle
        nsubr: number of subintervals for the distribution
        ngaur: number of Gauss points used in the calculations

        nlaysMAX: maximal number of layers
        ncoefsMAX: max number of coefs
        nfouMAX: max number of Fourier coefs
        nmuMAX:
        nsupMAX:
    OUTPUT:
        coefin : mie expansion coefficients for all layers in the model
        ncoefin : number of non-zero elements in coefin for each layer
        Also generates output file for the different aerosols types if output==True.

    NOTES:
    """

    # ---------------------------
    # Mie calculations parameters
    # ---------------------------
    nsupMAX = nmuMAX * nmatMAX

    par3 = 0.25  # this last parameter is only used for some PSD
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
    specie = aerosols.typ
    idis = aerosols.psd  # index of the particle size dist

    # Getting extrema for radii
    rmin, rmax = mie.rminmax(idis, r_eff, v_eff, par3, weight2, cutoff)
    # -----
    # LOOP ON WAVELENGTHS
    # -----

    print('Beginning of Mie program')

    for i, wav in enumerate(wavelengths):
        print('Wavelength {:06.3f}'.format(wav))
        nr = aerosols.nr[i]
        ni = aerosols.ni[i]
        m = nr - 1j * ni

        scfile_name = specie + '.sc.' + '{:06.3f}'.format(wav)

        # calculation of the scattering matrix
        u, wg, F, miec, nangle = mie.scatmat(m, wav, idis, nsubr, ngaur, rmin,
                                             rmax, r_eff, v_eff, par3, weight2,
                                             delta)

        ncoefs = nangle
        #expansion of the matrix
        coefs = matrix_expansion(ncoefs, nangle, u, wg, F)

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

    print('End of Mie program')


def mie_shell(aerosols, wavelengths, output=False, delta=1e-8, cutoff=1e-8, thmin=0, thmax=180,
              step=1, nsubr=20, ngaur=20, nlaysMAX=50, ncoefsMAX=4001,
              nfouMAX=4001, nmuMAX=201, nmatMAX=4):
    """ Takes an input Model object and computes the Mie expansion coefficients
    for layered spheres.
    Requires the module_mie module.
    INPUT:
        aerosols : an input aerosol type model containing all the modeling parameters
    OPTIONAL INPUT:
        delta: truncation of the Mie sum
        cutoff: cutoff value for the particle size distribution
        thmin: minimal value of phase angle (in degrees)
        thmax: maximal phase angle
        nsubr:
        ngaur: number of Gauss points used in the calculations

        nlaysMAX: maximal number of layers
        ncoefsMAX: max number of coefs
        nfouMAX: max number of Fourier coefs
        nmuMAX:
        nsupMAX:
    OUTPUT:
        coefin : mie expansion coefficients for all layers in the model
        ncoefin : number of non-zero elements in coefin for each layer
        Also generates output file for the different aerosols types if output==1.

    NOTES:
    """

    # ---------------------------
    # Mie calculations parameters
    # ---------------------------
    nsupMAX = nmuMAX * nmatMAX

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
        #expansion of the matrix
        coefs = matrix_expansion(ncoefs, nangle, u, wg, F)

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

    print('End of Mie program')
    return u,F, miec


def matrix_expansion(ncoefs, nangle, u, wg, F):
    """ Expansion of a scattering matrix in fourier coefficients.
    INPUT:
        ncoefs: the number of coefficients to be used in the expansion
        nangle: number of scattering angles
        u: cosine of scattering angles
        wg: gaussian weights associated with u
        F: scattering matrix
    OUTPUT:
        returns the expansion coefficients PRECISE SHAPE!
    """
    coefs = mie.devel(ncoefs, nangle, u, wg, F)
    return coefs


def read_mie_output(filename, full_output=False, nameout='stuff.dat'):
    """ This function reads the output of the Mie code (cf. function mie_code
    above)
    INPUTS:
        filename: the name of the files containing the expansion coefficients
    OUTPUT:
        theta: scattering angle
        Pl: degree of linear polarization
    if full_output=True, the output is:
        theta: scattering angle
        F : a (6,nangles) array with F11, F22, F33, F44, F12, F34 elements of
        the scattering matrix and degree of polarisation

    Also produces an output file 'stuff.dat'
    """

    theta, F = readmie.readmieoutput(filename, nameout)

    Pl = - F[4,:] / F[0,:]

    if full_output is True:
        return theta, F
    else:
        return theta, Pl


def dap_code(model, rename=False, output_name='modelA',
             path_output='./dap_database/', step=10,
             nlaysMAX=50, ncoefsMAX=4001, nfouMAX=4001, ngeosMAX=200,
             nmuMAX=201, nmatMAX=4, nmat=4, nmug=20):
    """ This function launches the DAP code to calculate the supermatrices
    produced by the doubling-adding code. Reads input from a model class.
    Requires the module module_dap.
    INPUT:
        model : a Model object
    OPTIONAL INPUT:
        rename: if True, the user set output file is used for the resulting coefficient files
        output_file: name of the output if rename==True
    OUTPUT:
        produces an output file
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

    pres = np.zeros(nlaysMAX, order='F')  # pressure levels for each layer
    ncoefs = np.zeros(nlaysMAX, order='F')
    wavels = model.wvl_list
    nwavels = len(wavels)
    taus = np.zeros(nlaysMAX, order='F')  # all values of tau
    taus_g = np.zeros(nlaysMAX, order='F')  # all values of tau_g at wvl z
    laylevel = np.zeros(nlaysMAX, order='F')  # level of layers

    # Creating array to receive the coefficients
    coefin = np.zeros((nmatMAX,nmatMAX,ncoefsMAX,nlaysMAX), order='F')
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

        print('Wavelength {:06.3f}'.format(wav))

        # Loop on layers
        l = 0  # layer number
        for layer_name, layer in vars(model.layers).items():

            #  reading information contained in all layers
            # reading the coefficients of each layer for the current wavelength
            taus[l] = layer.tau[z]  #value of tau at wvl z in layer l
            taus_g[l] = layer.tau_g[z]  #value of tau_g at wvl z in layer l
            pres[l] = layer.press
            laylevel[l] = layer.level  # vert position of layer

            # if the layer is transparent, then ignore the aerosols
            # coefficients
            if max(layer.tau)!=0:
                coefin[:,:,:,l] = layer.mixed_aerosols.coefs[z,:,:,:]
                ncoefin[l] = layer.mixed_aerosols.ncoefs[z]

            print('{}.sc.{:06.3f}'.format(layer.mixed_aerosols.typ, wav))
            l = l + 1

        laylevel = laylevel[laylevel!=0]  # getting useful layers positions
        layorder = np.argsort(laylevel)  # putting layers in order

        taus[:l] = taus[layorder]  # putting the taus in right order
        taus_g[:l] = taus_g[layorder]  # putting the taus_g in right order
        pres[:l] = pres[layorder]  # putting the pressures in right order

        coefin[:,:,:,:l] = coefin[:,:,:,layorder]
        ncoefin[:l] = ncoefin[layorder]

        # Calculate the molecular parameters of the atmosphere:
        bmsca, bmabs, coefsm = dap.bmolecules(wav, nlays, pres, dpol, mma, gravity)

        model.coefsm = coefsm
        #-----------------------------------------------------------------
        #     Calculate the aerosol scattering and absorption optical
        #     thicknesses, and the expansion coefficients:
        #-----------------------------------------------------------------
        #basca, baabs, coefs, coefsa, ncoefsa = dap.baerosols(nlays, nmat, taus,
        baabs, basca, coefs, coefsa, ncoefsa = dap.baerosols(nlays, nmat, taus,
                                                             coefin, ncoefin)
        # REM: why are baabs and bsca inverted? should be checked
        model.coefsa = coefsa

        # Storing the effective scattering and absorption opacities
        for layer_name, layer in vars(model.layers).items():
            lev = layer.level

            # force user-define rayleigh opacity
            if layer.rayscat==False:
                bmsca[lev-1] = layer.tau_ray[z]

            layer.bmsca[z] = bmsca[lev-1]
            layer.bmabs[z] = bmabs[lev-1]
            layer.basca[z] = basca[lev-1]
            layer.baabs[z] = baabs[lev-1]

        #---------------------------------------------------------------------
        #     Open the Fourier coefficients output file:
        #---------------------------------------------------------------------
        outputname = 'fou_' + '{:4.3f}'.format(wav) + '.dat'

        #---------------------------------------------------------------
        #     Calculate the combined expansion coefficients
        #---------------------------------------------------------------
       #for i in np.arange(nlays):
       #    ncoefs[i] = max(ncoefsa[i], 2)
       #    for j in np.arange(nmat):
       #        for k in np.arange(nmat):
       #            for m in np.arange(ncoefs[i]):
       #                com = bmsca[i] * coefsm[j,k,m]
       #                coa = basca[i] * coefsa[j,k,m,i]
       #                if ((bmsca[i]+basca[i]) < 1e-10):
       #                    coefs[j,k,m,i]= 0.
       #                else:
       #                    coefs[j,k,m,i]= (com+coa)/(bmsca[i]+basca[i])

        for i in np.arange(nlays):
            ncoefs[i] = max(ncoefsa[i], 2)
            # multiply all coefs by the associated optical thickness
            # in each layer
            com = bmsca[i,np.newaxis,np.newaxis] * coefsm
            coa = basca[i,np.newaxis,np.newaxis,np.newaxis] * coefsa[:,:,:,i]
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
        print('fou_{:4.3f}.dat'.format(wav))
        if rename is True:
            output_file = path_output + output_name + '_{:4.3f}.dat'.format(wav)
            output_file = os.path.normpath(output_file)
            model.name[z] = output_file
            os.rename('fou_{:4.3f}.dat'.format(wav),output_file)
        else:
            output_file = path_output + 'fou_{:4.3f}.dat'.format(wav)
            output_file = os.path.normpath(output_file)
            model.name[z] = output_file
            os.rename('fou_{:4.3f}.dat'.format(wav),output_file)
        print('End of DAP program')


def read_dap_output(phase, sza, emission, filename, beta=None, phi=None,
                    ngeosMAX=100000, nmuMAX=300, nfouMAX=2000, nmatMAX=4):
    """ This function takes a geometry and reads the supermatrices coefficients
    from the DAP code.
    Input:
        phase (deg)
        sza (deg)
        emission (deg)
        filename
    Output:
        Stokes vectors I,Q,U,V normalised with input flux unity (not Pi)
        Assuming normal reflection (i.e. multiply by cos(theta0) if you want real observed flux)
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


def compute_model(atm_model, force=False,
               path_input='./dap_database/', set_taus=False, rename=False,
               output_name='modelA', nmug_mie=20, nmug=20, nsubr=50, nmat=4):
    """
    Function to compute the Fourier files associated with a Model object.
    INPUTS:
        atm_model : a Model object with all the input parameters
    KEYWORDS:
        force : if 0, existing Fourier files are not overwritten; if 1
            existing files are replaced by newer versions
        path_input : path of the fourier DAP files
        set_taus: if True, will set opacities following scattering cross section and column density
        rename: if true, output_name is used
        output_name: custom name radical for the output files of the DAP code
        nmug: number of Gauss points for Mie and DAP calculations
        nmat: number of Stokes elements to compute
    OUTPUT: computes the Fourier files  related to the given model. Also stores
    their names in the model object.
    """

    # Get wvl list
    wvl = atm_model.wvl_list

    # If the model is not yet computed or is forced to
    if atm_model.name[0] == '' or force is True:
        # Execute Mie on all aerosols types on all layers
        for lay, layer in vars(atm_model.layers).items():
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
        force : if 0, existing Fourier files are not overwritten; if 1
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
    """ Calculates values of Q and U after rotation from local meridian
	plane to scattering plane with angle beta"""

    newQ = np.cos(2*beta)*Q + np.sin(2*beta)*U
    newU= -np.sin(2*beta)*Q + np.cos(2*beta)*U

    return newQ, newU


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
                  output_names=['modelA','modelB'], fixed_pattern=True,
                  input_pattern=None, cusp=False, thresh_lat=50., patchy=True,
                  xscale=0.1, yscale=0.01,
                  fclouds=[0.5,0.5], constant_fcloud=False, sscloud=False,
                  sigma_c=10., delta_c=[0.], nmug_mie=20, nmug=20., nsubr=50,
                  nmat=4, pixscaler=1, adaptive_pixels=False):
    """ Function to generate disk-resolved images of a planet according to model
    INPUT:
        atm_model: model to compute
        alpha: phase angles
        npix: number of pixels (total number of pixels will be npix**2)
        force: if True, will force recalculation of model
        set_taus: if True, will set opacities following scattering cross section and column density
        rename: if True, model output files will be renamed
        output_name: output name used for the DAP files if rename=True
        alternate_model: model to be used if inhomogeneous cloud cover. Defines the NON cloudy.
        fixed_pattern: if True, a cloud pattern is generated at start and then
            reused for all phase angle after.
        input_pattern: an existing pattern that can be used as a starter
            (caution: must have size npix*npix)
        cusp: if True, polar cusps are created
        thresh_lat: defines the latitude above which the cusps exist
        patchy: if True, patchy clouds are generated
        fcloud: fraction of the planet to be covered with clouds
        constant_fcloud: if True, the factor fcloud applies not to the whole
            planet but to the lit part of the planet
        sscloud: if True, a subsolar cloud is created
        sigma_c: extend in degrees of the subsolar cloud with respect to the
            SSP. Cloud exists for SZA<sigma_c
        delta_c: offset in degrees the position of the subsolar cloud with
            respect to subsolar point.
        nmug, nmug_mie: number of Gauss point for Mie and DAP codes
        nmat: number of Stokes elements to compute
        nsubr: number of divisions for size dist in Mie calculations
        adaptive_pixels: if True, npix increases with increasing phase angle
            (in sin**2 of alpha/2)
        pixscaler: factor used in combination with adaptive_pixels to set the
            rate of increase in pix number
        xscale: for patchy clouds gives the typical size on x-axis, as a function of npix
        yscale: for patchy clouds gives the typical size on y-axis, as a function of npix
    OUTPUT:
        produces two pdf files for intensity and degree of linear polarization
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
        ppI = PdfPages('I_phase_wvl{:1.3f}.pdf'.format(w))
        ppPl = PdfPages('Pl_phase_wvl{:1.3f}.pdf'.format(w))
        ppU = PdfPages('U_phase_wvl{:1.3f}.pdf'.format(w))
        ppV = PdfPages('V_phase_wvl{:1.3f}.pdf'.format(w))
        ppPt = PdfPages('Pt_phase_wvl{:1.3f}.pdf'.format(w))


        # Loop on phase angle
        # -------------------
        for A,alph in enumerate(alpha):

            if fixed_pattern is False:
                picture_full=None

            if input_pattern!=None:
                picture_full=input_pattern

            if adaptive_pixels is True:
                npix2 = np.ceil(npix * (1 + np.sin(np.radians(alph)/2.)**2))
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

        # PLOT
        font_size=14
        # BUG PIXEL SIZE?
        #pixsize = max(np.diff(lats))

        figsize = 850
        dpi = 90

           ## plot Pl
           #fig = mpl.figure(figsize=(figsize/dpi,figsize/dpi), dpi=dpi)
           #ax = fig.add_subplot(111, aspect=1)
           #ax.set_title('Polarization at phase angle: {:3.2f}'.format(alph))
           #circ = mpl.Circle((0,0),1,color='gray')
           #ax.add_patch(circ)
           #sc = ax.scatter(x, y, c=100*atm_model.P[j,A,:len(x)],lw=0,
           #                marker='s',s=(0.6*figsize/npix2)**2,
           #                cmap=mpl.cm.seismic, zorder=10,)
           #                #vmin=-10, vmax=10)
           #fig.tight_layout(pad=1.2)
           #cb = fig.colorbar(sc,pad=0.02, extend='both')
           #cb.set_label('Degree of linear polarization (%)',size=font_size)
           #ax.set_ylim(-1,1)
           #ax.set_xlim(-1,1)
           ##fig.savefig('Pl_phase_{:01.3f}_{:3.2f}.png'.format(w,alph))
           #ppPl.savefig(fig)
           #mpl.close(fig)

           ## plot U/I
           #fig = mpl.figure(figsize=(figsize/dpi,figsize/dpi), dpi=dpi)
           #ax = fig.add_subplot(111, aspect=1)
           #ax.set_title('Polarization at phase angle: {:3.2f}'.format(alph))
           #circ = mpl.Circle((0,0),1,color='gray')
           #ax.add_patch(circ)
           #sc = ax.scatter(x, y, c=100*(atm_model.U[j,A,:]/atm_model.I[j,A,:]),lw=0,
           #                marker='s',s=(0.6*figsize/npix2)**2,
           #                cmap=mpl.cm.seismic, zorder=10,
           #                vmin=-5, vmax=5)
           #fig.tight_layout(pad=1.2)
           #cb = fig.colorbar(sc,pad=0.02, extend='both')
           #cb.set_label('U/I (%)',size=font_size)
           #ax.set_ylim(-1,1)
           #ax.set_xlim(-1,1)
           ##fig.savefig('U_phase_{:01.3f}_{:3.2f}.png'.format(w,alph))
           #ppU.savefig(fig)
           #mpl.close(fig)

           ## plot I
           #fig = mpl.figure(figsize=(figsize/dpi,figsize/dpi), dpi=dpi)
           #ax = fig.add_subplot(111,aspect=1)
           #ax.set_title('Intensity at phase angle: {:3.2f}'.format(alph))
           #circ = mpl.Circle((0,0),1,color='gray')
           #ax.add_patch(circ)
           #sc = ax.scatter(x, y, c=atm_model.I[j,A,:],lw=0,marker='s',
           #                s=(0.6*figsize/npix2)**2,
           #                cmap=mpl.cm.YlOrRd, zorder=10,)
           #                #vmin=0, vmax=1)
           #fig.tight_layout(pad=1.2)
           #cb = fig.colorbar(sc,pad=0.02, extend='both')
           #cb.set_label('Intensity',size=font_size)
           #ax.set_ylim(-1,1)
           #ax.set_xlim(-1,1)
           ##fig.savefig('I_phase_{:01.3f}_{:3.2f}.png'.format(w,alph))
           #ppI.savefig(fig)
           #mpl.close(fig)

           ## plot V/I
           #fig = mpl.figure(figsize=(figsize/dpi,figsize/dpi), dpi=dpi)
           #ax = fig.add_subplot(111,aspect=1)
           #ax.set_title('Circular polarization V/I at alpha= {:3.2f}'.format(alph))
           #circ = mpl.Circle((0,0),1,color='gray')
           #ax.add_patch(circ)
           #sc = ax.scatter(x, y, c=100*(atm_model.V[j,A,:]/atm_model.I[j,A,:]),lw=0,marker='s',
           #                s=(0.6*figsize/npix2)**2,
           #                cmap=mpl.cm.seismic, zorder=10,
           #                vmin=-0.1, vmax=0.1)
           #fig.tight_layout(pad=1.2)
           #cb = fig.colorbar(sc,pad=0.02, extend='both')
           #cb.set_label('V/I (%)',size=font_size)
           #ax.set_ylim(-1,1)
           #ax.set_xlim(-1,1)
           ##fig.savefig('V_phase_{:01.3f}_{:3.2f}.png'.format(w,alph))
           #ppV.savefig(fig)
           #mpl.close(fig)

           ## plot Ptot
           #Ptot = np.sqrt((atm_model.Q**2+atm_model.U**2+atm_model.V**2)/atm_model.I)
           #fig = mpl.figure(figsize=(figsize/dpi,figsize/dpi), dpi=dpi)
           #ax = fig.add_subplot(111, aspect=1)
           #ax.set_title('Polarization at phase angle: {:3.2f}'.format(alph))
           #circ = mpl.Circle((0,0),1,color='gray')
           #ax.add_patch(circ)
           #sc = ax.scatter(x, y, c=100*Ptot[j,A,:],lw=0,
           #                marker='s',s=(0.6*figsize/npix2)**2,
           #                cmap=mpl.cm.YlOrRd, zorder=10,)
           #                #vmin=0, vmax=5)
           #fig.tight_layout(pad=1.2)
           #cb = fig.colorbar(sc,pad=0.02, extend='both')
           #cb.set_label('Degree of linear polarization (%)',size=font_size)
           #ax.set_ylim(-1,1)
           #ax.set_xlim(-1,1)
           ##fig.savefig('Pt_phase_{:01.3f}_{:3.2f}.png'.format(w,alph))
           #ppPt.savefig(fig)
           #mpl.close(fig)

        ppI.close()
        ppPl.close()
        ppU.close()
        ppV.close()
        ppPt.close()

    mpl.ion()

def plot_pixels(X,Y,Z, title='Polarization', cmap='YlOrRd',vmin=0,vmax=1, font_size=12, npix=20):
    """ Function to nicely plot a resolved planet """
    figsize = 850
    dpi = 90

    fig = mpl.figure(figsize=(figsize/dpi,figsize/dpi), dpi=dpi)
    ax = fig.add_subplot(111, aspect=1)
    ax.set_title(title)
    circ = mpl.Circle((0,0),1,color='gray')
    ax.add_patch(circ)
    sc = ax.scatter(X, Y, c=Z,lw=0, marker='s',
                    s=(0.6*figsize/npix)**2,
                    cmap=cmap, zorder=10,
                    vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=1.2)
    cb = fig.colorbar(sc,pad=0.02, extend='both')
    cb.set_label('Degree of linear polarization (%)',size=font_size)
    ax.set_xlim(-np.nanmax(X),np.nanmax(X))
    ax.set_ylim(-np.nanmax(Y),np.nanmax(Y))
    ax.set_aspect('equal')

    return fig,ax



def planet_integrated(models, alpha=[10], npix=15, force=False, set_taus=False,
                      rename=True, output_names=['modelA','modelB'], fixed_pattern=False,
                      input_pattern=None, cusp=False, thresh_lat=50., full_disk=False,
                      patchy=True, fclouds=[0.5,0.5], constant_fcloud=False,
                      xscale=0.1, yscale=0.01,
                      sscloud=False, sigma_c=10., delta_c=[0.], nmug_mie=20,
                      niter=1, nmug=20., nsubr=50, nmat=4,
                      adaptive_pixels=False):

    """ Function to generate disk-integrated images of a planet according to model
    INPUT:
        models: a list or tuples of models to compute
        alpha: phase angles
        npix: number of pixels (total number of pixels will be npix**2)
        force: if True, will force recalculation of models
        set_taus: if True, will set opacities of each layer in models following
            scattering cross section and column density
        rename: if True, model output files will be renamed
        fixed_pattern: if True, a cloud pattern is generated at start and then
            reused for all phase angles after.
        input_pattern: an existing pattern that can be used as a starter
            (caution: must have size npix*npix)
        cusp: if True, polar cusps are created
        thresh_lat: defines the latitude above which the cusps exist
        patchy: if True, patchy coverage is generated
        fclouds: fractions of the planet to be covered with each model
        constant_fcloud: if True, the factor fcloud applies not to the whole
            planet but to the lit part of the planet
        sscloud: if True, a subsolar cloud is created
        sigma_c: extend in degrees of the subsolar cloud with respect to the
            SSP. Cloud exists for SZA<sigma_c
        delta_c: offset in degrees the position of the subsolar cloud with
            respect to subsolar point.
        output_names: output names used for the DAP files if rename=True
        nmug_mie: number of Gauss point for Mie codes
        nmug: number of Gauss point for DAP codes
        nmat: number of Stokes elements to compute
        nsubr: number of divisions for size dist in Mie calculations
        niter: number of iterations for the coverage
        adaptive_pixels: if True, npix increases with increasing phase angle
            (in sin**2 of alpha/2)
        xscale: for patchy clouds gives the typical size on x-axis, as a function of npix
        yscale: for patchy clouds gives the typical size on y-axis, as a function of npix
    OUTPUT:
        I,Q,U,V: Stokes elements. I(alpha=0) being the geometric albedo
        P: -Q/I
        Pqmin,Pqmax: min and max values of -Q/I
        Plmin,Plmax: min and max values of Pl
        Ptmin,Plmax: min and max values of total polarization
        Imin,Imax: min and max values of intensity
        those parameters being stored in the first model object given as input.
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
    vecasym = np.zeros(len(alpha))

    if len(alpha)!=len(delta_c):
        delta_c = delta_c[0]*np.ones(len(alpha))

    # ===================
    # Loop on phase angle
    # ===================

    for a,alph in enumerate(alpha):

        if fixed_pattern is False:
            picture_full=None

        if input_pattern!=None:
            picture_full=input_pattern

        #Get geom of pixels
        if adaptive_pixels is True:
            npix2 = np.ceil(npix * (1 + np.sin(np.radians(alph)/2.)**2))
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

        for pixtype,model in enumerate(models): #for each pixel type
            for j,w in enumerate(wvl): # and each wvl
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
            # Generate a pixel mask
            picture, mask, picture_full, ncloud, asym = mask_planet(alpha=alph, npix=npix2,
                                                                    fixed_cover=picture_full,
                                                                    cusp=cusp,
                                                                    thresh_lat=thresh_lat,
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
                    vecasym[a] = asym
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
                vecasym[a] = asym
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


def mask_planet(alpha=15, npix=20, cusp=False, thresh_lat=50., patchy=True,
                 ntypes=2, full_disk=False, fclouds=[0.5,0.5], fixed_cover=None,
                 constant_fcloud=False, sscloud=False, sigma_c=10.,
                 delta_c=0., xscale=0.1, yscale=0.01):
    """ Generates a mask that can be used for inhomogeneous planetsi
    INPUTS:
        alpha: phase angle at which the calculation is made
        npix: number of pixels on each axis
        cusp: if True, polar cusps are added
        thresh_lat: latitude threshold above which the cusps extend
        patchy: generates a random patchy cloud cover
        ntypes: types of different pixels for patchy
        fcloud: fraction of the planet that should be covered with clouds
        fixed_cover: if None, a new cloud cover is generated. If a table is
            given, it will be used as the cloud cover
        constant_fcloud: if True, the cloud cover fraction is calculated for
            the given phase angle and not for the whole disk
        sscloud: if True, a subsolar cloud is generated
        sigma_c: extend in degrees of the subsolar cloud. Points between the
            subsolar point and the points with SZA=alpha+sigma are cloudy.
        delta_c: longitudinal offset for the cloud, in degrees
        xscale: for patchy clouds gives the typical size on x-axis, as a function of npix
        yscale: for patchy clouds gives the typical size on y-axis, as a function of npix

    OUTPUT:
        grid_lit: array corresponding to the points of the generated cloud
            cover that are lit
        grid_out: array with the points that actually will be used in an array form
        gird_full: cloud cover of the whole disk
        nb_cloud: fraction of the planet covered with clouds
        asym: asymetry parameter

    """

    # if no specific pattern, assume full_disk
    if patchy==False and cusp==False and sscloud==False:
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

    if full_disk==True:
        # remove outside of disk
        grid_full[xv*xv+yv*yv>1]=np.nan
        # validate pixels lit
        grid_lit[yidx,xidx] = 0. # validate those
        grid_full[yidx,xidx] = 0. # validate those
        # wARNING! arrays have shape (nlines, ncols), hence the grid[y,x]!


    if sscloud==True:
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
    if cusp==True:
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

    if patchy==True:
        #if no fixed cover is wanted
        # n types
        ntypes = len(fclouds)
        fclouds = np.array(fclouds).astype(float) #avoids issues if integers are given
        fclouds = fclouds/sum(fclouds) #renormalization

        if fixed_cover==None:
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
                    moy = (npr.randint(1,npix),npr.randint(1,npix))
                    cov = np.diag([npix*yscale,npix*xscale])
                    x,y = npr.multivariate_normal(moy,cov,50).T
                    # Warning: here x is N/S axis and y E/W axis
                    x = x.astype('int')
                    y = y.astype('int')
                    # if they go beyond the grid, wrap them around
                    x[abs(x)>=npix] = -1
                    y[abs(y)>=npix] = -1
                    # if a pixel is not already taken, give the value of the current type
                    #grid[x,y] = np.where(grid[x,y]==-1, T, grid[x,y])
                    grid_full[x,y] = np.where(grid_full[x,y]==-1, T, grid_full[x,y])
                    grid_lit[x,y] = np.where(grid_lit[x,y]==-1, T, grid_lit[x,y])
                    # wARNING! arrays have shape (nlines, ncols), hence the
                    # grid[x,y]!

                    # lit part of the planet and remove out-of-planet pixels
                    grid_full[xv*xv+yv*yv>1]=np.nan

                    # get current cloud coverage at given phase angle
                    if constant_fcloud==True:
                        cl = np.where(grid_lit>=0)[0].size
                        lit = np.where(~np.isnan(grid_lit))[0].size
                        nb_cloud = float(cl)/(lit)
                    else:
                        cl = np.where(grid_full>=0)[0].size
                        ondisk = np.where(~np.isnan(grid_full))[0].size
                        nb_cloud = float(cl)/(ondisk)

        else:
            # else take existing one
            grid_lit = np.copy(fixed_cover)
            grid_full = np.copy(fixed_cover)
            #grid_lit[xv*xv+yv*yv>1]=np.nan
            grid_full[xv*xv+yv*yv>1]=np.nan

    # get current cloud coverage at given phase angle
    cl = np.where(grid_lit>=0)[0].size
    lit = np.where(~np.isnan(grid_lit))[0].size
    nb_cloud = float(cl)/(lit+1) #+1 to avoid division by 0

    #compute asymmetry factor
    #difference north south hemispheres
    diffgrid = grid_lit[:npix/2.,:] - grid_lit[(npix/2.)-1::-1,:]
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


def surface_check(model,nmug, nmat=4, nmuMAX=201):
    """ A function to check what type of surface is defined by the user, and
    act accordingly for the DAP code"""

    if type(model.surface)==np.ndarray:
        mus, smf, Lfin = fourier_matrix(nmug=nmug, surf_mat=model.surface, nmat=nmat, nmuMAX=nmuMAX)
    elif model.surface==str:
        print('read fourier file for surface')

    return Lfin


def orthographic_projection(center=np.array([0,0]), npix=20):
    """ function to compute the orthographic proj given coordinates"""

    img = Image.open('./earth_contour.png')
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
