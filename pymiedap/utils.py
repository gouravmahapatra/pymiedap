# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

import numpy as np

# Small utilities
def calc_azimuth(phase, sza, emission, deg=True):
    """Compute the azimuth angle from the geometric
    data.

    Parameters
    ----------
    phase : float or array
        phase angle
    sza : float or array
        solar zenith angle
    emission : float or array
        emission angle
    deg : Boolean, optional
        If True, the output is given in degrees; in radians otherwise.
        By default True

    Returns
    -------
    azimuth : float or array
        the azimuthal angle
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
    """Calculate the rotation angle between the local meridian plane and the
    scattering plane.

    Parameters
    ----------
    PHA :
        phase angle (deg)
    SZA :
        solar zenith angle (deg)
    EMI :
        emission angle (deg)
    AZI :
        azimuthal angle (deg)

    Returns
    -------
    cb :
        cosine of angle beta (rad)

    """

    sgn = np.ones(len(AZI))
    #sgn[AZI<0.] = 1
    #sgn[AZI>0.] = -1
    #sgn = (-1.* AZI>=np.pi) + (1. * AZI<np.pi)
    #sgn = 1.
    num = np.cos(np.pi*SZA/180.)-np.cos(np.pi*EMI/180.)*np.cos(np.pi*PHA/180.)
    denom = (sgn*np.sin(np.pi*EMI/180.)*np.sin(np.pi*PHA/180.))
    cb = num/denom
    cb[denom==0] = 0.
    cb[cb>1.] = 1.
    cb[cb<-1.] = -1.

    return cb


def sunblackbody(L, Ts=5750):
    """ Compute the blackbody radiance of the Sun at given
    wavelength

    Parameters
    ----------
    L : float or array
        wavelength in metres
    Ts : float, optional
        Surface temperature of the star, in K.
        By default 5750 K

    Returns
    -------
    B : float or array
        Blackbody radiance in W/m2/sr/um

    """

    h = 6.63e-34  # Planck's constant
    c = 2.99e8  # Speed of light
    kb = 1.38e-23  # Boltzmann cst
    #Ts = 5750  # Sun's surf temp.

    A = (2*h*c*c) / (L**5)
    AA = 1. / (np.exp((h*c) / (L*kb*Ts)) - 1)
    B = A * AA  # units : W/m2/sr/m
    B = B * 1e-6  # W/m2/sr/um

    return B


