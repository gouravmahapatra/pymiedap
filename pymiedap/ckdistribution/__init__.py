"""
pymiedap.ckdistribution — Correlated k-distribution subpackage
==============================================================

Computes gaseous absorption cross-sections using the HITRAN database (via HAPI)
and converts them into band-mean absorption optical depths (``bmabs``) suitable
for PyMieDAP's radiative transfer engine.

The pipeline replaces the original Fortran ``salpha.f`` + ``spectrum.f``
workflow with pure-Python HAPI calls that support all HITRAN molecules.

Quickstart
----------
::

    import numpy as np
    from pymiedap.ckdistribution import (
        earth_standard,
        gauss_legendre_points,
        compute_bmabs,
    )

    atm  = earth_standard(nlev=30)
    wav  = np.arange(1.56, 1.65, 0.005)   # µm, CO₂ + CH₄ band
    gp, gw = gauss_legendre_points(10)

    bmabs = compute_bmabs(
        atmosphere=atm,
        molecule='CO2',
        wav=wav,
        gauss_points=gp,
        cache_dir='/tmp/hitran_cache',
        sigma_um=0.005,
    )
    # bmabs.shape == (nlayer, nwav, ngauss)

To compute the full reflected-light spectrum, see
:func:`~.integration.run_ckd_spectrum`.

Submodules
----------
constants       Physical / spectroscopic constants.
gases           HITRAN molecule registry (IDs, molar masses, isotopologue lists).
atmosphere      ``Atmosphere`` class and planet presets (Earth, Venus, Mars).
hitran          HAPI-based line-list fetch and cross-section computation.
kdistribution   CKD logic: slitfunction, ksort, kspec, compute_bmabs.
irf             Instrument response functions (Gaussian, box, SPICAV-IR).
integration     PyMieDAP interface: compute_reflected_spectrum, run_ckd_spectrum.
"""

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

# Constants
from .constants import (
    C2, LOSCH, K_B, N_A,
    T_REF, P_REF, BAR_TO_ATM,
    G_EARTH, G_VENUS, G_MARS,
)

# Molecule registry
from .gases import (
    GASES,
    get_gas,
    molar_mass,
    hitran_id,
    list_gases,
)

# Atmosphere
from .atmosphere import (
    Atmosphere,
    earth_standard,
    venus_ignatiev,
    mars_mgs,
)

# HITRAN / HAPI interface
from .hitran import (
    fetch_lines,
    compute_cross_section,
    compute_cross_section_subL,
    mean_cross_section,
)

# Sub-Lorentzian line profiles and χ-factors
from .lineprofiles import (
    chi_tonkov96,
    chi_perrin89,
    CHI_FACTORS,
    voigt_profile,
    doppler_hwhm,
)

# CKD core
from .kdistribution import (
    wvl2wvn,
    wvn2wvl,
    slitfunction,
    gauss_legendre_points,
    ksort,
    kspec_layer,
    compute_bmabs,
)

# IRF
from .irf import (
    gaussian_irf,
    box_irf,
    no_irf,
    spicav_ir_irf,
    spicav_ir_irf_at_center,
    convolve_spectrum,
)

# PyMieDAP integration
from .integration import (
    compute_reflected_spectrum,
    run_ckd_spectrum,
    save_bmabs,
    load_bmabs,
    write_bmabs_text,
)

__all__ = [
    # constants
    'C2', 'LOSCH', 'K_B', 'N_A', 'T_REF', 'P_REF', 'BAR_TO_ATM',
    'G_EARTH', 'G_VENUS', 'G_MARS',
    # gases
    'GASES', 'get_gas', 'molar_mass', 'hitran_id', 'list_gases',
    # atmosphere
    'Atmosphere', 'earth_standard', 'venus_ignatiev', 'mars_mgs',
    # hitran
    'fetch_lines', 'compute_cross_section', 'compute_cross_section_subL',
    'mean_cross_section',
    # lineprofiles
    'chi_tonkov96', 'chi_perrin89', 'CHI_FACTORS',
    'voigt_profile', 'doppler_hwhm',
    # kdistribution
    'wvl2wvn', 'wvn2wvl', 'slitfunction', 'gauss_legendre_points',
    'ksort', 'kspec_layer', 'compute_bmabs',
    # irf
    'gaussian_irf', 'box_irf', 'no_irf',
    'spicav_ir_irf', 'spicav_ir_irf_at_center', 'convolve_spectrum',
    # integration
    'compute_reflected_spectrum', 'run_ckd_spectrum',
    'save_bmabs', 'load_bmabs', 'write_bmabs_text',
]

__version__ = '0.1.0'
