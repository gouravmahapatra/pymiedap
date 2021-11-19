#!/usr/bin/env python
# coding: utf-8

# # Benchmark for PyMieDAP
# This script contains some tests that can tell you whether or not everything is
# working fine in your PyMieDAP install.

import pymiedap.pymiedap as pmd
import matplotlib.pyplot as mpl  # for plotting
import numpy as np 

# ## Test of the DAP code
# 
# in their paper, de Haan et al, have tested two atmosphere models:
# * one with an homogeneous atmosphere above a black surface, only with aerosols and opacity=1 (we'll call it modelA)
# * one with a Lambertian surface (A=0.1) and two layers: (modelB)
#     * an upper layer of molecules with opacity 0.1
#     * a lower layer with gas and aerosols mixed such as the molecular opt. thickness is 0.1 and aerosol opt. thickness is 0.4
#     
# Molecular depolarization is 0.0279. The aerosols are those of the type C defined above.

print(' \n \n')
print('Test of DAP CODE')
print('''Test of the DAP code
 
 in their paper, de Haan et al, have tested two atmosphere models:
 * one with an homogeneous atmosphere above a black surface, only with aerosols and opacity=1 (we'll call it modelA)
 * one with a Lambertian surface (A=0.1) and two layers: (modelB)
     * an upper layer of molecules with opacity 0.1
     * a lower layer with gas and aerosols mixed such as the molecular opt. thickness is 0.1 and aerosol opt. thickness is 0.4
     
 Molecular depolarization is 0.0279. The aerosols are those of the type C defined above.
''')

aerC = pmd.Aerosols(nr=[1.33], ni=[0], v_eff=0.07, r_eff=2, par3=0.5, psd='7')

modelB = pmd.Model()
modelB.wvl_list = [0.7]
del modelB.layers.gasbelow
del modelB.layers.haze
modelB.layers.gastop.rayscat=False
modelB.layers.gastop.tau_ray=[0.1]
modelB.layers.cloud.rayscat=False
modelB.layers.cloud.tau = [0.4]
modelB.layers.cloud.tau_ray = [0.1]
modelB.layers.cloud.aerosols = aerC
modelB.dpol=0.0279
modelB.surface[0,0] = 0.1


modelA = pmd.Model()
modelA.wvl_list = [0.7]
del modelA.layers.gastop
del modelA.layers.haze
del modelA.layers.cloud
modelA.layers.gasbelow.rayscat=False
modelA.layers.gasbelow.tau=[1.0]
modelA.layers.gasbelow.aerosols = aerC
modelA.surface[0,0] = 0.0

print('Tests are done for a series of mu, mu0 and phi, as in tables 5-12 from de Haan 87')

mus = np.array([0.1, 0.5, 1.0, 0.1, 0.5, 1.0, 0.1, 0.5, 1.0, 0.1, 0.5, 1.0])
emissions = np.degrees(np.arccos(mus))
mu0s = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
szas = np.degrees(np.arccos(mu0s))
dphi = np.radians([0, 0, 0, 30, 30, 30, 0, 0, 0, 30, 30, 30])
alphas = mus*mu0s  + np.sqrt(1-mus**2)*np.sqrt(1-mu0s**2)*np.cos(dphi) 
phases = np.degrees(np.arccos(alphas))
#print(phases)
#print(emissions)
#print(szas)
pmd.calc_azimuth(phases,szas,emissions, deg=True)


pmd.compute_model(modelA, output_name='modelA', rename=True, nmat=4, nmug=40)
pmd.compute_model(modelB, output_name='modelB', rename=True, nmat=4, nmug=40)


print('Comparing Stokes elements for models A and B, from PyMieDAP and de Haan')
Ia, Qa, Ua, Va = pmd.read_dap_output(phases, szas, emissions, modelA.name[0], beta=np.zeros_like(dphi), phi=np.degrees(dphi))
Ib, Qb, Ub, Vb = pmd.read_dap_output(phases, szas, emissions, modelB.name[0], beta=np.zeros_like(dphi), phi=np.degrees(dphi))

Ia_haan = [1.10269, 0.31943, 0.033033, 0.66414, 0.25209, 0.033033, 2.93214, 0.22054, 0.009287, 0.76910, 0.132828, 0.009287]
Qa_haan = [0.004604, -0.002881, -0.002979, 0.000303, -0.001444, -0.001489, 0.009900, 0.000976, -0.000815, -0.003758, 0.000220, -0.000408]
Ua_haan = [0, 0, 0, -0.002770, -0.004141, -0.002580, 0, 0, 0, 0.003124, -0.000525, -0.000706]
Va_haan = [0, 0, 0, 0.000038, 0.000017, 0., 0, 0, 0, 0.000012, 0.000007, 0.]


fig,ax = mpl.subplots(ncols=2, nrows=2)
ax[0][0].plot(Ia*mu0s*np.pi - Ia_haan)
ax[0][0].set_title('dI')
ax[0][1].plot(Qa*mu0s*np.pi - Qa_haan)
ax[0][1].set_title('dQ')
ax[1][0].plot(Ua*mu0s*np.pi - Ua_haan)
ax[1][0].set_title('dU')
ax[1][1].plot(Va*mu0s*np.pi - Va_haan)
ax[1][1].set_title('dV')
fig.tight_layout()
fig.suptitle('model A, all delta should be 0')
mpl.show()

Ib_haan = [0.53295, 0.20843, 0.093680, 0.41814, 0.18497, 0.093680, 0.52277, 0.106590, 0.026009, 0.27630, 0.083628, 0.026009]
Qb_haan = [-0.028340, -0.036299, -0.024156, -0.000058, -0.019649, -0.012078, 0.011506, -0.005186, -0.014984, 0.034368, 0.003839, -0.007492]
Ub_haan = [0, 0, 0, -0.073105, -0.041401, -0.020920, 0, 0, 0, -0.016042, -0.014492, -0.012976]
Vb_haan = [0, 0, 0, 0.000106, 0.000040, 0, 0, 0, 0, 0.000027, 0.000017, 0]


fig,ax = mpl.subplots(ncols=2, nrows=2)
ax[0][0].plot(Ia*mu0s*np.pi - Ia_haan)
ax[0][0].set_title('dI')
ax[0][1].plot(Qa*mu0s*np.pi - Qa_haan)
ax[0][1].set_title('dQ')
ax[1][0].plot(Ua*mu0s*np.pi - Ua_haan)
ax[1][0].set_title('dU')
ax[1][1].plot(Va*mu0s*np.pi - Va_haan)
ax[1][1].set_title('dV')
fig.tight_layout()
fig.suptitle('model A, all delta should be 0')
mpl.show()


print('End of test of DAP code')
