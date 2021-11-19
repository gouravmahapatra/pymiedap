#!/usr/bin/env python
# coding: utf-8

# # Benchmark for PyMieDAP
# This script contains some tests that can tell you whether or not everything is
# working fine in your PyMieDAP install.

import pymiedap.pymiedap as pmd
import matplotlib.pyplot as mpl  # for plotting
import numpy as np 

# Test of the Mie code
# ====================
print( "TEST of MIE CODE")

# We will first check that the Mie code runs fine. For that, we use the
# tabulated output from **De Rooij et al. 1984**.
# We define the same types of aerosols they use.

aerA = pmd.Aerosols(nr=[1.45], ni=[0], r_eff=0.23, v_eff=0.18)
aerB = pmd.Aerosols(nr=[1.44], ni=[0], r_eff=1.05, v_eff=0.07)
aerC = pmd.Aerosols(nr=[1.33], ni=[0], v_eff=0.07, r_eff=2, par3=0.5, psd='7')
aerD = pmd.Aerosols(nr=[1.33], ni=[0], r_eff=2.2, v_eff=0.07)

# then we compute the Mie scattering

coefsA = pmd.mie_code(aerA, [0.55], nsubr=40, ngaur=500) # wvl= 0.55 um
coefsB = pmd.mie_code(aerB, [0.55], nsubr=40, ngaur=500) 
coefsC = pmd.mie_code(aerC, [0.70], nsubr=40, ngaur=500) #wvl = 0.7 um
coefsD = pmd.mie_code(aerD, [0.70], nsubr=40, ngaur=500)

# First thing we can check is that the asymmetry parameter is correct.

print("Result for asymmetry parameter, compared against de Rooij et al. 1984")

text = "Aerosols {}, PyMieDAP {:1.8f}, Theory {:1.8f}, delta {:1.8f}"
theoryA = 0.721
theoryB = 0.71800
theoryC = 0.80420
theoryD = 0.80188

delta=aerA.asym-0.72100
print(text.format('A',float(aerA.asym), float(0.721), float(delta))) 
delta=aerB.asym-theoryB
print(text.format('B',aerB.asym[0], theoryB, float(delta))) 
delta=aerC.asym-theoryC
print(text.format('C',aerC.asym[0], theoryC, float(delta))) 
delta=aerD.asym-theoryD
print(text.format('D',aerD.asym[0], theoryD, float(delta)))

# Then we can check the coefficients on some lines
print(' \n \n')
print('\n Check some expansion coefficients')

coordsx = [0,1,2,3,0,2] #positions of alphas and betas in the coef table
coordsy = [0,1,2,3,1,3]
values_rooij = np.array([[1.0, 0.0, 0.0, 0.9290180850, 0., 0.],
    [0.0880294506, 0.1112167059, 0.1015169161, 0.0832370991, -0.0035893933,
        0.0261461577],
    [0.0022847331, 0.0021817248, 0.0022409582, 0.0025506376, 0.0003370957,
        0.0004637887],
    [0.0000055645, 0.0000068705, 0.0000024694, 0.0000020625, -0.0000007483,
        0.0000035651]])
print('Check of coefficients alpha and beta in table 2 of de Rooij84')

for i,l in enumerate([0, 10, 15, 24]):
    print('Line ',l)
    print('PyMieDAP gives:')
    print(aerA.coefs[0,coordsx,coordsy,l])
    print('Should be')
    print(values_rooij[i,:])
    print('delta',aerA.coefs[0,coordsx,coordsy,l]-values_rooij[i,:])

# And for aerosols B, C and D.
values_rooij = np.array([[1.0, 0.0, 0.0, 0.86462, 0., 0.],
    [5.06996, 5.33409, 5.31811, 5.12295, 0.03999, 0.02898],
    [4.29325, 4.22823, 4.16063, 4.35190, -0.04426, 0.31831],
    [1.48138, 1.66705, 1.50933, 1.42450, -0.04210, 0.24849]])
print('Check of coefficients alpha and beta in table 3 of de Rooij84')

for i,l in enumerate([0, 10, 15, 24]):
    print('Line ',l)
    print('PyMieDAP gives:')
    print(aerB.coefs[0,coordsx,coordsy,l])
    print('Should be')
    print(values_rooij[i,:])
    print('delta',aerB.coefs[0,coordsx,coordsy,l]-values_rooij[i,:])

values_rooij = np.array([[1.0, 0.0, 0.0, 0.9546147, 0., 0.],
    [1.1771468, 1.2847697, 1.2659658, 1.1741245, -0.0223080, 0.0949846],
    [0.4256280, 0.4257144, 0.4231548, 0.4298172, -0.0055369, 0.0426312],
    [0.0000075, 0.0000077, 0.0000072, 0.0000071, 0.0, 0.0000009]])
print('Check of coefficients alpha and beta in table 4 of de Rooij84')

for i,l in enumerate([0, 10, 15, 80]):
    print('Line ',l)
    print('PyMieDAP gives:')
    print(aerC.coefs[0,coordsx,coordsy,l])
    print('Should be')
    print(values_rooij[i,:])
    print('delta',aerC.coefs[0,coordsx,coordsy,l]-values_rooij[i,:])

values_rooij = np.array([[1.0, 0.0, 0.0, 0.91852, 0., 0.],
    [7.20704, 7.35296, 7.25694, 7.13469, -0.02996, 0.09669],
    [7.98540, 8.00874, 7.98885, 7.99970, 0.02916, 0.16475],
    [0.00304, 0.00318, 0.00283, 0.00277, 0.00003, 0.00056]])
print('Check of coefficients alpha and beta in table 3 of de Rooij84')

for i,l in enumerate([0, 10, 15, 80]):
    print('Line ',l)
    print('PyMieDAP gives:')
    print(aerD.coefs[0,coordsx,coordsy,l])
    print('Should be')
    print(values_rooij[i,:])
    print('delta',aerD.coefs[0,coordsx,coordsy,l]-values_rooij[i,:])

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
print(phases)
print(emissions)
print(szas)
pmd.calc_azimuth(phases,szas,emissions, deg=True)


pmd.compute_model(modelA, output_name='modelA', rename=True, nmat=4, nmug=40)
pmd.compute_model(modelB, output_name='modelB', rename=True, nmat=4, nmug=40)



print(modelA.bmsca)
print(modelA.basca)
print(modelB.bmsca)
print(modelB.basca)


# In[14]:

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
fig.suptitle('model A, all delta should be 0')
fig.tight_layout()
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
fig.suptitle('model A, all delta should be 0')
fig.tight_layout()
mpl.show()


print('''
 ## Testing disk integration

 Let's define a model that corresponds to a Lambertian surface with surface albedo 1.0.
 ''')

modelA = pmd.Model()
modelA.wvl_list = [0.7]
del modelA.layers.gastop
del modelA.layers.haze
del modelA.layers.cloud
modelA.layers.gasbelow.rayscat=False
modelA.layers.gasbelow.tau=[0.0]
modelA.surface[0,0] = 1.0


alphas = np.linspace(0,np.pi,80) # phase angles
alphas_deg = np.degrees(alphas)  # in degrees
theta = np.pi - alphas  # scattering angle
P = 2*(np.sin(theta) - theta*np.cos(theta))/(3.*np.pi)  # analytical Lambertian phase function


pmd.planet_integrated([modelA],npix=60, alpha=alphas_deg, force=True)

mpl.plot(modelA.phase, modelA.I[0,:]-P)  #comparing the two models
mpl.xlabel('Phase angle')
mpl.ylabel('Delta I')
mpl.title('Delta between PyMieDAP and lambertian analytic solution')
mpl.show()

