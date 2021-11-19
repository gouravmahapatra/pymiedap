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

print("Result for asymmetry parameter, compared against de Rooij et al. 1984,
        the difference delta should ideally be negligible \n")

text = "Aerosols {}, PyMieDAP {:1.8f}, Theory {:1.8f}, delta {:1.8f}"
theoryA = 0.721
theoryB = 0.71800
theoryC = 0.80420
theoryD = 0.80188

delta=aerA.asym-theoryA
print(text.format('A',float(aerA.asym), float(0.721), float(delta))) 
delta=aerB.asym-theoryB
print(text.format('B',aerB.asym[0], theoryB, float(delta))) 
delta=aerC.asym-theoryC
print(text.format('C',aerC.asym[0], theoryC, float(delta))) 
delta=aerD.asym-theoryD
print(text.format('D',aerD.asym[0], theoryD, float(delta)))

# Then we can check the coefficients on some lines
print(' \n \n')
print('\n Now we can check some expansion coefficients for the different types
        of aerosols. Comparison is made with de Rooij et al. 1984 ')

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

print("End of test of Mie code")
