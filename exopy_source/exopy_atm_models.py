from .. import pymiedap as _pmd
import numpy as _np
import os as _os

def Earth(alb=0.3, g=9.81, pc=0.7, reff=1.0, nr=1.33, wvlmin=0.4, wvlmax=0.41, dwvl=0.2, path_input='./dap_database/'):

   # Cloudy model

    cloudy = _pmd.Model()
    cloudy.mma = 29
    cloudy.gravity = g
    cloudy.dpol = 0.03
    del(cloudy.layers.haze)
    cloudy.wvl_list = _np.arange(wvlmin, wvlmax, dwvl)
    nwvl = len(cloudy.wvl_list)

    cloudy.layers.gasbelow.level = 1
    cloudy.layers.gasbelow.press = 1.

    cloudy.layers.cloud.level = 2
    cloudy.layers.cloud.press = pc+0.1
    cloudy.layers.cloud.tau = 6*_np.ones(nwvl)
    cloudy.layers.cloud.aerosols.nr = nr*_np.ones(nwvl)
    cloudy.layers.cloud.aerosols.r_eff = reff
    cloudy.layers.cloud.aerosols.v_eff = 0.10

    cloudy.layers.gastop.press = pc
    cloudy.layers.gastop.level = 3
    cloudy.surface[0,0] = alb
    
    cloudy.tag = 'cloudy_{:03.1f}_{:3.2f}_{:2.1f}'.format(g,alb,reff)

    check = _np.ones_like(cloudy.wvl_list)
    for i,wvl in enumerate(cloudy.wvl_list):
	    filename = cloudy.tag + '_{:4.3f}.dat'.format(wvl)

	    if _os.path.isfile(filename) and _os.access(filename, _os.R_OK):
		check[i] = 1
	    else:
		check[i] = 0

    if check.prod() ==1:
	    print('Files exist!')
	    cloudy.name = ['']*_np.size(cloudy.wvl_list)
	    for i, wvl in enumerate(cloudy.wvl_list):
	    	cloudy.name[i] = _os.path.normpath(cloudy.tag + '_{:4.3f}.dat'.format(wvl))


    # Clear sky model
    clear = _pmd.Model()
    clear.mma = 29 #for air
    clear.gravity = g
    clear.dpol = 0.03
    del(clear.layers.haze)
    clear.wvl_list = _np.arange(wvlmin, wvlmax, dwvl)
    nwvl = len(clear.wvl_list)

    clear.layers.gasbelow.level = 1
    clear.layers.gasbelow.press = 1.

    clear.layers.cloud.level = 2
    clear.layers.cloud.press = pc+0.1
    clear.layers.cloud.tau = 0*_np.ones(nwvl)
    clear.layers.cloud.aerosols.nr = nr*_np.ones(nwvl)
    clear.layers.cloud.aerosols.r_eff = reff
    clear.layers.cloud.aerosols.v_eff = 0.10

    clear.layers.gastop.level = 3
    clear.layers.gastop.press = pc
    clear.surface[0,0] = alb

    clear.tag = path_input + 'clear_{:03.1f}_{:3.2f}_{:2.1f}'.format(g,alb,reff)

    check = _np.ones_like(clear.wvl_list)
    for i,wvl in enumerate(clear.wvl_list):
	    filename = clear.tag + '_{:4.3f}.dat'.format(wvl)

	    if _os.path.isfile(filename) and _os.access(filename, _os.R_OK):
		check[i] = 1
	    else:
		check[i] = 0

    if check.prod() ==1:
	    print('Files exist!')
	    clear.name = ['']*_np.size(clear.wvl_list)
	    for i, wvl in enumerate(clear.wvl_list):
	    	clear.name[i] = _os.path.normpath(clear.tag + '_{:4.3f}.dat'.format(wvl))

    del i, wvl, check, nwvl

    return cloudy, clear

def Moon(alb=0.1, g=1, pc=0.7, reff=1.0, nr=1.33, wvlmin=0.4, wvlmax=0.41, dwvl=0.2, path_input = './dap_database/'):

   # Lambertian model

    lamb = _pmd.Model()
    lamb.mma = 29
    lamb.gravity = g
    lamb.dpol = 0.03
    del(lamb.layers.haze)
    lamb.wvl_list = _np.arange(wvlmin, wvlmax, dwvl)
    nwvl = len(lamb.wvl_list)

    lamb.layers.gasbelow.level = 1
    lamb.layers.gasbelow.press = 1.
    lamb.layers.gasbelow.rayscat = False
    lamb.layers.gasbelow.bmsca = _np.zeros(nwvl)

    lamb.layers.cloud.level = 2
    lamb.layers.cloud.press = pc+0.1
    lamb.layers.cloud.tau = 0*_np.ones(nwvl)
    lamb.layers.cloud.aerosols.nr = nr*_np.ones(nwvl)
    lamb.layers.cloud.aerosols.r_eff = reff
    lamb.layers.cloud.aerosols.v_eff = 0.10
    lamb.layers.cloud.rayscat = False
    lamb.layers.cloud.bmsca = _np.zeros(nwvl)

    lamb.layers.gastop.press = pc
    lamb.layers.gastop.level = 3
    lamb.layers.gastop.rayscat = False
    lamb.layers.gastop.bmsca = _np.zeros(nwvl)

    lamb.surface[0,0] = alb

    
    lamb.tag = 'lamb_{:05.3f}_{:4.3f}'.format(g,alb)

    check = _np.ones_like(lamb.wvl_list)
    for i,wvl in enumerate(lamb.wvl_list):
	    filename = lamb.tag + '_{:4.3f}.dat'.format(wvl)

	    if _os.path.isfile(filename) and _os.access(filename, _os.R_OK):
		check[i] = 1
	    else:
		check[i] = 0

    if check.prod() ==1:
	    print('Files exist!')
	    lamb.name = ['']*_np.size(lamb.wvl_list)
	    for i, wvl in enumerate(lamb.wvl_list):
	    	lamb.name[i] = _os.path.normpath(lamb.tag + '_{:4.3f}.dat'.format(wvl))

    del i, wvl, check, nwvl

    return lamb,
