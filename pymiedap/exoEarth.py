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
import matplotlib
matplotlib.use('Agg')
import pymiedap.pymiedap as pmd
#==============================================================================
# Customized by Ashwyn Groot to incorporate Earth as a model planet
#==============================================================================
# ==============
# IMPORT MODULES
# ==============
import numpy as np
import os, time, pickle, sys, os.path, rasterio, glob, h5py
from datetime import date, timedelta
from dateutil import relativedelta
import matplotlib.pyplot as mpl
from matplotlib import colors
from pyhdf.SD import SD, SDC
from time import sleep
import deepdish as dd

# ---------
# CLASSES DEFINITION
# ---------
if __name__ == "__main__":
      print(__doc__)


def model_generator(force=False, wvl=0.550, mma=29, dpol=0.03, tau=[0.0], tau_g=[0.0], rayscat=True, set_taus=False,
                    v_eff=0.1, psd='2', nr=[1.33], ni=[1e-08], layered=False, nmug_mie=100, nmug=150, nsubr=50, nmat=4,
                    dap_output_path=['./dap_database_0.550/'], data_output_path='./PyMieDAP_Data/', modelfile="modelcodes.txt",
                    bin_surfalb=np.array([0.0278,0.10936,0.98969,0.139387]),bin_surf=np.array([0,14,15,16]),
                    bin_ctpvals=np.array([500.,700.,850.]),
                    bin_cotvals=np.array([20.,10.,5.,0.]),
                    bin_cervals=np.array([17.5,15.,12.5,10.]),
                    custom_modelname=False):

    start=time.time()
    with open(data_output_path+modelfile, "rb") as fp:   # Unpickling
        model_codetot = pickle.load(fp)
    model_codeuniq=list(set(model_codetot))
    model_codesuniq_origi=model_codeuniq
    model_codeuniq_new=[]
    models_seq=[]
    if custom_modelname!=False:
        for cust_model in model_codeuniq:
            model_values=np.array([item for item in cust_model.split('_')])
            counter=0.
            for customname in custom_modelname[::2]:
                custom_values=np.array([item for item in customname.split('_')])
                idx_cusval=np.where(custom_values!='x')
                if len(model_values)==len(custom_values):
                    if all(model_values[idx_cusval]==custom_values[idx_cusval]):
                        model_codeuniq_new.extend([custom_modelname[int(counter+1.)]])
                        break
                    elif int(counter)==int(len(custom_modelname)-2):
                        model_codeuniq_new.extend([cust_model])
                counter+=2.
        model_codeuniq=model_codeuniq_new

    for model_code in list(model_codeuniq):
        while not 'model_'+model_code+'_{:4.3f}'.format(wvl) +'.dat' or 'model_'+model_code+'_{:4.3f}'.format(wvl) +'.hdf5' in os.listdir('./'):
            sleep(600)
        file_indir=[str('model_'+model_code) + '_{:4.7f}.dat'.format(wvl) in os.listdir(paths) for paths in dap_output_path]
        file_indirhdf5=[str('model_'+model_code) + '_{:4.7f}.hdf5'.format(wvl) in os.listdir(paths) for paths in dap_output_path]
        if (file_indir.count(True)>=1 or file_indirhdf5.count(True)>=1) and force==False:
            if file_indir.count(True)>1:
                print('Model seems to be computed twice')
            print(str('model_'+model_code)+' already computed.')
            models_seq.append(str('model_'+model_code))
        else:
            print(str('model_'+model_code)+' being computed.')
            model_values=np.array([int(item) for item in model_code.split('_')])
            if (1.<=bin_surf[model_values[0]] and 14.>=bin_surf[model_values[0]]) or bin_surf[model_values[0]]==16.:
                wvl_ref=[0.645]
            elif bin_surf[model_values[0]]==0.:
                wvl_ref=[0.858]
            elif bin_surf[model_values[0]]==15.:
                wvl_ref=[1.240]

            model_ref=pmd.Model(mma=mma,dpol=dpol,wvl_list=wvl_ref)
            model_ref.layers.cloud.aerosols=pmd.Aerosols(r_eff=bin_cervals[model_values[3]],v_eff=v_eff,psd=psd,nr=nr,ni=ni,layered=layered,typ='C')
            pmd.mie_code(model_ref.layers.cloud.aerosols, wvl_ref,ngaur=nmug_mie)
            sigma_ext_ref=model_ref.layers.cloud.aerosols.sext

            model=pmd.Model(mma=mma,dpol=dpol,wvl_list=[wvl])
            model.layers.cloud.aerosols=pmd.Aerosols(r_eff=bin_cervals[model_values[3]],v_eff=v_eff,psd=psd,nr=nr,ni=ni,layered=layered,typ='C')
            pmd.mie_code(model.layers.cloud.aerosols, [wvl],ngaur=nmug_mie)
            sigma_ext=model.layers.cloud.aerosols.sext
            
            Tau=bin_cotvals[model_values[2]]*(sigma_ext[0]/sigma_ext_ref[0])
            model=pmd.Model(mma=mma,dpol=dpol,wvl_list=[wvl])
            del(model.layers.cloud)
            del(model.layers.gasbelow)
            del(model.layers.gastop)
            del(model.layers.haze)
            model.layers.gasbelow=pmd.Layer(tau=tau,tau_g=tau_g,press=1.0,rayscat=rayscat)
            model.layers.gastop=pmd.Layer(tau=tau,tau_g=tau_g,press=bin_ctpvals[model_values[1]]/1000.,rayscat=rayscat)
            model.layers.cloud=pmd.Layer(tau=[Tau],tau_g=tau_g,press=model.layers.gastop.press+0.1,rayscat=rayscat)
            if model.layers.gasbelow.press==model.layers.cloud.press:
                del(model.layers.gasbelow)
            model.layers.cloud.aerosols=pmd.Aerosols(r_eff=bin_cervals[model_values[3]],v_eff=v_eff,psd=psd,nr=nr,ni=ni,layered=layered,typ='C')
            model.asurf=bin_surfalb[model_values[0]]

            pmd.compute_model(model, force=True, filetype=1, path_input=dap_output_path[0], rename=True,output_name=str('model_'+model_code),set_taus=set_taus, nmug_mie=nmug_mie, nmug=nmug, nsubr=nsubr, nmat=nmat)
            models_seq.append(str('model_'+model_code))
    print(time.time()-start)
    return models_seq,model_codesuniq_origi

def Observation_dates(obs_period,obs_interval,diurnal_var=False,diurnal_res='adaptive'): 
    f_obsday=obs_period[0].split('/')
    l_obsday=obs_period[1].split('/')
    f_date=date(int(f_obsday[0]),int(f_obsday[1]),int(f_obsday[2]))
    l_date=date(int(l_obsday[0]),int(l_obsday[1]),int(l_obsday[2]))+timedelta(days=1)
    f_yearday=date(f_date.year,1,1)
    l_yearday=date(f_date.year,12,31)
    f_delta=f_date-f_yearday
    l_delta=l_date-f_yearday
    full_year=l_yearday-f_yearday+timedelta(days=1)
    # assuming a circular orbit
    period_sunearth=full_year.days      #sidereal year in days 
    if type(obs_interval)!=str:
        obs_dates=np.array([])
        lon_var_delta=np.array([])
        for day in np.arange(f_delta.days+1,l_delta.days+1,obs_interval):
            if diurnal_res=='adaptive':
                if diurnal_var==False or ((f_yearday+timedelta(day-1)).strftime('%Y/%m/%d') not in diurnal_var and diurnal_var!='all'):
                    day_var=day
                    lon_var=-(day-int(day))*360.    
                elif (f_yearday+timedelta(day-1)).strftime('%Y/%m/%d') in diurnal_var or diurnal_var=='all':
                    alpha=(day-1)*360./full_year.days
                    day_resolution=np.ceil(12 * (1 + np.sin(np.radians(alpha)/2.)**2))
                    if day_resolution>30.:
                        day_resolution=30.
                    day_var=np.linspace(int(day),int(day+1),day_resolution+1)[:-1]
                    lon_var=np.linspace(0,-360,day_resolution+1)[:-1]
                lon_var_delta=np.append(lon_var_delta,lon_var)
                obs_dates=np.append(obs_dates,day_var)
            elif type(diurnal_res)!=str:
                if diurnal_var==False or ((f_yearday+timedelta(day-1)).strftime('%Y/%m/%d') not in diurnal_var and diurnal_var!='all'):
                    day_var=day
                    lon_var=-(day-int(day))*360.            
                elif (f_yearday+timedelta(day-1)).strftime('%Y/%m/%d') in diurnal_var or diurnal_var=='all':
                    day_resolution=diurnal_res
                    day_var=np.linspace(day,day+1,day_resolution+1)[:-1]
                    lon_var=np.linspace(0,-360,day_resolution+1)[:-1]
                lon_var_delta=np.append(lon_var_delta,lon_var)
                obs_dates=np.append(obs_dates,day_var)
        alpha=(obs_dates-1)*360./full_year.days
        lon_greenw_offset=-(obs_dates-1)*360./period_sunearth+lon_var_delta
    elif obs_interval=='M':
        obs_dates=np.array([])
        lon_var_delta=np.array([])
        rel_date=relativedelta.relativedelta(l_date,f_date)
        l_month=l_date.month+rel_date.years*12+12
        for month in np.arange(f_date.month,l_month+1,1):
            if month>12:
                day=date(f_date.year+int(month-1)//12,month-12*(int(month-1)//12),1)
            else:
                day=date(f_date.year,month,1)
            day=(day-f_yearday).days+1
            if diurnal_res=='adaptive':
                if diurnal_var==False or ((f_yearday+timedelta(day-1)).strftime('%Y/%m/%d') not in diurnal_var and diurnal_var!='all'):
                    day_var=day
                    lon_var=-(day-int(day))*360.            
                elif (f_yearday+timedelta(day-1)).strftime('%Y/%m/%d') in diurnal_var or diurnal_var=='all':
                    alpha=(day-1)*360./full_year.days
                    day_resolution=np.ceil(12 * (1 + np.sin(np.radians(alpha)/2.)**2))
                    if day_resolution>30.:
                        day_resolution=30.
                    day_var=np.linspace(day,day+1,day_resolution+1)[:-1]
                    lon_var=np.linspace(0,-360,day_resolution+1)[:-1]
                lon_var_delta=np.append(lon_var_delta,lon_var)
                obs_dates=np.append(obs_dates,day_var)
            elif type(diurnal_res)!=str:
                if diurnal_var==False or ((f_yearday+timedelta(day-1)).strftime('%Y/%m/%d') not in diurnal_var and diurnal_var!='all'):
                    day_var=day
                    lon_var=-(day-int(day))*360.              
                elif (f_yearday+timedelta(day-1)).strftime('%Y/%m/%d') in diurnal_var or diurnal_var=='all':
                    day_resolution=diurnal_res
                    day_var=np.linspace(day,day+1,day_resolution+1)[:-1]
                    lon_var=np.linspace(0,-360,day_resolution+1)[:-1]
                lon_var_delta=np.append(lon_var_delta,lon_var)
                obs_dates=np.append(obs_dates,day_var)
        alpha=(obs_dates-1)*360./full_year.days
        lon_greenw_offset=-(obs_dates-1)*360./period_sunearth+lon_var_delta
    return alpha, lon_greenw_offset, obs_dates, f_date,full_year

def asterspectrum_read(spectrum_file,x_int,xyvalues=False):
    f=open(spectrum_file,'r')
    lines=f.readlines()
    index=0
    x=np.array([])
    y=np.array([])
    for line in lines:
        if 'Name' in line:
            name=line.strip('Name:').strip('\r\n').upper()
        if 'Type' in line:
            Type=line.strip('Type:').strip(' ').strip('\r\n').upper()
        if 'X Units' in line:
            xunits=line.strip('X Units:').strip(' ').strip('\r\n')
        if 'Y Units' in line:
            yunits=line.strip('Y Units:').strip(' ').strip('\r\n')
        if 'First Column' in line:
            fcol=line.strip('First Column:').strip(' ').strip('\r\n')
        if 'Second Column' in line:
            scol=line.strip('Second Column:').strip(' ').strip('\r\n')
        if 'First X Value' in line:
            x0=float(line.strip('First X Value:').strip(' ').strip('\r\n'))
        if 'Last X Value' in line:
            xend=float(line.strip('Last X Value:').strip(' ').strip('\r\n'))
        if 'Number of X Values' in line:
            xlen=float(line.strip('Number of X Values:').strip(' ').strip('\r\n'))
    for line in lines:
        if str(x0) in line and 'Value' not in line:
            x0index=index
            break
        else:
            index+=1
    if fcol.lower().strip(' ')=='x' and scol.lower().strip(' ')=='y':
        row_ind1=0
        row_ind2=1
    elif fcol.lower().strip(' ')=='y' and scol.lower().strip(' ')=='x':
        row_ind1=1
        row_ind2=0
    else:
        print('Data structured in a nonlogical way')
    for data in np.arange(x0index,x0index+int(xlen),1):
        row=lines[data].strip('\r\n').split('\t')
        x=np.append(x,float(row[row_ind1]))
        y=np.append(y,float(row[row_ind2]))
    if x[-1]!=xend:
        print('Data not properly extracted')
    else:
        idx_x=np.argsort(x)
        x=x[idx_x]
        y=y[idx_x]
    if x_int<x[0] or x_int>x[-1]:
        print('Data out of range for interpolation')
        return
    else:
        y_inter=pmd.geos.splint(x,y,pmd.geos.spline(x,y,len(x)),len(x),x_int)
    if xyvalues==True:
        return y_inter/100,x,y
    else:
        return y_inter/100

def findModel(searchList, elem):
    idxx=np.array([])
    idxy=np.array([])
    for e in elem:
        for i,x in enumerate(searchList):
            if e in x:
                idxy=np.append(idxy,np.where(x==e)[0])
                idxx=np.append(idxx,np.ones((len(np.where(x==e)[0])))*i)
    return (idxx.astype(int),idxy.astype(int))

def planet_Earth(models, alpha=[0], npix=20, force=False, set_taus=False, rename=True,filetype=1,
                  dap_output_path=['./dap_database_0.550/'],stored_pixvals=True,
                  obs_path='./observation_database/',data_output_path='./PyMieDAP_Data/',
                  windobs=False,plot_obs=True,nmug_mie=100, nmug=150, nsubr=50,
                  nmat=3, diurnal_var='all',diurnal_res=18.,aerosols=False,
                  obliquity=[0.],rotation=[0.],long_alph_custom=False, custom_year=False,custom_modelname=False,
                  custom_glint=False,mask_data=False,
                  bin_surfalb=np.array([0.0278,0.10936,0.98969,0.139387]),bin_surf=np.array([0,14,15,16]), 
                  bin_aot=np.array([0.2625,0.1875,0.001,0.]),bin_aotvals=np.array([0.3,0.225,0.15,0.]),
                  bin_ctp=np.array([600.,800.,1e6]),bin_ctpvals=np.array([500.,700.,850.]),
                  bin_cot=np.array([15.,7.5,0.001,0.]), bin_cotvals=np.array([20.,10.,5.,0.]),
                  bin_cer=np.array([16.25,13.75,11.25,0.]), bin_cervals=np.array([17.5,15.,12.5,10.]),
                  obs_period=['2011/01/01','2011/12/31'], obs_interval=1.):
    start=time.strftime("%Y.%m.%d.%H.%M")
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
    nmug : int, optional
        number of Gauss point for DAP code
    nmug_mie : int, optional
        number of Gauss point for Mie code
    nmat : int, optional
        number of Stokes elements to compute
    nsubr : int, optional
        number of divisions for size dist in Mie calculations
    
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
    ngeosMax=1000000
    atm_model=models[0]
    atm_model.dpol=0.03
    atm_model.mma=29
    del(atm_model.layers.gasbelow)
    del(atm_model.layers.gastop)
    del(atm_model.layers.haze)
    atm_model.layers.gas=pmd.Layer(tau=[0.0],tau_g=[0.0],press=1.0,rayscat=True)
    atm_model.layers.cloud.aerosols=pmd.Aerosols(v_eff=0.1,psd='2',nr=[1.33],ni=[1e-08],layered=False)
    wvl = atm_model.wvl_list
    nwvl = len(atm_model.wvl_list)
    mma=atm_model.mma
    dpol=atm_model.dpol
    tau=atm_model.layers.gas.tau
    tau_g=atm_model.layers.gas.tau_g
    rayscat=atm_model.layers.gas.rayscat
    v_eff=atm_model.layers.cloud.aerosols.v_eff
    psd=atm_model.layers.cloud.aerosols.psd
    nr=atm_model.layers.cloud.aerosols.nr
    ni=atm_model.layers.cloud.aerosols.ni
    layered=atm_model.layers.cloud.aerosols.layered
    mpl.ioff()
    # calculating alpha from observational period and defining filenames
    #-------------------------------------------------------------------
    alpha_obs,lon_greenw_offset, obs_dates,f_obsday_year,full_year=Observation_dates(obs_period,obs_interval,diurnal_var=diurnal_var,diurnal_res=diurnal_res)
    f_obsday_year=f_obsday_year.year
    if long_alph_custom==False:
        alpha=alpha_obs
    elif long_alph_custom!=False and long_alph_custom!=True:
            lon_greenw_offset=long_alph_custom
    elif long_alph_custom==True:
        if len(alpha)!=len(lon_greenw_offset):
            alpha=alpha[0]*np.ones(len(lon_greenw_offset))
        else:
            alpha=alpha
    nalpha = len(alpha)
    atm_model.fcloud = np.zeros(nalpha)
    atm_model.asym = []
    atm_model.mask=[]
    if type(obs_interval)!=bool:
        if float(obs_interval)==8.:
            obs_interval_product='E'
        else:
            obs_interval_product='D'
    if type(custom_year)==np.ndarray:
        obs_dates+=custom_year
        obs_dates[obs_dates>=full_year.days+1]-=full_year.days
    product='MYD08_'+str(obs_interval_product)+'3'
    filenames=[product+'.A'+str(int(f_obsday_year))+str(int(item)).zfill(3) for item in np.floor(obs_dates)]
    if long_alph_custom!=False:
        filenames=filenames*len(alpha)
        obs_dates=obs_dates*np.ones(len(alpha))
        
    # Preparing arrays
    # ------------------
    If = np.zeros((nwvl,nalpha,(npix)**2))
    Qf = np.zeros((nwvl,nalpha,(npix)**2))
    Uf = np.zeros((nwvl,nalpha,(npix)**2))
    Vf = np.zeros((nwvl,nalpha,(npix)**2))

    phaf = np.zeros((nwvl,nalpha,(npix)**2))
    szaf = np.zeros((nwvl,nalpha,(npix)**2))
    emif = np.zeros((nwvl,nalpha,(npix)**2))
    azif = np.zeros((nwvl,nalpha,(npix)**2))
    betf = np.zeros((nwvl,nalpha,(npix)**2))
    xf = np.zeros((nwvl,nalpha,(npix)**2))
    yf = np.zeros((nwvl,nalpha,(npix)**2))

    if len(alpha)!=len(obliquity):
        obliquity = obliquity[0]*np.ones(len(alpha))
    if len(alpha)!=len(rotation):
        rotation = rotation[0]*np.ones(len(alpha))

    # Loop on wvl
    # -------------
    for j,w in enumerate(wvl):

    # Loop on phase angle
    # -------------------
    # Create table to store output
        It = np.zeros((nalpha,(npix)**2))
        Qt = np.zeros((nalpha,(npix)**2))
        Ut = np.zeros((nalpha,(npix)**2))
        Vt = np.zeros((nalpha,(npix)**2))
        
        Modelnames=[[]]*len(alpha)
        Modelnames_orig=[[]]*len(alpha)
        masks=[[]]*len(alpha)
        Patchy_clouds_mask=[[]]*len(alpha)
        
        pixelval_filename='Pixelvals_dres'+str(int(diurnal_res))+'_obsint'+str(int(obs_interval))+'_'+str(obs_period[0]).replace("/","-")+'_'+str(int(obs_dates[0]))+'_'+str(obs_period[1]).replace("/","-")+'_'+str(int(obs_dates[-1]))+'_wind'+str(windobs)+'_binsurf'+str(len(bin_surf))+'_npix'+str(npix)+'.pix'
        if stored_pixvals==True and os.path.exists(data_output_path+pixelval_filename) and mask_data==False:
            with open(data_output_path+pixelval_filename, "rb") as fp:   # Unpickling
                [masks,Patchy_clouds_mask,Unique_models,Unique_models_orig] = pickle.load(fp)
                
            for A,alph in enumerate(alpha):
                #Get geom
                ngeos, apix, theta0, theta, phi, beta, lats, longs, xs, ys = pmd.geos.getgeos(alph, npix)
        
                theta0 = theta0[:ngeos]
                theta = theta[:ngeos]
                phi = phi[:ngeos]
                beta = beta[:ngeos]
                lats = lats[:ngeos]
                longs = longs[:ngeos]        
        
                phaf[j,A,:ngeos] = alph*np.ones(ngeos)
                szaf[j,A,:ngeos] = theta0[:ngeos]
                emif[j,A,:ngeos] = theta[:ngeos]
                azif[j,A,:ngeos] = phi[:ngeos]
                betf[j,A,:ngeos] = beta[:ngeos]
                xf[j,A,:ngeos] = xs[:ngeos]
                yf[j,A,:ngeos] = ys[:ngeos]
                
                print('Computing only geometries for filename:',filenames[A])
                
        else:
            for A,alph in enumerate(alpha):        
                #Get geom
                ngeos, apix, theta0, theta, phi, beta, lats, longs, xs, ys = pmd.geos.getgeos(alph, npix)
        
                theta0 = theta0[:ngeos]
                theta = theta[:ngeos]
                phi = phi[:ngeos]
                beta = beta[:ngeos]
                lats = lats[:ngeos]
                longs = longs[:ngeos]        
        
                phaf[j,A,:ngeos] = alph*np.ones(ngeos)
                szaf[j,A,:ngeos] = theta0[:ngeos]
                emif[j,A,:ngeos] = theta[:ngeos]
                azif[j,A,:ngeos] = phi[:ngeos]
                betf[j,A,:ngeos] = beta[:ngeos]
                xf[j,A,:ngeos] = xs[:ngeos]
                yf[j,A,:ngeos] = ys[:ngeos]
        
                print(filenames[A])

                # call mask_Earth
                mask, ncloud, asym, Data_means, Patchy_clouds = mask_Earth(alpha=alph, npix=npix,
                                                                        windobs=windobs,plot_obs=plot_obs,
                                                                        aerosols=aerosols,custom_glint=custom_glint,
                                                                        filename=filenames[A],obs_day=obs_dates[A],
                                                                        longitudinalpos=lon_greenw_offset[A],
                                                                        obliquity=obliquity[A],rotation=rotation[A],nobins=False,
                                                                        bin_surf=bin_surf,bin_ctp=bin_ctp,
                                                                        bin_ctpvals=bin_ctpvals,bin_cotvals=bin_cotvals,bin_cervals=bin_cervals,
                                                                        bin_cot=bin_cot,bin_cer=bin_cer,start=start,obs_path=obs_path, data_output_path=data_output_path)
                masks[A]=mask
                Patchy_clouds_mask[A]=Patchy_clouds
                
                atm_model.fcloud[A]=ncloud
                atm_model.asym.append(asym)
                atm_model.asym.append(Data_means)
                atm_model.mask.append([mask,Patchy_clouds])
                models,model_codesuniq_origi=model_generator(force=force, wvl=w, mma=mma, dpol=dpol, tau=tau, tau_g=tau_g, rayscat=rayscat,
                                                                   v_eff=v_eff, psd=psd, nr=nr, ni=ni, layered=layered,bin_surfalb=bin_surfalb,bin_surf=bin_surf,
                                                                   nmug_mie=nmug_mie, nmug=nmug, nsubr=nsubr, nmat=nmat,set_taus=set_taus,bin_ctpvals=bin_ctpvals,
                                                                   bin_cotvals=bin_cotvals,bin_cervals=bin_cervals,dap_output_path=dap_output_path,
                                                                  data_output_path=data_output_path,custom_modelname=custom_modelname)
                
                Modelnames[A]=models
                Modelnames_orig[A]=model_codesuniq_origi
            Unique_models=list(set(sum(Modelnames, [])))
            Unique_models_orig=[]

            for model in Unique_models:
                idx_model=sum(Modelnames,[]).index(model)
                Unique_models_orig.append(sum(Modelnames_orig,[])[idx_model])

            with open(data_output_path+pixelval_filename, 'wb') as file:
                pickle.dump([masks,Patchy_clouds_mask,Unique_models,Unique_models_orig], file, protocol=2)
        
            if mask_data==True:
                return
        
        phase = phaf[j,:,:]
        theta0 = szaf[j,:,:]
        theta = emif[j,:,:]
        phi = azif[j,:,:]
        beta = betf[j,:,:]
        xs = xf[j,:,:]
        ys = yf[j,:,:]


        if filetype==1 or filetype==2:
              file_ext='.hdf5'
        elif filetype==3:
              file_ext='.dat'

        print('The following model files are read:')
        print(Unique_models)

        for pixtype,model in enumerate(Unique_models):
              if custom_modelname==False:
                    pixtype=model.split('model_')[-1]
                    pixtype=pixtype.encode()
                    print('Reading '+str(model) + '_{:4.7f}'.format(w)+file_ext)
              elif custom_modelname!=False:
                    pixtype=Unique_models_orig[pixtype]
                    pixtype=pixtype.encode()
                    print('Reading '+str(model) + '_{:4.7f}'.format(w)+file_ext+' instead of model_'+pixtype.decode()+'_{:4.7f}'.format(w)+file_ext)

              idx_pixels=findModel(masks,[pixtype])
              n=len(idx_pixels[0])
              geo_limit=np.append(np.arange(0,n,ngeosMax),np.array([n]))
              print('Computing  ',n,'  pixels with the current model type.')
              
              phaseB = phase[idx_pixels]
              theta0B = theta0[idx_pixels]
              thetaB = theta[idx_pixels]
              phiB = phi[idx_pixels]
              betaB = beta[idx_pixels]
              
              PCL=np.array([]).reshape(0,len(Patchy_clouds_mask[0][0]))
              for m,mask in enumerate(masks):
                  idx_pcl=np.where(mask==pixtype)
                  Patchy_clouds=Patchy_clouds_mask[m]
                  PCL=np.append(PCL,Patchy_clouds[idx_pcl],axis=0)

              if len(phaseB)==0:
                  print('No corresponding pixels with this model atmosphere.')
              else:
                  #==============================================================================
                  #           In case of aerosols, clouds and clear atmosphere in the same pixel
                  #==============================================================================
                  if Patchy_clouds.shape[1]>2 and aerosols==True:
                      if any(PCL[:,0]!=0.) or any(PCL[:,2]!=0.):
                          phaseB_PCL=phaseB
                          theta0B_PCL = theta0B
                          thetaB_PCL = thetaB
                          phiB_PCL = phiB
                          betaB_PCL = betaB
                          pixtype=pixtype.decode()
                          if pixtype.split('_')[2]=='3':
                              print('PCL Error')
                          model_cloud=pixtype.split('_')[:-1]
                          model_clear=pixtype.split('_')[:-1]
                          model_clear[2]='3'
                          model_cloud='model_'+'_'.join(model_cloud)
                          model_clear='model_'+'_'.join(model_clear)
                          if custom_modelname!=False:
                              model_values=np.array([item for item in model_clear.split('model_')[-1].split('_')])
                              counter=0.
                              for customname in custom_modelname[::2]:
                                  custom_values=np.array([item for item in customname.split('_')])
                                  idx_cusval=np.where(custom_values!='x')
                                  if len(model_values)==len(custom_values):
                                      if all(model_values[idx_cusval]==custom_values[idx_cusval]):
                                          model_clear='model_'+custom_modelname[int(counter+1.)]
                                          break
                                      elif int(counter)==int(len(custom_modelname)-2):
                                          model_clear=model_clear
                                  counter+=2.
                              modelcld_values=np.array([item for item in model_cloud.split('model_')[-1].split('_')])
                              counter=0.
                              for customname in custom_modelname[::2]:
                                  custom_values=np.array([item for item in customname.split('_')])
                                  idx_cusval=np.where(custom_values!='x')
                                  if len(modelcld_values)==len(custom_values):
                                      if all(modelcld_values[idx_cusval]==custom_values[idx_cusval]):
                                          model_cloud='model_'+custom_modelname[int(counter+1.)]
                                          break
                                      elif int(counter)==int(len(custom_modelname)-2):
                                          model_cloud=model_cloud
                                  counter+=2.
                          print('Reading a Patchy cloud with aerosol pixels, consisting of (cloudy,clear,aerosol): '+str(model_cloud)+file_ext+', '+str(model_clear)+file_ext+', '+str(model)+file_ext)
                          for output_path in dap_output_path:
                              if os.path.exists(output_path+str(model_clear) + '_{:4.7f}'.format(w)+file_ext):
                                  filepath_clear=output_path+str(model_clear) + '_{:4.7f}'.format(w)+file_ext
                                  break
                          for output_path in dap_output_path:
                              if os.path.exists(output_path+str(model_cloud) + '_{:4.7f}'.format(w)+file_ext):
                                  filepath_cloud=output_path+str(model_cloud) + '_{:4.7f}'.format(w)+file_ext
                                  break
                          for output_path in dap_output_path:
                              if os.path.exists(output_path+str(model) + '_{:4.7f}'.format(w)+file_ext):
                                  filepath=output_path+str(model) + '_{:4.7f}'.format(w)+file_ext
                                  break
                          IB_clear=np.zeros((n));QB_clear=np.zeros((n));UB_clear=np.zeros((n));VB_clear=np.zeros((n))
                          IB_cloud=np.zeros((n));QB_cloud=np.zeros((n));UB_cloud=np.zeros((n));VB_cloud=np.zeros((n))
                          IB_aer=np.zeros((n));QB_aer=np.zeros((n));UB_aer=np.zeros((n));VB_aer=np.zeros((n))
                          for g,geo_iter in enumerate(np.arange(0,n,ngeosMax)):
                              idxlow=geo_limit[g]
                              idxhigh=geo_limit[g+1]
                              IB_clear[idxlow:idxhigh],QB_clear[idxlow:idxhigh],UB_clear[idxlow:idxhigh],VB_clear[idxlow:idxhigh] = pmd.read_dap_output(phaseB_PCL[idxlow:idxhigh],theta0B_PCL[idxlow:idxhigh],thetaB_PCL[idxlow:idxhigh],
                                      filepath_clear,filetype=filetype,phi=phiB_PCL[idxlow:idxhigh], beta=betaB_PCL[idxlow:idxhigh])
                              IB_cloud[idxlow:idxhigh],QB_cloud[idxlow:idxhigh],UB_cloud[idxlow:idxhigh],VB_cloud[idxlow:idxhigh] = pmd.read_dap_output(phaseB_PCL[idxlow:idxhigh],theta0B_PCL[idxlow:idxhigh],thetaB_PCL[idxlow:idxhigh],
                                      filepath_cloud,filetype=filetype,phi=phiB_PCL[idxlow:idxhigh], beta=betaB_PCL[idxlow:idxhigh])
                              IB_aer[idxlow:idxhigh],QB_aer[idxlow:idxhigh],UB_aer[idxlow:idxhigh],VB_aer[idxlow:idxhigh] = pmd.read_dap_output(phaseB[idxlow:idxhigh],theta0B[idxlow:idxhigh],thetaB[idxlow:idxhigh],
                                      filepath,filetype=filetype,phi=phiB[idxlow:idxhigh], beta=betaB[idxlow:idxhigh])
                          if any((PCL[:,2]/PCL[:,1]+(1-(PCL[:,0]+PCL[:,2])/PCL[:,1])+(PCL[:,0]/PCL[:,1]))>1):
                              print('Fractions not computed correctly')
                          IB=IB_aer*(PCL[:,2]/PCL[:,1])+IB_clear*(1-(PCL[:,0]+PCL[:,2])/PCL[:,1])+IB_cloud*(PCL[:,0]/PCL[:,1])
                          QB=QB_aer*(PCL[:,2]/PCL[:,1])+QB_clear*(1-(PCL[:,0]+PCL[:,2])/PCL[:,1])+QB_cloud*(PCL[:,0]/PCL[:,1])
                          UB=UB_aer*(PCL[:,2]/PCL[:,1])+UB_clear*(1-(PCL[:,0]+PCL[:,2])/PCL[:,1])+UB_cloud*(PCL[:,0]/PCL[:,1])
                          VB=VB_aer*(PCL[:,2]/PCL[:,1])+VB_clear*(1-(PCL[:,0]+PCL[:,2])/PCL[:,1])+VB_cloud*(PCL[:,0]/PCL[:,1])
                          pixtype=pixtype.encode()
                      else:
                          for output_path in dap_output_path:
                              if os.path.exists(output_path+str(model) + '_{:4.7f}'.format(w))+file_ext:
                                  filepath=output_path+str(model) + '_{:4.7f}'.format(w)+file_ext
                                  break
                              
                          IB=np.zeros((n));QB=np.zeros((n));UB=np.zeros((n));VB=np.zeros((n))
                          for g,geo_iter in enumerate(np.arange(0,n,ngeosMax)):
                              idxlow=geo_limit[g]
                              idxhigh=geo_limit[g+1]
                              IB[idxlow:idxhigh],QB[idxlow:idxhigh],UB[idxlow:idxhigh],VB[idxlow:idxhigh] = pmd.read_dap_output(phaseB[idxlow:idxhigh],theta0B[idxlow:idxhigh],thetaB[idxlow:idxhigh],
                                      filepath,filetype=filetype,phi=phiB[idxlow:idxhigh], beta=betaB[idxlow:idxhigh])
                  
                  #==============================================================================
                  #        In case of clouds and clear pixels at the same pixel
                  #==============================================================================
                  else:
                      if any(PCL[:,0]!=0.):
                          phaseB_PCL=phaseB
                          theta0B_PCL = theta0B
                          thetaB_PCL = thetaB
                          phiB_PCL = phiB
                          betaB_PCL = betaB
                          pixtype=pixtype.decode()
                          if pixtype.split('_')[2]=='3':
                              print(PCL)
                              sys.exit("PCL Error")
                          model_clear=pixtype.split('_')
                          model_clear[2]='3'
                          model_clear='model_'+'_'.join(model_clear)
                          if custom_modelname!=False:
                              model_values=np.array([item for item in model_clear.split('model_')[-1].split('_')])
                              counter=0.
                              for customname in custom_modelname[::2]:
                                  custom_values=np.array([item for item in customname.split('_')])
                                  idx_cusval=np.where(custom_values!='x')
                                  if len(model_values)==len(custom_values):
                                      if all(model_values[idx_cusval]==custom_values[idx_cusval]):
                                          model_clear='model_'+custom_modelname[int(counter+1.)]
                                          break
                                      elif int(counter)==int(len(custom_modelname)-2):
                                          model_clear=model_clear
                                  counter+=2.
                          print('Reading a Patchy cloud, consisting of (cloudy,clear): '+str(model)+file_ext+', '+str(model_clear)+file_ext)
                          for output_path in dap_output_path:
                              if os.path.exists(output_path+str(model_clear) + '_{:4.7f}'.format(w)+file_ext):
                                  filepath_clear=output_path+str(model_clear) + '_{:4.7f}'.format(w)+file_ext
                                  break
                          for output_path in dap_output_path:
                              if os.path.exists(output_path+str(model) + '_{:4.7f}'.format(w)+file_ext):
                                  filepath=output_path+str(model) + '_{:4.7f}'.format(w)+file_ext
                                  break
                          IB_clear=np.zeros((n));QB_clear=np.zeros((n));UB_clear=np.zeros((n));VB_clear=np.zeros((n))
                          IB_cloud=np.zeros((n));QB_cloud=np.zeros((n));UB_cloud=np.zeros((n));VB_cloud=np.zeros((n))
                          for g,geo_iter in enumerate(np.arange(0,n,ngeosMax)):
                              idxlow=geo_limit[g]
                              idxhigh=geo_limit[g+1]
                              IB_clear[idxlow:idxhigh],QB_clear[idxlow:idxhigh],UB_clear[idxlow:idxhigh],VB_clear[idxlow:idxhigh] = pmd.read_dap_output(phaseB_PCL[idxlow:idxhigh],theta0B_PCL[idxlow:idxhigh],thetaB_PCL[idxlow:idxhigh],
                                      filepath_clear,filetype=filetype,phi=phiB_PCL[idxlow:idxhigh], beta=betaB_PCL[idxlow:idxhigh])
                              IB_cloud[idxlow:idxhigh],QB_cloud[idxlow:idxhigh],UB_cloud[idxlow:idxhigh],VB_cloud[idxlow:idxhigh] = pmd.read_dap_output(phaseB[idxlow:idxhigh],theta0B[idxlow:idxhigh],thetaB[idxlow:idxhigh],
                                      filepath,filetype=filetype,phi=phiB[idxlow:idxhigh], beta=betaB[idxlow:idxhigh])
                          if any(((1-PCL[:,0]/PCL[:,1])+(PCL[:,0]/PCL[:,1]))>1):
                              print('Fractions not computed correctly')
                          IB=IB_cloud*(PCL[:,0]/PCL[:,1])+IB_clear*(1-PCL[:,0]/PCL[:,1])
                          QB=QB_cloud*(PCL[:,0]/PCL[:,1])+QB_clear*(1-PCL[:,0]/PCL[:,1])
                          UB=UB_cloud*(PCL[:,0]/PCL[:,1])+UB_clear*(1-PCL[:,0]/PCL[:,1])
                          VB=VB_cloud*(PCL[:,0]/PCL[:,1])+VB_clear*(1-PCL[:,0]/PCL[:,1])
                          pixtype=pixtype.encode()
                      else:
                          for output_path in dap_output_path:
                              if os.path.exists(output_path+str(model) + '_{:4.7f}'.format(w)+file_ext):
                                  filepath=output_path+str(model) + '_{:4.7f}'.format(w)+file_ext
                                  break
                          IB=np.zeros((n));QB=np.zeros((n));UB=np.zeros((n));VB=np.zeros((n))
                          for g,geo_iter in enumerate(np.arange(0,n,ngeosMax)):
                              idxlow=geo_limit[g]
                              idxhigh=geo_limit[g+1]
                              IB[idxlow:idxhigh],QB[idxlow:idxhigh],UB[idxlow:idxhigh],VB[idxlow:idxhigh] = pmd.read_dap_output(phaseB[idxlow:idxhigh],theta0B[idxlow:idxhigh],thetaB[idxlow:idxhigh],
                                      filepath,filetype=filetype,phi=phiB[idxlow:idxhigh], beta=betaB[idxlow:idxhigh])
    
                  It[idx_pixels] = IB*np.cos(np.radians(theta0B))
                  Qt[idx_pixels] = QB*np.cos(np.radians(theta0B))
                  Ut[idx_pixels] = UB*np.cos(np.radians(theta0B))
                  Vt[idx_pixels] = VB*np.cos(np.radians(theta0B))
    
        If[j,:,:] = It
        Qf[j,:,:] = Qt
        Uf[j,:,:] = Ut
        Vf[j,:,:] = Vt

    # Compute adjusted radiance
    B = pmd.sunblackbody(np.array(atm_model.wvl_list)*1e-6, Ts=atm_model.Ts)
    Rs = 696342000.
    Dvs = 108208930000.0
    omegas = np.pi * (atm_model.Rs/atm_model.Dps)**2
    I2 = If * B[:,np.newaxis,np.newaxis] * omegas


    output_dir='Run_'+start+'_'+filenames[0].split('.')[0]
    try:
        os.makedirs(os.path.normpath(data_output_path+output_dir))
    except OSError:
        if not os.path.isdir(os.path.normpath(data_output_path+output_dir)):
            raise
    h5_filename = data_output_path+output_dir+'/'+"Model_"+str(int(f_obsday_year))+"_"+str(obs_interval)+"_"+str(npix)+"_"+time.strftime("%Y.%m.%d.%H.%M")

    dd.io.save(h5_filename, dict(model=atm_model))

    with h5py.File(h5_filename, 'a') as f:
        f.create_dataset("pha", data=phaf)
        f.create_dataset("x", data=xf)
        f.create_dataset("y", data=yf)
        f.create_dataset("npix", data=npix)
        f.create_dataset("I", data=If)
        f.create_dataset("Q", data=Qf)
        f.create_dataset("U", data=Uf)
        f.create_dataset("V", data=Vf)
    f.close()
    
    mpl.ion()


def plot_Earth(model, wvl_idx=0, display='grid', stokes='Ps', phase_idx=0,title='Polarization', cmap='YlOrRd',vmin=None,vmax=None,
                font_size=12, figsize=[8,8], Ylim=[-10,35,5], dpi=100, diskintegrate=False, data_scale=1.,data_only=False,
                data_direct=False,xin=None,yin=None,npix=None):

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
        which Stokes element to plot. Allowed are 'Ps' (-Q/I), 'I', 'Q', 'U',
        'V', 'Pt' (total polarization), 'Pl' total linear pol.
        Default is 'Ps'
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

    if data_only==False:
        fig = mpl.figure(figsize=(figsize), dpi=dpi)
        
    if data_direct==True:
        Z=np.nan*np.zeros(xin.shape[1])
        for phase_idx,alph in enumerate(xin[wvl_idx,:,0]):
            ngeos, apix, theta0, theta, phi, beta, lats, longs, xs, ys = pmd.geos.getgeos(alph, npix)
            Yx=yin[wvl_idx,phase_idx,:]*apix
            Z[phase_idx] = np.nansum(Yx)
        x=xin[wvl_idx,:,0]
        x=x[~np.isnan(Z)]
        Z=Z[~np.isnan(Z)]
        Z = data_scale * Z
        return x,Z
    elif diskintegrate==True:
        npix=model.npix
        Z=np.nan*np.zeros(model.phase.shape[1])
        for phase_idx,alph in enumerate(model.phase[wvl_idx,:,0]):
                ngeos, apix, theta0, theta, phi, beta, lats, longs, xs, ys = pmd.geos.getgeos(alph, npix)
                if stokes=='F':
                    Ix=model.I[wvl_idx,phase_idx,:]*apix
                    Z[phase_idx] = np.nansum(Ix)
                if stokes=='Q':
                    Qx=model.Q[wvl_idx,phase_idx,:]*apix
                    Z[phase_idx] = np.nansum(Qx)
                if stokes=='U':
                    Ux=model.U[wvl_idx,phase_idx,:]*apix
                    Z[phase_idx] = np.nansum(Ux)
                if stokes=='V':
                    Vx=model.V[wvl_idx,phase_idx,:]*apix
                    Z[phase_idx] = np.nansum(Vx)
                if stokes=='Ps':
                    Psx=-np.nansum(model.Q[wvl_idx,phase_idx,:])/(np.nansum(model.I[wvl_idx,phase_idx,:])+1e-40)
                    Z[phase_idx] = Psx
                if stokes=='Pl':
                    Plx=np.sqrt(np.nansum(model.Q[wvl_idx,phase_idx,:])**2 + np.nansum(model.U[wvl_idx,phase_idx,:])**2)/(np.nansum(model.I[wvl_idx,phase_idx,:])+1e-40)
                    Z[phase_idx] = Plx
                if stokes=='Pt':
                    Ptx=np.sqrt(np.nansum(model.Q[wvl_idx,phase_idx,:])**2 + np.nansum(model.U[wvl_idx,phase_idx,:])**2 + np.nansum(model.V[wvl_idx,phase_idx,:])**2)/(np.nansum(model.I[wvl_idx,phase_idx,:])+1e-40)
                    Z[phase_idx] = Ptx
                if stokes=='Pv':
                    Pvx=np.nansum(model.V[wvl_idx,phase_idx,:])/(np.nansum(model.I[wvl_idx,phase_idx,:])+1e-40)
                    Z[phase_idx] = Pvx
                if stokes=='Pu':
                    Pux=np.nansum(model.U[wvl_idx,phase_idx,:])/(np.nansum(model.I[wvl_idx,phase_idx,:])+1e-40)
                    Z[phase_idx] = Pux
        x=model.phase[wvl_idx,:,0]
        x=x[~np.isnan(Z)]
        Z=Z[~np.isnan(Z)]
        Z = data_scale * Z
        if data_only==False:
            ax = fig.add_subplot(111)
            ax.set_title(title)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
            ax.plot(x,Z,'k-')
            fig.tight_layout(pad=1.2)
            mpl.xlabel('Phase angle')
            mpl.ylabel(title)
            mpl.xticks(fontsize=14)
            mpl.yticks(np.round(np.arange(Ylim[0], Ylim[1], Ylim[2]),12), [str(xaxe) for xaxe in np.round(np.arange(Ylim[0], Ylim[1], Ylim[2]),12)], fontsize=14)
            for yaxe in np.round(np.arange(Ylim[0], Ylim[1], Ylim[2]),12):
                mpl.plot(x, [yaxe] * len(x), "--", lw=0.5, color="black", alpha=0.3)
            mpl.legend(stokes)
            mpl.ylim(Ylim[0],Ylim[1])
            mpl.autoscale(enable=True, axis='x', tight=True)
        return x,Z
    else:
        npix=model.npix
        stokes=stokes[0]
        alph=model.phase[wvl_idx,:,0]
        ngeos, apix, theta0, theta, phi, beta, lats, longs, xs, ys = pmd.geos.getgeos(alph, npix)
        ax = fig.add_subplot(111,aspect=0.5)
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
    
        if stokes=='F':
            Z = model.I[wvl_idx,phase_idx,:]*apix
        if stokes=='Q':
            Z = model.Q[wvl_idx,phase_idx,:]*apix
        if stokes=='U':
            Z = model.U[wvl_idx,phase_idx,:]*apix
        if stokes=='V':
            Z = model.V[wvl_idx,phase_idx,:]*apix
        if stokes=='Ps':
            Z = model.P[wvl_idx,phase_idx,:]
        if stokes=='Pl':
            Z = model.Pl[wvl_idx,phase_idx,:]
        if stokes=='Pt':
            Z = model.Pt[wvl_idx,phase_idx,:]
        if stokes=='IPl':
            Z1 = model.I[wvl_idx,phase_idx,:]*apix
            Z2 = model.Pl[wvl_idx,phase_idx,:]
            Z=Z1*Z2
    
    
        Z = data_scale * Z

        ax.set_title(title)
        sc = ax.scatter(X, Y, c=Z,lw=0, marker='s',s=area*scalingx*scalingy,cmap=cmap, zorder=10,vmin=vmin, vmax=vmax)
        fig.tight_layout(pad=1.2)
        cb = fig.colorbar(sc,pad=0.02, extend='both')
        cb.set_label(stokes,size=font_size)
        ax.set_aspect('equal')

    return fig,ax

def mask_Earth(alpha=0, npix=100, filename='MYD08_D3.A2011115', obs_day=115.,
                obs_path='./observation_database/',data_output_path='./PyMieDAP_Data/',
                longitudinalpos=0., obliquity=0., rotation=0.,custom_glint=False,
                nobins=False, iceclouds=False,windobs=False,plot_obs=True,aerosols=False,
                bin_surfalb=np.array([0.0278,0.10936,0.98969,0.139387]),
                bin_surf=np.array([1,1,1,1,1,2,4,2,2,2,4,2,4,2,3,4,0]),
                bin_aot=np.array([0.2625,0.1875,0.0001,0.]),bin_aotvals=np.array([0.3,0.225,0.15,0.]),
                bin_ctp=np.array([600.,800.,1e6]),bin_ctpvals=np.array([500.,700.,850.]),
                bin_cot=np.array([15.,7.5,0.001,0.]), bin_cotvals=np.array([20.,10.,5.,0.]),
                bin_cer=np.array([16.25,13.75,11.25,0.]), bin_cervals=np.array([17.5,15.,12.5,10.]),
                start=time.strftime("%Y.%m.%d")):
    starttime=time.time()             
    """ Generates a mask that can be used for inhomogeneous planets
    
    Parameters
    ----------
    alpha : int,optional
        phase angle at which the calculation is made
        default is 0
    npix : int, optional
        number of pixels on each axis
        default is 20
    
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
    
    # read the pixel geometries
    ngeos, apix, theta0, theta, phi, beta, lats, longs, xs, ys = pmd.geos.getgeos(alpha, npix)
    
    # prepare grids
    grid_lit = np.zeros((npix,npix))
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
    xs = xs.round(12)
    ys = ys.round(12)
    
    # X and Y axis
    step = 2./npix
    X = -1 + 0.5*step + np.arange(0,npix)*step
    Y = -1 + 0.5*step + np.arange(0,npix)*step
    X = X.round(12)
    Y = Y.round(12)
    xv, yv = np.meshgrid(X,Y)
 
    mpl.ioff()
    # open land cover file
    model_classtot=np.array([[],[],[],[],[]])
    surffile = rasterio.open(obs_path+'MCD12Q1.006_LC_Type1_doy2011_SAMP.tif')
    ds=surffile.read(1)
    surffile.close()
    for file in glob.glob(obs_path+filename+'*Cloud_Top_Pressure_Mean_INT_SAMP'+'*.tif'):
        CTPfile=rasterio.open(file)
        dsCTP=CTPfile.read(1)
        CTPfile.close()
    for file in glob.glob(obs_path+filename+'*Cloud_Optical_Thickness_Liquid_Log_Mean_INT_SAMP'+'*.tif'):
        COTfile=rasterio.open(file)
        dsCOT=COTfile.read(1)
        COTfile.close()
    for file in glob.glob(obs_path+filename+'*Cloud_Effective_Radius_Liquid_Mean_INT_SAMP'+'*.tif'):
        CERfile=rasterio.open(file)
        dsCER=CERfile.read(1)
        CERfile.close()
    for file in glob.glob(obs_path+filename+'*Cloud_Fraction_Mean_INT_SAMP'+'*.tif'):
        CFfile=rasterio.open(file)
        dsCF=CFfile.read(1)
        CFfile.close()
    if windobs==True:
        bin_wind=np.array([-100,6])
        bin_windvals=np.array([5,7])
        obs_hour=(obs_day-int(obs_day))*24
        obs_bin=np.arange(0,24,6)[::-1]
        obs_hour_idx=np.digitize(obs_hour,obs_bin,right=False)
        print('Reading '+filename+'.'+str(int(obs_bin[obs_hour_idx])).zfill(2)+' as a multiday windspeed observation')
        for file in glob.glob(obs_path+filename+'.'+str(int(obs_bin[obs_hour_idx])).zfill(2)+'*Ocean_Surface_Wind_Speed_10m_INT_SAMP.tif'):
            dsWdsp=rasterio.open(file)
            dsWdsp=dsWdsp.read(1)
        model_classtot=np.array([[],[],[],[],[],[]])
    if iceclouds==True:
        for file in glob.glob(obs_path+filename+'*Cloud_Effective_Radius_Ice_Mean_INT_SAMP'+'*.tif'):
            dsCERice=rasterio.open(file)
            dsCERice=dsCERice.read(1)
        for file in glob.glob(obs_path+filename+'*Cloud_Optical_Thickness_Ice_Log_Mean_INT_SAMP'+'*.tif'):
            dsCOTice=rasterio.open(file)
            dsCOTice=dsCOTice.read(1)
        model_classtot=np.array([[],[],[],[],[],[],[]])
    if aerosols==True:
        print('Reading an aerosol observation.')
        for file in glob.glob(obs_path+filename+'*AOD_550_Dark_Target_Deep_Blue_Combined_Mean_INT_SAMP'+'*.tif'):
            dsAOT=rasterio.open(file)
            dsAOT=dsAOT.read(1)
        if aerosols==True and windobs==False and iceclouds==False:
            model_classtot=np.array([[],[],[],[],[],[]])
        elif aerosols==True and windobs==True and iceclouds==False:
            model_classtot=np.array([[],[],[],[],[],[],[]])
        elif aerosols==True and windobs==True and iceclouds==True:
            model_classtot=np.array([[],[],[],[],[],[],[],[]])
    for file in glob.glob(obs_path+'MYD08_D3'+'*.hdf'):
        file = SD(file, SDC.READ)
        ctp_id = file.select('Cloud_Top_Pressure_Mean')
        cotL_id = file.select('Cloud_Optical_Thickness_Liquid_Log_Mean')
        cerL_id = file.select('Cloud_Effective_Radius_Liquid_Mean')
        cf_id = file.select('Cloud_Retrieval_Fraction_Liquid')
        if aerosols==True:
            aot_id = file.select('AOD_550_Dark_Target_Deep_Blue_Combined_Mean')
        if iceclouds==True:
            cotI_id = file.select('Cloud_Optical_Thickness_Ice_Log_Mean')
            cerI_id = file.select('Cloud_Effective_Radius_Ice_Mean')

    if custom_glint==True:
        idx_glint=glint_coordinates(clear_ocean_path='/home2/victor/eDAP2/Database_1layer_noclouds/model_0_0_3_1_7_0.8650000.dat',alpha=alpha,npix=npix,wvl=0.865,plot_glint=False)
    model_codetot=[]
    model_codestore=[]
    Data_means=[]
    Patchy_clouds=[]
    ssxidx=[]
    ssyidx=[]
    Y=Y[::-1]

    nlat=ds.shape[0]
    nlon=ds.shape[1]
    phi=np.ones((nlat,nlon))
    lamb=np.ones((nlat,nlon))
    lambsub=np.linspace(-180,180,nlon+1)[1:]
    phisub=np.linspace(90,-90,nlat)
    for i in np.arange(0,nlat,1):
          phi[i,:]=phi[i,:]*phisub[i]
    for i in np.arange(0,nlon,1):
          lamb[:,i]=lamb[:,i]*lambsub[i]
    mask=pmd.geos.obsmask(phi,lamb,obliquity,longitudinalpos,rotation,npix,alpha,nlon,nlat)
    
    for i in np.arange(0,ngeos,1):
        #==============================================================================
        #   The next set of lines extracts the observations
        #==============================================================================
        idx_obs=np.where(mask==i)
        datasetsurf=ds[idx_obs]
        if len(datasetsurf)==0:
            print('Pixel not in LC observations')
        datasetcf=(dsCF[idx_obs]-cf_id.attributes()['add_offset'])*cf_id.attributes()['scale_factor']
        if len(datasetcf)==0:
            print('Pixel not in CF observations')
        datasetctp=(dsCTP[idx_obs]-ctp_id.attributes()['add_offset'])*ctp_id.attributes()['scale_factor']
        if len(datasetctp)==0:
            print('Pixel not in CTP observations')
        datasetcot=(dsCOT[idx_obs]-cotL_id.attributes()['add_offset'])*cotL_id.attributes()['scale_factor']
        if len(datasetcot)==0:
            print('Pixel not in COT observations')
        datasetcer=(dsCER[idx_obs]-cerL_id.attributes()['add_offset'])*cerL_id.attributes()['scale_factor']
        if len(datasetcer)==0:
            print('Pixel not in CER observations')

        if windobs==True:
            datasetwdsp=dsWdsp[idx_obs]
            if len(datasetwdsp)==0:
                print('Pixel not in Wind observations')

        if aerosols==True:
            datasetaot=(dsAOT[idx_obs]-aot_id.attributes()['add_offset'])*aot_id.attributes()['scale_factor']
            if len(datasetaot)==0:
                print('Pixel not in Aerosol observations')

        if iceclouds==True:
            datasetcotice=(dsCOTice[idx_obs]-aot_id.attributes()['add_offset'])*aot_id.attributes()['scale_factor']
            if len(datasetcotice)==0:
                print('Pixel not in COT ice observations')
            
            datasetcerice=(dsCERice[idx_obs]-aot_id.attributes()['add_offset'])*aot_id.attributes()['scale_factor']
            if len(datasetcerice)==0:
                print('Pixel not in CER ice observations')
        
        datasetcf=np.round(datasetcf,5)
        datasetcf_backup=np.copy(datasetcf)
        if np.sum(datasetcf)==0.:
            ctp=700.
            cer=12.5
            cot=5.
        else:
            if len(datasetctp)==0:
                ctp=700.
            else:
                if len(datasetctp[datasetctp>0.])==0.:
                    ctp=700.
                else:
                    ctp=np.sum((datasetctp[datasetctp>0.]*datasetcf[datasetctp>0.]))/(np.sum(datasetcf[datasetctp>0.]))
            if len(datasetcer)==0:
                cer=12.5
            else:
                if len(datasetcer[datasetcer>0.])==0.:
                    cer=12.5
                else:
                    cer=np.sum((datasetcer[datasetcer>0.]*datasetcf[datasetcer>0.]))/(np.sum(datasetcf[datasetcer>0.]))
            if len(datasetcot)==0:
                cot=5.
            else:
                if len(datasetcot[datasetcot>-10.])==0.:
                    cot=5.
                else:
                    cot_log_mean=np.sum((datasetcot[datasetcot>-10.]*datasetcf[datasetcot>-10.]))/(np.sum(datasetcf[datasetcot>-10.]))
                    cot=10**cot_log_mean
        try:
            vals,counts=np.unique(datasetsurf,return_counts=True)
            IGBP=vals[np.argmax(counts)]
        except (IndexError,ValueError):
            IGBP=17
        if int(IGBP)==255 or int(IGBP)==-1: # Unclassified regions are set to ocean.
            IGBP=17
        ############
        if iceclouds==True:
            if len(datasetcerice)==0:
                cerice=30.
            else:
                if len(datasetcerice[datasetcerice>0.])==0.:
                    cerice=30.
                else:
                    cerice=np.sum((datasetcerice[datasetcerice>0.]*datasetcf[datasetcerice>0.]))/(np.sum(datasetcf[datasetcerice>0.]))
            if len(datasetcotice)==0:
                cotice=2.
            else:
                if len(datasetcot[datasetcot>-10.])==0.:
                    cotice=2.
                else:
                    cotice_log_mean=np.sum((datasetcotice[datasetcotice>-10.]*datasetcf[datasetcotice>-10.]))/(np.sum(datasetcf[datasetcotice>-10.]))
                    cotice=10**cotice_log_mean
        if windobs==True:
            if len(datasetwdsp)==0. or int(IGBP)!=17. or len(datasetwdsp[datasetwdsp>0.])==0:
                wdsp=0.
            else:
                wdsp=np.mean(datasetwdsp[datasetwdsp>0.])
        if aerosols==True:
            if len(datasetAOT)==0. or len(datasetAOT[datasetAOT>0.])==0.:
                Aot=0.
                AeroFrac=0.
            else:
                Aot=np.mean(datasetAOT[datasetAOT>0.])
                AeroFrac=len(datasetAOT[datasetAOT!=0.])/len(datasetAOT)
        ############
        #there always has to be a zero cot or cer in the bins for the no data!!!!!!!!!

        idxlclass = bin_surf[IGBP-1] # remove 1 because python array index starts at 0!

        if np.isnan(ctp) or np.isnan(cot) or np.isnan(cer):
            idxctp=np.argmin(bin_ctpvals)
            idxcot=np.where(bin_cotvals==0)[0][0]
            idxcer=np.argmin(bin_cervals)
        else:
            idxctp=np.digitize(ctp,bin_ctp,right=True)
            idxcot=np.digitize(cot,bin_cot,right=False)
            idxcer=np.digitize(cer,bin_cer,right=False)
            if aerosols==True:
                idxaot=np.digitize(Aot,bin_aot,right=False)
            if windobs==True:
                idxwnd=np.digitize(wdsp,bin_wind,right=True)

        if idxctp==len(bin_ctp):
            idxctp=idxctp-1
        if idxcot==len(bin_cot):
            idxcot=idxcot-1
        if idxcer==len(bin_cer):
            idxcer=idxcer-1
        if len(datasetcf)==0 and i!=0:
            datasetcf=np.array([Data_means[-1][3]])
        # plot index on grid
        newxs = xs[i]
        newys = ys[i]
        ssx = np.where(X==newxs)[0][0]
        ssxidx=np.append(ssxidx,ssx).astype(int)
        ssy = np.where(Y==newys)[0][0]
        ssyidx=np.append(ssyidx,ssy).astype(int)
        datasetcf_backup=datasetcf_backup[datasetcf_backup>=0.]
        if len(datasetcf_backup)==0 and i!=0:
            datasetcf_backup=np.array([Data_means[-1][3]])
        elif len(datasetcf_backup)==0 and i==0:
            datasetcf_backup=np.array([0.])
        if custom_glint==True:
            if IGBP<=0:
                if i in idx_glint[0]:
                    datasetcf=np.zeros(len(datasetcf))
                    datasetcf_backup=np.zeros(len(datasetcf_backup))
        if nobins==False:
            if iceclouds==True and windobs==False and aerosols==False:
                model_classtot=np.append(model_classtot,[[IGBP],[bin_ctpvals[idxctp]],[bin_cotvals[idxcot]],[bin_cervals[idxcer]],[np.mean(datasetcf_backup)],[cotice],[cerice]],axis=1)
            
            if windobs==True and aerosols==False and iceclouds==False:
                model_code=str(idxlclass)+str(idxwnd)+'_'+str(idxctp)+'_'+str(idxcot)+'_'+str(idxcer)
                model_codestore.append(str(idxlclass)+str(idxwnd)+'_'+str(idxctp)+'_'+'3'+'_'+str(idxcer))
                model_classtot=np.append(model_classtot,[[IGBP],[bin_ctpvals[idxctp]],[bin_cotvals[idxcot]],[bin_cervals[idxcer]],[np.mean(datasetcf_backup)],[bin_windvals[idxwnd-1]]],axis=1)
            
            if aerosols==True and windobs==True and iceclouds==False:
                model_code=str(idxlclass)+str(idxwnd)+'_'+str(idxctp)+'_'+str(idxcot)+'_'+str(idxcer)+'_'+str(idxaot)
                model_classtot=np.append(model_classtot,[[IGBP],[bin_ctpvals[idxctp]],[bin_cotvals[idxcot]],[bin_cervals[idxcer]],[np.mean(datasetcf_backup)],[bin_windvals[idxwnd-1]],[bin_aotvals[idxaot]]],axis=1)
                model_codestore.append(str(idxlclass)+str(idxwnd)+'_'+str(idxctp)+'_'+str(idxcot)+'_'+str(idxcer))
                model_codestore.append(str(idxlclass)+str(idxwnd)+'_'+str(idxctp)+'_'+'3'+'_'+str(idxcer))
            
            if iceclouds==False and windobs==False and aerosols==False:
                model_code=str(idxlclass)+'_'+str(idxctp)+'_'+str(idxcot)+'_'+str(idxcer)
                model_codestore.append(str(idxlclass)+'_'+str(idxctp)+'_'+'3'+'_'+str(idxcer))
                model_classtot=np.append(model_classtot,[[IGBP],[bin_ctpvals[idxctp]],[bin_cotvals[idxcot]],[bin_cervals[idxcer]],[np.mean(datasetcf_backup)]],axis=1)
        
        elif nobins==True:
            if iceclouds==True and windobs==False and aerosols==False:
                model_code=str(idxlclass)+'_'+str(idxctp)+'_'+str(idxcot)+'_'+str(idxcer)
                model_classtot=np.append(model_classtot,[[IGBP],[ctp],[cot],[cer],[np.mean(datasetcf_backup)],[cotice],[cerice]],axis=1)
            
            if windobs==True and aerosols==False and iceclouds==False:
                model_code=str(idxlclass)+str(idxwnd)+'_'+str(idxctp)+'_'+str(idxcot)+'_'+str(idxcer)
                if sum(datasetcf_backup)!=len(datasetcf_backup) and sum(datasetcf_backup)!=0.:
                    model_codestore.append(str(idxlclass)+str(idxwnd)+'_'+str(idxctp)+'_'+'3'+'_'+str(idxcer))
                model_classtot=np.append(model_classtot,[[IGBP],[ctp],[cot],[cer],[np.mean(datasetcf_backup)],[wdsp]],axis=1)
            
            if aerosols==True and windobs==True and iceclouds==False:
                model_code=str(idxlclass)+str(idxwnd)+'_'+str(idxctp)+'_'+str(idxcot)+'_'+str(idxcer)+'_'+str(idxaot)
                model_classtot=np.append(model_classtot,[[IGBP],[ctp],[cot],[cer],[np.mean(datasetcf_backup)],[Aot],[wdsp]],axis=1)
            
            if iceclouds==False and windobs==False and aerosols==False:
                model_code=str(idxlclass)+'_'+str(idxctp)+'_'+str(idxcot)+'_'+str(idxcer)
                if sum(datasetcf_backup)!=len(datasetcf_backup) and sum(datasetcf_backup)!=0.:
                    model_codestore.append(str(idxlclass)+'_'+str(idxctp)+'_'+'3'+'_'+str(idxcer))
                model_classtot=np.append(model_classtot,[[IGBP],[ctp],[cot],[cer],[np.mean(datasetcf_backup)]],axis=1)
        if custom_glint==True:
            if IGBP<=0:
                if i in idx_glint[0]:
                    model_code=model_code.split('_')
                    model_code[2]='3'
                    model_code='_'.join(model_code)
        model_codetot.append(model_code)
        model_codestore.append(model_code)
        Data_means.append([ctp,cot,cer,np.mean(datasetcf_backup)])
        if aerosols==True:
            TotAF=(1-np.sum(datasetcf_backup)/len(datasetcf_backup))*AeroFrac*len(datasetcf_backup)
            Patchy_clouds.append([np.sum(datasetcf_backup),len(datasetcf_backup),TotAF])
        else:
            Patchy_clouds.append([np.sum(datasetcf_backup),len(datasetcf_backup)])
        sys.stdout.write('\r'+str("{:2.2f}".format(i*100./(ngeos-1)).zfill(5)+" %"))
        sys.stdout.flush()
    Data_means=np.mean(np.array(Data_means),axis=0)
    dataland=np.copy(grid_lit)
    dataland[ssyidx,ssxidx]=model_classtot[0]
    datapres=np.copy(grid_lit)
    datapres[ssyidx,ssxidx]=model_classtot[1]
    dataopac=np.copy(grid_lit)
    dataopac[ssyidx,ssxidx]=model_classtot[2]
    dataeffr=np.copy(grid_lit)
    dataeffr[ssyidx,ssxidx]=model_classtot[3]
    datacfrac=np.copy(grid_lit)
    datacfrac[ssyidx,ssxidx]=model_classtot[4]
    if windobs==True and aerosols==False:
        datawind=np.copy(grid_lit)
        datawind[ssyidx,ssxidx]=model_classtot[5]
    if aerosols==True and windobs==False:
        dataAot=np.copy(grid_lit)
        dataAot[ssyidx,ssxidx]=model_classtot[5]
    if aerosols==True and windobs==True:
        datawind=np.copy(grid_lit)
        datawind[ssyidx,ssxidx]=model_classtot[5]
        dataAot=np.copy(grid_lit)
        dataAot[ssyidx,ssxidx]=model_classtot[6]
    if iceclouds==True:
        dataopacice=np.copy(grid_lit)
        dataopacice[ssyidx,ssxidx]=model_classtot[5]
        dataeffrice=np.copy(grid_lit)
        dataeffrice[ssyidx,ssxidx]=model_classtot[6]
    Land_frac=np.array([])
    Surfvals,Surfcounts=np.unique(dataland,return_counts=True)
    for c,L_type in enumerate(Surfvals):
        frac=Surfcounts[c]
        if ngeos!=0:
            frac=frac/float(ngeos)
        else:
            frac=np.nan
        Land_frac=np.append(Land_frac,frac)
    Data_means=[Data_means,Land_frac]
    AllData_lit=[dataland,datapres,dataopac,dataeffr,datacfrac]
    output_dir='Run_'+start+'_'+filename.split('.')[0]
    try:
        os.makedirs(os.path.normpath(data_output_path+output_dir))
    except OSError:
        if not os.path.isdir(os.path.normpath(data_output_path+output_dir)):
            raise
    if plot_obs==True:
        cmap = colors.ListedColormap(['forestgreen','forestgreen','forestgreen','forestgreen','forestgreen','yellowgreen','y','khaki','gold','limegreen','lightseagreen','green','grey','yellowgreen','snow','moccasin','darkblue'])
        bounds = np.linspace(1,17,33).tolist()
        norm = colors.BoundaryNorm(bounds,cmap.N)
        if iceclouds==True:
            f, axarr = mpl.subplots(7, 1,figsize=(6,32))
        elif windobs==True and aerosols==False:
            f, axarr = mpl.subplots(6, 1,figsize=(6,32))
        elif aerosols==True and windobs==False:
            f, axarr = mpl.subplots(6, 1,figsize=(6,32))
        elif aerosols==True and windobs==True:
            f, axarr = mpl.subplots(7, 1,figsize=(6,32))
        else:
            f, axarr = mpl.subplots(5, 1,figsize=(6,32))
#            mpl.tight_layout()
        im=axarr[0].matshow(dataland, norm=norm, cmap=cmap)
        axarr[0].set_title('Land cover',fontsize=14, family='serif')
        f.colorbar(im,ax=axarr[0])
        im=axarr[1].matshow(datapres, cmap='YlOrRd')
        axarr[1].set_title('Cloud Top Pressure',fontsize=14, family='serif')
        f.colorbar(im,ax=axarr[1])
        im=axarr[2].matshow(dataopac, cmap='YlOrRd')
        axarr[2].set_title('Cloud Optical Thickness',fontsize=14, family='serif')
        f.colorbar(im,ax=axarr[2])
        im=axarr[3].matshow(dataeffr, cmap='YlOrRd')
        axarr[3].set_title('Cloud Effective Radius',fontsize=14, family='serif')
        f.colorbar(im,ax=axarr[3])
        im=axarr[4].matshow(datacfrac, cmap='YlOrRd')
        axarr[4].set_title('Cloud Fraction',fontsize=14, family='serif')
        f.colorbar(im,ax=axarr[4])
        if windobs==True and aerosols==False:
            im=axarr[5].matshow(datawind, cmap='YlOrRd')
            axarr[5].set_title('Wind speed (m/s)',fontsize=14, family='serif')
            f.colorbar(im,ax=axarr[5])
        if aerosols==True and windobs==False:
            im=axarr[5].matshow(dataAot, cmap='YlOrRd')
            axarr[5].set_title('Aerosol Optical Depth',fontsize=14, family='serif')
            f.colorbar(im,ax=axarr[5])
        if aerosols==True and windobs==True:
            im=axarr[5].matshow(datawind, cmap='YlOrRd')
            axarr[5].set_title('Wind speed (m/s)',fontsize=14, family='serif')
            f.colorbar(im,ax=axarr[5])
            im=axarr[6].matshow(dataAot, cmap='YlOrRd')
            axarr[6].set_title('Aerosol Optical Depth',fontsize=14, family='serif')
            f.colorbar(im,ax=axarr[6])
        if iceclouds==True:
            im=axarr[5].matshow(dataopacice, cmap='YlOrRd')
            axarr[5].set_title('Cloud Optical Thickness (Ice)',fontsize=14, family='serif')
            f.colorbar(im,ax=axarr[5])
            im=axarr[6].matshow(dataeffrice, cmap='YlOrRd')
            axarr[6].set_title('Cloud Effective Radius (Ice)',fontsize=14, family='serif')
            f.colorbar(im,ax=axarr[6])
        mpl.savefig(data_output_path+output_dir+'/'+str(npix)+'_'+str(np.round(alpha,2))+'_'+str(abs(np.round(longitudinalpos,2))) + '.png',format = 'png', transparent=True)
        mpl.close()
    try:
        len_bitstr=len(max(model_codetot, key=len))
    except ValueError:
        len_bitstr=0
    if aerosols==True:
        grid_lit_Patchies=np.ones((npix,npix,3))*np.nan
    else:
        grid_lit_Patchies=np.ones((npix,npix,2))*np.nan
    grid_lit=grid_lit.astype(dtype='|S'+str(len_bitstr))
    grid_lit[ssyidx,ssxidx]= np.array(model_codetot)
    if ngeos!=0:
        grid_lit_Patchies[ssyidx,ssxidx]=np.array(Patchy_clouds)
    with open(data_output_path+"modelcodes.txt", "wb") as fp:
        pickle.dump(model_codestore, fp, protocol=2)
    ds=None
    Y=Y[::-1]
    grid_lit=np.flipud(grid_lit)
    grid_lit_Patchies=np.flipud(grid_lit_Patchies)
    print(time.time()-starttime)
    mpl.ion()

        # get current cloud coverage at given phase angle
    cl=np.where(dataopac[np.where(~np.isnan(dataopac))]>0)[0].size
    lit = np.where(~np.isnan(dataopac))[0].size
    if ngeos!=0:
        nb_cloud=float(cl)/(lit)
    else:
        nb_cloud=np.nan
    asym=[]
    xv, yv = np.meshgrid(X,Y)
    asymmetries=4
    for grid in AllData_lit:
      diffgrid=[]
            # 0 degrees
      diffgrid.append(grid[:npix//2,:] - (grid[(npix//2):,:])[::-1,:])
            # 90 degrees
      diffgrid.append(grid[:,:npix//2] - (grid[:,(npix//2):])[:,::-1])
            #+45 degrees
      diffgrid.append(grid[np.where(yv<-xv)]-grid[(np.array([yv.shape[0]-1])-np.where(yv<-xv))[1],(np.array([yv.shape[0]-1])-np.where(yv<-xv))[0]])
            #-45 degrees
      diffgrid.append(grid[np.where(yv<xv)]-grid[np.where(yv<xv)[::-1]])
      diffgrid=np.array(diffgrid)
      asym_sub=[]
      for types in np.arange(0,asymmetries,1):
            asympix = np.where(diffgrid[types][~np.isnan(diffgrid[types])]!=0)[0]
            if np.size(diffgrid[types][~np.isnan(diffgrid[types])])==0:
                  asym_sub_sub = 0
            else:
                  asym_sub_sub = np.size(asympix)/np.float(np.size(diffgrid[types][~np.isnan(diffgrid[types])]))
            asym_sub.append(asym_sub_sub)
      asym.append(asym_sub)

    # Flatten the grid with only lit points
    grid_out=grid_lit[np.where(grid_lit!='nan'.encode())]
    grid_lit_Patchies=grid_lit_Patchies[np.where(grid_lit!='nan'.encode())]
    grid_out = grid_out.flatten()
    
    return grid_out, nb_cloud, asym, Data_means, grid_lit_Patchies
    
def glint_coordinates(clear_ocean_path='/home2/victor/eDAP2/Database_1layer_noclouds/model_0_0_3_1_7_0.8650000.dat',alpha=87,npix=50,wvl=0.865,threshold=2.5,plot_glint=False):
    model=Model()
    model.wvl_list=[wvl]
    model.name[0]=clear_ocean_path
    
    pmd.planet_pixels([model],alpha=[alpha],npix=npix,full_disk=True)
    if plot_glint==True:
        plot_Earth(model,stokes=['Q'],figsize=[10,10])
    
    ngeos, apix, theta0, theta, phi, beta, lats, longs, xs, ys = pmd.geos.getgeos(alpha, npix)
    
    Ps=model.Q[0][0][:ngeos]
    threshold=np.nanmean(Ps)-threshold*np.nanstd(Ps)
    idx_glint=np.where(Ps<threshold)
    
    if plot_glint==True:
        X=xs[:ngeos]
        Y=ys[:ngeos]
        Z=np.zeros(ngeos)
        Z[idx_glint]=1.
        W=[str(a) for a in np.arange(0,ngeos,1)]
        dpi=120
        fig=mpl.figure(figsize=([10,10]),dpi=dpi)
        ax=fig.add_subplot(111,aspect=0.5)
        circ=mpl.Circle((0,0),1,color='gray')
        ax.add_patch(circ)
        ax.set_xlim(-1,1)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_ylim(-1,1)
        bbox=ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        axh=bbox.height
        axw=bbox.width
        scalingx=np.ones_like(X)
        scalingy=np.ones_like(Y)
        area=(axh*dpi*axw*dpi)/(1.5*npix)**2
        ax.set_title('Masked planetary disk (20x20)')
        ax.scatter(X,Y,c=Z,lw=0,marker='s',s=area*scalingx*scalingy,cmap='YlOrRd',zorder=2)
        fig.tight_layout(pad=1.2)
        ax.set_aspect('equal')
        for i,txt in enumerate(W):
            ax.annotate(txt,(X[i],Y[i]),size=7.,ha='center',va='center',rotation=45,color='black')
        mpl.xlabel('X coordinates')
        mpl.ylabel('Y coordinates')
    return idx_glint

# -----------
# MAIN
# -----------

if __name__ == '__main__':

    print('ok')
