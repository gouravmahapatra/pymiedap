import sys, warnings, numpy as np, time
warnings.filterwarnings('ignore')
sys.path.insert(0, '/sessions/eloquent-admiring-ptolemy/mnt/pymiedap')
import pymiedap.pymiedap as pmd

def _transparent():
    ma=pmd.Aerosols(); ma.typ='G'
    ma.coefs=np.zeros((1,4,4,1)); ma.ncoefs=np.ones(1)
    ma.ssalb=np.zeros(1); ma.sext=np.zeros(1); ma.ssca=np.zeros(1)
    ma.col_dens=0.0; return ma

NMUG=50
t_start=time.time()
aero=pmd.Aerosols(nr=[1.33],ni=[0.],r_eff=10.0,v_eff=0.1,psd='2',typ='C')
pmd.mie_code(aero,[0.550],ngaur=100,nsubr=50)
print(f"Mie done: ncoefs={int(aero.ncoefs[0])}", flush=True)
m=pmd.Model(); m.wvl_list=[0.550]
del m.layers.gasbelow, m.layers.haze
m.layers.gastop.rayscat=False; m.layers.gastop.tau=[0.0]
m.layers.gastop.tau_g=[0.0]; m.layers.gastop.tau_ray=[0.05]
m.layers.cloud.rayscat=False; m.layers.cloud.tau=[4.926]
m.layers.cloud.tau_g=[0.0]; m.layers.cloud.tau_ray=[0.0]
m.dpol=0.03; m.surface[0,0]=0.0
m.layers.cloud.aerosols=aero
m.layers.gastop.mixed_aerosols=_transparent()
m.layers.cloud.mix_aerosols()
m.name=[""]
print(f"Running DAP nmug={NMUG}...", flush=True)
pmd.dap_code(m,rename=True,output_name=f'r10_nmug{NMUG}',
             path_output='dap_database/',nmug=NMUG,nmat=4)
dt=time.time()-t_start
I,Q,_,_=pmd.read_dap_output(
    np.array([80.]),np.array([40.]),np.array([30.]),
    m.name[0],phi=np.array([10.]),beta=np.array([0.]))
print(f"RESULT nmug={NMUG}: DAP_time={dt:.0f}s  I={I[0]:.6f}  NaN={np.isnan(I[0])}")
print(f"STATUS: {'STABLE - valid result!' if not np.isnan(I[0]) else 'UNSTABLE - still NaN'}")
