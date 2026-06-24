import sys, warnings, numpy as np, time
warnings.filterwarnings('ignore')
sys.path.insert(0, '/sessions/eloquent-admiring-ptolemy/mnt/pymiedap')
import pymiedap.pymiedap as pmd
import os

os.chdir('/sessions/eloquent-admiring-ptolemy/mnt/pymiedap')

LOG = '/sessions/eloquent-admiring-ptolemy/mnt/pymiedap/nmug50_result.txt'

def _transparent(nwav=1):
    ma = pmd.Aerosols(); ma.typ = 'G'
    ma.coefs = np.zeros((nwav,4,4,1)); ma.ncoefs = np.ones(nwav)
    ma.ssalb = np.zeros(nwav); ma.sext = np.zeros(nwav); ma.ssca = np.zeros(nwav)
    ma.col_dens = 0.0; return ma

NMUG = 50
with open(LOG, 'w') as f:
    f.write(f'[{time.strftime("%H:%M:%S")}] nmug={NMUG} starting Mie\n'); f.flush()

aero = pmd.Aerosols(nr=[1.33], ni=[0.], r_eff=10.0, v_eff=0.1, psd='2', typ='C')
pmd.mie_code(aero, [0.550], ngaur=100, nsubr=50)

with open(LOG, 'a') as f:
    f.write(f'[{time.strftime("%H:%M:%S")}] Mie done, ncoefs={aero.ncoefs}\n'); f.flush()

m = pmd.Model(); m.wvl_list = [0.550]
del m.layers.gasbelow, m.layers.haze
m.layers.gastop.rayscat = False; m.layers.gastop.tau = [0.0]
m.layers.gastop.tau_g = [0.0]; m.layers.gastop.tau_ray = [0.05]
m.layers.cloud.rayscat = False; m.layers.cloud.tau = [4.926]
m.layers.cloud.tau_g = [0.0]; m.layers.cloud.tau_ray = [0.0]
m.dpol = 0.03; m.surface[0,0] = 0.0
m.layers.cloud.aerosols = aero
m.layers.gastop.mixed_aerosols = _transparent()
m.layers.cloud.mix_aerosols()
m.name = ['']

with open(LOG, 'a') as f:
    f.write(f'[{time.strftime("%H:%M:%S")}] Starting dap_code nmug={NMUG}\n'); f.flush()

t0 = time.time()
pmd.dap_code(m, rename=True, output_name=f'r10_nmug{NMUG}',
             path_output='dap_database/', nmug=NMUG, nmat=4)
dt = time.time() - t0

with open(LOG, 'a') as f:
    f.write(f'[{time.strftime("%H:%M:%S")}] dap_code done in {dt:.0f}s\n'); f.flush()

I, Q, _, _ = pmd.read_dap_output(
    np.array([80.]), np.array([40.]), np.array([30.]),
    m.name[0], phi=np.array([10.]), beta=np.array([0.]))

with open(LOG, 'a') as f:
    f.write(f'nmug={NMUG}: DAP={dt:.0f}s  I={I[0]:.5f}  NaN={np.isnan(I[0])}\n')
    f.write(f"=> {'STABLE' if not np.isnan(I[0]) else 'STILL NaN'}\n")
    f.flush()

print(f'nmug={NMUG}: DAP={dt:.0f}s  I={I[0]:.5f}  NaN={np.isnan(I[0])}')
print(f"=> {'STABLE' if not np.isnan(I[0]) else 'STILL NaN'}")
