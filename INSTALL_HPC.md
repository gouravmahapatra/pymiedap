# Installing PyMieDAP on a Linux HPC Cluster

This guide covers everything needed to run PyMieDAP — including the ocean surface model, CKD pipeline, and all figure reproduction scripts — on a university Linux HPC cluster with a Python virtual environment.

## Requirements

| Requirement | Version | Notes |
|-------------|---------|-------|
| Python | ≥ 3.9, **3.11 recommended** | Must be available via `module load` or system |
| gfortran | ≥ 7 | Part of GCC; almost always available |
| numpy | **< 2.0** | The Fortran f2py interface requires numpy 1.x |
| scipy | ≥ 1.9 | |
| matplotlib | ≥ 3.5 | |
| h5py | ≥ 3.7 | For HDF5 caches in CKD pipeline |
| hitran-api | ≥ 1.2 | For the `ckdistribution` absorption pipeline |

The `numpy < 2.0` constraint is critical. numpy 2.x changed the C-API and will produce an `ImportError` when loading the compiled Fortran modules.

---

## Quick start (automated)

```bash
# Clone or copy the repo to the cluster
git clone <repo-url> pymiedap    # or: scp -r pymiedap user@cluster:~/
cd pymiedap

# Load Python and GCC on your cluster first (see below)
module load python/3.11          # example — exact name varies by cluster
module load gcc/12               # provides gfortran

# Run the installer
bash install_hpc.sh
```

The script creates a virtual environment at `./venv/`, installs all dependencies, compiles all five Fortran modules with f2py, and runs a smoke test.

To use a custom venv path or a specific Python interpreter:
```bash
bash install_hpc.sh --python python3.11 --venv /scratch/$USER/pymiedap_venv
```

---

## Step-by-step (manual)

### 1. Load HPC modules

Every cluster uses a different module naming convention. Use `module avail` to find the right names:

```bash
module avail python    # find Python modules
module avail gcc       # find GCC/gfortran modules

# Load them (examples — replace with your cluster's actual names):
module load python/3.11
module load gcc/12
```

If your cluster uses Anaconda/Conda instead:
```bash
module load anaconda3
conda create -n pymiedap python=3.11 -y
conda activate pymiedap
```

> **No `sudo` on the cluster?** That's expected — you don't need root for any
> of this. `gfortran` should come from a module (`module avail gcc` →
> `module load gcc/...`). If no module provides a Fortran compiler, install one
> into your own environment without admin rights via conda
> (`conda install -c conda-forge gfortran`) or [Spack](https://spack.io). The
> `sudo apt-get`/`sudo yum` commands shown later apply **only to machines where
> you are root** (e.g. your own workstation), never to a shared cluster.

### 2. Create and activate a virtual environment

```bash
cd ~/pymiedap                      # or wherever you placed the repo
python3 -m venv venv
source venv/bin/activate
```

Verify you have the right Python:
```bash
python --version                   # should be 3.9+
which python                       # should point inside venv/
```

### 3. Install Python dependencies

```bash
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
```

Verify numpy version:
```bash
python -c "import numpy; print(numpy.__version__)"   # must be 1.x
```

If numpy 2.x was installed despite the pin (can happen when another package pulls it in), force-downgrade:
```bash
pip install "numpy>=1.23,<2.0" --force-reinstall
```

### 4. Compile the Fortran extensions

All five extensions must be compiled on the cluster with f2py. The `.so` files already in the repository are for other platforms and should be ignored.

```bash
cd ~/pymiedap

# Mie scattering module
python -m numpy.f2py -c -m module_mie \
    mie_sig.pyf \
    mie_source/anbn.f mie_source/devel.f mie_source/fichid.f \
    mie_source/gauleg.f mie_source/pitau.f mie_source/rminmax.f \
    mie_source/scatmat.f mie_source/sizedis.f mie_source/writsc.f \
    --f77flags="-O2"

# Read Mie output module
python -m numpy.f2py -c -m module_readmie \
    readmie_sig.pyf \
    readmie_source/read_mie_output.f readmie_source/file2coefs.f \
    --f77flags="-O2"

# Coated-sphere Mie module
python -m numpy.f2py -c -m module_mieshell \
    mieshell_sig.pyf \
    mieshell_source/bhcoat.f mieshell_source/gauleg.f mieshell_source/pitau.f \
    mieshell_source/scatmat.f mieshell_source/sizedis.f \
    --f77flags="-O2"

# Adding-doubling radiative transfer module
python -m numpy.f2py -c -m module_dap \
    dap_sig.pyf \
    dap_source/adding.f dap_source/addlay.f dap_source/addsm.f \
    dap_source/assign.f dap_source/bmolecules.f dap_source/brack.f \
    dap_source/bstart.f dap_source/double.f dap_source/expbmu.f \
    dap_source/fillup.f dap_source/gauleg.f dap_source/init.f \
    dap_source/layer0.f dap_source/layerm.f dap_source/ldiapr.f \
    dap_source/newfou.f dap_source/nobot.f dap_source/notop.f \
    dap_source/ord1m.f dap_source/ord2m.f dap_source/prod.f \
    dap_source/rdiapr.f dap_source/renorm.f dap_source/scalzm.f \
    dap_source/setfou.f dap_source/setmu.f dap_source/setzm.f \
    dap_source/star.f dap_source/top2bot.f dap_source/trace.f \
    dap_source/transf.f dap_source/tstar.f \
    --f77flags="-O2"

# Disk geometry / DAP output reader module
python -m numpy.f2py -c -m module_geos \
    geos_sig.pyf \
    geos_source/bracks.f geos_source/getgeos.f geos_source/rdfous.f \
    geos_source/read_dap.f geos_source/spline.f geos_source/splint.f \
    geos_source/wrout.f \
    --f77flags="-O2"
```

Each command takes 10–60 seconds. Successful compilation produces a file like `module_mie.cpython-311-x86_64-linux-gnu.so` in the repo root.

### 5. Verify the installation

```bash
python -c "
import module_mie, module_dap, module_geos, module_readmie, module_mieshell
import pymiedap.pymiedap as pmd
print('All modules OK')
print(pmd.Model())
"
```

---

## Activating the environment in future sessions

Each time you log in, you need to load modules and activate the venv:

```bash
module load python/3.11
module load gcc/12
source ~/pymiedap/venv/bin/activate
```

You can add these lines to your `~/.bashrc` or a project-specific `activate.sh` script.

---

## Running on compute nodes (SLURM example)

For long computations, submit to the queue. Example job script:

```bash
#!/bin/bash
#SBATCH --job-name=pymiedap_fig7
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=fig7_%j.log

module load python/3.11
module load gcc/12
source ~/pymiedap/venv/bin/activate
cd ~/pymiedap

# Run each step sequentially
python examples/trees_stam_fig7.py --step dap
python examples/trees_stam_fig7.py --step images
python examples/trees_stam_fig7.py --step plot
```

Submit with:
```bash
sbatch run_fig7.slurm
```

For PBS/Torque clusters, replace `#SBATCH` directives with `#PBS` equivalents.

---

## Troubleshooting

**`ImportError: numpy.core.multiarray failed to import`**
The Fortran `.so` was compiled with numpy 1.x but numpy 2.x is active.
Fix: `pip install "numpy>=1.23,<2.0" --force-reinstall` then recompile.

**`ModuleNotFoundError: No module named 'module_mie'`**
The `.so` files are platform-specific. You must compile on the cluster, not copy from another machine.
Fix: Run the compilation commands in Step 4 above.

**`meson.build ERROR: Compiler ... cannot compile programs`** (Python ≥ 3.12 only)
numpy ≥ 2.0 with Python ≥ 3.12 uses the Meson backend and requires a working gfortran.
Fix: Use Python 3.11 with `numpy<2.0` (the combination that has been tested).

**`FileNotFoundError: [Errno 2] No such file or directory: 'meson'`**
You are on Python ≥ 3.12 and numpy ≥ 2.0. Install meson: `pip install meson ninja`.
Alternatively, switch to Python 3.11 + numpy < 2.0.

**NaN values in disk images / Fourier files**
Increase `N_MUG` (Gauss quadrature points) in the example script. The default is 40 which is the minimum for numerical stability. Try 50 or 60 if NaN persists at specific wavelengths.

**Slow computation on the login node**
Do not run heavy computations on the login node. Submit a job via SLURM/PBS (see above). The `--step dap`, `--step images`, and `--step plot` sub-commands in the example scripts are designed for this workflow.

---

## Directory structure

```
pymiedap/
├── install_hpc.sh          # automated HPC installer (this file's companion)
├── INSTALL_HPC.md          # this guide
├── requirements.txt        # pinned Python dependencies
├── setup.py                # setuptools config (alternative to manual f2py)
├── pymiedap/               # Python package
│   ├── pymiedap.py         # main Model / compute_model API
│   ├── ocean.py            # ocean surface model (Trees & Stam 2019)
│   └── ckdistribution/     # correlated-k distribution pipeline
├── mie_source/             # Fortran: Mie scattering
├── dap_source/             # Fortran: adding-doubling polarised RT
├── geos_source/            # Fortran: disk geometry / output reader
├── readmie_source/         # Fortran: read Mie output files
├── mieshell_source/        # Fortran: coated-sphere Mie
├── *_sig.pyf               # f2py signature files (do not modify)
├── dap_database/           # pre-computed Fourier coefficient files
└── examples/               # reproduction scripts for published figures
    ├── trees_stam_fig7.py  # Disk-resolved ocean planet RGB images
    ├── trees_stam_fig3.py  # Phase curves vs surface pressure
    ├── trees_stam_fig1_single_wl.py
    ├── fresnel_only_fig1.py
    └── recreate_fig8.py    # Venus CO₂ 1.4 µm band (Mahapatra et al. 2024)
```
