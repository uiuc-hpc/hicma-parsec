#!/bin/bash
#SBATCH --account=k1205
#SBATCH --job-name=plasma-install
#SBATCH --output=plasma-install.out
#SBATCH --error=plasma-install.err
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#module switch PrgEnv-cray PrgEnv-gnu

#copy this script to /lustre/project/k1205/akbudak/codes/plasma-installer_2.8.0/plasma-install.sh

env=intel
installdir=/lustre/project/k1205/akbudak/codes/plasma-installer_2.8.0/install-$env
module switch PrgEnv-cray PrgEnv-$env

srun ./setup.py --clean; 
#srun ./setup.py --cc="cc" --fc="ftn" --blaslib="-L/opt/intel/composer_xe_2015.2.164/mkl/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lpthread" --lapacklib="-L/opt/intel/composer_xe_2015.2.164/mkl/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lpthread" --ldflags_c=-fopenmp --ldflags_c=-fopenmp --downall
#srun ./setup.py --cc="cc -shared -fPIC" --fc="ftn -shared -fPIC" --blaslib=" " --lapacklib=" " --ldflags_c=-fopenmp --ldflags_c=-fopenmp --downall
srun ./setup.py --prefix=$installdir --cc="cc -fPIC" --fc="ftn -fPIC" --blaslib=" " --lapacklib=" " --ldflags_c=-fopenmp --ldflags_c=-fopenmp --downall --notesting
