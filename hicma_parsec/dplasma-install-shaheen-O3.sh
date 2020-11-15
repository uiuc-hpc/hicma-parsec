#!/bin/bash
#SBATCH --account=k1205
#SBATCH --job-name=dplasma-install
#SBATCH --output=dplasma-install.out
#SBATCH --error=dplasma-install.err
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --partition=debug

#module unload PrgEnv-intel
#module load PrgEnv-gnu
env=intel
type=O3
installdir=/lustre/project/k1205/lei/dplasma/build-$type/install
export PKG_CONFIG_PATH=/lustre/project/k1205/akbudak/codes/plasma-installer_2.8.0/install-$env/lib/pkgconfig:$PKG_CONFIG_PATH
#module unload cray-libsci/17.12.1
#module unload perftools-base/7.1.3
#module load papi
#module load cray-python/2.7.15.7

#source /lustre/project/k1205/peiy/env/bin/activate

#module load PrgEnv-intel/5.2.82
#ftn --version

export CRAYPE_LINK_TYPE=dynamic

# cray complier
srun cmake .. -DCMAKE_Fortran_FLAGS="-nofor-main" -DDPLASMA_PRECISIONS="s;d" -DBUILD_SHARED_LIBS=ON -DPARSEC_WITH_HEADER_FILES=ON -DCMAKE_INSTALL_PREFIX=$installdir -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DPARSEC_ATOMIC_USE_GCC_32_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_64_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_128_BUILTINS_EXITCODE=0 -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DCMAKE_BUILD_TYPE=Release -DBLA_VENDOR=Intel10_64lp_seq -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0 -DPARSEC_DIST_COLLECTIVES=OFF

make clean
srun make -j 10 install
