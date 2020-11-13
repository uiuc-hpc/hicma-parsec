#!/bin/bash
#SBATCH --account=k1339
##SBATCH --account=k1205
#SBATCH --job-name=dplasma-install
#SBATCH --output=dplasma-install.out
#SBATCH --error=dplasma-install.out
#SBATCH --nodes=1
##SBATCH --time=00:30:00
##SBATCH --partition=debug

#copy this script to /lustre/project/k1205/akbudak/codes/parsec/dplasma-install.sh

#module load PrgEnv-intel
#module load PrgEnv-gnu
env=intel
installdir=/lustre/project/k1205/akbudak/codes/parsec/build-$env/install
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/lustre/project/k1205/akbudak/codes/plasma-installer_2.8.0/install-$env/lib/pkgconfig
module swap PrgEnv-cray PrgEnv-$env
module unload cray-libsci/17.12.1

#module load PrgEnv-intel/5.2.82
#ftn --version

export CRAYPE_LINK_TYPE=dynamic
#srun cmake .. -DPARSEC_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=/project/k1124/ltaiefh/codes/astronomy/build-gnu/install -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn
#srun cmake .. -DPARSEC_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=/project/k1124/ltaiefh/codes/astronomy/build-gnu/install -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DPARSEC_ATOMIC_USE_GCC_32_BUILTINS=1 -DPARSEC_ATOMIC_USE_GCC_64_BUILTINS=1 -DPARSEC_ATOMIC_USE_GCC_128_BUILTINS=1

#srun cmake .. -DPARSEC_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=/project/k1124/ltaiefh/codes/astronomy/build-intel/install -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DMPI_INCLUDE_PATH=/opt/cray/mpt/7.2.6/gni/mpich-intel/14.0/include -DMPI_LIBRARIES=/opt/cray/mpt/7.2.6/gni/mpich-intel/14.0/lib -DMPI_Fortran_LIBRARIES=mpichf90_intel -DMPI_C_LIBRARIES=mpich_intel -DMPI_CXX_LIBRARIES=mpichcxx_intel -DMPI_Fortran_INCLUDE_PATH=/opt/cray/mpt/7.2.6/gni/mpich-intel/14.0/include



cd /lustre/project/k1205/akbudak/codes/parsec/build-$env
#srun cmake /lustre/project/k1205/akbudak/codes/parsec  -DCMAKE_LINKER=cc -DBUILD_SHARED_LIBS=OFF -DPARSEC_WITH_HEADER_FILES=ON -DCMAKE_C_FLAGS="-fPIC -craympich-mt" -DPARSEC_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=$installdir -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DPARSEC_ATOMIC_USE_GCC_32_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_64_BUILTINS_EXITCODE=0 \
#    -DPARSEC_ATOMIC_USE_GCC_128_BUILTINS_EXITCODE=0 -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0 

#-DBUILD_DPLASMA=OFF  -DBUILD_TESTING=OFF


srun make  install
