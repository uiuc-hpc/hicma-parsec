#!/bin/bash
#SBATCH --account=k1339
#SBATCH --job-name=hicma-parsec-install
#SBATCH --output=hicma-parsec-install.out
#SBATCH --error=hicma-parsec-install.out
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --partition=debug

# Intel env
module swap PrgEnv-cray PrgEnv-intel

# stars-h and gsl
export PKG_CONFIG_PATH=/lustre/project/k1205/lei/stars-h/build/install/lib/pkgconfig:/sw/xc40cle7/gsl/1.14/cle7_intel19.0.1/install/lib/pkgconfig:$PKG_CONFIG_PATH

# cmake 3.18 
export PATH=/lustre/project/k1205/lei/software/CMake/bin:$PATH

# cray link type
export CRAYPE_LINK_TYPE=dynamic

# cray complier
srun cmake .. -DCMAKE_Fortran_FLAGS="-nofor-main" -DDPLASMA_PRECISIONS="d" -DBUILD_SHARED_LIBS=ON -DPARSEC_WITH_HEADER_FILES=ON -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DPARSEC_ATOMIC_USE_GCC_32_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_64_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_128_BUILTINS_EXITCODE=0 -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DCMAKE_BUILD_TYPE=Release -DBLA_VENDOR=Intel10_64lp_seq -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0 -DPARSEC_DIST_COLLECTIVES=ON

make clean
srun make -j 10
