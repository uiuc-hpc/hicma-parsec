#!/bin/bash
#SBATCH --account=k1339
#SBATCH --job-name=compile-lorapo
#SBATCH --output=compile-lorapo.out
#SBATCH --error=compile-lorapo.out
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --partition=debug

module unload cray-libsci/17.12.1
typ=O3

export PKG_CONFIG_PATH=/lustre/project/k1205/lei/dplasma/build-$typ/install/lib/pkgconfig:/lustre/project/k1205/lei/dplasma/build-$typ/install/lib64/pkgconfig:/lustre/project/k1205/akbudak/codes/plasma-installer_2.8.0/install-intel/lib/pkgconfig:/lustre/project/k1205/akbudak/hicma-dev/stars-h/build/install/lib/pkgconfig:$PKG_CONFIG_PATH

export CRAYPE_LINK_TYPE=dynamic

make clean
srun make
