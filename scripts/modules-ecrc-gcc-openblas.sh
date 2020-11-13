module load ecrc-extras
module load cmake
module load gcc/7.2.0 openmpi/3.0.0-gcc-7.2.0
module load openblas/0.2.20-gcc-7.2.0
module load plasma/2.8.0-gcc-7.2.0-openblas hwloc/1.11.8-gcc-7.2.0
module load parsec/master-gcc-7.2.0-openblas-openmpi-plasma-2.8.0
#If ecrc modules are used, export path of DPLASMA inside PARSEC
export PKG_CONFIG_PATH=${PARSEC_ROOT}/dplasma/lib/pkgconfig:$PKG_CONFIG_PATH
