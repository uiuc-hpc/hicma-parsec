module switch PrgEnv-cray PrgEnv-intel
module unload cray-libsci/17.12.1
export PARSEC_ROOT=/lustre/project/k1205/akbudak/codes/parsec/build-intel/install
export PKG_CONFIG_PATH=/lustre/project/k1205/akbudak/codes/plasma-installer_2.8.0/install-intel/lib/pkgconfig/:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=/lustre/project/k1205/akbudak/codes/parsec/build-intel/install/lib/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=/lustre/project/k1205/akbudak/codes/parsec/build-intel/install/dplasma/lib/pkgconfig/:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=/lustre/project/k1205/akbudak/hicma-dev/stars-h/build/install/lib/pkgconfig:$PKG_CONFIG_PATH
