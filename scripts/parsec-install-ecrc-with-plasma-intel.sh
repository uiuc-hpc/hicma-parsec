#cmake /home/akbudak/parsec  -DBUILD_SHARED_LIBS=OFF -DPARSEC_WITH_HEADER_FILES=ON -DCMAKE_C_FLAGS="-fPIC -DADD_" -DPARSEC_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=/home/akbudak/lib/parsec -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=ifort -DPARSEC_PROF_GRAPHER=ON -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0 
#make -j install
#export PATH=/home/akbudak/lib/parsec/bin:$PATH
#export PKG_CONFIG_PATH=/home/akbudak/lib/parsec/lib/pkgconfig:$PKG_CONFIG_PATH

scriptdir=`dirname ${BASH_SOURCE[0]}`
INSTALLDIR=/home/akbudak/lib/parsec_plasma_intel
cdir=$PWD
cd $HOME
if [ ! -d parsec ]; then
    git clone git@bitbucket.org:icldistcomp/parsec.git
fi
if [ ! -d parsec ]; then
    echo "parsec repo cannot be cloned. Exiting"
    exit 1
fi
if [ -d $INSTALLDIR ]; then
    echo "removing previous installation"
    rm -rf $INSTALLDIR
fi
cd parsec
if [ -d build ]; then
    echo "removing previous build dir"
    rm -rf build
fi
mkdir build
cd build

module purge
module load ecrc-extras
module load cmake/3.11.1
module load intel/2018
#. /opt/share/intel/2018/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64

module load hwloc/1.11.8-intel-2018
module load plasma/2.8.0-intel-2018-mkl

cmake ..  -DBUILD_SHARED_LIBS=OFF -DPARSEC_WITH_HEADER_FILES=ON -DCMAKE_C_FLAGS="-fPIC -DADD_" -DPARSEC_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=$INSTALLDIR -DPARSEC_WITH_DEVEL_HEADERS=ON  -DPARSEC_PROF_GRAPHER=ON -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0 -DBUILD_TESTING=ON -DDPLASMA_MPI=ON -DBUILD_DPLASMA=ON
#EDUARDO 
#cmake ..  -DBUILD_SHARED_LIBS=OFF -DPARSEC_WITH_HEADER_FILES=ON -DCMAKE_INSTALL_PREFIX=$INSTALLDIR -DPARSEC_WITH_DEVEL_HEADERS=ON  
#cmake ..  -DBUILD_SHARED_LIBS=OFF -DPARSEC_WITH_HEADER_FILES=ON -DCMAKE_C_FLAGS="-fPIC -DADD_" -DPARSEC_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=$INSTALLDIR -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=ifort -DPARSEC_PROF_GRAPHER=ON -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0 -DBUILD_TESTING=ON -DDPLASMA_MPI=ON -DBUILD_DPLASMA=ON
make -j install
export PATH=$INSTALLDIR/bin:$PATH
export PARSEC_ROOT=$INSTALLDIR
export PKG_CONFIG_PATH=$INSTALLDIR/lib/pkgconfig:$PKG_CONFIG_PATH
cp $scriptdir/make.inc-plasma-intel $scriptdir/../make.inc 
cd $cdir
