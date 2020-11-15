# This file is a preliminary version of a script to setup and install plasma and parsec.
# Please do not run this script directly
exit 0
wget http://icl.cs.utk.edu/projectsfiles/plasma/pubs/plasma-installer_2.8.0.tar.gz

export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/home/akbudak/lib/plasma/lib/pkgconfig


/setup.py --clean;./setup.py --cc=icc --fc=ifort --blaslib=" -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -lirc -lsvml" --lapacklib="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -lirc -lsvml" --downall --notesting --prefix=/home/akbudak/lib/plasma

git clone git@bitbucket.org:icldistcomp/parsec.git
mkdir build
cd build

#this worked for generating dot file
cmake /home/akbudak/parsec  -DBUILD_SHARED_LIBS=OFF -DPARSEC_WITH_HEADER_FILES=ON -DCMAKE_C_FLAGS="-fPIC" -DPARSEC_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=/home/akbudak/lib/parsec -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=ifort -DPARSEC_PROF_GRAPHER=ON -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0


cmake /home/akbudak/parsec  -DBUILD_SHARED_LIBS=OFF -DPARSEC_WITH_HEADER_FILES=ON -DCMAKE_C_FLAGS="-fPIC" -DPARSEC_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=/home/akbudak/lib/parsec -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=ifort -DPARSEC_ATOMIC_USE_GCC_32_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_64_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_128_BUILTINS_EXITCODE=0 \
        -DPARSEC_PROF_GRAPHER=ON -DPARSEC_DIST_EAGER_LIMIT=0 \
            -DPARSEC_DIST_SHORT_LIMIT=0 -DPARSEC_PROF_TRACE=ON 


                -DPARSEC_PROF_TRACE_SYSTEM=OTF2

cmake ..  -DBUILD_SHARED_LIBS=OFF -DPARSEC_WITH_HEADER_FILES=ON -DCMAKE_C_FLAGS="-fPIC" -DPARSEC_GPU_WITH_CUDA=OFF -DCMAKE_INSTALL_PREFIX=/home/akbudak/lib/parsec -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort -DPARSEC_ATOMIC_USE_GCC_32_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_64_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_128_BUILTINS_EXITCODE=0 \
    -DPARSEC_PROF_GRAPHER=ON -DPARSEC_DIST_EAGER_LIMIT=0 \
    -DPARSEC_DIST_SHORT_LIMIT=0 -DPARSEC_PROF_TRACE=ON \
    -DPARSEC_PROF_TRACE_SYSTEM=OTF2


export PKG_CONFIG_PATH=/home/akbudak/lib/parsec/dplasma/lib/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=/home/akbudak/lib/parsec/lib/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=/home/akbudak/hicma-dev/stars-h/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
