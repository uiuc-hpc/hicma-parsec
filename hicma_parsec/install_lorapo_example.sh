#/bin/bash

git_dplasma=1
git_parsec=1
git_starsh=1
git_lorapo=1
install_dplasma=1
install_starsh=1
install_lorapo=1

root=$PWD
echo "root = $root"
if [ ! -d $root ]; then
    echo "Root folder $root does not exist."
    mkdir $root
fi
cd $root

module list -l

echo "Downloading lorapo"
cd  $root
if [ $git_lorapo -eq 1 ]; then git clone https://Qinglei_Cao@bitbucket.org/bosilca/lorapo.git; fi
if [ $install_lorapo -eq 1 ]; then
    cd lorapo
    git checkout lei_2flow

    #git checkout lei_reorder_gemm_testing 
    git reset --hard HEAD
    echo "\
            CC = mpicc
            LD = mpicc
            #NVCC = nvcc

            BUILD_TYPE = build

            LDFLAGS = -lmkl_intel_lp64 -lmkl_sequential -lgsl 

            DPLASMA_SRC = $root/dplasma
            PARSEC_SRC = \${DPLASMA_SRC}/parsec
            DPLASMA_INSTALL = \${DPLASMA_SRC}/\${BUILD_TYPE}/install" > make.inc
fi

cd $root

echo "Install DPLASMA"
if [ $git_dplasma -eq 1 ]; then git clone https://bitbucket.org/icldistcomp/dplasma.git; fi
cd dplasma
if [ $git_dplasma -eq 1 ]; then 
    git reset --hard 2fbef5d868f51f08300231045f1310e469821bb0
    sed -i '31i#include \"dplasma_complex.h\"' ./src/cores/core_zblas.h
fi
if [ $git_parsec -eq 1 ]; then 
    rm -rf parsec
    git clone  https://bitbucket.org/icldistcomp/parsec.git; 
fi
cd parsec
if [ $git_parsec -eq 1 ]; then 
    git reset --hard 42a44ebdc37ea789aff17ea5511ffb3a86a22fb8
    #echo "#define PARSEC_KCYCLIC_WITH_VIEW 1" >> parsec/include/parsec/parsec_config_bottom.h 
    cp $root/lorapo/lorapo_parsec.patch .
    git apply lorapo_parsec.patch
fi
if [ $install_dplasma -eq 1 ];then
    cd $root/dplasma
    mkdir -p build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/install -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0 -DPARSEC_DIST_COLLECTIVES=OFF -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DDPLASMA_PRECISIONS="d" -DPARSEC_WITH_DEVEL_HEADERS=ON
    make -j 8
    make -j 8 install
fi
echo "export PKG_CONFIG_PATH=$root/dplasma/build/install/lib/pkgconfig:\$PKG_CONFIG_PATH" > $root/env_lorapo.sh
echo "export PKG_CONFIG_PATH=$root/dplasma/build/install/lib64/pkgconfig:\$PKG_CONFIG_PATH" >> $root/env_lorapo.sh
echo "PaRSEC_ROOT=$root/dplasma/build/install" >> $root/env_lorapo.sh
export PKG_CONFIG_PATH=$root/dplasma/build/install/lib/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=$root/dplasma/build/install/lib64/pkgconfig:$PKG_CONFIG_PATH
export PaRSEC_ROOT=$root/dplasma/build/install

echo "Install starsh"
cd $root
if [ $git_starsh -eq 1 ]; then git clone https://github.com/ecrc/stars-h.git; fi
if [ $install_starsh -eq 1 ]; then
    cd stars-h/
    git submodule update --init
    mkdir -p build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/install
    make -j 8
    make -j 8 install
fi
echo "export PKG_CONFIG_PATH=$root/stars-h/build/install/lib/pkgconfig:\$PKG_CONFIG_PATH" >> $root/env_lorapo.sh
export PKG_CONFIG_PATH=$root/stars-h/build/install/lib/pkgconfig:$PKG_CONFIG_PATH

echo "Install lorapo"
cd  $root
if [ $install_lorapo -eq 1 ]; then 
    cd lorapo
    make
fi

# NOTE
# Environment parameters are written in env_lorapo.sh
# . env_lorapo.sh before running experiments

# Running Exemple
./testing_dpotrf -N 27000 -t 2700 -e 1e-8 -u 200 -j 27000 -v -c 19 -G 200 -U 200 -D 2 -z 30 -Z 10 -Y 1
# mpirun -np 4 --npernode 1 ./testing_dpotrf -N 108000 -t 2700 -e 1e-8 -u 1200 -j 108000 -v -P 2 -c 19 -G 1200 -U 1200 -D 4 -z 300 -Z 6 -Y 1 -I 0

# More details ./testing_dpotrf --help
