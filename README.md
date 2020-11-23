It depends on dplasma ans star-sh

First, install star-sh, following https://github.com/ecrc/stars-h. Make sure to export PKG_CONFIG_PATH.

git clone https://github.com/ecrc/hicma-x-dev.git

cd hicma-x-dev && mkdir -p build

git checkout band_tlr_pasc 

cd build && cmake .. -DPARSEC_GPU_WITH_CUDA=OFF -DDPLASMA_PRECISIONS="d" -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0

In addtion, (1) if gonna to run a more dense case, e.g., 3d-exp, add -DPARSEC_DIST_COLLECTIVES=OFF; (2) intel compiler, add -DCMAKE_Fortran_FLAGS="-nofor-main".

make -j 8

cd hicma_parsec

./testing_dpotrf --help

For instance,

mpirun -np 4 -npernode 1 ./testing_dpotrf_tlr -N 2700 -t 270 -e 1e-8 -u 130 -j 2700 -v -c 19 -G 130 -U 130 -D 2 -z 30 -Z 1 -Y 1 -E 0 -P 2 -x -W 1

mpirun -np 4 --npernode 1 ./testing_dpotrf_tlr -N 108000 -t 2700 -e 1e-8 -u 1200 -j 108000 -v -P 2 -c 19 -G 1200 -U 1200 -D 4 -z 300 -Z 6 -Y 1 -I 0
