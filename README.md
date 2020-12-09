It depends on dplasma ans star-sh

First, install star-sh, following https://github.com/ecrc/stars-h. Make sure to export PKG_CONFIG_PATH.

git clone https://github.com/ecrc/hicma-x-dev.git

cd hicma-x-dev && mkdir -p build

git checkout band_tlr_pasc 

cd build && cmake .. 

In addtion, if intel compiler, add -DCMAKE_Fortran_FLAGS="-nofor-main". For instance, configuration on Shaheen II: cmake .. -DCMAKE_Fortran_FLAGS="-nofor-main" -DDPLASMA_PRECISIONS="d" -DBUILD_SHARED_LIBS=ON -DPARSEC_WITH_HEADER_FILES=ON -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DPARSEC_ATOMIC_USE_GCC_32_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_64_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_128_BUILTINS_EXITCODE=0 -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DCMAKE_BUILD_TYPE=Release -DBLA_VENDOR=Intel10_64lp_seq -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0

make -j 8

cd hicma_parsec

Running examples:

statistics-2d-sqexp:

mpirun -np 4 -npernode 1 ./testing_dpotrf_tlr -N 2700 -t 270 -e 1e-8 -u 130 -D 2 -P 2 -v

statistics-3d-exp:

mpirun -np 4 --npernode 1 ./testing_dpotrf_tlr -N 108000 -t 2700 -e 1e-8 -u 1200 -D 4 -P 2 -v

-N: matrix size; required 

-t: tile size; required

-e: accuracy threshold; default: 1.0e-8

-u: maxrank threshold for compressed tiles; default: tile_size/2

-P: row process grid; default: number_of_nodes

-D: kind of problem: default: 2

-v: print more info

More information:

./testing_dpotrf_tlr --help

Additional PaRSEC flags: ./testing_dpotrf_tlr -- --help



Tips: 

(1) Set argument -c to number_of_cores - 1;

(2) Choose the process grid to be as square as possible with P < Q;

(3) in most cases,
    for -D 2 (statistics-2d-sqexp), set maxrank= 150;
    for -D 3 (statistics-3d-sqexp), set maxrank= 500;
    for -D 4 (statistics-3d-exp), set maxrank= tile_size / 2. 
