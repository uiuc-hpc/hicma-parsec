It depends on dplasma ans star-sh

First, install star-sh, following https://github.com/ecrc/stars-h. Make sure to export PKG_CONFIG_PATH.

git clone https://github.com/ecrc/hicma-x-dev.git

cd hicma-x-dev && mkdir -p build

git checkout band_tlr_pasc 

cd build && cmake .. 

In addtion, if intel compiler, add -DCMAKE_Fortran_FLAGS="-nofor-main". For instance, configartion on Shaheen II: cmake .. -DCMAKE_Fortran_FLAGS="-nofor-main" -DDPLASMA_PRECISIONS="d" -DBUILD_SHARED_LIBS=ON -DPARSEC_WITH_HEADER_FILES=ON -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DPARSEC_ATOMIC_USE_GCC_32_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_64_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_128_BUILTINS_EXITCODE=0 -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DCMAKE_BUILD_TYPE=Release -DBLA_VENDOR=Intel10_64lp_seq -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0

make -j 8

cd hicma_parsec

./testing_dpotrf_tlr --help

Additional PaRSEC flags: ./testing_dpotrf_tlr -- --help

Running examples:

mpirun -np 4 -npernode 1 ./testing_dpotrf_tlr -N 2700 -t 270 -e 1e-8 -u 130 -j 2700 -v -c 19 -G 130 -U 130 -D 2 -z 30 -Z 1 -Y 1 -E 0 -P 2 -x

mpirun -np 4 --npernode 1 ./testing_dpotrf_tlr -N 108000 -t 2700 -e 1e-8 -u 1200 -j 108000 -v -P 2 -c 19 -G 1200 -U 1200 -D 4 -z 300 -Z 1 -Y 1 -- -mca runtime_comm_coll_bcast 0

Note: if the problem is a little dense, i.e., band_size > 1 after auto-tuning, "-- -mca runtime_comm_coll_bcast 0" is needed for better performance. 
