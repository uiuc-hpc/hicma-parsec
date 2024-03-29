# HICMA with PaRSEC

## Source codes

Suppose you will put the source codes under $HOME directory.

```
cd $HOME
git clone --recursive https://github.com/ecrc/hicma-x-dev.git
```

## External Dependencies

The following modules can be used on ECRC systems:

```
mkl/2018-update-1   gcc/5.5.0    cmake/3.17.3   openmpi/3.0.0-gcc-5.5.0
```

DPLASMA as well as PaRSEC which is a submodule in DPLASMA, HCORE and STARS-H are required. 
These libraries are provided as submodules of this repository 
so use these submodules for installation.

`git submodule update --init --recursive` can be used to get the submodules.

HCORE and STARS-H are manually installed as mentioned in the following subsections.
But DPLASMA and HiCMA are installed together using a single command
as will be mentoned in the next section.

### HCORE

A sample installation of HCORE:

```
cd $HOME/hicma-x-dev
cd hcore && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=`pwd`/installdir 
make -j install
```

Set `$PKG_CONFIG_FILE`:

```
export PKG_CONFIG_PATH=$HOME/hicma-x-dev/hcore/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
```

### STARS-H

STARS-H can be installed following the instructions at https://github.com/ecrc/stars-h. Make sure to export `PKG_CONFIG_PATH`.

A sample installation of STARS-H:

```
cd $HOME/hicma-x-dev
cd stars-h && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=`pwd`/installdir -DMPI=OFF -DSTARPU=OFF -DCUDA=OFF -DOPENMP=OFF -DGSL=OFF
make -j install
```

Set `$PKG_CONFIG_FILE`:

```
export PKG_CONFIG_PATH=$HOME/hicma-x-dev/stars-h/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
```

---

## Installation

Submodules must be updated via `git submodule update --init --recursive`.

```
cd $HOME/hicma-x-dev
mkdir -p build && cd build 
cmake .. 
```

In addition, if Intel compiler is used, add `-DCMAKE_Fortran_FLAGS="-nofor-main"`. 

A sample configuration for Shaheen II: 
```
cmake .. -DCMAKE_Fortran_FLAGS="-nofor-main" -DDPLASMA_PRECISIONS="d" -DBUILD_SHARED_LIBS=ON -DPARSEC_WITH_HEADER_FILES=ON -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DPARSEC_ATOMIC_USE_GCC_32_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_64_BUILTINS_EXITCODE=0 -DPARSEC_ATOMIC_USE_GCC_128_BUILTINS_EXITCODE=0 -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DCMAKE_BUILD_TYPE=Release -DBLA_VENDOR=Intel10_64lp_seq -DPARSEC_DIST_EAGER_LIMIT=0 -DPARSEC_DIST_SHORT_LIMIT=0
```

Run make for compilation. If following command fails, try removing `-j 8`

```
make -j 8
```

## Quick Example Runs

Go to HiCMA folder:

```
cd hicma_parsec
```

Running examples:

statistics-2d-sqexp:

```
mpirun -np 4 -npernode 1 ./testing_dpotrf_tlr -N 2700 -t 270 -e 1e-8 -u 130 -D 2 -P 2 -v
```

numerical correctness checking enabled:

```
mpirun -n 4 -npernode 1 ./testing_dpotrf_tlr -N 2700 -t 270 -e 1e-8 -u 130 -D 2 -P 2 -v --check
```

statistics-3d-exp:

```
mpirun -np 4 --npernode 1 ./testing_dpotrf_tlr -N 108000 -t 2700 -e 1e-8 -u 1200 -D 4 -P 2 -v -- -mca runtime_comm_coll_bcast 0
```

3-flow version:

```
mpirun -n 4 ./testing_dpotrf_tlr -N 2700 -t 270 -e 1e-8 -u 130 -D 2 -P 2 -v -E 0 -Z 1 --check
```

## Program Parameters

```
-N: matrix size; required 
-t: tile size; required
-e: accuracy threshold; default: 1.0e-8
-u: maxrank threshold for compressed tiles; default: tile_size/2
-P: row process grid; default: number_of_nodes
-D: kind of problem: default: 2
-v: print more info
```

More information:

```
./testing_dpotrf_tlr --help
```

Additional PaRSEC flags:
```
./testing_dpotrf_tlr -- --help
```



## Tips 

(1) if the problem is a little dense, i.e., band_size > 1 after auto-tuning (e.g., in statistics-3d-sqexp application with accuracy threshold -e 1.0e-8), "-- -mca runtime_comm_coll_bcast 0" is needed for better performance;

(2) Set argument -c to number_of_cores - 1;

(3) Choose the process grid to be as square as possible with P < Q;

(4) in most cases,

```
    for -D 2 (statistics-2d-sqexp), set maxrank= 150;
    for -D 3 (statistics-3d-sqexp), set maxrank= 500;
    for -D 4 (statistics-3d-exp), set maxrank= tile_size / 2. 
```

References
-----------
1. Q. Cao, Y. Pei, T. Herault, K. Akbudak, A. Mikhalev, G. Bosilca, H. Ltaief, D. E. Keyes, and J. Dongarra, *Performance Analysis of Tile Low-Rank Cholesky Factorization Using PaRSEC Instrumentation Tools*, **2019 IEEE/ACM International Workshop on Programming and Performance Visualization Tools (ProTools)**, Denver, CO, USA, 2019, pp. 25-32.

2. Q. Cao, Y. Pei, K. Akbudak, A. Mikhalev, G. Bosilca, H. Ltaief, D. E. Keyes, and J. Dongarra, *Extreme-Scale Task-Based Cholesky Factorization Toward Climate and Weather Prediction Applications*, **The Platform for Advanced Scientific Computing (PASC 2020)**.

3. Q. Cao, Y. Pei, K. Akbudak, G. Bosilca, H. Ltaief, D. E. Keyes, and J. Dongarra, *Leveraging PaRSEC Runtime Support to Tackle Challenging 3D Data-Sparse Matrix Problems*, **IEEE International Parallel & Distributed Processing Symposium (IPDPS 2021)**.
