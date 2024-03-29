# @copyright (c) 2021     King Abdullah University of Science and Technology (KAUST).
# @copyright (c) 2021     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
#                         All rights reserved.
# @version 0.1.0
# @date 2021-01-24

cmake_policy(SET CMP0074 NEW)
include(FindPkgConfig)

include(ParsecCompilePTG)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set( SRC
     descprint.c
     HiCMA_dpotrf_L_2flow_wrapper.c
     HiCMA_dpotrf_L_3flow_wrapper.c
     hicma_internal.c
     hicma_init.c
     hicma_parsec.c
     hicma_parsec_distribution.c
   )

set( JDF
     Av_memory.jdf
     band_free.jdf
     band_regenerate.jdf
     band_size_auto_tuning.jdf
     HiCMA_dpotrf_L_2flow.jdf
     HiCMA_dpotrf_L_3flow.jdf
     rank_gather.jdf
     rank_check.jdf
     rank_print.jdf
     STARSH_check.jdf
     STARSH_gen.jdf
   )

#git rev-parse --short HEAD
set( CMAKE_C_FLAGS ${CMAKE_C_FLAGS} " -DGITHASH=8cdc2d9 -D_GNU_SOURCE ")

set( LD_LIB_HICMA_PARSEC
     ${PARSEC_LIBRARIES}
     ${DPLASMA_LIBRARIES}
     ${HICMA_DEP}
   )

# If MKL, make sure the sequential lib is used
if(CBLAS_DIR_FOUND MATCHES "mkl")
    set( LD_LIB_HICMA_PARSEC
         mkl_intel_lp64
         mkl_sequential
         ${LD_LIB_HICMA_PARSEC} 
       )
    message(STATUS "Linked with MKL Library: ${LD_LIB_HICMA_PARSEC}")
endif( )

add_executable(testing_dpotrf_tlr ${SRC} testing_dpotrf.c)
target_link_libraries(testing_dpotrf_tlr ${LD_LIB_HICMA_PARSEC})

target_ptg_sources(testing_dpotrf_tlr PRIVATE ${JDF} )

add_test(NAME "dpotrf-3flow"
    COMMAND make -j && mpirun -n 4 ./testing_dpotrf_tlr -N 2700 -t 270 -e 1e-8 -u 130 -D 2 -P 2 -v -E 0 -Z 1 --check)
add_test(NAME "dpotrf-2flow"
    COMMAND make -j && mpirun -n 4 ./testing_dpotrf_tlr -N 2700 -t 270 -e 1e-8 -u 130 -D 2 -P 2 -v --check)

include(GNUInstallDirs)
install(TARGETS testing_dpotrf_tlr
        DESTINATION ${CMAKE_INSTALL_LIBEXECDIR}/hicma)
