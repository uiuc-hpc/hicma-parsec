cmake_policy(SET CMP0074 NEW)
include(FindPkgConfig)

include(ParsecCompilePTG)
include_directories("${PROJECT_SOURCE_DIR}/dplasma/parsec")
include_directories("${PROJECT_SOURCE_DIR}/dplasma")
include_directories("${PROJECT_SOURCE_DIR}/dplasma/src")
include_directories(BEFORE "${PROJECT_SOURCE_DIR}/dplasma/src/include")
include_directories("${PROJECT_SOURCE_DIR}/dplasma/src/cores")
include_directories("${PROJECT_SOURCE_DIR}/dplasma/parsec")
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}/../dplasma")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/../dplasma")
link_directories(${STARSH_LIBRARY_DIRS})
link_directories(${HCORE_LIBRARY_DIRS})

set( SRC
     descprint.c
     HiCMA_dpotrf_L_2flow_wrapper.c
     HiCMA_dpotrf_L_3flow_wrapper.c
     hicma_hcore.c
     hicma_init.c
     hicma_parsec.c
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
     parsec
     dplasma
     ${HICMA_DEP}
   )

# If MKL, make sure the sequential lib is used
string(FIND "${CBLAS_DIR_FOUND}" "mkl" _match_mkl)
if( _match_mkl )
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
