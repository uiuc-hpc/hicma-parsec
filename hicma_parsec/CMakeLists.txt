cmake_policy(SET CMP0074 NEW)
include(FindPkgConfig)

OPTION(LINK_WITH_MKL "Link with MKL" ON)

include(ParsecCompilePTG)
include_directories("${PROJECT_SOURCE_DIR}/dplasma/parsec")
include_directories("${PROJECT_SOURCE_DIR}/dplasma")
include_directories("${PROJECT_SOURCE_DIR}/dplasma/src")
include_directories(BEFORE "${PROJECT_SOURCE_DIR}/dplasma/src/include")
include_directories("${PROJECT_SOURCE_DIR}/dplasma/src/cores")
include_directories("${PROJECT_SOURCE_DIR}/dplasma/parsec")
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}")
link_directories(${STARSH_LIBRARY_DIRS})

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

set( CMAKE_C_FLAGS 
  ${CMAKE_C_FLAGS}
  -D_GNU_SOURCE
   )

add_executable(testing_dpotrf_tlr ${SRC} testing_dpotrf.c)
if( LINK_WITH_MKL ) 
    target_link_libraries(testing_dpotrf_tlr mkl_intel_lp64 mkl_sequential parsec dplasma lapacke ${HICMA_DEP})
else( LINK_WITH_MKL )
    target_link_libraries(testing_dpotrf_tlr parsec dplasma lapacke ${HICMA_DEP})
endif( LINK_WITH_MKL )

target_ptg_sources(testing_dpotrf_tlr PRIVATE ${JDF} )
