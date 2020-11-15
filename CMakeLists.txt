cmake_minimum_required (VERSION 3.16)
project (DTRSM C)

include(CMakeDependentOption)
include(CMakePushCheckState)
include(CheckSymbolExists)
include(CheckCSourceCompiles)
include(CheckIncludeFiles)
include(GNUInstallDirs)

# The current version number
set (DTRSM_VERSION_MAJOR 2)
set (DTRSM_VERSION_MINOR 0)

SET(DTRSM_VERSION "${DTRSM_VERSION_MAJOR}.${DTRSM_VERSION_MINOR}")

set(CMAKE_NO_SYSTEM_FROM_IMPORTED True)

# CMake Policies Tuning
if(POLICY CMP0074)
  # CMP0074: Starting with CMake 3.12, all FIND_<something> use <something>_ROOT in the search path
  #          in addition to the specified paths
  cmake_policy(SET CMP0074 NEW)
ENDIF(POLICY CMP0074)
if(POLICY CMP0094)
  # CMP0094: Starting with CMake 3.16, all FIND_Python will use the first matching version instead of
  #          of searching for the largest available version number (which defeats our selection logic)
  cmake_policy(SET CMP0094 NEW)
ENDIF(POLICY CMP0094)


option(BUILD_64bits
  "Build 64 bits mode" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

option(BUILD_SHARED_LIBS
       "Enable building shared DPLASMA lib (default ON)" ON)


option(SUPPORT_C11
       "Enable support for C11 capabilities. Might not work with OpenMP support due to an OpenMP compilers restriction (default ON)." ON)
if( SUPPORT_C11 )
  set(CMAKE_C_STANDARD 11)
else (SUPPORT_C11)
  set(CMAKE_C_STANDARD 99)
endif( SUPPORT_C11 )
set(CMAKE_C_STANDARD_REQUIRED ON)


STRING(COMPARE EQUAL ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} DTRSM_BUILD_INPLACE)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/dplasma/cmake_modules/")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/dplasma/parsec/cmake_modules/")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/dplasma/tools/PrecisionGenerator/")

#set(BLA_VENDOR CACHE STRING "Possible valid values are: Intel10_64lp_seq, IBMESSL, Apple, or All, which will default to something like Goto or ATLAS")
# Goto, ATLAS, PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_64lp_seq, ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic
#set_property(CACHE BLA_VENDOR PROPERTY STRINGS Intel10_64lp_seq IBMESSL Apple Goto ATLAS)
#find_package(LAPACKE REQUIRED)

# we have to build it
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/dplasma/CMakeLists.txt")
    find_package(Git QUIET)
    if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
        # Update submodules as needed
        option(GIT_SUBMODULE "Check submodules during build" ON)
        if(GIT_SUBMODULE)
            message(STATUS "Submodule update")
            execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --force --init --recursive --remote
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                RESULT_VARIABLE GIT_SUBMOD_RESULT)
            if(NOT GIT_SUBMOD_RESULT EQUAL "0")
                message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
            endif()
        endif()
    endif()
endif(NOT EXISTS "${PROJECT_SOURCE_DIR}/dplasma/CMakeLists.txt")

#option(PARSEC_WITH_DEVEL_HEADERS "" ON)
add_subdirectory(dplasma)

# cmake config files
#include(CMakePackageConfigHelpers)


include(DplasmaCompilerFlags)
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}/include")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/dplasma/tests")
#Add_Subdirectory(include)
#Add_Subdirectory(cores)

include(ParsecCompilePTG)

### Generate .c files from .jdf for all required precisions
set(generated_jdf
    dtrsm_LLN.jdf dtrsm_LLT.jdf
    )

#precisions_rules_py(generated_jdf
#    "${JDF}"
#    PRECISIONS "${DPLASMA_PRECISIONS}")

set(generated_wrappers
    dtrsm_wrapper.c
    )

add_executable(dtrsm testing_dtrsm.c common.c common.h)
target_link_libraries(dtrsm PRIVATE dplasma )

target_ptg_sources(dtrsm PRIVATE
    ${generated_jdf})

target_include_directories(dtrsm
    INTERFACE
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/dplasma/src>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/dplasma/tests>
    ${PROJECT_SOURCE_DIR}/dplasma/tests
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tests>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/include>
    $<$<NOT:${DPLASMA_BUILD_INPLACE}>:
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src/include>>)

set_target_properties(dtrsm PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

target_link_libraries(dtrsm
    PRIVATE
    dplasma_cores
    PUBLIC
    PaRSEC::parsec
    LAPACKE::LAPACKE
    $<$<BOOL:${DPLASMA_HAVE_CUDA}>:${CUDA_CUBLAS_LIBRARIES}>)
#set_target_properties(dtrsm PROPERTIES VERSION ${DPLASMA_VERSION_MAJOR}.${DPLASMA_VERSION_MINOR}
#    SOVERSION ${DPLASMA_VERSION_MAJOR})

