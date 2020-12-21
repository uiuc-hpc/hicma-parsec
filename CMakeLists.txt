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
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")

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

include(GenPkgConfig)

# HICMA depends on STARSH
# -------------------------------
find_package( STARSH REQUIRED )
if( STARSH_FOUND )
    include_directories(${STARSH_INCLUDE_DIRS_DEP})
    if(STARSH_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${STARSH_LINKER_FLAGS}")
    endif()
    if(STARSH_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${STARSH_LIBRARY_DIRS}")
    endif()
    if (STARSH_LIBRARIES)
        # look for gsl
        find_library( _STARSH_LIB NAME starsh PATHS ${STARSH_LIBRARY_DIRS} )
        if( _STARSH_LIB AND NOT "${STARSH_LIBRARIES_DEP}" MATCHES "gsl" )
            execute_process( COMMAND nm ${_STARSH_LIB} COMMAND grep -q gsl RESULT_VARIABLE GSL_IN_STARSH )
            if ( ${GSL_IN_STARSH} EQUAL 0)
                message( STATUS "STARSH depends on gsl. Adding it to dependency list")
                find_package(GSL REQUIRED)
                if (STARSH_LIBRARIES_DEP)
                    list( APPEND STARSH_LIBRARIES_DEP ${GSL_LIBRARIES})
                else()
                    list( APPEND STARSH_LIBRARIES ${GSL_LIBRARIES})
                endif()
            endif()
        endif()
        # insert to dependencies
        if (STARSH_LIBRARIES_DEP)
            list(APPEND HICMA_DEP ${STARSH_LIBRARIES_DEP})
        else()
            list(APPEND HICMA_DEP ${STARSH_LIBRARIES})
        endif()
    endif()
endif()

# HICMA depends on HCORE
# -------------------------------
find_package( HCORE REQUIRED )
if( HCORE_FOUND )
    include_directories(${HCORE_INCLUDE_DIRS_DEP})
    if(HCORE_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${HCORE_LINKER_FLAGS}")
    endif()
    if(HCORE_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${HCORE_LIBRARY_DIRS}")
    endif()
    if (HCORE_LIBRARIES)
        find_library( _HCORE_LIB NAME hcore PATHS ${HCORE_LIBRARY_DIRS} )
        list(APPEND HICMA_DEP ${HCORE_LIBRARIES})
    endif()
endif()

add_subdirectory(hicma_parsec)
