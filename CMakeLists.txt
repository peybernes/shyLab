# -*- shell-script -*-
# CMake shy configuration script.
# http://www.cmake.org

CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

# Add rules to the makefile to generate assembly files.
PROJECT(ASSEMBLER C)

IF(NOT CMAKE_BUILD_TYPE)
SET(CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

project(Shy)
message(STATUS "Configuring Shy")

set(LD_LIBRARY_PATH "${LD_LIBRARY_PATH}" "${CMAKE_SOURCE_DIR}/external/muparser/lib " "${CMAKE_SOURCE_DIR}/external/iaca-lin32/lib")

set (CMAKE_C_COMPILER g++)
set (CMAKE_CXX_COMPILER g++)

set(CMAKE_VERBOSE_MAKEFILE on)

IF(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_C_FLAGS "-g -O3 -msse4.2 -ffast-math -fno-tree-vectorize -ftree-vectorizer-verbose=1 -O3 -fopenmp -std=c99")
set(CMAKE_CXX_FLAGS "-g -O3 -fopenmp")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

IF(${CMAKE_C_COMPILER} MATCHES "icpc")
#set(CMAKE_C_FLAGS "-g -O3 -xsse4.2 -fp-model fast=2 -no-prec-div -no-prec-sqrt -fimf-domain-exclusion=15 -vec-report")
#set(CMAKE_C_FLAGS "-O3 -xCORE-AVX2 -vec-report2")
set(CMAKE_C_FLAGS "-O3 -xAVX -fimf-domain-exclusion=15 -vec-report2 -openmp")
ENDIF( ${CMAKE_C_COMPILER} MATCHES "icpc" )

IF(${CMAKE_CXX_COMPILER} MATCHES "icpc")
#set(CMAKE_CXX_FLAGS "-g -xsse4.2 -O3")
#set(CMAKE_CXX_FLAGS "-O3 -xCORE-AVX2 -vec-report2")
set(CMAKE_CXX_FLAGS "-O3 -xAVX -openmp")
ENDIF( ${CMAKE_CXX_COMPILER} MATCHES "icpc" )


# INCLUDE(TestBigEndian)
# TEST_BIG_ENDIAN(WORDS_BIGENDIAN)
# IF (WORDS_BIGENDIAN)
#  ADD_DEFINITIONS("-DBIG_ENDIAN")
# ELSE (WORDS_BIGENDIAN)
#  ADD_DEFINITIONS("-DLOW_ENDIAN")
# ENDIF (WORDS_BIGENDIAN)

ADD_DEFINITIONS(" -DRealType=double -I${CMAKE_SOURCE_DIR}/external/boost -I${CMAKE_SOURCE_DIR}/external/iaca-lin32/include")

add_executable(
    shy
    src/array_io.cpp
    src/cell.cpp
    src/clock.cpp
    src/event.cpp
    src/geom.cpp
    src/grid.cpp
    src/init.cpp
    src/kernel_advection_fv_uw_1_2d.c
    src/kernel_advection_fv_uw_1_3d.c
    src/kernel_advection_fv_uw_2_2d.c
    src/kernel_advection_fv_uw_2_3d.c
    src/kernel_advection_fv_uw_kappa_2d.c
    src/kernel_advection_fv_uw_kappa_3d.c
    src/kernel_advection_fv_central_compact_explicit_2d.c
    src/kernel_compressible_euler_fv_uw_kappa_2d.c
    src/math_parser.cpp
    src/mesh.cpp
    src/mpi_boundary_condition.cpp
    src/numerical_core.c
    src/numerical_params.cpp 
    src/output.cpp
    src/physical_params.cpp
    src/euler_riemann_analytic_solver.cpp
    src/shy.cpp 
    src/simulation.cpp 
    src/structured_grid.cpp
    src/timetable.cpp 
    src/variable_database.cpp
    src/variable_store.cpp)

# Telling cmake where our custom cmake modules are.
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set_property(SOURCE shy.cpp PROPERTY COMPILE_FLAGS "-O2 -g3")

# Eigen (for high performance array class, coming with template
# metaprogramming goodies, automatic vectorization etc.). This is a
# header only library.
#
# see : http://eigen.tuxfamily.org/
#
# CMake Eigen detection module is cmake/Modules/FindEigen.cmake
find_package(Eigen REQUIRED)

if (EIGEN_FOUND)
  include_directories(${EIGEN_INCLUDE_DIRS})
endif (EIGEN_FOUND)

# muParser -- mathematical expression parser.
#
# see : http://muparser.sourceforge.net/
#
# CMake muParser detection module is cmake/Modules/FindMuParser.cmake
find_package(MuParser REQUIRED)

if (MUPARSER_FOUND)
  include_directories(${MUPARSER_INCLUDE_DIR})
#  target_link_libraries(shy ${MUPARSER_LIBRARIES})
endif (MUPARSER_FOUND)

set(MUPARSER_LIBRARIES "${CMAKE_SOURCE_DIR}/external/muparser/lib/libmuparser.a")
target_link_libraries(shy ${MUPARSER_LIBRARIES})

set(REALTIME_LIBRARY "rt")
target_link_libraries(shy ${REALTIME_LIBRARY})

add_subdirectory(src)

# Copy files from source directory to destination directory, substituting any
# variables.  Create destination directory if it does not exist.

macro(configure_files srcDir destDir)
    make_directory(${destDir})

    file(GLOB templateFiles RELATIVE ${srcDir} ${srcDir}/*)
    foreach(templateFile ${templateFiles})
        set(srcTemplatePath ${srcDir}/${templateFile})
        if(NOT IS_DIRECTORY ${srcTemplatePath})
	    GET_FILENAME_COMPONENT(extension ${srcTemplatePath} EXT)

	    #Mesh binary files should be left untouched.
	    if("${extension}" STREQUAL ".msh")
	    	configure_file(
                    ${srcTemplatePath}
                    ${destDir}/${templateFile}
	    	    COPYONLY)
	    else("${extension}" STREQUAL ".msh")
	    	configure_file(
                    ${srcTemplatePath}
                    ${destDir}/${templateFile}
	    	    @COPY)
	    endif("${extension}" STREQUAL ".msh")
        endif(NOT IS_DIRECTORY ${srcTemplatePath})

    endforeach(templateFile)
endmacro(configure_files)

# Copy example from source directory to build directory.
message(STATUS "Copying examples")

set(exampleRoot "${CMAKE_CURRENT_SOURCE_DIR}/examples")
file(GLOB exampleDirs RELATIVE ${exampleRoot} "${exampleRoot}/*")

foreach(exampleDir ${exampleDirs})
set(examplePath ${exampleRoot}/${exampleDir})
if(IS_DIRECTORY ${examplePath})
message(STATUS "  - copying example simulation : ${exampleDir}")
configure_files(${examplePath} "${CMAKE_CURRENT_BINARY_DIR}/examples/${exampleDir}")
endif(IS_DIRECTORY ${examplePath})
endforeach(exampleDir ${exampleDirs})

#include(UseDoxygen OPTIONAL)

# # Go in the muparser directory and compile it.
# ExternalProject_Add(MuParser
#     PREFIX "${CMAKE_CURRENT_BUILD_DIR}"
#     SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/muparser"
#     BUILD_IN_SOURCE 0
#     CONFIGURE_COMMAND
#     "${CMAKE_CURRENT_SOURCE_DIR}/external/muparser/configure"
#     INSTALL_COMMAND ""
# )
