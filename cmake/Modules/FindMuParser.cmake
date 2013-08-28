# -*- shell-script -*- (for emacs users)

# - Find muParser
# Find the muParser library
#
#  This module defines the following variables:
#     MUPARSER_FOUND        - True if MUPARSER_INCLUDE_DIR & MUPARSER_LIBRARY are found
#     MUPARSER_LIBRARIES    - Set when MUPARSER_LIBRARY is found
#     MUPARSER_INCLUDE_DIRS - Set when MUPARSER_INCLUDE_DIR is found
#
#     MUPARSER_INCLUDE_DIR  - where to find muParser.h etc.
#     MUPARSER_LIBRARY      - the muParser library
#
#
# This file based on FindALSA.cmake
#
#=============================================================================
# Copyright 2009 Kitware, Inc.
# Copyright 2009 Philip Lowman <philip@yhbt.com>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)
find_path(MUPARSER_INCLUDE_DIR muParser.h
  #PATHS 
  ${CMAKE_SOURCE_DIR}/external/muparser/include
)

find_library(MUPARSER_LIBRARIES NAMES libmuparser.so
  #PATH 
    ${CMAKE_SOURCE_DIR}/external/muparser/lib
)

if (MUPARSER_LIBRARIES)
  MESSAGE(STATUS "Found muParser shared lib : ${MUPARSER_LIBRARIES}")
endif (MUPARSER_LIBRARIES)

if (MUPARSER_INCLUDE_DIR)
  set(MUPARSER_FOUND TRUE)
  #INCLUDE_DIRECTORIES(${MUPARSER_INCLUDE_DIRS})
  MESSAGE(STATUS "Found muParser headers in ${MUPARSER_INCLUDE_DIR}")
endif(MUPARSER_INCLUDE_DIR)
