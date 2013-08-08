# -*- shell-script -*- (for emacs users)


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
