# - Try to find Eigen
# Once done, this will define
#
#  Eigen_FOUND - system has Eigen
#  Eigen_INCLUDE_DIR - the Eigen include directory

find_path(EIGEN_INCLUDE_DIRS Eigen/Core
  ${CMAKE_SOURCE_DIR}/external/eigen
  #PATHS
  #/opt/eigen
  #/usr/local/eigen
  #${INCLUDE_INSTALL_DIR}
  #PATH_SUFFIXES eigen	
  )

if (EIGEN_INCLUDE_DIRS)
  set(EIGEN_FOUND TRUE)
  MESSAGE(STATUS "Found Eigen headers in ${EIGEN_INCLUDE_DIRS}")
endif(EIGEN_INCLUDE_DIRS)
