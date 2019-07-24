# - Find CGAL
# Find the CGAL includes and client library
# This module defines
#  CGAL_INCLUDE_DIR, where to find CGAL.h
#  CGAL_LIBRARIES, the libraries needed to use CGAL.
#  CGAL_FOUND, If false, do not try to use CGAL.

if(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
   set(CGAL_FOUND TRUE)

else(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
  if ((IS_DIRECTORY ${CGAL_ROOT}) OR (IS_DIRECTORY $ENV{CGAL_DIR}))
    message(STATUS "Preferring the environment CGAL settings over system installation...")
    find_path(CGAL_INCLUDE_DIR CGAL/basic.h
              HINTS
              $ENV{CGAL_DIR}/include
              $ENV{CGAL_ROOT}/include)

    find_library(CGAL_LIBRARIES NAMES CGAL libCGAL
                 HINTS
                 $ENV{CGAL_DIR}/lib
                 $ENV{CGAL_ROOT}/lib
                 $ENV{CGAL_DIR}/lib64
                 $ENV{CGAL_ROOT}/lib64)
    else()
    find_path(CGAL_INCLUDE_DIR CGAL/basic.h
              /usr/local/include
              /usr/include)

    find_library(CGAL_LIBRARIES NAMES CGAL libCGAL
                 PATHS
                 /usr/lib
                 /usr/lib64
                 /usr/local/lib
                 /usr/local/lib64)

  endif((IS_DIRECTORY ${CGAL_ROOT}) OR (IS_DIRECTORY $ENV{CGAL_DIR}))

  mark_as_advanced(CGAL_INCLUDE_DIR CGAL_LIBRARIES)

endif(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)

if(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
  message(STATUS "CGAL_INCLUDE_DIR=${CGAL_INCLUDE_DIR}")
  message(STATUS "CGAL_LIBRARIES=${CGAL_LIBRARIES}")
  set(CGAL_FOUND TRUE)
else(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
  set(CGAL_FOUND FALSE)
  message(STATUS "CGAL not found.")
endif(CGAL_INCLUDE_DIR AND CGAL_LIBRARIES)
