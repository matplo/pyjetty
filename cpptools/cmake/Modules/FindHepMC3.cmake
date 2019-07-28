# - Locate HepMC library
# Defines:
#
#  HEPMC3_FOUND
#  HEPMC3_INCLUDE_DIR
#  HEPMC3_INCLUDE_DIRS (not cached)
#  HEPMC3_LIBRARIES
#  adopted from https://gitlab.cern.ch/sft/lcgcmake/tree/master/cmake/modules

find_path(HEPMC3_INCLUDE_DIR HepMC/GenEvent.h 
          HINTS ${HEPMC3_ROOT_DIR}/include $ENV{HEPMC3_ROOT_DIR}/include $ENV{HEPMC3_DIR}/include)
find_library(HEPMC3_LIBRARY NAMES HepMC 
             HINTS ${HEPMC3_ROOT_DIR}/lib $ENV{HEPMC3_ROOT_DIR}/lib $ENV{HEPMC3_DIR}/lib)

set(HEPMC3_INCLUDE_DIRS ${HEPMC3_INCLUDE_DIR})
set(HEPMC3_LIBRARIES ${HEPMC3_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set HEPMC3_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HepMC3 DEFAULT_MSG HEPMC3_INCLUDE_DIR HEPMC3_LIBRARY)

mark_as_advanced(HEPMC3_FOUND HEPMC3_INCLUDE_DIR HEPMC3_LIBRARY)
