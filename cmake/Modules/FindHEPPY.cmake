# - Locate Heppy Setup
# Defines:
# HEPPY_FOUND
# HEPPY_DIR

message(STATUS "HEPPY_DIR $ENV{HEPPY_DIR}")
find_path(HEPPY_CMAKE_DIR common_heppy.cmake HINTS ${HEPPY_DIR}/cmake $ENV{HEPPY_DIR}/cmake)
if (HEPPY_CMAKE_DIR)
	set(HEPPY_FOUND TRUE)
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HEPPY DEFAULT_MSG HEPPY_FOUND HEPPY_CMAKE_DIR)
