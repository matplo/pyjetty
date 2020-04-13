# This file is just an orchestration
cmake_minimum_required(VERSION 3.13.5)
# SWIG: use SWIG_MODULE_NAME property.
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif(POLICY CMP0086)
# target_sources: use absolute path for INTERFACE_SOURCES.
if(POLICY CMP0076)
  cmake_policy(SET CMP0076 NEW)
endif(POLICY CMP0076)
# SWIG: use standard target name.
if(POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif(POLICY CMP0078)

project(pyjetty VERSION 1.0 LANGUAGES CXX)
message(STATUS "version: ${PROJECT_VERSION}")

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)
set(PYJETTY_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
list(APPEND CMAKE_MODULE_PATH "${PYJETTY_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${PYJETTY_DIR}/cmake/Modules")
message(STATUS "module path ${CMAKE_MODULE_PATH}")
find_package(HEPPY REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${HEPPY_CMAKE_DIR}/../cmake")
message(STATUS "module path ${CMAKE_MODULE_PATH}")
include(common_heppy)
include(common_heppy_finds)
include(swigify)
include(install_targets_custom)

# Layout build dir like install dir
if(UNIX)
  option(BUILD_SHARED_LIBS "Build shared libraries(.so or .dyld)." ON)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_MACOSX_RPATH ON)
  # note here we are installing with RPATH explicit
  # set(CMAKE_INSTALL_RPATH ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  if (APPLE)
    # something for SIP - always FULL RPATH
    # https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
    # use, i.e. don't skip the full RPATH for the build tree
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)
    # when building, don't use the install RPATH already
    # (but later on when installing)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    # the RPATH to be used when installing, but only if it's not a system directory
    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    if("${isSystemDir}" STREQUAL "-1")
       set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    endif("${isSystemDir}" STREQUAL "-1")
    # something for brew
    include_directories("/usr/local/include")
    link_directories("/usr/local/lib")
  endif(APPLE)
  # for multi-config builds (e.g. msvc,xcode)
  foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/lib)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/lib)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/bin)
  endforeach()
else()
  # Windows only support static build.
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  # for multi-config builds (e.g. msvc,xcode)
  foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/bin)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/bin)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/bin)
  endforeach()
  # Currently Only support static build for windows
  #include(swig)
endif()

message( STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# Force default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug, Release (default), RelWithDebInfo and MinSizeRel."
    FORCE)
endif(NOT CMAKE_BUILD_TYPE)

add_subdirectory(src)

# for debugging only
# include(list_vars)
# foreach(pack PROJECT_DIR HEPPY_DIR Python_EXECUTABLE PYTHIA8_DIR ROOT_HEPPY_PREFIX HEPMC3_DIR HEPMC_DIR LHAPDF6_DIR FASTJET_DIR CMAKE_MODULE_PATH )
#   dump_cmake_variables(${pack})
# endforeach()

if (MAKE_MODULE)
  include(make_module)
  set(heppy_module_file ${PYJETTY_DIR}/modules/pyjetty/1.0)
  make_module()

  module_append_command("prepend-path" "LD_LIBRARY_PATH" "${PYJETTY_DIR}" "${PYJETTY_DIR}/cpptools/lib" TRUE)
  module_append_command("prepend-path" "LD_LIBRARY_PATH" "${PYJETTY_DIR}" "${PYJETTY_DIR}/cpptools/lib64" TRUE)

  module_append_command("prepend-path" "DYLD_LIBRARY_PATH" "${PYJETTY_DIR}" "${PYJETTY_DIR}/cpptools/lib" TRUE)
  module_append_command("prepend-path" "DYLD_LIBRARY_PATH" "${PYJETTY_DIR}" "${PYJETTY_DIR}/cpptools/lib64" TRUE)

  module_append_command("prepend-path" "PYTHONPATH" "${PYJETTY_DIR}" "${PYJETTY_DIR}/cpptools/lib" TRUE)
  module_append_command("prepend-path" "PYTHONPATH" "${PYJETTY_DIR}" "${PYJETTY_DIR}/cpptools/lib64" TRUE)
  module_append_command("prepend-path" "PYTHONPATH" "${PYJETTY_DIR}" "${PYJETTY_DIR}")

  module_append_command_no_prop("setenv" "PYJETTY_DIR" "${PYJETTY_DIR}" TRUE)
  module_append_command_no_prop("setenv" "PYJETTY_PYTHON_EXECUTABLE" "${Python_EXECUTABLE}" TRUE)
  module_append_command_no_prop("set-alias" "pyjettypython" "\"${Python_EXECUTABLE}\"" FALSE)
  module_append_command_no_prop("set-alias" "pyjetty_cd" "\"cd ${PYJETTY_DIR}\"" FALSE)
endif(MAKE_MODULE)
