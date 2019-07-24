# This file is just an orchestration
cmake_minimum_required(VERSION 3.14.2)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

# Apple: Don't modify install_name when touching RPATH.
if(POLICY CMP0068)
  cmake_policy(SET CMP0068 NEW)
endif()
# target_sources: use absolute path for INTERFACE_SOURCES.
if(POLICY CMP0076)
  cmake_policy(SET CMP0076 NEW)
endif()
# SWIG: use standard target name.
if(POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif()
# SWIG: use SWIG_MODULE_NAME property.
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif()

project(pyjetty VERSION 1.0 LANGUAGES CXX)
message(STATUS "version: ${PROJECT_VERSION}")

find_program ( FASTJETCONFIG fastjet-config )
message(STATUS "Using fjconfig at ${FASTJETCONFIG}")
if (EXISTS ${FASTJETCONFIG})
  execute_process ( COMMAND ${FASTJETCONFIG} --prefix WORKING_DIRECTORY /tmp OUTPUT_VARIABLE FASTJET_DIR OUTPUT_STRIP_TRAILING_WHITESPACE )
  execute_process ( COMMAND ${FASTJETCONFIG} --cxxflags WORKING_DIRECTORY /tmp OUTPUT_VARIABLE FASTJET_CXXFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE )
  execute_process ( COMMAND ${FASTJETCONFIG} --libs --plugins WORKING_DIRECTORY /tmp OUTPUT_VARIABLE FASTJET_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE )
  message( STATUS "FASTJET_DIR: ${FASTJET_DIR}")
  message( STATUS "FASTJET_CXXFLAGS: ${FASTJET_CXXFLAGS}")
  message( STATUS "FASTJET_LIBS: ${FASTJET_LIBS}")
  add_compile_options(${FASTJET_CXXFLAGS})
else()
  message(FATAL_ERROR "This really requires fastjet-config... - must be in your \$PATH")
endif()

find_program ( PYTHIA8CONFIG pythia8-config )
message(STATUS "Using pythia8-config at ${PYTHIA8CONFIG}")
if (EXISTS ${PYTHIA8CONFIG})
  execute_process ( COMMAND ${PYTHIA8CONFIG} --prefix WORKING_DIRECTORY /tmp OUTPUT_VARIABLE PYTHIA8_DIR OUTPUT_STRIP_TRAILING_WHITESPACE )
  execute_process ( COMMAND ${PYTHIA8CONFIG} --cxxflags WORKING_DIRECTORY /tmp OUTPUT_VARIABLE PYTHIA8_CXXFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE )
  execute_process ( COMMAND ${PYTHIA8CONFIG} --ldflags WORKING_DIRECTORY /tmp OUTPUT_VARIABLE PYTHIA8_LDFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE )
  message( STATUS "PYTHIA8_DIR: ${PYTHIA8_DIR}")
  message( STATUS "PYTHIA8_CXXFLAGS: ${PYTHIA8_CXXFLAGS}")
  message( STATUS "PYTHIA8_LDFLAGS: ${PYTHIA8_LDFLAGS}")
  add_compile_options(${PYTHIA8_CXXFLAGS})
else()
  message(FATAL_ERROR "Some of this requires pythia8-config... - must be in your \$PATH")
endif()

# Layout build dir like install dir
if(UNIX)
  option(BUILD_SHARED_LIBS "Build shared libraries(.so or .dyld)." ON)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  include_directories("/usr/local/include")
  link_directories("/usr/local/lib")
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

# now actual build instructions

set(MODULE_NAME cpptools)
set(MODULE ${PROJECT_NAME}_${MODULE_NAME})
set(NAME_LIB ${MODULE})

message(STATUS "SOURCES: ${PROJECT_SOURCE_DIR}")
include_directories(${PROJECT_SOURCE_DIR}/src)
file(GLOB SOURCES_LIB ${PROJECT_SOURCE_DIR}/src/*.cxx)
message(STATUS "SOURCES: ${SOURCES_LIB}")
string(REPLACE ".cxx" ".h" HEADERS_LIB "${SOURCES_LIB}")
string(REPLACE ".cxx" "_wrap.c" SWIG_HEADERS_LIB "${SOURCES_LIB}")

add_library(${NAME_LIB} SHARED ${SOURCES_LIB})
target_include_directories(${NAME_LIB} PUBLIC ${PYTHIA8_DIR}/include ${ROOT_INCLUDE_DIR} ${FASTJET_DIR}/include)
target_link_libraries(${NAME_LIB} ${PYTHIA8_LDFLAGS} ${ROOT_LIBRARIES} ${FASTJET_LIBS} ${Boost_LIBRARIES})

################
# swigify...
################
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# Add subdirectories for each language if desired
option(BUILD_PYTHON "Build Python SWIG module" ON)
if(BUILD_PYTHON)

  # Include python
  find_package(PythonLibs REQUIRED)
  include_directories(${PYTHON_INCLUDE_PATH})

  set(CMAKE_SWIG_FLAGS "")
  set_source_files_properties(src/${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
  set_property(SOURCE src/${PROJECT_NAME}.i PROPERTY SWIG_MODULE_NAME ${PROJECT_NAME})

  # Add swig module
  swig_add_library(${PROJECT_NAME} TYPE SHARED LANGUAGE python 
                  SOURCES src/${PROJECT_NAME}.i ${SOURCES_LIB})
  swig_link_libraries(${PROJECT_NAME} ${NAME_LIB} ${PYTHON_LIBRARIES})

  # Files to install with Python
  set(PYTHON_INSTALL_FILES
          ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.py
          ${CMAKE_CURRENT_BINARY_DIR}/lib/_${PROJECT_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX})
endif(BUILD_PYTHON)

#install_targets_custom()
install(TARGETS ${NAME_LIB} LIBRARY DESTINATION lib)
install(FILES ${PYTHON_INSTALL_FILES} DESTINATION lib)
install(CODE "execute_process(COMMAND ln -s ${CMAKE_CURRENT_BINARY_DIR}/lib/_${PROJECT_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX} ${CMAKE_CURRENT_BINARY_DIR}/lib/_${PROJECT_NAME}.so)")
