add_compile_options(${PYTHIA8_CXXFLAGS})

set(MODULE_NAME pythiahepmc3)
set(MODULE ${PROJECT_NAME}_${MODULE_NAME})
set(NAME_LIB ${MODULE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SOURCES_LIB ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
#if (APPLE)
#	# this is something to be understood - on linux it is NOT needed... additional link needed?
#	list(APPEND SOURCES_LIB ${HEPMC3_INCLUDE_DIRS}/../interfaces/pythia8/src/Pythia8ToHepMC3.cc)
#endif(APPLE)
list(APPEND SOURCES_LIB ${HEPMC3_INCLUDE_DIRS}/../interfaces/pythia8/src/Pythia8ToHepMC3.cc)
message(STATUS "SOURCES: ${SOURCES_LIB}")
string(REPLACE ".cxx" ".hh" HEADERS_LIB "${SOURCES_LIB}")
string(REPLACE ".cxx" "_wrap.c" SWIG_HEADERS_LIB "${SOURCES_LIB}")

add_library(${NAME_LIB} SHARED ${SOURCES_LIB})
target_include_directories(${NAME_LIB} PUBLIC ${HEPMC3_INCLUDE_DIRS} ${HEPMC3_INCLUDE_DIRS}/../interfaces/pythia8/include)
target_link_libraries(${NAME_LIB} PUBLIC ${PYTHIA8_LDFLAGS} ${HEPMC3_LIBRARIES})

set(SWIG_TARGET_LINK_LIBRARIES ${PYTHIA8_LDFLAGS} ${HEPMC3_LIBRARIES})

swigify()

install_targets_custom()
