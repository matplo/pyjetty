add_compile_options(${PYTHIA8_CXXFLAGS})

set(MODULE_NAME pythiaext)
set(MODULE ${PROJECT_NAME}_${MODULE_NAME})
set(NAME_LIB ${MODULE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SOURCES_LIB ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
message(STATUS "SOURCES: ${SOURCES_LIB}")
string(REPLACE ".cxx" ".hh" HEADERS_LIB "${SOURCES_LIB}")
string(REPLACE ".cxx" "_wrap.c" SWIG_HEADERS_LIB "${SOURCES_LIB}")

add_library(${NAME_LIB} SHARED ${SOURCES_LIB})
target_include_directories(${NAME_LIB} PUBLIC ${HEPMC_INCLUDE_DIRS})
target_link_libraries(${NAME_LIB} PUBLIC ${PYTHIA8_LDFLAGS} ${HEPMC_LIBRARIES})

set(SWIG_TARGET_LINK_LIBRARIES ${PYTHIA8_LDFLAGS} ${HEPMC_LIBRARIES})

swigify()

install_targets_custom()
