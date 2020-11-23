add_definitions(${ROOT_CXX_FLAGS})

set(MODULE_NAME fjtools)
set(MODULE ${PROJECT_NAME}_${MODULE_NAME})
set(NAME_LIB ${MODULE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SOURCES_LIB ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
message(STATUS "SOURCES: ${SOURCES_LIB}")
string(REPLACE ".cxx" ".hh" HEADERS_LIB "${SOURCES_LIB}")
string(REPLACE ".cxx" "_wrap.c" SWIG_HEADERS_LIB "${SOURCES_LIB}")

add_library(${NAME_LIB} SHARED ${SOURCES_LIB})

# Define preprocessor command for .hh include statement
# And set include/link directories appropriately
if(${RooUnfold_FOUND})
  add_compile_definitions(USE_ROOUNFOLD=ON)
  target_include_directories(${NAME_LIB} PUBLIC ${FASTJET_DIR}/include ${ROOT_INCLUDE_DIRS} ${ROOUNFOLD_INCLUDE_DIR})
  target_link_libraries(${NAME_LIB} PUBLIC ${FASTJET_LIBS} ${ROOT_LIBRARIES} ${ROOUNFOLD_LIBRARIES})
else(${RooUnfold_FOUND})
  target_include_directories(${NAME_LIB} PUBLIC ${FASTJET_DIR}/include ${ROOT_INCLUDE_DIRS})
  target_link_libraries(${NAME_LIB} PUBLIC ${FASTJET_LIBS} ${ROOT_LIBRARIES} )
endif(${RooUnfold_FOUND})

set(SWIG_TARGET_LINK_LIBRARIES "${FASTJET_LIBS} ${ROOT_LIBRARIES}")

swigify()

install_targets_custom()
