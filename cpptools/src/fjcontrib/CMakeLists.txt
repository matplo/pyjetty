message( STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
execute_process ( COMMAND
                 ${CMAKE_CURRENT_SOURCE_DIR}/buildtools/get_fj_contrib.sh
                 ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/buildtools
                 OUTPUT_VARIABLE GET_FJ_CONTRIB
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
message( STATUS "GET_FJ_CONTRIB: ${GET_FJ_CONTRIB}")

add_compile_options(${FASTJET_CXXFLAGS})

set(MODULE_NAME fjcontrib)
set(MODULE ${PROJECT_NAME}_${MODULE_NAME})
set(NAME_LIB ${MODULE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB_RECURSE SOURCES_LIB ${CMAKE_CURRENT_SOURCE_DIR}/fjcontrib-1.041/*.cc)
message(STATUS "SOURCES: ${SOURCES_LIB}")
string(REPLACE ".cc" ".hh" HEADERS_LIB "${SOURCES_LIB}")
string(REPLACE ".cc" "_wrap.c" SWIG_HEADERS_LIB "${SOURCES_LIB}")
add_library(${NAME_LIB} SHARED ${SOURCES_LIB})
target_include_directories(${NAME_LIB} PUBLIC ${FASTJET_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/fjcontrib-1.041)
target_link_libraries(${NAME_LIB} PUBLIC ${FASTJET_LIBS})

set(SWIG_INTERFACE_FILE interface/${MODULE_NAME}.i)
set(SWIG_TARGET_LINK_LIBRARIES ${FASTJET_LIBS})
swigify()
install_targets_custom()
