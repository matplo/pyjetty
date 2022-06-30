# add_compile_options(${FASTJET_CXXFLAGS})
include(GNUInstallDirs)

set(MODULE_NAME rutilext)
set(MODULE ${PROJECT_NAME}_${MODULE_NAME})
set(NAME_LIB ${MODULE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${FASTJET_DIR}/include)
add_definitions(${ROOT_CXX_FLAGS})
file(GLOB SOURCES_LIB ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
message(STATUS "SOURCES: ${SOURCES_LIB}")
string(REPLACE ".cxx" ".hh" HEADERS_LIB "${SOURCES_LIB}")
# string(REPLACE ".cxx" "_wrap.c" SWIG_HEADERS_LIB "${SOURCES_LIB}")

# generate ROOT dictionaries
set(DICT ${NAME_LIB}_Dict )
set(LINKDEFS "${NAME_LIB}_LinkDef.h" )
ROOT_GENERATE_DICTIONARY(${DICT} ${HEADERS_LIB} LINKDEF ${LINKDEFS})
add_library(${NAME_LIB} SHARED ${SOURCES_LIB} ${DICT})

target_link_libraries(${NAME_LIB} PUBLIC ${FASTJET_LIBS} ${ROOT_LIBRARIES})

# set(SWIG_TARGET_LINK_LIBRARIES ${FASTJET_LIBS})
# swigify()
# install_targets_custom()

# message(STATUS "cmake install libdir=${CMAKE_INSTALL_LIBDIR}")
install(TARGETS ${NAME_LIB} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT library)
# install(FILES ${HEADERS_LIB} DESTINATION include)
# postfix MP - on 64 bit archs these files land in lib64 only while SO in lib - not needed using ${CMAKE_INSTALL_LIBDIR} instead of lib
# install(FILES ${PROJECT_BINARY_DIR}/lib/${libprefix}${DICT}_rdict.pcm DESTINATION lib)
# install(FILES ${PROJECT_BINARY_DIR}/lib/${libprefix}${DICT}.rootmap DESTINATION lib)
