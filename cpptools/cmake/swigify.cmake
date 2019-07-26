if (BUILD_PYTHON)
################
# swigify...
################
# Add subdirectories for each language if desired
# option(BUILD_PYTHON "Build Python SWIG module" ON)
# if(BUILD_PYTHON)
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# Include python
# find_package(PythonLibs REQUIRED)
# include_directories(${PYTHON_INCLUDE_PATH})
find_package(Python3 REQUIRED COMPONENTS Interpreter NumPy)
message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${Python3_NumPy_INCLUDE_DIRS})

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(${MODULE_NAME}.i PROPERTIES CPLUSPLUS ON)
set_property(SOURCE ${MODULE_NAME}.i PROPERTY SWIG_MODULE_NAME ${MODULE_NAME})

# Add swig module
swig_add_library(${MODULE_NAME} TYPE SHARED LANGUAGE python 
                SOURCES ${MODULE_NAME}.i ${SOURCES_LIB})
swig_link_libraries(${MODULE_NAME} ${NAME_LIB} ${Python3_LIBRARIES})

# Files to install with Python
set(PYTHON_INSTALL_FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.py
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/_${MODULE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX})
endif(BUILD_PYTHON)
