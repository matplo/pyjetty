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
if (Python_User)
	message(STATUS "Using Python paths specified by the user...")
	set(Python3_FOUND)
else()
	set (Python3_FIND_ABI "ANY" "ANY" "ANY")
	find_package(Python3 REQUIRED COMPONENTS Interpreter NumPy)
	message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")
endif()

message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")
message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_NumPy_INCLUDE_DIRS: ${Python3_NumPy_INCLUDE_DIRS}")

include_directories(${Python3_INCLUDE_DIRS})
include_directories(${Python3_NumPy_INCLUDE_DIRS})

if (SWIG_INTERFACE_FILE)
		message(STATUS "Using swig file ${SWIG_INTERFACE_FILE}")
	else()
		set(SWIG_INTERFACE_FILE ${MODULE_NAME}.i)
		message(STATUS "Using swig file ${SWIG_INTERFACE_FILE} - from MODULE_NAME := ${MODULE_NAME}")
endif()

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(${SWIG_INTERFACE_FILE} PROPERTIES CPLUSPLUS ON)
set_property(SOURCE ${SWIG_INTERFACE_FILE} PROPERTY SWIG_MODULE_NAME ${MODULE_NAME})

# Add swig module
swig_add_library(${MODULE_NAME} TYPE SHARED LANGUAGE python SOURCES ${SWIG_INTERFACE_FILE})
swig_link_libraries(${MODULE_NAME} ${NAME_LIB} ${Python3_LIBRARIES} ${SWIG_MODULE_LINK_LIBRARIES})

# Files to install with Python
set(PYTHON_INSTALL_FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.py
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/_${MODULE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX})
endif(BUILD_PYTHON)
