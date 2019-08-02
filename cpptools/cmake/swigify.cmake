function(SWIGIFY)
if (BUILD_PYTHON)
################
# swigify...
################
# Add subdirectories for each language if desired
# option(BUILD_PYTHON "Build Python SWIG module" ON)
# if(BUILD_PYTHON)
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
	message(STATUS "${Yellow}Swigify${ColourReset}")
	message(STATUS "${Yellow} ... for module ${MODULE_NAME} ${ColourReset}")
	# Include python
	# find_package(PythonLibs REQUIRED)
	# include_directories(${PYTHON_INCLUDE_PATH})
	if (Python_User)
			message(STATUS "Using Python environment specified by the user...")
			set(Python_NumPy_INCLUDE_DIRS $ENV{PYJETTY_PYTHON_NUMPY_INCLUDE_DIR})
			set(Python_INCLUDE_DIRS $ENV{PYJETTY_PYTHON_INCLUDE_DIR})
			set(Python_LIBRARIES $ENV{PYJETTY_PYTHON_CONFIG_LDFLAGS})
			set(Python_EXECUTABLE $ENV{PYJETTY_PYTHON_EXECUTABLE})
			set(Python_FOUND 1)
		else(Python_User)
			message(STATUS "No python environment specified by the user... - using findPython ...")
			set (Python_FIND_ABI "ANY" "ANY" "ANY")
			find_package(Python REQUIRED COMPONENTS Interpreter NumPy)
			message(STATUS "Python libraries: ${Python_LIBRARIES}")
	endif(Python_User)
	if (Python_FOUND)
		message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")
		message(STATUS "Python_INCLUDE_DIRS: ${Python_INCLUDE_DIRS}")
		message(STATUS "Python_NumPy_INCLUDE_DIRS: ${Python_NumPy_INCLUDE_DIRS}")
		message(STATUS "Python_LIBRARIES: ${Python_LIBRARIES}")

		include_directories(${Python_INCLUDE_DIRS})
		include_directories(${Python_NumPy_INCLUDE_DIRS})

		if (SWIG_INTERFACE_FILE)
				message(STATUS "Using swig file ${SWIG_INTERFACE_FILE}")
			else(SWIG_INTERFACE_FILE)
				set(SWIG_INTERFACE_FILE ${MODULE_NAME}.i)
				message(STATUS "Using swig file ${SWIG_INTERFACE_FILE} - from MODULE_NAME := ${MODULE_NAME}")
		endif(SWIG_INTERFACE_FILE)

		set(CMAKE_SWIG_FLAGS "")
		set_source_files_properties(${SWIG_INTERFACE_FILE} PROPERTIES CPLUSPLUS ON)
		set_property(SOURCE ${SWIG_INTERFACE_FILE} PROPERTY SWIG_MODULE_NAME ${MODULE_NAME})

		# Add swig module
		swig_add_library(${MODULE_NAME} TYPE SHARED LANGUAGE python SOURCES ${SWIG_INTERFACE_FILE})
		swig_link_libraries(${MODULE_NAME} ${NAME_LIB} ${Python_LIBRARIES} ${SWIG_MODULE_LINK_LIBRARIES})

		# Files to install with Python
		list(APPEND PYTHON_INSTALL_FILES
		        ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.py
		        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/_${MODULE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX})
		message(STATUS "${Green} ... built: ${PYTHON_INSTALL_FILES} ${ColourReset}")
		set(PYTHON_INSTALL_FILES ${PYTHON_INSTALL_FILES} PARENT_SCOPE)
	else(Python_FOUND)
		message(warning "Missing Python *** NO PYTHON INTERFACE WILL BE BUILD ***")
	endif(Python_FOUND)
endif(BUILD_PYTHON)
	message(STATUS "${Yellow} ... swigify done with ${MODULE_NAME} ${ColourReset}")
endfunction(SWIGIFY)