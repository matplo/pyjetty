if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
string(TOUPPER "${CMAKE_BUILD_TYPE}" cmake_build_type_toupper)

if(NOT cmake_build_type_tolower STREQUAL "debug" AND
   NOT cmake_build_type_tolower STREQUAL "release" AND
   NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
    message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" c11flag)
#message( STATUS "c11flag is ${c11flag}")
if (c11flag)
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11" )
else()
	CHECK_CXX_COMPILER_FLAG("-std=gnu++11" gnu11flag)
	#message( STATUS "gnu11flag is ${gnu11flag}")
	if(gnu11flag)
		SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=gnu++11" )
	else()
		message( FATAL_ERROR "Compiler does not support -std=c++11 nor -std=gnu++11" )
	endif()
endif()
message ( STATUS "CMAKE_CXX_FLAGS =${CMAKE_CXX_FLAGS}")
