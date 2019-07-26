message( STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
execute_process ( COMMAND
                 ${CMAKE_CURRENT_SOURCE_DIR}/buildtools/get_fj_contrib.sh
                 ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/buildtools
                 OUTPUT_VARIABLE GET_FJ_CONTRIB
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
message( STATUS "GET_FJ_CONTRIB: ${GET_FJ_CONTRIB}")
