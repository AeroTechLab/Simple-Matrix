cmake_minimum_required( VERSION 2.8 )
project( Simple-Matrix )

set( CMAKE_C_STANDARD 99 )
set( CMAKE_C_STANDARD_REQUIRED ON )

set( LIBRARY_DIR ${CMAKE_CURRENT_LIST_DIR} CACHE PATH "Relative or absolute path to directory where built shared libraries will be placed" )

find_package( BLAS REQUIRED )
find_package( LAPACK REQUIRED )

find_library( has_matrix Matrix PATHS ${LIBRARY_DIR} )
if( NOT has_matrix )
  add_library( Matrix SHARED ${CMAKE_CURRENT_LIST_DIR}/matrix.c )
  set_target_properties( Matrix PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_DIR} )
  target_include_directories( Matrix PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
  target_compile_definitions( Matrix PUBLIC -DDEBUG )
  target_link_libraries( Matrix -lm ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )
endif()
