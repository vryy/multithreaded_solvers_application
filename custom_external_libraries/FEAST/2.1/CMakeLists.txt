set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )

set( FEAST_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/f90_functions_wrapper.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/src/feast_nompi.f90
#        ${CMAKE_CURRENT_SOURCE_DIR}/src/feast_sparse.f90
    )

set( FEAST_SPARSE_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/feast_sparse.f90
    )

###############################################################
add_library(feast STATIC ${FEAST_SRC})
add_library(feast_sparse STATIC ${FEAST_SPARSE_SRC})
#target_link_libraries(feast ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
# install(TARGETS feast DESTINATION libs )
