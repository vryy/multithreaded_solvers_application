set(CMAKE_INCLUDE_CURRENT_DIR ON)

#include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/SRC )

set( UTIL_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/UTIL/icnteq.f
        ${CMAKE_CURRENT_SOURCE_DIR}/UTIL/icopy.f
        ${CMAKE_CURRENT_SOURCE_DIR}/UTIL/iset.f
        ${CMAKE_CURRENT_SOURCE_DIR}/UTIL/iswap.f
        ${CMAKE_CURRENT_SOURCE_DIR}/UTIL/ivout.f
        ${CMAKE_CURRENT_SOURCE_DIR}/UTIL/second.f
        ${CMAKE_CURRENT_SOURCE_DIR}/UTIL/dvout.f
        ${CMAKE_CURRENT_SOURCE_DIR}/UTIL/dmout.f
     )

set( SRC_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgetv0.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlaqrb.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dstqrb.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsortc.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsortr.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dstatn.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dstats.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dnaitr.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dnapps.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dnaup2.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dnaupd.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dnconv.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dneigh.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dngets.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsaitr.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsapps.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsaup2.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsaupd.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsconv.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dseigt.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsgets.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dneupd.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dseupd.f
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsesrt.f
     )

###############################################################
add_library(arpack STATIC ${UTIL_SRC} ${SRC_SRC})
target_link_libraries(arpack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
# install(TARGETS arpack DESTINATION libs )
