set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/mv/include )

set( MV_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvvd.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvmd.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvvf.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvmf.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvvdio.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvblasd.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvblasf.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvblasi.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvvi.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvmi.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvvc.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvmc.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvvcio.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/mv/src/mvblasc.cc
     )

set( SPBLAS_SOURCES     
        ${CMAKE_CURRENT_SOURCE_DIR}/spblas/spmm.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/spblas/spsm.cc
     )
add_library(spblas STATIC ${SPBLAS_SOURCES})

set(SPBLAS_LIBRARY spblas)
#set(SPBLAS_LIBRARY "$ENV{MKLROOT}/lib/intel64/libmkl_intel_ilp64.so;$ENV{MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a;$ENV{MKLROOT}/lib/intel64/libmkl_core.so;$ENV{MKLROOT}/../../lib/intel64/libiomp5.so;$ENV{MKLROOT}/../../lib/intel64/libifcore.so")

set( SPARSEpp_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/compcol_double.cc
    	${CMAKE_CURRENT_SOURCE_DIR}/src/comprow_double.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/coord_double.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/iotext_double.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/diagpre_double.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/icpre_double.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/ilupre_double.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/ilutpre_double.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/qsort_double.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/qsort_int.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/iohb_double.cc
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/iohb.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/tools.f
	    ${CMAKE_CURRENT_SOURCE_DIR}/src/pre.cc
    )

add_definitions( -DCOMPLEX=std::complex<double> )

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")


###############################################################
add_library(Sparse++ STATIC ${SPARSEpp_SOURCES} ${MV_SOURCES})
target_link_libraries(Sparse++ ${SPBLAS_LIBRARY})

