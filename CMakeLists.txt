set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultithreadedSolversApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/ )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/external_includes )
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SparseLib++/1.7/include )
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SparseLib++/1.7/mv/include )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SuperLU_MT_2.2 )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/HSL )

#add_definitions( -DCOMPLEX=std::complex<double> )
#add_subdirectory(custom_external_libraries/SparseLib++/1.7)

add_subdirectory(custom_external_libraries/SuperLU_MT_2.2)
add_subdirectory(custom_external_libraries/HSL)

set(MULTITHREADED_SOLVERS_APP_USE_ARPACK TRUE)
if(${MULTITHREADED_SOLVERS_APP_USE_ARPACK} MATCHES TRUE)
    add_definitions(-DMULTITHREADED_SOLVERS_APP_USE_ARPACK)
    include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/arpack++/include )
    add_subdirectory(custom_external_libraries/ARPACK)
    #REMARKS: if compiler complains about DefineParameters, just put this-> before it
endif()

if(${MULTITHREADED_SOLVERS_APPLICATION_ENABLE_FEAST} MATCHES ON) # enable this in configure script
    set(MULTITHREADED_SOLVERS_APP_USE_FEAST TRUE)
endif()

if(${MULTITHREADED_SOLVERS_APP_USE_FEAST} MATCHES TRUE)
    add_definitions(-DMULTITHREADED_SOLVERS_APP_USE_FEAST)
    include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/FEAST/2.1/include )
    add_subdirectory(custom_external_libraries/FEAST/2.1)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/custom_cmake_modules")
if(UMFPACK_DIR)
    message("UMFPACK_DIR: " "${UMFPACK_DIR}")
    find_package(UMFPACK REQUIRED)
    message("UMFPACK_INCLUDE_DIRS: " "${UMFPACK_INCLUDE_DIRS}")
    message("UMFPACK_LIBRARIES: " "${UMFPACK_LIBRARIES}")
    include_directories( ${UMFPACK_INCLUDE_DIRS} )
    add_definitions(-DMULTITHREADED_SOLVERS_APPLICATION_USE_UMFPACK)
    add_definitions(-DMULTITHREADED_SOLVERS_APPLICATION_USE_AMD)
endif(UMFPACK_DIR)

if(DEFINED MULTITHREADED_SOLVERS_APP_METIS_DIR)
    INCLUDE_DIRECTORIES(${MULTITHREADED_SOLVERS_APP_METIS_DIR}/include)
    FIND_LIBRARY(METIS_LIB_FOR_MULTITHREADED_SOLVERS_APP metis "${MULTITHREADED_SOLVERS_APP_METIS_DIR}/lib" NO_DEFAULT_PATH)
    add_definitions(-DMULTITHREADED_SOLVERS_APPLICATION_USE_METIS)
    message("METIS_LIB_FOR_MULTITHREADED_SOLVERS_APP: " "${METIS_LIB_FOR_MULTITHREADED_SOLVERS_APP}")
endif()

## generate variables with the sources
set( KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/multithreaded_solvers_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/external_includes/r8lib.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/external_includes/condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mesh_rcm.cpp
)

if(${MULTITHREADED_SOLVERS_APP_USE_FEAST} MATCHES TRUE)
set( KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SOURCES
    ${KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_eigen_solvers/feast_solver.cpp
)
endif()

if(${MULTITHREADED_SOLVERS_APPLICATION_ENABLE_ARPACK} MATCHES ON) # enable this in configure script
    set(MULTITHREADED_SOLVERS_APP_USE_ARPACK TRUE)
endif()

if(${MULTITHREADED_SOLVERS_APP_USE_ARPACK} MATCHES TRUE)
set( KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SOURCES
    ${KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_eigen_solvers/arpack_solver.cpp
)
endif()

set( KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SOURCES
    ${KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SOURCES}
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_eigen_solvers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tools.f
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/multithreaded_solvers_python_application.cpp
)

#find_library(UMFPACK_LIB NAMES umfpack PATHS /usr/local/lib)
#find_library(AMD_LIB NAMES amd PATHS /usr/local/lib)
#find_library(CHOLMOD_LIB NAMES cholmod PATHS /usr/local/lib)
#find_library(COLAMD_LIB NAMES colamd PATHS /usr/local/lib)
#find_library(BLAS_LIB NAMES blas PATHS /usr/local/lib)
##find_library(BOOST_REGEX_LIB NAMES boost_regex PATHS ${BOOST_ROOT}/lib)
##find_library(BOOST_TIMER_LIB NAMES boost_timer PATHS ${BOOST_ROOT}/lib)
##find_library(BOOST_SYSTEM_LIB NAMES boost_system PATHS ${BOOST_ROOT}/lib)
##find_library(BOOST_RANDOM_LIB NAMES boost_random PATHS ${BOOST_ROOT}/lib)

if(PARDISO_LIBRARY)
    message("PARDISO_LIBRARY: " "${PARDISO_LIBRARY}")
    add_definitions(-DMULTITHREADED_SOLVERS_APPLICATION_USE_PARDISO)
endif(PARDISO_LIBRARY)

if(${USE_INTEL_GREATER_THAN_13} MATCHES TRUE)
    message("***************************************************************************")
    message("****** kratos configured for using the mkl shipped with intel >= 2013 *****")
    message("***************************************************************************")
  find_package(MKLSOLVER13 REQUIRED)
elseif(${USE_INTEL_GREATER_THAN_12} MATCHES TRUE)
    message("***************************************************************************")
    message("****** kratos configured for using the mkl shipped with intel >= 2012 *****")
    message("****** set USE_INTEL_GREATER_THAN_13 to TRUE to use a more modern version *")
    message("***************************************************************************")
    find_package(MKLSOLVER12 REQUIRED)
elseif(${USE_INTEL_GREATER_THAN_15} MATCHES TRUE)
    message("***************************************************************************")
    message("****** kratos configured for using the mkl shipped with intel >= 2015 *****")
    message("****** set USE_INTEL_GREATER_THAN_13 to TRUE to use a more modern version *")
    message("***************************************************************************")
    find_package(MKLSOLVER15 REQUIRED)
elseif(${USE_INTEL_GREATER_THAN_17} MATCHES TRUE)
    message("***************************************************************************")
    message("****** kratos configured for using the mkl shipped with intel >= 2017 *****")
    message("****** set USE_INTEL_GREATER_THAN_13 to TRUE to use a more modern version *")
    message("***************************************************************************")
    find_package(MKLSOLVER17 REQUIRED)
elseif(${USE_INTEL_MIC} MATCHES TRUE)
    message("***************************************************************************")
    message("****** kratos configured for using the mkl for Intel Xeon Phi         *****")
    message("***************************************************************************")
    find_package(MKLSOLVERmic REQUIRED)
    else(${USE_INTEL_GREATER_THAN_13} MATCHES TRUE)
    message("***************************************************************************")
    message("****** kratos configured for using the mkl shipped with intel < 2012*******")
    message("****** set USE_INTEL_GREATER_THAN_13 to TRUE to use a more modern version *")
    message("***************************************************************************")
    find_package(MKLSOLVER REQUIRED)
endif()

if (${MKLSOLVER_FOUND} MATCHES "YES")
    set(MULTITHREADED_SOLVERS_APP_USE_MKL TRUE)
    add_definitions(-DMULTITHREADED_SOLVERS_APP_USE_MKL)
endif()

########ADD SPECTRA LIBRARY########
FIND_PACKAGE(Eigen)
IF(${SPECTRA_FOUND} MATCHES TRUE)
    MESSAGE("Eigen is found at " ${EIGEN_DIR})
ENDIF()

FIND_PACKAGE(Spectra)
IF(${SPECTRA_FOUND} MATCHES TRUE)
    MESSAGE("Spectra is found at " ${SPECTRA_DIR})
    ADD_DEFINITIONS(-DMULTITHREADED_SOLVERS_APP_USE_SPECTRA)
ELSE()
    MESSAGE("Spetra is not found. The eigenvalues solver is not available.")
ENDIF()
#####################

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosMultithreadedSolversApplication SHARED ${KRATOS_MULTITHREADED_SOLVERS_APPLICATION_SOURCES})
#target_link_libraries(KratosMultithreadedSolversApplication KratosCore ${SPARSE_LIB} ${MV_LIB} ${SPBLAS_LIB} ${UMFPACK_LIB} ${AMD_LIB} ${CHOLMOD_LIB} ${COLAMD_LIB} ${BLAS_LIB})
#target_link_libraries(KratosMultithreadedSolversApplication KratosCore Sparse++ superlu_mt ${UMFPACK_LIBRARIES} rt ${PARDISO_LIBRARY} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(KratosMultithreadedSolversApplication KratosCore)
target_link_libraries(KratosMultithreadedSolversApplication superlu_mt)
target_link_libraries(KratosMultithreadedSolversApplication ${UMFPACK_LIBRARIES})
target_link_libraries(KratosMultithreadedSolversApplication rt)
target_link_libraries(KratosMultithreadedSolversApplication ${METIS_LIB_FOR_MULTITHREADED_SOLVERS_APP})
target_link_libraries(KratosMultithreadedSolversApplication ${PARDISO_LIBRARY})

if(${MULTITHREADED_SOLVERS_APP_USE_ARPACK} MATCHES TRUE)
    target_link_libraries(KratosMultithreadedSolversApplication arpack)
endif()

if(${MULTITHREADED_SOLVERS_APP_USE_FEAST} MATCHES TRUE)
    # REMARKS: Feast uses various subroutines of MKL, so we have to link manually here
    target_link_libraries(KratosMultithreadedSolversApplication feast_sparse)
    target_link_libraries(KratosMultithreadedSolversApplication feast)
    target_link_libraries(KratosMultithreadedSolversApplication mkl_intel_lp64)
    target_link_libraries(KratosMultithreadedSolversApplication mkl_intel_thread)
    target_link_libraries(KratosMultithreadedSolversApplication mkl_core)
    target_link_libraries(KratosMultithreadedSolversApplication imf)
    target_link_libraries(KratosMultithreadedSolversApplication svml)
    target_link_libraries(KratosMultithreadedSolversApplication irc)
    target_link_libraries(KratosMultithreadedSolversApplication iomp5)
    target_link_libraries(KratosMultithreadedSolversApplication pthread)
    target_link_libraries(KratosMultithreadedSolversApplication m)
endif()

target_link_libraries(KratosMultithreadedSolversApplication ${BLAS_LIBRARIES})
target_link_libraries(KratosMultithreadedSolversApplication ${LAPACK_LIBRARIES})
target_link_libraries(KratosMultithreadedSolversApplication ${MKLSOLVER_LIBRARIES})
target_link_libraries(KratosMultithreadedSolversApplication hsl)
target_link_libraries(KratosMultithreadedSolversApplication gfortran)
set_target_properties(KratosMultithreadedSolversApplication PROPERTIES COMPILE_DEFINITIONS "KRATOS_CORE=IMPORT;MULTITHREADED_SOLVERS_APPLICATION=EXPORT;")
set_target_properties(KratosMultithreadedSolversApplication PROPERTIES PREFIX "")
install(TARGETS KratosMultithreadedSolversApplication DESTINATION libs )
install(FILES ${PARDISO_LIBRARY} DESTINATION libs )

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMultithreadedSolversApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MultithreadedSolversApplication.py" DESTINATION KratosMultiphysics )
