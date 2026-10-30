# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.

find_path(SuperLUMT_INCLUDE_DIRS
  NAMES
  slu_mt_Cnames.h
  PATHS
  ${SuperLUMT_DIR}/include
  ${SuperLUMT_DIR}/include/superlu_mt
  NO_DEFAULT_PATH
)

message("SuperLUMT_DIR: " "${SuperLUMT_DIR}")
message("SuperLUMT_INCLUDE_DIRS: " "${SuperLUMT_INCLUDE_DIRS}")

find_library(SuperLUMT_LIBRARIES
  NAMES
  superlu_mt
  superlu_mt_OPENMP
  PATHS
  ${SuperLUMT_DIR}/lib
  NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SuperLUMT DEFAULT_MSG SuperLUMT_INCLUDE_DIRS SuperLUMT_LIBRARIES)

mark_as_advanced(SuperLUMT_INCLUDE_DIRS SuperLUMT_LIBRARIES)
