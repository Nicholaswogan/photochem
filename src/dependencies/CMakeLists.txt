
# fortran yaml
add_subdirectory(fortran-yaml-c EXCLUDE_FROM_ALL)

# binning
add_subdirectory(futils EXCLUDE_FROM_ALL)

# CVODE
set(BUILD_SHARED_LIBS OFF)
set(EXAMPLES_ENABLE_C OFF)
set(EXAMPLES_ENABLE_F77 OFF)
set(EXAMPLES_ENABLE_F90 OFF)
set(EXAMPLES_ENABLE_F2003 OFF)
set(BUILD_FORTRAN_MODULE_INTERFACE ON)
set(BUILD_ARKODE OFF)
set(BUILD_CVODES OFF)
set(BUILD_IDA OFF)
set(BUILD_KINSOL OFF)
add_subdirectory(cvode-5.7.0 EXCLUDE_FROM_ALL)

add_library(minpack 
minpack/src/minpack.f90 
minpack/src/minpack_kinds.f90)

