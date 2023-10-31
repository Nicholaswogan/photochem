
CPMAddPackage(
  NAME futils
  VERSION 0.1.7
  GITHUB_REPOSITORY "Nicholaswogan/futils"
  GIT_TAG "v0.1.7"
  EXCLUDE_FROM_ALL ON
)

# CPMAddPackage(
#   NAME finterp
#   VERSION 1.3.0
#   GITHUB_REPOSITORY "jacobwilliams/finterp"
#   GIT_TAG "1.3.0"
#   DOWNLOAD_ONLY ON
# )

# add_library(finterp ${finterp_SOURCE_DIR}/src/linear_interpolation_module.F90)

CPMAddPackage(
  NAME fortran-yaml-c
  VERSION 0.2.4
  GITHUB_REPOSITORY "Nicholaswogan/fortran-yaml-c"
  GIT_TAG "v0.2.4"
  EXCLUDE_FROM_ALL ON
)

CPMAddPackage(
  NAME cvode 
  VERSION 5.7.0
  OPTIONS 
    "BUILD_SHARED_LIBS OFF" 
    "BUILD_FORTRAN_MODULE_INTERFACE ON"
  URL "https://github.com/LLNL/sundials/releases/download/v5.7.0/cvode-5.7.0.tar.gz"
  EXCLUDE_FROM_ALL ON
)

# CPMAddPackage(
#   NAME minpack
#   VERSION 2.0.0
#   GITHUB_REPOSITORY "fortran-lang/minpack"
#   GIT_TAG "v2.0.0-rc.1"
#   DOWNLOAD_ONLY ON
# )

# add_library(minpack
#   ${minpack_SOURCE_DIR}/src/minpack.f90
# )

CPMAddPackage(
  NAME clima 
  VERSION 0.4.2
  OPTIONS 
    "BUILD_EXECUTABLE OFF"
    "BUILD_WITH_OPENMP ${BUILD_WITH_OPENMP}"
    "SKBUILD ${SKBUILD}"
    "BUILD_PYTHON_CLIMA ${BUILD_PYTHON_PHOTOCHEM}"
    "PYTHON_CLIMA_DESTINATION photochem"
  GITHUB_REPOSITORY "Nicholaswogan/clima"
  GIT_TAG "v0.4.2"
  EXCLUDE_FROM_ALL OFF
)