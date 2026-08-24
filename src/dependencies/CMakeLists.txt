
CPMAddPackage(
  NAME futils
  VERSION 0.1.15
  GITHUB_REPOSITORY "Nicholaswogan/futils"
  GIT_TAG "v0.1.15"
  EXCLUDE_FROM_ALL ON
)

CPMAddPackage(
  NAME fortran-yaml-c
  VERSION 0.2.6
  GITHUB_REPOSITORY "Nicholaswogan/fortran-yaml-c"
  GIT_TAG "v0.2.6"
  EXCLUDE_FROM_ALL ON
)

CPMAddPackage(
  NAME cvode 
  VERSION 5.7.0
  OPTIONS 
    "BUILD_SHARED_LIBS OFF" 
    "BUILD_FORTRAN_MODULE_INTERFACE ON"
    # Its legacy F2003 probe cannot configure with CMake 4. Supported
    # Photochem compilers already provide Fortran 2003.
    "F2003_FOUND TRUE"
  URL "https://github.com/LLNL/sundials/releases/download/v5.7.0/cvode-5.7.0.tar.gz"
  EXCLUDE_FROM_ALL ON
)

CPMAddPackage(
  NAME differentia
  VERSION 0.1.4
  OPTIONS 
    "BUILD_EXECUTABLE OFF"
  GITHUB_REPOSITORY "Nicholaswogan/Differentia"
  GIT_TAG "v0.1.4"
  EXCLUDE_FROM_ALL ON
)

# This fork keeps conda library directories available to HDF5 discovery.
CPMAddPackage(
  NAME h5fortran
  VERSION 4.5.0_c55f15e94300be746b22b4b840a6c4051d9fb323
  OPTIONS
    "CMAKE_POSITION_INDEPENDENT_CODE ON"
  GITHUB_REPOSITORY "Nicholaswogan/h5fortran"
  GIT_TAG "c55f15e94300be746b22b4b840a6c4051d9fb323"
  EXCLUDE_FROM_ALL ON
)

fypp_f90("" "linear_interpolation_module.F90" finterp_sources)
add_library(finterp ${finterp_sources})
target_link_libraries(finterp differentia)

CPMAddPackage(
  NAME dop853
  VERSION 1.2.0
  GITHUB_REPOSITORY "jacobwilliams/dop853"
  GIT_TAG "1.2.0"
  DOWNLOAD_ONLY ON
)
add_library(dop853
  ${dop853_SOURCE_DIR}/src/dop853_constants.f90
  ${dop853_SOURCE_DIR}/src/dop853_module.f90)

CPMAddPackage(
  NAME minpack
  VERSION 2.0.0
  GITHUB_REPOSITORY "fortran-lang/minpack"
  GIT_TAG "v2.0.0-rc.1"
  DOWNLOAD_ONLY ON
)
add_library(minpack
  ${minpack_SOURCE_DIR}/src/minpack.f90)
