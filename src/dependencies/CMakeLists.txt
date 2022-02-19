
add_subdirectory(futils)

CPMAddPackage(
  NAME cminpack-f
  VERSION 0.1.2
  GITHUB_REPOSITORY "Nicholaswogan/cminpack-f"
  GIT_TAG "v0.1.2"
  EXCLUDE_FROM_ALL ON
)

CPMAddPackage(
  NAME fortran-yaml-c
  VERSION 0.1.2
  GITHUB_REPOSITORY "Nicholaswogan/fortran-yaml-c"
  GIT_TAG "v0.1.2"
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
