cmake_minimum_required(VERSION "3.12")
project(photochem LANGUAGES Fortran C CXX)

include(FortranCInterface)
FortranCInterface_VERIFY()
FortranCInterface_VERIFY(CXX)

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

option(BUILD_STATIC_PHOTOCHEM "if ON, then will build the
static library libphotochem.a" ON)

option(BUILD_PYTHON_PHOTOCHEM "if ON, then will build a python
version with f2py" OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" 
   AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wunused")
endif()

option(BUILD_WITH_OPENMP "Compile with muli-threading" OFF)
if (BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
endif()
message (STATUS "Building with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")

add_subdirectory(src)

if (BUILD_STATIC_PHOTOCHEM)

  # the main executable
  add_executable(photorun ${CMAKE_SOURCE_DIR}/src/main.f90)
  target_link_libraries(photorun photochem)
  target_include_directories(photorun PRIVATE ${CMAKE_BINARY_DIR}/modules)

  #  tests
  add_executable(testrun ${CMAKE_SOURCE_DIR}/tests/test.f90)
  target_link_libraries(testrun photochem)
  target_include_directories(testrun PRIVATE ${CMAKE_BINARY_DIR}/modules)

  enable_testing()
  add_test(test_photochem testrun)
endif()

if (BUILD_PYTHON_PHOTOCHEM)
  if (NOT SKBUILD)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
  endif()
  
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Cython REQUIRED)
  
  add_subdirectory(Photochem)
  
  if (SKBUILD)
    install(TARGETS _Photochem  DESTINATION Photochem)
  else()
    install(TARGETS _Photochem  DESTINATION ${CMAKE_SOURCE_DIR}/Photochem)
  endif()

                    
endif()



