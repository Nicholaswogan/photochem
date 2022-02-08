cmake_minimum_required(VERSION "3.14")
project(Photochem LANGUAGES Fortran C VERSION "0.1.11")

include(FortranCInterface)
FortranCInterface_VERIFY()

# includes
include(cmake/CPM.cmake)

# options
option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

option(TRAVIS "Should be ON if it is a TravisCI build" OFF)

option(BUILD_EXECUTABLES "if ON, then will build the
Fortran executables" ON)

option(BUILD_PYTHON_PHOTOCHEM "if ON, then will build a python
version via Cython" OFF)

# version
configure_file(src/photochem_version.f90.in photochem_version.f90)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wunused")
  endif()
  if (TRAVIS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage")
  endif()
endif()

option(BUILD_WITH_OPENMP "Compile with muli-threading" OFF)
if (BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if (OpenMP_Fortran_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  endif()
endif()
message (STATUS "Building with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")

add_subdirectory(src)

if (BUILD_EXECUTABLES)

  # the main executable
  add_executable(photorun ${CMAKE_SOURCE_DIR}/src/main.f90)
  target_link_libraries(photorun photochem)
  target_include_directories(photorun PRIVATE ${CMAKE_BINARY_DIR}/modules)

  #  tests
  add_executable(testrun ${CMAKE_SOURCE_DIR}/tests/test.f90)
  target_link_libraries(testrun photochem)
  target_include_directories(testrun PRIVATE ${CMAKE_BINARY_DIR}/modules)
  
  add_executable(memtest ${CMAKE_SOURCE_DIR}/tests/memtest.f90)
  target_link_libraries(memtest photochem)
  target_include_directories(memtest PRIVATE ${CMAKE_BINARY_DIR}/modules)
  
endif()

if (BUILD_PYTHON_PHOTOCHEM)
  if (NOT SKBUILD)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
  endif()
  
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Cython REQUIRED)
  
  add_subdirectory(photochem)
  
  if (SKBUILD)
    install(TARGETS _photochem  DESTINATION photochem)
  else()
    install(TARGETS _photochem  DESTINATION ${CMAKE_SOURCE_DIR}/photochem)
  endif()
                
endif()

