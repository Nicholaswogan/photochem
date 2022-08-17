cmake_minimum_required(VERSION "3.14")
project(Photochem LANGUAGES Fortran C VERSION "0.3.2")

include(FortranCInterface)
FortranCInterface_VERIFY()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" AND 
    "${CMAKE_SYSTEM_PROCESSOR}" MATCHES "arm64" AND
    "${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
  if (${CMAKE_Fortran_COMPILER_VERSION} LESS "11.3.0")
    message (FATAL_ERROR "Photochem will only work with gfortran >= 11.3.0 for MacOS M1")
  endif()
endif()

# includes
include(cmake/CPM.cmake)

# options
option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

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

add_subdirectory(src)

if (BUILD_EXECUTABLES)

  #  tests
  add_executable(testrun ${CMAKE_SOURCE_DIR}/tests/testrun.f90)
  target_link_libraries(testrun photochem)
  target_include_directories(testrun PRIVATE ${CMAKE_BINARY_DIR}/modules)

  add_executable(testevo ${CMAKE_SOURCE_DIR}/tests/testevo.f90)
  target_link_libraries(testevo photochem)
  target_include_directories(testevo PRIVATE ${CMAKE_BINARY_DIR}/modules)
  
  add_executable(memtest ${CMAKE_SOURCE_DIR}/tests/memtest.f90)
  target_link_libraries(memtest photochem)
  target_include_directories(memtest PRIVATE ${CMAKE_BINARY_DIR}/modules)
  
endif()

if (BUILD_PYTHON_PHOTOCHEM)
  # set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  if (NOT SKBUILD)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
  endif()
  
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Cython REQUIRED)
  
  add_subdirectory(photochem)
  
  if (SKBUILD)
    install(TARGETS _photochem DESTINATION photochem)
  else()
    install(TARGETS _photochem DESTINATION ${CMAKE_SOURCE_DIR}/photochem)
  endif()
                
endif()

