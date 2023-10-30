cmake_minimum_required(VERSION "3.14")
cmake_policy(SET CMP0148 OLD)
project(Photochem LANGUAGES Fortran C VERSION "0.4.4")

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
message(STATUS "The project is built using scikit-build: ${SKBUILD}")

option(BUILD_EXECUTABLES "if ON, then will build the
Fortran executables" ON)

option(BUILD_PYTHON_PHOTOCHEM "if ON, then will build a python
version via Cython" OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

add_subdirectory(src)
add_subdirectory(tests)

if (BUILD_PYTHON_PHOTOCHEM)
  if (NOT SKBUILD)
    if (NOT DEFINED SKBUILD_CMAKE_MODULE_DIR)
      # Here, we try to find scikit-build cmake modules
      find_package(Python COMPONENTS Development)
      set(SKBUILD_CMAKE_MODULE_DIR "${Python_LIBRARY_DIRS}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/skbuild/resources/cmake")
    endif()
    if (EXISTS ${SKBUILD_CMAKE_MODULE_DIR})
      message(STATUS "Scikit-build CMake modules: ${SKBUILD_CMAKE_MODULE_DIR}")
    else()
      message(FATAL_ERROR "Failed to find scikit-build CMake modules in directory: ${SKBUILD_CMAKE_MODULE_DIR}")
    endif()
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SKBUILD_CMAKE_MODULE_DIR})
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

