cmake_minimum_required(VERSION "3.14")
cmake_policy(SET CMP0148 OLD)

project(Photochem LANGUAGES Fortran C VERSION "0.8.4")
set(PHOTOCHEM_CLIMA_DATA_VERSION "0.3.1")

include(FortranCInterface)
FortranCInterface_VERIFY()

# Ensure the compiler will work, and add some necessary global compiler and linker flags
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  if (${CMAKE_Fortran_COMPILER_VERSION} LESS "14.0.0")
    message(FATAL_ERROR "Photochem will only work with gfortran >= 14.0.0.")
  endif()
  # Avoid executable-stack trampolines from internal procedure callbacks.
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ftrampoline-impl=heap>)
endif()
if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
  # Linux/ELF link-time hardening: mark outputs as not requiring an executable stack.
  add_link_options(
    $<$<AND:$<LINK_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:-Wl,-z,noexecstack>
    $<$<AND:$<LINK_LANGUAGE:C>,$<C_COMPILER_ID:GNU>>:-Wl,-z,noexecstack>
  )
endif()

# includes
include(cmake/CPM.cmake)
include(cmake/fypp.cmake)

# Need fypp
find_program(FYPP fypp REQUIRED)

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
