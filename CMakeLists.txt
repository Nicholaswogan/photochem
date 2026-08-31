cmake_minimum_required(VERSION "3.15")
# Allow older vendored dependencies to configure with CMake 4.
set(CMAKE_POLICY_VERSION_MINIMUM "3.5")
cmake_policy(SET CMP0148 OLD)

project(Photochem LANGUAGES Fortran C VERSION "0.9.0")
set(PHOTOCHEM_CLIMA_DATA_VERSION "0.3.2")

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
if (SKBUILD)
  message(STATUS "The project is built using scikit-build-core ${SKBUILD_CORE_VERSION}")
else()
  message(STATUS "The project is built directly with CMake")
endif()

option(BUILD_EXECUTABLES "if ON, then will build the
Fortran executables" ON)

option(BUILD_PYTHON_PHOTOCHEM "if ON, then will build a python
version via Cython" OFF)

option(BUILD_WITH_OPENMP "Compile with multi-threading" OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

add_subdirectory(src)
add_subdirectory(tests)

if (BUILD_PYTHON_PHOTOCHEM)
  find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
  include(cmake/cython.cmake)
  add_subdirectory(bindings)
endif()
