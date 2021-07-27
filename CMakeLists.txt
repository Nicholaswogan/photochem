cmake_minimum_required(VERSION "3.12")
project(photochem Fortran)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
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
file(COPY "${CMAKE_SOURCE_DIR}/data" DESTINATION "${CMAKE_BINARY_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/templates" DESTINATION "${CMAKE_BINARY_DIR}")

add_executable(photo.run ${CMAKE_SOURCE_DIR}/src/main.f90)
target_link_libraries(photo.run PUBLIC photochem)
target_include_directories(photo.run PRIVATE ${CMAKE_BINARY_DIR}/modules)

#  tests
add_executable(test.run ${CMAKE_SOURCE_DIR}/tests/test.f90)
target_link_libraries(test.run PUBLIC photochem)
target_include_directories(test.run PRIVATE ${CMAKE_BINARY_DIR}/modules)

enable_testing()
add_test(test_photochem test.run)