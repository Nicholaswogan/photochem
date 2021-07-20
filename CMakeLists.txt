cmake_minimum_required(VERSION "3.19.1")
project(photochem LANGUAGES Fortran)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
add_subdirectory(src)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

file(COPY "${CMAKE_SOURCE_DIR}/data" DESTINATION "${CMAKE_BINARY_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/templates" DESTINATION "${CMAKE_BINARY_DIR}")

#  tests
add_executable(test.run ${CMAKE_SOURCE_DIR}/tests/test.f90)
target_link_libraries(test.run PUBLIC photochem binning yaml)
target_include_directories(test.run PRIVATE ${CMAKE_BINARY_DIR}/modules)

enable_testing()
add_test(test_photochem test.run)