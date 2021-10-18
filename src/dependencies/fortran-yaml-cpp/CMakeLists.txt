cmake_minimum_required(VERSION 3.4)

project(fortran-yaml-cpp VERSION 0.0.2 LANGUAGES CXX C Fortran)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/yaml-cpp/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

add_subdirectory(yaml-cpp EXCLUDE_FROM_ALL)

set(CMAKE_CXX_STANDARD 11)

add_library(fortran-yaml-cpp yaml_types.f90 yaml.f90 yaml.cpp)
target_link_libraries(fortran-yaml-cpp yaml-cpp)
target_include_directories(fortran-yaml-cpp PUBLIC yaml-cpp/include)

add_executable(test_yaml test_yaml.f90)
target_link_libraries(test_yaml fortran-yaml-cpp)

add_executable(example example.f90)
target_link_libraries(example fortran-yaml-cpp)

