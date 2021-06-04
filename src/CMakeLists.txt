cmake_minimum_required(VERSION "3.19.6")


project(photochem LANGUAGES Fortran)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# fortran yaml
add_subdirectory(fortran-yaml)
# set_property(TARGET yaml PROPERTY Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# binning
add_subdirectory(binning)

# CVODE
set(BUILD_SHARED_LIBS OFF)
set(EXAMPLES_ENABLE_C OFF)
set(EXAMPLES_ENABLE_F77 OFF)
set(EXAMPLES_ENABLE_F90 OFF)
set(EXAMPLES_ENABLE_F2003 OFF)
set(BUILD_FORTRAN_MODULE_INTERFACE ON)
set(BUILD_ARKODE OFF)
set(BUILD_CVODES OFF)
set(BUILD_IDA OFF)
set(BUILD_KINSOL OFF)
add_subdirectory(cvode-5.7.0)


add_library(photochem 
            photochem_const.f90 
        		photochem_types.f90 
        		photochem_input.f90 
        		photochem_wrk.f90 
        		photochem_vars.f90 
        		photochem_data.f90 
            photochem_setup.f90
        		photochem_radtran.f90
            photochem.f90 
)

# set_property(TARGET photochem PROPERTY Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(photochem PUBLIC yaml)
target_include_directories(photochem PUBLIC
                          "${CMAKE_BINARY_DIR}/modules"
                          )

add_executable(photo.run main.f90)
target_link_libraries(photo.run PUBLIC photochem binning yaml)
target_include_directories(photo.run PRIVATE ${CMAKE_BINARY_DIR}/modules)