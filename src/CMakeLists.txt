
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

option(BUILD_MODERN_FORTRAN_EXAMPLES OFF)
add_subdirectory(cminpack)

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

target_link_libraries(photochem PUBLIC yaml)
target_link_libraries(photochem PUBLIC sundials_fcvode_mod)
target_link_libraries(photochem PUBLIC sundials_cvode)
target_link_libraries(photochem PUBLIC cminpack)
target_link_libraries(photochem PUBLIC binning)
target_include_directories(photochem PUBLIC
                          "${CMAKE_BINARY_DIR}/modules"
                          )

