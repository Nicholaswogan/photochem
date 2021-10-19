
add_subdirectory(dependencies)

set(PHOTOCHEM_MODULES
  photochem_const.f90
  photochem_radtran.f90 
  photochem_eqns.f90 
  photochem_types.f90)

set(PHOTOCHEM_INPUT
    input/photochem_input.f90
    input/photochem_input_after_read.f90
    input/photochem_input_read.f90)

set(PHOTOCHEM_ATMOSPHERE
    atmosphere/photochem_atmosphere.f90
    atmosphere/photochem_atmosphere_init.f90
    atmosphere/photochem_atmosphere_rhs.f90
    atmosphere/photochem_atmosphere_integrate.f90
    atmosphere/photochem_atmosphere_utils.f90)

set(PHOTOCHEM_SOURCES 
  ${PHOTOCHEM_MODULES}
  ${PHOTOCHEM_INPUT}
  ${PHOTOCHEM_ATMOSPHERE}
  photochem.f90)

# libphotochem.a
add_library(photochem 
            ${PHOTOCHEM_SOURCES})
target_link_libraries(photochem fortran-yaml-cpp sundials_fcvode_mod sundials_cvode cminpack binning)
target_include_directories(photochem PUBLIC
                          "${CMAKE_BINARY_DIR}/modules")
