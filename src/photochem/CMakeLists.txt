set(fppFiles
  photochem_eqns.f90
  photochem_chemistry.f90
  photochem_evoatmosphere_rhs.f90)

fypp_f90("" "${fppFiles}" outFiles)

set(PHOTOCHEM_MODULES
  photochem_const.f90
  photochem_radtran.f90
  photochem_enum.f90
  photochem_settings.f90
  photochem_data.f90
  photochem_vars.f90
  photochem_wrk.f90)

set(PHOTOCHEM_EVOATMOSPHERE
  photochem_evoatmosphere.f90
  photochem_evoatmosphere_boundary.f90
  photochem_evoatmosphere_init.f90
  photochem_evoatmosphere_integrate.f90
  photochem_evoatmosphere_output.f90
  photochem_evoatmosphere_profiles.f90
  photochem_evoatmosphere_grid.f90)

set(PHOTOCHEM_SOURCES
  ${PHOTOCHEM_MODULES}
  ${PHOTOCHEM_EVOATMOSPHERE}
  ${outFiles}
  ${CMAKE_CURRENT_SOURCE_DIR}/../photochem.f90)

add_library(photochem
  ${PHOTOCHEM_VERSION_SOURCE}
  ${PHOTOCHEM_SOURCES})
target_link_libraries(photochem
  fortran-yaml-c
  sundials_fcvode_mod
  sundials_cvode
  futils
  clima
  finterp
  minpack
  differentia)
target_include_directories(photochem PUBLIC
  "${CMAKE_BINARY_DIR}/modules")

set_property(TARGET photochem PROPERTY POSITION_INDEPENDENT_CODE ON)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(photochem PRIVATE -Wunused -Wimplicit-interface -fimplicit-none)
endif()

option(BUILD_WITH_OPENMP "Compile with muli-threading" OFF)
if (BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if (OpenMP_Fortran_FOUND)
    target_compile_options(photochem PRIVATE ${OpenMP_Fortran_FLAGS})
    target_link_options(photochem PUBLIC ${OpenMP_Fortran_FLAGS})
  endif()
endif()
message(STATUS "Building Photochem with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")
