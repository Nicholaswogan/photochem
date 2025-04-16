
add_subdirectory(dependencies)

# version
configure_file(photochem_version.f90.in ${CMAKE_CURRENT_BINARY_DIR}/photochem_version.f90)

# Adding subdirectories just to create this folder in the build dir
add_subdirectory(evoatmosphere)

set(fppFiles 
  photochem_eqns.f90
  photochem_common.f90
  evoatmosphere/photochem_evoatmosphere_rhs.f90)

fypp_f90("" "${fppFiles}" outFiles)

set(PHOTOCHEM_MODULES
  photochem_const.f90
  photochem_radtran.f90 
  photochem_enum.f90
  photochem_types.f90
  photochem_types_create.f90)

set(PHOTOCHEM_INPUT
    input/photochem_input.f90
    input/photochem_input_after_read.f90
    input/photochem_input_read.f90)

set(PHOTOCHEM_EVOATMOSPHERE
    evoatmosphere/photochem_evoatmosphere.f90
    evoatmosphere/photochem_evoatmosphere_init.f90
    evoatmosphere/photochem_evoatmosphere_rhs_climate.f90
    evoatmosphere/photochem_evoatmosphere_integrate.f90
    evoatmosphere/photochem_evoatmosphere_utils.f90)

set(PHOTOCHEM_SOURCES 
  ${PHOTOCHEM_MODULES}
  ${PHOTOCHEM_INPUT}
  ${PHOTOCHEM_EVOATMOSPHERE}
  ${outFiles}
  photochem.f90)

# libphotochem.a
add_library(photochem 
            ${CMAKE_CURRENT_BINARY_DIR}/photochem_version.f90
            ${PHOTOCHEM_SOURCES})
target_link_libraries(photochem 
  fortran-yaml-c 
  sundials_fcvode_mod 
  sundials_cvode 
  futils
  clima
  finterp 
  minpack
  differentia
)
target_include_directories(photochem PUBLIC
                          "${CMAKE_BINARY_DIR}/modules")

# This is important for making Python wrapper
set_property(TARGET photochem PROPERTY POSITION_INDEPENDENT_CODE ON)

# Add special compiler flags
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
message (STATUS "Building Photochem with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")