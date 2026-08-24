set(clima_fpp_files
  clima_saturationdata.f90)
fypp_f90("" "${clima_fpp_files}" clima_fpp_sources)

set(CLIMA_BASE_SOURCES
  clima_const.f90
  clima_useful.f90
  clima_eqns.f90
  clima_eqns_water.f90
  clima_types.f90
  clima_types_create.f90
  clima_ptc.f90)

set(CLIMA_RADTRAN_SOURCES
  clima_radtran_types.f90
  clima_radtran_types_create.f90
  clima_radtran_twostream.f90
  clima_radtran_radiate.f90
  clima_radtran.f90)

set(CLIMA_ADIABAT_SOURCES
  clima_adiabat_general.f90
  clima_adiabat_rc.f90
  clima_adiabat_solve.f90
  clima_adiabat_dry.f90
  clima_adiabat.f90
  clima_adiabat_altitude.f90)

find_package(LAPACK REQUIRED)

add_library(clima
  ${CLIMA_BASE_SOURCES}
  ${CLIMA_RADTRAN_SOURCES}
  ${CLIMA_ADIABAT_SOURCES}
  ${clima_fpp_sources}
  clima.f90)
target_link_libraries(clima
  fortran-yaml-c
  h5fortran
  futils
  finterp
  dop853
  minpack
  differentia
  ${LAPACK_LIBRARIES})
target_include_directories(clima PUBLIC
  "${CMAKE_BINARY_DIR}/modules")
set_property(TARGET clima PROPERTY POSITION_INDEPENDENT_CODE ON)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(clima PRIVATE
    -Wunused
    -Wimplicit-interface
    -fimplicit-none)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(clima PRIVATE -fcheck=all,no-array-temps)
  endif()
endif()

if (BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if (OpenMP_Fortran_FOUND)
    target_compile_options(clima PRIVATE ${OpenMP_Fortran_FLAGS})
    target_link_options(clima PUBLIC ${OpenMP_Fortran_FLAGS})
  endif()
endif()
message(STATUS "Building Clima with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")
