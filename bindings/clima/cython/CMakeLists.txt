if (WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMS_WIN64")
endif()

set(_clima_deps
  ${CMAKE_CURRENT_SOURCE_DIR}/futils.pyx
  ${CMAKE_CURRENT_SOURCE_DIR}/futils_pxd.pxd
  ${CMAKE_CURRENT_SOURCE_DIR}/RTChannel.pyx
  ${CMAKE_CURRENT_SOURCE_DIR}/RTChannel_pxd.pxd
  ${CMAKE_CURRENT_SOURCE_DIR}/ClimaRadtranWrk.pyx
  ${CMAKE_CURRENT_SOURCE_DIR}/ClimaRadtranWrk_pxd.pxd
  ${CMAKE_CURRENT_SOURCE_DIR}/Radtran.pyx
  ${CMAKE_CURRENT_SOURCE_DIR}/Radtran_pxd.pxd
  ${CMAKE_CURRENT_SOURCE_DIR}/AdiabatClimate.pyx
  ${CMAKE_CURRENT_SOURCE_DIR}/AdiabatClimate_pxd.pxd)
add_cython_module(_clima _clima.pyx ${_clima_deps})

target_link_libraries(_clima PRIVATE clima_c_api)
target_include_directories(_clima PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET _clima PROPERTY LINKER_LANGUAGE Fortran)
