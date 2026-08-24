add_cython_target(_clima _clima.pyx)

if (WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMS_WIN64")
endif()

add_library(_clima MODULE ${_clima})
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
set_source_files_properties(_clima.pyx PROPERTIES
  OBJECT_DEPENDS "${_clima_deps}")

python_extension_module(_clima)

target_link_libraries(_clima clima_c_api)
target_include_directories(_clima PUBLIC
  ${NumPy_INCLUDE_DIRS}
  "${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET _clima PROPERTY LINKER_LANGUAGE Fortran)
