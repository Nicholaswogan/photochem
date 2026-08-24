
add_cython_target(_photochem _photochem.pyx)

if(WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMS_WIN64")
endif()

add_library(_photochem MODULE ${_photochem})
set(_photochem_deps
  AtomConservation.pyx
  AtomConservation_pxd.pxd
  EvoAtmosphere.pyx
  EvoAtmosphere_pxd.pxd
  PhotochemData.pyx
  PhotochemData_pxd.pxd
  PhotochemVars.pyx
  PhotochemVars_pxd.pxd
  PhotochemWrk.pyx
  PhotochemWrk_pxd.pxd
  PhotochemWrk.h
  ProductionLoss.pyx
  ProductionLoss_pxd.pxd
)
set_source_files_properties(_photochem.pyx PROPERTIES
  OBJECT_DEPENDS "${_photochem_deps}"
)

python_extension_module(_photochem)

target_link_libraries(_photochem photochem_wrapper)
target_include_directories(_photochem PUBLIC
  ${NumPy_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
)
set_property(TARGET _photochem PROPERTY LINKER_LANGUAGE Fortran)
