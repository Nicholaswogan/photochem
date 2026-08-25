
if(WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMS_WIN64")
endif()

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
add_cython_module(_photochem _photochem.pyx ${_photochem_deps})

target_link_libraries(_photochem PRIVATE photochem_wrapper)
target_include_directories(_photochem PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)
set_property(TARGET _photochem PROPERTY LINKER_LANGUAGE Fortran)
