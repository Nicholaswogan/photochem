
add_cython_target(_photochem _photochem.pyx)

if(WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMS_WIN64")
endif()

add_library(_photochem MODULE ${_photochem})
set(_photochem_deps
  EvoAtmosphere.pyx
  PhotochemData.pyx 
  PhotochemVars.pyx
  ProductionLoss.pyx 
  AtomConservation.pyx
)
foreach(PDEP ${_photochem_deps})
  set_source_files_properties(_photochem.pyx PROPERTIES OBJECT_DEPENDS ${PDEP})
endforeach()

python_extension_module(_photochem)

target_link_libraries(_photochem photochem_wrapper)
target_include_directories(_photochem PUBLIC ${NumPy_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET _photochem PROPERTY LINKER_LANGUAGE Fortran)