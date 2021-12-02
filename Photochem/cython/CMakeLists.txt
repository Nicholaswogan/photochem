
add_cython_target(_Photochem _Photochem.pyx)

if(WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMS_WIN64")
endif()

add_library(_Photochem MODULE ${_Photochem})
set_source_files_properties(_Photochem.pyx 
PROPERTIES OBJECT_DEPENDS Atmosphere.pyx PhotochemData.pyx PhotochemVars.pyx
                          ProductionLoss.pyx AtomConservation.pyx)
python_extension_module(_Photochem)

target_link_libraries(_Photochem photochem_wrapper)
target_include_directories(_Photochem PUBLIC ${NumPy_INCLUDE_DIRS})