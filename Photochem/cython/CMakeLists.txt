
add_cython_target(_Photochem _Photochem.pyx)
add_library(_Photochem MODULE ${_Photochem})
set_source_files_properties(_Photochem.pyx
                            PROPERTIES OBJECT_DEPENDS Atmosphere.pyx PhotochemData.pyx PhotochemVars.pyx)
python_extension_module(_Photochem)

target_link_libraries(_Photochem photochem_wrapper)
target_include_directories(_Photochem PUBLIC ${NumPy_INCLUDE_DIRS})