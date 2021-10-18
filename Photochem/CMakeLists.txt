add_library(photochem_wrapper wrapper_utils.f90 
                              Atmosphere_wrapper.f90 
                              PhotochemData_wrapper.f90
                              PhotochemVars_wrapper.f90)
target_link_libraries(photochem_wrapper photochem)

add_cython_target(_Photochem _Photochem.pyx)
add_library(_Photochem MODULE ${_Photochem})
set_source_files_properties(_Photochem.pyx
                            PROPERTIES OBJECT_DEPENDS Atmosphere.pyx PhotochemData.pyx PhotochemVars.pyx)
python_extension_module(_Photochem)

target_link_libraries(_Photochem photochem_wrapper)
target_include_directories(_Photochem PUBLIC ${NumPy_INCLUDE_DIRS}
${CMAKE_CURRENT_SOURCE_DIR})