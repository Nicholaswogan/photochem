
add_library(photochem_wrapper wrapper_utils.f90 
                              Atmosphere_wrapper.f90 
                              PhotochemData_wrapper.f90
                              PhotochemVars_wrapper.f90
                              PhotochemWrk_wrapper.f90
                              ProductionLoss_wrapper.f90)
target_link_libraries(photochem_wrapper photochem)
