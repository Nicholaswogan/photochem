add_library(clima_c_api clima_c_api.f90)
set(clima_c_api_deps
  ${CMAKE_CURRENT_SOURCE_DIR}/futils.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/RTChannel.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/ClimaRadtranWrk.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/Radtran.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/AdiabatClimate.f90)
set_source_files_properties(clima_c_api.f90 PROPERTIES
  OBJECT_DEPENDS "${clima_c_api_deps}")
target_include_directories(clima_c_api PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(clima_c_api clima)
