
add_library(binning interp.f90 sorting.f90)
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/interp.f90
PROPERTIES OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/binning.f90)
