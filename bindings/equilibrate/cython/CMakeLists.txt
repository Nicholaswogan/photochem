if (WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMS_WIN64")
endif()

add_cython_module(
  _equilibrate
  _equilibrate.pyx
  "${CMAKE_CURRENT_SOURCE_DIR}/ChemEquiAnalysis_pxd.pxd"
  "${CMAKE_CURRENT_SOURCE_DIR}/_equilibrate.h"
)

target_link_libraries(_equilibrate PRIVATE equilibrate_c_api)
target_include_directories(_equilibrate PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}"
)
set_property(TARGET _equilibrate PROPERTY LINKER_LANGUAGE Fortran)
