add_library(equilibrate
  equilibrate_yaml.f90
  equilibrate_const.f90
  equilibrate_cea.f90
  equilibrate.f90
)
target_link_libraries(equilibrate fortran-yaml-c)
target_include_directories(equilibrate PUBLIC
  "${CMAKE_BINARY_DIR}/modules"
)
set_property(TARGET equilibrate PROPERTY POSITION_INDEPENDENT_CODE ON)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(equilibrate PRIVATE
    -Wunused
    -Wimplicit-interface
    -fimplicit-none
  )
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(equilibrate PRIVATE -fcheck=all,no-array-temps)
  endif()
endif()
