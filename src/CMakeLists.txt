
add_subdirectory(dependencies)

# version
configure_file(photochem_version.f90.in ${CMAKE_CURRENT_BINARY_DIR}/photochem_version.f90)
add_library(photochem_version STATIC ${CMAKE_CURRENT_BINARY_DIR}/photochem_version.f90)
target_include_directories(photochem_version PUBLIC
  "${CMAKE_BINARY_DIR}/modules"
)
set_property(TARGET photochem_version PROPERTY POSITION_INDEPENDENT_CODE ON)

add_subdirectory(clima)
add_subdirectory(equilibrate)
add_subdirectory(photochem)
