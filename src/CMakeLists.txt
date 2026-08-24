
add_subdirectory(dependencies)

# version
configure_file(photochem_version.f90.in ${CMAKE_CURRENT_BINARY_DIR}/photochem_version.f90)
set(PHOTOCHEM_VERSION_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/photochem_version.f90")

add_subdirectory(clima)
add_subdirectory(equilibrate)
add_subdirectory(photochem)
