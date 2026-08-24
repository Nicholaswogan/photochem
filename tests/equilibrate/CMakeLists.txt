file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/equilibrate_test_paths.f90.in"
  "${CMAKE_CURRENT_BINARY_DIR}/equilibrate_test_paths.f90"
  @ONLY
)
add_library(equilibrate_test_support STATIC
  "${CMAKE_CURRENT_BINARY_DIR}/equilibrate_test_paths.f90"
)
set(EQUILIBRATE_TEST_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/modules")
set_target_properties(equilibrate_test_support PROPERTIES
  Fortran_MODULE_DIRECTORY "${EQUILIBRATE_TEST_MODULE_DIR}"
)
target_include_directories(equilibrate_test_support PUBLIC
  "${EQUILIBRATE_TEST_MODULE_DIR}"
)

add_executable(equilibrate_test_equilibrate test_equilibrate.f90)
set_target_properties(equilibrate_test_equilibrate PROPERTIES
  OUTPUT_NAME test_equilibrate
)
target_link_libraries(equilibrate_test_equilibrate
  equilibrate
  equilibrate_test_support
)
target_include_directories(equilibrate_test_equilibrate PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)
