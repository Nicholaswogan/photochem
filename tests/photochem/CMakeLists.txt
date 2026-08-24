set(PHOTOCHEM_TESTS test_input test_api test_jacobian test_production_loss
                    test_memory test_evolution)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/photochem_test_paths.f90.in"
  "${CMAKE_CURRENT_BINARY_DIR}/photochem_test_paths.f90"
  @ONLY
)
add_library(photochem_test_support STATIC
  "${CMAKE_CURRENT_BINARY_DIR}/photochem_test_paths.f90"
)
set(PHOTOCHEM_TEST_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/modules")
set_target_properties(photochem_test_support PROPERTIES
  Fortran_MODULE_DIRECTORY "${PHOTOCHEM_TEST_MODULE_DIR}"
)
target_include_directories(photochem_test_support PUBLIC
  "${PHOTOCHEM_TEST_MODULE_DIR}"
)

foreach(PHOTOCHEM_TEST ${PHOTOCHEM_TESTS})
  set(PHOTOCHEM_TEST_TARGET "photochem_${PHOTOCHEM_TEST}")
  add_executable(${PHOTOCHEM_TEST_TARGET} "${PHOTOCHEM_TEST}.f90")
  set_target_properties(${PHOTOCHEM_TEST_TARGET} PROPERTIES
    OUTPUT_NAME "${PHOTOCHEM_TEST}"
  )
  target_link_libraries(${PHOTOCHEM_TEST_TARGET}
    photochem
    photochem_test_support
  )
  target_include_directories(${PHOTOCHEM_TEST_TARGET} PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
  )
endforeach()
