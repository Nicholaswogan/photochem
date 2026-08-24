set(CLIMA_TESTS test_radtran test_adiabat)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/clima_test_paths.f90.in"
  "${CMAKE_CURRENT_BINARY_DIR}/clima_test_paths.f90"
  @ONLY
)
add_library(clima_test_support STATIC
  "${CMAKE_CURRENT_BINARY_DIR}/clima_test_paths.f90"
)
set(CLIMA_TEST_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/modules")
set_target_properties(clima_test_support PROPERTIES
  Fortran_MODULE_DIRECTORY "${CLIMA_TEST_MODULE_DIR}"
)
target_include_directories(clima_test_support PUBLIC
  "${CLIMA_TEST_MODULE_DIR}"
)

foreach(CLIMA_TEST ${CLIMA_TESTS})
  set(CLIMA_TEST_TARGET "clima_${CLIMA_TEST}")
  add_executable(${CLIMA_TEST_TARGET} "${CLIMA_TEST}.f90")
  set_target_properties(${CLIMA_TEST_TARGET} PROPERTIES
    OUTPUT_NAME "${CLIMA_TEST}"
  )
  target_link_libraries(${CLIMA_TEST_TARGET}
    clima
    clima_test_support
  )
  target_include_directories(${CLIMA_TEST_TARGET} PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
  )
endforeach()
