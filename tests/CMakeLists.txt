# First download data, if needed
if (NOT EXISTS "${CMAKE_SOURCE_DIR}/data")
  file(DOWNLOAD 
    "https://github.com/Nicholaswogan/photochem_clima_data/archive/v${PHOTOCHEM_CLIMA_DATA_VERSION}.tar.gz"
    "photochem_clima_data_v${PHOTOCHEM_CLIMA_DATA_VERSION}.tar.gz"
  )
  file(ARCHIVE_EXTRACT 
    INPUT "photochem_clima_data_v${PHOTOCHEM_CLIMA_DATA_VERSION}.tar.gz"
    DESTINATION ${CMAKE_BINARY_DIR}
  )
  file(RENAME 
    "${CMAKE_BINARY_DIR}/photochem_clima_data-${PHOTOCHEM_CLIMA_DATA_VERSION}/photochem_clima_data/data"
    "${CMAKE_SOURCE_DIR}/data"
  )
endif()

set(PHOTOCHEM_TESTS test_input test_api test_jacobian test_production_loss
                    test_memory test_evolution)

if (BUILD_EXECUTABLES)
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
    add_executable(${PHOTOCHEM_TEST} "${PHOTOCHEM_TEST}.f90")
    target_link_libraries(${PHOTOCHEM_TEST} photochem photochem_test_support)
    target_include_directories(${PHOTOCHEM_TEST} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
  endforeach()
endif()
