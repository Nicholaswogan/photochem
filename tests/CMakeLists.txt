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

set(PHOTOCHEM_TESTS test_input test_api test_memory test_evolution test_block_thomas)

if (BUILD_EXECUTABLES)
  foreach(PHOTOCHEM_TEST ${PHOTOCHEM_TESTS})
    add_executable(${PHOTOCHEM_TEST} "${PHOTOCHEM_TEST}.f90")
    target_link_libraries(${PHOTOCHEM_TEST} photochem)
    target_include_directories(${PHOTOCHEM_TEST} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
  endforeach()

  add_executable(benchmark_block_thomas benchmark_block_thomas.f90)
  target_link_libraries(benchmark_block_thomas photochem)
  target_include_directories(benchmark_block_thomas PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
endif()
