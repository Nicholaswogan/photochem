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

if (BUILD_EXECUTABLES)
  add_subdirectory(clima)
  add_subdirectory(equilibrate)
  add_subdirectory(photochem)
endif()
