add_subdirectory(c_api)
add_subdirectory(cython)

if (SKBUILD)
  install(TARGETS _equilibrate DESTINATION photochem)
else()
  install(TARGETS _equilibrate DESTINATION "${CMAKE_SOURCE_DIR}/photochem")
endif()
