add_subdirectory(c_api)
add_subdirectory(cython)

if (SKBUILD)
  install(TARGETS _clima DESTINATION photochem)
else()
  install(TARGETS _clima DESTINATION "${CMAKE_SOURCE_DIR}/photochem")
endif()
