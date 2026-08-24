add_subdirectory(c_api)
add_subdirectory(cython)

if (SKBUILD)
  install(TARGETS _photochem DESTINATION photochem)
else()
  install(TARGETS _photochem DESTINATION "${CMAKE_SOURCE_DIR}/photochem")
endif()
