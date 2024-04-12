
set(PHOTOCHEM_TESTS testrun testevo memtest memtest_evo)

if (BUILD_EXECUTABLES)
  foreach(PHOTOCHEM_TEST ${PHOTOCHEM_TESTS})
    add_executable(${PHOTOCHEM_TEST} "${PHOTOCHEM_TEST}.f90")
    target_link_libraries(${PHOTOCHEM_TEST} photochem)
    target_include_directories(${PHOTOCHEM_TEST} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
  endforeach()
endif()