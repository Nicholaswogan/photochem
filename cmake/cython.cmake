function(add_cython_module target source)
  set(cython_source "${CMAKE_CURRENT_BINARY_DIR}/${target}.c")

  add_custom_command(
    OUTPUT "${cython_source}"
    COMMAND
      "${Python_EXECUTABLE}" -m cython --3
      --directive embedsignature=True
      --directive embedsignature.format=clinic
      --module-name "photochem.${target}"
      --output-file "${cython_source}"
      "${CMAKE_CURRENT_SOURCE_DIR}/${source}"
    DEPENDS "${source}" ${ARGN}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Cythonizing ${source}"
    VERBATIM
  )

  Python_add_library(${target} MODULE WITH_SOABI "${cython_source}")
  target_link_libraries(${target} PRIVATE Python::NumPy)
endfunction()
