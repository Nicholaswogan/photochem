cmake_minimum_required(VERSION "3.12")
project(photochem LANGUAGES Fortran C CXX)

include(FortranCInterface)
FortranCInterface_VERIFY()
FortranCInterface_VERIFY(CXX)


option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

option(BUILD_STATIC_PHOTOCHEM "if ON, then will build the
static library libphotochem.a" ON)

option(BUILD_F2PY_PHOTOCHEM "if ON, then will build a python
version with f2py" OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" 
   AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wunused")
endif()

option(BUILD_WITH_OPENMP "Compile with muli-threading" OFF)
if (BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if (OpenMP_Fortran_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  endif()
endif()
message (STATUS "Building with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")

# adding src build dependencies
add_subdirectory(src)

set(PHOTOCHEM_SOURCES ${CMAKE_SOURCE_DIR}/src/photochem_const.f90
                      ${CMAKE_SOURCE_DIR}/src/photochem_radtran.f90 
                      ${CMAKE_SOURCE_DIR}/src/photochem_eqns.f90 
                      ${CMAKE_SOURCE_DIR}/src/photochem_types.f90
                      ${CMAKE_SOURCE_DIR}/src/photochem_mie.f90  
                      ${CMAKE_SOURCE_DIR}/src/photochem_input.f90
                      ${CMAKE_SOURCE_DIR}/src/photochem_object.f90
                      ${CMAKE_SOURCE_DIR}/src/photochem_object_init.f90
                      ${CMAKE_SOURCE_DIR}/src/photochem_object_rhs.f90
                      ${CMAKE_SOURCE_DIR}/src/atmos.f90)
                      
set(F2PY_SOURCES ${CMAKE_SOURCE_DIR}/src/photochem_const.f90 
                 ${CMAKE_SOURCE_DIR}/src/photochem_vars.f90 
                 ${CMAKE_SOURCE_DIR}/src/photochem_data.f90 
                 ${CMAKE_SOURCE_DIR}/src/photochem_setup.f90
                 ${CMAKE_SOURCE_DIR}/src/photochem.f90)

if (BUILD_STATIC_PHOTOCHEM)
  
  # libphotochem.a
  add_library(photochem 
              ${PHOTOCHEM_SOURCES})
  target_link_libraries(photochem fortran-yaml-cpp sundials_fcvode_mod sundials_cvode cminpack binning)
  target_include_directories(photochem PUBLIC
                            "${CMAKE_BINARY_DIR}/modules"
                            )

  # the main executable
  add_executable(photorun ${CMAKE_SOURCE_DIR}/src/main.f90)
  target_link_libraries(photorun photochem)
  target_include_directories(photorun PRIVATE ${CMAKE_BINARY_DIR}/modules)

  #  tests
  # add_executable(testrun ${CMAKE_SOURCE_DIR}/tests/test.f90)
  # target_link_libraries(testrun photochem)
  # target_include_directories(testrun PRIVATE ${CMAKE_BINARY_DIR}/modules)

  # enable_testing()
  # add_test(test_photochem testrun)
endif()

if (BUILD_F2PY_PHOTOCHEM)
  if (NOT SKBUILD)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
  endif()
  
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy)
  find_package(F2PY REQUIRED)
  
  message(STATUS "F2PY included from: ${F2PY_INCLUDE_DIR}")
  
  set(f2py_module_name "_photochem")
  set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
  
  # f2py -m photochem photochem_setup.f90 photochem_vars.f90 
  # photochem_data.f90 photochem_const.f90 photochem.f90 
  # only: setup photo_equilibrium compute_surface_fluxes out2atmosphere_txt :
  
  add_custom_target(${f2py_module_name} ALL
  DEPENDS ${F2PY_SOURCES}
  )
  
  add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    ${F2PY_SOURCES}
    --lower
    only: setup photo_equilibrium output_surface_fluxes out2in out2atmosphere_txt
          get_species_names get_atoms_names get_output_density get_output_pressure :
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${F2PY_SOURCES}
  )
  
  add_library(${generated_module_file} MODULE
              "${f2py_module_name}module.c"
              "${F2PY_INCLUDE_DIR}/fortranobject.c"
              "${f2py_module_name}-f2pywrappers2.f90"
               ${PHOTOCHEM_SOURCES})
    
  target_link_libraries(${generated_module_file} fortran-yaml-cpp sundials_fcvode_mod sundials_cvode cminpack binning)
  target_include_directories(${generated_module_file} PUBLIC
                            "${CMAKE_BINARY_DIR}/modules"
                            ${F2PY_INCLUDE_DIRS}
                            ${PYTHON_INCLUDE_DIRS})
                            
  set_target_properties(${generated_module_file} PROPERTIES SUFFIX "")
  set_target_properties(${generated_module_file} PROPERTIES PREFIX "")    
  
  if (UNIX)
    if (APPLE)
      set_target_properties(${generated_module_file} PROPERTIES
      LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
    else()
      set_target_properties(${generated_module_file} PROPERTIES
    LINK_FLAGS  '-Wl,--allow-shlib-undefined')
    endif()
  endif()
  
  
  if (SKBUILD)
    install(TARGETS ${generated_module_file} DESTINATION Photochem)
  else()
    install(TARGETS ${generated_module_file} DESTINATION ${CMAKE_SOURCE_DIR}/Photochem)
  endif()

                    
endif()



