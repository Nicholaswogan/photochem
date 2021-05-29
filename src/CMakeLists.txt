cmake_minimum_required(VERSION "3.19.6")


project(photochem LANGUAGES Fortran)

add_library(photochem 
            photochem_const.f90 
        		photochem_types.f90 
        		photochem_input.f90 
        		photochem_wrk.f90 
        		photochem_vars.f90 
        		photochem_data.f90 
        		photochem_radtran.f90 
)

set_property(TARGET photochem PROPERTY Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

target_include_directories(photochem
   PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../dependencies/modules>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
)

target_link_directories(photochem
   PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../dependencies/lib>
)

target_link_libraries(photochem PRIVATE yaml binning stringifor)

add_executable(photo.run main.f90)

target_include_directories(photo.run PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_link_libraries(photo.run PRIVATE photochem)