cmake_minimum_required(VERSION 3.18...3.27)
project(couple_mcdc_g4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Geant4 11.3.2 REQUIRED ui_all vis_all)

# Prefer Python from the currently activated conda environment.
if(NOT Python_ROOT_DIR AND DEFINED ENV{CONDA_PREFIX})
  set(Python_ROOT_DIR "$ENV{CONDA_PREFIX}")
endif()
set(Python_FIND_VIRTUALENV FIRST)
set(Python_FIND_STRATEGY LOCATION)
find_package(Python COMPONENTS Interpreter Development.Module Development.Embed REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

set(B1_DIR ${CMAKE_CURRENT_SOURCE_DIR}/examples/B1)

file(GLOB B1_SOURCES CONFIGURE_DEPENDS ${B1_DIR}/src/*.cc)
file(GLOB B1_HEADERS CONFIGURE_DEPENDS ${B1_DIR}/include/*.hh)

add_executable(exampleB1 ${B1_DIR}/exampleB1.cc ${B1_SOURCES} ${B1_HEADERS})
target_include_directories(exampleB1 PRIVATE ${B1_DIR}/include)
target_link_libraries(exampleB1 PRIVATE ${Geant4_LIBRARIES})

add_library(core
  bridge/detector_construction.cpp
  bridge/session.cpp
  bridge/primary_bank.cpp
  bridge/source_distribution.cpp
  bridge/bank_primary_generator_action.cpp
  bridge/results.cpp
)
set_target_properties(core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/bridge
  PRIVATE
    ${Geant4_INCLUDE_DIRS}
)
target_link_libraries(core PUBLIC pybind11::headers Python::Module PRIVATE ${Geant4_LIBRARIES})

Python_add_library(geant4_bridge MODULE WITH_SOABI
  bridge/bindings.cpp
)
target_include_directories(geant4_bridge PRIVATE ${Geant4_INCLUDE_DIRS})
target_link_libraries(geant4_bridge PRIVATE core Python::Module pybind11::headers)

enable_testing()

add_executable(source_distribution_sample_test
  tests/source_distribution_sample_test.cpp
)
target_link_libraries(source_distribution_sample_test PRIVATE core pybind11::embed Python::Python)
add_test(NAME source_distribution_sample_test COMMAND source_distribution_sample_test)

add_executable(device_geometry_config_test
  tests/device_geometry_config_test.cpp
)
target_link_libraries(device_geometry_config_test PRIVATE core pybind11::embed Python::Python)
add_test(NAME device_geometry_config_test COMMAND device_geometry_config_test)

add_executable(device_geometry_mass_test
  tests/device_geometry_mass_test.cpp
)
target_include_directories(device_geometry_mass_test PRIVATE ${Geant4_INCLUDE_DIRS})
target_link_libraries(device_geometry_mass_test PRIVATE core pybind11::embed Python::Python ${Geant4_LIBRARIES})
add_test(NAME device_geometry_mass_test COMMAND device_geometry_mass_test)

add_test(NAME seu_scoring_test
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/seu_scoring_test.py
          ${CMAKE_CURRENT_BINARY_DIR})

set(B1_SCRIPTS
  exampleB1.in
  exampleB1.out
  init_vis.mac
  run1.mac
  run2.mac
  vis.mac
)

foreach(_script ${B1_SCRIPTS})
  configure_file(${B1_DIR}/${_script} ${CMAKE_CURRENT_BINARY_DIR}/${_script} COPYONLY)
endforeach()
