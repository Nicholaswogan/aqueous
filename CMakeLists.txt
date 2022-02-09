cmake_minimum_required(VERSION 3.14)

project(gibbs_minimizer VERSION 0.1.0 LANGUAGES Fortran C)

# options
option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

option(BUILD_PYTHON_GIBBS "if ON, then will build a python
version via Cython" OFF)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

include(cmake/CPM.cmake)

CPMAddPackage("gh:Nicholaswogan/fortran-yaml-c@0.1.0")
CPMAddPackage("gh:grimme-lab/nlopt-f@0.1.0")

add_subdirectory(src)

add_executable(main test/test.f90)
target_link_libraries(main gibbs)
target_include_directories(main PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

if (BUILD_PYTHON_GIBBS)
  
  if (NOT SKBUILD)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
  endif()

  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Cython REQUIRED)

  add_subdirectory(gibbs)

  if (SKBUILD)
    install(TARGETS _gibbs  DESTINATION gibbs)
  else()
    install(TARGETS _gibbs  DESTINATION ${CMAKE_SOURCE_DIR}/gibbs)
  endif()
  
endif ()