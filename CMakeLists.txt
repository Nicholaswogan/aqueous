cmake_minimum_required(VERSION 3.14)

project(AQUEOUS LANGUAGES Fortran C VERSION "0.1.0")

# options
option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

option(BUILD_PYTHON_AQUEOUS "if ON, then will build a python
version via Cython" OFF)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wunused")
  endif()
endif()


include(cmake/CPM.cmake)
CPMAddPackage("gh:Nicholaswogan/fortran-yaml-c@0.1.2")

CPMAddPackage(
  NAME nlopt-f
  VERSION 0.1.1
  URL "https://github.com/grimme-lab/nlopt-f/archive/8e36c15a6d83e11d615702ae5865c7ebe447e4de.zip"
  EXCLUDE_FROM_ALL ON
)

add_subdirectory(src)

add_executable(test_aqueous test/test_aqueous.f90)
target_link_libraries(test_aqueous aqueous)
target_include_directories(test_aqueous PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

if (BUILD_PYTHON_AQUEOUS)
  
  if (NOT SKBUILD)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
  endif()

  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Cython REQUIRED)

  add_subdirectory(aqueous)

  if (SKBUILD)
    install(TARGETS _aqueous  DESTINATION aqueous)
  else()
    install(TARGETS _aqueous  DESTINATION ${CMAKE_SOURCE_DIR}/aqueous)
  endif()
  
endif ()