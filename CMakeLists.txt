cmake_minimum_required(VERSION 3.14)

project(gibbs_minimizer VERSION 0.1.0 LANGUAGES Fortran C)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

include(cmake/CPM.cmake)

CPMAddPackage("gh:Nicholaswogan/fortran-yaml-c@0.1.0")
CPMAddPackage("gh:grimme-lab/nlopt-f@0.1.0")

add_subdirectory(src)

add_executable(main test/test.f90)
target_link_libraries(main gibbs)
target_include_directories(main PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")

find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)
find_package(Cython REQUIRED)

add_subdirectory(gibbs)

install(TARGETS _gibbs  DESTINATION ${CMAKE_SOURCE_DIR}/gibbs)