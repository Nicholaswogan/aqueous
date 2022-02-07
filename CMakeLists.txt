cmake_minimum_required(VERSION 3.14)

project(gibbs_minimizer VERSION 0.1.0 LANGUAGES Fortran C)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

include(cmake/CPM.cmake)

CPMAddPackage("gh:Nicholaswogan/fortran-yaml-c@0.1.0")
CPMAddPackage("gh:grimme-lab/nlopt-f@0.1.0")


add_library(gibbs src/gibbs.f90 src/gibbs_types.f90)
target_link_libraries(gibbs nlopt-f-lib fortran-yaml-c)

add_executable(main test/test.f90)
target_link_libraries(main gibbs)
target_include_directories(main PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})


