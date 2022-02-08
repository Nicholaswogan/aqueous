
add_library(gibbs_wrapper gibbs_wrapper.f90)
target_link_libraries(gibbs_wrapper gibbs)

add_cython_target(gibbs_pyx _gibbs.pyx)
add_library(_gibbs MODULE ${gibbs_pyx})
python_extension_module(_gibbs)

target_link_libraries(_gibbs gibbs_wrapper)
target_include_directories(_gibbs PUBLIC ${NumPy_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
