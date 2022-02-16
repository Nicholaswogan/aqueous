
add_library(aqueous_wrapper aqueous_wrapper.f90)
target_link_libraries(aqueous_wrapper aqueous)

add_cython_target(aqueous_pyx _aqueous.pyx)
add_library(_aqueous MODULE ${aqueous_pyx})
python_extension_module(_aqueous)

target_link_libraries(_aqueous aqueous_wrapper)
target_include_directories(_aqueous PUBLIC ${NumPy_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
