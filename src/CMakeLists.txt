add_library(gibbs 
  gibbs_constants.f90
  gibbs_types.f90
  gibbs_database.f90
  gibbs_equilibrium.f90
  gibbs.f90 )
target_link_libraries(gibbs nlopt-f-lib fortran-yaml-c)