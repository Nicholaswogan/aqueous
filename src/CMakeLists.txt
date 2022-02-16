add_library(aqueous
  aqueous_constants.f90
  aqueous_types.f90
  aqueous_database.f90
  aqueous_equilibrium.f90
  aqueous.f90 )
target_link_libraries(aqueous nlopt-f-lib fortran-yaml-c)