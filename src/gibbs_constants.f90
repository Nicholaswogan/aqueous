module gibbs_constants
  use iso_fortran_env, only: dp => real64
  implicit none
  integer, parameter :: STR_LEN = 20
  real(dp), parameter :: Rgas = 8.31446261815324e0_dp
  real(dp), parameter :: mu_H2O = 0.01801528_dp
end module