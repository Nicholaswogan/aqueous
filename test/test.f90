program test
  use gibbs, only: gibbs_energy, load_spronsbl, dp
  implicit none
  
  character(len=:), allocatable :: err
  real(dp) :: G
  
  call load_spronsbl("../gibbs/data/")
  G = gibbs_energy("CO2,aq", 298.0_dp, 1.0_dp, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  print*,G
  
end program