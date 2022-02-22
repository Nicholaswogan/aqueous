program test
  use aqueous, only: gibbs_energy, load_spronsbl, dp, AqueousSolution, STR_LEN
  implicit none
  
  character(len=:), allocatable :: err
  real(dp) :: G
  type(AqueousSolution) :: a
  character(len=:), allocatable :: species(:)
  real(dp), allocatable :: m(:)
  
  call load_spronsbl("../aqueous/data/")
  
  G = gibbs_energy("NH3,aq", 298.0_dp, 1.0_dp, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  species = [character(len=STR_LEN) :: "H+", "OH-"]
  allocate(m(2))
  
  a = AqueousSolution(species, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  m = [1.0e-40_dp,1.0e-40_dp]
  call a%equilibrate(m, 298.0_dp, 1.0_dp, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  print*,m
  
end program