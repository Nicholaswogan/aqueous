program test
  use gibbs, only: gibbs_energy, load_spronsbl, dp, AqueousSolution, STR_LEN
  implicit none
  
  character(len=:), allocatable :: err
  real(dp) :: G
  type(AqueousSolution) :: a
  character(len=STR_LEN), allocatable :: species(:)
  real(dp), allocatable :: m(:)
  
  call load_spronsbl("../gibbs/data/")
  G = gibbs_energy("H2O", 350.0_dp, 1.0_dp, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  allocate(species(2))
  allocate(m(2))
  species(1) = "H+"
  species(2) = "OH-"
  
  call a%init(species, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  ! m = [1.0e-20_dp, 1.0e-20_dp]
  m = [1.0e-2_dp, 1.0e-2_dp]
  call a%equilibrate(m, 298.0_dp, 1.0_dp, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  print*,a%m_opt
  
end program