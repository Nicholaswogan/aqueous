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
  
  allocate(species(7))
  allocate(m(7))
  species(1) = "H+"
  species(2) = "OH-"
  species(3) = "CO2,aq"
  species(4) = "CH4"
  species(5) = "CO,AQ"
  species(6) = "CO3-2"
  species(7) = "HCO3-"
  
  call a%init(species, err)
  if (allocated(err)) then
    print*,err
    stop 1
  endif
  
  ! m = [1.0e-50_dp, 1.0e-50_dp]
  m = 1.0e-10_dp
  call a%equilibrate(m, 298.0_dp, 1.0_dp, err)
  if (allocated(err)) then
    print*,err
    ! stop 1
  endif
  
  print*,a%m_init
  print*,a%m_opt
  
end program