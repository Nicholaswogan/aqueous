module gibbs_types
  use gibbs_constants, only: dp, STR_LEN
  implicit none
  
  type :: ThermodynamicData
    integer :: dtype
    integer :: ntemps
    real(dp), allocatable :: temps(:)
    real(dp), allocatable :: data(:,:)
  end type
  
  type :: GibbsData
    integer :: natoms
    integer :: nsp
    character(len=STR_LEN), allocatable :: species_names(:)
    character(len=STR_LEN), allocatable :: atoms_names(:)
    integer, allocatable :: species_atoms(:,:)
    type(ThermodynamicData), allocatable :: thermo(:)
  end type
  
  type, extends(GibbsData) :: AllData
    character(len=STR_LEN), allocatable :: alt_names(:)
  end type
  
  type(AllData) :: as
  
end module