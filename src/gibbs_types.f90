module gibbs_types
  use iso_fortran_env, only: dp => real64
  implicit none
  
  type :: GibbsData
    integer :: natoms
    integer :: nsp
    character(len=20), allocatable :: species_names(:)
    character(len=20), allocatable :: atoms_names(:)
    integer, allocatable :: species_atoms(:,:)
    real(dp), allocatable :: coeffs(:,:)
  end type
  
  type, extends(GibbsData) :: AllData
    character(len=20), allocatable :: alt_names(:)
  end type
  
  type(AllData) :: as
  
end module