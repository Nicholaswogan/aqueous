module gibbs_equilibrium
  use gibbs_constants, only: dp, STR_LEN
  use gibbs_types, only: GibbsData
  use nlopt_wrap, only : nlopt_opt
  implicit none
  
  type :: constraint_data
    real(dp) :: d(2)
  end type
  
  type :: AqueousSolution
    type(GibbsData) :: d
    real(dp), allocatable :: x(:)
    type(nlopt_opt) :: opt
  contains
    procedure :: init => AqueousSolution_init
  end type

contains

  function AqueousSolution_obj(x, gradient, func_data) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout), optional :: gradient(:)
    class(*), intent(in), optional :: func_data
    real(dp) :: f
    
    type(AqueousSolution), pointer :: s
    
    select type(func_data)
    type is(AqueousSolution)
      s => func_data
    end select
    
    if (present(gradient)) then
      
    endif
    
    f = 1.0_dp
    
  end function
  
  subroutine AqueousSolution_init(self, infile, err)
    use nlopt_wrap, only : create, nlopt_func
    use nlopt_enum, only : algorithm_from_string
    use fortran_yaml_c, only: parse, error_length
    use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                          type_list, type_list_item, type_scalar
    
    class(AqueousSolution), intent(inout), target :: self
    character(len=*), intent(in) :: infile
    character(len=:), allocatable, intent(out) :: err

    class(type_node), pointer :: root
    character(len=error_length) :: error
    
    real(dp) :: lb(2),x(2), minf
    integer :: stat, i
    real(dp), parameter :: xtol = 1.0e-4_dp
    
    root => parse(infile, error = error)
    if (error/='') then
      print*,trim(error)
      stop 1
    endif
  
    select type (root)
    class is (type_list)
      call process_datafile(root, self%d, self%x, err)
    end select
  
    call root%finalize()
    deallocate(root)
    
    ! setup nlopt 
    ! call create(self%opt, algorithm_from_string("LD_MMA"), self%d%nsp)
    
    ! set constrains and functions
    
  end subroutine
  
  subroutine process_datafile(root, dat, x, err)
    use gibbs_database, only: find_species_ind
    use gibbs_types, only: as
    use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                          type_list, type_list_item, type_scalar, type_key_value_pair
                    
    type(type_list), intent(in) :: root
    type(GibbsData), intent(out) :: dat
    real(dp), allocatable, intent(out) :: x(:)
    character(len=:), allocatable, intent(out) :: err
    
    character(len=:), allocatable :: dups
    type(type_list_item), pointer :: item
    type (type_error), pointer :: io_err
    integer :: i, j, ind, k
    integer, allocatable :: inds(:)
    logical, allocatable :: existing_atoms(:)
    
    dat%nsp = root%size()
    allocate(dat%species_names(dat%nsp))
    allocate(x(dat%nsp))
    
    j = 1
    item => root%first
    do while(associated(item))
      select type (dict => item%node)
      class is (type_dictionary)
        dat%species_names(j) = dict%get_string("name",error=io_err)
        if (associated(io_err)) then; err = trim(io_err%message); return; endif
          
        x(j) = dict%get_real("concentration",default=0.0_dp,error=io_err)
      class default
        err = "Each entry in the data file must be a dictionary."
        return
      end select
      item => item%next
      j = j + 1
    enddo
    
    allocate(inds(dat%nsp))
    
    do j = 1,dat%nsp
      ind = find_species_ind(dat%species_names(j), as%species_names, as%alt_names, err)
      if (allocated(err)) return
      inds(j) = ind
    enddo
    
    ! Check that all inds are unique
    do i = 1,dat%nsp
      do j = i+1,dat%nsp
        if (inds(i) == inds(j)) then
          err = 'These two species are the same: "'//&
                trim(dat%species_names(j))//'" and "'//&
                trim(dat%species_names(j))//'"'
          return
        endif
      enddo
    enddo
    
    allocate(dat%coeffs(10,dat%nsp))
    do j = 1,dat%nsp
      dat%coeffs(:,j) = as%coeffs(:,inds(j))
    enddo
    
    allocate(existing_atoms(as%natoms))
    existing_atoms = .false.
    
    do j = 1,dat%nsp
      do i = 1,as%natoms
        if (as%species_atoms(i,inds(j)) /= 0) then
          existing_atoms(i) = .true.
        endif
      enddo
    enddo
    
    ! count number of true
    dat%natoms = 0
    do i = 1,as%natoms
      if (existing_atoms(i)) dat%natoms = dat%natoms + 1
    enddo
    
    allocate(dat%atoms_names(dat%natoms))
    allocate(dat%species_atoms(dat%natoms,dat%nsp))
    k = 1
    do i = 1,as%natoms
      if (existing_atoms(i)) then
        dat%atoms_names(k) = as%atoms_names(i)
        do j = 1,dat%nsp
          dat%species_atoms(k,j) = as%species_atoms(i,inds(j))
        enddo
        k = k + 1
      endif
    enddo
    
  end subroutine

  subroutine tester1()
    use nlopt_wrap, only : nlopt_func, create, destroy
    use nlopt_enum, only : NLOPT_SUCCESS, algorithm_from_string
    implicit none
    type(nlopt_opt) :: opt
    real(dp) :: lb(2), x(2), minf
    integer :: stat, i
    type(constraint_data), target :: d1, d2
    real(dp), parameter :: xtol = 1.0e-4_dp

    call create(opt, algorithm_from_string("LD_MMA"), 2)

    call opt%get_lower_bounds(lb)
    lb(2) = 0.0_dp
    call opt%set_lower_bounds(lb)

    d1%d = [+2.0_dp, +0.0_dp]
    d2%d = [-1.0_dp, +1.0_dp]
    associate(&
        & f => nlopt_func(myfunc), &
        & fc1 => nlopt_func(myconstraint, d1), &
        & fc2 => nlopt_func(myconstraint, d2))
      call opt%set_min_objective(f)
      call opt%add_inequality_constraint(fc1, 1.0e-8_dp)
      call opt%add_inequality_constraint(fc2, 1.0e-8_dp)
    end associate
    
    call opt%set_xtol_rel(xtol)
    x = [1.234_dp, 5.678_dp]
    call opt%optimize(x, minf, stat)

    if (stat < NLOPT_SUCCESS) then
      write(*, '(a)') "NLopt failed!"
      stop 1
    endif

    write(*, '(a, *(1x, g0))') "Found minimum at", x
    write(*, '(a, *(1x, g0))') "Minimum value is", minf

    call destroy(opt)
  end subroutine
  
  function myfunc(x, gradient, func_data) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout), optional :: gradient(:)
    class(*), intent(in), optional :: func_data
    real(dp) :: f

    if (present(gradient)) then
      gradient(1) = 0.0_dp
      gradient(2) = 0.5_dp / sqrt(x(2))
    endif
    f = sqrt(x(2))
    
  end function myfunc
  
  function myconstraint(x, gradient, func_data) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout), optional :: gradient(:)
    class(*), intent(in), optional :: func_data
    real(dp) :: f

    select type(func_data)
    type is(constraint_data)
      associate(a => func_data%d(1), b => func_data%d(2))
        if (present(gradient)) then
          gradient(1) = 3.0_dp * a * (a*x(1) + b)**2
          gradient(2) = -1.0_dp
        endif
        f = (a*x(1) + b)**3 - x(2)
      end associate
    end select
  end function myconstraint
  
end module