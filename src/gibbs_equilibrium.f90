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
    
    real(dp), allocatable :: DG(:)
    real(dp) :: G_init
    real(dp), allocatable :: n_init(:)
    real(dp), allocatable :: m_init(:)
    real(dp) :: G_opt
    real(dp), allocatable :: n_opt(:)
    real(dp), allocatable :: m_opt(:)
    real(dp), allocatable :: atoms_init(:)
    real(dp), allocatable :: T
    real(dp), allocatable :: P
    
    real(dp) :: xtol = 1.0e-6_dp
    real(dp) :: conserv_tol = 1.0e-9_dp
    real(dp) :: lb = 1.0e-50_dp
    real(dp) :: ub = 5.0_dp
    real(dp) :: maxtime = 5.0_dp
    character(len=STR_LEN) :: algorithm = "LD_MMA"
      
  contains
    procedure :: init => AqueousSolution_init
    procedure :: equilibrate => AqueousSolution_equilibrate
  end type

contains
  
  subroutine AqueousSolution_equilibrate(self, m, T, P, err)
    use gibbs_constants, only: Rgas, mu_H2O
    use gibbs_database, only: gibbs_energy_eval
    ! use nlopt_enum, only : NLOPT_SUCCESS
    
    use nlopt_wrap, only : create, nlopt_func, nlopt_mfunc
    use nlopt_enum, only : NLOPT_SUCCESS, algorithm_from_string, result_to_string
    
    class(AqueousSolution), intent(inout), target :: self
    real(dp), intent(inout) :: m(:)
    real(dp), intent(in) :: T
    real(dp), intent(in) :: P
    character(len=:), allocatable, intent(out) :: err
    
    real(dp), parameter :: kg_H2O = 1.0_dp
    real(dp), parameter :: mol_H2O = kg_H2O/mu_H2O
    
    integer :: i, stat, algorithm
    logical :: found
    real(dp) :: minf
    type(nlopt_opt) :: opt, opt_other
    
    real(dp), allocatable :: tol(:), lb(:), ub(:), w(:)
    type(nlopt_func) :: f
    type(nlopt_mfunc) :: fc
    
    if (size(m) /= self%d%nsp-1) then
      err = 'Input "m" has the wrong size.'
      return
    endif
  
    self%m_init = m
    self%n_init(1) = mol_H2O
    self%n_init(2:) = m*kg_H2O
    
    do i = 1,self%d%natoms
      self%atoms_init(i) = sum(self%d%species_atoms(i,:)*self%n_init)
    enddo
    
    self%T = T
    self%P = P
    
    do i = 1,self%d%nsp
      call gibbs_energy_eval(self%d%thermo(i), T, P, found, self%DG(i))
      if (.not. found) then
        err = 'Species "'//trim(self%d%species_names(i))//'" has no thermodynamic data for input temperature.'
        return
      endif
    enddo
    
    ! setup nlopt 
    call create(opt, algorithm_from_string("AUGLAG_EQ"), self%d%nsp)
    
    allocate(w(self%d%nsp))
    w = 1.0_dp
    w(1) = 0.0_dp
    call opt%set_xtol_rel(self%xtol)
    call opt%set_x_weights(w)
    call opt%set_maxtime(self%maxtime)
    
    ! Set Optimizer
    algorithm = algorithm_from_string(trim(self%algorithm))
    if (algorithm == -1) then
      err = "Not a valid algorithm: "//trim(self%algorithm)
      return
    endif
    call create(opt_other, algorithm, self%d%nsp)
    call opt_other%set_xtol_rel(self%xtol)
    call opt_other%set_x_weights(w)
    call opt_other%set_maxtime(self%maxtime)
    
    call opt%set_local_optimizer(opt_other, stat)
    if (stat /= NLOPT_SUCCESS) then
      err = "NLOPT setup failed: "//result_to_string(stat)
      return
    endif
    
    f = nlopt_func(AqueousSolution_obj, self)
    fc = nlopt_mfunc(AqueousSolution_con, self)
    
    call opt%set_min_objective(f, stat)
    if (stat /= NLOPT_SUCCESS) then
      err = "NLOPT setup failed: "//result_to_string(stat)
      return
    endif
    
    allocate(lb(self%d%nsp))
    lb(1) = 30.0_dp
    lb(2:) = self%lb
    call opt%set_lower_bounds(lb, stat)
    if (stat /= NLOPT_SUCCESS) then
      err = "NLOPT setup failed: "//result_to_string(stat)
      return
    endif
    
    allocate(ub(self%d%nsp))
    ub(1) = 70.0_dp
    ub(2:) = self%ub
    call opt%set_upper_bounds(ub, stat)
    if (stat /= NLOPT_SUCCESS) then
      err = "NLOPT setup failed: "//result_to_string(stat)
      return
    endif
    
    allocate(tol(self%d%natoms))
    tol = self%conserv_tol
    call opt%add_equality_mconstraint(self%d%natoms, fc, tol, stat)
    if (stat /= NLOPT_SUCCESS) then
      err = "NLOPT setup failed: "//result_to_string(stat)
      return
    endif
    
    self%G_init = AqueousSolution_obj(x=self%n_init, func_data=self)
    self%n_opt = self%n_init
    call opt%optimize(self%n_opt, minf, stat)
    if (stat < NLOPT_SUCCESS .or. stat == 6) then
      err = "NLOPT optimization failed: "//result_to_string(stat)
      return
    endif
    
    self%G_opt = minf
    self%m_opt = self%n_opt(2:)/(self%n_opt(1)*mu_H2O)
    m = self%m_opt
    
  end subroutine

  function AqueousSolution_obj(x, gradient, func_data) result(f)
    use gibbs_constants, only: Rgas, mu_H2O
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout), optional :: gradient(:)
    class(*), intent(in), optional :: func_data
    real(dp) :: f
    
    type(AqueousSolution), pointer :: s
    integer :: i
    
    select type(func_data)
    type is(AqueousSolution)
      s => func_data
    end select
    
    if (present(gradient)) then
      gradient(1) = s%DG(1) + (Rgas*s%T/x(1))*sum(-x(2:))
      do i = 2,size(gradient)
        gradient(i) = s%DG(i) + Rgas*s%T*log(x(i)/x(1)/mu_H2O)
      enddo
    endif
    
    f = x(1)*(s%DG(1) - (sum(x(2:))/x(1))*Rgas*s%T) + sum(x(2:)*(s%DG(2:) + Rgas*s%T*log(x(2:)/x(1)/mu_H2O)))
    
  end function
  
  subroutine AqueousSolution_con(result, x, gradient, func_data)
      real(dp), intent(inout) :: result(:)
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout), optional :: gradient(:,:)
      class(*), intent(in), optional :: func_data
      
      type(AqueousSolution), pointer :: s
      integer :: i, j, n, m
      
      select type(func_data)
      type is(AqueousSolution)
        s => func_data
      end select
      
      m = s%d%natoms
      n = s%d%nsp
      
      if (present(gradient)) then
        do j = 1,m
          do i = 1,n
            gradient(i,j) = s%d%species_atoms(j,i)
          enddo
        enddo
      endif
      
      do i = 1,m
        result(i) = sum(s%d%species_atoms(i,:)*x)
      enddo
      result = result - s%atoms_init
      

    end subroutine
  
  subroutine AqueousSolution_init(self, species, err)
    use nlopt_wrap, only : create, nlopt_func, nlopt_mfunc
    use nlopt_enum, only : algorithm_from_string, NLOPT_SUCCESS
    
    class(AqueousSolution), intent(inout), target :: self
    character(len=STR_LEN), intent(in) :: species(:)
    character(len=:), allocatable, intent(out) :: err
    
    integer :: stat, i
    real(dp), allocatable :: tol(:), lb(:), ub(:)
    type(nlopt_func) :: f
    type(nlopt_mfunc) :: fc
    
    call process_species(species, self%d, err)
    if (allocated(err)) return
    
    ! other allocations
    allocate(self%DG(self%d%nsp))
    allocate(self%n_init(self%d%nsp))
    allocate(self%m_init(self%d%nsp-1))
    allocate(self%n_opt(self%d%nsp))
    allocate(self%m_opt(self%d%nsp-1))
    allocate(self%atoms_init(self%d%natoms))

  end subroutine
  
  subroutine process_species(species, dat, err)
    use gibbs_database, only: find_species_ind
    use gibbs_types, only: as
    
    character(len=STR_LEN), intent(in) :: species(:)
    type(GibbsData), intent(out) :: dat
    character(len=:), allocatable, intent(out) :: err
    
    character(len=:), allocatable :: dups
    integer :: i, j, ind, k
    integer, allocatable :: inds(:)
    logical, allocatable :: existing_atoms(:)
    
    dat%nsp = size(species) + 1
    allocate(dat%species_names(dat%nsp))
    dat%species_names(1) = "H2O"
    dat%species_names(2:) = species
    
    allocate(inds(dat%nsp))
    
    do j = 1,dat%nsp
      ind = find_species_ind(dat%species_names(j), as%species_names, as%alt_names, err)
      if (allocated(err)) return
      inds(j) = ind
    enddo
    
    ! Check that all inds are unique
    ind = findloc(species, "H2O", 1)
    if (ind /= 0) then
      err = 'Do not include "H2O" as a species. It is added automatically.'
      return
    endif
    
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
    
    allocate(dat%thermo(dat%nsp))
    do j = 1,dat%nsp
      dat%thermo(j) = as%thermo(inds(j))
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

end module