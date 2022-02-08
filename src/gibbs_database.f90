module gibbs_database
  use gibbs_constants, only: dp, STR_LEN
  implicit none
  
contains

  function gibbs_energy_coeff(coef, T, P) result(G)
    real(dp), intent(in) :: coef(10), T, P
    real(dp) :: G
    
    real(dp), parameter :: Tr = 298.15_dp &
                         , Pr = 1.0_dp &
                         , Psi = 2600.0_dp &
                         , Theta = 228.0_dp &
                         , Y = -5.81e-5_dp
                         
    real(dp) :: Gr, Hr, Sr, a1, a2, a3, a4, c1, c2, w
    real(dp) :: G1, G2, G3, G4, G5, G6, G7, G8
    real(dp) :: h1, h2, h3, h4
    
    Gr = coef(1)
    Hr = coef(2)
    Sr = coef(3)
    a1 = coef(4)/10.0_dp
    a2 = coef(5)*100.0_dp
    a3 = coef(6)
    a4 = coef(7)*10000.0_dp
    c1 = coef(8)
    c2 = coef(9)*10000.0_dp
    w = coef(10)*100000.0_dp
    
    G1 = Gr
    G2 = -1.0_dp*Sr*(T-Tr)
    G3 = -1.0_dp*c1*(T*log(T/Tr)-T+Tr)
    G4 = a1*(P-Pr)
    
    h1 = log( (Psi+P) / (Psi + Pr))

    G5 = a2*h1

    h2 = (1.0_dp/(T-Theta))-(1.0_dp/(Tr-Theta))
    h3 = (Theta-T)/Theta
    h4 = log(Tr*(T-Theta)/(T*(Tr-Theta)))

    G6 = -1.0_dp*c2*(h2*h3-T/(Theta*Theta)*h4)
    G7 = (1.0_dp/(T-Theta))*(a3*(P-Pr)+a4*h1)
    G8 = w*Y*(T-Tr)

    G = 4.184_dp*(G1+G2+G3+G4+G5+G6+G7+G8)
  
  end function
  
  function gibbs_energy(species, T, P, err) result(G)
    use gibbs_types, only: as
    character(len=*), intent(in) :: species
    real(dp), intent(in) :: T, P
    character(len=:), allocatable, intent(out) :: err
    
    real(dp) :: G
    integer :: ind
    
    ind = find_species_ind(species, as%species_names, as%alt_names, err)
    G = gibbs_energy_coeff(as%coeffs(:,ind), T, P)
  end function
  
  function find_species_ind(species, species_names, alt_names, err) result(ind)
    character(len=*), intent(in) :: species
    character(len=STR_LEN), intent(in) :: species_names(:)
    character(len=STR_LEN), intent(in) :: alt_names(:)
    character(len=:), allocatable, intent(out) :: err
    
    integer :: ind
    
    character(len=:), allocatable :: dups
    integer :: i
    
    ind = findloc(species_names, species, 1)
    if (ind == 0) then
      ind = findloc(alt_names, species, 1)
      if (ind == 0) then
        err = '"'//trim(species)//'" is not in the data base.'
        return 
      else
        ! check that entry is unique
        do i = 1,size(species_names)
          if (i == ind) cycle
          if (alt_names(ind) == alt_names(i)) then
            if (.not. allocated(dups)) then
              dups = '{unique-name: "'//trim(species_names(ind))//'", alt-name: "'//trim(alt_names(ind))//'"}'
            endif
            dups = dups//', {unique-name: "'//trim(species_names(i))//'", alt-name: "'//trim(alt_names(i))//'"}'
          endif
        enddo
        if (allocated(dups)) then
          err = '"'//trim(species)//'" is not a unique entry in the database. '//dups
          return
        endif
      endif
    endif
    
  end function
  
  subroutine load_spronsbl(path)
    use gibbs_types, only: as
    use fortran_yaml_c, only: parse, error_length
    use yaml_types, only: type_node, type_dictionary
    
    character(len=*), intent(in) :: path
    
    class(type_node), pointer :: root
    character(len=error_length) :: error
    character(len=:), allocatable :: err

    root => parse(trim(path)//"sprons96.yaml", error = error)
    if (error/='') then
      print*,trim(error)
      stop 1
    endif
    
    select type (root)
    class is (type_dictionary)
      call process_spronsbl(root, as, err)
    end select
    
    call root%finalize()
    deallocate(root)
    
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    
  end subroutine
  
  subroutine process_spronsbl(root, as, err)
    use gibbs_types, only: AllData
    use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                          type_list, type_list_item, type_scalar, type_key_value_pair
                    
    type(type_dictionary), intent(in) :: root
    type(AllData), intent(out) :: as
    character(len=:), allocatable, intent(out) :: err
    
    type(type_key_value_pair), pointer :: key_value_pair
    type(type_list_item), pointer :: item, item1
    type(type_list), pointer :: atoms, species, list1
    type(type_dictionary), pointer :: dict
    type (type_error), pointer :: io_err
    logical :: success
    integer :: j, i

    atoms => root%get_list("atoms",required=.true.,error = io_err)
    if (associated(io_err)) then; err = trim(io_err%message); return; endif
    species => root%get_list("species",required=.true.,error = io_err)
    if (associated(io_err)) then; err = trim(io_err%message); return; endif
    
    as%natoms = atoms%size()
    allocate(as%atoms_names(as%natoms))
    
    j = 1
    item => atoms%first
    do while(associated(item))
      select type (scal => item%node)
      class is (type_scalar)
        as%atoms_names(j) = trim(scal%string)
      end select
      item => item%next
      j = j + 1
    enddo
    
    as%nsp = species%size()
    allocate(as%species_names(as%nsp))
    allocate(as%alt_names(as%nsp))
    allocate(as%species_atoms(as%natoms, as%nsp))
    allocate(as%coeffs(10, as%nsp))
    
    j = 1
    item => species%first
    do while(associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        as%species_names(j) = trim(element%get_string("name",error = io_err))
        if (associated(io_err)) then; err = trim(io_err%message); return; endif
        as%alt_names(j) = trim(element%get_string("alt-name",error = io_err))
        if (associated(io_err)) then; err = trim(io_err%message); return; endif
        
        dict => element%get_dictionary("composition",required=.true.,error = io_err)
        if (associated(io_err)) then; err = trim(io_err%message); return; endif
        do i = 1, as%natoms
          as%species_atoms(i,j) = dict%get_integer(trim(as%atoms_names(i)), default=0, error = io_err)
          if (associated(io_err)) then; err = trim(io_err%message); return; endif
        enddo
        
        list1 => element%get_list("thermo",required=.true.,error = io_err)
        if (associated(io_err)) then; err = trim(io_err%message); return; endif
        i = 1
        item1 => list1%first
        do while(associated(item1))
          select type (scal => item1%node)
          class is (type_scalar)
            as%coeffs(i,j) = scal%to_real(default=1.0_dp,success=success)
            if (.not. success) then
              err = "Issue converting string or scalar!"
              return
            endif
          end select
          item1 => item1%next
          i = i + 1
        enddo
      end select
      item => item%next
      j = j + 1
    enddo
    
  end subroutine
  
  ! useless!
  pure function dielectric(T_in, P_in) result(di)
    real(dp), intent(in) :: T_in, P_in
    real(dp) :: di
    
    real(dp), parameter :: a(8) = [-22.5713_dp, -0.032066_dp, -0.00028568_dp,&
                                   0.0011832_dp, 0.000027895_dp, -0.00000001476_dp,&
                                   2300.64_dp, -0.13476_dp]
    real(dp), parameter :: D0 = 4.476150
    real(dp) :: T, P
    
    T = T_in - 273.15_dp
    P = P_in - 1.0_dp
    
    di = exp(D0 &
            + 2_dp*a(1)*P &
            + 2_dp*a(2)*P*T &
            + 2_dp*a(3)*P*T**2_dp &
            + 2_dp*a(4)*P**2_dp &
            + 2_dp*a(5)*P**2_dp*T &
            + 2_dp*a(6)*P**2_dp*T**2_dp &
            + 2_dp*a(7)*T &
            + 2_dp*a(8)*T**2_dp)

  end function
  
end module