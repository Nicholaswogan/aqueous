module aqueous_database
  use aqueous_constants, only: dp, STR_LEN
  implicit none
  
  public
  
contains
    
  function gibbs_energy(species, T, P, err) result(G)
    use aqueous_types, only: as
    character(len=*), intent(in) :: species
    real(dp), intent(in) :: T, P
    character(len=:), allocatable, intent(out) :: err
    
    real(dp) :: G
    integer :: ind
    logical :: found
    
    ind = find_species_ind(species, as%species_names, as%alt_names, err)
    if (allocated(err)) return 
    
    call gibbs_energy_eval(as%thermo(ind), T, P, found, G)
    if (.not. found) then
      err = 'Species "'//trim(species)//'" has no thermodynamic data for input temperature.'
      return
    endif
    
  end function
  
  pure subroutine gibbs_energy_eval(thermo, T, P, found, gibbs_energy)
    use aqueous_types, only: ThermodynamicData
    
    type(ThermodynamicData), intent(in) :: thermo
    real(dp), intent(in) :: T, P
    logical, intent(out) :: found
    real(dp), intent(out) :: gibbs_energy
    
    integer :: k
    
    found = .false.
    do k = 1,thermo%ntemps
      if (T >= thermo%temps(k) .and. &
          T <  thermo%temps(k+1)) then
          
        found = .true.
        if (thermo%dtype == 1) then ! aqueous
          ! check to see if liquid phase is possible
          gibbs_energy = gibbs_energy_sprons96(thermo%data(:,k), T, P)
        elseif (thermo%dtype == 2) then
          gibbs_energy = gibbs_energy_nasa9(thermo%data(:,k), T)
        endif
        
        exit
        
      endif
    enddo

  end subroutine
  
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
  
  pure function gibbs_energy_nasa9(coeffs, T) result(gibbs)
    use aqueous_constants, only: Rgas
    real(dp), intent(in) :: coeffs(9)
    real(dp), intent(in) :: T
    real(dp) :: gibbs
    
    real(dp) :: enthalpy, entropy
    
    enthalpy = (- coeffs(1)*T**(-2.0_dp) + coeffs(2)*log(T)/T &
                + coeffs(3) + coeffs(4)*T/2.0_dp + coeffs(5)*T**(2.0_dp)/3.0_dp &
                + coeffs(6)*T**(3.0_dp)/4.0_dp + coeffs(7)*T**(4.0_dp)/5.0_dp &
                + coeffs(8)/T)*T*Rgas
             
    entropy = (- coeffs(1)*T**(-2.0_dp)/2.0_dp - coeffs(2)*T**(-1.0_dp) &
               + coeffs(3)*log(T) + coeffs(4)*T + coeffs(5)*T**(2.0_dp)/2.0_dp &
               + coeffs(6)*T**(3.0_dp)/3.0_dp + coeffs(7)*T**(4.0_dp)/4.0_dp &
               + coeffs(9))*Rgas
               
    gibbs = enthalpy - T*entropy
  end function
  
  pure function gibbs_energy_sprons96(coef, T, P) result(G)
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
  
  subroutine load_spronsbl(path)
    use aqueous_types, only: as
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
    use aqueous_types, only: AllData
    use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                          type_list, type_list_item, type_scalar, type_key_value_pair
                    
    type(type_dictionary), intent(in) :: root
    type(AllData), intent(inout) :: as
    character(len=:), allocatable, intent(out) :: err
    
    type(type_key_value_pair), pointer :: key_value_pair
    type(type_list_item), pointer :: item, item1
    type(type_list), pointer :: atoms, species, list1
    type(type_dictionary), pointer :: dict
    type (type_error), pointer :: io_err
    character (len=:), allocatable :: tmp
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
    allocate(as%thermo(as%nsp))
    
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
        
        dict => element%get_dictionary("thermo",required=.true.,error = io_err)
        if (associated(io_err)) then; err = trim(io_err%message); return; endif
        
        tmp = trim(dict%get_string("model",error = io_err))
        if (associated(io_err)) then; err = trim(io_err%message); return; endif
        
        if (tmp == "sprons96") then
          call get_sprons96_thermo(dict, as%thermo(j), err)
          if (allocated(err)) return
        elseif (tmp == "nasa9-sprons96") then
          call get_nasa9_thermo(dict, as%thermo(j), err)
          if (allocated(err)) return
        else
          err = "Unknown thermodynamic type: "//tmp
          return
        endif
        
      end select
      item => item%next
      j = j + 1
    enddo
    
  end subroutine
  
  subroutine get_sprons96_thermo(dict, thermo, err)
    use aqueous_types, only: ThermodynamicData
    use yaml_types, only: type_dictionary, type_error, &
                          type_list, type_list_item, type_scalar
    
    type(type_dictionary), intent(in) :: dict
    type(ThermodynamicData), intent(out) :: thermo
    character(len=:), allocatable, intent(out) :: err
    
    type(type_list_item), pointer :: item
    type(type_list), pointer :: list
    type (type_error), pointer :: io_err
    logical :: success
    integer :: i
    
    thermo%dtype = 1
    thermo%ntemps = 1
    allocate(thermo%temps(2))
    allocate(thermo%data(10,1))
    
    list = dict%get_list("temperature-ranges",required=.true.,error = io_err)
    if (associated(io_err)) then; err = trim(io_err%message); return; endif
    i = 1
    item => list%first
    do while(associated(item))
      select type (scal => item%node)
      class is (type_scalar)
        thermo%temps(i) = scal%to_real(default=1.0_dp,success=success)
        if (.not. success) then
          err = "Issue converting string or scalar!"
          return
        endif
      end select
      item => item%next
      i = i + 1
    enddo
    
    list = dict%get_list("data",required=.true.,error = io_err)
    
    i = 1
    item => list%first
    do while(associated(item))
      select type (scal => item%node)
      class is (type_scalar)
        thermo%data(i,1) = scal%to_real(default=1.0_dp,success=success)
        if (.not. success) then
          err = "Issue converting string or scalar!"
          return
        endif
      end select
      item => item%next
      i = i + 1
    enddo

  end subroutine
  
  subroutine get_nasa9_thermo(dict, thermo, err)
    use aqueous_types, only: ThermodynamicData
    use yaml_types, only: type_dictionary, type_error, &
                          type_list, type_list_item, type_scalar
    
    type(type_dictionary), intent(in) :: dict
    type(ThermodynamicData), intent(out) :: thermo
    character(len=:), allocatable, intent(out) :: err
    
    type(type_list_item), pointer :: item, item1
    type(type_list), pointer :: list
    type (type_error), pointer :: io_err
    logical :: success
    integer :: j, i
    
    thermo%dtype = 2
    
    list = dict%get_list("temperature-ranges",required=.true.,error = io_err)
    if (associated(io_err)) then; err = trim(io_err%message); return; endif
      
    thermo%ntemps = list%size() - 1
    allocate(thermo%temps(thermo%ntemps + 1))
      
    i = 1
    item => list%first
    do while(associated(item))
      select type (scal => item%node)
      class is (type_scalar)
        thermo%temps(i) = scal%to_real(default=1.0_dp,success=success)
        if (.not. success) then
          err = "Issue converting string or scalar!"
          return
        endif
      end select
      item => item%next
      i = i + 1
    enddo
    
    list = dict%get_list("data",required=.true.,error = io_err)
  
    allocate(thermo%data(9, thermo%ntemps))
    
    i = 1
    item => list%first
    do while(associated(item))
      select type (list1 => item%node)
      class is (type_list)
        
        j = 1
        item1 => list1%first
        do while(associated(item1))
          select type (scal => item1%node)
          class is (type_scalar)
            thermo%data(j,i) = scal%to_real(default=1.0_dp,success=success)
            if (.not. success) then
              err = "Issue converting string or scalar!"
              return
            endif
          end select
          item1 => item1%next
          j = j + 1
        enddo
      end select
      item => item%next
      i = i + 1
    enddo

  end subroutine
  
end module