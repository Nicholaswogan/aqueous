module gibbs_wrapper
  use iso_c_binding
  implicit none

contains
  
  subroutine gibbs_gibbs_energy(species_len, species, T, P, err_len, err, G) bind(c)
    use gibbs, only: gibbs_energy
    integer(c_int64_t), intent(in) :: species_len
    character(kind=c_char), intent(in) :: species(species_len)
    real(c_double), intent(in) :: T
    real(c_double), intent(in) :: P
    integer(c_int64_t), intent(out) :: err_len
    type(c_ptr), intent(out) :: err
    real(c_double), intent(out) :: G
  
    character(len=:), allocatable, target :: err_f
    character(len=:), pointer :: err_p
    character(len=species_len) :: species_f
    
    call copy_string_ctof(species, species_f)
    
    G = gibbs_energy(species_f, T, P, err_f)
    if (allocated(err_f)) then
      err_len = len(err_f)
      allocate(character(len=err_len)::err_p)
      err_p = err_f
      err = c_loc(err_p)
    else
      err_len = 0
      err = c_null_ptr
    endif

  end subroutine
  
  subroutine gibbs_gibbs_energy_err(err_len, err_cp, err) bind(c)
    integer(c_int64_t), intent(in) :: err_len
    type(c_ptr), intent(in), value :: err_cp
    character(kind=c_char), intent(out) :: err(err_len)
    
    character(kind=c_char), pointer :: err_p(:)
    
    call c_f_pointer(err_cp, err_p, [err_len])
    err = err_p
    deallocate(err_p)
    
  end subroutine
  
  subroutine gibbs_load_spronsbl(path_len, path) bind(c)
    use gibbs, only: load_spronsbl
    integer(c_int64_t), intent(in) :: path_len
    character(kind=c_char), intent(in) :: path(path_len)
    
    character(len=path_len) :: path_f
    
    call copy_string_ctof(path, path_f)
    
    call load_spronsbl(path_f)
  end subroutine
  
  subroutine copy_string_ctof(stringc,stringf)
    ! utility function to convert c string to fortran string
    character(len=*), intent(out) :: stringf
    character(c_char), intent(in) :: stringc(:)
    integer j
    stringf = ''
    char_loop: do j=1,min(size(stringc),len(stringf))
       if (stringc(j)==c_null_char) exit char_loop
       stringf(j:j) = stringc(j)
    end do char_loop
  end subroutine copy_string_ctof

  subroutine copy_string_ftoc(stringf,stringc)
    ! utility function to convert c string to fortran string
    character(len=*), intent(in) :: stringf
    character(c_char), intent(out) :: stringc(:)
    integer j,n
    n = len_trim(stringf)   
    do j=1,n    
      stringc(j) = stringf(j:j)   
    end do
    stringc(n+1) = c_null_char
  end subroutine copy_string_ftoc
  
end module