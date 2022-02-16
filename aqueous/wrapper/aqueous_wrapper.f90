module aqueous_wrapper
  use iso_c_binding
  implicit none

contains
  
  subroutine aqueous_alloc_aqueoussolution(ptr) bind(c)
    use aqueous, only: dtype => aqueoussolution
    type(c_ptr), intent(out) :: ptr
    type(dtype), pointer :: t
    allocate(t)
    ptr = c_loc(t)
  end subroutine
  
  subroutine aqueous_dealloc_aqueoussolution(ptr) bind(c)
    use aqueous, only: dtype => aqueoussolution
    type(c_ptr), intent(in) :: ptr
    type(dtype), pointer :: t
    call c_f_pointer(ptr, t)
    deallocate(t)
  end subroutine
  
  subroutine aqueous_aqueoussolution_init(ptr, species_dim, species, err_len, err) bind(c)
    use aqueous, only: dtype => aqueoussolution
    use aqueous, only: STR_LEN
    type(c_ptr), intent(in) :: ptr
    integer(c_int64_t), intent(in) :: species_dim
    character(kind=c_char), intent(in) :: species(species_dim*STR_LEN)
    integer(c_int64_t), intent(out) :: err_len
    type(c_ptr), intent(out) :: err
    
    character(len=STR_LEN) :: species_f(species_dim)
    character(len=:), allocatable, target :: err_f
    character(len=:), pointer :: err_p
    type(dtype), pointer :: t
    integer :: i, k
    
    call c_f_pointer(ptr, t)
    
    do i =1,species_dim
      k = (i-1)*STR_LEN+1
      call copy_string_ctof(species(k:k+STR_LEN), species_f(i)(1:STR_LEN))
    enddo

    call t%init(species_f, err_f)
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
  
  subroutine aqueous_aqueoussolution_equilibrate(ptr, m_dim, m, T, P, err_len, err) bind(c)
    use aqueous, only: dtype => aqueoussolution
    type(c_ptr), intent(in) :: ptr
    integer(c_int64_t), intent(in) :: m_dim
    real(c_double), intent(inout) :: m(m_dim)
    real(c_double), intent(in) :: T
    real(c_double), intent(in) :: P
    integer(c_int64_t), intent(out) :: err_len
    type(c_ptr), intent(out) :: err
    
    character(len=:), allocatable, target :: err_f
    character(len=:), pointer :: err_p
    type(dtype), pointer :: tt
    
    call c_f_pointer(ptr, tt)

    call tt%equilibrate(m, T, P, err_f)
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
  
  subroutine aqueous_aqueoussolution_xtol_get(ptr, val) bind(c)
    use aqueous, only: dtype => aqueoussolution
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    
    type(dtype), pointer :: t
    call c_f_pointer(ptr, t)
    
    val = t%xtol
  end subroutine
  
  subroutine aqueous_aqueoussolution_xtol_set(ptr, val) bind(c)
    use aqueous, only: dtype => aqueoussolution
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    
    type(dtype), pointer :: t
    call c_f_pointer(ptr, t)
    
    t%xtol = val
  end subroutine
  
  subroutine aqueous_aqueoussolution_conserv_tol_get(ptr, val) bind(c)
    use aqueous, only: dtype => aqueoussolution
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    
    type(dtype), pointer :: t
    call c_f_pointer(ptr, t)
    
    val = t%conserv_tol
  end subroutine
  
  subroutine aqueous_aqueoussolution_conserv_tol_set(ptr, val) bind(c)
    use aqueous, only: dtype => aqueoussolution
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(in) :: val
    
    type(dtype), pointer :: t
    call c_f_pointer(ptr, t)
    
    t%conserv_tol = val
  end subroutine
  
  subroutine aqueous_aqueoussolution_g_init_get(ptr, val) bind(c)
    use aqueous, only: dtype => aqueoussolution
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    
    type(dtype), pointer :: t
    call c_f_pointer(ptr, t)
    
    val = t%g_init
  end subroutine
  
  subroutine aqueous_aqueoussolution_g_opt_get(ptr, val) bind(c)
    use aqueous, only: dtype => aqueoussolution
    type(c_ptr), intent(in) :: ptr
    real(c_double), intent(out) :: val
    
    type(dtype), pointer :: t
    call c_f_pointer(ptr, t)
    
    val = t%g_opt
  end subroutine
  
  subroutine aqueous_aqueoussolution_algorithm_get(ptr, val) bind(c)
    use aqueous, only: dtype => aqueoussolution
    use aqueous, only: STR_LEN
    type(c_ptr), intent(in) :: ptr
    character(kind=c_char), intent(out) :: val(STR_LEN)
    
    type(dtype), pointer :: t
    call c_f_pointer(ptr, t)
    
    call copy_string_ftoc(t%algorithm, val)
  end subroutine
  
  subroutine aqueous_aqueoussolution_algorithm_set(ptr, val_len, val) bind(c)
    use aqueous, only: dtype => aqueoussolution
    type(c_ptr), intent(in) :: ptr
    integer(c_int64_t), intent(in) :: val_len
    character(kind=c_char), intent(in) :: val(val_len)
    
    type(dtype), pointer :: t
    call c_f_pointer(ptr, t)
    
    call copy_string_ctof(val, t%algorithm)
  end subroutine
  
  subroutine aqueous_gibbs_energy(species_len, species, T, P, err_len, err, G) bind(c)
    use aqueous, only: gibbs_energy
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
  
  subroutine aqueous_err(err_len, err_cp, err) bind(c)
    integer(c_int64_t), intent(in) :: err_len
    type(c_ptr), intent(in), value :: err_cp
    character(kind=c_char), intent(out) :: err(err_len)
    
    character(kind=c_char), pointer :: err_p(:)
    
    call c_f_pointer(err_cp, err_p, [err_len])
    err = err_p
    deallocate(err_p)
    
  end subroutine
  
  subroutine aqueous_load_spronsbl(path_len, path) bind(c)
    use aqueous, only: load_spronsbl
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