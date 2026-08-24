module equilibrate_c_api
  use iso_c_binding
  implicit none

  integer, parameter :: err_len = 1024

contains

  !~~ Allocator and destroyer ~~!

  function allocate_chemequianalysis() result(ptr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr) :: ptr
    type(ChemEquiAnalysis), pointer :: cea
    allocate(cea)
    ptr = c_loc(cea)
  end function

  subroutine deallocate_chemequianalysis(ptr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    deallocate(cea)
  end subroutine

  !~~ Subroutine wrappers ~~!

  subroutine chemequianalysis_create_wrapper(ptr, thermofile, &
                                             atoms_present, atoms_dim, atoms, &
                                             species_present, species_dim, species, err) bind(c)
    use equilibrate, only: ChemEquiAnalysis, s_str_len
    type(c_ptr), value, intent(in) :: ptr
    character(kind=c_char), intent(in) :: thermofile(*)
    logical(c_bool), intent(in) :: atoms_present
    integer(c_int), intent(in) :: atoms_dim
    character(kind=c_char), intent(in) :: atoms(atoms_dim*s_str_len+1)
    logical(c_bool), intent(in) :: species_present
    integer(c_int), intent(in) :: species_dim
    character(kind=c_char), intent(in) :: species(species_dim*s_str_len+1)
    character(kind=c_char), intent(out) :: err(err_len+1)

    character(len=:), allocatable :: thermofile_f
    character(s_str_len), allocatable :: atoms_f(:), species_f(:)
    character(:), allocatable :: err_f
    integer :: i, j, k
    type(ChemEquiAnalysis), pointer :: cea

    call c_f_pointer(ptr, cea)

    allocate(character(len=len_cstring(thermofile))::thermofile_f)
    call copy_string_ctof(thermofile, thermofile_f)

    allocate(atoms_f(atoms_dim))
    do i = 1,atoms_dim
      do j = 1,s_str_len
        k = j + (i - 1) * s_str_len
        atoms_f(i)(j:j) = atoms(k)
      enddo
    enddo

    allocate(species_f(species_dim))
    do i = 1,species_dim
      do j = 1,s_str_len
        k = j + (i - 1) * s_str_len
        species_f(i)(j:j) = species(k)
      enddo
    enddo

    if (atoms_present .and. species_present) then
      cea = ChemEquiAnalysis(thermofile_f, atoms=atoms_f, species=species_f, err=err_f)
    elseif (atoms_present .and. .not.species_present) then
      cea = ChemEquiAnalysis(thermofile_f, atoms=atoms_f, err=err_f)
    elseif (.not.atoms_present .and. species_present) then
      cea = ChemEquiAnalysis(thermofile_f, species=species_f, err=err_f)
    elseif (.not.atoms_present .and. .not.species_present) then
      cea = ChemEquiAnalysis(thermofile_f, err=err_f)
    endif

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif

  end subroutine

  subroutine chemequianalysis_solve_wrapper(ptr, P, T, &
                                            molfracs_atoms_present, molfracs_atoms_dim, molfracs_atoms, &
                                            molfracs_species_present, molfracs_species_dim, molfracs_species, &
                                            converged, err) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(in) :: P
    real(c_double), intent(in) :: T
    logical(c_bool), intent(in) :: molfracs_atoms_present
    integer(c_int), intent(in) :: molfracs_atoms_dim
    real(c_double), intent(in) :: molfracs_atoms(molfracs_atoms_dim)
    logical(c_bool), intent(in) :: molfracs_species_present
    integer(c_int), intent(in) :: molfracs_species_dim
    real(c_double), intent(in) :: molfracs_species(molfracs_species_dim)
    logical(c_bool), intent(out) :: converged
    character(kind=c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(ChemEquiAnalysis), pointer :: cea

    call c_f_pointer(ptr, cea)

    if (molfracs_atoms_present .and. molfracs_species_present) then
      converged = cea%solve(P, T, molfracs_atoms=molfracs_atoms, molfracs_species=molfracs_species, err=err_f)
    elseif (molfracs_atoms_present .and. .not.molfracs_species_present) then
      converged = cea%solve(P, T, molfracs_atoms=molfracs_atoms, err=err_f)
    elseif (.not.molfracs_atoms_present .and. molfracs_species_present) then
      converged = cea%solve(P, T, molfracs_species=molfracs_species, err=err_f)
    elseif (.not.molfracs_atoms_present .and. .not.molfracs_species_present) then
      converged = cea%solve(P, T, err=err_f)
    endif

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif

  end subroutine

  subroutine chemequianalysis_solve_metallicity_wrapper(ptr, P, T, &
                                                        metallicity, CtoO_present, CtoO, &
                                                        converged, err) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(in) :: P
    real(c_double), intent(in) :: T
    real(c_double), intent(in) :: metallicity
    logical(c_bool), intent(in) :: CtoO_present
    real(c_double), intent(in) :: CtoO
    logical(c_bool), intent(out) :: converged
    character(kind=c_char), intent(out) :: err(err_len+1)

    character(:), allocatable :: err_f
    type(ChemEquiAnalysis), pointer :: cea

    call c_f_pointer(ptr, cea)

    if (CtoO_present) then
      converged = cea%solve_metallicity(P, T, metallicity, CtoO=CtoO, err=err_f)
    else
      converged = cea%solve_metallicity(P, T, metallicity, err=err_f)
    endif

    err(1) = c_null_char
    if (allocated(err_f)) then
      call copy_string_ftoc(err_f, err)
    endif

  end subroutine

  !~~ Getters and setters ~~!

  subroutine chemequianalysis_atoms_names_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%atoms_names)
  end subroutine
  
  subroutine chemequianalysis_atoms_names_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis, s_str_len
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    character(kind=c_char), intent(out) :: arr(dim1*s_str_len+1)
    type(ChemEquiAnalysis), pointer :: cea
    
    integer :: i, j, k
    
    call c_f_pointer(ptr, cea)
    do i = 1,dim1
      do j = 1,s_str_len
        k = j + (i - 1) * s_str_len
        arr(k) = cea%atoms_names(i)(j:j)
      enddo
    enddo
    arr(dim1*s_str_len+1) = c_null_char
    
  end subroutine

  subroutine chemequianalysis_species_names_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%species_names)
  end subroutine
  
  subroutine chemequianalysis_species_names_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis, s_str_len
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    character(kind=c_char), intent(out) :: arr(dim1*s_str_len+1)
    type(ChemEquiAnalysis), pointer :: cea
    
    integer :: i, j, k
    
    call c_f_pointer(ptr, cea)
    do i = 1,dim1
      do j = 1,s_str_len
        k = j + (i - 1) * s_str_len
        arr(k) = cea%species_names(i)(j:j)
      enddo
    enddo
    arr(dim1*s_str_len+1) = c_null_char
    
  end subroutine

  subroutine chemequianalysis_gas_names_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%gas_names)
  end subroutine
  
  subroutine chemequianalysis_gas_names_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis, s_str_len
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    character(kind=c_char), intent(out) :: arr(dim1*s_str_len+1)
    type(ChemEquiAnalysis), pointer :: cea
    
    integer :: i, j, k
    
    call c_f_pointer(ptr, cea)
    do i = 1,dim1
      do j = 1,s_str_len
        k = j + (i - 1) * s_str_len
        arr(k) = cea%gas_names(i)(j:j)
      enddo
    enddo
    arr(dim1*s_str_len+1) = c_null_char
    
  end subroutine

  subroutine chemequianalysis_condensate_names_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%condensate_names)
  end subroutine
  
  subroutine chemequianalysis_condensate_names_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis, s_str_len
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    character(kind=c_char), intent(out) :: arr(dim1*s_str_len+1)
    type(ChemEquiAnalysis), pointer :: cea
    
    integer :: i, j, k
    
    call c_f_pointer(ptr, cea)
    do i = 1,dim1
      do j = 1,s_str_len
        k = j + (i - 1) * s_str_len
        arr(k) = cea%condensate_names(i)(j:j)
      enddo
    enddo
    arr(dim1*s_str_len+1) = c_null_char
    
  end subroutine

  subroutine chemequianalysis_molfracs_atoms_sun_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%molfracs_atoms_sun)
  end subroutine
  
  subroutine chemequianalysis_molfracs_atoms_sun_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    arr = cea%molfracs_atoms_sun
  end subroutine

  subroutine chemequianalysis_molfracs_atoms_sun_set(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    cea%molfracs_atoms_sun = arr
  end subroutine

  subroutine chemequianalysis_molfracs_atoms_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%molfracs_atoms)
  end subroutine
  
  subroutine chemequianalysis_molfracs_atoms_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    arr = cea%molfracs_atoms
  end subroutine

  subroutine chemequianalysis_molfracs_species_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%molfracs_species)
  end subroutine
  
  subroutine chemequianalysis_molfracs_species_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    arr = cea%molfracs_species
  end subroutine

  subroutine chemequianalysis_massfracs_species_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%massfracs_species)
  end subroutine
  
  subroutine chemequianalysis_massfracs_species_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    arr = cea%massfracs_species
  end subroutine

  subroutine chemequianalysis_molfracs_atoms_gas_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%molfracs_atoms_gas)
  end subroutine
  
  subroutine chemequianalysis_molfracs_atoms_gas_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    arr = cea%molfracs_atoms_gas
  end subroutine

  subroutine chemequianalysis_molfracs_species_gas_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%molfracs_species_gas)
  end subroutine
  
  subroutine chemequianalysis_molfracs_species_gas_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    arr = cea%molfracs_species_gas
  end subroutine

  subroutine chemequianalysis_molfracs_atoms_condensate_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%molfracs_atoms_condensate)
  end subroutine
  
  subroutine chemequianalysis_molfracs_atoms_condensate_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    arr = cea%molfracs_atoms_condensate
  end subroutine

  subroutine chemequianalysis_molfracs_species_condensate_get_size(ptr, dim1) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(out) :: dim1
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    dim1 = size(cea%molfracs_species_condensate)
  end subroutine
  
  subroutine chemequianalysis_molfracs_species_condensate_get(ptr, dim1, arr) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), intent(in) :: dim1
    real(c_double), intent(out) :: arr(dim1)
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    arr = cea%molfracs_species_condensate
  end subroutine

  subroutine chemequianalysis_mubar_get(ptr, val) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    val = cea%mubar
  end subroutine

  subroutine chemequianalysis_verbose_get(ptr, val) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    logical(c_bool), intent(out) :: val
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    val = cea%verbose
  end subroutine
  
  subroutine chemequianalysis_verbose_set(ptr, val) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    logical(c_bool), intent(in) :: val
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    cea%verbose = val
  end subroutine

  subroutine chemequianalysis_use_prev_guess_get(ptr, val) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    logical(c_bool), intent(out) :: val
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    val = cea%use_prev_guess
  end subroutine
  
  subroutine chemequianalysis_use_prev_guess_set(ptr, val) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    logical(c_bool), intent(in) :: val
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    cea%use_prev_guess = val
  end subroutine

  subroutine chemequianalysis_mass_tol_get(ptr, val) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(out) :: val
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    val = cea%mass_tol
  end subroutine
  
  subroutine chemequianalysis_mass_tol_set(ptr, val) bind(c)
    use equilibrate, only: ChemEquiAnalysis
    type(c_ptr), value, intent(in) :: ptr
    real(c_double), intent(in) :: val
    type(ChemEquiAnalysis), pointer :: cea
    call c_f_pointer(ptr, cea)
    cea%mass_tol = val
  end subroutine

  !~~ Utilities  ~~!
    
  function len_cstring(stringc) result (length)
    ! DOES NOT include the null character terminating c string
    character(kind=c_char), intent(in) :: stringc(*)
    integer(c_int) :: length
    integer, parameter :: max_len = 10000
    integer :: j  
    j = 1
    do
      if (stringc(j)==c_null_char) then
        length = j - 1
        exit
      endif
      if (j == max_len) then
        print*,"'len_cstring' tried to determine the length of an invalid C string"
        stop 1
      endif
      j = j + 1
    end do
  end function

  subroutine copy_string_ctof(stringc,stringf)
    ! utility function to convert c string to fortran string
    character(len=*), intent(out) :: stringf
    character(c_char), intent(in) :: stringc(*)
    integer j
    stringf = ''
    char_loop: do j=1,len(stringf)
      if (stringc(j)==c_null_char) exit char_loop
      stringf(j:j) = stringc(j)
    end do char_loop
  end subroutine copy_string_ctof

  subroutine copy_string_ftoc(stringf,stringc)
    ! utility function to convert c string to fortran string
    character(len=*), intent(in) :: stringf
    character(c_char), intent(out) :: stringc(:)
    integer j, n, n1, n2
    n1 = len_trim(stringf)  
    n2 = size(stringc) - 1
    n = min(n1, n2)
    do j=1,n    
      stringc(j) = stringf(j:j)   
    end do
    stringc(n+1) = c_null_char
  end subroutine copy_string_ftoc

end module
