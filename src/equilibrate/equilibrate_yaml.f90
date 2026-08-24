module equilibrate_yaml
  use equilibrate_const, only: dp, s_str_len
  implicit none
  private

  public :: ReactantList, ShomatePolynomial, Nasa9Polynomial, Nasa7Polynomial, check_for_duplicates

  enum, bind(c)
  enumerator :: &
    ShomatePolynomial = 1, &
    Nasa9Polynomial = 2, &
    Nasa7Polynomial = 3
  end enum
  
  type :: ThermodynamicData
    integer :: dtype
    integer :: ntemps
    real(dp), allocatable :: temps(:)
    real(dp), allocatable :: data(:,:)
  end type

  type :: Reactant
    character(:), allocatable :: name
    integer, allocatable :: composition(:) ! (natoms)
    real(dp) :: mass
    logical :: condensate
    
    ! thermodynamics
    type(ThermodynamicData) :: thermo
    
  end type

  type :: ReactantList
    integer :: natoms
    character(s_str_len), allocatable :: atoms_names(:)
    real(dp), allocatable :: atoms_mass(:)

    integer :: nr
    type(Reactant), allocatable :: r(:) ! (nr)
  end type
  interface ReactantList
    module procedure create_ReactantList
  end interface

contains

  function create_ReactantList(filename, err) result(rl)
    use fortran_yaml_c, only: YamlFile, type_dictionary
    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: err 
    
    type(ReactantList) :: rl
    
    type(YamlFile) :: file
    
    call file%parse(filename, err)
    if (allocated(err)) return

    select type (root => file%root)
      class is (type_dictionary)
        call unpack_speciesfile(root, filename, rl, err)
      class default
        err = 'yaml file "'//filename//'" must have dictionaries at the root level.'
        return
    end select
    call file%finalize()
    
  end function
  
  subroutine unpack_speciesfile(root, filename, sp, err)
    use fortran_yaml_c_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
    type(type_dictionary), pointer :: root
    character(*), intent(in) :: filename
    type(ReactantList), intent(inout) :: sp
    character(:), allocatable, intent(out) :: err
    
    type(type_list), pointer :: tmp_list
    type(type_dictionary), pointer :: dict
    type(type_key_value_pair), pointer :: key_value_pair
    type(type_list_item), pointer :: item
    type(type_error), allocatable :: io_err
    
    character(:), allocatable :: tmp_str
    character(s_str_len) :: tmp_str1
    integer :: i, j, ind

    !!! atoms !!!
    tmp_list => root%get_list('atoms',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    sp%natoms = tmp_list%size()
    allocate(sp%atoms_names(sp%natoms))
    allocate(sp%atoms_mass(sp%natoms))
    
    j = 1
    item => tmp_list%first
    do while(associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        sp%atoms_names(j) = element%get_string("name",error = io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        sp%atoms_mass(j) = element%get_real("mass",error = io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      class default
        err = '"atoms" in "'//trim(filename)//'" must made of dictionaries.'
        return 
      end select
      j = j + 1
      item => item%next
    enddo

    ! check for duplicates
    ind = check_for_duplicates(sp%atoms_names)
    if (ind /= 0) then
      err = '"'//trim(sp%atoms_names(ind))//'" is a duplicate atom in "'//filename//'"'
      return
    endif

    !!! done with atoms !!!
    
    !!! species !!!
    tmp_list => root%get_list('species',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
    sp%nr = tmp_list%size()
    allocate(sp%r(sp%nr))
    do j = 1,sp%nr
      allocate(sp%r(j)%composition(sp%natoms))
      sp%r(j)%composition = 0
    enddo
    
    j = 1
    item => tmp_list%first
    do while(associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        ! name
        tmp_str = element%get_string("name",error = io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        sp%r(j)%name = trim(tmp_str)

        ! composition
        dict => element%get_dictionary("composition",.true.,error = io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        key_value_pair => dict%first
        do while (associated(key_value_pair))
          tmp_str1 = trim(key_value_pair%key)
          ind = findloc(sp%atoms_names,tmp_str1, 1)
          if (ind == 0) then
            err = 'The atom "'// trim(key_value_pair%key)//'" in species "'// &
                  sp%r(j)%name//'" is not in the list of atoms.'
            return
          endif
          key_value_pair =>key_value_pair%next
        enddo
        do i=1,sp%natoms
          sp%r(j)%composition(i) =  &
              dict%get_integer(sp%atoms_names(i), default = 0, error = io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        enddo
        
        ! mass
        sp%r(j)%mass = sum(sp%r(j)%composition*sp%atoms_mass)

        ! condensible?
        sp%r(j)%condensate = element%get_logical('condensate', default=.false., error=io_err)
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        
        ! thermodynamics
        call unpack_thermo(element, sp%r(j)%name, filename, sp%r(j)%thermo, err)
        if (allocated(err)) return
        
      class default
        err = '"species" in "'//trim(filename)//'" must made of dictionaries.'
        return 
      end select
      j = j + 1
      item => item%next
    enddo

    ! Check that at least one species has every atom
    block
      integer, allocatable :: composition(:)
      allocate(composition(sp%natoms))
      composition = 0
      do i = 1,sp%nr
        composition = composition + abs(sp%r(i)%composition)
      enddo
      ind = findloc(composition, 0, 1)
      if (ind /= 0) then
        err = 'Atom '//trim(sp%atoms_names(ind))//' is not in any species'
        return
      endif
    end block
    
    ! check for duplicates
    block
      character(s_str_len), allocatable :: tmp_str_list(:)
      allocate(tmp_str_list(sp%nr))
      do i = 1,sp%nr
        tmp_str_list(i) = sp%r(i)%name
      enddo
      ind = check_for_duplicates(tmp_str_list)
      if (ind /= 0) then
        err = '"'//sp%r(ind)%name//'" is a duplicate species in "'//filename//'"'
        return
      endif
    end block

    !!! done with species !!!

  end subroutine
  
  subroutine unpack_thermo(molecule, molecule_name, infile, thermo, err)
    use fortran_yaml_c_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
    class(type_dictionary), intent(in) :: molecule
    character(len=*), intent(in) :: molecule_name
    character(len=*), intent(in) :: infile
    
    type(ThermodynamicData), intent(out) :: thermo
    character(:), allocatable, intent(out) :: err
    
    type (type_error), allocatable :: io_err
    class(type_dictionary), pointer :: tmpdict
    class(type_list), pointer :: tmplist
    class(type_list_item), pointer :: item, item1
    character(len=:), allocatable :: model
    logical :: success
    
    integer :: j, k
    
    tmpdict => molecule%get_dictionary("thermo",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    ! check thermodynamic model
    model = tmpdict%get_string("model",error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (model == "Shomate") then
      thermo%dtype = ShomatePolynomial
    elseif (model == "NASA9") then
      thermo%dtype = Nasa9Polynomial
    elseif (model == "NASA7") then
      thermo%dtype = Nasa7Polynomial
    else
      err = "Thermodynamic data must be in Shomate, NASA9 or NASA7 format for "//trim(molecule_name)
      return
    endif
    
    ! get temperature ranges
    tmplist =>tmpdict%get_list("temperature-ranges",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    thermo%ntemps = tmplist%size() - 1
    if (thermo%ntemps < 1) then
      err = "Problem reading thermodynamic data for "//trim(molecule_name)
      return
    endif
    allocate(thermo%temps(thermo%ntemps + 1))
    
    j = 1
    item => tmplist%first
    do while (associated(item))
      select type (listitem => item%node)
      class is (type_scalar)
        thermo%temps(j) = listitem%to_real(-1.0_dp,success)
        if (.not. success) then
          err = "Problem reading thermodynamic data for "//trim(molecule_name)
          return
        endif
      class default
        err = "Problem reading thermodynamic data for "//trim(molecule_name)
        return
      end select
      item => item%next
      j = j + 1
    enddo
    
    ! get data
    tmplist => tmpdict%get_list("data",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    if (tmplist%size() /= thermo%ntemps) then
      err = "Problem reading thermodynamic data for "//trim(molecule_name)
      return
    endif
    
    if (thermo%dtype == ShomatePolynomial) then
      ! Shomate
      allocate(thermo%data(7,thermo%ntemps))
    elseif (thermo%dtype == Nasa9Polynomial) then
      ! NASA9
      allocate(thermo%data(9,thermo%ntemps))
    elseif (thermo%dtype == Nasa7Polynomial) then
      allocate(thermo%data(7,thermo%ntemps))
    endif
    
    k = 1
    item => tmplist%first
    do while (associated(item))
      select type (listitem => item%node)
      class is (type_list)
        
        if (listitem%size() /= size(thermo%data,1)) then
          err = "Too much or too little thermodynamic data for "//trim(molecule_name)
          return
        endif
        
        j = 1
        item1 => listitem%first
        do while (associated(item1)) 
          select type (listitem1 => item1%node)
          class is (type_scalar)

            thermo%data(j, k) = listitem1%to_real(-1.0_dp,success)
            if (.not.success) then
              err = "Problem reading thermodynamic data for "//trim(molecule_name)
              return
            endif
          class default
            err = "Problem reading thermodynamic data for "//trim(molecule_name)
            return
          end select
        item1 => item1%next
        j = j + 1
        enddo
      class default
        err = "IOError: Problem reading thermodynamic data for "//trim(molecule_name)
        return
      end select
      item => item%next
      k = k + 1
    enddo          
                            
  end subroutine

  pure function check_for_duplicates(str_list) result(ind)
    character(*), intent(in) :: str_list(:)
    integer :: ind
    integer :: i, j
    ind = 0
    do i = 1,size(str_list)-1
      do j = i+1,size(str_list)
        if (str_list(i) == str_list(j)) then
          ind = i
          exit
        endif
      enddo
    enddo
  end function

end module