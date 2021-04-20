
module photochem_types ! make a giant IO object
  implicit none
  
  private
  integer,parameter :: real_kind = kind(1.0d0)
  
  type, public :: yamldata
    integer :: nsp, nr, natoms
    character(len=8), allocatable :: atoms_names(:)
    character(len=8), allocatable :: species_names(:)
    integer, allocatable :: species_composition(:,:)
    integer, allocatable :: lowerboundcond(:)
    real(real_kind), allocatable :: surface_vdep(:)
    real(real_kind), allocatable :: surface_flux(:)
    real(real_kind), allocatable :: surface_distributed_flux(:)
    real(real_kind), allocatable :: surface_fixed_mr(:)
    integer, allocatable :: upperboundcond(:)
    real(real_kind), allocatable :: top_veff(:)
    real(real_kind), allocatable :: top_flux(:)
    real(real_kind), allocatable :: thermo_data(:,:,:)
    real(real_kind), allocatable :: thermo_temps(:,:)
    character(len=8), allocatable :: reactions_names(:,:)
    integer, allocatable :: reactions_indices(:,:)
    character(len=15), allocatable :: rxtypes(:)
    real(real_kind), allocatable :: rateparams(:,:)
  end type
  
end module


module photochem_vars ! unpack the IO object to photochem variable 
  use photochem_types, only : yamldata
  implicit none
  
  integer,parameter :: real_kind = kind(1.0d0)
  
  integer :: nsp, nr, natoms
  character(len=8), allocatable :: atoms_names(:)
  character(len=8), allocatable :: species_names(:)
  integer, allocatable :: species_composition(:,:)
  integer, allocatable :: lowerboundcond(:)
  real(real_kind), allocatable :: surface_vdep(:)
  real(real_kind), allocatable :: surface_flux(:)
  real(real_kind), allocatable :: surface_distributed_flux(:)
  real(real_kind), allocatable :: surface_fixed_mr(:)
  integer, allocatable :: upperboundcond(:)
  real(real_kind), allocatable :: top_veff(:)
  real(real_kind), allocatable :: top_flux(:)
  real(real_kind), allocatable :: thermo_data(:,:,:)
  real(real_kind), allocatable :: thermo_temps(:,:)
  character(len=8), allocatable :: reactions_names(:,:)
  integer, allocatable :: reactions_indices(:,:)
  character(len=15), allocatable :: rxtypes(:)
  real(real_kind), allocatable :: rateparams(:,:)
  
contains

end module


module photochem_io
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar
  use stringifor, only : string
  implicit none

  private 
  
  public tester
  
  integer,parameter :: real_kind = kind(1.0d0)
  
contains
  
  subroutine tester(infile) ! read yaml and make giant IO object
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    character(error_length) :: error
    class (type_node), pointer :: root
    
    root => parse(infile,unit=100,error=error)
    if (error/='') then
      write (*,*) 'PARSE ERROR: '//trim(error)
      stop
    end if
    
    select type (root)
      class is (type_dictionary)
        call parserootdict(root)
        call root%finalize()
        deallocate(root)
      class default
        print*,"yaml file must have dictionaries at root level"
        stop
    end select
  end subroutine
  
  subroutine parserootdict(mapping)
    class (type_dictionary), intent(in), pointer :: mapping
    class (type_dictionary), pointer :: settings, planet
    class (type_list), pointer :: atoms, species, reactions
    type (type_error), pointer :: config_error
    class (type_list_item), pointer :: item
    class (type_dictionary), pointer :: dict

    ! temporary work
    type(string) :: tmp
    type(string), allocatable :: tmps(:)
    character(len=8) :: outstr(5)
    integer :: outarr(5)
    integer :: i, j
    
    ! useful things to output
    integer :: nsp, nr, natoms
    character(len=8), allocatable :: atoms_names(:)
    character(len=8), allocatable :: species_names(:)
    integer, allocatable :: species_composition(:,:)
    integer, allocatable :: lowerboundcond(:)
    real(real_kind), allocatable :: surface_vdep(:)
    real(real_kind), allocatable :: surface_flux(:)
    real(real_kind), allocatable :: surface_distributed_flux(:)
    real(real_kind), allocatable :: surface_fixed_mr(:)
    integer, allocatable :: upperboundcond(:)
    real(real_kind), allocatable :: top_veff(:)
    real(real_kind), allocatable :: top_flux(:)
    real(real_kind), allocatable :: thermo_data(:,:,:)
    real(real_kind), allocatable :: thermo_temps(:,:)
    character(len=8), allocatable :: reactions_names(:,:)
    integer, allocatable :: reactions_indices(:,:)
    character(len=15), allocatable :: rxtypes(:)
    real(real_kind), allocatable :: rateparams(:,:)
    ! useful things above
    
    settings => mapping%get_dictionary('settings',.true.,error = config_error)
    planet => mapping%get_dictionary('planet',.true.,error = config_error)
    atoms => mapping%get_list('atoms',.true.,error = config_error)
    species => mapping%get_list('species',.true.,error = config_error)
    reactions => mapping%get_list('reactions',.true.,error = config_error) 
    
    ! first do atoms
    item => atoms%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_scalar)
        tmp = tmp//" "//trim(element%string)
      end select
      item => item%next
    enddo
    call tmp%split(tokens=tmps) ! list of atoms
    natoms = size(tmps)
    allocate(atoms_names(natoms))
    do i = 1,natoms
      atoms_names(i) = tmps(i)%chars()
    enddo
    ! done with atoms
    
    ! now do species
    nsp = 0 ! count number of species
    item => species%first
    do while (associated(item))
      item => item%next
      nsp = nsp + 1
    enddo
        
    allocate(species_composition(natoms,nsp))
    allocate(species_names(nsp))
    allocate(surface_vdep(nsp))
    allocate(surface_flux(nsp))
    allocate(surface_distributed_flux(nsp))
    allocate(surface_fixed_mr(nsp))
    allocate(top_veff(nsp))
    allocate(top_flux(nsp))
    allocate(thermo_data(2,7,nsp))
    allocate(thermo_temps(3,nsp))
    j = 1
    item => species%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        dict => element%get_dictionary("composition",.true.,error = config_error)
        do i=1,natoms
          species_composition(i,j) = dict%get_integer(atoms_names(i),0,error = config_error)
        enddo
        species_names(j) = trim(element%get_string("name","",error = config_error))
      class default
        print*,"Problem with species number ", j,"  in the input file"
        stop
      end select
      item => item%next
      j = j + 1
    enddo
    
    ! reactions
    nr = 0 ! count equations
    item => reactions%first
    do while (associated(item))
      item => item%next
      nr = nr + 1
    enddo
    
    allocate(reactions_names(5,nr))
    allocate(reactions_indices(5,nr))
    allocate(rateparams(6,nr))
    allocate(rxtypes(nr))
    j = 1
    item => reactions%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmp = trim(element%get_string("equation","",error = config_error))
        call parse_equation(tmp,outstr)
        reactions_names(:,j) = outstr
        call species_name2number(tmp ,outstr, species_names, species_composition, natoms,nsp, outarr)
        reactions_indices(:,j) = outarr
        call find_rateparams(element, rxtypes(j), rateparams(:,j))
      class default
        print*,"IOError: Problem with reaction number ",j," in the input file."
        stop
      end select
      item => item%next
      j = j + 1
    enddo

  end subroutine
  
  subroutine find_rateparams(reaction, rxtype, rateparam)
    class(type_dictionary), intent(in) :: reaction
    character(len=15), intent(out) :: rxtype
    real(real_kind), intent(out) :: rateparam(6)
    
    type (type_error), pointer :: config_error
    class(type_dictionary), pointer :: tmpdict
    class(type_scalar), pointer :: tmpscalar
    
    rateparam = 0.d0
    
    rxtype= reaction%get_string("type","",error = config_error)
    if (trim(rxtype) == '') rxtype = "elementary"
    
    ! get params
    if ((trim(rxtype) == 'elementary') .or. (trim(rxtype) == 'three-body')) then
      tmpdict => reaction%get_dictionary('rate-constant',.true.,error = config_error)
      rateparam(1) = tmpdict%get_real('A',0.d0,error = config_error)
      rateparam(2) = tmpdict%get_real('b',0.d0,error = config_error)
      rateparam(3) = tmpdict%get_real('Ea',0.d0,error = config_error)
    elseif (trim(rxtype) == 'falloff') then
      tmpdict => reaction%get_dictionary('low-P-rate-constant',.true.,error = config_error)
      rateparam(1) = tmpdict%get_real('A',0.d0,error = config_error)
      rateparam(2) = tmpdict%get_real('b',0.d0,error = config_error)
      rateparam(3) = tmpdict%get_real('Ea',0.d0,error = config_error)
      tmpdict => reaction%get_dictionary('high-P-rate-constant',.true.,error = config_error)
      rateparam(4) = tmpdict%get_real('A',0.d0,error = config_error)
      rateparam(5) = tmpdict%get_real('b',0.d0,error = config_error)
      rateparam(6) = tmpdict%get_real('Ea',0.d0,error = config_error)
    elseif (trim(rxtype) == 'photolysis') then
      ! nothing
    else
      print*,'IOError: reaction type ',trim(rxtype),' is not a valid reaction type.'
      stop
    endif
    
  end subroutine
  
  subroutine parse_equation(instring,outstring)
    type(string), intent(in) :: instring
    character(len=8), intent(out) :: outstring(5)
    type(string) :: string1, string2, string3, string4
    type(string), allocatable :: eq1(:), eqr(:), eqp(:)
    integer i
    string1 = instring%replace(old='+', new=' ')
    string2 = string1%replace(old='M', new=' ')
    string3 = string2%replace(old='(', new=' ')
    string1 = string3%replace(old='hv', new=' ')
    string4 = string1%replace(old=')', new=' ')
    if (index(instring%chars(), "<=>") /= 0) then
      call string4%split(eq1, sep="<=>")
    elseif (index(instring%chars(), " =>") /= 0) then
      call string4%split(eq1, sep=" =>")
    else
      ! problem
    endif
    call eq1(1)%split(eqr, sep=" ")
    call eq1(2)%split(eqp, sep=" ")
    outstring = ''
    do i=1,size(eqr)
      outstring(i) = eqr(i)%chars()
    enddo
    do i=1,size(eqp)
      outstring(2+i) = eqp(i)%chars()
    enddo
  end subroutine
  
  subroutine species_name2number(reaction, instring, species_names, species_composition, natoms, nsp, outarr)
    type(string) :: reaction
    character(len=8), intent(in) :: instring(5)
    character(len=8), intent(in) :: species_names(nsp)
    integer, intent(in) :: species_composition(natoms,nsp)
    integer, intent(in) :: nsp, natoms
    integer, intent(out) :: outarr(5)
    integer :: i, ind(1)
    integer :: reactant_atoms(natoms), product_atoms(natoms)
    reactant_atoms = 0
    product_atoms = 0
    do i = 1,5
      ind = findloc(species_names,instring(i))
      outarr(i) = ind(1)
      if ((instring(i) /= '') .and. (ind(1) == 0)) then
        print*,"IOError: ", & 
               "Species ",trim(instring(i))," in reaction ",reaction, &
               "is not in the list of species."
        stop
      endif
    enddo
    do i=1,2
      if (outarr(i) /= 0) then
        reactant_atoms = reactant_atoms + species_composition(:,outarr(i))
      endif
    enddo
    do i=3,5
      if (outarr(i) /= 0) then
        product_atoms = product_atoms + species_composition(:,outarr(i))
      endif
    enddo
    if (.not. all(reactant_atoms == product_atoms)) then
      print*,"IOError: ", & 
             "Bad mass balance in reaction ",reaction
      stop
    endif
  end subroutine
  
  subroutine handleerror(message)
    character(len=*), intent(in) :: message
    print*,message
    stop
  end subroutine

end module

program main
  use photochem_io
  implicit none

  call tester("../zahnle.yaml")
end program





