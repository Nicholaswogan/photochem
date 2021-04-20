
module photochem_types ! make a giant IO object
  implicit none
  
  private
  integer,parameter :: real_kind = kind(1.0d0)
  
  public PhotoChemData, PhotoSettings, PhotoPlanet
  
  type :: PhotoSettings
    real(real_kind) :: bottom_atmosphere
    real(real_kind) :: top_atmosphere 
    integer :: nz = 0
    
    real(real_kind) :: lower_wavelength
    real(real_kind) :: upper_wavelength
    integer :: nw !number of bins
  end type
  
  type :: PhotoPlanet
    real(real_kind) :: gravity
    real(real_kind) :: surface_pressure
    real(real_kind) :: planet_radius
    real(real_kind) :: surface_albedo
    real(real_kind) :: trop_alt
    logical :: lightning
    real(real_kind) :: lightning_NO_production
    logical :: rainout
    real(real_kind) :: rainout_multiplier
  end type
  
  type :: PhotoChemData
    type(PhotoSettings) :: settings
    type(PhotoPlanet) :: planet
    integer :: nsp, nr, natoms
    character(len=8), allocatable :: atoms_names(:)
    character(len=8), allocatable :: species_names(:)
    integer, allocatable :: species_composition(:,:)
    integer, allocatable :: lowerboundcond(:)
    real(real_kind), allocatable :: lower_vdep(:)
    real(real_kind), allocatable :: lower_flux(:)
    real(real_kind), allocatable :: lower_distributed_height(:)
    real(real_kind), allocatable :: lower_fixed_mr(:)
    integer, allocatable :: upperboundcond(:)
    real(real_kind), allocatable :: upper_veff(:)
    real(real_kind), allocatable :: upper_flux(:)
    real(real_kind), allocatable :: thermo_data(:,:,:)
    real(real_kind), allocatable :: thermo_temps(:,:)
    character(len=8), allocatable :: reactions_names(:,:)
    integer, allocatable :: reactions_indices(:,:)
    character(len=15), allocatable :: rxtypes(:)
    real(real_kind), allocatable :: rateparams(:,:)
  end type
  
end module


module photochem_vars ! unpack the IO object to plain variables 
  use photochem_types, only : PhotoChemData
  implicit none
  integer,parameter :: real_kind = kind(1.0d0)
end module


module photochem_io
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar
  use stringifor, only : string
  use photochem_types, only : PhotoChemData, PhotoSettings, PhotoPlanet
  implicit none

  private 
  
  public tester
  
  integer,parameter :: real_kind = kind(1.0d0)
  
contains
  
  subroutine tester(infile,photodata) ! read yaml and make giant IO object
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    type(PhotoChemData), intent(out) :: photodata
    
    character(error_length) :: error
    class (type_node), pointer :: root
    
    root => parse(infile,unit=100,error=error)
    if (error/='') then
      write (*,*) 'PARSE ERROR: '//trim(error)
      stop
    end if
    
    select type (root)
      class is (type_dictionary)
        call parserootdict(root,infile,photodata)
        call root%finalize()
        deallocate(root)
      class default
        print*,"yaml file must have dictionaries at root level"
        stop
    end select
  end subroutine
  
  subroutine parserootdict(mapping, infile,photodata)
    class (type_dictionary), intent(in), pointer :: mapping
    character(len=*), intent(in) :: infile
    type(PhotoChemData), intent(out) :: photodata
    
    
    class (type_dictionary), pointer :: settings, planet
    class (type_list), pointer :: atoms, species, reactions
    type (type_error), pointer :: config_error
    class (type_list_item), pointer :: item
    class (type_dictionary), pointer :: dict

    ! temporary work variables
    type(string) :: tmp
    type(string), allocatable :: tmps(:)
    character(len=8) :: outstr(5)
    integer :: outarr(5)
    integer :: i, j
    
    settings => mapping%get_dictionary('settings',.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    planet => mapping%get_dictionary('planet',.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    atoms => mapping%get_list('atoms',.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    species => mapping%get_list('species',.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    reactions => mapping%get_list('reactions',.true.,error = config_error) 
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    ! settings
    call get_settings(settings, infile, photodata%settings)
    
    ! planet
    call get_planet(planet, infile, photodata%planet)
    
    ! atmoms
    item => atoms%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_scalar)
        tmp = tmp//" "//trim(element%string)
      class default
        print*,"IOError: Problem reading in atoms."
        stop
      end select
      item => item%next
    enddo
    call tmp%split(tokens=tmps) ! list of atoms
    photodata%natoms = size(tmps)
    allocate(photodata%atoms_names(photodata%natoms))
    do i = 1,photodata%natoms
      photodata%atoms_names(i) = tmps(i)%chars()
    enddo
    ! done with atoms
    
    ! now do species
    photodata%nsp = 0 ! count number of species
    item => species%first
    do while (associated(item))
      item => item%next
      photodata%nsp = photodata%nsp + 1
    enddo
        
    allocate(photodata%species_composition(photodata%natoms,photodata%nsp))
    allocate(photodata%species_names(photodata%nsp))
    allocate(photodata%lowerboundcond(photodata%nsp))
    allocate(photodata%lower_vdep(photodata%nsp))
    allocate(photodata%lower_flux(photodata%nsp))
    allocate(photodata%lower_distributed_height(photodata%nsp))
    allocate(photodata%lower_fixed_mr(photodata%nsp))
    allocate(photodata%upperboundcond(photodata%nsp))
    allocate(photodata%upper_veff(photodata%nsp))
    allocate(photodata%upper_flux(photodata%nsp))
    allocate(photodata%thermo_data(7,2,photodata%nsp))
    allocate(photodata%thermo_temps(3,photodata%nsp))
    j = 1
    item => species%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        dict => element%get_dictionary("composition",.true.,error = config_error)  ! get composition
        if (associated(config_error)) call handleerror(config_error%message,infile)
        do i=1,photodata%natoms
          photodata%species_composition(i,j) = dict%get_integer(photodata%atoms_names(i),0,error = config_error) ! no error possible.
        enddo
        photodata%species_names(j) = trim(element%get_string("name",error = config_error)) ! get name
        if (associated(config_error)) call handleerror(config_error%message,infile)
        call get_boundaryconds(element,photodata%species_names(j), infile, &
                               photodata%lowerboundcond(j), photodata%lower_vdep(j), &
                               photodata%lower_flux(j), photodata%lower_distributed_height(j), &
                               photodata%lower_fixed_mr(j), &
                               photodata%upperboundcond(j), photodata%upper_veff(j), photodata%upper_flux(j))! get boundary conditions
        
        call get_thermodata(element,photodata%species_names(j), infile,photodata%thermo_temps(:,j),photodata%thermo_data(:,:,j)) ! get thermodynamic data
        
      class default
        print*,"IOError: Problem with species number ", j,"  in the input file"
        stop
      end select
      item => item%next
      j = j + 1
    enddo
    
    ! reactions
    photodata%nr = 0 ! count equations
    item => reactions%first
    do while (associated(item))
      item => item%next
      photodata%nr = photodata%nr + 1
    enddo
    
    allocate(photodata%reactions_names(5,photodata%nr))
    allocate(photodata%reactions_indices(5,photodata%nr))
    allocate(photodata%rateparams(6,photodata%nr))
    allocate(photodata%rxtypes(photodata%nr))
    j = 1
    item => reactions%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmp = trim(element%get_string("equation","",error = config_error))
        call parse_equation(tmp,outstr)
        photodata%reactions_names(:,j) = outstr
        call species_name2number(tmp ,outstr, photodata%species_names, &
                                 photodata%species_composition, photodata%natoms, &
                                 photodata%nsp, outarr)
        photodata%reactions_indices(:,j) = outarr
        call get_rateparams(element, photodata%rxtypes(j), photodata%rateparams(:,j))
      class default
        print*,"IOError: Problem with reaction number ",j," in the input file."
        stop
      end select
      item => item%next
      j = j + 1
    enddo

  end subroutine
  
  
  subroutine get_settings(filesettings, infile, outsettings)
    class(type_dictionary), intent(in) :: filesettings
    character(len=*), intent(in) :: infile
    type(PhotoSettings), intent(out) :: outsettings
    
    type (type_error), pointer :: config_error
    class(type_dictionary), pointer :: tmpdict
    
    tmpdict => filesettings%get_dictionary("atmosphere-grid",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    outsettings%bottom_atmosphere = tmpdict%get_real("bottom",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outsettings%top_atmosphere = tmpdict%get_real("top",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outsettings%nz = tmpdict%get_integer("number-of-layers",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    tmpdict => filesettings%get_dictionary("photo-grid",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    outsettings%lower_wavelength = tmpdict%get_real("lower-wavelength",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outsettings%upper_wavelength = tmpdict%get_real("upper-wavelength",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outsettings%nw = tmpdict%get_integer("number-of-bins",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
  
  end subroutine
  
  subroutine get_planet(fileplanet, infile, outplanet)
    class(type_dictionary), intent(in) :: fileplanet
    character(len=*), intent(in) :: infile
    type(PhotoPlanet), intent(out) :: outplanet
    
    type (type_error), pointer :: config_error
    class(type_dictionary), pointer :: tmpdict
    
    outplanet%gravity = fileplanet%get_real("gravity",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outplanet%surface_pressure = fileplanet%get_real("surface-pressure",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    outplanet%planet_radius = fileplanet%get_real("planet-radius",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)  
    outplanet%surface_albedo = fileplanet%get_real("surface-albedo",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)  
    outplanet%trop_alt = fileplanet%get_real("tropopause-altitude",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    tmpdict => fileplanet%get_dictionary("lightning",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    outplanet%lightning = tmpdict%get_logical("on-off",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    if (.not. outplanet%lightning) then
      outplanet%lightning_NO_production = 0.d0
    else
      outplanet%lightning_NO_production = tmpdict%get_real("NO-production",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    endif
    tmpdict => fileplanet%get_dictionary("rainout",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    outplanet%rainout= tmpdict%get_logical("on-off",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    if (.not.outplanet%rainout) then
      outplanet%rainout_multiplier = 1.d0
    else
      outplanet%rainout_multiplier = tmpdict%get_real("multiplier",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    endif
  end subroutine
  
  subroutine get_boundaryconds(molecule, molecule_name, infile, &
                               lowercond, Lvdep, Lflux, LdistH, Lmr, &
                               uppercond, Uveff, Uflux)
    class(type_dictionary), intent(in) :: molecule
    character(len=*), intent(in) :: molecule_name
    character(len=*), intent(in) :: infile
    
    integer, intent(out) :: lowercond
    real(real_kind), intent(out) :: Lvdep, Lflux, LdistH, Lmr
    integer, intent(out) :: uppercond
    real(real_kind), intent(out) :: Uveff, Uflux
    
    type (type_error), pointer :: config_error
    class(type_dictionary), pointer :: tmpdict
    character(len=:), allocatable :: bctype
    
    tmpdict => molecule%get_dictionary("lower-boundary",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    bctype = tmpdict%get_string("type",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    ! constant deposition velocity = vdep
    ! constant mixing-ratio = mixing-ratio
    ! constant flux = flux
    ! deposition velocity + distributed flux = vdep + flux + distributed-height
    if (bctype == "constant deposition velocity") then
      lowercond = 0
      Lvdep = tmpdict%get_real("vdep",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      Lflux = 0.d0
      LdistH = 0.d0
      Lmr = 0.d0 
    elseif (bctype == "constant mixing-ratio") then
      lowercond = 1
      Lvdep = 0.d0
      Lflux = 0.d0
      LdistH = 0.d0
      Lmr = tmpdict%get_real("mixing-ratio",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    elseif (bctype == "constant flux") then
      lowercond = 2
      Lvdep = 0.d0
      Lflux = tmpdict%get_real("flux",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      LdistH = 0.d0
      Lmr = 0.d0 
    elseif (bctype == "deposition velocity + distributed flux") then
      lowercond = 3
      Lvdep = tmpdict%get_real("vdep",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      Lflux = tmpdict%get_real("flux",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      LdistH = tmpdict%get_real("distributed-height",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      Lmr = 0.d0 
    else
      print*,'IOError: "',trim(bctype),'" is not a valid lower boundary condition for ',trim(molecule_name)
      stop
    endif
    
    tmpdict => molecule%get_dictionary("upper-boundary",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    bctype = tmpdict%get_string("type",error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    ! constant effusion velocity = veff
    ! constant flux = flux
    if (bctype == "constant effusion velocity") then
      uppercond = 0
      Uveff = tmpdict%get_real("veff",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
      Uflux = 0.d0
    elseif (bctype == "constant flux") then
      uppercond = 2
      Uveff = 0.d0
      Uflux = tmpdict%get_real("flux",error = config_error)
      if (associated(config_error)) call handleerror(config_error%message,infile)
    else
      print*,'IOError: "',trim(bctype),'" is not a valid upper boundary condition for ',trim(molecule_name)
      stop
    endif
    
  end subroutine
    
  
  subroutine get_thermodata(molecule, molecule_name, infile, thermo_temps_entry, thermo_data_entry)
    class(type_dictionary), intent(in) :: molecule
    character(len=*), intent(in) :: molecule_name
    character(len=*), intent(in) :: infile
    
    real(real_kind), intent(out) :: thermo_temps_entry(3)
    real(real_kind), intent(out) :: thermo_data_entry(7,2)
    
    type (type_error), pointer :: config_error
    class(type_dictionary), pointer :: tmpdict, tmpdict1
    class(type_list), pointer :: tmplist
    class(type_list_item), pointer :: item
    logical :: success
    
    integer :: j, i
    
    thermo_temps_entry = -1.d0
    thermo_data_entry = -1.d0
    
    tmpdict => molecule%get_dictionary("thermo",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    ! check thermodynamic model
    if (tmpdict%get_string("model",error = config_error) /= "Shomate") then
      print*,"IOError: Thermodynamic data must be in Shomate format for ",trim(molecule_name)
      stop
    endif
    if (associated(config_error)) call handleerror(config_error%message,infile)
    
    ! get temperature ranges
    tmplist =>tmpdict%get_list("temperature-ranges",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    j = 1
    item => tmplist%first
    do while (associated(item))
      select type (listitem => item%node)
      class is (type_scalar)
        if (j > 3) then
          print*,"IOError: Too many temperature ranges for ",trim(molecule_name)
          stop
        endif
        thermo_temps_entry(j) = listitem%to_real(-1.d0,success)
        if (.not. success) then
          print*,"IOError: Problem reading thermodynamic data for  ",trim(molecule_name)
          stop
        endif
      class default
        print*,"IOError: Problem reading thermodynamic data for ",trim(molecule_name)
        stop
      end select
      item => item%next
      j = j + 1
    enddo
    
    ! check amount of thermodynamic data
    if (thermo_temps_entry(3) == -1.d0) then
      i = 1
    else
      i = 2
    endif
    
    ! get data
    tmpdict1 => tmpdict%get_dictionary("data",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    tmplist =>tmpdict1%get_list("P1",.true.,error = config_error)
    if (associated(config_error)) call handleerror(config_error%message,infile)
    j = 1
    item => tmplist%first
    do while (associated(item)) 
      select type (listitem => item%node)
      class is (type_scalar)
        if (j > 7) then
          print*,"IOError: Too much thermodynamic data for ",trim(molecule_name)
          stop
        endif
        thermo_data_entry(j, 1) = listitem%to_real(-1.d0,success)
        if (.not. success) then
          print*,"IOError: Problem reading thermodynamic data for  ",trim(molecule_name)
          stop
        endif
      class default
        print*,"IOError: Problem reading thermodynamic data for ",trim(molecule_name)
        stop
      end select
      item => item%next
      j = j + 1
    enddo

    tmplist =>tmpdict1%get_list("P2",.false.,error = config_error)
    if ((.not.associated(tmplist)) .and. (i == 1)) then
      ! nothing happens
    elseif ((.not.associated(tmplist)) .and. (i == 2)) then
      print*,'IOError: More temperature ranges than thermodynamic data for ',trim(molecule_name)
      stop
    else
      j = 1
      item => tmplist%first
      do while (associated(item)) 
        select type (listitem => item%node)
        class is (type_scalar)
          if (j > 7) then
            print*,"IOError: Too much thermodynamic data for ",trim(molecule_name)
            stop
          endif
          thermo_data_entry(j, 2) = listitem%to_real(-1.d0,success)
          if (.not. success) then
            print*,"IOError: Problem reading thermodynamic data for  ",trim(molecule_name)
            stop
          endif
        class default
          print*,"IOError: Problem reading thermodynamic data for ",trim(molecule_name)
          stop
        end select
        item => item%next
        j = j + 1
      enddo
    endif
    
  end subroutine
  
  subroutine get_rateparams(reaction, rxtype, rateparam)
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
             "Bad mass balance in reaction",reaction
      print*,"You could have messed up how many atoms one of the species has."
      stop
    endif
  end subroutine
  
  subroutine handleerror(message,infile)
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: infile
    print*,"IOError: ",trim(infile),trim(message)
    stop
  end subroutine

end module

program main
  use photochem_io, only: tester
  use photochem_types, only: PhotoChemData
  implicit none
  type(PhotoChemData) :: photodata

  call tester("../zahnle.yaml", photodata)
  
  print*,photodata%settings%top_atmosphere
  
end program





