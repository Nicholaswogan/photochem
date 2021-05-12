
module photochem_io
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  use stringifor, only : string
  use photochem_types, only : PhotoMechanism, PhotoSettings, PhotoRadTran
  implicit none
  private 
  integer,parameter :: real_kind = kind(1.0d0)
  integer, parameter :: err_len = 1000
  integer, parameter :: str_len = 1000

  public get_photomech, get_photorad, get_photoset, reaction_string
    
contains
  
  subroutine get_photomech(infile, photomech, err) 
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    type(PhotoMechanism), intent(out) :: photomech
    character(len=err_len), intent(out) :: err
    
    character(error_length) :: error
    class (type_node), pointer :: root
    err = ''
    
    ! parse yaml file
    root => parse(infile,unit=100,error=error)
    if (error /= '') then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        call get_rxmechanism(root, infile, photomech, err)
        if (len_trim(err) > 0) return
        call root%finalize()
        deallocate(root)
      class default
        err = "yaml file must have dictionaries at root level"
        return
    end select
     
  end subroutine
  
  subroutine get_rxmechanism(mapping, infile, photomech, err)
    class (type_dictionary), intent(in), pointer :: mapping
    character(len=*), intent(in) :: infile
    type(PhotoMechanism), intent(out) :: photomech
    character(len=err_len), intent(out) :: err
    
    class (type_dictionary), pointer :: atoms
    class (type_list), pointer :: species, reactions
    type (type_error), pointer :: io_err
    class (type_list_item), pointer :: item
    class (type_dictionary), pointer :: dict
    class (type_key_value_pair), pointer :: key_value_pair

    ! temporary work variables
    type(string) :: tmp
    type(string), allocatable :: eqr(:), eqp(:)
    integer :: i, j, k, kk, l, ind(1)
    logical :: reverse
    err = ''

    atoms => mapping%get_dictionary('atoms',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    species => mapping%get_list('species',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    reactions => mapping%get_list('reactions',.true.,error = io_err) 
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    ! should i reverse reactions?
    photomech%reverse = mapping%get_logical('reverse-reactions',.true.,error = io_err)
    
    !!! atoms !!
    photomech%natoms = 0
    key_value_pair => atoms%first
    do while (associated(key_value_pair))
      photomech%natoms = photomech%natoms +1
      key_value_pair => key_value_pair%next
    enddo
    allocate(photomech%atoms_names(photomech%natoms))
    allocate(photomech%atoms_mass(photomech%natoms))
    j = 1
    key_value_pair => atoms%first
    do while (associated(key_value_pair))
      photomech%atoms_names(j) = trim(key_value_pair%key)
      photomech%atoms_mass(j) = atoms%get_real(trim(key_value_pair%key),error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      key_value_pair => key_value_pair%next
      j = j + 1
    enddo
    !!! done with atoms !!!
    
    !!! species !!!
    photomech%nsp = 0 ! count number of species
    item => species%first
    do while (associated(item))
      item => item%next
      photomech%nsp = photomech%nsp + 1
    enddo
    
    allocate(photomech%species_mass(photomech%nsp))
    allocate(photomech%species_composition(photomech%natoms,photomech%nsp+2))
    photomech%species_composition = 0
    allocate(photomech%species_names(photomech%nsp+2))
    photomech%species_names(photomech%nsp+1) = "hv" ! always add these guys
    photomech%species_names(photomech%nsp+2) = "M"
    if (photomech%reverse) then
      allocate(photomech%thermo_data(7,2,photomech%nsp))
      allocate(photomech%thermo_temps(3,photomech%nsp))
    endif
    j = 1
    item => species%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        dict => element%get_dictionary("composition",.true.,error = io_err)  ! get composition
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        key_value_pair => dict%first ! dont allow unspecified atoms
        do while (associated(key_value_pair))
          ind = findloc(photomech%atoms_names,trim(key_value_pair%key))
          if (ind(1) == 0) then
            err = 'IOError: The atom "'// trim(key_value_pair%key)// '" is not in the list of atoms.'
            return
          endif
          key_value_pair =>key_value_pair%next
        enddo
        
        do i=1,photomech%natoms
          photomech%species_composition(i,j) =  &
              dict%get_integer(photomech%atoms_names(i),0,error = io_err) ! no error possible.
        enddo
        photomech%species_mass(j) = sum(photomech%species_composition(:,j) * photomech%atoms_mass)
        photomech%species_names(j) = trim(element%get_string("name",error = io_err)) ! get name
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        if (photomech%reverse) then
          call get_thermodata(element,photomech%species_names(j), infile,photomech%thermo_temps(:,j), &
                              photomech%thermo_data(:,:,j), err) ! get thermodynamic data
          if (len_trim(err) > 0) return
        endif
      class default
        err = "IOError: Problem with species number "//char(j)//"  in the input file"
        return
      end select
      item => item%next
      j = j + 1
    enddo
    !!! done with species !!!
    
    !!! reactions !!!
    photomech%nrF = 0 ! count forward reactions
    item => reactions%first
    do while (associated(item))
      item => item%next
      photomech%nrF = photomech%nrF + 1
    enddo
    
    allocate(photomech%rateparams(10,photomech%nrF))
    allocate(photomech%rxtypes(photomech%nrF))
    allocate(photomech%falloff_type(photomech%nrF))
    allocate(photomech%num_efficient(photomech%nrF))
    photomech%falloff_type = - huge(1)
    ! determine which reactions to reverse. 
    ! Determine maximum number of reactants, and productants
    photomech%max_num_reactants = 1
    photomech%max_num_products = 1
    photomech%nrR = 0
    photomech%kj = 0
    j = 1
    item => reactions%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmp = trim(element%get_string("equation",error = io_err))
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        call get_rateparams(element, infile, photomech%rxtypes(j), &
                          photomech%falloff_type(j), photomech%rateparams(:,j), err)
        if (len_trim(err) > 0) return
        
        call parse_reaction(tmp, reverse, eqr, eqp, err)
        if (len_trim(err) > 0) return
        if (reverse) then
          if (.not.photomech%reverse) then
            err = 'IOError: reaction file '//trim(infile)//' contains reverse reaction '//tmp// &
                  ', which is incompatible with "reverse-reactions: false"'
            return
          endif
          photomech%nrR = photomech%nrR + 1
          if (size(eqr) > photomech%max_num_products) photomech%max_num_products = size(eqr)
          if (size(eqp) > photomech%max_num_reactants) photomech%max_num_reactants = size(eqp)
        endif
        if (size(eqr) > photomech%max_num_reactants) photomech%max_num_reactants = size(eqr)
        if (size(eqp) > photomech%max_num_products) photomech%max_num_products = size(eqp)
        
        call count_efficiencies(element, photomech%num_efficient(j))
        
        if (photomech%rxtypes(j) == 0) then ! if photolysis reaction
          photomech%kj = photomech%kj + 1
        endif
        call compare_rxtype_string(tmp, eqr, eqp, photomech%rxtypes(j),err)
        if (len_trim(err) > 0) return
      class default
        err = "IOError: Problem with reaction number "//char(j)//" in the input file."
        return
      end select
      item => item%next
      j = j+1
    enddo
    photomech%nrT = photomech%nrR + photomech%nrF

    ! allocate stuff and loop through reactions again
    allocate(photomech%nreactants(photomech%nrT))
    allocate(photomech%nproducts(photomech%nrT))
    allocate(photomech%reactants_sp_inds(photomech%max_num_reactants,photomech%nrT))
    allocate(photomech%products_sp_inds(photomech%max_num_products,photomech%nrT))
    if (photomech%reverse) then
      allocate(photomech%reverse_info(photomech%nrT))
    endif
    photomech%reverse_info = 0 ! initialize
    allocate(photomech%reactants_names(photomech%max_num_reactants,photomech%nrF))
    allocate(photomech%products_names(photomech%max_num_products,photomech%nrF))
    ! efficiency stuff
    allocate(photomech%efficiencies(maxval(photomech%num_efficient),photomech%nrF))
    allocate(photomech%eff_sp_inds(maxval(photomech%num_efficient),photomech%nrF))
    allocate(photomech%def_eff(photomech%nrF))
    photomech%efficiencies = -huge(1.d0) ! so everything blows up if we make a mistake
    photomech%eff_sp_inds = -huge(0)
    photomech%def_eff = 1.d0 ! default is 1
    
    j = 1
    k = 1
    item => reactions%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmp = trim(element%get_string("equation",error = io_err))
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        call get_efficient(element, j, infile, photomech, err)
        if (len_trim(err)>0) return
        
        call get_reaction_chars(tmp, photomech%max_num_reactants, photomech%max_num_products, &
                            photomech%nreactants(j), photomech%nproducts(j), &
                            photomech%reactants_names(:,j), photomech%products_names(:,j), reverse, err)
        if (len_trim(err)>0) return
                            
        call get_reaction_sp_nums(tmp, photomech%max_num_reactants, photomech%max_num_products, &
                                 photomech%reactants_names(:,j), photomech%products_names(:,j), &
                                 photomech%species_names, photomech%species_composition, &
                                 photomech%natoms, photomech%nsp, &
                                 photomech%reactants_sp_inds(:,j), photomech%products_sp_inds(:,j), err)
        if (len_trim(err)>0) return
        if (reverse) then
          ! reaction has a reverse
          i = photomech%nrF + k
          photomech%reverse_info(j) = i ! the reaction number of reversed reaction
          photomech%reverse_info(i) = j ! the reaction number of the forward reaction
          photomech%nreactants(i) = photomech%nproducts(j)
          photomech%nproducts(i) = photomech%nreactants(j)
          photomech%reactants_sp_inds(:,i) = photomech%products_sp_inds(:,j)
          photomech%products_sp_inds(:,i) = photomech%reactants_sp_inds(:,j)
          k = k + 1
        endif
      class default
        err = "IOError: Problem with reaction number "//char(j)//" in the input file."
        return
      end select
      item => item%next
      j = j + 1
    enddo
    
    ! nump, numl, iprod, iloss
    ! first find nump and numl then allocate iprod and iloss
    allocate(photomech%nump(photomech%nsp))
    allocate(photomech%numl(photomech%nsp))
    photomech%numl = 0
    photomech%nump = 0
    do j = 1,photomech%nrT
      k = photomech%nreactants(j)
      do i = 1,k
        kk = photomech%reactants_sp_inds(i,j)
        if ((kk /= photomech%nsp+1) .and. (kk /= photomech%nsp+2)) then
          photomech%numl(kk) = photomech%numl(kk) + 1
        endif
      enddo
      k = photomech%nproducts(j)
      do i = 1,k
        kk = photomech%products_sp_inds(i,j)
        if ((kk /= photomech%nsp+1) .and. (kk /= photomech%nsp+2)) then
          photomech%nump(kk) = photomech%nump(kk) + 1
        endif
      enddo
    enddo
    allocate(photomech%iprod(maxval(photomech%nump),photomech%nsp))
    allocate(photomech%iloss(maxval(photomech%numl),photomech%nsp))
    photomech%iprod = 0
    photomech%iloss = 0
    photomech%numl = 0
    photomech%nump = 0
    ! loop again and get iprod and iloss
    do j = 1,photomech%nrT
      k = photomech%nreactants(j)
      do i = 1,k
        kk = photomech%reactants_sp_inds(i,j)
        if ((kk /= photomech%nsp+1) .and. (kk /= photomech%nsp+2)) then
          photomech%numl(kk) = photomech%numl(kk) + 1
          l = photomech%numl(kk)
          photomech%iloss(l,kk) = j
        endif
      enddo
      k = photomech%nproducts(j)
      do i = 1,k
        kk = photomech%products_sp_inds(i,j)
        if ((kk /= photomech%nsp+1) .and. (kk /= photomech%nsp+2)) then
          photomech%nump(kk) = photomech%nump(kk) + 1
          l = photomech%nump(kk)
          photomech%iprod(l,kk) = j
        endif
      enddo
    enddo
    
    ! photolysis
    allocate(photomech%photonums(photomech%kj))
    j = 1
    do i = 1, photomech%nrF
      if (photomech%rxtypes(i) == 0) then
        photomech%photonums(j) = i
        j = j + 1
      endif
    enddo
    
    !!! end reactions !!!
    
    call check_for_duplicates(photomech,err)
    if (len(trim(err)) > 0) return
    
  end subroutine
  
  subroutine get_efficient(reaction, rxn, infile, photomech, err)
    class(type_dictionary), intent(in) :: reaction
    integer, intent(in) :: rxn
    character(len=*), intent(in) :: infile
    type(PhotoMechanism), intent(inout) :: photomech
    character(len=*), intent(out) :: err
    
    class(type_dictionary), pointer :: tmpdict
    class (type_key_value_pair), pointer :: key_value_pair
    type (type_error), pointer :: io_err
    character(len=20) :: rxn_str
    integer :: ind(1), j
    
    tmpdict => reaction%get_dictionary("efficiencies",.false.,error = io_err)
    
    if (associated(tmpdict)) then
      j = 1
      key_value_pair => tmpdict%first
      do while (associated(key_value_pair))
        
        ind = findloc(photomech%species_names,trim(key_value_pair%key))
        if (ind(1) == 0) then
          write(rxn_str,*) rxn
          err = 'IOError: Reaction number '//trim(adjustl(rxn_str))//' has efficiencies for species that are'// &
          ' not in the list of species'
        endif
        photomech%eff_sp_inds(j,rxn) = ind(1)
        photomech%efficiencies(j,rxn) = tmpdict%get_real(trim(key_value_pair%key),error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif

        key_value_pair => key_value_pair%next
        j = j + 1
      enddo
    endif
    
    photomech%def_eff(rxn) = reaction%get_real("default-efficiency",1.d0,error = io_err)

  end subroutine
  
  
  subroutine count_efficiencies(reaction, numeff)
    class(type_dictionary), intent(in) :: reaction
    integer, intent(out) :: numeff
    
    class(type_dictionary), pointer :: tmpdict
    class (type_key_value_pair), pointer :: key_value_pair
    type (type_error), pointer :: io_err
    
    tmpdict => reaction%get_dictionary("efficiencies",.false.,error = io_err)
    
    numeff = 0
    if (associated(tmpdict)) then
      ! how many?
      key_value_pair => tmpdict%first
      do while (associated(key_value_pair))
        numeff = numeff + 1
        key_value_pair => key_value_pair%next
      enddo
    endif
    
  end subroutine
  
  
  subroutine check_for_duplicates(photomech,err)
    type(PhotoMechanism), intent(in) :: photomech
    character(len=err_len), intent(out) :: err
    character(len=:), allocatable :: rxstring
    
    integer i, ii
    logical l, m
    err = ''
    
    do i = 1,photomech%nrT-1
      do ii = i+1,photomech%nrT
        l = all(photomech%reactants_sp_inds(:,i) == photomech%reactants_sp_inds(:,ii))
        m = all(photomech%products_sp_inds(:,i) == photomech%products_sp_inds(:,ii))
        if ((m) .and. (l)) then
          err = "IOError: This reaction is a duplicate: "
          call reaction_string(photomech,i,rxstring)
          err(len_trim(err)+2:) = rxstring
          return
        endif
      enddo
    enddo
    
  end subroutine
  
  subroutine reaction_string(photomech,rxn,rxstring)
    type(PhotoMechanism), intent(in) :: photomech
    integer, intent(in) :: rxn
    character(len=:), allocatable, intent(out) :: rxstring
    integer j, k
    rxstring = ''
    do j = 1,photomech%nreactants(rxn)-1
      k = photomech%reactants_sp_inds(j,rxn)
      rxstring = rxstring //(trim(photomech%species_names(k))//' + ')
    enddo
    k = photomech%reactants_sp_inds(photomech%nreactants(rxn),rxn)
    rxstring = rxstring // trim(photomech%species_names(k))//' => '
    
    do j = 1,photomech%nproducts(rxn)-1
      k = photomech%products_sp_inds(j,rxn)
      rxstring = rxstring // trim(photomech%species_names(k))//' + '
    enddo
    k = photomech%products_sp_inds(photomech%nproducts(rxn),rxn)
    rxstring = rxstring // trim(photomech%species_names(k))
  end subroutine
  
  subroutine compare_rxtype_string(tmp, eqr, eqp, rxtype_int, err)
    type(string), allocatable, intent(in) :: eqr(:), eqp(:)
    type(string), intent(in) :: tmp
    integer, intent(in) :: rxtype_int
    character(len=err_len), intent(out) :: err
    character(len=15) :: rxtype
    integer i
    logical k, j, m, l, kk, jj
    l = .false.
    m = .false.
    k = .false.
    kk = .false.
    j = .false.
    jj = .false.
    err = ''
    if (rxtype_int == 0) then
      rxtype = 'photolysis'
    elseif (rxtype_int == 1) then
      rxtype = 'elementary'
    elseif (rxtype_int == 2) then
      rxtype = 'three-body'
    elseif (rxtype_int == 3) then
      rxtype = 'falloff'
    endif
  
    if ((trim(rxtype) == 'three-body') .or. (trim(rxtype) == 'falloff')) then
      do i = 1,size(eqr)
        if ((trim(eqr(i)%chars()) == 'M').and.j) jj = .true.
        if (trim(eqr(i)%chars()) == 'M') j = .true.
        if (trim(eqr(i)%chars()) == 'hv') m = .true.
      enddo
      do i = 1,size(eqp)
        if ((trim(eqp(i)%chars()) == 'M').and.k) kk = .true.
        if (trim(eqp(i)%chars()) == 'M') k = .true.
        if (trim(eqp(i)%chars()) == 'hv') l = .true.
      enddo
      if (trim(eqr(size(eqr))%chars()) /= 'M') then
        err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
                ' must have "M" as the last reactant'
        return
      endif
      if (trim(eqp(size(eqp))%chars()) /= 'M') then
        err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
                ' must have "M" as the last product'
        return
      endif
      if ((j).and.(k)) then
        ! good
      else
        err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
                ' must have "M" on both sides'
        return
      endif
      if ((jj).or.(kk)) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
                ' can only have one "M" on either side'
        return
      endif
      if ((m).or.(l)) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
                ' can not contain "hv". Only photolysis reactions can.'
        return
      endif
    elseif (trim(rxtype) == 'elementary') then
      do i = 1,size(eqr)
        if (trim(eqr(i)%chars()) == 'M') j = .true.
        if (trim(eqr(i)%chars()) == 'hv') m = .true.
      enddo
      do i = 1,size(eqp)
        if (trim(eqp(i)%chars()) == 'M') k = .true.
        if (trim(eqp(i)%chars()) == 'hv') l = .true.
      enddo
      if ((j).or.(k)) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
                ' can not contain "M".'
        return
      endif
      if ((m).or.(l)) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
                ' can not contain "hv". Only photolysis reactions can.'
        return
      endif
    elseif (trim(rxtype) == 'photolysis') then
      do i = 1,size(eqr)
        if (trim(eqr(i)%chars()) == 'M') j = .true.
        if ((trim(eqr(i)%chars()) == 'hv').and.(m)) jj = .true.
        if (trim(eqr(i)%chars()) == 'hv') m = .true.
      enddo
      do i = 1,size(eqp)
        if (trim(eqp(i)%chars()) == 'M') k = .true.
        if (trim(eqp(i)%chars()) == 'hv') l = .true.
      enddo
      if ((j).or.(k)) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
                ' can not contain "M".'
        return
      endif
      if (jj) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
                ' can only have one "hv" on the left side.'
        return
      endif
      if ((m).and..not.(l)) then
        ! good
      else
        err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
                ' must have "hv" on the left and no "hv" on the right.'
        return
      endif
    endif
  end subroutine
    
  
  subroutine get_thermodata(molecule, molecule_name, infile, &
                            thermo_temps_entry, thermo_data_entry, err)
    class(type_dictionary), intent(in) :: molecule
    character(len=*), intent(in) :: molecule_name
    character(len=*), intent(in) :: infile
    
    real(real_kind), intent(out) :: thermo_temps_entry(3)
    real(real_kind), intent(out) :: thermo_data_entry(7,2)
    character(len=err_len), intent(out) :: err
    
    type (type_error), pointer :: io_err
    class(type_dictionary), pointer :: tmpdict
    class(type_list), pointer :: tmplist
    class(type_list_item), pointer :: item, item1
    logical :: success
    
    integer :: j, i, k
    
    err = ''
    thermo_temps_entry = -1.d0
    thermo_data_entry = -1.d0
    
    tmpdict => molecule%get_dictionary("thermo",.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    
    ! check thermodynamic model
    if (tmpdict%get_string("model",error = io_err) /= "Shomate") then
      err = "IOError: Thermodynamic data must be in Shomate format for "//trim(molecule_name)
      return
    endif
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    
    ! get temperature ranges
    tmplist =>tmpdict%get_list("temperature-ranges",.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    j = 1
    item => tmplist%first
    do while (associated(item))
      select type (listitem => item%node)
      class is (type_scalar)
        if (j > 3) then
          err = "IOError: Too many temperature ranges for "//trim(molecule_name)
          return
        endif
        thermo_temps_entry(j) = listitem%to_real(-1.d0,success)
        if (.not. success) then
          err = "IOError: Problem reading thermodynamic data for  "//trim(molecule_name)
          return
        endif
      class default
        err = "IOError: Problem reading thermodynamic data for "//trim(molecule_name)
        return
      end select
      item => item%next
      j = j + 1
    enddo
    
    ! check amount of thermodynamic data
    if (thermo_temps_entry(3) < -0.5d0) then
      i = 1
    else
      i = 2
    endif
    
    ! get data
    tmplist => tmpdict%get_list("data",.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    k = 1
    item => tmplist%first
    do while (associated(item))
      if (k > i) then
        err = "IOError: Too much thermodynamic data for "//trim(molecule_name)
        return
      endif
      select type (listitem => item%node)
      class is (type_list)
        j = 1
        item1 => listitem%first
        do while (associated(item1)) 
          select type (listitem1 => item1%node)
          class is (type_scalar)
            if (j > 7) then
              err = "IOError: Too much thermodynamic data for "//trim(molecule_name)
              return
            endif
            thermo_data_entry(j, k) = listitem1%to_real(-1.d0,success)
            if (.not.success) then
              err = "IOError: Problem reading thermodynamic data for "//trim(molecule_name)
              return
            endif
          class default
            err = "IOError: Problem reading thermodynamic data for "//trim(molecule_name)
            return
          end select
        item1 => item1%next
        j = j + 1
        enddo
        if (j-1 /= 7) then
          err = "IOError: Missing thermodynamic data for "//trim(molecule_name)
          return
        endif
      class default
        err = "IOError: Problem reading thermodynamic data for "//trim(molecule_name)
        return
      end select
      item => item%next
      k = k + 1
    enddo
    if (k - 1 /= i) then
      err = "IOError: More temperature ranges than thermodynamic data for "//trim(molecule_name)
      return
    endif
    
  end subroutine
  
  subroutine get_rateparams(reaction, infile, rxtype, falloff_type, rateparam, err)
    class(type_dictionary), intent(in) :: reaction
    character(len=*), intent(in) :: infile
    integer, intent(out) :: rxtype
    integer, intent(out) :: falloff_type
    real(real_kind), intent(out) :: rateparam(10)
    character(len=err_len), intent(out) :: err
    
    type (type_error), pointer :: io_err
    class(type_dictionary), pointer :: tmpdict
    character(len=15) :: rxtype_str
    err = ''
    rateparam = 0.d0
    ! no error possible
    rxtype_str = reaction%get_string("type","",error = io_err) 
    if (rxtype_str == 'photolysis') then
      rxtype = 0
    elseif ((rxtype_str == 'elementary') .or. (rxtype_str == '')) then
      rxtype = 1
    elseif (rxtype_str == 'three-body') then
      rxtype = 2
    elseif (rxtype_str == 'falloff') then
      rxtype = 3
    else
      err = 'IOError: reaction type '//trim(rxtype_str)//' is not a valid reaction type.'
      return
    endif
    falloff_type = -huge(1)
    
    ! get params
    if ((rxtype_str == 'elementary') .or. (rxtype_str == 'three-body')) then
      tmpdict => reaction%get_dictionary('rate-constant',.true.,error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      rateparam(1) = tmpdict%get_real('A',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      rateparam(2) = tmpdict%get_real('b',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      rateparam(3) = tmpdict%get_real('Ea',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    elseif (rxtype_str == 'falloff') then
      tmpdict => reaction%get_dictionary('low-P-rate-constant',.true.,error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      rateparam(1) = tmpdict%get_real('A',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      rateparam(2) = tmpdict%get_real('b',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      rateparam(3) = tmpdict%get_real('Ea',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      tmpdict => reaction%get_dictionary('high-P-rate-constant',.true.,error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      rateparam(4) = tmpdict%get_real('A',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      rateparam(5) = tmpdict%get_real('b',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      rateparam(6) = tmpdict%get_real('Ea',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
      nullify(tmpdict)
      tmpdict => reaction%get_dictionary('Troe',.false.,error = io_err)
      if (associated(tmpdict)) then
        rateparam(7) = tmpdict%get_real('A',error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        rateparam(8) = tmpdict%get_real('T1',error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif

        rateparam(9) = tmpdict%get_real('T2',error = io_err)
        if (associated(io_err)) then ! T2 is not there
          falloff_type = 1
          nullify(io_err)
        else ! T2 is there
          falloff_type = 2
        endif
        
        rateparam(10) = tmpdict%get_real('T3',error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      else
        falloff_type = 0 ! no falloff function
      endif
        
    endif
    
  end subroutine
  
  subroutine parse_reaction(instring, reverse, eqr, eqp, err)
    type(string), intent(in) :: instring
    logical, intent(out) :: reverse
    type(string), allocatable, intent(out) :: eqr(:), eqp(:)
    character(len=err_len), intent(out) :: err
    
    type(string) :: string1, string2, string3
    type(string), allocatable :: eq1(:), eq2(:)
    type(string), allocatable :: eqr1(:), eqp1(:)
    integer i
    string1 = instring%replace(old='(', new=' ')
    string2 = string1%replace(old=')', new=' ')
    string3 = string2%replace(old='+', new=' ')
    if (index(instring%chars(), "<=>") /= 0) then
      call string2%split(eq1, sep="<=>")
      call string3%split(eq2, sep="<=>")
      reverse = .true.
    elseif (index(instring%chars(), " =>") /= 0) then
      call string2%split(eq1, sep="=>")
      call string3%split(eq2, sep="=>")
      reverse = .false.
    else
      err = "IOError: Invalid reaction arrow in reaction "//instring// &
            '. Note, forward reactions must have a space before the arrow, like " =>"'
      return
    endif
    
    call eq1(1)%split(eqr1, sep="+")
    call eq1(2)%split(eqp1, sep="+")
    ! remove white space
    allocate(eqr(size(eqr1)))
    allocate(eqp(size(eqp1)))
    do i=1,size(eqr1)
      eqr(i) = eqr1(i)%replace(old=' ', new='')
    enddo
    do i=1,size(eqp1)
      eqp(i) = eqp1(i)%replace(old=' ', new='')
    enddo
    
    call eq2(1)%split(eqr1, sep=" ")
    call eq2(2)%split(eqp1, sep=" ")
    
    if ((size(eqr1) /= size(eqr)) .or. (size(eqp1) /= size(eqp))) then
      err = 'IOError: Missing "+" sign(s) in reaction '//instring
      return
    endif
    
    if (size(eqr) + size(eqp) /= instring%count('+') + 2) then
      err = 'IOError: Too many "+" signs in reaction '//instring
      return
    endif
  end subroutine
  
  subroutine get_reaction_chars(instring, max_num_react, max_num_prod, numr, nump, &
                                outreact, outprod, reverse, err)
    type(string), intent(in) :: instring
    integer, intent(in) :: max_num_react, max_num_prod
    
    integer, intent(out) :: numr, nump
    character(len=8), intent(out) :: outreact(max_num_react), outprod(max_num_prod)
    logical, intent(out) :: reverse
    character(len=err_len), intent(out) :: err
    
    type(string), allocatable :: eqr(:), eqp(:)
    integer :: i
    
    call parse_reaction(instring, reverse, eqr, eqp, err)
    if (len_trim(err) > 0) return
    
    numr = size(eqr)
    nump = size(eqp)
    
    outreact = ''
    outprod = ''
    do i=1,numr
      outreact(i) = eqr(i)%chars()
    enddo
    do i=1,nump
      outprod(i) = eqp(i)%chars()
    enddo
  end subroutine
  
  subroutine get_reaction_sp_nums(reaction, max_num_react, max_num_prod, reacts, prods, &
                                 species_names, species_composition, natoms, nsp, &
                                 react_sp_nums, prod_sp_nums, err)
    type(string), intent(in) :: reaction
    integer, intent(in) :: max_num_react, max_num_prod
    character(len=8), intent(in) :: reacts(max_num_react)
    character(len=8), intent(in) :: prods(max_num_prod)
    integer, intent(in) :: nsp, natoms
    character(len=8), intent(in) :: species_names(nsp+2)
    integer, intent(in) :: species_composition(natoms,nsp+2)
    
    integer, intent(out) :: react_sp_nums(max_num_react)
    integer, intent(out) :: prod_sp_nums(max_num_prod)
    character(len=err_len), intent(out) :: err
    
    integer :: i, ind(1)
    integer :: reactant_atoms(natoms), product_atoms(natoms)
    err = ''
    reactant_atoms = 0
    product_atoms = 0
    
    do i = 1,max_num_react
      ind = findloc(species_names,reacts(i))
      react_sp_nums(i) = ind(1)
      if ((reacts(i) /= '') .and. (ind(1) == 0)) then
        err = "IOError: "// & 
               "Species "//trim(reacts(i))//" in reaction "//reaction// &
               "is not in the list of species."
        return
      endif
    enddo
    
    do i = 1,max_num_prod
      ind = findloc(species_names,prods(i))
      prod_sp_nums(i) = ind(1)
      if ((prods(i) /= '') .and. (ind(1) == 0)) then
        err = "IOError: "// & 
               "Species "//trim(prods(i))//" in reaction "//reaction// &
               "is not in the list of species."
        return
      endif
    enddo
  
    do i=1,max_num_react
      if (react_sp_nums(i) /= 0) then
        reactant_atoms = reactant_atoms + species_composition(:,react_sp_nums(i))
      endif
    enddo
    do i=1,max_num_prod
      if (prod_sp_nums(i) /= 0) then
        product_atoms = product_atoms + species_composition(:,prod_sp_nums(i))
      endif
    enddo
    if (.not. all(reactant_atoms == product_atoms)) then
      err = "IOError: "//& 
             'Bad mass balance in reaction "'//reaction// &
             '". You could have messed up how many atoms one of the species has.'
      return
    endif
  end subroutine
  
  subroutine get_photoset(infile, photomech, photoset, err)
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoSettings), intent(out) :: photoset
    character(len=err_len), intent(out) :: err
  
    character(error_length) :: error
    class (type_node), pointer :: root
    err = ''
    
    root => parse(infile,unit=100,error=error)
    if (error /= '') then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        call unpack_settings(root, infile, photomech, photoset, err)
        call root%finalize()
        deallocate(root)
      class default
        err = trim(infile)//" file must have dictionaries at root level"
        return
    end select
    ! tmpdict => photoset%get_dictionary("atmosphere-grid",.true.,error = io_err)
    ! if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif

  end subroutine
  
  subroutine unpack_settings(mapping, infile, photomech, photoset, err)
    class (type_dictionary), intent(in), pointer :: mapping
    character(len=*), intent(in) :: infile
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoSettings), intent(out) :: photoset
    character(len=err_len), intent(out) :: err
    
    class (type_dictionary), pointer :: tmp1, tmp2
    class (type_list), pointer :: bcs
    class (type_list_item), pointer :: item
    type (type_error), pointer :: io_err
    
    character(len=str_len) :: background_gas, spec_type, spec_name
    integer :: j, i, ind(1), ind1(1), ll, sl
    character(len=20), allocatable :: dups(:)
    real(real_kind) :: grav
    
    err = ''
    
    ! photolysis grid
    tmp1 => mapping%get_dictionary('photolysis-grid',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photoset%regular_grid = tmp1%get_logical('regular-grid',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    if (photoset%regular_grid) then
      photoset%lower_wavelength = tmp1%get_real('lower-wavelength',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      photoset%upper_wavelength = tmp1%get_real('upper-wavelength',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      photoset%nw = tmp1%get_real('number-of-bins',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      if (photoset%nw < 1) then 
        err = 'Number of photolysis bins must be >= 1 in '//trim(infile)
        return
      endif
      if (photoset%lower_wavelength > photoset%upper_wavelength) then
        err = 'lower-wavelength must be smaller than upper-wavelength in '//trim(infile)
        return
      endif
    else
      photoset%grid_file = tmp1%get_string('input-file',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    endif
    
    ! atmosphere grid
    tmp1 => mapping%get_dictionary('atmosphere-grid',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photoset%bottom_atmos = tmp1%get_real('bottom',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photoset%top_atmos = tmp1%get_real('top',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photoset%nz = tmp1%get_real('number-of-layers',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    ! Planet
    tmp1 => mapping%get_dictionary('planet',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    photoset%back_gas = tmp1%get_logical('use-background-gas',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    if (photoset%back_gas) then
      background_gas = tmp1%get_string('background-gas',error = io_err)
      ind = findloc(photomech%species_names,trim(background_gas))
      if (ind(1) == 0) then
        err = 'IOError: Background gas "'//trim(background_gas)// &
              '" is not one of the species in the reaction mechanism.'
        return
      else
        photoset%back_gas_ind = ind(1)
      endif
    else
      photoset%back_gas_ind = -1
      ! err = "Currently, the model requires there to be a background gas."// &
      !       " You must set 'use-background-gas: true'"
      ! return
    endif
    
    
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photoset%surface_pressure = tmp1%get_real('surface-pressure',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photoset%surface_pressure < 0.d0) then
      err = 'IOError: Planet surface pressure must be greater than zero.'
      return
    endif
    photoset%planet_mass = tmp1%get_real('planet-mass',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photoset%planet_mass < 0.d0) then
      err = 'IOError: Planet mass must be greater than zero.'
      return
    endif
    photoset%planet_radius = tmp1%get_real('planet-radius',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photoset%planet_radius < 0.d0) then
      err = 'IOError: Planet radius must be greater than zero.'
      return
    endif
    photoset%surface_albedo = tmp1%get_real('surface-albedo',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photoset%surface_albedo < 0.d0) then
      err = 'IOError: Surface albedo must be greater than zero.'
      return
    endif
    tmp2 => tmp1%get_dictionary('water',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photoset%water_sat_trop = tmp2%get_logical('water-saturated-troposphere',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    ind = findloc(photomech%species_names,'H2O')
    if (ind(1) == 0) then
      err = 'IOError: H2O must be a species if water-saturated-troposhere = True.'
      return
    endif
    if (photoset%water_sat_trop) then  
      photoset%trop_alt = tmp2%get_real('tropopause-altitude',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      if ((photoset%trop_alt < photoset%bottom_atmos) .or. &
          (photoset%trop_alt > photoset%top_atmos)) then
          err = 'IOError: tropopause-altitude must be between the top and bottom of the atmosphere'
          return
      endif
    endif
    
    ! boundary conditions and species types
    bcs => mapping%get_list('boundary-conditions',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    allocate(photoset%species_type(photomech%nsp))
    photoset%species_type = 1 ! long lived by default
    if (photoset%back_gas) then
      photoset%species_type(photoset%back_gas_ind) = 0 ! background gas
    endif
      
    ! determine number of short lived species. 
    allocate(dups(photomech%nsp))
    j = 1
    photoset%nsl = 0
    item => bcs%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        if (j > photomech%nsp) then
          err = "IOError: Too many boundary condition entries in settings file."
          return
        endif
        dups(j) = element%get_string('name',error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        ! check for duplicates
        if (j > 1) then
          do i = 1,j-1
            if (dups(j) == dups(i)) then
              err = "IOError: Species "//trim(dups(i))// &
              " has more than one boundary conditions entry in the settings file."
              return
            endif
          enddo
        endif
        ! check if in rxmech
        ind = findloc(photomech%species_names,dups(j))
        if (ind(1) == 0) then
          err = "IOError: Species "//trim(dups(j))// &
          ' in settings file is not in the reaction mechanism file.'
          return 
        endif
        if ((ind(1) == photoset%back_gas_ind) .and. (photoset%back_gas)) then ! can't be background gas
          err = "IOError: Species "//trim(dups(j))// &
          ' in settings file is the background gas, and can not have boundary conditions.'
          return
        endif
        
        spec_type = element%get_string('type','long lived',error = io_err)
        if (spec_type == 'short lived') then
          photoset%nsl = photoset%nsl + 1
          photoset%species_type(ind(1)) = 2
        elseif (spec_type == 'long lived') then
          ! do nothing
        else
          err = 'IOError: species type '//trim(spec_type)//' is not a valid.' 
          return
        endif
      class default
        err = "IOError: Boundary conditions must be a list of dictionaries."
        return
      end select 
      item => item%next
      j = j + 1
    enddo
    
    ! allocate book keeping
    if (photoset%back_gas) then
      photoset%nq = photomech%nsp - photoset%nsl  - 1 ! minus one for background gas
    else
      photoset%nq = photomech%nsp - photoset%nsl 
    endif
    allocate(photoset%LL_inds(photoset%nq))
    allocate(photoset%SL_inds(photoset%nsl))
    sl = 0
    ll = 0
    do i = 1,photomech%nsp
      if (photoset%species_type(i) == 0) then
        ! nothing. background gas
      elseif (photoset%species_type(i) == 1) then
        ll = ll + 1
        photoset%LL_inds(ll) = i
      elseif (photoset%species_type(i) == 2) then
        sl = sl + 1
        photoset%SL_inds(sl) = i
      else
        err = "IOError: Problem reading in species types in settings file."
        return
      endif
    enddo
    
    ! allocate boundary conditions
    allocate(photoset%lowerboundcond(photoset%nq))
    allocate(photoset%lower_vdep(photoset%nq))
    allocate(photoset%lower_flux(photoset%nq))
    allocate(photoset%lower_dist_height(photoset%nq))
    allocate(photoset%lower_fix_mr(photoset%nq))
    allocate(photoset%upperboundcond(photoset%nq))
    allocate(photoset%upper_veff(photoset%nq))
    allocate(photoset%upper_flux(photoset%nq))
    ! default boundary conditions
    photoset%lowerboundcond = 0
    photoset%lower_vdep = 0.d0
    photoset%upperboundcond = 0
    photoset%upper_veff = 0.d0
  
    item => bcs%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        spec_name = element%get_string('name',error = io_err)
        ind = findloc(photomech%species_names,spec_name)
    
        spec_type = element%get_string('type','long lived',error = io_err)
        if (spec_type == 'long lived') then
          ! find proper index
          ind1 = findloc(photoset%LL_inds,ind(1))
          i = ind1(1)
          ! get boundary condition
          call get_boundaryconds(element, spec_name, infile, &
                                 photoset%lowerboundcond(i), photoset%lower_vdep(i), &
                                 photoset%lower_flux(i), photoset%lower_dist_height(i), &
                                 photoset%lower_fix_mr(i), &
                                 photoset%upperboundcond(i), photoset%upper_veff(i), &
                                 photoset%upper_flux(i), err)
          if (len_trim(err) /= 0) return
        endif
      end select 
      item => item%next
    enddo    
    
    
    ! check for SL nonlinearities
    call check_sl(photomech, photoset, err)
    if (len_trim(err) /= 0) return

    ! ! set up the atmosphere grid
    ! allocate(photoset%z(photoset%nz))
    ! allocate(photoset%dz(photoset%nz))
    ! call vertical_grid(photoset%bottom_atmos,photoset%top_atmos,photoset%nz, &
    !                    photoset%z, photoset%dz)
    ! 
    ! ! compute the gravity
    ! allocate(photoset%grav(photoset%nz))
    ! do i =1,photoset%nz
    !   call compute_gravity((photoset%planet_radius + photoset%z(i))/1.d2, &
    !                         photoset%planet_mass/1.d3, grav)
    !   photoset%grav(i) = grav*1.d2 ! convert to cgs
    ! enddo
    
  end subroutine
  
  subroutine check_sl(photomech, photoset, err)
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoSettings), intent(in) :: photoset
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, l, k, kk, m, mm, n, ind(1), counter
    character(len=:), allocatable :: reaction
    err = ''
    
    do i = 1, photoset%nsl
      j = photoset%SL_inds(i)
      ! can not be an efficiency.
      do k = 1,photomech%nrF
        if ((photomech%rxtypes(k) == 2) .or. (photomech%rxtypes(k) == 3)) then ! if three body or falloff
          ind = findloc(photomech%eff_sp_inds(:,k),j)
          if (ind(1) /= 0) then
            call reaction_string(photomech,k,reaction)
            err = 'IOError: Reaction "'//reaction//'" has short-lived species collision efficiencies.' // &
            ' This is not allowed. Either remove the efficiencies, or change the species to long lived.'
            return
          endif
        endif
      enddo
      
      
      l = photomech%nump(j)
      do k = 1,l
        kk = photomech%iprod(k,j)
        call reaction_string(photomech,kk,reaction)
        m = photomech%nreactants(kk)
        do mm = 1, m
          ! are SL species produced by other SL species?
          ind = findloc(photoset%SL_inds,photomech%reactants_sp_inds(mm,kk))
          if (ind(1) /= 0) then
            err = 'IOError: Reaction "'//reaction//'" has short-lived species as reactants'// &
            ' and products. This is not allowed. Change one or both of the species to long-lived.'
            return
          endif
        enddo
      enddo

      l = photomech%numl(j)
      do k = 1,l
        kk = photomech%iloss(k,j)
        call reaction_string(photomech,kk,reaction)
        m = photomech%nreactants(kk)
        counter = 0
        do mm = 1, m
          n = photomech%reactants_sp_inds(mm,kk)
          ind = findloc(photoset%SL_inds,n)
          if ((ind(1) /= 0) .and. (n == j)) then
            counter = counter + 1
            if (counter > 1) then
              err = 'IOError: Reaction "'//reaction//'" short lived species react'// &
              ' with themselves. This is not allowed. Change the species to long lived.'
              return
            endif
          elseif ((ind(1) /= 0) .and. (n /= j)) then
            err = 'IOError: Reaction "'//reaction//'" short lived species react'// &
            ' with other short lived species. This is not allowed.'
            return
          elseif (photomech%species_names(n) == 'hv') then
            err = 'IOError: Photolysis reaction "'//reaction//'" can not have short lived species.'
            return
          endif
        enddo
      enddo
    enddo
    
  end subroutine
  
  
  subroutine compute_gravity(radius, mass, grav)
    use photochem_const, only: G_grav
    real(real_kind), intent(in) :: radius, mass
    real(real_kind), intent(out) :: grav
    grav = G_grav * mass / radius**2.d0
  end subroutine

  subroutine vertical_grid(bottom, top, nz, z, dz)
    real(real_kind), intent(in) :: bottom, top
    integer, intent(in) :: nz
    real(real_kind), intent(out) :: z(nz), dz(nz)
    
    integer :: i
  
    dz = (top - bottom)/nz
    z(1) = dz(1)/2.d0
    do i = 2,nz
      z(i) = z(i-1) + dz(i)
    enddo
  end subroutine
  
  subroutine get_boundaryconds(molecule, molecule_name, infile, &
                               lowercond, Lvdep, Lflux, LdistH, Lmr, &
                               uppercond, Uveff, Uflux, err)
    class(type_dictionary), intent(in) :: molecule
    character(len=*), intent(in) :: molecule_name
    character(len=*), intent(in) :: infile
    
    integer, intent(out) :: lowercond
    real(real_kind), intent(out) :: Lvdep, Lflux, LdistH, Lmr
    integer, intent(out) :: uppercond
    real(real_kind), intent(out) :: Uveff, Uflux
    character(len=err_len), intent(out) :: err
    
    type (type_error), pointer :: io_err
    class(type_dictionary), pointer :: tmpdict
    character(len=:), allocatable :: bctype
    err = ''

    tmpdict => molecule%get_dictionary("lower-boundary",.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    bctype = tmpdict%get_string("type",error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    
    ! constant deposition velocity = vdep
    ! constant mixing-ratio = mixing-ratio
    ! constant flux = flux
    ! deposition velocity + distributed flux = vdep + flux + distributed-height
    if (bctype == "constant deposition velocity") then
      lowercond = 0
      Lvdep = tmpdict%get_real("vdep",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      Lflux = 0.d0
      LdistH = 0.d0
      Lmr = 0.d0 
    elseif (bctype == "constant mixing-ratio") then
      lowercond = 1
      Lvdep = 0.d0
      Lflux = 0.d0
      LdistH = 0.d0
      Lmr = tmpdict%get_real("mixing-ratio",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    elseif (bctype == "constant flux") then
      lowercond = 2
      Lvdep = 0.d0
      Lflux = tmpdict%get_real("flux",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      LdistH = 0.d0
      Lmr = 0.d0 
    elseif (bctype == "deposition velocity + distributed flux") then
      lowercond = 3
      Lvdep = tmpdict%get_real("vdep",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      Lflux = tmpdict%get_real("flux",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      LdistH = tmpdict%get_real("distributed-height",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      Lmr = 0.d0 
    else
      err = 'IOError: "'//trim(bctype)//'" is not a valid lower boundary condition for '//trim(molecule_name)
      return
    endif
    
    tmpdict => molecule%get_dictionary("upper-boundary",.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    bctype = tmpdict%get_string("type",error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    
    ! constant effusion velocity = veff
    ! constant flux = flux
    if (bctype == "constant effusion velocity") then
      uppercond = 0
      Uveff = tmpdict%get_real("veff",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      Uflux = 0.d0
    elseif (bctype == "constant flux") then
      uppercond = 2
      Uveff = 0.d0
      Uflux = tmpdict%get_real("flux",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    else
      err = 'IOError: "'//trim(bctype)//'" is not a valid upper boundary condition for '//trim(molecule_name)
      return
    endif
    
  end subroutine
    
  subroutine get_photorad(photomech, photoset, photorad, err)
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoSettings), intent(in) :: photoset
    type(PhotoRadTran), intent(out) :: photorad
    character(len=err_len), intent(out) :: err
    
    integer :: i
    err = ''
    
    ! compute wavelength grid
    if (photoset%regular_grid) then
      photorad%nw = photoset%nw
      allocate(photorad%wavl(photoset%nw+1))
      photorad%wavl(1) = photoset%lower_wavelength
      do i = 2,photoset%nw+1
        photorad%wavl(i) = photorad%wavl(i-1) + &
                           (photoset%upper_wavelength - photoset%lower_wavelength)/photoset%nw
      enddo
    else
      ! read file
      err = 'Still need to add support for reading in wavelength grid from file'
      return
    endif
    
    ! get rayleigh
    call get_rayleigh(photomech, photorad, err)
    if (len_trim(err) /= 0) return
    
    ! get photolysis xsections data
    call get_photolysis_xs(photomech, photorad, err)
    if (len_trim(err) /= 0) return
    
  end subroutine
  
  subroutine get_photolysis_xs(photomech, photorad, err)
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoRadTran), intent(inout) :: photorad ! inout!
    character(len=err_len), intent(out) :: err
    
    integer, parameter :: maxcols = 200
    character(len=:), allocatable :: xsroot
    character(len=:), allocatable :: filename, xsfilename, reaction
    character(len=str_len) :: line
    character(len=100) :: tmp(maxcols), tmp1
    real(real_kind), allocatable :: file_xs(:,:), file_qy(:,:), file_wav(:), file_line(:)
    real(real_kind), allocatable :: dumby(:,:)
    real(real_kind), parameter :: rdelta = 1.d-4
    
    integer :: i, j, k, l, m, mm, io, kk, ierr
    err = ''
    
    xsroot = "../data/xsections/"
    
    ! count temperature columns
    allocate(photorad%num_temp_cols(photomech%kj))
    allocate(photorad%sum_temp_cols(photomech%kj))
    do i = 1,photomech%kj
      filename = ''
      j = photomech%photonums(i)
      call reaction_string(photomech,j,reaction)
      filename = reaction
      do k = 1,len(filename)-1
        if (filename(k:k) == ' ') then
          filename(k:k) = '_'
        endif
      enddo
      filename = filename//'.txt'

      k = photomech%reactants_sp_inds(1,j)
      filename = trim(photomech%species_names(k))//'/'//filename
      xsfilename = trim(photomech%species_names(k))//'/'//trim(photomech%species_names(k))//'_xs.txt'
      open(101, file=xsroot//filename,status='old',iostat=io)
      if (io /= 0) then
        err = 'The photolysis reaction '//reaction//' does not have quantum yield data'
        return
      endif
      read(101,*)
      read(101,'(A)') line
      do k=1,maxcols
        read(line,*,iostat=io) tmp(1:k)
        if (io /= 0) exit
      enddo
      if (k == maxcols+1) then
        err = 'More cross section temperature data than allowed for reaction '//reaction
        return
      endif
      photorad%num_temp_cols(i) = k - 2
      close(101)
    enddo
    photorad%sum_temp_cols(1) = 0
    do i = 2,photomech%kj
      photorad%sum_temp_cols(i) = photorad%sum_temp_cols(i-1) + photorad%num_temp_cols(i-1)
    enddo
    
    ! allocate
    allocate(photorad%xs_data(sum(photorad%num_temp_cols)*photorad%nw))
    allocate(photorad%xs_data_temps(maxval(photorad%num_temp_cols),photomech%kj))
    photorad%xs_data_temps = 0.d0

    ! read in data
    do i = 1,photomech%kj
      filename = ''
      j = photomech%photonums(i)
      call reaction_string(photomech,j,reaction)
      filename = reaction
      do k = 1,len(filename)-1
        if (filename(k:k) == ' ') then
          filename(k:k) = '_'
        endif
      enddo
      filename = filename//'.txt'

      k = photomech%reactants_sp_inds(1,j)
      filename = trim(photomech%species_names(k))//'/'//filename
      xsfilename = trim(photomech%species_names(k))//'/'//trim(photomech%species_names(k))//'_xs.txt'
      open(101, file=xsroot//filename,status='old',iostat=io)
      read(101,*)
      read(101,'(A)') line
       
      do k=1,photorad%num_temp_cols(i) + 1
        read(line,*,iostat=io) tmp(1:k)
      enddo
      do k=1,photorad%num_temp_cols(i)
        tmp1 = tmp(k+1)
        read(tmp1(1:index(tmp1,'K')-1),*,iostat=io) photorad%xs_data_temps(k,i)
        if (io /= 0) then
          err = 'Problem reading in cross sections for reaction '//reaction
          return
        endif
      enddo
      
      ! count lines
      k = 0
      do while(io == 0)
        read(101,*, iostat=io)
        if (io == 0) k = k + 1
      enddo
      allocate(file_wav(k+4))
      allocate(file_qy(k+4,photorad%num_temp_cols(i)))
      allocate(file_line(photorad%num_temp_cols(i)+1))
      allocate(dumby(photorad%nw,photorad%num_temp_cols(i)))

      rewind(101)
      read(101,*)
      read(101,*)
      do l = 1, k
        read(101,'(A)',iostat=io) line
        read(line,*) file_wav(l)
        do m = 1,photorad%num_temp_cols(i)+1
          read(line,*) file_line(1:m)
        enddo
        file_qy(l,:) = file_line(2:)
      enddo
      
      ! interpolate to grid. save in photorad%xs_data
      ierr = 0
      do l = 1, photorad%num_temp_cols(i)
        kk = k
        call addpnt(file_wav, file_qy(:,l), kk+4, k, file_wav(1)*(1.d0-rdelta), 0.d0, ierr)
        call addpnt(file_wav, file_qy(:,l), kk+4, k, 0.d0, 0.d0, ierr)
        call addpnt(file_wav, file_qy(:,l), kk+4, k, file_wav(k)*(1.d0+rdelta), 0.d0, ierr)
        call addpnt(file_wav, file_qy(:,l), kk+4, k, huge(rdelta), 0.d0, ierr)
        if (ierr /= 0) then
          err = 'IOError: Problem interpolating quantum yield data to photolysis grid for reaction '// &
                trim(reaction)
          return
        endif
        m = ((l-1)*photorad%nw + 1) + (photorad%sum_temp_cols(i)*photorad%nw)
        mm = (photorad%nw*l) + (photorad%sum_temp_cols(i)*photorad%nw)
        call inter2(photorad%nw+1,photorad%wavl,dumby(:,l), &
                    kk+4,file_wav,file_qy(:,l),ierr)
        photorad%xs_data(m:mm) = dumby(:,l)
        k = kk
        if (ierr /= 0) then
          err = 'IOError: Problem interpolating quantum yield data to photolysis grid for reaction '// &
                trim(reaction)
          return
        endif
      enddo
      
      close(101)
      deallocate(file_wav,file_qy)
      
      ! now do xs
      open(102, file=xsroot//xsfilename,status='old',iostat=io)
      if (io /= 0) then
        err = 'The photolysis reaction '//reaction//' does not have cross section data'
        return
      endif
      ! count lines
      k = -2 ! skip first two lines
      do while(io == 0)
        read(102,*, iostat=io)
        if (io == 0) k = k + 1
      enddo
      allocate(file_wav(k+4))
      allocate(file_xs(k+4,photorad%num_temp_cols(i)))
      
      rewind(102)
      read(102,*)
      read(102,*)
      do l = 1, k
        read(102,'(A)',iostat=io) line
        read(line,*) file_wav(l)
        do m = 1,photorad%num_temp_cols(i)+1
          read(line,*) file_line(1:m)
        enddo
        file_xs(l,:) = file_line(2:)
      enddo

      ierr = 0
      do l = 1, photorad%num_temp_cols(i)
        kk = k
        call addpnt(file_wav, file_xs(:,l), kk+4, k, file_wav(1)*(1.d0-rdelta), 0.d0,ierr)
        call addpnt(file_wav, file_xs(:,l), kk+4, k, 0.d0, 0.d0,ierr)
        call addpnt(file_wav, file_xs(:,l), kk+4, k, file_wav(k)*(1.d0+rdelta), 0.d0,ierr)
        call addpnt(file_wav, file_xs(:,l), kk+4, k, huge(rdelta), 0.d0,ierr)
        if (ierr /= 0) then
          err = 'IOError: Problem interpolating xs data to photolysis grid for reaction '// &
                trim(reaction)
          return
        endif
        m = ((l-1)*photorad%nw + 1) + (photorad%sum_temp_cols(i)*photorad%nw)
        mm =(photorad%nw*l) + (photorad%sum_temp_cols(i)*photorad%nw)

        call inter2(photorad%nw+1,photorad%wavl,dumby(:,l), &
                    kk+4,file_wav,file_xs(:,l),ierr)
        photorad%xs_data(m:mm) = photorad%xs_data(m:mm) * dumby(:,l)
        k = kk
        
        if (ierr /= 0) then
          err = 'IOError: Problem interpolating xs data to photolysis grid for reaction '// &
                trim(reaction)
          return
        endif
      enddo

      close(102)
      deallocate(file_xs, file_wav, file_line, dumby)

    enddo
    
  end subroutine
  
  
  subroutine get_rayleigh(photomech, photorad, err)
    use yaml, only : parse, error_length
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoRadTran), intent(inout) :: photorad ! inout!
    character(len=err_len), intent(out) :: err
    
    real(real_kind), allocatable :: A(:), B(:), Delta(:)
    character(len=str_len) :: rayleigh_file

    character(error_length) :: error
    class (type_node), pointer :: root
    integer :: i, j
    err = ''
    
    rayleigh_file = "../data/xsections/rayleigh.yaml"
    
    ! parse yaml file
    root => parse(rayleigh_file,unit=100,error=error)
    if (error /= '') then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        call rayleigh_params(root,photomech,trim(rayleigh_file),err, &
                             photorad%raynums, A, B, Delta)
        if (len_trim(err) /= 0) return
        call root%finalize()
        deallocate(root)
    end select
    
    ! compute cross sections
    allocate(photorad%sigray(size(A),photorad%nw))
    do i = 1,photorad%nw
      do j = 1,size(A)
        call rayleigh_vardavas(A(j), B(j), Delta(j), photorad%wavl(i), &
                               photorad%sigray(j, i))
      enddo
    enddo
    deallocate(A,B,Delta)
  end subroutine
  
  subroutine rayleigh_params(mapping,photomech,infile,err, raynums, A, B, Delta)
    class (type_dictionary), intent(in), pointer :: mapping
    type(PhotoMechanism), intent(in) :: photomech
    character(len=*), intent(in) :: infile
    character(len=err_len), intent(out) :: err
    
    class (type_key_value_pair), pointer :: key_value_pair
    class (type_dictionary), pointer :: tmp1, tmp2
    type (type_error), pointer :: io_err
    real(real_kind), allocatable, intent(out) :: A(:), B(:), Delta(:)
    integer, allocatable, intent(out) :: raynums(:)
    
    integer :: j, ind(1)  
    
    j = 0
    key_value_pair => mapping%first
    do while (associated(key_value_pair))
      ind = findloc(photomech%species_names,trim(key_value_pair%key))
      if (ind(1) /= 0) then 
        j = j + 1
      endif
      key_value_pair => key_value_pair%next
    enddo
    allocate(raynums(j))
    allocate(A(j))
    allocate(B(j))
    allocate(Delta(j))
    j = 1
    key_value_pair => mapping%first
    do while (associated(key_value_pair))
      ind = findloc(photomech%species_names,trim(key_value_pair%key))
      if (ind(1) /= 0) then 
        raynums(j) = ind(1)
        
        tmp1 => mapping%get_dictionary(key_value_pair%key,.true.,error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        if (trim(tmp1%get_string('formalism',error=io_err)) /= 'vardavas') then
          err = "Unknown formalism for Rayleigh cross section for "//trim(key_value_pair%key)
          return
        endif
        tmp2 => tmp1%get_dictionary("data",.true.,error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        Delta(j) = tmp2%get_real("Delta",error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        A(j) = tmp2%get_real("A",error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        B(j) = tmp2%get_real("B",error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        j = j + 1
      endif
      key_value_pair => key_value_pair%next
    enddo
  end subroutine
  
  subroutine rayleigh_vardavas(A, B, Delta, lambda, sigray)
    real(real_kind), intent(in) :: A, B, Delta, lambda
    real(real_kind), intent(out) :: sigray
    
    sigray = 4.577d-21*((6.d0+3.d0*Delta)/(6.d0-7.d0*Delta)) * &
            (A*(1.d0+B/(lambda*1.d-3)**2.d0))**2.d0 * &
            (1.d0/(lambda*1.d-3)**4.d0)

  end subroutine

end module



