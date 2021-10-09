
module photochem_input
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  use photochem_types, only : PhotoMechanism, PhotoSettings, PhotoRadTran, PhotoInitAtm
  implicit none
  private 
  integer,parameter :: real_kind = kind(1.0d0)
  integer, parameter :: err_len = 1024
  integer, parameter :: str_len = 1024

  public :: get_photomech, get_photorad, get_photoset, reaction_string, read_all_files
    
contains
  
  subroutine read_all_files(mechanism_file, settings_file, flux_file, atmosphere_txt, &
                            photomech, photoset, photorad, photoinit, err)
    
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    type(PhotoMechanism), intent(out) :: photomech
    type(PhotoSettings), intent(out) :: photoset
    type(PhotoRadTran), intent(out) :: photorad
    type(PhotoInitAtm), intent(out) :: photoinit
    character(len=err_len), intent(out) :: err
    
    ! first get SL and background species from settings
    call get_SL_and_background(settings_file, photoset, err)
    if (len(trim(err)) /= 0) return
    
    call get_photomech(mechanism_file, photoset, photomech, err)
    if (len(trim(err)) /= 0) return
    
    call get_photoset(settings_file, photomech, photoset, err)
    if (len(trim(err)) /= 0) return
    
    call get_photorad(photomech, photoset, photorad, err)
    if (len(trim(err)) /= 0) return
    
    ! stelar flux
    allocate(photorad%photon_flux(photorad%nw))
    call read_stellar_flux(flux_file, photorad%nw, photorad%wavl, photorad%photon_flux, err)
    if (len(trim(err)) /= 0) return
    
    ! initial atmosphere
    call read_atmosphere_file(atmosphere_txt, photomech, photoset, photoinit,err)
    if (len(trim(err)) /= 0) return
    
  end subroutine
  
  subroutine get_photomech(infile, photoset, photomech, err) 
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    type(PhotoSettings), intent(inout) :: photoset
    type(PhotoMechanism), intent(out) :: photomech
    character(len=err_len), intent(out) :: err
    
    character(error_length) :: error
    class (type_node), pointer :: root
    err = ''
    
    ! parse yaml file
    root => parse(infile,error=error)
    if (error /= '') then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        call get_rxmechanism(root, infile, photoset, photomech, err)
        if (len_trim(err) > 0) return
        call root%finalize()
        deallocate(root)
      class default
        err = "yaml file must have dictionaries at root level"
        return
    end select
     
  end subroutine
  
  subroutine get_rxmechanism(mapping, infile, photoset, photomech, err)
    class (type_dictionary), intent(in), pointer :: mapping
    character(len=*), intent(in) :: infile
    type(PhotoSettings), intent(inout) :: photoset
    type(PhotoMechanism), intent(out) :: photomech
    character(len=err_len), intent(out) :: err
    
    class (type_dictionary), pointer :: atoms, sat_params
    class (type_list), pointer :: species, reactions
    class (type_list), pointer :: particles
    type (type_error), pointer :: io_err
    class (type_list_item), pointer :: item, next
    class (type_dictionary), pointer :: dict
    class (type_key_value_pair), pointer :: key_value_pair

    ! temporary work variables
    character(len=str_len) :: tmpchar
    character(len=str_len) :: tmp
    character(len=20), allocatable :: eqr(:), eqp(:)
    integer :: i, ii, j, k, kk, l, ind(1), size_eqr, size_eqp
    logical :: reverse
    ! all_species causes a small memory leak. Not sure how to free the memory properly
    type(type_list) :: all_species, all_reactions ! will include particles

    err = ''

    atoms => mapping%get_dictionary('atoms',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    species => mapping%get_list('species',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    reactions => mapping%get_list('reactions',.true.,error = io_err) 
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    ! should i reverse reactions?
    photomech%reverse = mapping%get_logical('reverse-reactions',.true.,error = io_err)
    
    !!! atoms !!!
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
    
    !!! particles !!!
    ! get particles.
    particles => mapping%get_list('particles',.false.,error = io_err) 
    if (associated(particles)) then
      ! there are particles
      photomech%there_are_particles = .true.
      photomech%np = 0
      item => particles%first
      do while (associated(item))
        photomech%np = photomech%np + 1
        item => item%next
      enddo
      
      allocate(photomech%particle_names(photomech%np))
      allocate(photomech%particle_formation_method(photomech%np))
      allocate(photomech%particle_density(photomech%np))
      allocate(photomech%particle_sat_params(3,photomech%np))
      allocate(photomech%particle_gas_phase(photomech%np))
      allocate(photomech%particle_optical_prop(photomech%np))
      allocate(photomech%particle_optical_type(photomech%np))
      
      item => particles%first
      j = 1
      do while (associated(item))
        call all_species%append(item%node)
        select type (element => item%node)
        class is (type_dictionary)
          photomech%particle_names(j) = element%get_string("name",error = io_err) ! get name
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          
          tmpchar = element%get_string("formation",error = io_err) 
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          if (trim(tmpchar) == 'saturation') then
            photomech%particle_formation_method(j) = 1
          elseif (trim(tmpchar) == 'reaction') then
            photomech%particle_formation_method(j) = 2
          else
            err = "IOError: the only formation mechanism for particles is 'saturation'"
            return
          endif
          photomech%particle_density(j) = element%get_real("density",error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          photomech%particle_optical_prop(j) = element%get_string("optical-properties",error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          tmpchar = element%get_string("optical-type",error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          if (trim(tmpchar) == "mie") then
            photomech%particle_optical_type(j) = 0
          elseif  (trim(tmpchar) == "fractal") then
            err = "IOError: 'fractal' is not an optional optical type for "//trim(photomech%particle_names(j))
            return
          else
            err = "IOError: "//trim(tmpchar)//" is not an optional optical type for "//trim(photomech%particle_names(j))
            return
          endif
  
          if (photomech%particle_formation_method(j) == 1) then
            ! there should be saturation vapor pressure information
            sat_params => element%get_dictionary('saturation-parameters',.true.,error = io_err) 
            if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
            i = 0
            key_value_pair => sat_params%first
            do while (associated(key_value_pair))
              tmpchar = trim(key_value_pair%key)
              
              if (trim(tmpchar) == "A") then
                photomech%particle_sat_params(1,j) = sat_params%get_real(trim(tmpchar),error = io_err)
                if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
              elseif (trim(tmpchar) == "B") then
                photomech%particle_sat_params(2,j) = sat_params%get_real(trim(tmpchar),error = io_err)
                if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
              elseif (trim(tmpchar) == "C") then
                photomech%particle_sat_params(3,j) = sat_params%get_real(trim(tmpchar),error = io_err)
                if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
              else
                err = "Particle "//trim(photomech%particle_names(j))//" saturation parameters "//&
                      "can only be 'A', 'B', or 'C'"
                return
              endif                
              key_value_pair => key_value_pair%next
              i = i + 1
            enddo
            if (i /= 3) then
              err = "IOError: Missing or two many saturation parameters for "//trim(photomech%particle_names(j))
              return 
            endif
            
            ! gas phase
            photomech%particle_gas_phase(j) = element%get_string("gas-phase",error = io_err) 
            if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          elseif (photomech%particle_formation_method(j) == 2) then
            ! add the reaction to the list of reactions
            call all_reactions%append(item%node)
          endif
        
        class default
          err = "IOError: Problem with particle number "//char(j)//"  in the input file"
          return
        end select
        item => item%next
        j = j + 1
      enddo
    else ! there are no particles
      photomech%there_are_particles = .false.
      photomech%np = 0
    endif
    
    ! for now number particle equations will be the same 
    ! as number of particles
    photomech%npq = photomech%np
    
    !!! done with particles !!!
    
    !!! species !!!
    photomech%ng = 0 ! count number of gas phase species
    item => species%first
    do while (associated(item))
      item => item%next
      photomech%ng = photomech%ng + 1
    enddo
    
    ! get number of sl from photoset
    photomech%nsl = photoset%nsl
    
    if (photoset%back_gas) then
      photomech%nll = photomech%ng - photomech%nsl - 1 ! minus 1 for background
    else
      photomech%nll = photomech%ng - photomech%nsl
    endif
    
    photomech%ng_1 = photomech%npq + 1 ! the long lived gas index
    ! photomech%nq is the last ll gas index
    
    ! now we now nq, the number of PDEs
    photomech%nq = photomech%npq + photomech%nll
    photoset%nq = photomech%nq
    
    ! we also now nsp, the index of the backgorund gas 
    photomech%nsp = photomech%npq + photomech%ng
    
    ! species_mass, species_composition, and species_names
    ! will include the particles, thus we allocate nsp
    allocate(photomech%species_mass(photomech%nsp))
    allocate(photomech%species_composition(photomech%natoms,photomech%nsp+2))
    photomech%species_composition = 0
    allocate(photomech%species_names(photomech%nsp+2))
    photomech%species_names(photomech%nsp+1) = "hv" ! always add these guys
    photomech%species_names(photomech%nsp+2) = "M"
    ! we will not include particles in thermodynamic data.
    if (photomech%reverse) then
      allocate(photomech%thermo_data(7,2,photomech%ng))
      allocate(photomech%thermo_temps(3,photomech%ng))
    endif
    
    ! Append the species to the end of a list
    ! which has particles in the beginning
    item => species%first
    do while (associated(item))
      call all_species%append(item%node)
      item => item%next
    enddo

    ! Loop through particles and gases
    kk = photomech%ng_1
    l = 1
    ii = 1 ! overall counter
    item => all_species%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmpchar = trim(element%get_string("name",error = io_err)) ! get name
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        
        if (ii < photomech%ng_1) then
          ! we are dealing with particles
          j = ii
        else
          ! we are dealing with gases
          ind = findloc(photoset%SL_names,tmpchar)
          if (ind(1) /= 0) then ! short lived species
            j = photoset%nq + l 
            l = l + 1
          elseif (tmpchar == photoset%back_gas_name) then ! background gas
            j = photomech%nsp
          else ! long lived species
            j = kk
            kk = kk + 1
          endif
        endif
                  
        photomech%species_names(j) = tmpchar
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
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        if (photomech%reverse .and. (ii >= photomech%ng_1)) then
          call get_thermodata(element,photomech%species_names(j), infile,photomech%thermo_temps(:,j-photomech%npq), &
                              photomech%thermo_data(:,:,j-photomech%npq), err) ! get thermodynamic data
          if (len_trim(err) > 0) return
        endif
      class default
        err = "IOError: Problem with species number "//char(j)//"  in the input file"
        return
      end select
      ii = ii + 1
      item => item%next
    enddo
    
    if (l-1 /= photoset%nsl) then
      err = 'IOError: One of the short lived species is not in the file '//trim(infile)
      return
    endif
    if (photoset%back_gas) then
      ind = findloc(photomech%species_names,photoset%back_gas_name)
      if (ind(1) == 0) then
        err = 'IOError: The specified background gas is not in '//trim(infile)
        return
      endif
    endif
    item => all_species%first
    do while (associated(item))
       next => item%next
       deallocate(item)
       item => next
    end do
    nullify(all_species%first)
    !!! done with species !!!
    
    if (photomech%there_are_particles) then
      ! get indexes of gas phase condensing species
      allocate(photomech%particle_gas_phase_ind(photomech%np))
      do i = 1,photomech%np
        if (photomech%particle_formation_method(i) == 1) then
          ! if a condensing molecule
          ind = findloc(photomech%species_names,photomech%particle_gas_phase(i))
          if (ind(1) /= 0) then
            photomech%particle_gas_phase_ind(i) = ind(1)
          else
            err = "IOError: particle "//trim(photomech%particle_names(i))// &
                  " can not be made from "//trim(photomech%particle_gas_phase(i))// &
                  " because "//trim(photomech%particle_gas_phase(i))//" is not a gas"// &
                  " in the model."
            return
          endif
        endif
      enddo
    endif
    
    !!! reactions !!!
    item => reactions%first
    do while (associated(item))
      call all_reactions%append(item%node)
      item => item%next
    enddo

    photomech%nrF = 0 ! count forward reactions
    item => all_reactions%first
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
    item => all_reactions%first
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
        
        call compare_rxtype_string(tmp, eqr, eqp, reverse, photomech%rxtypes(j),err)
        if (len_trim(err) > 0) return
        
        size_eqp = size(eqp)
        size_eqr = size(eqr)
        if ((photomech%rxtypes(j) == 2) .or. (photomech%rxtypes(j) == 3)) then ! if threebody or falloff
          size_eqr = size_eqr - 1 ! remove the M
          size_eqp = size_eqp - 1
        endif
        
        if (reverse) then
          if (.not.photomech%reverse) then
            err = 'IOError: reaction file '//trim(infile)//' contains reverse reaction '//tmp// &
                  ', which is incompatible with "reverse-reactions: false"'
            return
          endif
          photomech%nrR = photomech%nrR + 1
          if (size_eqr > photomech%max_num_products) photomech%max_num_products = size_eqr
          if (size_eqp > photomech%max_num_reactants) photomech%max_num_reactants = size_eqp
        endif
        if (size_eqr > photomech%max_num_reactants) photomech%max_num_reactants = size_eqr
        if (size_eqp > photomech%max_num_products) photomech%max_num_products = size_eqp
        
        call count_efficiencies(element, photomech%num_efficient(j))
        
        if (photomech%rxtypes(j) == 0) then ! if photolysis reaction
          photomech%kj = photomech%kj + 1
        endif
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
    photomech%reactants_sp_inds = huge(1)
    photomech%products_sp_inds = huge(1)
    if (photomech%reverse) then
      allocate(photomech%reverse_info(photomech%nrT))
      photomech%reverse_info = 0 ! initialize
    endif
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
    item => all_reactions%first
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

          kk = photomech%nreactants(i)
          l = photomech%nreactants(j)
          
          photomech%reactants_sp_inds(1:kk,i) = photomech%products_sp_inds(1:kk,j)
          photomech%products_sp_inds(1:l,i) = photomech%reactants_sp_inds(1:l,j)
          k = k + 1
        endif
      class default
        err = "IOError: Problem with reaction number "//char(j)//" in the input file."
        return
      end select
      item => item%next
      j = j + 1
    enddo
    
    item => all_reactions%first
    do while (associated(item))
       next => item%next
       deallocate(item)
       item => next
    end do
    nullify(all_reactions%first)
    

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
    
    call check_for_duplicates(photomech,err)
    if (len(trim(err)) > 0) return
    !!! end reactions !!!
    
    
    !!! henrys law !!!
    call get_henry(photomech, err)
    if (len(trim(err)) > 0) return
    !!! end henrys law !!!
    
  end subroutine
  
  subroutine get_henry(photomech, err)
    use yaml, only : parse, error_length
    use photochem_vars, only: data_dir
    type(PhotoMechanism), intent(inout) :: photomech
    character(len=*), intent(out) :: err
    
    character(error_length) :: error
    class (type_node), pointer :: root
    type (type_error), pointer :: io_err
    class (type_list_item), pointer :: item
    integer :: j, ind(1), i
    
    character(len=15), allocatable :: henry_names(:)
    real(real_kind), allocatable :: henry_data(:,:)
    
    err = ''

    ! parse yaml file
    root => parse(trim(data_dir)//"/henry/henry.yaml",error=error)
    if (error /= '') then
      err = trim(error)
      return
    end if
    select type (root)
    class is (type_list)
      
      j = 0
      item => root%first
      do while(associated(item))
        j = j + 1
        item => item%next
      enddo
      
      allocate(henry_names(j))
      allocate(henry_data(2,j))
      j = 1
      item => root%first
      do while(associated(item))
        select type (element => item%node)
        class is (type_dictionary)
          henry_names(j) = element%get_string('name',error = io_err)
          if (associated(io_err)) then; err = trim(io_err%message); return; endif
          henry_data(1,j) = element%get_real('A',error = io_err)
          if (associated(io_err)) then; err = trim(io_err%message); return; endif
          henry_data(2,j) = element%get_real('B',error = io_err)
          if (associated(io_err)) then; err = trim(io_err%message); return; endif
        
        j = j + 1
        item => item%next
        end select
      enddo
      
      
      call root%finalize()
      deallocate(root)
    class default
      err = "yaml file must have dictionaries at root level"
      return
    end select
    
    allocate(photomech%henry_data(2,photomech%nsp))
    photomech%henry_data = 0.d0
    do j = 1,size(henry_names)
      ind = findloc(photomech%species_names,henry_names(j))
      if (ind(1) /= 0) then
        i = ind(1)
        photomech%henry_data(:,i) = henry_data(:,j)
      endif
    enddo
    ! set particle solubility to super high number
    photomech%henry_data(1,1:photomech%npq) = 1.d11
    
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
    
    integer i, ii, nr, j, np, k, matches, matches_old, kk, jj, jjj
    logical l, m
    
    integer, allocatable :: tmp_arr(:), tmp_arr_save(:)
    
    err = ''
    
    if (photomech%max_num_reactants > photomech%max_num_products) then
      allocate(tmp_arr(photomech%max_num_reactants))
      allocate(tmp_arr_save(photomech%max_num_reactants))
    else
      allocate(tmp_arr(photomech%max_num_products))
      allocate(tmp_arr_save(photomech%max_num_products))
    endif
    
    do i = 1,photomech%nrT-1
      do ii = i+1,photomech%nrT
        
        ! check the same num of reactants and products
        nr = photomech%nreactants(i)
        np = photomech%nproducts(i)
        if (nr == photomech%nreactants(ii) .and. np == photomech%nproducts(i)) then
          
          !!! check if all the reactants are the same
          tmp_arr(1:nr) = photomech%reactants_sp_inds(1:nr,ii)
          matches = 0
          matches_old = 0
          do j = 1,nr
            do k = 1,nr-matches
              if (photomech%reactants_sp_inds(j,i) == tmp_arr(k)) then
                ! we have found a match. We new delete the match
                ! from tmp_arr and compare next element of
                ! photomech%reactants_sp_inds(j,i)
                jjj = nr-matches
                jj = 1
                do kk = 1,jjj
                  if (kk == k) then
                    cycle
                  else
                    tmp_arr_save(jj) = tmp_arr(kk)
                    jj = jj + 1
                  endif
                enddo
                tmp_arr(1:jjj-1) = tmp_arr_save(1:jjj-1)
                
                matches_old = matches
                matches = matches + 1
                exit
              endif
            enddo
            if (matches == matches_old) then
              ! this means there wasn't a match for one species
              ! so we can exit
              exit
            endif
          enddo
          
          if (nr == matches) then
            m = .true.
          else
            m = .false.
          endif
          
          !!! check if the products are the same
          tmp_arr(1:np) = photomech%products_sp_inds(1:np,ii)
          matches = 0
          matches_old = 0
          do j = 1,np
            do k = 1,np-matches
              if (photomech%products_sp_inds(j,i) == tmp_arr(k)) then
                
                jjj = np-matches
                jj = 1
                do kk = 1,jjj
                  if (kk == k) then
                    cycle
                  else
                    tmp_arr_save(jj) = tmp_arr(kk)
                    jj = jj + 1
                  endif
                enddo
                tmp_arr(1:jjj-1) = tmp_arr_save(1:jjj-1)
                
                matches_old = matches
                matches = matches + 1
                exit
              endif
            enddo
            if (matches == matches_old) then
              exit
            endif
          enddo
          
          if (np == matches) then
            l = .true.
          else
            l = .false.
          endif
          
          if (m .and. l) then
            err = "IOError: This reaction is a duplicate: "
            call reaction_string(photomech,i,rxstring)
            err(len_trim(err)+2:) = rxstring
            return
          endif
          
        endif
      enddo
    enddo
    
  end subroutine
  
  subroutine reaction_string(photomech,rxn,rxstring)
    type(PhotoMechanism), intent(in) :: photomech
    integer, intent(in) :: rxn
    character(len=:), allocatable, intent(out) :: rxstring
    integer j, k, i
    rxstring = ''
    if (rxn > photomech%nrF) then
      i = photomech%reverse_info(rxn)
    else
      i = rxn
    endif
    do j = 1,photomech%nreactants(rxn)-1
      k = photomech%reactants_sp_inds(j,rxn)
      rxstring = rxstring //(trim(photomech%species_names(k))//' + ')
    enddo
    
    k = photomech%reactants_sp_inds(photomech%nreactants(rxn),rxn)
    rxstring = rxstring // trim(photomech%species_names(k))//' => '
    
    if ((photomech%rxtypes(i) == 2) .or.(photomech%rxtypes(i) == 3)) then
      rxstring = rxstring(1:len(rxstring)-4) //(' + M'//' => ')
    endif
    
    do j = 1,photomech%nproducts(rxn)-1
      k = photomech%products_sp_inds(j,rxn)
      rxstring = rxstring // trim(photomech%species_names(k))//' + '
    enddo
    k = photomech%products_sp_inds(photomech%nproducts(rxn),rxn)
    rxstring = rxstring // trim(photomech%species_names(k))
    
    if ((photomech%rxtypes(i) == 2) .or.(photomech%rxtypes(i) == 3)) then
      rxstring = rxstring //' + M'
    endif
  end subroutine
  
  subroutine compare_rxtype_string(tmp, eqr, eqp, reverse, rxtype_int, err)
    character(len=*), intent(in) :: tmp
    character(len=20), allocatable, intent(in) :: eqr(:), eqp(:)
    logical, intent(in) :: reverse
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
        if ((trim(eqr(i)) == 'M').and.j) jj = .true.
        if (trim(eqr(i)) == 'M') j = .true.
        if (trim(eqr(i)) == 'hv') m = .true.
      enddo
      do i = 1,size(eqp)
        if ((trim(eqp(i)) == 'M').and.k) kk = .true.
        if (trim(eqp(i)) == 'M') k = .true.
        if (trim(eqp(i)) == 'hv') l = .true.
      enddo
      if (trim(eqr(size(eqr))) /= 'M') then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' must have "M" as the last reactant'
        return
      endif
      if (trim(eqp(size(eqp))) /= 'M') then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' must have "M" as the last product'
        return
      endif
      if ((j).and.(k)) then
        ! good
      else
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' must have "M" on both sides'
        return
      endif
      if ((jj).or.(kk)) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can only have one "M" on either side'
        return
      endif
      if ((m).or.(l)) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can not contain "hv". Only photolysis reactions can.'
        return
      endif
    elseif (trim(rxtype) == 'elementary') then
      do i = 1,size(eqr)
        if (trim(eqr(i)) == 'M') j = .true.
        if (trim(eqr(i)) == 'hv') m = .true.
      enddo
      do i = 1,size(eqp)
        if (trim(eqp(i)) == 'M') k = .true.
        if (trim(eqp(i)) == 'hv') l = .true.
      enddo
      if ((j).or.(k)) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can not contain "M".'
        return
      endif
      if ((m).or.(l)) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can not contain "hv". Only photolysis reactions can.'
        return
      endif
    elseif (trim(rxtype) == 'photolysis') then
      if (reverse) then
        err = 'IOError: Photolysis reaction '//trim(tmp)//' can not be reversed.'
      endif
      
      do i = 1,size(eqr)
        if (trim(eqr(i)) == 'M') j = .true.
        if ((trim(eqr(i)) == 'hv').and.(m)) jj = .true.
        if (trim(eqr(i)) == 'hv') m = .true.
      enddo
      do i = 1,size(eqp)
        if (trim(eqp(i)) == 'M') k = .true.
        if (trim(eqp(i)) == 'hv') l = .true.
      enddo
      if ((j).or.(k)) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can not contain "M".'
        return
      endif
      if (jj) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can only have one "hv" on the left side.'
        return
      endif
      if ((m).and..not.(l)) then
        ! good
      else
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
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
    character(len=str_len) :: rxtype_str
    err = ''
    rateparam = 0.d0
    ! no error possible
    rxtype_str = reaction%get_string("type"," ",error = io_err) 
    if (rxtype_str == 'photolysis') then
      rxtype = 0
    elseif ((rxtype_str == 'elementary') .or. (len_trim(rxtype_str) == 0)) then
      rxtype = 1
      rxtype_str = 'elementary'
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
    character(len=*), intent(in) :: instring
    logical, intent(out) :: reverse
    character(len=20), allocatable, intent(out) :: eqr(:), eqp(:)
    character(len=err_len), intent(out) :: err
    
    character(len=:), allocatable :: string2, string3
    character(len=:), allocatable :: eqr1, eqp1, eqr2, eqp2
    integer :: i, nreac1, nprod1, nreac2, nprod2, start
    logical :: switch
    
    string2 = ''
    string3 = ''
    do i = 1,len_trim(instring)
      if (instring(i:i) /= '(' .and. instring(i:i) /= ')') then
        string2 = string2//instring(i:i)
        if (instring(i:i) == '+') then
          string3 = string3//' '
        else
          string3 = string3//instring(i:i)
        endif
      endif
    enddo
    
    if (index(instring, "<=>") /= 0) then
      i = index(string2, "<=>")
      eqr1 = string2(1:i-1)
      eqp1 = string2(i+3:)
      i = index(string3, "<=>")
      eqr2 = string3(1:i-1)
      eqp2 = string3(i+3:)
      reverse = .true.
    elseif (index(instring, " =>") /= 0) then
      i = index(string2, " =>")
      eqr1 = string2(1:i-1)
      eqp1 = string2(i+3:)
      i = index(string3, " =>")
      eqr2 = string3(1:i-1)
      eqp2 = string3(i+3:)
      reverse = .false.
    else
      err = "IOError: Invalid reaction arrow in reaction "//instring// &
            '. Note, forward reactions must have a space before the arrow, like " =>"'
      return
    endif
    
    ! trim whitespace
    call trim_whitespace(eqr1)
    call trim_whitespace(eqp1)
    call trim_whitespace(eqp2)
    call trim_whitespace(eqr2)
    
    ! count number of + signs
    nreac1 = 1
    do i = 1,len(eqr1)
      if (eqr1(i:i) == '+') nreac1 = nreac1 + 1
    enddo
    nprod1 = 1
    do i = 1,len(eqp1)
      if (eqp1(i:i) == '+') nprod1 = nprod1 + 1
    enddo
    
    ! count number of gaps
    nreac2 = 1
    switch = .true.
    do i = 1,len(eqr2)
      if (eqr2(i:i) == ' ' .and. switch) then
        nreac2 = nreac2 + 1
        switch = .false.
      endif
      if (eqr2(i:i) /= ' ' .and. .not.switch) switch = .true.
    enddo
    nprod2 = 1
    switch = .true.
    do i = 1,len(eqp2)
      if (eqp2(i:i) == ' ' .and. switch) then
        nprod2 = nprod2 + 1
        switch = .false.
      endif
      if (eqp2(i:i) /= ' ' .and. .not.switch) switch = .true.
    enddo
    
    if (nreac1 /= nreac2 .or. nprod1 /= nprod2) then
      err = 'IOError: Missing or too many "+" sign(s) in reaction '//instring
      return
    endif
    
    if (allocated(eqr)) then
      deallocate(eqr, eqp)
    endif
    allocate(eqr(nreac1), eqp(nprod1))
    eqr = ''
    eqp = ''
    
    nreac2 = 1
    start = 1
    switch = .true.
    do i = 1,len(eqr2)
      if (eqr2(i:i) == ' ' .and. switch) then
        eqr(nreac2) = eqr2(start:i-1)
        nreac2 = nreac2 + 1
        switch = .false.
      endif
      if (eqr2(i:i) /= ' ' .and. .not.switch) then
        start = i
        switch = .true.
      endif
    enddo
    eqr(nreac2) = eqr2(start:)
    
    nprod2 = 1
    start = 1
    switch = .true.
    do i = 1,len(eqp2)
      if (eqp2(i:i) == ' ' .and. switch) then
        eqp(nprod2) = eqp2(start:i-1)
        nprod2 = nprod2 + 1
        switch = .false.
      endif
      if (eqp2(i:i) /= ' ' .and. .not.switch) then
        start = i
        switch = .true.
      endif
    enddo
    eqp(nprod2) = eqp2(start:)
    
  end subroutine
  
  subroutine trim_whitespace(str)
    character(len=:), allocatable, intent(inout) :: str
    
    integer :: i
    do i = 1,len(str)
      if (str(i:i) /= ' ') exit
    enddo
    str= str(i:)
    do i = len(str),1,-1
      if (str(i:i) /= ' ') exit
    enddo
    str = str(1:i)
  end subroutine
    
  subroutine get_reaction_chars(instring, max_num_react, max_num_prod, numr, nump, &
                                    outreact, outprod, reverse, err)
    character(len=*), intent(in) :: instring
    integer, intent(in) :: max_num_react, max_num_prod
    
    integer, intent(out) :: numr, nump
    character(len=*), intent(out) :: outreact(max_num_react), outprod(max_num_prod)
    logical, intent(out) :: reverse
    character(len=err_len), intent(out) :: err
    
    character(len=20), allocatable :: eqr(:), eqp(:)
    integer :: i
    
    call parse_reaction(instring, reverse, eqr, eqp, err)
    if (len_trim(err) > 0) return
    
    numr = size(eqr)
    nump = size(eqp)

    if (eqr(size(eqr)) == 'M') then
      numr = numr - 1
      nump = nump - 1
    endif
    
    outreact = ''
    outprod = ''
    do i=1,numr
      outreact(i) = eqr(i)
    enddo
    do i=1,nump
      outprod(i) = eqp(i)
    enddo
  end subroutine
    
  subroutine get_reaction_sp_nums(reaction, max_num_react, max_num_prod, reacts, prods, &
                                  species_names, species_composition, natoms, nsp, &
                                  react_sp_nums, prod_sp_nums, err)
    character(len=*), intent(in) :: reaction
    integer, intent(in) :: max_num_react, max_num_prod
    character(len=*), intent(in) :: reacts(max_num_react)
    character(len=*), intent(in) :: prods(max_num_prod)
    integer, intent(in) :: nsp, natoms
    character(len=*), intent(in) :: species_names(nsp+2)
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
               "Species "//trim(reacts(i))//" in reaction "//trim(reaction)// &
               " is not in the list of species."
        return
      endif
    enddo
    
    do i = 1,max_num_prod
      ind = findloc(species_names,prods(i))
      prod_sp_nums(i) = ind(1)
      if ((prods(i) /= '') .and. (ind(1) == 0)) then
        err = "IOError: "// & 
               "Species "//trim(prods(i))//" in reaction "//trim(reaction)// &
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
             'Bad mass balance in reaction "'//trim(reaction)// &
             '". You could have messed up how many atoms one of the species has.'
      return
    endif
  end subroutine
  
  subroutine get_SL_and_background(infile, photoset, err)
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    type(PhotoSettings), intent(inout) :: photoset
    character(len=err_len), intent(out) :: err
  
    character(error_length) :: error
    class (type_node), pointer :: root
    class (type_dictionary), pointer :: tmp1
    type (type_error), pointer :: io_err
    class (type_list), pointer :: bcs
    class (type_list_item), pointer :: item
    character(len=str_len) :: spec_type
    integer :: i, j
    err = ''
    
    root => parse(infile,error=error)
    if (error /= '') then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        
        ! get background species
        tmp1 => root%get_dictionary('planet',.true.,error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        photoset%back_gas = tmp1%get_logical('use-background-gas',error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        if (photoset%back_gas) then
          photoset%back_gas_name = tmp1%get_string('background-gas',error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        else
          photoset%back_gas_ind = -1
          err = "Currently, the model requires there to be a background gas."// &
                " You must set 'use-background-gas: true'"
          return
        endif
        
        bcs => root%get_list('boundary-conditions',.true.,error = io_err)
        
        photoset%nsl = 0
        item => bcs%first
        do while (associated(item))
          select type (element => item%node)
          class is (type_dictionary)
            spec_type = element%get_string('type','long lived',error = io_err)
            if (spec_type == 'short lived') then
              photoset%nsl = photoset%nsl + 1
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
        enddo
        allocate(photoset%SL_names(photoset%nsl))
        
        photoset%nsl = 0
        item => bcs%first
        do while (associated(item))
          select type (element => item%node)
          class is (type_dictionary)
            spec_type = element%get_string('type','long lived',error = io_err)
            if (spec_type == 'short lived') then
              photoset%nsl = photoset%nsl + 1
              photoset%SL_names(photoset%nsl) = element%get_string('name',error = io_err)
              if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
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
        enddo

        call root%finalize()
        deallocate(root)
      class default
        err = trim(infile)//" file must have dictionaries at root level"
        return
    end select
    
    ! check for duplicates
    do i = 1,photoset%nsl-1
      do j = i+1,photoset%nsl
        if (photoset%SL_names(i) == photoset%SL_names(j)) then
          err = 'IOError: Short lived species '//trim(photoset%SL_names(i))// &
                ' is a duplicate.'
          return
        endif
      enddo
    enddo
  
  end subroutine
  
  subroutine get_photoset(infile, photomech, photoset, err)
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoSettings), intent(inout) :: photoset
    character(len=err_len), intent(out) :: err
  
    character(error_length) :: error
    class (type_node), pointer :: root
    err = ''
    
    root => parse(infile,error=error)
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

  end subroutine
  
  subroutine unpack_settings(mapping, infile, photomech, photoset, err)
    class (type_dictionary), intent(in), pointer :: mapping
    character(len=*), intent(in) :: infile
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoSettings), intent(inout) :: photoset
    character(len=err_len), intent(out) :: err
    
    class (type_dictionary), pointer :: tmp1, tmp2, tmp3
    class (type_list), pointer :: bcs, particles
    class (type_list_item), pointer :: item
    type (type_error), pointer :: io_err
    
    character(len=str_len) :: spec_type
    integer :: j, i, ind(1), io
    character(len=str_len), allocatable :: dups(:)
    character(len=30) :: temp_char
    integer :: default_lowerboundcond
    logical, allocatable :: particle_checklist(:)
    
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
      photoset%nw = tmp1%get_integer('number-of-bins',error = io_err)
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
    ! scale factor for photon flux. Its optional
    photoset%photon_scale_factor = tmp1%get_real('photon-scale-factor', 1.d0,error = io_err)
    
    ! atmosphere grid
    tmp1 => mapping%get_dictionary('atmosphere-grid',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photoset%bottom_atmos = tmp1%get_real('bottom',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photoset%top_atmos = tmp1%get_real('top',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photoset%nz = tmp1%get_integer('number-of-layers',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    ! Planet
    tmp1 => mapping%get_dictionary('planet',.true.,error = io_err)
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
    photoset%diurnal_fac = tmp1%get_real('diurnal-averaging-factor',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photoset%diurnal_fac < 0.d0 .or. photoset%diurnal_fac > 1.d0) then
      err = 'IOError: diurnal-averaging-factor must be between 0 and 1.'
      return
    endif
    
    photoset%solar_zenith = tmp1%get_real('solar-zenith-angle',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photoset%solar_zenith < 0.d0 .or. photoset%solar_zenith > 90.d0) then
      err = 'IOError: solar zenith must be between 0 and 90.'
      return
    endif
    
    ! H2 escape
    photoset%diff_H_escape = tmp1%get_logical('diff-lim-hydrogen-escape',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    ind = findloc(photomech%species_names,'H2')
    photoset%LH2 = ind(1)
    if (ind(1) == 0 .and. photoset%diff_H_escape) then
      err = 'IOError: H2 must be a species if diff-lim-hydrogen-escape = True.'
      return
    endif
    ind = findloc(photomech%species_names,'H')
    photoset%LH = ind(1)
    if (ind(1) == 0 .and. photoset%diff_H_escape) then
      err = 'IOError: H must be a species if diff-lim-hydrogen-escape = True.'
      return
    endif
    
    ! default lower boundary
    temp_char = tmp1%get_string('default-lower-boundary',"deposition velocity",error = io_err)
    if (trim(temp_char) == 'deposition velocity') then
      default_lowerboundcond = 0
    elseif (trim(temp_char) == 'Moses') then
      default_lowerboundcond = -1
    else
      err = "IOError: Only 'deposition velocity' or 'Moses' can be default boundary conditions."
      return
    endif
    
    tmp2 => tmp1%get_dictionary('water',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photoset%fix_water_in_trop = tmp2%get_logical('fix-water-in-troposphere',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    ind = findloc(photomech%species_names,'H2O')
    photoset%LH2O = ind(1)
    if (ind(1) == 0 .and. photoset%fix_water_in_trop) then
      err = 'IOError: H2O must be a species if water-saturated-troposhere = True.'
      return
    endif
    if (photoset%fix_water_in_trop) then  
      photoset%trop_alt = tmp2%get_real('tropopause-altitude',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      if ((photoset%trop_alt < photoset%bottom_atmos) .or. &
          (photoset%trop_alt > photoset%top_atmos)) then
          err = 'IOError: tropopause-altitude must be between the top and bottom of the atmosphere'
          return
      endif
      
      temp_char = tmp2%get_string('relative-humidity',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
      read(temp_char,*,iostat = io) photoset%relative_humidity
      
      if (io /= 0) then
        ! it isn't a float
        if (trim(temp_char) == "manabe") then
          photoset%use_manabe = .true.
        else
          err = '"relative-humidity" can only be a number between 0 and 1, or "manabe". See '//trim(infile)
          return 
        endif
      else
        photoset%use_manabe = .false.
      endif
      
      photoset%gas_rainout = tmp2%get_logical('gas-rainout',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      photoset%stratospheric_cond = tmp2%get_logical('stratospheric-condensation',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      if (photoset%stratospheric_cond) then
        photoset%relative_humidity_cold_trap = tmp2%get_real('cold-trap-relative-humitity',error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        
        tmp3 => tmp2%get_dictionary('condensation-rate',.true.,error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        photoset%H2O_condensation_rate(1) = tmp3%get_real('A',error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        photoset%H2O_condensation_rate(2) = tmp3%get_real('rh0',error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        if (photoset%H2O_condensation_rate(2) <= 1) then
          err = 'IOError: Rate constant "rh0" for H2O condensation must be > 1. See '//trim(infile)
          return
        endif
      endif
    else
      photoset%gas_rainout = .false.
        
    endif
    
    ! particle parameters
    if (photomech%there_are_particles) then
      particles => mapping%get_list('particles',.true.,error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
      allocate(particle_checklist(photomech%np))
      allocate(photoset%condensation_rate(2,photomech%np))
      particle_checklist = .false.
      
      item => particles%first
      do while (associated(item))
        select type (element => item%node)
        class is (type_dictionary)
          temp_char = element%get_string('name',error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          ind = findloc(photomech%particle_names,trim(temp_char))
          if (particle_checklist(ind(1))) then
            err = "IOError: particle "//trim(temp_char)//" in the settings"// &
                  " file is listed more than once"
            return
          endif
          if (ind(1) == 0) then
            err = "IOError: particle "//trim(temp_char)//" in the settings"// &
                  " file isn't in the list of particles in the reaction mechanism file"
            return
          else
            particle_checklist(ind(1)) = .true.
          endif
          
          tmp1 => element%get_dictionary('condensation-rate',.true.,error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          
          photoset%condensation_rate(1,ind(1)) = tmp1%get_real('A',error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          photoset%condensation_rate(2,ind(1)) = tmp1%get_real('rh0',error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          
        class default
          err = "IOError: Particle settings must be a list of dictionaries."
          return
        end select
        item => item%next
      end do
      
      do i = 1,photomech%np
        if (photomech%particle_formation_method(i) == 1 .and. .not. particle_checklist(i)) then
          err = 'IOError: Particle '//trim(photomech%particle_names(i))// &
                ' does not have any condensation rate data in the file '//trim(infile)
          return
        endif
      enddo
        
    endif
      
    ! boundary conditions and species types
    bcs => mapping%get_list('boundary-conditions',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    if (photoset%back_gas) then
      ind = findloc(photomech%species_names,trim(photoset%back_gas_name))
      photoset%back_gas_ind = ind(1)
      photoset%back_gas_mu = photomech%species_mass(ind(1))
    endif
    
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
    photoset%lowerboundcond = default_lowerboundcond
    photoset%lower_vdep = 0.d0
    photoset%upperboundcond = 0
    photoset%upper_veff = 0.d0
      
    ! determine number of short lived species. 
    allocate(dups(photomech%nsp))
    j = 1
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
        ! make sure it isn't water
        if (photoset%fix_water_in_trop .and. dups(j) == "H2O") then
          err = "IOError: H2O can not have a specified boundary condition"// &
                " if water-saturated-troposphere = true in the settings file."
          return
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
          ! nothing
          
        elseif (spec_type == 'long lived') then
          i = ind(1)
          ! get boundary condition
          call get_boundaryconds(element, dups(j), infile, &
                                 photoset%lowerboundcond(i), photoset%lower_vdep(i), &
                                 photoset%lower_flux(i), photoset%lower_dist_height(i), &
                                 photoset%lower_fix_mr(i), &
                                 photoset%upperboundcond(i), photoset%upper_veff(i), &
                                 photoset%upper_flux(i), err)
          if (len_trim(err) /= 0) return
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
    
    ! Make sure that upper boundary condition for H and H2 are
    ! effusion velocities, if diffusion limited escape
    if (photoset%diff_H_escape) then
      if (photoset%back_gas_name /= "H2") then
        if (photoset%upperboundcond(photoset%LH2) /= 0) then
          err = "IOError: H2 must have a have a effusion velocity upper boundary"// &
                " if diff-lim-hydrogen-escape = True"
          return
        endif
      endif
      if (photoset%upperboundcond(photoset%LH) /= 0) then
        err = "IOError: H must have a have a effusion velocity upper boundary"// &
              " if diff-lim-hydrogen-escape = True"
        return
      endif
    endif

    ! check for SL nonlinearities
    call check_sl(photomech, photoset, err)
    if (len_trim(err) /= 0) return

    
  end subroutine
  
  subroutine check_sl(photomech, photoset, err)
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoSettings), intent(in) :: photoset
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, l, k, kk, m, mm, n, nn, ind(1), counter
    character(len=:), allocatable :: reaction
    err = ''
    
    do i = 1, photoset%nsl
      j = photoset%nq + i
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
          do n = photoset%nq+1,photoset%nq + photoset%nsl
            if (n == photomech%reactants_sp_inds(mm,kk)) then
              err = 'IOError: Reaction "'//reaction//'" has short-lived species as reactants'// &
              ' and products. This is not allowed. Change one or both of the species to long-lived.'
              return
            endif
          enddo
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
          do nn = photoset%nq+1,photoset%nq + photoset%nsl
            if ((nn == n) .and. (n == j)) then
              counter = counter + 1
              if (counter > 1) then
                err = 'IOError: Reaction "'//reaction//'" short lived species react'// &
                ' with themselves. This is not allowed. Change the species to long lived.'
                return
              endif
            elseif ((nn == n) .and. (n /= j)) then
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
    
    
    ! deposition velocity = vdep
    ! mixing-ratio = mixing-ratio
    ! flux = flux
    ! deposition velocity + distributed flux = vdep + flux + distributed-height
    if (bctype == "deposition velocity") then
      lowercond = 0
      Lvdep = tmpdict%get_real("vdep",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      Lflux = 0.d0
      LdistH = 0.d0
      Lmr = 0.d0 
    elseif (bctype == "mixing-ratio") then
      lowercond = 1
      Lvdep = 0.d0
      Lflux = 0.d0
      LdistH = 0.d0
      Lmr = tmpdict%get_real("mixing-ratio",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    elseif (bctype == "flux") then
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
    elseif (bctype == "Moses") then
      lowercond = -1
    else
      err = 'IOError: "'//trim(bctype)//'" is not a valid lower boundary condition for '//trim(molecule_name)
      return
    endif
    
    tmpdict => molecule%get_dictionary("upper-boundary",.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    bctype = tmpdict%get_string("type",error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    
    ! effusion velocity = veff
    ! flux = flux
    if (bctype == "effusion velocity") then
      uppercond = 0
      Uveff = tmpdict%get_real("veff",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      Uflux = 0.d0
    elseif (bctype == "flux") then
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
    
    if (photomech%there_are_particles) then
      call get_aerosol_xs(photomech, photorad, err)
      if (len_trim(err) /= 0) return
    endif
    
  end subroutine
  
  subroutine get_aerosol_xs(photomech, photorad, err)
    use photochem_mie, only: read_mie_data_file
    use photochem_vars, only: data_dir
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoRadTran), intent(inout) :: photorad
    character(len=err_len), intent(out) :: err
    
    integer :: nrad
    integer, parameter :: nrad_fixed = 50
    real(real_kind), allocatable :: radii(:)
    real(real_kind), allocatable :: w0(:,:), qext(:,:), g(:,:)
    character(len=:), allocatable :: xsroot
    character(len=:), allocatable :: filename
    integer :: i
    
    xsroot = trim(data_dir)//"/aerosol_xsections/"
    
    allocate(photorad%radii_file(nrad_fixed,photomech%np))
    allocate(photorad%w0_file(nrad_fixed, photomech%np, photorad%nw))
    allocate(photorad%qext_file(nrad_fixed, photomech%np, photorad%nw))
    allocate(photorad%g_file(nrad_fixed, photomech%np, photorad%nw))
    photorad%nrad_file = nrad_fixed
    
    do i = 1,photomech%np
      
      if (photomech%particle_optical_type(i) == 0) then
        filename = xsroot//trim(photomech%particle_optical_prop(i))// &
                  "/mie_"//trim(photomech%particle_optical_prop(i))//".dat"
      elseif (photomech%particle_optical_type(i) == 1) then
        filename = xsroot//trim(photomech%particle_optical_prop(i))// &
                  "/frac_"//trim(photomech%particle_optical_prop(i))//".dat"
      endif
      
      if (allocated(radii)) then
        deallocate(radii, w0, qext, g)
      endif
      call read_mie_data_file(filename, photorad%nw, photorad%wavl, &
                              nrad, radii, w0, qext, g, err) 
      if (len_trim(err) /= 0) return
      if (nrad /= nrad_fixed) then
        err = "IOError: Aerosol data file "//filename// &
              "must have 50 radii bins."
        return
      endif
      
      photorad%radii_file(:,i) = radii/1.d4 ! convert from micron to cm
      photorad%w0_file(:,i,:) = w0
      photorad%qext_file(:,i,:) = qext
      photorad%g_file(:,i,:) = g
      
    enddo
    
  end subroutine
  
  subroutine get_photolysis_xs(photomech, photorad, err)
    use interp_tools, only: inter2, addpnt
    use photochem_vars, only: data_dir, xs_folder_name
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
    
    xsroot = trim(data_dir)//"/"//trim(xs_folder_name)//"/"
    
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
    use photochem_vars, only: data_dir
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
    
    rayleigh_file = trim(data_dir)//"/rayleigh/rayleigh.yaml"
    
    ! parse yaml file
    root => parse(rayleigh_file,error=error)
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
    photorad%nray = size(A)
    allocate(photorad%sigray(photorad%nray,photorad%nw))
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
  
  subroutine read_stellar_flux(star_file, nw, wavl, photon_flux, err)
    use interp_tools, only: inter2, addpnt
    use photochem_const, only: c_light, plank
    
    character(len=*), intent(in) :: star_file
    integer, intent(in) :: nw
    real(real_kind), intent(in) :: wavl(nw+1)
    real(real_kind), intent(out) :: photon_flux(nw)
    character(len=err_len), intent(out) :: err
    
    real(real_kind), allocatable :: file_wav(:), file_flux(:)
    real(real_kind) :: flux(nw)
    real(real_kind) :: dum1, dum2
    integer :: io, i, n, ierr
    real(real_kind), parameter :: rdelta = 1.d-4
    
    open(1,file=star_file,status='old',iostat=io)
    if (io /= 0) then
      err = "The input file "//star_file//' does not exist.'
      return
    endif
    
    ! count lines
    n = -1 
    read(1,*)
    do while (io == 0)
      read(1,*,iostat=io) dum1, dum2
      n = n + 1
    enddo
    
    allocate(file_wav(n+4), file_flux(n+4))
    
    ! read data
    rewind(1)
    read(1,*)
    do i = 1,n
      read(1,*,iostat=io) file_wav(i), file_flux(i)
      if (io /= 0) then
        err = "Problem reading "//star_file
        return
      endif
    enddo
    close(1)
    
    i = n
    ! interpolate 
    call addpnt(file_wav, file_flux, n+4, i, file_wav(1)*(1.d0-rdelta), 0.d0, ierr)
    call addpnt(file_wav, file_flux, n+4, i, 0.d0, 0.d0, ierr)
    call addpnt(file_wav, file_flux, n+4, i, file_wav(i)*(1.d0+rdelta), 0.d0,ierr)
    call addpnt(file_wav, file_flux, n+4, i, huge(rdelta), 0.d0,ierr)
    if (ierr /= 0) then
      err = "Problem interpolating "//trim(star_file)
      return
    endif

    call inter2(nw+1, wavl, flux, n+4, file_wav, file_flux, ierr)
    if (ierr /= 0) then
      err = "Problem interpolating "//trim(star_file)
      return
    endif
    
    ! now convert to photons/cm2/s
    do i = 1,nw
      photon_flux(i) = (1/(plank*c_light*1.d16))*flux(i)*(wavl(i+1)-wavl(i))* &
                       ((wavl(i+1)+wavl(i))/2.d0)
    enddo
    
  end subroutine
  
  
  subroutine read_atmosphere_file(atmosphere_txt, photomech, photoset, &
                                  photoinit, err)
    character(len=*), intent(in) :: atmosphere_txt
    type(PhotoMechanism), intent(in) :: photomech
    type(PhotoSettings), intent(in) :: photoset
    type(PhotoInitAtm), intent(out) :: photoinit
    character(len=err_len), intent(out) :: err
    
    character(len=10000) :: line
    character(len=:), allocatable :: message
    character(len=8) :: arr1(1000)
    character(len=24) :: arr11(1000)
    character(len=24),allocatable, dimension(:) :: labels
    integer :: ind(1)
    real(real_kind), allocatable :: temp(:)
    integer :: io, i, n, nn, iii, k, j, ii
    
    err = ''
    open(4, file=trim(atmosphere_txt),status='old',iostat=io)
    if (io /= 0) then
      err = 'Can not open file '//trim(atmosphere_txt)
      return
    endif
    read(4,'(A)') line
    
    photoinit%nzf = -1
    io = 0
    do while (io == 0)
      read(4,*,iostat=io)
      photoinit%nzf = photoinit%nzf + 1
    enddo
    
    allocate(photoinit%z_file(photoinit%nzf))
    allocate(photoinit%T_file(photoinit%nzf))
    allocate(photoinit%edd_file(photoinit%nzf))
    allocate(photoinit%usol_file(photomech%nq, photoinit%nzf))
    photoinit%z_file = 0.d0
    photoinit%T_file = 0.d0
    photoinit%edd_file = 0.d0
    photoinit%usol_file = 1.d-40
    if (photomech%there_are_particles) then
      allocate(photoinit%particle_radius_file(photomech%npq, photoinit%nzf))
    endif
    
    rewind(4)
    read(4,'(A)') line
    n = 0
    nn = 0
    do i=1,1000
      read(line,*,iostat=io) arr1(1:i)
      if (io==-1) exit
      n = n+1
    enddo
    read(4,'(A)') line
    do i=1,1000
      read(line,*,iostat=io) arr11(1:i)
      if (io==-1) exit
      nn = nn+1
    enddo
    if (n /= nn) then
      err = 'There is a missing column label in the file '//trim(atmosphere_txt)
      return
    endif
    
    allocate(labels(n))
    allocate(temp(n))
    rewind(4)
    read(4,'(A)') line
    read(line,*) (labels(i),i=1,n)
    
    ! reads in mixing ratios
    iii = 0
    do i=1,photomech%nq
      do j=1,n
        if (labels(j).eq.photomech%species_names(i)) then
          iii = iii+1
          do k = 1,photoinit%nzf
            read(4,*,iostat=io) (temp(ii),ii=1,n)
            if (io /= 0) then
              err = 'Problem reading in initial atmosphere in '//trim(atmosphere_txt)
              return
            endif
            photoinit%usol_file(i,k) = temp(j)
          enddo
          rewind(4) ! rewind!
          read(4,*) ! skip first line
          exit
        endif
      enddo
    enddo
    
    photoinit%no_water_profile = .false.
    if (photoset%fix_water_in_trop) then
      ind = findloc(labels,'H2O')
      if (ind(1) == 0) then
        photoinit%no_water_profile = .true. 
      endif
    endif
    
    if (iii.ne.photomech%nq) then
      message = 'Warning: Did not find initial data for some species in '// &
                trim(atmosphere_txt)//' . The program will assume initial mixing ratios of 1.0e-40'
      if (photoinit%no_water_profile) then
        message = message // " except H2O, which will be set to saturation in troposphere with constant "//&
                              "extrapolation above the tropopause."
      endif
      print*,message
    endif
    
    if (photomech%there_are_particles) then
      ! reads in particles radius
      iii = 1
      ! particle names
      do i = 1,photomech%npq
        ind = findloc(labels,trim(photomech%species_names(i))//"_r")
        if (ind(1) /= 0) then
          rewind(4)
          read(4,*)
          do k=1,photoinit%nzf
            read(4,*,iostat = io) (temp(ii),ii=1,n)
            if (io /= 0) then
              err = 'Problem reading in particle radius in '//trim(atmosphere_txt)
              return
            endif
            photoinit%particle_radius_file(i,k) = temp(ind(1))
          enddo
        else
          ! did not find the data
          ! will set to 0.1 micron
          photoinit%particle_radius_file(i,:) = 1.d-7
          iii = 0
        endif
      enddo
      
      if (iii == 0) then
        print*,'Warning: Did not find particle radii for some species in '//&
                trim(atmosphere_txt)//' . The program will assume 0.1 micron raddii.'
      endif
      
    endif
    
    rewind(4)
    read(4,*)
    ! reads in temperature
    ind = findloc(labels,'temp')
    if (ind(1) /= 0) then
      do k=1,photoinit%nzf
        read(4,*,iostat = io) (temp(ii),ii=1,n)
        if (io /= 0) then
          err = 'Problem reading in temperature in '//trim(atmosphere_txt)
          return
        endif
        photoinit%T_file(k) = temp(ind(1))
      enddo
    else
      err = 'temp was not found in input file '//trim(atmosphere_txt)
      return
    endif
    
    rewind(4)
    read(4,*)
    ! reads in alt
    ind = findloc(labels,'alt')
    if (ind(1) /= 0) then
      do k=1,photoinit%nzf
        read(4,*,iostat = io) (temp(ii),ii=1,n)
        if (io /= 0) then
          err = 'Problem reading in altitude in '//trim(atmosphere_txt)
          return
        endif
        photoinit%z_file(k) = temp(ind(1))*1.d5 ! conver to cm
      enddo
    else
      err = '"alt" was not found in input file '//trim(atmosphere_txt)
      return
    endif

    rewind(4)
    read(4,*)
    ! reads in eddy diffusion?
    ind = findloc(labels,'eddy')
    if (ind(1) /= 0) then
      do k=1,photoinit%nzf
        read(4,*,iostat = io) (temp(ii),ii=1,n)
        if (io /= 0) then
          err = 'Problem reading in eddy diffusion in '//trim(atmosphere_txt)
          return
        endif
        photoinit%edd_file(k) = temp(ind(1))
      enddo
    else
      err = 'eddy was not found in input file '//trim(atmosphere_txt)
      return
    endif
    
    close(4)

  end subroutine
  
end module



