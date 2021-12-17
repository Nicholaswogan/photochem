submodule (photochem_input) photochem_input_read
  implicit none
  
contains
  
  module subroutine read_all_files(mechanism_file, settings_file, flux_file, atmosphere_txt, &
                                   photodata, photovars, err)
    
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(inout) :: photovars
    character(len=err_len), intent(out) :: err
    
    err = ""
    
    ! first get SL and background species from settings
    call get_SL_and_background(settings_file, photodata, err)
    if (len(trim(err)) /= 0) return
    
    call get_photomech(mechanism_file, photodata, photovars, err)
    if (len(trim(err)) /= 0) return
    
    call get_photoset(settings_file, photodata, photovars, err)
    if (len(trim(err)) /= 0) return
    
    call get_photorad(photodata, photovars, err)
    if (len(trim(err)) /= 0) return
    
    ! stelar flux
    allocate(photovars%photon_flux(photodata%nw))
    call read_stellar_flux(flux_file, photodata%nw, photodata%wavl, photovars%photon_flux, err)
    if (len(trim(err)) /= 0) return
    
    ! initial atmosphere
    call read_atmosphere_file(atmosphere_txt, photodata, photovars, err)
    if (len(trim(err)) /= 0) return
    
  end subroutine
  
  subroutine get_photomech(infile, photodata, photovars, err) 
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(in) :: photovars
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
        call get_rxmechanism(root, infile, photodata, photovars, err)
      class default
        err = "yaml file must have dictionaries at root level"
    end select
    call root%finalize()
    deallocate(root)  
    if (len_trim(err) > 0) return
     
  end subroutine
  
  subroutine get_rxmechanism(mapping, infile, photodata, photovars, err)
    class (type_dictionary), intent(in), pointer :: mapping
    character(len=*), intent(in) :: infile
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(in) :: photovars
    character(len=err_len), intent(out) :: err
    
    class (type_dictionary), pointer :: sat_params
    class (type_list), pointer :: species, reactions, atoms
    class (type_list), pointer :: particles
    type (type_error), pointer :: io_err
    class (type_list_item), pointer :: item, next
    class (type_dictionary), pointer :: dict
    class (type_key_value_pair), pointer :: key_value_pair

    ! temporary work variables
    character(len=str_len) :: tmpchar
    character(len=str_len) :: tmp
    character(len=:), allocatable :: rxstring
    character(len=s_str_len), allocatable :: eqr(:), eqp(:)
    integer :: i, ii, j, k, kk, l, ind(1), size_eqr, size_eqp
    logical :: reverse
    ! all_species causes a small memory leak. Not sure how to free the memory properly
    type(type_list) :: all_species, all_reactions ! will include particles
    logical, allocatable :: duplicate(:)

    err = ''

    atoms => mapping%get_list('atoms',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    species => mapping%get_list('species',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    reactions => mapping%get_list('reactions',.true.,error = io_err) 
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    ! should i reverse reactions?
    photodata%reverse = mapping%get_logical('reverse-reactions',.true.,error = io_err)
    
    !!! atoms !!!
    photodata%natoms = atoms%size()
    allocate(photodata%atoms_names(photodata%natoms))
    allocate(photodata%atoms_mass(photodata%natoms))
    allocate(photodata%atoms_redox(photodata%natoms))
    
    j = 1
    item => atoms%first
    do while(associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        photodata%atoms_names(j) = element%get_string("name",error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        photodata%atoms_mass(j) = element%get_real("mass",error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        photodata%atoms_redox(j) = element%get_real("redox",error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      class default
        err = "atoms in "//trim(infile)//" must made of dictionaries"
        return 
      end select
      j = j + 1
      item => item%next
    enddo
    !!! done with atoms !!!
    
    !!! particles !!!
    ! get particles.
    particles => mapping%get_list('particles',.false.,error = io_err) 
    if (associated(particles)) then
      ! there are particles
      photodata%there_are_particles = .true.
      photodata%np = 0
      item => particles%first
      do while (associated(item))
        photodata%np = photodata%np + 1
        item => item%next
      enddo
      
      allocate(photodata%particle_names(photodata%np))
      allocate(photodata%particle_formation_method(photodata%np))
      allocate(photodata%particle_density(photodata%np))
      allocate(photodata%particle_sat_type(photodata%np))
      allocate(photodata%particle_sat_params(3,photodata%np))
      allocate(photodata%particle_gas_phase(photodata%np))
      allocate(photodata%particle_optical_prop(photodata%np))
      allocate(photodata%particle_optical_type(photodata%np))
      
      item => particles%first
      j = 1
      do while (associated(item))
        call all_species%append(item%node)
        select type (element => item%node)
        class is (type_dictionary)
          photodata%particle_names(j) = element%get_string("name",error = io_err) ! get name
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          
          tmpchar = element%get_string("formation",error = io_err) 
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          if (trim(tmpchar) == 'saturation') then
            photodata%particle_formation_method(j) = 1
          elseif (trim(tmpchar) == 'reaction') then
            photodata%particle_formation_method(j) = 2
          else
            err = "IOError: the only formation mechanism for particles is 'saturation'"
            return
          endif
          photodata%particle_density(j) = element%get_real("density",error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          photodata%particle_optical_prop(j) = element%get_string("optical-properties",error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          ! only require optical type if not "none"
          if (photodata%particle_optical_prop(j) /= 'none') then
            tmpchar = element%get_string("optical-type",error = io_err)
            if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
            if (trim(tmpchar) == "mie") then
              photodata%particle_optical_type(j) = 0
            elseif  (trim(tmpchar) == "fractal") then
              err = "IOError: 'fractal' is not an optional optical type for "// &
                    trim(photodata%particle_names(j))
              return
            else
              err = "IOError: "//trim(tmpchar)//" is not an optional optical type for "// &
                    trim(photodata%particle_names(j))
              return
            endif
          endif
  
          if (photodata%particle_formation_method(j) == 1) then
            ! there should be saturation vapor pressure information
            tmpchar = element%get_string("saturation-type",default="arrhenius",error = io_err)
            if (tmpchar == 'arrhenius') then
              photodata%particle_sat_type(j) = 1
              sat_params => element%get_dictionary('saturation-parameters',.true.,error = io_err) 
              if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
              i = 0
              key_value_pair => sat_params%first
              do while (associated(key_value_pair))
                tmpchar = trim(key_value_pair%key)
                
                if (trim(tmpchar) == "A") then
                  photodata%particle_sat_params(1,j) = sat_params%get_real(trim(tmpchar),error = io_err)
                  if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
                elseif (trim(tmpchar) == "B") then
                  photodata%particle_sat_params(2,j) = sat_params%get_real(trim(tmpchar),error = io_err)
                  if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
                elseif (trim(tmpchar) == "C") then
                  photodata%particle_sat_params(3,j) = sat_params%get_real(trim(tmpchar),error = io_err)
                  if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
                else
                  err = "Particle "//trim(photodata%particle_names(j))//" saturation parameters "//&
                        "can only be 'A', 'B', or 'C'"
                  return
                endif                
                key_value_pair => key_value_pair%next
                i = i + 1
              enddo
              if (i /= 3) then
                err = "IOError: Missing or two many saturation parameters for "//trim(photodata%particle_names(j))
                return 
              endif
            elseif (tmpchar == 'H2SO4') then
              photodata%particle_sat_type(j) = 2
              ! make a H2SO4 interpolator
              call H2SO4_interpolator(photovars, photodata%H2SO4_sat, err)
              if (len_trim(err) /= 0) return
            else
              err = "Saturation type '"//trim(tmpchar)//"' is not a valid type."
              return
            endif
            
            ! gas phase
            photodata%particle_gas_phase(j) = element%get_string("gas-phase",error = io_err) 
            if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          elseif (photodata%particle_formation_method(j) == 2) then
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
      photodata%there_are_particles = .false.
      photodata%np = 0
    endif
    
    ! for now number particle equations will be the same 
    ! as number of particles
    photodata%npq = photodata%np
    
    !!! done with particles !!!
    
    !!! species !!!
    photodata%ng = 0 ! count number of gas phase species
    item => species%first
    do while (associated(item))
      item => item%next
      photodata%ng = photodata%ng + 1
    enddo
    
    if (photodata%back_gas) then
      photodata%nll = photodata%ng - photodata%nsl - 1 ! minus 1 for background
    else
      photodata%nll = photodata%ng - photodata%nsl
    endif
    
    photodata%ng_1 = photodata%npq + 1 ! the long lived gas index
    ! photodata%nq is the last ll gas index
    
    ! now we now nq, the number of PDEs
    photodata%nq = photodata%npq + photodata%nll
    
    ! we also now nsp, the index of the backgorund gas 
    photodata%nsp = photodata%npq + photodata%ng
    
    ! species_mass, species_composition, and species_names
    ! will include the particles, thus we allocate nsp
    allocate(photodata%species_redox(photodata%nsp))
    allocate(photodata%species_mass(photodata%nsp))
    allocate(photodata%species_composition(photodata%natoms,photodata%nsp+2))
    photodata%species_composition = 0
    allocate(photodata%species_names(photodata%nsp+2))
    photodata%species_names(photodata%nsp+1) = "hv" ! always add these guys
    photodata%species_names(photodata%nsp+2) = "M"
    ! we will not include particles in thermodynamic data.
    if (photodata%reverse) then
      allocate(photodata%thermo_data(photodata%ng))
    endif
    
    ! Append the species to the end of a list
    ! which has particles in the beginning
    item => species%first
    do while (associated(item))
      call all_species%append(item%node)
      item => item%next
    enddo

    ! Loop through particles and gases
    kk = photodata%ng_1
    l = 1
    ii = 1 ! overall counter
    item => all_species%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmpchar = trim(element%get_string("name",error = io_err)) ! get name
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        
        if (ii < photodata%ng_1) then
          ! we are dealing with particles
          j = ii
        else
          ! we are dealing with gases
          ind = findloc(photodata%SL_names,tmpchar)
          if (ind(1) /= 0) then ! short lived species
            j = photodata%nq + l 
            l = l + 1
          elseif (tmpchar == photodata%back_gas_name) then ! background gas
            j = photodata%nsp
          else ! long lived species
            j = kk
            kk = kk + 1
          endif
        endif
                  
        photodata%species_names(j) = tmpchar
        dict => element%get_dictionary("composition",.true.,error = io_err)  ! get composition
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        key_value_pair => dict%first ! dont allow unspecified atoms
        do while (associated(key_value_pair))
          ind = findloc(photodata%atoms_names,trim(key_value_pair%key))
          if (ind(1) == 0) then
            err = 'IOError: The atom "'// trim(key_value_pair%key)// '" is not in the list of atoms.'
            return
          endif
          key_value_pair =>key_value_pair%next
        enddo
        
        do i=1,photodata%natoms
          photodata%species_composition(i,j) =  &
              dict%get_integer(photodata%atoms_names(i),0,error = io_err) ! no error possible.
        enddo
        photodata%species_mass(j) = sum(photodata%species_composition(:,j) * photodata%atoms_mass)
        photodata%species_redox(j) = sum(photodata%species_composition(:,j) * photodata%atoms_redox)
        
        if (photodata%reverse .and. (ii >= photodata%ng_1)) then
          call get_thermodata(element,photodata%species_names(j), infile, &
                              photodata%thermo_data(j-photodata%npq), err)
          if (len_trim(err) > 0) return
        endif
      class default
        err = "IOError: Problem with species number "//char(j)//"  in the input file"
        return
      end select
      ii = ii + 1
      item => item%next
    enddo
    
    if (l-1 /= photodata%nsl) then
      err = 'IOError: One of the short lived species is not in the file '//trim(infile)
      return
    endif
    if (photodata%back_gas) then
      ind = findloc(photodata%species_names,photodata%back_gas_name)
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
    
    if (photodata%there_are_particles) then
      ! get indexes of gas phase condensing species
      allocate(photodata%particle_gas_phase_ind(photodata%np))
      do i = 1,photodata%np
        if (photodata%particle_formation_method(i) == 1) then
          ! if a condensing molecule
          ind = findloc(photodata%species_names,photodata%particle_gas_phase(i))
          if (ind(1) /= 0) then
            photodata%particle_gas_phase_ind(i) = ind(1)
          else
            err = "IOError: particle "//trim(photodata%particle_names(i))// &
                  " can not be made from "//trim(photodata%particle_gas_phase(i))// &
                  " because "//trim(photodata%particle_gas_phase(i))//" is not a gas"// &
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

    photodata%nrF = 0 ! count forward reactions
    item => all_reactions%first
    do while (associated(item))
      item => item%next
      photodata%nrF = photodata%nrF + 1
    enddo
    
    allocate(photodata%rateparams(10,photodata%nrF))
    allocate(photodata%rxtypes(photodata%nrF))
    allocate(photodata%falloff_type(photodata%nrF))
    allocate(photodata%num_efficient(photodata%nrF))
    photodata%falloff_type = - huge(1)
    ! determine which reactions to reverse. 
    ! Determine maximum number of reactants, and productants
    photodata%max_num_reactants = 1
    photodata%max_num_products = 1
    photodata%nrR = 0
    photodata%kj = 0
    j = 1
    item => all_reactions%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmp = trim(element%get_string("equation",error = io_err))
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        call get_rateparams(element, infile, photodata%rxtypes(j), &
                            photodata%falloff_type(j), photodata%rateparams(:,j), err)
        if (len_trim(err) > 0) return        
        
        call parse_reaction(tmp, reverse, eqr, eqp, err)
        if (len_trim(err) > 0) return
        
        call compare_rxtype_string(tmp, eqr, eqp, reverse, photodata%rxtypes(j),err)
        if (len_trim(err) > 0) return
        
        size_eqp = size(eqp)
        size_eqr = size(eqr)
        if ((photodata%rxtypes(j) == 2) .or. (photodata%rxtypes(j) == 3)) then ! if threebody or falloff
          size_eqr = size_eqr - 1 ! remove the M
          size_eqp = size_eqp - 1
        endif
        
        if (reverse) then
          if (.not.photodata%reverse) then
            err = 'IOError: reaction file '//trim(infile)//' contains reverse reaction '//tmp// &
                  ', which is incompatible with "reverse-reactions: false"'
            return
          endif
          photodata%nrR = photodata%nrR + 1
          if (size_eqr > photodata%max_num_products) photodata%max_num_products = size_eqr
          if (size_eqp > photodata%max_num_reactants) photodata%max_num_reactants = size_eqp
        endif
        if (size_eqr > photodata%max_num_reactants) photodata%max_num_reactants = size_eqr
        if (size_eqp > photodata%max_num_products) photodata%max_num_products = size_eqp
        
        call count_efficiencies(element, photodata%num_efficient(j))
        
        if (photodata%rxtypes(j) == 0) then ! if photolysis reaction
          photodata%kj = photodata%kj + 1
        endif
      class default
        err = "IOError: Problem with reaction number "//char(j)//" in the input file."
        return
      end select
      item => item%next
      j = j+1
    enddo
    photodata%nrT = photodata%nrR + photodata%nrF
    
    ! allocate stuff and loop through reactions again
    allocate(duplicate(photodata%nrT))
    allocate(photodata%nreactants(photodata%nrT))
    allocate(photodata%nproducts(photodata%nrT))
    allocate(photodata%reactants_sp_inds(photodata%max_num_reactants,photodata%nrT))
    allocate(photodata%products_sp_inds(photodata%max_num_products,photodata%nrT))
    photodata%reactants_sp_inds = huge(1)
    photodata%products_sp_inds = huge(1)
    if (photodata%reverse) then
      allocate(photodata%reverse_info(photodata%nrT))
      photodata%reverse_info = 0 ! initialize
    endif
    allocate(photodata%reactants_names(photodata%max_num_reactants,photodata%nrF))
    allocate(photodata%products_names(photodata%max_num_products,photodata%nrF))
    ! efficiency stuff
    allocate(photodata%efficiencies(maxval(photodata%num_efficient),photodata%nrF))
    allocate(photodata%eff_sp_inds(maxval(photodata%num_efficient),photodata%nrF))
    allocate(photodata%def_eff(photodata%nrF))
    photodata%efficiencies = -huge(1.d0) ! so everything blows up if we make a mistake
    photodata%eff_sp_inds = -huge(0)
    photodata%def_eff = 1.d0 ! default is 1
    
    j = 1
    k = 1
    item => all_reactions%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmp = trim(element%get_string("equation",error = io_err))
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        call get_efficient(element, j, infile, photodata, err)
        if (len_trim(err)>0) return
        
        call get_reaction_chars(tmp, photodata%max_num_reactants, photodata%max_num_products, &
                            photodata%nreactants(j), photodata%nproducts(j), &
                            photodata%reactants_names(:,j), photodata%products_names(:,j), reverse, err)
        if (len_trim(err)>0) return
        
        call get_reaction_sp_nums(tmp, photodata%max_num_reactants, photodata%max_num_products, &
                                 photodata%reactants_names(:,j), photodata%products_names(:,j), &
                                 photodata%species_names, photodata%species_composition, &
                                 photodata%natoms, photodata%nsp, &
                                 photodata%reactants_sp_inds(:,j), photodata%products_sp_inds(:,j), err)
        if (len_trim(err)>0) return
        
        ! check if duplicate
        duplicate(j) = element%get_logical("duplicate",default=.false.,error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        if (reverse) then
          ! reaction has a reverse
          i = photodata%nrF + k
          photodata%reverse_info(j) = i ! the reaction number of reversed reaction
          photodata%reverse_info(i) = j ! the reaction number of the forward reaction
          photodata%nreactants(i) = photodata%nproducts(j)
          photodata%nproducts(i) = photodata%nreactants(j)
          
          duplicate(i) = duplicate(j)

          kk = photodata%nreactants(i)
          l = photodata%nreactants(j)
          
          photodata%reactants_sp_inds(1:kk,i) = photodata%products_sp_inds(1:kk,j)
          photodata%products_sp_inds(1:l,i) = photodata%reactants_sp_inds(1:l,j)
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
    allocate(photodata%nump(photodata%nsp))
    allocate(photodata%numl(photodata%nsp))
    photodata%numl = 0
    photodata%nump = 0
    do j = 1,photodata%nrT
      k = photodata%nreactants(j)
      do i = 1,k
        kk = photodata%reactants_sp_inds(i,j)
        if ((kk /= photodata%nsp+1) .and. (kk /= photodata%nsp+2)) then
          photodata%numl(kk) = photodata%numl(kk) + 1
        endif
      enddo
      k = photodata%nproducts(j)
      do i = 1,k
        kk = photodata%products_sp_inds(i,j)
        if ((kk /= photodata%nsp+1) .and. (kk /= photodata%nsp+2)) then
          photodata%nump(kk) = photodata%nump(kk) + 1
        endif
      enddo
    enddo
    allocate(photodata%iprod(maxval(photodata%nump),photodata%nsp))
    allocate(photodata%iloss(maxval(photodata%numl),photodata%nsp))
    photodata%iprod = 0
    photodata%iloss = 0
    photodata%numl = 0
    photodata%nump = 0
    ! loop again and get iprod and iloss
    do j = 1,photodata%nrT
      k = photodata%nreactants(j)
      do i = 1,k
        kk = photodata%reactants_sp_inds(i,j)
        if ((kk /= photodata%nsp+1) .and. (kk /= photodata%nsp+2)) then
          photodata%numl(kk) = photodata%numl(kk) + 1
          l = photodata%numl(kk)
          photodata%iloss(l,kk) = j
        endif
      enddo
      k = photodata%nproducts(j)
      do i = 1,k
        kk = photodata%products_sp_inds(i,j)
        if ((kk /= photodata%nsp+1) .and. (kk /= photodata%nsp+2)) then
          photodata%nump(kk) = photodata%nump(kk) + 1
          l = photodata%nump(kk)
          photodata%iprod(l,kk) = j
        endif
      enddo
    enddo
    
    ! photolysis
    allocate(photodata%photonums(photodata%kj))
    j = 1
    do i = 1, photodata%nrF
      if (photodata%rxtypes(i) == 0) then
        photodata%photonums(j) = i
        j = j + 1
      endif
    enddo
    
    ! save reaction names
    allocate(photodata%reaction_equations(photodata%nrT))
    do i = 1,photodata%nrT
      call reaction_string(photodata,i,rxstring)
      photodata%reaction_equations(i) = rxstring
    enddo
    
    call check_for_duplicates(photodata, duplicate, err)
    if (len(trim(err)) > 0) return
    !!! end reactions !!!
    
    !!! henrys law !!!
    call get_henry(photodata, photovars, err)
    if (len(trim(err)) > 0) return
    !!! end henrys law !!!
    
  end subroutine
  
  subroutine H2SO4_interpolator(photovars, s2, err)
    use linear_interpolation_module, only: linear_interp_2d
    type(PhotochemVars), intent(in) :: photovars
    type(linear_interp_2d), intent(out) :: s2
    character(len=err_len), intent(out) :: err
    
    real(real_kind), allocatable :: H2O(:)
    real(real_kind), allocatable :: Temp(:)
    real(real_kind), allocatable :: H2SO4(:,:)
    integer :: io
    integer :: nT, nH2O
    character(len=:), allocatable :: filename
    
    err = ""
    
    filename = trim(photovars%data_dir)//"/misc/H2SO4.dat"
    open(unit=1,file=filename, status='old',iostat=io,form='unformatted')
    if (io /= 0) then
      err = "Could not open "//trim(filename)
      return
    endif
    read(1) nT
    read(1) nH2O
    allocate(Temp(nT))
    allocate(H2O(nH2O))
    allocate(H2SO4(nT,nH2O))
    read(1) Temp
    read(1) H2O
    read(1) H2SO4
    close(1)
    
    call s2%initialize(Temp, H2O, H2SO4, io)
    if (io /= 0) then
      err = "Failed to initialize H2SO4 interpolator."
      return
    endif
    
  end subroutine
  
  subroutine get_henry_parse(root, photodata, photovars, henry_names, henry_data, err)
    class (type_list), intent(in) :: root
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(in) :: photovars
    character(len=s_str_len), allocatable, intent(out) :: henry_names(:)
    real(real_kind), allocatable, intent(out) :: henry_data(:,:)
    character(len=*), intent(out) :: err
  
    type (type_error), pointer :: io_err
    class (type_list_item), pointer :: item
    integer :: j  
    
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
  end subroutine
  
  subroutine get_henry(photodata, photovars, err)
    use yaml, only : parse, error_length
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(in) :: photovars
    character(len=*), intent(out) :: err
    
    character(error_length) :: error
    class (type_node), pointer :: root
    integer :: j, ind(1), i
    
    character(len=s_str_len), allocatable :: henry_names(:)
    real(real_kind), allocatable :: henry_data(:,:)
    
    err = ''

    ! parse yaml file
    root => parse(trim(photovars%data_dir)//"/henry/henry.yaml",error=error)
    if (error /= '') then
      err = trim(error)
      return
    end if
    select type (root)
    class is (type_list)
      call get_henry_parse(root, photodata, photovars, henry_names, henry_data, err)
    class default
      err = "yaml file must have dictionaries at root level"
    end select
    call root%finalize()
    deallocate(root)
    if (len_trim(err) /= 0) return
    
    allocate(photodata%henry_data(2,photodata%nsp))
    photodata%henry_data = 0.d0
    do j = 1,size(henry_names)
      ind = findloc(photodata%species_names,henry_names(j))
      if (ind(1) /= 0) then
        i = ind(1)
        photodata%henry_data(:,i) = henry_data(:,j)
      endif
    enddo
    ! set particle solubility to super high number
    photodata%henry_data(1,1:photodata%npq) = 7.d11
    
  end subroutine
  
  subroutine get_efficient(reaction, rxn, infile, photodata, err)
    class(type_dictionary), intent(in) :: reaction
    integer, intent(in) :: rxn
    character(len=*), intent(in) :: infile
    type(PhotochemData), intent(inout) :: photodata
    character(len=*), intent(out) :: err
    
    class(type_dictionary), pointer :: tmpdict
    class (type_key_value_pair), pointer :: key_value_pair
    type (type_error), pointer :: io_err
    character(len=s_str_len) :: rxn_str
    integer :: ind(1), j
    
    tmpdict => reaction%get_dictionary("efficiencies",.false.,error = io_err)
    
    if (associated(tmpdict)) then
      j = 1
      key_value_pair => tmpdict%first
      do while (associated(key_value_pair))
        
        ind = findloc(photodata%species_names,trim(key_value_pair%key))
        if (ind(1) == 0) then
          write(rxn_str,*) rxn
          err = 'IOError: Reaction number '//trim(adjustl(rxn_str))//' has efficiencies for species that are'// &
          ' not in the list of species'
        endif
        photodata%eff_sp_inds(j,rxn) = ind(1)
        photodata%efficiencies(j,rxn) = tmpdict%get_real(trim(key_value_pair%key),error = io_err)
        if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif

        key_value_pair => key_value_pair%next
        j = j + 1
      enddo
    endif
    
    photodata%def_eff(rxn) = reaction%get_real("default-efficiency",1.d0,error = io_err)

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
  
  
  subroutine check_for_duplicates(photodata, duplicate, err)
    use sorting, only: sort
    type(PhotochemData), intent(in) :: photodata
    logical, intent(in) :: duplicate(:)
    character(len=err_len), intent(out) :: err
    character(len=:), allocatable :: rxstring
    
    integer i, ii, j, jj, nr, np, rxt, rxt_ii
    logical l, m
    
    integer, allocatable :: tmp_arr1(:), tmp_arr2(:)
    
    err = ''
    
    if (photodata%max_num_reactants > photodata%max_num_products) then
      allocate(tmp_arr1(photodata%max_num_reactants))
      allocate(tmp_arr2(photodata%max_num_reactants))
    else
      allocate(tmp_arr1(photodata%max_num_products))
      allocate(tmp_arr2(photodata%max_num_products))
    endif
    
    do i = 1,photodata%nrT-1
      do ii = i+1,photodata%nrT
        
        ! if not designated as duplicates then check if they are
        ! duplicates
        if (.not.(duplicate(i) .and. duplicate(ii))) then
        
        ! check the same num of reactants and products
        nr = photodata%nreactants(i)
        np = photodata%nproducts(i)
        if (i > photodata%nrF) then
          j = photodata%reverse_info(i)
        else
          j = i
        endif
        if (ii > photodata%nrF) then
          jj = photodata%reverse_info(ii)
        else
          jj = ii
        endif
        rxt = photodata%rxtypes(j)
        rxt_ii = photodata%rxtypes(jj)
        if (nr == photodata%nreactants(ii) .and. np == photodata%nproducts(ii) &
            .and. rxt == rxt_ii) then
          tmp_arr1(1:nr) = photodata%reactants_sp_inds(1:nr,ii)
          tmp_arr2(1:nr) = photodata%reactants_sp_inds(1:nr,i)
          
          call sort(tmp_arr1(1:nr))
          call sort(tmp_arr2(1:nr))
          
          m = all(tmp_arr1(1:nr) == tmp_arr2(1:nr))
          if (m) then
          
            tmp_arr1(1:np) = photodata%products_sp_inds(1:np,ii)
            tmp_arr2(1:np) = photodata%products_sp_inds(1:np,i)
            
            call sort(tmp_arr1(1:np))
            call sort(tmp_arr2(1:np))
          
            l = all(tmp_arr1(1:np) == tmp_arr2(1:np))
          
            if (l) then
              err = "IOError: This reaction is a duplicate: "
              call reaction_string(photodata,i,rxstring)
              err(len_trim(err)+2:) = rxstring
              return
            endif
          endif
        endif
        
        endif
      enddo
    enddo
    
  end subroutine
  
  subroutine reaction_string(photodata,rxn,rxstring)
    type(PhotochemData), intent(in) :: photodata
    integer, intent(in) :: rxn
    character(len=:), allocatable, intent(out) :: rxstring
    integer j, k, i
    rxstring = ''
    if (rxn > photodata%nrF) then
      i = photodata%reverse_info(rxn)
    else
      i = rxn
    endif
    do j = 1,photodata%nreactants(rxn)-1
      k = photodata%reactants_sp_inds(j,rxn)
      rxstring = rxstring //(trim(photodata%species_names(k))//' + ')
    enddo
    
    k = photodata%reactants_sp_inds(photodata%nreactants(rxn),rxn)
    rxstring = rxstring // trim(photodata%species_names(k))//' => '
    
    if ((photodata%rxtypes(i) == 2) .or.(photodata%rxtypes(i) == 3)) then
      rxstring = rxstring(1:len(rxstring)-4) //(' + M'//' => ')
    endif
    
    do j = 1,photodata%nproducts(rxn)-1
      k = photodata%products_sp_inds(j,rxn)
      rxstring = rxstring // trim(photodata%species_names(k))//' + '
    enddo
    k = photodata%products_sp_inds(photodata%nproducts(rxn),rxn)
    rxstring = rxstring // trim(photodata%species_names(k))
    
    if ((photodata%rxtypes(i) == 2) .or.(photodata%rxtypes(i) == 3)) then
      rxstring = rxstring //' + M'
    endif
  end subroutine
  
  subroutine compare_rxtype_string(tmp, eqr, eqp, reverse, rxtype_int, err)
    character(len=*), intent(in) :: tmp
    character(len=s_str_len), allocatable, intent(in) :: eqr(:), eqp(:)
    logical, intent(in) :: reverse
    integer, intent(in) :: rxtype_int
    character(len=err_len), intent(out) :: err
    character(len=s_str_len) :: rxtype
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
                            thermo, err)
    use photochem_types, only: ThermodynamicData
    class(type_dictionary), intent(in) :: molecule
    character(len=*), intent(in) :: molecule_name
    character(len=*), intent(in) :: infile
    
    type(ThermodynamicData), intent(out) :: thermo
    character(len=err_len), intent(out) :: err
    
    type (type_error), pointer :: io_err
    class(type_dictionary), pointer :: tmpdict
    class(type_list), pointer :: tmplist
    class(type_list_item), pointer :: item, item1
    character(len=:), allocatable :: model
    logical :: success
    
    integer :: j, k
    
    err = ''
    
    tmpdict => molecule%get_dictionary("thermo",.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    ! check thermodynamic model
    model = tmpdict%get_string("model",error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (model == "Shomate") then
      thermo%dtype = 1
    elseif (model == "NASA9") then
      thermo%dtype = 2
    else
      err = "IOError: Thermodynamic data must be in Shomate format for "//trim(molecule_name)
      return
    endif
    
    ! get temperature ranges
    tmplist =>tmpdict%get_list("temperature-ranges",.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    thermo%ntemps = tmplist%size() - 1
    if (thermo%ntemps /= 1 .and. thermo%ntemps /= 2) then
      err = "IOError: Problem reading thermodynamic data for  "//trim(molecule_name)
      return
    endif
    allocate(thermo%temps(thermo%ntemps + 1))
    
    j = 1
    item => tmplist%first
    do while (associated(item))
      select type (listitem => item%node)
      class is (type_scalar)
        thermo%temps(j) = listitem%to_real(-1.d0,success)
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
    
    ! get data
    tmplist => tmpdict%get_list("data",.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    if (tmplist%size() /= thermo%ntemps) then
      err = "IOError: Problem reading thermodynamic data for "//trim(molecule_name)
      return
    endif
    
    if (thermo%dtype == 1) then
      ! Shomate
      allocate(thermo%data(7,thermo%ntemps))
    elseif (thermo%dtype == 2) then
      ! NASA9
      allocate(thermo%data(9,thermo%ntemps))
    endif
    
    k = 1
    item => tmplist%first
    do while (associated(item))
      select type (listitem => item%node)
      class is (type_list)
        
        if (listitem%size() /= size(thermo%data,1)) then
          err = "IOError: Too much or too little thermodynamic data for "//trim(molecule_name)
          return
        endif
        
        j = 1
        item1 => listitem%first
        do while (associated(item1)) 
          select type (listitem1 => item1%node)
          class is (type_scalar)

            thermo%data(j, k) = listitem1%to_real(-1.d0,success)
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
      class default
        err = "IOError: Problem reading thermodynamic data for "//trim(molecule_name)
        return
      end select
      item => item%next
      k = k + 1
    enddo          
                            
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
    logical :: use_jpl
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
        
      use_jpl = reaction%get_logical('JPL',default=.false.,error = io_err)
      if (use_jpl) then
        falloff_type = 3
      else
        falloff_type = 0
      endif

      nullify(tmpdict)
      tmpdict => reaction%get_dictionary('Troe',.false.,error = io_err)
      if (associated(tmpdict)) then
        if (use_jpl) then
          err = "Both 'Troe' and 'JPL' falloff types are specified for reaction "// &
                trim(reaction%get_string("equation",error = io_err))//". Only one is allowed"
          return
        endif
        
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
    character(len=s_str_len), allocatable, intent(out) :: eqr(:), eqp(:)
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
    
    character(len=s_str_len), allocatable :: eqr(:), eqp(:)
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
  
  subroutine SL_and_background(root, infile, photodata, err)
    class (type_dictionary), intent(in) :: root
    character(len=*), intent(in) :: infile
    type(PhotochemData), intent(inout) :: photodata
    character(len=err_len), intent(out) :: err
  
    class (type_dictionary), pointer :: tmp1
    type (type_error), pointer :: io_err
    class (type_list), pointer :: bcs
    class (type_list_item), pointer :: item
    character(len=str_len) :: spec_type
    
    ! get background species
    tmp1 => root%get_dictionary('planet',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    photodata%back_gas = tmp1%get_logical('use-background-gas',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    if (photodata%back_gas) then
      photodata%back_gas_name = tmp1%get_string('background-gas',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    else
      photodata%back_gas_ind = -1
      err = "Currently, the model requires there to be a background gas."// &
            " You must set 'use-background-gas: true'"
      return
    endif
    
    bcs => root%get_list('boundary-conditions',.true.,error = io_err)
    
    photodata%nsl = 0
    item => bcs%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        spec_type = element%get_string('type','long lived',error = io_err)
        if (spec_type == 'short lived') then
          photodata%nsl = photodata%nsl + 1
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
    allocate(photodata%SL_names(photodata%nsl))
    
    photodata%nsl = 0
    item => bcs%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        spec_type = element%get_string('type','long lived',error = io_err)
        if (spec_type == 'short lived') then
          photodata%nsl = photodata%nsl + 1
          photodata%SL_names(photodata%nsl) = element%get_string('name',error = io_err)
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
    
  end subroutine
  
  subroutine get_SL_and_background(infile, photodata, err)
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    type(PhotochemData), intent(inout) :: photodata
    character(len=err_len), intent(out) :: err
  
    character(error_length) :: error
    class (type_node), pointer :: root
    integer :: i, j
    err = ''
    
    root => parse(infile,error=error)
    if (error /= '') then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        call SL_and_background(root, infile, photodata, err)
      class default
        err = trim(infile)//" file must have dictionaries at root level"
    end select
    call root%finalize()
    deallocate(root)
    if (len_trim(err) /= 0) return
    
    ! check for duplicates
    do i = 1,photodata%nsl-1
      do j = i+1,photodata%nsl
        if (photodata%SL_names(i) == photodata%SL_names(j)) then
          err = 'IOError: Short lived species '//trim(photodata%SL_names(i))// &
                ' is a duplicate.'
          return
        endif
      enddo
    enddo
  
  end subroutine
  
  subroutine get_photoset(infile, photodata, photovars, err)
    use yaml, only : parse, error_length
    character(len=*), intent(in) :: infile
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(inout) :: photovars
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
        call unpack_settings(root, infile, photodata, photovars, err)
      class default
        err = trim(infile)//" file must have dictionaries at root level"
    end select
    call root%finalize()
    deallocate(root)
    if (len_trim(err) /= 0) return

  end subroutine
  
  subroutine unpack_settings(mapping, infile, photodata, photovars, err)
    class (type_dictionary), intent(in), pointer :: mapping
    character(len=*), intent(in) :: infile
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(inout) :: photovars
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
    photodata%regular_grid = tmp1%get_logical('regular-grid',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    if (photodata%regular_grid) then
      photodata%lower_wavelength = tmp1%get_real('lower-wavelength',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      photodata%upper_wavelength = tmp1%get_real('upper-wavelength',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      photodata%nw = tmp1%get_integer('number-of-bins',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      if (photodata%nw < 1) then 
        err = 'Number of photolysis bins must be >= 1 in '//trim(infile)
        return
      endif
      if (photodata%lower_wavelength > photodata%upper_wavelength) then
        err = 'lower-wavelength must be smaller than upper-wavelength in '//trim(infile)
        return
      endif
    else
      photodata%grid_file = tmp1%get_string('input-file',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    endif
    ! scale factor for photon flux. Its optional
    photovars%photon_scale_factor = tmp1%get_real('photon-scale-factor', 1.d0,error = io_err)
    
    ! atmosphere grid
    tmp1 => mapping%get_dictionary('atmosphere-grid',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photovars%bottom_atmos = tmp1%get_real('bottom',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photovars%top_atmos = tmp1%get_real('top',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photovars%nz = tmp1%get_integer('number-of-layers',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    ! Planet
    tmp1 => mapping%get_dictionary('planet',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    photovars%surface_pressure = tmp1%get_real('surface-pressure',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photovars%surface_pressure <= 0.d0) then
      err = 'IOError: Planet surface pressure must be greater than zero.'
      return
    endif
    photodata%planet_mass = tmp1%get_real('planet-mass',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photodata%planet_mass < 0.d0) then
      err = 'IOError: Planet mass must be greater than zero.'
      return
    endif
    photodata%planet_radius = tmp1%get_real('planet-radius',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photodata%planet_radius < 0.d0) then
      err = 'IOError: Planet radius must be greater than zero.'
      return
    endif
    photovars%surface_albedo = tmp1%get_real('surface-albedo',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photovars%surface_albedo < 0.d0) then
      err = 'IOError: Surface albedo must be greater than zero.'
      return
    endif
    photovars%diurnal_fac = tmp1%get_real('diurnal-averaging-factor',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photovars%diurnal_fac < 0.d0 .or. photovars%diurnal_fac > 1.d0) then
      err = 'IOError: diurnal-averaging-factor must be between 0 and 1.'
      return
    endif
    
    photovars%solar_zenith = tmp1%get_real('solar-zenith-angle',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    if (photovars%solar_zenith < 0.d0 .or. photovars%solar_zenith > 90.d0) then
      err = 'IOError: solar zenith must be between 0 and 90.'
      return
    endif
    
    ! H2 escape
    photodata%diff_H_escape = tmp1%get_logical('diff-lim-hydrogen-escape',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    ind = findloc(photodata%species_names,'H2')
    photodata%LH2 = ind(1)
    if (ind(1) == 0 .and. photodata%diff_H_escape) then
      err = 'IOError: H2 must be a species if diff-lim-hydrogen-escape = True.'
      return
    endif
    ind = findloc(photodata%species_names,'H')
    photodata%LH = ind(1)
    if (ind(1) == 0 .and. photodata%diff_H_escape) then
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
    photodata%fix_water_in_trop = tmp2%get_logical('fix-water-in-troposphere',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photodata%water_cond = tmp2%get_logical('water-condensation',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    photodata%gas_rainout = tmp2%get_logical('gas-rainout',error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    ind = findloc(photodata%species_names,'H2O')
    photodata%LH2O = ind(1)
    if (ind(1) == 0 .and. photodata%fix_water_in_trop) then
      err = 'IOError: H2O must be a species if water-saturated-troposhere = True.'
      return
    elseif (ind(1) == 0 .and. photodata%water_cond) then
      err = 'IOError: H2O must be a species if water-condensation = True.'
      return
    elseif (ind(1) == 0 .and. photodata%gas_rainout) then
      err = 'IOError: H2O must be a species if gas-rainout = True.'
      return
    elseif (ind(1) == 0 .and. photodata%there_are_particles) then
      if (any(photodata%particle_sat_type == 2)) then
        err = 'IOError: H2O must be a species if H2SO4 condensation is on.'
        return
      endif
    endif
    
    if (photodata%fix_water_in_trop) then  
      
      temp_char = tmp2%get_string('relative-humidity',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
      read(temp_char,*,iostat = io) photovars%relative_humidity
      
      if (io /= 0) then
        ! it isn't a float
        if (trim(temp_char) == "manabe") then
          photovars%use_manabe = .true.
        else
          err = '"relative-humidity" can only be a number between 0 and 1, or "manabe". See '//trim(infile)
          return 
        endif
      else
        photovars%use_manabe = .false.
      endif
      
    endif
    
    if (photodata%gas_rainout) then
      photovars%rainfall_rate = tmp2%get_real('rainfall-rate',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    endif
    
    if (photodata%fix_water_in_trop .or. photodata%gas_rainout) then
      ! we need a tropopause altitude
      photovars%trop_alt = tmp2%get_real('tropopause-altitude',error = io_err)
      if (associated(io_err)) then
        err = "tropopause-altitude must be specified if fix-water-in-troposphere = true, or gas-rainout = true"
        return
      endif
      if ((photovars%trop_alt < photovars%bottom_atmos) .or. &
          (photovars%trop_alt > photovars%top_atmos)) then
          err = 'IOError: tropopause-altitude must be between the top and bottom of the atmosphere'
          return
      endif
      
    endif
    
    if (photodata%water_cond) then        
      
      tmp3 => tmp2%get_dictionary('condensation-rate',.true.,error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      photovars%H2O_condensation_rate(1) = tmp3%get_real('A',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      photovars%H2O_condensation_rate(2) = tmp3%get_real('rhc',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      photovars%H2O_condensation_rate(3) = tmp3%get_real('rh0',error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      if (photovars%H2O_condensation_rate(3) <= photovars%H2O_condensation_rate(2)) then
        err = 'IOError: Rate constant "rh0" for H2O condensation must be > "rhc". See '//trim(infile)
        return
      endif
    endif
    
    ! particle parameters
    if (photodata%there_are_particles) then
      particles => mapping%get_list('particles',.true.,error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
      allocate(particle_checklist(photodata%np))
      allocate(photovars%condensation_rate(3,photodata%np))
      particle_checklist = .false.
      
      item => particles%first
      do while (associated(item))
        select type (element => item%node)
        class is (type_dictionary)
          temp_char = element%get_string('name',error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          ind = findloc(photodata%particle_names,trim(temp_char))
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
          
          photovars%condensation_rate(1,ind(1)) = tmp1%get_real('A',error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          photovars%condensation_rate(2,ind(1)) = tmp1%get_real('rhc',error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          photovars%condensation_rate(3,ind(1)) = tmp1%get_real('rh0',error = io_err)
          if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          if (photovars%condensation_rate(3,ind(1)) <= photovars%condensation_rate(2,ind(1))) then
            err = 'IOError: Rate constant "rh0" for '//trim(temp_char) &
                  //' condensation must be > "rhc". See '//trim(infile)
            return
          endif
          
        class default
          err = "IOError: Particle settings must be a list of dictionaries."
          return
        end select
        item => item%next
      end do
      
      do i = 1,photodata%np
        if (photodata%particle_formation_method(i) == 1 .and. .not. particle_checklist(i)) then
          err = 'IOError: Particle '//trim(photodata%particle_names(i))// &
                ' does not have any condensation rate data in the file '//trim(infile)
          return
        endif
      enddo
        
    endif
      
    ! boundary conditions and species types
    bcs => mapping%get_list('boundary-conditions',.true.,error = io_err)
    if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    if (photodata%back_gas) then
      ind = findloc(photodata%species_names,trim(photodata%back_gas_name))
      photodata%back_gas_ind = ind(1)
      photodata%back_gas_mu = photodata%species_mass(ind(1))
    endif
    
    ! allocate boundary conditions
    allocate(photovars%lowerboundcond(photodata%nq))
    allocate(photovars%lower_vdep(photodata%nq))
    allocate(photovars%lower_flux(photodata%nq))
    allocate(photovars%lower_dist_height(photodata%nq))
    allocate(photovars%lower_fix_mr(photodata%nq))
    allocate(photovars%upperboundcond(photodata%nq))
    allocate(photovars%upper_veff(photodata%nq))
    allocate(photovars%upper_flux(photodata%nq))
    ! default boundary conditions
    photovars%lowerboundcond = default_lowerboundcond
    photovars%lower_vdep = 0.d0
    photovars%upperboundcond = 0
    photovars%upper_veff = 0.d0
      
    ! determine number of short lived species. 
    allocate(dups(photodata%nsp))
    j = 1
    item => bcs%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        if (j > photodata%nsp) then
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
        if (photodata%fix_water_in_trop .and. dups(j) == "H2O") then
          err = "IOError: H2O can not have a specified boundary condition"// &
                " if water-saturated-troposphere = true in the settings file."
          return
        endif
        ! check if in rxmech
        ind = findloc(photodata%species_names,dups(j))
        if (ind(1) == 0) then
          err = "IOError: Species "//trim(dups(j))// &
          ' in settings file is not in the reaction mechanism file.'
          return 
        endif
        if ((ind(1) == photodata%back_gas_ind) .and. (photodata%back_gas)) then ! can't be background gas
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
                                 photovars%lowerboundcond(i), photovars%lower_vdep(i), &
                                 photovars%lower_flux(i), photovars%lower_dist_height(i), &
                                 photovars%lower_fix_mr(i), &
                                 photovars%upperboundcond(i), photovars%upper_veff(i), &
                                 photovars%upper_flux(i), err)
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
    if (photodata%diff_H_escape) then
      if (photodata%back_gas_name /= "H2") then
        if (photovars%upperboundcond(photodata%LH2) /= 0) then
          err = "IOError: H2 must have a have a effusion velocity upper boundary"// &
                " if diff-lim-hydrogen-escape = True"
          return
        endif
      endif
      if (photovars%upperboundcond(photodata%LH) /= 0) then
        err = "IOError: H must have a have a effusion velocity upper boundary"// &
              " if diff-lim-hydrogen-escape = True"
        return
      endif
    endif

    ! check for SL nonlinearities
    call check_sl(photodata, err)
    if (len_trim(err) /= 0) return

    
  end subroutine
  
  subroutine check_sl(photodata, err)
    type(PhotochemData), intent(in) :: photodata
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, l, k, kk, m, mm, n, nn, ind(1), counter
    character(len=:), allocatable :: reaction
    err = ''
    
    do i = 1, photodata%nsl
      j = photodata%nq + i
      ! can not be an efficiency.
      do k = 1,photodata%nrF
        if ((photodata%rxtypes(k) == 2) .or. (photodata%rxtypes(k) == 3)) then ! if three body or falloff
          ind = findloc(photodata%eff_sp_inds(:,k),j)
          if (ind(1) /= 0) then
            call reaction_string(photodata,k,reaction)
            err = 'IOError: Reaction "'//reaction//'" has short-lived species collision efficiencies.' // &
            ' This is not allowed. Either remove the efficiencies, or change the species to long lived.'
            return
          endif
        endif
      enddo
      
      l = photodata%nump(j)
      do k = 1,l
        kk = photodata%iprod(k,j)
        call reaction_string(photodata,kk,reaction)
        m = photodata%nreactants(kk)
        do mm = 1, m
          ! are SL species produced by other SL species?
          do n = photodata%nq+1,photodata%nq + photodata%nsl
            if (n == photodata%reactants_sp_inds(mm,kk)) then
              err = 'IOError: Reaction "'//reaction//'" has short-lived species as reactants'// &
              ' and products. This is not allowed. Change one or both of the species to long-lived.'
              return
            endif
          enddo
        enddo
      enddo

      l = photodata%numl(j)
      do k = 1,l
        kk = photodata%iloss(k,j)
        call reaction_string(photodata,kk,reaction)
        m = photodata%nreactants(kk)
        
        counter = 0
        do mm = 1, m
          n = photodata%reactants_sp_inds(mm,kk)
          do nn = photodata%nq+1,photodata%nq + photodata%nsl
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
            elseif (photodata%species_names(n) == 'hv') then
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
    if (bctype == "vdep") then
      lowercond = 0
      Lvdep = tmpdict%get_real("vdep",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      Lflux = 0.d0
      LdistH = 0.d0
      Lmr = 0.d0 
    elseif (bctype == "mix") then
      lowercond = 1
      Lvdep = 0.d0
      Lflux = 0.d0
      LdistH = 0.d0
      Lmr = tmpdict%get_real("mix",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    elseif (bctype == "flux") then
      lowercond = 2
      Lvdep = 0.d0
      Lflux = tmpdict%get_real("flux",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      LdistH = 0.d0
      Lmr = 0.d0 
    elseif (bctype == "vdep + dist flux") then
      lowercond = 3
      Lvdep = tmpdict%get_real("vdep",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      Lflux = tmpdict%get_real("flux",error = io_err)
      if (associated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      LdistH = tmpdict%get_real("height",error = io_err)
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
    if (bctype == "veff") then
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
    
  subroutine get_photorad(photodata, photovars, err)
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(in) :: photovars
    character(len=err_len), intent(out) :: err
    
    integer :: i
    err = ''
    
    ! compute wavelength grid
    if (photodata%regular_grid) then
      allocate(photodata%wavl(photodata%nw+1))
      photodata%wavl(1) = photodata%lower_wavelength
      do i = 2,photodata%nw+1
        photodata%wavl(i) = photodata%wavl(i-1) + &
                           (photodata%upper_wavelength - photodata%lower_wavelength)/photodata%nw
      enddo
    else
      ! read file
      err = 'Still need to add support for reading in wavelength grid from file'
      return
    endif
    
    ! get rayleigh
    call get_rayleigh(photodata, photovars, err)
    if (len_trim(err) /= 0) return
    
    ! get photolysis xsections data
    call get_photolysis_xs(photodata, photovars, err)
    if (len_trim(err) /= 0) return
    
    if (photodata%there_are_particles) then
      call get_aerosol_xs(photodata, photovars, err)
      if (len_trim(err) /= 0) return
    endif
    
  end subroutine
  
  ! Reads mie binary data file, then interpolates optical data to wavelength grid.
  ! also returns the radii of particles for that file.
  subroutine read_mie_data_file(filename, nw, wavl, &
                                 nrad_file, radii_file, w0_file, qext_file, g_file, err)
    use futils, only: addpnt, inter2
    
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nw
    real(real_kind), intent(in) :: wavl(nw+1)
    
    integer, intent(out) :: nrad_file
    real(real_kind), allocatable, intent(out) :: radii_file(:)
    real(real_kind), allocatable, intent(out) :: w0_file(:,:), qext_file(:,:), g_file(:,:)
    character(len=err_len), intent(out) :: err
    
    real(real_kind), allocatable :: wavl_tmp(:)
    real(real_kind), allocatable :: w0_tmp(:,:), qext_tmp(:,:), g_tmp(:,:)
    real(real_kind), allocatable :: temp_data(:), temp_wavelength(:)
    
    integer :: nw_tmp
    real(real_kind) :: dum
    integer :: i, j, io, ierr
    
    err = ''
    
    open(2,file=filename,form="unformatted",status='old',iostat=io)
    if (io /= 0) then
      err = "Was unable to open mie data file "//trim(filename)
      return
    endif
    
    read(2, iostat=io) nw_tmp
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    allocate(wavl_tmp(nw_tmp))
    read(2, iostat=io) wavl_tmp
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    
    read(2, iostat=io) nrad_file
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    allocate(radii_file(nrad_file))
    read(2, iostat=io) radii_file
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    
    allocate(w0_tmp(nw_tmp,nrad_file))
    allocate(qext_tmp(nw_tmp,nrad_file))
    allocate(g_tmp(nw_tmp,nrad_file))
    
    read(2, iostat=io) w0_tmp
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    read(2, iostat=io) qext_tmp
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    read(2, iostat=io) g_tmp
    if (io /= 0) then
      err = "Problem reading mie data file "//trim(filename)
      return
    endif
    
    ! this next read should be usuccessful
    read(2, iostat=io) dum
    if (io == 0) then
      err = "Problem reading mie data file "//trim(filename)// &
            ". Should have reached the end of the file, but did not."
      return
    endif
    close(2)
    
    ! now lets interpolate a bunch
    allocate(w0_file(nrad_file,nw))
    allocate(qext_file(nrad_file,nw))
    allocate(g_file(nrad_file,nw))
    
    allocate(temp_data(nw_tmp+2)) ! for interpolation
    allocate(temp_wavelength(nw_tmp+2))
    
    ! We do a constant extrapolation
    ! beyond the aerosol data, if it is necessary.
    ! If extrapolation happens, we print a warning
    ! if (wavl(1) < wavl_tmp(1)) then
    !   print*,"Warning: Constantly extrapolating aerosol data file"//trim(filename)// &
    !           " to lower wavelengths than there is data."
    ! endif
    ! if (wavl(nw+1) > wavl_tmp(nw_tmp)) then
    !   print*,"Warning: Constantly extrapolating aerosol data file"//trim(filename)// &
    !           " to larger wavelengths than there is data."
    ! endif
    do i = 1, nrad_file
      
      ! w0 (single scattering albedo)
      temp_data(1:nw_tmp) = w0_tmp(1:nw_tmp,i)
      temp_wavelength(1:nw_tmp) =  wavl_tmp(1:nw_tmp)
      j = nw_tmp
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, 0.d0, w0_tmp(1,i), ierr)
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, huge(1.d0), w0_tmp(nw_tmp,i), ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      call inter2(nw+1, wavl, w0_file(i,:), &
                  nw_tmp+2, temp_wavelength, temp_data, ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      
      ! qext (The extinction efficiency)
      temp_data(1:nw_tmp) = qext_tmp(1:nw_tmp,i)
      temp_wavelength(1:nw_tmp) =  wavl_tmp(1:nw_tmp)
      j = nw_tmp
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, 0.d0, qext_tmp(1,i), ierr)
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, huge(1.d0), qext_tmp(nw_tmp,i), ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      call inter2(nw+1, wavl, qext_file(i,:), &
                  nw_tmp+2, temp_wavelength, temp_data, ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      
      ! g (The scattering anisotropy or asymmetry factor)
      temp_data(1:nw_tmp) = g_tmp(1:nw_tmp,i)
      temp_wavelength(1:nw_tmp) =  wavl_tmp(1:nw_tmp)
      j = nw_tmp
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, 0.d0, g_tmp(1,i), ierr)
      call addpnt(temp_wavelength, temp_data, nw_tmp+2, j, huge(1.d0), g_tmp(nw_tmp,i), ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      call inter2(nw+1, wavl, g_file(i,:), &
                  nw_tmp+2, temp_wavelength, temp_data, ierr)
      if (ierr /= 0) then
        err = "Problems interpolating mie data from file "//trim(filename)// &
              " to the wavelength grid"
        return
      endif
      
    enddo
  
  end subroutine
  
  subroutine get_aerosol_xs(photodata, photovars, err)
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(in) :: photovars
    character(len=err_len), intent(out) :: err
    
    integer :: nrad
    integer, parameter :: nrad_fixed = 50
    real(real_kind), allocatable :: radii(:)
    real(real_kind), allocatable :: w0(:,:), qext(:,:), g(:,:)
    character(len=:), allocatable :: xsroot
    character(len=:), allocatable :: filename
    integer :: i
    
    xsroot = trim(photovars%data_dir)//"/aerosol_xsections/"
    
    allocate(photodata%radii_file(nrad_fixed,photodata%np))
    allocate(photodata%part_xs_file(photodata%np))
    
    photodata%nrad_file = nrad_fixed
    
    do i = 1,photodata%np
      
      if (photodata%particle_optical_prop(i) == 'none') then
        ! there is no optical data, so we skip
        photodata%part_xs_file(i)%ThereIsData = .false.
        cycle
      else
        ! there is optical data, so we allocate and get it
        photodata%part_xs_file(i)%ThereIsData = .true.
        allocate(photodata%part_xs_file(i)%w0(nrad_fixed,photodata%nw))
        allocate(photodata%part_xs_file(i)%qext(nrad_fixed,photodata%nw))
        allocate(photodata%part_xs_file(i)%gt(nrad_fixed,photodata%nw))
      endif
      
      if (photodata%particle_optical_type(i) == 0) then
        filename = xsroot//trim(photodata%particle_optical_prop(i))// &
                  "/mie_"//trim(photodata%particle_optical_prop(i))//".dat"
      elseif (photodata%particle_optical_type(i) == 1) then
        filename = xsroot//trim(photodata%particle_optical_prop(i))// &
                  "/frac_"//trim(photodata%particle_optical_prop(i))//".dat"
      endif
      
      if (allocated(radii)) then
        deallocate(radii, w0, qext, g)
      endif
      call read_mie_data_file(filename, photodata%nw, photodata%wavl, &
                              nrad, radii, w0, qext, g, err) 
      if (len_trim(err) /= 0) return
      if (nrad /= nrad_fixed) then
        err = "IOError: Aerosol data file "//filename// &
              "must have 50 radii bins."
        return
      endif
      
      photodata%radii_file(:,i) = radii/1.d4 ! convert from micron to cm
      photodata%part_xs_file(i)%w0 = w0
      photodata%part_xs_file(i)%qext = qext
      photodata%part_xs_file(i)%gt = g

    enddo
    
  end subroutine
  
  subroutine get_photolysis_xs(photodata, photovars, err)
    use futils, only: inter2, addpnt, replaceStr
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(in) :: photovars
    character(len=err_len), intent(out) :: err
    
    integer, parameter :: maxcols = 200
    character(len=:), allocatable :: xsroot
    character(len=:), allocatable :: filename, xsfilename, reaction
    character(len=str_len) :: line
    character(len=100) :: tmp(maxcols), tmp1
    real(real_kind), allocatable :: file_xs(:,:), file_qy(:,:), file_wav(:), file_wav_save(:), file_line(:)
    real(real_kind), allocatable :: dumby(:,:)
    real(real_kind), parameter :: rdelta = 1.d-4
    
    integer :: i, j, k, l, m, io, kk, ierr
    err = ''
    
    xsroot = trim(photovars%data_dir)//"/"//trim(photovars%xs_folder_name)//"/"
    
    allocate(photodata%xs_data(photodata%kj))
    
    do i = 1,photodata%kj
      filename = ''
      j = photodata%photonums(i)
      call reaction_string(photodata,j,reaction)
      filename = reaction
      filename = replaceStr(filename, ' ', '_')
      filename = replaceStr(filename, '>', '')
      filename = filename//'.txt'

      k = photodata%reactants_sp_inds(1,j)
      filename = trim(photodata%species_names(k))//'/'//filename
      xsfilename = trim(photodata%species_names(k))//'/'//trim(photodata%species_names(k))//'_xs.txt'
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
      photodata%xs_data(i)%n_temps = k - 2
      allocate(photodata%xs_data(i)%xs(k - 2, photodata%nw))
      allocate(photodata%xs_data(i)%xs_temps(k - 2))
      
      do k=1,photodata%xs_data(i)%n_temps
        tmp1 = tmp(k+1)
        read(tmp1(1:index(tmp1,'K')-1),*,iostat=io) photodata%xs_data(i)%xs_temps(k)
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
      allocate(file_wav_save(k+4))
      allocate(file_qy(k+4,photodata%xs_data(i)%n_temps))
      allocate(file_line(photodata%xs_data(i)%n_temps+1))
      allocate(dumby(photodata%nw,photodata%xs_data(i)%n_temps))

      rewind(101)
      read(101,*)
      read(101,*)
      do l = 1, k
        read(101,'(A)',iostat=io) line
        read(line,*) file_wav(l)
        do m = 1,photodata%xs_data(i)%n_temps+1
          read(line,*) file_line(1:m)
        enddo
        file_qy(l,:) = file_line(2:)
      enddo
      file_wav_save = file_wav
      
      ! interpolate to grid
      ierr = 0
      do l = 1, photodata%xs_data(i)%n_temps
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
        call inter2(photodata%nw+1,photodata%wavl,dumby(:,l), &
                    kk+4,file_wav,file_qy(:,l),ierr)
        if (ierr /= 0) then
          err = 'IOError: Problem interpolating quantum yield data to photolysis grid for reaction '// &
                trim(reaction)
          return
        endif
        photodata%xs_data(i)%xs(l,:) = dumby(:,l)
        k = kk
        file_wav = file_wav_save
      enddo
      
      close(101)
      deallocate(file_wav,file_wav_save,file_qy)
      
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
      allocate(file_wav_save(k+4))
      allocate(file_xs(k+4,photodata%xs_data(i)%n_temps))
      
      rewind(102)
      read(102,*)
      read(102,*)
      do l = 1, k
        read(102,'(A)',iostat=io) line
        read(line,*) file_wav(l)
        do m = 1,photodata%xs_data(i)%n_temps+1
          read(line,*) file_line(1:m)
        enddo
        file_xs(l,:) = file_line(2:)
      enddo
      file_wav_save = file_wav
      
      ierr = 0
      do l = 1, photodata%xs_data(i)%n_temps
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
        call inter2(photodata%nw+1,photodata%wavl,dumby(:,l), &
                    kk+4,file_wav,file_xs(:,l),ierr)  
        if (ierr /= 0) then
          err = 'IOError: Problem interpolating xs data to photolysis grid for reaction '// &
                trim(reaction)
          return
        endif
        photodata%xs_data(i)%xs(l,:) = photodata%xs_data(i)%xs(l,:)*dumby(:,l)
        k = kk
        file_wav = file_wav_save
      enddo

      close(102)
      deallocate(file_xs, file_wav, file_wav_save, file_line, dumby)
    enddo
    
  end subroutine
  
  
  subroutine get_rayleigh(photodata, photovars, err)
    use yaml, only : parse, error_length
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(in) :: photovars
    character(len=err_len), intent(out) :: err
    
    real(real_kind), allocatable :: A(:), B(:), Delta(:)
    character(len=str_len) :: rayleigh_file

    character(error_length) :: error
    class (type_node), pointer :: root
    integer :: i, j
    err = ''
    
    rayleigh_file = trim(photovars%data_dir)//"/rayleigh/rayleigh.yaml"
    
    ! parse yaml file
    root => parse(rayleigh_file,error=error)
    if (error /= '') then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        call rayleigh_params(root,photodata,trim(rayleigh_file),err, &
                             photodata%raynums, A, B, Delta)
    end select
    call root%finalize()
    deallocate(root)
    if (len_trim(err) /= 0) return
    
    ! compute cross sections
    photodata%nray = size(A)
    allocate(photodata%sigray(photodata%nray,photodata%nw))
    do i = 1,photodata%nw
      do j = 1,size(A)
        call rayleigh_vardavas(A(j), B(j), Delta(j), photodata%wavl(i), &
                               photodata%sigray(j, i))
      enddo
    enddo
    deallocate(A,B,Delta)
  end subroutine
  
  subroutine rayleigh_params(mapping,photodata,infile,err, raynums, A, B, Delta)
    class (type_dictionary), intent(in), pointer :: mapping
    type(PhotochemData), intent(in) :: photodata
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
      ind = findloc(photodata%species_names,trim(key_value_pair%key))
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
      ind = findloc(photodata%species_names,trim(key_value_pair%key))
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
    use futils, only: inter2, addpnt
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
  
  
  subroutine read_atmosphere_file(atmosphere_txt, photodata, photovars, err)
    character(len=*), intent(in) :: atmosphere_txt
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(inout) :: photovars
    character(len=err_len), intent(out) :: err
    
    character(len=10000) :: line
    character(len=:), allocatable :: message
    character(len=s_str_len) :: arr1(1000)
    character(len=s_str_len) :: arr11(1000)
    character(len=s_str_len),allocatable, dimension(:) :: labels
    integer :: ind(1)
    real(real_kind), allocatable :: temp(:,:)
    integer :: io, i, n, nn, ii
    logical :: missing
    
    err = ''
    open(4, file=trim(atmosphere_txt),status='old',iostat=io)
    if (io /= 0) then
      err = 'Can not open file '//trim(atmosphere_txt)
      return
    endif
    read(4,'(A)') line
    
    photodata%nzf = -1
    io = 0
    do while (io == 0)
      read(4,*,iostat=io)
      photodata%nzf = photodata%nzf + 1
    enddo
    
    allocate(photodata%z_file(photodata%nzf))
    allocate(photodata%T_file(photodata%nzf))
    allocate(photodata%edd_file(photodata%nzf))
    allocate(photodata%usol_file(photodata%nq, photodata%nzf))
    photodata%z_file = 0.d0
    photodata%T_file = 0.d0
    photodata%edd_file = 0.d0
    photodata%usol_file = 1.d-40
    if (photodata%there_are_particles) then
      allocate(photodata%particle_radius_file(photodata%npq, photodata%nzf))
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
    
    ! allocate memory
    allocate(labels(n))
    allocate(temp(n,photodata%nzf))
    rewind(4)
    read(4,'(A)') line
    read(line,*) (labels(i),i=1,n)
    
    ! First read in all the data into big array
    do i = 1,photodata%nzf
      read(4,*,iostat=io) (temp(ii,i),ii=1,n)
      if (io /= 0) then
        err = 'Problem reading in initial atmosphere in '//trim(atmosphere_txt)
        return
      endif
    enddo
    close(4)
    
    ! reads in mixing ratios
    missing = .false.
    do i=1,photodata%nq
      ind = findloc(labels,photodata%species_names(i))
      if (ind(1) /= 0) then
        photodata%usol_file(i,:) = temp(ind(1),:)
      else
        missing = .true.
      endif
    enddo
    
    photovars%no_water_profile = .false.
    if (photodata%fix_water_in_trop) then
      ind = findloc(labels,'H2O')
      if (ind(1) == 0) then
        photovars%no_water_profile = .true. 
      endif
    endif
    
    if (missing) then
      message = 'Warning: Did not find initial data for some species in '// &
                trim(atmosphere_txt)//' . The program will assume initial mixing ratios of 1.0e-40'
      if (photovars%no_water_profile) then
        message = message // " except H2O, which will be set to saturation in troposphere with constant "//&
                              "extrapolation above the tropopause."
      endif
      print*,message
    endif
    
    if (photodata%there_are_particles) then
      missing = .false.
      do i=1,photodata%npq
        ind = findloc(labels,trim(photodata%species_names(i))//"_r")
        if (ind(1) /= 0) then
          photodata%particle_radius_file(i,:) = temp(ind(1),:)
        else
          ! did not find the data
          ! will set to 0.1 micron
          photodata%particle_radius_file(i,:) = 1.d-5
          missing = .true.
        endif
      enddo
      
      if (missing) then
        print*,'Warning: Did not find particle radii for some species in '//&
                trim(atmosphere_txt)//' . The program will assume 0.1 micron raddii.'
      endif
    endif
    
    ! reads in temperature
    ind = findloc(labels,'temp')
    if (ind(1) /= 0) then
      photodata%T_file(:) = temp(ind(1),:)
    else
      err = '"temp" was not found in input file '//trim(atmosphere_txt)
      return
    endif
    
    ! reads in alt
    ind = findloc(labels,'alt')
    if (ind(1) /= 0) then
      photodata%z_file(:) = temp(ind(1),:)*1.d5 ! conver to cm
    else
      err = '"alt" was not found in input file '//trim(atmosphere_txt)
      return
    endif
    
    ! reads in eddy diffusion
    ind = findloc(labels,'eddy')
    if (ind(1) /= 0) then
      photodata%edd_file(:) = temp(ind(1),:)
    else
      err = '"eddy" was not found in input file '//trim(atmosphere_txt)
      return
    endif

  end subroutine
  
end submodule