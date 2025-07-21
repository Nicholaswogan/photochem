submodule (photochem_input) photochem_input_read
  implicit none
  
contains
  
  module subroutine read_all_files(mechanism_file, s, flux_file, atmosphere_txt, dat, var, err)
    character(len=*), intent(in) :: mechanism_file
    type(PhotoSettings), intent(in) :: s
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    ! stuff dat needs before entering get_photomech
    dat%nsl = s%nsl
    dat%SL_names = s%SL_names

    call get_photomech(mechanism_file, dat, var, err)
    if (allocated(err)) return

    ! initial atmosphere
    call read_atmosphere_file(atmosphere_txt, dat, var, err)
    if (allocated(err)) return
    
    call unpack_settings(s%filename, s, dat, var, err)
    if (allocated(err)) return

    !!! henrys law !!!
    if (dat%gas_rainout) then
      call get_henry(dat, var, s, err)
      if (allocated(err)) return
    endif
    
    call get_photorad(dat, var, err)
    if (allocated(err)) return
    
    ! stellar flux
    allocate(var%photon_flux(dat%nw))
    call read_stellar_flux(flux_file, dat%nw, dat%wavl, var%photon_flux, err)
    if (allocated(err)) return
    
  end subroutine
  
  subroutine get_photomech(infile, dat, var, err) 
    use fortran_yaml_c, only : YamlFile
    character(len=*), intent(in) :: infile
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    character(:), allocatable, intent(out) :: err
    
    type(YamlFile) :: file

    ! parse yaml file
    call file%parse(infile, err)
    if (allocated(err)) return
    select type (root => file%root)
      class is (type_dictionary)
        call get_rxmechanism(root, infile, dat, var, err)
      class default
        err = "yaml file must have dictionaries at root level"
    end select
    call file%finalize()
    if (allocated(err)) return
     
  end subroutine
  
  subroutine get_rxmechanism(mapping, infile, dat, var, err)
    use photochem_types, only: ReverseRate
    use photochem_enum, only: CondensingParticle, ReactionParticle
    use photochem_enum, only: MieParticle, FractalParticle 
    use photochem_enum, only: PhotolysisRateType, ReverseRateType
    use clima_saturationdata, only: SaturationData
    class (type_dictionary), intent(in), pointer :: mapping
    character(len=*), intent(in) :: infile
    type(PhotochemData), target, intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    character(:), allocatable, intent(out) :: err
    
    class (type_dictionary), pointer :: sat_params
    class (type_list), pointer :: species, reactions, atoms
    class (type_list), pointer :: particles
    type (type_error), allocatable :: io_err
    class (type_list_item), pointer :: item
    class (type_dictionary), pointer :: dict
    class (type_key_value_pair), pointer :: key_value_pair

    ! temporary work variables
    character(len=str_len) :: tmpchar
    character(len=str_len) :: tmp
    character(len=:), allocatable :: rxstring
    integer :: i, ii, j, k, kk, l, ind(1)
    logical :: reverse
    ! all_species causes a small memory leak. Not sure how to free the memory properly
    type(type_list_tmp) :: all_species, all_reactions ! will include particles
    logical, allocatable :: duplicate(:)


    atoms => mapping%get_list('atoms',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    species => mapping%get_list('species',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    reactions => mapping%get_list('reactions',.true.,error = io_err) 
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
    
    ! should i reverse reactions?
    dat%reverse = mapping%get_logical('reverse-reactions',.true.,error = io_err)
    
    !!! atoms !!!
    dat%natoms = atoms%size()
    allocate(dat%atoms_names(dat%natoms))
    allocate(dat%atoms_mass(dat%natoms))
    allocate(dat%atoms_redox(dat%natoms))
    
    j = 1
    item => atoms%first
    do while(associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        dat%atoms_names(j) = element%get_string("name",error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        dat%atoms_mass(j) = element%get_real("mass",error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        dat%atoms_redox(j) = element%get_real("redox",error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
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
      dat%there_are_particles = .true.
      dat%np = 0
      item => particles%first
      do while (associated(item))
        dat%np = dat%np + 1
        item => item%next
      enddo
      
      allocate(dat%particle_names(dat%np))
      allocate(dat%particle_formation_method(dat%np))
      allocate(dat%particle_density(dat%np))
      allocate(dat%particle_sat(dat%np))
      allocate(dat%particle_gas_phase(dat%np))
      allocate(dat%particle_optical_prop(dat%np))
      allocate(dat%particle_optical_type(dat%np))
      
      item => particles%first
      j = 1
      do while (associated(item))
        call all_species%append(item%node)
        select type (element => item%node)
        class is (type_dictionary)
          dat%particle_names(j) = element%get_string("name",error = io_err) ! get name
          if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          
          tmpchar = element%get_string("formation",error = io_err) 
          if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          if (trim(tmpchar) == 'saturation') then
            dat%particle_formation_method(j) = CondensingParticle
          elseif (trim(tmpchar) == 'reaction') then
            dat%particle_formation_method(j) = ReactionParticle
          else
            err = "IOError: the only formation mechanism for particles is 'saturation'"
            return
          endif
          dat%particle_density(j) = element%get_real("density",error = io_err)
          if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          dat%particle_optical_prop(j) = element%get_string("optical-properties",error = io_err)
          if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
          ! only require optical type if not "none"
          if (dat%particle_optical_prop(j) /= 'none') then
            tmpchar = element%get_string("optical-type",error = io_err)
            if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
            if (trim(tmpchar) == "mie") then
              dat%particle_optical_type(j) = MieParticle
            elseif  (trim(tmpchar) == "fractal") then
              err = "IOError: 'fractal' is not an optional optical type for "// &
                    trim(dat%particle_names(j))
              return
            else
              err = "IOError: "//trim(tmpchar)//" is not an optional optical type for "// &
                    trim(dat%particle_names(j))
              return
            endif
          endif
  
          if (dat%particle_formation_method(j) == CondensingParticle) then

            dat%particle_gas_phase(j) = element%get_string("gas-phase",error = io_err) 
            if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif

            sat_params => element%get_dictionary('saturation',.true.,error=io_err)
            if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif

            dat%particle_sat(j) = SaturationData(sat_params, trim(dat%particle_names(j)), trim(infile), err)
            if (allocated(err)) return
            
          elseif (dat%particle_formation_method(j) == ReactionParticle) then
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

      block
        character(s_str_len), allocatable :: str_list(:), str_list1(:)
        ! Check for duplicate gas phases species
        str_list = ['']
        str_list1 = ['']
        do j = 1,dat%np
          if (dat%particle_formation_method(j) == CondensingParticle) then
            str_list = [str_list, dat%particle_gas_phase(j)]
            str_list1 = [str_list1, dat%particle_names(j)]
          endif
        enddo
        str_list = str_list(2:size(str_list))
        str_list1 = str_list1(2:size(str_list1))
        i = check_for_duplicates(str_list)
        if (i /= 0) then
          err = 'Particle "'//trim(str_list1(i))//'" has gas phase species "'//trim(str_list(i))// &
                '"  which is duplicate with another particle.'
          return
        endif
      endblock

    else ! there are no particles
      dat%there_are_particles = .false.
      dat%np = 0
    endif
    
    ! for now number particle equations will be the same 
    ! as number of particles
    dat%npq = dat%np
    
    !!! done with particles !!!
    
    !!! species !!!
    dat%ng = 0 ! count number of gas phase species
    item => species%first
    do while (associated(item))
      item => item%next
      dat%ng = dat%ng + 1
    enddo
    
    dat%nll = dat%ng - dat%nsl
    
    dat%ng_1 = dat%npq + 1 ! the long lived gas index
    ! dat%nq is the last ll gas index
    
    ! now we now nq, the number of PDEs
    dat%nq = dat%npq + dat%nll
    
    ! we also now nsp, the index of the backgorund gas 
    dat%nsp = dat%npq + dat%ng
    
    ! species_mass, species_composition, and species_names
    ! will include the particles, thus we allocate nsp
    allocate(dat%species_redox(dat%nsp))
    allocate(dat%species_mass(dat%nsp))
    allocate(dat%species_composition(dat%natoms,dat%nsp+2))
    dat%species_composition = 0
    allocate(dat%species_names(dat%nsp+2))
    dat%species_names(dat%nsp+1) = "hv" ! always add these guys
    dat%species_names(dat%nsp+2) = "M"
    ! we will not include particles in thermodynamic data.
    if (dat%reverse) then
      allocate(dat%thermo_data(dat%ng))
    endif
    
    ! Append the species to the end of a list
    ! which has particles in the beginning
    item => species%first
    do while (associated(item))
      call all_species%append(item%node)
      item => item%next
    enddo

    ! Loop through particles and gases
    kk = dat%ng_1
    l = 1
    ii = 1 ! overall counter
    item => all_species%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmpchar = trim(element%get_string("name",error = io_err)) ! get name
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        
        if (ii < dat%ng_1) then
          ! we are dealing with particles
          j = ii
        else
          ! we are dealing with gases
          ind = findloc(dat%SL_names,tmpchar)
          if (ind(1) /= 0) then ! short lived species
            j = dat%nq + l 
            l = l + 1
          else ! long lived species
            j = kk
            kk = kk + 1
          endif
        endif
                  
        dat%species_names(j) = tmpchar
        dict => element%get_dictionary("composition",.true.,error = io_err)  ! get composition
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        key_value_pair => dict%first ! dont allow unspecified atoms
        do while (associated(key_value_pair))
          ind = findloc(dat%atoms_names,trim(key_value_pair%key))
          if (ind(1) == 0) then
            err = 'IOError: The atom "'// trim(key_value_pair%key)// '" is not in the list of atoms.'
            return
          endif
          key_value_pair =>key_value_pair%next
        enddo
        
        do i=1,dat%natoms
          dat%species_composition(i,j) =  &
              dict%get_integer(dat%atoms_names(i),0,error = io_err)
          if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        enddo
        dat%species_mass(j) = sum(dat%species_composition(:,j) * dat%atoms_mass)
        dat%species_redox(j) = sum(dat%species_composition(:,j) * dat%atoms_redox)
        
        if (dat%reverse .and. (ii >= dat%ng_1)) then
          call get_thermodata(element,dat%species_names(j), infile, &
                              dat%thermo_data(j-dat%npq), err)
          if (allocated(err)) return
        endif
      class default
        err = "IOError: Problem with species number "//char(j)//"  in the input file"
        return
      end select
      ii = ii + 1
      item => item%next
    enddo
    
    if (l-1 /= dat%nsl) then
      err = 'IOError: One of the short lived species is not in the file '//trim(infile)
      return
    endif
    i = check_for_duplicates(dat%species_names)
    if (i /= 0) then
      err = 'Species "'//trim(dat%species_names(i))//'" is a duplicate in '//trim(infile)
      return
    endif
    !!! done with species !!!

    allocate(dat%gas_particle_ind(dat%nq))
    dat%gas_particle_ind = 0
    if (dat%there_are_particles) then
      ! get indexes of gas phase condensing species
      allocate(dat%particle_gas_phase_ind(dat%np))
      do i = 1,dat%np
        if (dat%particle_formation_method(i) == CondensingParticle) then
          ! if a condensing molecule
          ind = findloc(dat%species_names,dat%particle_gas_phase(i))
          if (ind(1) == 0) then
            err = "IOError: particle "//trim(dat%particle_names(i))// &
                  " can not be made from "//trim(dat%particle_gas_phase(i))// &
                  " because "//trim(dat%particle_gas_phase(i))//" is not a gas"// &
                  " in the model."
            return
          elseif (ind(1) > dat%nq .or. ind(1) < dat%ng_1) then
            err = "Particle "//trim(dat%particle_names(i))// &
                  " can not be made from "//trim(dat%particle_gas_phase(i))// &
                  " because "//trim(dat%particle_gas_phase(i))//" is"// &
                  " is a particle or is short-lived species."
            return
          else
            dat%particle_gas_phase_ind(i) = ind(1)
          endif

          ! Index of particle
          dat%gas_particle_ind(ind(1)) = i
        endif
      enddo
    endif
    
    !!! reactions !!!
    item => reactions%first
    do while (associated(item))
      call all_reactions%append(item%node)
      item => item%next
    enddo

    dat%nrF = all_reactions%size()
    dat%nrR = 0
    item => all_reactions%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmp = trim(element%get_string("equation",error = io_err))
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        reverse = is_rx_reverse(tmp, err)
        if (allocated(err)) return 
        
        if (reverse) then
          if (.not.dat%reverse) then
            err = 'IOError: reaction file '//trim(infile)//' contains reverse reaction '//tmp// &
                  ', which is incompatible with "reverse-reactions: false"'
            return
          endif
          dat%nrR = dat%nrR + 1
        endif

      class default
        err = "IOError: Problem with reaction number "//char(j)//" in the input file."
        return
      end select
      item => item%next
    enddo
    dat%nrT = dat%nrR + dat%nrF
    
    ! allocate stuff and loop through reactions again
    allocate(duplicate(dat%nrT))
    allocate(dat%rx(dat%nrT))
    
    j = 1
    k = 1
    item => all_reactions%first
    do while (associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        tmp = trim(element%get_string("equation",error = io_err))
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        call get_rateparams(dat, element, infile, dat%rx(j), err)
        if (allocated(err)) return
        call get_reaction_sp_nums(dat, tmp, dat%rx(j), reverse, err)
        if (allocated(err)) return
        
        ! check if duplicate
        duplicate(j) = element%get_logical("duplicate",default=.false.,error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        if (reverse) then
          ! reaction has a reverse
          i = dat%nrF + k
          duplicate(i) = duplicate(j)
          allocate(dat%rx(j)%reverse_info)
          allocate(dat%rx(i)%reverse_info)
          dat%rx(j)%reverse_info = i
          dat%rx(i)%reverse_info = j
          dat%rx(i)%nreact = dat%rx(j)%nprod
          dat%rx(i)%nprod = dat%rx(j)%nreact
          allocate(dat%rx(i)%react_sp_inds(dat%rx(i)%nreact))
          allocate(dat%rx(i)%prod_sp_inds(dat%rx(i)%nprod))
          dat%rx(i)%react_sp_inds = dat%rx(j)%prod_sp_inds
          dat%rx(i)%prod_sp_inds = dat%rx(j)%react_sp_inds
          
          allocate(ReverseRate::dat%rx(i)%rp)
          dat%rx(i)%rp%rxtype = ReverseRateType

          k = k + 1
        endif
      class default
        err = "IOError: Problem with reaction number "//char(j)//" in the input file."
        return
      end select
      item => item%next
      j = j + 1
    enddo
    
    ! production and loss mechanisms for each species
    allocate(dat%pl(dat%nsp))
    do i = 1,dat%nsp
      dat%pl(i)%nump = 0
      dat%pl(i)%numl = 0
    enddo
    do j = 1,dat%nrT
      do i = 1,dat%rx(j)%nreact
        kk = dat%rx(j)%react_sp_inds(i)
        if (kk <= dat%nsp) then
          dat%pl(kk)%numl = dat%pl(kk)%numl + 1
        endif
      enddo
      do i = 1,dat%rx(j)%nprod
        kk = dat%rx(j)%prod_sp_inds(i)
        if (kk <= dat%nsp) then
          dat%pl(kk)%nump = dat%pl(kk)%nump + 1
        endif
      enddo
    enddo
    do i = 1,dat%nsp
      allocate(dat%pl(i)%iprod(dat%pl(i)%nump))
      allocate(dat%pl(i)%iloss(dat%pl(i)%numl))
      dat%pl(i)%nump = 0
      dat%pl(i)%numl = 0
    enddo

    do j = 1,dat%nrT
      do i = 1,dat%rx(j)%nreact
        kk = dat%rx(j)%react_sp_inds(i)
        if (kk <= dat%nsp) then
          dat%pl(kk)%numl = dat%pl(kk)%numl + 1
          l = dat%pl(kk)%numl
          dat%pl(kk)%iloss(l) = j
        endif
      enddo
      do i = 1,dat%rx(j)%nprod
        kk = dat%rx(j)%prod_sp_inds(i)
        if (kk <= dat%nsp) then
          dat%pl(kk)%nump = dat%pl(kk)%nump + 1
          l = dat%pl(kk)%nump
          dat%pl(kk)%iprod(l) = j
        endif
      enddo
    enddo
    
    ! photolysis
    dat%kj = 0
    do i = 1, dat%nrF
      if (dat%rx(i)%rp%rxtype == PhotolysisRateType) then
        dat%kj = dat%kj + 1
      endif
    enddo
    allocate(dat%photonums(dat%kj))
    j = 1
    do i = 1, dat%nrF
      if (dat%rx(i)%rp%rxtype == PhotolysisRateType) then
        dat%photonums(j) = i
        j = j + 1
      endif
    enddo
    
    ! save reaction names
    allocate(dat%reaction_equations(dat%nrT))
    do i = 1,dat%nrT
      call reaction_string(dat,i,rxstring)
      dat%reaction_equations(i) = rxstring
    enddo
    
    call check_for_reaction_duplicates(dat, duplicate, err)
    if (allocated(err)) return
    
    ! Make sure particles are not being destroyed from reactions
    do i = 1,dat%np
      if (dat%pl(i)%numl /= 0) then
        err = 'Particle "'//trim(dat%species_names(i))//'" is destroyed by reactions. '// &
              'Particles can not be destroyed by reactions! Remove these reaction(s).'
        return
      endif
    enddo
    
    !!! end reactions !!!
    
  end subroutine
  
  subroutine unpack_settings(infile, s, dat, var, err)
    use photochem_enum, only: CondensingParticle
    use photochem_enum, only: VelocityBC, DensityBC, PressureBC
    use photochem_enum, only: DiffusionLimHydrogenEscape, ZahnleHydrogenEscape, NoHydrogenEscape
    use photochem_types, only: PhotoSettings
    character(len=*), intent(in) :: infile
    type(PhotoSettings), intent(in) :: s
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    integer :: j, i, ind(1), io
    
    !!!!!!!!!!!!!!!!!!!!!!!
    !!! atmosphere-grid !!!
    !!!!!!!!!!!!!!!!!!!!!!!
    var%bottom_atmos = s%bottom

    read(s%top,*,iostat = io) var%top_atmos
    if (io /= 0) then
      ! it isn't a float
      if (trim(s%top) == "atmospherefile") then
        ! Use the TOA in the initialization file
        var%top_atmos = dat%z_file(dat%nzf) + 0.5_dp*(dat%z_file(2) - dat%z_file(1))
      else
        err = '"top" can only be a number or "atmospherefile". See '//trim(infile)
        return 
      endif
    endif
    if (var%top_atmos < var%bottom_atmos) then
      err = 'The top of the atmosphere must be bigger than the bottom'
      return
    endif
    var%nz = s%nz
    
    !!!!!!!!!!!!!!
    !!! planet !!!
    !!!!!!!!!!!!!!
    dat%planet_mass = s%planet_mass
    dat%planet_radius = s%planet_radius
    var%surface_albedo = s%surface_albedo
    var%photon_scale_factor = s%photon_scale_factor
    var%solar_zenith = s%solar_zenith
    dat%H_escape_type = s%H_escape_type
    if (dat%H_escape_type /= NoHydrogenEscape) then

      ! make sure H2 and H are in the model
      ind = findloc(dat%species_names,'H2')
      dat%LH2 = ind(1)
      if (ind(1) == 0) then
        err = 'IOError: H2 must be a species if hydrogen-escape is on.'
        return
      endif
      ind = findloc(dat%species_names,'H')
      dat%LH = ind(1)
      if (ind(1) == 0) then
        err = 'IOError: H must be a species if hydrogen-escape is on.'
        return
      endif

      ! make sure H2 and H are not short-lived
      ind = findloc(dat%species_names(dat%nq+1:dat%nq+dat%nsl),'H2')
      if (ind(1) /= 0) then
        err = 'H2 must not be short lived if hydrogen escape is on'
        return
      endif
      ind = findloc(dat%species_names(dat%nq+1:dat%nq+dat%nsl),'H')
      if (ind(1) /= 0) then
        err = 'H must not be short lived if hydrogen escape is on'
        return
      endif

      if (dat%H_escape_type == ZahnleHydrogenEscape) then; block
        use photochem_eqns, only: zahnle_Hescape_coeff

        allocate(dat%H_escape_coeff)
        dat%H_escape_coeff = zahnle_Hescape_coeff(s%H_escape_S1)

      endblock; endif

    endif
    
    ! default-gas-lower-boundary already applied to PhotoSettings

    ! conserving initialization?
    dat%conserving_init = s%conserving_init
    
    ! water
    dat%fix_water_in_trop = s%fix_water_in_trop
    dat%water_cond = s%water_cond
    dat%gas_rainout = s%gas_rainout
    ind = findloc(dat%species_names,'H2O')
    dat%LH2O = ind(1)
    if (ind(1) == 0 .and. dat%fix_water_in_trop) then
      err = 'IOError: H2O must be a species if fix-water-in-troposphere = True.'
      return
    elseif (ind(1) == 0 .and. dat%water_cond) then
      err = 'IOError: H2O must be a species if water-condensation = True.'
      return
    elseif (ind(1) == 0 .and. dat%gas_rainout) then
      err = 'IOError: H2O must be a species if gas-rainout = True.'
      return
    endif
    if ((dat%fix_water_in_trop .or. dat%water_cond) .and. dat%there_are_particles) then
      ! Make sure there isn't already H2O condensation implemented by a particle.
      do i = 1,dat%np
        if (dat%particle_formation_method(i) == CondensingParticle) then
          if (dat%particle_gas_phase(i) == 'H2O') then
            err = 'IOError: Either "fix-water-in-troposphere" or "water-condensation" is turned on'// &
            ' in the settings file, but the reaction mechanism already implements H2O condensation'// &
            ' via a particle.'
            return
          endif
        endif
      enddo
    endif
    
    if (dat%fix_water_in_trop) then  
      
      read(s%relative_humidity,*,iostat = io) var%relative_humidity
      
      if (io /= 0) then
        ! it isn't a float
        if (trim(s%relative_humidity) == "manabe") then
          var%use_manabe = .true.
          var%relative_humidity = -1.0
        else
          err = '"relative-humidity" can only be a number between 0 and 1, or "manabe". See '//trim(infile)
          return 
        endif
      else
        var%use_manabe = .false.
      endif
      
    endif
    
    if (dat%gas_rainout) then
      var%rainfall_rate = s%rainfall_rate
    endif
    
    if (dat%fix_water_in_trop .or. dat%gas_rainout) then
      ! we need a tropopause altitude
      var%trop_alt = s%trop_alt
      if (var%trop_alt > var%top_atmos) then
        err = 'IOError: tropopause-altitude must be between the top and bottom of the atmosphere'
        return
      endif
    endif
    
    !!!!!!!!!!!!!!!!!
    !!! particles !!!
    !!!!!!!!!!!!!!!!!
    ! Condensation rate parameters. size is zero when there are no particles
    allocate(var%cond_params(dat%np))
    ! Replace default values with values from settings file, if needed
    if (allocated(s%particles) .and. dat%there_are_particles) then
      do i = 1,size(s%particles)
        ind = findloc(dat%species_names(1:dat%np),s%particles(i)%name)
        if (ind(1) == 0) then
          err = 'Particle '//s%particles(i)%name//' in the settings file is not a particle in the reaction file.'
          return
        endif
        var%cond_params(ind(1)) = s%particles(i)%params
      enddo
    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! boundary-conditions !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(var%lowerboundcond(dat%nq))
    allocate(var%lower_vdep(dat%nq))
    allocate(var%lower_flux(dat%nq))
    allocate(var%lower_dist_height(dat%nq))
    allocate(var%lower_fix_den(dat%nq))
    allocate(var%lower_fix_press(dat%nq))
    allocate(var%upperboundcond(dat%nq))
    allocate(var%upper_veff(dat%nq))
    allocate(var%upper_flux(dat%nq))

    allocate(var%only_eddy(dat%nq))
    allocate(var%rate_fcns(dat%nq))
    ! default boundary conditions
    var%lowerboundcond(:dat%np) = VelocityBC ! default particle BC is alway velocity
    var%lowerboundcond(dat%ng_1:) = s%default_lowerboundcond ! can be -1 (Moses) or 0 (velocity)
    var%lower_vdep = 0.0_dp
    var%upperboundcond = VelocityBC
    var%upper_veff = 0.0_dp
    var%only_eddy = .false.
    
    do j = 1,size(s%ubcs)
      ! make sure it isn't water
      if (dat%fix_water_in_trop .and. s%sp_names(j) == "H2O") then
        err = "IOError: H2O can not have a specified boundary condition"// &
              " if water-saturated-troposphere = true in the settings file."
        return
      endif
      ! check if in rxmech
      ind = findloc(dat%species_names,s%sp_names(j))
      if (ind(1) == 0) then
        err = "IOError: Species "//trim(s%sp_names(j))// &
        ' in settings file is not in the reaction mechanism file.'
        return 
      endif
      
      if (s%sp_types(j) == 'long lived') then
      
        var%lowerboundcond(ind(1)) = s%lbcs(j)%bc_type
        var%lower_vdep(ind(1)) = s%lbcs(j)%vel
        var%lower_flux(ind(1)) = s%lbcs(j)%flux
        var%lower_dist_height(ind(1)) = s%lbcs(j)%height
        var%lower_fix_den(ind(1)) = s%lbcs(j)%den
        var%lower_fix_press(ind(1)) = s%lbcs(j)%press
        
        var%upperboundcond(ind(1)) = s%ubcs(j)%bc_type
        var%upper_veff(ind(1)) = s%ubcs(j)%vel
        var%upper_flux(ind(1)) = s%ubcs(j)%flux
        
        var%only_eddy(ind(1)) = s%only_eddy(j)
        
      endif
      
    enddo
    
    ! Make sure that upper boundary condition for H and H2 are
    ! effusion velocities, if diffusion limited escape
    if (dat%H_escape_type == DiffusionLimHydrogenEscape) then
      if (var%upperboundcond(dat%LH2) /= VelocityBC) then
        err = "IOError: H2 must have a have a effusion velocity upper boundary"// &
              " for diffusion limited hydrogen escape"
        return
      endif
      if (var%upperboundcond(dat%LH) /= VelocityBC) then
        err = "IOError: H must have a have a effusion velocity upper boundary"// &
              " for diffusion limited hydrogen escape"
        return
      endif
    endif

    ! Make sure all lower boundary conditions for particles are deposition
    ! velocities, so that particles actually fall out of the model.
    if (dat%there_are_particles) then
      do i = 1,dat%npq
        if (var%lowerboundcond(i) /= VelocityBC) then
          err = 'Particle "'//trim(dat%species_names(i))//'" must have deposition velocity '// &
                'lower boundary condition.'
          return
        endif
      enddo
    endif

    ! check for SL nonlinearities
    call check_sl(dat, err)
    if (allocated(err)) return
    
  end subroutine
  
  subroutine get_henry_parse(root, dat, var, henry_names, henry_data, err)
    class (type_list), intent(in) :: root
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    character(len=s_str_len), allocatable, intent(out) :: henry_names(:)
    real(dp), allocatable, intent(out) :: henry_data(:,:)
    character(:), allocatable :: err
  
    type (type_error), allocatable :: io_err
    class (type_list_item), pointer :: item
    integer :: j  
    
    j = root%size()
    
    allocate(henry_names(j))
    allocate(henry_data(2,j))
    j = 1
    item => root%first
    do while(associated(item))
      select type (element => item%node)
      class is (type_dictionary)
        henry_names(j) = element%get_string('name',error = io_err)
        if (allocated(io_err)) then; err = trim(io_err%message); return; endif
        henry_data(1,j) = element%get_real('A',error = io_err)
        if (allocated(io_err)) then; err = trim(io_err%message); return; endif
        henry_data(2,j) = element%get_real('B',error = io_err)
        if (allocated(io_err)) then; err = trim(io_err%message); return; endif
      
      j = j + 1
      item => item%next
      end select
    enddo
  end subroutine
  
  subroutine get_henry(dat, var, s, err)
    use photochem_types, only: PhotoSettings
    use fortran_yaml_c, only : YamlFile
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    type(PhotoSettings), intent(in) :: s
    character(:), allocatable :: err
    
    type(YamlFile) :: file
    integer :: j, ind, ind1, i
    
    character(len=s_str_len), allocatable :: henry_names(:)
    real(dp), allocatable :: henry_data(:,:)

    ! parse yaml file
    call file%parse(trim(var%data_dir)//"/henry/henry.yaml", err)
    if (allocated(err)) return
    select type (root => file%root)
    class is (type_list)
      call get_henry_parse(root, dat, var, henry_names, henry_data, err)
    class default
      err = "yaml file must have dictionaries at root level"
    end select
    call file%finalize()
    if (allocated(err)) return

    allocate(dat%henry_data(2,dat%nsp))
    dat%henry_data = 0.0_dp

    if (allocated(s%rainout_species)) then
      do j = 1,size(s%rainout_species)
        ind = findloc(dat%species_names, s%rainout_species(j), 1)
        if (ind == 0) then
          err = 'Rainout species "'//trim(s%rainout_species(j))//'" is not in the list of species.'
          return
        endif

        ! particle or gas?
        if (ind < dat%ng_1) then
          ! particle
          dat%henry_data(1,ind) = 7.e11_dp
        else
          ! gas
          ! look for the gas in the henry data
          ind1 = findloc(henry_names, s%rainout_species(j), 1)
          if (ind1 == 0) then
            err = 'No solubility data exits for rainout species "'//trim(s%rainout_species(j))//'"'
            return
          endif

          dat%henry_data(:,ind) = henry_data(:,ind1)
        endif

      enddo
    else
      ! we do all possible rainout species
      do j = 1,size(henry_names)
        ind = findloc(dat%species_names,henry_names(j), 1)
        if (ind /= 0) then
          i = ind
          dat%henry_data(:,i) = henry_data(:,j)
        endif
      enddo
      ! set particle solubility to super high number
      dat%henry_data(1,1:dat%npq) = 7.e11_dp
    endif

  end subroutine  
  
  subroutine check_for_reaction_duplicates(dat, duplicate, err)
    use futils, only: sort
    type(PhotochemData), intent(in) :: dat
    logical, intent(in) :: duplicate(:)
    character(:), allocatable, intent(out) :: err
    character(len=:), allocatable :: rxstring
    
    integer i, ii, j, jj, nr, np, rxt, rxt_ii
    logical l, m
    
    integer, allocatable :: tmp_arr1(:), tmp_arr2(:)
    
    
    do i = 1,dat%nrT-1
      do ii = i+1,dat%nrT
        
        ! if not designated as duplicates then check if they are
        ! duplicates
        if (.not.(duplicate(i) .and. duplicate(ii))) then
        
        ! check the same num of reactants and products
        nr = dat%rx(i)%nreact
        np = dat%rx(i)%nprod
        if (i > dat%nrF) then
          j = dat%rx(i)%reverse_info
        else
          j = i
        endif
        if (ii > dat%nrF) then
          jj = dat%rx(ii)%reverse_info
        else
          jj = ii
        endif
        rxt = dat%rx(j)%rp%rxtype
        rxt_ii = dat%rx(jj)%rp%rxtype
        if (nr == dat%rx(ii)%nreact .and. np == dat%rx(ii)%nprod &
            .and. rxt == rxt_ii) then
          if (allocated(tmp_arr1)) then
            deallocate(tmp_arr1, tmp_arr2)
          endif
          allocate(tmp_arr1(nr), tmp_arr2(nr))
          
          tmp_arr1 = dat%rx(ii)%react_sp_inds
          tmp_arr2 = dat%rx(i)%react_sp_inds
          
          call sort(tmp_arr1)
          call sort(tmp_arr2)
          
          m = all(tmp_arr1 == tmp_arr2)
          if (m) then
            if (allocated(tmp_arr1)) then
              deallocate(tmp_arr1, tmp_arr2)
            endif
            allocate(tmp_arr1(np), tmp_arr2(np))
          
            tmp_arr1 = dat%rx(ii)%prod_sp_inds
            tmp_arr2 = dat%rx(i)%prod_sp_inds
            
            call sort(tmp_arr1)
            call sort(tmp_arr2)
          
            l = all(tmp_arr1 == tmp_arr2)
          
            if (l) then
              err = "IOError: This reaction is a duplicate: "
              call reaction_string(dat, i, rxstring)
              err = err//rxstring
              return
            endif
          endif
        endif
        
        endif
      enddo
    enddo
    
  end subroutine
  
  subroutine reaction_string(dat,rxn,rxstring)
    use photochem_enum, only: ThreeBodyRateType, FalloffRateType
    type(PhotochemData), intent(in) :: dat
    integer, intent(in) :: rxn
    character(len=:), allocatable, intent(out) :: rxstring
    integer j, k, i
    rxstring = ''
    if (rxn > dat%nrF) then
      i = dat%rx(rxn)%reverse_info
    else
      i = rxn
    endif
    do j = 1,dat%rx(rxn)%nreact-1
      k = dat%rx(rxn)%react_sp_inds(j)
      rxstring = rxstring //(trim(dat%species_names(k))//' + ')
    enddo
    
    k = dat%rx(rxn)%react_sp_inds(dat%rx(rxn)%nreact)
    rxstring = rxstring // trim(dat%species_names(k))//' => '
    
    if (dat%rx(i)%rp%rxtype == ThreeBodyRateType .or. dat%rx(i)%rp%rxtype == FalloffRateType) then
      rxstring = rxstring(1:len(rxstring)-4) //(' + M'//' => ')
    endif
    
    do j = 1,dat%rx(rxn)%nprod-1
      k = dat%rx(rxn)%prod_sp_inds(j)
      rxstring = rxstring // trim(dat%species_names(k))//' + '
    enddo
    k = dat%rx(rxn)%prod_sp_inds(dat%rx(rxn)%nprod)
    rxstring = rxstring // trim(dat%species_names(k))
    
    if (dat%rx(i)%rp%rxtype == ThreeBodyRateType .or. dat%rx(i)%rp%rxtype == FalloffRateType) then
      rxstring = rxstring //' + M'
    endif
  end subroutine
  
  subroutine compare_rxtype_string(tmp, eqr, eqp, reverse, rxtype_int, err)
    use photochem_enum, only: PhotolysisRateType, ElementaryRateType, ThreeBodyRateType, FalloffRateType, PressDependentRateType
    character(len=*), intent(in) :: tmp
    character(len=s_str_len), allocatable, intent(in) :: eqr(:), eqp(:)
    logical, intent(in) :: reverse
    integer, intent(in) :: rxtype_int
    character(:), allocatable, intent(out) :: err
    character(:), allocatable :: rxtype
    integer i
    logical k, j, m, l, kk, jj
    l = .false.
    m = .false.
    k = .false.
    kk = .false.
    j = .false.
    jj = .false.
    if (rxtype_int == PhotolysisRateType) then
      rxtype = 'photolysis'
    elseif (rxtype_int == ElementaryRateType) then
      rxtype = 'elementary'
    elseif (rxtype_int == ThreeBodyRateType) then
      rxtype = 'three-body'
    elseif (rxtype_int == FalloffRateType) then
      rxtype = 'falloff'
    elseif (rxtype_int == PressDependentRateType) then
      rxtype = 'pressure-dependent-Arrhenius'
    else
      err = 'Internal error in Photochem involving "compare_rxtype_string". Report this bug!'
      return
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
    elseif (trim(rxtype) == 'elementary' .or. trim(rxtype) == 'pressure-dependent-Arrhenius') then
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
    use photochem_enum, only: ShomatePolynomial, Nasa9Polynomial, Nasa7Polynomial
    use photochem_types, only: ThermodynamicData
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
      err = "IOError: Thermodynamic data must be in Shomate format for "//trim(molecule_name)
      return
    endif
    
    ! get temperature ranges
    tmplist =>tmpdict%get_list("temperature-ranges",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    thermo%ntemps = tmplist%size() - 1
    if (thermo%ntemps < 1) then
      err = "IOError: Problem reading thermodynamic data for  "//trim(molecule_name)
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
    if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
    if (tmplist%size() /= thermo%ntemps) then
      err = "IOError: Problem reading thermodynamic data for "//trim(molecule_name)
      return
    endif
    
    if (thermo%dtype == ShomatePolynomial) then
      allocate(thermo%data(7,thermo%ntemps))
    elseif (thermo%dtype == Nasa9Polynomial) then
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
          err = "IOError: Too much or too little thermodynamic data for "//trim(molecule_name)
          return
        endif
        
        j = 1
        item1 => listitem%first
        do while (associated(item1)) 
          select type (listitem1 => item1%node)
          class is (type_scalar)

            thermo%data(j, k) = listitem1%to_real(-1.0_dp,success)
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
  
  subroutine get_reaction_sp_nums(dat, rx_str, rx, reverse, err)
    use photochem_types, only: Reaction
    use photochem_enum, only: ThreeBodyRateType, FalloffRateType
    type(PhotochemData), intent(in) :: dat
    character(len=*), intent(in) :: rx_str
    type(Reaction), intent(inout) :: rx ! already has rate parameters
    logical, intent(out) :: reverse
    character(:), allocatable, intent(out) :: err
    
    character(len=s_str_len), allocatable :: eqr(:), eqp(:), eqr1(:), eqp1(:)
    integer :: reactant_atoms(dat%natoms), product_atoms(dat%natoms)
    integer :: ind(1), i
    
    
    call parse_reaction(rx_str, reverse, eqr1, eqp1, err)
    if (allocated(err)) return
    call compare_rxtype_string(rx_str, eqr1, eqp1, reverse, rx%rp%rxtype, err)
    if (allocated(err)) return
    
    if (rx%rp%rxtype == ThreeBodyRateType .or. rx%rp%rxtype == FalloffRateType) then
      ! remove the M
      eqr = eqr1(1:size(eqr1)-1)
      eqp = eqp1(1:size(eqp1)-1)
    else
      eqr = eqr1
      eqp = eqp1
    endif
    
    rx%nreact = size(eqr)
    rx%nprod = size(eqp)
    allocate(rx%react_sp_inds(rx%nreact))
    allocate(rx%prod_sp_inds(rx%nprod))
    
    do i = 1,rx%nreact
      ind = findloc(dat%species_names,eqr(i))
      rx%react_sp_inds(i) = ind(1)
      if (ind(1) == 0) then
        err = "IOError: "// & 
               "Species "//trim(eqr(i))//" in reaction "//trim(rx_str)// &
               " is not in the list of species."
        return
      endif
    enddo
    
    do i = 1,rx%nprod
      ind = findloc(dat%species_names,eqp(i))
      rx%prod_sp_inds(i) = ind(1)
      if (ind(1) == 0) then
        err = "IOError: "// & 
               "Species "//trim(eqp(i))//" in reaction "//trim(rx_str)// &
               " is not in the list of species."
        return
      endif
    enddo
    
    reactant_atoms = 0
    product_atoms = 0
    do i=1,rx%nreact
      reactant_atoms = reactant_atoms + dat%species_composition(:,rx%react_sp_inds(i))
    enddo
    do i=1,rx%nprod
      product_atoms = product_atoms + dat%species_composition(:,rx%prod_sp_inds(i))
    enddo
    if (.not. all(reactant_atoms == product_atoms)) then
      err = "IOError: "//& 
             'Bad mass balance in reaction "'//trim(rx_str)// &
             '". You could have messed up how many atoms one of the species has.'
      return
    endif
    
  end subroutine
  
  function is_rx_reverse(rx_string, err) result(reverse)
    character(len=*), intent(in) :: rx_string
    logical :: reverse
    character(:), allocatable, intent(out) :: err
    
    
    if (index(rx_string, "<=>") /= 0) then
      reverse = .true.
    elseif (index(rx_string, " =>") /= 0) then
      reverse = .false.
    else
      err = "IOError: Invalid reaction arrow in reaction "//trim(rx_string)// &
            '. Note, forward reactions must have a space before the arrow, like " =>"'
      return
    endif
    
  end function
  
  subroutine get_rateparams(dat, reaction_d, infile, rx, err)
    use photochem_enum, only: NoFalloff, TroeWithoutT2Falloff, TroeWithT2Falloff, JPLFalloff
    use photochem_enum, only: PhotolysisRateType, ElementaryRateType, ThreeBodyRateType, FalloffRateType, PressDependentRateType
    use photochem_types, only: Reaction, BaseRate, ElementaryRate, ThreeBodyRate, FalloffRate, PhotolysisRate, PressDependentRate
    type(PhotochemData), intent(in) :: dat
    class(type_dictionary), intent(in) :: reaction_d
    character(len=*), intent(in) :: infile
    type(Reaction), target, intent(out) :: rx
    character(:), allocatable, intent(out) :: err
    
    class(BaseRate), pointer :: rp
    type (type_error), allocatable :: io_err
    character(len=str_len) :: rxtype_str
    type(type_dictionary), pointer :: dict
    type(type_list), pointer :: list
    class(type_list_item), pointer :: item
    logical :: use_jpl
    real(dp) :: T2
        
    rxtype_str = reaction_d%get_string("type", default="elementary", error = io_err) 
    
    if (rxtype_str == 'photolysis') then
      allocate(PhotolysisRate::rx%rp)
      rx%rp%rxtype = PhotolysisRateType
    elseif (rxtype_str == 'elementary') then
      allocate(ElementaryRate::rx%rp)
      rx%rp%rxtype = ElementaryRateType
    elseif (rxtype_str == 'three-body') then
      allocate(ThreeBodyRate::rx%rp)
      rx%rp%rxtype = ThreeBodyRateType
    elseif (rxtype_str == 'falloff') then
      allocate(FalloffRate::rx%rp)
      rx%rp%rxtype = FalloffRateType
    elseif (rxtype_str == 'pressure-dependent-Arrhenius') then
      allocate(PressDependentRate::rx%rp)
      rx%rp%rxtype = PressDependentRateType
    else
      err = 'IOError: reaction type '//trim(rxtype_str)//' is not a valid reaction type.'
      return
    endif
    
    rp => rx%rp
    
    select type (rp)
    class is (PhotolysisRate)
      ! No rate info
    class is (ElementaryRate)
      dict => reaction_d%get_dictionary('rate-constant',.true.,error = io_err)
      if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      call get_arrhenius(dict, rp%A, rp%b, rp%Ea, err)
      if (allocated(err)) return
        
    class is (ThreeBodyRate)
      dict => reaction_d%get_dictionary('rate-constant',.true.,error = io_err)
      if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      call get_arrhenius(dict, rp%A, rp%b, rp%Ea, err)
      if (allocated(err)) return
      
      call get_efficiencies(dat, reaction_d, rp%eff, err)
      if (allocated(err)) return
      
    class is (FalloffRate)
      
      dict => reaction_d%get_dictionary('low-P-rate-constant',.true.,error = io_err)
      if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      call get_arrhenius(dict, rp%A0, rp%b0, rp%Ea0, err)
      if (allocated(err)) return
        
      dict => reaction_d%get_dictionary('high-P-rate-constant',.true.,error = io_err)
      if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      
      call get_arrhenius(dict, rp%Ainf, rp%binf, rp%Eainf, err)
      if (allocated(err)) return
      
      call get_efficiencies(dat, reaction_d, rp%eff, err)
      if (allocated(err)) return
      
      ! get falloff stuff
      use_jpl = reaction_d%get_logical('JPL',default=.false.,error = io_err)
      nullify(dict)
      dict => reaction_d%get_dictionary('Troe',required=.false.,error = io_err)
      if (associated(dict) .and. use_jpl) then
        err = "Both 'Troe' and 'JPL' falloff types are specified for reaction "// &
              trim(reaction_d%get_string("equation",error = io_err))//". Only one is allowed"
        return
      endif
      
      rp%falloff_type = NoFalloff
      
      if (use_jpl) then
        rp%falloff_type = JPLFalloff
      endif
      
      if (associated(dict)) then
        allocate(rp%A_T)
        allocate(rp%T1)
        allocate(rp%T3)
        
        rp%A_T = dict%get_real('A',error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        rp%T1 = dict%get_real('T1',error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        
        T2 = dict%get_real('T2',error = io_err)
        if (allocated(io_err)) then ! T2 is not there
          rp%falloff_type = TroeWithoutT2Falloff
          deallocate(io_err)
        else ! T2 is there
          rp%falloff_type = TroeWithT2Falloff
          allocate(rp%T2)
          rp%T2 = T2
        endif
        
        rp%T3 = dict%get_real('T3',error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
      endif

    class is (PressDependentRate); block
      use futils, only: argsort
      use photochem_types, only: MultiArrheniusRate
      type(MultiArrheniusRate) :: rate
      real(dp), allocatable :: P(:), A(:), b(:), Ea(:)
      integer, allocatable :: inds(:)
      integer :: i, j

      list => reaction_d%get_list('rate-constants',.true.,error=io_err)
      if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif

      i = list%size()
      allocate(P(i), A(i), b(i), Ea(i), inds(i))

      i = 1
      item => list%first
      do while (associated(item))
        select type (listitem => item%node)
        class is (type_dictionary)

          P(i) = listitem%get_real('P',error = io_err)
          if (allocated(io_err)) then; err = trim(io_err%message); return; endif

          if (P(i) <= 0.0_dp) then
            err = 'Pressures must be positive in pressure-dependent-Arrhenius reaction '// &
            trim(reaction_d%get_string("equation",error=io_err))
            return
          endif
          
          call get_arrhenius(listitem, A(i), b(i), Ea(i), err)
          if (allocated(err)) return

        class default
          err = '"rate-constants" must be a list of dictionaries for reaction '// &
          trim(reaction_d%get_string("equation",error=io_err))
          return
        end select
        i = i + 1
        item => item%next
      enddo

      ! sort
      inds = argsort(P)
      P = P(inds)
      A = A(inds)
      b = b(inds)
      Ea = Ea(inds)

      allocate(rp%logP(0))
      allocate(rp%rate(0))
      i = 1
      do while (i <= size(P))
        rp%logP = [rp%logP, log(P(i))]
        rate = MultiArrheniusRate([A(i)], [b(i)], [Ea(i)])  
        j = i + 1
        do
          if (j > size(P)) then
            i = j
            exit
          endif
          if (P(j) == P(i)) then
            rate%A = [rate%A, A(j)]
            rate%b = [rate%b, b(j)]
            rate%Ea = [rate%Ea, Ea(j)]
          else
            i = j
            exit
          endif
          j = j + 1
        enddo
        rp%rate = [rp%rate, rate]
      enddo

      if (size(rp%logP) < 2) then
        err = 'pressure-dependent-Arrhenius reaction '// &
        trim(reaction_d%get_string("equation",error=io_err))// &
        ' must have > 1 unique pressures.'
        return
      endif
      
    endblock; end select
    
  end subroutine
  
  subroutine get_efficiencies(dat, rx, eff, err)
    use photochem_types, only: Efficiencies
    type(PhotochemData), intent(in) :: dat
    type(type_dictionary), intent(in) :: rx
    type(Efficiencies), intent(out) :: eff
    character(:), allocatable, intent(out) :: err
    
    type (type_error), allocatable :: io_err
    type (type_key_value_pair), pointer :: key_value_pair
    type(type_dictionary), pointer :: d
    integer :: ind(1), j
    
    
    eff%def_eff = rx%get_real("default-efficiency",default=1.0_dp,error = io_err)
  
    eff%n_eff = 0
    d => rx%get_dictionary("efficiencies",required=.false.,error = io_err)
    
    if (associated(d)) then
      eff%n_eff = d%size()
      allocate(eff%efficiencies(eff%n_eff))
      allocate(eff%eff_sp_inds(eff%n_eff))
      
      j = 1
      key_value_pair => d%first
      do while (associated(key_value_pair))
        ind = findloc(dat%species_names,trim(key_value_pair%key))
        if (ind(1) == 0) then
          err = 'IOError: Reaction '//trim(rx%get_string("equation",error = io_err))//&
                ' has efficiencies for species that are'// &
                ' not in the list of species'
          return 
        endif
        
        eff%eff_sp_inds(j) = ind(1)
        eff%efficiencies(j) = d%get_real(trim(key_value_pair%key),error = io_err)
        if (allocated(io_err)) then; err = trim(io_err%message); return; endif
        
        key_value_pair => key_value_pair%next
        j = j + 1
      enddo
      
    endif
    
  end subroutine
  
  subroutine get_arrhenius(d, A, b, Ea, err)
    type(type_dictionary), intent(in) :: d
    real(dp), intent(out) :: A, b, Ea
    character(:), allocatable, intent(out) :: err
    
    type (type_error), allocatable :: io_err
    
    
    A = d%get_real('A',error = io_err)
    if (allocated(io_err)) then; err = trim(io_err%message); return; endif
    
    b = d%get_real('b',error = io_err)
    if (allocated(io_err)) then; err = trim(io_err%message); return; endif
    
    Ea = d%get_real('Ea',error = io_err)
    if (allocated(io_err)) then; err = trim(io_err%message); return; endif
    
  end subroutine

  subroutine parse_reaction(instring, reverse, eqr, eqp, err)
    character(len=*), intent(in) :: instring
    logical, intent(out) :: reverse
    character(len=s_str_len), allocatable, intent(out) :: eqr(:), eqp(:)
    character(:), allocatable, intent(out) :: err
    
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
  
  subroutine check_sl(dat, err)
    use photochem_types, only: ThreeBodyRate, FalloffRate
    type(PhotochemData), intent(in) :: dat
    character(:), allocatable, intent(out) :: err
    
    integer :: i, j, l, k, kk, mm, n, nn, ind(1), counter
    character(len=:), allocatable :: reaction
    
    do i = 1, dat%nsl
      
      j = dat%nq + i
      ! can not be an efficiency.
      do k = 1,dat%nrF
        select type(rp => dat%rx(k)%rp)        
        class is (ThreeBodyRate)
          if (rp%eff%n_eff > 0) then
            ind = findloc(rp%eff%eff_sp_inds,j)
            if (ind(1) /= 0) then
              call reaction_string(dat,k,reaction)
              err = 'IOError: Reaction "'//reaction//'" has short-lived species collision efficiencies.' // &
              ' This is not allowed. Either remove the efficiencies, or change the species to long lived.'
              return
            endif
          endif
        class is (FalloffRate)
          if (rp%eff%n_eff > 0) then
            ind = findloc(rp%eff%eff_sp_inds, j)
            if (ind(1) /= 0) then
              call reaction_string(dat,k,reaction)
              err = 'IOError: Reaction "'//reaction//'" has short-lived species collision efficiencies.' // &
              ' This is not allowed. Either remove the efficiencies, or change the species to long lived.'
              return
            endif
          endif
        end select
      enddo
      
      l = dat%pl(j)%nump
      do k = 1,l
        kk = dat%pl(j)%iprod(k)
        do mm = 1, dat%rx(kk)%nreact
          ! are SL species produced by other SL species?
          do n = dat%nq+1,dat%nq + dat%nsl
            if (n == dat%rx(kk)%react_sp_inds(mm)) then
              call reaction_string(dat, kk, reaction)
              err = 'IOError: Reaction "'//reaction//'" has short-lived species as reactants'// &
              ' and products. This is not allowed. Change one or both of the species to long-lived.'
              return
            endif
          enddo
        enddo
      enddo

      l = dat%pl(j)%numl
      do k = 1,l
        kk = dat%pl(j)%iloss(k)     
        counter = 0
        do mm = 1, dat%rx(kk)%nreact
          n = dat%rx(kk)%react_sp_inds(mm)
          do nn = dat%nq+1,dat%nq + dat%nsl
            if (nn == n .and. n == j) then
              counter = counter + 1
              if (counter > 1) then
                call reaction_string(dat, kk, reaction)   
                err = 'IOError: Reaction "'//reaction//'" short lived species react'// &
                ' with themselves. This is not allowed. Change the species to long lived.'
                return
              endif
            elseif (nn == n .and. n /= j) then
              call reaction_string(dat, kk, reaction)   
              err = 'IOError: Reaction "'//reaction//'" short lived species react'// &
              ' with other short lived species. This is not allowed.'
              return
            elseif (dat%species_names(n) == 'hv') then
              call reaction_string(dat, kk, reaction)   
              err = 'IOError: Photolysis reaction "'//reaction//'" can not have short lived species.'
              return
            endif
          enddo
        enddo
      enddo
      
    enddo
    
  end subroutine
  
  subroutine get_photorad(dat, var, err)
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    character(:), allocatable, intent(out) :: err
    
    ! Read wavelength grid
    call get_wavl(dat, var, err)
    if (allocated(err)) return
    
    ! get rayleigh
    call get_rayleigh(dat, var, err)
    if (allocated(err)) return
    
    ! get photolysis xsections data
    call get_photolysis_xs(dat, var, err)
    if (allocated(err)) return
    
    if (dat%there_are_particles) then
      call get_aerosol_xs(dat, var, err)
      if (allocated(err)) return
    endif
    
  end subroutine

  subroutine get_wavl(dat, var, err)
    use h5fortran
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    character(:), allocatable, intent(out) :: err

    type(hdf5_file) :: h
    integer(HSIZE_T), allocatable :: dims(:)
    character(:), allocatable :: xsroot, filename

    ! Root directory
    xsroot = trim(var%data_dir)//"/xsections/"
    filename = xsroot//'bins.h5'

    if (.not.is_hdf5(filename)) then
      err = 'Failed to read in photolysis grid from "'//filename//'"'
      return
    endif

    call h%open(filename,'r')

    call check_h5_dataset(h, 'wavl', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close()
      return
    endif
    call h%shape("wavl", dims)
    allocate(dat%wavl(dims(1)))
    call h%read("wavl", dat%wavl)
    dat%nw = size(dat%wavl) - 1

    call h%close()

  end subroutine

  subroutine read_mie_data_file(filename, nw, wavl, &
                                 nrad_file, radii_file, w0_file, qext_file, g_file, err)
    use h5fortran
    use futils, only: interp_discrete_to_bins
    use clima_useful, only: hdf5_file_closer
    character(*), intent(in) :: filename
    integer, intent(in) :: nw
    real(dp), intent(in) :: wavl(:)
    
    integer, intent(out) :: nrad_file
    real(dp), allocatable, intent(out) :: radii_file(:)
    real(dp), allocatable, intent(out) :: w0_file(:,:), qext_file(:,:), g_file(:,:)
    character(:), allocatable, intent(out) :: err

    type(hdf5_file), target :: h
    integer(HSIZE_T), allocatable :: dims(:)
    real(dp), allocatable :: wv_tmp(:)
    real(dp), allocatable :: w0_tmp(:,:), qext_tmp(:,:), g0_tmp(:,:)
    integer :: i
    
    if (.not.is_hdf5(filename)) then
      err = "Was unable to open mie data file "//trim(filename)
      return
    endif
    
    block
    type(hdf5_file_closer) :: h_closer

    ! Open file
    call h%open(filename,'r')
    h_closer%h => h

    ! Wavelengths
    call check_h5_dataset(h, 'wavelengths', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("wavelengths", dims)
    allocate(wv_tmp(dims(1)))
    call h%read("wavelengths", wv_tmp)

    ! Radii
    call check_h5_dataset(h, 'radii', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("radii", dims)
    nrad_file = dims(1)
    allocate(radii_file(dims(1)))
    call h%read("radii", radii_file)

    ! w0
    call check_h5_dataset(h, 'w0', 2, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("w0", dims)
    allocate(w0_tmp(dims(1),dims(2)))
    call h%read("w0", w0_tmp)

    ! qext
    call check_h5_dataset(h, 'qext', 2, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("qext", dims)
    allocate(qext_tmp(dims(1),dims(2)))
    call h%read("qext", qext_tmp)

    ! g0
    call check_h5_dataset(h, 'g0', 2, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("g0", dims)
    allocate(g0_tmp(dims(1),dims(2)))
    call h%read("g0", g0_tmp)

    if (size(w0_tmp,1) /= nrad_file .or. size(w0_tmp,2) /= size(wv_tmp)) then
      err = '"w0" has the wrong shape in "'//filename//'"'
      return
    endif
    if (size(qext_tmp,1) /= nrad_file.or. size(qext_tmp,2) /= size(wv_tmp)) then
      err = '"qext" has the wrong shape in "'//filename//'"'
      return
    endif
    if (size(g0_tmp,1) /= nrad_file .or. size(g0_tmp,2) /= size(wv_tmp)) then
      err = '"g0" has the wrong shape in "'//filename//'"'
      return
    endif

    end block

    ! Now lets interpolate a bunch
    allocate(w0_file(nrad_file,size(wavl)-1))
    allocate(qext_file(nrad_file,size(wavl)-1))
    allocate(g_file(nrad_file,size(wavl)-1))

    do i = 1, nrad_file

      ! w0
      call interp_discrete_to_bins(wavl, wv_tmp, w0_tmp(i,:), w0_file(i,:), 'Constant', err=err)
      if (allocated(err)) return

      ! ext
      call interp_discrete_to_bins(wavl, wv_tmp, qext_tmp(i,:), qext_file(i,:), 'Constant', err=err)
      if (allocated(err)) return

      ! g0
      call interp_discrete_to_bins(wavl, wv_tmp, g0_tmp(i,:), g_file(i,:), 'Constant', err=err)
      if (allocated(err)) return

    enddo

  end subroutine
  
  subroutine get_aerosol_xs(dat, var, err)
    use photochem_enum, only: MieParticle, FractalParticle
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    character(:), allocatable, intent(out) :: err
    
    integer :: nrad
    integer, parameter :: nrad_fixed = 50
    real(dp), allocatable :: radii(:)
    real(dp), allocatable :: w0(:,:), qext(:,:), g(:,:)
    character(len=:), allocatable :: xsroot
    character(len=:), allocatable :: filename
    integer :: i
    
    xsroot = trim(var%data_dir)//"/aerosol_xsections/"
    
    allocate(dat%radii_file(nrad_fixed,dat%np))
    allocate(dat%part_xs_file(dat%np))
    
    dat%nrad_file = nrad_fixed
    
    do i = 1,dat%np
      
      if (dat%particle_optical_prop(i) == 'none') then
        ! there is no optical data, so we skip
        dat%part_xs_file(i)%ThereIsData = .false.
        cycle
      else
        ! there is optical data, so we allocate and get it
        dat%part_xs_file(i)%ThereIsData = .true.
        allocate(dat%part_xs_file(i)%w0(nrad_fixed,dat%nw))
        allocate(dat%part_xs_file(i)%qext(nrad_fixed,dat%nw))
        allocate(dat%part_xs_file(i)%gt(nrad_fixed,dat%nw))
      endif
      
      if (dat%particle_optical_type(i) == MieParticle) then
        filename = xsroot//trim(dat%particle_optical_prop(i))// &
                  "/mie_"//trim(dat%particle_optical_prop(i))//".h5"
      elseif (dat%particle_optical_type(i) == FractalParticle) then
        filename = xsroot//trim(dat%particle_optical_prop(i))// &
                  "/frac_"//trim(dat%particle_optical_prop(i))//".h5"
      endif
      
      if (allocated(radii)) then
        deallocate(radii, w0, qext, g)
      endif
      call read_mie_data_file(filename, dat%nw, dat%wavl, &
                              nrad, radii, w0, qext, g, err) 
      if (allocated(err)) return
      if (nrad /= nrad_fixed) then
        err = "IOError: Aerosol data file "//filename// &
              "must have 50 radii bins."
        return
      endif
      
      dat%radii_file(:,i) = radii/1.e4_dp ! convert from micron to cm
      dat%part_xs_file(i)%w0 = w0
      dat%part_xs_file(i)%qext = qext
      dat%part_xs_file(i)%gt = g

    enddo

  end subroutine
  
  subroutine get_photolysis_xs(dat, var, err)
    use h5fortran, only: is_hdf5
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    character(:), allocatable, intent(out) :: err
    
    character(:), allocatable :: xsroot, species, filename, reaction
    integer :: i, j, k

    ! Root directory
    xsroot = trim(var%data_dir)//"/xsections/"

    !~~ First, we get total absorption ~~!

    ! Count Species with opacities
    k = 0
    do i = 1,size(dat%species_names(dat%ng_1:dat%nq))
      j = i + dat%np
      species = trim(dat%species_names(j))
      filename = xsroot//species//'.h5'

      ! Skip the gases without opacities
      if (.not. is_hdf5(filename)) cycle

      k = k + 1
    enddo

    ! Allocate
    allocate(dat%absorp_xs(k))

    ! Get the opacities
    k = 1
    do i = 1,size(dat%species_names(dat%ng_1:dat%nq))
      j = i + dat%np
      species = trim(dat%species_names(j))
      filename = xsroot//species//'.h5'

      ! Skip the gases without opacities
      if (.not. is_hdf5(filename)) cycle

      dat%absorp_xs(k) = create_XsectionData(filename, j, dat%wavl, err)
      if (allocated(err)) return

      k = k + 1
    enddo

    !~~ Now, get the cross sections for specific reactions ~~!

    allocate(dat%photolysis_xs(dat%kj,dat%nw))
    do i = 1,dat%kj

      ! Get the reaction name
      j = dat%photonums(i)
      call reaction_string(dat, j, reaction)

      ! Get the absorbing species
      k = dat%rx(j)%react_sp_inds(1)
      species = trim(dat%species_names(k))

      ! Filename is based on the species
      filename = xsroot//species//'.h5'

      call get_photolysis_xs_for_reaction(filename, reaction, dat%wavl, dat%photolysis_xs(i,:), err)
      if (allocated(err)) return

    enddo
    
  end subroutine

  subroutine get_photolysis_xs_for_reaction(filename, reaction, wavl, xs, err)
    use futils, only: interp
    use h5fortran
    use futils, only: interp_discrete_to_bins
    use clima_useful, only: hdf5_file_closer
    character(*), intent(in) :: filename
    character(*), intent(in) :: reaction
    real(dp), intent(in) :: wavl(:)
    real(dp), intent(out) :: xs(:)
    character(:), allocatable, intent(out) :: err
    
    type(hdf5_file), target :: h
    integer(HSIZE_T), allocatable :: dims(:)
    real(dp), allocatable :: wv_tmp(:), xs_tmp(:), wv1(:), qy1(:), qy2(:), wv_av(:)
    integer :: ierr
    character(:), allocatable :: err_base

    err_base = 'Failed to retrieve photolysis cross section data for the reaction "'//reaction//'". '

    block
    type(hdf5_file_closer) :: h_closer

    if (.not. is_hdf5(filename)) then
      err = err_base//'The file "'//filename//'" does not exist.'
      return
    endif

    ! Open file
    call h%open(filename,'r')
    h_closer%h => h

    ! Wavelengths
    call check_h5_dataset(h, 'wavelengths', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("wavelengths", dims)
    allocate(wv_tmp(dims(1)))
    call h%read("wavelengths", wv_tmp)

    ! Photodissociation
    call check_h5_dataset(h, 'photodissociation', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("photodissociation", dims)
    allocate(xs_tmp(dims(1)))
    call h%read("photodissociation", xs_tmp)

    ! Interpolate
    call interp_discrete_to_bins(wavl, wv_tmp, xs_tmp, xs, 'FillValue', tiny(1.0_dp), err)
    if (allocated(err)) return

    ! Now do the quantum yields
    if (.not. h%exists("photodissociation-qy")) then
      err = err_base//filename//': dataset "photodissociation-qy" does not exist'
      return
    endif

    call h%open_group("photodissociation-qy")
    
    ! Wavelengths
    call check_h5_dataset(h, 'wavelengths', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close_group()
      return
    endif
    call h%shape("wavelengths", dims)
    allocate(wv1(dims(1)))
    call h%read("wavelengths", wv1)

    ! Quantum yields
    call check_h5_dataset(h, reaction, 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close_group()
      return
    endif
    call h%shape(reaction, dims)
    allocate(qy1(dims(1)))
    call h%read(reaction, qy1)

    ! Interpolate yield to the grid
    wv_av = (wavl(2:) + wavl(1:size(wavl)-1))/2.0_dp
    allocate(qy2(size(wv_av)))
    call interp(wv_av, wv1, qy1, qy2, ierr=ierr)
    if (ierr /= 0) then
      call h%close_group()
      err = err_base//'Subroutine "interp" returned an error.'
      return
    endif

    call h%close_group()

    ! Multiply the photolysis cross section by the yield.
    xs = xs*qy2
    
    end block

  end subroutine

  function create_XsectionData(filename, sp_ind, wavl, err) result(xs)
    use h5fortran
    use clima_useful, only: hdf5_file_closer
    use futils, only: interp_discrete_to_bins
    use photochem_types, only: XsectionData
    character(*), intent(in) :: filename
    integer, intent(in) :: sp_ind
    real(dp), intent(in) :: wavl(:)
    character(:), allocatable, intent(out) :: err
    type(XsectionData) :: xs

    type(hdf5_file), target :: h
    integer(HSIZE_T), allocatable :: dims(:)
    real(dp), allocatable :: wv_tmp(:), xs_tmp(:)

    ! Set the species index
    xs%sp_ind = sp_ind

    ! Block to close hdf5 file when we leave this function
    block
    type(hdf5_file_closer) :: h_closer

    ! Open file
    call h%open(filename,'r')
    h_closer%h => h

    ! Wavelengths
    call check_h5_dataset(h, 'wavelengths', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("wavelengths", dims)
    allocate(wv_tmp(dims(1)))
    call h%read("wavelengths", wv_tmp)

    ! Photoabsorption
    call check_h5_dataset(h, 'photoabsorption', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("photoabsorption", dims)
    allocate(xs_tmp(dims(1)))
    call h%read("photoabsorption", xs_tmp)

    ! Interpolate
    allocate(xs%xs(size(wavl)-1))
    call interp_discrete_to_bins(wavl, wv_tmp, xs_tmp, xs%xs, 'FillValue', tiny(1.0_dp), err)
    if (allocated(err)) return

    end block

  end function

  subroutine check_h5_dataset(h, dataset, ndims, dtype, prefix, err)
    use h5fortran, only: hdf5_file
    
    type(hdf5_file), intent(in) :: h
    character(*), intent(in) :: dataset
    integer, intent(in) :: ndims
    integer, intent(in) :: dtype
    character(*), intent(in) :: prefix
    character(:), allocatable, intent(out) :: err
    
    if (.not. h%exists(dataset)) then
      err = prefix//': dataset "'//dataset//'" does not exist'
      return
    endif
    
    if (h%ndims(dataset) /= ndims) then
      err = prefix//': dataset "'//dataset//'" has wrong number of dimensions'
      return
    endif
    
    if (h%class(dataset) /= dtype) then
      err = prefix//': dataset "'//dataset//'" has the wrong type'
      return
    endif
    
  end subroutine
  
  subroutine get_rayleigh(dat, var, err)
    use fortran_yaml_c, only : YamlFile
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(in) :: var
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: A(:), B(:), Delta(:)
    character(len=str_len) :: rayleigh_file

    type(YamlFile) :: file
    integer :: i, j
    
    rayleigh_file = trim(var%data_dir)//"/rayleigh/rayleigh.yaml"
    
    ! parse yaml file
    call file%parse(rayleigh_file, err)
    if (allocated(err)) return
    select type (root => file%root)
      class is (type_dictionary)
        call rayleigh_params(root,dat,trim(rayleigh_file),err, &
                             dat%raynums, A, B, Delta)
    end select
    call file%finalize()
    if (allocated(err)) return
    
    ! compute cross sections
    dat%nray = size(A)
    allocate(dat%sigray(dat%nray,dat%nw))
    do i = 1,dat%nw
      do j = 1,size(A)
        call rayleigh_vardavas(A(j), B(j), Delta(j), dat%wavl(i), &
                               dat%sigray(j, i))
      enddo
    enddo
    deallocate(A,B,Delta)
  end subroutine
  
  subroutine rayleigh_params(mapping,dat,infile,err, raynums, A, B, Delta)
    class (type_dictionary), intent(in), pointer :: mapping
    type(PhotochemData), intent(in) :: dat
    character(len=*), intent(in) :: infile
    character(:), allocatable, intent(out) :: err
    
    class (type_key_value_pair), pointer :: key_value_pair
    class (type_dictionary), pointer :: tmp1, tmp2
    type (type_error), allocatable :: io_err
    real(dp), allocatable, intent(out) :: A(:), B(:), Delta(:)
    integer, allocatable, intent(out) :: raynums(:)
    
    integer :: j, ind(1)  
    
    j = 0
    key_value_pair => mapping%first
    do while (associated(key_value_pair))
      ind = findloc(dat%species_names,trim(key_value_pair%key))
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
      ind = findloc(dat%species_names,trim(key_value_pair%key))
      if (ind(1) /= 0) then 
        raynums(j) = ind(1)
        
        tmp1 => mapping%get_dictionary(key_value_pair%key,.true.,error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        if (trim(tmp1%get_string('formalism',error=io_err)) /= 'vardavas') then
          err = "Unknown formalism for Rayleigh cross section for "//trim(key_value_pair%key)
          return
        endif
        tmp2 => tmp1%get_dictionary("data",.true.,error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        Delta(j) = tmp2%get_real("Delta",error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        A(j) = tmp2%get_real("A",error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        B(j) = tmp2%get_real("B",error = io_err)
        if (allocated(io_err)) then; err = trim(infile)//trim(io_err%message); return; endif
        j = j + 1
      endif
      key_value_pair => key_value_pair%next
    enddo
  end subroutine
  
  subroutine rayleigh_vardavas(A, B, Delta, lambda, sigray)
    real(dp), intent(in) :: A, B, Delta, lambda
    real(dp), intent(out) :: sigray
    
    sigray = 4.577e-21_dp*((6.0_dp+3.0_dp*Delta)/(6.0_dp-7.0_dp*Delta)) * &
            (A*(1.0_dp+B/(lambda*1.0e-3_dp)**2.0_dp))**2.0_dp * &
            (1.0_dp/(lambda*1.0e-3_dp)**4.0_dp)

  end subroutine
  
  subroutine read_stellar_flux(star_file, nw, wavl, photon_flux, err)
    use futils, only: inter2, addpnt, FileCloser
    use photochem_const, only: c_light, plank
    
    character(len=*), intent(in) :: star_file
    integer, intent(in) :: nw
    real(dp), intent(in) :: wavl(nw+1)
    real(dp), intent(out) :: photon_flux(nw)
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: file_wav(:), file_flux(:)
    real(dp) :: flux(nw)
    real(dp) :: dum1, dum2
    integer :: io, i, n, ierr
    real(dp), parameter :: rdelta = 1.0e-4_dp
    type(FileCloser) :: file
    
    open(1,file=star_file,status='old',iostat=io)
    file%unit = 1
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
    
    i = n
    ! interpolate 
    call addpnt(file_wav, file_flux, n+4, i, file_wav(1)*(1.0_dp-rdelta), 0.0_dp, ierr)
    call addpnt(file_wav, file_flux, n+4, i, 0.0_dp, 0.0_dp, ierr)
    call addpnt(file_wav, file_flux, n+4, i, file_wav(i)*(1.0_dp+rdelta), 0.0_dp,ierr)
    call addpnt(file_wav, file_flux, n+4, i, huge(rdelta), 0.0_dp,ierr)
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
      photon_flux(i) = (1/(plank*c_light*1.e16_dp))*flux(i)*(wavl(i+1)-wavl(i))* &
                       ((wavl(i+1)+wavl(i))/2.0_dp)
    enddo
    
  end subroutine
  
  subroutine read_atmosphere_file(atmosphere_txt, dat, var, err)
    use futils, only: FileCloser
    character(len=*), intent(in) :: atmosphere_txt
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    character(len=10000) :: line
    character(len=s_str_len) :: arr1(1000)
    character(len=s_str_len) :: arr11(1000)
    character(len=s_str_len),allocatable, dimension(:) :: labels
    integer :: ind(1)
    real(dp), allocatable :: temp(:,:)
    integer :: io, i, n, nn, ii
    type(FileCloser) :: file
    
    open(4, file=trim(atmosphere_txt),status='old',iostat=io)
    file%unit = 4
    if (io /= 0) then
      err = 'Can not open file '//trim(atmosphere_txt)
      return
    endif
    read(4,'(A)') line
    
    dat%nzf = -1
    io = 0
    do while (io == 0)
      read(4,*,iostat=io)
      dat%nzf = dat%nzf + 1
    enddo
    
    allocate(dat%z_file(dat%nzf))
    allocate(dat%T_file(dat%nzf))
    allocate(dat%edd_file(dat%nzf))
    allocate(dat%den_file(dat%nzf))
    allocate(dat%mix_file(dat%nq, dat%nzf))
    dat%z_file = 0.0_dp
    dat%T_file = 0.0_dp
    dat%edd_file = 0.0_dp
    dat%mix_file = 1.0e-40_dp
    if (dat%there_are_particles) then
      allocate(dat%particle_radius_file(dat%npq, dat%nzf))
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
    allocate(temp(n,dat%nzf))
    rewind(4)
    read(4,'(A)') line
    read(line,*) (labels(i),i=1,n)
    
    ! First read in all the data into big array
    do i = 1,dat%nzf
      read(4,*,iostat=io) (temp(ii,i),ii=1,n)
      if (io /= 0) then
        err = 'Problem reading in initial atmosphere in '//trim(atmosphere_txt)
        return
      endif
    enddo
    
    ! reads in mixing ratios
    do i=1,dat%nq
      ind = findloc(labels,dat%species_names(i))
      if (ind(1) /= 0) then
        dat%mix_file(i,:) = temp(ind(1),:)
      endif
    enddo
    
    if (dat%there_are_particles) then
      do i=1,dat%npq
        ind = findloc(labels,trim(dat%species_names(i))//"_r")
        if (ind(1) /= 0) then
          dat%particle_radius_file(i,:) = temp(ind(1),:)
        else
          ! did not find the data
          ! will set to 0.1 micron
          dat%particle_radius_file(i,:) = 1.0e-5_dp
        endif
      enddo
    endif
    
    ! reads in temperature
    ind = findloc(labels,'temp')
    if (ind(1) /= 0) then
      dat%T_file(:) = temp(ind(1),:)
    else
      err = '"temp" was not found in input file '//trim(atmosphere_txt)
      return
    endif
    
    ! reads in alt
    ind = findloc(labels,'alt')
    if (ind(1) /= 0) then
      dat%z_file(:) = temp(ind(1),:)*1.e5_dp ! conver to cm
    else
      err = '"alt" was not found in input file '//trim(atmosphere_txt)
      return
    endif
    
    ! reads in eddy diffusion
    ind = findloc(labels,'eddy')
    if (ind(1) /= 0) then
      dat%edd_file(:) = temp(ind(1),:)
    else
      err = '"eddy" was not found in input file '//trim(atmosphere_txt)
      return
    endif

    ! reads in density.
    ind = findloc(labels,'den')
    if (ind(1) /= 0) then
      dat%den_file(:) = temp(ind(1),:)
    else
      err = '"den" was not found in input file '//trim(atmosphere_txt)
      return
    endif

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
  
end submodule