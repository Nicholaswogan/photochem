module photochem_data
  use fortran_yaml_c_types, only: type_node, type_dictionary, type_list, type_error, &
                                 type_list_item, type_scalar, type_key_value_pair
  use photochem_const, only: dp, str_len, s_str_len, m_str_len
  use photochem_settings, only: PhotoSettings
  use clima_saturationdata, only: SaturationData
  implicit none
  private

  public :: PhotochemData
  public :: XsectionData, ParticleXsections, ThermodynamicData
  public :: Reaction, Efficiencies, BaseRate, ReverseRate, PhotolysisRate
  public :: ElementaryRate, ThreeBodyRate, FalloffRate, PressDependentRate
  public :: MultiArrheniusRate, ProdLoss
  public :: read_photochem_mechanism, read_photochem_supporting_data
  public :: parse_reaction

  type, extends(type_list) :: type_list_tmp
    ! Temporary list for accessing all reactions and species in a row.
  contains
    final :: list_destroy
  end type

  type :: XsectionData
    integer :: sp_ind !! Species index
    real(dp), allocatable :: xs(:) !! The cross section in cm^2/molecule (nw)
  end type

  type :: ParticleXsections
    logical :: ThereIsData
    real(dp), allocatable :: w0(:,:) ! (nz,nw) or (nrad_file, nw)
    real(dp), allocatable :: qext(:,:)
    real(dp), allocatable :: gt(:,:)
  end type

  type :: ThermodynamicData
    integer :: dtype ! shomate = 1
    integer :: ntemps
    real(dp), allocatable :: temps(:)
    real(dp), allocatable :: data(:,:)
  end type

  type :: Efficiencies
    integer :: n_eff !! number of efficiencies
    real(dp) :: def_eff !! default efficiency
    real(dp), allocatable :: efficiencies(:) !! 3-body efficiencies
    integer, allocatable :: eff_sp_inds(:) !! species indices for each efficiency
  end type

  type, abstract :: BaseRate
    integer :: rxtype !! Specifies the types of rate (see enums)
  end type

  type, extends(BaseRate) :: ReverseRate
    ! no rate parameters.
  end type

  type, extends(BaseRate) :: PhotolysisRate
    ! no rate parameters. calculated via radiative transfer
  end type

  type, extends(BaseRate) :: ElementaryRate
    ! rate = A*T^b*exp(-Ea/T)
    real(dp) :: A !! pre-exponential factor (various units)
    real(dp) :: b !! temperature exponent (unitless)
    real(dp) :: Ea !! Activate energy (T)
  end type

  type, extends(BaseRate) :: ThreeBodyRate
    real(dp) :: A !! pre-exponential factor (various units)
    real(dp) :: b !! temperature exponent (unitless)
    real(dp) :: Ea !! Activate energy (T)
    type(Efficiencies) :: eff
  end type

  type, extends(BaseRate) :: FalloffRate
    real(dp) :: A0 !! For low-P rate constant
    real(dp) :: b0 !! For low-P rate constant
    real(dp) :: Ea0 !! For low-P rate constant
    real(dp) :: Ainf !! For high-P rate constant
    real(dp) :: binf !! For high-P rate constant
    real(dp) :: Eainf !! For high-P rate constant

    integer :: falloff_type !! Type of falloff
    real(dp), allocatable :: A_T !! Troe falloff parameter
    real(dp), allocatable :: T1 !! Troe falloff parameter
    real(dp), allocatable :: T2 !! Troe falloff parameter
    real(dp), allocatable :: T3 !! Troe falloff parameter

    type(Efficiencies) :: eff
  end type

  type :: MultiArrheniusRate
    real(dp), allocatable :: A(:)
    real(dp), allocatable :: b(:)
    real(dp), allocatable :: Ea(:)
  end type

  type, extends(BaseRate) :: PressDependentRate
    real(dp), allocatable :: logP(:)
    type(MultiArrheniusRate), allocatable :: rate(:)
  end type

  type :: Reaction
    integer :: nreact !! number of reactants
    integer :: nprod !! number of products
    integer, allocatable :: react_sp_inds(:) !! (nreact) species indexes for reactants
    integer, allocatable :: prod_sp_inds(:) !! (nprod) species indexes for products
    integer, allocatable :: reverse_info !! if a reversed reaction,
                                         !! then it is the index of the forward
    class(BaseRate), allocatable :: rp !! rate parameters
  end type

  type :: ProdLoss
    integer :: nump ! (iprod) number off production processes
    integer :: numl ! (iloss) number off loss processes
    integer, allocatable :: iprod(:) ! reaction #s of production mechanism for sp
    integer, allocatable :: iloss(:) ! reaction #s of loss mechanism for sp
  end type

  type :: PhotochemData
    ! PhotochemData contains information that is never changed
    ! after file read-in

    integer :: natoms !! number of atoms
    character(len=s_str_len), allocatable :: atoms_names(:) !! (natoms)
    real(dp), allocatable :: atoms_mass(:) !! g/mol (natoms)
    real(dp), allocatable :: atoms_redox(:) !!  (natoms)

    ! species
    ! Organization is as follows
    ! [     nsp     ]
    ! [   nq   + nsl]
    ! [np + ng + nsl]
    ! |_______|
    !     |
    ! Only np + ng = nq evolve through time. nsl are assumed to be in equilibrium.
    integer :: nq !! number of gases + particles which evolve over time from integration
    integer :: ng_1 !! index of first gas
    integer :: nll !! number of long-lived gas molecules
    integer :: nsl !! number of short-lived gas molecules
    integer :: ng  !! number of gases
    integer :: nsp !! total number of species (nq + nsl + 1)
    integer :: kd, kl, ku !! not read in. It is nq + nq + 1 (diagonal width of jacobian)
    integer :: lda !! not read in. It is nq + nq + nq + 1. leading dimension of array which stores jacobian
    character(len=s_str_len), allocatable :: SL_names(:) !! (nsl)
    character(len=s_str_len), allocatable :: species_names(:) !! (nsp+2) + 2 for hv and M
    integer, allocatable :: species_composition(:,:) !! (natoms, nsp+2)
    real(dp), allocatable :: species_mass(:) !! (nsp)
    real(dp), allocatable :: species_redox(:) !! (nsp)
    type(ThermodynamicData), allocatable :: thermo_data(:) !! (ng)
    real(dp), allocatable :: henry_data(:,:) !! (2, nsp).
    ! henry_data(:,i) = [A, B], and [mol/(kg * Pa)] = A*exp(B*(1.0_dp/298.15e0_dp - 1.0_dp/T))

    ! particles
    logical :: there_are_particles
    integer :: np !! number of particles
    integer :: npq !! number of particle equations. for now nq = npq.
    character(len=s_str_len), allocatable :: particle_names(:) !! np
    integer, allocatable :: particle_formation_method(:) !! np. 1 == saturation, 2 == reaction
    real(dp), allocatable :: particle_density(:) !! np (g/cm3)
    type(SaturationData), allocatable :: particle_sat(:) !! (np)
    character(len=s_str_len), allocatable :: particle_gas_phase(:) !! (np). gas phase of particle.
    ! Only for saturation particles
    integer, allocatable :: particle_gas_phase_ind(:) !! np. index of gas phase of particle
    integer, allocatable :: gas_particle_ind(:) !! (nq). Index of particle phase of gas
    character(len=s_str_len), allocatable :: particle_optical_prop(:) !! (np)
    integer, allocatable :: particle_optical_type(:) !! (np) 1 == mie, 2 == fractal

    ! reactions
    logical :: reverse !! True if there are reverse reactions
    integer :: nrF !! number of forward reactions
    integer :: nrR !! number of reverse reactions
    integer :: nrT !! number of total reactions
    type(Reaction), allocatable :: rx(:) !! (nrT) array of reaction objects
    character(len=m_str_len), allocatable :: reaction_equations(:) !! (nrT)
    type(ProdLoss), allocatable :: pl(:) !! (nsp) reactions producing and destroying each species
    integer :: kj !! number of photolysis reactions
    integer, allocatable :: photonums(:) !! (kj) the reaction number of each photolysis reaction

    ! radiative transfer
    integer :: nw !! number of wavelength bins
    real(dp), allocatable :: wavl(:) !! (nw+1) wavelength bins in nm
    type(XsectionData), allocatable :: absorp_xs(:) !! Contains photoabsorption species
    real(dp), allocatable :: photolysis_xs(:,:) !! (kj,nw) The xs time qy for each photolysis reaction
    integer :: nray !! number of species with rayleigh scattering
    real(dp), allocatable :: sigray(:,:) !! (len(raynums), nw)
    integer, allocatable :: raynums(:) !! species number of rayleigh species

    ! particle radiative transfer
    integer :: nrad_file
    real(dp), allocatable  :: radii_file(:,:) !! particle radii in optical data files
    ! We use array of types for particle xs because we want the option
    ! to exclude optical properties, but not take up a ton of useless memory.
    ! So some elements of this array have nothing in it.
    type(ParticleXsections), allocatable :: part_xs_file(:) !! np in length

    ! settings
    real(dp) :: planet_mass !! grams
    real(dp) :: planet_radius !! cm
    integer :: LH2O !! index of H2O; used by gas rainout
    logical :: gas_rainout !! True if gas rains out
    integer :: H_escape_type !! Diffusion-limited, Zahnle, or None
    real(dp) :: H_escape_coeff ! Coefficient for zahnle hydrogen escape
    integer :: LH2 !! H2 index
    integer :: LH !! H index
  end type

contains

  subroutine read_photochem_mechanism(mechanism_file, s, dat, err)
    character(len=*), intent(in) :: mechanism_file
    type(PhotoSettings), intent(in) :: s
    type(PhotochemData), intent(inout) :: dat
    character(:), allocatable, intent(out) :: err

    dat%nsl = s%nsl
    dat%SL_names = s%SL_names

    call get_photomech(mechanism_file, dat, err)
  end subroutine

  subroutine read_photochem_supporting_data(dat, data_dir, s, err)
    type(PhotochemData), intent(inout) :: dat
    character(len=*), intent(in) :: data_dir
    type(PhotoSettings), intent(in) :: s
    character(:), allocatable, intent(out) :: err

    call check_sl(dat, err)
    if (allocated(err)) return

    if (dat%gas_rainout) then
      call get_henry(dat, data_dir, s, err)
      if (allocated(err)) return
    endif

    call get_photorad(dat, data_dir, err)
  end subroutine

  subroutine get_photomech(infile, dat, err)
    use fortran_yaml_c, only : YamlFile
    character(len=*), intent(in) :: infile
    type(PhotochemData), intent(inout) :: dat
    character(:), allocatable, intent(out) :: err

    type(YamlFile) :: file

    ! parse yaml file
    call file%parse(infile, err)
    if (allocated(err)) return
    select type (root => file%root)
      class is (type_dictionary)
        call get_rxmechanism(root, infile, dat, err)
      class default
        err = "yaml file must have dictionaries at root level"
    end select
    call file%finalize()
    if (allocated(err)) return

  end subroutine

  subroutine get_rxmechanism(mapping, infile, dat, err)
    use photochem_enum, only: CondensingParticle, ReactionParticle
    use photochem_enum, only: MieParticle, FractalParticle
    use photochem_enum, only: PhotolysisRateType, ReverseRateType
    use clima_saturationdata, only: SaturationData
    class (type_dictionary), intent(in), pointer :: mapping
    character(len=*), intent(in) :: infile
    type(PhotochemData), target, intent(inout) :: dat
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

  subroutine get_henry_parse(root, henry_names, henry_data, err)
    class (type_list), intent(in) :: root
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

  subroutine get_henry(dat, data_dir, s, err)
    use photochem_settings, only: PhotoSettings
    use fortran_yaml_c, only : YamlFile
    type(PhotochemData), intent(inout) :: dat
    character(len=*), intent(in) :: data_dir
    type(PhotoSettings), intent(in) :: s
    character(:), allocatable :: err

    type(YamlFile) :: file
    integer :: j, ind, ind1, i

    character(len=s_str_len), allocatable :: henry_names(:)
    real(dp), allocatable :: henry_data(:,:)

    ! parse yaml file
    call file%parse(trim(data_dir)//"/henry/henry.yaml", err)
    if (allocated(err)) return
    select type (root => file%root)
    class is (type_list)
      call get_henry_parse(root, henry_names, henry_data, err)
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

  pure integer function count_reaction_tokens(tokens, token)
    character(len=*), intent(in) :: tokens(:)
    character(len=*), intent(in) :: token
    integer :: i

    count_reaction_tokens = 0
    do i = 1,size(tokens)
      if (trim(tokens(i)) == trim(token)) count_reaction_tokens = count_reaction_tokens + 1
    enddo
  end function

  subroutine compare_rxtype_string(tmp, eqr, eqp, reverse, rxtype_int, err)
    use photochem_enum, only: PhotolysisRateType, ElementaryRateType, ThreeBodyRateType, FalloffRateType, PressDependentRateType
    character(len=*), intent(in) :: tmp
    character(len=s_str_len), allocatable, intent(in) :: eqr(:), eqp(:)
    logical, intent(in) :: reverse
    integer, intent(in) :: rxtype_int
    character(:), allocatable, intent(out) :: err
    character(:), allocatable :: rxtype
    integer :: n_m_react, n_m_prod, n_hv_react, n_hv_prod
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

    n_m_react = count_reaction_tokens(eqr, 'M')
    n_m_prod = count_reaction_tokens(eqp, 'M')
    n_hv_react = count_reaction_tokens(eqr, 'hv')
    n_hv_prod = count_reaction_tokens(eqp, 'hv')

    if (trim(rxtype) == 'three-body' .or. trim(rxtype) == 'falloff') then
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
      if (n_m_react /= 1 .or. n_m_prod /= 1) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can only have one "M" on either side'
        return
      endif
      if (n_hv_react > 0 .or. n_hv_prod > 0) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can not contain "hv". Only photolysis reactions can.'
        return
      endif
    elseif (trim(rxtype) == 'elementary' .or. &
            trim(rxtype) == 'pressure-dependent-Arrhenius') then
      if (n_m_react > 0 .or. n_m_prod > 0) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can not contain "M".'
        return
      endif
      if (n_hv_react > 0 .or. n_hv_prod > 0) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can not contain "hv". Only photolysis reactions can.'
        return
      endif
    elseif (trim(rxtype) == 'photolysis') then
      if (reverse) then
        err = 'IOError: Photolysis reaction '//trim(tmp)//' can not be reversed.'
        return
      endif

      if (n_m_react > 0 .or. n_m_prod > 0) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can not contain "M".'
        return
      endif
      if (n_hv_react > 1) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' can only have one "hv" on the left side.'
        return
      endif
      if (n_hv_react /= 1 .or. n_hv_prod /= 0) then
        err = 'IOError: '//trim(rxtype)// ' reaction '//trim(tmp)// &
                ' must have "hv" on the left and no "hv" on the right.'
        return
      endif
    endif
  end subroutine

  subroutine get_thermodata(molecule, molecule_name, infile, &
                            thermo, err)
    use photochem_enum, only: ShomatePolynomial, Nasa9Polynomial, Nasa7Polynomial
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
    use photochem_enum, only: ThreeBodyRateType, FalloffRateType
    type(PhotochemData), intent(in) :: dat
    character(len=*), intent(in) :: rx_str
    type(Reaction), intent(inout) :: rx ! already has rate parameters
    logical, intent(out) :: reverse
    character(:), allocatable, intent(out) :: err

    character(len=s_str_len), allocatable :: eqr(:), eqp(:)
    integer :: reactant_atoms(dat%natoms), product_atoms(dat%natoms)
    integer :: ind(1), i


    call parse_reaction(rx_str, reverse, eqr, eqp, err)
    if (allocated(err)) return
    call compare_rxtype_string(rx_str, eqr, eqp, reverse, rx%rp%rxtype, err)
    if (allocated(err)) return

    if (rx%rp%rxtype == ThreeBodyRateType .or. rx%rp%rxtype == FalloffRateType) then
      ! remove the M
      eqr = eqr(1:size(eqr)-1)
      eqp = eqp(1:size(eqp)-1)
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
    character(:), allocatable :: lhs, rhs

    call split_reaction_arrow(rx_string, reverse, lhs, rhs, err)

  end function

  subroutine get_rateparams(dat, reaction_d, infile, rx, err)
    use photochem_enum, only: NoFalloff, TroeWithoutT2Falloff, TroeWithT2Falloff, JPLFalloff
    use photochem_enum, only: PhotolysisRateType, ElementaryRateType, ThreeBodyRateType, FalloffRateType, PressDependentRateType
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
    character(:), allocatable :: lhs, rhs

    call split_reaction_arrow(instring, reverse, lhs, rhs, err)
    if (allocated(err)) return

    call tokenize_reaction_side(lhs, eqr, err)
    if (allocated(err)) return

    call tokenize_reaction_side(rhs, eqp, err)
  end subroutine

  subroutine split_reaction_arrow(instring, reverse, lhs, rhs, err)
    character(len=*), intent(in) :: instring
    logical, intent(out) :: reverse
    character(:), allocatable, intent(out) :: lhs, rhs
    character(:), allocatable, intent(out) :: err

    integer :: i, n, arrow_start, arrow_length, arrow_count

    reverse = .false.
    lhs = ''
    rhs = ''
    arrow_start = 0
    arrow_length = 0
    arrow_count = 0
    n = len_trim(instring)

    if (n == 0) then
      err = 'IOError: Reaction equation is empty.'
      return
    endif

    i = 1
    do while (i <= n)
      if (i <= n-2) then
        if (instring(i:i+2) == '<=>') then
          arrow_count = arrow_count + 1
          arrow_start = i
          arrow_length = 3
          reverse = .true.
          i = i + 3
          cycle
        endif
      endif

      if (i <= n-1) then
        if (instring(i:i+1) == '=>') then
          arrow_count = arrow_count + 1
          arrow_start = i
          arrow_length = 2
          reverse = .false.
          i = i + 2
          cycle
        endif
      endif

      i = i + 1
    enddo

    if (arrow_count == 0) then
      err = 'IOError: Reaction "'//trim(instring)//'" has no valid arrow. '// &
            'Expected "=>" or "<=>".'
      return
    elseif (arrow_count > 1) then
      err = 'IOError: Reaction "'//trim(instring)//'" has more than one reaction arrow.'
      return
    endif

    if (arrow_start > 1) lhs = trim(instring(:arrow_start-1))
    if (arrow_start + arrow_length <= n) then
      rhs = trim(instring(arrow_start+arrow_length:n))
    endif

    if (len_trim(lhs) == 0) then
      err = 'IOError: Reaction "'//trim(instring)//'" has no reactants.'
      return
    endif
    if (len_trim(rhs) == 0) then
      err = 'IOError: Reaction "'//trim(instring)//'" has no products.'
      return
    endif
  end subroutine

  subroutine tokenize_reaction_side(side, tokens, err)
    character(len=*), intent(in) :: side
    character(len=s_str_len), allocatable, intent(out) :: tokens(:)
    character(:), allocatable, intent(out) :: err

    integer :: i, n, token_start, parentheses_depth
    logical :: term_complete, expect_term
    character :: c

    allocate(tokens(0))
    n = len_trim(side)
    token_start = 0
    parentheses_depth = 0
    term_complete = .false.
    expect_term = .false.

    do i = 1,n
      c = side(i:i)

      if (c == '+' .and. is_reaction_plus_marker(side, i)) then
        ! In parenthesized third-body notation, '+' is a marker, not a
        ! separator: both (+M) and (+ M) represent the species M.
        expect_term = .true.
        term_complete = .false.
        cycle
      endif

      if (c == '+' .and. is_reaction_plus_separator(side, i)) then
        if (.not. term_complete) then
          err = 'IOError: Reaction side "'//trim(side)// &
                '" has a plus sign without a species on its left.'
          return
        endif
        expect_term = .true.
        term_complete = .false.
        cycle
      endif

      if (c == '(') parentheses_depth = parentheses_depth + 1
      if (c == ')') then
        if (parentheses_depth == 0) then
          err = 'IOError: Reaction side "'//trim(side)//'" has an unmatched closing parenthesis.'
          return
        endif
        parentheses_depth = parentheses_depth - 1
      endif

      if (is_reaction_space(c) .or. c == '(' .or. c == ')') then
        if (token_start /= 0) then
          call append_reaction_token(tokens, side, token_start, i-1, err)
          if (allocated(err)) return
          token_start = 0
          term_complete = .true.
        endif
        cycle
      endif

      if (token_start == 0) then
        if (term_complete .and. .not. expect_term) then
          err = 'IOError: Reaction side "'//trim(side)// &
                '" has adjacent species without a "+" separator.'
          return
        endif
        token_start = i
        term_complete = .false.
        expect_term = .false.
      endif
    enddo

    if (token_start /= 0) then
      call append_reaction_token(tokens, side, token_start, n, err)
      if (allocated(err)) return
      term_complete = .true.
    endif

    if (parentheses_depth /= 0) then
      err = 'IOError: Reaction side "'//trim(side)//'" has an unmatched opening parenthesis.'
      return
    endif

    if (expect_term) then
      err = 'IOError: Reaction side "'//trim(side)//'" ends with a plus sign.'
      return
    endif
    if (size(tokens) == 0) then
      err = 'IOError: Reaction side "'//trim(side)//'" contains no species.'
    endif
  end subroutine

  subroutine append_reaction_token(tokens, side, first, last, err)
    character(len=s_str_len), allocatable, intent(inout) :: tokens(:)
    character(len=*), intent(in) :: side
    integer, intent(in) :: first, last
    character(:), allocatable, intent(out) :: err

    character(:), allocatable :: token

    if (last < first) return
    token = trim(side(first:last))
    if (len_trim(token) == 0) return
    if (len(token) > s_str_len) then
      err = 'IOError: Species token "'//trim(token)//'" in reaction side "'// &
            trim(side)//'" is longer than the supported species-name length.'
      return
    endif

    if (size(tokens) == 0) then
      tokens = [repeat(' ', s_str_len)]
    else
      tokens = [tokens, repeat(' ', s_str_len)]
    endif
    tokens(size(tokens)) = token
  end subroutine

  pure logical function is_reaction_space(c)
    character, intent(in) :: c
    is_reaction_space = c == ' ' .or. c == achar(9)
  end function

  pure logical function is_reaction_plus_separator(side, position)
    character(len=*), intent(in) :: side
    integer, intent(in) :: position
    integer :: n
    logical :: left_space, right_space

    n = len_trim(side)
    if (position == 1) then
      left_space = .true.
    else
      left_space = is_reaction_space(side(position-1:position-1))
    endif
    if (position == n) then
      right_space = .true.
    else
      right_space = is_reaction_space(side(position+1:position+1))
    endif

    ! Parenthesized third-body notation commonly uses (+M), (+ M), or ( + M).
    is_reaction_plus_separator = left_space .and. right_space
    if (position > 1) then
      is_reaction_plus_separator = is_reaction_plus_separator .or. &
                                   side(position-1:position-1) == '('
    endif
    if (position < n) then
      is_reaction_plus_separator = is_reaction_plus_separator .or. &
                                   side(position+1:position+1) == ')'
    endif
  end function

  pure logical function is_reaction_plus_marker(side, position)
    character(len=*), intent(in) :: side
    integer, intent(in) :: position
    integer :: i

    i = position - 1
    do while (i >= 1)
      if (.not. is_reaction_space(side(i:i))) exit
      i = i - 1
    enddo
    is_reaction_plus_marker = .false.
    if (i >= 1) is_reaction_plus_marker = side(i:i) == '('
  end function

  subroutine check_sl(dat, err)
    type(PhotochemData), intent(in) :: dat
    character(:), allocatable, intent(out) :: err

    integer :: i, j, l, k, kk, mm, n, nn, ind(1), counter
    character(len=:), allocatable :: reaction_equation

    do i = 1, dat%nsl

      j = dat%nq + i
      ! can not be an efficiency.
      do k = 1,dat%nrF
        select type(rp => dat%rx(k)%rp)
        class is (ThreeBodyRate)
          if (rp%eff%n_eff > 0) then
            ind = findloc(rp%eff%eff_sp_inds,j)
            if (ind(1) /= 0) then
              call reaction_string(dat,k,reaction_equation)
              err = 'IOError: Reaction "'//reaction_equation//'" has short-lived species collision efficiencies.' // &
              ' This is not allowed. Either remove the efficiencies, or change the species to long lived.'
              return
            endif
          endif
        class is (FalloffRate)
          if (rp%eff%n_eff > 0) then
            ind = findloc(rp%eff%eff_sp_inds, j)
            if (ind(1) /= 0) then
              call reaction_string(dat,k,reaction_equation)
              err = 'IOError: Reaction "'//reaction_equation//'" has short-lived species collision efficiencies.' // &
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
              call reaction_string(dat, kk, reaction_equation)
              err = 'IOError: Reaction "'//reaction_equation//'" has short-lived species as reactants'// &
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
                call reaction_string(dat, kk, reaction_equation)
                err = 'IOError: Reaction "'//reaction_equation//'" short lived species react'// &
                ' with themselves. This is not allowed. Change the species to long lived.'
                return
              endif
            elseif (nn == n .and. n /= j) then
              call reaction_string(dat, kk, reaction_equation)
              err = 'IOError: Reaction "'//reaction_equation//'" short lived species react'// &
              ' with other short lived species. This is not allowed.'
              return
            elseif (dat%species_names(n) == 'hv') then
              call reaction_string(dat, kk, reaction_equation)
              err = 'IOError: Photolysis reaction "'//reaction_equation//'" can not have short lived species.'
              return
            endif
          enddo
        enddo
      enddo

    enddo

  end subroutine

  subroutine get_photorad(dat, data_dir, err)
    type(PhotochemData), intent(inout) :: dat
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err

    ! Read wavelength grid
    call get_wavl(dat, data_dir, err)
    if (allocated(err)) return

    ! get rayleigh
    call get_rayleigh(dat, data_dir, err)
    if (allocated(err)) return

    ! get photolysis xsections data
    call get_photolysis_xs(dat, data_dir, err)
    if (allocated(err)) return

    if (dat%there_are_particles) then
      call get_aerosol_xs(dat, data_dir, err)
      if (allocated(err)) return
    endif

  end subroutine

  subroutine get_wavl(dat, data_dir, err)
    use h5fortran
    type(PhotochemData), intent(inout) :: dat
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err

    type(hdf5_file) :: h
    integer(HSIZE_T), allocatable :: dims(:)
    character(:), allocatable :: xsroot, filename

    ! Root directory
    xsroot = trim(data_dir)//"/xsections/"
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
    use ieee_arithmetic, only: ieee_is_finite
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

    if (nrad_file < 2) then
      err = '"radii" must contain at least two values in "'//filename//'"'
      return
    endif
    if (.not. all(ieee_is_finite(radii_file))) then
      err = '"radii" contains a non-finite value in "'//filename//'"'
      return
    endif
    if (any(radii_file <= 0.0_dp)) then
      err = '"radii" must contain only positive values in "'//filename//'"'
      return
    endif
    do i = 1,nrad_file-1
      if (radii_file(i+1) <= radii_file(i)) then
        err = '"radii" must be strictly increasing in "'//filename//'"'
        return
      endif
    enddo

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

    call validate_mie_dataset(w0_tmp, 'w0', filename, 0.0_dp, 1.0_dp, err)
    if (allocated(err)) return
    call validate_mie_dataset(qext_tmp, 'qext', filename, 0.0_dp, err=err)
    if (allocated(err)) return
    call validate_mie_dataset(g0_tmp, 'g0', filename, -1.0_dp, 1.0_dp, err)
    if (allocated(err)) return

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

  subroutine get_aerosol_xs(dat, data_dir, err)
    use photochem_enum, only: MieParticle, FractalParticle
    type(PhotochemData), intent(inout) :: dat
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err

    integer :: nrad
    integer, parameter :: nrad_fixed = 50
    real(dp), allocatable :: radii(:)
    real(dp), allocatable :: w0(:,:), qext(:,:), g(:,:)
    character(len=:), allocatable :: xsroot
    character(len=:), allocatable :: filename
    integer :: i

    xsroot = trim(data_dir)//"/aerosol_xsections/"

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

  subroutine get_photolysis_xs(dat, data_dir, err)
    use h5fortran, only: is_hdf5
    type(PhotochemData), intent(inout) :: dat
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err

    character(:), allocatable :: xsroot, species, filename, reaction_equation
    integer :: i, j, k

    ! Root directory
    xsroot = trim(data_dir)//"/xsections/"

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
      call reaction_string(dat, j, reaction_equation)

      ! Get the absorbing species
      k = dat%rx(j)%react_sp_inds(1)
      species = trim(dat%species_names(k))

      ! Filename is based on the species
      filename = xsroot//species//'.h5'

      call get_photolysis_xs_for_reaction(filename, reaction_equation, dat%wavl, dat%photolysis_xs(i,:), err)
      if (allocated(err)) return

    enddo

  end subroutine

  subroutine get_photolysis_xs_for_reaction(filename, reaction_equation, wavl, xs, err)
    use futils, only: interp
    use h5fortran
    use futils, only: interp_discrete_to_bins
    use clima_useful, only: hdf5_file_closer
    character(*), intent(in) :: filename
    character(*), intent(in) :: reaction_equation
    real(dp), intent(in) :: wavl(:)
    real(dp), intent(out) :: xs(:)
    character(:), allocatable, intent(out) :: err

    type(hdf5_file), target :: h
    integer(HSIZE_T), allocatable :: dims(:)
    real(dp), allocatable :: wv_tmp(:), xs_tmp(:), wv1(:), qy1(:), qy2(:), wv_av(:)
    integer :: ierr
    character(:), allocatable :: err_base, photo_type

    err_base = 'Failed to retrieve photolysis cross section data for the reaction "'//reaction_equation//'". '

    block
    type(hdf5_file_closer) :: h_closer

    if (.not. is_hdf5(filename)) then
      err = err_base//'The file "'//filename//'" does not exist.'
      return
    endif

    ! Open file
    call h%open(filename,'r')
    h_closer%h => h

    ! Find if reaction is photolysis or ionization
    if (h%exists("photodissociation-qy")) then
      call h%open_group("photodissociation-qy")
      if (h%exists(reaction_equation)) photo_type = 'photodissociation'
      call h%close_group()
    endif
    if (.not.allocated(photo_type) .and. h%exists("photoionisation-qy")) then
      call h%open_group("photoionisation-qy")
      if (h%exists(reaction_equation)) photo_type = 'photoionisation'
      call h%close_group()
    endif
    if (.not.allocated(photo_type)) then
      err = err_base//filename//': could not find quantum yield information.'
      return
    endif

    ! Wavelengths
    call check_h5_dataset(h, 'wavelengths', 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape("wavelengths", dims)
    allocate(wv_tmp(dims(1)))
    call h%read("wavelengths", wv_tmp)

    ! Photodissociation
    call check_h5_dataset(h, photo_type, 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) return
    call h%shape(photo_type, dims)
    allocate(xs_tmp(dims(1)))
    call h%read(photo_type, xs_tmp)

    ! Interpolate
    call interp_discrete_to_bins(wavl, wv_tmp, xs_tmp, xs, 'FillValue', tiny(1.0_dp), err)
    if (allocated(err)) return

    ! Now do the quantum yields
    if (.not. h%exists(photo_type//"-qy")) then
      err = err_base//filename//': dataset "'//photo_type//'-qy" does not exist'
      return
    endif

    call h%open_group(photo_type//"-qy")

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
    call check_h5_dataset(h, reaction_equation, 1, H5T_FLOAT_F, filename, err)
    if (allocated(err)) then
      call h%close_group()
      return
    endif
    call h%shape(reaction_equation, dims)
    allocate(qy1(dims(1)))
    call h%read(reaction_equation, qy1)

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

  subroutine validate_mie_dataset(values, dataset, filename, lower_bound, &
                                  upper_bound, err)
    use ieee_arithmetic, only: ieee_is_finite

    real(dp), intent(in) :: values(:,:)
    character(*), intent(in) :: dataset, filename
    real(dp), intent(in) :: lower_bound
    real(dp), intent(in), optional :: upper_bound
    character(:), allocatable, intent(out) :: err

    if (.not. all(ieee_is_finite(values))) then
      err = '"'//trim(dataset)//'" contains a non-finite value in "'// &
            trim(filename)//'"'
      return
    endif

    if (any(values < lower_bound)) then
      err = '"'//trim(dataset)//'" contains a value below the allowed '// &
            'range in "'//trim(filename)//'"'
      return
    endif

    if (present(upper_bound)) then
      if (any(values > upper_bound)) then
        err = '"'//trim(dataset)//'" contains a value above the allowed '// &
              'range in "'//trim(filename)//'"'
        return
      endif
    endif
  end subroutine

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

  subroutine get_rayleigh(dat, data_dir, err)
    use fortran_yaml_c, only : YamlFile
    type(PhotochemData), intent(inout) :: dat
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: A(:), B(:), Delta(:)
    character(len=str_len) :: rayleigh_file

    type(YamlFile) :: file
    integer :: i, j

    rayleigh_file = trim(data_dir)//"/rayleigh/rayleigh.yaml"

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
  subroutine list_destroy(self)
    type(type_list_tmp), intent(inout) :: self

    type(type_list_item), pointer :: item, next

    item => self%first
    do while (associated(item))
      next => item%next
      deallocate(item)
      item => next
    enddo
    nullify(self%first)
  end subroutine

end module
