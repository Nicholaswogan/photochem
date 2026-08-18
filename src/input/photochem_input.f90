
module photochem_input
  use fortran_yaml_c_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  use photochem_types, only : PhotochemData, PhotochemVars, PhotoSettings, &
                              AtmosphereState, AtmosphereStateDerived
  use photochem_const, only: dp, str_len, s_str_len
  implicit none
  private 

  public :: setup_static, setup_atmosphere_from_file, map_atmosphere_z_to_grid
  public :: map_atmosphere_p_to_grid
  public :: finalize_atmosphere_initialization
  public :: finalize_atmosphere_state
  public :: refresh_temperature_dependent_state
  public :: parse_reaction
  
  type, extends(type_list) :: type_list_tmp
  ! temporary list for accessing all reactions and
  ! species in a row.
  contains
    final :: list_destroy
  end type

  ! Raw atmosphere-file data is initialization-only state. Keep it in a
  ! short-lived profile object rather than retaining it in PhotochemData.
  type :: AtmosphereFileProfile
    integer :: nlayer = 0
    real(dp), allocatable :: z(:) !! Layer-center altitude (cm)
    real(dp), allocatable :: temperature(:) !! Layer-center temperature (K)
    real(dp), allocatable :: edd(:) !! Layer-center eddy diffusion (cm^2/s)
    real(dp), allocatable :: density(:) !! Layer-center total density (molecules/cm^3)
    real(dp), allocatable :: mix(:,:) !! Evolved-species mixing ratios
    real(dp), allocatable :: particle_radius(:,:) !! Particle radii (cm)
  end type
  
  interface
    !> Complete derived atmosphere setup after a profile has been mapped onto
    !! the model grid and the static model state has been loaded.
    module subroutine finalize_atmosphere_initialization(dat, var, err)
      type(PhotochemData), intent(inout) :: dat
      type(PhotochemVars), intent(inout) :: var
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Complete derived atmosphere setup on a standalone atmospheric state.
    !! This is the transaction-friendly form used by initialization paths that
    !! must validate all atmospheric changes before committing them to `var`.
    module subroutine finalize_atmosphere_state(dat, state, derived, err)
      type(PhotochemData), intent(in) :: dat
      type(AtmosphereState), intent(inout) :: state
      type(AtmosphereStateDerived), intent(inout) :: derived
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Refresh state derived from temperature and the current tropopause.
    !! `var%temperature` must already contain the proposed temperature. If
    !! `trop_alt` is supplied, it replaces `var%trop_alt` when gas rainout is
    !! enabled.
    module subroutine refresh_temperature_dependent_state(dat, var, trop_alt, err)
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
      real(dp), optional, intent(in) :: trop_alt
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Map a temporary legacy atmosphere-file profile onto the common altitude
    !! grid. File density remains an explicit legacy policy; it is not replaced
    !! by hydrostatic reconstruction in this adapter.
    module subroutine map_atmosphere_file_to_grid(dat, var, profile, usol, err)
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
      type(AtmosphereFileProfile), intent(in) :: profile
      real(dp), intent(out) :: usol(:,:)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    !> Read files and settings that do not require an initialized atmosphere.
    module subroutine read_static_files(mechanism_file, s, flux_file, dat, var, err)
      character(len=*), intent(in) :: mechanism_file
      type(PhotoSettings), intent(in) :: s
      character(len=*), intent(in) :: flux_file
      type(PhotochemData), intent(inout) :: dat
      type(PhotochemVars), intent(inout) :: var
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Read raw profiles from a legacy atmosphere text file.
    module subroutine read_atmosphere_file(atmosphere_txt, dat, profile, err)
      character(len=*), intent(in) :: atmosphere_txt
      type(PhotochemData), intent(in) :: dat
      type(AtmosphereFileProfile), intent(out) :: profile
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Resolve and validate settings that depend on raw atmosphere-file data.
    module subroutine resolve_atmosphere_settings(profile, dat, var, err)
      type(AtmosphereFileProfile), intent(in) :: profile
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
      character(:), allocatable, intent(out) :: err
    end subroutine

    module subroutine parse_reaction(instring, reverse, eqr, eqp, err)
      character(len=*), intent(in) :: instring
      logical, intent(out) :: reverse
      character(len=s_str_len), allocatable, intent(out) :: eqr(:), eqp(:)
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Map altitude-based atmospheric inputs onto the model grid and construct
    !! a hydrostatic initial state. This is an internal initialization kernel;
    !! species-name handling and transactional model updates are performed by
    !! EvoAtmosphere.
    module subroutine map_atmosphere_z_to_grid(dat, nz, trop_alt, z, temperature, &
                                               edd, surface_pressure, mix, &
                                               particle_radius, state, derived, err)
      type(PhotochemData), intent(in) :: dat
      integer, intent(in) :: nz
      real(dp), intent(in) :: trop_alt
      real(dp), intent(in) :: z(:) !! Altitude profile knots (cm), including both domain edges.
      real(dp), intent(in) :: temperature(:) !! Temperature at `z` (K).
      real(dp), intent(in) :: edd(:) !! Eddy diffusion at `z` (cm^2/s).
      real(dp), intent(in) :: surface_pressure !! Pressure at the lower domain edge (dyn/cm^2).
      !> Mixing ratios in evolved-species order (`1:dat%nq`) at `z`.
      !! Particle species precede evolved gases, matching the first dimension
      !! of `usol`.
      real(dp), intent(in) :: mix(:,:)
      !> Particle radii in mechanism order (`1:dat%npq`) at `z` (cm).
      real(dp), intent(in) :: particle_radius(:,:)
      type(AtmosphereState), intent(inout) :: state
      type(AtmosphereStateDerived), intent(inout) :: derived
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Map pressure-based atmospheric inputs onto the model altitude grid and
    !! construct a hydrostatic initial state. The first and last pressure
    !! points define the lower and upper boundaries of the model domain. This
    !! is an internal initialization kernel; transactional model updates are
    !! performed by EvoAtmosphere. When supplied, `trop_p` is converted to the
    !! corresponding altitude before the common altitude mapper validates and
    !! finalizes the candidate state.
    module subroutine map_atmosphere_p_to_grid(dat, var, profile_pressure, &
                                               temperature, edd, mix, &
                                               particle_radius, trop_p, &
                                               pressure, density, mubar, usol, err)
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
      !> Strictly decreasing pressure profile knots (dyn/cm^2), including both domain edges.
      real(dp), intent(in) :: profile_pressure(:)
      real(dp), intent(in) :: temperature(:) !! Temperature at `profile_pressure` (K).
      real(dp), intent(in) :: edd(:) !! Eddy diffusion at `profile_pressure` (cm^2/s).
      !> Mixing ratios in evolved-species order (`1:dat%nq`) at `profile_pressure`.
      real(dp), intent(in) :: mix(:,:)
      !> Particle radii in mechanism order (`1:dat%npq`) at `profile_pressure` (cm).
      real(dp), intent(in) :: particle_radius(:,:)
      !> Optional tropopause pressure (dyn/cm^2), used when gas rainout is enabled.
      real(dp), optional, intent(in) :: trop_p
      real(dp), intent(out) :: pressure(:) !! Hydrostatic pressure at model centers (dyn/cm^2).
      real(dp), intent(out) :: density(:) !! Total gas number density at model centers (molecules/cm^3).
      real(dp), intent(out) :: mubar(:) !! Gas mean molecular weight at model centers (g/mol).
      real(dp), intent(out) :: usol(:,:) !! Initial number densities (molecules/cm^3).
      character(:), allocatable, intent(out) :: err
    end subroutine
    
  end interface
    
contains

  !> Read and allocate model state that is independent of atmospheric profiles.
  !!
  !! This establishes mechanism and radiative dimensions, reads settings and
  !! stellar data, and allocates persistent arrays whose shapes depend only on
  !! those dimensions and `number-of-layers`. It does not construct a vertical
  !! grid, initialize atmospheric profiles, or prepare an RHS state.
  subroutine setup_static(mechanism_file, s, flux_file, dat, var, err)

    character(len=*), intent(in) :: mechanism_file
    type(PhotoSettings), intent(in) :: s
    character(len=*), intent(in) :: flux_file
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err

    call read_static_files(mechanism_file, s, flux_file, dat, var, err)
    if (allocated(err)) return

    dat%kd = 2*dat%nq + 1
    dat%kl = dat%kd + dat%nq
    dat%ku = dat%kd - dat%nq
    dat%lda = 3*dat%nq + 1

    call allocate_nz_vars(dat, var)

  end subroutine

  subroutine allocate_nz_vars(dat, var)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    
    integer :: i
    
    var%neqs = dat%nq*var%nz

    allocate(var%temperature(var%nz))
    allocate(var%z(var%nz))
    allocate(var%dz(var%nz))
    allocate(var%edd(var%nz))
    allocate(var%grav(var%nz))
    allocate(var%particle_radius(dat%npq,var%nz))
    allocate(var%xs_x_qy(var%nz,dat%kj,dat%nw))
    
    allocate(var%particle_xs(dat%np))
    do i = 1,dat%np
      ! only allocate space if there is data
      if (dat%part_xs_file(i)%ThereIsData) then
        var%particle_xs(i)%ThereIsData = .true.
        allocate(var%particle_xs(i)%w0(var%nz,dat%nw))
        allocate(var%particle_xs(i)%qext(var%nz,dat%nw))
        allocate(var%particle_xs(i)%gt(var%nz,dat%nw))
      else
        var%particle_xs(i)%ThereIsData = .false.
      endif
    enddo
    
    if (dat%reverse) then
      allocate(var%gibbs_energy(var%nz,dat%ng))
    endif

    allocate(var%tauc(var%nz,dat%nw))
    var%tauc = 0.0_dp
    allocate(var%w0c(var%nz,dat%nw))
    var%w0c = 0.0_dp
    allocate(var%g0c(var%nz,dat%nw))
    var%g0c = 0.0_dp

  end subroutine

  !> Read and map a legacy atmosphere text file onto the model grid.
  !!
  !! [[setup_static]] must have completed successfully before this routine is
  !! called. This routine reads the raw atmosphere, resolves grid-dependent
  !! settings, constructs the grid and gravity, and interpolates atmospheric
  !! profiles and particle radii. The caller must invoke
  !! `finalize_atmosphere_initialization` before preparing the RHS work arrays.
  subroutine setup_atmosphere_from_file(atmosphere_txt, dat, var, usol, err)

    character(len=*), intent(in) :: atmosphere_txt
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    real(dp), intent(out) :: usol(:,:)
    character(:), allocatable, intent(out) :: err

    type(AtmosphereFileProfile) :: profile

    call read_atmosphere_file(atmosphere_txt, dat, profile, err)
    if (allocated(err)) return

    call map_atmosphere_file_to_grid(dat, var, profile, usol, err)
    if (allocated(err)) return

  end subroutine
  
  subroutine list_destroy(self)
    type(type_list_tmp), intent(inout) :: self
    
    type (type_list_item),pointer :: item, next
    
    item => self%first
    do while (associated(item))
       next => item%next
       deallocate(item)
       item => next
    end do
    nullify(self%first)
  end subroutine
  
end module
