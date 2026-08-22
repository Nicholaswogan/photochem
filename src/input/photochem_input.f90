
module photochem_input
  use photochem_data, only: PhotochemData, read_photochem_mechanism, &
                            read_photochem_supporting_data, parse_reaction
  use photochem_settings, only: PhotoSettings
  use photochem_vars, only: PhotochemVars, apply_vars_settings, &
                            allocate_model_grid, read_stellar_flux
  use photochem_types, only: AtmosphereState
  use photochem_const, only: dp, s_str_len
  implicit none
  private

  public :: setup_static, map_atmosphere_file_to_grid, map_atmosphere_z_to_grid
  public :: map_atmosphere_p_to_grid
  public :: finalize_atmosphere_initialization
  public :: finalize_atmosphere_state
  public :: refresh_temperature_dependent_vars
  public :: parse_reaction

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
    module subroutine finalize_atmosphere_state(dat, state, err)
      type(PhotochemData), intent(in) :: dat
      type(AtmosphereState), intent(inout) :: state
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Refresh temperature-dependent quantities without modifying
    !! `PhotochemVars` directly.
    module subroutine refresh_temperature_dependent_vars(dat, temperature, z, &
                                                         bottom_atmos, top_atmos, &
                                                         trop_alt_new, xs_x_qy, &
                                                         gibbs_energy, trop_alt, &
                                                         trop_ind, err)
      type(PhotochemData), intent(in) :: dat
      real(dp), intent(in) :: temperature(:), z(:)
      real(dp), intent(in) :: bottom_atmos, top_atmos
      real(dp), optional, intent(in) :: trop_alt_new
      real(dp), intent(inout) :: xs_x_qy(:,:,:)
      real(dp), allocatable, intent(inout) :: gibbs_energy(:,:)
      real(dp), intent(inout) :: trop_alt
      integer, intent(inout) :: trop_ind
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Read a legacy atmosphere file and map it onto the common altitude grid.
    !! File density remains an explicit legacy policy; it is not replaced by
    !! hydrostatic reconstruction in this adapter.
    module subroutine map_atmosphere_file_to_grid(dat, atmosphere_txt, trop_alt_default, state, err)
      type(PhotochemData), intent(in) :: dat
      character(len=*), intent(in) :: atmosphere_txt
      real(dp), intent(in) :: trop_alt_default
      type(AtmosphereState), intent(inout) :: state
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Read files and settings that do not require an initialized atmosphere.
    module subroutine read_static_files(mechanism_file, s, flux_file, data_dir, &
                                        dat, var, err)
      character(len=*), intent(in) :: mechanism_file
      type(PhotoSettings), intent(in) :: s
      character(len=*), intent(in) :: flux_file
      character(len=*), intent(in) :: data_dir
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

    !> Map altitude-based atmospheric inputs onto the model grid and construct
    !! a hydrostatic initial state. This is an internal initialization kernel;
    !! species-name handling and transactional model updates are performed by
    !! EvoAtmosphere.
    module subroutine map_atmosphere_z_to_grid(dat, nz, trop_alt, z, temperature, &
                                               edd, surface_pressure, mix, &
                                               particle_radius, state, err)
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
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Map pressure-based atmospheric inputs onto the model altitude grid and
    !! construct a hydrostatic initial state. The first and last pressure
    !! points define the lower and upper boundaries of the model domain. This
    !! is an internal initialization kernel; transactional model updates are
    !! performed by EvoAtmosphere. When supplied, `trop_p` is converted to the
    !! corresponding altitude before the common altitude mapper validates and
    !! finalizes the supplied `AtmosphereState`.
    module subroutine map_atmosphere_p_to_grid(dat, nz, trop_alt_default, &
                                               profile_pressure, temperature, &
                                               edd, mix, particle_radius, trop_p, &
                                               state, err)
      type(PhotochemData), intent(in) :: dat
      integer, intent(in) :: nz
      real(dp), intent(in) :: trop_alt_default
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
      type(AtmosphereState), intent(inout) :: state
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
  subroutine setup_static(mechanism_file, s, flux_file, data_dir, dat, var, err)

    character(len=*), intent(in) :: mechanism_file
    type(PhotoSettings), intent(in) :: s
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: data_dir
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err

    call read_static_files(mechanism_file, s, flux_file, data_dir, dat, var, err)
    if (allocated(err)) return

    dat%kd = 2*dat%nq + 1
    dat%kl = dat%kd + dat%nq
    dat%ku = dat%kd - dat%nq
    dat%lda = 3*dat%nq + 1

    call allocate_model_grid(dat, var)

  end subroutine

end module
