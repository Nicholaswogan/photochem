
module photochem_input
  use fortran_yaml_c_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  use photochem_types, only : PhotochemData, PhotochemVars, PhotoSettings
  use photochem_const, only: dp, str_len, s_str_len
  implicit none
  private 

  public :: setup_static, setup_atmosphere_from_file
  public :: interp2xsdata, compute_gibbs_energy, interp2particlexsdata, parse_reaction
  
  type, extends(type_list) :: type_list_tmp
  ! temporary list for accessing all reactions and
  ! species in a row.
  contains
    final :: list_destroy
  end type
  
  interface
    !> Construct and prepare atmosphere-dependent state after raw atmospheric
    !! profiles and the static model state have been loaded.
    module subroutine after_read_setup(dat, var, err)
      use photochem_eqns, only: vertical_grid, gravity
      type(PhotochemData), intent(inout) :: dat
      type(PhotochemVars), intent(inout) :: var
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
    module subroutine read_atmosphere_file(atmosphere_txt, dat, var, err)
      character(len=*), intent(in) :: atmosphere_txt
      type(PhotochemData), intent(inout) :: dat
      type(PhotochemVars), intent(inout) :: var
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Resolve and validate settings that depend on raw atmosphere-file data.
    module subroutine resolve_atmosphere_settings(dat, var, err)
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
      character(:), allocatable, intent(out) :: err
    end subroutine

    !> Allocate persistent arrays whose dimensions depend on the mechanism and nz.
    module subroutine allocate_nz_vars(dat, var)
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
    end subroutine

    module subroutine parse_reaction(instring, reverse, eqr, eqp, err)
      character(len=*), intent(in) :: instring
      logical, intent(out) :: reverse
      character(len=s_str_len), allocatable, intent(out) :: eqr(:), eqp(:)
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    module subroutine interp2xsdata(dat, var, err)
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    module subroutine compute_gibbs_energy(dat, var, err)
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
      character(:), allocatable, intent(out) :: err
    end subroutine

    module subroutine interp2particlexsdata(dat, var, err)
      use photochem_const, only: smaller_real
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
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

  !> Initialize atmosphere-dependent state from a legacy atmosphere text file.
  !!
  !! [[setup_static]] must have completed successfully before this routine is
  !! called. This routine reads the raw atmosphere, resolves grid-dependent
  !! settings, constructs the grid and gravity, interpolates atmospheric and
  !! optical properties, computes temperature-dependent data, and validates the
  !! tropopause. It does not allocate or prepare EvoAtmosphere RHS work arrays.
  subroutine setup_atmosphere_from_file(atmosphere_txt, dat, var, err)

    character(len=*), intent(in) :: atmosphere_txt
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err

    call read_atmosphere_file(atmosphere_txt, dat, var, err)
    if (allocated(err)) return

    call resolve_atmosphere_settings(dat, var, err)
    if (allocated(err)) return

    call after_read_setup(dat, var, err)
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
