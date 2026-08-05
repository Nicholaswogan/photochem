submodule(photochem_evoatmosphere) photochem_evoatmosphere_init
  implicit none
  
  ! Contains the Constructor for the EvoAtmosphere derived type.
  
contains
  
  module function create_EvoAtmosphere(mechanism_file, settings_file, flux_file, atmosphere_txt, data_dir, err) result(self)
    use photochem_input, only: setup_static
    use photochem_types, only: PhotoSettings
    
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err
    type(EvoAtmosphere) :: self

    type(PhotoSettings) :: s

    s = PhotoSettings(settings_file, err)
    if (allocated(err)) return
    
    if (allocated(self%dat)) then
      deallocate(self%dat)
      deallocate(self%var)
      deallocate(self%wrk)
    endif
    
    allocate(self%dat)
    allocate(self%var)
    allocate(self%wrk)
    
    self%var%data_dir = data_dir
    call setup_static(mechanism_file, s, flux_file, self%dat, self%var, err)
    if (allocated(err)) return
    
    call self%wrk%init(self%dat%nsp, self%dat%np, self%dat%nq, &
                       self%var%nz, self%dat%nrT, self%dat%kj, &
                       self%dat%nw)

    call self%initialize_from_atmosphere_file(atmosphere_txt, err)
    if (allocated(err)) return

  end function

  module subroutine initialize_from_atmosphere_file(self, atmosphere_txt, err)
    use photochem_input, only: setup_atmosphere_from_file

    class(EvoAtmosphere), intent(inout) :: self
    character(len=*), intent(in) :: atmosphere_txt
    character(:), allocatable, intent(out) :: err

    self%atmosphere_initialized = .false.

    call setup_atmosphere_from_file(atmosphere_txt, self%dat, self%var, err)
    if (allocated(err)) return

    call self%prep_atmosphere(self%var%usol_init, err)
    if (allocated(err)) return

    self%atmosphere_initialized = .true.

  end subroutine

  module subroutine require_atmosphere_initialized(self, operation, err)
    class(EvoAtmosphere), intent(in) :: self
    character(len=*), intent(in) :: operation
    character(:), allocatable, intent(out) :: err

    if (.not. self%atmosphere_initialized) then
      err = 'EvoAtmosphere atmosphere is not initialized. Initialize the atmosphere before calling "'// &
            trim(operation)//'".'
    endif

  end subroutine
  
end submodule
