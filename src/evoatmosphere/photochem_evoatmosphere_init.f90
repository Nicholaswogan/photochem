submodule(photochem_evoatmosphere) photochem_evoatmosphere_init
  implicit none
  
  ! Contains the Constructor for the EvoAtmosphere derived type.
  
contains
  
  module function create_EvoAtmosphere(mechanism_file, settings_file, flux_file, atmosphere_txt, data_dir, err) result(self)
    use iso_c_binding, only : c_associated
    use photochem_input, only: setup
    use photochem_types, only: PhotoSettings
    use photochem_enum, only: PressureBC
    use clima_types, only: ClimaSettings
    
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    character(len=*), intent(in) :: data_dir
    character(:), allocatable, intent(out) :: err
    type(EvoAtmosphere) :: self

    type(PhotoSettings) :: s
    type(ClimaSettings) :: cs

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
    call setup(mechanism_file, s, flux_file, atmosphere_txt, self%dat, self%var, err)
    if (allocated(err)) return 
    
    call self%wrk%init(self%dat%nsp, self%dat%np, self%dat%nq, &
                       self%var%nz, self%dat%nrT, self%dat%kj, &
                       self%dat%nw)
    
    if (s%evolve_climate) then
      self%evolve_climate = .true.

      ! Make sure there are no pressure boundary conditions.
      if (any(self%var%lowerboundcond == PressureBC)) then
        err = 'Fixed pressure boundary conditions are not allowed for class "EvoAtmosphere" '// &
              'when evolve_climate is true.'
        return
      endif

      ! allocate
      if (allocated(self%rad)) deallocate(self%rad)
      allocate(self%rad)

      ! create the settings object
      cs = ClimaSettings(settings_file, err)
      if (allocated(err)) return

      ! create radtran
      block
      integer :: num_zenith_angles
      num_zenith_angles = 1
      self%rad = Radtran(self%dat%species_names(self%dat%ng_1:self%dat%nq), &
                         self%dat%species_names(1:self%dat%np), &
                         cs, flux_file, &
                         num_zenith_angles, s%surface_albedo, self%var%nz, data_dir, err)
      if (allocated(err)) return
      end block

    else
      self%evolve_climate = .false.
    endif
    ! The initial guess for T_surf
    self%T_surf = self%var%temperature(1)
    
    call self%prep_atmosphere(self%var%usol_init, err)
    if (allocated(err)) return   

  end function
  
end submodule