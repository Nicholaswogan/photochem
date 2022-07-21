submodule(photochem_evoatmosphere) photochem_evoatmosphere_init
  implicit none
  
  ! Contains the Constructor for the EvoAtmosphere derived type.
  
contains
  
  module subroutine EvoAtmosphere_init(self, data_dir, mechanism_file, settings_file, flux_file, atmosphere_txt, err)
    use iso_c_binding, only : c_associated
    use photochem_input, only: setup
    use photochem_types, only: PhotoSettings
    
    class(EvoAtmosphere), intent(inout) :: self
    character(len=*), intent(in) :: data_dir
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    character(:), allocatable, intent(out) :: err

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
    call setup(mechanism_file, s, flux_file, atmosphere_txt, .false., &
               self%dat, self%var, err)
    if (allocated(err)) return 
    
    call self%wrk%init(self%dat%nsp, self%dat%np, self%dat%nq, &
                       self%var%nz, self%dat%nrT, self%dat%kj, &
                       self%dat%nw, self%var%trop_ind)
    call self%prep_atmosphere(self%var%usol_init, err)
    if (allocated(err)) return            
  end subroutine
  
end submodule