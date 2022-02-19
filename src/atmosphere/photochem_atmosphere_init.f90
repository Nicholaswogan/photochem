submodule(photochem_atmosphere) photochem_atmosphere_init
  implicit none
  
  ! Contains the Constructor for the Atmosphere derived type.
  ! This constructor reads all data files, setting up atmosphere
  ! for photochemical calculations.
  
contains
  
  module subroutine Atmosphere_init(self, data_dir, mechanism_file, settings_file, flux_file, atmosphere_txt, err)
    use iso_c_binding, only : c_associated
    use photochem_input, only: setup
    
    class(Atmosphere), intent(inout) :: self
    character(len=*), intent(in) :: data_dir
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    character(:), allocatable, intent(out) :: err
    
    if (allocated(self%dat)) then
      if (c_associated(self%wrk%cvode_mem)) then
        call self%destroy_stepper(err)
        if (allocated(err)) return 
      endif
      deallocate(self%dat)
      deallocate(self%var)
      deallocate(self%wrk)
    endif
    
    allocate(self%dat)
    allocate(self%var)
    allocate(self%wrk)
    
    self%var%data_dir = data_dir
    call setup(mechanism_file, settings_file, flux_file, atmosphere_txt, &
               self%dat, self%var, err)
    if (allocated(err)) return 
    
    call self%wrk%init(self%dat%nsp, self%dat%np, self%dat%nq, &
                       self%var%nz, self%dat%nrT, self%dat%kj, &
                       self%dat%nw, self%var%trop_ind)
                       
    call self%prep_atmosphere(self%var%usol_init, err)
    if (allocated(err)) return 
  end subroutine
  
end submodule