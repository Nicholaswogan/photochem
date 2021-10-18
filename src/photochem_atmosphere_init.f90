submodule(photochem_atmosphere) photochem_atmosphere_init
  implicit none
  
contains
  
  subroutine Atmosphere_init(self, data_dir, mechanism_file, settings_file, flux_file, atmosphere_txt, err)
    use photochem_input, only: setup
    
    class(Atmosphere), intent(inout) :: self
    character(len=*), intent(in) :: data_dir
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    character(len=err_len), intent(out) :: err
    err = ""
    
    if (allocated(self%dat)) then
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
    if (len_trim(err) /= 0) return 
    
    call self%wrk%init(self%dat%nsp, self%dat%np, self%dat%nq, &
                       self%var%nz, self%dat%nrT, self%dat%kj, &
                       self%dat%nw, self%var%trop_ind)
    
  end subroutine
  
end submodule