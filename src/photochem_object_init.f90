submodule(photochem_object) photochem_object_init
  implicit none
  
contains
  
  subroutine Photochem_init(self, data_dir, mechanism_file, settings_file, flux_file, atmosphere_txt, err)
    use photochem_input, only: setup
    
    class(Photochem), intent(inout) :: self
    character(len=*), intent(in) :: data_dir
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    character(len=err_len), intent(out) :: err
    err = ""
    self%var%data_dir = data_dir
    call setup(mechanism_file, settings_file, flux_file, atmosphere_txt, &
               self%dat, self%var, err)
    if (len_trim(err) /= 0) return 
    
  end subroutine
  
end submodule