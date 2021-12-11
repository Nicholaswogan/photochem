
module photochem_input
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  use photochem_types, only : PhotochemData, PhotochemVars
  use photochem_const, only: real_kind, str_len, err_len, s_str_len
  implicit none
  private 

  public :: setup, interp2xsdata, compute_gibbs_energy
  
  interface
    module subroutine after_read_setup(photodata, photovars, err)
      use photochem_eqns, only: vertical_grid, gravity
      type(PhotochemData), intent(inout) :: photodata
      type(PhotochemVars), intent(inout) :: photovars
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine read_all_files(mechanism_file, settings_file, flux_file, atmosphere_txt, &
                              photodata, photovars, err)
      character(len=*), intent(in) :: mechanism_file
      character(len=*), intent(in) :: settings_file
      character(len=*), intent(in) :: flux_file
      character(len=*), intent(in) :: atmosphere_txt
      type(PhotochemData), intent(inout) :: photodata
      type(PhotochemVars), intent(inout) :: photovars
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine interp2xsdata(dat, var, err)
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
      character(len=err_len), intent(out) :: err
    end subroutine
    
    module subroutine compute_gibbs_energy(dat, var, err)
      type(PhotochemData), intent(in) :: dat
      type(PhotochemVars), intent(inout) :: var
      character(len=err_len), intent(out) :: err
    end subroutine
    
  end interface
    
contains
  
  subroutine setup(mechanism_file, settings_file, flux_file, atmosphere_txt, &
                   photodata, photovars, err)
    
    character(len=*), intent(in) :: mechanism_file
    character(len=*), intent(in) :: settings_file
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    type(PhotochemData), intent(inout) :: photodata
    type(PhotochemVars), intent(inout) :: photovars
    character(len=err_len), intent(out) :: err
    
    err = ""
    
    call read_all_files(mechanism_file, settings_file, flux_file, atmosphere_txt, &
                        photodata, photovars, err)
    if (len_trim(err) /= 0) return     
                 
    call after_read_setup(photodata, photovars, err)
    if (len_trim(err) /= 0) return
    
  end subroutine
  
end module



