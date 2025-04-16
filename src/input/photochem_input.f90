
module photochem_input
  use fortran_yaml_c_types, only : type_node, type_dictionary, type_list, type_error, &
                         type_list_item, type_scalar, type_key_value_pair
  use photochem_types, only : PhotochemData, PhotochemVars, PhotoSettings
  use photochem_const, only: dp, str_len, s_str_len
  implicit none
  private 

  public :: setup, interp2xsdata, compute_gibbs_energy, interp2particlexsdata
  
  type, extends(type_list) :: type_list_tmp
  ! temporary list for accessing all reactions and
  ! species in a row.
  contains
    final :: list_destroy
  end type
  
  interface
    module subroutine after_read_setup(dat, var, err)
      use photochem_eqns, only: vertical_grid, gravity
      type(PhotochemData), intent(inout) :: dat
      type(PhotochemVars), intent(inout) :: var
      character(:), allocatable, intent(out) :: err
    end subroutine
    
    module subroutine read_all_files(mechanism_file, s, flux_file, atmosphere_txt, dat, var, err)
      character(len=*), intent(in) :: mechanism_file
      type(PhotoSettings), intent(in) :: s
      character(len=*), intent(in) :: flux_file
      character(len=*), intent(in) :: atmosphere_txt
      type(PhotochemData), intent(inout) :: dat
      type(PhotochemVars), intent(inout) :: var
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
  
  subroutine setup(mechanism_file, s, flux_file, atmosphere_txt, dat, var, err)
    
    character(len=*), intent(in) :: mechanism_file
    type(PhotoSettings), intent(in) :: s
    character(len=*), intent(in) :: flux_file
    character(len=*), intent(in) :: atmosphere_txt
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    call read_all_files(mechanism_file, s, flux_file, atmosphere_txt, dat, var, err)
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



