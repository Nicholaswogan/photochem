
module photochem_data
  implicit none  
  public ! public but protected
  integer, private, parameter :: real_kind = kind(1.0d0)
  ! Data that doesn't change during model runs.
  ! No writing to it!

  ! molecules
  integer, protected :: nsp
  integer, protected :: natoms
  character(len=8), allocatable, protected :: atoms_names(:) 
  character(len=8), allocatable, protected :: species_names(:)
  integer, allocatable, protected :: species_composition(:,:)
  integer, allocatable, protected :: lowerboundcond(:)
  real(real_kind), allocatable, protected :: lower_vdep(:)
  real(real_kind), allocatable, protected :: lower_flux(:)
  real(real_kind), allocatable, protected :: lower_distributed_height(:)
  real(real_kind), allocatable, protected :: lower_fixed_mr(:)
  integer, allocatable, protected :: upperboundcond(:)
  real(real_kind), allocatable, protected :: upper_veff(:)
  real(real_kind), allocatable, protected :: upper_flux(:)
  real(real_kind), allocatable, protected :: thermo_data(:,:,:)
  real(real_kind), allocatable, protected :: thermo_temps(:,:)
  
  ! reaction
  integer, protected :: nrF
  integer, protected :: nrR
  integer, protected :: nrT
  integer, protected :: max_num_reactants
  integer, protected :: max_num_products
  character(len=8), allocatable, protected :: reactants_names(:,:) ! not really needed.
  character(len=8), allocatable, protected :: products_names(:,:)
  integer, allocatable, protected :: reactants_sp_inds(:,:)
  integer, allocatable, protected :: products_sp_inds(:,:)
  integer, allocatable, protected :: nreactants(:)
  integer, allocatable, protected :: nproducts(:)
  integer, allocatable, protected :: reverse_info(:) ! all for calculating rates
  character(len=15), allocatable, protected :: rxtypes(:)
  real(real_kind), allocatable, protected :: rateparams(:,:)
  
  integer, allocatable, protected :: nump(:) ! length nsp. number of 
  integer, allocatable, protected :: numl(:)
  integer, allocatable, protected :: iprod(:,:)
  integer, allocatable, protected :: iloss(:,:)

contains
  
  subroutine get_data(yaml_file,)
    use photochem_types, only: PhotoMechanism
    use photochem_io, only: get_photomech
    use photochem_vars, only: ...  
    
    character(len=*), intent(in) :: 
    
    type(PhotoMechanism) :: photomech
    
    ! if allocated, deallocate
    if (allocated(atoms_names)) then
      deallocate(atoms_names)
      deallocate(species_names)
      .
      .
      .
      ! deallocate work variables

    endif
    
    call get_photomech(yaml_file,)
    
    ! allocate
    
    ! allocate work variables
  end subroutine

  
end module