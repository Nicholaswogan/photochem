
program main
  use Photochem, only: Atmosphere, err_len, real_kind
  use photochem_types, only: ThermodynamicData
  use photochem_thermo
  implicit none
  character(len=err_len) :: err
  type(Atmosphere) :: pc
  
  type(EquilibriumSolver) :: eq
  integer :: n_species, i
  integer, allocatable :: species_composition(:,:)
  real(real_kind), allocatable :: mole_fractions(:)
  type(ThermodynamicData), allocatable :: thermo_data(:)
  
  logical :: success
  
  err = ""
  call pc%init("../Photochem/data", &
               "../Photochem/data/reaction_mechanisms/zahnle_earth.yaml", &
               "../templates/ModernEarth/settings_ModernEarth.yaml", &
               "../templates/ModernEarth/Sun_now.txt", &
               "../templates/ModernEarth/atmosphere_ModernEarth.txt", &
               err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  print*,pc%dat%atoms_names
  n_species = pc%dat%nll+1
  allocate(species_composition(pc%dat%natoms,n_species))
  allocate(thermo_data(n_species))
  allocate(mole_fractions(n_species))
  
  species_composition(:,1:n_species-1) = pc%dat%species_composition(:, pc%dat%ng_1:pc%dat%nq)
  species_composition(:,n_species) = pc%dat%species_composition(:, pc%dat%nsp)
  
  thermo_data(1:n_species-1) = pc%dat%thermo_data(1:pc%dat%nll)
  thermo_data(n_species) = pc%dat%thermo_data(pc%dat%ng)

  ! call eq%init(pc%dat%species_composition(:,pc%dat%ng_1:pc%dat%nsp), pc%dat%thermo_data, err)
  call eq%init(species_composition, thermo_data, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  mole_fractions(1:n_species-1) = pc%var%usol_init(pc%dat%ng_1:pc%dat%nq,1)
  mole_fractions(n_species) = 1.d0 - sum(pc%var%usol_init(pc%dat%ng_1:pc%dat%nq,1))
  
  mole_fractions = 1.d0/n_species
  call eq%equilibrate(300.d0, 1.d0, mole_fractions, err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop 1
  endif
  
  do i = 1,n_species-1
    print*,pc%dat%species_names(pc%dat%npq+i),mole_fractions(i)
  enddo
  print*,pc%dat%species_names(pc%dat%nsp),mole_fractions(n_species)
  
  

end program
