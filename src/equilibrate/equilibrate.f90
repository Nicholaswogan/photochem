module equilibrate
  use equilibrate_const, only: dp, atom_str_len, reac_str_len, s_str_len
  use equilibrate_cea, only: CEAData
  implicit none
  private

  public :: ChemEquiAnalysis, dp, s_str_len

  !> Chemical-equilibrium solver and the state from its most recent solve.
  !>
  !> Name and mass arrays are established during construction. Composition and
  !> thermodynamic outputs are updated by [[ChemEquiAnalysis:solve]] and
  !> [[ChemEquiAnalysis:solve_metallicity]]; check the returned convergence flag
  !> before using them. Array ordering follows the corresponding name array.
  type :: ChemEquiAnalysis
    character(atom_str_len), allocatable :: atoms_names(:) !! Atom names, shape `(na)`.
    character(reac_str_len), allocatable :: species_names(:) !! All species names, shape `(ns)`.
    real(dp), allocatable :: species_mass(:) !! Species molar masses, shape `(ns)` (g/mol).
    character(reac_str_len), allocatable :: gas_names(:) !! Gas species names, shape `(ng)`.
    real(dp), allocatable :: gas_mass(:) !! Gas molar masses, shape `(ng)` (g/mol).
    character(reac_str_len), allocatable :: condensate_names(:) !! Condensed species names, shape `(nc)`.
    real(dp), allocatable :: condensate_mass(:) !! Condensate molar masses, shape `(nc)` (g/mol).

    !> Describes the composition of each species (i.e., how many atoms)
    real(dp), allocatable, private :: species_composition(:,:)

    !> Reference solar atom mole fractions in `atoms_names` order, shape `(na)`.
    real(dp), allocatable :: molfracs_atoms_sun(:)

    ! Results from the most recent equilibrium solve.

    real(dp), allocatable :: molfracs_atoms(:) !! Atom mole fractions across all phases, shape `(na)`.
    real(dp), allocatable :: molfracs_species(:) !! Species mole fractions across all phases, shape `(ns)`.
    real(dp), allocatable :: massfracs_species(:) !! Species mass fractions across all phases, shape `(ns)`.

    real(dp), allocatable :: molfracs_atoms_gas(:) !! Gas-phase atom mole fractions, shape `(na)`.
    real(dp), allocatable :: molfracs_species_gas(:) !! Gas-phase species mole fractions, shape `(ng)`.

    real(dp), allocatable :: molfracs_atoms_condensate(:) !! Condensed-phase atom mole fractions, shape `(na)`.
    real(dp), allocatable :: molfracs_species_condensate(:) !! Condensed species mole fractions, shape `(nc)`.

    real(dp) :: nabla_ad !! Adiabatic logarithmic temperature gradient, `dln(T)/dln(P)`.
    real(dp) :: gamma2 !! Second adiabatic exponent, `1/(1 - nabla_ad)`.
    real(dp) :: rho !! Equilibrium gas mass density (g/cm^3).
    real(dp) :: c_pe !! Equilibrium specific heat at constant pressure (erg/(g K)).

    real(dp) :: mubar !! Mean molecular weight of the gas phase (g/mol).

    logical :: verbose = .false. !! Whether the equilibrium solver prints iteration details.
    real(dp) :: mass_tol = 1.0e-6_dp !! Dimensionless mass-balance convergence tolerance.
    logical :: use_prev_guess = .false. !! Whether to initialize from the previous converged solution.

    !> Driver class
    type(CEAData), allocatable :: dat

  contains
    procedure :: solve
    procedure :: solve_metallicity
  end type
  interface ChemEquiAnalysis
    module procedure :: create_ChemEquiAnalysis
  end interface

contains
    
  !> Initializes a chemical-equilibrium solver from thermodynamic data.
  !>
  !> Legacy `.inp` files require both `atoms` and `species`. For YAML files:
  !>
  !> - If both lists are supplied, exactly those atoms and species are used.
  !> - If only `atoms` is supplied, every compatible species is used.
  !> - If only `species` is supplied, every atom in those species is used.
  !> - If neither is supplied, every atom and species in the file is used.
  function create_ChemEquiAnalysis(thermopath, atoms, species, err) result(cea)
    character(*), intent(in) :: thermopath !! Path to a `.yaml` or legacy `.inp` thermodynamic file.
    character(*), optional, intent(in) :: atoms(:) !! Nonempty atom-name selection.
    character(*), optional, intent(in) :: species(:) !! Nonempty species-name selection.
    character(:), allocatable, intent(out) :: err !! Error description; unallocated on success.
    type(ChemEquiAnalysis) :: cea

    integer :: i, j, jj, k

    if (present(atoms)) then
      if (size(atoms) < 1) then
        err = '"atoms" must not be empty.'
        return
      endif
      if (len(atoms(1)) > atom_str_len) then
        err = 'Entries in "atoms" exceed the supported string length.'
        return
      endif
    endif
    if (present(species)) then
      if (size(species) < 1) then
        err = '"species" must not be empty.'
        return
      endif 
      if (len(species(1)) > reac_str_len) then
        err = 'Entries in "species" exceed the supported string length.'
        return
      endif
    endif

    if (allocated(cea%dat)) deallocate(cea%dat)
    allocate(cea%dat)

    call cea%dat%set_data(thermopath, atoms, species)
    if (cea%dat%error) then
      err = trim(cea%dat%err_msg)
      return
    endif

    cea%atoms_names = cea%dat%names_atoms(1:size(cea%dat%names_atoms)-1) ! skip E
    cea%species_names = cea%dat%names_reactants_orig
    allocate(cea%species_mass(size(cea%dat%mol_weight)))
    do i = 1,size(cea%species_mass)
      cea%species_mass(i) = cea%dat%mol_weight(cea%dat%id_reactants(i,1))
    enddo

    allocate(cea%condensate_names(cea%dat%N_cond))
    allocate(cea%gas_names(cea%dat%N_gas))
    allocate(cea%condensate_mass(cea%dat%N_cond))
    allocate(cea%gas_mass(cea%dat%N_gas))
    j = 1
    jj = 1
    do i = 1,size(cea%species_names)
      if (cea%dat%reac_condensed(i)) then
        cea%condensate_names(j) = cea%species_names(i)
        cea%condensate_mass(j) = cea%species_mass(i)
        j = j + 1
      else
        cea%gas_names(jj) = cea%species_names(i)
        cea%gas_mass(jj) = cea%species_mass(i)
        jj = jj + 1
      endif
    enddo

    if (allocated(cea%molfracs_species)) then
      deallocate(cea%species_composition)
      deallocate(cea%molfracs_atoms)
      deallocate(cea%molfracs_atoms_gas)
      deallocate(cea%molfracs_atoms_condensate)
      deallocate(cea%molfracs_species)
      deallocate(cea%molfracs_species_gas)
      deallocate(cea%molfracs_species_condensate)
      deallocate(cea%massfracs_species)
    endif
    allocate(cea%species_composition(size(cea%atoms_names),size(cea%species_names)))
    allocate(cea%molfracs_atoms(size(cea%atoms_names)))
    allocate(cea%molfracs_atoms_gas(size(cea%atoms_names)))
    allocate(cea%molfracs_atoms_condensate(size(cea%atoms_names)))
    allocate(cea%molfracs_species(size(cea%species_names)))
    allocate(cea%molfracs_species_gas(size(cea%gas_names)))
    allocate(cea%molfracs_species_condensate(size(cea%condensate_names)))
    allocate(cea%massfracs_species(size(cea%species_names)))

    ! Get composition of each species
    cea%species_composition = 0.0_dp
    do j = 1,size(cea%species_names)
      jj = cea%dat%id_reactants(j,1)
      do k = 1,size(cea%dat%reac_atoms_id,1)
        if (cea%dat%reac_atoms_id(k,jj) > 0) then
          i = findloc(cea%dat%id_atoms, cea%dat%reac_atoms_id(k,jj), 1)
          if (i == 0 .or. i > size(cea%atoms_names) + 1) then
            err = 'Indexing error during initialization.'
            return
          endif
          if (i == size(cea%atoms_names) + 1) then
            ! Electron so we skip
            cycle
          endif
          cea%species_composition(i,j) = cea%dat%reac_stoich(k,jj)
        endif
      enddo
    enddo

    block
      use equilibrate_const, only: atoms_sun, molfracs_atoms_sun
      integer :: ind
      ! Default composition of Sun
      allocate(cea%molfracs_atoms_sun(size(cea%atoms_names)))
      cea%molfracs_atoms_sun = 1.0e-50_dp
      do i = 1,size(cea%atoms_names)
        ind = findloc(atoms_sun, cea%atoms_names(i), 1)
        if (ind /= 0) then
          cea%molfracs_atoms_sun(i) = molfracs_atoms_sun(ind)
        endif
      enddo
      cea%molfracs_atoms_sun = cea%molfracs_atoms_sun/sum(cea%molfracs_atoms_sun)
    endblock

  end function

  !> Computes chemical equilibrium from atom or species mole fractions.
  !>
  !> Supply exactly one composition array. The input is normalized internally.
  !> On a completed solve, the returned flag reports convergence and the result
  !> fields of `self` contain the latest equilibrium state.
  function solve(self, P, T, molfracs_atoms, molfracs_species, err) result(converged)
    class(ChemEquiAnalysis), intent(inout) :: self
    real(dp), intent(in) :: P !! Pressure (dyn/cm^2).
    real(dp), intent(in) :: T !! Temperature (K).
    !> Nonnegative atom mole fractions in `atoms_names` order, shape `(na)`.
    real(dp), optional, intent(in) :: molfracs_atoms(:)
    !> Nonnegative species mole fractions in `species_names` order, shape `(ns)`.
    real(dp), optional, intent(in) :: molfracs_species(:)
    character(:), allocatable, intent(out) :: err !! Error description; unallocated on a completed solve.
    logical :: converged !! Whether the equilibrium iteration met its tolerances.

    real(dp), allocatable :: molfracs_atoms_(:)
    real(dp) :: P_bars
    integer :: i, j, jj

    converged = .false.
    P_bars = P/1.0e6_dp

    if (present(molfracs_atoms) .and. present(molfracs_species)) then
      err = 'Both "molfracs_atoms" and "molfracs_species" are inputs, but only one is allowed.'
      return
    endif
    if (.not.present(molfracs_atoms) .and. .not.present(molfracs_species)) then
      err = 'Neither "molfracs_atoms" nor "molfracs_species" was provided; exactly one is required.'
      return
    endif

    ! Input is an array of atom mole fractions
    if (present(molfracs_atoms)) then
      if (size(molfracs_atoms) /= size(self%atoms_names)) then
        err = 'Input "molfracs_atoms" has the wrong size'
        return
      endif
      if (any(molfracs_atoms < 0)) then
        err = 'Input "molfracs_atoms" cannot contain negative values.'
        return
      endif
      molfracs_atoms_ = molfracs_atoms/max(sum(molfracs_atoms),tiny(1.0_dp))
    endif

    ! Input is an array of species mole fractions
    if (present(molfracs_species)) then
      if (size(molfracs_species) /= size(self%species_names)) then
        err = 'Input "molfracs_species" has the wrong size'
        return
      endif
      if (any(molfracs_species < 0)) then
        err = 'Input "molfracs_species" cannot contain negative values.'
        return
      endif

      allocate(molfracs_atoms_(size(self%atoms_names)))
      molfracs_atoms_ = 0.0_dp
      do j = 1,size(self%species_names)
        molfracs_atoms_(:) = molfracs_atoms_(:) + self%species_composition(:,j)*molfracs_species(j)
      enddo
      molfracs_atoms_ = molfracs_atoms_/max(sum(molfracs_atoms_),tiny(1.0_dp))
    endif

    ! Compute chemical equilibrium
    call self%dat%solve(mode='q', verbo='  ', verbose2=self%verbose, N_atoms_in=size(self%atoms_names), &
                           N_reactants_in=size(self%species_names), molfracs_atoms=molfracs_atoms_, &
                           mass_tol=self%mass_tol, &
                           molfracs_reactants=self%molfracs_species, &
                           massfracs_reactants=self%massfracs_species, &
                           temp=T, press=P_bars, &
                           nabla_ad=self%nabla_ad, gamma2=self%gamma2, MMW=self%mubar, rho=self%rho, c_pe=self%c_pe, &
                           use_prev_guess=self%use_prev_guess)
    if (self%dat%error) then
      err = trim(self%dat%err_msg)
      return
    endif
    converged = self%dat%converged

    ! Compute mole fractions of gases and condensates in the solution
    j = 1
    jj = 1
    do i = 1,size(self%species_names)
      if (self%dat%reac_condensed(i)) then
        self%molfracs_species_condensate(j) = self%molfracs_species(i)
        j = j + 1
      else
        self%molfracs_species_gas(jj) = self%molfracs_species(i)
        jj = jj + 1
      endif
    enddo
    if (size(self%condensate_names) > 0) then
      self%molfracs_species_condensate = self%molfracs_species_condensate/max(sum(self%molfracs_species_condensate),tiny(1.0_dp))
    endif
    if (size(self%gas_names) > 0) then
      self%molfracs_species_gas = self%molfracs_species_gas/max(sum(self%molfracs_species_gas),tiny(1.0_dp))
    endif

    ! Compute the mole fractions of the atoms in the solution
    self%molfracs_atoms = 0.0_dp
    self%molfracs_atoms_gas = 0.0_dp
    self%molfracs_atoms_condensate = 0.0_dp
    do j = 1,size(self%species_names)
      self%molfracs_atoms(:) = self%molfracs_atoms(:) + self%species_composition(:,j)*self%molfracs_species(j)
      if (self%dat%reac_condensed(j)) then
        self%molfracs_atoms_condensate(:) = self%molfracs_atoms_condensate(:) &
                                            + self%species_composition(:,j)*self%molfracs_species(j)
      else
        self%molfracs_atoms_gas(:) = self%molfracs_atoms_gas(:) &
                                     + self%species_composition(:,j)*self%molfracs_species(j)
      endif
    enddo
    self%molfracs_atoms = self%molfracs_atoms/max(sum(self%molfracs_atoms),tiny(1.0_dp))
    self%molfracs_atoms_condensate = self%molfracs_atoms_condensate/max(sum(self%molfracs_atoms_condensate),tiny(1.0_dp))
    self%molfracs_atoms_gas = self%molfracs_atoms_gas/max(sum(self%molfracs_atoms_gas),tiny(1.0_dp))

    ! Compute mean molecular weight in gas
    self%mubar = sum(self%molfracs_species_gas*self%gas_mass)

  end function

  !> Computes chemical equilibrium from metallicity and an optional C/O ratio.
  !>
  !> Heavy-element abundances in `molfracs_atoms_sun` are scaled relative to H
  !> and He. If `CtoO` is present, carbon and oxygen are redistributed while
  !> preserving their combined abundance. Result fields are updated as in
  !> [[ChemEquiAnalysis:solve]].
  function solve_metallicity(self, P, T, metallicity, CtoO, err) result(converged)
    class(ChemEquiAnalysis), intent(inout) :: self
    real(dp), intent(in) :: P !! Pressure (dyn/cm^2).
    real(dp), intent(in) :: T !! Temperature (K).
    real(dp), intent(in) :: metallicity !! Positive metallicity relative to the reference solar composition.
    !> Positive C/O ratio relative to the reference solar value; `1` preserves solar C/O.
    real(dp), optional, intent(in) :: CtoO
    character(:), allocatable, intent(out) :: err !! Error description; unallocated on a completed solve.
    logical :: converged !! Whether the equilibrium iteration met its tolerances.

    real(dp), allocatable :: molfracs_atoms(:)
    integer :: i, indC, indO
    real(dp) :: x, a

    if (metallicity <= 0) then
      err = '"metallicity" must be greater than 0.'
      return
    endif

    ! Compute atoms based on the input metallicity.
    molfracs_atoms = self%molfracs_atoms_sun
    do i = 1,size(molfracs_atoms)
      if (self%atoms_names(i) /= 'H' .and. self%atoms_names(i) /= 'He') then
        molfracs_atoms(i) = self%molfracs_atoms_sun(i)*metallicity
      endif
    enddo
    molfracs_atoms = molfracs_atoms/sum(molfracs_atoms)

    ! Adjust C/O ratio, if specified
    if (present(CtoO)) then
      if (CtoO <= 0) then
        err = '"CtoO" must be greater than 0.'
        return
      endif

      indC = findloc(self%atoms_names, 'C', 1)
      if (indC == 0) then
        err = 'C must be an atom if CtoO is specified.'
        return
      endif

      indO = findloc(self%atoms_names, 'O', 1)
      if (indO == 0) then
        err = 'O must be an atom if CtoO is specified.'
        return
      endif

      x = CtoO*(molfracs_atoms(indC)/molfracs_atoms(indO))
      a = (x*molfracs_atoms(indO) - molfracs_atoms(indC))/(1 + x)
      molfracs_atoms(indC) = molfracs_atoms(indC) + a
      molfracs_atoms(indO) = molfracs_atoms(indO) - a
    endif

    converged = self%solve(P, T, molfracs_atoms=molfracs_atoms, err=err)
    if (allocated(err)) return

  end function

end module
