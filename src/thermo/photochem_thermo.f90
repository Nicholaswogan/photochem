module photochem_thermo
  use photochem_const, only: real_kind, err_len
  use photochem_types, only: ThermodynamicData
  implicit none
    

type :: EquilibriumSolver
  
  ! data
  integer :: n_species
  integer :: n_atoms
  type(ThermodynamicData), allocatable :: thermo_data(:)
  real(real_kind), allocatable :: atoms(:,:)
  real(real_kind), allocatable :: atoms_copy(:,:)
  real(real_kind), allocatable :: lambda_guess(:)
  
  ! work
  real(real_kind), allocatable :: lambda(:)
  real(real_kind), allocatable :: gibbs(:)
  real(real_kind), allocatable :: g_tilde(:)
  real(real_kind), allocatable :: alpha(:)
  real(real_kind), allocatable :: Z_vec(:)
  
  real(real_kind), allocatable :: B_dist(:) ! maybe get rid of

  ! lapack work arrays
  real(real_kind), allocatable :: S_vals(:)
  integer :: rank
  real(real_kind), allocatable :: work(:)
  integer :: lwork
  integer, allocatable :: iwork(:)
  integer :: liwork
  ! end lapack work arrays
  
  ! hybrd1 work arrays
  real(real_kind) :: tol = 1.d-8
  integer :: lwa
  real(real_kind), allocatable :: wa(:)
  real(real_kind), allocatable :: fvec(:)
  ! end hybrd1
  
contains
  
  procedure :: init => EquilibriumSolver_init
  procedure :: equilibrate
    
end type

contains
  
  subroutine dumby()
    
    
  end subroutine
  
  subroutine EquilibriumSolver_init(eq, species_composition, thermo_data, err)
    class(EquilibriumSolver), intent(inout) :: eq
  
    integer, intent(in) :: species_composition(:,:)
    type(ThermodynamicData), intent(in) :: thermo_data(:)
    character(len=err_len), intent(out) :: err
    
    integer :: i, j, info
    
    err = ""
  
    eq%n_species = size(species_composition,2)
    eq%n_atoms = size(species_composition,1)
    
    if (eq%n_species /= size(thermo_data)) then
      err = "Dimensions issue in EquilibriumSolver_init"
      return
    endif
  
    !!! allocate memory
    allocate(eq%atoms(eq%n_species,eq%n_atoms))
    allocate(eq%atoms_copy(eq%n_species,eq%n_atoms))
    allocate(eq%thermo_data(eq%n_species))
    allocate(eq%gibbs(eq%n_species))
    allocate(eq%g_tilde(eq%n_species))
    allocate(eq%Z_vec(eq%n_species))
    allocate(eq%B_dist(eq%n_species))
    allocate(eq%lambda(eq%n_atoms))
    allocate(eq%lambda_guess(eq%n_atoms))
    allocate(eq%alpha(eq%n_atoms))
  
    ! lapack work arrays
    allocate(eq%S_vals(min(eq%n_species,eq%n_atoms)))
    ! use dgelsd to get length of work arrays
    eq%lwork = -1
    allocate(eq%iwork(1))
    allocate(eq%work(1))
    call dgelsd(eq%n_species, eq%n_atoms, 1, eq%atoms, eq%n_species, &
                eq%B_dist, eq%n_species, eq%S_vals, -1.d0, eq%rank, eq%work, &
                eq%lwork, eq%iwork, info)
    eq%lwork = eq%work(1)
    eq%liwork = eq%iwork(1)
    deallocate(eq%work, eq%iwork)
    allocate(eq%work(eq%lwork))
    allocate(eq%iwork(eq%liwork))
    ! end lapack work arrays
    ! minpack
    eq%lwa = (eq%n_atoms*(3*eq%n_atoms+13))/2 + 2
    allocate(eq%wa(eq%lwa))
    allocate(eq%fvec(eq%n_atoms))
    ! end minpack
    !!! end allocating memory
  
    !!! data for all calculations
    do i = 1,eq%n_species
      do j = 1,eq%n_atoms
        eq%atoms(i,j) = species_composition(j,i)
      enddo
    enddo
    eq%thermo_data = thermo_data
    !!! end data for all calculations
  
    call InitialGuess(eq, err)
    if (len_trim(err) /= 0) return
  
  end subroutine
  
  subroutine InitialGuess(eq, err)
    class(EquilibriumSolver), intent(inout) :: eq
    character(len=err_len), intent(out) :: err
    
    integer :: info
    
    err = ""
    
    ! Camberos+ 2001 recommends uniform distriubtion
    eq%B_dist = 1.d0/eq%n_species

    eq%atoms_copy = eq%atoms
    call dgelsd(eq%n_species, eq%n_atoms, 1, eq%atoms, eq%n_species, &
                eq%B_dist, eq%n_species, eq%S_vals, -1.d0, eq%rank, eq%work, &
                eq%lwork, eq%iwork, info)
    if (info /= 0) then
      err  = "Failed to find initial guess for EquilibriumSolver"
      return
    endif
    ! restore atoms to original state
    eq%atoms = eq%atoms_copy
    ! solution. Initial guess for lambda
    eq%lambda_guess = eq%B_dist(1:eq%n_atoms)
    
  end subroutine
  
  subroutine Equilibrate(eq, T, P, mole_fractions, err)
    use iso_c_binding, only: c_ptr, c_loc
    use photochem_eqns, only: gibbs_energy_eval
    use cminpack2fort, only: hybrd1
    use photochem_const, only: Rgas
    
    class(EquilibriumSolver), target, intent(inout) :: eq
    
    real(real_kind), intent(in) :: T, P
    real(real_kind), intent(inout) :: mole_fractions(:)
    character(len=err_len), intent(out) :: err
    
    integer :: i, info
    real(real_kind) :: error, sum_mix
    logical :: found
    type(c_ptr) :: ptr
    type(EquilibriumSolver), pointer :: eq_ptr
    
    err = ""
    
    if (size(mole_fractions) /= eq%n_species) then
      err = "mole_fractions has the wrong dimensions in Equilibrate"
      return
    endif
    
    sum_mix = sum(mole_fractions)
    if (sum_mix < 1.d0-1.d-5 .or. sum_mix > 1.d0+1.d-5) then
      err = "mole_fractions in Equilibrate does not sum to 1"
      return
    endif
    
    do i = 1,eq%n_species
      call gibbs_energy_eval(eq%thermo_data(i), T, found, eq%gibbs(i))
      if (.not. found) then
        err = "The temperature is not within the ranges "// &
                "given for the thermodynamic data."
        return
      endif
    enddo
    
    eq%g_tilde = eq%gibbs/(Rgas*T) + log(P/1.d0) ! Here we assume P0 is 1 bar
    
    do i = 1,eq%n_atoms
      ! here we assume total moles
      eq%alpha(i) = sum(eq%atoms(:,i)*mole_fractions(:))*1.d7
    enddo
    print*,eq%alpha
    
    ! now we are ready for optimization
    eq_ptr => eq
    ptr = c_loc(eq_ptr)
    eq%lambda = eq%lambda_guess
    
    call hybrd1(fcn_gibbs, ptr, eq%n_atoms, eq%lambda, eq%fvec, eq%tol, info, eq%wa, eq%lwa)
    
    error = sum(sqrt(eq%fvec**2.d0))
    print*,info,error,eq%fvec
    if (info /= 1 .or. error > eq%tol) then  
      err = "Non-linear solved failed in Equilibrate."
      return 
    endif
    
    ! compute solution
    do i = 1,eq%n_species
      mole_fractions(i) = exp(-eq%g_tilde(i) + sum(eq%lambda(:)*eq%atoms(i,:)))
    enddo
    
  end subroutine

  function fcn_gibbs(ptr, n_atoms, lambda, fvec, iflag) result(res) bind(c)
    use iso_c_binding, only: c_ptr, c_f_pointer
    
    type(c_ptr) :: ptr ! void pointer. Be careful!
    integer, value :: n_atoms, iflag ! n == n_atoms
    real(real_kind), intent(in) :: lambda(n_atoms) 
    real(real_kind), intent(out) :: fvec(n_atoms) ! fvec == residual
    integer :: res
    
    integer :: i, j
    real(real_kind) :: val1, val2
    type(EquilibriumSolver), pointer :: eq

    call c_f_pointer(ptr, eq)

    ! compute Z_vec
    do i = 1,eq%n_species
      eq%Z_vec(i) = exp(-eq%g_tilde(i) + sum(lambda(:)*eq%atoms(i,:)))
    enddo

    ! first equation
    fvec(1) = log10(sum(eq%Z_vec)) - log10(1.d0)
    ! other equations
    val2 = sum(eq%atoms(:,eq%n_atoms)*eq%Z_vec(:))
    do j = 1,eq%n_atoms-1
      val1 = sum(eq%atoms(:,j)*eq%Z_vec(:))
      fvec(j+1) = log10(1 + eq%alpha(eq%n_atoms)*val1) - log10(1 + eq%alpha(j)*val2)
    enddo
    
    ! stop
    res = 0
  end function


end module