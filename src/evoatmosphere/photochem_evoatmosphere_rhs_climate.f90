submodule(photochem_evoatmosphere:photochem_evoatmosphere_rhs) photochem_evoatmosphere_rhs_climate
  use futils, only: brent_class
  implicit none

  type :: EvoClimateData

    ! inputs
    real(dp), pointer :: usol_den(:,:)
    real(dp), pointer :: T_trop
    
    ! calculated from inputs
    real(dp), allocatable :: pdensities(:,:)
    real(dp), allocatable :: radii(:,:)

    ! outputs
    real(dp), pointer :: T_surf
    real(dp), pointer :: T(:)
    real(dp), pointer :: z_trop

    ! work
    real(dp), allocatable :: densities(:,:)
    real(dp), allocatable :: P(:)

  end type

  type, extends(brent_class) :: brent_custom
    type(EvoAtmosphere), pointer :: self
    type(EvoClimateData), pointer :: d
  end type

contains

  module subroutine equilibrium_climate(self, usol_den, molecules_per_particle, T_trop, T_surf_guess, &
                                        T_surf, T, z_trop, err)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), target, intent(in) :: usol_den(:,:)
    real(dp), intent(in) :: molecules_per_particle(:,:)
    real(dp), target, intent(in) :: T_trop, T_surf_guess
    real(dp), target, intent(out) :: T_surf, T(:)
    real(dp), target, intent(out) :: z_trop
    character(:), allocatable, intent(out) :: err

    integer :: i
    type(MinpackHybrd1Vars) :: mv

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(EvoClimateData), target :: d

    dat => self%dat
    var => self%var

    ! compute particle densities (particles/cm3)
    ! and the radii
    allocate(d%pdensities(var%nz,dat%np))
    allocate(d%radii(var%nz,dat%np))
    do i = 1,dat%np
      d%pdensities(:,i) = usol_den(i,:)*(1.0_dp/molecules_per_particle(i,:))
      d%radii(:,i) = var%particle_radius(i,:)
    enddo

    ! allocate work space
    allocate(d%densities(var%nz,dat%nll))
    allocate(d%P(var%nz))
    ! associate
    d%usol_den => usol_den
    d%T_trop => T_trop
    d%T_surf => T_surf
    d%T => T
    d%z_trop => z_trop

    mv = MinpackHybrd1Vars(1, tol=1.0e-8_dp)
    mv%x(1) = log10(T_surf_guess)
    call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
    if (mv%info /= 1) then; block
      type(brent_custom) :: brent
      real(dp) :: T_surf_small, T_surf_big
      ! hybrd1 failed to find the root.
      ! Fallback on brent's method. 
      call brent%set_function(fcn_brent)
      brent%self => self
      brent%d => d
      
      T_surf_small = T_surf_guess*0.99_dp
      T_surf_big = T_surf_guess*1.01_dp
      call brent%find_zero(log10(T_surf_small), log10(T_surf_big), mv%tol, mv%x(1), mv%fvec(1), mv%info)
      if (mv%info /= 0) then
        err = 'Non-linear solve for temperature structure failed.'
        return
      endif

    endblock; endif

    ! output T_surf
    T_surf = 10.0_dp**mv%x(1)
    ! output T(:) and z_trop from calling make_profile one more time
    ! with the solution T_surf.
    call make_profile_discrete(self, usol_den, T_surf, T_trop, &
                               T, d%P, d%densities, &
                               z_trop, err)
    if (allocated(err)) return

  contains
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_

      character(:), allocatable :: err_

      ! self and d are passed in by scope
      d%T_surf = 10.0_dp**x_(1)
      if (associated(self%albedo_fcn)) then
        self%rad%surface_albedo = self%albedo_fcn(d%T_surf)
      endif
      fvec_(1) = equilibrium_climate_equation(self, d, err_)
      if (allocated(err_)) then
        iflag_ = -1
        return
      endif
      
    end subroutine
  end subroutine

  function fcn_brent(me, x) result(f)
    class(brent_class),intent(inout) :: me
    real(dp),intent(in) :: x
    real(dp) :: f

    character(:), allocatable :: err

    select type (me)
    class is (brent_custom)
        me%d%T_surf = 10.0_dp**x
        f = equilibrium_climate_equation(me%self, me%d, err)
    end select

  end function

  function equilibrium_climate_equation(self, d, err) result(F_in)
    type(EvoAtmosphere), intent(inout) :: self
    type(EvoClimateData), intent(inout) :: d
    character(:), allocatable, intent(out) :: err
    real(dp) :: F_in

    ! make a T profile
    call make_profile_discrete(self, d%usol_den, d%T_surf, d%T_trop, &
                               d%T, d%P, d%densities, &
                               d%z_trop, err)
    if (allocated(err)) return

    ! do raditative transfer
    call self%rad%radiate(d%T_surf, d%T, d%P, d%densities, self%var%dz, d%pdensities, d%radii, err=err)
    if (allocated(err)) return

    ! compute the energy going into the atmosphere
    F_in = self%rad%wrk_ir%fdn_n(self%var%nz+1) - self%rad%wrk_ir%fup_n(self%var%nz+1) + &
           self%rad%wrk_sol%fdn_n(self%var%nz+1) - self%rad%wrk_sol%fup_n(self%var%nz+1)

  end function

  subroutine make_profile_discrete(self, usol_den, T_surf, T_trop, &
                                   T, P, densities, &
                                   z_trop, err)
    use photochem_const, only: k_boltz, T_crit_H2O
    use clima_eqns_water, only: sat_pressure_H2O
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_den(:,:)
    real(dp), intent(in) :: T_surf, T_trop
    real(dp), intent(out) :: T(:), P(:), densities(:,:)
    real(dp), intent(out) :: z_trop
    character(:), allocatable, intent(out) :: err

    real(dp) :: n_dry(self%var%nz), mu_dry(self%var%nz)
    real(dp) :: f_H2O_trop, P_H2O, n_H2O(self%var%nz)

    integer :: i, j, k, trop_ind

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! check inputs
    if (T_surf < T_trop) then
      err = "Surface temperature is lower then the tropopause temperture"
      return
    endif

    if (T_surf > T_crit_H2O) then
      err = "Surface temperature is greater than the critical temperature for H2O"
      return
    endif

    ! compute dry atmospheric properties
    do j = 1,var%nz
      n_dry(j) = 0.0_dp
      do i = dat%ng_1,dat%nq
        if (i /= dat%LH2O) then
          n_dry(j) = n_dry(j) + usol_den(i,j)
        endif
      enddo
      mu_dry(j) = 0.0_dp
      do i = dat%ng_1,dat%nq
        if (i /= dat%LH2O) then
          mu_dry(j) = mu_dry(j) + (usol_den(i,j)/n_dry(j))*dat%species_mass(i)
        endif
      enddo
    enddo

    call integrate_moist(self, usol_den, n_dry, mu_dry, T_surf, T_trop, T, trop_ind, z_trop, err)
    if (allocated(err)) return

    ! For raditave transfer we use an estimate of the H2O density profile (NOT the real one)
    ! H2O is fixed to some relative humidity in the troposphere
    ! H2O is constant above the troposphere
    do j = 1,trop_ind
      P_H2O = var%relative_humidity*sat_pressure_H2O(T(j))
      n_H2O(j) = P_H2O/(k_boltz*T(j))
    enddo
    ! Guess mixing ratio of H2O at the tropopause
    f_H2O_trop = n_H2O(trop_ind)/(n_dry(trop_ind) + n_H2O(trop_ind))
    ! assume H2O mixing ratio is constant above the troposphere
    do j = trop_ind+1,var%nz
      n_H2O(j) = n_dry(j)*(f_H2O_trop/(1.0_dp - f_H2O_trop))
    enddo

    do i = 1,dat%nll
      do j = 1,var%nz
        k = i + dat%npq
        densities(j,i) = usol_den(k,j)
      enddo
    enddo
    k = dat%LH2O - dat%npq
    densities(:,k) = n_H2O

    do j = 1,var%nz
      ! pressure in bars
      P(j) = sum(densities(j,:))*k_boltz*T(j)/1.0e6_dp
    enddo

  end subroutine
  
  subroutine integrate_moist(self, usol_den, n_dry, mu_dry, T_surf, T_trop, temperature, trop_ind, z_trop, err)
    use futils, only: interp
    use photochem_const, only: k_boltz
    use photochem_eqns, only: heat_capacity_eval
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_den(:,:)
    real(dp), intent(in) :: n_dry(:), mu_dry(:)
    real(dp), intent(in) :: T_surf, T_trop
    real(dp), intent(out) :: temperature(:)
    integer, intent(out) :: trop_ind
    real(dp), intent(out) :: z_trop
    character(:), allocatable, intent(out) :: err

    real(dp) :: cp_dry
    real(dp) :: T_edge, dT_dz
    real(dp), allocatable :: z_inv(:), T_inv(:)
    real(dp) :: z_trop_(1)

    integer :: j, ierr
    logical :: trop_hit

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! integrate from surface to center of first grid cell
    cp_dry = dry_heat_capacity_of_layer(dat, T_surf, usol_den(:,1), n_dry(1), err)
    if (allocated(err)) return
    dT_dz = dT_dz_moist(T_surf, var%relative_humidity, n_dry(1), mu_dry(1), cp_dry, var%grav(1))    
    temperature(1) = T_surf + dT_dz*(0.5_dp*var%dz(1))

    trop_hit = .false. 
    do j = 1,var%nz-1

      ! Integrate to top edge of layer j
      cp_dry = dry_heat_capacity_of_layer(dat, temperature(j), usol_den(:,j), n_dry(j), err)
      if (allocated(err)) return
      dT_dz = dT_dz_moist(temperature(j), var%relative_humidity, n_dry(j), mu_dry(j), cp_dry, var%grav(j))    
      T_edge = temperature(j) + dT_dz*(0.5_dp*var%dz(j))
      
      ! Integrate from top edge of j to middle of j+1
      cp_dry = dry_heat_capacity_of_layer(dat, T_edge, usol_den(:,j+1), n_dry(j+1), err)
      if (allocated(err)) return
      dT_dz = dT_dz_moist(T_edge, var%relative_humidity, n_dry(j+1), mu_dry(j+1), cp_dry, var%grav(j+1))
      temperature(j+1) = T_edge + dT_dz*(0.5_dp*var%dz(j+1))

      if (temperature(j+1) <= T_trop) then
        trop_hit = .true.
        trop_ind = j+1
        exit
      endif
    enddo

    if (.not. trop_hit) then
      err = "Tropopause was never reached"
      return
    endif

    ! now we find the exact altitude of the tropopause
    allocate(z_inv(trop_ind), T_inv(trop_ind))

    z_inv = var%z(1:trop_ind)
    T_inv = temperature(1:trop_ind)
    z_inv = z_inv(trop_ind:1:-1)
    T_inv = T_inv(trop_ind:1:-1)

    call interp(1, size(z_inv), [T_trop], T_inv, z_inv, z_trop_, ierr=ierr)
    if (ierr /= 0) then
      err = "interpolation failed while making a T profile"
      return
    endif

    temperature(trop_ind:) = T_trop
    z_trop = z_trop_(1)

  end subroutine

  function dry_heat_capacity_of_layer(dat, T, usol_den_layer, n_dry, err) result(cp_dry)
    use photochem_eqns, only: heat_capacity_eval
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: T, usol_den_layer(:)
    real(dp), intent(in) :: n_dry
    character(:), allocatable, intent(out) :: err
    real(dp) :: cp_dry

    real(dp) :: cp, mubar_dry
    logical :: found
    integer :: i

    cp_dry = 0.0_dp
    mubar_dry = 0.0_dp
    do i = dat%ng_1,dat%nq
      if (i /= dat%LH2O) then

        ! J/(mol*K)
        call heat_capacity_eval(dat%thermo_data(i-dat%npq), T, found, cp)
        if (.not. found) then
          err = 'Could not compute heat capacity for '//trim(dat%species_names(i))
          return
        endif

        cp_dry = cp_dry + (usol_den_layer(i)/n_dry)*cp ! J/(mol*K)
        mubar_dry = mubar_dry + (usol_den_layer(i)/n_dry)*dat%species_mass(i) ! g/mol
      endif
    enddo
    ! J/(mol*K) * (kg/mol) = J/(kg*K)
    cp_dry = cp_dry*(1.0_dp/(mubar_dry*1.0e-3_dp)) ! J/(kg*K)
    ! convert J/(kg*K) to erg/(g*K)
    cp_dry = cp_dry*1.0e4_dp

  end function

  pure function dP_dry_dT(T, P_dry, mu_dry, cp_dry, P_H2O, mu_H2O, cp_H2O, L_H2O) result(dPdT)
    use clima_eqns_water, only: Rgas
    real(dp), intent(in) :: T
    real(dp), intent(in) :: P_dry, mu_dry, cp_dry
    real(dp), intent(in) :: P_H2O, mu_H2O, cp_H2O, L_H2O
    real(dp) :: dPdT

    dPdT = (P_dry*mu_dry*cp_dry)/(Rgas*T)*&
           (P_dry + (cp_H2O/cp_dry + (L_H2O*mu_H2O/(Rgas*T) - 1.0_dp)&
                  *(L_H2O/(cp_dry*T)))*((mu_H2O*P_H2O)/(mu_dry)))/ &
           (P_dry + (L_H2O*mu_dry)/(Rgas*T)*(mu_H2O*P_H2O)/(mu_dry))
  end function
  
  pure function dP_H2O_dT(T, P_H2O, mu_H2O, L_H2O) result(dPdT)
    use clima_eqns_water, only: Rgas
    real(dp), intent(in) :: T
    real(dp), intent(in) :: P_H2O, mu_H2O, L_H2O
    real(dp) :: dPdT
    dPdT = (L_H2O*mu_H2O*P_H2O)/(Rgas*T**2.0_dp)
  end function

  function dT_dz_moist(T, rh, n_dry, mu_dry, cp_dry, grav) result(dTdz)
    use clima_eqns_water, only: latent_heat_H2O, sat_pressure_H2O, mu_H2O
    use photochem_const, only: k_boltz, N_avo
    real(dp), intent(in) :: T, rh
    real(dp), intent(in) :: n_dry, mu_dry, cp_dry, grav
    real(dp) :: dTdz

    real(dp) :: dPdT_dry, dPdT_H2O
    real(dp) :: P_H2O, n_H2O
    real(dp) :: cp_H2O
    real(dp) :: P_dry
    real(dp) :: n_tot, mubar

    !! latent heat of H2O vaporization/condensation (erg/g)
    real(dp) :: L_H2O

    L_H2O = latent_heat_H2O(T)
    P_H2O = rh*sat_pressure_H2O(T)
    n_H2O = P_H2O/(k_boltz*T)
    cp_H2O = heat_capacity_H2O(T) ! J/(mol*K)
    cp_H2O = cp_H2O*(1.0_dp/(mu_H2O*1.0e-3_dp)) ! J/(kg*K)
    cp_H2O = cp_H2O*1.0e4_dp ! convert to erg/(g*K)

    ! dry pressure
    P_dry = n_dry*k_boltz*T

    ! combine dry and moist atmosphere
    n_tot = n_dry + n_H2O
    mubar = (n_dry/n_tot)*mu_dry + (n_H2O/n_tot)*mu_H2O

    dPdT_dry = dP_dry_dT(T, P_dry, mu_dry, cp_dry, P_H2O, mu_H2O, cp_H2O, L_H2O)
    dPdT_H2O = dP_H2O_dT(T, P_H2O, mu_H2O, L_H2O) 
    
    dTdz = (-grav*mubar*n_tot/N_avo) * (dPdT_dry + dPdT_H2O)**(-1.0_dp)

  end function

  pure function heat_capacity_H2O(T) result(cp)
    use photochem_eqns, only: heat_capacity_shomate
    real(dp), intent(in) :: T !! K
    real(dp) :: cp !! J/(mol*K)
    
    real(dp), parameter :: coeffs_1(7) = &
          [30.092, 6.832514, 6.793435, -2.53448, 0.082139, -250.881, 223.3967]
    real(dp), parameter :: coeffs_2(7) = &
          [41.96426, 8.622053, -1.49978, 0.098119, -11.15764, -272.1797, 219.7809]
    
    if (T > 0.0_dp .and. T <= 1700.0_dp) then
      cp = heat_capacity_shomate(coeffs_1, T)
    elseif (T > 1700.0_dp .and. T <= 6000.0_dp) then
      cp = heat_capacity_shomate(coeffs_2, T)
    endif

  end function

end submodule