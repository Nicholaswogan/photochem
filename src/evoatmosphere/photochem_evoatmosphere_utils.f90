
submodule(photochem_evoatmosphere) photochem_evoatmosphere_utils
  implicit none

contains

module subroutine out2atmosphere_txt(self, filename, number_of_decimals, overwrite, clip, err)
    use photochem_common, only: out2atmosphere_txt_base
    class(EvoAtmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: number_of_decimals
    logical, intent(in) :: overwrite
    logical, intent(in) :: clip
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: rhs(self%var%neqs)  
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    ! update wrk variables
    call self%right_hand_side_chem(wrk%usol, rhs, err)
    if (allocated(err)) return

    call out2atmosphere_txt_base(dat, var, &
                                 wrk%pressure, wrk%density, wrk%densities, wrk%molecules_per_particle, &
                                 filename, number_of_decimals, overwrite, clip, err)
    if (allocated(err)) return

  end subroutine

  module subroutine gas_fluxes(self, surf_fluxes, top_fluxes, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(out) :: surf_fluxes(:)
    real(dp), intent(out) :: top_fluxes(:)
    character(:), allocatable, intent(out) :: err
  
    real(dp) :: rhs(self%var%neqs)  
    real(dp) :: diffusive_production
    real(dp) :: chemical_production
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
  
    integer :: i
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (size(surf_fluxes) /= dat%nq .or. size(top_fluxes) /= dat%nq) then
      err = "Input fluxes to gas_fluxes has the wrong dimensions"
      return
    endif
  
    call self%right_hand_side_chem(wrk%usol, rhs, err)
    if (allocated(err)) return
    
    ! surface flux is molecules required to sustain the lower boundary
    ! chemical production + diffusion production = total change in lower cell    
    do i = 1,dat%nq
      diffusive_production = (wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                            + wrk%DD(i,1)*wrk%usol(i,1) + wrk%ADD(i,1)*wrk%usol(i,1)) &
                              *var%dz(1)
      chemical_production = rhs(i)*var%dz(1)
      surf_fluxes(i) = -(diffusive_production + chemical_production)
    enddo
    
    ! fluxes going into or out of the top of the atmosphere.
    do i = 1,dat%nq
      diffusive_production = &
         (wrk%DD(i,var%nz)*wrk%usol(i,var%nz) + wrk%ADD(i,var%nz)*wrk%usol(i,var%nz) &
        + wrk%DL(i,var%nz)*wrk%usol(i,var%nz-1) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz-1)) &
          *var%dz(var%nz)
    
      chemical_production = rhs(i + (var%nz-1)*dat%nq)*var%dz(var%nz)
      top_fluxes(i) = diffusive_production + chemical_production
    enddo
    
  end subroutine

  module subroutine set_lower_bc(self, species, bc_type, vdep, den, press, flux, height, err)
    use photochem_enum, only: MosesBC, VelocityBC, DensityBC, PressureBC, FluxBC, VelocityDistributedFluxBC
    class(EvoAtmosphere), intent(inout) :: self
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: bc_type
    real(dp), optional, intent(in) :: vdep
    real(dp), optional, intent(in) :: den
    real(dp), optional, intent(in) :: press
    real(dp), optional, intent(in) :: flux
    real(dp), optional, intent(in) :: height
    character(:), allocatable, intent(out) :: err
    
    integer :: ind
    
    ind = findloc(self%dat%species_names(1:self%dat%nq), trim(species), 1)
    if (ind == 0) then
      err = "Can not change boundary condition of '"//trim(species)// &
            "' because it is not in the list of species"
      return
    endif
    if (self%dat%there_are_particles) then
      if (ind <= self%dat%np .and. bc_type /= 'vdep') then
        err = 'Particles must have a deposition velocity lower boundary condition.'
        return
      endif
    endif
    
    if (bc_type == 'vdep') then
      if (.not. present(vdep)) then
        err = "To change boundary condition to deposition"// &
              " velocity must supply the 'vdep' argument"
        return
      endif
      self%var%lowerboundcond(ind) = VelocityBC
      self%var%lower_vdep(ind) = vdep
      
    elseif (bc_type == 'den') then
      if (.not. present(den)) then
        err = "To change boundary condition to fixed density"// &
              " you must supply the 'den' argument"
        return
      endif
      self%var%lowerboundcond(ind) = DensityBC
      self%var%lower_fix_den(ind) = den

    elseif (bc_type == 'press') then
      if (.not. present(press)) then
        err = "To change boundary condition to fixed pressure"// &
              " you must supply the 'press' argument"
        return
      endif
      self%var%lowerboundcond(ind) = PressureBC
      self%var%lower_fix_press(ind) = press
      
    elseif (bc_type == 'flux') then
      if (.not. present(flux)) then
        err = "To change boundary condition to a surface flux"// &
              " must supply the 'flux' argument"
        return
      endif
      self%var%lowerboundcond(ind) = FluxBC
      self%var%lower_flux(ind) = flux

    elseif (bc_type == 'vdep + dist flux') then
      if (.not.present(vdep) .or. .not.present(flux) .or. .not.present(height)) then
        err = "To change boundary condition to deposition velocity with"// &
              " a distributed flux, must supply the 'vdep', 'flux', and 'height' arguments"
        return
      endif
      self%var%lowerboundcond(ind) = VelocityDistributedFluxBC
      self%var%lower_vdep(ind) = vdep
      self%var%lower_flux(ind) = flux
      self%var%lower_dist_height(ind) = height
      
    elseif (bc_type == 'Moses') then
      self%var%lowerboundcond(ind) = MosesBC
    else
      err = "Boundary condition type '"//trim(bc_type)//"' is not a valid"// &
            " boundary condition type"
      return
    endif
    
  end subroutine

  module subroutine set_upper_bc(self, species, bc_type, veff, flux, err)
    use photochem_enum, only: VelocityBC, FluxBC
    use photochem_enum, only: DiffusionLimHydrogenEscape
    class(EvoAtmosphere), intent(inout) :: self
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: bc_type
    real(dp), optional, intent(in) :: veff
    real(dp), optional, intent(in) :: flux
    character(:), allocatable, intent(out) :: err
    
    integer :: ind(1)
    
    ind = findloc(self%dat%species_names(1:self%dat%nq), trim(species))
    if (ind(1) == 0) then
      err = "Can not change boundary condition of '"//trim(species)// &
            "' because it is not in the list of species"
      return
    endif
    
    if (self%dat%H_escape_type == DiffusionLimHydrogenEscape) then
      if (species == "H2") then
        err = "You can not change the boundary condition for H2 because"// &
              " diffusion limited H2 escape is on."
        return
      endif
      if (species == "H") then
        err = "You can not change the boundary condition for H because"// &
              " diffusion limited H escape is on."
        return
      endif
    endif
    
    if (bc_type == 'veff') then
      if (.not. present(veff)) then
        err = "To change boundary condition to effusion"// &
              " velocity must supply the 'veff' argument"
        return
      endif
      self%var%upperboundcond(ind(1)) = VelocityBC
      self%var%upper_veff(ind(1)) = veff

    elseif (bc_type == 'flux') then
      if (.not. present(flux)) then
        err = "To change boundary condition to a flux"// &
              " must supply the 'flux' argument"
        return
      endif
      self%var%upperboundcond(ind(1)) = FluxBC
      self%var%upper_flux(ind(1)) = flux

    else
      err = "Boundary condition type '"//trim(bc_type)//"' is not a valid"// &
            " upper boundary condition type"
      return
    endif
  
  end subroutine

  module subroutine set_rate_fcn(self, species, fcn, err)
    use photochem_types, only: time_dependent_rate_fcn
    class(EvoAtmosphere), target, intent(inout) :: self
    character(*), intent(in) :: species
    procedure(time_dependent_rate_fcn), pointer :: fcn
    character(:), allocatable, intent(inout) :: err
    
    integer :: ind

    ind = findloc(self%dat%species_names(1:self%dat%nq), species, 1)
    if (ind == 0) then
      err = 'Species "'//species//'" is not in the list of species, '// &
            'or is a short-lived species.'
      return
    endif

    self%var%rate_fcns(ind)%fcn => fcn

  end subroutine

  module subroutine set_temperature(self, temperature, trop_alt, err)
    use photochem_input, only: interp2xsdata, compute_gibbs_energy
    
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: temperature(:)
    real(dp), optional, intent(in) :: trop_alt
    character(:), allocatable, intent(out) :: err
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemVars) :: var_save
    
    dat => self%dat
    var => self%var
    
    if (size(temperature) /= var%nz) then
      err = "temperature has the wrong input dimension"
      return
    endif
    ! save in case there is an issue
    var_save = var
    
    var%temperature = temperature
    
    ! xsections and gibbs energy needs updating
    call interp2xsdata(dat, var, err)
    if (allocated(err)) then
      var = var_save
      return
    endif
    if (dat%reverse) then
      call compute_gibbs_energy(dat, var, err)
      if (allocated(err)) then
        var = var_save
        return
      endif
    endif
    
    ! If gas rainout is enabled and trop_alt is present, then we need to
    ! change trop_ind, reallocate some stuff
    ! in wrk, then we will re-prep the atmosphere
    if (dat%gas_rainout .and. present(trop_alt)) then
      if (trop_alt < var%bottom_atmos .or. trop_alt > var%top_atmos) then
        var = var_save
        err = "trop_alt is above or bellow the atmosphere!"
        return
      endif
      
      var%trop_alt = trop_alt
      var%trop_ind = max(minloc(abs(var%z - var%trop_alt), 1) - 1, 1)

      if (var%trop_ind < 3) then
        var = var_save
        err = 'Tropopause is too low.'
        return
      elseif (var%trop_ind > var%nz-2) then
        var = var_save
        err = 'Tropopause is too high.'
        return
      endif

    endif

    ! Fill wrk with new values
    call self%prep_atmosphere(self%wrk%usol, err)
    if (allocated(err)) return
    
  end subroutine

  module subroutine set_press_temp_edd(self, P, T, edd, trop_p, hydro_pressure, err)
    use futils, only: interp, brent_class
    use ieee_arithmetic, only: ieee_is_finite
    use photochem_const, only: small_real, k_boltz, N_avo
    use photochem_enum, only: DensityBC, PressureBC
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: P(:)
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: edd(:)
    real(dp), optional, intent(in) :: trop_p
    logical, optional, intent(in) :: hydro_pressure
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: log10P_in(:), T_in(:), log10edd_in(:)
    real(dp), allocatable :: P_wrk(:), log10P_wrk(:), T_new(:), log10edd_new(:), edd_save(:)
    real(dp), allocatable :: usol_base(:,:), density_base(:), mubar_base(:)
    real(dp) :: xzero, Psurf_initial, Psurf_final
    real(dp) :: log10P_previous, temperature_previous
    logical :: hydro_pressure_
    integer :: ierr, i, j, residual_layer
    character(32) :: layer_string

    real(dp), parameter :: log10e = log10(exp(1.0_dp))
    real(dp), parameter :: root_tol = 1.0e-10_dp
    integer, parameter :: max_bracket_expansions = 60

    type(brent_class) :: root_solver

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var
    
    if (size(P) /= size(T)) then
      err = '"P" and "T" not the same size'
      return
    endif

    if (size(P) /= size(edd)) then
      err = '"P" and "edd" not the same size'
      return
    endif

    if (size(P) < 2) then
      err = 'size(P) must be >= 2'
      return
    endif

    if (.not. all(ieee_is_finite(P)) .or. .not. all(ieee_is_finite(T)) .or. &
        .not. all(ieee_is_finite(edd))) then
      err = 'All elements of "P", "T" and "edd" must be finite'
      return
    endif

    if (any(P <= 0.0_dp) .or. any(T <= 0.0_dp) .or. any(edd <= 0.0_dp)) then
      err = 'All elements of "P", "T" and "edd" must be positive'
      return
    endif

    if (.not. all(P(2:) < P(:size(P)-1))) then
      err = '"P" must be strictly decreasing'
      return
    endif

    if (present(trop_p)) then
      if (.not. ieee_is_finite(trop_p) .or. trop_p <= 0.0_dp) then
        err = '"trop_p" must be finite and positive'
        return
      endif
    endif

    if (dat%gas_rainout .and. .not.present(trop_p)) then
      err = '"trop_p" is a required input.'
      return
    endif

    ! optional arguments
    if (present(hydro_pressure)) then
      hydro_pressure_ = hydro_pressure
    else
      hydro_pressure_ = .true. ! default is True
    endif

    allocate(P_wrk(var%nz), log10P_wrk(var%nz), T_new(var%nz))
    allocate(log10edd_new(var%nz), edd_save(var%nz))
    allocate(usol_base(dat%nq,var%nz), density_base(var%nz), mubar_base(var%nz))

    ! Copy and clip the current state in the same way as prep_atm_evo_gas.
    usol_base = self%wrk%usol
    where (usol_base >= 0.0_dp)
      usol_base = max(usol_base, small_real)
    elsewhere
      usol_base = min(usol_base, -small_real)
    endwhere

    ! Above the bottom layer, density and mean molecular weight do not
    ! depend on temperature. The bottom layer is handled separately because
    ! fixed-pressure boundary conditions depend on temperature.
    do j = 2,var%nz
      density_base(j) = sum(usol_base(dat%ng_1:,j))
      if (.not. ieee_is_finite(density_base(j)) .or. density_base(j) <= 0.0_dp) then
        write(layer_string,'(i0)') j
        err = 'The gas density is not finite and positive in layer '//trim(layer_string)
        return
      endif
      mubar_base(j) = sum(dat%species_mass(dat%ng_1:dat%nq)* &
                          usol_base(dat%ng_1:dat%nq,j))/density_base(j)
      if (.not. ieee_is_finite(mubar_base(j)) .or. mubar_base(j) <= 0.0_dp) then
        write(layer_string,'(i0)') j
        err = 'The mean molecular weight is not finite and positive in layer '//trim(layer_string)
        return
      endif
    enddo

    call bottom_column_state(var%temperature(1), density_base(1), mubar_base(1), Psurf_initial)
    if (.not. ieee_is_finite(density_base(1)) .or. density_base(1) <= 0.0_dp .or. &
        .not. ieee_is_finite(mubar_base(1)) .or. mubar_base(1) <= 0.0_dp .or. &
        .not. ieee_is_finite(Psurf_initial) .or. Psurf_initial <= 0.0_dp) then
      err = 'Could not compute a finite, positive bottom-layer state'
      return
    endif

    ! copy over inputs, and covert to log10 space
    allocate(log10P_in(size(P)),T_in(size(P)),log10edd_in(size(P)))
    log10P_in = log10(P)
    T_in = T
    log10edd_in = log10(edd)

    ! if the P-T-edd profile does not extend to the surface,
    ! then we log-linearly extrapolate to surface
    if (Psurf_initial > P(1)) then; block
      real(dp) :: slope, intercept, P_surf, T_surf, edd_surf

      ! log10 surface pressure in dynes/cm^2
      P_surf = log10(Psurf_initial)

      slope = (T_in(2) - T_in(1))/(log10P_in(2) - log10P_in(1))
      intercept = T_in(1) - slope*log10P_in(1)
      T_surf = slope*P_surf + intercept
      if (.not. ieee_is_finite(T_surf) .or. T_surf <= 0.0_dp) then
        err = 'Extrapolating the input P-T profile to the surface produced a non-positive temperature'
        return
      endif

      slope = (log10edd_in(2) - log10edd_in(1))/(log10P_in(2) - log10P_in(1))
      intercept = log10edd_in(1) - slope*log10P_in(1)
      edd_surf = slope*P_surf + intercept
      if (.not. ieee_is_finite(edd_surf)) then
        err = 'Extrapolating the input P-Kzz profile to the surface failed'
        return
      endif

      log10P_in = [P_surf, log10P_in]
      T_in = [T_surf, T_in]
      log10edd_in = [edd_surf, log10edd_in]
      
    endblock; endif

    ! Flip order for interpolation purposes
    log10P_in = log10P_in(size(log10P_in):1:-1)
    T_in = T_in(size(log10P_in):1:-1)
    log10edd_in = log10edd_in(size(log10P_in):1:-1)

    call root_solver%set_function(pressure_residual)

    ! Solve one scalar pressure equation per layer. In hydrostatic mode the
    ! recurrence is triangular, so each solved layer provides the lower
    ! pressure and temperature boundary for the next layer.
    do i = 1,var%nz
      call solve_layer_pressure(i, xzero, err)
      if (allocated(err)) return

      log10P_wrk(i) = xzero
      P_wrk(i) = 10.0_dp**xzero
      T_new(i) = profile_value(xzero, log10P_in, T_in)
      log10edd_new(i) = profile_value(xzero, log10P_in, log10edd_in)

      if (.not. ieee_is_finite(P_wrk(i)) .or. P_wrk(i) <= 0.0_dp) then
        write(layer_string,'(i0)') i
        err = 'The pressure solve produced a non-finite or non-positive pressure in layer '// &
              trim(layer_string)
        return
      endif
      if (hydro_pressure_ .and. i > 1) then
        if (P_wrk(i) >= P_wrk(i-1)) then
          write(layer_string,'(i0)') i
          err = 'The hydrostatic pressure does not decrease upward at layer '//trim(layer_string)
          return
        endif
      endif

      if (i == 1) then
        call bottom_column_state(T_new(1), density_base(1), mubar_base(1), Psurf_final)
        if (hydro_pressure_ .and. P_wrk(1) >= Psurf_final) then
          err = 'The solved bottom-layer pressure is not below the surface pressure'
          return
        endif
      endif

      log10P_previous = xzero
      temperature_previous = T_new(i)
    enddo

    ! Commit point: everything above only reads self and works with local
    ! arrays. Updating Kzz and calling set_temperature below are the only
    ! operations in this routine that mutate self.
    if (present(trop_p)) then; block
      real(dp) :: trop_alt(1)

      call interp([log10(trop_p)], log10P_wrk(var%nz:1:-1), &
                  var%z(var%nz:1:-1), trop_alt, ierr=ierr)
      if (ierr /= 0) then
        err = 'Subroutine interp returned an error.'
        return
      endif

      edd_save = var%edd
      var%edd = 10.0_dp**log10edd_new
      call self%set_temperature(T_new, trop_alt(1), err)
      if (allocated(err)) var%edd = edd_save
      if (allocated(err)) return
    endblock; else
      edd_save = var%edd
      var%edd = 10.0_dp**log10edd_new
      call self%set_temperature(T_new, err=err)
      if (allocated(err)) var%edd = edd_save
      if (allocated(err)) return
    endif

  contains
    subroutine solve_layer_pressure(layer, xzero, err)
      integer, intent(in) :: layer
      real(dp), intent(out) :: xzero
      character(:), allocatable, intent(out) :: err

      real(dp) :: fzero, xcenter, xlower, xupper, flower, fupper, width
      real(dp) :: reference_density
      integer :: iflag, nexpand

      ! pressure_residual has the callback signature required by brent_class.
      ! residual_layer selects the layer; all other callback context is
      ! read-only and consists of the input profiles, atmospheric state, and
      ! the previously solved pressure and temperature in hydrostatic mode.
      residual_layer = layer

      if (hydro_pressure_) then
        if (layer == 1) then
          xcenter = log10(max(self%wrk%pressure_hydro(1), small_real))
          width = max(abs(log10(Psurf_initial) - xcenter), 1.0e-3_dp)
        else
          xcenter = log10P_previous - &
                    (mubar_base(layer)*var%grav(layer)*var%dz(layer)*log10e)/ &
                    (N_avo*k_boltz*max(temperature_previous,small_real))
          width = max(abs(log10P_previous - xcenter), 1.0e-3_dp)
        endif
      else
        reference_density = density_base(layer)
        xcenter = log10(reference_density*k_boltz*max(var%temperature(layer),small_real))
        width = 1.0e-3_dp
      endif

      xlower = xcenter - width
      xupper = xcenter + width
      flower = pressure_residual(root_solver, xlower)
      fupper = pressure_residual(root_solver, xupper)
      nexpand = 0
      do while (.not. opposite_signs(flower, fupper))
        width = 2.0_dp*width
        xlower = xcenter - width
        xupper = xcenter + width
        flower = pressure_residual(root_solver, xlower)
        fupper = pressure_residual(root_solver, xupper)
        nexpand = nexpand + 1
        if (nexpand >= max_bracket_expansions) then
          write(layer_string,'(i0)') layer
          err = 'Could not bracket the pressure root in layer '//trim(layer_string)
          return
        endif
      enddo

      call root_solver%find_zero(xlower, xupper, root_tol, xzero, fzero, iflag, flower, fupper)
      if (iflag /= 0 .or. .not. ieee_is_finite(xzero)) then
        write(layer_string,'(i0)') layer
        err = 'The scalar pressure solve failed in layer '//trim(layer_string)
        return
      endif

    end subroutine

    function pressure_residual(me, x) result(residual)
      class(brent_class), intent(inout) :: me
      real(dp), intent(in) :: x
      real(dp) :: residual

      real(dp) :: temperature_trial, density_trial, mubar_trial, Psurf_trial

      temperature_trial = profile_value(x, log10P_in, T_in)

      if (hydro_pressure_) then
        if (residual_layer == 1) then
          call bottom_column_state(temperature_trial, density_trial, mubar_trial, Psurf_trial)
          residual = x - log10(Psurf_trial) + &
                     (mubar_trial*var%grav(1)*0.5_dp*var%dz(1)*log10e)/ &
                     (N_avo*k_boltz*temperature_trial)
        else
          residual = x - log10P_previous + &
                     (mubar_base(residual_layer)*var%grav(residual_layer)*var%dz(residual_layer)*log10e)/ &
                     (N_avo*k_boltz*0.5_dp*(temperature_previous + temperature_trial))
        endif
      else
        if (residual_layer == 1) then
          call bottom_column_state(temperature_trial, density_trial, mubar_trial, Psurf_trial)
        else
          density_trial = density_base(residual_layer)
        endif
        residual = x - log10(density_trial*k_boltz*temperature_trial)
      endif

    end function

    subroutine bottom_column_state(temperature, density_bottom, mubar_bottom, surface_pressure)
      real(dp), intent(in) :: temperature
      real(dp), intent(out) :: density_bottom, mubar_bottom, surface_pressure

      real(dp) :: usol_bottom(dat%nq), Psat, column_mass
      integer :: gas_ind, particle_ind

      ! This helper reads the saved composition and boundary conditions but
      ! applies them to a local bottom-layer copy; it never modifies self.
      usol_bottom = usol_base(:,1)
      do gas_ind = 1,dat%nq
        if (var%lowerboundcond(gas_ind) == DensityBC) then
          usol_bottom(gas_ind) = var%lower_fix_den(gas_ind)
        elseif (var%lowerboundcond(gas_ind) == PressureBC) then
          Psat = huge(1.0_dp)
          if (dat%gas_particle_ind(gas_ind) /= 0) then
            particle_ind = dat%gas_particle_ind(gas_ind)
            Psat = dat%particle_sat(particle_ind)%sat_pressure(temperature)* &
                   var%cond_params(particle_ind)%RHc
          endif
          usol_bottom(gas_ind) = min(var%lower_fix_press(gas_ind),Psat)/ &
                                   (k_boltz*temperature)
        endif
      enddo

      density_bottom = sum(usol_bottom(dat%ng_1:))
      mubar_bottom = sum(dat%species_mass(dat%ng_1:dat%nq)* &
                         usol_bottom(dat%ng_1:dat%nq))/density_bottom

      column_mass = density_bottom*mubar_bottom*var%grav(1)*var%dz(1)
      do gas_ind = 2,var%nz
        column_mass = column_mass + density_base(gas_ind)*mubar_base(gas_ind)* &
                                    var%grav(gas_ind)*var%dz(gas_ind)
      enddo
      surface_pressure = column_mass/N_avo

    end subroutine

    pure function profile_value(x, pressure_grid, values) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: pressure_grid(:)
      real(dp), intent(in) :: values(:)
      real(dp) :: value

      real(dp) :: fraction
      integer :: low, high, mid

      if (x <= pressure_grid(1)) then
        value = values(1)
      elseif (x >= pressure_grid(size(pressure_grid))) then
        value = values(size(values))
      else
        low = 1
        high = size(pressure_grid)
        do while (high-low > 1)
          mid = (low+high)/2
          if (x >= pressure_grid(mid)) then
            low = mid
          else
            high = mid
          endif
        enddo
        fraction = (x-pressure_grid(low))/(pressure_grid(high)-pressure_grid(low))
        value = values(low) + fraction*(values(high)-values(low))
      endif

    end function

    pure function opposite_signs(a, b) result(opposite)
      real(dp), intent(in) :: a, b
      logical :: opposite

      opposite = (a <= 0.0_dp .and. b >= 0.0_dp) .or. &
                 (a >= 0.0_dp .and. b <= 0.0_dp)

    end function
  end subroutine

  function TOA_at_pressure(self, usol, TOA_pressure, err) result(top_atmos)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: TOA_pressure !! dynes/cm^2
    character(:), allocatable, intent(out) :: err
    real(dp) :: top_atmos !! cm

    type(MinpackHybrd1Vars) :: mv

    mv = MinpackHybrd1Vars(1,tol=1.0e-5_dp)
    mv%x(1) = self%var%z(self%var%nz)
    call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
    if (mv%info == 0 .or. mv%info > 1) then
      err = 'hybrd1 root solve failed in TOA_at_pressure.'
      return
    elseif (mv%info < 0) then
      err = 'hybrd1 root solve failed in TOA_at_pressure: '//err
      return
    endif

    top_atmos = mv%x(1)

  contains
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_
      real(dp) :: TOA_pressure_
      TOA_pressure_ = pressure_at_TOA(self, usol, x_(1), err)
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      fvec_(1) = log10(TOA_pressure_) - log10(TOA_pressure)
    end subroutine
  end function

  function pressure_at_TOA(self, usol, top_atmos_new, err) result(TOA_pressure)
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: top_atmos_new !! cm
    character(:), allocatable, intent(out) :: err
    real(dp) :: TOA_pressure !! dynes/cm^2

    real(dp), allocatable :: usol_new(:,:)
    real(dp), allocatable :: z_new(:), dz_new(:), grav_new(:)
    real(dp), allocatable :: temperature_new(:), edd_new(:), particle_radius_new(:,:)
    real(dp), allocatable :: pressure_new(:)

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! Allocate
    allocate(usol_new(dat%nq,var%nz))
    allocate(z_new(var%nz),dz_new(var%nz),grav_new(var%nz))
    allocate(temperature_new(var%nz),edd_new(var%nz), particle_radius_new(dat%np,var%nz))
    allocate(pressure_new(var%nz))

    call properties_for_new_TOA(self, usol, top_atmos_new, &
      z_new, dz_new, grav_new, temperature_new, edd_new, usol_new, &
      particle_radius_new, pressure_new, err)
    if (allocated(err)) return

    TOA_pressure = pressure_new(var%nz)

  end function

  subroutine properties_for_new_TOA(self, usol, top_atmos_new, &
                              z_new, dz_new, grav_new, temperature_new, edd_new, usol_new, &
                              particle_radius_new, pressure_new, err)
    use photochem_enum, only: DensityBC, PressureBC
    use futils, only: interp
    use photochem_eqns, only: vertical_grid, press_and_den, gravity
    use photochem_const, only: small_real, k_boltz
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: top_atmos_new !! cm
    real(dp), intent(out) :: z_new(:)
    real(dp), intent(out) :: dz_new(:)
    real(dp), intent(out) :: grav_new(:)
    real(dp), intent(out) :: temperature_new(:)
    real(dp), intent(out) :: edd_new(:)
    real(dp), intent(out) :: usol_new(:,:)
    real(dp), intent(out) :: particle_radius_new(:,:)
    real(dp), intent(out) :: pressure_new(:)
    character(:), allocatable, intent(out) :: err

    real(dp) :: Psat
    real(dp), allocatable :: mix(:,:), mix_new(:,:)
    real(dp), allocatable :: density(:), density_new(:)
    integer :: i, j, ierr

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! Allocate
    allocate(mix(dat%nq,var%nz), mix_new(dat%nq,var%nz))
    allocate(density(var%nz),density_new(var%nz))

    ! Remake the vertical grid and gravity
    call vertical_grid(var%bottom_atmos, top_atmos_new, &
                      var%nz, z_new, dz_new)
    call gravity(dat%planet_radius, dat%planet_mass, &
                var%nz, z_new, grav_new)

    ! Temperature
    call interp(var%nz, var%nz, z_new, var%z, var%temperature, temperature_new, ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
      return
    endif

    ! Eddy diffusion
    call interp(var%nz, var%nz, z_new, var%z, log10(max(var%edd,1.0e-40_dp)), edd_new, ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
      return
    endif
    edd_new = 10.0_dp**edd_new

    ! determine mixing ratios and density
    do j = 1,var%nz
      density(j) = sum(usol(dat%ng_1:,j))
      mix(:,j) = usol(:,j)/density(j) ! mixing ratios
    enddo
    ! Interpolate mixing ratios, with constant extrapolation
    do i = 1,dat%nq
      call interp(z_new, var%z, log10(max(mix(i,:),small_real)), mix_new(i,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Subroutine interp returned an error.'
        return
      endif
    enddo
    mix_new = 10.0_dp**mix_new

    ! Interpolate density, with linear extrapolation
    call interp(z_new, var%z, log10(density), density_new, linear_extrap=.true., ierr=ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
      return
    endif
    density_new = 10.0_dp**density_new

    ! Compute usol_new with mixing ratios and densities
    do i = 1,var%nz
      usol_new(:,i) = mix_new(:,i)*density_new(i)
    enddo 

    ! Particle radii
    if (dat%there_are_particles) then
      do i = 1,dat%npq
        call interp(var%nz, var%nz, z_new, var%z, &
                    log10(max(var%particle_radius(i,:),small_real)), particle_radius_new(i,:), ierr)
        if (ierr /= 0) then
          err = 'Subroutine interp returned an error.'
          return
        endif
      enddo
      particle_radius_new = 10.0_dp**particle_radius_new
    endif

    ! Account for fixed surface mixing ratios
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == DensityBC) then
        usol_new(i,1) = var%lower_fix_den(i)
      elseif (var%lowerboundcond(i) == PressureBC) then
        Psat = huge(1.0_dp)
        if (dat%gas_particle_ind(i) /= 0) then
          j = dat%gas_particle_ind(i)
          Psat = dat%particle_sat(j)%sat_pressure(var%temperature(1))*var%cond_params(j)%RHc
        endif
        usol_new(i,1) = min(var%lower_fix_press(i), Psat)/(k_boltz*temperature_new(1))
      endif
    enddo

    pressure_new = density_new*k_boltz*temperature_new

  end subroutine


  module subroutine update_vertical_grid(self, TOA_alt, TOA_pressure, err)
    use photochem_input, only: interp2particlexsdata
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), optional, intent(in) :: TOA_alt !! cm
    real(dp), optional, intent(in) :: TOA_pressure !! dynes/cm^2
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: usol_new(:,:)
    real(dp), allocatable :: z_new(:), dz_new(:), grav_new(:)
    real(dp), allocatable :: temperature_new(:), edd_new(:), particle_radius_new(:,:)
    real(dp), allocatable :: pressure_new(:)

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk

    dat => self%dat
    var => self%var
    wrk => self%wrk

    if (present(TOA_alt) .and. present(TOA_pressure)) then
      err = 'Both "TOA_alt" and "TOA_pressure" can not be specified'
      return
    endif
    if (.not.present(TOA_alt) .and. .not.present(TOA_pressure)) then
      err = 'Either "TOA_alt" and "TOA_pressure" must be specified'
      return
    endif

    if (present(TOA_alt)) then
      if (TOA_alt < 0.0_dp) then
        err = '"TOA_alt" must be positive.'
        return
      endif
    endif

    if (present(TOA_pressure)) then
      if (TOA_pressure < 0.0_dp) then
        err = '"TOA_pressure" must be positive.'
        return
      endif
    endif

    ! Allocate
    allocate(usol_new(dat%nq,var%nz))
    allocate(z_new(var%nz),dz_new(var%nz),grav_new(var%nz))
    allocate(temperature_new(var%nz),edd_new(var%nz), particle_radius_new(dat%np,var%nz))
    allocate(pressure_new(var%nz))

    if (present(TOA_alt)) then
      var%top_atmos = TOA_alt
    endif
    if (present(TOA_pressure)) then
      ! Compute new TOA. We use wrk%usol
      var%top_atmos = TOA_at_pressure(self, wrk%usol, TOA_pressure, err)
      if (allocated(err)) return
    endif

    ! Compute properties associated with new TOA
    call properties_for_new_TOA(self, wrk%usol, var%top_atmos, &
        z_new, dz_new, grav_new, temperature_new, edd_new, usol_new, &
        particle_radius_new, pressure_new, err)
    if (allocated(err)) return

    var%z = z_new
    var%dz = dz_new
    var%grav = grav_new
    var%temperature = temperature_new
    var%edd = edd_new
    wrk%usol = usol_new
    var%usol_init = usol_new
    var%particle_radius = particle_radius_new

    ! Get new optical properties associated with new particle radii
    call interp2particlexsdata(dat, var, err)
    if (allocated(err)) return

    ! Update variables that depend on temperature
    if (dat%gas_rainout) then
      call self%set_temperature(var%temperature, var%trop_alt, err)
    else
      call self%set_temperature(var%temperature, err=err)
    endif
    if (allocated(err)) return

  end subroutine

  ! Below is mostly stuff needed in `evolve` routine

  module subroutine rebin_update_vertical_grid(self, usol_old, top_atmos, usol_new, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_old(:,:)
    real(dp), intent(in) :: top_atmos
    real(dp), intent(out) :: usol_new(:,:)
    character(:), allocatable, intent(out) :: err

    call rebin_densities(self, usol_old, top_atmos, usol_new, err)
    if(allocated(err)) return
    call update_vertical_grid_file(self, usol_new, top_atmos, err)
    if(allocated(err)) return

  end subroutine

  subroutine rebin_densities(self, usol_old, top_atmos, usol_new, err)
    use photochem_eqns, only: vertical_grid
    use futils, only: conserving_rebin
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol_old(:,:)
    real(dp), intent(in) :: top_atmos
    real(dp), intent(out) :: usol_new(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: ze_old(:), ze_new(:)
    real(dp), allocatable :: z_old(:), dz_old(:), z_new(:), dz_new(:)
    integer :: i, ierr

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! save z and dz
    allocate(z_old(var%nz), dz_old(var%nz))
    allocate(z_new(var%nz), dz_new(var%nz))
    z_old = var%z
    dz_old = var%dz

    ! compute the new grid
    call vertical_grid(var%bottom_atmos, top_atmos, &
                       var%nz, z_new, dz_new)

    ! We do a conserving rebin of the densities
    allocate(ze_old(var%nz+1))
    allocate(ze_new(var%nz+1))

    ze_old(1) = z_old(1) - 0.5_dp*dz_old(1)
    do i = 1,size(z_old)
      ze_old(i+1) = z_old(i) + 0.5_dp*dz_old(i)
    enddo
    ze_new(1) = z_new(1) - 0.5_dp*dz_new(1)
    do i = 1,var%nz
      ze_new(i+1) = z_new(i) + 0.5_dp*dz_new(i)
    enddo

    do i = 1,dat%nq
      call conserving_rebin(ze_old, usol_old(i,:), ze_new, usol_new(i,:), ierr)
      if (ierr /= 0) then
        err = 'subroutine conserving_rebin returned an error'
        return
      endif
    enddo

  end subroutine

  subroutine update_vertical_grid_file(self, usol_new, top_atmos, err)
    use photochem_input, only: interp2particlexsdata, interp2xsdata, compute_gibbs_energy
    use photochem_eqns, only: gravity, vertical_grid
    use futils, only: interp
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_new(:,:)
    real(dp), intent(in) :: top_atmos
    character(:), allocatable, intent(out) :: err

    integer :: i, ierr

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! set the new top of the  atmosphere
    var%top_atmos = top_atmos

    ! remake the vertical grid and gravity
    call vertical_grid(var%bottom_atmos, var%top_atmos, &
                      var%nz, var%z, var%dz)
    call gravity(dat%planet_radius, dat%planet_mass, &
                var%nz, var%z, var%grav)

    ! We will assume Temperature, eddy diffusion, and particle radius from the file
    call interp(var%nz, dat%nzf, var%z, dat%z_file, dat%T_file, var%Temperature, ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
      return
    endif

    call interp(var%nz, dat%nzf, var%z, dat%z_file, log10(abs(dat%edd_file)), var%edd, ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
      return
    endif
    var%edd = 10.0_dp**var%edd

    if (dat%there_are_particles) then
      do i = 1,dat%npq
        call interp(var%nz, dat%nzf, var%z, dat%z_file, &
                    log10(abs(dat%particle_radius_file(i,:))), var%particle_radius(i,:), ierr)
        if (ierr /= 0) then
          err = 'Subroutine interp returned an error.'
          return
        endif
      enddo
      var%particle_radius = 10.0_dp**var%particle_radius
    endif

    call interp2particlexsdata(dat, var, err)
    if (allocated(err)) return

    ! all below depends on Temperature
    call interp2xsdata(dat, var, err)
    if (allocated(err)) return
    
    if (dat%reverse) then
      call compute_gibbs_energy(dat, var, err)
      if (allocated(err)) return
    endif
    
    if (dat%gas_rainout) then
      call self%set_trop_ind(usol_new, err)
    else
      var%trop_ind = 1
    endif

  end subroutine

  module subroutine regrid_prep_atmosphere(self, usol_new, top_atmos, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_new(:,:)
    real(dp), intent(in) :: top_atmos
    character(:), allocatable, intent(out) :: err

    call update_vertical_grid_file(self, usol_new, top_atmos, err)
    if(allocated(err)) return
    call self%prep_atmosphere(usol_new, err)
    if(allocated(err)) return

  end subroutine

end submodule
