
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
    
    call self%require_atmosphere_initialized('out2atmosphere_txt', err)
    if (allocated(err)) return

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
    
    call self%require_atmosphere_initialized('gas_fluxes', err)
    if (allocated(err)) return

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
    use iso_c_binding, only: c_associated
    use photochem_input, only: refresh_temperature_dependent_state
    
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: temperature(:)
    real(dp), optional, intent(in) :: trop_alt
    character(:), allocatable, intent(out) :: err
    
    type(PhotochemVars), pointer :: var
    type(PhotochemVars) :: var_save
    real(dp), allocatable :: usol_start(:,:)
    character(:), allocatable :: original_err, rollback_err

    var => self%var
    
    call self%require_atmosphere_initialized('set_temperature', err)
    if (allocated(err)) return

    if (c_associated(self%wrk%sun%cvode_mem)) then
      err = "Can not change the temperature profile while a CVODE stepper "// &
            "is initialized. Call 'destroy_stepper' first."
      return
    endif

    if (var%press_temp_edd_profile%enabled) then
      err = "The persistent pressure-temperature-eddy profile is enabled. "// &
            "Call 'clear_press_temp_edd_profile' before 'set_temperature'."
      return
    endif
    
    if (size(temperature) /= var%nz) then
      err = "temperature has the wrong input dimension"
      return
    endif

    ! Save in case there is an issue
    var_save = var
    
    ! Commit the new temperature
    var%temperature = temperature

    ! Update temperature-dependent state in `var`
    call refresh_temperature_dependent_state(self%dat, var, trop_alt, err)
    if (allocated(err)) then
      var = var_save
      return
    endif

    ! Fill wrk with new values. Keep the input separate from wrk%usol because
    ! prep_atmosphere has distinct input and output arrays.
    usol_start = self%wrk%usol
    call self%prep_atmosphere(usol_start, err)
    if (allocated(err)) then
      ! Restore both the configuration and derived work state if preparation
      ! failed after making partial changes.
      original_err = err
      var = var_save
      call self%prep_atmosphere(usol_start, rollback_err)
      if (allocated(rollback_err)) then
        err = original_err//' Rollback failed: '//rollback_err
      else
        err = original_err
      endif
      return
    endif
    
  end subroutine

  module subroutine set_press_temp_edd(self, P, T, edd, trop_p, hydro_pressure, err)
    use iso_c_binding, only: c_associated
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: P(:)
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: edd(:)
    real(dp), optional, intent(in) :: trop_p
    logical, optional, intent(in) :: hydro_pressure
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: T_grid(:), edd_grid(:), log10P_grid(:)
    real(dp), allocatable :: edd_save(:)
    real(dp) :: trop_alt
    real(dp) :: trop_p_
    logical :: has_trop_p, hydro_pressure_

    call self%require_atmosphere_initialized('set_press_temp_edd', err)
    if (allocated(err)) return

    if (c_associated(self%wrk%sun%cvode_mem)) then
      err = "Can not set a pressure-temperature-eddy profile while a CVODE "// &
            "stepper is initialized. Call 'destroy_stepper' first."
      return
    endif

    if (self%var%press_temp_edd_profile%enabled) then
      err = "The persistent pressure-temperature-eddy profile is enabled. "// &
            "Call 'clear_press_temp_edd_profile' before 'set_press_temp_edd'."
      return
    endif

    trop_p_ = -1.0_dp
    if (present(trop_p)) trop_p_ = trop_p
    has_trop_p = trop_p_ > 0.0_dp
    if (present(hydro_pressure)) then
      hydro_pressure_ = hydro_pressure
    else
      hydro_pressure_ = .true.
    endif

    ! Allocate only after the lifecycle checks above. In particular, an
    ! uninitialized object must return the normal API error without using an
    ! undefined `var%nz` as an automatic-array extent.
    allocate(T_grid(self%var%nz), edd_grid(self%var%nz), log10P_grid(self%var%nz), &
             edd_save(self%var%nz))

    ! First compute the mapping without changing model state. This kernel is
    ! also suitable for applying a persistent pressure-based profile to an
    ! arbitrary trial composition during a future RHS evaluation.
    call map_press_temp_edd( &
      self, &
      self%wrk%usol, &
      P, &
      T, &
      edd, &
      trop_p=trop_p_, &
      hydro_pressure=hydro_pressure_, &
      grid_z=self%var%z, &
      grid_dz=self%var%dz, &
      grid_grav=self%var%grav, &
      temperature_reference=self%var%temperature, &
      pressure_reference=self%wrk%pressure_hydro, &
      T_grid=T_grid, &
      edd_grid=edd_grid, &
      log10P_grid=log10P_grid, &
      trop_alt=trop_alt, &
      err=err &
    )
    if (allocated(err)) return

    ! Commit point: the mapping above only reads self. Updating Kzz and
    ! calling set_temperature below are the only operations that mutate it.
    edd_save = self%var%edd
    self%var%edd = edd_grid
    if (has_trop_p) then
      call self%set_temperature(T_grid, trop_alt, err)
    else
      call self%set_temperature(T_grid, err=err)
    endif
    if (allocated(err)) self%var%edd = edd_save

  end subroutine

  module subroutine set_press_temp_edd_profile(self, P, T, edd, trop_p, hydro_pressure, &
                                               maintain_toa_pressure, target_pressure, err)
    use iso_c_binding, only: c_associated
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: P(:)
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: edd(:)
    real(dp), optional, intent(in) :: trop_p
    logical, optional, intent(in) :: hydro_pressure
    logical, optional, intent(in) :: maintain_toa_pressure
    real(dp), optional, intent(in) :: target_pressure
    character(:), allocatable, intent(out) :: err

    type(PhotochemVars) :: var_save
    real(dp), allocatable :: usol_start(:,:)
    logical :: maintain_toa_pressure_
    character(:), allocatable :: original_err, rollback_err

    call self%require_atmosphere_initialized('set_press_temp_edd_profile', err)
    if (allocated(err)) return

    if (c_associated(self%wrk%sun%cvode_mem)) then
      err = "Can not set a persistent pressure-temperature-eddy profile while "// &
            "a stepper is initialized. Call 'destroy_stepper' first."
      return
    endif

    maintain_toa_pressure_ = .true.
    if (present(maintain_toa_pressure)) maintain_toa_pressure_ = maintain_toa_pressure
    if (present(target_pressure)) then
      if (.not.ieee_is_finite(target_pressure) .or. target_pressure <= 0.0_dp) then
        err = '`target_pressure` must be finite and positive.'
        return
      endif
    endif

    var_save = self%var

    self%var%press_temp_edd_profile%pressure = P
    self%var%press_temp_edd_profile%temperature = T
    self%var%press_temp_edd_profile%edd = edd
    if (present(trop_p)) then
      self%var%press_temp_edd_profile%trop_p = trop_p
    else
      self%var%press_temp_edd_profile%trop_p = -1.0_dp
    endif
    if (present(hydro_pressure)) then
      self%var%press_temp_edd_profile%hydro_pressure = hydro_pressure
    else
      self%var%press_temp_edd_profile%hydro_pressure = .true.
    endif
    self%var%press_temp_edd_profile%enabled = .true.

    ! prep_atmosphere validates and applies the stored profile through the
    ! same path used by RHS and Jacobian evaluations.
    ! Keep the input separate from wrk%usol because preparation writes the
    ! canonical working state through a distinct output argument.
    usol_start = self%wrk%usol
    call self%prep_atmosphere(usol_start, err)
    if (allocated(err)) then
      original_err = err
      self%var = var_save
      call self%prep_atmosphere(usol_start, rollback_err)
      if (allocated(rollback_err)) then
        err = original_err//' Rollback failed: '//rollback_err
      else
        err = original_err
      endif
      return
    endif

    ! The profile is now committed. Couple TOA maintenance to the same
    ! successful profile installation, using the requested target or the
    ! standard default when the caller did not provide one.
    self%var%toa_pressure_maintenance%enabled = maintain_toa_pressure_
    if (maintain_toa_pressure_) then
      if (present(target_pressure)) then
        self%var%toa_pressure_maintenance%target_pressure = target_pressure
      else
        self%var%toa_pressure_maintenance%target_pressure = 0.1_dp
      endif
    endif

  end subroutine

  module subroutine clear_press_temp_edd_profile(self, err)
    use iso_c_binding, only: c_associated
    class(EvoAtmosphere), intent(inout) :: self
    character(:), allocatable, intent(out) :: err

    if (c_associated(self%wrk%sun%cvode_mem)) then
      err = "Can not clear the persistent pressure-temperature-eddy profile while "// &
            "a stepper is initialized. Call 'destroy_stepper' first."
      return
    endif

    call reset_press_temp_edd_profile(self%var)

  end subroutine

  module subroutine reset_press_temp_edd_profile(var)
    type(PhotochemVars), intent(inout) :: var

    var%press_temp_edd_profile%enabled = .false.
    var%press_temp_edd_profile%hydro_pressure = .true.
    var%press_temp_edd_profile%trop_p = -1.0_dp
    ! A cleared persistent profile cannot support automatic TOA maintenance.
    ! Preserve the maintenance tuning values, but clear its active state and
    ! any target that was coupled to the previous profile.
    var%toa_pressure_maintenance%enabled = .false.
    var%toa_pressure_maintenance%target_pressure = 0.0_dp
    if (allocated(var%press_temp_edd_profile%pressure)) &
      deallocate(var%press_temp_edd_profile%pressure)
    if (allocated(var%press_temp_edd_profile%temperature)) &
      deallocate(var%press_temp_edd_profile%temperature)
    if (allocated(var%press_temp_edd_profile%edd)) &
      deallocate(var%press_temp_edd_profile%edd)

  end subroutine

  module subroutine apply_press_temp_edd_profile(self, usol_in, err)
    use photochem_input, only: refresh_temperature_dependent_state
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_in(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp) :: T_grid(self%var%nz), edd_grid(self%var%nz)
    real(dp) :: log10P_grid(self%var%nz), trop_alt

    if (.not. self%var%press_temp_edd_profile%enabled) return

    call map_press_temp_edd( &
      self, &
      usol_in, &
      self%var%press_temp_edd_profile%pressure, &
      self%var%press_temp_edd_profile%temperature, &
      self%var%press_temp_edd_profile%edd, &
      trop_p=self%var%press_temp_edd_profile%trop_p, &
      hydro_pressure=self%var%press_temp_edd_profile%hydro_pressure, &
      grid_z=self%var%z, &
      grid_dz=self%var%dz, &
      grid_grav=self%var%grav, &
      temperature_reference=self%var%temperature, &
      pressure_reference=self%wrk%pressure_hydro, &
      T_grid=T_grid, &
      edd_grid=edd_grid, &
      log10P_grid=log10P_grid, &
      trop_alt=trop_alt, &
      err=err &
    )
    if (allocated(err)) return

    ! Commit only after the profile has mapped successfully and the
    ! tropopause has been validated. Do not call set_temperature here: that
    ! routine calls prep_atmosphere and would recurse back into this routine.
    self%var%temperature = T_grid
    self%var%edd = edd_grid
    call refresh_temperature_dependent_state(self%dat, self%var, trop_alt, err)

  end subroutine

  module subroutine map_press_temp_edd(self, usol, P, T, edd, trop_p, hydro_pressure, &
                                       grid_z, grid_dz, grid_grav, temperature_reference, &
                                       pressure_reference, &
                                       T_grid, edd_grid, log10P_grid, trop_alt, err)
    use futils, only: interp, brent_class
    use ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use photochem_const, only: small_real, k_boltz, N_avo
    use photochem_enum, only: DensityBC, PressureBC
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: P(:)
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: edd(:)
    real(dp), intent(in) :: trop_p
    logical, intent(in) :: hydro_pressure
    real(dp), intent(in) :: grid_z(:), grid_dz(:), grid_grav(:)
    real(dp), intent(in) :: temperature_reference(:), pressure_reference(:)
    real(dp), intent(out) :: T_grid(:)
    real(dp), intent(out) :: edd_grid(:)
    real(dp), intent(out) :: log10P_grid(:)
    real(dp), intent(out) :: trop_alt
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: log10P_in(:), T_in(:), log10edd_in(:)
    real(dp), allocatable :: P_wrk(:), log10edd_new(:)
    real(dp), allocatable :: usol_base(:,:), density_base(:), mubar_base(:)
    real(dp), allocatable :: usol_bottom(:)
    real(dp) :: xzero, Psurf_initial, Psurf_final
    real(dp) :: log10P_previous, temperature_previous
    real(dp) :: trop_alt_array(1)
    logical :: has_trop_p
    integer :: ierr, i, j, residual_layer
    character(32) :: layer_string
    character(:), allocatable :: residual_error

    real(dp), parameter :: log10e = log10(exp(1.0_dp))
    real(dp), parameter :: root_tol = 1.0e-10_dp
    integer, parameter :: max_bracket_expansions = 60

    type(brent_class) :: root_solver

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var

    ! Check inputs
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz) then
      err = '"usol" has the wrong dimensions'
      return
    endif
    if (size(T_grid) /= var%nz .or. size(edd_grid) /= var%nz .or. &
        size(log10P_grid) /= var%nz) then
      err = 'Pressure-profile mapping output arrays have the wrong dimensions'
      return
    endif
    if (size(grid_z) /= var%nz .or. size(grid_dz) /= var%nz .or. &
        size(grid_grav) /= var%nz .or. &
        size(temperature_reference) /= var%nz .or. &
        size(pressure_reference) /= var%nz) then
      err = 'Pressure-profile mapping grid arrays have the wrong dimensions'
      return
    endif
    if (.not.all(ieee_is_finite(grid_z)) .or. &
        .not.all(ieee_is_finite(grid_dz)) .or. any(grid_dz <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(grid_grav)) .or. any(grid_grav <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(temperature_reference)) .or. &
        any(temperature_reference <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(pressure_reference)) .or. &
        any(pressure_reference <= 0.0_dp)) then
      err = 'Pressure-profile mapping requires finite, positive grid reference state'
      return
    endif
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
    if (.not. ieee_is_finite(trop_p)) then
      err = '"trop_p" must be finite; use a non-positive value to omit it'
      return
    endif
    has_trop_p = trop_p > 0.0_dp
    if (has_trop_p .and. .not.dat%gas_rainout) then
      err = '"trop_p" can only be supplied when gas rainout is enabled.'
      return
    endif
    if (dat%gas_rainout .and. .not.has_trop_p) then
      err = '"trop_p" is a required input.'
      return
    endif

    ! Allocate work space
    allocate(P_wrk(var%nz), log10edd_new(var%nz))
    allocate(usol_base(dat%nq,var%nz), density_base(var%nz), mubar_base(var%nz))
    allocate(usol_bottom(dat%nq))

    ! Copy and clip the current state in the same way as prep_atm_evo_gas.
    call clip_usol(usol, usol_base)

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

    call bottom_column_state(temperature_reference(1), density_base(1), &
                             mubar_base(1), Psurf_initial, err)
    if (allocated(err)) return

    ! Copy the input profiles into log-pressure space. The arrays remain at
    ! the supplied pressure points; values below the deepest point are
    ! extrapolated directly from the deepest two points by evaluate_profile.
    allocate(log10P_in(size(P)),T_in(size(P)),log10edd_in(size(P)))
    log10P_in = log10(P)
    T_in = T
    log10edd_in = log10(edd)
    ! Flip order for interpolation purposes
    log10P_in = log10P_in(size(log10P_in):1:-1)
    T_in = T_in(size(log10P_in):1:-1)
    log10edd_in = log10edd_in(size(log10P_in):1:-1)

    call root_solver%set_function(pressure_residual)

    ! Solve one scalar pressure equation per layer. In hydrostatic mode the
    ! recurrence is triangular, so each solved layer provides the lower
    ! pressure and temperature boundary for the next layer.
    do i = 1,var%nz

      ! Perform the solve
      call solve_layer_pressure(i, xzero, err)
      if (allocated(err)) return

      ! Save the result, and evaluate T and Kzz at the root.
      log10P_grid(i) = xzero
      P_wrk(i) = 10.0_dp**xzero
      call evaluate_profile(xzero, log10P_in, T_in, T_grid(i), 'temperature', .true., err)
      if (allocated(err)) return

      call evaluate_profile(xzero, log10P_in, log10edd_in, log10edd_new(i), &
                            'eddy-diffusion', .false., err)
      if (allocated(err)) return

      ! Check for problems
      if (.not. ieee_is_finite(P_wrk(i)) .or. P_wrk(i) <= 0.0_dp) then
        write(layer_string,'(i0)') i
        err = 'The pressure solve produced a non-finite or non-positive pressure in layer '// &
              trim(layer_string)
        return
      endif
      if (i > 1) then
        if (P_wrk(i) >= P_wrk(i-1)) then
          if (hydro_pressure) then
            write(layer_string,'(i0)') i
            err = 'The hydrostatic pressure does not decrease upward at layer '//trim(layer_string)
            return
          elseif (has_trop_p) then
            write(layer_string,'(i0," and ",i0)') i-1, i
            err = 'Cannot determine the tropopause altitude because actual pressure is not '// &
                  'strictly decreasing between layers '//trim(layer_string)
            return
          endif
        endif
      endif

      if (i == 1) then; block
        real(dp) :: density_final, mubar_final

        ! Little sanity check to make sure pressure is decreasing.

        call bottom_column_state(T_grid(1), density_final, mubar_final, Psurf_final, err)
        if (allocated(err)) return
        if (hydro_pressure .and. P_wrk(1) >= Psurf_final) then
          err = 'The solved bottom-layer pressure is not below the surface pressure'
          return
        endif
      end block; endif

      ! Save the previous solve. This is needed for initial guesses and for the residuals on
      ! subsequent layers.
      log10P_previous = xzero
      temperature_previous = T_grid(i)
    enddo

    ! Ensure eddy diffusion is valid
    if (.not. all(ieee_is_finite(log10edd_new))) then
      err = 'The mapped eddy-diffusion profile is not finite'
      return
    endif
    edd_grid = 10.0_dp**log10edd_new
    if (.not. all(ieee_is_finite(edd_grid)) .or. any(edd_grid <= 0.0_dp)) then
      err = 'The mapped eddy-diffusion profile is not finite and positive'
      return
    endif

    ! Compute the tropopause altitude
    trop_alt = 0.0_dp
    if (has_trop_p) then
      call interp([log10(trop_p)], log10P_grid(var%nz:1:-1), &
                  grid_z(var%nz:1:-1), trop_alt_array, ierr=ierr)
      if (ierr /= 0) then
        err = 'Subroutine interp returned an error.'
        return
      endif
      trop_alt = trop_alt_array(1)
    endif

  contains
    subroutine solve_layer_pressure(layer, xzero_, err_)
      integer, intent(in) :: layer
      real(dp), intent(out) :: xzero_
      character(:), allocatable, intent(out) :: err_

      real(dp) :: fzero, xcenter, xlower, xupper, flower, fupper, width
      real(dp) :: reference_density
      integer :: iflag, nexpand

      ! pressure_residual has the callback signature required by brent_class.
      ! residual_layer selects the layer; all other callback context is
      ! read-only and consists of the input profiles, atmospheric state, and
      ! the previously solved pressure and temperature in hydrostatic mode.
      residual_layer = layer
      if (allocated(residual_error)) deallocate(residual_error)

      ! `xcenter` is the root estimate, while `width` is how far above and below
      ! `xcenter` we will search for a bracket that we can use for the solve.
      if (hydro_pressure) then
        if (layer == 1) then
          ! Existing bottom pressure is the root estimate.
          xcenter = log10(max(pressure_reference(1), small_real))
          width = max(abs(log10(Psurf_initial) - xcenter), 1.0e-3_dp)
        else
          ! Subtract the expected hydrostatic drop from the previous pressure;
          ! use that drop for the initial bracket.
          xcenter = log10P_previous - &
                    (mubar_base(layer)*grid_grav(layer)*grid_dz(layer)*log10e)/ &
                    (N_avo*k_boltz*max(temperature_previous,small_real))
          width = max(abs(log10P_previous - xcenter), 1.0e-3_dp)
        endif
      else
        ! Use the ideal-gas pressure from the current density and temperature.
        reference_density = density_base(layer)
        xcenter = log10(reference_density*k_boltz* &
                        max(temperature_reference(layer),small_real))
        ! Start narrow; the bracket expands if the residual has no sign change.
        width = 1.0e-3_dp
      endif

      ! Now find a bracket that we can solve on.
      xlower = xcenter - width
      xupper = xcenter + width
      flower = pressure_residual(root_solver, xlower)
      if (allocated(residual_error)) then
        err_ = residual_error
        return
      endif
      fupper = pressure_residual(root_solver, xupper)
      if (allocated(residual_error)) then
        err_ = residual_error
        return
      endif
      nexpand = 0
      do while (.not. opposite_signs(flower, fupper))
        width = 2.0_dp*width
        xlower = xcenter - width
        xupper = xcenter + width
        flower = pressure_residual(root_solver, xlower)
        if (allocated(residual_error)) then
          err_ = residual_error
          return
        endif
        fupper = pressure_residual(root_solver, xupper)
        if (allocated(residual_error)) then
          err_ = residual_error
          return
        endif
        nexpand = nexpand + 1
        if (nexpand >= max_bracket_expansions) then
          write(layer_string,'(i0)') layer
          err_ = 'Could not bracket the pressure root in layer '//trim(layer_string)
          return
        endif
      enddo

      ! Perform the root solve.
      call root_solver%find_zero(xlower, xupper, root_tol, xzero_, fzero, iflag, flower, fupper)
      if (allocated(residual_error)) then
        err_ = residual_error
        return
      endif
      if (iflag /= 0 .or. .not. ieee_is_finite(xzero_)) then
        write(layer_string,'(i0)') layer
        err_ = 'The scalar pressure solve failed in layer '//trim(layer_string)
        return
      endif

    end subroutine

    function pressure_residual(me, x) result(residual)
      class(brent_class), intent(inout) :: me
      real(dp), intent(in) :: x
      real(dp) :: residual

      real(dp) :: temperature_trial, density_trial, mubar_trial, Psurf_trial

      ! Preserve the first invalid evaluation while Brent is unwinding. This
      ! prevents a later callback from replacing the useful diagnostic.
      if (allocated(residual_error)) then
        residual = ieee_value(0.0_dp, ieee_quiet_nan)
        return
      endif

      ! For out input log10 pressure `x`, compute the `temperature_trial`
      call evaluate_profile(x, log10P_in, T_in, temperature_trial, 'temperature', .true., residual_error)
      if (allocated(residual_error)) then
        residual = ieee_value(0.0_dp, ieee_quiet_nan)
        return
      endif

      if (hydro_pressure) then
        ! x is trial log10 pressure. Since the input temperature is T(P),
        ! each trial pressure determines the temperature used in the residual.
        if (residual_layer == 1) then
          ! Recompute the surface pressure because lower boundary conditions
          ! can depend on temperature. Enforce balance to the layer center,
          ! hence the half-layer thickness.
          call bottom_column_state(temperature_trial, density_trial, mubar_trial, &
                                   Psurf_trial, residual_error)
          if (allocated(residual_error)) then
            residual = ieee_value(0.0_dp, ieee_quiet_nan)
            return
          endif
          residual = x - log10(Psurf_trial) + &
                     (mubar_trial*grid_grav(1)*0.5_dp*grid_dz(1)*log10e)/ &
                     (N_avo*k_boltz*temperature_trial)
        else
          ! For higher layers, enforce the hydrostatic pressure drop from the
          ! previous layer using the average of the two temperatures.
          residual = x - log10P_previous + &
                     (mubar_base(residual_layer)*grid_grav(residual_layer)* &
                      grid_dz(residual_layer)*log10e)/ &
                     (N_avo*k_boltz*0.5_dp*(temperature_previous + temperature_trial))
        endif
      else
        ! Without hydrostatic coupling, enforce the local ideal-gas relation
        ! P = density*k_B*T at the trial pressure.
        if (residual_layer == 1) then
          call bottom_column_state(temperature_trial, density_trial, mubar_trial, &
                                   Psurf_trial, residual_error)
          if (allocated(residual_error)) then
            residual = ieee_value(0.0_dp, ieee_quiet_nan)
            return
          endif
        else
          density_trial = density_base(residual_layer)
        endif
        residual = x - log10(density_trial*k_boltz*temperature_trial)
      endif

      if (.not. ieee_is_finite(residual)) then
        residual_error = 'The pressure residual became non-finite during the pressure solve'
        residual = ieee_value(0.0_dp, ieee_quiet_nan)
      endif

    end function

    subroutine bottom_column_state(temperature, density_bottom, mubar_bottom, surface_pressure, err_)
      real(dp), intent(in) :: temperature
      real(dp), intent(out) :: density_bottom, mubar_bottom, surface_pressure
      character(:), allocatable, intent(out) :: err_

      real(dp) :: column_mass
      integer :: gas_ind

      ! Apply boundary conditions to a local bottom-layer copy; self is never
      ! modified by this calculation.
      usol_bottom = usol_base(:,1)
      call self%apply_lower_boundary_conditions(temperature, usol_bottom, err_)
      if (allocated(err_)) return

      density_bottom = sum(usol_bottom(dat%ng_1:))
      if (.not. ieee_is_finite(density_bottom) .or. density_bottom <= 0.0_dp) then
        err_ = 'The bottom-layer number density is not finite and positive'
        return
      endif
      mubar_bottom = sum(dat%species_mass(dat%ng_1:dat%nq)* &
                         usol_bottom(dat%ng_1:dat%nq))/density_bottom
      if (.not. ieee_is_finite(mubar_bottom) .or. mubar_bottom <= 0.0_dp) then
        err_ = 'The bottom-layer mean molecular weight is not finite and positive'
        return
      endif

      column_mass = density_bottom*mubar_bottom*grid_grav(1)*grid_dz(1)
      do gas_ind = 2,var%nz
        column_mass = column_mass + density_base(gas_ind)*mubar_base(gas_ind)* &
                                    grid_grav(gas_ind)*grid_dz(gas_ind)
      enddo
      surface_pressure = column_mass/N_avo
      if (.not. ieee_is_finite(surface_pressure) .or. surface_pressure <= 0.0_dp) then
        err_ = 'The bottom-layer surface pressure is not finite and positive'
        return
      endif

    end subroutine

    subroutine evaluate_profile(x, pressure_grid, values, value, profile_name, require_positive, err_)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: pressure_grid(:)
      real(dp), intent(in) :: values(:)
      real(dp), intent(out) :: value
      character(*), intent(in) :: profile_name
      logical, intent(in) :: require_positive
      character(:), allocatable, intent(out) :: err_

      real(dp) :: fraction
      integer :: low, high, mid

      if (x <= pressure_grid(1)) then
        ! Above the supplied top of the atmosphere, retain the top value.
        value = values(1)
      elseif (x >= pressure_grid(size(pressure_grid))) then
        ! Below the supplied deepest point, continue the deepest two-point
        ! slope in log-pressure rather than attaching a moving surface knot.
        fraction = (x-pressure_grid(size(pressure_grid))) / &
                   (pressure_grid(size(pressure_grid))-pressure_grid(size(pressure_grid)-1))
        value = values(size(values)) + fraction * &
                (values(size(values))-values(size(values)-1))
      else
        ! Otherwise, we perform linear interpolation
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

      if (.not. ieee_is_finite(value) .or. (require_positive .and. value <= 0.0_dp)) then
        if (require_positive) then
          err_ = 'Linear extrapolation of the input '//trim(profile_name)// &
                 ' profile produced a non-finite or non-positive value'
        else
          err_ = 'Linear extrapolation of the input '//trim(profile_name)// &
                 ' profile produced a non-finite value'
        endif
      endif

    end subroutine

    pure function opposite_signs(a, b) result(opposite)
      real(dp), intent(in) :: a, b
      logical :: opposite

      opposite = (a <= 0.0_dp .and. b >= 0.0_dp) .or. &
                 (a >= 0.0_dp .and. b <= 0.0_dp)

    end function
  end subroutine

  function TOA_at_pressure(self, usol, TOA_pressure, candidate, err) result(top_atmos)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use futils, only: brent_class
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: TOA_pressure !! dynes/cm^2
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err
    real(dp) :: top_atmos !! cm

    integer, parameter :: max_bracket_expansions = 60
    real(dp), parameter :: pressure_residual_tolerance = 1.0e-12_dp
    real(dp) :: bottom, current_top, minimum_top, span, search_distance
    real(dp) :: lower, upper, flower, fupper, fcurrent, fnew, fzero
    real(dp) :: altitude_tolerance
    integer :: expansion, iflag
    type(brent_class) :: root_solver

    bottom = self%var%bottom_atmos
    current_top = self%var%top_atmos
    span = current_top - bottom
    ! Lowest trial model top allowed during downward bracket search.
    minimum_top = bottom + max(100.0_dp*epsilon(1.0_dp)* &
                               max(abs(bottom),span,1.0_dp), tiny(1.0_dp))
    ! Absolute altitude tolerance passed to the Brent root solver.
    altitude_tolerance = max(1.0e-10_dp*span, &
                             100.0_dp*epsilon(1.0_dp)*max(abs(current_top),1.0_dp))
    ! Initial altitude increment used to expand the pressure bracket.
    search_distance = max(0.05_dp*span, maxval(self%var%dz))

    ! Set the residual function in the root solver.
    call root_solver%set_function(altitude_residual)
    fcurrent = altitude_residual(root_solver, current_top)
    if (allocated(err)) then
      err = 'Could not evaluate TOA pressure at the current model top: '//err
      return
    endif
    if (abs(fcurrent) <= pressure_residual_tolerance) then
      ! We are at the right pressure, so we can just return
      top_atmos = current_top
      return
    endif

    ! Now determine the bracket for which to do the solve
    if (fcurrent > 0.0_dp) then
      ! The requested pressure is lower than the current TOA pressure, so the
      ! model top must move upward.
      lower = current_top
      flower = fcurrent
      fnew = flower
      do expansion = 1,max_bracket_expansions
        upper = current_top+search_distance
        fupper = altitude_residual(root_solver, upper)
        if (allocated(err)) then
          err = 'Could not bracket the requested TOA pressure while raising the model top: '//err
          return
        endif
        if (fupper >= fnew) then
          err = 'TOA pressure did not decrease while raising the model top; '// &
                'a monotonic pressure bracket could not be constructed.'
          return
        endif
        if (fupper <= 0.0_dp) exit
        fnew = fupper
        search_distance = 2.0_dp*search_distance
      enddo
      if (fupper > 0.0_dp) then
        err = 'Could not bracket the requested TOA pressure above the current model top.'
        return
      endif
    else
      ! The requested pressure is higher than the current TOA pressure, so the
      ! model top must move downward without crossing the model bottom.
      upper = current_top
      fupper = fcurrent
      fnew = fupper
      do expansion = 1,max_bracket_expansions
        lower = max(current_top-search_distance,minimum_top)
        flower = altitude_residual(root_solver, lower)
        if (allocated(err)) then
          err = 'Could not bracket the requested TOA pressure while lowering the model top: '//err
          return
        endif
        if (flower <= fnew) then
          err = 'TOA pressure did not increase while lowering the model top; '// &
                'a monotonic pressure bracket could not be constructed.'
          return
        endif
        if (flower >= 0.0_dp) exit
        if (lower == minimum_top) then
          err = 'The requested TOA pressure exceeds the maximum reachable pressure '// &
                'above the model bottom.'
          return
        endif
        fnew = flower
        search_distance = 2.0_dp*search_distance
      enddo
      if (flower < 0.0_dp) then
        err = 'Could not bracket the requested TOA pressure below the current model top.'
        return
      endif
    endif

    ! Call the root solver for the determined bracket
    call root_solver%find_zero(lower, upper, altitude_tolerance, top_atmos, &
                               fzero, iflag, flower, fupper)
    if (allocated(err)) then
      err = 'Evaluating the bracketed TOA-pressure solve failed: '//err
      return
    endif
    if (iflag /= 0 .or. .not.ieee_is_finite(top_atmos) .or. &
        top_atmos <= bottom .or. .not.ieee_is_finite(fzero)) then
      err = 'The bracketed TOA-pressure solve failed.'
      return
    endif

    ! Successful if reached the end of the function.

  contains
    function altitude_residual(me, altitude) result(residual)
      class(brent_class), intent(inout) :: me
      real(dp), intent(in) :: altitude
      real(dp) :: residual

      real(dp) :: pressure

      pressure = pressure_at_TOA(self, usol, altitude, candidate, err)
      if (allocated(err)) then
        residual = ieee_value(0.0_dp,ieee_quiet_nan)
        return
      endif
      if (.not.ieee_is_finite(pressure) .or. pressure <= 0.0_dp) then
        err = 'Candidate TOA pressure was not finite and positive.'
        residual = ieee_value(0.0_dp,ieee_quiet_nan)
        return
      endif
      residual = log10(pressure)-log10(TOA_pressure)
    end function
  end function

  function pressure_at_TOA(self, usol, top_atmos_new, candidate, err) result(TOA_pressure)
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: top_atmos_new !! cm
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err
    real(dp) :: TOA_pressure !! dynes/cm^2

    call build_vertical_grid_candidate(self, usol, top_atmos_new, candidate, err)
    if (allocated(err)) return

    TOA_pressure = candidate%pressure(self%var%nz)

  end function

  subroutine build_vertical_grid_candidate(self, usol, top_atmos_new, candidate, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: top_atmos_new !! cm
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    if (.not.ieee_is_finite(top_atmos_new) .or. top_atmos_new <= var%bottom_atmos) then
      err = 'The candidate TOA altitude must be finite and above the model bottom.'
      return
    endif
    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz) then
      err = 'The candidate regrid state has the wrong dimensions.'
      return
    endif
    if (.not.all(ieee_is_finite(usol))) then
      err = 'The candidate regrid state must be finite.'
      return
    endif

    if (var%press_temp_edd_profile%enabled) then
      call build_vertical_grid_candidate_p(self, usol, top_atmos_new, candidate, err)
    else
      call build_vertical_grid_candidate_z(self, usol, top_atmos_new, candidate, err)
    endif

  end subroutine

  subroutine initialize_candidate_grid(self, top_atmos_new, candidate, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use photochem_eqns, only: vertical_grid, gravity
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: top_atmos_new !! cm
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    candidate%top_atmos = top_atmos_new
    candidate%trop_alt = var%trop_alt

    call vertical_grid(var%bottom_atmos, top_atmos_new, var%nz, &
                       candidate%z, candidate%dz)
    call gravity(dat%planet_radius, dat%planet_mass, var%nz, &
                 candidate%z, candidate%grav)
    if (.not.all(ieee_is_finite(candidate%z)) .or. &
        .not.all(ieee_is_finite(candidate%dz)) .or. &
        .not.all(ieee_is_finite(candidate%grav)) .or. &
        any(candidate%dz <= 0.0_dp) .or. any(candidate%grav <= 0.0_dp)) then
      err = 'The candidate TOA altitude produced invalid grid geometry or gravity.'
      return
    endif

  end subroutine

  subroutine seed_candidate_temperature_edd(self, candidate, err)
    use futils, only: interp
    class(EvoAtmosphere), target, intent(in) :: self
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    type(PhotochemVars), pointer :: var
    integer :: ierr

    var => self%var

    ! Map the altitude-based profiles. This is the final profile mapping for
    ! the altitude path and an initial guess for the pressure-profile path.
    call interp(candidate%z, var%z, var%temperature, candidate%temperature, ierr=ierr)
    if (ierr /= 0) then
      err = 'Temperature interpolation failed while constructing a candidate vertical grid.'
      return
    endif
    call interp(candidate%z, var%z, log10(max(var%edd,1.0e-40_dp)), &
                candidate%edd, ierr=ierr)
    if (ierr /= 0) then
      err = 'Eddy-diffusion interpolation failed while constructing a candidate vertical grid.'
      return
    endif
    candidate%edd = 10.0_dp**candidate%edd
    call validate_candidate_temperature_edd(candidate, err)
    if (allocated(err)) return

  end subroutine

  subroutine validate_candidate_temperature_edd(candidate, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(VerticalGridCandidate), intent(in) :: candidate
    character(:), allocatable, intent(out) :: err

    if (.not.all(ieee_is_finite(candidate%temperature)) .or. &
        any(candidate%temperature <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(candidate%edd)) .or. &
        any(candidate%edd <= 0.0_dp)) then
      err = 'The candidate vertical grid produced invalid temperature or eddy diffusion.'
      return
    endif

  end subroutine

  subroutine prepare_candidate_composition(self, usol, candidate, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_const, only: small_real
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol(:,:)
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    real(dp) :: gas_mix_total
    integer :: i, j, ierr

    dat => self%dat
    var => self%var

    ! Convert the input atmosphere to mixing ratios and total gas density.
    do j = 1,var%nz
      candidate%density(j) = sum(usol(dat%ng_1:,j))
      if (.not.ieee_is_finite(candidate%density(j)) .or. &
          candidate%density(j) <= 0.0_dp) then
        err = 'Gas density must be finite and positive to construct a candidate vertical grid.'
        return
      endif
      candidate%mix(:,j) = usol(:,j)/candidate%density(j)
    enddo

    ! Interpolate mixing ratios in log-space. Gas mixing ratios are
    ! normalized below; particle abundances remain ratios to total gas.
    do i = 1,dat%nq
      call interp(candidate%z, var%z, &
                  log10(max(candidate%mix(i,:),small_real)), &
                  candidate%mix_new(i,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Mixing-ratio interpolation failed while constructing a candidate vertical grid.'
        return
      endif
    enddo
    candidate%mix_new = 10.0_dp**candidate%mix_new
    if (.not.all(ieee_is_finite(candidate%mix_new))) then
      err = 'The candidate vertical grid produced nonfinite mixing ratios.'
      return
    endif

    ! Normalize gas mixing ratios so they sum to one. Particle abundances
    ! remain ratios relative to total gas density.
    do j = 1,var%nz
      gas_mix_total = sum(candidate%mix_new(dat%ng_1:dat%nq,j))
      if (.not.ieee_is_finite(gas_mix_total) .or. gas_mix_total <= 0.0_dp) then
        err = 'The candidate vertical grid produced invalid gas mixing ratios.'
        return
      endif
      candidate%mix_new(dat%ng_1:dat%nq,j) = &
          candidate%mix_new(dat%ng_1:dat%nq,j)/gas_mix_total
    enddo

  end subroutine

  subroutine extend_candidate_density(self, candidate, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_const, only: k_boltz, N_avo
    class(EvoAtmosphere), target, intent(in) :: self
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    real(dp) :: pressure_previous, temperature_previous
    real(dp) :: delta_z, mubar
    integer :: j, ierr, first_extended

    dat => self%dat
    var => self%var

    ! Interpolate density inside the old grid. Above the old top, replace the
    ! interpolation with a hydrostatic continuation using candidate T and z.
    call interp(candidate%z, var%z, log10(candidate%density), &
                candidate%density_new, ierr=ierr)
    if (ierr /= 0) then
      err = 'Gas-density interpolation failed while constructing a candidate vertical grid.'
      return
    endif
    candidate%density_new = 10.0_dp**candidate%density_new
    if (.not.all(ieee_is_finite(candidate%density_new)) .or. &
        any(candidate%density_new <= 0.0_dp)) then
      err = 'The candidate vertical grid produced invalid gas density.'
      return
    endif

    ! If the model was extended upward, replace the extrapolation above the
    ! old top with a hydrostatic continuation.
    first_extended = var%nz + 1
    do j = 1,var%nz
      if (candidate%z(j) > var%z(var%nz)) then
        first_extended = j
        exit
      endif
    enddo
    if (first_extended > var%nz) return

    pressure_previous = candidate%density(var%nz)*k_boltz*var%temperature(var%nz)
    temperature_previous = var%temperature(var%nz)
    delta_z = candidate%z(first_extended)-var%z(var%nz)
    do j = first_extended,var%nz
      if (j > first_extended) delta_z = candidate%z(j)-candidate%z(j-1)
      mubar = sum(dat%species_mass(dat%ng_1:dat%nq)* &
                  candidate%mix_new(dat%ng_1:dat%nq,j))
      pressure_previous = pressure_previous*exp( &
          -(mubar*candidate%grav(j)*delta_z)/ &
           (N_avo*k_boltz*0.5_dp* &
            (temperature_previous+candidate%temperature(j))))
      candidate%density_new(j) = pressure_previous/(k_boltz*candidate%temperature(j))
      temperature_previous = candidate%temperature(j)
    enddo

  end subroutine

  subroutine map_candidate_particle_radii(self, candidate, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_const, only: small_real
    class(EvoAtmosphere), target, intent(in) :: self
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    integer :: i, ierr

    dat => self%dat
    var => self%var

    if (.not.dat%there_are_particles) return

    do i = 1,dat%npq
      call interp(candidate%z, var%z, &
                  log10(max(var%particle_radius(i,:),small_real)), &
                  candidate%particle_radius(i,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Particle-radius interpolation failed while constructing a candidate vertical grid.'
        return
      endif
    enddo
    candidate%particle_radius = 10.0_dp**candidate%particle_radius
    if (.not.all(ieee_is_finite(candidate%particle_radius)) .or. &
        any(candidate%particle_radius <= 0.0_dp)) then
      err = 'The candidate vertical grid produced invalid particle radii.'
      return
    endif

  end subroutine

  subroutine assemble_candidate_atmosphere(self, candidate, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use photochem_const, only: k_boltz
    class(EvoAtmosphere), target, intent(in) :: self
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    integer :: j

    dat => self%dat
    var => self%var

    ! Construct the candidate atmospheric state and enforce the lower BCs.
    do j = 1,var%nz
      candidate%usol(:,j) = candidate%mix_new(:,j)*candidate%density_new(j)
    enddo
    call self%apply_lower_boundary_conditions(candidate%temperature(1), &
                                              candidate%usol(:,1), err)
    if (allocated(err)) return

    do j = 1,var%nz
      candidate%pressure(j) = sum(candidate%usol(dat%ng_1:,j))* &
                              k_boltz*candidate%temperature(j)
    enddo
    if (.not.all(ieee_is_finite(candidate%usol)) .or. &
        .not.all(ieee_is_finite(candidate%pressure)) .or. &
        any(candidate%pressure <= 0.0_dp)) then
      err = 'The candidate vertical grid produced invalid composition or pressure.'
      return
    endif

  end subroutine

  subroutine build_vertical_grid_candidate_z(self, usol, top_atmos_new, candidate, err)
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: top_atmos_new !! cm
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    call initialize_candidate_grid(self, top_atmos_new, candidate, err)
    if (allocated(err)) return

    call seed_candidate_temperature_edd(self, candidate, err)
    if (allocated(err)) return

    call prepare_candidate_composition(self, usol, candidate, err)
    if (allocated(err)) return

    call extend_candidate_density(self, candidate, err)
    if (allocated(err)) return

    call map_candidate_particle_radii(self, candidate, err)
    if (allocated(err)) return

    call assemble_candidate_atmosphere(self, candidate, err)

  end subroutine

  subroutine build_vertical_grid_candidate_p(self, usol, top_atmos_new, candidate, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: top_atmos_new !! cm
    type(VerticalGridCandidate), intent(inout) :: candidate
    character(:), allocatable, intent(out) :: err

    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
    real(dp), parameter :: persistent_tolerance = 1.0e-10_dp
    integer, parameter :: persistent_max_iterations = 50
    real(dp) :: temperature_change
    integer :: iteration
    logical :: converged

    var => self%var
    wrk => self%wrk

    ! Build an initial candidate using the shared altitude-based seed. The
    ! pressure-profile mapping replaces that seed during the fixed-point loop.
    call initialize_candidate_grid(self, top_atmos_new, candidate, err)
    if (allocated(err)) return
    call seed_candidate_temperature_edd(self, candidate, err)
    if (allocated(err)) return
    call prepare_candidate_composition(self, usol, candidate, err)
    if (allocated(err)) return
    call extend_candidate_density(self, candidate, err)
    if (allocated(err)) return
    call assemble_candidate_atmosphere(self, candidate, err)
    if (allocated(err)) return

    ! Use the current hydrostatic pressure as the initial reference for the
    ! pressure-profile solve. Each later iteration uses the pressure profile
    ! produced by the preceding mapping pass.
    candidate%pressure_reference = wrk%pressure_hydro
    converged = .false.
    do iteration = 1,persistent_max_iterations
      candidate%temperature_reference = candidate%temperature
      call map_press_temp_edd( &
        self, &
        candidate%usol, &
        var%press_temp_edd_profile%pressure, &
        var%press_temp_edd_profile%temperature, &
        var%press_temp_edd_profile%edd, &
        trop_p=var%press_temp_edd_profile%trop_p, &
        hydro_pressure=var%press_temp_edd_profile%hydro_pressure, &
        grid_z=candidate%z, &
        grid_dz=candidate%dz, &
        grid_grav=candidate%grav, &
        temperature_reference=candidate%temperature_reference, &
        pressure_reference=candidate%pressure_reference, &
        T_grid=candidate%temperature, &
        edd_grid=candidate%edd, &
        log10P_grid=candidate%log10P, &
        trop_alt=candidate%trop_alt, &
        err=err &
      )
      if (allocated(err)) return

      call validate_candidate_temperature_edd(candidate, err)
      if (allocated(err)) return
      if (.not.all(ieee_is_finite(candidate%log10P))) then
        err = 'The pressure-profile mapping produced nonfinite log-pressure values.'
        return
      endif

      temperature_change = maxval(abs(candidate%temperature - &
                                      candidate%temperature_reference)/ &
                                  max(candidate%temperature,1.0_dp))
      candidate%pressure_reference = 10.0_dp**candidate%log10P
      if (.not.all(ieee_is_finite(candidate%pressure_reference)) .or. &
          any(candidate%pressure_reference <= 0.0_dp)) then
        err = 'The pressure-profile mapping produced invalid pressure values.'
        return
      endif

      if (temperature_change <= persistent_tolerance) then
        converged = .true.
        exit
      endif
      if (iteration == persistent_max_iterations) exit

      ! Rebuild the state used by the next mapping pass from the updated
      ! pressure-profile temperature.
      call extend_candidate_density(self, candidate, err)
      if (allocated(err)) return
      call assemble_candidate_atmosphere(self, candidate, err)
      if (allocated(err)) return
    enddo

    if (.not.converged) then
      err = 'The persistent pressure profile did not converge on the candidate vertical grid.'
      return
    endif

    ! Rebuild the final candidate state using the converged temperature.
    call extend_candidate_density(self, candidate, err)
    if (allocated(err)) return
    call map_candidate_particle_radii(self, candidate, err)
    if (allocated(err)) return
    call assemble_candidate_atmosphere(self, candidate, err)

  end subroutine

  subroutine VerticalGridCandidate_allocate(candidate, np, nq, nz)
    type(VerticalGridCandidate), intent(inout) :: candidate
    integer, intent(in) :: np, nq, nz

    allocate(candidate%usol(nq,nz))
    allocate(candidate%pressure(nz))
    allocate(candidate%z(nz))
    allocate(candidate%dz(nz))
    allocate(candidate%grav(nz))
    allocate(candidate%temperature(nz))
    allocate(candidate%edd(nz))
    allocate(candidate%particle_radius(np, nz))
    allocate(candidate%mix(nq, nz))
    allocate(candidate%mix_new(nq, nz))
    allocate(candidate%density(nz))
    allocate(candidate%density_new(nz))
    allocate(candidate%temperature_reference(nz))
    allocate(candidate%pressure_reference(nz))
    allocate(candidate%log10P(nz))

  end subroutine

  module subroutine update_vertical_grid(self, TOA_alt, TOA_pressure, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use photochem_input, only: finalize_atmosphere_initialization
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), optional, intent(in) :: TOA_alt !! cm
    real(dp), optional, intent(in) :: TOA_pressure !! dynes/cm^2
    character(:), allocatable, intent(out) :: err

    real(dp) :: top_atmos_new
    type(VerticalGridCandidate) :: candidate, candidate_copy

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
    character(:), allocatable :: original_err

    dat => self%dat
    var => self%var
    wrk => self%wrk

    ! Check inputs
    call self%require_atmosphere_initialized('update_vertical_grid', err)
    if (allocated(err)) return

    if (present(TOA_alt) .and. present(TOA_pressure)) then
      err = 'Both "TOA_alt" and "TOA_pressure" can not be specified'
      return
    endif
    if (.not.present(TOA_alt) .and. .not.present(TOA_pressure)) then
      err = 'Either "TOA_alt" and "TOA_pressure" must be specified'
      return
    endif

    if (present(TOA_alt)) then
      if (.not.ieee_is_finite(TOA_alt) .or. TOA_alt <= var%bottom_atmos) then
        err = '"TOA_alt" must be finite and greater than the model bottom.'
        return
      endif
    endif

    if (present(TOA_pressure)) then
      if (.not.ieee_is_finite(TOA_pressure) .or. TOA_pressure <= 0.0_dp) then
        err = '"TOA_pressure" must be finite and positive.'
        return
      endif
    endif

    call VerticalGridCandidate_allocate(candidate, dat%np, dat%nq, var%nz)
    call VerticalGridCandidate_allocate(candidate_copy, dat%np, dat%nq, var%nz)

    ! Snapshot the committed state before refreshing any derived work arrays.
    ! The hydrostatic pressure used by candidate construction must be derived
    ! from the canonical atmospheric state, not a possibly stale work cache.
    call copy_var_to_candidate(var, candidate_copy)
    candidate_copy%usol = wrk%usol

    ! Refresh all atmospheric work state from the committed composition. When
    ! enabled, the persistent profile is applied before pressure and
    ! hydrostatic quantities are recomputed.
    call self%prep_atm_evo_gas(candidate_copy%usol, wrk%usol, &
                               wrk%molecules_per_particle, wrk%pressure, &
                               wrk%density, wrk%mix, wrk%mubar, &
                               wrk%pressure_hydro, wrk%density_hydro, &
                               var%press_temp_edd_profile%enabled, err)
    if (allocated(err)) then
      original_err = err
      call rollback_vertical_grid(self, candidate_copy, original_err, err)
      return
    endif

    if (present(TOA_alt)) then
      top_atmos_new = TOA_alt
    else
      top_atmos_new = TOA_at_pressure(self, self%wrk%usol, TOA_pressure, candidate, err)
      if (allocated(err)) then
        original_err = err
        call rollback_vertical_grid(self, candidate_copy, original_err, err)
        return
      endif
    endif

    ! Compute properties associated with new TOA
    call build_vertical_grid_candidate(self, self%wrk%usol, top_atmos_new, candidate, err)
    if (allocated(err)) then
      original_err = err
      call rollback_vertical_grid(self, candidate_copy, original_err, err)
      return
    endif

    ! Now commit
    call copy_candidate_to_var(candidate, var)
    call finalize_atmosphere_initialization(dat, var, err)
    if (allocated(err)) then
      original_err = err
      call rollback_vertical_grid(self, candidate_copy, original_err, err)
      return
    endif
    call self%prep_atmosphere(candidate%usol, err)
    if (allocated(err)) then
      original_err = err
      call rollback_vertical_grid(self, candidate_copy, original_err, err)
      return
    endif
    call self%destroy_stepper(err)
    if (allocated(err)) then
      original_err = err
      call rollback_vertical_grid(self, candidate_copy, original_err, err)
      return
    endif

  end subroutine

  subroutine copy_candidate_to_var(candidate, var)
    type(VerticalGridCandidate), intent(in) :: candidate
    type(PhotochemVars), intent(inout) :: var

    var%trop_alt = candidate%trop_alt
    var%top_atmos = candidate%top_atmos
    var%z = candidate%z
    var%dz = candidate%dz
    var%grav = candidate%grav
    var%temperature = candidate%temperature
    var%edd = candidate%edd
    var%particle_radius = candidate%particle_radius

  end subroutine

  subroutine copy_var_to_candidate(var, candidate)
    type(PhotochemVars), intent(in) :: var
    type(VerticalGridCandidate), intent(inout) :: candidate

    candidate%trop_alt = var%trop_alt
    candidate%top_atmos = var%top_atmos
    candidate%z = var%z
    candidate%dz = var%dz
    candidate%grav = var%grav
    candidate%temperature = var%temperature
    candidate%edd = var%edd
    candidate%particle_radius = var%particle_radius

  end subroutine

  subroutine rollback_vertical_grid(self, candidate, original_err, err)
    use photochem_input, only: finalize_atmosphere_initialization
    class(EvoAtmosphere), target, intent(inout) :: self
    type(VerticalGridCandidate), intent(in) :: candidate
    character(:), allocatable, intent(in) :: original_err
    character(:), allocatable, intent(out) :: err

    character(:), allocatable :: rollback_err

    call copy_candidate_to_var(candidate, self%var)

    call finalize_atmosphere_initialization(self%dat, self%var, rollback_err)
    if (.not.allocated(rollback_err)) then
      call self%prep_atmosphere(candidate%usol, rollback_err)
    endif

    if (allocated(rollback_err)) then
      err = original_err//' Rollback failed: '//rollback_err
    else
      err = original_err
    endif

  end subroutine

end submodule
