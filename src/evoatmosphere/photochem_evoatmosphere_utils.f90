
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
    if (allocated(err)) return
    
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

    real(dp) :: T_new(self%var%nz), edd_new(self%var%nz)
    real(dp) :: log10P_wrk(self%var%nz), trop_alt
    real(dp) :: edd_save(self%var%nz)
    real(dp) :: trop_p_ = 0.0_dp
    logical :: has_trop_p

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

    has_trop_p = present(trop_p)
    if (has_trop_p) trop_p_ = trop_p

    ! First compute the mapping without changing model state. This kernel is
    ! also suitable for applying a persistent pressure-based profile to an
    ! arbitrary trial composition during a future RHS evaluation.
    call map_press_temp_edd(self, self%wrk%usol, P, T, edd, &
                            trop_p=trop_p_, has_trop_p=has_trop_p, &
                            hydro_pressure=hydro_pressure, &
                            grid_z=self%var%z, grid_dz=self%var%dz, &
                            grid_grav=self%var%grav, &
                            temperature_reference=self%var%temperature, &
                            pressure_reference=self%wrk%pressure_hydro, &
                            T_new=T_new, edd_new=edd_new, &
                            log10P_wrk=log10P_wrk, trop_alt=trop_alt, err=err)
    if (allocated(err)) return

    ! Commit point: the mapping above only reads self. Updating Kzz and
    ! calling set_temperature below are the only operations that mutate it.
    edd_save = self%var%edd
    self%var%edd = edd_new
    if (present(trop_p)) then
      call self%set_temperature(T_new, trop_alt, err)
    else
      call self%set_temperature(T_new, err=err)
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
    self%var%press_temp_edd_profile%has_trop_p = present(trop_p)
    if (present(trop_p)) then
      self%var%press_temp_edd_profile%trop_p = trop_p
    else
      self%var%press_temp_edd_profile%trop_p = 0.0_dp
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
      self%var = var_save
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
    var%press_temp_edd_profile%has_trop_p = .false.
    var%press_temp_edd_profile%trop_p = 0.0_dp
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

    real(dp) :: T_new(self%var%nz), edd_new(self%var%nz)
    real(dp) :: log10P_wrk(self%var%nz), trop_alt

    if (.not. self%var%press_temp_edd_profile%enabled) return

    call map_press_temp_edd(self, usol_in, &
         self%var%press_temp_edd_profile%pressure, &
         self%var%press_temp_edd_profile%temperature, &
         self%var%press_temp_edd_profile%edd, &
         trop_p=self%var%press_temp_edd_profile%trop_p, &
         has_trop_p=self%var%press_temp_edd_profile%has_trop_p, &
         hydro_pressure=self%var%press_temp_edd_profile%hydro_pressure, &
         grid_z=self%var%z, grid_dz=self%var%dz, grid_grav=self%var%grav, &
         temperature_reference=self%var%temperature, &
         pressure_reference=self%wrk%pressure_hydro, &
         T_new=T_new, edd_new=edd_new, log10P_wrk=log10P_wrk, &
         trop_alt=trop_alt, err=err)
    if (allocated(err)) return

    ! Commit only after the profile has mapped successfully and the
    ! tropopause has been validated. Do not call set_temperature here: that
    ! routine calls prep_atmosphere and would recurse back into this routine.
    self%var%temperature = T_new
    self%var%edd = edd_new
    call refresh_temperature_dependent_state(self%dat, self%var, trop_alt=trop_alt, err=err)

  end subroutine

  module subroutine map_press_temp_edd(self, usol_in, P, T, edd, trop_p, has_trop_p, hydro_pressure, &
                                       grid_z, grid_dz, grid_grav, temperature_reference, &
                                       pressure_reference, &
                                       T_new, edd_new, log10P_wrk, trop_alt, err)
    use futils, only: interp, brent_class
    use ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use photochem_const, only: small_real, k_boltz, N_avo
    use photochem_enum, only: DensityBC, PressureBC
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol_in(:,:)
    real(dp), intent(in) :: P(:)
    real(dp), intent(in) :: T(:)
    real(dp), intent(in) :: edd(:)
    real(dp), intent(in) :: trop_p
    logical, intent(in) :: has_trop_p
    logical, optional, intent(in) :: hydro_pressure
    real(dp), intent(in) :: grid_z(:), grid_dz(:), grid_grav(:)
    real(dp), intent(in) :: temperature_reference(:), pressure_reference(:)
    real(dp), intent(out) :: T_new(:)
    real(dp), intent(out) :: edd_new(:)
    real(dp), intent(out) :: log10P_wrk(:)
    real(dp), intent(out) :: trop_alt
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: log10P_in(:), T_in(:), log10edd_in(:)
    real(dp), allocatable :: P_wrk(:), log10edd_new(:)
    real(dp), allocatable :: usol_base(:,:), density_base(:), mubar_base(:)
    real(dp) :: xzero, Psurf_initial, Psurf_final
    real(dp) :: log10P_previous, temperature_previous
    real(dp) :: trop_alt_array(1)
    logical :: hydro_pressure_
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
    if (size(usol_in,1) /= dat%nq .or. size(usol_in,2) /= var%nz) then
      err = '"usol_in" has the wrong dimensions'
      return
    endif
    if (size(T_new) /= var%nz .or. size(edd_new) /= var%nz .or. &
        size(log10P_wrk) /= var%nz) then
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
    if (has_trop_p) then
      if (.not. ieee_is_finite(trop_p) .or. trop_p <= 0.0_dp) then
        err = '"trop_p" must be finite and positive'
        return
      endif
    endif
    if (dat%gas_rainout .and. .not.has_trop_p) then
      err = '"trop_p" is a required input.'
      return
    endif

    ! optional arguments
    if (present(hydro_pressure)) then
      hydro_pressure_ = hydro_pressure
    else
      hydro_pressure_ = .true. ! default is True
    endif

    ! Allocate work space
    allocate(P_wrk(var%nz), log10edd_new(var%nz))
    allocate(usol_base(dat%nq,var%nz), density_base(var%nz), mubar_base(var%nz))

    ! Copy and clip the current state in the same way as prep_atm_evo_gas.
    call clip_usol(usol_in, usol_base)

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
      log10P_wrk(i) = xzero
      P_wrk(i) = 10.0_dp**xzero
      call evaluate_profile(xzero, log10P_in, T_in, T_new(i), 'temperature', .true., err)
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
          if (hydro_pressure_) then
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

        call bottom_column_state(T_new(1), density_final, mubar_final, Psurf_final, err)
        if (allocated(err)) return
        if (hydro_pressure_ .and. P_wrk(1) >= Psurf_final) then
          err = 'The solved bottom-layer pressure is not below the surface pressure'
          return
        endif
      end block; endif

      ! Save the previous solve. This is needed for initial guesses and for the residuals on
      ! subsequent layers.
      log10P_previous = xzero
      temperature_previous = T_new(i)
    enddo

    ! Ensure eddy diffusion is valid
    if (.not. all(ieee_is_finite(log10edd_new))) then
      err = 'The mapped eddy-diffusion profile is not finite'
      return
    endif
    edd_new = 10.0_dp**log10edd_new
    if (.not. all(ieee_is_finite(edd_new)) .or. any(edd_new <= 0.0_dp)) then
      err = 'The mapped eddy-diffusion profile is not finite and positive'
      return
    endif

    ! Compute the tropopause altitude
    trop_alt = 0.0_dp
    if (has_trop_p) then
      call interp([log10(trop_p)], log10P_wrk(var%nz:1:-1), &
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
      if (hydro_pressure_) then
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

      if (hydro_pressure_) then
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

      real(dp) :: usol_bottom(dat%nq), column_mass
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

  function TOA_at_pressure(self, usol, TOA_pressure, err) result(top_atmos)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use futils, only: brent_class
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: TOA_pressure !! dynes/cm^2
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

      pressure = pressure_at_TOA(self, usol, altitude, err)
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

  function pressure_at_TOA(self, usol, top_atmos_new, err) result(TOA_pressure)
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: top_atmos_new !! cm
    character(:), allocatable, intent(out) :: err
    real(dp) :: TOA_pressure !! dynes/cm^2

    type(VerticalGridCandidate) :: candidate

    call build_vertical_grid_candidate(self, usol, top_atmos_new, candidate, err)
    if (allocated(err)) return

    TOA_pressure = candidate%pressure(self%var%nz)

  end function

  subroutine build_vertical_grid_candidate(self, usol, top_atmos_new, candidate, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_eqns, only: vertical_grid, gravity
    use photochem_const, only: small_real, k_boltz, N_avo
    class(EvoAtmosphere), target, intent(in) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: top_atmos_new !! cm
    type(VerticalGridCandidate), intent(out) :: candidate
    character(:), allocatable, intent(out) :: err

    real(dp), parameter :: persistent_tolerance = 1.0e-10_dp
    integer, parameter :: persistent_max_iterations = 50
    real(dp) :: gas_mix_total, pressure_previous, temperature_previous
    real(dp) :: delta_z, mubar
    real(dp), allocatable :: mix(:,:), mix_new(:,:)
    real(dp), allocatable :: density(:), density_new(:)
    real(dp), allocatable :: temperature_mapped(:), edd_mapped(:), log10P_mapped(:)
    real(dp), allocatable :: pressure_reference_candidate(:)
    real(dp) :: trop_alt, temperature_change
    integer :: i, j, ierr, iteration, first_extended

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

    candidate%top_atmos = top_atmos_new
    candidate%surface_pressure = var%surface_pressure
    candidate%trop_alt = var%trop_alt
    candidate%trop_ind = var%trop_ind
    allocate(candidate%usol(dat%nq,var%nz))
    allocate(candidate%z(var%nz),candidate%dz(var%nz),candidate%grav(var%nz))
    allocate(candidate%temperature(var%nz),candidate%edd(var%nz))
    allocate(candidate%particle_radius(dat%np,var%nz))
    allocate(candidate%pressure(var%nz))
    allocate(mix(dat%nq,var%nz), mix_new(dat%nq,var%nz))
    allocate(density(var%nz),density_new(var%nz))

    ! Remake the vertical grid and gravity
    call vertical_grid(var%bottom_atmos, top_atmos_new, &
                      var%nz, candidate%z, candidate%dz)
    call gravity(dat%planet_radius, dat%planet_mass, &
                var%nz, candidate%z, candidate%grav)
    if (.not.all(ieee_is_finite(candidate%z)) .or. &
        .not.all(ieee_is_finite(candidate%dz)) .or. &
        .not.all(ieee_is_finite(candidate%grav)) .or. &
        any(candidate%dz <= 0.0_dp) .or. any(candidate%grav <= 0.0_dp)) then
      err = 'The candidate TOA altitude produced invalid grid geometry or gravity.'
      return
    endif

    ! Temperature
    call interp(candidate%z, var%z, var%temperature, candidate%temperature, ierr=ierr)
    if (ierr /= 0) then
      err = 'Temperature interpolation failed while constructing a candidate vertical grid.'
      return
    endif

    ! Eddy diffusion
    call interp(candidate%z, var%z, log10(max(var%edd,1.0e-40_dp)), &
                candidate%edd, ierr=ierr)
    if (ierr /= 0) then
      err = 'Eddy-diffusion interpolation failed while constructing a candidate vertical grid.'
      return
    endif
    candidate%edd = 10.0_dp**candidate%edd
    if (.not.all(ieee_is_finite(candidate%temperature)) .or. &
        any(candidate%temperature <= 0.0_dp) .or. &
        .not.all(ieee_is_finite(candidate%edd)) .or. any(candidate%edd <= 0.0_dp)) then
      err = 'The candidate vertical grid produced invalid temperature or eddy diffusion.'
      return
    endif

    ! determine mixing ratios and density
    do j = 1,var%nz
      density(j) = sum(usol(dat%ng_1:,j))
      if (.not.ieee_is_finite(density(j)) .or. density(j) <= 0.0_dp) then
        err = 'Gas density must be finite and positive to construct a candidate vertical grid.'
        return
      endif
      mix(:,j) = usol(:,j)/density(j) ! mixing ratios
    enddo
    ! Interpolate mixing ratios, with constant extrapolation
    do i = 1,dat%nq
      call interp(candidate%z, var%z, log10(max(mix(i,:),small_real)), &
                  mix_new(i,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Mixing-ratio interpolation failed while constructing a candidate vertical grid.'
        return
      endif
    enddo
    mix_new = 10.0_dp**mix_new
    if (.not.all(ieee_is_finite(mix_new))) then
      err = 'The candidate vertical grid produced nonfinite mixing ratios.'
      return
    endif

    ! Gas mixing ratios define the mean molecular weight used by the
    ! hydrostatic extension and must sum to exactly one. Particle abundances
    ! remain ratios relative to total gas density. Beyond the old top model
    ! center, both are held at their normalized old-top values explicitly.
    do j = 1,var%nz
      gas_mix_total = sum(mix_new(dat%ng_1:dat%nq,j))
      if (.not.ieee_is_finite(gas_mix_total) .or. gas_mix_total <= 0.0_dp) then
        err = 'The candidate vertical grid produced invalid gas mixing ratios.'
        return
      endif
      mix_new(dat%ng_1:dat%nq,j) = mix_new(dat%ng_1:dat%nq,j)/gas_mix_total
    enddo
    first_extended = var%nz + 1
    do j = 1,var%nz
      if (candidate%z(j) > var%z(var%nz)) then
        first_extended = j
        exit
      endif
    enddo
    if (first_extended <= var%nz) then
      gas_mix_total = sum(mix(dat%ng_1:dat%nq,var%nz))
      mix_new(dat%ng_1:dat%nq,first_extended:) = spread( &
          mix(dat%ng_1:dat%nq,var%nz)/gas_mix_total, 2, var%nz-first_extended+1)
      if (dat%npq > 0) then
        mix_new(1:dat%npq,first_extended:) = spread( &
            max(mix(1:dat%npq,var%nz),small_real), 2, var%nz-first_extended+1)
      endif
    endif

    ! Interpolate density inside the old grid. Values above the old top model
    ! center are replaced below by a hydrostatic continuation.
    call interp(candidate%z, var%z, log10(density), density_new, ierr=ierr)
    if (ierr /= 0) then
      err = 'Gas-density interpolation failed while constructing a candidate vertical grid.'
      return
    endif
    density_new = 10.0_dp**density_new
    if (.not.all(ieee_is_finite(density_new)) .or. any(density_new <= 0.0_dp)) then
      err = 'The candidate vertical grid produced invalid gas density.'
      return
    endif

    ! Particle radii
    if (dat%there_are_particles) then
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
      if (first_extended <= var%nz) then
        candidate%particle_radius(:,first_extended:) = spread( &
            var%particle_radius(:,var%nz), 2, var%nz-first_extended+1)
      endif
      if (.not.all(ieee_is_finite(candidate%particle_radius)) .or. &
          any(candidate%particle_radius <= 0.0_dp)) then
        err = 'The candidate vertical grid produced invalid particle radii.'
        return
      endif
    endif

    ! Reconcile a persistent pressure-based T-Kzz profile with the candidate
    ! composition before committing anything to the live model. Mapping the
    ! profile changes temperature, while the hydrostatic continuation depends
    ! on that temperature, so iterate the two operations to consistency.
    if (var%press_temp_edd_profile%enabled) then
      allocate(temperature_mapped(var%nz), edd_mapped(var%nz), &
               log10P_mapped(var%nz), pressure_reference_candidate(var%nz))
      pressure_reference_candidate = self%wrk%pressure_hydro

      do iteration = 1,persistent_max_iterations
        call extend_density_hydrostatically()
        call fill_candidate_usol()

        call map_press_temp_edd(self, candidate%usol, &
             var%press_temp_edd_profile%pressure, &
             var%press_temp_edd_profile%temperature, &
             var%press_temp_edd_profile%edd, &
             trop_p=var%press_temp_edd_profile%trop_p, &
             has_trop_p=var%press_temp_edd_profile%has_trop_p, &
             hydro_pressure=var%press_temp_edd_profile%hydro_pressure, &
             grid_z=candidate%z, grid_dz=candidate%dz, &
             grid_grav=candidate%grav, &
             temperature_reference=candidate%temperature, &
             pressure_reference=pressure_reference_candidate, &
             T_new=temperature_mapped, edd_new=edd_mapped, &
             log10P_wrk=log10P_mapped, trop_alt=trop_alt, err=err)
        if (allocated(err)) return

        temperature_change = maxval(abs(temperature_mapped-candidate%temperature)/ &
                                    max(temperature_mapped,1.0_dp))
        candidate%temperature = temperature_mapped
        candidate%edd = edd_mapped
        pressure_reference_candidate = 10.0_dp**log10P_mapped
        if (temperature_change <= persistent_tolerance) exit
      enddo
      if (iteration > persistent_max_iterations) then
        err = 'The persistent pressure profile did not converge on the candidate vertical grid.'
        return
      endif
    endif

    ! One final continuation uses the converged persistent temperature, or
    ! the altitude-interpolated temperature when persistence is disabled.
    call extend_density_hydrostatically()
    call fill_candidate_usol()

    ! Account for fixed surface mixing ratios.
    call self%apply_lower_boundary_conditions(candidate%temperature(1), candidate%usol(:,1), err)
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

    if (dat%gas_rainout) then
      if (var%press_temp_edd_profile%enabled) candidate%trop_alt = trop_alt
      candidate%trop_ind = max(minloc(abs(candidate%z-candidate%trop_alt), 1)-1, 1)
      if (candidate%trop_ind < 3) then
        err = 'The candidate vertical grid places the tropopause too low.'
        return
      elseif (candidate%trop_ind > var%nz-2) then
        err = 'The candidate vertical grid places the tropopause too high.'
        return
      endif
    else
      candidate%trop_ind = 1
    endif

  contains

    subroutine fill_candidate_usol()
      integer :: layer

      do layer = 1,var%nz
        candidate%usol(:,layer) = mix_new(:,layer)*density_new(layer)
      enddo
    end subroutine

    subroutine extend_density_hydrostatically()
      integer :: layer

      if (first_extended > var%nz) return

      pressure_previous = density(var%nz)*k_boltz*var%temperature(var%nz)
      temperature_previous = var%temperature(var%nz)
      delta_z = candidate%z(first_extended)-var%z(var%nz)
      do layer = first_extended,var%nz
        if (layer > first_extended) delta_z = candidate%z(layer)-candidate%z(layer-1)
        mubar = sum(dat%species_mass(dat%ng_1:dat%nq)* &
                    mix_new(dat%ng_1:dat%nq,layer))
        pressure_previous = pressure_previous*exp( &
            -(mubar*candidate%grav(layer)*delta_z)/ &
             (N_avo*k_boltz*0.5_dp* &
              (temperature_previous+candidate%temperature(layer))))
        density_new(layer) = pressure_previous/(k_boltz*candidate%temperature(layer))
        temperature_previous = candidate%temperature(layer)
      enddo
    end subroutine

  end subroutine


  module subroutine update_vertical_grid(self, TOA_alt, TOA_pressure, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use photochem_input, only: interp2particlexsdata, refresh_temperature_dependent_state
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), optional, intent(in) :: TOA_alt !! cm
    real(dp), optional, intent(in) :: TOA_pressure !! dynes/cm^2
    character(:), allocatable, intent(out) :: err

    real(dp) :: top_atmos_new
    type(VerticalGridCandidate) :: candidate
    ! The candidate owns proposed grid-dependent arrays.  The work structures
    ! below deliberately remain separate: original_wrk owns the currently
    ! committed state (and any live CVODE resources), while prepared_wrk owns
    ! the replacement state until the commit boundary is crossed.
    type(PhotochemWrkEvo), allocatable :: original_wrk, prepared_wrk

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

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

    if (present(TOA_alt)) then
      top_atmos_new = TOA_alt
    else
      top_atmos_new = TOA_at_pressure(self, self%wrk%usol, TOA_pressure, err)
      if (allocated(err)) return
    endif

    ! Compute properties associated with new TOA
    call build_vertical_grid_candidate(self, self%wrk%usol, top_atmos_new, candidate, err)
    if (allocated(err)) return

    ! Prepare the candidate using a fresh work structure. The original work
    ! structure, including all CVODE-owned pointers, is parked without copying.
    call initialize_vertical_grid_candidate_derived(dat, var, candidate)
    allocate(prepared_wrk)
    call prepared_wrk%init(dat%nsp, dat%np, dat%nq, var%nz, dat%nrT, dat%kj, dat%nw)
    prepared_wrk%tn = self%wrk%tn

    call swap_vertical_grid_candidate(var, candidate)
    call move_alloc(self%wrk, original_wrk)
    call move_alloc(prepared_wrk, self%wrk)

    call interp2particlexsdata(dat, var, err)
    if (.not.allocated(err)) call refresh_temperature_dependent_state(dat, var, err=err)
    if (.not.allocated(err)) then
      ! The persistent profile was reconciled during candidate construction.
      call prep_all_evo_gas(self, candidate%usol, &
                            apply_persistent_profile=.false., err=err)
    endif

    if (allocated(err)) then
      ! Restore the exact committed state. original_wrk still owns the active
      ! stepper, so a failed candidate does not disturb integration.
      call restore_vertical_grid_candidate(self, var, candidate, original_wrk, prepared_wrk)
      return
    endif

    ! Restore the committed state before crossing the ownership boundary.
    ! A successful regrid always invalidates the old stepper before commit.
    call commit_vertical_grid_candidate(self, var, candidate, original_wrk, prepared_wrk, err)

  end subroutine

  subroutine restore_vertical_grid_candidate(self, var, candidate, original_wrk, prepared_wrk)
    ! Return ownership to the pre-regrid state.  This helper is also used on
    ! candidate-preparation failure, so no candidate allocation or work-array
    ! copy can leak into the committed model.
    class(EvoAtmosphere), target, intent(inout) :: self
    type(PhotochemVars), intent(inout) :: var
    type(VerticalGridCandidate), intent(inout) :: candidate
    type(PhotochemWrkEvo), allocatable, intent(inout) :: original_wrk, prepared_wrk

    call move_alloc(self%wrk, prepared_wrk)
    call move_alloc(original_wrk, self%wrk)
    call swap_vertical_grid_candidate(var, candidate)

  end subroutine

  subroutine commit_vertical_grid_candidate(self, var, candidate, original_wrk, prepared_wrk, err)
    ! Commit a fully prepared candidate without ever exposing its partial work
    ! arrays as the live state.  The old CVODE state is destroyed only after
    ! candidate preparation succeeds; if destruction reports an error, the
    ! original state remains installed and the caller receives that error.
    class(EvoAtmosphere), target, intent(inout) :: self
    type(PhotochemVars), intent(inout) :: var
    type(VerticalGridCandidate), intent(inout) :: candidate
    type(PhotochemWrkEvo), allocatable, intent(inout) :: original_wrk, prepared_wrk
    character(:), allocatable, intent(out) :: err

    call restore_vertical_grid_candidate(self, var, candidate, original_wrk, prepared_wrk)
    call self%destroy_stepper(err)
    if (allocated(err)) return

    ! The old work structure is now finalized and can be discarded. Reinstall
    ! the candidate grid and promote prepared_wrk%usol to the sole committed
    ! atmospheric state.
    call move_alloc(self%wrk, original_wrk)
    call swap_vertical_grid_candidate(var, candidate)
    call move_alloc(prepared_wrk, self%wrk)

  end subroutine

  subroutine initialize_vertical_grid_candidate_derived(dat, var, candidate)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(in) :: var
    type(VerticalGridCandidate), intent(inout) :: candidate

    integer :: i

    allocate(candidate%xs_x_qy(var%nz,dat%kj,dat%nw))
    allocate(candidate%particle_xs(dat%np))
    do i = 1,dat%np
      if (dat%part_xs_file(i)%ThereIsData) then
        allocate(candidate%particle_xs(i)%w0(var%nz,dat%nw))
        allocate(candidate%particle_xs(i)%qext(var%nz,dat%nw))
        allocate(candidate%particle_xs(i)%gt(var%nz,dat%nw))
      endif
    enddo
    if (dat%reverse) allocate(candidate%gibbs_energy(var%nz,dat%ng))
    candidate%photon_flux = var%photon_flux

  end subroutine

  subroutine swap_vertical_grid_candidate(var, candidate)
    type(PhotochemVars), intent(inout) :: var
    type(VerticalGridCandidate), intent(inout) :: candidate

    real(dp), allocatable :: tmp_1d(:), tmp_2d(:,:), tmp_3d(:,:,:)
    real(dp) :: tmp_real
    integer :: tmp_integer, i

    tmp_real = var%top_atmos
    var%top_atmos = candidate%top_atmos
    candidate%top_atmos = tmp_real
    tmp_real = var%surface_pressure
    var%surface_pressure = candidate%surface_pressure
    candidate%surface_pressure = tmp_real
    tmp_real = var%trop_alt
    var%trop_alt = candidate%trop_alt
    candidate%trop_alt = tmp_real
    tmp_integer = var%trop_ind
    var%trop_ind = candidate%trop_ind
    candidate%trop_ind = tmp_integer

    call swap_real_1d(var%z, candidate%z)
    call swap_real_1d(var%dz, candidate%dz)
    call swap_real_1d(var%grav, candidate%grav)
    call swap_real_1d(var%temperature, candidate%temperature)
    call swap_real_1d(var%edd, candidate%edd)
    call swap_real_1d(var%photon_flux, candidate%photon_flux)

    call move_alloc(var%particle_radius, tmp_2d)
    call move_alloc(candidate%particle_radius, var%particle_radius)
    call move_alloc(tmp_2d, candidate%particle_radius)
    if (allocated(var%gibbs_energy) .or. allocated(candidate%gibbs_energy)) then
      call move_alloc(var%gibbs_energy, tmp_2d)
      call move_alloc(candidate%gibbs_energy, var%gibbs_energy)
      call move_alloc(tmp_2d, candidate%gibbs_energy)
    endif

    call move_alloc(var%xs_x_qy, tmp_3d)
    call move_alloc(candidate%xs_x_qy, var%xs_x_qy)
    call move_alloc(tmp_3d, candidate%xs_x_qy)
    do i = 1,size(candidate%particle_xs)
      if (allocated(candidate%particle_xs(i)%w0)) then
        call swap_real_2d(var%particle_xs(i)%w0, candidate%particle_xs(i)%w0)
        call swap_real_2d(var%particle_xs(i)%qext, candidate%particle_xs(i)%qext)
        call swap_real_2d(var%particle_xs(i)%gt, candidate%particle_xs(i)%gt)
      endif
    enddo

  contains

    subroutine swap_real_1d(a, b)
      real(dp), allocatable, intent(inout) :: a(:), b(:)

      call move_alloc(a, tmp_1d)
      call move_alloc(b, a)
      call move_alloc(tmp_1d, b)
    end subroutine

    subroutine swap_real_2d(a, b)
      real(dp), allocatable, intent(inout) :: a(:,:), b(:,:)

      call move_alloc(a, tmp_2d)
      call move_alloc(b, a)
      call move_alloc(tmp_2d, b)
    end subroutine

  end subroutine

end submodule
