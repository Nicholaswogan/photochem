
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
    
    if (self%dat%fix_water_in_trop) then
      if (species == "H2O") then
        err = "You can not change the boundary condition for H2O because"// &
              " you have water fixed in the troposphere."
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
      if (self%evolve_climate) then
        err = "You can not set a fixed surface pressure when"// &
              " evolving climate."
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
    if (self%evolve_climate) then
      err = "You can not set the temperature when evolving climate"
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
    
    ! if water is fixed in troposhere or gas rainout, and trop_alt present
    ! then we need to change trop_ind, reallocate some stuff
    ! in wrk, then we will re-prep the atmosphere
    if ((dat%fix_water_in_trop .or. dat%gas_rainout) .and. present(trop_alt)) then
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
    use futils, only: interp
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars

    use photochem_const, only: small_real
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: P(:) !! dynes/cm^2
    real(dp), intent(in) :: T(:) !! K
    real(dp), intent(in) :: edd(:) !! cm^2/s
    real(dp), optional, intent(in) :: trop_p !! dynes/cm^2
    logical, optional, intent(in) :: hydro_pressure
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: log10P_in(:), T_in(:), log10edd_in(:)
    real(dp), allocatable :: P_wrk(:), log10P_wrk(:), T_in_interp(:), T_save(:), log10edd_new(:)
    logical :: hydro_pressure_
    integer :: ierr

    type(MinpackHybrd1Vars) :: mv

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    
    dat => self%dat
    var => self%var
    
    if (self%evolve_climate) then
      err = "You can not set the temperature when evolving climate"
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

    if (size(P) <= 2) then
      err = 'size(P) must be >= 2'
      return
    endif

    if (any(P < 0.0_dp) .or. any(T < 0.0_dp) .or. any(edd < 0.0_dp)) then
      err = 'All elements of "P", "T" and "edd" must be positive'
      return
    endif

    if ((dat%fix_water_in_trop .or. dat%gas_rainout) .and. .not.present(trop_p)) then
      err = '"trop_p" is a required input.'
      return
    endif

    ! optional arguments
    if (present(hydro_pressure)) then
      hydro_pressure_ = hydro_pressure
    else
      hydro_pressure_ = .true. ! default is True
    endif

    ! Work
    allocate(P_wrk(var%nz),log10P_wrk(var%nz))
    allocate(T_in_interp(var%nz), T_save(var%nz))
    allocate(log10edd_new(var%nz))

    ! copy over inputs, and covert to log10 space
    allocate(log10P_in(size(P)),T_in(size(P)),log10edd_in(size(P)))
    log10P_in = log10(P)
    T_in = T
    log10edd_in = log10(edd)

    ! if the P-T-edd profile does not extend to the surface,
    ! then we log-linearly extrapolate to surface
    if (self%var%surface_pressure*1.0e6_dp > P(1)) then; block
      real(dp) :: slope, intercept, P_surf, T_surf, edd_surf

      ! log10 surface pressure in dynes/cm^2
      P_surf = log10(self%var%surface_pressure*1.0e6_dp)

      slope = (T_in(2) - T_in(1))/(log10P_in(2) - log10P_in(1))
      intercept = T_in(1) - slope*log10P_in(1)
      T_surf = slope*P_surf + intercept
      T_surf = max(T_surf, small_real)

      slope = (log10edd_in(2) - log10edd_in(1))/(log10P_in(2) - log10P_in(1))
      intercept = log10edd_in(1) - slope*log10P_in(1)
      edd_surf = slope*P_surf + intercept

      log10P_in = [P_surf, log10P_in]
      T_in = [T_surf, T_in]
      log10edd_in = [edd_surf, log10edd_in]
      
    endblock; endif

    ! Flip order for interpolation purposes
    log10P_in = log10P_in(size(log10P_in):1:-1)
    T_in = T_in(size(log10P_in):1:-1)
    log10edd_in = log10edd_in(size(log10P_in):1:-1)

    ! Do non-linear solve for the T profile that matches
    ! input P-T profile
    mv = MinpackHybrd1Vars(var%nz,tol=1.0e-5_dp)
    mv%x = var%temperature
    call hybrd1(fcn, mv%n, mv%x, mv%fvec, mv%tol, mv%info, mv%wa, mv%lwa)
    if (mv%info == 0 .or. mv%info > 1) then
      err = 'hybrd1 root solve failed in set_press_temp_edd.'
      return
    elseif (mv%info < 0) then
      err = 'hybrd1 root solve failed in set_press_temp_edd: '//err
      return
    endif

    ! set the temperature
    if (present(trop_p)) then; block
      real(dp) :: trop_alt(1)
      real(dp), allocatable :: z_(:)

      allocate(z_(var%nz))
      z_ = var%z
      z_ = z_(var%nz:1:-1)
      
      call interp(1, var%nz, [log10(trop_p)], log10P_wrk, z_, trop_alt, ierr)
      if (ierr /= 0) then
        err = 'Subroutine interp returned an error.'
        return
      endif

      call self%set_temperature(mv%x, trop_alt(1), err)
      if (allocated(err)) return
    endblock; else
      call self%set_temperature(mv%x, err=err)
      if (allocated(err)) return
    endif

    ! finally, interpolate input eddy to new grid
    call interp(var%nz, size(log10P_in), log10P_wrk, log10P_in, log10edd_in, log10edd_new, ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
      return
    endif
    log10edd_new = log10edd_new(var%nz:1:-1)
    var%edd = 10.0_dp**log10edd_new

  contains
    subroutine fcn(n_, x_, fvec_, iflag_)
      integer, intent(in) :: n_
      real(dp), intent(in) :: x_(n_)
      real(dp), intent(out) :: fvec_(n_)
      integer, intent(inout) :: iflag_

      ! x_ is temperature. 
      ! First save the temperature in var.
      T_save = var%temperature
      ! Next, set var%temperature to the x_ temperature, to pass the temperature into `prep_atm_evo_gas`
      var%temperature = x_
      call self%prep_atm_evo_gas(self%wrk%usol, self%wrk%usol, &
                          self%wrk%molecules_per_particle, self%wrk%pressure, &
                          self%wrk%density, self%wrk%mix, self%wrk%mubar, &
                          self%wrk%pressure_hydro, self%wrk%density_hydro, err)
      ! Return var%temperature to its original value, before error checking.
      var%temperature = T_save
      if (allocated(err)) then
        iflag_ = -1
        return
      endif
      ! There are two pressures we could use. The hydrostatic pressure
      ! and the true pressure
      if (hydro_pressure_) then
        P_wrk = self%wrk%pressure_hydro
      else
        P_wrk = self%wrk%pressure
      endif
      log10P_wrk = log10(P_wrk) ! log10
      log10P_wrk = log10P_wrk(size(log10P_wrk):1:-1) ! flip order

      ! Then we have P-T relation. Interpolate input T to this P grid.
      ! Assume constant extrapolation of input T above top of atmosphere.
      ! Code log-linearly extrapolates to surface_pressure
      call interp(var%nz, size(log10P_in), log10P_wrk, log10P_in, T_in, T_in_interp, ierr)
      if (ierr /= 0) then
        err = 'Subroutine interp returned an error.'
        iflag_ = -1
        return
      endif

      ! Flip
      T_in_interp = T_in_interp(size(T_in_interp):1:-1)

      fvec_ = T_in_interp - x_

    end subroutine
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

    if (self%evolve_climate) then
      err = 'You can not call "update_vertical_grid" when evolving climate'
      return
    endif

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
    call self%set_temperature(var%temperature, var%trop_alt, err)
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
    
    if (dat%fix_water_in_trop .or. dat%gas_rainout) then
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