submodule (photochem_input) photochem_input_after_read
  implicit none
  
contains

  module subroutine map_atmosphere_z_to_grid(dat, var, z, temperature, &
                                             edd, surface_pressure, mix, &
                                             particle_radius, &
                                             pressure, density, mubar, usol, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_eqns, only: press_and_den

    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    real(dp), intent(in) :: z(:), temperature(:), edd(:)
    real(dp), intent(in) :: surface_pressure
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    real(dp), intent(out) :: pressure(:), density(:), mubar(:)
    real(dp), intent(out) :: usol(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), parameter :: mixing_ratio_floor = 1.0e-40_dp
    real(dp), allocatable :: gas_mix_normalized(:,:), gas_mix_model(:,:), particle_mix_model(:,:)
    real(dp), allocatable :: interpolation_input(:), interpolation_output(:)
    real(dp) :: gas_total, z_tolerance
    integer :: i, j, ierr, nprofile, ngas

    nprofile = size(z)
    ngas = dat%nq - dat%ng_1 + 1

    if (nprofile < 2) then
      err = 'Altitude initialization requires at least two profile points.'
      return
    endif
    if (size(temperature) /= nprofile .or. size(edd) /= nprofile) then
      err = 'Altitude, temperature, and eddy-diffusion profiles must have the same length.'
      return
    endif
    if (size(mix,1) /= dat%nq .or. size(mix,2) /= nprofile) then
      err = 'The mixing-ratio array has the wrong shape.'
      return
    endif
    if (size(particle_radius,1) /= dat%npq .or. size(particle_radius,2) /= nprofile) then
      err = 'The particle-radius array has the wrong shape.'
      return
    endif

    if (.not. all(ieee_is_finite(z))) then
      err = 'Altitude profile contains a nonfinite value.'
      return
    endif
    if (.not. all(ieee_is_finite(temperature)) .or. any(temperature <= 0.0_dp)) then
      err = 'Temperature profile must contain only finite, positive values.'
      return
    endif
    if (.not. all(ieee_is_finite(edd)) .or. any(edd <= 0.0_dp)) then
      err = 'Eddy-diffusion profile must contain only finite, positive values.'
      return
    endif
    if (.not. ieee_is_finite(surface_pressure) .or. surface_pressure <= 0.0_dp) then
      err = 'Surface pressure must be finite and positive.'
      return
    endif
    if (.not. all(ieee_is_finite(mix)) .or. any(mix < 0.0_dp)) then
      err = 'Mixing ratios must contain only finite, nonnegative values.'
      return
    endif
    if (.not. all(ieee_is_finite(particle_radius)) .or. any(particle_radius <= 0.0_dp)) then
      err = 'Particle radii must contain only finite, positive values.'
      return
    endif
    if (any(z(2:) <= z(:nprofile-1))) then
      err = 'Altitude profile must be strictly increasing.'
      return
    endif
    z_tolerance = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(z(nprofile)))
    if (abs(z(1)) > z_tolerance) then
      err = 'The first altitude profile point must be zero.'
      return
    endif
    if (dat%gas_rainout .and. var%trop_alt >= z(nprofile)) then
      err = 'Tropopause altitude must be below the model-top altitude.'
      return
    endif

    allocate(gas_mix_normalized(ngas,nprofile))
    do j = 1,nprofile
      gas_total = sum(mix(dat%ng_1:dat%nq,j))
      if (.not. ieee_is_finite(gas_total) .or. gas_total <= 0.0_dp) then
        err = 'At least one gas mixing ratio must be positive at every altitude.'
        return
      endif
      gas_mix_normalized(:,j) = max(mix(dat%ng_1:dat%nq,j)/gas_total, mixing_ratio_floor)
      gas_mix_normalized(:,j) = gas_mix_normalized(:,j)/sum(gas_mix_normalized(:,j))
    enddo

    call initialize_altitude_grid(dat, var, z(nprofile))
    call interpolate_temperature_edd(var, z, temperature, edd, err)
    if (allocated(err)) return

    allocate(interpolation_input(nprofile), interpolation_output(var%nz))

    allocate(gas_mix_model(ngas,var%nz))
    do i = 1,ngas
      interpolation_input = log10(max(gas_mix_normalized(i,:), mixing_ratio_floor))
      call interp(var%z, z, interpolation_input, interpolation_output, ierr=ierr)
      if (ierr /= 0) then
        err = 'Unable to interpolate gas mixing ratios onto the model grid.'
        return
      endif
      gas_mix_model(i,:) = 10.0_dp**interpolation_output
    enddo
    do j = 1,var%nz
      gas_mix_model(:,j) = gas_mix_model(:,j)/sum(gas_mix_model(:,j))
      mubar(j) = sum(gas_mix_model(:,j)*dat%species_mass(dat%ng_1:dat%nq))
    enddo
    call press_and_den(var%nz, var%temperature, var%grav, surface_pressure, &
                       var%dz, mubar, pressure, density)
    if (.not. all(ieee_is_finite(pressure)) .or. any(pressure <= 0.0_dp) .or. &
        .not. all(ieee_is_finite(density)) .or. any(density <= 0.0_dp)) then
      err = 'Hydrostatic integration produced an invalid pressure or density.'
      return
    endif

    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz) then
      err = 'The initial atmospheric-state array has the wrong shape.'
      return
    endif
    usol = 0.0_dp
    do i = 1,ngas
      usol(dat%ng_1+i-1,:) = gas_mix_model(i,:)*density
    enddo

    if (dat%npq > 0) then
      allocate(particle_mix_model(dat%npq,var%nz))
      do i = 1,dat%npq
        interpolation_input = log10(max(mix(i,:), mixing_ratio_floor))
        call interp(var%z, z, interpolation_input, interpolation_output, ierr=ierr)
        if (ierr /= 0) then
          err = 'Unable to interpolate particle mixing ratios onto the model grid.'
          return
        endif
        particle_mix_model(i,:) = 10.0_dp**interpolation_output
        usol(i,:) = particle_mix_model(i,:)*density

        interpolation_input = log10(particle_radius(i,:))
        call interp(var%z, z, interpolation_input, interpolation_output, ierr=ierr)
        if (ierr /= 0) then
          err = 'Unable to interpolate particle radii onto the model grid.'
          return
        endif
        var%particle_radius(i,:) = 10.0_dp**interpolation_output
      enddo
    endif

    var%surface_pressure = surface_pressure/1.0e6_dp

  end subroutine

  module subroutine map_atmosphere_p_to_grid(dat, var, profile_pressure, &
                                             temperature, edd, mix, &
                                             particle_radius, trop_p, &
                                             pressure, density, mubar, usol, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_const, only: k_boltz, N_avo
    use photochem_eqns, only: gravity

    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    real(dp), intent(in) :: profile_pressure(:), temperature(:), edd(:)
    real(dp), intent(in) :: mix(:,:), particle_radius(:,:)
    real(dp), optional, intent(in) :: trop_p
    real(dp), intent(out) :: pressure(:), density(:), mubar(:)
    real(dp), intent(out) :: usol(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: z(:), profile_mubar(:)
    real(dp) :: gas_total, inverse_radius, inverse_radius_new
    real(dp) :: inverse_radius_factor, delta_log_pressure
    real(dp) :: surface_z(1), surface_gravity(1)
    real(dp) :: trop_alt_array(1)
    integer :: i, nprofile, ierr

    nprofile = size(profile_pressure)

    ! These checks are needed before pressure and composition can safely be
    ! converted to altitude. The altitude mapper performs the remaining common
    ! profile validation before changing var.
    if (nprofile < 2) then
      err = 'Pressure initialization requires at least two profile points.'
      return
    endif
    if (size(temperature) /= nprofile .or. size(edd) /= nprofile) then
      err = 'Pressure, temperature, and eddy-diffusion profiles must have the same length.'
      return
    endif
    if (size(mix,1) /= dat%nq .or. size(mix,2) /= nprofile) then
      err = 'The mixing-ratio array has the wrong shape.'
      return
    endif
    if (.not. all(ieee_is_finite(profile_pressure)) .or. &
        any(profile_pressure <= 0.0_dp)) then
      err = 'Pressure profile must contain only finite, positive values.'
      return
    endif
    if (any(profile_pressure(2:) >= profile_pressure(:nprofile-1))) then
      err = 'Pressure profile must be strictly decreasing.'
      return
    endif
    if (.not. all(ieee_is_finite(temperature)) .or. any(temperature <= 0.0_dp)) then
      err = 'Temperature profile must contain only finite, positive values.'
      return
    endif
    if (.not. all(ieee_is_finite(mix)) .or. any(mix < 0.0_dp)) then
      err = 'Mixing ratios must contain only finite, nonnegative values.'
      return
    endif

    allocate(z(nprofile), profile_mubar(nprofile))
    do i = 1,nprofile
      gas_total = sum(mix(dat%ng_1:dat%nq,i))
      if (.not. ieee_is_finite(gas_total) .or. gas_total <= 0.0_dp) then
        err = 'At least one gas mixing ratio must be positive at every pressure.'
        return
      endif
      profile_mubar(i) = sum(mix(dat%ng_1:dat%nq,i)* &
                             dat%species_mass(dat%ng_1:dat%nq))/gas_total
      if (.not. ieee_is_finite(profile_mubar(i)) .or. profile_mubar(i) <= 0.0_dp) then
        err = 'Gas composition produced a nonfinite or nonpositive mean molecular weight.'
        return
      endif
    enddo

    ! For spherical gravity, hydrostatic balance can be integrated in inverse
    ! radius without a nonlinear solve:
    !
    !   d(1/r)/d(log(P)) = N_A*k_B*T/(mubar*g_surface*R_surface^2).
    !
    ! Apply the trapezoid rule to T/mubar between pressure knots. The configured
    ! planet radius is the radius at the lower pressure boundary.
    surface_z = 0.0_dp
    call gravity(dat%planet_radius, dat%planet_mass, 1, surface_z, surface_gravity)
    if (.not. ieee_is_finite(surface_gravity(1)) .or. surface_gravity(1) <= 0.0_dp) then
      err = 'Could not compute finite, positive gravity at the lower boundary.'
      return
    endif

    inverse_radius_factor = N_avo*k_boltz/ &
                            (surface_gravity(1)*dat%planet_radius**2)
    inverse_radius = 1.0_dp/dat%planet_radius
    z(1) = 0.0_dp
    do i = 2,nprofile
      delta_log_pressure = log(profile_pressure(i))-log(profile_pressure(i-1))
      inverse_radius_new = inverse_radius + inverse_radius_factor*0.5_dp* &
                           (temperature(i-1)/profile_mubar(i-1) + &
                            temperature(i)/profile_mubar(i))*delta_log_pressure
      if (.not. ieee_is_finite(inverse_radius_new) .or. inverse_radius_new <= 0.0_dp) then
        err = 'Pressure profile extends beyond a finite hydrostatic altitude.'
        return
      endif

      z(i) = 1.0_dp/inverse_radius_new-dat%planet_radius
      if (.not. ieee_is_finite(z(i)) .or. z(i) <= z(i-1)) then
        err = 'Hydrostatic pressure integration did not produce increasing altitude.'
        return
      endif
      inverse_radius = inverse_radius_new
    enddo

    if (present(trop_p)) then
      if (.not. dat%gas_rainout) then
        err = '"trop_p" can only be supplied when gas rainout is enabled.'
        return
      endif
      if (.not. ieee_is_finite(trop_p) .or. trop_p <= 0.0_dp) then
        err = '"trop_p" must be finite and positive.'
        return
      endif
      if (trop_p > profile_pressure(1) .or. trop_p < profile_pressure(nprofile)) then
        err = '"trop_p" must lie within the supplied pressure profile.'
        return
      endif

      ! The pressure profile is descending while the interpolation abscissa
      ! must be ascending. Install the resulting altitude before the common
      ! altitude mapper validates the candidate tropopause.
      call interp([log(trop_p)], log(profile_pressure(nprofile:1:-1)), &
                  z(nprofile:1:-1), trop_alt_array, ierr=ierr)
      if (ierr /= 0 .or. .not. ieee_is_finite(trop_alt_array(1))) then
        err = 'Unable to determine the tropopause altitude from "trop_p".'
        return
      endif
      var%trop_alt = trop_alt_array(1)
    endif

    call map_atmosphere_z_to_grid(dat, var, z, temperature, edd, &
                                  profile_pressure(1), mix, particle_radius, &
                                  pressure, density, mubar, usol, err)

  end subroutine
  
  subroutine initialize_altitude_grid(dat, var, top_atmos)
    use photochem_eqns, only: vertical_grid, gravity

    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    real(dp), intent(in) :: top_atmos

    var%bottom_atmos = 0.0_dp
    var%top_atmos = top_atmos
    call vertical_grid(var%bottom_atmos, var%top_atmos, var%nz, var%z, var%dz)
    call gravity(dat%planet_radius, dat%planet_mass, var%nz, var%z, var%grav)

  end subroutine

  module subroutine map_atmosphere_file_to_grid(dat, var, profile, usol, err)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    type(AtmosphereFileProfile), intent(in) :: profile
    real(dp), intent(out) :: usol(:,:)
    character(:), allocatable, intent(out) :: err

    call resolve_atmosphere_settings(profile, dat, var, err)
    if (allocated(err)) return

    call initialize_altitude_grid(dat, var, var%top_atmos)
    call interp2atmosfile(dat, var, profile, usol, err)
    if (allocated(err)) return

  end subroutine

  module subroutine resolve_atmosphere_settings(profile, dat, var, err)
    type(AtmosphereFileProfile), intent(in) :: profile
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err

    if (profile%nlayer < 2) then
      err = 'Atmosphere file must contain at least two data rows.'
      return
    endif

    var%bottom_atmos = 0.0_dp
    var%top_atmos = profile%z(profile%nlayer) + &
                    0.5_dp*(profile%z(profile%nlayer) - profile%z(profile%nlayer-1))

    if (var%top_atmos < var%bottom_atmos) then
      err = 'The top of the atmosphere must be bigger than the bottom'
      return
    endif

    if (dat%gas_rainout .and. var%trop_alt > var%top_atmos) then
      err = 'IOError: tropopause-altitude must be between the top and bottom of the atmosphere'
      return
    endif

  end subroutine

  subroutine interpolate_temperature_edd(var, z, temperature, edd, err, absolute_eddy)
    use futils, only: interp

    type(PhotochemVars), intent(inout) :: var
    real(dp), intent(in) :: z(:), temperature(:), edd(:)
    character(:), allocatable, intent(out) :: err
    logical, intent(in), optional :: absolute_eddy

    real(dp), allocatable :: interpolation_input(:), interpolation_output(:)
    logical :: use_absolute_eddy
    integer :: ierr

    use_absolute_eddy = .false.
    if (present(absolute_eddy)) use_absolute_eddy = absolute_eddy

    call interp(var%z, z, temperature, var%temperature, ierr=ierr)
    if (ierr /= 0) then
      err = 'Unable to interpolate temperature onto the model grid.'
      return
    endif

    allocate(interpolation_input(size(z)), interpolation_output(var%nz))
    if (use_absolute_eddy) then
      interpolation_input = log10(abs(edd))
    else
      interpolation_input = log10(edd)
    endif
    call interp(var%z, z, interpolation_input, interpolation_output, ierr=ierr)
    if (ierr /= 0) then
      err = 'Unable to interpolate eddy diffusion onto the model grid.'
      return
    endif
    var%edd = 10.0_dp**interpolation_output

  end subroutine

  subroutine interp2atmosfile(dat, var, profile, usol, err)
    use futils, only: interp
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    type(AtmosphereFileProfile), intent(in) :: profile
    real(dp), intent(out) :: usol(:,:)
    character(:), allocatable, intent(out) :: err
    
    integer :: i, ierr

    ! Preserve the file initializer's historical absolute-value treatment of
    ! eddy diffusion while sharing the common temperature/Kzz interpolation
    ! kernel used by altitude-based initialization.
    call interpolate_temperature_edd(var, profile%z, profile%temperature, &
                                     profile%edd, err, absolute_eddy=.true.)
    if (allocated(err)) return

    call interp2atmosfile_mix(dat, var, profile, usol, err)
    if (allocated(err)) return
    
    if (dat%there_are_particles) then
      do i = 1,dat%npq
        call interp(var%nz, profile%nlayer, var%z, profile%z, &
                    log10(abs(profile%particle_radius(i,:))), var%particle_radius(i,:), ierr)
        if (ierr /= 0) then
          err = 'Subroutine interp returned an error.'
          return
        endif
      enddo
      var%particle_radius = 10.0_dp**var%particle_radius
    endif
    
  end subroutine

  subroutine interp2atmosfile_mix(dat, var, profile, usol, err)
    use futils, only: interp
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    type(AtmosphereFileProfile), intent(in) :: profile
    real(dp), intent(out) :: usol(:,:)
    character(:), allocatable, intent(out) :: err

    integer :: i, ierr
    real(dp), allocatable :: density(:)

    allocate(density(var%nz))

    ! Interpolate file density to model grid
    call interp(var%z, profile%z, log10(profile%density), density, linear_extrap=.true., ierr=ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
      return
    endif
    density = 10.0_dp**density

    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz) then
      err = 'The initial atmospheric-state array has the wrong shape.'
      return
    endif

    do i = 1,dat%nq
      call interp(var%nz, profile%nlayer, var%z, profile%z,&
                  log10(abs(profile%mix(i,:))), usol(i,:), ierr)
      if (ierr /= 0) then
        err = 'Subroutine interp returned an error.'
        return
      endif
    enddo
    usol = 10.0_dp**usol

    do i = 1,var%nz
      usol(:,i) = usol(:,i)*density(i)
    enddo 

  end subroutine

  !~~ The routines below finalized setup after grid construction.

  module subroutine finalize_atmosphere_initialization(dat, var, err)
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err

    call interp2particlexsdata(dat, var%particle_radius, var%particle_xs, err)
    if (allocated(err)) return

    call refresh_temperature_dependent_state(dat, var, err=err)
    if (allocated(err)) return

  end subroutine

  module subroutine refresh_temperature_dependent_state(dat, var, trop_alt, err)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    real(dp), optional, intent(in) :: trop_alt
    character(:), allocatable, intent(out) :: err

    call interp2xsdata(dat, var%xs_x_qy, err)
    if (allocated(err)) return
    
    call compute_gibbs_energy(dat, var%temperature, var%gibbs_energy, err)
    if (allocated(err)) return

    call set_tropopause(dat, var%z, var%bottom_atmos, var%top_atmos, trop_alt, var%trop_alt, var%trop_ind, err)
    if (allocated(err)) return

  end subroutine

  subroutine set_tropopause(dat, z, bottom_atmos, top_atmos, trop_alt_new, trop_alt, trop_ind, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: z(:), bottom_atmos, top_atmos
    real(dp), optional, intent(in) :: trop_alt_new
    real(dp), intent(inout) :: trop_alt
    integer, intent(inout) :: trop_ind
    character(:), allocatable, intent(out) :: err

    trop_ind = 1
    if (.not. dat%gas_rainout) return

    if (present(trop_alt_new)) trop_alt = trop_alt_new

    if (.not. ieee_is_finite(trop_alt) .or. &
        trop_alt < bottom_atmos .or. trop_alt > top_atmos) then
      err = 'trop_alt is above or bellow the atmosphere!'
      return
    endif

    trop_ind = max(minloc(abs(z - trop_alt), 1) - 1, 1)
    if (trop_ind < 3) then
      err = 'Tropopause is too low.'
      return
    elseif (trop_ind > size(z) - 2) then
      err = 'Tropopause is too high.'
      return
    endif

  end subroutine

  subroutine compute_gibbs_energy(dat, temperature, gibbs_energy, err)
    use photochem_eqns, only: gibbs_energy_eval
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: temperature(:)
    real(dp), intent(inout) :: gibbs_energy(:,:)
    character(:), allocatable, intent(out) :: err
    
    integer :: i, j
    logical :: found

    if (.not. dat%reverse) return
    
    do i = 1,dat%ng
      do j = 1,size(temperature)
        call gibbs_energy_eval(dat%thermo_data(i), temperature(j), &
                               found, gibbs_energy(j,i))
        if (.not. found) then
          err = 'The temperature is not within the ranges '// &
                'given for the thermodynamic data for '//trim(dat%species_names(i+dat%npq))
          return
        endif
      enddo
    enddo

  end subroutine
  
  subroutine interp2xsdata(dat, xs_x_qy, err)
    use photochem_const, only: smaller_real
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(inout) :: xs_x_qy(:,:,:)
    character(:), allocatable, intent(out) :: err
    
    integer :: i, k

    ! No temperature dependence, so we just copy over
    do k = 1, dat%nw
      do i = 1,dat%kj
        xs_x_qy(:,i,k) = abs(dat%photolysis_xs(i,k)) + smaller_real
      enddo
    enddo

  end subroutine

  subroutine interp2particlexsdata(dat, particle_radius, particle_xs, err)
    use photochem_types, only: ParticleXsections
    type(PhotochemData), intent(in) :: dat
    real(dp), intent(in) :: particle_radius(:,:)
    type(ParticleXsections), intent(inout) :: particle_xs(:)
    character(:), allocatable, intent(out) :: err
    
    integer :: i, j, k, jj, nz
    real(dp) :: dr, slope, intercept

    if (.not.dat%there_are_particles) return

    nz = size(particle_radius, 2)
    
    do j = 1,nz
      do k = 1,dat%np
        ! if there is optical data, then check that the
        ! data covers the particle radii in the atmosphere.
        if (dat%part_xs_file(k)%ThereIsData) then
        
          if (particle_radius(k,j) <= dat%radii_file(1,k)) then
            err = "There is not any optical data for the "// &
                  "particle radii specified in the atmosphere."
            return
          endif
          if (particle_radius(k,j) >= dat%radii_file(dat%nrad_file,k)) then
            err = "There is not any optical data for the "// &
                  "particle radii specified in the atmosphere."
            return
          endif
          
        endif
        
      enddo
    enddo
    do i = 1,dat%nw
      do j = 1,nz
        do k = 1,dat%np
          
        ! if there is particle optical data, then linearly interpolate
        ! it to to the particle radii in the atmosphere.
        if (dat%part_xs_file(k)%ThereIsData) then
          
          do jj = 1,dat%nrad_file-1
            if (particle_radius(k,j) >= dat%radii_file(jj,k) .and. &
                particle_radius(k,j) < dat%radii_file(jj+1,k)) then
  
              dr = dat%radii_file(jj+1,k) - dat%radii_file(jj,k)
  
              slope = (dat%part_xs_file(k)%w0(jj+1,i) - dat%part_xs_file(k)%w0(jj,i))/dr
              intercept = dat%part_xs_file(k)%w0(jj,i) - dat%radii_file(jj,k)*slope
              particle_xs(k)%w0(j,i) = slope*particle_radius(k,j) + intercept
              
              slope = (dat%part_xs_file(k)%qext(jj+1,i) - dat%part_xs_file(k)%qext(jj,i))/dr
              intercept = dat%part_xs_file(k)%qext(jj,i) - dat%radii_file(jj,k)*slope
              particle_xs(k)%qext(j,i) = slope*particle_radius(k,j) + intercept
  
              slope = (dat%part_xs_file(k)%gt(jj+1,i) - dat%part_xs_file(k)%gt(jj,i))/dr
              intercept = dat%part_xs_file(k)%gt(jj,i) - dat%radii_file(jj,k)*slope
              particle_xs(k)%gt(j,i) = slope*particle_radius(k,j) + intercept
            endif
          enddo
        
        endif 
          
        enddo
      enddo
    enddo
    
  end subroutine

end submodule
