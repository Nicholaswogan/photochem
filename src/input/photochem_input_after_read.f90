submodule (photochem_input) photochem_input_after_read
  implicit none
  
contains

  module subroutine map_atmosphere_z_to_grid(dat, var, z, temperature, &
                                             edd, surface_pressure, gas_mix, &
                                             particle_mix, particle_radius, &
                                             pressure, density, mubar, err)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use futils, only: interp
    use photochem_eqns, only: vertical_grid, gravity, press_and_den

    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    real(dp), intent(in) :: z(:), temperature(:), edd(:)
    real(dp), intent(in) :: surface_pressure
    real(dp), intent(in) :: gas_mix(:,:), particle_mix(:,:), particle_radius(:,:)
    real(dp), intent(out) :: pressure(:), density(:), mubar(:)
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
    if (size(gas_mix,1) /= ngas .or. size(gas_mix,2) /= nprofile) then
      err = 'The gas mixing-ratio array has the wrong shape.'
      return
    endif
    if (size(particle_mix,1) /= dat%npq .or. size(particle_mix,2) /= nprofile) then
      err = 'The particle mixing-ratio array has the wrong shape.'
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
    if (.not. all(ieee_is_finite(gas_mix)) .or. any(gas_mix < 0.0_dp)) then
      err = 'Gas mixing ratios must contain only finite, nonnegative values.'
      return
    endif
    if (.not. all(ieee_is_finite(particle_mix)) .or. any(particle_mix < 0.0_dp)) then
      err = 'Particle mixing ratios must contain only finite, nonnegative values.'
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
      gas_total = sum(gas_mix(:,j))
      if (.not. ieee_is_finite(gas_total) .or. gas_total <= 0.0_dp) then
        err = 'At least one gas mixing ratio must be positive at every altitude.'
        return
      endif
      gas_mix_normalized(:,j) = max(gas_mix(:,j)/gas_total, mixing_ratio_floor)
      gas_mix_normalized(:,j) = gas_mix_normalized(:,j)/sum(gas_mix_normalized(:,j))
    enddo

    var%bottom_atmos = 0.0_dp
    var%top_atmos = z(nprofile)
    call vertical_grid(var%bottom_atmos, var%top_atmos, var%nz, var%z, var%dz)
    call gravity(dat%planet_radius, dat%planet_mass, var%nz, var%z, var%grav)

    call interp(var%z, z, temperature, var%temperature, ierr=ierr)
    if (ierr /= 0) then
      err = 'Unable to interpolate temperature onto the model grid.'
      return
    endif

    allocate(interpolation_input(nprofile), interpolation_output(var%nz))
    interpolation_input = log10(edd)
    call interp(var%z, z, interpolation_input, interpolation_output, ierr=ierr)
    if (ierr /= 0) then
      err = 'Unable to interpolate eddy diffusion onto the model grid.'
      return
    endif
    var%edd = 10.0_dp**interpolation_output

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

    var%usol_init = 0.0_dp
    do i = 1,ngas
      var%usol_init(dat%ng_1+i-1,:) = gas_mix_model(i,:)*density
    enddo

    if (dat%npq > 0) then
      allocate(particle_mix_model(dat%npq,var%nz))
      do i = 1,dat%npq
        interpolation_input = log10(max(particle_mix(i,:), mixing_ratio_floor))
        call interp(var%z, z, interpolation_input, interpolation_output, ierr=ierr)
        if (ierr /= 0) then
          err = 'Unable to interpolate particle mixing ratios onto the model grid.'
          return
        endif
        particle_mix_model(i,:) = 10.0_dp**interpolation_output
        var%usol_init(i,:) = particle_mix_model(i,:)*density

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
  
  module subroutine after_read_setup(dat, var, err)
    use photochem_eqns, only: vertical_grid, gravity
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    ! set up the atmosphere grid
    call vertical_grid(var%bottom_atmos, var%top_atmos, &
                       var%nz, var%z, var%dz)
    call gravity(dat%planet_radius, dat%planet_mass, &
                 var%nz, var%z, var%grav)
    call interp2atmosfile(dat, var, err)
    if (allocated(err)) return

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
      ! we have a tropopause
      var%trop_ind = max(minloc(abs(var%z - var%trop_alt), 1) - 1, 1)

      if (var%trop_ind < 3) then
        err = 'Tropopause is too low.'
        return
      elseif (var%trop_ind > var%nz-2) then
        err = 'Tropopause is too high.'
        return
      endif
    else
      var%trop_ind = 1
    endif
    
  end subroutine
  
  subroutine interp2atmosfile(dat, var, err)
    use futils, only: interp, conserving_rebin
    use photochem_const, only: small_real
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    integer :: i, ierr
    
    call interp(var%nz, dat%nzf, var%z, dat%z_file, dat%T_file, var%Temperature, ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
      return
    endif
    
    call interp(var%nz, dat%nzf, var%z, dat%z_file, log10(dabs(dat%edd_file)), var%edd, ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
      return
    endif
    var%edd = 10.0_dp**var%edd

    if (dat%conserving_init) then
      call interp2atmosfile_mixconserving(dat, var, err)
      if (allocated(err)) return
    else
      call interp2atmosfile_mix(dat, var, err)
      if (allocated(err)) return
    endif
    
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
    
  end subroutine

  subroutine interp2atmosfile_mixconserving(dat, var, err)
    use futils, only: interp, conserving_rebin
    use photochem_const, only: small_real
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err

    integer :: i, ierr
    real(dp) :: dz_file
    real(dp), allocatable :: densities_file(:,:) ! molecules/cm3
    real(dp), allocatable :: ze_file(:), ze(:)

    dz_file = dat%z_file(2)-dat%z_file(1)

    allocate(densities_file(dat%nq,dat%nzf))
    allocate(ze_file(dat%nzf+1))
    allocate(ze(var%nz+1))

    do i = 1,dat%nq
      densities_file(i,:) = dat%mix_file(i,:)*dat%den_file
    enddo

    ze_file(1) = dat%z_file(1) - 0.5_dp*dz_file
    do i = 1,dat%nzf
      ze_file(i+1) = dat%z_file(i) + 0.5_dp*dz_file
    enddo
    ze = var%z(1) - 0.5_dp*var%dz(1)
    do i = 1,var%nz
      ze(i+1) = var%z(i) + 0.5_dp*var%dz(i)
    enddo

    do i = 1,dat%nq
      call conserving_rebin(ze_file, densities_file(i,:), ze, var%usol_init(i,:), ierr)
      if (ierr /= 0) then
        err = 'subroutine conserving_rebin returned an error'
        return
      endif
    enddo

  end subroutine

  subroutine interp2atmosfile_mix(dat, var, err)
    use futils, only: interp
    use photochem_const, only: small_real
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err

    integer :: i, ierr
    real(dp), allocatable :: density(:)

    allocate(density(var%nz))

    ! Interpolate file density to model grid
    call interp(var%z, dat%z_file, log10(dat%den_file), density, linear_extrap=.true., ierr=ierr)
    if (ierr /= 0) then
      err = 'Subroutine interp returned an error.'
      return
    endif
    density = 10.0_dp**density

    do i = 1,dat%nq
      call interp(var%nz, dat%nzf, var%z, dat%z_file,&
                  log10(abs(dat%mix_file(i,:))), var%usol_init(i,:), ierr)
      if (ierr /= 0) then
        err = 'Subroutine interp returned an error.'
        return
      endif
    enddo
    var%usol_init = 10.0_dp**var%usol_init

    do i = 1,var%nz
      var%usol_init(:,i) = var%usol_init(:,i)*density(i)
    enddo 

  end subroutine

  module subroutine compute_gibbs_energy(dat, var, err)
    use photochem_eqns, only: gibbs_energy_eval
    
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    integer :: i, j
    logical :: found

    
    do i = 1,dat%ng
      do j = 1,var%nz
        call gibbs_energy_eval(dat%thermo_data(i), var%temperature(j), &
                               found, var%gibbs_energy(j,i))
        if (.not. found) then
          err = 'The temperature is not within the ranges '// &
                'given for the thermodynamic data for '//trim(dat%species_names(i+dat%npq))
          return
        endif
      enddo
    enddo

  end subroutine
  
  module subroutine interp2xsdata(dat, var, err)
    use futils, only: interp
    use photochem_const, only: smaller_real
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var

    character(:), allocatable, intent(out) :: err
    
    integer :: i, k

    ! No temperature dependence, so we just copy over
    do k = 1, dat%nw
      do i = 1,dat%kj
        var%xs_x_qy(:,i,k) = abs(dat%photolysis_xs(i,k))+smaller_real
      enddo
    enddo

  end subroutine

  module subroutine interp2particlexsdata(dat, var, err)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var

    character(:), allocatable, intent(out) :: err
    
    integer :: i, j, k, jj
    real(dp) :: dr, slope, intercept
    
    ! particles
    if (dat%there_are_particles) then
      do j = 1,var%nz
        do k = 1,dat%np
          ! if there is optical data, then check that the
          ! data covers the particle radii in the atmosphere.
          if (dat%part_xs_file(k)%ThereIsData) then
          
            if (var%particle_radius(k,j) <= dat%radii_file(1,k)) then
              err = "There is not any optical data for the "// &
                    "particle radii specified in the atmosphere."
              return
            endif
            if (var%particle_radius(k,j) >= dat%radii_file(dat%nrad_file,k)) then
              err = "There is not any optical data for the "// &
                    "particle radii specified in the atmosphere."
              return
            endif
            
          endif
          
        enddo
      enddo
      do i = 1,dat%nw
        do j = 1,var%nz
          do k = 1,dat%np
            
          ! if there is particle optical data, then linearly interpolate
          ! it to to the particle radii in the atmosphere.
          if (dat%part_xs_file(k)%ThereIsData) then
            
            do jj = 1,dat%nrad_file-1
              if (var%particle_radius(k,j) >= dat%radii_file(jj,k) .and. &
                  var%particle_radius(k,j) < dat%radii_file(jj+1,k)) then
    
                dr = dat%radii_file(jj+1,k) - dat%radii_file(jj,k)
    
                slope = (dat%part_xs_file(k)%w0(jj+1,i) - dat%part_xs_file(k)%w0(jj,i))/dr
                intercept = dat%part_xs_file(k)%w0(jj,i) - dat%radii_file(jj,k)*slope
                var%particle_xs(k)%w0(j,i) = slope*var%particle_radius(k,j) + intercept
                
                slope = (dat%part_xs_file(k)%qext(jj+1,i) - dat%part_xs_file(k)%qext(jj,i))/dr
                intercept = dat%part_xs_file(k)%qext(jj,i) - dat%radii_file(jj,k)*slope
                var%particle_xs(k)%qext(j,i) = slope*var%particle_radius(k,j) + intercept
    
                slope = (dat%part_xs_file(k)%gt(jj+1,i) - dat%part_xs_file(k)%gt(jj,i))/dr
                intercept = dat%part_xs_file(k)%gt(jj,i) - dat%radii_file(jj,k)*slope
                var%particle_xs(k)%gt(j,i) = slope*var%particle_radius(k,j) + intercept
              endif
            enddo
          
          endif 
            
          enddo
        enddo
      enddo
    endif
    
  end subroutine
  
  module subroutine allocate_nz_vars(dat, var)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    
    integer :: i
    
    var%neqs = dat%nq*var%nz

    allocate(var%temperature(var%nz))
    allocate(var%z(var%nz))
    allocate(var%dz(var%nz))
    allocate(var%edd(var%nz))
    allocate(var%grav(var%nz))
    allocate(var%usol_init(dat%nq,var%nz))
    allocate(var%particle_radius(dat%npq,var%nz))
    allocate(var%xs_x_qy(var%nz,dat%kj,dat%nw))
    allocate(var%usol_out(dat%nq,var%nz))
    
    allocate(var%particle_xs(dat%np))
    do i = 1,dat%np
      ! only allocate space if there is data
      if (dat%part_xs_file(i)%ThereIsData) then
        var%particle_xs(i)%ThereIsData = .true.
        allocate(var%particle_xs(i)%w0(var%nz,dat%nw))
        allocate(var%particle_xs(i)%qext(var%nz,dat%nw))
        allocate(var%particle_xs(i)%gt(var%nz,dat%nw))
      else
        var%particle_xs(i)%ThereIsData = .false.
      endif
    enddo
    
    if (dat%reverse) then
      allocate(var%gibbs_energy(var%nz,dat%ng))
    endif

    allocate(var%tauc(var%nz,dat%nw))
    var%tauc = 0.0_dp
    allocate(var%w0c(var%nz,dat%nw))
    var%w0c = 0.0_dp
    allocate(var%g0c(var%nz,dat%nw))
    var%g0c = 0.0_dp

  end subroutine

end submodule
