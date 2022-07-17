submodule (photochem_input) photochem_input_after_read
  implicit none
  
contains
  
  module subroutine after_read_setup(dat, var, err)
    use photochem_eqns, only: vertical_grid, gravity
    type(PhotochemData), intent(inout) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    dat%kd = 2*dat%nq + 1
    dat%kl = dat%kd + dat%nq
    dat%ku = dat%kd - dat%nq
    dat%lda = 3*dat%nq + 1
    
    call allocate_nz_vars(dat, var)
    ! set up the atmosphere grid
    call vertical_grid(var%bottom_atmos, var%top_atmos, &
                       var%nz, var%z, var%dz)
    call gravity(dat%planet_radius, dat%planet_mass, &
                 var%nz, var%z, var%grav)
    call interp2atmosfile(dat, var, err)
    if (allocated(err)) return
    
    ! all below depends on Temperature
    call interp2xsdata(dat, var, err)
    if (allocated(err)) return
    
    if (dat%reverse) then
      call compute_gibbs_energy(dat, var, err)
      if (allocated(err)) return
    endif
    
    if (dat%fix_water_in_trop .or. dat%gas_rainout) then
      ! we have a tropopause
      var%trop_ind = minloc(var%z,1, &
                          var%z .ge. var%trop_alt) - 1
    else
      var%trop_ind = 0
    endif
    
  end subroutine
  
  subroutine interp2atmosfile(dat, var, err)
    use futils, only: interp, rebin, is_close
    use photochem_const, only: small_real
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: var
    character(:), allocatable, intent(out) :: err
    
    integer :: i, j, k, ierr
    
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

    if (dat%back_gas) then
      do i = 1,dat%nq
        call interp(var%nz, dat%nzf, var%z, dat%z_file,&
                    log10(abs(dat%mix_file(i,:))), var%usol_init(i,:), ierr)
        if (ierr /= 0) then
          err = 'Subroutine interp returned an error.'
          return
        endif
      enddo
      var%usol_init = 10.0_dp**var%usol_init

    else
    block
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
        call rebin(ze_file, densities_file(i,:), ze, var%usol_init(i,:))
      enddo

      if (ze(var%nz+1) >= ze_file(dat%nzf+1) .and. ze(var%nz) < ze_file(dat%nzf+1)) then
        ! The model grid empasses the file grid. We are OK. No mass is lost.
        ! <--down--|___|___|____|___| (grid)
        ! <--down--|___|___|___|___| (file)
      elseif (ze(var%nz) > ze_file(dat%nzf+1)) then
        ! upper cell of grid is completely above upper cell in file.
        ! We must re-distribute molecules to prevent program from crashing
        ! <--down--|___|___|___|___| (grid)
        ! <--down--|___|___|         (file)
        do k = 1,var%nz
          j = var%nz+1-k
          if (all(var%usol_init(:,j) /= 0.0_dp)) then
            exit
          endif
        enddo
        ! distribute molecules in the upper most grid cell among the higher altitude cells
        do i = 1,dat%nq
          var%usol_init(i,j:) = var%usol_init(i,j)*var%dz(j)/sum(var%dz(j:))
        enddo

      elseif (ze(var%nz+1) < ze_file(dat%nzf+1)) then
        ! we have lost mass
        ! <--down--|___|___|___|___| (grid)
        ! <--down--|___|___|____|____| (file)
        ! err = "highest edge of the grid must be >= the highest edge in the 'atmosphere.txt' file"
        ! return

      endif
      
      if (ze(1) /= ze_file(1)) then
        err = "Lowest edge of grid must be the same as the lowest edge in the 'atmosphere.txt' file"
        return
      endif

    end block
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
    
    if (var%z(1) < dat%z_file(1)) then
      print*,'Warning: vertical grid is being extrapolated below where there is input data.'
    endif
    
    if (var%z(var%nz) > dat%z_file(dat%nzf)) then
      print*,'Warning: vertical grid is being extrapolated above where there is input data.'
    endif

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
    
    integer :: i, j, k, jj
    real(dp) :: val(1), T_temp(1)
    real(dp) :: dr, slope, intercept
    

    do k = 1, dat%nw
      do i = 1,dat%kj
        if (dat%xs_data(i)%n_temps > 1) then
          do j = 1,var%nz
            T_temp(1) = var%temperature(j)
      
            call interp(1, dat%xs_data(i)%n_temps, T_temp, dat%xs_data(i)%xs_temps, &
                        log10(abs(dat%xs_data(i)%xs(:,k))+smaller_real), val)
                        
            var%xs_x_qy(j,i,k) = 10.0_dp**val(1)
          enddo
        else
          do j = 1,var%nz
            var%xs_x_qy(j,i,k) = 10.0_dp**log10(abs(dat%xs_data(i)%xs(1,k))+smaller_real)   
          enddo
        endif
      enddo
    enddo
    
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
  
  subroutine allocate_nz_vars(dat, vars)
    type(PhotochemData), intent(in) :: dat
    type(PhotochemVars), intent(inout) :: vars
    
    integer :: i
    
    vars%neqs = dat%nq*vars%nz

    allocate(vars%temperature(vars%nz))
    allocate(vars%z(vars%nz))
    allocate(vars%dz(vars%nz))
    allocate(vars%edd(vars%nz))
    allocate(vars%grav(vars%nz))
    allocate(vars%usol_init(dat%nq,vars%nz))
    allocate(vars%particle_radius(dat%npq,vars%nz))
    allocate(vars%xs_x_qy(vars%nz,dat%kj,dat%nw))
    allocate(vars%usol_out(dat%nq,vars%nz))
    
    allocate(vars%particle_xs(dat%np))
    do i = 1,dat%np
      ! only allocate space if there is data
      if (dat%part_xs_file(i)%ThereIsData) then
        vars%particle_xs(i)%ThereIsData = .true.
        allocate(vars%particle_xs(i)%w0(vars%nz,dat%nw))
        allocate(vars%particle_xs(i)%qext(vars%nz,dat%nw))
        allocate(vars%particle_xs(i)%gt(vars%nz,dat%nw))
      else
        vars%particle_xs(i)%ThereIsData = .false.
      endif
    enddo
    
    if (dat%reverse) then
      allocate(vars%gibbs_energy(vars%nz,dat%ng))
    endif

  end subroutine

end submodule