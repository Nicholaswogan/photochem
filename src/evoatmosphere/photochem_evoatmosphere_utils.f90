
submodule(photochem_evoatmosphere) photochem_evoatmosphere_utils
  implicit none

contains

  module subroutine update_vertical_grid(self, usol_old, top_atmos, usol_new, err)
    use photochem_input, only: interp2particlexsdata, interp2xsdata, compute_gibbs_energy
    use photochem_eqns, only: vertical_grid, gravity
    use futils, only: interp, conserving_rebin
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_old(:,:)
    real(dp), intent(in) :: top_atmos
    real(dp), intent(out) :: usol_new(:,:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: ze_old(:), ze(:)
    real(dp), allocatable :: z_old(:), dz_old(:)
    integer :: i, ierr

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! save z and dz
    allocate(z_old(var%nz), dz_old(var%nz))
    z_old = var%z
    dz_old = var%dz

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

    ! We do a conserving rebin of the densities
    allocate(ze_old(size(z_old)+1))
    allocate(ze(var%nz+1))

    ze_old(1) = z_old(1) - 0.5_dp*dz_old(1)
    do i = 1,size(z_old)
      ze_old(i+1) = z_old(i) + 0.5_dp*dz_old(i)
    enddo
    ze = var%z(1) - 0.5_dp*var%dz(1)
    do i = 1,var%nz
      ze(i+1) = var%z(i) + 0.5_dp*var%dz(i)
    enddo

    do i = 1,dat%nq
      call conserving_rebin(ze_old, usol_old(i,:), ze, usol_new(i,:), ierr)
      if (ierr /= 0) then
        err = 'subroutine conserving_rebin returned an error'
        return
      endif
    enddo

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
      var%trop_ind = 0
    endif

  end subroutine

end submodule