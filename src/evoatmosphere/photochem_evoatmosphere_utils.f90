
submodule(photochem_evoatmosphere) photochem_evoatmosphere_utils
  implicit none

contains

module subroutine out2atmosphere_txt(self, filename, overwrite, clip, err)
    use photochem_common, only: out2atmosphere_txt_base
    class(EvoAtmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    logical, intent(in) :: overwrite, clip
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
                                 filename, overwrite, clip, err)
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

  module subroutine rebin_update_vertical_grid(self, usol_old, top_atmos, usol_new, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol_old(:,:)
    real(dp), intent(in) :: top_atmos
    real(dp), intent(out) :: usol_new(:,:)
    character(:), allocatable, intent(out) :: err

    call rebin_densities(self, usol_old, top_atmos, usol_new, err)
    if(allocated(err)) return
    call update_vertical_grid(self, usol_new, top_atmos, err)
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

  subroutine update_vertical_grid(self, usol_new, top_atmos, err)
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

    call update_vertical_grid(self, usol_new, top_atmos, err)
    if(allocated(err)) return
    call self%prep_atmosphere(usol_new, err)
    if(allocated(err)) return

  end subroutine

end submodule