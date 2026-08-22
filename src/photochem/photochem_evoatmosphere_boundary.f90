submodule(photochem_evoatmosphere) photochem_evoatmosphere_boundary
  implicit none

contains

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
    type(PhotochemWrk), pointer :: wrk

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

    call self%chemistry_right_hand_side(wrk%usol, rhs, err)
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
    use photochem_vars, only: time_dependent_rate_fcn
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

end submodule
