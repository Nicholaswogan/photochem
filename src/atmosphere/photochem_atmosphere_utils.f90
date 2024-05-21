submodule(photochem_atmosphere) photochem_atmosphere_utils
  implicit none
  
  ! Contains routines utility routines for returning or saving 
  ! model output
  
contains
  
  module subroutine out2atmosphere_txt(self, filename, overwrite, clip, err)
    use photochem_common, only: out2atmosphere_txt_base
    class(Atmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    logical, intent(in) :: overwrite, clip
    character(:), allocatable, intent(out) :: err
    
    real(dp) :: rhs(self%var%neqs)  
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
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
  
  module subroutine out2in(self, err)
    class(Atmosphere), intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    
    if (self%var%at_photo_equilibrium) then
      self%var%usol_init = self%var%usol_out
    else
      err = "Can not set output to input without first converging to photochemical equilibrium."
      return
    endif
  end subroutine
  
  module subroutine gas_fluxes(self, surf_fluxes, top_fluxes, err)
    class(Atmosphere), target, intent(inout) :: self
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
                              *wrk%density(1)*var%dz(1)
      chemical_production = rhs(i)*wrk%density(1)*var%dz(1)
      surf_fluxes(i) = -(diffusive_production + chemical_production)
    enddo
    
    ! fluxes going into or out of the top of the atmosphere.
    do i = 1,dat%nq
      diffusive_production = &
         (wrk%DD(i,var%nz)*wrk%usol(i,var%nz) + wrk%ADD(i,var%nz)*wrk%usol(i,var%nz) &
        + wrk%DL(i,var%nz)*wrk%usol(i,var%nz-1) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz-1)) &
          *wrk%density(var%nz)*var%dz(var%nz)
    
      chemical_production = rhs(i + (var%nz-1)*dat%nq)*wrk%density(var%nz)*var%dz(var%nz)
      top_fluxes(i) = diffusive_production + chemical_production
    enddo
    
  end subroutine
  
  module function atom_conservation(self, atom, err) result(con)
    use photochem_enum, only: VelocityDistributedFluxBC
    use photochem_eqns, only: damp_condensation_rate
    use photochem_types, only: AtomConservation
    use photochem_const, only: small_real
    class(Atmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: atom
    character(:), allocatable, intent(out) :: err
    type(AtomConservation) :: con
    
    real(dp) :: surf_fluxes(self%dat%nq)
    real(dp) :: top_fluxes(self%dat%nq)
    real(dp) :: integrated_rainout(self%dat%nq)
    real(dp) :: rh, df_gas_dt, cond_rate0, cond_rate, con_evap_rate
    
    integer :: ind(1), i, j, kk
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    ind = findloc(dat%atoms_names,atom)
    kk = ind(1)
    if (ind(1) == 0) then
      err = "Atom "//trim(atom)//" is not in the list of atoms."
      return
    endif
    if (dat%species_composition(ind(1),dat%nsp) /= 0) then
      err = "Atom "//trim(atom)//" makes up the background gas"// &
            " which is not conserved."
      return 
    endif
    
    call self%gas_fluxes(surf_fluxes, top_fluxes, err)
    if (allocated(err)) return
    
    con%in_surf = 0
    con%in_top = 0
    con%in_dist = 0
    con%in_other = 0
    con%out_surf = 0
    con%out_top = 0
    con%out_rain = 0
    con%out_other = 0
    
    ! Upper and lower boundary
    do i = 1,dat%nq
      if (surf_fluxes(i) > 0) then
        con%in_surf = con%in_surf + surf_fluxes(i)*dat%species_composition(kk,i)
      else
        con%out_surf = con%out_surf + (-1.0_dp)*surf_fluxes(i)*dat%species_composition(kk,i)
      endif
      if (top_fluxes(i) > 0) then
        con%out_top = con%out_top + top_fluxes(i)*dat%species_composition(kk,i)
      else
        con%in_top = con%in_top + (-1.0_dp)*top_fluxes(i)*dat%species_composition(kk,i)
      endif
    enddo
    
    ! distributed fluxes
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == VelocityDistributedFluxBC) then
        con%in_dist = con%in_dist + var%lower_flux(i)*dat%species_composition(kk,i)
      endif
    enddo
    
    ! rainout
    if (dat%gas_rainout) then
      integrated_rainout = 0.0_dp
      do j = 1,var%trop_ind
        do i = 1,dat%nq
          integrated_rainout(i) = integrated_rainout(i) + &
                wrk%rainout_rates(i,j)*wrk%usol(i,j)*wrk%density(j)*var%dz(j)
        enddo
      enddo
      
      do i = 1,dat%nq
        con%out_rain = con%out_rain + integrated_rainout(i)*dat%species_composition(kk,i)
      enddo
    endif

    if (dat%fix_water_in_trop) then
      do j = 1,var%trop_ind
        con_evap_rate = var%fast_arbitrary_rate*(wrk%H2O_sat_mix(j)*wrk%H2O_rh(j) - wrk%usol(dat%LH2O,j)) &
                        *wrk%density(j)*var%dz(j)*dat%species_composition(kk,dat%LH2O)
        if (con_evap_rate > 0.0_dp) then
          con%in_other = con%in_other + con_evap_rate
        else
          con%out_other = con%out_other + (-1.0_dp)*con_evap_rate
        endif
      enddo
    endif
    
    if (dat%water_cond) then
      if (dat%fix_water_in_trop) then
        i = var%trop_ind+1
      else
        i = 1
      endif
      do j = i,var%nz

        rh = max(wrk%usol(dat%LH2O,j)/wrk%H2O_sat_mix(j),small_real)

        if (rh > var%H2O_cond_params%RHc) then

          cond_rate0 = var%H2O_cond_params%k_cond*(var%edd(j)/wrk%scale_height(j)**2.0_dp)
          cond_rate = damp_condensation_rate(cond_rate0, &
                                             var%H2O_cond_params%RHc, &
                                             (1.0_dp + var%H2O_cond_params%smooth_factor)*var%H2O_cond_params%RHc, &
                                             rh)
                 
          ! Rate H2O gas is destroyed (1/s)
          df_gas_dt = - cond_rate*wrk%usol(dat%LH2O,j)
          con%out_other = con%out_other - df_gas_dt*wrk%density(j)*var%dz(j)*dat%species_composition(kk,dat%LH2O)
           
        endif
      enddo
    endif

    ! custom rate functions
    ! NOTE, This might be wrong. The model does not seem to conserve well
    ! for these custom functions.
    do i = 1,dat%nq
      if (associated(var%rate_fcns(i)%fcn)) then; block
        real(dp) :: tmp_rate
        ! note that we use the time set in wrk%tn
        call var%rate_fcns(i)%fcn(wrk%tn, var%nz, wrk%xp) ! using wrk%xp space.
        tmp_rate = sum(wrk%xp*var%dz)*dat%species_composition(kk,i) ! atoms/cm^2/s
        if (tmp_rate > 0.0_dp) then
          con%in_other = con%in_other + tmp_rate
        else
          con%out_other = con%out_other + (-1.0_dp)*tmp_rate
        endif
      endblock; endif
    enddo

    con%net = con%in_surf + con%in_top + con%in_dist + con%in_other &
            - con%out_surf - con%out_top - con%out_rain - con%out_other
    
    con%factor = abs(con%net/maxval([con%in_surf, con%in_top, con%in_dist, con%in_other, &
                                     con%out_surf, con%out_top, con%out_rain, con%out_other]))
    
  end function
  
  module function redox_conservation(self, err) result(redox_factor)
    use photochem_enum, only: VelocityDistributedFluxBC
    class(Atmosphere), target, intent(inout) :: self
    character(:), allocatable, intent(out) :: err
    real(dp) :: redox_factor
    
    real(dp) :: surf_fluxes(self%dat%nq)
    real(dp) :: top_fluxes(self%dat%nq)
    real(dp) :: integrated_rainout(self%dat%nq)
    
    real(dp) :: oxi_in_surf, oxi_out_surf
    real(dp) :: red_in_surf, red_out_surf
    real(dp) :: oxi_in_top, oxi_out_top
    real(dp) :: red_in_top, red_out_top
    real(dp) :: oxi_in_dist
    real(dp) :: red_in_dist
    real(dp) :: oxi_out_rain
    real(dp) :: red_out_rain
    real(dp) :: oxi_in_other, oxi_out_other
    real(dp) :: red_in_other, red_out_other
    
    real(dp) :: oxi_in, oxi_out, red_in, red_out
    real(dp) :: net_redox
    
    integer :: i, j
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    call self%gas_fluxes(surf_fluxes, top_fluxes, err)
    if (allocated(err)) return
    
    oxi_in_surf = 0.0_dp
    oxi_out_surf = 0.0_dp
    red_in_surf = 0.0_dp
    red_out_surf = 0.0_dp
    oxi_in_top = 0.0_dp
    oxi_out_top = 0.0_dp
    red_in_top = 0.0_dp
    red_out_top = 0.0_dp
    oxi_in_dist = 0.0_dp
    red_in_dist = 0.0_dp
    oxi_out_rain = 0.0_dp
    red_out_rain = 0.0_dp
    oxi_in_other = 0.0_dp
    oxi_out_other = 0.0_dp
    red_in_other = 0.0_dp
    red_out_other = 0.0_dp
    
    ! All numbers will be treated as positive.
    ! Later on, we implement signs
    
    ! boundary fluxes
    do i = 1,dat%nq
      if (dat%species_redox(i) > 0) then
        
        if (surf_fluxes(i) > 0) then
          oxi_in_surf = oxi_in_surf + surf_fluxes(i)*dat%species_redox(i)
        else
          oxi_out_surf = oxi_out_surf + (-1.0_dp)*surf_fluxes(i)*dat%species_redox(i)
        endif
        
        if (top_fluxes(i) > 0) then
          oxi_out_top = oxi_out_top + top_fluxes(i)*dat%species_redox(i)
        else
          oxi_in_top = oxi_in_top + (-1.0_dp)*top_fluxes(i)*dat%species_redox(i)
        endif
        
      elseif (dat%species_redox(i) < 0) then
        
        if (surf_fluxes(i) > 0) then
          red_in_surf = red_in_surf + surf_fluxes(i)*dat%species_redox(i)*(-1.0_dp)
        else
          red_out_surf = red_out_surf + (-1.0_dp)*surf_fluxes(i)*dat%species_redox(i)*(-1.0_dp)
        endif
        
        if (top_fluxes(i) > 0) then
          red_out_top = red_out_top + top_fluxes(i)*dat%species_redox(i)*(-1.0_dp)
        else
          red_in_top = red_in_top + (-1.0_dp)*top_fluxes(i)*dat%species_redox(i)*(-1.0_dp)
        endif
        
      endif
    enddo
    
    ! distributed fluxes
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == VelocityDistributedFluxBC) then
        if (dat%species_redox(i) > 0) then
          oxi_in_dist = oxi_in_dist + var%lower_flux(i)*dat%species_redox(i)
        elseif (dat%species_redox(i) < 0) then
          red_in_dist = red_in_dist + var%lower_flux(i)*dat%species_redox(i)*(-1.0_dp)
        endif
      endif
    enddo
    
    ! rainout
    if (dat%gas_rainout) then
      integrated_rainout = 0.0_dp
      ! rhs_chem already got layer 1. So
      ! we start at layer 2
      do j = 1,var%trop_ind
        do i = 1,dat%nq
          integrated_rainout(i) = integrated_rainout(i) + &
                wrk%rainout_rates(i,j)*wrk%usol(i,j)*wrk%density(j)*var%dz(j)
        enddo
      enddo
      
      do i = 1,dat%nq
        if (dat%species_redox(i) > 0) then
          oxi_out_rain = oxi_out_rain + integrated_rainout(i)*dat%species_redox(i)
        elseif (dat%species_redox(i) < 0) then
          red_out_rain = red_out_rain + integrated_rainout(i)*dat%species_redox(i)*(-1.0_dp)
        endif
      enddo
    endif

    ! custom rate functions
    ! NOTE, This might be wrong. The model does not seem to conserve well
    ! for these custom functions.
    do i = 1,dat%nq
      if (associated(var%rate_fcns(i)%fcn)) then; block
        real(dp) :: tmp_rate
        ! note that we use the time set in wrk%tn
        call var%rate_fcns(i)%fcn(wrk%tn, var%nz, wrk%xp) ! using wrk%xp space.
        tmp_rate = sum(wrk%xp*var%dz)

        if (dat%species_redox(i) > 0) then
          if (tmp_rate > 0) then
            oxi_in_other = oxi_in_other + tmp_rate*dat%species_redox(i)
          else
            oxi_out_other = oxi_out_other + (-1.0_dp)*tmp_rate*dat%species_redox(i)
          endif
        elseif (dat%species_redox(i) < 0) then
          if (tmp_rate > 0) then
            red_in_other = red_in_other + tmp_rate*dat%species_redox(i)*(-1.0_dp)
          else
            red_out_other = red_out_other + (-1.0_dp)*tmp_rate*dat%species_redox(i)*(-1.0_dp)
          endif
        endif
      endblock; endif
    enddo
    
    ! total fluxes going in and out
    oxi_in = oxi_in_surf + oxi_in_top + oxi_in_dist + oxi_in_other
    red_in = red_in_surf + red_in_top + red_in_dist + red_in_other
    oxi_out = oxi_out_surf + oxi_out_top + oxi_out_rain + oxi_out_other
    red_out = red_out_surf + red_out_top + red_out_rain + red_out_other
    ! Net redox. Should be close to zero
    net_redox = oxi_in - red_in - oxi_out + red_out
    ! compute how close net_redox is to zero, relative to redox fluxes going in and out
    redox_factor = abs(net_redox/maxval([oxi_in, red_in, oxi_out, red_out]))

  end function
  
  module subroutine set_lower_bc(self, species, bc_type, vdep, mix, flux, height, err)
    use photochem_enum, only: MosesBC, VelocityBC, MixingRatioBC, FluxBC, VelocityDistributedFluxBC
    class(Atmosphere), intent(inout) :: self
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: bc_type
    real(dp), optional, intent(in) :: vdep
    real(dp), optional, intent(in) :: mix
    real(dp), optional, intent(in) :: flux
    real(dp), optional, intent(in) :: height
    character(:), allocatable, intent(out) :: err
    
    integer :: ind(1)
    
    
    ind = findloc(self%dat%species_names(1:self%dat%nq), trim(species))
    if (ind(1) == 0) then
      err = "Can not change boundary condntion of '"//trim(species)// &
            "' because it is not in the list of species"
      return
    endif
    if (self%dat%there_are_particles) then
      if (ind(1) <= self%dat%np .and. bc_type /= 'vdep') then
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
      self%var%lowerboundcond(ind(1)) = VelocityBC
      self%var%lower_vdep(ind(1)) = vdep
      
    elseif (bc_type == 'mix') then
      if (.not. present(mix)) then
        err = "To change boundary condition to fixed mixing"// &
              " ratio must supply the 'mix' argument"
        return
      endif
      self%var%lowerboundcond(ind(1)) = MixingRatioBC
      self%var%lower_fix_mr(ind(1)) = mix
      
    elseif (bc_type == 'flux') then
      if (.not. present(flux)) then
        err = "To change boundary condition to a surface flux"// &
              " must supply the 'flux' argument"
        return
      endif
      self%var%lowerboundcond(ind(1)) = FluxBC
      self%var%lower_flux(ind(1)) = flux

    elseif (bc_type == 'vdep + dist flux') then
      if (.not.present(vdep) .or. .not.present(flux) .or. .not.present(height)) then
        err = "To change boundary condition to deposition velocity with"// &
              " a distributed flux, must supply the 'vdep', 'flux', and 'height' arguments"
        return
      endif
      self%var%lowerboundcond(ind(1)) = VelocityDistributedFluxBC
      self%var%lower_vdep(ind(1)) = vdep
      self%var%lower_flux(ind(1)) = flux
      self%var%lower_dist_height(ind(1)) = height
      
    elseif (bc_type == 'Moses') then
      self%var%lowerboundcond(ind(1)) = MosesBC
    else
      err = "Boundary condition type '"//trim(bc_type)//"' is not a valid"// &
            " boundary condition type"
      return
    endif
    
  end subroutine
    
  module subroutine set_upper_bc(self, species, bc_type, veff, flux, err)
    use photochem_enum, only: VelocityBC, FluxBC
    use photochem_enum, only: DiffusionLimHydrogenEscape
    class(Atmosphere), intent(inout) :: self
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
  
  module subroutine set_temperature(self, temperature, trop_alt, err)
    use photochem_input, only: interp2xsdata, compute_gibbs_energy
    
    class(Atmosphere), target, intent(inout) :: self
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

  module subroutine set_press_temp_edd(self, P, T, edd, trop_p, err)
    use futils, only: interp
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars

    use photochem_const, only: small_real
    class(Atmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: P(:) !! dynes/cm^2
    real(dp), intent(in) :: T(:) !! K
    real(dp), intent(in) :: edd(:) !! cm^2/s
    real(dp), optional, intent(in) :: trop_p !! dynes/cm^2
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: log10P_in(:), T_in(:), log10edd_in(:)
    real(dp), allocatable :: P_wrk(:), log10P_wrk(:), T_in_interp(:), log10edd_new(:)

    integer :: ierr

    type(MinpackHybrd1Vars) :: mv

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

    ! Work
    allocate(P_wrk(var%nz),log10P_wrk(var%nz))
    allocate(T_in_interp(var%nz))
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

      ! x_ is temperature
      call pressure_for_T(self, self%wrk%usol, x_, P_wrk, err)
      if (allocated(err)) then
        iflag_ = -1
        return
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

  ! Give a T profile on the grid and a usol, calculate pressure
  subroutine pressure_for_T(self, usol, temperature, pressure, err)
    use photochem_enum, only: MixingRatioBC
    use photochem_eqns, only: molar_weight, press_and_den
    class(Atmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: usol(:,:)
    real(dp), intent(in) :: temperature(:)
    real(dp), intent(out) :: pressure(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: usol_new(:,:)
    real(dp), allocatable :: sum_usol_new(:), mubar_new(:)
    real(dp), allocatable :: density_new(:)

    integer :: i

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    if (size(usol,1) /= dat%nq .or. size(usol,2) /= var%nz) then
      err = '"usol" is the wrong shape'
      return
    endif

    if (size(temperature) /= var%nz) then
      err = '"temperature" is the wrong shape'
      return
    endif

    if (size(pressure) /= var%nz) then
      err = '"pressure" is the wrong shape'
      return
    endif

    allocate(usol_new(dat%nq,var%nz))
    allocate(sum_usol_new(var%nz),mubar_new(var%nz))
    allocate(density_new(var%nz))

    ! make copy of mixing ratio
    usol_new = usol

    do i = 1,dat%nq
      if (var%lowerboundcond(i) == MixingRatioBC) then
        usol_new(i,1) = var%lower_fix_mr(i)
      endif
    enddo

    !!! pressure, density and mean molcular weight
    do i = 1,var%nz
      sum_usol_new(i) = sum(usol_new(dat%ng_1:,i))
      if (sum_usol_new(i) > 1.0e0_dp) then
        err = 'Mixing ratios sum to >1.0 at some altitude (should be <=1).' // &
              ' The atmosphere is probably in a run-away state'
        return
      endif
    enddo

    do i = 1,var%nz
      call molar_weight(dat%nll, usol_new(dat%ng_1:,i), sum_usol_new(i), dat%species_mass(dat%ng_1:), dat%back_gas_mu, mubar_new(i))
    enddo
    
    call press_and_den(var%nz, temperature, var%grav, var%surface_pressure*1.0e6_dp, var%dz, &
                       mubar_new, pressure, density_new)

  end subroutine

  module subroutine set_rate_fcn(self, species, fcn, err)
    use photochem_types, only: time_dependent_rate_fcn
    class(Atmosphere), target, intent(inout) :: self
    character(*), intent(in) :: species
    procedure(time_dependent_rate_fcn), pointer :: fcn
    character(:), allocatable, intent(inout) :: err
    
    integer :: ind

    ind = findloc(self%dat%species_names(1:self%dat%nq), species, 1)
    if (ind == 0) then
      err = 'Species "'//species//'" is not in the list of species, '// &
            'or is a background or short-lived species.'
      return
    endif

    self%var%rate_fcns(ind)%fcn => fcn

  end subroutine

  function TOA_at_pressure(self, usol, TOA_pressure, err) result(top_atmos)
    use minpack_module, only: hybrd1
    use clima_useful, only: MinpackHybrd1Vars
    class(Atmosphere), target, intent(inout) :: self
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
    use photochem_enum, only: MixingRatioBC
    use futils, only: interp
    use photochem_eqns, only: vertical_grid, molar_weight, press_and_den, gravity
    use photochem_const, only: small_real
    class(Atmosphere), target, intent(inout) :: self
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
    use photochem_enum, only: MixingRatioBC
    use futils, only: interp
    use photochem_eqns, only: vertical_grid, molar_weight, press_and_den, gravity
    use photochem_const, only: small_real
    class(Atmosphere), target, intent(inout) :: self
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

    real(dp), allocatable :: sum_usol_new(:), mubar_new(:)
    real(dp), allocatable :: density_new(:)
    integer :: i, ierr

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var

    dat => self%dat
    var => self%var

    ! Allocate
    allocate(sum_usol_new(var%nz),mubar_new(var%nz))
    allocate(density_new(var%nz))

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

    ! Mixing ratios
    do i = 1,dat%nq
      call interp(var%nz, var%nz, z_new, var%z, &
                  log10(max(usol(i,:),small_real)), usol_new(i,:), ierr)
      if (ierr /= 0) then
        err = 'Subroutine interp returned an error.'
        return
      endif
    enddo
    usol_new = 10.0_dp**usol_new

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
      if (var%lowerboundcond(i) == MixingRatioBC) then
        usol_new(i,1) = var%lower_fix_mr(i)
      endif
    enddo

    ! Pressure, density and mean molecular weight
    do i = 1,var%nz
      sum_usol_new(i) = sum(usol_new(dat%ng_1:,i))
      if (sum_usol_new(i) > 1.0e0_dp) then
        err = 'Mixing ratios sum to >1.0 at some altitude (should be <=1).' // &
              ' The atmosphere is probably in a run-away state'
        return
      endif
    enddo

    do i = 1,var%nz
      call molar_weight(dat%nll, usol_new(dat%ng_1:,i), sum_usol_new(i), dat%species_mass(dat%ng_1:), dat%back_gas_mu, mubar_new(i))
    enddo
    
    call press_and_den(var%nz, temperature_new, grav_new, var%surface_pressure*1.e6_dp, dz_new, &
                       mubar_new, pressure_new, density_new)

  end subroutine

  module subroutine update_vertical_grid(self, TOA_alt, TOA_pressure, err)
    use photochem_const, only: small_real
    use futils, only: interp
    use photochem_eqns, only: vertical_grid, gravity
    use photochem_input, only: interp2particlexsdata, interp2xsdata, compute_gibbs_energy
    class(Atmosphere), target, intent(inout) :: self
    real(dp), optional, intent(in) :: TOA_alt !! cm
    real(dp), optional, intent(in) :: TOA_pressure !! dynes/cm^2
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: usol_new(:,:)
    real(dp), allocatable :: z_new(:), dz_new(:), grav_new(:)
    real(dp), allocatable :: temperature_new(:), edd_new(:), particle_radius_new(:,:)
    real(dp), allocatable :: pressure_new(:)

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

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
    
    ! Set new properties
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
  
end submodule