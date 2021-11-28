submodule(photochem_atmosphere) photochem_atmosphere_utils
  implicit none
  
  ! Contains routines utility routines for returning or saving 
  ! model output
  
contains
  
  module subroutine out2atmosphere_txt(self, filename, overwrite, clip, err)
    class(Atmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    logical, intent(in) :: overwrite, clip
    character(len=1024), intent(out) :: err
    
    character(len=100) :: tmp
    integer :: io, i, j
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    err = ""
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    ! update wrk variables
    call self%prep_atmosphere(wrk%usol, err)
    if (len_trim(err) /= 0) return
    
    if (overwrite) then
      open(1, file=filename, form='formatted', status='replace', iostat=io)
      if (io /= 0) then
        err = "Unable to overwrite file "//trim(filename)
        return
      endif
    else
      open(1, file=filename, form='formatted', status='new', iostat=io)
      if (io /= 0) then
        err = "Unable to create file "//trim(filename)//" because it already exists"
        return
      endif
    endif
    
    tmp = 'alt'
    write(unit=1,fmt="(3x,a27)",advance='no') tmp
    tmp = 'press'
    write(unit=1,fmt="(a27)",advance='no') tmp
    tmp = 'den'
    write(unit=1,fmt="(a27)",advance='no') tmp
    tmp = 'temp'
    write(unit=1,fmt="(a27)",advance='no') tmp
    tmp = 'eddy'
    write(unit=1,fmt="(a27)",advance='no') tmp
    do j = 1,dat%nq
      tmp = dat%species_names(j)
      write(unit=1,fmt="(a27)",advance='no') tmp
    enddo
    if (dat%there_are_particles) then
      do j = 1,dat%npq
        tmp = trim(dat%species_names(j))//"_r"
        write(unit=1,fmt="(a27)",advance='no') tmp
      enddo
    endif
    
    do i = 1,var%nz
      write(1,*)
      write(unit=1,fmt="(es27.17e3)",advance='no') var%z(i)/1.d5
      write(unit=1,fmt="(es27.17e3)",advance='no') wrk%pressure(i)/1.d6
      write(unit=1,fmt="(es27.17e3)",advance='no') wrk%density(i)
      write(unit=1,fmt="(es27.17e3)",advance='no') var%temperature(i)
      write(unit=1,fmt="(es27.17e3)",advance='no') var%edd(i)
      do j = 1,dat%nq
        if (clip) then
          write(unit=1,fmt="(es27.17e3)",advance='no') max(wrk%usol(j,i),1.d-40)
        else
          write(unit=1,fmt="(es27.17e3)",advance='no') wrk%usol(j,i)
        endif
      enddo
      if (dat%there_are_particles) then
        do j = 1,dat%npq
          write(unit=1,fmt="(es27.17e3)",advance='no') var%particle_radius(j,i)
        enddo
      endif
    enddo
    
    close(1)
    
  end subroutine
  
  module subroutine out2in(self, err)
    class(Atmosphere), intent(inout) :: self
    character(len=err_len), intent(out) :: err
    err = ''
    
    if (self%var%at_photo_equilibrium) then
      self%var%usol_init = self%var%usol_out
      self%var%no_water_profile = .false.
    else
      err = "Can not set output to input without first converging to photochemical equilibrium."
      return
    endif
  end subroutine
  
  module subroutine gas_fluxes(self, surf_fluxes, top_fluxes, err)
    class(Atmosphere), target, intent(inout) :: self
    real(real_kind), intent(out) :: surf_fluxes(:)
    real(real_kind), intent(out) :: top_fluxes(:)
    character(len=err_len), intent(out) :: err
  
    real(real_kind) :: rhs(self%var%neqs)  
    real(real_kind) :: diffusive_production
    real(real_kind) :: chemical_production
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
  
    integer :: i
    
    err = ""
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    if (size(surf_fluxes) /= dat%nq .or. size(top_fluxes) /= dat%nq) then
      err = "Input fluxes to gas_fluxes has the wrong dimensions"
      return
    endif
  
    call self%right_hand_side_chem(wrk%usol, rhs, err)
    if (len_trim(err) /= 0) return
    
    ! surface flux is molecules required to sustain the lower boundary
    ! chemical production + diffusion production = total change in lower cell    
    do i = 1,dat%nq
      diffusive_production = (wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                              - wrk%DU(i,1)*wrk%usol(i,1) + wrk%ADU(i,1)*wrk%usol(i,1)) &
                              *wrk%density(1)*var%dz(1)
      chemical_production = rhs(i)*wrk%density(1)*var%dz(1)
      surf_fluxes(i) = -(diffusive_production + chemical_production)
      ! We don't count chemical production for water
      if (dat%fix_water_in_trop) then
        if (i == dat%LH2O) then
          surf_fluxes(i) = - diffusive_production
        endif
      endif
    enddo
    
    ! fluxes going into or out of the top of the atmosphere.
    do i = 1,dat%nq
      diffusive_production = &
         (- wrk%DL(i,var%nz)*wrk%usol(i,var%nz) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz) &
          + wrk%DL(i,var%nz)*wrk%usol(i,var%nz-1) + wrk%ADL(i,var%nz)*wrk%usol(i,var%nz-1)) &
          *wrk%density(var%nz)*var%dz(var%nz)
    
      chemical_production = rhs(i + (var%nz-1)*dat%nq)*wrk%density(var%nz)*var%dz(var%nz)
      top_fluxes(i) = diffusive_production + chemical_production
    enddo
    
  end subroutine
  
  module function atom_conservation(self, atom, err) result(con)
    use photochem_types, only: AtomConservation
    class(Atmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: atom
    character(len=err_len), intent(out) :: err
    type(AtomConservation) :: con
    
    real(real_kind) :: surf_fluxes(self%dat%nq)
    real(real_kind) :: top_fluxes(self%dat%nq)
    real(real_kind) :: integrated_rainout(self%dat%nq)
    
    integer :: ind(1), i, j, kk
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    err = ""
    
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
    
    if (dat%fix_water_in_trop) then
      if (atom == "H" .or. atom == "O") then
        err = "Atom "//trim(atom)//" is not conserved in this model because"// &
              " H2O is fixed in the troposphere"
        return
      endif
    endif
    
    call self%gas_fluxes(surf_fluxes, top_fluxes, err)
    if (len_trim(err) /= 0) return
    
    con%in_surf = 0
    con%in_top = 0
    con%in_dist = 0
    con%out_surf = 0
    con%out_top = 0
    con%out_rain = 0
    
    ! Upper and lower boundary
    do i = 1,dat%nq
      if (surf_fluxes(i) > 0) then
        con%in_surf = con%in_surf + surf_fluxes(i)*dat%species_composition(kk,i)
      else
        con%out_surf = con%out_surf + (-1.d0)*surf_fluxes(i)*dat%species_composition(kk,i)
      endif
      if (top_fluxes(i) > 0) then
        con%out_top = con%out_top + top_fluxes(i)*dat%species_composition(kk,i)
      else
        con%in_top = con%in_top + (-1.d0)*top_fluxes(i)*dat%species_composition(kk,i)
      endif
    enddo
    
    ! distributed fluxes
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == 3) then
        con%in_dist = con%in_dist + var%lower_flux(i)*dat%species_composition(kk,i)
      endif
    enddo
    
    ! rainout
    if (var%gas_rainout) then
      integrated_rainout = 0.d0
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
    
    con%net = con%in_surf + con%in_top + con%in_dist &
            - con%out_surf - con%out_top - con%out_rain
    
    con%factor = abs(con%net/maxval([con%in_surf, con%in_top, con%in_dist, &
                                     con%out_surf, con%out_top, con%out_rain]))
    
  end function
  
  module function redox_conservation(self, err) result(redox_factor)
    class(Atmosphere), target, intent(inout) :: self
    character(len=err_len), intent(out) :: err
    real(real_kind) :: redox_factor
    
    real(real_kind) :: surf_fluxes(self%dat%nq)
    real(real_kind) :: top_fluxes(self%dat%nq)
    real(real_kind) :: integrated_rainout(self%dat%nq)
    
    real(real_kind) :: oxi_in_surf, oxi_out_surf
    real(real_kind) :: red_in_surf, red_out_surf
    real(real_kind) :: oxi_in_top, oxi_out_top
    real(real_kind) :: red_in_top, red_out_top
    real(real_kind) :: oxi_in_dist
    real(real_kind) :: red_in_dist
    real(real_kind) :: oxi_out_rain
    real(real_kind) :: red_out_rain
    
    real(real_kind) :: oxi_in, oxi_out, red_in, red_out
    real(real_kind) :: net_redox
    
    integer :: i, j
    
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk
    
    err = ""
    
    dat => self%dat
    var => self%var
    wrk => self%wrk
    
    call self%gas_fluxes(surf_fluxes, top_fluxes, err)
    if (len_trim(err) /= 0) return
    
    oxi_in_surf = 0.d0
    oxi_out_surf = 0.d0
    red_in_surf = 0.d0
    red_out_surf = 0.d0
    oxi_in_top = 0.d0
    oxi_out_top = 0.d0
    red_in_top = 0.d0
    red_out_top = 0.d0
    oxi_in_dist = 0.d0
    red_in_dist = 0.d0
    oxi_out_rain = 0.d0
    red_out_rain = 0.d0
    
    ! All numbers will be treated as positive.
    ! Later on, we implement signs
    
    ! boundary fluxes
    do i = 1,dat%nq
      if (dat%species_redox(i) > 0) then
        
        if (surf_fluxes(i) > 0) then
          oxi_in_surf = oxi_in_surf + surf_fluxes(i)*dat%species_redox(i)
        else
          oxi_out_surf = oxi_out_surf + (-1.d0)*surf_fluxes(i)*dat%species_redox(i)
        endif
        
        if (top_fluxes(i) > 0) then
          oxi_out_top = oxi_out_top + top_fluxes(i)*dat%species_redox(i)
        else
          oxi_in_top = oxi_in_top + (-1.d0)*top_fluxes(i)*dat%species_redox(i)
        endif
        
      elseif (dat%species_redox(i) < 0) then
        
        if (surf_fluxes(i) > 0) then
          red_in_surf = red_in_surf + surf_fluxes(i)*dat%species_redox(i)*(-1.d0)
        else
          red_out_surf = red_out_surf + (-1.d0)*surf_fluxes(i)*dat%species_redox(i)*(-1.d0)
        endif
        
        if (top_fluxes(i) > 0) then
          red_out_top = red_out_top + top_fluxes(i)*dat%species_redox(i)*(-1.d0)
        else
          red_in_top = red_in_top + (-1.d0)*top_fluxes(i)*dat%species_redox(i)*(-1.d0)
        endif
        
      endif
    enddo
    
    ! distributed fluxes
    do i = 1,dat%nq
      if (var%lowerboundcond(i) == 3) then
        if (dat%species_redox(i) > 0) then
          oxi_in_dist = oxi_in_dist + var%lower_flux(i)*dat%species_redox(i)
        elseif (dat%species_redox(i) < 0) then
          red_in_dist = red_in_dist + var%lower_flux(i)*dat%species_redox(i)*(-1.d0)
        endif
      endif
    enddo
    
    ! rainout
    if (var%gas_rainout) then
      integrated_rainout = 0.d0
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
          red_out_rain = red_out_rain + integrated_rainout(i)*dat%species_redox(i)*(-1.d0)
        endif
      enddo
    endif
    
    ! total fluxes going in and out
    oxi_in = oxi_in_surf + oxi_in_top + oxi_in_dist
    red_in = red_in_surf + red_in_top + red_in_dist
    oxi_out = oxi_out_surf + oxi_out_top + oxi_out_rain
    red_out = red_out_surf + red_out_top + red_out_rain
    ! Net redox. Should be close to zero
    net_redox = oxi_in - red_in - oxi_out + red_out
    ! compute how close net_redox is to zero, relative to redox fluxes going in and out
    redox_factor = abs(net_redox/maxval([oxi_in, red_in, oxi_out, red_out]))

  end function
  
  module subroutine set_lower_bc(self, species, bc_type, vdep, mix, flux, height, err)
    class(Atmosphere), intent(inout) :: self
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: bc_type
    real(real_kind), optional, intent(in) :: vdep
    real(real_kind), optional, intent(in) :: mix
    real(real_kind), optional, intent(in) :: flux
    real(real_kind), optional, intent(in) :: height
    character(len=err_len), intent(out) :: err
    
    integer :: ind(1)
    
    err = ""
    
    ind = findloc(self%dat%species_names(1:self%dat%nq), trim(species))
    if (ind(1) == 0) then
      err = "Can not change boundary condntion of '"//trim(species)// &
            "' because it is not in the list of species"
      return
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
      self%var%lowerboundcond(ind(1)) = 0
      self%var%lower_vdep(ind(1)) = vdep
      
    elseif (bc_type == 'mix') then
      if (.not. present(mix)) then
        err = "To change boundary condition to fixed mixing"// &
              " ratio must supply the 'mix' argument"
        return
      endif
      self%var%lowerboundcond(ind(1)) = 1
      self%var%lower_fix_mr(ind(1)) = mix
      
    elseif (bc_type == 'flux') then
      if (.not. present(flux)) then
        err = "To change boundary condition to a surface flux"// &
              " must supply the 'flux' argument"
        return
      endif
      self%var%lowerboundcond(ind(1)) = 2
      self%var%lower_flux(ind(1)) = flux

    elseif (bc_type == 'vdep + dist flux') then
      if (.not.present(vdep) .or. .not.present(flux) .or. .not.present(height)) then
        err = "To change boundary condition to deposition velocity with"// &
              " a distributed flux, must supply the 'vdep', 'flux', and 'height' arguments"
        return
      endif
      self%var%lowerboundcond(ind(1)) = 3
      self%var%lower_vdep(ind(1)) = vdep
      self%var%lower_flux(ind(1)) = flux
      self%var%lower_dist_height(ind(1)) = height
      
    elseif (bc_type == 'Moses') then
      self%var%lowerboundcond(ind(1)) = -1
    else
      err = "Boundary condition type '"//trim(bc_type)//"' is not a valid"// &
            " boundary condition type"
      return
    endif
    
  end subroutine
    
  module subroutine set_upper_bc(self, species, bc_type, veff, flux, err)
    class(Atmosphere), intent(inout) :: self
    character(len=*), intent(in) :: species
    character(len=*), intent(in) :: bc_type
    real(real_kind), optional, intent(in) :: veff
    real(real_kind), optional, intent(in) :: flux
    character(len=err_len), intent(out) :: err
    
    integer :: ind(1)
    
    err = ""
    
    ind = findloc(self%dat%species_names(1:self%dat%nq), trim(species))
    if (ind(1) == 0) then
      err = "Can not change boundary condition of '"//trim(species)// &
            "' because it is not in the list of species"
      return
    endif
    
    if (self%dat%diff_H_escape) then
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
      self%var%lowerboundcond(ind(1)) = 0
      self%var%lower_vdep(ind(1)) = veff

    elseif (bc_type == 'flux') then
      if (.not. present(flux)) then
        err = "To change boundary condition to a flux"// &
              " must supply the 'flux' argument"
        return
      endif
      self%var%lowerboundcond(ind(1)) = 2
      self%var%lower_flux(ind(1)) = flux

    else
      err = "Boundary condition type '"//trim(bc_type)//"' is not a valid"// &
            " upper boundary condition type"
      return
    endif
  
  end subroutine
  
end submodule