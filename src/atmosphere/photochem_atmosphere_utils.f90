submodule(photochem_atmosphere) photochem_atmosphere_utils
  implicit none
  
  ! Contains routines utility routines for returning or saving 
  ! model output
  
contains
  
  subroutine out2atmosphere_txt(self, filename, overwrite, clip, err)
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
    
    if (.not.var%at_photo_equilibrium) then
      err = "Can not write an output atmosphere until photochemical equilibrium is achieved."
      return
    endif
    
    ! update wrk variables
    call self%prep_atmosphere(var%usol_out, err)
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
          write(unit=1,fmt="(es27.17e3)",advance='no') max(var%usol_out(j,i),1.d-40)
        else
          write(unit=1,fmt="(es27.17e3)",advance='no') var%usol_out(j,i)
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
  
  subroutine out2in(self, err)
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
  
  module subroutine surface_fluxes(self, fluxes, err)
    class(Atmosphere), target, intent(inout) :: self
    real(real_kind), intent(out) :: fluxes(:)
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
    
    if (size(fluxes) /= dat%nq) then
      err = "Input fluxes to surface_fluxes has the wrong dimensions"
      return
    endif
    
    if (.not. var%at_photo_equilibrium) then
      err = "Must integrate to photochemical equilibrium before calculating surface fluxes"
      return
    endif
  
    call self%right_hand_side_chem(var%usol_out, rhs, err)
    if (len_trim(err) /= 0) return
  
    ! surface flux is molecules required to sustain the lower boundary
    ! chemical production + diffusion production = total change in lower cell    
    do i = 1,dat%nq
      diffusive_production = (wrk%DU(i,1)*wrk%usol(i,2) + wrk%ADU(i,1)*wrk%usol(i,2) &
                              - wrk%DU(i,1)*wrk%usol(i,1)) &
                              *wrk%density(1)*var%dz(1)
      chemical_production = rhs(i)*wrk%density(1)*var%dz(1)
      fluxes(i) = -(diffusive_production + chemical_production)
      ! We don't count chemical production for water
      if (dat%fix_water_in_trop) then
        if (i == dat%LH2O) then
          fluxes(i) = - diffusive_production
        endif
      endif
    enddo
  
  end subroutine
  
  module subroutine change_lower_bc(self, species, bc_type, vdep, mix, flux, height, err)
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
    
  module subroutine change_upper_bc(self, species, bc_type, veff, flux, err)
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