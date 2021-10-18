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
  
end submodule