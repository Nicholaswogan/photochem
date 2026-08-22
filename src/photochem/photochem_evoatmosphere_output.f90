submodule(photochem_evoatmosphere) photochem_evoatmosphere_output
  implicit none

contains

  module subroutine out2atmosphere_txt(self, filename, number_of_decimals, overwrite, clip, err)
    class(EvoAtmosphere), target, intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: number_of_decimals
    logical, intent(in) :: overwrite
    logical, intent(in) :: clip
    character(:), allocatable, intent(out) :: err

    real(dp) :: rhs(self%var%neqs)
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrk), pointer :: wrk

    call self%require_atmosphere_initialized('out2atmosphere_txt', err)
    if (allocated(err)) return

    dat => self%dat
    var => self%var
    wrk => self%wrk

    ! Update work variables needed for output.
    call self%chemistry_right_hand_side(wrk%usol, rhs, err)
    if (allocated(err)) return

    call out2atmosphere_txt_base(dat, var, &
                                 wrk%pressure, wrk%density, wrk%densities, &
                                 wrk%molecules_per_particle, filename, &
                                 number_of_decimals, overwrite, clip, err)
    if (allocated(err)) return

  end subroutine

  subroutine out2atmosphere_txt_base(dat, var, &
                                     pressure, density, densities, molecules_per_particle, &
                                     filename, number_of_decimals, overwrite, clip, err)
    use futils, only: FileCloser
    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    real(dp), intent(in) :: pressure(:), density(:), densities(:,:), molecules_per_particle(:,:)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: number_of_decimals
    logical, intent(in) :: overwrite, clip
    character(:), allocatable, intent(out) :: err

    character(len=100) :: tmp
    integer :: io, i, j
    integer :: number_of_spaces
    character(len=10) :: number_of_decimals_str, number_of_spaces_str
    character(:), allocatable :: fmt_label, fmt_number
    real(dp) :: val, clip_value
    type(FileCloser) :: file

    if (clip) then
      clip_value = 1.0e-40_dp
    else
      clip_value = - huge(1.0_dp)
    endif

    if (overwrite) then
      open(1, file=filename, form='formatted', status='replace', iostat=io)
      file%unit = 1
      if (io /= 0) then
        err = "Unable to overwrite file "//trim(filename)
        return
      endif
    else
      open(1, file=filename, form='formatted', status='new', iostat=io)
      file%unit = 1
      if (io /= 0) then
        err = "Unable to create file "//trim(filename)//" because it already exists"
        return
      endif
    endif

    ! number of decimals must be reasonable
    if (number_of_decimals < 2 .or. number_of_decimals > 17) then
      err = '"number_of_decimals" should be between 1 and 17.'
      return
    endif
    number_of_spaces = number_of_decimals + 9
    ! make sure number of spaces works with the length of species names
    do i=1,dat%nsp
      number_of_spaces = max(number_of_spaces,len_trim(dat%species_names(i)) + 3)
    enddo
    write(number_of_decimals_str,'(i10)') number_of_decimals
    write(number_of_spaces_str,'(i10)') number_of_spaces

    fmt_label = "(a"//trim(adjustl(number_of_spaces_str))//")"
    fmt_number = "(es"//trim(adjustl(number_of_spaces_str))//"."//trim(adjustl(number_of_decimals_str))//"e3)"

    tmp = 'alt'
    write(unit=1,fmt=fmt_label,advance='no') tmp
    tmp = 'press'
    write(unit=1,fmt=fmt_label,advance='no') tmp
    tmp = 'den'
    write(unit=1,fmt=fmt_label,advance='no') tmp
    tmp = 'temp'
    write(unit=1,fmt=fmt_label,advance='no') tmp
    tmp = 'eddy'
    write(unit=1,fmt=fmt_label,advance='no') tmp
    do j = 1,dat%nsp
      tmp = dat%species_names(j)
      write(unit=1,fmt=fmt_label,advance='no') tmp
    enddo
    if (dat%there_are_particles) then
      do j = 1,dat%npq
        tmp = trim(dat%species_names(j))//"_r"
        write(unit=1,fmt=fmt_label,advance='no') tmp
      enddo
    endif

    do i = 1,var%nz
      write(1,*)
      write(tmp,fmt=fmt_number) var%z(i)/1.e5_dp
      write(unit=1,fmt=fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt=fmt_number) pressure(i)/1.e6_dp
      write(unit=1,fmt=fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt=fmt_number) density(i)
      write(unit=1,fmt=fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt=fmt_number) var%temperature(i)
      write(unit=1,fmt=fmt_label,advance='no') adjustl(tmp)

      write(tmp,fmt=fmt_number) var%edd(i)
      write(unit=1,fmt=fmt_label,advance='no') adjustl(tmp)
      ! particles
      if (dat%there_are_particles) then
        do j = 1,dat%npq
          val = densities(j,i)*(molecules_per_particle(j,i)/density(i)) ! mixing ratio of particle
          write(tmp,fmt=fmt_number) max(val, clip_value)
          write(unit=1,fmt=fmt_label,advance='no') adjustl(tmp)
        enddo
      endif
      ! gases
      do j = dat%ng_1,dat%nsp
        val = densities(j,i)/density(i)
        write(tmp,fmt=fmt_number) max(val, clip_value)
        write(unit=1,fmt=fmt_label,advance='no') adjustl(tmp)
      enddo
      ! particle radii
      if (dat%there_are_particles) then
        do j = 1,dat%npq
          write(tmp,fmt=fmt_number) var%particle_radius(j,i)
          write(unit=1,fmt=fmt_label,advance='no') adjustl(tmp)
        enddo
      endif
    enddo

  end subroutine

end submodule
