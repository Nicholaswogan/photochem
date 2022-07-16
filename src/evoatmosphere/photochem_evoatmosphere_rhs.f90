
submodule(photochem_evoatmosphere) photochem_evoatmosphere_rhs
  implicit none


contains

  module subroutine prep_all_evo_gas(self, dsol_in, err)

    use photochem_enum, only: MixingRatioBC
    use photochem_const, only: pi, k_boltz, N_avo, small_real

    class(EvoAtmosphere), target, intent(inout) :: self
    real(dp), intent(in) :: dsol_in(:,:)
    character(:), allocatable, intent(out) :: err

    type(PhotochemData), pointer :: dat
    type(PhotochemVars), pointer :: var
    type(PhotochemWrkEvo), pointer :: wrk
    integer :: i, j, k

    dat => self%dat
    var => self%var
    wrk => self%wrk

    ! make a copy of the densities.
    do j = 1,var%nz
      do i = 1,dat%nq
        if (dsol_in(i,j) < 0.0_dp) then
          wrk%dsol(i,j) = min(dsol_in(i,j),-small_real)
        else
          wrk%dsol(i,j) = max(dsol_in(i,j), small_real)
        endif
      enddo
    enddo

    wrk%upper_veff_copy = var%upper_veff
    wrk%lower_vdep_copy = var%lower_vdep

    ! total density and mixing ratio
    do j = 1,var%nz
      wrk%density(j) = sum(wrk%dsol(:,j))
      wrk%usol(:,j) = wrk%dsol(:,j)/wrk%density(j)
    enddo

    do i = 1,dat%nq
      if (var%lowerboundcond(i) == MixingRatioBC) then
        wrk%dsol(i,1) = var%lower_fix_mr(i)*wrk%density(1)
      endif
    enddo
    

  end subroutine


end submodule


