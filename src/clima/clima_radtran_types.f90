
module clima_radtran_types
  use iso_c_binding
  use clima_const, only: dp, s_str_len
  use linear_interpolation_module, only: linear_interp_1d, linear_interp_2d
  implicit none
  public

  real(dp), parameter :: max_w0 = 0.99999_dp
  real(dp), parameter :: max_gt = 0.999999_dp
  real(dp), parameter :: tau_min = 1.0e-20_dp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Optical Properties  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Two main different types: `Ktable` and `Xsection`.
  ! **`Ktable`** is a k-distribution table representative of line opacities.
  ! **`Xsection`** is a continuum opacity source There are 4 types (CIA, Rayleigh,
  ! Photolysis, Absorption), and can be grids, which are interpolated over 
  ! constant (no interpolation), or T or T and log10P.
  
  type :: Ktable
    integer :: sp_ind
    integer :: ngauss
    integer :: npress
    integer :: ntemp
    integer :: nwav
    real(dp), allocatable :: weight_e(:) ! (ngauss+1) bin edges
    real(dp), allocatable :: weights(:) ! (ngauss)
    real(dp), allocatable :: log10P(:) ! (nwav) log10(bars)
    real(dp) :: log10P_min
    real(dp) :: log10P_max
    real(dp), allocatable :: temp(:) ! (ntemp) kelvin
    real(dp) :: T_min
    real(dp) :: T_max
    type(linear_interp_2d), allocatable :: log10k(:,:) ! (ngauss, nwav)
  end type
  
  enum, bind(c)
    enumerator :: CIAXsection, RayleighXsection, AbsorptionXsection, PhotolysisXsection
  end enum
  
  type :: Xsection
    integer :: xs_type ! see enum above.
    integer :: dim ! 0, 1
    integer, allocatable :: sp_ind(:) ! (1 or 2). species indexes
    integer, allocatable :: rxn ! photolysis reaction number (only for PhotolysisXsection)
    integer, allocatable :: ntemp ! number of temperatures (only for xs_dim = 1)
    real(dp), allocatable :: temp(:) ! (ntemp) Kelvin
    real(dp), allocatable :: T_min
    real(dp), allocatable :: T_max
    real(dp), allocatable :: xs_0d(:) ! (nw) 
    type(linear_interp_1d), allocatable :: log10_xs_1d(:) ! (nw) 
  end type

  type :: ParticleXsection
    character(:), allocatable :: dat_name
    integer :: p_ind
    integer :: nrad
    real(dp), allocatable :: radii(:) ! cm
    real(dp) :: r_min, r_max
    type(linear_interp_1d), allocatable :: w0(:) ! (nw)
    type(linear_interp_1d), allocatable :: qext(:) ! (nw)
    type(linear_interp_1d), allocatable :: gt(:) ! (nw)
  end type
  
  type :: WaterContinuum
    character(:), allocatable :: model
    integer :: LH2O ! index of H2O
    integer :: ntemp ! number of temperatures
    real(dp), allocatable :: temp(:) ! (ntemp) Kelvin
    real(dp) :: T_min
    real(dp) :: T_max
    type(linear_interp_1d), allocatable :: log10_xs_H2O(:) ! (nw) 
    type(linear_interp_1d), allocatable :: log10_xs_foreign(:) ! (nw) 
  end type
  
  ! integer :: k_method 
  enum, bind(c)
    enumerator :: k_RandomOverlapResortRebin, k_AdaptiveEquivalentExtinction
  end enum
  
  type :: Ksettings
    ! approach to combining k-distributions
    character(:), allocatable :: k_method_name ! name
    integer :: k_method ! enum (see above)
    ! if k_method == k_RandomOverlapResortRebin
    ! then these are the weights we are re-binning to.
    integer :: nbin = -1
    real(dp), allocatable :: wxy(:) ! (nbin*nbin)
    real(dp), allocatable :: wbin(:) ! (nbin)
    real(dp), allocatable :: wbin_e(:) ! (nbin+1)
  end type
  
  type :: OpticalProperties
    integer :: nw
    real(dp), allocatable :: wavl(:)
    real(dp), allocatable :: freq(:)

    ! Copy of species and particle names
    character(s_str_len), allocatable :: species_names(:)
    character(s_str_len), allocatable :: particle_names(:)
    
    ! K-distributions (e.g. H2O)
    integer :: nk
    type(Ksettings) :: kset
    type(Ktable), allocatable :: k(:)
    ! T-dependent CIA coefficients (e.g. H2-H2)
    integer :: ncia
    type(Xsection), allocatable :: cia(:)
    ! Rayleigh Scattering
    integer :: nray
    type(Xsection), allocatable :: ray(:)
    ! Absorption cross section
    integer :: naxs
    type(Xsection), allocatable :: axs(:)
    ! Photolysis cross sections (e.g. O3 photolysis)
    integer :: npxs
    type(Xsection), allocatable :: pxs(:)

    ! Particles
    integer :: npart
    type(ParticleXsection), allocatable :: part(:)
    
    ! Water Continuum absorption. Can only be a thing if H2O is a species
    ! and if there are other gases in the atmosphere.
    type(WaterContinuum), allocatable :: cont

    ! Custom optical properties
    type(linear_interp_1d), private, allocatable :: dtau_dz_interp(:) ! nw
    type(linear_interp_1d), private, allocatable :: w0_interp(:) ! nw
    type(linear_interp_1d), private, allocatable :: g0_interp(:) ! nw

  contains
    procedure :: compute_opacity => OpticalProperties_compute_opacity
    procedure :: opacities2yaml => OpticalProperties_opacities2yaml
    procedure :: set_custom_optical_properties => OpticalProperties_set_custom_optical_properties
    procedure :: unset_custom_optical_properties => OpticalProperties_unset_custom_optical_properties
    procedure :: custom_optical_properties => OpticalProperties_custom_optical_properties
  end type
  interface
    module function create_OpticalProperties(datadir, species_names, &
                                             particle_names, sop, err) result(op)
      use clima_types, only: SettingsOpacity
      character(*), intent(in) :: datadir
      character(*), intent(in) :: species_names(:)
      character(*), intent(in) :: particle_names(:)
      type(SettingsOpacity), intent(in) :: sop
      character(:), allocatable, intent(out) :: err
      type(OpticalProperties) :: op
    end function
  end interface
  interface OpticalProperties
    module procedure :: create_OpticalProperties
  end interface
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Work arrays for optical properties  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type :: Kcoefficients
    real(dp), allocatable :: k(:,:) ! (nz, ngauss)
  end type
  
  type :: OpacityBinWork

    ! Results of interpolation
    type(Kcoefficients), allocatable :: ks(:)
    real(dp), allocatable :: cia(:,:)
    real(dp), allocatable :: pxs(:,:)
    real(dp), allocatable :: H2O(:), foreign(:)
    real(dp), allocatable :: w0p(:,:)
    real(dp), allocatable :: qextp(:,:)
    real(dp), allocatable :: gtp(:,:)

    ! Intermediate opacities
    real(dp), allocatable :: tausg(:) ! Rayleigh scattering tau
    real(dp), allocatable :: taua(:) ! Continuum absorption tau
    real(dp), allocatable :: tauc(:), tausc(:), w0c(:), g0c(:) ! Custom
    real(dp), allocatable :: tausp(:)
    real(dp), allocatable :: tausp_1(:,:)
    real(dp), allocatable :: taup(:)
    real(dp), allocatable :: taua_1(:)
    real(dp), allocatable :: tau(:)
    real(dp), allocatable :: w0(:)
    real(dp), allocatable :: gt(:)
    
    ! work arrays that are needed only if
    ! k_method == k_RandomOverlapResortRebin
    real(dp), allocatable :: tau_k(:,:) ! (nz,nbin)
    real(dp), allocatable :: tau_k1(:) ! (nbin)
    real(dp), allocatable :: tau_xy(:,:) ! (nz,nbin*nbin)
    real(dp), allocatable :: tau_xy1(:) ! (nbin*nbin)
    real(dp), allocatable :: tau_xy2(:) ! (nbin*nbin)
    real(dp), allocatable :: wxy1(:) ! (nbin*nbin)
    real(dp), allocatable :: wxy_e(:) ! (nbin*nbin+1)
    integer, allocatable :: inds(:) ! (nbin*nbin)
    ! end work arrays that are needed for k_RandomOverlapResortRebin

    ! work arrays that are only needed for if 
    ! k_method == k_AdapativeEquivalentExtinction
    real(dp), allocatable :: tau_grey(:,:) ! (nz,op%nk)
    real(dp), allocatable :: tau_grey_sum(:) ! (nz)
    integer, allocatable :: ind_major(:) ! (nz)
  end type
  interface
    module function create_OpacityBinWork(op, nz) result(wrk)
      type(OpticalProperties), target, intent(in) :: op
      integer, intent(in) :: nz
      type(OpacityBinWork) :: wrk
    end function
  end interface
  interface OpacityBinWork
    module procedure :: create_OpacityBinWork
  end interface

  type :: OpticalPropertiesWork
    ! Indicates errors (if /= 0, then error)
    integer, allocatable :: ierrs(:) ! (nw)

    ! Pressures and columns
    real(dp), allocatable :: log10P(:) ! (nz)
    real(dp), allocatable :: log10P_cgs(:) ! (nz)
    real(dp), allocatable :: cols(:,:) ! (nz,ng)
    real(dp), allocatable :: foreign_col(:) ! (nz)
    logical, allocatable :: pair_reuse(:) ! (nz) doubled-grid paired-layer cache mask

    ! Work space for each wavelength bin
    type(OpacityBinWork), allocatable :: bins(:)
  end type
  interface
    module function create_OpticalPropertiesWork(op, nz) result(opw)
      type(OpticalProperties), target, intent(in) :: op
      integer, intent(in) :: nz
      type(OpticalPropertiesWork) :: opw
    end function
  end interface
  interface OpticalPropertiesWork
    module procedure :: create_OpticalPropertiesWork
  end interface

  type :: OpticalPropertiesResult
    real(dp), allocatable :: tau(:,:,:) ! (nz, ngauss, nw)
    real(dp), allocatable :: tau_band(:,:) ! (nz, nw)
    real(dp), allocatable :: w0(:,:,:) ! (nz, ngauss, nw)
    real(dp), allocatable :: g(:,:) ! (nz,nw)
  end type
  interface
    module function create_OpticalPropertiesResult(op, nz) result(res)
      type(OpticalProperties), target, intent(in) :: op
      integer, intent(in) :: nz
      type(OpticalPropertiesResult) :: res
    end function
  end interface
  interface OpticalPropertiesResult
    module procedure :: create_OpticalPropertiesResult
  end interface

  enum, bind(c)
    enumerator :: SolarChannel, IRChannel
  end enum

  type :: RTChannel
    integer :: channel_type
    integer :: ind_start, ind_end
    integer :: nw
    real(dp), allocatable :: wavl(:)
    real(dp), allocatable :: freq(:)
  end type
  interface
    module function create_RTChannel(datadir, channel_type, wavelength_bins_file, op, err) result(rtc)
      character(*), intent(in) :: datadir
      integer, intent(in) :: channel_type
      character(:), allocatable, intent(in) :: wavelength_bins_file
      type(OpticalProperties), intent(in) :: op
      character(:), allocatable, intent(out) :: err
      type(RTChannel) :: rtc
    end function
  end interface
  interface RTChannel
    module procedure :: create_RTChannel
  end interface

  type :: RadiateBinWork
    real(dp), allocatable :: bplanck(:) ! (nz+1)
    ! All (nz+1) below
    real(dp), allocatable :: fup0(:), fup1(:), fup2(:)
    real(dp), allocatable :: fdn0(:), fdn1(:), fdn2(:)
    real(dp), allocatable :: amean0(:), amean1(:), amean2(:)
  endtype
  interface
    module function create_RadiateBinWork(nz) result(rbw)
      integer, intent(in) :: nz
      type(RadiateBinWork) :: rbw
    end function
  end interface
  interface RadiateBinWork
    module procedure :: create_RadiateBinWork
  end interface

  type :: RadiateWork
    type(RadiateBinWork), allocatable :: rbw(:) ! nw
  end type
  interface
    module function create_RadiateWork(nw, nz) result(rw)
      integer, intent(in) :: nw
      integer, intent(in) :: nz
      type(RadiateWork) :: rw
    end function
  end interface
  interface RadiateWork
    module procedure :: create_RadiateWork
  end interface
  
  ! for reading the stellar flux
  interface
    module subroutine read_stellar_flux(star_file, nw, wavl, photon_flux, err)
      character(len=*), intent(in) :: star_file
      integer, intent(in) :: nw
      real(dp), intent(in) :: wavl(nw+1)
      real(dp), intent(out) :: photon_flux(nw)
      character(:), allocatable, intent(out) :: err
    end subroutine
  end interface

contains
  
  function OpticalProperties_opacities2yaml(self) result(out)
    use clima_const, only: 
    class(OpticalProperties), intent(inout) :: self
    character(:), allocatable :: out

    character(:), allocatable :: line
    integer :: i

    out = ''

    line = '  '
    line = line//'k-method: '//self%kset%k_method_name
    out = out//line

    out = out//new_line('(a)')
    line = '  '
    line = line//'opacities:'
    out = out//line

    if (allocated(self%k)) then
      out = out//new_line('(a)')
      line = '    '
      line = line//'k-distributions: ['
      do i = 1,self%nk
        line = line//trim(self%species_names(self%k(i)%sp_ind))
        if (i /= self%nk) then
          line = line//', '
        endif
      enddo
      line = line//']'
      out = out//line
    endif

    if (allocated(self%cia)) then
      out = out//new_line('(a)')
      line = '    '
      line = line//'CIA: ['
      do i = 1,self%ncia
        line = line//trim(self%species_names(self%cia(i)%sp_ind(1)))// &
                '-'//trim(self%species_names(self%cia(i)%sp_ind(2)))
        if (i /= self%ncia) then
          line = line//', '
        endif
      enddo
      line = line//']'
      out = out//line
    endif

    if (allocated(self%ray)) then
      out = out//new_line('(a)')
      line = '    '
      line = line//'rayleigh: ['
      do i = 1,self%nray
        line = line//trim(self%species_names(self%ray(i)%sp_ind(1)))
        if (i /= self%nray) then
          line = line//', '
        endif
      enddo
      line = line//']'
      out = out//line
    endif

    if (allocated(self%pxs)) then
      out = out//new_line('(a)')
      line = '    '
      line = line//'photolysis-xs: ['
      do i = 1,self%npxs
        line = line//trim(self%species_names(self%pxs(i)%sp_ind(1)))
        if (i /= self%npxs) then
          line = line//', '
        endif
      enddo
      line = line//']'
      out = out//line
    endif

    if (allocated(self%cont)) then
      out = out//new_line('(a)')
      line = '    '
      line = line//'water-continuum: '//self%cont%model
      out = out//line
    endif

    if (allocated(self%part)) then
      out = out//new_line('(a)')
      line = '    '
      line = line//'particle-xs: ['
      do i = 1,self%npart
        line = line//'{name: '//trim(self%particle_names(self%part(i)%p_ind))
        line = line//', data: '//trim(self%part(i)%dat_name)//'}'
        if (i /= self%npart) then
          line = line//', '
        endif
      enddo
      line = line//']'
      out = out//line
    endif

  end function

  !> Sets custom optical properties
  subroutine OpticalProperties_set_custom_optical_properties(self, wv, P, dtau_dz, w0, g0, err)
    use futils, only: interp
    class(OpticalProperties), intent(inout) :: self
    real(dp), intent(in) :: wv(:) !! Array of of wavelengths in nm
    real(dp), intent(in) :: P(:) !! Array of pressures in dynes/cm^2. Must be decreasing.
    real(dp), intent(in) :: dtau_dz(:,:) !! (size(P),size(wv)), Optical depth per altitude (1/cm).
    real(dp), intent(in) :: w0(:,:) !! (size(P),size(wv)), Single scattering albedo
    real(dp), intent(in) :: g0(:,:) !! (size(P),size(wv)), Asymmetry parameter
    character(:), allocatable, intent(out) :: err
    
    real(dp), allocatable :: wv1(:), log10P(:), dtau_dz_tmp(:,:), w0_tmp(:,:), g0_tmp(:,:)
    integer :: i, j, ierr

    ! Check inputs
    if (any(wv <= 0.0_dp)) then
      err = 'All elements of `wv` must be larger than zero'
      return
    endif
    if (any(P <= 0.0_dp)) then
      err = 'All elements of `P` must be larger than zero'
      return
    endif
    if (size(P) /= size(dtau_dz,1)) then
      err = '`P` and `dtau_dz` have incompatible shapes'
      return
    endif
    if (size(wv) /= size(dtau_dz,2)) then
      err = '`wv` and `dtau_dz` have incompatible shapes'
      return
    endif
    if (size(P) /= size(w0,1)) then
      err = '`P` and `w0` have incompatible shapes'
      return
    endif
    if (size(wv) /= size(w0,2)) then
      err = '`wv` and `w0` have incompatible shapes'
      return
    endif
    if (size(P) /= size(g0,1)) then
      err = '`P` and `g0` have incompatible shapes'
      return
    endif
    if (size(wv) /= size(g0,2)) then
      err = '`wv` and `g0` have incompatible shapes'
      return
    endif

    allocate(wv1(size(self%wavl)-1))
    allocate(dtau_dz_tmp(size(P),size(self%wavl)-1))
    allocate(w0_tmp(size(P),size(self%wavl)-1))
    allocate(g0_tmp(size(P),size(self%wavl)-1))

    ! Median wavelength
    wv1 = 0.5_dp*(self%wavl(2:) + self%wavl(1:size(self%wavl)-1))

    ! At each pressure, interpolate to the wavelength grid
    do i = 1,size(P)
      j = size(P) + 1 - i
      call interp(wv1, wv, dtau_dz(i,:), dtau_dz_tmp(j,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Interpolation error in `set_custom_optical_properties`'
        return
      endif
      call interp(wv1, wv, w0(i,:), w0_tmp(j,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Interpolation error in `set_custom_optical_properties`'
        return
      endif
      call interp(wv1, wv, g0(i,:), g0_tmp(j,:), ierr=ierr)
      if (ierr /= 0) then
        err = 'Interpolation error in `set_custom_optical_properties`'
        return
      endif
    enddo

    ! Consider log10 pressure
    log10P = log10(P)
    log10P = log10P(size(log10P):1:-1)

    if (.not.allocated(self%dtau_dz_interp)) then
      allocate(self%dtau_dz_interp(self%nw))
      allocate(self%w0_interp(self%nw))
      allocate(self%g0_interp(self%nw))
    endif
   
    ! At each wavelength bin, create an interpolator for altitude
    do i = 1,size(wv1)
      call self%dtau_dz_interp(i)%initialize(log10P, dtau_dz_tmp(:,i), ierr)
      if (ierr /= 0) then
        err = 'Interpolation initialization error in `set_custom_optical_properties`'
        return
      endif
      call self%w0_interp(i)%initialize(log10P, w0_tmp(:,i), ierr)
      if (ierr /= 0) then
        err = 'Interpolation initialization error in `set_custom_optical_properties`'
        return
      endif
      call self%g0_interp(i)%initialize(log10P, g0_tmp(:,i), ierr)
      if (ierr /= 0) then
        err = 'Interpolation initialization error in `set_custom_optical_properties`'
        return
      endif
    enddo

  end subroutine

  !> Unsets custom optical properties set with `set_custom_optical_properties`.
  subroutine OpticalProperties_unset_custom_optical_properties(self)
    class(OpticalProperties), intent(inout) :: self
    if (allocated(self%dtau_dz_interp)) then
      deallocate(self%dtau_dz_interp)
      deallocate(self%w0_interp)
      deallocate(self%g0_interp)
    endif
  end subroutine

  !> Interpolates the custom optical properties
  subroutine OpticalProperties_custom_optical_properties(self, log10P, dz, l, tau, w0, g0)
    use clima_eqns, only: ten2power
    class(OpticalProperties), intent(inout) :: self
    real(dp), intent(in) :: log10P(:) !! Pressure in dynes/cm^2
    real(dp), intent(in) :: dz(:) !! Layer thickness in cm
    integer, intent(in) :: l !! Wavelength index
    real(dp), intent(out) :: tau(:) !! (nz), Optical depth.
    real(dp), intent(out) :: w0(:) !! (nz), Single scattering albedo
    real(dp), intent(out) :: g0(:) !! (nz), Asymetry parameter

    integer :: j

    if (.not.allocated(self%dtau_dz_interp)) then
      tau = tiny(0.0_dp)
      w0 = tiny(0.0_dp)
      g0 = tiny(0.0_dp)
      return
    endif

    do j = 1,size(log10P)
      call self%dtau_dz_interp(l)%evaluate(log10P(j), tau(j))
      tau(j) = tau(j)*dz(j)
      call self%w0_interp(l)%evaluate(log10P(j), w0(j))
      call self%g0_interp(l)%evaluate(log10P(j), g0(j))
    enddo

  end subroutine

  subroutine OpticalProperties_compute_opacity(self, P, T, densities, dz, pdensities, radii, opw, res, err)
    use clima_const, only: pi
    use clima_eqns, only: ten2power
    use futils, only: is_close
    class(OpticalProperties), intent(inout) :: self

    real(dp), intent(in) :: P(:) !! (nz) Pressure (bars)
    real(dp), intent(in) :: T(:) !! (nz) Temperature (K) 
    real(dp), intent(in) :: densities(:,:) !! (nz,ng) number density of each 
                                           !! molecule in each layer (molcules/cm3)
    real(dp), intent(in) :: dz(:) !! (nz) thickness of each layer (cm)
    real(dp), optional, intent(in) :: pdensities(:,:) !! (nz,np) particles/cm3
    real(dp), optional, intent(in) :: radii(:,:) !! (nz,np) cm
    class(OpticalPropertiesWork), target, intent(inout) :: opw !! Work space
    class(OpticalPropertiesResult), intent(inout) :: res !! Result
    character(:), allocatable, intent(out) :: err

    integer :: nz, ng
    integer :: l
    ! Ensure
    if (self%nk == 0) then
      err = 'There are no k-distributions, yet there must be some to compute total opacity.'
      return
    endif

    ! Pre-compute log10P and columns.
    block
    integer :: i, j
    real(dp), parameter :: pair_tol = 1.0e-12_dp
      nz = size(dz)
      ng = size(densities, 2)
      opw%log10P = log10(P)
      opw%log10P_cgs = log10(P*1.0e6_dp)
      do i = 1,ng
        opw%cols(:,i) = densities(:,i)*dz(:)
      enddo
      if (allocated(self%cont)) then
        do j = 1,nz
          opw%foreign_col(j) = 0.0_dp
          do i = 1,ng
            if (i /= self%cont%LH2O) then
              opw%foreign_col(j) = opw%foreign_col(j) + opw%cols(j,i)
            endif
          enddo
        enddo
      endif

      opw%pair_reuse = .false.
      if (mod(nz,2) == 0) then
        do j = 2,nz,2
          opw%pair_reuse(j) = .true.
          opw%pair_reuse(j) = opw%pair_reuse(j) .and. is_close(P(j), P(j-1), tol=pair_tol)
          opw%pair_reuse(j) = opw%pair_reuse(j) .and. is_close(T(j), T(j-1), tol=pair_tol)
          opw%pair_reuse(j) = opw%pair_reuse(j) .and. all(is_close(opw%cols(j,:), opw%cols(j-1,:), tol=pair_tol))
          if (present(radii) .and. self%npart > 0) then
            opw%pair_reuse(j) = opw%pair_reuse(j) .and. all(is_close(radii(j,:), radii(j-1,:), tol=pair_tol))
          endif
        enddo
      endif
    end block

    ! Zero out error indicator
    opw%ierrs = 0

    !$omp parallel private(l)
    !$omp do
    do l = 1,self%nw; block
      integer :: i, j, k, n, jj, ierr
      real(dp) :: TT, log10PP, taup_1
      type(OpacityBinWork), pointer :: wrk

      ! Get work space for a given bin
      wrk => opw%bins(l)

      ! Interpolate k-distributions
      do i = 1,self%nk
        do k = 1,self%k(i)%ngauss
          do j = 1,nz
            if (opw%pair_reuse(j)) then
              wrk%ks(i)%k(j,k) = wrk%ks(i)%k(j-1,k)
            else
              TT = min(max(T(j), self%k(i)%T_min), self%k(i)%T_max)
              log10PP = min(max(opw%log10P(j), self%k(i)%log10P_min), self%k(i)%log10P_max)
              call self%k(i)%log10k(k,l)%evaluate(log10PP, TT, wrk%ks(i)%k(j,k))
              wrk%ks(i)%k(j,k) = ten2power(wrk%ks(i)%k(j,k))
            endif
          enddo
        enddo
      enddo

      ! Interpolate CIA
      do i = 1,self%ncia
        call interpolate_Xsection(self%cia(i), l, T, wrk%cia(:,i), opw%pair_reuse)
      enddo

      ! Interpolate photolysis
      do i = 1,self%npxs
        call interpolate_Xsection(self%pxs(i), l, T, wrk%pxs(:,i), opw%pair_reuse)
      enddo

      ! Interpolate water continuum
      if (allocated(self%cont)) then
        call interpolate_WaterContinuum(self%cont, l, T, wrk%H2O, wrk%foreign, opw%pair_reuse)
      endif

      ! Interpolate particles
      do i = 1,self%npart
        call interpolate_Particle(self%part(i), l, radii, wrk%w0p(:,i), wrk%qextp(:,i), wrk%gtp(:,i), ierr, opw%pair_reuse)
        opw%ierrs(l) = opw%ierrs(l) + ierr
      enddo

      ! Rayleigh scattering
      wrk%tausg(:) = 0.0_dp
      do i = 1,self%nray
        j = self%ray(i)%sp_ind(1)
        do k = 1,nz
          n = nz+1-k
          wrk%tausg(n) = wrk%tausg(n) + self%ray(i)%xs_0d(l)*opw%cols(k,j)
        enddo
      enddo
      
      ! CIA
      wrk%taua(:) = 0.0_dp
      do i = 1,self%ncia
        j = self%cia(i)%sp_ind(1)
        jj = self%cia(i)%sp_ind(2)
        do k = 1,nz
          n = nz+1-k
          wrk%taua(n) = wrk%taua(n) + wrk%cia(k,i)*densities(k,j)*densities(k,jj)*dz(k)
        enddo
      enddo
      
      ! Photolysis/Absorption
      do i = 1,self%npxs
        j = self%pxs(i)%sp_ind(1)
        do k = 1,nz
          n = nz+1-k
          wrk%taua(n) = wrk%taua(n) + wrk%pxs(k,i)*opw%cols(k,j)
        enddo
      enddo

      ! Continuum absorption
      if (allocated(self%cont)) then
        do k = 1,nz
          n = nz+1-k
          wrk%taua(n) = wrk%taua(n) &
                       + wrk%H2O(k)*densities(k,self%cont%LH2O)*opw%cols(k,self%cont%LH2O) &
                       + wrk%foreign(k)*densities(k,self%cont%LH2O)*opw%foreign_col(k)
        enddo
      endif

      ! Custom opacity
      call self%custom_optical_properties(opw%log10P_cgs, dz, l, wrk%tauc, wrk%w0c, wrk%g0c)
      wrk%tauc = wrk%tauc(nz:1:-1)
      wrk%w0c = wrk%w0c(nz:1:-1)
      wrk%g0c = wrk%g0c(nz:1:-1)
      wrk%tausc = wrk%w0c*wrk%tauc

      ! Particles
      wrk%tausp(:) = 0.0_dp
      wrk%taup(:) = 0.0_dp
      do i = 1,self%npart
        j = self%part(i)%p_ind
        do k = 1,nz
          n = nz+1-k
          taup_1 = wrk%qextp(k,i)*pi*radii(k,j)**2.0_dp*pdensities(k,j)*dz(k)
          wrk%taup(n) = wrk%taup(n) + taup_1
          wrk%tausp_1(n,i) = wrk%w0p(k,i)*taup_1
          wrk%tausp(n) = wrk%tausp(n) + wrk%tausp_1(n,i)
        enddo
      enddo

      wrk%gt(:) = 0.0_dp
      do i = 1,self%npart
        do k = 1,nz
          n = nz+1-k
          wrk%gt(n) = wrk%gt(n) + wrk%gtp(k,i)*wrk%tausp_1(n,i)/max(tau_min, (wrk%tausp(n) + wrk%tausg(n) + wrk%tausc(n)))
        enddo
      enddo
      wrk%gt = wrk%gt + wrk%g0c*wrk%tausc/max(tau_min, (wrk%tausp + wrk%tausg + wrk%tausc))
      do k = 1,nz
        n = nz+1-k
        wrk%gt(n) = min(wrk%gt(n), max_gt)
      enddo

      if (self%kset%k_method == k_RandomOverlapResortRebin) then
        call k_rorr(self, l, opw, res)
      elseif (self%kset%k_method == k_AdaptiveEquivalentExtinction) then
        ! Will implement eventually
        opw%ierrs(l) = 1
      else
        ! This should never happen
        opw%ierrs(l) = 1
      endif

    end block; enddo
    !$omp enddo
    !$omp end parallel

    if (any(opw%ierrs /= 0)) then
      err = 'Opacity computation failed in one or more wavelength bins.'
      return
    endif

  end subroutine

  subroutine k_rorr(op, l, opw, res)
    use futils_mrgrnk, only: mrgrnk
    use futils, only: rebin
    use clima_eqns, only: weights_to_bins
    type(OpticalProperties), target, intent(in) :: op
    integer, intent(in) :: l
    type(OpticalPropertiesWork), target, intent(inout) :: opw
    type(OpticalPropertiesResult), intent(inout) :: res

    type(OpacityBinWork), pointer :: wrk
    type(Ksettings), pointer :: kset
    real(dp), pointer :: tau_k(:,:)
    real(dp), pointer :: tau_k1(:)
    real(dp), pointer :: tau_xy(:,:)
    real(dp), pointer :: tau_xy1(:)
    real(dp), pointer :: tau_xy2(:)
    real(dp), pointer :: wxy1(:)
    real(dp), pointer :: wxy_e(:)
    integer, pointer :: inds(:)

    integer :: i, j, jj, j2, j1, k, n
    integer :: nz, ngauss

    kset => op%kset
    wrk => opw%bins(l)
    tau_k => wrk%tau_k
    tau_k1 => wrk%tau_k1
    tau_xy => wrk%tau_xy
    tau_xy1 => wrk%tau_xy1
    tau_xy2 => wrk%tau_xy2
    wxy1 => wrk%wxy1
    wxy_e => wrk%wxy_e
    inds => wrk%inds
    nz = size(opw%cols, 1)

    ! First k-distribution
    j1 = op%k(1)%sp_ind
    do i = 1,kset%nbin
      tau_k(:,i) = wrk%ks(1)%k(:,i)*opw%cols(:,j1) ! cm^2/molecule * molecules/cm^2
    enddo

    ! Mix rest of k-coeff species with the first species
    ngauss = kset%nbin
    do jj = 2,op%nk
      j2 = op%k(jj)%sp_ind

      do i = 1,ngauss
        do j = 1,ngauss
          tau_xy(:,j + (i-1)*ngauss) = tau_k(:,i) + wrk%ks(jj)%k(:,j)*opw%cols(:,j2)
        enddo
      enddo

      do i = 1,nz
        if (opw%pair_reuse(i)) then
          tau_k(i,:) = tau_k(i-1,:)
        else
          ! Sort tau_xy and the weights
          do j = 1,size(tau_xy,2)
            tau_xy1(j) = tau_xy(i,j)
          enddo
          call mrgrnk(tau_xy1, inds)
          do j = 1,size(inds)
            tau_xy2(j) = tau_xy1(inds(j))
            wxy1(j) = kset%wxy(inds(j))
          enddo
          ! Rebin to smaller grid
          call weights_to_bins(wxy1, wxy_e)
          call rebin(wxy_e, tau_xy2, kset%wbin_e, tau_k1)
          do j = 1,size(tau_k1)
            tau_k(i,j) = tau_k1(j)
          enddo
        endif
      enddo

    enddo

    res%tau_band(:,l) = 0.0_dp
    do i = 1,kset%nbin

      ! tau_k(:,i) is optical depth of ith mixed k-coeff.
      ! tau_k(1,i) is lowest level. Need to reorder so that
      ! first element in array is top of the atmosphere.
      do k = 1,nz
        n = nz+1-k
        wrk%taua_1(n) = tau_k(k,i)
      enddo

      ! Sum all optical depths
      ! total = gas scattering + continumm opacities + particle absorption + k-coeff + custom opacity
      wrk%tau(:) = wrk%tausg(:) + wrk%taua(:) + wrk%taup(:) + wrk%taua_1(:) + wrk%tauc(:)
      do j = 1,nz
        if (wrk%tau(j) <= tau_min) then
          wrk%w0(j) = 0.0_dp
        else
          wrk%w0(j) = min(max_w0,(wrk%tausg(j) + wrk%tausp(j) + wrk%tausc(j))/wrk%tau(j))
        endif
      enddo

      ! Save results
      res%tau(:,i,l) = wrk%tau
      res%tau_band(:,l) = res%tau_band(:,l) + wrk%tau*kset%wbin(i)
      res%w0(:,i,l) = wrk%w0

    enddo

    ! Asymmetry parameter
    res%g(:,l) = wrk%gt
    
  end subroutine

  subroutine interpolate_Xsection(xs, l, T, res, pair_reuse)
    use clima_eqns, only: ten2power
    type(Xsection), intent(inout) :: xs
    integer, intent(in) :: l
    real(dp), intent(in) :: T(:)
    real(dp), intent(out) :: res(:)
    logical, intent(in) :: pair_reuse(:)
    
    integer :: j
    real(dp) :: val, TT

    if (xs%dim == 0) then
      do j = 1,size(T)
        res(j) = xs%xs_0d(l)
      enddo
    elseif (xs%dim == 1) then
      do j = 1,size(T)
        if (pair_reuse(j)) then
          res(j) = res(j-1)
        else
          TT = min(max(T(j), xs%T_min), xs%T_max)
          call xs%log10_xs_1d(l)%evaluate(TT, val)
          res(j) = ten2power(val)
        endif
      enddo
    endif
    
  end subroutine
  
  subroutine interpolate_WaterContinuum(cont, l, T, H2O, foreign, pair_reuse)
    use clima_eqns, only: ten2power
    
    type(WaterContinuum), intent(inout) :: cont
    integer, intent(in) :: l
    real(dp), intent(in) :: T(:)
    real(dp), intent(out) :: H2O(:)
    real(dp), intent(out) :: foreign(:)
    logical, intent(in) :: pair_reuse(:)
    
    integer :: j
    real(dp) :: val, TT
    
    do j = 1,size(T)
      if (pair_reuse(j)) then
        H2O(j) = H2O(j-1)
        foreign(j) = foreign(j-1)
      else
        TT = min(max(T(j), cont%T_min), cont%T_max)
        call cont%log10_xs_H2O(l)%evaluate(TT, val)
        H2O(j) = ten2power(val)
        call cont%log10_xs_foreign(l)%evaluate(TT, val)
        foreign(j) = ten2power(val)
      endif
    enddo
    
  end subroutine

  subroutine interpolate_Particle(part, l, radii, w0, qext, gt, ierr, pair_reuse)
    type(ParticleXsection), intent(inout) :: part
    integer, intent(in) :: l
    real(dp), intent(in) :: radii(:,:)
    real(dp), intent(out) :: w0(:)
    real(dp), intent(out) :: qext(:)
    real(dp), intent(out) :: gt(:)
    integer, intent(out) :: ierr
    logical, intent(in) :: pair_reuse(:)

    integer :: j
    real(dp) :: rp

    ierr = 0

    do j = 1,size(radii,1)
      if (pair_reuse(j)) then
        w0(j) = w0(j-1)
        qext(j) = qext(j-1)
        gt(j) = gt(j-1)
        cycle
      endif

      rp = radii(j,part%p_ind)

      ! Error indication if radius is outside bounds
      if (rp < part%r_min .or. rp > part%r_max) then
        rp = min(max(rp, part%r_min), part%r_max)
        ierr = 1
      endif

      call part%w0(l)%evaluate(rp, w0(j))
      call part%qext(l)%evaluate(rp, qext(j))
      call part%gt(l)%evaluate(rp, gt(j))
    enddo

  end subroutine

end module
