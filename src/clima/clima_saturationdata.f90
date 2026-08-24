#:set TYPES = ['real(dp)', 'type(dual)']
#:set NAMES = ['real', 'dual']
#:set TYPES_NAMES = list(zip(TYPES, NAMES))

module clima_saturationdata
  use clima_const, only: dp, s_str_len
  implicit none
  private
  
  public :: SaturationData

  !! This module specifies the SaturationData class, which is used
  !! to compute latents heats of sublimiation/vaporization, as well as 
  !! saturation vapor pressures above solids or liquids

  type :: SaturationData
    real(dp) :: mu !! molar mass (g/mol)
    real(dp) :: T_ref !! Reference Temperature (K)
    real(dp) :: P_ref !! Reference pressure (dynes/cm^2)
    real(dp) :: T_triple !! Triple point temperature
    real(dp) :: T_critical !! Critical temperature (K)

    !> Latent heat fit parameters for vaporization (no units)
    real(dp) :: a_v, b_v
    !> Latent heat fit parameters for sublimation (no units)
    real(dp) :: a_s, b_s
    !> Non-physical latent heat fit parameters for super-critical gas (no units)
    real(dp) :: a_c, b_c

  contains
    procedure :: latent_heat_crit
    procedure :: latent_heat_vap
    procedure :: latent_heat_sub
    procedure :: latent_heat
    procedure :: sat_pressure_crit => sat_pressure_crit_real
    procedure :: sat_pressure_vap => sat_pressure_vap_real
    procedure :: sat_pressure_sub => sat_pressure_sub_real
    procedure :: sat_pressure => sat_pressure_real
    procedure :: sat_pressure_derivative
  end type
  interface SaturationData
    procedure create_SaturationData
  end interface

  #:for NAME in ['sat_pressure_vap','sat_pressure_sub','sat_pressure','integral_fcn','sat_pressure_crit']
  interface ${NAME}$
    module procedure :: ${NAME}$_real, ${NAME}$_dual
  end interface
  #:endfor

contains

  !> Latent heat of condensation above the critical point.
  !> This is non-physical. Its to approximate transitions in atmospheres
  !> between super-critical and sub-critical gases.
  function latent_heat_crit(self, T) result(L)
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: L !! erg/g
    L = self%a_c + self%b_c*T
  end function  

  !> Latent heat of vaporization
  function latent_heat_vap(self, T) result(L)
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: L !! erg/g
    L = self%a_v + self%b_v*T
  end function

  !> Latent heat of sublimation
  function latent_heat_sub(self, T) result(L)
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: L !! erg/g
    L = self%a_s + self%b_s*T
  end function

  !> Latent heat of vaporization or sublimation
  function latent_heat(self, T) result(L)
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: L !! erg/g
    if (T >= self%T_critical) then
      L = self%latent_heat_crit(T) ! non-physical approximation
    elseif (T > self%T_triple .and. T < self%T_critical) then
      L = self%latent_heat_vap(T)
    else ! (T <= T_triple) then
      L = self%latent_heat_sub(T)
    endif
  end function

  #:for TYPE1, NAME in TYPES_NAMES
  !> Saturation pressure of a super-critical gas.
  !> This is non-physical. Its to approximate transitions in atmospheres
  !> between super-critical and sub-critical gases.
  function sat_pressure_crit_${NAME}$(self, T) result(p_sat)
    #:if NAME == 'dual'
    use differentia
    #:endif
    use clima_const, only: Rgas
    class(SaturationData), intent(inout) :: self
    ${TYPE1}$, intent(in) :: T !! K
    ${TYPE1}$ :: p_sat !! dynes/cm2
    ${TYPE1}$ :: tmp
    tmp = (integral_fcn(self%a_v, self%b_v, self%T_critical) - integral_fcn(self%a_v, self%b_v, self%T_ref)) + &
          (integral_fcn(self%a_c, self%b_c, T) - integral_fcn(self%a_c, self%b_c, self%T_critical))
    p_sat = self%P_ref*exp((self%mu/Rgas)*(tmp))
  end function
  
  !> Saturation pressure over liquid
  function sat_pressure_vap_${NAME}$(self, T) result(p_sat)
    #:if NAME == 'dual'
    use differentia
    #:endif
    use clima_const, only: Rgas
    class(SaturationData), intent(inout) :: self
    ${TYPE1}$, intent(in) :: T !! K
    ${TYPE1}$ :: p_sat !! dynes/cm2
    ${TYPE1}$ :: tmp
    tmp = integral_fcn(self%a_v, self%b_v, T) - integral_fcn(self%a_v, self%b_v, self%T_ref)
    p_sat = self%P_ref*exp((self%mu/Rgas)*(tmp))
  end function

  !> Saturation pressure over solid
  function sat_pressure_sub_${NAME}$(self, T) result(p_sat)
    #:if NAME == 'dual'
    use differentia
    #:endif
    use clima_const, only: Rgas
    class(SaturationData), intent(inout) :: self
    ${TYPE1}$, intent(in) :: T !! K
    ${TYPE1}$ :: p_sat !! dynes/cm2
    ${TYPE1}$ :: tmp
    tmp = (integral_fcn(self%a_v, self%b_v, self%T_triple) - integral_fcn(self%a_v, self%b_v, self%T_ref)) + &
          (integral_fcn(self%a_s, self%b_s, T) - integral_fcn(self%a_s, self%b_s, self%T_triple))
    p_sat = self%P_ref*exp((self%mu/Rgas)*(tmp))
  end function

  !> Saturation pressure over liquid or solid
  function sat_pressure_${NAME}$(self, T) result(p_sat)
    #:if NAME == 'dual'
    use differentia
    #:endif
    class(SaturationData), intent(inout) :: self
    ${TYPE1}$, intent(in) :: T !! K
    ${TYPE1}$ :: p_sat !! dynes/cm2
    if (T >= self%T_critical) then
      p_sat = sat_pressure_crit(self, T) ! non-physical
    elseif (T > self%T_triple .and. T < self%T_critical) then
      p_sat = sat_pressure_vap(self, T)
    else ! (T <= T_triple) then
      p_sat = sat_pressure_sub(self, T)
    endif
  end function

  !> This is $\int L/T^2 dT$
  function integral_fcn_${NAME}$(A, B, T) result(res)
    #:if NAME == 'dual'
    use differentia
    #:endif
    real(dp), intent(in) :: A, B
    ${TYPE1}$, intent(in) :: T !! K
    ${TYPE1}$ :: res
    res = -A/T + B*log(T)
  end function
  #:endfor

  !> Compute the derivative of the SVP function.
  function sat_pressure_derivative(self, T) result(dPdT)
    use differentia, only: derivative
    class(SaturationData), intent(inout) :: self
    real(dp), intent(in) :: T !! K
    real(dp) :: dPdT !! dP/dT where P is dynes/cm^2 and T is K
    real(dp) :: p_sat
    call derivative(fcn, T, p_sat, dPdT)
  contains
    function fcn(x_) result(res_)
      use differentia, only: dual
      type(dual), intent(in) :: x_
      type(dual) :: res_
      res_ = sat_pressure_dual(self, x_)
    end function
  end function

  function create_SaturationData(s, name, filename, err) result(sat)
    use fortran_yaml_c_types, only: type_dictionary, type_error
    type(type_dictionary), intent(in) :: s
    character(*), intent(in) :: name
    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: err
    type(SaturationData) :: sat

    type(type_error), allocatable :: io_err
    class(type_dictionary), pointer :: tmpdict
    character(len=:), allocatable :: model
    
    model = s%get_string("model",error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (model /= "LinearLatentHeat") then
      err = 'Saturation "model" must be "LinearLatentHeat" for species "'//name//'" in '//filename
      if (model == 'Wogan') then
        err = err // '. Saturation model "Wogan" is no longer supported. *DO NOT* try to use'// &
                     ' "Wogan" model A, B, and C fitting parameters.'
      endif
      return
    endif

    tmpdict => s%get_dictionary("parameters",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%mu = tmpdict%get_real('mu',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (sat%mu <= 0.0_dp) then
      err = 'Saturation "mu" must be positive for species "'//name//'" in "'//filename//'"'
      return
    endif

    sat%T_ref = tmpdict%get_real('T-ref',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (sat%T_ref <= 0.0_dp) then
      err = 'Saturation "T-ref" must be positive for species "'//name//'" in "'//filename//'"'
      return
    endif

    sat%P_ref = tmpdict%get_real('P-ref',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (sat%P_ref <= 0.0_dp) then
      err = 'Saturation "P-ref" must be positive for species "'//name//'" in "'//filename//'"'
      return
    endif

    sat%T_triple = tmpdict%get_real('T-triple',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (sat%T_triple <= 0.0_dp) then
      err = 'Saturation "T_triple" must be positive for species "'//name//'" in "'//filename//'"'
      return
    endif
    if (sat%T_ref <= sat%T_triple) then
      err = 'Saturation "T-ref" must be bigger than "T-triple" for species "'//name//'" in "'//filename//'"'
      return
    endif

    sat%T_critical = tmpdict%get_real('T-critical',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (sat%T_critical <= 0.0_dp) then
      err = 'Saturation "T_critical" must be positive for species "'//name//'" in "'//filename//'"'
      return
    endif
    if (sat%T_ref >= sat%T_critical) then
      err = 'Saturation "T-ref" must be less than "T-critical" for species "'//name//'" in "'//filename//'"'
      return
    endif

    tmpdict => s%get_dictionary("vaporization",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%a_v = tmpdict%get_real('a',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%b_v = tmpdict%get_real('b',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    tmpdict => s%get_dictionary("sublimation",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%a_s = tmpdict%get_real('a',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%b_s = tmpdict%get_real('b',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    ! non-physical parameters to approximate atmospheric transitions from super-critical
    ! to sub-critical states.
    tmpdict => s%get_dictionary("super-critical",.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%a_c = tmpdict%get_real('a',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    sat%b_c = tmpdict%get_real('b',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

    ! test evaluations
    block
      real(dp) :: L(3), P(3)
      L(1) = sat%latent_heat(sat%T_triple-1.0e-8_dp)
      L(2) = sat%latent_heat(sat%T_triple+1.0e-8_dp)
      L(3) = sat%latent_heat(sat%T_critical-1.0e-8_dp)
      if (any(L <= 0.0_dp) .or. any(L /= L)) then
        err = 'Problem with saturation data for species "'//name//'" in "'//filename//'". '// &
              'Computed latent heats are negative or nan.'
        return
      endif

      P(1) = sat%sat_pressure(sat%T_triple-1.0e-8_dp)
      P(2) = sat%sat_pressure(sat%T_triple+1.0e-8_dp)
      P(3) = sat%sat_pressure(sat%T_critical-1.0e-8_dp)
      if (any(P <= 0.0_dp) .or. any(P /= P)) then
        err = 'Problem with saturation data for species "'//name//'" in "'//filename//'". '// &
              'Computed saturation vapor pressures are negative or nan.'
        return
      endif

      P(2) = sat%sat_pressure_derivative(sat%T_triple+1.0e-8_dp)

    endblock

  end function

end module