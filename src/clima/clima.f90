
module clima
  ! constants that matter
  use clima_const, only: dp, s_str_len 

  ! Callback functions
  use clima_eqns, only: temp_dependent_albedo_fcn, ocean_solubility_fcn

  ! Radiative transfer classes
  use clima_radtran, only: Radtran, TAU_MODE_DZ, TAU_MODE_DLOG10P

  ! Adiabatic climate model
  use clima_adiabat, only: AdiabatClimate 

end module
