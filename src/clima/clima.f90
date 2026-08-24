
module clima
  ! version
  use clima_version, only: version

  ! constants that matter
  use clima_const, only: dp, s_str_len 

  ! Callback functions
  use clima_eqns, only: temp_dependent_albedo_fcn, ocean_solubility_fcn

  ! Radiative transfer classes
  use clima_radtran, only: Radtran

  ! Adiabatic climate model
  use clima_adiabat, only: AdiabatClimate 

end module
