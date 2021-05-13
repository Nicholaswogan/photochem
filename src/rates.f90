


function normal_rate(A, b, Ea, T) result(k)
  real(real_kind), intent(in) :: A, b, Ea, T
  real(real_kind) :: k
  k = A * T**b * dexp(-Ea/T)
end function

function falloff_rate(kinf, Pr, F) result(k)
  real(real_kind), intent(in) :: kinf, Pr, F
  real(real_kind) :: k
  
  k = kinf * (Pr / (1.d0 + Pr)) * F
end function

function Troe_noT2(A, T1, T3, T, Pr) result(F)
  real(real_kind), intent(in) :: A, T1, T3, T, Pr
  real(real_kind) :: F
  
  real(real_kind) :: log10Fcent, f1, C, N
  
  log10Fcent = dlog10((1.d0-A)*dexp(-T/T3) + A*dexp(-T/T1))
  C = -0.4d0 - 0.67d0*log10Fcent
  N = 0.75d0 - 1.27d0*log10Fcent
  f1 = (dlog10(Pr) + C)/(N - 0.14d0*(dlog10(Pr + C)))
  F = 10.d0**((log10Fcent)/(1.d0 + f1**2.d0))
end function

function Troe_withT2(A, T1, T2, T3, T, Pr) result(F)
  real(real_kind), intent(in) :: A, T1, T2, T3, T, Pr
  real(real_kind) :: F
  
  real(real_kind) :: log10Fcent, f1, C, N
  
  log10Fcent = dlog10((1.d0-A)*dexp(-T/T3) + A*dexp(-T/T1) + dexp(-T2/T))
  C = -0.4d0 - 0.67d0*log10Fcent
  N = 0.75d0 - 1.27d0*log10Fcent
  f1 = (dlog10(Pr) + C)/(N - 0.14d0*(dlog10(Pr + C)))
  F = 10.d0**((log10Fcent)/(1.d0 + f1**2.d0))

end function


