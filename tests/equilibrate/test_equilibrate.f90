program main
  use equilibrate, only: ChemEquiAnalysis, dp, s_str_len
  use equilibrate_test_paths, only: fixture_file
  implicit none

  call test()
  call test_sonora()
  call test_memory()
  call test_nasa7()

contains

  subroutine test()
    type(ChemEquiAnalysis) :: cea, cea2
    character(:), allocatable :: err
    character(s_str_len), allocatable :: species(:)
    character(s_str_len), allocatable :: atoms(:)
    real(dp), allocatable :: X(:), correct_answer(:)
    logical :: converged
    integer :: i

    species = [ &
      'H              ', &
      'H2             ', &
      'He             ', &
      'O              ', &
      'C              ', &
      'N              ', &
      'Mg             ', &
      'Si             ', &
      'Fe             ', &
      'S              ', &
      'Al             ', &
      'Ca             ', &
      'Na             ', &
      'Ni             ', &
      'P              ', &
      'K              ', &
      'Ti             ', &
      'CO             ', &
      'OH             ', &
      'SH             ', &
      'N2             ', &
      'O2             ', &
      'SiO            ', &
      'TiO            ', &
      'SiS            ', &
      'H2O            ', &
      'C2             ', &
      'CH             ', &
      'CN             ', &
      'CS             ', &
      'SiC            ', &
      'NH             ', &
      'SiH            ', &
      'NO             ', &
      'SN             ', &
      'SiN            ', &
      'SO             ', &
      'S2             ', &
      'C2H            ', &
      'HCN            ', &
      'C2H2,acetylene ', &
      'CH4            ', &
      'AlH            ', &
      'AlOH           ', &
      'Al2O           ', &
      'CaOH           ', &
      'MgH            ', &
      'MgOH           ', &
      'PH3            ', &
      'CO2            ', &
      'TiO2           ', &
      'Si2C           ', &
      'SiO2           ', &
      'FeO            ', &
      'NH2            ', &
      'NH3            ', &
      'CH2            ', &
      'CH3            ', &
      'H2S            ', &
      'VO             ', &
      'VO2            ', &
      'NaCl           ', &
      'KCl            ', &
      'e-             ', &
      'H+             ', &
      'H-             ', &
      'Na+            ', &
      'K+             ', &
      'PH2            ', &
      'P2             ', &
      'PS             ', &
      'PO             ', &
      'P4O6           ', &
      'PH             ', &
      'V              ', &
      'VO(c)          ', &
      'VO(L)          ', &
      'MgSiO3(c)      ', &
      'SiC(c)         ', &
      'Fe(c)          ', &
      'Al2O3(c)       ', &
      'Na2S(c)        ', &
      'KCl(c)         ', &
      'Fe(L)          ', &
      'SiC(L)         ', &
      'MgSiO3(L)      ', &
      'H2O(L)         ', &
      'H2O(c)         ', &
      'TiO(c)         ', &
      'TiO(L)         ', &
      'TiO2(c)        ', &
      'TiO2(L)        ', &
      'H3PO4(c)       ', &
      'H3PO4(L)       ' &
    ]

    atoms = [ &
      'H ', &
      'He', &
      'C ', &
      'N ', &
      'O ', &
      'Na', &
      'Mg', &
      'Al', &
      'Si', &
      'P ', &
      'S ', &
      'Cl', &
      'K ', &
      'Ca', &
      'Ti', &
      'V ', &
      'Fe', &
      'Ni' &
    ]

    X = [ &
      9.207539305000000e-01_dp, &
      7.836886940000000e-02_dp, &
      2.478241000000000e-04_dp, &
      6.225060569498810e-05_dp, &
      4.509658000000000e-04_dp, &
      1.600086943532050e-06_dp, &
      3.665587420553620e-05_dp, &
      2.595000000000000e-06_dp, &
      2.979500000000000e-05_dp, &
      2.366702019976680e-07_dp, &
      1.213790073460400e-05_dp, &
      2.911679584995890e-07_dp, &
      9.866056119256769e-08_dp, &
      2.014390114292550e-06_dp, &
      8.206228043663590e-08_dp, &
      7.836886940899920e-09_dp, &
      2.911679584995890e-05_dp, &
      1.528071168062810e-06_dp &
    ]

    correct_answer = [ &
      2.097572751458766e-09_dp, &
      8.531731613887943e-01_dp, &
      1.455015322134050e-01_dp, &
      1.042780525932260e-23_dp, &
      3.324518992641168e-32_dp, &
      2.215657306846010e-24_dp, &
      1.272590434869170e-05_dp, &
      1.407156949719381e-39_dp, &
      1.501096396346585e-14_dp, &
      1.279732322882177e-15_dp, &
      8.000386979842814e-27_dp, &
      5.942520535535404e-07_dp, &
      2.536655340883834e-06_dp, &
      2.837053768754091e-06_dp, &
      6.063950848321747e-15_dp, &
      7.669076761338264e-08_dp, &
      3.779282843519088e-31_dp, &
      1.230756622872010e-05_dp, &
      4.517417332085708e-15_dp, &
      1.018993933154473e-10_dp, &
      5.606388106792407e-05_dp, &
      4.395887817992532e-27_dp, &
      2.665884349002063e-27_dp, &
      6.406957289164702e-25_dp, &
      2.916480301529460e-28_dp, &
      6.482812282847260e-04_dp, &
      1.933633849489186e-38_dp, &
      2.835207940323098e-28_dp, &
      1.228293692087729e-22_dp, &
      4.463306559593289e-16_dp, &
      0.000000000000000e+00_dp, &
      1.729011724690401e-20_dp, &
      2.376605339725764e-37_dp, &
      3.817089469303068e-20_dp, &
      1.580758530397593e-19_dp, &
      8.630873323290326e-42_dp, &
      2.933670533266407e-17_dp, &
      3.477277095869979e-14_dp, &
      1.672217504051772e-27_dp, &
      2.901964760389726e-10_dp, &
      4.009508311969189e-14_dp, &
      4.477949124882182e-04_dp, &
      5.802425487208889e-25_dp, &
      1.507070706230124e-18_dp, &
      3.117988834775885e-33_dp, &
      3.145713280112587e-06_dp, &
      3.640291080449329e-11_dp, &
      1.212369761144567e-08_dp, &
      4.330783027573506e-07_dp, &
      1.340999161150335e-08_dp, &
      1.637887633641227e-23_dp, &
      0.000000000000000e+00_dp, &
      3.272522560297460e-33_dp, &
      1.822570028256660e-21_dp, &
      1.163861617512337e-14_dp, &
      3.447922902061464e-06_dp, &
      1.750125291174481e-20_dp, &
      1.952232981780049e-11_dp, &
      2.253541113131996e-05_dp, &
      3.869277213853278e-21_dp, &
      3.627438028270097e-18_dp, &
      4.341046527015537e-07_dp, &
      1.064847984094941e-07_dp, &
      1.038902490444878e-14_dp, &
      0.000000000000000e+00_dp, &
      3.271205495198492e-21_dp, &
      3.242755947534449e-17_dp, &
      1.035660061617894e-14_dp, &
      2.700056103682119e-09_dp, &
      1.809450315782121e-09_dp, &
      6.017508591220193e-12_dp, &
      1.809876338439503e-12_dp, &
      4.411013966372366e-34_dp, &
      2.509637985492699e-12_dp, &
      3.301416573163102e-24_dp, &
      1.455015321938003e-08_dp, &
      0.000000000000000e+00_dp, &
      5.531811528594553e-05_dp, &
      0.000000000000000e+00_dp, &
      5.405894509609382e-05_dp, &
      2.408969779610401e-06_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      1.523588081832998e-07_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp &
    ]

    cea = ChemEquiAnalysis(fixture_file('thermo_easy_chem_simp_own.yaml'), atoms, species, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    cea2 = ChemEquiAnalysis(fixture_file('thermo_easy_chem_simp_own.inp'), atoms, species, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    converged = cea%solve(1.0e6_dp, 1000.0_dp, molfracs_atoms=X, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    if (.not.converged) stop 1

    converged = cea2%solve(1.0e6_dp, 1000.0_dp, molfracs_atoms=X, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    if (.not.converged) stop 1

    do i = 1,size(cea%molfracs_species)
      if (.not.is_close(cea%molfracs_species(i),correct_answer(i),tol=1.0e-4_dp) .and. cea%molfracs_species(i) > 1.0e-50_dp) then
        print*,cea%molfracs_species(i),correct_answer(i)
        print*,'ChemEquiAnalysis failed to compute the right equilibrium.'
        stop 1
      endif
    enddo

    do i = 1,size(cea2%molfracs_species)
      if (.not.is_close(cea2%molfracs_species(i),correct_answer(i),tol=1.0e-4_dp) .and. cea2%molfracs_species(i) > 1.0e-50_dp) then
        print*,cea2%molfracs_species(i),correct_answer(i)
        print*,'ChemEquiAnalysis failed to compute the right equilibrium.'
        stop 1
      endif
    enddo

    converged = cea%solve(1.0e6_dp, 1000.0_dp, molfracs_species=cea%molfracs_species, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    if (.not.converged) stop 1

    do i = 1,size(cea%molfracs_species)
      if (.not.is_close(cea%molfracs_species(i),correct_answer(i),tol=1.0e-4_dp) .and. cea%molfracs_species(i) > 1.0e-50_dp) then
        print*,cea%molfracs_species(i),correct_answer(i)
        print*,'ChemEquiAnalysis failed to compute the right equilibrium.'
        stop 1
      endif
    enddo

    converged = cea%solve_metallicity(1.0e6_dp, 1000.0_dp, 1.0_dp, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    if (.not.converged) stop 1

    do i = 1,size(cea%molfracs_species)
      if (.not.is_close(cea%molfracs_species(i),correct_answer(i),tol=1.0e-4_dp) .and. cea%molfracs_species(i) > 1.0e-50_dp) then
        print*,cea%molfracs_species(i),correct_answer(i)
        print*,'ChemEquiAnalysis failed to compute the right equilibrium.'
        stop 1
      endif
    enddo

    print*,'test passed.'

  end subroutine

  subroutine test_sonora()
    type(ChemEquiAnalysis) :: cea
    character(:), allocatable :: err
    character(s_str_len), allocatable :: species(:)
    real(dp), allocatable :: molfracs_atoms_sun(:)
    logical :: converged
    integer, parameter :: nT = 10, nP = 10
    integer :: i, iP, iT
    real(dp) :: P, T

    species = [ &
      'e-                  ', 'H2                  ', 'H                   ', 'H+                  ', &
      'H-                  ', 'H2-                 ', 'H2+                 ', 'H3+                 ', &
      'He                  ', 'H2O                 ', 'CH4                 ', 'CO                  ', &
      'NH3                 ', 'N2                  ', 'PH3                 ', 'H2S                 ', &
      'TiO                 ', 'VO                  ', 'Fe                  ', 'FeH                 ', &
      'CrH                 ', 'Na                  ', 'K                   ', 'Rb                  ', &
      'atCs                ', 'CO2                 ', 'HCN                 ', 'C2H2                ', &
      'C2H4                ', 'C2H6                ', 'COS                 ', 'SiO                 ', &
      'MgH                 ', 'Li                  ', 'LiOH                ', 'LiH                 ', &
      'LiCl                ', 'Li+                 ', 'LiF                 ', 'OH                  ', &
      'C-gr                ', 'Mg                  ', 'Mg+                 ', 'Si                  ', &
      'Fe+                 ', 'Ti                  ', 'Ti+                 ', 'C                   ', &
      'O                   ', 'C+                  ', 'O+                  ', 'He+                 ', &
      'C2                  ', 'CH                  ', 'CN                  ', 'CS                  ', &
      'C2H                 ', 'CH2                 ', 'CH3                 ', 'C3H8                ', &
      'HCHO                ', 'CH2OH               ', 'CH3OH               ', 'CH3O                ', &
      'N                   ', 'NH                  ', 'NH2                 ', 'NO                  ', &
      'N2H2                ', 'N2H4                ', 'O2                  ', 'H2O2                ', &
      'P                   ', 'PH2                 ', 'P2                  ', 'PO                  ', &
      'PH                  ', 'P4O6(Gurvich)       ', 'S                   ', 'SH                  ', &
      'SN                  ', 'SO                  ', 'S2                  ', 'SO2                 ', &
      'S-                  ', 'SH-                 ', 'Cr                  ', 'Cr+                 ', &
      'CrO                 ', 'CrO2                ', 'FeO                 ', 'FeOH                ', &
      'FeS                 ', 'Fe(OH)2             ', 'FeCl                ', 'MgO                 ', &
      'MgOH                ', 'MgS                 ', 'Mg(OH)2             ', 'Si+                 ', &
      'SiS                 ', 'SiH                 ', 'SiO2                ', 'SiH2                ', &
      'SiH3                ', 'SiH4                ', 'Na+                 ', 'NaCl                ', &
      'NaOH                ', 'NaH                 ', 'K+                  ', 'KCl                 ', &
      'KH                  ', 'KOH                 ', 'V                   ', 'V+                  ', &
      'VO2                 ', 'TiO2                ', 'Cl-                 ', 'Cl                  ', &
      'HCl                 ', 'RbCl                ', 'Rb+                 ', 'RbH                 ', &
      'RbO                 ', 'RbOH                ', 'RbF                 ', 'CsCl                ', &
      'Cs+                 ', 'CsH                 ', 'F                   ', 'F-                  ', &
      'HF                  ', 'NaF                 ', 'NH4H2PO4(c)         ', 'VO(c)               ', &
      'VO(L)               ', 'TiO2(c)             ', 'TiO2(L)             ', 'MgO(c)              ', &
      'MgO(L)              ', 'SiO2(c)             ', 'SiO2(L)             ', 'Cr(c)               ', &
      'Cr(L)               ', 'Fe(c)               ', 'Fe(L)               ', 'H2O(L)              ', &
      'H2O(c)              ', 'Na2S(c)             ', 'KCl(c)              ', 'RbCl(c)             ', &
      'CsCl(c)             ', 'Li2S(c)             ', 'LiF(cr)             ' &
    ]

    cea = ChemEquiAnalysis(fixture_file('thermo_sonora_component.yaml'), species=species, err=err)
    if (allocated(err)) then
      print*, err
      stop 1
    endif

    cea%mass_tol = 1.0e-2_dp

    allocate(molfracs_atoms_sun(size(cea%atoms_names)))
    do i = 1,size(cea%atoms_names)
      select case (trim(cea%atoms_names(i)))
      case ('H')
        molfracs_atoms_sun(i) = 9.082387e-01_dp
      case ('He')
        molfracs_atoms_sun(i) = 9.046346e-02_dp
      case ('Li')
        molfracs_atoms_sun(i) = 2.050745e-09_dp
      case ('C')
        molfracs_atoms_sun(i) = 3.286959e-04_dp
      case ('N')
        molfracs_atoms_sun(i) = 7.893027e-05_dp
      case ('O')
        molfracs_atoms_sun(i) = 5.982842e-04_dp
      case ('F')
        molfracs_atoms_sun(i) = 4.577235e-08_dp
      case ('Na')
        molfracs_atoms_sun(i) = 2.083182e-06_dp
      case ('Mg')
        molfracs_atoms_sun(i) = 3.712245e-05_dp
      case ('Si')
        molfracs_atoms_sun(i) = 3.604122e-05_dp
      case ('P')
        molfracs_atoms_sun(i) = 2.977005e-07_dp
      case ('S')
        molfracs_atoms_sun(i) = 1.575001e-05_dp
      case ('Cl')
        molfracs_atoms_sun(i) = 1.906580e-07_dp
      case ('K')
        molfracs_atoms_sun(i) = 1.301448e-07_dp
      case ('Ti')
        molfracs_atoms_sun(i) = 8.862536e-08_dp
      case ('V')
        molfracs_atoms_sun(i) = 9.911335e-09_dp
      case ('Cr')
        molfracs_atoms_sun(i) = 4.732212e-07_dp
      case ('Fe')
        molfracs_atoms_sun(i) = 3.142794e-05_dp
      case ('Rb')
        molfracs_atoms_sun(i) = 2.584155e-10_dp
      case ('Cs')
        molfracs_atoms_sun(i) = 1.326317e-11_dp
      case default
        print*, 'Unexpected Sonora atom in test_sonora: ', trim(cea%atoms_names(i))
        stop 1
      end select
    enddo
    molfracs_atoms_sun = molfracs_atoms_sun/sum(molfracs_atoms_sun)
    cea%molfracs_atoms_sun = molfracs_atoms_sun

    do iT = 1,nT
      if (nT == 1) then
        T = 5.0e2_dp
      else
        T = 5.0e2_dp + (1.0e3_dp - 5.0e2_dp)*real(iT - 1, dp)/real(nT - 1, dp)
      endif
      do iP = 1,nP
        if (nP == 1) then
          P = 1.0e9_dp
        else
          P = 1.0e1_dp**(9.0_dp - 9.0_dp*real(iP - 1, dp)/real(nP - 1, dp))
        endif

        converged = cea%solve_metallicity(P, T, 10.0_dp**3.5_dp, CtoO=0.5_dp, err=err)
        if (allocated(err)) then
          print*, 'test_sonora failed at P, T = ', P, T
          print*, err
          stop 1
        endif
        if (.not.converged) then
          print*, 'test_sonora did not converge at P, T = ', P, T
          stop 1
        endif

        if (size(cea%molfracs_species) /= size(species)) then
          print*, 'Sonora test returned the wrong number of species.'
          print*, 'P, T = ', P, T
          stop 1
        endif
        if (size(cea%molfracs_species_gas) /= size(cea%gas_names)) then
          print*, 'Sonora test returned the wrong number of gas species.'
          print*, 'P, T = ', P, T
          stop 1
        endif
        if (size(cea%molfracs_species_condensate) /= size(cea%condensate_names)) then
          print*, 'Sonora test returned the wrong number of condensates.'
          print*, 'P, T = ', P, T
          stop 1
        endif
        if (any(.not.(cea%molfracs_species >= 0.0_dp))) then
          print*, 'Sonora test produced negative species mole fractions.'
          print*, 'P, T = ', P, T
          stop 1
        endif
      enddo
    enddo

    print*, 'test_sonora passed.'
  end subroutine

  subroutine test_memory()
    type(ChemEquiAnalysis), pointer :: cea
    character(:), allocatable :: err

    allocate(cea)

    cea = ChemEquiAnalysis(fixture_file('thermo_easy_chem_simp_own.yaml'), err=err)
    if (allocated(err)) then
      deallocate(cea)
      print*,err
      stop 1
    endif

    cea = ChemEquiAnalysis(fixture_file('thermo_easy_chem_simp_own.yaml'), err=err)
    if (allocated(err)) then
      deallocate(cea)
      print*,err
      stop 1
    endif

    deallocate(cea)

    print*,'test_memory passed.'

  end subroutine

  subroutine test_nasa7()
    use equilibrate_cea, only: entropy_nasa7, enthalpy_nasa7, heat_capacity_nasa7
    type(ChemEquiAnalysis) :: cea
    character(:), allocatable :: err
    logical :: converged
    real(dp) :: entropy, enthalpy, heat_capacity
    real(dp), allocatable :: coeffs(:)

    cea = ChemEquiAnalysis(fixture_file('photochem_thermo.yaml'), err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    converged = cea%solve_metallicity(1.0e6_dp, 1000.0_dp, 1.0_dp, err=err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif
    if (.not.converged) stop 1

    coeffs = [-18.42562_dp, 0.1279622_dp, -0.000241671_dp, 2.073899e-07_dp, -6.684715e-11_dp, -107350.5_dp, 81.90459_dp]
    entropy = entropy_nasa7(coeffs, 900.0_dp)
    if (.not. is_close(entropy, 110.471033963557_dp)) stop 1
    enthalpy = enthalpy_nasa7(coeffs, 900.0_dp)
    if (.not. is_close(enthalpy, -870626.2619052551_dp)) stop 1
    heat_capacity = heat_capacity_nasa7(coeffs, 900.0_dp)
    if (.not. is_close(heat_capacity, 69.14032042926877_dp)) stop 1

    print*,'test_nasa7 passed.'

  end subroutine

  !> coppied from fortran stdlib v0.2.0
  elemental function is_close(a, b, tol, abs_tol, equal_nan) result(close)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: tol, abs_tol
    logical, intent(in), optional :: equal_nan
    logical :: close

    real(dp) :: rel_tol_, abs_tol_
    logical :: equal_nan_

    if (present(tol)) then
      rel_tol_ = tol
    else
      rel_tol_ = 1.0e-5_dp
    endif

    if (present(abs_tol)) then
      abs_tol_ = abs_tol
    else
      abs_tol_ = 0.0_dp
    endif

    if (present(equal_nan)) then
      equal_nan_ = equal_nan
    else
      equal_nan_ = .false.
    endif

    if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
        close = merge(.true., .false., equal_nan_ .and. ieee_is_nan(a) .and. ieee_is_nan(b))
    else
        close = abs(a - b) <= max(abs(rel_tol_*max(abs(a), abs(b))), abs(abs_tol_))
    end if     

  end function


end program
