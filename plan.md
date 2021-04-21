
Problem with PhotochemPy
- So much global data. Need to split globals into multiple modules, and pass more data to subroutines
- xsections suck
- mysterious if statements
- haze particles


Modules
- Photochem



I need to deal with more than 2 reactants and 3 products.

I think I understand the production rates. But i'm still worried about double counting reactions for stuff like H + H => H2. haha no double counting makes sense!

```fortran
np = nump(k) ! k is a species
! np is number of reactions that produce species k
do i=1,np
  m = iprod(i,k) ! m is reaction number
  l = numreactants(m) ! l is the number of reactants
  do ii = 1,l
    sp_inds(ii) = jchem(ii,k)
  enddo
  do iii = 1,nz
    DD = 1.d0
    do iiii = 1,l
      DD = DD * D(sp_inds(iiii),iii)
    enddo
    xp(iii) = xp(iii) + A(iii,k) * DD
  enddo
enddo
```

What about loss rates? We can mirror above. This will make `xl` in units of molecules/cm3/s. This is different than PhotochemPy

```fortran
nl = numl(k) ! k is a species
! nl is number of reactions that destroy species k
do i=1,np
  m = iloss(i,k) ! This will JUST be reaction number
  l = numreactants(m) ! number of reactants
  do ii = 1,l
    sp_inds(ii) = jchem(ii,k)
  enddo
  do iii = 1,nz
    DD = 1.d0
    do iiii = 1,l
      DD = DD * D(sp_inds(iiii),iii)
    enddo
    xl(iii) = xl(iii) + A(iii,k) * DD
  enddo
enddo
```

