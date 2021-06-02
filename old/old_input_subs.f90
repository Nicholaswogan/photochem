subroutine compare_rxtype_string(tmp, eqr, eqp, reverse, rxtype_int, err)
  type(string), allocatable, intent(in) :: eqr(:), eqp(:)
  type(string), intent(in) :: tmp    
  logical, intent(in) :: reverse
  integer, intent(in) :: rxtype_int
  character(len=err_len), intent(out) :: err
  character(len=15) :: rxtype
  integer i
  logical k, j, m, l, kk, jj
  l = .false.
  m = .false.
  k = .false.
  kk = .false.
  j = .false.
  jj = .false.
  err = ''
  if (rxtype_int == 0) then
    rxtype = 'photolysis'
  elseif (rxtype_int == 1) then
    rxtype = 'elementary'
  elseif (rxtype_int == 2) then
    rxtype = 'three-body'
  elseif (rxtype_int == 3) then
    rxtype = 'falloff'
  endif

  if ((trim(rxtype) == 'three-body') .or. (trim(rxtype) == 'falloff')) then
    do i = 1,size(eqr)
      if ((trim(eqr(i)%chars()) == 'M').and.j) jj = .true.
      if (trim(eqr(i)%chars()) == 'M') j = .true.
      if (trim(eqr(i)%chars()) == 'hv') m = .true.
    enddo
    do i = 1,size(eqp)
      if ((trim(eqp(i)%chars()) == 'M').and.k) kk = .true.
      if (trim(eqp(i)%chars()) == 'M') k = .true.
      if (trim(eqp(i)%chars()) == 'hv') l = .true.
    enddo
    if (trim(eqr(size(eqr))%chars()) /= 'M') then
      err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
              ' must have "M" as the last reactant'
      return
    endif
    if (trim(eqp(size(eqp))%chars()) /= 'M') then
      err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
              ' must have "M" as the last product'
      return
    endif
    if ((j).and.(k)) then
      ! good
    else
      err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
              ' must have "M" on both sides'
      return
    endif
    if ((jj).or.(kk)) then
      err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
              ' can only have one "M" on either side'
      return
    endif
    if ((m).or.(l)) then
      err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
              ' can not contain "hv". Only photolysis reactions can.'
      return
    endif
  elseif (trim(rxtype) == 'elementary') then
    do i = 1,size(eqr)
      if (trim(eqr(i)%chars()) == 'M') j = .true.
      if (trim(eqr(i)%chars()) == 'hv') m = .true.
    enddo
    do i = 1,size(eqp)
      if (trim(eqp(i)%chars()) == 'M') k = .true.
      if (trim(eqp(i)%chars()) == 'hv') l = .true.
    enddo
    if ((j).or.(k)) then
      err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
              ' can not contain "M".'
      return
    endif
    if ((m).or.(l)) then
      err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
              ' can not contain "hv". Only photolysis reactions can.'
      return
    endif
  elseif (trim(rxtype) == 'photolysis') then
    if (reverse) then
      err = 'IOError: Photolysis reaction '//tmp//' can not be reversed.'
    endif
    
    do i = 1,size(eqr)
      if (trim(eqr(i)%chars()) == 'M') j = .true.
      if ((trim(eqr(i)%chars()) == 'hv').and.(m)) jj = .true.
      if (trim(eqr(i)%chars()) == 'hv') m = .true.
    enddo
    do i = 1,size(eqp)
      if (trim(eqp(i)%chars()) == 'M') k = .true.
      if (trim(eqp(i)%chars()) == 'hv') l = .true.
    enddo
    if ((j).or.(k)) then
      err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
              ' can not contain "M".'
      return
    endif
    if (jj) then
      err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
              ' can only have one "hv" on the left side.'
      return
    endif
    if ((m).and..not.(l)) then
      ! good
    else
      err = 'IOError: '//trim(rxtype)// ' reaction '//tmp// &
              ' must have "hv" on the left and no "hv" on the right.'
      return
    endif
  endif
end subroutine



subroutine parse_reaction(instring, reverse, eqr, eqp, err)
  type(string), intent(in) :: instring
  logical, intent(out) :: reverse
  type(string), allocatable, intent(out) :: eqr(:), eqp(:)
  character(len=err_len), intent(out) :: err
  
  type(string) :: string1, string2, string3
  type(string), allocatable :: eq1(:), eq2(:)
  type(string), allocatable :: eqr1(:), eqp1(:)
  integer i
  string1 = instring%replace(old='(', new=' ')
  string2 = string1%replace(old=')', new=' ')
  string3 = string2%replace(old='+', new=' ')
  if (index(instring%chars(), "<=>") /= 0) then
    call string2%split(eq1, sep="<=>")
    call string3%split(eq2, sep="<=>")
    reverse = .true.
  elseif (index(instring%chars(), " =>") /= 0) then
    call string2%split(eq1, sep="=>")
    call string3%split(eq2, sep="=>")
    reverse = .false.
  else
    err = "IOError: Invalid reaction arrow in reaction "//instring// &
          '. Note, forward reactions must have a space before the arrow, like " =>"'
    return
  endif
  
  call eq1(1)%split(eqr1, sep="+")
  call eq1(2)%split(eqp1, sep="+")
  ! remove white space
  allocate(eqr(size(eqr1)))
  allocate(eqp(size(eqp1)))
  do i=1,size(eqr1)
    eqr(i) = eqr1(i)%replace(old=' ', new='')
  enddo
  do i=1,size(eqp1)
    eqp(i) = eqp1(i)%replace(old=' ', new='')
  enddo
  
  call eq2(1)%split(eqr1, sep=" ")
  call eq2(2)%split(eqp1, sep=" ")
  
  if ((size(eqr1) /= size(eqr)) .or. (size(eqp1) /= size(eqp))) then
    err = 'IOError: Missing "+" sign(s) in reaction '//instring
    return
  endif
  
  if (size(eqr) + size(eqp) /= instring%count('+') + 2) then
    err = 'IOError: Too many "+" signs in reaction '//instring
    return
  endif
end subroutine


subroutine get_reaction_chars(instring, max_num_react, max_num_prod, numr, nump, &
                              outreact, outprod, reverse, err)
  type(string), intent(in) :: instring
  integer, intent(in) :: max_num_react, max_num_prod
  
  integer, intent(out) :: numr, nump
  character(len=8), intent(out) :: outreact(max_num_react), outprod(max_num_prod)
  logical, intent(out) :: reverse
  character(len=err_len), intent(out) :: err
  
  type(string), allocatable :: eqr(:), eqp(:)
  integer :: i
  
  call parse_reaction(instring, reverse, eqr, eqp, err)
  if (len_trim(err) > 0) return
  
  if (eqr(size(eqr)) == 'M') then
    eqr = eqr(1:size(eqr)-1)
    eqp = eqp(1:size(eqp)-1)
  endif
  
  numr = size(eqr)
  nump = size(eqp)
  
  outreact = ''
  outprod = ''
  do i=1,numr
    outreact(i) = eqr(i)%chars()
  enddo
  do i=1,nump
    outprod(i) = eqp(i)%chars()
  enddo
end subroutine


