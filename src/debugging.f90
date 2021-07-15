
program main
  
  use photochem_setup, only: setup
  use photochem
  
  use photochem_const, only: pi
  use photochem_data, only: nq, nsp, nsl, nrT, kj, nw, LH2O, water_sat_trop, photonums, &
                            reactants_names, species_names
  use photochem_vars, only: nz, z, dz, trop_ind, &
                            usol_init, data_dir, neqs, z, use_fast_jacobian

  implicit none
  character(len=1024) :: err
  integer, parameter :: real_kind = kind(1.0d0)
  real(real_kind), allocatable :: mubar(:), pressure(:)
  real(real_kind), allocatable :: density(:), fH2O(:)
  real(real_kind), allocatable :: densities(:,:)
  real(real_kind), allocatable :: rx_rates(:,:)
  real(real_kind), allocatable :: prates(:,:), surf_radiance(:)
  real(real_kind), allocatable :: xp(:), xl(:)
  real(real_kind), allocatable :: DU(:,:), DD(:,:), DL(:,:), ADU(:,:), ADL(:,:)
  real(real_kind), allocatable :: fH2O_save(:), lower_fix_mr_save(:)
  
  real(real_kind), pointer :: usol_flat(:)
  real(real_kind), allocatable :: rhs(:)
  
  real(real_kind) :: t_eval(1)
  real(real_kind), allocatable, target :: solution(:,:,:)
  logical :: success
  
  
  character(len=:), allocatable :: rxstring
  character(len=50) :: rxstring1
  
  real(real_kind) :: disth, ztop, ztop1    
  integer :: i, k, j, jdisth

  data_dir = "../data"

  call setup("../zahnle_earth.yaml", "../settings_Hadean.yaml", "../Sun_4.0Ga.txt", "../atmosphere_Hadean.txt", err)
  if (len(trim(err)) > 0) then
    print*,trim(err)
    stop
  endif
  
  allocate(mubar(nz), pressure(nz))
  allocate(density(nz), fH2O(nz))
  allocate(densities(nsp+1,nz))
  allocate(rx_rates(nz,nrT))
  allocate(prates(nz,kj), surf_radiance(nw))
  allocate(xp(nz), xl(nz))
  allocate(DU(nq,nz), DD(nq,nz), DL(nq,nz), ADU(nq,nz), ADL(nq,nz))
  allocate(fH2O_save(trop_ind), lower_fix_mr_save(nz))
  
  allocate(solution(nq,nz,1))
  allocate(rhs(neqs))
  ! use_fast_jacobian = .false.
  t_eval = 1.d10
  ! max_order = 2
  call evolve_background_atm(0.d0, nq, nz, usol_init, 1, t_eval, 1.d-3, 1.d-30, &
                             100000, solution, success, err)

  usol_init = solution(:,:,1)
  call prep_all_background_gas(nsp, nq, nz, nrT, kj, nw, trop_ind, usol_init, densities, &
                               density, rx_rates, mubar, pressure, fH2O, fH2O_save, &
                               prates, surf_radiance, &
                               DU, DD, DL, ADU, ADL, err)
  ! do i = 1,nq                         
    ! print*,species_names(i),DU(i,22),DD(i,22),DL(i,22),ADU(i,22),ADL(i,22)
  ! enddo
                               
  usol_flat(1:neqs) => usol_init            
  call rhs_background_gas(neqs, usol_flat, rhs, err)
  ! usol_init = solution(:,:,1)  
  ! do i = 1,kj
  !   k = photonums(i)
  !   call reaction_string(k,rxstring)
  !   print*,rxstring,prates(nz,i)
  ! enddo
  ! do i = 1,nz
  !   k = photonums(1)
  !   call reaction_string(k,rxstring)
  !   print*,z(i)/1.d5,i,rxstring,prates(i,1),usol_init(3,i)
  ! enddo
  
  
  ! do i = 1,nrT
  !   call reaction_string(i,rxstring)
  !   rxstring1 = ''
  !   rxstring1 = rxstring
  !   print*,i,rxstring1,rx_rates(1, i)
  ! enddo
  
  ! do i = 1,nq
  !   k = i + (1-1) * nq
  ! 
  !   print*,species_names(i),usol_flat(k),rhs(k)
  ! 
  ! enddo
  ! stop
  
  i = 1
  do j = 1,nz
    k = i + (j-1) * nq
    
    print*,usol_flat(k)
     
  enddo
  stop
  
  
  
  do j = 1,nz
    do k = 1,nq
      densities(k,j) = usol_init(k,j)*density(j)
    enddo
    densities(nsp,j) = (1.d0-sum(usol_init(:,j)))*density(j) ! background gas
    densities(nsp+1,j) = 1.d0 ! for hv
  enddo
  
  ! short lived
  do k = nq+1,nq+nsl
    call chempl(nz, nsp, nrT, densities, rx_rates, k, 2, xp, xl) 
    densities(k,:) = xp/xl
  enddo
  
  ! long lived              
  ! do i = 1,nq


  

  

end program
