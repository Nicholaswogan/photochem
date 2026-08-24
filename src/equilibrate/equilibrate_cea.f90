
module equilibrate_cea
   use equilibrate_const, only: dp, atom_str_len, reac_str_len
   implicit none
   private

   public :: CEAData
   public :: entropy_nasa7, enthalpy_nasa7, heat_capacity_nasa7

   type :: CEAData

      real(dp) :: mass_tol = 1.0e-2_dp !! Degree to which mass will be balanced. 
                                       !! Gordon & McBride's default is 1.0e-6, but
                                       !! this seems too strict.
      logical :: verbose2 !! actual verbosity variable
      logical :: converged !! for passing convergence

      ! Run variables
      integer     :: iter_max
      integer     :: N_atoms = 0, N_reactants = 0, N_gas = 0, N_cond = 0, N_ions = 0
      logical     :: verbose, verbose_cond, quick, ions, remove_ions, error
      character(len=500)   :: err_msg

      ! List of atoms
      character(len=atom_str_len), allocatable   :: names_atoms(:)  !(N_atoms)
      integer, allocatable            :: id_atoms(:)  !(N_atoms)

      ! List of reactants reordered with condensates at the end
      character(len=reac_str_len), allocatable  :: names_reactants(:), names_reactants_orig(:)  !(N_reac)
      integer, allocatable            :: id_reactants(:,:)  !(N_reac,2)

      ! Atomic data for each reactant
      character(len=atom_str_len), allocatable   :: reac_atoms_names(:,:)  !(5,N_reac)
      integer, allocatable            :: reac_atoms_id(:,:)  !(5,N_reac)
      real(dp), allocatable   :: reac_stoich(:,:)  !(5,N_reac)

      ! Nature of reactant
      logical, allocatable    :: reac_condensed(:), reac_ion(:)  !(N_reac)

      ! Thermodynamic data arrays
      ! integer :: N_coeffs = 10, N_temps = 10
      integer, allocatable    :: thermo_data_n_coeffs(:,:)  !(N_temps,N_reac)
      integer, allocatable    :: thermo_data_n_intervs(:)  !(N_reac)
      integer, allocatable    :: thermo_data_kind(:) !(N_reac)
      real(dp), allocatable   :: thermo_data(:,:,:)  !(N_coeffs,N_temps,N_reac)
      real(dp), allocatable   :: thermo_data_temps(:,:,:)  !(2,N_temps,N_reac)
      real(dp), allocatable   :: thermo_data_T_exps(:,:,:)  !(8,N_temps,N_reac)
      real(dp), allocatable   :: form_heat_Jmol_298_15_K(:)  !(N_reac)
      real(dp), allocatable   :: H_0_298_15_K_m_H_0_0_K(:,:)  !(N_temps, N_reac)
      real(dp), allocatable   :: mol_weight(:)  !(N_reac)
      real(dp), allocatable   :: a_stoich(:,:)  !(N_reac,N_atoms)
      real(dp), allocatable   :: b_0_cache(:)  !(N_atoms)
      real(dp)                :: b_0_norm_cache = 0d0
      logical                 :: b_0_cache_valid = .false.
      integer                 :: b_0_cache_atoms_use = 0
      logical                 :: prev_guess_valid = .false.
      real(dp), allocatable   :: prev_n_spec(:)  !(N_reac)
      real(dp)                :: prev_n = 0d0
   
   contains
      procedure :: set_data
      procedure :: solve
   end type

contains

   !> INITIALIZE ALL DATA
   subroutine SET_DATA(self, fpath, atoms_names, reactants_names)
      use equilibrate_const, only: N_coeffs, N_temps
      use equilibrate_yaml, only: check_for_duplicates, ReactantList
      class(CEAData), intent(inout) :: self
      character(*), intent(in) :: fpath
      character(atom_str_len), optional, intent(in) :: atoms_names(:)
      character(reac_str_len), optional, intent(in) :: reactants_names(:)

      type(ReactantList) :: rl
      character(:), allocatable :: err
      integer :: N_atoms_in, N_reactants_in
      character(atom_str_len), allocatable :: atoms_names_(:)
      character(reac_str_len), allocatable :: reactants_names_(:)
      character(atom_str_len) :: atom_tmp 
      character(reac_str_len) :: reactant_tmp
      integer :: ind, i, j
      logical :: keep, found

      self%error = .false.

      if (file_extension(fpath) == 'yaml') then
         ! We parse the yaml file
         rl = ReactantList(fpath, err)
         if (allocated(err)) then
            self%error = .true.
            self%err_msg = err
            return
         endif
      endif

      if (.not.present(atoms_names) .or. .not.present(reactants_names)) then
         if (.not.(file_extension(fpath) == 'yaml')) then
            self%error = .true.
            self%err_msg = 'The thermo file must be a yaml file, if you do not include '// &
                           'a list of reactants or atoms.'
            return
         endif
      endif

      ! 4 cases to consider
      if (present(atoms_names) .and. present(reactants_names)) then
         ! Copy over, but exclude "E" if it present.
         allocate(atoms_names_(1))
         do i = 1,size(atoms_names)
            if (trim(atoms_names(i)) /= 'E') then
               atom_tmp = atoms_names(i)
               atoms_names_ = [atoms_names_, atom_tmp]
            endif
         enddo
         atoms_names_ = atoms_names_(2:size(atoms_names_))
         reactants_names_ = reactants_names 
      elseif (present(atoms_names) .and. .not.present(reactants_names)) then
         ! Copy over, but exclude "E" if it present.
         allocate(atoms_names_(1))
         do i = 1,size(atoms_names)
            if (trim(atoms_names(i)) /= 'E') then
               atom_tmp = atoms_names(i)
               atoms_names_ = [atoms_names_, atom_tmp]
            endif
         enddo
         atoms_names_ = atoms_names_(2:size(atoms_names_))

         ! Make sure that input atoms are in yaml file
         do i = 1,size(atoms_names)
            ind = findloc(rl%atoms_names, atoms_names(i), 1)
            if (ind == 0) then
               self%error = .true.
               self%err_msg = 'Input atom '//trim(atoms_names(i))//' is not in '//trim(fpath)
               return
            endif
         enddo

         ! Determine reactants from Yaml file
         allocate(reactants_names_(1))
         do i = 1,rl%nr
            keep = .true.
            do j = 1,rl%natoms
               if (rl%r(i)%composition(j) /= 0) then
                  ind = findloc(atoms_names, rl%atoms_names(j), 1)
                  if (ind == 0) then
                     keep = .false.
                     exit
                  endif
               endif
            enddo
            if (keep) then
               reactant_tmp = rl%r(i)%name
               reactants_names_ = [reactants_names_, reactant_tmp]
            endif
         enddo
         reactants_names_ = reactants_names_(2:size(reactants_names_))

      elseif (.not.present(atoms_names) .and. present(reactants_names)) then; block
         integer, allocatable :: composition(:)
         reactants_names_ = reactants_names

         allocate(composition(rl%natoms))
         composition = 0
         do i = 1,size(reactants_names_)
            found = .false.
            do j = 1,rl%nr
               if (trim(reactants_names_(i)) == trim(rl%r(j)%name)) then
                  found = .true.
                  composition = composition + abs(rl%r(j)%composition)
                  exit
               endif
            enddo
            if (.not.found) then
               self%error = .true.
               self%err_msg = 'Input reactant '//trim(reactants_names_(i))//' is not in '//trim(fpath)
               return
            endif
         enddo
         allocate(atoms_names_(1))
         do i = 1,size(composition)
            if (composition(i) /= 0 .and. trim(rl%atoms_names(i)) /= 'E') then
               atom_tmp = rl%atoms_names(i)
               atoms_names_ = [atoms_names_, atom_tmp]
            endif
         enddo
         atoms_names_ = atoms_names_(2:size(atoms_names_))

      endblock; elseif (.not.present(atoms_names) .and. .not.present(reactants_names)) then

         allocate(reactants_names_(1))
         do i = 1,rl%nr
            reactant_tmp = rl%r(i)%name
            reactants_names_ = [reactants_names_, reactant_tmp]
         enddo
         reactants_names_ = reactants_names_(2:size(reactants_names_))
         allocate(atoms_names_(1))
         do i = 1,rl%natoms
            if (trim(rl%atoms_names(i)) /= 'E') then
               atom_tmp = rl%atoms_names(i)
               atoms_names_ = [atoms_names_, atom_tmp]
            endif
         enddo
         atoms_names_ = atoms_names_(2:size(atoms_names_))
         
      endif

      N_atoms_in = size(atoms_names_)
      N_reactants_in = size(reactants_names_)

      ! REACTANTS
      if (N_reactants_in /= self%N_reactants .and. allocated(self%names_reactants_orig)) then
         ! Deallocate everything if the number of reactants changed
         deallocate(self%names_reactants_orig)
         deallocate(self%names_reactants)
         deallocate(self%id_reactants)
         deallocate(self%reac_atoms_names)
         deallocate(self%reac_atoms_id)
         deallocate(self%reac_stoich)
         deallocate(self%reac_condensed)
         deallocate(self%reac_ion)
         deallocate(self%thermo_data_n_coeffs)
         deallocate(self%thermo_data_n_intervs)
         deallocate(self%thermo_data_kind)
         deallocate(self%thermo_data)
         deallocate(self%thermo_data_temps)
         deallocate(self%thermo_data_T_exps)
         deallocate(self%form_heat_Jmol_298_15_K)
         deallocate(self%H_0_298_15_K_m_H_0_0_K)
         deallocate(self%mol_weight)
         if (allocated(self%a_stoich)) then
            deallocate(self%a_stoich)
         end if
         if (allocated(self%b_0_cache)) then
            deallocate(self%b_0_cache)
         end if
         if (allocated(self%prev_n_spec)) then
            deallocate(self%prev_n_spec)
         end if
         self%b_0_cache_valid = .false.
         self%b_0_cache_atoms_use = 0
         self%prev_guess_valid = .false.
         self%prev_n = 0d0
      end if
      if (.not. allocated(self%names_reactants_orig)) then
         ! Allocate reactant-related arrays
         self%N_reactants = N_reactants_in
         allocate(self%names_reactants_orig(self%N_reactants))
         allocate(self%names_reactants(self%N_reactants))
         allocate(self%id_reactants(self%N_reactants,2))
         allocate(self%reac_atoms_names(5,self%N_reactants))
         allocate(self%reac_atoms_id(5,self%N_reactants))
         allocate(self%reac_stoich(5,self%N_reactants))
         allocate(self%reac_condensed(self%N_reactants))
         allocate(self%reac_ion(self%N_reactants))

         allocate(self%thermo_data_n_coeffs(N_temps,self%N_reactants))
         allocate(self%thermo_data_n_intervs(self%N_reactants))
         allocate(self%thermo_data_kind(self%N_reactants))
         allocate(self%thermo_data(N_coeffs,N_temps,self%N_reactants))
         allocate(self%thermo_data_temps(2,N_temps,self%N_reactants))
         allocate(self%thermo_data_T_exps(8,N_temps,self%N_reactants))
         allocate(self%form_heat_Jmol_298_15_K(self%N_reactants))
         allocate(self%H_0_298_15_K_m_H_0_0_K(N_temps,self%N_reactants))
         allocate(self%mol_weight(self%N_reactants))
         allocate(self%b_0_cache(self%N_atoms))
         allocate(self%prev_n_spec(self%N_reactants))
      end if

      ! Set self%names_reactants_orig with the given list of reactants
      ! call da_CH2STR(reac_char, self%names_reactants_orig)
      self%names_reactants_orig = reactants_names_
      ! check for duplicates
      ind = check_for_duplicates(self%names_reactants_orig)
      if (ind /= 0) then
         self%error = .true.
         self%err_msg = '"'//trim(self%names_reactants_orig(ind))//'" is a duplicate reactant input.'
         return
      endif

      if (file_extension(fpath) == 'yaml') then

      call read_thermo_yaml(self, rl, fpath)
      if (self%error) return

      else

      ! Set all thermodynamic data
      call da_READ_THERMO(self, fpath)
      if (self%error) RETURN

      endif

      ! ATOMS
      if (N_atoms_in+1 /= self%N_atoms .and. allocated(self%names_atoms)) then
         deallocate(self%names_atoms)
         deallocate(self%id_atoms)
      end if
      if (.not. allocated(self%names_atoms)) then
         self%N_atoms = N_atoms_in + 1
         allocate(self%names_atoms(self%N_atoms))
         allocate(self%id_atoms(self%N_atoms))
      end if
      if (allocated(self%b_0_cache)) then
         if (size(self%b_0_cache) /= self%N_atoms) then
            deallocate(self%b_0_cache)
         end if
      end if
      if (.not. allocated(self%b_0_cache)) then
         allocate(self%b_0_cache(self%N_atoms))
      end if
      self%b_0_cache_valid = .false.
      self%b_0_cache_atoms_use = 0
      if (allocated(self%prev_n_spec)) then
         deallocate(self%prev_n_spec)
      end if
      allocate(self%prev_n_spec(self%N_reactants))
      self%prev_guess_valid = .false.
      self%prev_n = 0d0

      ! call da_CH2STR(atoms_char, self%names_atoms(1:N_atoms_in))
      self%names_atoms(1:N_atoms_in) = atoms_names_
      self%names_atoms(self%N_atoms) = 'E'
      ! check for duplicates
      ind = check_for_duplicates(self%names_atoms)
      if (ind /= 0) then
         self%error = .true.
         self%err_msg = '"'//trim(self%names_atoms(ind))//'" is a duplicate atom input. '// &
         'Note that "E" should not be an atom input.'
         return
      endif

      call da_ATOMS_ID(self)
      if (self%error) RETURN

      call da_REAC_ATOMS_ID(self)
      if (self%error) RETURN

      call da_BUILD_A_STOICH(self)
      if (self%error) RETURN

      self%prev_guess_valid = .false.

   end subroutine SET_DATA

   subroutine read_thermo_yaml(self, rl, fpath)
      use equilibrate_yaml, only: ReactantList, ShomatePolynomial, Nasa9Polynomial, Nasa7Polynomial
      class(CEAData), intent(inout) :: self
      type(ReactantList), intent(in) :: rl
      character(*), intent(in) :: fpath

      character(len=reac_str_len) :: name1, name2
      logical :: found
      integer :: i, j, k, kk

      ! Not set and not used in code later
      self%form_heat_Jmol_298_15_K = 0.0_dp
      self%thermo_data_n_coeffs = 0
      self%thermo_data_T_exps = 0.0_dp
      self%H_0_298_15_K_m_H_0_0_K = 0.0_dp

      self%N_gas = 0
      ! first loop sets self%N_gas and self%reac_condensed
      do i = 1,self%N_reactants
         call uppercase(self%names_reactants_orig(i), name1)

         found = .false.
         do j = 1,rl%nr
            call uppercase(rl%r(j)%name, name2)
            if (trim(name1) == trim(name2)) then
               found = .true.

               self%reac_condensed(i) = rl%r(j)%condensate
               if (.not.self%reac_condensed(i)) then
                  self%N_gas = self%N_gas + 1
               endif

               exit
            endif
         enddo

         if (.not.found) then
            self%error = .true.
            self%err_msg = 'Reactant '//trim(self%names_reactants_orig(i))//' was not found in '// &
                           trim(fpath)
            return
         endif
      enddo

      self%N_cond = self%N_reactants - self%N_gas

      call da_REORDER_SPECS(self)

      self%thermo_data_n_intervs = -1
      self%N_ions = 0
      self%reac_ion = .false.
      self%ions = .false.
      self%reac_atoms_names = ''
      self%reac_stoich = 0

      do i = 1,self%N_reactants
         call uppercase(self%names_reactants(i), name1)

         do j = 1,rl%nr
            call uppercase(rl%r(j)%name, name2)
            if (trim(name1) == trim(name2)) then

               ! Number of temperature intervals
               self%thermo_data_n_intervs(i) = rl%r(j)%thermo%ntemps

               ! Atoms on species
               kk = 1
               do k = 1,rl%natoms
                  if (rl%r(j)%composition(k) /= 0) then
                     if (kk > size(self%reac_atoms_names,1)) then
                        self%error = .true.
                        self%err_msg = 'The maximum number of atoms for a species is 5. '// &
                                       'The species '//trim(self%names_reactants(i))//' has more than this.'
                        return
                     endif

                     self%reac_atoms_names(kk,i) = trim(rl%atoms_names(k))
                     self%reac_stoich(kk,i) = rl%r(j)%composition(k)
                     if (trim(rl%atoms_names(k)) == 'E') then
                        self%ions = .true.
                        self%reac_ion(i) = .true.
                        self%N_ions = self%N_ions + 1
                     endif
                     kk = kk + 1
                  endif
               enddo

               self%mol_weight(i) = rl%r(j)%mass

               self%thermo_data_kind(i) = rl%r(j)%thermo%dtype

               ! Thermodynamic data
               do k = 1,rl%r(j)%thermo%ntemps
                  if (k > size(self%thermo_data_temps,2)) then
                     self%error = .true.
                     self%err_msg = 'The maximum number of thermodynamic intervals is 10. '// &
                                    'The species '//trim(self%names_reactants(i))//' has more than this.'
                     return
                  endif

                  self%thermo_data_temps(1,k,i) = rl%r(j)%thermo%temps(k)
                  self%thermo_data_temps(2,k,i) = rl%r(j)%thermo%temps(k+1)

                  if (rl%r(j)%thermo%dtype == ShomatePolynomial) then
                     self%thermo_data(1:7,k,i) = rl%r(j)%thermo%data(:,k)
                  elseif (rl%r(j)%thermo%dtype == Nasa9Polynomial) then
                     self%thermo_data(1:9,k,i) = rl%r(j)%thermo%data(:,k)
                     self%thermo_data(9:10,k,i) = self%thermo_data(8:9,k,i) 
                     self%thermo_data(8,k,i) = 0.0_dp
                  elseif (rl%r(j)%thermo%dtype == Nasa7Polynomial) then
                     self%thermo_data(1:7,k,i) = rl%r(j)%thermo%data(:,k)
                  endif
               enddo

               exit
            endif
         enddo
      enddo

   end subroutine

   function file_extension(filename) result(ext)
      character(*), intent(in) :: filename
      character(:), allocatable :: ext

      character(:), allocatable :: filename_trim
      integer :: i, j

      filename_trim = trim(filename)

      do i = len(filename_trim),1,-1
         if (filename_trim(i:i) == '.') exit
      enddo
      ext = ''
      do j = i+1,len(filename_trim)
         ext = ext//filename_trim(j:j)
      enddo 

   end function

   !> Sets id_atoms, where the i-th cell corresponds to names_atoms(i) and contains the index of the same atom in names_atoms_save
   !> names_atoms_save(id_atoms(i)) = names_atoms(i)
   subroutine da_ATOMS_ID(self)
      use equilibrate_const, only: N_atoms_save, names_atoms_save
      class(CEAData), intent(inout) :: self
      integer           :: i_atom, i_atom_save
      logical           :: change
      character(len=atom_str_len)  :: atom_upper, atom_upper_save
      character(len=3)  :: num

      do i_atom = 1, self%N_atoms
         if (trim(self%names_atoms(i_atom)) == '') then
            self%id_atoms(i_atom) = -1
            CYCLE
         end if
         change = .false.
         call uppercase(self%names_atoms(i_atom), atom_upper)
         do i_atom_save = 1, N_atoms_save
            call uppercase(names_atoms_save(i_atom_save), atom_upper_save)
            if (atom_upper == atom_upper_save) then
               self%id_atoms(i_atom) = i_atom_save
               change = .true.
               exit
            end if
         end do
         if (.not. change) then
            self%id_atoms(i_atom) = -1
            self%error = .true.
            write(num, '(i3.3)') i_atom
            self%err_msg = "READ DATA ERROR: the atom #"//num//" '"//self%names_atoms(i_atom)//"' was not recognised."
            RETURN
         end if
      end do
   end subroutine da_ATOMS_ID

   subroutine da_REAC_ATOMS_ID(self)
      use equilibrate_const, only: N_atoms_save, names_atoms_save
      class(CEAData), intent(inout) :: self
      integer           :: i_reac, i_atom, i_atom_save
      logical           :: change
      character(len=atom_str_len)  :: atom_upper, atom_upper_save
      character(len=3)  :: num

      do i_reac = 1, self%N_reactants
         do i_atom = 1, 5
            if (trim(self%reac_atoms_names(i_atom,i_reac))=='') then
               self%reac_atoms_id(i_atom, i_reac) = -1
               CYCLE
            end if
            change = .false.
            call uppercase(self%reac_atoms_names(i_atom,i_reac), atom_upper)
            do i_atom_save = 1, N_atoms_save
               call uppercase(names_atoms_save(i_atom_save), atom_upper_save)
               if (atom_upper == atom_upper_save) then
                  self%reac_atoms_id(i_atom,i_reac) = i_atom_save
                  change = .true.
                  EXIT
               end if
            end do
            if (.not. change) then
               self%reac_atoms_id(i_atom, i_reac) = -1
               self%error = .true.
               write(num, '(i3.3)') i_atom
               self%err_msg = "READ DATA ERROR: the atom #"//num//" '"//self%names_atoms(i_atom)//&
               "' in reactant '"//self%names_reactants(i_reac)//"' was not recognised."
               RETURN
            end if
         end do
      end do
   end subroutine da_REAC_ATOMS_ID

   subroutine da_BUILD_A_STOICH(self)
      use equilibrate_const, only: mol
      class(CEAData), intent(inout) :: self
      integer :: i_reac, i_ratom, i_atom, atom_id

      if (allocated(self%a_stoich)) then
         if (size(self%a_stoich,1) /= self%N_reactants .or. &
             size(self%a_stoich,2) /= self%N_atoms) then
            deallocate(self%a_stoich)
         end if
      end if
      if (.not. allocated(self%a_stoich)) then
         allocate(self%a_stoich(self%N_reactants,self%N_atoms))
      end if
      self%a_stoich = 0d0

      do i_reac = 1, self%N_reactants
         do i_ratom = 1, 5
            atom_id = self%reac_atoms_id(i_ratom,i_reac)
            if (atom_id > 0) then
               i_atom = findloc(self%id_atoms, atom_id, 1)
               if (i_atom == 0) then
                  self%error = .true.
                  self%err_msg = 'Stoichiometry mapping failed for reactant '// &
                                 trim(self%names_reactants(i_reac))
                  return
               end if
               self%a_stoich(i_reac,i_atom) = self%reac_stoich(i_ratom,i_reac)*mol
            end if
         end do
      end do
   end subroutine da_BUILD_A_STOICH

   !> Read in provided file all thermodynamic data
   subroutine da_READ_THERMO(self, fpath)
      use equilibrate_yaml, only: Nasa9Polynomial
      class(CEAData), intent(inout) :: self
      character(*), intent(in) :: fpath
      character(len=80)             :: file_line, file_line_up
      character(len=reac_str_len)             :: name_reac_up
      integer                       :: i_reac, n_interv, i_interv, i_stoich, stoich_start, reac_found

      reac_found = 0
      self%thermo_data_n_intervs = -1
      self%N_gas = 0
      self%N_ions = 0

      self%reac_ion = .FALSE.
      self%ions = .FALSE.

      self%thermo_data_kind = Nasa9Polynomial

      open(unit=17,file=fpath, action="READ", status="OLD")
      do while (1>0)
         read(17,'(A80)',end=122) file_line
         call uppercase(file_line, file_line_up)
         do i_reac = 1, self%N_reactants
            call uppercase(self%names_reactants_orig(i_reac), name_reac_up)

            if (trim(adjustl(file_line_up(1:18))) == trim(adjustl(name_reac_up))) then
               ! Recognized a reactant in file
               reac_found = reac_found+1

               ! Mark it as found
               read(17,'(A80)',end=122) file_line
               read(file_line(1:3),'(I2)') n_interv
               self%thermo_data_n_intervs(i_reac) = n_interv

               ! Gas or condensate ?
               if (file_line(52:52) == '0') then
                  self%reac_condensed(i_reac) = .FALSE.
                  self%N_gas = self%N_gas + 1
               else
                  self%reac_condensed(i_reac) = .TRUE.
               end if
            end if

         end do
      end do
      122 close(17)

      ! If not all reactants were found
      if (reac_found /= self%N_reactants) then
         self%error = .true.
         self%err_msg = 'READ DATA ERROR: For the following species no thermodynamical data was found:'
         n_interv = 0
         do i_reac = 1, self%N_reactants
            if (self%thermo_data_n_intervs(i_reac) .EQ. -1) then
               if (n_interv == 0) then
                  self%err_msg = trim(self%err_msg)//' '//trim(adjustl(self%names_reactants_orig(i_reac)))
               else
                  self%err_msg = trim(self%err_msg)//', '//trim(adjustl(self%names_reactants_orig(i_reac)))
               end if
            end if
         end do
         RETURN
      end if

      self%N_cond = self%N_reactants - self%N_gas

      ! Puts in self%names_reactants the ordered list of reactants, sets self%id_reactants as a link between the two
      call da_REORDER_SPECS(self)

      ! BASED ON THE ~RIGHT DESCRIPTION GIVEN IN GORDON 1996, page 73 AND THE APPEARANCE OF THERMO.INP.
      open(unit=17,file=fpath, action="READ", status="OLD")
      do while (1>0)
         read(17,'(A80)',end=123) file_line
         call uppercase(file_line, file_line_up)
         do i_reac = 1, self%N_reactants
            call uppercase(self%names_reactants(i_reac), name_reac_up)

            if (trim(adjustl(file_line_up(1:18))) == trim(adjustl(name_reac_up))) then
               ! Recognized a reactant in file
               read(17,'(A80)',end=123) file_line
               read(file_line(1:3),'(I2)') n_interv
               self%thermo_data_n_intervs(i_reac) = n_interv

               stoich_start = 11
               do i_stoich = 1, 5
                  ! Set the atoms for each reactant
                  self%reac_atoms_names(i_stoich,i_reac) = file_line(stoich_start:stoich_start+1)
                  ! Are there self%ions to be treated?
                  if (trim(adjustl(self%reac_atoms_names(i_stoich,i_reac))) == 'E') then
                     self%ions = .TRUE.
                     self%reac_ion(i_reac) = .TRUE.
                     self%N_ions = self%N_ions + 1
                  end if
                  read(file_line(stoich_start+2:stoich_start+7),'(F6.2)') self%reac_stoich(i_stoich,i_reac)
                  stoich_start = stoich_start+8
               end do

               read(file_line(53:65),'(F13.5)') self%mol_weight(i_reac)
               read(file_line(66:80),'(F13.5)') self%form_heat_Jmol_298_15_K(i_reac)

               do i_interv = 1, 3*n_interv
                  read(17,'(A80)',end=123) file_line
                  if (MOD(i_interv,3) == 1) then
                     read(file_line(2:22),'(F10.3,X,F10.3)') self%thermo_data_temps(1,i_interv/3+1,i_reac), &
                     self%thermo_data_temps(2,i_interv/3+1,i_reac)
                     read(file_line(23:23),'(I1)') self%thermo_data_n_coeffs(i_interv/3+1,i_reac)
                     read(file_line(24:63),'(8F5.1)') self%thermo_data_T_exps(1,i_interv/3+1,i_reac), &
                     self%thermo_data_T_exps(2,i_interv/3+1,i_reac), self%thermo_data_T_exps(3,i_interv/3+1,i_reac), &
                     self%thermo_data_T_exps(4,i_interv/3+1,i_reac), self%thermo_data_T_exps(5,i_interv/3+1,i_reac), &
                     self%thermo_data_T_exps(6,i_interv/3+1,i_reac), self%thermo_data_T_exps(7,i_interv/3+1,i_reac), &
                     self%thermo_data_T_exps(8,i_interv/3+1,i_reac)
                     read(file_line(66:80),'(F15.3)') self%H_0_298_15_K_m_H_0_0_K(i_interv/3+1,i_reac)
                  end if
                  if (MOD(i_interv,3) == 2) then
                     read(file_line(1:80),'(5D16.8)') self%thermo_data(1,i_interv/3+1,i_reac), &
                     self%thermo_data(2,i_interv/3+1,i_reac), self%thermo_data(3,i_interv/3+1,i_reac) , &
                     self%thermo_data(4,i_interv/3+1,i_reac), self%thermo_data(5,i_interv/3+1,i_reac)
                  end if
                  if (MOD(i_interv,3) == 0) then
                     read(file_line(1:80),'(5D16.8)') self%thermo_data(6,i_interv/3,i_reac), &
                     self%thermo_data(7,i_interv/3,i_reac), self%thermo_data(8,i_interv/3,i_reac) , &
                     self%thermo_data(9,i_interv/3,i_reac), self%thermo_data(10,i_interv/3,i_reac)
                  end if
               end do
            end if
         end do
      end do
      123 CLOSE(17)
   end subroutine da_READ_THERMO

   !> Puts every letter in strin in uppercase
   subroutine uppercase(strin,strout)

      character(len=*)  :: strin,strout
      integer           :: i

      strout = strin
      do i = 1, len(strout)
         select case(strout(i:i))
         case("a":"z")
            strout(i:i) = achar(iachar(strout(i:i))-32)
         end select
      end do
   end subroutine uppercase

   !> Puts in self%names_reactants the ordered list of reactants, sets self%id_reactants as a link between the two
   subroutine da_REORDER_SPECS(self)
      class(CEAData), intent(inout) :: self
      integer  :: i_reac, gas_offset, cond_offset

      gas_offset = 1
      cond_offset = 1
      do i_reac = 1, self%N_reactants
         if (self%reac_condensed(i_reac)) then
            self%names_reactants(self%N_gas+cond_offset) = self%names_reactants_orig(i_reac)
            ! 1st row = new index of the original reactant at i_reac
            self%id_reactants(i_reac,1) = self%N_gas + cond_offset
            ! 2nd row = original index of reactants in ordered list
            self%id_reactants(self%N_gas+cond_offset,2) = i_reac
            cond_offset = cond_offset + 1
         else
            self%names_reactants(gas_offset) = self%names_reactants_orig(i_reac)
            self%id_reactants(i_reac,1) = gas_offset
            self%id_reactants(gas_offset,2) = i_reac
            gas_offset = gas_offset + 1
         end if
      end do
   end subroutine da_REORDER_SPECS

   !> MAIN SUBROUTINE
   subroutine solve(self,mode,verbo,verbose2,N_atoms_in,N_reactants_in,molfracs_atoms,mass_tol, &
      molfracs_reactants,massfracs_reactants,temp,press,nabla_ad,gamma2,MMW,rho,c_pe, use_prev_guess)

      !! I/O:
      class(CEAData), intent(inout) :: self
      character, intent(in)            :: mode
      character(len=2), intent(in)     :: verbo
      logical, intent(in) :: verbose2
      real(dp), intent(in)     :: molfracs_atoms(N_atoms_in)
      real(dp), intent(in) :: mass_tol
      real(dp), intent(out)    :: molfracs_reactants(N_reactants_in), massfracs_reactants(N_reactants_in)
      integer, intent(in)              :: N_atoms_in, N_reactants_in
      real(dp), intent(in)     :: temp, press
      real(dp), intent(out)    :: nabla_ad,gamma2,MMW,rho,c_pe
      logical, intent(in) :: use_prev_guess

      !! Internal:
      real(dp)                 :: C_P_0(self%N_reactants), H_0(self%N_reactants), S_0(self%N_reactants)
      real(dp)                 :: molfracs_atoms_ions(N_atoms_in+1), temp_use
      integer                          :: N_atoms_use

      self%error = .false.

      if (self%N_atoms /= N_atoms_in+1 .or. self%N_reactants /= N_reactants_in) then
         self%error = .true.
         self%err_msg = "VALUE ERROR: The initialized and given arrays are not of the same size..."
         RETURN
      end if

      self%converged = .false.
      self%verbose = .FALSE.
      self%verbose_cond = (verbo == 'vy')
      self%verbose2 = verbose2
      self%quick = (mode /= 's')
      self%mass_tol = mass_tol
      self%remove_ions = .FALSE.

      ! call INIT_RAND_SEED()
      ! call random_init(repeatable=.true., image_distinct=.true.)

      molfracs_atoms_ions(1:N_atoms_in) = molfracs_atoms
      molfracs_atoms_ions(self%N_atoms) = 0d0

      if (self%ions) then
         if (temp > 750) then
            N_atoms_use = self%N_atoms
         else
            self%remove_ions = .true.
            N_atoms_use = N_atoms_in
         end if
      else
         N_atoms_use = N_atoms_in
      end if

      call ec_b_0(self,N_atoms_use, molfracs_atoms_ions(1:N_atoms_use), self%b_0_norm_cache, &
                 self%b_0_cache(1:N_atoms_use))
      self%b_0_cache_valid = .true.
      self%b_0_cache_atoms_use = N_atoms_use

      ! CALCULATION BEGINS

      call ec_COMP_THERMO_QUANTS(self,temp,self%N_reactants,C_P_0, H_0, S_0)
      gamma2 = 0d0
      temp_use = temp
      if (use_prev_guess .and. self%prev_guess_valid) then
         call ec_COMP_EQU_CHEM(self,N_atoms_use, self%N_reactants,  molfracs_atoms_ions(1:N_atoms_use), &
         molfracs_reactants, massfracs_reactants, &
         temp_use, press, C_P_0, H_0, S_0, &
         nabla_ad, gamma2, MMW, rho, c_pe, use_prev_guess=.true.)
         if (self%error) RETURN
         if (.not. self%converged) then
            call ec_COMP_EQU_CHEM(self,N_atoms_use, self%N_reactants,  molfracs_atoms_ions(1:N_atoms_use), &
            molfracs_reactants, massfracs_reactants, &
            temp_use, press, C_P_0, H_0, S_0, &
            nabla_ad, gamma2, MMW, rho, c_pe, use_prev_guess=.false.)
            if (self%error) RETURN
         end if
      else
         call ec_COMP_EQU_CHEM(self,N_atoms_use, self%N_reactants,  molfracs_atoms_ions(1:N_atoms_use), &
         molfracs_reactants, massfracs_reactants, &
         temp_use, press, C_P_0, H_0, S_0, &
         nabla_ad, gamma2, MMW, rho, c_pe, use_prev_guess=.false.)
         if (self%error) RETURN
      end if

      c_pe = c_pe*1d7 ! J/(g K) to erg/(g K)

   end subroutine

   !> Computes the values of C_P_0, H_0 and S_0
   subroutine ec_COMP_THERMO_QUANTS(self,temp,N_reac,C_P_0, H_0, S_0)
      use equilibrate_const, only: R
      use equilibrate_yaml, only: ShomatePolynomial, Nasa9Polynomial, Nasa7Polynomial
      !! I/O
      class(CEAData), intent(inout) :: self
      real(dp), intent(in)  :: temp
      integer, intent(in)           :: N_reac
      real(dp), intent(out) :: C_P_0(N_reac), H_0(N_reac), &
      S_0(N_reac)
      !! internal
      integer                       :: i_reac, i_tempinv, tempinv_ind

      do i_reac = 1, self%N_reactants

         ! Get temperature interpolation range
         if (temp < self%thermo_data_temps(1,1,i_reac)) then
            tempinv_ind = 1
         else if (temp >= self%thermo_data_temps(2,self%thermo_data_n_intervs(i_reac),i_reac)) then
            tempinv_ind = self%thermo_data_n_intervs(i_reac)
         else
            do i_tempinv = 1, self%thermo_data_n_intervs(i_reac)
               if ((temp >= self%thermo_data_temps(1,i_tempinv,i_reac)) .AND. &
               (temp < self%thermo_data_temps(2,i_tempinv,i_reac))) then
                  tempinv_ind = i_tempinv
                  EXIT
               end if
            end do
         end if

         if (self%thermo_data_kind(i_reac) == Nasa9Polynomial) then
         ! Calculate thermodynamic quantities as explained in Gordon 1996, page 74
         C_P_0(i_reac) = (self%thermo_data(1,tempinv_ind,i_reac)*temp**(-2)+ &
         self%thermo_data(2,tempinv_ind,i_reac)*temp**(-1)+ &
         self%thermo_data(3,tempinv_ind,i_reac)+self%thermo_data(4,tempinv_ind,i_reac)* &
         temp**(1)+self%thermo_data(5,tempinv_ind,i_reac)*temp**(2)+ &
         self%thermo_data(6,tempinv_ind,i_reac)*temp**(3)+self%thermo_data(7,tempinv_ind,i_reac)* &
         temp**(4))*R ! J/(K*mol)
         H_0(i_reac) = (-self%thermo_data(1,tempinv_ind,i_reac)*temp**(-2)+ &
         self%thermo_data(2,tempinv_ind,i_reac)*temp**(-1)*log(temp)+ &
         self%thermo_data(3,tempinv_ind,i_reac)+self%thermo_data(4,tempinv_ind,i_reac)*temp**(1)/2d0+ &
         self%thermo_data(5,tempinv_ind,i_reac)*temp**(2)/3d0+ &
         self%thermo_data(6,tempinv_ind,i_reac)*temp**(3)/4d0+self%thermo_data(7,tempinv_ind,i_reac)* &
         temp**(4)/5d0+self%thermo_data(9,tempinv_ind,i_reac)/temp)* &
         R*temp ! J/mol
         S_0(i_reac) = (-self%thermo_data(1,tempinv_ind,i_reac)*temp**(-2)/2d0- &
         self%thermo_data(2,tempinv_ind,i_reac)*temp**(-1)+ &
         self%thermo_data(3,tempinv_ind,i_reac)*log(temp)+ &
         self%thermo_data(4,tempinv_ind,i_reac)*temp**(1)+ &
         self%thermo_data(5,tempinv_ind,i_reac)*temp**(2)/2d0+ &
         self%thermo_data(6,tempinv_ind,i_reac)*temp**(3)/3d0+self%thermo_data(7,tempinv_ind,i_reac)* &
         temp**(4)/4d0+self%thermo_data(10,tempinv_ind,i_reac))*R ! J/(K*mol)

         elseif (self%thermo_data_kind(i_reac) == ShomatePolynomial) then

         C_P_0(i_reac) = heat_capacity_shomate(self%thermo_data(1:7,tempinv_ind,i_reac), temp) ! J/(K*mol)
         H_0(i_reac) = enthalpy_shomate(self%thermo_data(1:7,tempinv_ind,i_reac), temp) ! J/mol
         S_0(i_reac) = entropy_shomate(self%thermo_data(1:7,tempinv_ind,i_reac), temp) ! J/(K*mol)

         elseif (self%thermo_data_kind(i_reac) == Nasa7Polynomial) then

         C_P_0(i_reac) = heat_capacity_nasa7(self%thermo_data(1:7,tempinv_ind,i_reac), temp) ! J/(K*mol)
         H_0(i_reac) = enthalpy_nasa7(self%thermo_data(1:7,tempinv_ind,i_reac), temp) ! J/mol
         S_0(i_reac) = entropy_nasa7(self%thermo_data(1:7,tempinv_ind,i_reac), temp) ! J/(K*mol)

         endif

      end do

   end subroutine ec_COMP_THERMO_QUANTS

   pure function heat_capacity_shomate(coeffs, T) result(cp)
      real(dp), intent(in) :: coeffs(7)
      real(dp), intent(in) :: T
      real(dp) :: cp !! J/(mol K)
      
      real(dp) :: TT
      
      TT = T/1000.0_dp
      cp = coeffs(1) + coeffs(2)*TT + coeffs(3)*TT**2 + &
            coeffs(4)*TT**3 + coeffs(5)/TT**2
   end function

   pure function enthalpy_shomate(coeffs, T) result(enthalpy)
      real(dp), intent(in) :: coeffs(7)
      real(dp), intent(in) :: T
      real(dp) :: enthalpy !! J/mol

      real(dp) :: TT

      TT = T/1000.0_dp
      enthalpy = coeffs(1)*TT + (coeffs(2)*TT**2)/2.0_dp &
               + (coeffs(3)*TT**3)/3.0_dp  + (coeffs(4)*TT**4)/4.0_dp &
               - coeffs(5)/TT + coeffs(6)
      enthalpy = enthalpy*1000.0_dp
   end function

   pure function entropy_shomate(coeffs, T) result(entropy)
      real(dp), intent(in) :: coeffs(7)
      real(dp), intent(in) :: T
      real(dp) :: entropy !! J/(mol K)
      
      real(dp) :: TT
      
      TT = T/1000.0_dp
      entropy = coeffs(1)*log(TT) + coeffs(2)*TT &
               + (coeffs(3)*TT**2)/2.0_dp + (coeffs(4)*TT**3)/3.0_dp &
               - coeffs(5)/(2.0_dp * TT**2) + coeffs(7)
   end function

   pure function heat_capacity_nasa7(coeffs, T) result(cp)
      use equilibrate_const, only: R
      real(dp), intent(in) :: coeffs(:)
      real(dp), intent(in) :: T
      real(dp) :: cp !! J/(mol K)

      real(dp) :: a0, a1, a2, a3, a4

      a0 = coeffs(1)
      a1 = coeffs(2)
      a2 = coeffs(3)
      a3 = coeffs(4)
      a4 = coeffs(5)

      cp = a0 + a1*T + a2*T**2.0_dp + a3*T**3.0_dp + a4*T**4.0_dp
      cp = cp*R

   end function

   pure function enthalpy_nasa7(coeffs, T) result(enthalpy)
      use equilibrate_const, only: R
      real(dp), intent(in) :: coeffs(:)
      real(dp), intent(in) :: T
      real(dp) :: enthalpy !! J/mol

      real(dp) :: a0, a1, a2, a3, a4, a5

      a0 = coeffs(1)
      a1 = coeffs(2)
      a2 = coeffs(3)
      a3 = coeffs(4)
      a4 = coeffs(5)
      a5 = coeffs(6)

      enthalpy = &
        a0 + &
        a1*T/2.0_dp + &
        a2*T**2.0_dp/3.0_dp + &
        a3*T**3.0_dp/4.0_dp + &
        a4*T**4.0_dp/5.0_dp + &
        a5/T
      enthalpy = enthalpy*R*T
  
   end function

   pure function entropy_nasa7(coeffs, T) result(entropy)
      use equilibrate_const, only: R
      real(dp), intent(in) :: coeffs(:)
      real(dp), intent(in) :: T
      real(dp) :: entropy !! J/(mol K)

      real(dp) :: a0, a1, a2, a3, a4, a6

      a0 = coeffs(1)
      a1 = coeffs(2)
      a2 = coeffs(3)
      a3 = coeffs(4)
      a4 = coeffs(5)
      a6 = coeffs(7)
      
      entropy = &
        a0*log(T) + &
        a1*T + &
        a2*T**2.0_dp/2.0_dp + &
        a3*T**3.0_dp/3.0_dp + &
        a4*T**4.0_dp/4.0_dp + &
        a6
      entropy = entropy*R

   end function

   !> Computes the specie abundances (molar and mass)
   recursive subroutine ec_COMP_EQU_CHEM(self, N_atoms_use, N_reac, molfracs_atoms, &
      molfracs_reactants, massfracs_reactants, &
      temp, press, C_P_0, H_0, S_0, nabla_ad, gamma2, MMW, rho, c_pe, use_prev_guess)
      use equilibrate_const, only: amu, kB, masses_atoms_save

      !! I/O:
      class(CEAData), intent(inout) :: self
      integer, intent(in)              :: N_atoms_use, N_reac
      real(dp), intent(in)     :: molfracs_atoms(N_atoms_use)
      real(dp), intent(inout)  :: molfracs_reactants(N_reac), massfracs_reactants(N_reac)
      real(dp), intent(in)     :: temp, press
      real(dp), intent(in)     :: C_P_0(N_reac), H_0(N_reac), S_0(N_reac)
      real(dp), intent(out)    :: nabla_ad, gamma2, MMW, rho, c_pe
      logical, intent(in) :: use_prev_guess

      !! CEA McBride 1994 style variables:
      real(dp)  :: n ! Moles of gas particles per total mass of mixture in kg
      real(dp)  :: n_spec(self%N_reactants) ! Moles of species per total mass of mixture in kg
      real(dp)  :: n_spec_old(self%N_reactants) ! Moles of species per total mass of mixture in kg of previous iteration
      real(dp)  :: pi_atom(N_atoms_use) ! Lagrangian multipliers for the atomic species divided by (R*T)
      real(dp)  :: matrix(self%N_reactants+N_atoms_use+1,self%N_reactants+N_atoms_use+1)
      ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for condensed species, the pis and the delta log(n)
      real(dp)  :: vector(self%N_reactants+N_atoms_use+1), solution_vector(self%N_reactants+N_atoms_use+1)

      !! Internal:
      INTEGER           :: i_iter, i_reac, inc_next, current_solids_number, N_spec_eff, buffer_ind, i_atom
      LOGICAL           :: converged, remove_cond, slowed
      LOGICAL           :: solid_inclu(self%N_cond), neg_cond(self%N_cond)
      real(dp)  :: dgdnj(self%N_cond)
      INTEGER           :: solid_indices(self%N_cond), solid_indices_buff(self%N_cond)
      real(dp)  :: nsum, mu_gas(self%N_gas), a_gas(self%N_gas,N_atoms_use), mass_species, atom_mass, msum

      converged = .FALSE.
      slowed = .FALSE.
      call ec_INIT_ALL_VALS(self,N_atoms_use,self%N_reactants,n,n_spec,pi_atom,use_prev_guess)

      self%iter_max = 50 + self%N_reactants/2
      current_solids_number = 0

      MMW = 0d0

      n_spec_old = n_spec

      ! FIRST: DO GAS ONLY!
      DO i_iter = 1, self%iter_max

         ! IF (self%quick) THEN
         call ec_PREP_MATRIX_SHORT(self,N_atoms_use,self%N_reactants, molfracs_atoms,self%N_gas,press,temp,&
         H_0,S_0,n,n_spec,matrix(1:N_atoms_use+1,1:N_atoms_use+1),vector(1:N_atoms_use+1),&
         (/1,1,1,1,1/),5, mu_gas,a_gas)
         call ec_INVERT_MATRIX_SHORT(self,N_atoms_use+1, &
         matrix(1:N_atoms_use+1,1:N_atoms_use+1),vector(1:N_atoms_use+1), &
         solution_vector(1:N_atoms_use+1))
         if (self%error) RETURN

         call ec_UPDATE_ABUNDS_SHORT(self,N_atoms_use,self%N_reactants,self%N_gas,&
         solution_vector(1:N_atoms_use+1), n_spec,pi_atom,n,converged,&
         (/1,1,1,1,1/),5,mu_gas,a_gas,temp,molfracs_atoms,n_spec_old)
         ! ELSE
         !    call ec_PREP_MATRIX_LONG(N_atoms,self%id_atoms,molfracs_atoms,self%N_gas,&
         !    press,temp,H_0,S_0,n,n_spec,&
         !    matrix(1:self%N_gas+N_atoms+1,1:self%N_gas+N_atoms+1),vector(1:self%N_gas+N_atoms+1),&
         !    (/1,1,1,1,1/),self%names_reactants,N_reactants,5)
         !    call ec_INVERT_MATRIX_LONG(N_atoms+self%N_gas+1, &
         !    matrix(1:self%N_gas+N_atoms+1,1:self%N_gas+N_atoms+1),vector(1:self%N_gas+N_atoms+1), &
         !    solution_vector(1:self%N_gas+N_atoms+1))
         !    call ec_UPDATE_ABUNDS_LONG(N_atoms,self%N_gas,solution_vector(1:self%N_gas+N_atoms+1), &
         !    n_spec,pi_atom,n,converged,(/1,1,1,1,1/),5,self%id_atoms,molfracs_atoms,N_reactants,&
         !    n_spec_old)
         ! END IF

         n_spec_old = n_spec

         IF (self%verbose) THEN
            write(*,*)
            write(*,*)
            write(*,*) i_iter
            DO i_reac = 1, self%N_reactants
               write(*,*) self%names_reactants(i_reac), n_spec(i_reac)/SUM(n_spec)
            END DO
         END IF

         IF (converged) THEN
            EXIT
         END IF
      END DO

      self%converged = converged
      IF (.NOT. converged) THEN
        if (self%verbose2) then
          print*,'One or more convergence criteria not satisfied (gas only).'
        endif
      END IF

      converged = .FALSE.
      remove_cond = .FALSE.

      IF (self%N_gas .EQ. self%N_reactants) THEN
         call ec_COMP_ADIABATIC_GRAD(self,N_atoms_use, self%N_reactants, self%N_gas,  n_spec, &
         n, H_0, C_P_0, (/ 1,1,1,1,1 /), 5, temp, nabla_ad, gamma2, c_pe)
         if (self%error) RETURN
      END IF

      ! THEN: INCLUDE CONDENSATES!
      IF (self%N_cond > 0) THEN
         solid_inclu = .FALSE.
         solid_indices = 0
         inc_next = 0
         neg_cond = .FALSE.

         N_spec_eff = self%N_gas

         DO WHILE (inc_next /= -1)

            call ec_INCLUDE_WHICH_SOLID(self,N_atoms_use,N_reac,pi_atom,H_0,S_0,temp, &
            n_spec,solid_inclu,neg_cond,dgdnj,remove_cond,inc_next)

            if (inc_next/=-1) then

               IF (remove_cond) THEN
                  solid_indices_buff = 0
                  buffer_ind = 1
                  do i_reac = 1, current_solids_number
                     IF (solid_indices(i_reac) .NE. inc_next) THEN
                        solid_indices_buff(buffer_ind) = solid_indices(i_reac)
                        buffer_ind = buffer_ind + 1
                     END IF
                  END DO
                  current_solids_number = buffer_ind - 1
                  solid_indices = solid_indices_buff
                  solid_inclu(inc_next-self%N_gas) = .FALSE.
                  neg_cond(inc_next-self%N_gas) = .TRUE.
                  if (self%verbose_cond) then
                     print *, '-   ', self%names_reactants(inc_next), dgdnj(inc_next-self%N_gas), n_spec(inc_next)
                  end if
                  n_spec(inc_next) = 0d0
               ELSE
                  current_solids_number = current_solids_number + 1
                  solid_indices(current_solids_number) = inc_next
                  solid_inclu(inc_next-self%N_gas) = .TRUE.
                  call ec_INIT_COND_VALS(self,N_atoms_use, self%N_reactants, molfracs_atoms, inc_next, n_spec)
                  if (self%error) return
                  
                  if (self%verbose_cond) then
                     print *, '+   ', self%names_reactants(inc_next), dgdnj(inc_next-self%N_gas), n_spec(inc_next)
                  end if
               END IF

               N_spec_eff = self%N_gas+current_solids_number
               DO i_iter = 1, self%iter_max

                  IF (self%quick) THEN
                     call ec_PREP_MATRIX_SHORT(self,N_atoms_use, self%N_reactants, molfracs_atoms,N_spec_eff, &
                     press, temp, H_0, S_0, n, n_spec, &
                     matrix(1:N_atoms_use+1+N_spec_eff-self%N_gas,1:N_atoms_use+1+N_spec_eff-self%N_gas), &
                     vector(1:N_atoms_use+1+N_spec_eff-self%N_gas), solid_indices, &
                     N_spec_eff-self%N_gas, mu_gas, a_gas)
                     call ec_INVERT_MATRIX_SHORT(self,N_atoms_use+1+N_spec_eff-self%N_gas, &
                     matrix(1:N_atoms_use+1+N_spec_eff-self%N_gas,1:N_atoms_use+1+N_spec_eff-self%N_gas), &
                     vector(1:N_atoms_use+1+N_spec_eff-self%N_gas), &
                     solution_vector(1:N_atoms_use+1+N_spec_eff-self%N_gas))
                     if (self%error) RETURN

                     call ec_UPDATE_ABUNDS_SHORT(self,N_atoms_use,self%N_reactants,N_spec_eff,&
                     solution_vector(1:N_atoms_use+1+N_spec_eff-self%N_gas), &
                     n_spec,pi_atom,n,converged,&
                     solid_indices,N_spec_eff-self%N_gas,mu_gas,a_gas,temp,molfracs_atoms, &
                     n_spec_old)
                  ELSE
                     call ec_PREP_MATRIX_LONG(self,N_atoms_use,self%N_reactants,molfracs_atoms,N_spec_eff,&
                     press,temp,H_0,S_0,n,n_spec,&
                     matrix(1:N_spec_eff+N_atoms_use+1,1:N_spec_eff+N_atoms_use+1),&
                     vector(1:N_spec_eff+N_atoms_use+1),&
                     solid_indices,N_spec_eff-self%N_gas)
                     call ec_INVERT_MATRIX_LONG(self,N_atoms_use+N_spec_eff+1, &
                     matrix(1:N_spec_eff+N_atoms_use+1,1:N_spec_eff+N_atoms_use+1),&
                     vector(1:N_spec_eff+N_atoms_use+1), &
                     solution_vector(1:N_spec_eff+N_atoms_use+1))
                     if (self%error) RETURN

                     call ec_UPDATE_ABUNDS_LONG(self,N_atoms_use, N_reac, N_spec_eff, &
                     solution_vector(1:N_spec_eff+N_atoms_use+1), &
                     n_spec, pi_atom, n, converged, &
                     solid_indices, N_spec_eff-self%N_gas, molfracs_atoms, n_spec_old)
                  END IF

                  ! call writetxtall(self%N_reactants, n_spec)
                  n_spec_old = n_spec

                  IF (self%verbose) THEN
                     write(*,*)
                     write(*,*)
                     write(*,*) i_iter
                     DO i_reac = 1, self%N_reactants
                        write(*,*) self%names_reactants(i_reac), n_spec(i_reac)/SUM(n_spec)
                     END DO
                  END IF

                  DO i_reac = self%N_gas+1, self%N_reactants
                     IF ((n_spec(i_reac) < 0d0) .AND. (i_iter > 30)) THEN
                        converged = .TRUE.
                        EXIT
                     END IF
                  END DO

                  IF (converged) THEN
                     EXIT
                  END IF

               END DO

               self%converged = converged
               IF (.NOT. converged) THEN
                  IF (self%quick) THEN
                     self%quick = .FALSE.
                     if (self%verbose2) then
                        print*,'SLOW!'
                     endif
                     call ec_COMP_EQU_CHEM(self,N_atoms_use, self%N_reactants, molfracs_atoms, &
                     molfracs_reactants, massfracs_reactants, &
                     temp, press, C_P_0, H_0, S_0, &
                     nabla_ad,gamma2,MMW,rho,c_pe, use_prev_guess)
                     if (self%error) return
                     self%quick = .TRUE.
                     slowed = .TRUE.
                     EXIT
                  ELSE
                    if (self%verbose2) then
                       print*,'One or more convergence criteria not satisfied (gas & condensates).'
                    endif
                  END IF
               END IF

               converged = .FALSE.

            END IF

            remove_cond = .FALSE.

         END DO

         ! Calc. nabla_ad
         IF (.NOT. slowed) THEN
            call ec_COMP_ADIABATIC_GRAD(self,N_atoms_use, self%N_reactants, N_spec_eff, n_spec, &
            n,H_0,C_P_0,solid_indices,N_spec_eff-self%N_gas,temp, nabla_ad,gamma2,c_pe)
            if (self%error) RETURN
         END IF

         if (self%verbose_cond) then
            print *
            print *, 'Solids included:'
            do i_reac = 1, self%N_reactants-self%N_gas
               if (solid_inclu(i_reac)) then
                  print *, ' ', self%names_reactants(i_reac+self%N_gas), n_spec(i_reac+self%N_gas)
               end if
            end do
         end if

      END IF

      ! PREPARE FINAL OUTPUT

      IF (.NOT. slowed) THEN
        !  nsum = SUM(n_spec)
        !  DO i_reac = 1, self%N_reactants
        !     IF (n_spec(i_reac)/nsum < 1d-50) THEN
        !        n_spec(i_reac) = 0d0
        !     END IF
        !  END DO

         nsum = SUM(n_spec)
         do i_reac = 1, self%N_reactants
            molfracs_reactants(i_reac) = n_spec(self%id_reactants(i_reac,1))/nsum
         end do

         msum = 0d0
         DO i_reac = 1, self%N_reactants
            mass_species = 0d0
            DO i_atom = 1, 5
               if (self%reac_atoms_id(i_atom,i_reac)>0) then
                  atom_mass = masses_atoms_save(self%reac_atoms_id(i_atom,i_reac))
                  mass_species = mass_species+atom_mass*DBLE(self%reac_stoich(i_atom,i_reac))
               END IF
            END DO
            massfracs_reactants(self%id_reactants(i_reac,2)) = n_spec(i_reac) * mass_species
            if (i_reac <= self%N_gas) then
               MMW = MMW + massfracs_reactants(self%id_reactants(i_reac,2))/mass_species
               msum = msum + massfracs_reactants(self%id_reactants(i_reac,2))
            end if
         END DO
         massfracs_reactants = massfracs_reactants / SUM(massfracs_reactants)

         ! Mean molecular weight is only calculated as mean molecular gas weight per all gas species
         MMW = MMW/msum
         MMW = 1d0/MMW

         rho = (press*1d6)/kB/temp*MMW
         MMW = MMW/amu

      END IF

      if (self%converged) then
         if (.not. allocated(self%prev_n_spec)) then
            allocate(self%prev_n_spec(self%N_reactants))
         end if
         self%prev_n_spec = n_spec
         self%prev_n = n
         self%prev_guess_valid = .true.
      end if

   end subroutine ec_COMP_EQU_CHEM

   !> Initialize all abundances with uniform abundances for gas and 0 for condensates
   subroutine ec_INIT_ALL_VALS(self,N_atoms_use,N_reac,n,n_spec,pi_atom,use_prev_guess)
      !! I/O:
      class(CEAData), intent(inout) :: self
      integer, intent(in)           :: N_atoms_use, N_reac
      real(dp), intent(out) :: n ! Moles of gas particles per total mass of mixture in kg
      real(dp), intent(out) :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg
      real(dp), intent(out) :: pi_atom(N_atoms_use) ! Lagrangian multipliers for the atomic species divided
      ! by (R*T)
      logical, intent(in) :: use_prev_guess

      !! Internal:
      INTEGER                      :: i_reac

      if (use_prev_guess .and. self%prev_guess_valid .and. &
          allocated(self%prev_n_spec) .and. size(self%prev_n_spec) == N_reac) then
         n_spec = self%prev_n_spec
         n = self%prev_n
         if (self%remove_ions) then
            DO i_reac = 1, self%N_gas
               if (self%reac_ion(i_reac)) then
                  n_spec(i_reac) = 0d0
               end if
            END DO
            n = max(sum(n_spec(1:self%N_gas)), tiny(1.0_dp))
         else
            n = max(n, sum(n_spec(1:self%N_gas)), tiny(1.0_dp))
         end if
      else
         n = 0.1d0
         n_spec = 0d0
         pi_atom = 0d0
         DO i_reac = 1, self%N_gas
            n_spec(i_reac) = n/DBLE(self%N_gas)
            IF (self%remove_ions) THEN
               IF(self%reac_ion(i_reac)) THEN
                  n_spec(i_reac) = 0d0
               END IF
            END IF
         END DO
      end if
      pi_atom = 0d0

   end subroutine ec_INIT_ALL_VALS

   !> Selects which solid to include next
   subroutine ec_INCLUDE_WHICH_SOLID(self,N_atoms_use,N_reac,pi_atom,H_0,S_0,temp, &
      n_spec,solid_inclu,neg_cond,dgdnj,remove_cond,inc_next)
      use equilibrate_const, only: mol, R
      !! I/O
      class(CEAData), intent(inout) :: self
      integer, intent(in)           :: N_atoms_use, N_reac
      real(dp), intent(in)  :: pi_atom(N_atoms_use)
      real(dp), intent(in)  :: H_0(N_reac), S_0(N_reac), temp
      LOGICAL, intent(in)           :: solid_inclu(N_reac-self%N_gas), neg_cond(N_reac-self%N_gas)

      real(dp), intent(inout)  :: n_spec(N_reac), dgdnj(N_reac-self%N_gas)
      logical, intent(out)          :: remove_cond
      integer, intent(out)          :: inc_next

      !! Internal:
      real(dp)             :: a(self%N_reactants,N_atoms_use)
      real(dp)             :: mu(self%N_reactants), minval_inc
      INTEGER                      :: i_atom, i_reac, i_ratom, remove_count, remove_id

      !f2py integer, intent(aux) :: self%N_gas

      remove_count = 0
      remove_id = 0
      remove_cond = .false.
      DO i_reac = self%N_gas+1, self%N_reactants
         IF (n_spec(i_reac) < 0d0) THEN
            ! n_spec(i_reac) = 0d0
            remove_count = remove_count + 1
            remove_id = i_reac
            ! remove_cond = .TRUE.
         END IF
      END DO

      if (remove_count >= 2) then
         inc_next = remove_id
         remove_cond = .true.
      else

         inc_next = -1
         minval_inc = 2d0

         ! Set up a_ij
         a = 0d0
         DO i_atom = 1, N_atoms_use
            DO i_reac = 1, self%N_reactants
               IF (self%remove_ions) THEN
                  IF (self%reac_ion(i_reac)) THEN
                     CYCLE
                  END IF
               END IF
               DO i_ratom = 1, 5
                  IF (self%reac_atoms_id(i_ratom, i_reac)>0 .and. self%id_atoms(i_atom) == self%reac_atoms_id(i_ratom, i_reac)) then
                     a(i_reac,i_atom) = self%reac_stoich(i_ratom,i_reac)*mol
                  END IF
               END DO
            END DO
         END DO

         ! EVAL Eq. 3.7 in McBride Manual
         DO i_reac = self%N_gas+1, self%N_reactants
            mu(i_reac) = H_0(i_reac) - temp*S_0(i_reac)
            dgdnj(i_reac-self%N_gas) = mu(i_reac)/R/temp - SUM(a(i_reac,1:N_atoms_use)*pi_atom)
         END DO

         DO i_reac = self%N_gas+1, self%N_reactants
            IF ((dgdnj(i_reac-self%N_gas) < 0d0) .AND. (.NOT. solid_inclu(i_reac-self%N_gas))) THEN
               IF (((dgdnj(i_reac-self%N_gas) < minval_inc) .AND. (.NOT. neg_cond(i_reac-self%N_gas)) .AND. &
               temp <= self%thermo_data_temps(2,self%thermo_data_n_intervs(i_reac),i_reac)) .AND. &
               (temp >= self%thermo_data_temps(1,1,i_reac))) THEN
                  minval_inc = dgdnj(i_reac-self%N_gas)
                  inc_next = i_reac
               END IF
            END IF
         END DO

         if (inc_next==-1 .and. remove_count==1) then
            inc_next = remove_id
            remove_cond = .true.
         end if

      end if

   end subroutine ec_INCLUDE_WHICH_SOLID

   !> Initialize one condensate abundance, as if the most molecules condensed
   subroutine ec_INIT_COND_VALS(self,N_atoms_use, N_reac, molfracs_atoms, i_cond, n_spec)
      class(CEAData), intent(inout) :: self
      integer, intent(in)              :: N_atoms_use, i_cond, N_reac
      real(dp), intent(in)     :: molfracs_atoms(N_atoms_use)
      real(dp), intent(inout)  :: n_spec(N_reac)

      integer           :: i_ratom, i_atom
      real(dp)  :: min_molfrac, stoich_molfrac

      min_molfrac = -1
      do i_ratom = 1, 5
         if (self%reac_atoms_id(i_ratom,i_cond) > 0) then
            do i_atom = 1, N_atoms_use
               if (self%id_atoms(i_atom) == self%reac_atoms_id(i_ratom,i_cond)) then
                  stoich_molfrac = molfracs_atoms(i_atom) / self%reac_stoich(i_ratom,i_cond)
                  if (min_molfrac==-1 .or. stoich_molfrac < min_molfrac) then
                     min_molfrac = stoich_molfrac
                  end if
               end if
            end do
         end if
      end do

      if (min_molfrac >= 0) then
         n_spec(i_cond) = min_molfrac
      else
         self%error = .true.
         self%err_msg = 'No data found for the atoms of the given condensate.'
         return
      end if

   end subroutine ec_INIT_COND_VALS

   !> Build the small matrix
   subroutine ec_PREP_MATRIX_SHORT(self,N_atoms_use, N_reac, molfracs_atoms, N_species, press, temp, &
   H_0, S_0, n, n_spec, matrix, vector, solid_indices, N_solids, mu_gas, a_gas)
      use equilibrate_const, only: R
      !! I/O:
      class(CEAData), intent(inout) :: self
      INTEGER, intent(in)          :: N_atoms_use, N_reac, N_species, N_solids
      INTEGER, intent(in)          :: solid_indices(N_solids)
      real(dp), intent(in) :: molfracs_atoms(N_atoms_use), press, temp
      real(dp), intent(in) ::  H_0(N_reac), S_0(N_reac)
      real(dp), intent(in) :: n ! Moles of gas particles per total mass of mixture in kg
      real(dp), intent(inout)  :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg

      real(dp), intent(out):: matrix(N_atoms_use+1+(N_species-self%N_gas),N_atoms_use+1+(N_species-self%N_gas))
      ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
      ! condensed species, the pis and the delta log(n)
      real(dp), intent(out):: vector(N_atoms_use+1+(N_species-self%N_gas))
      real(dp), intent(out):: mu_gas(self%N_gas), a_gas(self%N_gas,N_atoms_use)

      !! Internal:
      real(dp)             :: b_0(N_atoms_use), b_0_norm, b(N_atoms_use)
      real(dp)             :: mu_solids(max(1,N_solids))
      real(dp)             :: inv_RT, logp, sum_n_gas, sum_n_mu, a_ia
      INTEGER                      :: i_atom, i_reac, i_atom2, i_solid, solid_idx

      !f2py integer, intent(aux) :: self%N_gas

      ! print *, "START ec_PREP_MATRIX_SHORT"

      ! Set up b0
      ! b_0_norm = 0d0
      ! DO i_atom = 1, N_atoms_use
      !    ! call ec_ATOM_MASS(self%names_atoms(i_atom),mass_atom)
      !    mass_atom = masses_atoms_save(self%id_atoms(i_atom))
      !    b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
      ! END DO
      ! b_0 = molfracs_atoms/b_0_norm
      if (self%b_0_cache_valid .and. self%b_0_cache_atoms_use == N_atoms_use) then
         b_0 = self%b_0_cache(1:N_atoms_use)
         b_0_norm = self%b_0_norm_cache
      else
         call ec_b_0(self,N_atoms_use, molfracs_atoms, b_0_norm, b_0)
      end if

      ! Set up mu_j
      inv_RT = 1d0/(R*temp)
      logp = log(press)
      mu_gas = 0d0
      if (N_solids > 0) then
         mu_solids = 0d0
      end if
      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         ! Taken from Venot et al. (2012), in comparison with McBride 1996.
         mu_gas(i_reac) = H_0(i_reac) - temp*S_0(i_reac)

         IF (n_spec(i_reac) > 1d-290) THEN
            mu_gas(i_reac) = mu_gas(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*logp
         ELSE
            IF (self%verbose) THEN
               write(*,*) 'n_spec(i_reac) == 0 for '//trim(adjustl(self%names_reactants(i_reac)))// &
               ' set to 1d-13 and try again.'
            END IF
            n_spec(i_reac) = 1d-13*(1d0 + 1d-6*DBLE(i_reac))
            mu_gas(i_reac) = mu_gas(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*logp
         END IF
      END DO
      if (N_solids > 0) then
         DO i_solid = 1, N_solids
            solid_idx = solid_indices(i_solid)
            mu_solids(i_solid) = H_0(solid_idx) - temp*S_0(solid_idx)
         END DO
      end if

      a_gas = self%a_stoich(1:self%N_gas,1:N_atoms_use)

      ! MATRIX SETUP
      matrix = 0d0
      b = 0d0
      sum_n_gas = 0d0
      sum_n_mu = 0d0

      ! Set up the matrix for the N_atoms equations (Eq. 2.24)
      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         sum_n_gas = sum_n_gas + n_spec(i_reac)
         sum_n_mu = sum_n_mu + n_spec(i_reac)*mu_gas(i_reac)*inv_RT
         DO i_atom = 1, N_atoms_use
            a_ia = self%a_stoich(i_reac,i_atom)
            b(i_atom) = b(i_atom) + a_ia*n_spec(i_reac)
            matrix(i_atom,N_atoms_use+1) = matrix(i_atom,N_atoms_use+1) + a_ia*n_spec(i_reac)
            DO i_atom2 = 1, N_atoms_use
               matrix(i_atom,i_atom2) = matrix(i_atom,i_atom2) + a_ia*self%a_stoich(i_reac,i_atom2)*n_spec(i_reac)
            END DO
         END DO
      END DO
      matrix(N_atoms_use+1,1:N_atoms_use) = matrix(1:N_atoms_use,N_atoms_use+1)
      matrix(N_atoms_use+1,N_atoms_use+1) = sum_n_gas - n

      IF (self%N_gas < N_species) THEN
         DO i_solid = 1, N_species-self%N_gas
            solid_idx = solid_indices(i_solid)
            DO i_atom = 1, N_atoms_use
               a_ia = self%a_stoich(solid_idx,i_atom)
               b(i_atom) = b(i_atom) + a_ia*n_spec(solid_idx)
               matrix(i_atom,N_atoms_use+1+i_solid) = a_ia
               matrix(N_atoms_use+1+i_solid,i_atom) = a_ia
            END DO
         END DO
      END IF

      ! Set up the matrix for the (self%N_reactants-self%N_gas) equations (Eq. 2.25)

      ! VECTOR SETUP
      !vector(N_atoms+1+(self%N_reactants-self%N_gas))
      vector = 0d0

      ! (Eq. 2.25)
      IF (self%N_gas < N_species) THEN
         vector(N_atoms_use+2:N_atoms_use+1+(N_species-self%N_gas)) = mu_solids(1:N_species-self%N_gas)*inv_RT
      END IF

      ! (Eq. 2.24)
      vector(1:N_atoms_use) = b_0 - b
      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         vector(1:N_atoms_use) = vector(1:N_atoms_use) + &
         self%a_stoich(i_reac,1:N_atoms_use)*n_spec(i_reac)*mu_gas(i_reac)*inv_RT
      END DO

      ! (Eq. 2.26)
      vector(N_atoms_use+1) = n - sum_n_gas + sum_n_mu

   end subroutine ec_PREP_MATRIX_SHORT



   !> Return the result of one step of computation with small matrix(problem: AX=B)
   subroutine ec_UPDATE_ABUNDS_SHORT(self,N_atoms_use,N_reac,N_species,solution_vector,n_spec,pi_atom,&
   n,converged,solid_indices,N_solids,mu_gas,a_gas,temp,molfracs_atoms,n_spec_old)
      use equilibrate_const, only: mol, R
      !! I/O:
      class(CEAData), intent(inout) :: self
      INTEGER, intent(in)          :: N_atoms_use, N_reac, N_species, N_solids
      INTEGER, intent(in)          :: solid_indices(N_solids)
      real(dp), intent(in) :: solution_vector(N_atoms_use+1+(N_species-self%N_gas))
      real(dp), intent(inout)  :: n ! Moles of gas particles per total mass of mixture in kg
      real(dp), intent(inout)  :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg
      real(dp), intent(in) :: n_spec_old(N_reac) ! Moles of species per total mass of mixture in kg
      real(dp), intent(inout)  :: pi_atom(N_atoms_use) ! Lagrangian multipliers for the atomic species divided
      ! by (R*T)
      LOGICAL, intent(out)          :: converged
      real(dp), intent(in) :: mu_gas(self%N_gas), a_gas(self%N_gas,N_atoms_use), temp

      !! Internal:
      INTEGER                      :: i_reac
      INTEGER, save                :: n_done = 0
      real(dp)             :: lambda, lambda1, lambda2
      real(dp), parameter  :: SIZE = 18.420681
      LOGICAL                      :: gas_good, solids_good, total_good
      real(dp)             :: delta_n_gas(self%N_gas)

      ! IONS
      INTEGER                      :: i_ion, i_stoich
      real(dp)             :: pi_ion, pi_ion_norm, mass

      ! MASS BALANCE CHECKS
      real(dp)             :: b_0(N_atoms_use), b_0_norm, pi_atom_old(N_atoms_use)
      real(dp)             :: mval_mass_good
      INTEGER                      :: i_atom
      LOGICAL                      :: mass_good, pi_good
      real(dp)             :: molfracs_atoms(N_atoms_use)
      real(dp)             :: change

      !f2py integer, intent(aux) :: self%N_gas

      ! print *, "START ec_UPDATE_ABUNDS_SHORT"

      ! Get delta_self%N_gas, following Eq. 2.18:
      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         delta_n_gas(i_reac) = SUM(a_gas(i_reac,1:N_atoms_use)*solution_vector(1:N_atoms_use)) + &
         solution_vector(N_atoms_use+1) - mu_gas(i_reac)/R/temp
      END DO

      ! Calculate correction factors as described in Section 3.3 of the McBride Manual
      lambda1 = 9d99
      lambda2 = 9d99
      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         IF (LOG(n_spec(i_reac)/n) > -SIZE) THEN
            lambda1 = MIN(lambda1,2d0/(MAX(5d0*ABS(solution_vector(N_atoms_use+1)), &
            ABS(delta_n_gas(i_reac)))))
         ELSE IF ((LOG(n_spec(i_reac)/n) <= -SIZE) .AND. (delta_n_gas(i_reac) >= 0d0)) THEN
            lambda2 = MIN(lambda2,ABS((-LOG(n_spec(i_reac)/n)-9.2103404)/ &
            (delta_n_gas(i_reac)-solution_vector(N_atoms_use+1))))
         END IF
      END DO
      lambda = MIN(1d0,lambda1,lambda2)

      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         n_spec(i_reac) = n_spec(i_reac)*exp(lambda*delta_n_gas(i_reac))
      END DO

      IF (self%N_gas < N_species) THEN
         DO i_reac = self%N_gas+1, N_species
            change = lambda*solution_vector(N_atoms_use+1+i_reac-self%N_gas)
            ! if (2*abs(change) < n_spec(solid_indices(i_reac-self%N_gas)) .OR. n_spec(solid_indices(i_reac-self%N_gas)) < tiny(0d0)) then
            !    n_spec(solid_indices(i_reac-self%N_gas)) = n_spec(solid_indices(i_reac-self%N_gas)) + change
            ! else
            !    n_spec(solid_indices(i_reac-self%N_gas)) = (1d0 + sign(0.5d0,change)) * n_spec(solid_indices(i_reac-self%N_gas))
            ! end if
            n_spec(solid_indices(i_reac-self%N_gas)) = n_spec(solid_indices(i_reac-self%N_gas)) + change
         END DO
      END IF
      pi_atom_old = pi_atom
      pi_atom = solution_vector(1:N_atoms_use)
      n = n*exp(lambda*solution_vector(N_atoms_use+1))

      gas_good = .TRUE.
      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         IF (n_spec(i_reac)*ABS(delta_n_gas(i_reac))/SUM(n_spec) > 0.5d-5) THEN
            gas_good = .FALSE.
         END IF
      END DO
      solids_good = .TRUE.
      IF (self%N_gas < N_species) THEN
         DO i_reac = self%N_gas+1, N_species
            IF (ABS(solution_vector(N_atoms_use+1+i_reac-self%N_gas))/SUM(n_spec) > 0.5d-5) THEN
               solids_good = .FALSE.
            END IF
         END DO
      END IF
      total_good = .TRUE.
      IF (n*ABS(solution_vector(N_atoms_use+1))/SUM(n_spec) > 0.5d-5) THEN
         total_good = .FALSE.
      END IF

      !!!-----------------------

      mass_good = .TRUE.
      pi_good = .TRUE.

      if (self%b_0_cache_valid .and. self%b_0_cache_atoms_use == N_atoms_use) then
         b_0 = self%b_0_cache(1:N_atoms_use)
         b_0_norm = self%b_0_norm_cache
      else
         call ec_b_0(self,N_atoms_use, molfracs_atoms, b_0_norm, b_0)
      end if

      mval_mass_good = MAXVAL(b_0)*self%mass_tol
      DO i_atom = 1, N_atoms_use
        mass = 0.0_dp
        do i_reac = 1,self%N_gas
          if (.not. self%remove_ions .or. .not. self%reac_ion(i_reac)) then
             mass = mass + self%a_stoich(i_reac,i_atom)*n_spec(i_reac)
          end if
        enddo
        do i_reac = self%N_gas+1, N_species
          mass = mass + self%a_stoich(solid_indices(i_reac-self%N_gas),i_atom)* &
                 n_spec(solid_indices(i_reac-self%N_gas))
        enddo
        IF ((abs(b_0(i_atom) - mass) > mval_mass_good) .AND. (b_0(i_atom) > 1d-6)) THEN
          mass_good = .FALSE.
        END IF
      END DO

      DO i_atom = 1, N_atoms_use
         IF (abs((pi_atom_old(i_atom)-pi_atom(i_atom))/pi_atom(i_atom)) > 1d-3) THEN
            pi_good = .FALSE.
         END IF
      END DO

      ! IF ((.NOT. mass_good) .OR. (.NOT. pi_good)) THEN
      !    mass_good = .TRUE.
      !    pi_good = .TRUE.
      !    DO i_reac = 1, self%N_reactants
      !       IF (ABS(n_spec(i_reac)-n_spec_old(i_reac)) > 1d-10) THEN
      !          mass_good = .FALSE.
      !          pi_good = .FALSE.
      !       END IF
      !    END DO
      ! END IF

      !!!-------------------

      ! ION CONVERGENCE?

      IF (self%ions .AND. (.NOT. self%remove_ions)) THEN

         ! DO THE MAGIC THEY DO IN SECT. 3.7 in McBride

         pi_ion = 0d0
         pi_ion_norm = 0d0
         DO i_reac = 1, N_species
            DO i_stoich = 1, 5
               ! IF (trim(adjustl(self%reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
               IF (self%reac_atoms_id(i_stoich,i_reac) == 1) THEN
                  pi_ion = pi_ion - n_spec(i_reac)*self%reac_stoich(i_stoich,i_reac)
                  pi_ion_norm = pi_ion_norm + n_spec(i_reac)*self%reac_stoich(i_stoich,i_reac)**2d0
                  EXIT
               END IF
            END DO
         END DO

         pi_ion = pi_ion / pi_ion_norm

         IF (ABS(pi_ion) > 1d-4) THEN
            DO i_ion = 1, 80
               DO i_reac = 1, N_species
                  DO i_stoich = 1, 5
                     ! IF (trim(adjustl(self%reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                     IF (self%reac_atoms_id(i_stoich,i_reac) == 1) THEN
                        n_spec(i_reac) = n_spec(i_reac)*exp(self%reac_stoich(i_stoich,i_reac)*pi_ion)
                        EXIT
                     END IF
                  END DO
               END DO

               pi_ion = 0d0
               pi_ion_norm = 0d0
               DO i_reac = 1, N_species
                  DO i_stoich = 1, 5
                     ! IF (trim(adjustl(self%reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                     IF (self%reac_atoms_id(i_stoich,i_reac) == 1) THEN
                        pi_ion = pi_ion - n_spec(i_reac)*self%reac_stoich(i_stoich,i_reac)
                        pi_ion_norm = pi_ion_norm + n_spec(i_reac)*self%reac_stoich(i_stoich,i_reac)**2d0
                        EXIT
                     END IF
                  END DO
               END DO

            END DO
         END IF

         IF (((gas_good .AND. solids_good) .AND. total_good) .AND. (ABS(pi_ion) <= 1d-4)) THEN
            IF (mass_good .AND. pi_good) THEN
               converged = .TRUE.
            END IF
         END IF

      ELSE

         IF ((gas_good .AND. solids_good) .AND. total_good) THEN
            IF (mass_good .AND. pi_good) THEN
               converged = .TRUE.
            END IF
         END IF

      END IF

      n_done = n_done + 1

   end subroutine ec_UPDATE_ABUNDS_SHORT

   !> Build the big matrix
   subroutine ec_PREP_MATRIX_LONG(self,N_atoms_use, N_reac, molfracs_atoms, N_species,press,temp, &
   H_0, S_0, n, n_spec, matrix, vector, solid_indices, N_solids)
      use equilibrate_const, only: mol, R
      !! I/O:
      class(CEAData), intent(inout) :: self
      INTEGER, intent(in)          :: N_atoms_use, N_reac, N_species, N_solids
      INTEGER, intent(in)          :: solid_indices(N_solids)
      real(dp), intent(in) :: molfracs_atoms(N_atoms_use), press, temp
      real(dp), intent(in) :: H_0(N_reac), S_0(N_reac)
      real(dp), intent(in) :: n ! Moles of gas particles per total mass of mixture in kg
      real(dp), intent(inout)  :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg
      ! real(dp), intent(in) :: pi_atom(N_atoms) ! Lagrangian multipliers for the atomic species divided
      ! ! by (R*T)
      real(dp), intent(out):: matrix(N_species+N_atoms_use+1,N_species+N_atoms_use+1)
      ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
      ! condensed species, the pis and the delta log(n)
      real(dp), intent(out):: vector(N_species+N_atoms_use+1)

      !! Internal:
      real(dp)             :: b_0(N_atoms_use), b_0_norm, b(N_atoms_use)
      real(dp)             :: a(N_species,N_atoms_use), mu(N_species)
      INTEGER                      :: i_atom, i_reac, i_ratom

      ! Set up b0
      ! b_0_norm = 0d0
      ! DO i_atom = 1, N_atoms_use
      !    ! call ec_ATOM_MASS(self%names_atoms(i_atom),mass_atom)
      !    if (self%id_atoms(i_atom) > 0) then
      !       mass_atom = masses_atoms_save(self%id_atoms(i_atom))
      !       b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
      !    end if
      ! END DO
      ! b_0 = molfracs_atoms/b_0_norm
      call ec_b_0(self,N_atoms_use, molfracs_atoms, b_0_norm, b_0)

      ! Set up a_ij
      a = 0d0
      DO i_atom = 1, N_atoms_use
         ! call uppercase(self%names_atoms(i_atom),upper_atom_name)
         DO i_reac = 1, self%N_gas
            IF (self%remove_ions) THEN
               IF (self%reac_ion(i_reac)) THEN
                  a(i_reac,1:N_atoms_use) = 0d0
                  CYCLE
               END IF
            END IF
            DO i_ratom = 1, 5
               ! call uppercase(self%reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (self%reac_atoms_id(i_ratom,i_reac)>0 .and. self%id_atoms(i_atom) == self%reac_atoms_id(i_ratom, i_reac)) then
                  a(i_reac,i_atom) = self%reac_stoich(i_ratom,i_reac)*mol
               END IF
            END DO
         END DO
         DO i_reac = self%N_gas+1, N_species
            DO i_ratom = 1, 5
               ! call uppercase(self%reac_atoms_names(i_ratom,solid_indices(i_reac-self%N_gas)),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (self%reac_atoms_id(i_ratom, solid_indices(i_reac - self%N_gas))>0 .and.&
               self%id_atoms(i_atom) == self%reac_atoms_id(i_ratom, solid_indices(i_reac - self%N_gas))) then
                  a(i_reac,i_atom) = self%reac_stoich(i_ratom,solid_indices(i_reac-self%N_gas))*mol
               END IF
            END DO
         END DO
      END DO

      ! Set up mu_j
      DO i_reac = 1, N_species
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               mu(i_reac) = 0d0
               CYCLE
            END IF
         END IF
         ! Taken from Venot et al. (2012), in comparison with McBride 1996.
         IF (i_reac <= self%N_gas) THEN
            mu(i_reac) = H_0(i_reac) - temp*S_0(i_reac)

            IF (n_spec(i_reac) > 1d-290) THEN
               mu(i_reac) = mu(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*log(press)
            ELSE
               IF (self%verbose) THEN
                  write(*,*) 'n_spec(i_reac) == 0 for '//trim(adjustl(self%names_reactants(i_reac)))// &
                  ' set to 1d-13 and try again.'
               END IF
               n_spec(i_reac) = 1d-13*(1d0 + 1d-6*DBLE(i_reac))
               mu(i_reac) = mu(i_reac) + R*temp*log(n_spec(i_reac)/n)+R*temp*log(press)
            END IF

         ELSE
            mu(i_reac) = H_0(solid_indices(i_reac-self%N_gas)) - temp*S_0(solid_indices(i_reac-self%N_gas))
         END IF
      END DO

      ! MATRIX SETUP
      matrix = 0d0
      ! Set up the matrix for the self%N_gas equations (Eq. 2.18)
      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         matrix(i_reac,i_reac) = 1d0
         DO i_atom = 1, N_atoms_use
            matrix(i_reac,N_species+i_atom) = -a(i_reac,i_atom)
         END DO
         matrix(i_reac,N_species+N_atoms_use+1) = -1d0
      END DO

      ! Set up the matrix for the self%N_reactants-self%N_gas equations (Eq. 2.19)
      IF (self%N_gas < N_species) THEN
         DO i_reac = self%N_gas+1, N_species
            DO i_atom = 1, N_atoms_use
               matrix(i_reac,N_species+i_atom) = -a(i_reac,i_atom)
            END DO
         END DO
      END IF

      ! Set up the matrix for the N_atom equations (Eq. 2.20)
      DO i_atom = 1, N_atoms_use
         DO i_reac = 1, self%N_gas
            IF (self%remove_ions) THEN
               IF (self%reac_ion(i_reac)) THEN
                  CYCLE
               END IF
            END IF
            matrix(N_species+i_atom,i_reac) = a(i_reac,i_atom)*n_spec(i_reac)
         END DO
         IF (self%N_gas < N_species) THEN
            DO i_reac = self%N_gas+1, N_species
               matrix(N_species+i_atom,i_reac) = a(i_reac,i_atom)
            END DO
         END IF
      END DO

      ! Set up the matrix for the last equation (Eq. 2.21)
      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         matrix(N_species+N_atoms_use+1,i_reac) = n_spec(i_reac)
      END DO
      matrix(N_species+N_atoms_use+1,N_species+N_atoms_use+1) = -n

      ! VECTOR SETUP
      !vector(self%N_reactants+N_atoms+1)
      vector = 0d0

      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         vector(i_reac) = -mu(i_reac)/R/temp ! (Eq. 2.18)
      END DO

      IF (self%N_gas < N_species) THEN
         vector(self%N_gas+1:N_species) = -mu(self%N_gas+1:N_species)/R/temp ! (Eq. 2.19)
      END IF

      b = 0d0
      DO i_atom = 1, N_atoms_use
         DO i_reac = 1, self%N_gas
            IF (self%remove_ions) THEN
               IF (self%reac_ion(i_reac)) THEN
                  CYCLE
               END IF
            END IF
            b(i_atom) = b(i_atom) + a(i_reac,i_atom)*n_spec(i_reac)
         END DO
         DO i_reac = self%N_gas+1, N_species
            b(i_atom) = b(i_atom) + a(i_reac,i_atom)*n_spec(solid_indices(i_reac-self%N_gas))
         END DO
      END DO
      vector(N_species+1:N_species+N_atoms_use) = b_0 - b ! (Eq. 2.20)

      vector(N_species+N_atoms_use+1) = n - SUM(n_spec(1:self%N_gas)) ! (Eq. 2.21)

   end subroutine ec_PREP_MATRIX_LONG

   !> Return the result of one step of computation with big matrix(problem: AX=B)
   subroutine ec_UPDATE_ABUNDS_LONG(self,N_atoms_use,N_reac,N_species,solution_vector,n_spec,pi_atom,&
   n,converged,solid_indices,N_solids,molfracs_atoms,n_spec_old)
      use equilibrate_const, only: mol
      !! I/O:
      class(CEAData), intent(inout) :: self
      INTEGER, intent(in)          :: N_atoms_use, N_reac, N_species, N_solids
      INTEGER , intent(in)         :: solid_indices(N_solids)
      real(dp), intent(in) :: solution_vector(N_species+N_atoms_use+1)
      real(dp), intent(inout)  :: n ! Moles of gas particles per total mass of mixture in kg
      real(dp), intent(inout)  :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg
      real(dp), intent(in)     :: n_spec_old(N_reac) ! Moles of species per total mass of mixture in kg
      real(dp), intent(inout)  :: pi_atom(N_atoms_use) ! Lagrangian multipliers for the atomic species divided
      ! by (R*T)
      LOGICAL , intent(out)         :: converged

      !! Internal:
      INTEGER                      :: i_reac
      INTEGER, save                :: n_done = 0
      real(dp)             :: lambda, lambda1, lambda2
      real(dp), parameter  :: SIZE = 18.420681
      LOGICAL                      :: gas_good, solids_good, total_good

      ! IONS
      INTEGER                      :: i_ion, i_stoich
      real(dp)             :: pi_ion, pi_ion_norm, mass

      ! MASS BALANCE CHECKS
      real(dp)             :: b_0(N_atoms_use), b_0_norm, pi_atom_old(N_atoms_use)
      real(dp)             :: a(N_species,N_atoms_use), mval_mass_good
      INTEGER                      :: i_atom, i_ratom
      LOGICAL                      :: mass_good, pi_good
      real(dp)             :: molfracs_atoms(N_atoms_use)
      real(dp)             :: change


      ! Calculate correction factors as described in Section 3.3 of the McBride Manual
      lambda1 = 9d99
      lambda2 = 9d99
      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         IF (LOG(n_spec(i_reac)/n) > -SIZE) THEN
            lambda1 = MIN(lambda1,2d0/(MAX(5d0*ABS(solution_vector(N_species+N_atoms_use+1)), &
            ABS(solution_vector(i_reac)))))
         ELSE IF ((LOG(n_spec(i_reac)/n) <= -SIZE) .AND. (solution_vector(i_reac) >= 0d0)) THEN
            lambda2 = MIN(lambda2,ABS((-LOG(n_spec(i_reac)/n)-9.2103404)/ &
            (solution_vector(i_reac)-solution_vector(N_species+N_atoms_use+1))))
         END IF
      END DO
      lambda = MIN(1d0,lambda1,lambda2)

      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         n_spec(i_reac) = n_spec(i_reac)*exp(lambda*solution_vector(i_reac))
      END DO

      IF (self%N_gas < N_species) THEN
         DO i_reac = self%N_gas+1, N_species
            change = lambda*solution_vector(N_atoms_use+1+i_reac-self%N_gas)
            ! if (2*abs(change) < n_spec(solid_indices(i_reac-self%N_gas)) .OR. n_spec(solid_indices(i_reac-self%N_gas)) < tiny(0d0)) then
            !    n_spec(solid_indices(i_reac-self%N_gas)) = n_spec(solid_indices(i_reac-self%N_gas)) + changepi_atom
            ! else
            !    n_spec(solid_indices(i_reac-self%N_gas)) = (1d0 + sign(0.5d0,change)) * n_spec(solid_indices(i_reac-self%N_gas))
            ! end if
            n_spec(solid_indices(i_reac-self%N_gas)) = n_spec(solid_indices(i_reac-self%N_gas)) + change
         END DO
      END IF
      pi_atom_old = pi_atom
      pi_atom = solution_vector(N_species+1:N_species+N_atoms_use)
      n = n*exp(lambda*solution_vector(N_species+N_atoms_use+1))

      gas_good = .TRUE.
      DO i_reac = 1, self%N_gas
         IF (self%remove_ions) THEN
            IF (self%reac_ion(i_reac)) THEN
               CYCLE
            END IF
         END IF
         IF (n_spec(i_reac)*ABS(solution_vector(i_reac))/SUM(n_spec) > 0.5d-5) THEN
            gas_good = .FALSE.
         END IF
      END DO
      solids_good = .TRUE.
      IF (self%N_gas < N_species) THEN
         DO i_reac = self%N_gas+1, N_species
            IF (ABS(solution_vector(i_reac))/SUM(n_spec) > 0.5d-5) THEN
               solids_good = .FALSE.
            END IF
         END DO
      END IF
      total_good = .TRUE.
      IF (n*ABS(solution_vector(N_species+N_atoms_use+1))/SUM(n_spec) > 0.5d-5) THEN
         total_good = .FALSE.
      END IF

      !!!-----------------------

      mass_good = .TRUE.
      pi_good = .TRUE.


      ! Set up b0
      ! b_0_norm = 0d0
      ! DO i_atom = 1, N_atoms_use
      !    ! call ec_ATOM_MASS(self%names_atoms(i_atom),mass_atom)
      !    mass_atom = masses_atoms_save(self%id_atoms(i_atom))
      !    b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
      ! END DO
      ! b_0 = molfracs_atoms/b_0_norm
      call ec_b_0(self,N_atoms_use, molfracs_atoms, b_0_norm, b_0)

      ! Set up a_ij
      a = 0d0
      DO i_atom = 1, N_atoms_use
         ! call uppercase(self%names_atoms(i_atom),upper_atom_name)
         DO i_reac = 1, self%N_gas
            IF (self%remove_ions) THEN
               IF (self%reac_ion(i_reac)) THEN
                  a(i_reac,1:N_atoms_use) = 0d0
                  CYCLE
               END IF
            END IF
            DO i_ratom = 1, 5
               ! call uppercase(self%reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (self%reac_atoms_id(i_ratom, i_reac)>0 .and. self%id_atoms(i_atom) == self%reac_atoms_id(i_ratom, i_reac)) then
                  a(i_reac,i_atom) = self%reac_stoich(i_ratom,i_reac)*mol
               END IF
            END DO
         END DO
         DO i_reac = self%N_gas+1, N_species
            DO i_ratom = 1, 5
               ! call uppercase(self%reac_atoms_names(i_ratom,solid_indices(i_reac-self%N_gas)),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (self%reac_atoms_id(i_ratom, solid_indices(i_reac - self%N_gas))>0 .and. &
               self%id_atoms(i_atom) == self%reac_atoms_id(i_ratom, solid_indices(i_reac - self%N_gas))) then
                  a(i_reac,i_atom) = self%reac_stoich(i_ratom,solid_indices(i_reac-self%N_gas))*mol
               END IF
            END DO
         END DO
      END DO

      mval_mass_good = MAXVAL(b_0)*self%mass_tol
      DO i_atom = 1, N_atoms_use
        mass = 0.0_dp
        do i_reac = 1,self%N_gas
          mass = mass + a(i_reac,i_atom)*n_spec(i_reac)
        enddo
        do i_reac = self%N_gas+1, N_species
          mass = mass + a(i_reac,i_atom)*n_spec(solid_indices(i_reac-self%N_gas))
        enddo
        IF ((abs(b_0(i_atom) - mass) > mval_mass_good) .AND. (b_0(i_atom) > 1d-6)) THEN
          mass_good = .FALSE.
        END IF
      END DO

      DO i_atom = 1, N_atoms_use
         IF (abs((pi_atom_old(i_atom)-pi_atom(i_atom))/pi_atom(i_atom)) > 1d-3) THEN
            pi_good = .FALSE.
         END IF
      END DO

      ! IF ((.NOT. mass_good) .OR. (.NOT. pi_good)) THEN
      !    mass_good = .TRUE.
      !    pi_good = .TRUE.
      !    DO i_reac = 1, self%N_reactants
      !       IF (ABS(n_spec(i_reac)-n_spec_old(i_reac)) > 1d-10) THEN
      !          mass_good = .FALSE.
      !          pi_good = .FALSE.
      !       END IF
      !    END DO
      ! END IF

      !!!-------------------

      ! ION CONVERGENCE?

      IF (self%ions .AND. (.NOT. self%remove_ions)) THEN

         ! DO THE MAGIC THEY DO IN SECT. 3.7 in McBride

         pi_ion = 0d0
         pi_ion_norm = 0d0
         DO i_reac = 1, N_species
            DO i_stoich = 1, 5
               ! IF (trim(adjustl(self%reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
               IF (self%reac_atoms_id(i_stoich,i_reac) == 1) THEN
                  pi_ion = pi_ion - n_spec(i_reac)*self%reac_stoich(i_stoich,i_reac)
                  pi_ion_norm = pi_ion_norm + n_spec(i_reac)*self%reac_stoich(i_stoich,i_reac)**2d0
                  EXIT
               END IF
            END DO
         END DO

         pi_ion = pi_ion / pi_ion_norm

         IF (ABS(pi_ion) > 1d-4) THEN
            DO i_ion = 1, 80
               DO i_reac = 1, N_species
                  DO i_stoich = 1, 5
                     ! IF (trim(adjustl(self%reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                     IF (self%reac_atoms_id(i_stoich,i_reac) == 1) THEN
                        n_spec(i_reac) = n_spec(i_reac)*exp(self%reac_stoich(i_stoich,i_reac)*pi_ion)
                        EXIT
                     END IF
                  END DO
               END DO

               pi_ion = 0d0
               pi_ion_norm = 0d0
               DO i_reac = 1, N_species
                  DO i_stoich = 1, 5
                     ! IF (trim(adjustl(self%reac_atoms_names(i_stoich,i_reac))) .EQ. 'E') THEN
                     IF (self%reac_atoms_id(i_stoich,i_reac) == 1) THEN
                        pi_ion = pi_ion - n_spec(i_reac)*self%reac_stoich(i_stoich,i_reac)
                        pi_ion_norm = pi_ion_norm + n_spec(i_reac)*self%reac_stoich(i_stoich,i_reac)**2d0
                        EXIT
                     END IF
                  END DO
               END DO
            END DO
         END IF

         IF (((gas_good .AND. solids_good) .AND. total_good) .AND. (ABS(pi_ion) <= 1d-4)) THEN
            IF (mass_good .AND. pi_good) THEN
               converged = .TRUE.
            END IF
         END IF

      ELSE

         IF ((gas_good .AND. solids_good) .AND. total_good) THEN
            IF (mass_good .AND. pi_good) THEN
               converged = .TRUE.
            END IF
         END IF

      END IF

      n_done = n_done + 1

   end subroutine ec_UPDATE_ABUNDS_LONG

   !> Computes the adiabatic gradient
   subroutine ec_COMP_ADIABATIC_GRAD(self,N_atoms_use,N_reac,N_spec_eff,n_spec, &
   n,H_0,C_P_0,solid_indices,N_solids,temp,nabla_ad,gamma2,c_pe)
      use equilibrate_const, only: mol, R

      !! I/O:
      class(CEAData), intent(inout) :: self
      INTEGER, intent(in)          :: N_atoms_use, N_reac, N_spec_eff, N_solids
      INTEGER, intent(in)          :: solid_indices(N_solids)
      real(dp), intent(in) :: temp
      real(dp), intent(in) :: C_P_0(N_reac), H_0(N_reac)
      real(dp), intent(in) :: n ! Moles of gas particles per total mass of mixture in kg
      real(dp), intent(in) :: n_spec(N_reac) ! Moles of species per total mass of mixture in kg
      real(dp), intent(out):: nabla_ad, gamma2, c_pe

      !! Internal:
      real(dp)             :: matrix(N_atoms_use+1+(N_spec_eff-self%N_gas),N_atoms_use+1+(N_spec_eff-self%N_gas))
      ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
      ! condensed species, the pis and the delta log(n)
      real(dp)             :: vector(N_atoms_use+1+(N_spec_eff-self%N_gas)), &
      solution_vector(N_atoms_use+1+(N_spec_eff-self%N_gas))
      real(dp)             :: a(N_spec_eff,N_atoms_use)
      INTEGER                      :: i_atom, i_reac, i_ratom, i_atom2

      ! Set up a_ij
      a = 0d0
      DO i_atom = 1, N_atoms_use
         ! call uppercase(self%names_atoms(i_atom),upper_atom_name)
         DO i_reac = 1, self%N_gas
            IF (self%remove_ions) THEN
               IF (self%reac_ion(i_reac)) THEN
                  a(i_reac,i_atom) = 0d0
                  CYCLE
               END IF
            END IF
            DO i_ratom = 1, 5
               ! call uppercase(self%reac_atoms_names(i_ratom,i_reac),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (self%reac_atoms_id(i_ratom, i_reac)>0 .and. self%id_atoms(i_atom) == self%reac_atoms_id(i_ratom, i_reac)) then
                  a(i_reac,i_atom) = self%reac_stoich(i_ratom,i_reac)*mol
               END IF
            END DO
         END DO
         DO i_reac = self%N_gas+1, N_spec_eff
            DO i_ratom = 1, 5
               ! call uppercase(self%reac_atoms_names(i_ratom,solid_indices(i_reac-self%N_gas)),upper_ratom_name)
               ! IF (trim(adjustl(upper_atom_name)) .EQ. trim(adjustl(upper_ratom_name))) THEN
               IF (self%reac_atoms_id(i_ratom, solid_indices(i_reac - self%N_gas))>0 .and. &
               self%id_atoms(i_atom) == self%reac_atoms_id(i_ratom, solid_indices(i_reac - self%N_gas))) then
                  a(i_reac,i_atom) = self%reac_stoich(i_ratom,solid_indices(i_reac-self%N_gas))*mol
               END IF
            END DO
         END DO
      END DO

      matrix = 0d0
      ! Setup matrix, following Eq. 2.56
      DO i_atom = 1, N_atoms_use
         ! First term, LHS
         DO i_atom2 = 1, N_atoms_use
            DO i_reac = 1, self%N_gas
               IF (self%remove_ions) THEN
                  IF (self%reac_ion(i_reac)) THEN
                     CYCLE
                  END IF
               END IF
               matrix(i_atom,i_atom2) = matrix(i_atom,i_atom2) + &
               n_spec(i_reac)*a(i_reac,i_atom2)*a(i_reac,i_atom)
            END DO
         END DO
         ! Second term, LHS
         DO i_reac = self%N_gas+1, N_spec_eff
            matrix(i_atom,N_atoms_use+1+i_reac-self%N_gas) = a(i_reac,i_atom)
         END DO
         ! Third term, LHS
         DO i_reac = 1, self%N_gas
            IF (self%remove_ions) THEN
               IF (self%reac_ion(i_reac)) THEN
                  CYCLE
               END IF
            END IF
            matrix(i_atom,N_atoms_use+1) = matrix(i_atom,N_atoms_use+1) + &
            a(i_reac,i_atom)*n_spec(i_reac)
         END DO
      END DO

      ! Setup matrix, following Eq. 2.58
      DO i_atom = 1, N_atoms_use
         DO i_reac = 1, self%N_gas
            IF (self%remove_ions) THEN
               IF (self%reac_ion(i_reac)) THEN
                  CYCLE
               END IF
            END IF
            matrix(N_atoms_use+1,i_atom) = matrix(N_atoms_use+1,i_atom) + &
            a(i_reac,i_atom)*n_spec(i_reac)
         END DO
      END DO

      ! Setup matrix, following Eq. 2.57
      DO i_reac = self%N_gas+1, N_spec_eff
         DO i_atom = 1, N_atoms_use
            matrix(N_atoms_use+1+i_reac-self%N_gas,i_atom) = a(i_reac,i_atom)
         END DO
      END DO

      vector = 0d0
      ! Setup the vector, following Eq. 2.56
      DO i_atom = 1, N_atoms_use
         vector(i_atom) = -SUM(a(1:self%N_gas,i_atom)*n_spec(1:self%N_gas)*H_0(1:self%N_gas)) &
         /R/temp
      END DO

      ! Setup the vector, following Eq. 2.58
      vector(N_atoms_use+1) = -SUM(n_spec(1:self%N_gas)*H_0(1:self%N_gas))/R/temp

      ! Setup the vector, following Eq. 2.57
      DO i_reac = self%N_gas+1, N_spec_eff
         vector(N_atoms_use+1+i_reac-self%N_gas) = -H_0(solid_indices(i_reac-self%N_gas))/R/temp
      END DO

      ! Solve the system
      call ec_INVERT_MATRIX_SHORT(self,N_atoms_use+1+N_spec_eff-self%N_gas,matrix,vector,solution_vector)
      if (self%error) RETURN

      ! Calculate c_pe, following Eq. 2.59
      c_pe = 0d0
      DO i_atom = 1, N_atoms_use
         c_pe = c_pe + SUM(a(1:self%N_gas,i_atom)*n_spec(1:self%N_gas)*H_0(1:self%N_gas)/R/temp) * &
         solution_vector(i_atom)
      END DO
      DO i_reac = self%N_gas+1, N_spec_eff
         c_pe = c_pe + H_0(solid_indices(i_reac-self%N_gas))/R/temp* &
         solution_vector(N_atoms_use+1+i_reac-self%N_gas) + &
         n_spec(solid_indices(i_reac-self%N_gas))* &
         C_P_0(solid_indices(i_reac-self%N_gas))/R
      END DO
      c_pe = c_pe + SUM(n_spec(1:self%N_gas)*C_P_0(1:self%N_gas)/R)
      c_pe = c_pe + SUM(n_spec(1:self%N_gas)*H_0(1:self%N_gas)/R/temp)* &
      solution_vector(N_atoms_use+1)
      c_pe = c_pe + SUM(n_spec(1:self%N_gas)*(H_0(1:self%N_gas)/R/temp)**2d0)
      c_pe = c_pe*R

      ! Calculate nabla_ad, using Eq. 2.50 and Eq. 2.75
      nabla_ad = c_pe/n/R/(1d0+solution_vector(N_atoms_use+1))
      nabla_ad = 1/nabla_ad
      gamma2 = 1d0/(1d0-nabla_ad)

   end subroutine ec_COMP_ADIABATIC_GRAD

   subroutine ec_b_0(self,N_atoms_use, molfracs_atoms, b_0_norm, b_0)
      use equilibrate_const, only: masses_atoms_save
      class(CEAData), intent(inout) :: self
      integer, intent(in)           :: N_atoms_use
      real(dp), intent(in)  :: molfracs_atoms(N_atoms_use)
      real(dp), intent(out) :: b_0_norm, b_0(N_atoms_use)

      integer                       :: i_atom
      real(dp)              :: mass_atom

      ! Set up b0
      b_0_norm = 0d0
      DO i_atom = 1, N_atoms_use
         ! call ec_ATOM_MASS(self%names_atoms(i_atom),mass_atom)
         mass_atom = masses_atoms_save(self%id_atoms(i_atom))
         b_0_norm = b_0_norm + mass_atom*molfracs_atoms(i_atom)
      END DO
      b_0 = molfracs_atoms/b_0_norm
   end subroutine ec_b_0

   !> Invert the small matrix
   subroutine ec_INVERT_MATRIX_SHORT(self, lens,matrix,vector,solution_vector)
      ! use data_block, only: error
      ! implicit none
      !! I/O:
      class(CEAData), intent(inout) :: self
      INTEGER, intent(in)           :: lens
      real(dp), intent(inout)  :: matrix(lens,lens)
      ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
      ! condensed species, the pis and the delta log(n)
      real(dp), intent(in)  :: vector(lens)
      real(dp), intent(out) :: solution_vector(lens)
      !! Internal:
      INTEGER                      :: index(lens)
      real(dp)             :: indic

      solution_vector = vector

      call ec_LUDCMP(self,matrix,index,indic)
      if (self%error) RETURN
      call ec_LUBKSB(self,matrix,index,solution_vector)

   end subroutine ec_INVERT_MATRIX_SHORT

   !> Invert the big matrix
   subroutine ec_INVERT_MATRIX_LONG(self,lens,matrix,vector,solution_vector)
      ! use data_block, only: N_gas, N_ions, remove_ions, reac_ion, error
      ! implicit none
      !! I/O:
      class(CEAData), intent(inout) :: self
      INTEGER, intent(in)           :: lens
      real(dp), intent(inout)  :: matrix(lens,lens)
      real(dp), intent(in)  :: vector(lens)
      real(dp), intent(out) :: solution_vector(lens)
      ! So the solution vector will contain the delta log(n_j) for gas, the delta n_j for
      ! condensed species, the pis and the delta log(n)

      real(dp)              :: matrix_nions(lens-self%N_ions,lens-self%N_ions)
      real(dp)              :: vector_nions(lens-self%N_ions), &
      solution_vector_nions(lens-self%N_ions)

      !! Internal:
      INTEGER                       :: index(lens), corrf_i, corrf_j, index_nions(lens-self%N_ions)
      real(dp)              :: indic
      INTEGER                       :: i_mat, j_mat

      solution_vector = vector

      !matrix(self%N_reactants+N_atoms+1,self%N_reactants+N_atoms+1)
      IF (self%remove_ions) THEN
         vector_nions = 0d0
         matrix_nions = 0d0
         corrf_i = 0
         DO i_mat = 1, lens
            corrf_j = 0
            IF (i_mat <= self%N_gas) THEN
               IF (self%reac_ion(i_mat)) THEN
                  corrf_i = corrf_i + 1
                  cycle
               END IF
            END IF
            DO j_mat = 1, lens
               IF (j_mat <= self%N_gas) THEN
                  IF (self%reac_ion(j_mat)) THEN
                     corrf_j = corrf_j + 1
                     cycle
                  END IF
               END IF
               matrix_nions(j_mat-corrf_j,i_mat-corrf_i) = matrix(j_mat,i_mat)
            END DO
            vector_nions(i_mat-corrf_i) = vector(i_mat)
         END DO
         solution_vector_nions = vector_nions

         call ec_LUDCMP(self,matrix_nions,index_nions,indic)
         if (self%error) RETURN
         call ec_LUBKSB(self,matrix_nions,index_nions,solution_vector_nions)
         if (self%error) RETURN

         corrf_i = 0
         DO i_mat = 1, lens
            IF (i_mat <= self%N_gas) THEN
               IF (self%reac_ion(i_mat)) THEN
                  corrf_i = corrf_i + 1
                  cycle
               END IF
            END IF
            solution_vector(i_mat) = solution_vector_nions(i_mat-corrf_i)
         END DO
      ELSE
         call ec_LUDCMP(self,matrix,index,indic)
         if (self%error) RETURN
         call ec_LUBKSB(self,matrix,index,solution_vector)
         if (self%error) RETURN
      END IF

   end subroutine ec_INVERT_MATRIX_LONG

   !> LU decomposition, from numerical recipes
   subroutine ec_LUDCMP(self,a,indx,d)
      ! use data_block, only: error, err_msg
      ! implicit none
      class(CEAData), intent(inout) :: self
      real(dp), intent(inout)  :: a(:,:)
      integer, intent(out)             :: indx(:)
      real(dp), intent(out)    :: d

      real(dp)                 :: vv(size(a,1))
      real(dp), parameter      :: TINY = 1.0e-20
      integer                          :: j, n, imax

      n = lu_asserteq(self,size(a,1),size(a,2),size(indx),'ec_LUDCMP')
      if (self%error) RETURN

      d = 1.0
      vv = maxval(abs(a),dim=2)
      if (any(vv == 0.0)) then
         self%error = .true.
         self%err_msg = "ERROR: singular matrix in LUDCMP"
         RETURN
      end if
      vv = 1.0/vv
      do j = 1,n
         imax = (j-1)+maxloc(vv(j:n)*abs(a(j:n,j)), dim=1)
         if (j /= imax) then
            call lu_swap(self,n,a(imax,:),a(j,:))
            d = -d
            vv(imax) = vv(j)
         end if
         indx(j) = imax
         if (a(j,j) == 0.0) a(j,j) = TINY
         a(j+1:n,j) = a(j+1:n,j) / a(j,j)
         ! a(j+1:n,j+1:n) = a(j+1:n,j+1:n)-lu_outerprod(n-j,a(j+1:n,j),n-j,a(j,j+1:n))
         a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - &
         (spread(a(j+1:n,j), dim=2, ncopies=n-j) * spread(a(j,j+1:n), dim=1, ncopies=n-j))
      end do
   END SUBROUTINE ec_LUDCMP

   !> LU back substitution
   SUBROUTINE ec_LUBKSB(self,a,indx,b)
      ! use data_block, only: error
      ! implicit none
      class(CEAData), intent(inout) :: self
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: indx(:)
      real(dp), intent(inout) :: b(:)
      integer :: i,n,ii,ll
      real(dp) :: summ

      n=lu_asserteq(self,size(a,1),size(a,2),size(indx),'ec_LUBKSB')
      if (self%error) RETURN

      ii=0
      do i=1,n
         ll=indx(i)
         summ=b(ll)
         b(ll)=b(i)
         if (ii /= 0) then
            summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
         else if (summ /= 0.0) then
            ii=i
         end if
         b(i)=summ
      end do
      do i=n,1,-1
         b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
      end do
   END SUBROUTINE ec_LUBKSB

   subroutine lu_swap(self,len,arr1,arr2)
      class(CEAData), intent(inout) :: self
      integer, intent(in)              :: len
      real(dp), intent(inout)  :: arr1(len), arr2(len)
      real(dp)                 :: tamp(len)

      tamp = arr1
      arr1 = arr2
      arr2 = tamp
   end subroutine lu_swap

   function lu_asserteq(self,n1,n2,n3,label) result(m)
      ! use data_block, only: error, err_msg
      ! implicit none
      class(CEAData), intent(inout) :: self
      integer, intent(in)        :: n1, n2, n3
      character(len=9), intent(in)   :: label
      integer                    :: m

      if (n1==n2 .and. n2==n3) then
         m = n1
      else
         self%error = .true.
         self%err_msg = 'ERROR in '//trim(adjustl(label))//": the sizes of the matrix' don't add up..."
         ! print *, 'nrerror: assert_eq failed in: ', label
         ! STOP 'Program terminated by lu_asserteq'
         RETURN
      end if
   end function lu_asserteq

end module
