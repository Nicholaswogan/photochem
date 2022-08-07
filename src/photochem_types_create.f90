submodule (photochem_types) photochem_types_create
  use yaml_types, only : type_node, type_dictionary, type_list, type_error, &
                          type_list_item, type_scalar, type_key_value_pair
  implicit none

contains

  module function create_PhotoSettings(filename, err) result(s)
    use fortran_yaml_c, only : parse, error_length
    use photochem_types, only: PhotoSettings
    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: err
    
    type(PhotoSettings) :: s
    
    character(error_length) :: error
    class (type_node), pointer :: root
    
    ! parse yaml file
    root => parse(filename,error=error)
    if (error /= '') then
      err = trim(error)
      return
    end if
    select type (root)
      class is (type_dictionary)
        s = unpack_PhotoSettings(root, filename, err)
      class default
        err = "yaml file must have dictionaries at root level"
    end select
    call root%finalize()
    deallocate(root)  
    if (allocated(err)) return
      
  end function
  
  function unpack_PhotoSettings(root, filename, err) result(s)
    use photochem_enum, only: VelocityBC, MosesBC
    use photochem_enum, only: DiffusionLimHydrogenEscape, ZahnleHydrogenEscape, NoHydrogenEscape
    use photochem_types, only: PhotoSettings
    type(type_dictionary), intent(in) :: root
    character(*), intent(in) :: filename
    character(:), allocatable, intent(out) :: err
    
    type(PhotoSettings) :: s
    
    type(type_dictionary), pointer :: dict, tmp2, tmp3
    type(type_list), pointer :: list, bcs
    type(type_list_item), pointer :: item
    type(type_scalar), pointer :: scalar
    type(type_error), allocatable :: io_err
    character(:), allocatable :: temp_char
    logical :: success
    integer :: i, j

    ! filename
    s%filename = filename
    
    !!!!!!!!!!!!!!!!!!!!!!!
    !!! atmosphere-grid !!!
    !!!!!!!!!!!!!!!!!!!!!!!
    dict => root%get_dictionary('atmosphere-grid',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    s%bottom = dict%get_real('bottom',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    s%top = dict%get_real('top',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    s%nz = dict%get_integer('number-of-layers',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
  
    !!!!!!!!!!!!!!!!!!!!!!!
    !!! photolysis-grid !!!
    !!!!!!!!!!!!!!!!!!!!!!!
    dict => root%get_dictionary('photolysis-grid',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    s%regular_grid = dict%get_logical('regular-grid',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    if (s%regular_grid) then
      s%lower_wv = dict%get_real('lower-wavelength',error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      s%upper_wv = dict%get_real('upper-wavelength',error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      s%nw = dict%get_integer('number-of-bins',error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      if (s%nw < 1) then 
        err = 'Number of photolysis bins must be >= 1 in '//trim(filename)
        return
      endif
      if (s%lower_wv > s%upper_wv) then
        err = 'lower-wavelength must be smaller than upper-wavelength in '//trim(filename)
        return
      endif
    else
      s%grid_file = dict%get_string('input-file',error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    endif
    ! scale factor for photon flux. Its optional
    s%photon_scale_factor = dict%get_real('photon-scale-factor', 1.0_dp,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    !!!!!!!!!!!!!!
    !!! planet !!!
    !!!!!!!!!!!!!!
    dict => root%get_dictionary('planet',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    scalar => dict%get_scalar('background-gas',required=.false.,error = io_err)
    if (associated(scalar)) then
      s%back_gas_name = trim(scalar%string)
    endif

    scalar => dict%get_scalar('surface-pressure',required=.false.,error = io_err)
    if (associated(scalar)) then
      s%P_surf = scalar%to_real(-1.0_dp, success=success)
      if (.not. success) then
        err = trim(filename)//trim(scalar%path)//' can not be converted to a real number.'
        return
      endif
      if (s%P_surf <= 0.0_dp) then
        err = 'IOError: Planet surface pressure must be greater than zero.'
        return
      endif
    endif

    s%planet_mass = dict%get_real('planet-mass',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%planet_mass <= 0.0_dp) then
      err = 'IOError: Planet mass must be greater than zero.'
      return
    endif
    s%planet_radius = dict%get_real('planet-radius',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%planet_radius <= 0.0_dp) then
      err = 'IOError: Planet radius must be greater than zero.'
      return
    endif
    s%surface_albedo = dict%get_real('surface-albedo',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%surface_albedo < 0.0_dp) then
      err = 'IOError: Surface albedo must be greater than zero.'
      return
    endif
    s%diurnal_fac = dict%get_real('diurnal-averaging-factor',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%diurnal_fac < 0.0_dp .or. s%diurnal_fac > 1.0_dp) then
      err = 'IOError: diurnal-averaging-factor must be between 0 and 1.'
      return
    endif
    
    s%solar_zenith = dict%get_real('solar-zenith-angle',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%solar_zenith < 0.0_dp .or. s%solar_zenith > 90.0_dp) then
      err = 'IOError: solar zenith must be between 0 and 90.'
      return
    endif
    
    ! H2 escape
    tmp2 => dict%get_dictionary('hydrogen-escape',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    temp_char = tmp2%get_string('type',error=io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    temp_char = trim(temp_char)
    if (temp_char == 'diffusion limited') then
      s%H_escape_type = DiffusionLimHydrogenEscape
    elseif (temp_char == 'zahnle') then
      s%H_escape_type = ZahnleHydrogenEscape
      allocate(s%H_escape_S1)
      s%H_escape_S1 = tmp2%get_real('S1',error=io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    elseif (temp_char == 'none') then
      s%H_escape_type = NoHydrogenEscape
    else
      err = '"'//temp_char//'" is not an a valid hydrogen escape type in '//trim(filename)
      return
    endif
    
    ! default lower boundary
    temp_char = trim(dict%get_string('default-gas-lower-boundary',"deposition velocity",error = io_err))
    if (trim(temp_char) == 'deposition velocity') then
      s%default_lowerboundcond = VelocityBC
    elseif (trim(temp_char) == 'Moses') then
      s%default_lowerboundcond = MosesBC
    else
      err = "IOError: Only 'deposition velocity' or 'Moses' can be default boundary conditions."
      return
    endif

    ! climate
    s%evolve_climate = dict%get_logical('evolve-climate',default=.false.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    tmp2 => dict%get_dictionary('water',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    s%fix_water_in_trop = tmp2%get_logical('fix-water-in-troposphere',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%evolve_climate .and. .not. s%fix_water_in_trop) then
      err = 'fix-water-in-troposphere must be true if evolve-climate is true in '//filename
      return
    endif
    s%water_cond = tmp2%get_logical('water-condensation',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    s%gas_rainout = tmp2%get_logical('gas-rainout',error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    if (s%fix_water_in_trop) then  
    
      s%relative_humidity = trim(tmp2%get_string('relative-humidity',error = io_err))
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif

      if (s%evolve_climate .and. s%relative_humidity == 'manabe') then
        err = 'relative-humidity must not be "manabe" if evolve-climate is true in '//filename
        return
      endif
      
    endif
    
    if (s%gas_rainout) then
      s%rainfall_rate = tmp2%get_real('rainfall-rate',error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    endif

    if (s%gas_rainout) then
      nullify(list)
      list => tmp2%get_list('rainout-species', required = .false., error = io_err)
      if (associated(list)) then
        call unpack_string_list(list, s%rainout_species, err)
        if (allocated(err)) return
        i = check_for_duplicates(s%rainout_species)
        if (i /= 0) then
          err = '"'//trim(s%rainout_species(i))//'" is a duplicate in '//trim(list%path)
          return
        endif
      endif
    endif
    
    if (s%fix_water_in_trop .or. s%gas_rainout) then
      ! we need a tropopause altitude
      s%trop_alt = tmp2%get_real('tropopause-altitude',error = io_err)
      if (allocated(io_err)) then
        err = "tropopause-altitude must be specified if fix-water-in-troposphere = true, or gas-rainout = true"
        return
      endif
      if ((s%trop_alt < s%bottom) .or. &
          (s%trop_alt > s%top)) then
          err = 'IOError: tropopause-altitude must be between the top and bottom of the atmosphere'
          return
      endif
    
    endif
    
    if (s%water_cond) then        
    
      tmp3 => tmp2%get_dictionary('condensation-rate',.true.,error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
      s%H2O_condensation_rate(1) = tmp3%get_real('A',error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      s%H2O_condensation_rate(2) = tmp3%get_real('rhc',error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      s%H2O_condensation_rate(3) = tmp3%get_real('rh0',error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      if (s%H2O_condensation_rate(3) <= s%H2O_condensation_rate(2)) then
        err = 'IOError: Rate constant "rh0" for H2O condensation must be > "rhc". See '//trim(filename)
        return
      endif
      
    endif
    
    !!!!!!!!!!!!!!!!!
    !!! particles !!!
    !!!!!!!!!!!!!!!!!
    
    list => root%get_list('particles', .false., error = io_err)
      
    if (associated(list)) then 
      allocate(s%con_names(list%size()))
      allocate(s%con(size(s%con_names)))
      
      j = 1
      item => list%first
      do while (associated(item))
        select type (e => item%node)
        class is (type_dictionary)
          s%con_names(j) = trim(e%get_string('name',error = io_err))
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
          tmp2 => e%get_dictionary('condensation-rate',.true.,error = io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
          s%con(j)%A = tmp2%get_real('A',error = io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
          s%con(j)%rhc = tmp2%get_real('rhc',error = io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
          s%con(j)%rh0 = tmp2%get_real('rh0',error = io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
          if (s%con(j)%rh0 <= s%con(j)%rhc) then
            err = 'IOError: Rate constant "rh0" for '//trim(s%con_names(j)) &
                  //' condensation must be > "rhc". See '//trim(filename)
            return
          endif
      
        class default
          err = "IOError: Particle settings must be a list of dictionaries."
          return
        end select
        j = j + 1
        item => item%next
      end do

    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! boundary-conditions !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    bcs => root%get_list('boundary-conditions',.true.,error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
   
    ! allocate boundary conditions
    allocate(s%ubcs(bcs%size()))
    allocate(s%lbcs(size(s%ubcs)))
    allocate(s%sp_types(size(s%ubcs)))
    allocate(s%sp_names(size(s%ubcs)))
    allocate(s%only_eddy(size(s%ubcs)))
   
    ! default boundary conditions
    do j = 1,size(s%ubcs)
     s%lbcs(j)%bc_type = s%default_lowerboundcond
     s%lbcs(j)%vel = 0.0_dp
     s%ubcs(j)%bc_type = VelocityBC
     s%ubcs(j)%vel = 0.0_dp
     s%only_eddy(j) = .false.
    enddo
     
    s%nsl = 0
    j = 1
    item => bcs%first
    do while (associated(item))
      select type (e => item%node)
      class is (type_dictionary)
        s%sp_names(j) = trim(e%get_string('name',error = io_err))
        if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        ! check for duplicates
        if (j > 1) then
          do i = 1,j-1
            if (s%sp_names(j) == s%sp_names(i)) then
              err = "IOError: Species "//trim(s%sp_names(i))// &
              " has more than one boundary conditions entry in the settings file."
              return
            endif
          enddo
        endif
        
        s%sp_types(j) = trim(e%get_string('type','long lived',error = io_err))
        if (s%sp_types(j) == 'short lived') then
          s%nsl = s%nsl + 1
        
        elseif (s%sp_types(j) == 'long lived') then
          ! get boundary condition
          dict => e%get_dictionary("upper-boundary",.true.,error = io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
          call unpack_SettingsBC(dict, "upper", s%sp_names(j), filename, s%ubcs(j), err)
          if (allocated(err)) return
          
          dict => e%get_dictionary("lower-boundary",.true.,error = io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
          call unpack_SettingsBC(dict, "lower", s%sp_names(j), filename, s%lbcs(j), err)
          if (allocated(err)) return
        
          s%only_eddy(j) = e%get_logical("only-eddy",default=.false., error = io_err)
          if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
        else
          err = 'IOError: species type '//s%sp_types(j)//' is not a valid.' 
          return
        endif
      class default
        err = "IOError: Boundary conditions must be a list of dictionaries."
        return
      end select 
      j = j + 1
      item => item%next
    enddo
    
    allocate(s%SL_names(s%nsl))
    i = 1
    do j = 1,size(s%sp_names)
      if (s%sp_types(j) == "short lived") then
        s%SL_names(i) = s%sp_names(j)
        i = i + 1
      endif
    enddo   
    
  end function
  
  subroutine unpack_SettingsBC(bc, bc_kind, sp_name, filename, sbc, err)
    use photochem_types, only: SettingsBC
    use photochem_enum, only: MosesBC, VelocityBC, MixingRatioBC, FluxBC
    use photochem_enum, only: VelocityDistributedFluxBC, DensityBC
    type(type_dictionary), intent(in) :: bc
    character(*), intent(in) :: bc_kind
    character(*), intent(in) :: sp_name
    character(*), intent(in) :: filename
    type(SettingsBC), intent(inout) :: sbc
    character(:), allocatable, intent(out) :: err
    
    character(:), allocatable :: vel, bctype
    type(type_error), allocatable :: io_err

    if (bc_kind == "upper") then
      vel = 'veff'
    elseif (bc_kind == "lower") then
      vel = 'vdep'
    endif
    
    bctype = bc%get_string("type",error = io_err)
    if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
    if (bctype == vel) then
      sbc%bc_type = VelocityBC
      sbc%vel = bc%get_real(vel,error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
      sbc%mix = -huge(1.0_dp)
      sbc%flux = -huge(1.0_dp)
      sbc%height = -huge(1.0_dp)
      sbc%den = -huge(1.0_dp)

      if (sbc%vel < 0.0_dp) then
        err = 'Velocity '//trim(bc_kind)//' boundary condition for '//trim(sp_name)// &
              ' must be positive.'
        return
      endif
    elseif (bctype == "mix") then
      if (bc_kind == "upper") then
        err = 'Upper boundary conditions can not be "mix" for '//trim(sp_name)
        return
      endif
      
      sbc%bc_type = MixingRatioBC
      sbc%mix = bc%get_real("mix",error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
      sbc%vel = -huge(1.0_dp)
      sbc%flux = -huge(1.0_dp)
      sbc%height = -huge(1.0_dp)
      sbc%den = -huge(1.0_dp)

      if (sbc%mix < 0.0_dp .or. sbc%mix > 1.0_dp) then
        err = 'Fixed '//trim(bc_kind)//' boundary condition for '//trim(sp_name)// &
              ' must be between 0 and 1.'
        return
      endif
    elseif (bctype == "flux") then
      sbc%bc_type = FluxBC
      sbc%flux = bc%get_real("flux",error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
    
      sbc%vel = -huge(1.0_dp)
      sbc%mix = -huge(1.0_dp)
      sbc%height = -huge(1.0_dp)
      sbc%den = -huge(1.0_dp)
    elseif (bctype == "vdep + dist flux") then
      if (bc_kind == "upper") then
        err = 'Upper boundary conditions can not be "vdep + dist flux" for '//trim(sp_name)
        return
      endif
      
      sbc%bc_type = VelocityDistributedFluxBC
      sbc%vel = bc%get_real(vel,error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
      sbc%flux = bc%get_real("flux",error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
      sbc%height = bc%get_real("height",error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
      sbc%mix = -huge(1.0_dp)
      sbc%den = -huge(1.0_dp)

      if (sbc%vel < 0.0_dp) then
        err = 'Velocity '//trim(bc_kind)//' boundary condition for '//trim(sp_name)// &
              ' must be positive.'
        return
      endif
      if (sbc%flux < 0.0_dp) then
        err = 'Distributed flux in '//trim(bc_kind)//' boundary condition for '//trim(sp_name)// &
              ' must be positive.'
        return
      endif
      if (sbc%height < 0.0_dp) then
        err = 'Distributed height in '//trim(bc_kind)//' boundary condition for '//trim(sp_name)// &
              ' must be positive.'
        return
      endif
    elseif (bctype == "den") then
      if (bc_kind == "upper") then
        err = 'Upper boundary conditions can not be "den" for '//trim(sp_name)
        return
      endif
      
      sbc%bc_type = DensityBC
      sbc%den = bc%get_real("den",error = io_err)
      if (allocated(io_err)) then; err = trim(filename)//trim(io_err%message); return; endif
      
      sbc%vel = -huge(1.0_dp)
      sbc%mix = -huge(1.0_dp)
      sbc%flux = -huge(1.0_dp)
      sbc%height = -huge(1.0_dp)

      if (sbc%den < 0.0_dp) then
        err = 'Fixed density '//trim(bc_kind)//' boundary condition for '//trim(sp_name)// &
              ' must be greater than 1.'
        return
      endif
    elseif (bctype == "Moses") then
      sbc%bc_type = MosesBC
      
      sbc%vel = -huge(1.0_dp)
      sbc%mix = -huge(1.0_dp)
      sbc%flux = -huge(1.0_dp)
      sbc%height = -huge(1.0_dp)
      sbc%den = -huge(1.0_dp)
    else
      err = 'IOError: "'//trim(bctype)//'" is not a valid lower boundary condition for '//trim(sp_name)
      return
    endif
    
  end subroutine

  pure function check_for_duplicates(str_list) result(ind)
    character(*), intent(in) :: str_list(:)
    integer :: ind
    integer :: i, j
    ind = 0
    do i = 1,size(str_list)-1
      do j = i+1,size(str_list)
        if (str_list(i) == str_list(j)) then
          ind = i
          exit
        endif
      enddo
    enddo
  end function

  subroutine unpack_string_list(list, str_list, err)
    type(type_list), intent(in) :: list
    character(*), allocatable, intent(out) :: str_list(:)
    character(:), allocatable, intent(out) :: err
    
    integer :: i
    type(type_list_item), pointer :: item
    
    allocate(str_list(list%size()))
    i = 1
    item => list%first
    do while (associated(item))
      select type (it => item%node)
      class is (type_scalar)
        str_list(i) = trim(it%string)
      class default
        err = '"'//trim(it%path)//'" must be a scalar.'
        return
      end select
      i = i + 1
      item => item%next
    enddo
    
  end subroutine

end submodule