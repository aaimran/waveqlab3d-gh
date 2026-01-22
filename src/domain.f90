!> @brief domain module holds multiple blocks and interfaces, and is used by time-stepping routines

!> @todo Hari, I've never been careful about public vs private in modules, maybe you can help with that?

!> domain = set of blocks and interfaces, with fields defined at time t
!> blocks and interfaces are allocatable to allow for any number of them,
!> but initialization is specific to two blocks and one interface (which I think
!> should be the initial focus since that is sufficient to do many interesting problems),
!> but we should generalize this to an arbitrary number of blocks and interfaces
!> with an unstructured connectivity pattern (presumably then the domain_type would also
!> contain some type that would handle the block/interface layout)

module domain

  use mpi
  use common, only : wp
  use datatypes, only : domain_type, block_type, iface_type, mms_type,fault_type,moment_tensor
  use block, only : block_temp_parameters, eval_block_mms, norm_fields_block
  use input_utils, only : locate_block_list_file
  use slice_output, only : init_slice_output, end_slice_output
  use seismogram
  use plane_output, only : init_plane_output, write_plane_output, end_plane_output

  implicit none

  logical :: in_block_comm(2), in_fault_comm(2), in_hslice_comm(2), in_vslice_comm(2)

contains


  subroutine init_domain(D, infile, input_filename)

    !> initialize the domain at t=0

    use block, only : init_block, block_time_step
    use iface, only : init_iface
    use fault_output, only : init_fault_output
    use mms, only : init_mms
    use interface_condition, only : init_vel_state
    use mpi3dbasic, only: is_master, rank, nprocs, new_communicator, warning, print_run_info
    use mpi3d_interface
    use, intrinsic :: iso_c_binding, only: c_int

    implicit none

    interface
      subroutine c_exit(status) bind(C, name="exit")
        import :: c_int
        integer(c_int), value :: status
      end subroutine c_exit
    end interface

    type(domain_type),intent(out) :: D
    integer, intent(in) :: infile
    character(*), intent(in), optional :: input_filename


    integer :: i, cart_size(3)
    integer :: n, ierr, cart_rank, comm_cart, coord(3), normal(3), &
               block_comms(2), fault_comms(2), hslice_comms(2), &
               vslice_comms(2), my_block_comm
    logical :: periodic(3) = (/ .false., .false., .false. /), reorder=.true.
    logical :: use_block_list_file
    type(interface3d) :: II

    ! Variables for parsing the input data file
    integer :: stat
    character(256) :: name, problem,response,plastic_model,fd_type
    character(256) :: response_norm
    logical :: invalid_response
    integer :: nt, nblocks
    logical :: output_fault_topo, w_fault, interpol, use_topography, mollify_source
    integer :: w_stride, ny, nz, order
    real(kind = wp) :: CFL, t_final, topo
    ! interface conditions
    character(64) :: coupling !< locked, slip-weakening_friction, linear_friction
    character(64) :: mesh_source, type_of_mesh, material_source !< cartesian or curvilinear

    ! Moment tensor info summary (optional)
    integer :: show_info_local, show_info_global
    integer :: use_mt_local, use_mt_global
    integer :: pref_local, pref_global
    logical :: show_moment_tensor_info
    integer :: total_point_sources, n_u, n_v, n_union
    logical :: have_union_list
    integer :: block_rank_local
    integer :: b1_assigned_local, b2_assigned_local
    integer :: b1_assigned, b2_assigned
    integer :: common_points
    integer, allocatable :: bbox1_local(:), bbox2_local(:), bbox1(:), bbox2(:)
    character(256) :: line

    real(kind = wp) :: dt, dtmin, dt1, dt2, spatialsetep1(3), spatialsetep2(3)
    
    real(kind = wp) :: meshvolume1, meshvolume2, totalvolume, ratio1, ratio2
    integer :: nprocs_1, nprocs_2
    integer :: block_list_unit

    type(block_temp_parameters) :: btp(2)
    character(len=512) :: block_list_file_path


    ! namelists to read in parameters

    namelist /problem_list/ name, problem,response,plastic_model,nblocks, &
                            nt, CFL, coupling, fd_type, order, t_final, &
                            mesh_source, type_of_mesh,material_source, &
                            interpol, w_stride, w_fault, use_topography, topo, mollify_source

    namelist /block_list/ btp


    !---------------------------------------------------------------------------
    !                       END OF VARIABLES
    !---------------------------------------------------------------------------
    ! set default parameters and read problem list

    name = 'default'
    problem = 'TPV5'
    response = 'elastic'
    plastic_model = 'default'
    nt = 0
    CFL = 0.5_wp
    coupling = 'locked'
    t_final = 0.0_wp
    w_fault = .true.
    output_fault_topo = .true.
    mesh_source = 'compute'
    material_source = 'hardcode'
    w_stride = 1
    interpol = .false.
    use_topography = .false.
    mollify_source = .false.
    topo = 1.0_wp
    fd_type = 'traditional'
    order = 5

    

    read(infile,nml=problem_list,iostat=stat)
    if (stat>0) stop 'error reading namelist problem_list'

     ! Validate `response` only on the master rank, then broadcast to all ranks.
     invalid_response = .false.
     if (is_master()) then
       response_norm = trim(adjustl(response))
       select case (response_norm)
       case ('elastic','plastic','anelastic')
         ! ok
       case default
         invalid_response = .true.
         write(*,*) ':: Error: invalid response option in &problem_list: "' // trim(response_norm) // '"'
         write(*,*) ':: Supported response options are: elastic, plastic, anelastic'
       end select
     end if

     call MPI_BCAST(invalid_response, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(response_norm, len(response_norm), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
     if (invalid_response) then
       call MPI_FINALIZE(ierr)
       call c_exit(1_c_int)
     end if

    if (present(input_filename)) then
      call locate_block_list_file(infile, input_filename = input_filename, block_list_file = block_list_file_path, has_file = use_block_list_file)
    else
      call locate_block_list_file(infile, block_list_file = block_list_file_path, has_file = use_block_list_file)
    end if

    if (is_master() .and. present(input_filename)) then
      call print_run_info(input_filename, name, problem)
      write(*,'(/,/)')
    end if

    D%name = name
    D%problem = problem
    D%response = response_norm
    D%plastic_model = plastic_model
    D%coupling = coupling
    D%w_fault = w_fault
    D%output_fault_topo = output_fault_topo
    D%t_final = t_final
    D%type_of_mesh = type_of_mesh
    D%mesh_source = mesh_source
    D%material_source = material_source
    D%CFL = CFL
    D%w_stride = w_stride
    D%interpol = interpol
    D%use_topography = use_topography
    D%topo = topo
    D%mollify_source = mollify_source 
    D%fd_type = fd_type
    D%order = order




    if (use_block_list_file) then
      open(newunit = block_list_unit, file = trim(block_list_file_path), status = 'old', action = 'read', iostat = stat)
      if (stat /= 0) stop 'error opening block list file '//trim(block_list_file_path)
      read(block_list_unit, nml = block_list, iostat = stat)
      if (stat > 0) then
        close(block_list_unit)
        stop 'error reading namelist block_list from '//trim(block_list_file_path)
      end if
      close(block_list_unit)
    else
      rewind(infile)
      read(infile,nml=block_list,iostat=stat)
      if (stat>0) stop 'error reading namelist block_list'
    end if


    
       

    !> set initial time

    D%t = 0.0_wp

    ! compute global time-step
    spatialsetep1(1) = (btp(1)%bqrs(1)-btp(1)%aqrs(1))/real(btp(1)%nqrs(1)-1)
    spatialsetep1(2) = (btp(1)%bqrs(2)-btp(1)%aqrs(2))/real(btp(1)%nqrs(2)-1)
    spatialsetep1(3) = (btp(1)%bqrs(3)-btp(1)%aqrs(3))/real(btp(1)%nqrs(3)-1)
    
    spatialsetep2(1) = (btp(2)%bqrs(1)-btp(2)%aqrs(1))/real(btp(2)%nqrs(1)-1)
    spatialsetep2(2) = (btp(2)%bqrs(2)-btp(2)%aqrs(2))/real(btp(2)%nqrs(2)-1)
    spatialsetep2(3) = (btp(2)%bqrs(3)-btp(2)%aqrs(3))/real(btp(2)%nqrs(3)-1)
       
    dt1 =   block_time_step(spatialsetep1, CFL, btp(1)%rho_s_p)
    dt2 =   block_time_step(spatialsetep2, CFL, btp(2)%rho_s_p)
    
    dtmin = min(dt1, dt2)
    
    D%dt = dtmin

        if (is_master()) call warning('warning: current method for setting time step does not use material properties; ' // &
          'mesh information from files is ignored (only scalars from btp%... in input file are used)', 'init_domain')
    
    !> set number of blocks and interfaces
    !> (for now, just two blocks and one interface)

    D%nblocks = 2
    D%nifaces = 1

    ! compute the volume of work in each block
    meshvolume1 = real(btp(1)%nqrs(1)*btp(1)%nqrs(2)*btp(1)%nqrs(3))
    meshvolume2 = real(btp(2)%nqrs(1)*btp(2)%nqrs(2)*btp(2)%nqrs(3))

    ! compute the work ratio
    totalvolume =  meshvolume1 +  meshvolume2
    ratio1 =  meshvolume1/totalvolume
    ratio2 =  meshvolume2/totalvolume

    ! assign processes using the work ratio
    nprocs_1 =  int((nprocs + 1)*ratio1)
    nprocs_2 =  int((nprocs + 1)*ratio2)

 
    ! For now, just two blocks.
    ! Can be generalized to N blocks by dividing by N instead of 2.
    in_block_comm(1) = (rank < (nprocs + 1)/2)
    in_block_comm(2) = (rank >= nprocs/2)

    !in_block_comm(1) = (rank .le. nprocs_1)
    !in_block_comm(2) = (rank > nprocs_1)

    call new_communicator(in_block_comm(1), block_comms(1))
    call new_communicator(in_block_comm(2), block_comms(2))

    if (in_block_comm(1)) my_block_comm = block_comms(1)
    if (in_block_comm(2)) my_block_comm = block_comms(2)

    allocate(D%B(D%nblocks))
    allocate(D%I(D%nifaces))
    allocate(D%seismometers(D%nblocks))
    allocate(D%plane_outputs(D%nblocks))

    !> initialize blocks
    do i = 1,D%nblocks
       ny = btp(i)%nqrs(2)
       nz = btp(i)%nqrs(3)

    if((btp(i)%profile_type == 'read_from_memomry_fractal') .and. &
         (btp(i)%topography_type == 'read_topo_from_memory')) then

       ny = btp(i)%faultsize(1)
       nz = btp(i)%faultsize(2)

    end if
       
      if (.not.in_block_comm(i)) cycle
      call init_block(D%mesh_source, D%type_of_mesh, D%material_source,&
           D%response, D%fd_type,  D%order, D%interpol, D%use_topography, topo, D%B(i), &
           problem, btp(i), block_comms(i),infile,i, ny, nz)      
  
      call MPI_Comm_size(block_comms(i),n,ierr)
      cart_size = (/0,0,0/)
      call MPI_Dims_create(n,3,cart_size,ierr)      

      call MPI_Cart_create(block_comms(i),3,cart_size,periodic,reorder,comm_cart,ierr)
      call MPI_Comm_rank(comm_cart,cart_rank,ierr)
      call MPI_Cart_coords(comm_cart,cart_rank,3,coord,ierr)

      
    end do

    ! Determine whether to print the moment tensor summary based on per-block settings.
    show_info_local = 0
    use_mt_local = 0
    pref_local = 1

    if (in_block_comm(1)) then
      call MPI_Comm_rank(block_comms(1), block_rank_local, ierr)
      if (block_rank_local == 0) then
        if (D%B(1)%MT%show_moment_tensor_info) show_info_local = 1
        if (D%B(1)%MT%use_moment_tensor) use_mt_local = 1
        pref_local = D%B(1)%MT%tensor_block_preference
      end if
    end if
    if (in_block_comm(2)) then
      call MPI_Comm_rank(block_comms(2), block_rank_local, ierr)
      if (block_rank_local == 0) then
        if (D%B(2)%MT%show_moment_tensor_info) show_info_local = 1
        if (D%B(2)%MT%use_moment_tensor) use_mt_local = 1
        pref_local = max(pref_local, D%B(2)%MT%tensor_block_preference)
      end if
    end if

    call MPI_Allreduce(show_info_local, show_info_global, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(use_mt_local, use_mt_global, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(pref_local, pref_global, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    show_moment_tensor_info = (show_info_global == 1)

    ! Print moment tensor source summary once (MPI-safe).
    if (show_moment_tensor_info .and. (use_mt_global == 1)) then
      have_union_list = .false.
      n_union = 0
      n_u = 0
      n_v = 0

      if (is_master()) then
        rewind(infile)
        do
          read(infile,'(a)', iostat=stat) line
          if (stat /= 0) exit
          if (trim(line) == '!---begin:tensor_list---') then
            have_union_list = .true.
            exit
          end if
        end do

        if (have_union_list) then
          rewind(infile)
          do
            read(infile,'(a)', iostat=stat) line
            if (stat /= 0) exit
            if (trim(line) == '!---begin:tensor_list---') exit
          end do
          do
            read(infile,'(a)', iostat=stat) line
            if (stat /= 0) exit
            if (trim(line) == '!---end:tensor_list---') exit
            line = adjustl(trim(line))
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '!') cycle
            n_union = n_union + 1
          end do
          total_point_sources = n_union
        else
          rewind(infile)
          do
            read(infile,'(a)', iostat=stat) line
            if (stat /= 0) exit
            if (trim(line) == '!---begin:tensor_listU---') exit
          end do
          do
            read(infile,'(a)', iostat=stat) line
            if (stat /= 0) exit
            if (trim(line) == '!---end:tensor_listU---') exit
            line = adjustl(trim(line))
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '!') cycle
            n_u = n_u + 1
          end do

          rewind(infile)
          do
            read(infile,'(a)', iostat=stat) line
            if (stat /= 0) exit
            if (trim(line) == '!---begin:tensor_listV---') exit
          end do
          do
            read(infile,'(a)', iostat=stat) line
            if (stat /= 0) exit
            if (trim(line) == '!---end:tensor_listV---') exit
            line = adjustl(trim(line))
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '!') cycle
            n_v = n_v + 1
          end do
          total_point_sources = n_u + n_v
        end if
      end if

      call MPI_BCAST(have_union_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(total_point_sources, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      b1_assigned_local = 0
      b2_assigned_local = 0

      if (in_block_comm(1)) then
        call MPI_Comm_rank(block_comms(1), block_rank_local, ierr)
        if (block_rank_local == 0) then
          if (allocated(D%B(1)%MT%near_x)) then
            b1_assigned_local = count((D%B(1)%MT%near_x /= -1) .and. (D%B(1)%MT%near_y /= -1) .and. (D%B(1)%MT%near_z /= -1))
          end if
        end if
      end if

      if (in_block_comm(2)) then
        call MPI_Comm_rank(block_comms(2), block_rank_local, ierr)
        if (block_rank_local == 0) then
          if (allocated(D%B(2)%MT%near_x)) then
            b2_assigned_local = count((D%B(2)%MT%near_x /= -1) .and. (D%B(2)%MT%near_y /= -1) .and. (D%B(2)%MT%near_z /= -1))
          end if
        end if
      end if

      call MPI_Allreduce(b1_assigned_local, b1_assigned, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(b2_assigned_local, b2_assigned, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

      common_points = 0
      if (have_union_list .and. total_point_sources > 0) then
        allocate(bbox1_local(total_point_sources), bbox2_local(total_point_sources))
        allocate(bbox1(total_point_sources), bbox2(total_point_sources))
        bbox1_local = 0
        bbox2_local = 0

        if (in_block_comm(1)) then
          call MPI_Comm_rank(block_comms(1), block_rank_local, ierr)
          if (block_rank_local == 0) then
            if (allocated(D%B(1)%MT%in_bbox)) then
              bbox1_local = merge(1, 0, D%B(1)%MT%in_bbox)
            end if
          end if
        end if

        if (in_block_comm(2)) then
          call MPI_Comm_rank(block_comms(2), block_rank_local, ierr)
          if (block_rank_local == 0) then
            if (allocated(D%B(2)%MT%in_bbox)) then
              bbox2_local = merge(1, 0, D%B(2)%MT%in_bbox)
            end if
          end if
        end if

        call MPI_Allreduce(bbox1_local, bbox1, total_point_sources, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
        call MPI_Allreduce(bbox2_local, bbox2, total_point_sources, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
        common_points = count((bbox1 == 1) .and. (bbox2 == 1))

        deallocate(bbox1_local, bbox2_local, bbox1, bbox2)
      end if

      if (is_master()) then
        write(*,'(a)') '+----------------------------------------------------------------------------+'
        write(*,'(a)') '| Moment Tensor Info                                                        |'
        write(*,'(a)') '+----------------------------------------------------------------------------+'

        write(line,'(a,i0)') 'Total point sources: ', total_point_sources
        write(*,'(a)') '| ' // trim(line) // repeat(' ', max(0, 74 - len_trim(line))) // ' |'

        write(line,'(a,i0)') 'In block 1: ', b1_assigned
        write(*,'(a)') '| ' // trim(line) // repeat(' ', max(0, 74 - len_trim(line))) // ' |'

        write(line,'(a,i0)') 'In block 2: ', b2_assigned
        write(*,'(a)') '| ' // trim(line) // repeat(' ', max(0, 74 - len_trim(line))) // ' |'

        write(line,'(a,i0)') 'Common points: ', common_points
        write(*,'(a)') '| ' // trim(line) // repeat(' ', max(0, 74 - len_trim(line))) // ' |'

        if (pref_global == 2) then
          line = 'Preference: block 2'
        else
          line = 'Preference: block 1'
        end if
        write(*,'(a)') '| ' // trim(line) // repeat(' ', max(0, 74 - len_trim(line))) // ' |'

        write(*,'(a)') '+----------------------------------------------------------------------------+'
      end if
    end if


    D%nt = floor(D%t_final/dtmin)
   
    if (is_master()) then
      write (*,*) name  
      write (*,*) "Domain time parameters :"
      write (*,*) "        Number of time steps = ", D%nt
      write (*,*) "        Initial dt = ", D%dt
      write (*,*) "        Final time = ", D%t_final
    end if

    ! Setup communicator for interface partitioning.


    if (in_block_comm(1)) then
      normal = [1, 0, 0]
      call new_interface(coord, cart_size, normal, MPI_COMM_WORLD, D%B(1)%G%C, II)
    end if
    if (in_block_comm(2).and. .not.in_block_comm(1)) then
      normal = [-1, 0, 0]
      call new_interface(coord, cart_size, normal, MPI_COMM_WORLD, D%B(2)%G%C, II)
    end if

    !> initialize interfaces

    do i = 1,D%nifaces
       !> create an interface that joins two blocks
       !> first define which blocks will be coupled: block 1 to block 2
       D%I(i)%im = 1
       D%I(i)%ip = 2
       !> and in which direction they will be coupled: q
       D%I(i)%direction = 'q'
       !> then initialize interface
       if (in_block_comm(1)) call init_iface(D%I(i),D%B(1), II)
       if (in_block_comm(2) .and. .not.in_block_comm(1)) call init_iface(D%I(i),D%B(2), II)
       if (in_block_comm(1)) call init_vel_state(D%problem, D%I(i), D%B(1), -1.0_wp, 2)
       if (in_block_comm(2) .and. .not.in_block_comm(1)) call init_vel_state(D%problem, D%I(i), D%B(2), 1.0_wp, 1)
    end do

    call exchange_materials_interface(D)

    in_fault_comm(1) = (in_block_comm(1) .and. D%I(1)%II%on_interface)
    in_fault_comm(2) = (in_block_comm(2) .and. D%I(1)%II%on_interface)
    call new_communicator(in_fault_comm(1), fault_comms(1))
    call new_communicator(in_fault_comm(2), fault_comms(2))

    !if ( D%w_fault .eqv.  .true.) then
    if (in_fault_comm(1)) call init_fault_output(D%w_fault,D%name, D%fault, D%B(1)%G%C, fault_comms(1))
    if (in_fault_comm(2) .and. .not.in_fault_comm(1)) call init_fault_output(D%w_fault,D%name, D%fault, D%B(2)%G%C, fault_comms(2))
    !end if
!     if (in_block_comm(1)) call init_slice_output(infile, D%name, D%slicer, D%B(1)%G%C)
!     if (in_block_comm(2)) call init_slice_output(infile, D%name, D%slicer, D%B(2)%G%C)

    call init_mms(infile, D%mms_vars)

    D%seismometers(1)%block_num = 1
    D%seismometers(2)%block_num = 2

    if (in_block_comm(1)) call init_seismogram(infile,D%seismometers(1),D%name, D%B(1)%G)
    if (in_block_comm(2)) call init_seismogram(infile,D%seismometers(2),D%name, D%B(2)%G)

    if (in_block_comm(1)) call init_plane_output(infile, D%name, D%plane_outputs(1), D%B(1)%G, 1)
    if (in_block_comm(2)) call init_plane_output(infile, D%name, D%plane_outputs(2), D%B(2)%G, 2)
    
    !> above assumes that blocks contain fields defined on a structured mesh,
    !> which is a cube in the computational domain; the computational coordinates
    !> are q,r,s, and there is a mapping to a curvilinear mesh in physical coordinates x,y,z
    !> I recommend that we use q,r,s in all references to directions, sides, etc., because
    !> it is possible that a boundary in q is mapped to a boundary in y or z rather than x
    !> (this will be a different naming convention than in original version of 3D code)

  end subroutine init_domain


  subroutine close_domain(D)

    use fault_output, only: destroy_fault

    type(domain_type),intent(inout) :: D

    if ( D%w_fault .eqv.  .true.) then
       if (in_fault_comm(1)) call destroy_fault(D%fault)
       if (in_fault_comm(2) .and. .not.in_fault_comm(1)) call destroy_fault(D%fault)
    end if
    !call end_slice_output(D%slicer)
    if (in_block_comm(1)) call destroy_seismogram(D%seismometers(1))
    if (in_block_comm(2)) call destroy_seismogram(D%seismometers(2))

    if (in_block_comm(1)) call end_plane_output(D%plane_outputs(1))
    if (in_block_comm(2)) call end_plane_output(D%plane_outputs(2))

  end subroutine close_domain

  subroutine enforce_bound_iface_conditions(D, stage)

    !> enforce boundary and interface conditions

    use block, only : enforce_bound_conditions

    implicit none

    type(domain_type),intent(inout) :: D
    integer, intent(in) :: stage

    !> enforce boundary conditions on external sides of each block
    !> enforce interface conditions to couple blocks

      if (in_block_comm(1)) then
        call enforce_bound_conditions(D%B(1), D%mms_vars, D%t)
        call enforce_iface_conditions(D%problem, D%coupling, D%I(1), &
          D%B(1),2,D%t, stage, D%mms_vars, D%fault)
      end if
      if (in_block_comm(2)) then
        call enforce_bound_conditions(D%B(2), D%mms_vars, D%t)
        call enforce_iface_conditions(D%problem, D%coupling, D%I(1), &
          D%B(2),1,D%t, stage, D%mms_vars, D%fault)
      end if

  end subroutine enforce_bound_iface_conditions

  subroutine enforce_iface_conditions(problem, coupling, I,B, ib, t, stage, mms_vars, handles)

    use datatypes, only : block_type, iface_type, mms_type, fault_type, block_boundary
    use RHS_Interior, only : Impose_Interface_Condition
    use mpi3dbasic, only : rank

    implicit none

    character(*), intent(in) :: problem, coupling
    type(iface_type),intent(inout) :: I
    type(block_type),intent(inout) :: B
    type(mms_type), intent(inout) :: mms_vars
    type(fault_type), intent(inout) :: handles
    real(kind = wp),intent(in) :: t
    integer, intent(in) :: stage, ib

    select case(I%direction)

    case('q')

!!! coupling: side 2 of block 1, B(2) <==> side 1 of block 2, Bp(1)

!!! this solves interface conditions for hat variables, constructs SAT forcing terms,
!!! and adds SAT forcing terms to rates

      if (B%boundary_vars%Rx == 0 .or. B%boundary_vars%Lx == 0) then
        call Impose_Interface_Condition(problem, coupling, I, B, ib, &
                              t, stage, mms_vars, handles)
      end if

    case('r','s')

       stop 'interfaces in r and s direction not implemented in enforce_iface_conditions'

    end select

  end subroutine enforce_iface_conditions

  subroutine exchange_fields(D)

    use block, only : exchange_fields_block

    implicit none

    type(domain_type),intent(inout) :: D
    integer :: i

    if (in_block_comm(1)) call exchange_fields_block(D%B(1))
    if (in_block_comm(2)) call exchange_fields_block(D%B(2))

  end subroutine exchange_fields

  subroutine exchange_fields_interface(D)

    use mpi3dbasic, only : nprocs
    use block, only : copy_fields_to_boundary
    use boundary, only : exchange_fields_across_interface

    implicit none
    type(domain_type),intent(inout) :: D

      if (nprocs == 1) then
        call copy_fields_to_boundary(D%B(1))
        call copy_fields_to_boundary(D%B(2))
        D%B(1)%B(2)%Fopp(:,:,:) = D%B(2)%B(1)%F(:,:,:)
        D%B(2)%B(1)%Fopp(:,:,:) = D%B(1)%B(2)%F(:,:,:)
      else
        if (in_block_comm(1)) then
          call copy_fields_to_boundary(D%B(1))
          call exchange_fields_across_interface(D%B(1)%B(2), D%B(1)%G%C, D%I(1)%II)
        end if
        if (in_block_comm(2)) then
          call copy_fields_to_boundary(D%B(2))
          call exchange_fields_across_interface(D%B(2)%B(1), D%B(2)%G%C, D%I(1)%II)
        end if
      end if

  end subroutine exchange_fields_interface

  subroutine exchange_materials_interface(D)

    use mpi3dbasic, only : nprocs
    use block, only : copy_fields_to_boundary
    use boundary, only : exchange_materials_across_interface

    implicit none
    type(domain_type),intent(inout) :: D

      if (nprocs == 1) then
        D%B(1)%B(2)%Mopp = D%B(2)%B(1)%M
        D%B(2)%B(1)%Mopp = D%B(1)%B(2)%M
      else
        if (in_block_comm(1)) then
          call exchange_materials_across_interface(D%B(1)%B(2), D%B(1)%G%C, D%I(1)%II)
        end if
        if (in_block_comm(2)) then
          call exchange_materials_across_interface(D%B(2)%B(1), D%B(2)%G%C, D%I(1)%II)
        end if
      end if

  end subroutine exchange_materials_interface

  subroutine write_output(D)

    !> @brief write fields (and, in some cases, rates) at time t
    !> note that the way it is written exposes the details of how fields are stored
    !> in blocks and on interfaces; a potentially better way is to use subroutines
    !> that retrieve certain fields and return them in an output array

    use,intrinsic :: iso_fortran_env, only : output_unit

    implicit none

    type(domain_type),intent(inout) :: D

    if ( D%w_fault .eqv.  .true.) then
       call write_fault_output(D)
       !call write_slice_output(D)
       call write_hat_output(D)
    end if

     if (in_block_comm(1)) call write_plane_output(D%plane_outputs(1), D%t, D%B(1)%F%F, D%B(1)%G)
     if (in_block_comm(2)) call write_plane_output(D%plane_outputs(2), D%t, D%B(2)%F%F, D%B(2)%G)
    
    call write_seismogram_output(D)

  end subroutine write_output

  subroutine write_fault_output(D)

    use fault_output
    use mpi3dbasic, only : rank
    implicit none

    type(domain_type), intent(in) :: D
    integer :: mq1, mr1, ms1, pq1, pr1, ps1, mq2, mr2, ms2, pq2, pr2, ps2

    if ( D%w_fault .eqv.  .true.) then
       if (in_fault_comm(1)) then
          
          mq1 = D%B(1)%G%C%mq
          mr1 = D%B(1)%G%C%mr
          ms1 = D%B(1)%G%C%ms
          pq1 = D%B(1)%G%C%pq
          pr1 = D%B(1)%G%C%pr
          ps1 = D%B(1)%G%C%ps
          
          call write_fault(D%B(1)%F%F(pq1,mr1:pr1,ms1:ps1,:),  &
               D%I(1)%S(mr1:pr1,ms1:ps1,:), &
               D%I(1)%W(mr1:pr1,ms1:ps1,:), D%fault)
          
       end if
       if (in_fault_comm(2)) then
          
          mq2 = D%B(2)%G%C%mq
          mr2 = D%B(2)%G%C%mr
          ms2 = D%B(2)%G%C%ms
          pq2 = D%B(2)%G%C%pq
          pr2 = D%B(2)%G%C%pr
          ps2 = D%B(2)%G%C%ps
          
          call write_file_distributed(D%fault%handles(2), D%B(2)%F%F(mq2,mr2:pr2,ms2:ps2,:))
       end if
    end if


  end subroutine write_fault_output

  subroutine write_hat_output(D)
    use fault_output
    implicit none

    type(domain_type), intent(in) :: D
    integer :: mq1, mr1, ms1, pq1, pr1, ps1

    if ( D%w_fault .eqv.  .true.) then

       if (.not.in_fault_comm(1)) return
       
       mq1 = D%B(1)%G%C%mq
       mr1 = D%B(1)%G%C%mr
       ms1 = D%B(1)%G%C%ms
       pq1 = D%B(1)%G%C%pq
       pr1 = D%B(1)%G%C%pr
       ps1 = D%B(1)%G%C%ps
       
       call write_hats(D%fault%Uhat_pluspres(mr1:pr1,ms1:ps1,:), &
            D%fault%Vhat_pluspres(mr1:pr1,ms1:ps1,:), &
            D%fault%Uhat_pluspres(mr1:pr1,ms1:ps1,1:3), &
            D%fault%time_rup(mr1:pr1,ms1:ps1,1),D%fault)
       
    end if
       
     end subroutine write_hat_output

  subroutine write_slice_output(D)

    use slice_output
    implicit none

    type(domain_type), intent(in) :: D

    if (in_block_comm(1)) call write_slice(D%B(1)%F%F,D%B(1)%G%C, D%slicer)

  end subroutine write_slice_output

  subroutine write_seismogram_output(D)

    use slice_output
    implicit none

    type(domain_type), intent(in) :: D
    integer :: mq1, mr1, ms1, pq1, pr1, ps1, mq2, mr2, ms2, pq2, pr2, ps2

    if(in_block_comm(1)) then
        mq1 = D%B(1)%G%C%mq
        mr1 = D%B(1)%G%C%mr
        ms1 = D%B(1)%G%C%ms
        pq1 = D%B(1)%G%C%pq
        pr1 = D%B(1)%G%C%pr
        ps1 = D%B(1)%G%C%ps

        call write_seismogram(D%seismometers(1), D%t, D%B(1)%F%F)

    end if

    if(in_block_comm(2)) then
        mq2 = D%B(2)%G%C%mq
        mr2 = D%B(2)%G%C%mr
        ms2 = D%B(2)%G%C%ms
        pq2 = D%B(2)%G%C%pq
        pr2 = D%B(2)%G%C%pr
        ps2 = D%B(2)%G%C%ps

        call write_seismogram(D%seismometers(2), D%t, D%B(2)%F%F) 

    end if


  end subroutine write_seismogram_output

  subroutine eval_mms(D)

    implicit none
    type(domain_type),intent(inout) :: D

    integer :: i

    do i = 1, D%nblocks
      if (.not.in_block_comm(i)) cycle
       call eval_block_mms(D%B(i), D%t, D%mms_vars)
    end do

  end subroutine eval_mms

  subroutine norm_fields(D)

    implicit none
    type(domain_type),intent(inout) :: D

    integer :: i

     if (in_block_comm(1)) call norm_fields_block(D%B(1))
     if (in_block_comm(2)) call norm_fields_block(D%B(2))

  end subroutine norm_fields

  subroutine scale_rates(D,A)

    !> @brief multiply all rates by RK coefficient A

    use block, only : scale_rates_block
    use iface, only : scale_rates_iface

    implicit none

    type(domain_type),intent(inout) :: D
    real(kind = wp),intent(in) :: A

    integer :: i

    ! first within blocks and on their boundaries

      if (in_block_comm(1)) call scale_rates_block(D%B(1),A)
      if (in_block_comm(2)) call scale_rates_block(D%B(2),A)


    ! and then on interfaces

    do i = 1,D%nifaces
       call scale_rates_iface(D%I(i),A)
    end do

  end subroutine scale_rates

  subroutine set_rates(D)

    !> @brief set rates using the PDE
    !> @details set rates using the PDE (in a low storage RK method, the new rates
    !> at the current stage are added to the old rates, instead of overwriting
    !> the rates array)

    use block, only : set_rates_block
    use moment_tensor, only : set_moment_tensor, moment_tensor_body_force, set_moment_tensor_smooth

    implicit none

    type(domain_type),intent(inout) :: D

       if (in_block_comm(1)) call set_rates_block(D%B(1), D%type_of_mesh)
       if (in_block_comm(2)) call set_rates_block(D%B(2), D%type_of_mesh)


    if (D%B(1)%MT%use_moment_tensor) then
       !call set_moment_tensor(D%B(1),D%t)
       if(D%mollify_source) then
          
          call set_moment_tensor_smooth(D%B(1),D%t)

       else
          
          call set_moment_tensor(D%B(1),D%t)
       end if
       
       !call moment_tensor_body_force(D%B(1),D%t) 
    end if

    if (D%B(2)%MT%use_moment_tensor) then
       !call set_moment_tensor(D%B(2),D%t)

       if(D%mollify_source) then
          
          call set_moment_tensor_smooth(D%B(2),D%t)

       else
          
          call set_moment_tensor(D%B(2),D%t)
       end if
       
       !call set_moment_tensor2(D%B(2),D%t)
       
       !call moment_tensor_body_force(D%B(2),D%t)
    end if
  end subroutine set_rates


  subroutine update_fields(D,dt,stage,RKstage)

    !> @brief use rates to update fields

    use block, only : update_fields_block
    use iface, only : update_fields_iface
    use plastic, only : update_fields_plastic
    implicit none

    type(domain_type),intent(inout) :: D
    real(kind = wp),intent(in) :: dt
    integer,intent(in) :: stage,RKstage

    integer :: i

    ! first within blocks and on their boundaries

      if (in_block_comm(1)) call update_fields_block(D%B(1),dt)
      if (in_block_comm(2)) call update_fields_block(D%B(2),dt)


    ! and then on interfaces

    do i = 1,D%nifaces
       call update_fields_iface(D%I(i),dt)
    end do

    if (stage==RKstage) then
       !if (in_block_comm(1)) call update_fields_plastic(D%B(1),D%B(1)%P,D%B(1)%G,D%B(1)%M,D%dt,D%t,D%problem,D%response)
       !if (in_block_comm(2)) call update_fields_plastic(D%B(2),D%B(2)%P,D%B(2)%G,D%B(2)%M,D%dt,D%t,D%problem,D%response)
       if (in_block_comm(1)) call update_fields_plastic(D%B(1),D%B(1)%P,D%B(1)%G,D%B(1)%M,D%dt,D%t,D%problem,D%response,&
            D%plastic_model)
       if (in_block_comm(2)) call update_fields_plastic(D%B(2),D%B(2)%P,D%B(2)%G,D%B(2)%M,D%dt,D%t,D%problem,D%response,&
            D%plastic_model)
    end if
  end subroutine update_fields

end module domain
