
!> This module contains all the derived datatypes used in WaveQLab.
!> Initially, these were distributed amongst the different source
!> files but that was causing circular dependencies. To maintain
!> source compatability with KD3D while killing the circular
!> dependencies, this approach was chosen.
module datatypes

   use common, only : wp
   use mpi3dcomm, only : cartesian3d_t
   use mpi3d_interface, only : interface3d
   use mpi3dio, only : file_distributed
 
   implicit none
 
   !> block_indices datatype originally created by Eric.
   !> Its functionality is redundant and should be deleted at some point.
   type :: block_indices
      integer :: nq,nr,ns !< number of grid points in each direction in computational coordinates q,r,s
      real(kind = wp) :: hq, hr, hs !< grid spacing along the 3 computational axes
   end type block_indices
 
   !> block_grid datatype that contains the grid positions, derivatives, and jacobians
   type :: block_grid_t
      type(cartesian3d_t) :: C
      real(kind = wp) :: hq, hr, hs, bhq, bhr, bhs      !< Grid spacing in unit cube (h*) and block (bh*)
      real(kind = wp),dimension(:,:,:,:), allocatable :: x !< coordinates
      real(kind = wp),dimension(:,:,:,:), allocatable :: metricx , &
                                         metricy , &
                                         metricz !< metric derivative
      real(kind = wp),dimension(:,:,:), allocatable :: J !> Jacobian
   end type block_grid_t
 
   !> block_material datatype to hold the material properties
   type :: block_material
      real(kind = wp),dimension(:,:,:,:), allocatable :: M !< material properties

      ! --- Anelastic (Q) / attenuation state (used only when response == 'anelastic')
      logical :: anelastic = .false.
      real(kind = wp) :: dt = 0.0_wp
      real(kind = wp) :: c = 1.0_wp
      real(kind = wp) :: weight_exp = 0.0_wp
      real(kind = wp) :: fref = 1.0_wp

      real(kind = wp), dimension(:,:,:), allocatable :: Qp_inv, Qs_inv
      real(kind = wp), dimension(4) :: tau = 0.0_wp
      real(kind = wp), dimension(4) :: weight = 0.0_wp

      ! Six stress-component memory variables (4 mechanisms each)
      real(kind = wp), dimension(:,:,:,:), allocatable :: eta4,eta5,eta6,eta7,eta8,eta9
      real(kind = wp), dimension(:,:,:,:), allocatable :: Deta4,Deta5,Deta6,Deta7,Deta8,Deta9
   end type block_material
 
   !> block_plastic datatype to hold Drucker-Prager plasticity variables
   type :: block_plastic
      real(kind = wp) :: mu_beta_eta(3) !< internal friction: mu, plastic dilatancy: beta and viscosity: eta
      real(kind = wp),dimension(:,:,:,:), allocatable :: P !< plastic array (plaststic strain, and strain rate) 
   end type block_plastic
 
   !> block_fields datatype holds the different velocity and stress fields,
   !> and their temporal derivatives
   type :: block_fields
      real(kind = wp),dimension(:,:,:,:), allocatable :: F !< fields array (velocities, stresses)
      real(kind = wp),dimension(:,:,:,:), allocatable :: DF !< rates array (stores time derivative of fields, dF/dt)     
   end type block_fields
   
   !> block_boundary datatype holds the different components of the 6 2D block boundaries.
   type :: block_boundary
      real(kind = wp),dimension(:,:,:), allocatable :: X !< grid values
      real(kind = wp),dimension(:,:,:), allocatable :: M !< material properties
      real(kind = wp),dimension(:,:,:), allocatable :: n_l , &
           n_m , &
           n_n !< unit normal vector
      real(kind = wp),dimension(:,:,:), allocatable :: F !< fields on boundary (copied from block_fields)
      real(kind = wp),dimension(:,:,:), allocatable :: DF !< rates on boundary (copied to block_rates)
      real(kind = wp),dimension(:,:,:), allocatable :: Fopp !< field values from other side
      real(kind = wp),dimension(:,:,:), allocatable :: Mopp !< material values from other side
      real(kind = wp),dimension(:,:,:), allocatable :: U , &
           DU !< displacements (and associated rate array)
   end type block_boundary
   
   !> boundary_type datatype created by Kenneth. Its functionality should be
   !> folded into the block_boundary above
   type :: boundary_type
      integer :: block_num
      integer :: Rx, Ry, Rz, Lx, Ly, Lz
   end type boundary_type
 
   type :: moment_tensor
       logical :: use_moment_tensor
   ! When a source lies on/near the shared interface plane between blocks,
   ! it may appear inside both blocks' bounding boxes. In that ambiguous case,
   ! we assign it to block 1 by default. Set this to 2 to prefer block 2 instead.
   integer :: tensor_block_preference = 1
       character(64), dimension(:), allocatable :: source_type
       real(kind = wp), dimension(:), allocatable :: mXX,mXY,mXZ,mYY,mYZ,mZZ
       real(kind = wp), dimension(:), allocatable :: location_x,location_y, location_z
       integer, dimension(:), allocatable :: near_x, near_y, near_z, alpha
       real(kind = wp), dimension(:), allocatable :: near_phys_x, near_phys_y, near_phys_z
       real(kind = wp), dimension(:,:,:,:), allocatable :: exact
       real(kind = wp), dimension(:), allocatable :: duration,t_init
       integer :: order,block_id,num_tensor
   end type moment_tensor  
 
   !> block_pml datatype holds the different components of the pml auxiliary variables and rates.
   type :: block_pml
      logical :: pml !< logical variable (true, false) denoting if PML needs to implemented
      integer :: N_pml !< number of grid points inside the PML
      type(cartesian3d_t) :: C
      real(kind = wp),dimension(:,:,:,:), allocatable :: Q !< PML auxiliary fields 
      real(kind = wp),dimension(:,:,:,:), allocatable :: DQ !< PML auxiliary rates 
   end type block_pml
   
 
   !> block datatype holds all the different pieces above.
   type :: block_type
      integer :: id !< ID of the block, is equal to 2^N where N is the block number
      character(64) :: physics !< physics (e.g., elastic or acoustic medium)
      character(64) :: fd_type
      integer :: order, nb
      type(block_indices) :: I !< indices and such for defining fields on a structured mesh
      type(block_grid_t) :: G !< spatial grid (physical coordinates and metric derivatives of mapping)
      type(block_material) :: M !< material properties
      type(block_fields) :: F !< fields (on mesh within block)
      type(moment_tensor) :: MT   !<moment tensor sources within block
      type(block_plastic) :: P !< plasticity variables 
      type(block_boundary) :: B(6) !< six sides on a block in 3D, with various arrays defined on them
      type(boundary_type) :: boundary_vars
      type(block_pml) :: PMLB(6) !< six sides on a block in 3D, with various PML auxiliary arrays defined on them
      real(kind = wp) :: tau0
      real(kind = wp) :: rho_s_p(3)
      real(kind = wp) :: sum
   end type block_type
 
   !> iface datatype, which includes various fields defined on the interface as well as
   !> information like indices of blocks on either side of interface
 
   type :: iface_type
      integer :: id !< ID of the interface, is equal to the sum of the two block IDs
      integer :: im,ip !< indices of minus and plus blocks
      character(1) :: direction !< direction of interface (q,r,s)
      real(kind = wp),dimension(:,:,:), allocatable :: V , &
                                       DV !< discontinuity in velocity vector across interface
      real(kind = wp),dimension(:,:,:), allocatable :: T !< traction vector on interface
      real(kind = wp),dimension(:,:,:), allocatable :: S , &
                                       DS !< displacement discontinuity across interface and its rates array
      real(kind = wp),dimension(:,:,:), allocatable :: W , &
                                       DW !< state variable and its rate
      real(kind = wp),dimension(:,:,:), allocatable :: Svel !< slip velocity across the interface
      real(kind = wp),dimension(:,:,:), allocatable :: trup !< slip velocity across the interface
      type(interface3d) :: II !< Interface communicator
   end type iface_type
 
   !> fault datatype for outputting the different values at the interface
   type fault_type
      type(file_distributed), dimension(8) :: handles
      integer :: array_s
      real(kind = wp),  allocatable, dimension(:, :, :) :: time_rup
      real(kind = wp),  allocatable, dimension(:, :, :) :: W
      real(kind = wp),  allocatable, dimension(:, :, :) :: slip , &
                                           Svel
      real(kind = wp),  allocatable, dimension(:, :, :) :: U_pluspres , &
                                           V_pluspres , &
                                           Uhat_pluspres , &
                                           Vhat_pluspres
   end type fault_type
 
   !> mms datatype
   type :: mms_type
      logical :: use_mms
      real(kind = wp) :: nx, ny, nz, nt ! spatial wavenumbers and frequency
   end type mms_type
 
   !> seismogram datatype
   type :: seismogram_type
      logical :: output_exact_moment
      logical :: output_seismograms,output_fields_block1,output_fields_block2
        logical :: output_station_info
      logical :: station_xyz_index
      integer :: stride_fields,file_unit_block1(9),file_unit_block2(9)
      integer :: nstations,block_num
      integer,dimension(:),allocatable :: i,j,k,file_unit
      real(kind = wp), dimension(:), allocatable :: i_phys,j_phys,k_phys
   end type seismogram_type
 
   !> slice datatype
   type slice_type
      type(file_distributed) :: Hslice,Vslice
      logical :: horz_slice,vert_slice
      integer :: locationH,locationV
      integer :: array_s
   end type slice_type

   !> plane output datatype (per block)
   type :: plane_output_plane
      logical :: enabled = .false.
      logical :: active = .false.
      character(len=2) :: axis = ''
      real(kind = wp) :: coord = 0.0_wp
      character(len=64) :: name = ''
      character(len=256) :: file_directory = '.'
      integer :: stride = 1
      integer :: step_counter = 0

      integer :: fixed_dir = 0
      integer :: fixed_index_global = 0

      integer :: plane_comm = -1
      integer :: plane_rank = -1
      integer :: plane_size = 0

      integer :: n1 = 0
      integer :: n2 = 0

      integer :: file_unit = -1
      logical :: header_written = .false.

      real(kind = wp), allocatable :: vx(:,:), vy(:,:), vz(:,:)
   end type plane_output_plane

   type :: plane_output_type
      logical :: enabled = .false.
      integer :: nplanes = 0
      type(plane_output_plane), allocatable :: P(:)
   end type plane_output_type
 
 
   !> domain datatype contains the different blocks and interfaces.
   !> For now, the output datatypes are also held here but that
   !> should be cleaned out.
 
   type :: domain_type
      integer :: nblocks,nifaces !< number of blocks and interfaces
      integer :: w_stride !< output stride
      real(kind = wp) :: dt, t, t_final !< time
      real(kind = wp) :: CFL !< CFL
      real(kind = wp) :: topo !< topography prefactor
      integer :: nt !< number of time steps
      type(block_type),allocatable :: B(:) !< blocks
      type(iface_type),allocatable :: I(:) !< interfaces
      type(fault_type) :: fault
      type(mms_type) :: mms_vars
      character(256) :: name, problem,response,plastic_model
      character(64) :: coupling, fd_type
      character(64) :: type_of_mesh, mesh_source,material_source
      logical :: w_fault, output_fault_topo, interpol,use_topography,mollify_source
      integer :: order !< order of accuracy of fd operator
      type(seismogram_type), allocatable :: seismometers(:)   !< seismograms
      type(slice_type) :: slicer
      type(plane_output_type), allocatable :: plane_outputs(:)
   end type domain_type
 
   type:: block_temp_parameters
      integer :: nqrs(3)
      real(kind = wp) :: aqrs(3), bqrs(3)
      real(kind = wp) :: lc, rc
      real(kind = wp) :: rho_s_p(3)
      real(kind = wp) :: mu_beta_eta(3)
      integer :: lqrs(3), rqrs(3)
      character(len=256) :: profile_type, profile_path,material_path(3)
      character(len=256) :: topography_type, topography_path
      logical :: pml_lqrs(3), pml_rqrs(3)
      integer :: npml, faultsize(2)
   end type block_temp_parameters
 
 end module datatypes
 
