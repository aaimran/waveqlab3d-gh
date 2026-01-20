module material

!> material module defines material properties on a block
  use common, only : wp
  use datatypes, only : block_material,block_grid_t,block_indices
  implicit none

contains


   subroutine init_anelastic_properties(M, G, infile)

      ! Initializes GSLS attenuation parameters and allocates Q and memory-variable arrays.
      ! Must be called only when response == 'anelastic'.

      use mpi3dcomm, only : allocate_array_body

      implicit none

      type(block_material), intent(inout) :: M
      type(block_grid_t), intent(in) :: G
      integer, intent(in) :: infile

      real(kind = wp) :: c, weight_exp, fref
      real(kind = wp) :: FAC, taumin, taumax, wref, facex
      real(kind = wp) :: val_S, val_P, denom_S, denom_P, vs, vp, mu_unrelax_S, mu_unrelax_P
      integer :: stat, i, l, j, k, N
      real(kind = wp), parameter :: pi = 3.141592653589793_wp

      namelist /anelastic_list/ c, weight_exp, fref

      c = 1.0_wp
      weight_exp = 0.0_wp
      fref = 1.0_wp
      M%anelastic = .true.

      rewind(infile)
      read(infile, nml=anelastic_list, iostat=stat)
      if (stat > 0) stop 'error reading namelist anelastic_list'

      M%c = c
      M%weight_exp = weight_exp
      M%fref = fref

      call allocate_array_body(M%Qp_inv, G%C, ghost_nodes=.true.)
      M%Qp_inv(:,:,:) = 0.0_wp
      call allocate_array_body(M%Qs_inv, G%C, ghost_nodes=.true.)
      M%Qs_inv(:,:,:) = 0.0_wp

      N = 4

      call allocate_array_body(M%eta4,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta4, G%C, N, ghost_nodes=.true.)
      M%eta4(:,:,:,:) = 0.0_wp
      M%Deta4(:,:,:,:) = 0.0_wp
      call allocate_array_body(M%eta5,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta5, G%C, N, ghost_nodes=.true.)
      M%eta5(:,:,:,:) = 0.0_wp
      M%Deta5(:,:,:,:) = 0.0_wp
      call allocate_array_body(M%eta6,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta6, G%C, N, ghost_nodes=.true.)
      M%eta6(:,:,:,:) = 0.0_wp
      M%Deta6(:,:,:,:) = 0.0_wp
      call allocate_array_body(M%eta7,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta7, G%C, N, ghost_nodes=.true.)
      M%eta7(:,:,:,:) = 0.0_wp
      M%Deta7(:,:,:,:) = 0.0_wp
      call allocate_array_body(M%eta8,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta8, G%C, N, ghost_nodes=.true.)
      M%eta8(:,:,:,:) = 0.0_wp
      M%Deta8(:,:,:,:) = 0.0_wp
      call allocate_array_body(M%eta9,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta9, G%C, N, ghost_nodes=.true.)
      M%eta9(:,:,:,:) = 0.0_wp
      M%Deta9(:,:,:,:) = 0.0_wp

      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps
               M%Qs_inv(i,j,k) = 1.0_wp/(c*sqrt(M%M(i,j,k,2)/M%M(i,j,k,3)))
               M%Qp_inv(i,j,k) = 0.5_wp*M%Qs_inv(i,j,k)
            end do
         end do
      end do

      if (M%weight_exp < 0.01_wp) then
         FAC = 1.0_wp
         taumin = 1.0_wp/(2.0_wp*pi*15.0_wp)*1.0_wp*FAC
         taumax = 1.0_wp/(2.0_wp*pi*0.08_wp)*200.0_wp*FAC
         do k = 1, N
            M%tau(k) = exp(log(taumin) + (2.0_wp*k-1.0_wp)/16.0_wp*(log(taumax) - log(taumin)))
         end do
         M%weight(1) = 1.6126_wp
         M%weight(2) = 0.6255_wp
         M%weight(3) = 0.6382_wp
         M%weight(4) = 1.5969_wp
      end if

      if (M%weight_exp > 0.59_wp .AND. M%weight_exp < 0.61_wp) then
         FAC = 1.0_wp
         taumin = 1.0_wp/(2.0_wp*pi*15.0_wp)*1.0_wp*FAC
         taumax = 1.0_wp/(2.0_wp*pi*0.08_wp)*200.0_wp*FAC
         do k = 1, N
            M%tau(k) = exp(log(taumin) + (2.0_wp*k-1.0_wp)/16.0_wp*(log(taumax) - log(taumin)))
         end do
         M%weight(1) = 0.0336_wp
         M%weight(2) = 0.6873_wp
         M%weight(3) = 0.8767_wp
         M%weight(4) = 1.5202_wp
      end if

      wref = 2.0_wp*pi*fref
      facex = 1.0_wp

      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps

               val_S = 0.0_wp
               val_P = 0.0_wp
               do l = 1, N
                  denom_S = (wref**2.0_wp*M%tau(l)**2.0_wp + 1.0_wp) * (1.0_wp/M%Qs_inv(i,j,k)) * facex
                  denom_P = (wref**2.0_wp*M%tau(l)**2.0_wp + 1.0_wp) * (1.0_wp/M%Qp_inv(i,j,k)) * facex
                  val_S = val_S + M%weight(l)/denom_S
                  val_P = val_P + M%weight(l)/denom_P
               end do

               vs = sqrt(M%M(i,j,k,2)/M%M(i,j,k,3))
               vp = sqrt((M%M(i,j,k,1) + 2.0_wp*M%M(i,j,k,2))/M%M(i,j,k,3))
               mu_unrelax_S = M%M(i,j,k,3)*vs**2.0_wp/(1.0_wp - val_S)
               mu_unrelax_P = M%M(i,j,k,3)*vp**2.0_wp/(1.0_wp - val_P)
               vs = sqrt(mu_unrelax_S/M%M(i,j,k,3))
               vp = sqrt(mu_unrelax_P/M%M(i,j,k,3))

               M%M(i,j,k,2) = vs**2.0_wp*M%M(i,j,k,3)
               M%M(i,j,k,1) = vp**2.0_wp*M%M(i,j,k,3) - 2.0_wp*M%M(i,j,k,2)

            end do
         end do
      end do

   end subroutine init_anelastic_properties


  !> initialize material properties
  subroutine init_material(M, G, I, physics,problem, rho_s_p, nb)

    use mpi3dcomm

    implicit none

    type(block_material),intent(out) :: M
    type(block_grid_t),intent(in) :: G
    type(block_indices),intent(in) :: I
    real(kind = wp), intent(in) :: rho_s_p(3)
    character(*),intent(in) :: physics,problem
    integer, intent(in) :: nb
    integer :: ii,l,j,k, np
    real(kind = wp) :: Vp,Vs,rho

    !> the number of material properties will be specific to the physics
    !> (I included the acoustic case to illustrate this, but I don't foresee
    !> it being a priority to implement all the necessary routines to solve
    !> the acoustic wave equation)

   
    select case(physics)

    case default

       stop 'invalid block physics in init_material'

    case('elastic')

       !> note use of the array allocation routine, which makes it easier
       !> to allocate arrays without having to make sure that indices are being
       !> typed correctly

       select case(problem)

       case default

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.)

          M%M(:,:,:,1) = rho_s_p(1)*(rho_s_p(3)**2 - 2.0_wp*rho_s_p(2)**2) ! lambda (Lame's first parameter)
          M%M(:,:,:,2) = rho_s_p(1)*rho_s_p(2)**2 ! mu (shear modulus)
          M%M(:,:,:,3) = rho_s_p(1) ! rho (density)

       case('LOH1')

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.) 

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,1) < 0.9999_wp) then

                      Vp = 4.0_wp ! P-wave speed
                      Vs  = 2.0_wp ! S-wave speed
                      rho  = 2.6_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) >= 0.9999_wp .AND. G%x(l,j,k,1) <= 1.0001_wp) then

                      Vp = 5.0_wp ! P-wave speed
                      Vs  = 2.732_wp ! S-wave speed
                      rho  = 2.65_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) > 1.0001_wp) then

                      Vp = 6_wp    ! P-wave speed 
                      Vs = 3.464_wp    ! S-wave speed 
                      rho = 2.7_wp    ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   end if
                   !print *, l, j, k, M%M(l,j,k,1)
                   !if (l == 25 .and. j == 94 .and. k == 1) print *, l, j, k, M%M(l,j,k,1)
                end do
             end do
          end do

       case('LOH1_Harmonic')

          call allocate_array_body(M%M,G%C,3, ghost_nodes=.true.)

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,1) < 0.9999_wp) then

                      Vp = 4_wp ! P-wave speed
                      Vs  = 2_wp ! S-wave speed
                      rho  = 2.6_wp ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2)   ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2                    ! mu (shear modulus)
                      M%M(l,j,k,3) = rho                          ! rho (density)

                   else if (G%x(l,j,k,1) >= 0.9999_wp .AND. G%x(l,j,k,1) <= 1.0001_wp) then

                      Vp   = 2.0_wp*4.0_wp*6.0_wp/(4.0_wp + 6.0_wp)           ! P-wave speed
                      Vs   = 2.0_wp*2.0_wp*3.464_wp/(2.0_wp + 3.464_wp)       ! S-wave speed
                      rho  = 2.0_wp*2.7_wp*2.6_wp/(2.7_wp + 2.6_wp)           ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2)   ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2                    ! mu (shear modulus)
                      M%M(l,j,k,3) = rho                          ! rho (density)

                   else if (G%x(l,j,k,1) > 1.0001_wp) then

                      Vp = 6_wp                                   ! P-wave speed
                      Vs = 3.464_wp                               ! S-wave speed
                      rho = 2.7_wp                                ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   end if
                end do
             end do
          end do

       case('TPV31')

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.)

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,2) < 2.3999_wp) then

                      Vp = 4.05_wp ! P-wave speed
                      Vs  = 2.25_wp ! S-wave speed
                      rho  = 2.58_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,2) >= 2.3999_wp .AND. G%x(l,j,k,2) <= 2.4001_wp) then

                      Vp = 4.25_wp ! P-wave speed
                      Vs  = 2.4_wp ! S-wave speed
                      rho  = 2.59_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,2) > 2.4001_wp .AND. G%x(l,j,k,2) < 4.9999_wp) then

                      Vp  = 4.45_wp + 0.75_wp*(G%x(l,j,k,2)-2.4_wp)/2.6_wp  ! P-wave speed
                      Vs  = 2.55_wp + 0.5_wp*(G%x(l,j,k,2)-2.4_wp)/2.6_wp ! S-wave speed 
                      rho  = 2.6_wp + 0.02_wp*(G%x(l,j,k,2)-2.4_wp)/2.6_wp  ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,2) >= 4.9999_wp .AND. G%x(l,j,k,2) <= 5.0001_wp) then

                      Vp = 5.475_wp ! P-wave speed
                      Vs  = 3.25_wp ! S-wave speed
                      rho  = 2.67_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,2) > 5.0001_wp .AND. G%x(l,j,k,2) < 9.9999_wp) then

                      Vp = 5.75_wp   ! P-wave speed 
                      Vs = 3.45_wp   ! S-wave speed 
                      rho = 2.72_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,2) >= 9.9999_wp .AND. G%x(l,j,k,2) <= 10.0001_wp) then

                      Vp = 6.125_wp ! P-wave speed
                      Vs  = 3.625_wp ! S-wave speed
                      rho  = 2.76_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,2) > 10.0001_wp) then

                      Vp = 6.5_wp    ! P-wave speed 
                      Vs = 3.8_wp    ! S-wave speed 
                      rho = 3.0_wp    ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   end if

                end do
             end do
          end do


       case('TPV32')

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.)

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr 
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,2) <= 0.5_wp) then

                      Vp = 2.2_wp + 0.8_wp*(G%x(l,j,k,2))/0.5_wp      ! P-wave speed 
                      Vs = 1.05_wp + 0.35_wp*(G%x(l,j,k,2))/0.5_wp    ! S-wave speed 
                      rho = 2.2_wp + 0.25_wp*(G%x(l,j,k,2))/0.5_wp     ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 0.5_wp .AND. G%x(l,j,k,2) <= 1.0_wp) then

                      Vp = 3.0_wp + 0.6_wp*(G%x(l,j,k,2)-0.5_wp)/0.5_wp    ! P-wave speed 
                      Vs = 1.4_wp + 0.55_wp*(G%x(l,j,k,2)-0.5_wp)/0.5_wp   ! S-wave speed 
                      rho = 2.45_wp + 0.1_wp*(G%x(l,j,k,2)-0.5_wp)/0.5_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,2) > 1.0_wp .AND. G%x(l,j,k,2) <= 1.6_wp) then

                      Vp = 3.6_wp + 0.8_wp*(G%x(l,j,k,2)-1.0_wp)/0.6_wp     ! P-wave speed 
                      Vs = 1.95_wp + 0.55_wp*(G%x(l,j,k,2)-1.0_wp)/0.6_wp   ! S-wave speed 
                      rho = 2.55_wp + 0.05_wp*(G%x(l,j,k,2)-1.0_wp)/0.6_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 1.6_wp .AND. G%x(l,j,k,2) <= 2.4_wp) then

                      Vp = 4.4_wp + 0.4_wp*(G%x(l,j,k,2)-1.6_wp)/0.8_wp   ! P-wave speed 
                      Vs = 2.5_wp + 0.3_wp*(G%x(l,j,k,2)-1.6_wp)/0.8_wp   ! S-wave speed 
                      rho = 2.6_wp                                         ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 2.4_wp .AND. G%x(l,j,k,2) <= 3.6_wp) then

                      Vp = 4.8_wp + 0.45_wp*(G%x(l,j,k,2)-2.4_wp)/1.2_wp  ! P-wave speed 
                      Vs = 2.8_wp + 0.3_wp*(G%x(l,j,k,2)-2.4_wp)/1.2_wp   ! S-wave speed 
                      rho = 2.6_wp + 0.02_wp*(G%x(l,j,k,2)-2.4_wp)/1.2_wp  ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 3.6_wp .AND. G%x(l,j,k,2) <= 5.0_wp) then

                      Vp = 5.25_wp + 0.25_wp*(G%x(l,j,k,2)-3.6_wp)/1.4_wp  ! P-wave speed 
                      Vs = 3.1_wp + 0.15_wp*(G%x(l,j,k,2)-3.6_wp)/1.4_wp   ! S-wave speed 
                      rho = 2.62_wp + 0.03_wp*(G%x(l,j,k,2)-3.6_wp)/1.4_wp  ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 5.0_wp .AND. G%x(l,j,k,2) <= 9.0_wp) then

                      Vp = 5.5_wp + 0.25_wp*(G%x(l,j,k,2)-5.0_wp)/4.0_wp   ! P-wave speed 
                      Vs = 3.25_wp + 0.2_wp*(G%x(l,j,k,2)-5.0_wp)/4.0_wp   ! S-wave speed 
                      rho = 2.65_wp + 0.07_wp*(G%x(l,j,k,2)-5.0_wp)/4.0_wp  ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 9.0_wp .AND. G%x(l,j,k,2) <= 11.0_wp) then

                      Vp = 5.75_wp + 0.35_wp*(G%x(l,j,k,2)-9.0_wp)/2.0_wp   ! P-wave speed 
                      Vs = 3.45_wp + 0.25_wp*(G%x(l,j,k,2)-9.0_wp)/2.0_wp   ! S-wave speed 
                      rho = 2.72_wp + 0.03_wp*(G%x(l,j,k,2)-9.0_wp)/2.0_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 11.0_wp .AND. G%x(l,j,k,2) <= 15.0_wp) then

                      Vp = 6.1_wp + 0.2_wp*(G%x(l,j,k,2)-11.0_wp)/4.0_wp     ! P-wave speed 
                      Vs = 3.6_wp + 0.1_wp*(G%x(l,j,k,2)-11.0_wp)/4.0_wp     ! S-wave speed 
                      rho = 2.75_wp + 0.15_wp*(G%x(l,j,k,2)-11.0_wp)/4.0_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 15.0_wp) then

                      Vp = 6.3_wp     ! P-wave speed 
                      Vs = 3.7_wp     ! S-wave speed 
                      rho = 2.9_wp     ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   end if

                end do
             end do
          end do

          ! Low velocity zone benchmark
       case('TPV33')

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.)

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr 
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,1) <= -0.8_wp) then

                      Vp = 5.626_wp      ! P-wave speed 
                      Vs = 3.248_wp      ! S-wave speed 
                      rho = 2.67_wp      ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,1) >= 0.8_wp) then

                      Vp = 6.0_wp        ! P-wave speed 
                      Vs = 3.464_wp      ! S-wave speed 
                      rho = 2.67_wp      ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > -0.8_wp .AND. G%x(l,j,k,1) < 0.8_wp) then

                      Vp = 3.75_wp       ! P-wave speed 
                      Vs = 2.165_wp      ! S-wave speed 
                      rho = 2.67_wp      ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   end if

                end do
             end do
          end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case('OKLAHOMA')

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.)

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,1) < 1.4999_wp) then

                      Vp = 2.77_wp ! P-wave speed
                      Vs  = 1.49_wp ! S-wave speed
                      rho  = 2.169_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) >= 1.4999_wp .AND. G%x(l,j,k,1) <= 1.5001_wp) then

                      Vp = 3.74_wp ! P-wave speed
                      Vs  = 2.01_wp ! S-wave speed
                      rho  = 2.418_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) > 1.5001_wp .AND. G%x(l,j,k,1) < 2.2999_wp) then

                      Vp  = 5.76_wp  ! P-wave speed
                      Vs  = 3.1_wp  ! S-wave speed 
                      rho  = 2.667_wp  ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 2.2999_wp .AND. G%x(l,j,k,1) <= 2.3001_wp) then

                      Vp = 5.755_wp ! P-wave speed
                      Vs  = 3.26_wp ! S-wave speed
                      rho  = 2.666_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) > 2.3001_wp .AND. G%x(l,j,k,1) < 5.2999_wp) then

                      Vp = 5.75_wp   ! P-wave speed 
                      Vs = 3.423_wp   ! S-wave speed 
                      rho = 2.665_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 5.2999_wp .AND. G%x(l,j,k,1) <= 5.3001_wp) then

                      Vp = 5.96_wp ! P-wave speed
                      Vs  = 3.54_wp ! S-wave speed
                      rho  = 2.711_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > 5.3001_wp .AND. G%x(l,j,k,1) < 8.2999_wp) then

                      Vp = 6.18_wp   ! P-wave speed 
                      Vs = 3.65_wp   ! S-wave speed 
                      rho = 2.756_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 8.2999_wp .AND. G%x(l,j,k,1) <= 8.3001_wp) then

                      Vp = 6.20_wp ! P-wave speed
                      Vs  = 3.63_wp ! S-wave speed
                      rho  = 2.762_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > 8.3001_wp .AND. G%x(l,j,k,1) < 11.2999_wp) then

                      Vp = 6.23_wp   ! P-wave speed 
                      Vs = 3.60_wp   ! S-wave speed 
                      rho = 2.767_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 11.2999_wp .AND. G%x(l,j,k,1) <= 11.3001_wp) then

                      Vp = 6.26_wp ! P-wave speed
                      Vs  = 3.64_wp ! S-wave speed
                      rho  = 2.776_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > 11.3001_wp .AND. G%x(l,j,k,1) < 20.2999_wp) then

                      Vp = 6.30_wp   ! P-wave speed 
                      Vs = 3.67_wp   ! S-wave speed 
                      rho = 2.784_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 20.2999_wp .AND. G%x(l,j,k,1) <= 20.3001_wp) then

                      Vp = 6.54_wp ! P-wave speed
                      Vs  = 3.80_wp ! S-wave speed
                      rho  = 2.848_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > 20.3001_wp .AND. G%x(l,j,k,1) < 39.9999_wp) then

                      Vp = 6.80_wp   ! P-wave speed 
                      Vs = 3.93_wp   ! S-wave speed 
                      rho = 2.911_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 39.9999_wp .AND. G%x(l,j,k,1) <= 40.0001_wp) then

                      Vp = 7.39_wp ! P-wave speed
                      Vs  = 4.25_wp ! S-wave speed
                      rho  = 3.119_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > 40.0001_wp) then

                      Vp = 8.10_wp   ! P-wave speed 
                      Vs = 4.62_wp   ! S-wave speed 
                      rho = 3.326_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   end if

                end do
             end do
          end do


       end select
       do ii = 1,3
          call exchange_all_neighbors(G%C,M%M(:,:,:,ii))
       end do


    case('acoustic')

       call allocate_array_body(M%M,G%C,2, ghost_nodes=.true.)

       M%M(:,:,:,1) = 1.0_wp ! rho (density)
       M%M(:,:,:,2) = 1.0_wp ! K (bulk modulus)

    end select

  end subroutine init_material


  subroutine init_material_from_file(id,nb,M,G,material_path)

    use mpi3dbasic, only : pw, ps
    use mpi3dio, only: file_distributed, open_file_distributed, read_file_distributed, close_file_distributed
    use mpi3dcomm, only : allocate_array_body, exchange_all_neighbors 

    implicit none

    integer, intent(in) :: id, nb
    type(block_grid_t), intent(in) :: G
    type(block_material),intent(out) :: M
    type(file_distributed) :: fids(3)
    character(len=256),intent(in) :: material_path(3)
    character(len=256) :: name
    integer :: i,j,k


    ! allocate memory for array

    call allocate_array_body(M%M,G%C,3,ghost_nodes=.true.)
    M%M(:,:,:,:) = 0_wp

    do i = 1,3
        name = material_path(i)
        call open_file_distributed(fids(i),name,"read",G%C%comm,G%C%array_w,pw)
        call read_file_distributed(fids(i),M%M(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps,i))
        call close_file_distributed(fids(i))
    end do

    do i = 1,3
        call exchange_all_neighbors(G%C,M%M(:,:,:,i))
    end do

  !  do k = G%C%mq,G%C%pq
  !    print*, k
  !    print*,'lambda = ',M%M(k,G%C%mr,G%C%ms,1)
  !    print*,'mu = ' ,M%M(k,G%C%mr,G%C%ms,2)
  !    print*,'rho = ',M%M(k,G%C%mr,G%C%ms,3)
  !    print*,' '
  !  end do

  end subroutine init_material_from_file

end module material
