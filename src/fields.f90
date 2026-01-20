module fields

  !> fields modules handles fields defined within blocks

  use common, only : wp
  use datatypes, only : block_type, block_fields
  implicit none

  !> first few indices are for spatial dimensions, last index is for field component

contains


  subroutine init_fields(F, C, physics)

!!! initialize fields in a block

    use mpi3dcomm

    implicit none

    type(block_fields),intent(out) :: F
    type(cartesian3d_t),intent(in) :: C
    character(*),intent(in) :: physics

!!! fields will be specific to block physics

    select case(physics)

    case default

       stop 'invalid block physics in init_fields'

    case('elastic')

       call allocate_array_body(F%F ,C,9, ghost_nodes=.true., Fval = 0.0_wp)
       call allocate_array_body(F%DF,C,9, ghost_nodes=.true., Fval = 1.0e40_wp)

    case('acoustic')

       call allocate_array_body(F%F , C, 4, ghost_nodes=.true.)
       call allocate_array_body(F%DF, C, 4, ghost_nodes=.true.)

    end select

  end subroutine init_fields


  subroutine scale_rates_interior(F,A)

!!! multiply rates by RK coefficient A

    implicit none

    type(block_type),intent(inout) :: F
    real(kind = wp),intent(in) :: A
    

    F%F%DF = A*F%F%DF

      if (F%M%anelastic) then
         F%M%Deta4 = A*F%M%Deta4
         F%M%Deta5 = A*F%M%Deta5
         F%M%Deta6 = A*F%M%Deta6
         F%M%Deta7 = A*F%M%Deta7
         F%M%Deta8 = A*F%M%Deta8
         F%M%Deta9 = A*F%M%Deta9
      end if
    
    if( F%PMLB(1)%pml .EQV. .TRUE.) then
       if(F%G%C%mq .le.  F%PMLB(1)%N_pml) then
          F%PMLB(1)%DQ = A*F%PMLB(1)%DQ
       end if
    end if

    
    if( F%PMLB(2)%pml .EQV. .TRUE.) then
        if(F%G%C%pq .ge. (F%G%C%nq- F%PMLB(2)%N_pml+1)) then
           F%PMLB(2)%DQ = A*F%PMLB(2)%DQ
        end if
     end if
     
     if( F%PMLB(3)%pml .EQV. .TRUE.) then
        if(F%G%C%mr .le. F%PMLB(3)%N_pml) then
           F%PMLB(3)%DQ = A*F%PMLB(3)%DQ
        end if
     end if
     
     if( F%PMLB(4)%pml .EQV. .TRUE.) then
         if(F%G%C%pr .ge. (F%G%C%nr- F%PMLB(4)%N_pml+1)) then
            F%PMLB(4)%DQ = A*F%PMLB(4)%DQ
         end if
     end if
     
     if( F%PMLB(5)%pml .EQV. .TRUE.) then
        if(F%G%C%ms .le.  F%PMLB(5)%N_pml) then
           F%PMLB(5)%DQ = A*F%PMLB(5)%DQ
        end if
     end if
     
     if( F%PMLB(6)%pml .EQV. .TRUE.) then
         if(F%G%C%ps .ge. (F%G%C%ns- F%PMLB(6)%N_pml+1)) then
            F%PMLB(6)%DQ = A*F%PMLB(6)%DQ
         end if
     end if
     

  end subroutine scale_rates_interior


  subroutine update_fields_interior(F,dt)

!!! update fields using rates

    implicit none

    type(block_type),intent(inout) :: F
    real(kind = wp),intent(in) :: dt

    F%F%F = F%F%F + dt*F%F%DF

      if (F%M%anelastic) then
         F%M%eta4 = F%M%eta4 + dt*F%M%Deta4
         F%M%eta5 = F%M%eta5 + dt*F%M%Deta5
         F%M%eta6 = F%M%eta6 + dt*F%M%Deta6
         F%M%eta7 = F%M%eta7 + dt*F%M%Deta7
         F%M%eta8 = F%M%eta8 + dt*F%M%Deta8
         F%M%eta9 = F%M%eta9 + dt*F%M%Deta9
      end if

    if( F%PMLB(1)%pml .EQV. .TRUE.) then
       if(F%G%C%mq .le. F%PMLB(1)%N_pml) then
         F%PMLB(1)%Q = F%PMLB(1)%Q + dt*F%PMLB(1)%DQ
       end if
    end if

    
    if( F%PMLB(2)%pml .EQV. .TRUE.) then
        if(F%G%C%pq .ge. (F%G%C%nq- F%PMLB(2)%N_pml+1)) then
          F%PMLB(2)%Q = F%PMLB(2)%Q + dt*F%PMLB(2)%DQ  
        end if
     end if
     
     if( F%PMLB(3)%pml .EQV. .TRUE.) then
        if(F%G%C%mr .le.  F%PMLB(3)%N_pml) then
           F%PMLB(3)%Q = F%PMLB(3)%Q + dt*F%PMLB(3)%DQ
        end if
     end if
     
     if( F%PMLB(4)%pml .EQV. .TRUE.) then
         if(F%G%C%pr .ge. (F%G%C%nr- F%PMLB(4)%N_pml+1)) then
            F%PMLB(4)%Q = F%PMLB(4)%Q + dt*F%PMLB(4)%DQ 
         end if
     end if
     
     if( F%PMLB(5)%pml .EQV. .TRUE.) then
        if(F%G%C%ms .le.  F%PMLB(5)%N_pml) then
           F%PMLB(5)%Q = F%PMLB(5)%Q + dt*F%PMLB(5)%DQ
        end if
     end if
     
     if( F%PMLB(6)%pml .EQV. .TRUE.) then
         if(F%G%C%ps .ge. (F%G%C%ns- F%PMLB(6)%N_pml+1)) then
            F%PMLB(6)%Q = F%PMLB(6)%Q + dt*F%PMLB(6)%DQ
         end if
     end if
     
    
         

  end subroutine update_fields_interior


  subroutine exchange_all_fields(F, C)

    use mpi3dcomm
    use mpi3dbasic, only: rank

    implicit none

    type(block_fields), intent(inout) :: F
    type(cartesian3d_t),intent(in) :: C
    integer :: i

    do i = 1, size(F%F, 4)
       call exchange_all_neighbors(C, F%F(:,:,:,i))
    end do

  end subroutine exchange_all_fields
!!! note that there is no set_rates_interior routine in this module;
!!! such a routine, which is quite involved, is stored separately in other modules,
!!! like elastic in elastic.f90, with one for each physics

  subroutine norm_fields(F,mx,my,mz,px, py, pz, n,sum)

    implicit none

    integer, intent(in) :: mx, my, mz, px, py, pz, n
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(in) :: F
    real(kind = wp), intent(inout) :: sum

    integer :: i,j,k,l

    sum = 0.0_wp

    do l = 1,n
       do k = mz, pz
          do j = my, py
             do i = mx, px
                sum = sum + F(i,j,k,l)**2
             end do
          end do
       end do
    end do

    sum = sqrt(sum)

  end subroutine norm_fields


end module fields
