!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-20 (S.Nishizawa)   [new] split from dynamical core
!! @li      2018-02-12 (S.Shima)       [mod] for heating of Khain et al. (2004)
!!
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  use scale_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_setup0
  public :: USER_setup
  public :: USER_step
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, save        :: HEAT_TIME = 0.0_RP      ! heating duration [s] 
  real(RP), private, save        :: HEAT_RATE = 0.0_RP      ! heating rate [K/h] 
  real(RP), private, save        :: HEAT_CX   = 0.0_RP      ! x coordinate of the center of the heating domain [m] 
  real(RP), private, save        :: HEAT_WX   = 5000.0_RP   ! width along x-axis of the heating domain [m] 
  real(RP), private, save, allocatable :: HEAT_WEIGHT(:,:,:)      ! weight of heating HEAT_AREA(KA,IA,JA)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup0
  subroutine USER_setup0
  end subroutine USER_setup0
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
         PRC_MPIstop
    use scale_les_process, only: &
         PRC_NUM_X
    use scale_grid, only: &
         GRID_CZ,  &
         GRID_CX,  &
         GRID_CXG

    namelist / PARAM_USER / &
       HEAT_TIME,HEAT_RATE,HEAT_CX,HEAT_WX

    real(RP) :: xaxis_weight, zaxis_weight
    integer :: i,j,k
    integer :: ierr
    !---------------------------------------------------------------------------
    !--- read namelist
    HEAT_CX = (GRID_CXG(IMAX*PRC_NUM_X)-GRID_CXG(1))/2.0_RP + GRID_CXG(1) ! set the default position of the heating area as the domain center 
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)

    ! allocate variables
    allocate( HEAT_WEIGHT(KA,IA,JA) )
    
    ! define the heating area
    HEAT_WEIGHT(:,:,:) = 0.0_RP    
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       xaxis_weight = (-4.0_RP/(HEAT_WX**2)) *                &
            & (GRID_CX(i)-(HEAT_CX-HEAT_WX/2.0_RP))*(GRID_CX(i)-(HEAT_CX+HEAT_WX/2.0_RP)) 
       xaxis_weight = max(0.0_RP,xaxis_weight)

       zaxis_weight = exp(-2.0_RP*(GRID_CZ(k)/1000.0_RP - 0.5_RP))
       
       HEAT_WEIGHT(k,i,j) = xaxis_weight * zaxis_weight 
    enddo
    enddo
    enddo

    call USER_step

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> User step
  !! This subroutine is called every TIME_DT 
  subroutine USER_step
    use mod_atmos_vars, only: &
         DENS,    &
         RHOT_tp
    use scale_time, only: &
         TIME_NOWSEC

    implicit none
    !---------------------------------------------------------------------------
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** USER_step'

    ! Add tendency of RHOT due to the heating
    if( TIME_NOWSEC <= HEAT_TIME)then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + DENS(k,i,j)*(HEAT_RATE/3600.0_RP)*HEAT_WEIGHT(k,i,j)
       end do
       end do
       end do
    end if

    return
  end subroutine USER_step

end module mod_user
