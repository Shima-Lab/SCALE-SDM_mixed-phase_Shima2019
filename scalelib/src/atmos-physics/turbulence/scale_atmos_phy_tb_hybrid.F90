!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Boundary layer turbulence model
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-09-18 (S.Nishizawa) [new]
!!
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_tb_hybrid
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

#ifdef DEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_TB_hybrid_setup
  public :: ATMOS_PHY_TB_hybrid

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
  abstract interface
     subroutine tb( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_rhoq,                & ! (out)
       tke,                                         & ! (inout)
       tke_t, nu_C, Ri, Pr, N2,                     & ! (out) diagnostic variables
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,          & ! (in)
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV, & ! (in)
       GSQRT, J13G, J23G, J33G, MAPF, dt            ) ! (in)
       use scale_precision
       use scale_grid_index
       use scale_tracer
       implicit none
       real(RP), intent(out) :: qflx_sgs_momz(KA,IA,JA,3)
       real(RP), intent(out) :: qflx_sgs_momx(KA,IA,JA,3)
       real(RP), intent(out) :: qflx_sgs_momy(KA,IA,JA,3)
       real(RP), intent(out) :: qflx_sgs_rhot(KA,IA,JA,3)
       real(RP), intent(out) :: qflx_sgs_rhoq(KA,IA,JA,3,QA)

       real(RP), intent(inout) :: tke (KA,IA,JA) ! TKE
       real(RP), intent(out)   :: tke_t(KA,IA,JA) ! TKE
       real(RP), intent(out)   :: nu_C(KA,IA,JA) ! eddy viscosity (center)
       real(RP), intent(out)   :: Ri  (KA,IA,JA) ! Richardson number
       real(RP), intent(out)   :: Pr  (KA,IA,JA) ! Prantle number
       real(RP), intent(out)   :: N2  (KA,IA,JA) ! squared Brunt-Vaisala frequency

       real(RP), intent(in)  :: MOMZ(KA,IA,JA)
       real(RP), intent(in)  :: MOMX(KA,IA,JA)
       real(RP), intent(in)  :: MOMY(KA,IA,JA)
       real(RP), intent(in)  :: RHOT(KA,IA,JA)
       real(RP), intent(in)  :: DENS(KA,IA,JA)
       real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)

       real(RP), intent(in)  :: SFLX_MW(IA,JA)
       real(RP), intent(in)  :: SFLX_MU(IA,JA)
       real(RP), intent(in)  :: SFLX_MV(IA,JA)
       real(RP), intent(in)  :: SFLX_SH(IA,JA)
       real(RP), intent(in)  :: SFLX_QV(IA,JA)

       real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
       real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
       real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
       real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
       real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor
       real(DP), intent(in)  :: dt
     end subroutine tb
  end interface
  procedure(tb), pointer :: SGS_TB => NULL()
  procedure(tb), pointer :: PBL_TB => NULL()
  real(RP), allocatable :: frac_sgs(:,:)
  real(RP), allocatable :: frac_pbl(:,:)
  real(RP), allocatable :: frac_sgs_tke(:,:)
  real(RP), allocatable :: frac_pbl_tke(:,:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_hybrid_setup( &
       TB_TYPE,       &
       CDZ, CDX, CDY, &
       CZ             )
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_phy_tb_smg, only: &
       ATMOS_PHY_TB_smg_setup, &
       ATMOS_PHY_TB_smg
    use scale_atmos_phy_tb_mynn, only: &
       ATMOS_PHY_TB_mynn_setup, &
       ATMOS_PHY_TB_mynn
    implicit none

    character(len=*), intent(in) :: TB_TYPE

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: CZ (KA,IA,JA)

    real(RP) :: ATMOS_PHY_TB_HYBRID_SGS_DX = 100.0_RP                      !< horizontal resolution for SGS model
    real(RP) :: ATMOS_PHY_TB_HYBRID_PBL_DX = 500.0_RP                      !< horizontal resolution for PBL model
    character(len=H_SHORT) :: ATMOS_PHY_TB_HYBRID_SGS_TYPE = 'SMAGORINSKY' !< scheme type for SGS
    character(len=H_SHORT) :: ATMOS_PHY_TB_HYBRID_PBL_TYPE = 'MYNN'        !< scheme type for turbulent parametarization

    character(len=H_SHORT) :: ATMOS_PHY_TB_HYBRID_TKE_TYPE = 'PBL'         !< SGS, PBL, or MIXED

    NAMELIST / PARAM_ATMOS_PHY_TB_HYBRID / &
         ATMOS_PHY_TB_HYBRID_SGS_DX, &
         ATMOS_PHY_TB_HYBRID_PBL_DX, &
         ATMOS_PHY_TB_HYBRID_SGS_TYPE, &
         ATMOS_PHY_TB_HYBRID_PBL_TYPE, &
         ATMOS_PHY_TB_HYBRID_TKE_TYPE


    real(RP) :: dxy

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[TURBULENCE] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ SGS-parameterization hybrid Model'

    if ( TB_TYPE /= 'HYBRID' ) then
       write(*,*) 'xxx ATMOS_PHY_TB_TYPE is not HYBRID. Check!'
       call PRC_MPIstop
    endif


    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_HYBRID,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_TB_HYBRID. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_TB_HYBRID)

    select case ( ATMOS_PHY_TB_HYBRID_SGS_TYPE )
    case ('SMAGORINSKY')
       call ATMOS_PHY_TB_SMG_setup( &
            ATMOS_PHY_TB_HYBRID_SGS_TYPE, &
            CDZ, CDX, CDY, &
            CZ             )
       SGS_TB => ATMOS_PHY_TB_SMG
    case default
       write(*,*) 'xxx ATMOS_PHY_TB_HYBRID_SGS_TYPE is invalid'
       call PRC_MPIstop
    end select

    select case ( ATMOS_PHY_TB_HYBRID_PBL_TYPE )
    case ('MYNN')
       call ATMOS_PHY_TB_mynn_setup( &
            ATMOS_PHY_TB_HYBRID_PBL_TYPE, &
            CDZ, CDX, CDY, &
            CZ             )
       PBL_TB => ATMOS_PHY_TB_mynn
    case default
       write(*,*) 'xxx ATMOS_PHY_TB_HYBRID_PBL_TYPE is invalid'
       call PRC_MPIstop
    end select

    allocate( frac_sgs(IA,JA) )
    allocate( frac_pbl(IA,JA) )
    allocate( frac_sgs_tke(IA,JA) )
    allocate( frac_pbl_tke(IA,JA) )

    do j = 1, JA
    do i = 1, IA
       dxy = sqrt( ( CDX(i)**2 + CDY(j)**2 )*0.5_RP )
       frac_pbl(i,j) = &
               min( 1.0_RP, &
               max( 0.0_RP, &
                    ( dxy - ATMOS_PHY_TB_HYBRID_SGS_DX ) &
                  / ( ATMOS_PHY_TB_HYBRID_PBL_DX - ATMOS_PHY_TB_HYBRID_SGS_DX ) ) )
       frac_sgs(i,j) = 1.0_RP - frac_pbl(i,j)
    end do
    end do

    select case ( ATMOS_PHY_TB_HYBRID_TKE_TYPE )
    case ('SGS')
       frac_pbl_tke(:,:) = 0.0_RP
       frac_sgs_tke(:,:) = 1.0_RP
    case ('PBL')
       frac_pbl_tke(:,:) = 1.0_RP
       frac_sgs_tke(:,:) = 0.0_RP
    case ('MIXED')
       frac_pbl_tke(:,:) = frac_pbl(:,:)
       frac_sgs_tke(:,:) = frac_sgs(:,:)
    case default
       write(*,*) 'xxx ATMOS_PHY_TB_HYBRID_TKE_TYPE is invalid'
       call PRC_MPIstop
    end select

    return
  end subroutine ATMOS_PHY_TB_hybrid_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_hybrid( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_rhoq,                & ! (out)
       tke,                                         & ! (inout)
       tke_t, Nu, Ri, Pr, N2,                       & ! (out) diagnostic variables
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,          & ! (in)
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV, & ! (in)
       GSQRT, J13G, J23G, J33G, MAPF, dt            ) ! (in)
    use scale_precision
    use scale_grid_index
    use scale_tracer
    use scale_const, only: &
       GRAV   => CONST_GRAV, &
       EPS    => CONST_EPS
    use scale_grid, only: &
       CDZ  => GRID_CDZ,  &
       FDZ  => GRID_FDZ,  &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ
    implicit none

    real(RP), intent(out) :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhoq(KA,IA,JA,3,QA)

    real(RP), intent(inout) :: tke (KA,IA,JA) ! TKE
    real(RP), intent(out) :: tke_t(KA,IA,JA) ! tendency TKE
    real(RP), intent(out) :: Nu(KA,IA,JA) ! eddy viscosity (center)
    real(RP), intent(out) :: Pr(KA,IA,JA) ! Plandtle number
    real(RP), intent(out) :: Ri(KA,IA,JA) ! Richardson number
    real(RP), intent(out) :: N2(KA,IA,JA) ! squared Brunt-Vaisala frequency

    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)

    real(RP), intent(in)  :: SFLX_MW(IA,JA)
    real(RP), intent(in)  :: SFLX_MU(IA,JA)
    real(RP), intent(in)  :: SFLX_MV(IA,JA)
    real(RP), intent(in)  :: SFLX_SH(IA,JA)
    real(RP), intent(in)  :: SFLX_QV(IA,JA)

    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor
    real(DP), intent(in)  :: dt

    real(RP) :: w_qflx_sgs_momz(KA,IA,JA,3,2)
    real(RP) :: w_qflx_sgs_momx(KA,IA,JA,3,2)
    real(RP) :: w_qflx_sgs_momy(KA,IA,JA,3,2)
    real(RP) :: w_qflx_sgs_rhot(KA,IA,JA,3,2)
    real(RP) :: w_qflx_sgs_rhoq(KA,IA,JA,3,QA,2)

    real(RP) :: w_tke(KA,IA,JA,2)
    real(RP) :: w_tke_t(KA,IA,JA,2)
    real(RP) :: w_Nu(KA,IA,JA,2)
    real(RP) :: w_Ri(KA,IA,JA,2)
    real(RP) :: w_Pr(KA,IA,JA,2)
    real(RP) :: w_N2(KA,IA,JA,2)

    integer :: k, i, j, iq

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       w_tke(k,i,j,1) = tke(k,i,j)
       w_tke(k,i,j,2) = tke(k,i,j)
    end do
    end do
    end do

    call SGS_TB( &
         w_qflx_sgs_momz(:,:,:,:,1), w_qflx_sgs_momx(:,:,:,:,1), & ! (out)
         w_qflx_sgs_momy(:,:,:,:,1), w_qflx_sgs_rhot(:,:,:,:,1), & ! (out
         w_qflx_sgs_rhoq(:,:,:,:,:,1),                           & ! (out)
         w_tke(:,:,:,1),                                         & ! (inout)
         w_tke_t(:,:,:,1),                                       & ! (out)
         w_Nu(:,:,:,1), w_Ri(:,:,:,1), w_Pr(:,:,:,1), w_N2(:,:,:,1), & ! (out)
         MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,                     & ! (in)
         SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV,            & ! (in)
         GSQRT, J13G, J23G, J33G, MAPF, dt                       ) ! (in)

    call PBL_TB( &
         w_qflx_sgs_momz(:,:,:,:,2), w_qflx_sgs_momx(:,:,:,:,2), & ! (out)
         w_qflx_sgs_momy(:,:,:,:,2), w_qflx_sgs_rhot(:,:,:,:,2), & ! (out
         w_qflx_sgs_rhoq(:,:,:,:,:,2),                           & ! (out)
         w_tke(:,:,:,2),                                         & ! (inout)
         w_tke_t(:,:,:,2),                                       & ! (out)
         w_Nu(:,:,:,2), w_Ri(:,:,:,2), w_Pr(:,:,:,2), w_N2(:,:,:,2), & ! (out)
         MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,                     & ! (in)
         SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV,            & ! (in)
         GSQRT, J13G, J23G, J33G, MAPF, dt                       ) ! (in)

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       qflx_sgs_momz(k,i,j,ZDIR) = w_qflx_sgs_momz(k,i,j,ZDIR,1) * frac_sgs(i,j) &
                                 + w_qflx_sgs_momz(k,i,j,ZDIR,2) * frac_pbl(i,j)
       qflx_sgs_momz(k,i,j,XDIR) = w_qflx_sgs_momz(k,i,j,XDIR,1)
       qflx_sgs_momz(k,i,j,YDIR) = w_qflx_sgs_momz(k,i,j,YDIR,1)
    end do
    end do
    end do

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       qflx_sgs_momx(k,i,j,ZDIR) = w_qflx_sgs_momx(k,i,j,ZDIR,1) * frac_sgs(i,j) &
                                 + w_qflx_sgs_momx(k,i,j,ZDIR,2) * frac_pbl(i,j)
       qflx_sgs_momx(k,i,j,XDIR) = w_qflx_sgs_momx(k,i,j,XDIR,1)
       qflx_sgs_momx(k,i,j,YDIR) = w_qflx_sgs_momx(k,i,j,YDIR,1)
    end do
    end do
    end do

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       qflx_sgs_momy(k,i,j,ZDIR) = w_qflx_sgs_momy(k,i,j,ZDIR,1) * frac_sgs(i,j) &
                                 + w_qflx_sgs_momy(k,i,j,ZDIR,2) * frac_pbl(i,j)
       qflx_sgs_momy(k,i,j,XDIR) = w_qflx_sgs_momy(k,i,j,XDIR,1)
       qflx_sgs_momy(k,i,j,YDIR) = w_qflx_sgs_momy(k,i,j,YDIR,1)
    end do
    end do
    end do

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       qflx_sgs_rhot(k,i,j,ZDIR) = w_qflx_sgs_rhot(k,i,j,ZDIR,1) * frac_sgs(i,j) &
                                 + w_qflx_sgs_rhot(k,i,j,ZDIR,2) * frac_pbl(i,j)
       qflx_sgs_rhot(k,i,j,XDIR) = w_qflx_sgs_rhot(k,i,j,XDIR,1)
       qflx_sgs_rhot(k,i,j,YDIR) = w_qflx_sgs_rhot(k,i,j,YDIR,1)
    end do
    end do
    end do

    do iq = 1, QA
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       qflx_sgs_rhoq(k,i,j,ZDIR,iq) = w_qflx_sgs_rhoq(k,i,j,ZDIR,iq,1) * frac_sgs(i,j) &
                                    + w_qflx_sgs_rhoq(k,i,j,ZDIR,iq,2) * frac_pbl(i,j)
       qflx_sgs_rhoq(k,i,j,XDIR,iq) = w_qflx_sgs_rhoq(k,i,j,XDIR,iq,1)
       qflx_sgs_rhoq(k,i,j,YDIR,iq) = w_qflx_sgs_rhoq(k,i,j,YDIR,iq,1)
    end do
    end do
    end do
    end do

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       tke(k,i,j) = w_tke(k,i,j,1) * frac_sgs_tke(i,j) &
                  + w_tke(k,i,j,2) * frac_pbl_tke(i,j)
    end do
    end do
    end do

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       tke_t(k,i,j) = w_tke_t(k,i,j,1) * frac_sgs_tke(i,j) &
                    + w_tke_t(k,i,j,2) * frac_pbl_tke(i,j)
    end do
    end do
    end do

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       Nu(k,i,j) = w_Nu(k,i,j,1) * frac_sgs(i,j) &
                 + w_Nu(k,i,j,2) * frac_pbl(i,j)
    end do
    end do
    end do

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       Ri(k,i,j) = w_Ri(k,i,j,1) * frac_sgs(i,j) &
                 + w_Ri(k,i,j,2) * frac_pbl(i,j)
    end do
    end do
    end do

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       Pr(k,i,j) = w_Pr(k,i,j,1) * frac_sgs(i,j) &
                 + w_Pr(k,i,j,2) * frac_pbl(i,j)
    end do
    end do
    end do

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       N2(k,i,j) = w_N2(k,i,j,1) * frac_sgs(i,j) &
                 + w_N2(k,i,j,2) * frac_pbl(i,j)
    end do
    end do
    end do

    return
  end subroutine ATMOS_PHY_TB_hybrid

end module scale_atmos_phy_tb_hybrid
