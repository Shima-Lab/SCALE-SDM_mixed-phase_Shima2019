!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics - Common
!!
!! @par Description
!!          Common module for Cloud Microphysics
!!          Sedimentation/Precipitation and Saturation adjustment
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-23 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_mp_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  use scale_const, only: &
     UNDEF => CONST_UNDEF
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_negative_fixer
  public :: ATMOS_PHY_MP_saturation_adjustment
  public :: ATMOS_PHY_MP_precipitation

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: moist_conversion_liq
  private :: moist_conversion_all

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Negative fixer
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_negative_fixer( &
       DENS, &
       RHOT, &
       QTRC  )
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)

    real(RP) :: diffq

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_filter', 3)

    !$omp parallel do private(i,j,diffq) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE

       diffq = 0.0_RP
       do iq = QQS+1, QQE
          ! total hydrometeor (before correction)
          diffq = diffq + QTRC(k,i,j,iq)
          ! remove negative value of hydrometeors (mass)
          QTRC(k,i,j,iq) = max( QTRC(k,i,j,iq), 0.0_RP )
       enddo

       do iq = QQS+1, QQE
          ! difference between before and after correction
          diffq = diffq - QTRC(k,i,j,iq)
       enddo

       ! Compensate for the lack of hydrometeors by the water vapor
       QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) + diffq

       ! TODO: We have to consider energy conservation (but very small)

       ! remove negative value of water vapor (mass)
       diffq = QTRC(k,i,j,I_QV)
       QTRC(k,i,j,I_QV) = max( QTRC(k,i,j,I_QV), 0.0_RP )
       diffq = diffq - QTRC(k,i,j,I_QV)

       ! Apply correction to total density
       ! TODO: We have to consider energy conservation (but very small)
       DENS(k,i,j) = DENS(k,i,j) * ( 1.0_RP - diffq ) ! diffq is negative
       RHOT(k,i,j) = RHOT(k,i,j) * ( 1.0_RP - diffq )

    enddo
    enddo
    enddo

    call PROF_rapend('MP_filter', 3)

    return
  end subroutine ATMOS_PHY_MP_negative_fixer

  !-----------------------------------------------------------------------------
  !> Saturation adjustment
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_saturation_adjustment( &
       RHOE_t, &
       QTRC_t, &
       RHOE0,  &
       QTRC0,  &
       DENS0   )
    use scale_const, only: &
       LHV => CONST_LHV, &
       LHF => CONST_LHF
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd          => ATMOS_THERMODYN_qd,         &
       THERMODYN_cv          => ATMOS_THERMODYN_cv,         &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       SATURATION_dens2qsat_all => ATMOS_SATURATION_dens2qsat_all
    implicit none

    real(RP), intent(inout) :: RHOE_t(KA,IA,JA)    ! tendency rhoe             [J/m3/s]
    real(RP), intent(inout) :: QTRC_t(KA,IA,JA,QA) ! tendency tracer           [kg/kg/s]
    real(RP), intent(inout) :: RHOE0 (KA,IA,JA)    ! density * internal energy [J/m3]
    real(RP), intent(inout) :: QTRC0 (KA,IA,JA,QA) ! mass concentration        [kg/kg]
    real(RP), intent(in)    :: DENS0 (KA,IA,JA)    ! density                   [kg/m3]

    ! working
    real(RP) :: TEMP0 (KA,IA,JA)
    real(RP) :: PRES0 (KA,IA,JA)
    real(RP) :: QDRY0 (KA,IA,JA)
    real(RP) :: CVtot (KA,IA,JA)

    real(RP) :: Emoist(KA,IA,JA) ! moist internal energy
    real(RP) :: QSUM1 (KA,IA,JA) ! QV+QC+QI
    real(RP) :: TEMP1 (KA,IA,JA)

    real(RP) :: RHOE1 (KA,IA,JA)
    real(RP) :: QTRC1 (KA,IA,JA,QA)
    real(RP) :: rdt

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

#ifndef DRY

    call PROF_rapstart('MP_Saturation_adjustment', 3)

    rdt = 1.0_RP / dt

    !$omp parallel do private(i,j,k,iq) OMP_SCHEDULE_ collapse(4)
    do iq = QQS, QQE
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       QTRC1(k,i,j,iq) = QTRC0(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    call THERMODYN_temp_pres_E( TEMP0(:,:,:),  & ! [OUT]
                                PRES0(:,:,:),  & ! [OUT]
                                DENS0(:,:,:),  & ! [IN]
                                RHOE0(:,:,:),  & ! [IN]
                                QTRC0(:,:,:,:) ) ! [IN]

    ! qdry dont change through the process
    call THERMODYN_qd( QDRY0(:,:,:),   & ! [OUT]
                       QTRC0(:,:,:,:)  ) ! [IN]

    call THERMODYN_cv( CVtot(:,:,:),   & ! [OUT]
                       QTRC0(:,:,:,:), & ! [IN]
                       QDRY0(:,:,:)    ) ! [IN]

    if ( I_QI <= 0 ) then ! warm rain

       ! Turn QC into QV with consistency of moist internal energy
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          Emoist(k,i,j) = TEMP0(k,i,j) * CVtot(k,i,j) &
                        + QTRC1(k,i,j,I_QV) * LHV

          QSUM1(k,i,j) = QTRC1(k,i,j,I_QV) &
                       + QTRC1(k,i,j,I_QC)

          QTRC1(k,i,j,I_QV) = QSUM1(k,i,j)
          QTRC1(k,i,j,I_QC) = 0.0_RP
       enddo
       enddo
       enddo

       call THERMODYN_cv( CVtot(:,:,:),   & ! [OUT]
                          QTRC1(:,:,:,:), & ! [IN]
                          QDRY0(:,:,:)    ) ! [IN]

       ! new temperature (after QC evaporation)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          TEMP1(k,i,j) = ( Emoist(k,i,j) - QTRC1(k,i,j,I_QV) * LHV ) / CVtot(k,i,j)
       enddo
       enddo
       enddo

       call moist_conversion_liq( TEMP1 (:,:,:),   & ! [INOUT]
                                  QTRC1 (:,:,:,:), & ! [INOUT]
                                  DENS0 (:,:,:),   & ! [IN]
                                  QSUM1 (:,:,:),   & ! [IN]
                                  QDRY0 (:,:,:),   & ! [IN]
                                  Emoist(:,:,:)    ) ! [IN]

    else ! cold rain

       ! Turn QC & QI into QV with consistency of moist internal energy
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          Emoist(k,i,j) = TEMP0(k,i,j) * CVtot(k,i,j) &
                        + QTRC1(k,i,j,I_QV) * LHV     &
                        - QTRC1(k,i,j,I_QI) * LHF

          QSUM1(k,i,j) = QTRC1(k,i,j,I_QV) &
                       + QTRC1(k,i,j,I_QC) &
                       + QTRC1(k,i,j,I_QI)

          QTRC1(k,i,j,I_QV) = QSUM1(k,i,j)
          QTRC1(k,i,j,I_QC) = 0.0_RP
          QTRC1(k,i,j,I_QI) = 0.0_RP
       enddo
       enddo
       enddo

       call THERMODYN_cv( CVtot(:,:,:),   & ! [OUT]
                          QTRC1(:,:,:,:), & ! [IN]
                          QDRY0(:,:,:)    ) ! [IN]

       ! new temperature (after QC & QI evaporation)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          TEMP1(k,i,j) = ( Emoist(k,i,j) - QTRC1(k,i,j,I_QV) * LHV ) / CVtot(k,i,j)
       enddo
       enddo
       enddo

       call moist_conversion_all( TEMP1 (:,:,:),   & ! [INOUT]
                                  QTRC1 (:,:,:,:), & ! [INOUT]
                                  DENS0 (:,:,:),   & ! [IN]
                                  QSUM1 (:,:,:),   & ! [IN]
                                  QDRY0 (:,:,:),   & ! [IN]
                                  Emoist(:,:,:)    ) ! [IN]

    endif

    call THERMODYN_cv( CVtot(:,:,:),   & ! [OUT]
                       QTRC1(:,:,:,:), & ! [IN]
                       QDRY0(:,:,:)    ) ! [IN]

    ! mass & energy update
    !$omp parallel do private(i,j,k,iq) OMP_SCHEDULE_ collapse(4)
    do iq = QQS, QQE
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC_t(k,i,j,iq) = QTRC_t(k,i,j,iq) + ( QTRC1(k,i,j,iq) - QTRC0(k,i,j,iq) ) * rdt

       QTRC0(k,i,j,iq) = QTRC1(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOE1(k,i,j) = DENS0(k,i,j) * TEMP1(k,i,j) * CVtot(k,i,j)

       RHOE_t(k,i,j) = RHOE_t(k,i,j) + ( RHOE1(k,i,j) - RHOE0(k,i,j) ) * rdt

       RHOE0(k,i,j) = RHOE1(k,i,j)
    enddo
    enddo
    enddo

    call PROF_rapend  ('MP_Saturation_adjustment', 3)

#else
    RHOE_t = UNDEF
    QTRC_t = UNDEF
    RHOE0  = UNDEF
    QTRC0  = UNDEF
#endif
    return
  end subroutine ATMOS_PHY_MP_saturation_adjustment

  !-----------------------------------------------------------------------------
  !> Iterative moist conversion for warm rain
  !-----------------------------------------------------------------------------
  subroutine moist_conversion_liq( &
       TEMP1, &
       QTRC1, &
       DENS0, &
       QSUM1, &
       QDRY0, &
       Emoist )
    use scale_const, only: &
       LHV => CONST_LHV
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_thermodyn, only: &
       THERMODYN_cv => ATMOS_THERMODYN_cv, &
       CVw => AQ_CV
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       CVovR_liq, &
       LovR_liq
    implicit none

    real(RP), intent(inout) :: TEMP1 (KA,IA,JA)
    real(RP), intent(inout) :: QTRC1 (KA,IA,JA,QA)
    real(RP), intent(in)    :: DENS0 (KA,IA,JA)
    real(RP), intent(in)    :: QSUM1 (KA,IA,JA)
    real(RP), intent(in)    :: QDRY0 (KA,IA,JA)
    real(RP), intent(in)    :: Emoist(KA,IA,JA)

    real(RP) :: QSAT(KA,IA,JA) ! saturated water vapor

    ! working
    real(RP) :: temp
    real(RP) :: q(QA)
    real(RP) :: CVtot
    real(RP) :: qsatl_new
    real(RP) :: Emoist_new ! moist internal energy

    ! d(X)/dT
    real(RP) :: dqsatl_dT
    real(RP) :: dqc_dT
    real(RP) :: dCVtot_dT
    real(RP) :: dEmoist_dT
    real(RP) :: dtemp

    integer  :: ijk_sat
    integer  :: index_sat(KA*IA*JA,3) ! list vector

    integer, parameter :: itelim = 100
    real(RP) :: dtemp_criteria
    logical  :: converged
    integer  :: k, i, j, ijk, iq, ite
    !---------------------------------------------------------------------------

#ifndef DRY

    dtemp_criteria = 0.1_RP**(2+RP/2)

    call SATURATION_dens2qsat_liq( QSAT (:,:,:), & ! [OUT]
                                   TEMP1(:,:,:), & ! [IN]
                                   DENS0(:,:,:)  ) ! [IN]

    ijk_sat = 0
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( QSUM1(k,i,j) > QSAT(k,i,j) ) then
          !$omp critical
          ijk_sat = ijk_sat + 1
          index_sat(ijk_sat,1) = k
          index_sat(ijk_sat,2) = i
          index_sat(ijk_sat,3) = j
          !$omp end critical
       endif
    enddo
    enddo
    enddo

    do ijk = 1, ijk_sat
       k = index_sat(ijk,1)
       i = index_sat(ijk,2)
       j = index_sat(ijk,3)

       ! store to work
       temp = TEMP1(k,i,j)
       do iq = QQS, QQE
          q(iq) = QTRC1(k,i,j,iq)
       enddo

       converged = .false.
       do ite = 1, itelim

          call SATURATION_dens2qsat_liq( qsatl_new,   & ! [OUT]
                                         temp,        & ! [IN]
                                         DENS0(k,i,j) ) ! [IN]

          ! Separation
          q(I_QV) = qsatl_new
          q(I_QC) = QSUM1(k,i,j) - qsatl_new

          call THERMODYN_cv( CVtot,       & ! [OUT]
                             q(:),        & ! [IN]
                             QDRY0(k,i,j) ) ! [IN]

          Emoist_new = temp * CVtot + qsatl_new * LHV

          ! dX/dT
          dqsatl_dT = ( LovR_liq / ( temp*temp ) + CVovR_liq / temp ) * qsatl_new

          dqc_dT = - dqsatl_dT

          dCVtot_dT = dqsatl_dT * CVw(I_QV) &
                    + dqc_dT    * CVw(I_QC)

          dEmoist_dT = qsatl_new * dCVtot_dT + CVtot + dqsatl_dT * LHV

          dtemp = ( Emoist_new - Emoist(k,i,j) ) / dEmoist_dT
          temp  = temp - dtemp

          if ( abs(dtemp) < dtemp_criteria ) then
             converged = .true.
             exit
          endif

          if( temp*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) 'xxx [moist_conversion] not converged! dtemp=', dtemp,k,i,j,ite
          call PRC_MPIstop
       endif

       TEMP1(k,i,j) = temp
       do iq = QQS, QQE
          QTRC1(k,i,j,iq) = q(iq)
       enddo
    enddo

#else
    TEMP1 = UNDEF
    QTRC1 = UNDEF
#endif

    return
  end subroutine moist_conversion_liq

  !-----------------------------------------------------------------------------
  !> Iterative moist conversion (liquid/ice mixture)
  !-----------------------------------------------------------------------------
  subroutine moist_conversion_all( &
       TEMP1, &
       QTRC1, &
       DENS0, &
       QSUM1, &
       QDRY0, &
       Emoist )
    use scale_const, only: &
       LHV => CONST_LHV,  &
       LHF => CONST_LHF
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_thermodyn, only: &
       THERMODYN_cv => ATMOS_THERMODYN_cv, &
       CVw => AQ_CV
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_all => ATMOS_SATURATION_dens2qsat_all, &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq, &
       SATURATION_dens2qsat_ice => ATMOS_SATURATION_dens2qsat_ice, &
       SATURATION_alpha         => ATMOS_SATURATION_alpha,         &
       SATURATION_dalphadT      => ATMOS_SATURATION_dalphadT,      &
       CVovR_liq, &
       CVovR_ice, &
       LovR_liq,  &
       LovR_ice
    implicit none

    real(RP), intent(inout) :: TEMP1 (KA,IA,JA)
    real(RP), intent(inout) :: QTRC1 (KA,IA,JA,QA)
    real(RP), intent(in)    :: DENS0 (KA,IA,JA)
    real(RP), intent(in)    :: QSUM1 (KA,IA,JA)
    real(RP), intent(in)    :: QDRY0 (KA,IA,JA)
    real(RP), intent(in)    :: Emoist(KA,IA,JA)

    real(RP) :: QSAT(KA,IA,JA) ! saturated water vapor

    ! working
    real(RP) :: temp
    real(RP) :: q(QA)
    real(RP) :: CVtot
    real(RP) :: alpha
    real(RP) :: qsat_new, qsatl_new, qsati_new
    real(RP) :: Emoist_new ! moist internal energy

    ! d(X)/dT
    real(RP) :: dalpha_dT
    real(RP) :: dqsat_dT, dqsatl_dT, dqsati_dT
    real(RP) :: dqc_dT, dqi_dT
    real(RP) :: dCVtot_dT
    real(RP) :: dEmoist_dT
    real(RP) :: dtemp

    integer  :: ijk_sat
    integer  :: index_sat(KA*IA*JA,3) ! list vector

    integer, parameter :: itelim = 100
    real(RP) :: dtemp_criteria
    logical  :: converged
    integer  :: k, i, j, ijk, iq, ite
    !---------------------------------------------------------------------------

#ifndef DRY
    dtemp_criteria = 0.1_RP**(2+RP/2)

    call SATURATION_dens2qsat_all( QSAT (:,:,:), & ! [OUT]
                                   TEMP1(:,:,:), & ! [IN]
                                   DENS0(:,:,:)  ) ! [IN]

    ijk_sat = 0
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( QSUM1(k,i,j) > QSAT(k,i,j) ) then
          !$omp critical
          ijk_sat = ijk_sat + 1
          index_sat(ijk_sat,1) = k
          index_sat(ijk_sat,2) = i
          index_sat(ijk_sat,3) = j
          !$omp end critical
       endif
    enddo
    enddo
    enddo

    do ijk = 1, ijk_sat
       k = index_sat(ijk,1)
       i = index_sat(ijk,2)
       j = index_sat(ijk,3)

       ! store to work
       temp = TEMP1(k,i,j)
       do iq = QQS, QQE
          q(iq) = QTRC1(k,i,j,iq)
       enddo

       converged = .false.
       do ite = 1, itelim

          ! liquid/ice separation factor
          call SATURATION_alpha( alpha, temp )
          ! Saturation
          call SATURATION_dens2qsat_all( qsat_new,  temp, DENS0(k,i,j) )
          call SATURATION_dens2qsat_liq( qsatl_new, temp, DENS0(k,i,j) )
          call SATURATION_dens2qsat_ice( qsati_new, temp, DENS0(k,i,j) )

          ! Separation
          q(I_QV) = qsat_new
          q(I_QC) = ( QSUM1(k,i,j)-qsat_new ) * (        alpha )
          q(I_QI) = ( QSUM1(k,i,j)-qsat_new ) * ( 1.0_RP-alpha )

          call THERMODYN_cv( CVtot,       & ! [OUT]
                             q(:),        & ! [IN]
                             QDRY0(k,i,j) ) ! [IN]

          Emoist_new = temp * CVtot + qsat_new * LHV - q(I_QI) * LHF

          ! dX/dT
          call SATURATION_dalphadT( dalpha_dT, temp )

          dqsatl_dT = ( LovR_liq / ( temp*temp ) + CVovR_liq / temp ) * qsatl_new
          dqsati_dT = ( LovR_ice / ( temp*temp ) + CVovR_ice / temp ) * qsati_new

          dqsat_dT  = qsatl_new * dalpha_dT + dqsatl_dT * (        alpha ) &
                    - qsati_new * dalpha_dT + dqsati_dT * ( 1.0_RP-alpha )

          dqc_dT =  ( QSUM1(k,i,j)-qsat_new ) * dalpha_dT - dqsat_dT * (        alpha )
          dqi_dT = -( QSUM1(k,i,j)-qsat_new ) * dalpha_dT - dqsat_dT * ( 1.0_RP-alpha )

          dCVtot_dT = dqsat_dT * CVw(I_QV) &
                    + dqc_dT   * CVw(I_QC) &
                    + dqi_dT   * CVw(I_QI)

          dEmoist_dT = temp * dCVtot_dT + CVtot + dqsat_dT * LHV - dqi_dT * LHF

          dtemp = ( Emoist_new - Emoist(k,i,j) ) / dEmoist_dT
          temp  = temp - dtemp

          if ( abs(dtemp) < dtemp_criteria ) then
             converged = .true.
             exit
          endif

          if( temp*0.0_RP /= 0.0_RP) exit
       enddo

       if ( .NOT. converged ) then
          write(*,*) TEMP1(k,i,j), dens0(k,i,j),q,QTRC1(k,i,j,QQS:QQE)
          write(*,*) 'xxx [moist_conversion] not converged! dtemp=', dtemp, k,i,j,ite
          call PRC_MPIstop
       endif

       TEMP1(k,i,j) = temp
       do iq = QQS, QQE
          QTRC1(k,i,j,iq) = q(iq)
       enddo
    enddo

#else
    TEMP1 = UNDEF
    QTRC1 = UNDEF
#endif

    return
  end subroutine moist_conversion_all

  !-----------------------------------------------------------------------------
  !> precipitation transport
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_precipitation( &
       flux_rain, &
       flux_snow, &
       DENS,      &
       MOMZ,      &
       MOMX,      &
       MOMY,      &
       RHOE,      &
       QTRC,      &
       vterm,     &
       temp,      &
       dt         )
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_FZ
    use scale_gridtrans, only: &
       J33G => GTRANS_J33G
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_thermodyn, only: &
       CVw => AQ_CV
    implicit none

    real(RP), intent(out)   :: flux_rain(KA,IA,JA)
    real(RP), intent(out)   :: flux_snow(KA,IA,JA)
    real(RP), intent(inout) :: DENS     (KA,IA,JA)
    real(RP), intent(inout) :: MOMZ     (KA,IA,JA)
    real(RP), intent(inout) :: MOMX     (KA,IA,JA)
    real(RP), intent(inout) :: MOMY     (KA,IA,JA)
    real(RP), intent(inout) :: RHOE     (KA,IA,JA)
    real(RP), intent(inout) :: QTRC     (KA,IA,JA,QA)
    real(RP), intent(inout) :: vterm    (KA,IA,JA,QA) ! terminal velocity of cloud mass
    real(RP), intent(in)    :: temp     (KA,IA,JA)
    real(DP), intent(in)    :: dt

    real(RP) :: rhoq(KA,IA,JA,QA) ! rho * q before precipitation
    real(RP) :: qflx(KA,IA,JA,QA)
    real(RP) :: eflx(KA,IA,JA)

    real(RP) :: rcdz  (KA,IA,JA)
    real(RP) :: rcdz_u(KA,IA,JA)
    real(RP) :: rcdz_v(KA,IA,JA)
    real(RP) :: rfdz  (KA,IA,JA)

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_Precipitation', 3)

    if ( TRACER_TYPE /= 'SUZUKI10' ) then
       do iq = I_QC, QQE
          call COMM_vars8( vterm(:,:,:,iq), iq )
       enddo
    endif
    do iq = I_QC, QQE
       call COMM_vars8( QTRC(:,:,:,iq), QQE+iq )
    enddo

!OCL XFILL
    flux_rain(:,:,:) = 0.0_RP
!OCL XFILL
    flux_snow(:,:,:) = 0.0_RP

    ! tracer/energy transport by falldown
    ! 1st order upwind, forward euler, velocity is always negative

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       rfdz(KS-1,i,j) = 1.0_RP / ( REAL_CZ(KS,i,j) - REAL_FZ(KS-1,i,j) )
       do k = KS, KE
          rcdz  (k,i,j) = 1.0_RP / ( REAL_FZ(k,i,j) - REAL_FZ(k-1,i,j) )
          rcdz_u(k,i,j) = 2.0_RP / ( ( REAL_FZ(k,i+1,j) - REAL_FZ(k-1,i+1,j) ) &
                                   + ( REAL_FZ(k,i  ,j) - REAL_FZ(k-1,i  ,j) ) )
          rcdz_v(k,i,j) = 2.0_RP / ( ( REAL_FZ(k,i,j+1) - REAL_FZ(k-1,i,j+1) ) &
                                   + ( REAL_FZ(k,i,j  ) - REAL_FZ(k-1,i,j  ) ) )
          rfdz(k,i,j) = 1.0_RP / ( REAL_CZ(k+1,i,j) - REAL_CZ(k,i,j) )
       enddo
    enddo
    enddo

    do iq = I_QC, QQE

       if ( TRACER_TYPE /= 'SUZUKI10' ) then
          call COMM_wait( vterm(:,:,:,iq), iq )
       endif
       call COMM_wait( QTRC(:,:,:,iq), QQE+iq )

       do j  = JS, JE
       do i  = IS, IE
          !--- mass flux for each mass tracer, upwind with vel < 0
          do k  = KS-1, KE-1
             qflx(k,i,j,iq) = vterm(k+1,i,j,iq) * DENS(k+1,i,j) * QTRC(k+1,i,j,iq) * J33G
          enddo
          qflx(KE,i,j,iq) = 0.0_RP

          !--- internal energy
          eflx(KS-1,i,j) = qflx(KS-1,i,j,iq) * temp(KS,i,j) * CVw(iq)
          do k  = KS, KE-1
             eflx(k,i,j) = qflx(k,i,j,iq) * temp(k+1,i,j) * CVw(iq)
             RHOE(k,i,j) = RHOE(k,i,j) - dt * ( eflx(k,i,j) - eflx(k-1,i,j) ) * rcdz(k,i,j)
          enddo
          eflx(KE,i,j) = 0.0_RP
          RHOE(KE,i,j) = RHOE(KE,i,j) - dt * ( - eflx(KE-1,i,j) ) * rcdz(KE,i,j)

          !--- potential energy
          eflx(KS-1,i,j) = qflx(KS-1,i,j,iq) * GRAV / rfdz(KS-1,i,j)
          do k  = KS, KE-1
             eflx(k,i,j) = qflx(k,i,j,iq) * GRAV / rfdz(k,i,j)
             RHOE(k,i,j) = RHOE(k,i,j) - dt * ( eflx(k,i,j) - eflx(k-1,i,j) ) * rcdz(k,i,j)
          enddo
          RHOE(KE,i,j) = RHOE(KE,i,j) - dt * ( - eflx(KE-1,i,j) ) * rcdz(KE,i,j)

          !--- momentum z (half level)
          do k  = KS-1, KE-1
             eflx(k,i,j) = 0.25_RP * ( vterm(k+1,i,j,iq) + vterm(k,i,j,iq) ) &
                                   * ( QTRC (k+1,i,j,iq) + QTRC (k,i,j,iq) ) &
                                   * MOMZ(k,i,j)
          enddo
          do k  = KS, KE-1
             MOMZ(k,i,j) = MOMZ(k,i,j) - dt * ( eflx(k+1,i,j) - eflx(k,i,j) ) * rfdz(k,i,j)
          enddo

          !--- momentum x
          eflx(KS-1,i,j) = 0.25_RP * ( vterm(KS,i,j,iq) + vterm(KS,i+1,j,iq) ) &
                                   * ( QTRC (KS,i,j,iq) + QTRC (KS,i+1,j,iq) ) &
                                   * MOMX(KS,i,j)
          do k  = KS, KE-1
             eflx(k,i,j) = 0.25_RP * ( vterm(k+1,i,j,iq) + vterm(k+1,i+1,j,iq) ) &
                                   * ( QTRC (k+1,i,j,iq) + QTRC (k+1,i+1,j,iq) ) &
                                   * MOMX(k+1,i,j)
             MOMX(k,i,j) = MOMX(k,i,j) - dt * ( eflx(k,i,j) - eflx(k-1,i,j) ) * rcdz_u(k,i,j)
          enddo
          MOMX(KE,i,j) = MOMX(KE,i,j) - dt * ( - eflx(KE-1,i,j) ) * rcdz_u(KE,i,j)

          !--- momentum y
          eflx(KS-1,i,j) = 0.25_RP * ( vterm(KS,i,j,iq) + vterm(KS,i,j+1,iq) ) &
                                   * ( QTRC (KS,i,j,iq) + QTRC (KS,i,j+1,iq) ) &
                                   * MOMY(KS,i,j)
          do k  = KS, KE-1
             eflx(k,i,j) = 0.25_RP * ( vterm(k+1,i,j,iq) + vterm(k+1,i,j+1,iq) ) &
                                   * ( QTRC (k+1,i,j,iq) + QTRC (k+1,i,j+1,iq) ) &
                                   * MOMY(k+1,i,j)
             MOMY(k,i,j) = MOMY(k,i,j) - dt * ( eflx(k,i,j) - eflx(k-1,i,j) ) * rcdz_v(k,i,j)
          enddo
          MOMY(KE,i,j) = MOMY(KE,i,j) - dt * ( - eflx(KE-1,i,j) ) * rcdz_v(KE,i,j)
       enddo
       enddo

    enddo ! mass tracer loop

    ! save previous value
    do iq = 1, QA
       do j = JS, JE
       do i = IS, IE
          rhoq(KS-1,i,j,iq) = DENS(KS,i,j) * QTRC(KS,i,j,iq)
          do k = KS, KE
             rhoq(k,i,j,iq) = DENS(k,i,j) * QTRC(k,i,j,iq)
          enddo
       enddo
       enddo
    enddo

    !--- mass flux for each tracer, upwind with vel < 0
    do iq = I_QC, QA

       do j  = JS, JE
       do i  = IS, IE
          do k  = KS-1, KE-1
             qflx(k,i,j,iq) = vterm(k+1,i,j,iq) * rhoq(k+1,i,j,iq)
          enddo
          qflx(KE,i,j,iq) = 0.0_RP
       enddo
       enddo
    enddo

    !--- update total density
    do iq = I_QC, QQE
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          DENS(k,i,j) = DENS(k,i,j) - dt * ( qflx(k,i,j,iq) - qflx(k-1,i,j,iq) ) * rcdz(k,i,j)
       enddo
       enddo
       enddo
    enddo ! mass tracer loop

    !--- update falling tracer (use updated total density)
    do iq = I_QC, QA
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          QTRC(k,i,j,iq) = ( rhoq(k,i,j,iq) - dt * ( qflx(k,i,j,iq) - qflx(k-1,i,j,iq) ) * rcdz(k,i,j) ) / DENS(k,i,j)
       enddo
       enddo
       enddo
    enddo

    !--- update no-falling tracer (use updated total density)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_QV) = rhoq(k,i,j,I_QV) / DENS(k,i,j)
    enddo
    enddo
    enddo

    !--- lowermost flux is saved for land process
    if ( QWS > 0 ) then
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS-1, KE
       do iq = QWS, QWE
          flux_rain(k,i,j) = flux_rain(k,i,j) - qflx(k,i,j,iq)
       enddo
       enddo
       enddo
       enddo
    endif
    if ( QIS > 0 ) then
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS-1, KE
       do iq = QIS, QIE
          flux_snow(k,i,j) = flux_snow(k,i,j) - qflx(k,i,j,iq)
       enddo
       enddo
       enddo
       enddo
    endif


    call PROF_rapend  ('MP_Precipitation', 3)

    return
  end subroutine ATMOS_PHY_MP_precipitation

end module scale_atmos_phy_mp_common
