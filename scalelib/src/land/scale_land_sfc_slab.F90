!-------------------------------------------------------------------------------
!> module LAND / Surface fluxes with slab land model
!!
!! @par Description
!!          Surface flux with slab land model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_land_sfc_slab
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_SFC_SLAB_setup
  public :: LAND_SFC_SLAB

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
  !-----------------------------------------------------------------------------
  logical,  private :: LST_UPDATE

  integer,  private :: LAND_SFC_SLAB_itr_max = 100 ! maximum iteration number

  real(RP), private :: LAND_SFC_SLAB_res_min = 1.0E+0_RP ! minimum number of residual
  real(RP), private :: LAND_SFC_SLAB_dTS_max = 5.0E-2_RP ! maximum delta surface temperature [K/s]

  logical, allocatable, private :: is_LND(:,:)

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_SFC_SLAB_setup( LAND_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_landuse, only: &
       LANDUSE_fact_land
    implicit none

    character(len=*), intent(in) :: LAND_TYPE

    NAMELIST / PARAM_LAND_SFC_SLAB / &
       LAND_SFC_SLAB_itr_max,  &
       LAND_SFC_SLAB_res_min,  &
       LAND_SFC_SLAB_dTS_max

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SLAB] / Categ[LAND SFC] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_SFC_SLAB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_SFC_SLAB. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LAND_SFC_SLAB)

    if( LAND_TYPE == 'CONST' ) then
       LST_UPDATE = .false.
    else if( LAND_TYPE == 'SLAB' ) then
       LST_UPDATE = .true.
    else
       write(*,*) 'xxx wrong LAND_TYPE. Check!'
       call PRC_MPIstop
    end if

    ! judge to run slab land model
    allocate( is_LND(IA,JA) )

    do j = JS, JE
    do i = IS, IE
      if( LANDUSE_fact_land(i,j) > 0.0_RP ) then
        is_LND(i,j) = .true.
      else
        is_LND(i,j) = .false.
      end if
    end do
    end do

    return
  end subroutine LAND_SFC_SLAB_setup

  !-----------------------------------------------------------------------------
  subroutine LAND_SFC_SLAB( &
        LST_t,      &
        ZMFLX,      &
        XMFLX,      &
        YMFLX,      &
        SHFLX,      &
        LHFLX,      &
        GHFLX,      &
        U10,        &
        V10,        &
        T2,         &
        Q2,         &
        TMPA,       &
        PRSA,       &
        WA,         &
        UA,         &
        VA,         &
        RHOA,       &
        QVA,        &
        Z1,         &
        PBL,        &
        PRSS,       &
        LWD,        &
        SWD,        &
        TG,         &
        LST,        &
        QVEF,       &
        ALB_LW,     &
        ALB_SW,     &
        DZG,        &
        TCS,        &
        Z0M,        &
        Z0H,        &
        Z0E,        &
        dt          )
    use scale_process, only: &
       PRC_myrank,  &
       PRC_MPIstop
    use scale_const, only: &
      Rdry  => CONST_Rdry, &
      CPdry => CONST_CPdry, &
      STB   => CONST_STB
    use scale_grid_index
    use scale_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_templhv
    use scale_bulkflux, only: &
      BULKFLUX
    implicit none

    ! parameters
    real(RP), parameter :: dTS0     = 1.0E-4_RP ! delta surface temp.
    real(RP), parameter :: dres_lim = 1.0E+2_RP ! limiter of d(residual)

    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0_RP    ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5_RP    ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1_RP    ! factor b in Tomita (2009)

    ! arguments
    real(RP), intent(out) :: LST_t(IA,JA) ! tendency of LST
    real(RP), intent(out) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SHFLX(IA,JA) ! sensible heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: LHFLX(IA,JA) ! latent heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: GHFLX(IA,JA) ! ground heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: U10  (IA,JA) ! velocity u at 10m [m/s]
    real(RP), intent(out) :: V10  (IA,JA) ! velocity v at 10m [m/s]
    real(RP), intent(out) :: T2   (IA,JA) ! temperature at 2m [K]
    real(RP), intent(out) :: Q2   (IA,JA) ! water vapor at 2m [kg/kg]

    real(RP), intent(in) :: TMPA(IA,JA) ! temperature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: PRSA(IA,JA) ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: WA  (IA,JA) ! velocity w at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: UA  (IA,JA) ! velocity u at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: VA  (IA,JA) ! velocity v at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: RHOA(IA,JA) ! density at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in) :: QVA (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: Z1  (IA,JA) ! cell center height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: PBL (IA,JA) ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: PRSS(IA,JA) ! pressure at the surface [Pa]
    real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface [J/m2/s]
    real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface [J/m2/s]

    real(RP), intent(in) :: TG    (IA,JA) ! soil temperature [K]
    real(RP), intent(in) :: LST   (IA,JA) ! land surface temperature [K]
    real(RP), intent(in) :: QVEF  (IA,JA) ! efficiency of evaporation [0-1]
    real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
    real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
    real(RP), intent(in) :: DZG   (IA,JA) ! soil depth [m]
    real(RP), intent(in) :: TCS   (IA,JA) ! thermal conductivity for soil [J/m/K/s]
    real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momemtum [m]
    real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]
    real(DP), intent(in) :: dt            ! delta time

    ! works
    real(RP) :: LST1(IA,JA)

    real(RP) :: res    ! residual
    real(RP) :: dres   ! d(residual)/dLST
    real(RP) :: oldres ! residual in previous step
    real(RP) :: redf   ! reduced factor

    real(RP) :: Ustar, dUstar ! friction velocity [m]
    real(RP) :: Tstar, dTstar ! friction potential temperature [K]
    real(RP) :: Qstar, dQstar ! friction water vapor mass ratio [kg/kg]
    real(RP) :: Uabs, dUabs   ! modified absolute velocity [m/s]
    real(RP) :: SQV, dSQV     ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: LHV(IA,JA)    ! latent heat for vaporization depending on temperature [J/kg]

    integer :: i, j, n
    !---------------------------------------------------------------------------

    ! copy land surfce temperature for iteration
    do j = JS, JE
    do i = IS, IE
      LST1(i,j) = LST(i,j)
    end do
    end do

    call ATMOS_THERMODYN_templhv( LHV, TMPA )

    ! update surface temperature
    if( LST_UPDATE ) then

      do j = JS, JE
      do i = IS, IE

        if( is_LND(i,j) ) then

          redf   = 1.0_RP
          oldres = 1.0E+10_RP

          ! modified Newton-Raphson method (Tomita 2009)
          do n = 1, LAND_SFC_SLAB_itr_max

            call qsat( SQV,       & ! [OUT]
                       LST1(i,j), & ! [IN]
                       PRSS(i,j)  ) ! [IN]
            call qsat( dSQV,           & ! [OUT]
                       LST1(i,j)+dTS0, & ! [IN]
                       PRSS(i,j)       ) ! [IN]

            call BULKFLUX( &
                Ustar,     & ! [OUT]
                Tstar,     & ! [OUT]
                Qstar,     & ! [OUT]
                Uabs,      & ! [OUT]
                TMPA(i,j), & ! [IN]
                LST1(i,j), & ! [IN]
                PRSA(i,j), & ! [IN]
                PRSS(i,j), & ! [IN]
                QVA (i,j), & ! [IN]
                SQV,       & ! [IN]
                UA  (i,j), & ! [IN]
                VA  (i,j), & ! [IN]
                Z1  (i,j), & ! [IN]
                PBL (i,j), & ! [IN]
                Z0M (i,j), & ! [IN]
                Z0H (i,j), & ! [IN]
                Z0E (i,j)  ) ! [IN]

            call BULKFLUX( &
                dUstar,         & ! [OUT]
                dTstar,         & ! [OUT]
                dQstar,         & ! [OUT]
                dUabs,          & ! [OUT]
                TMPA(i,j),      & ! [IN]
                LST1(i,j)+dTS0, & ! [IN]
                PRSA(i,j),      & ! [IN]
                PRSS(i,j),      & ! [IN]
                QVA (i,j),      & ! [IN]
                dSQV,           & ! [IN]
                UA  (i,j),      & ! [IN]
                VA  (i,j),      & ! [IN]
                Z1  (i,j),      & ! [IN]
                PBL (i,j),      & ! [IN]
                Z0M (i,j),      & ! [IN]
                Z0H (i,j),      & ! [IN]
                Z0E (i,j)       ) ! [IN]

            ! calculation for residual
            res = ( 1.0_RP - ALB_SW(i,j) ) * SWD(i,j) &
                + ( 1.0_RP - ALB_LW(i,j) ) * ( LWD(i,j) - STB * LST1(i,j)**4 ) &
                + CPdry    * RHOA(i,j) * Ustar * Tstar &
                + LHV(i,j) * RHOA(i,j) * Ustar * Qstar * QVEF(i,j) &
                - 2.0_RP * TCS(i,j) * ( LST1(i,j) - TG(i,j)  ) / DZG(i,j)

            if( abs(res) < LAND_SFC_SLAB_res_min ) then
              ! iteration converged
              exit
            end if

            ! calculation for d(residual)/dLST
            dres = -4.0_RP * ( 1.0_RP - ALB_LW(i,j) ) * STB * LST1(i,j)**3 &
                 + CPdry    * RHOA(i,j) * ( (dUstar-Ustar)/dTS0 * Tstar + Ustar * (dTstar-Tstar)/dTS0 ) &
                 + LHV(i,j) * RHOA(i,j) * ( (dUstar-Ustar)/dTS0 * Qstar + Ustar * (dQstar-Qstar)/dTS0 ) * QVEF(i,j) &
                 - 2.0_RP * TCS(i,j) / DZG(i,j)

            if( abs(dres)*dres_lim < abs(res) ) then
              ! stop iteration to prevent numerical error
              exit
            end if

            ! calculate reduced factor
            if( dres < 0.0_RP ) then
              if( abs(res) > abs(oldres) ) then
                redf = max( TFa*abs(redf), redf_min )
              else
                redf = min( TFb*abs(redf), redf_max )
              end if
            else
              redf = -1.0_RP
            end if

            ! estimate next surface temperature
            LST1(i,j) = LST1(i,j) - redf * res / dres

            ! save residual in this step
            oldres = res

          end do

          ! update land surface temperature with limitation
          LST1(i,j) = min( max( LST1(i,j), &
                                LST (i,j) - LAND_SFC_SLAB_dTS_max * dt ), &
                                LST (i,j) + LAND_SFC_SLAB_dTS_max * dt )

          if( n > LAND_SFC_SLAB_itr_max ) then
            ! check NaN
            if ( .NOT. ( res > -1.0_RP .OR. res < 1.0_RP ) ) then ! must be NaN
               write(*,*) 'xxx NaN is detected for land surface temperature in rank, i, j: ', PRC_myrank, i, j

               write(*,*) 'DEBUG Message --- Residual                            [J/m2/s]  :', res
               write(*,*) 'DEBUG Message --- delta Residual                      [J/m2/s]  :', dres
               write(*,*) ''
               write(*,*) 'DEBUG Message --- temperature                         [K]       :', TMPA  (i,j)
               write(*,*) 'DEBUG Message --- pressure                            [Pa]      :', PRSA  (i,j)
               write(*,*) 'DEBUG Message --- velocity                            [m/s]     :', WA    (i,j)
               write(*,*) 'DEBUG Message --- velocity u                          [m/s]     :', UA    (i,j)
               write(*,*) 'DEBUG Message --- velocity v                          [m/s]     :', VA    (i,j)
               write(*,*) 'DEBUG Message --- density                             [kg/m3]   :', RHOA  (i,j)
               write(*,*) 'DEBUG Message --- water vapor mass ratio              [kg/kg]   :', QVA   (i,j)
               write(*,*) 'DEBUG Message --- cell center height                  [m]       :', Z1    (i,j)
               write(*,*) 'DEBUG Message --- the top of atmospheric mixing layer [m]       :', PBL   (i,j)
               write(*,*) 'DEBUG Message --- pressure at the surface             [Pa]      :', PRSS  (i,j)
               write(*,*) 'DEBUG Message --- downward long-wave radiation flux   [J/m2/s]  :', LWD   (i,j)
               write(*,*) 'DEBUG Message --- downward short-wave radiation flux  [J/m2/s]  :', SWD   (i,j)
               write(*,*) ''
               write(*,*) 'DEBUG Message --- soil temperature                    [K]       :', TG    (i,j)
               write(*,*) 'DEBUG Message --- land surface temperature            [K]       :', LST   (i,j)
               write(*,*) 'DEBUG Message --- efficiency of evaporation           [0-1]     :', QVEF  (i,j)
               write(*,*) 'DEBUG Message --- surface albedo for LW               [0-1]     :', ALB_LW(i,j)
               write(*,*) 'DEBUG Message --- surface albedo for SW               [0-1]     :', ALB_SW(i,j)
               write(*,*) 'DEBUG Message --- soil depth                          [m]       :', DZG   (i,j)
               write(*,*) 'DEBUG Message --- thermal conductivity for soil       [J/m/K/s] :', TCS   (i,j)
               write(*,*) 'DEBUG Message --- roughness length for momemtum       [m]       :', Z0M   (i,j)
               write(*,*) 'DEBUG Message --- roughness length for heat           [m]       :', Z0H   (i,j)
               write(*,*) 'DEBUG Message --- roughness length for vapor          [m]       :', Z0E   (i,j)
               write(*,*) ''
               write(*,*) 'DEBUG Message --- friction velocity                   [m]       :', Ustar
               write(*,*) 'DEBUG Message --- friction potential temperature      [K]       :', Tstar
               write(*,*) 'DEBUG Message --- friction water vapor mass ratio     [kg/kg]   :', Qstar
               write(*,*) 'DEBUG Message --- modified absolute velocity          [m/s]     :', Uabs
               write(*,*) 'DEBUG Message --- next land surface temperature       [K]       :', LST1  (i,j)

               call PRC_MPIstop
            endif

            ! land surface temperature was not converged
            if( IO_L ) write(IO_FID_LOG,*) 'Warning: land surface tempearture was not converged.'

            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- Residual                            [J/m2/s]  :', res
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- delta Residual                      [J/m2/s]  :', dres
            if( IO_L ) write(IO_FID_LOG,*) ''
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- temperature                         [K]       :', TMPA  (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- pressure                            [Pa]      :', PRSA  (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- velocity                            [m/s]     :', WA    (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- velocity u                          [m/s]     :', UA    (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- velocity v                          [m/s]     :', VA    (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- density                             [kg/m3]   :', RHOA  (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- water vapor mass ratio              [kg/kg]   :', QVA   (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- cell center height                  [m]       :', Z1    (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- the top of atmospheric mixing layer [m]       :', PBL   (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- pressure at the surface             [Pa]      :', PRSS  (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- downward long-wave radiation flux   [J/m2/s]  :', LWD   (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- downward short-wave radiation flux  [J/m2/s]  :', SWD   (i,j)
            if( IO_L ) write(IO_FID_LOG,*) ''
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- soil temperature                    [K]       :', TG    (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- land surface temperature            [K]       :', LST   (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- efficiency of evaporation           [0-1]     :', QVEF  (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- surface albedo for LW               [0-1]     :', ALB_LW(i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- surface albedo for SW               [0-1]     :', ALB_SW(i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- soil depth                          [m]       :', DZG   (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- thermal conductivity for soil       [J/m/K/s] :', TCS   (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- roughness length for momemtum       [m]       :', Z0M   (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- roughness length for heat           [m]       :', Z0H   (i,j)
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- roughness length for vapor          [m]       :', Z0E   (i,j)
            if( IO_L ) write(IO_FID_LOG,*) ''
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- friction velocity                   [m]       :', Ustar
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- friction potential temperature      [K]       :', Tstar
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- friction water vapor mass ratio     [kg/kg]   :', Qstar
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- modified absolute velocity          [m/s]     :', Uabs
            if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- next land surface temperature       [K]       :', LST1  (i,j)
          end if

        end if

        ! calculate tendency
        LST_t(i,j) = ( LST1(i,j) - LST(i,j) ) / dt

      end do
      end do

    ! not update temperature
    else

      do j = JS, JE
      do i = IS, IE
         LST_t(i,j) = 0.0_RP
      end do
      end do

    end if


    ! calculate surface flux
    do j = JS, JE
    do i = IS, IE

      if( is_LND(i,j) ) then

        call qsat( SQV,       & ! [OUT]
                   LST1(i,j), & ! [IN]
                   PRSS(i,j)  ) ! [IN]

        call BULKFLUX( &
            Ustar,     & ! [OUT]
            Tstar,     & ! [OUT]
            Qstar,     & ! [OUT]
            Uabs,      & ! [OUT]
            TMPA(i,j), & ! [IN]
            LST1(i,j), & ! [IN]
            PRSA(i,j), & ! [IN]
            PRSS(i,j), & ! [IN]
            QVA (i,j), & ! [IN]
            SQV,       & ! [IN]
            UA  (i,j), & ! [IN]
            VA  (i,j), & ! [IN]
            Z1  (i,j), & ! [IN]
            PBL (i,j), & ! [IN]
            Z0M (i,j), & ! [IN]
            Z0H (i,j), & ! [IN]
            Z0E (i,j)  ) ! [IN]

        ZMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * WA(i,j)
        XMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * UA(i,j)
        YMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * VA(i,j)
        SHFLX(i,j) = -CPdry    * RHOA(i,j) * Ustar * Tstar
        LHFLX(i,j) = -LHV(i,j) * RHOA(i,j) * Ustar * Qstar * QVEF(i,j)
        GHFLX(i,j) = -2.0_RP * TCS(i,j) * ( LST1(i,j) - TG(i,j) ) / DZG(i,j)

        ! calculation for residual
        res = ( 1.0_RP - ALB_SW(i,j) ) * SWD(i,j) &
            + ( 1.0_RP - ALB_LW(i,j) ) * ( LWD(i,j) - STB * LST1(i,j)**4 ) &
            - SHFLX(i,j) - LHFLX(i,j) + GHFLX(i,j)

        ! put residual in ground heat flux
        GHFLX(i,j) = GHFLX(i,j) - res

        ! diagnostic variables
        U10(i,j) = UA  (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        V10(i,j) = VA  (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        T2 (i,j) = LST1(i,j) + ( TMPA(i,j) - LST1(i,j) ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0H(i,j) ) ) &
                                                         / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0H(i,j) ) )
        Q2 (i,j) = SQV       + (  QVA(i,j) - SQV       ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0E(i,j) ) ) &
                                                         / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0E(i,j) ) )
      else

        ! not calculate surface flux
        ZMFLX(i,j) = 0.0_RP
        XMFLX(i,j) = 0.0_RP
        YMFLX(i,j) = 0.0_RP
        SHFLX(i,j) = 0.0_RP
        LHFLX(i,j) = 0.0_RP
        GHFLX(i,j) = 0.0_RP
        U10  (i,j) = 0.0_RP
        V10  (i,j) = 0.0_RP
        T2   (i,j) = 0.0_RP
        Q2   (i,j) = 0.0_RP

      end if

    end do
    end do

    return
  end subroutine LAND_SFC_SLAB

end module scale_land_sfc_slab
