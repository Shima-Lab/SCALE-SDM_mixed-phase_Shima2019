!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Covert super-droplets to fluid variables
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-07-14 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2014-07-14 (S.Shima) [rev] sdm_sd2prec added
!! @li      2014-07-22 (Y.Sato ) [mod] Modify the definition of kl and ku for calculating drate
!! @li      2014-07-24 (Y.Sato ) [mod] Modify bugs accessing upper/lower boundary
!! @li      2015-06-27 (S.Shima) [add] Add fapp_start/stop calls to monitor the performance
!! @li      2015-06-27 (S.Shima) [mod] Working arrays are introduced to parallelize sd2rhow and sd2rhocr 
!! @li      2015-07-30 (Y.Sato)  [add] Add "ifdef" for fapp and fipp module
!! @li      2016-07-16 (S.Shima) [add] sdm_sd2qiqsqg added
!! @li      2016-07-16 (S.Shima) [mod] sdm_sd2qcqr, sdm_sd2rhocr, sdm_sd2rhow modified
!! @li      2016-07-16 (S.Shima) [add] sdm_sd2rhosol added
!! @li      2016-07-16 (S.Shima) [mod] sdm_sd2prec modified
!! @li      2016-07-22 (S.Shima) [fix] sdm_sd2rhoisg is fixed not to overcount graupel and snow 
!! @li      2018-07-01 (S.Shima) [mod] categorization of ice, snow, graupel, based on mrime and nmono
!! @li      2019-01-09 (S.Shima) [mod] difinition of ice, snow, graupel
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_sd2fluid
  use scale_precision

  implicit none
  private
  public :: sdm_sd2prec, sdm_sd2rhow, sdm_sd2rhocr, sdm_sd2qcqr, sdm_sd2qiqsqg, sdm_sd2rhosol

contains
  subroutine sdm_sd2prec(dtb_crs,                        &
       prec,sd_num,sd_n,sd_liqice,sd_x,sd_y,     &
       sd_r,sdi,sd_ri,sd_rj,sd_rk,ilist_s,ilist_r,pr_sdm)
    use scale_grid_index, only: &
         IA,JA,KA,IS,IE,JS,JE,KS,KE
    use scale_const, only: &
         rw => CONST_DWATR
    use m_sdm_common, only: &
         F_THRD, ONE_PI, dxiv_sdm, dyiv_sdm, VALID2INVALID, PREC2INVALID,INVALID, &
         sdm_cold, i2, sdicedef, STAT_LIQ, STAT_ICE
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    real(RP), intent(in) :: dtb_crs
    integer, intent(in) :: sd_num  ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    ! Input and output variables
    real(RP), intent(out) :: sd_ri(1:sd_num)   ! index-i(real) of super-droplets
    real(RP), intent(out) :: sd_rj(1:sd_num)   ! index-j(real) of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    real(RP), intent(inout) :: prec(IA,JA,1:6) ! precipitation rate [m/s] and accumlation [m]
                                               ! 1:rain rate, 2:rain accumulation
                                               ! 3:snow rate, 4:snow accumulation
                                               ! 5:precipitation rate, 6:precipitation accumulation
    ! Output variables
    real(RP), intent(out) :: pr_sdm(1:IA,1:JA) ! temporary buffer of CReSS dimension
    integer, intent(out) :: ilist_r(1:sd_num) ! buffer for list vectorization
    integer, intent(out) :: ilist_s(1:sd_num) ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(IA,JA) ! coef.
    real(RP) :: dtmp         ! temporary variables
    real(RP) :: dtbiv        ! 1.e0 / time step

    integer :: tlist_r, tlist_s              ! total list number
    integer :: cnt_r, cnt_s                ! counter

    integer :: i, j, m, n ! index
    !-------------------------------------------------------------------

    ! Initialize
    dtbiv = 1.0d0 / dtb_crs
    dcoef(1:IA,1:JA)=0.0d0
    do i = IS, IE
    do j = JS, JE
       dcoef(i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j)
    enddo
    enddo

    ! Get index list for compressing buffer.
    cnt_r=0
    cnt_s=0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID .and. &
            sd_rk(n)>PREC2INVALID ) then
          if( sd_liqice(n) .eq. STAT_LIQ) then
             cnt_r = cnt_r + 1
             ilist_r(cnt_r) = n
          else if( sd_liqice(n) .eq. STAT_ICE) then
             cnt_s = cnt_s + 1
             ilist_s(cnt_s) = n
          end if
       end if
    end do

    tlist_r = cnt_r
    tlist_s = cnt_s

    ! Get rain rate and accumulation
    if( tlist_r>0 ) then
       do j=1,JA
       do i=1,IA
          pr_sdm(i,j) = 0.d0
       end do
       end do

       !### get horizontal face index(real) of super-droplets ###!
       call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
       call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

       do m=1,tlist_r
          n=ilist_r(m)
          i=floor(sd_ri(n))+1
          j=floor(sd_rj(n))+1

          pr_sdm(i,j) = pr_sdm(i,j)                         &
               + sd_r(n) * sd_r(n) * sd_r(n)           &
               * real(sd_n(n),kind=RP)

          sd_rk(n) = INVALID     !! convert to invalid
       end do

       !### convert super-droplets to precipitation ###!

       do j=1,JA
       do i=1,IA

          dtmp = real( pr_sdm(i,j) * dcoef(i,j) )

          !! rain fall rate
          prec(i,j,1) = dtmp * dtbiv

          !! rain accumulation
          prec(i,j,2) = prec(i,j,2) + dtmp
          
       end do
       end do

    end if

    ! Get snow rate and accumulation
    if( sdm_cold .and. (tlist_s>0) ) then
       do j=1,JA
       do i=1,IA
          pr_sdm(i,j) = 0.d0
       end do
       end do

       !### get horizontal face index(real) of super-droplets ###!
       call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
       call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

       do m=1,tlist_s
          n=ilist_s(m)
          i=floor(sd_ri(n))+1
          j=floor(sd_rj(n))+1

          pr_sdm(i,j) = pr_sdm(i,j) &
               + sdi%re(n) * sdi%re(n) * sdi%rp(n) * sdi%rho(n) / rw &
               * real(sd_n(n),kind=RP)

          sd_rk(n) = INVALID     !! convert to invalid
       end do

       !### convert super-droplets to precipitation ###!

       do j=1,JA
       do i=1,IA

          dtmp = real( pr_sdm(i,j) * dcoef(i,j) )

          !! snow fall rate
          prec(i,j,3) = dtmp * dtbiv

          !! snow accumulation
          prec(i,j,4) = prec(i,j,4) + dtmp
          
       end do
       end do

    end if

    ! Get precipitation rate and accumulation
    do j=1,JA
    do i=1,IA

       !! precipitation rate
       prec(i,j,5) = prec(i,j,1) + prec(i,j,3) 

       !! precipitation accumulation
       prec(i,j,6) = prec(i,j,2) + prec(i,j,4) 
          
    end do
    end do

    return
  end subroutine sdm_sd2prec
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2qcqr(DENS,QC,QR,        &
       zph_crs,           &
       sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk, &
       sd_rkl,sd_rku,                           &
       rhoc_sdm,rhor_sdm,              &
       sd_itmp1,sd_itmp2,crs_dtmp1,crs_dtmp2    )
    use scale_grid_index, only: &
         IA,JA,KA,IS,IE,JS,JE,KS,KE
    use m_sdm_common, only: &
         sdm_rqc2qr,i2
    ! Input variables
    real(RP), intent(in) :: DENS(KA,IA,JA) ! Density     [kg/m3]
    real(RP), intent(out):: QC(KA,IA,JA)   ! rhoc/rho
    real(RP), intent(out):: QR(KA,IA,JA)   ! rhor/rho
    real(RP), intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinate
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num)      ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)      ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num)      ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)     ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)     ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)     ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA)  ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA)  ! upper boundary index[k/real] at scalar point in SDM calculation area
    ! Output variables
    real(RP), intent(out) :: rhoc_sdm(KA,IA,JA)   ! density of cloud water
    real(RP), intent(out) :: rhor_sdm(KA,IA,JA)   ! density of rain water
    real(RP), intent(out) :: crs_dtmp1(KA,IA,JA)    ! temporary buffer of CReSS dimension
    real(RP), intent(out) :: crs_dtmp2(KA,IA,JA)    ! temporary buffer of CReSS dimension
    integer, intent(out) :: sd_itmp1(1:sd_num)    ! temporary array of the size of the number of super-droplets.
    integer, intent(out) :: sd_itmp2(1:sd_num)    ! temporary array of the size of the number of super-droplets.
    ! Work variables
    integer :: i, j, k   ! index

    !-------------------------------------------------------------------

    ! Get density of water hydrometeor

    !### the case updating water hydrometeor by SDM ###!
    !! convert super-droplets to density of cloud water
    !! and rain water.
    
    call sdm_sd2rhocr(sdm_rqc2qr,                          &
         zph_crs,rhoc_sdm,rhor_sdm,           &
         sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,    &
         sd_rkl,sd_rku,                       &
         crs_dtmp1,crs_dtmp2,sd_itmp1,sd_itmp2)

    ! Convert water hydrometeor density to mixing ratio
    do k=KS,KE
    do j=JS,JE
    do i=IS,IE
       QC(k,i,j) = rhoc_sdm(k,i,j)/DENS(k,i,j)
       QR(k,i,j) = rhor_sdm(k,i,j)/DENS(k,i,j)
    end do
    end do
    end do

    return
  end subroutine sdm_sd2qcqr
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2qiqsqg(DENS,QI,QS,QG,        &
       zph_crs,           &
       sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk, &
       sd_rkl,sd_rku,                           &
       rhoi_sdm,rhos_sdm,rhog_sdm,              &
       sd_itmp1,sd_itmp2,sd_itmp3,crs_dtmp1,crs_dtmp2,crs_dtmp3    )
    use scale_grid_index, only: &
         IA,JA,KA,IS,IE,JS,JE,KS,KE
    use m_sdm_common, only: &
         sdm_rqi2qs, i2, sdicedef
    ! Input variables
    real(RP), intent(in) :: DENS(KA,IA,JA) ! Density     [kg/m3]
    real(RP), intent(out):: QI(KA,IA,JA)   ! rhoi/rho
    real(RP), intent(out):: QS(KA,IA,JA)   ! rhos/rho
    real(RP), intent(out):: QG(KA,IA,JA)   ! rhog/rho
    real(RP), intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinate
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num)      ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)      ! y-coordinate of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)     ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)     ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)     ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA)  ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA)  ! upper boundary index[k/real] at scalar point in SDM calculation area
    ! Output variables
    real(RP), intent(inout) :: rhoi_sdm(KA,IA,JA)   ! density of cloud ice
    real(RP), intent(inout) :: rhos_sdm(KA,IA,JA)   ! density of snow
    real(RP), intent(inout) :: rhog_sdm(KA,IA,JA)   ! density of graupel
    real(RP), intent(out) :: crs_dtmp1(KA,IA,JA)    ! temporary buffer of CReSS dimension
    real(RP), intent(out) :: crs_dtmp2(KA,IA,JA)    ! temporary buffer of CReSS dimension
    real(RP), intent(out) :: crs_dtmp3(KA,IA,JA)    ! temporary buffer of CReSS dimension
    integer, intent(out) :: sd_itmp1(1:sd_num)    ! temporary array of the size of the number of super-droplets.
    integer, intent(out) :: sd_itmp2(1:sd_num)    ! temporary array of the size of the number of super-droplets.
    integer, intent(out) :: sd_itmp3(1:sd_num)    ! temporary array of the size of the number of super-droplets.
    ! Work variables
    integer :: i, j, k   ! index

    !-------------------------------------------------------------------

    ! Get density of solid hydrometeor
    call sdm_sd2rhoisg(sdm_rqi2qs,                        &
         zph_crs,rhoi_sdm,rhos_sdm,rhog_sdm,              &
         sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk, &
         sd_rkl,sd_rku,                                   &
         crs_dtmp1,crs_dtmp2,crs_dtmp3,sd_itmp1,sd_itmp2,sd_itmp3)

    ! Convert water hydrometeor density to mixing ratio
    do k=KS,KE
    do j=JS,JE
    do i=IS,IE
       QI(k,i,j) = rhoi_sdm(k,i,j)/DENS(k,i,j)
       QS(k,i,j) = rhos_sdm(k,i,j)/DENS(k,i,j)
       QG(k,i,j) = rhog_sdm(k,i,j)/DENS(k,i,j)
    end do
    end do
    end do

    return
  end subroutine sdm_sd2qiqsqg
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2rhocr(sdm_rqc2qr,                        &
       zph_crs,rhoc_sdm,rhor_sdm,         &
       sd_num, sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk, &
       sd_rkl,sd_rku,                     &
       liqc_sdm,liqr_sdm,ilist_c,ilist_r)
    use scale_grid_index, only: &
         IA,JA,KA,IS,IE,JS,JE,KS,KE
    use scale_const, only: &
         rw => CONST_DWATR
    use m_sdm_common, only: &
         F_THRD, ONE_PI, dxiv_sdm, dyiv_sdm, VALID2INVALID, num_threads, i2, STAT_LIQ
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    real(RP),intent(in) :: sdm_rqc2qr          ! Threshould between qc and qr [m]
    real(RP),intent(in) :: zph_crs(KA,IA,JA)   ! z physical coordinate
    integer, intent(in) :: sd_num              ! Number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num)     ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)     ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num)     ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)    ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA)      ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA)      ! upper boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(out) :: rhoc_sdm(KA,IA,JA)  ! densitiy of cloud water
    real(RP), intent(out) :: rhor_sdm(KA,IA,JA)  ! densitiy of rain water
    real(RP), intent(out) :: liqc_sdm(KA,IA,JA)  ! liquid cloud water of super-droplets
    real(RP), intent(out) :: liqr_sdm(KA,IA,JA)  ! liquid rain water of super-droplets
    integer, intent(out) :: ilist_c(1:sd_num)  ! buffer for list vectorization
    integer, intent(out) :: ilist_r(1:sd_num)  ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(KA,IA,JA)   ! coef.
    real(RP) :: drate             ! temporary
    integer :: tlist_c            ! total list number for cloud
    integer :: tlist_r            ! total list number for rain
    integer :: ccnt               ! counter
    integer :: rcnt               ! counter
    integer :: i                  ! index
    integer :: j                  ! index
    integer :: k                  ! index
    integer :: kl                 ! index
    integer :: ku                 ! index
    integer :: m                  ! index
    integer :: n                  ! index
    integer :: i_threads
    real(RP) :: tmp_liq_sdm(num_threads,KA,IA,JA)
    !-----------------------------------------------------------------------------
    
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_sd2rhocr",1,1)
#endif
    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)
    
    ! Initialize
    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
       dcoef(k,i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j) &
            / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
    enddo
    enddo
    enddo

    tlist_c = 0
    tlist_r = 0

    do k=1,KA
    do j=1,JA
    do i=1,IA
       liqc_sdm(k,i,j) = 0.0_RP
       liqr_sdm(k,i,j) = 0.0_RP
       rhoc_sdm(k,i,j) = 0.0_RP
       rhor_sdm(k,i,j) = 0.0_RP
    end do
    end do
    end do

    ! Get index list for compressing buffer.
    ccnt = 0
    rcnt = 0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n)/=STAT_LIQ ) cycle

       if( sd_r(n)<real(sdm_rqc2qr,kind=RP) ) then

          !### cloud-water ###!

          ccnt = ccnt + 1
          ilist_c(ccnt) = n

       else
          
          !### rain-water ###!
          rcnt = rcnt + 1
          ilist_r(rcnt) = n

       end if

    end do

    tlist_c = ccnt
    tlist_r = rcnt

    ! Get density of cloud-water and rain-water.

    !### cloud-water ###!

    if( tlist_c>0 ) then

       tmp_liq_sdm(:,:,:,:) = 0.0_RP

       do i_threads=1,num_threads
       do m=i_threads,tlist_c,num_threads
          n = ilist_c(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          tmp_liq_sdm(i_threads,k,i,j) = tmp_liq_sdm(i_threads,k,i,j)                        &
               &                         + sd_r(n) * sd_r(n) * sd_r(n) * rw            &
               &                         * real(sd_n(n),kind=RP)
       end do
       end do

       do k=1,KA
       do j=1,JA
       do i=1,IA
       do i_threads=1,num_threads
          liqc_sdm(k,i,j)=liqc_sdm(k,i,j)+tmp_liq_sdm(i_threads,k,i,j)
       end do
       end do
       end do
       end do

!!$       do m=1,tlist_c
!!$          n = ilist_c(m)
!!$
!!$          i = floor(sd_ri(n))+1
!!$          j = floor(sd_rj(n))+1
!!$          k = floor(sd_rk(n))+1
!!$
!!$          liqc_sdm(k,i,j) = liqc_sdm(k,i,j)                     &
!!$               + sd_r(n) * sd_r(n) * sd_r(n)         &
!!$               * real(sd_n(n),kind=RP)
!!$       end do

       !=== adjust cloud-water in verical boundary. ===!

       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             liqc_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             liqc_sdm(kl,i,j) = liqc_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             liqc_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             liqc_sdm(ku,i,j) = liqc_sdm(ku,i,j)/drate
          end if

       end do
       end do

       !=== convert super-droplets to density of cloud-water. ===!

       do k=KS,KE
       do j=JS,JE
       do i=IS,IE
          rhoc_sdm(k,i,j) = liqc_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do

    end if


    !### rain-water ###!

    if( tlist_r>0 ) then

       tmp_liq_sdm(:,:,:,:) = 0.0_RP

       do i_threads=1,num_threads
       do m=i_threads,tlist_r,num_threads
          n = ilist_r(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          tmp_liq_sdm(i_threads,k,i,j) = tmp_liq_sdm(i_threads,k,i,j)                        &
               &                         + sd_r(n) * sd_r(n) * sd_r(n) * rw           &
               &                         * real(sd_n(n),kind=RP)
       end do
       end do

       do k=1,KA
       do j=1,JA
       do i=1,IA
       do i_threads=1,num_threads
          liqr_sdm(k,i,j)=liqr_sdm(k,i,j)+tmp_liq_sdm(i_threads,k,i,j)
       end do
       end do
       end do
       end do

!!$       do m=1,tlist_r
!!$          n = ilist_r(m)
!!$
!!$          i = floor(sd_ri(n))+1
!!$          j = floor(sd_rj(n))+1
!!$          k = floor(sd_rk(n))+1
!!$
!!$          liqr_sdm(k,i,j) = liqr_sdm(k,i,j)                     &
!!$               + sd_r(n) * sd_r(n) * sd_r(n)         &
!!$               * real(sd_n(n),kind=RP)
!!$       end do
       
       !=== adjust rain-water in verical boundary. ===!

       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)

          if( drate<0.5d0 ) then
             liqr_sdm(kl,i,j) = 0.e0           !! <50% in share
          else
             liqr_sdm(kl,i,j) = liqr_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.5d0 ) then
             liqr_sdm(ku,i,j) = 0.e0           !! <50% in share
          else
             liqr_sdm(ku,i,j) = liqr_sdm(ku,i,j)/drate
          end if

       end do
       end do

       !=== convert super-droplets to density of rain-water. ===!
       
       do k=KS,KE
       do j=JS,JE
       do i=IS,IE
          rhor_sdm(k,i,j) = liqr_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do

    end if

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_sd2rhocr",1,1)
#endif
    return
  end subroutine sdm_sd2rhocr
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2rhow(zph_crs,rhow_sdm,sd_num,sd_n,sd_liqice,  &
       sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,sd_rkl,sd_rku,      &
       liqw_sdm,ilist)
    use scale_grid_index, only: &
         IA,JA,KA,IS,IE,JS,JE,KS,KE
    use scale_const, only: &
         rw => CONST_DWATR
    use m_sdm_common, only: &
         F_THRD, ONE_PI, dxiv_sdm, dyiv_sdm, VALID2INVALID, num_threads, i2, STAT_LIQ
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    real(RP),intent(in) :: zph_crs(KA,IA,JA)    ! z physical coordinate
    integer, intent(in) :: sd_num               ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA) ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA) ! upper boundary index[k/real] at scalar point in SDM calculation area
    ! Output variables
    real(RP), intent(out) :: rhow_sdm(KA,IA,JA) ! densitiy of liquid water
    real(RP), intent(out) :: liqw_sdm(KA,IA,JA) ! liquid water of super-droplets
    integer, intent(out) :: ilist(1:sd_num)   ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(KA,IA,JA)    ! coef.
    real(RP) :: drate        ! temporary
    integer :: cnt                ! counter
    integer :: i, j, k, kl, ku, m, n    ! index
    integer :: tlist
    integer :: i_threads
    real(RP) :: tmp_liqw_sdm(num_threads,KA,IA,JA)
    !--------------------------------------------------------------------

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_sd2rhow",1,1)
#endif
    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

    ! Initialize
    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
       dcoef(k,i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j) &
            / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
    enddo
    enddo
    enddo

    do k=1,KA
    do j=1,JA
    do i=1,IA
       liqw_sdm(k,i,j) = 0.0_RP
       rhow_sdm(k,i,j) = 0.0_RP
    end do
    end do
    end do

    ! Get index list for compressing buffer.
    cnt =0
    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n)/=STAT_LIQ ) cycle

       cnt = cnt + 1
       ilist(cnt) = n
    end do

    tlist = cnt

    ! Get density of liquid-water.
    !### count voulme of super-droplets ###!

    tmp_liqw_sdm(:,:,:,:) = 0.0_RP

    do i_threads=1,num_threads
    do m=i_threads,tlist,num_threads
       n = ilist(m)

       i = floor(sd_ri(n))+1
       j = floor(sd_rj(n))+1
       k = floor(sd_rk(n))+1

       tmp_liqw_sdm(i_threads,k,i,j) = tmp_liqw_sdm(i_threads,k,i,j)                        &
            &                         + sd_r(n) * sd_r(n) * sd_r(n) * rw            &
            &                         * real(sd_n(n),kind=RP)
    end do
    end do

    do k=1,KA
    do j=1,JA
    do i=1,IA
    do i_threads=1,num_threads
       liqw_sdm(k,i,j)=liqw_sdm(k,i,j)+tmp_liqw_sdm(i_threads,k,i,j)
    end do
    end do
    end do
    end do
    
!!$    do m=1,tlist
!!$       n = ilist(m)
!!$
!!$       i = floor(sd_ri(n))+1
!!$       j = floor(sd_rj(n))+1
!!$       k = floor(sd_rk(n))+1
!!$
!!$       liqw_sdm(k,i,j) = liqw_sdm(k,i,j)                        &
!!$            &                         + sd_r(n) * sd_r(n) * sd_r(n)            &
!!$            &                         * real(sd_n(n),kind=RP)
!!$    end do

    ! Adjust liquid-water in verical boundary.
    do j=JS,JE
    do i=IS,IE
       !! at lower boundary
       kl    = floor(sd_rkl(i,j))+1
       drate = real(kl,kind=RP) - sd_rkl(i,j)
       if( drate<0.50_RP ) then
          liqw_sdm(kl,i,j) = 0.0_RP           !! <50% in share
       else
          liqw_sdm(kl,i,j) = liqw_sdm(kl,i,j)/drate
       end if

       !! at upper boundary
       ku    = floor(sd_rku(i,j))+1
       drate = sd_rku(i,j) - real(ku-1,kind=RP)
       if( drate<0.50_RP ) then
          liqw_sdm(ku,i,j) = 0.0_RP           !! <50% in share
       else
          liqw_sdm(ku,i,j) = liqw_sdm(ku,i,j)/drate
       end if
    end do
    end do

    ! Convert super-droplets to density of liquid-water.
    do k=KS,KE
    do j=JS,JE
    do i=IS,IE
       rhow_sdm(k,i,j) = liqw_sdm(k,i,j) * dcoef(k,i,j)
    end do
    end do
    end do

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_sd2rhow",1,1)
#endif
    return
  end subroutine sdm_sd2rhow
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2rhoisg(sdm_rqi2qs,                        &
       zph_crs,rhoi_sdm,rhos_sdm,rhog_sdm,         &
       sd_num, sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk, &
       sd_rkl,sd_rku,                     &
       soli_sdm,sols_sdm,solg_sdm,ilist_i,ilist_s,ilist_g)
    use scale_grid_index, only: &
         IA,JA,KA,IS,IE,JS,JE,KS,KE
    use m_sdm_common, only: &
         F_THRD, ONE_PI, dxiv_sdm, dyiv_sdm, VALID2INVALID, num_threads, &
         i2, STAT_ICE, sdicedef 
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    real(RP),intent(in) :: sdm_rqi2qs          ! Threshould between qi and qs [m]
    real(RP),intent(in) :: zph_crs(KA,IA,JA)   ! z physical coordinate
    integer, intent(in) :: sd_num              ! Number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num)     ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)     ! y-coordinate of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)    ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA)      ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA)      ! upper boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(out) :: rhoi_sdm(KA,IA,JA)  ! density of cloud ice
    real(RP), intent(out) :: rhos_sdm(KA,IA,JA)  ! density of snow
    real(RP), intent(out) :: rhog_sdm(KA,IA,JA)  ! density of graupel
    real(RP), intent(out) :: soli_sdm(KA,IA,JA)  ! cloud ice of super-droplets
    real(RP), intent(out) :: sols_sdm(KA,IA,JA)  ! snow of super-droplets
    real(RP), intent(out) :: solg_sdm(KA,IA,JA)  ! graupel of super-droplets
    integer, intent(out) :: ilist_i(1:sd_num)  ! buffer for list vectorization
    integer, intent(out) :: ilist_s(1:sd_num)  ! buffer for list vectorization
    integer, intent(out) :: ilist_g(1:sd_num)  ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(KA,IA,JA)   ! coef.
    real(RP) :: drate             ! temporary
    integer :: tlist_i            ! total list number for cloud
    integer :: tlist_s            ! total list number for rain
    integer :: tlist_g            ! total list number for rain
    integer :: icnt               ! counter
    integer :: scnt               ! counter
    integer :: gcnt               ! counter
    integer :: i                  ! index
    integer :: j                  ! index
    integer :: k                  ! index
    integer :: kl                 ! index
    integer :: ku                 ! index
    integer :: m                  ! index
    integer :: n                  ! index
    integer :: i_threads
    real(RP) :: tmp_sol_sdm(num_threads,KA,IA,JA)
    real(RP) :: sd_vol, sd_mass, sd_mrime
    integer  :: sd_nmono
    !-----------------------------------------------------------------------------
    
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_sd2rhoisg",1,1)
#endif
    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)
    
    ! Initialize
    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
       dcoef(k,i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j) &
            / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
    enddo
    enddo
    enddo

    tlist_i = 0
    tlist_s = 0
    tlist_g = 0

    do k=1,KA
    do j=1,JA
    do i=1,IA
       soli_sdm(k,i,j) = 0.0_RP
       sols_sdm(k,i,j) = 0.0_RP
       solg_sdm(k,i,j) = 0.0_RP
       rhoi_sdm(k,i,j) = 0.0_RP
       rhos_sdm(k,i,j) = 0.0_RP
       rhog_sdm(k,i,j) = 0.0_RP
    end do
    end do
    end do

    ! Get index list for compressing buffer.
    icnt = 0
    scnt = 0
    gcnt = 0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n)/=STAT_ICE ) cycle

       sd_vol   = F_THRD * ONE_PI * sdi%re(n)**2 * sdi%rp(n)
       sd_mass  = sd_vol * sdi%rho(n) ! mass of the particle [kg]
       sd_mrime = sdi%mrime(n)
       sd_nmono = sdi%nmono(n)

       if( (sd_mass*0.3d0) < sd_mrime )then
          !### graupel ###!
          gcnt = gcnt + 1
          ilist_g(gcnt) = n
       else if(sd_nmono>10)then
          !### snow flake ###!
          scnt = scnt + 1
          ilist_s(scnt) = n
       else
          !### ice crystal ###!
          icnt = icnt + 1
          ilist_i(icnt) = n
       end if

    end do

    tlist_i = icnt
    tlist_s = scnt
    tlist_g = gcnt

    ! Get density of cloud ice, snow, and graupel.

    !### cloud ice ###!

    if( tlist_i>0 ) then

       tmp_sol_sdm(:,:,:,:) = 0.0_RP

       do i_threads=1,num_threads
       do m=i_threads,tlist_i,num_threads
          n = ilist_i(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          tmp_sol_sdm(i_threads,k,i,j) = tmp_sol_sdm(i_threads,k,i,j)                        &
               &                         + sdi%re(n) * sdi%re(n) * sdi%rp(n) * sdi%rho(n)     &
               &                         * real(sd_n(n),kind=RP)
       end do
       end do

       do k=1,KA
       do j=1,JA
       do i=1,IA
       do i_threads=1,num_threads
          soli_sdm(k,i,j)=soli_sdm(k,i,j)+tmp_sol_sdm(i_threads,k,i,j)
       end do
       end do
       end do
       end do

       !=== adjust at verical boundary. ===!

       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             soli_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             soli_sdm(kl,i,j) = soli_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             soli_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             soli_sdm(ku,i,j) = soli_sdm(ku,i,j)/drate
          end if

       end do
       end do

       !=== convert super-droplets to density of cloud-water. ===!

       do k=KS,KE
       do j=JS,JE
       do i=IS,IE
          rhoi_sdm(k,i,j) = soli_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do

    end if

    !### snow ###!

    if( tlist_s>0 ) then

       tmp_sol_sdm(:,:,:,:) = 0.0_RP

       do i_threads=1,num_threads
       do m=i_threads,tlist_s,num_threads
          n = ilist_s(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          tmp_sol_sdm(i_threads,k,i,j) = tmp_sol_sdm(i_threads,k,i,j)                        &
               &                         + sdi%re(n) * sdi%re(n) * sdi%rp(n) *sdi%rho(n)     &
               &                         * real(sd_n(n),kind=RP)
       end do
       end do

       do k=1,KA
       do j=1,JA
       do i=1,IA
       do i_threads=1,num_threads
          sols_sdm(k,i,j)=sols_sdm(k,i,j)+tmp_sol_sdm(i_threads,k,i,j)
       end do
       end do
       end do
       end do

       !=== adjust at verical boundary. ===!

       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             sols_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             sols_sdm(kl,i,j) = sols_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             sols_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             sols_sdm(ku,i,j) = sols_sdm(ku,i,j)/drate
          end if

       end do
       end do

       !=== convert super-droplets to density of cloud-water. ===!

       do k=KS,KE
       do j=JS,JE
       do i=IS,IE
          rhos_sdm(k,i,j) = sols_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do

    end if

    !### graupel ###!

    if( tlist_g>0 ) then

       tmp_sol_sdm(:,:,:,:) = 0.0_RP

       do i_threads=1,num_threads
       do m=i_threads,tlist_g,num_threads
          n = ilist_g(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          tmp_sol_sdm(i_threads,k,i,j) = tmp_sol_sdm(i_threads,k,i,j)                        &
               &                         + sdi%re(n) * sdi%re(n) * sdi%rp(n) *sdi%rho(n)     &
               &                         * real(sd_n(n),kind=RP)
       end do
       end do

       do k=1,KA
       do j=1,JA
       do i=1,IA
       do i_threads=1,num_threads
          solg_sdm(k,i,j)=solg_sdm(k,i,j)+tmp_sol_sdm(i_threads,k,i,j)
       end do
       end do
       end do
       end do

       !=== adjust at verical boundary. ===!

       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             solg_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             solg_sdm(kl,i,j) = solg_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             solg_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             solg_sdm(ku,i,j) = solg_sdm(ku,i,j)/drate
          end if

       end do
       end do

       !=== convert super-droplets to density of cloud-water. ===!

       do k=KS,KE
       do j=JS,JE
       do i=IS,IE
          rhog_sdm(k,i,j) = solg_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do

    end if

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_sd2rhoisg",1,1)
#endif
    return
  end subroutine sdm_sd2rhoisg
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2rhosol(                        &
       zph_crs,rhosol_sdm,         &
       sd_num, sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk, &
       sd_rkl,sd_rku,                     &
       sol_sdm,ilist)
    use scale_grid_index, only: &
         IA,JA,KA,IS,IE,JS,JE,KS,KE
    use m_sdm_common, only: &
         F_THRD, ONE_PI, dxiv_sdm, dyiv_sdm, VALID2INVALID, num_threads, &
         i2, STAT_ICE, sdicedef 
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    real(RP),intent(in) :: zph_crs(KA,IA,JA)   ! z physical coordinate
    integer, intent(in) :: sd_num              ! Number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num)     ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)     ! y-coordinate of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)    ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA)      ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA)      ! upper boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(out) :: rhosol_sdm(KA,IA,JA)  ! density of cloud ice
    real(RP), intent(out) :: sol_sdm(KA,IA,JA)  ! cloud ice of super-droplets
    integer, intent(out) :: ilist(1:sd_num)  ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(KA,IA,JA)   ! coef.
    real(RP) :: drate             ! temporary
    integer :: tlist              ! total list number for cloud
    integer :: cnt               ! counter
    integer :: i                  ! index
    integer :: j                  ! index
    integer :: k                  ! index
    integer :: kl                 ! index
    integer :: ku                 ! index
    integer :: m                  ! index
    integer :: n                  ! index
    integer :: i_threads
    real(RP) :: tmp_sol_sdm(num_threads,KA,IA,JA)
    !-----------------------------------------------------------------------------
    
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_sd2rhosol",1,1)
#endif
    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)
    
    ! Initialize
    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
       dcoef(k,i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j) &
            / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
    enddo
    enddo
    enddo

    tlist = 0

    do k=1,KA
    do j=1,JA
    do i=1,IA
       sol_sdm(k,i,j) = 0.0_RP
       rhosol_sdm(k,i,j) = 0.0_RP
    end do
    end do
    end do

    ! Get index list for compressing buffer.
    cnt = 0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n)/=STAT_ICE ) cycle

       cnt = cnt + 1
       ilist(cnt) = n

    end do

    tlist = cnt

    ! Get density of solid phase water

    if( tlist>0 ) then

       tmp_sol_sdm(:,:,:,:) = 0.0_RP

       do i_threads=1,num_threads
       do m=i_threads,tlist,num_threads
          n = ilist(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          tmp_sol_sdm(i_threads,k,i,j) = tmp_sol_sdm(i_threads,k,i,j)                        &
               &                         + sdi%re(n) * sdi%re(n) * sdi%rp(n) * sdi%rho(n)     &
               &                         * real(sd_n(n),kind=RP)
       end do
       end do

       do k=1,KA
       do j=1,JA
       do i=1,IA
       do i_threads=1,num_threads
          sol_sdm(k,i,j)=sol_sdm(k,i,j)+tmp_sol_sdm(i_threads,k,i,j)
       end do
       end do
       end do
       end do

       !=== adjust at verical boundary. ===!

       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             sol_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             sol_sdm(kl,i,j) = sol_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             sol_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             sol_sdm(ku,i,j) = sol_sdm(ku,i,j)/drate
          end if

       end do
       end do

       !=== convert super-droplets to density of cloud-water. ===!

       do k=KS,KE
       do j=JS,JE
       do i=IS,IE
          rhosol_sdm(k,i,j) = sol_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do

    end if

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_sd2rhosol",1,1)
#endif
    return
  end subroutine sdm_sd2rhosol
  !-----------------------------------------------------------------------------
end module m_sdm_sd2fluid
