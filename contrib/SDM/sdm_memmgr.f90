!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Memory management of the SDM variables
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
!! @li      2014-06-23 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2014-07-14 (S.Shima) [rev] Removed unused variables
!! @li      2014-07-18 (Y.Sato)  [add] add QTRC_sdm
!! @li      2014-07-22 (Y.Sato)  [mod] modify the bug at the allocation of QTRC_sdm
!! @li      2016-07-11 (S.Shima) [mod] allocate sdice variables
!! @li      2016-07-13 (S.Shima) [mod] Modified for cold SDM
!! @li      2016-07-14 (S.Shima) [mod] Evaluation of bufsiz1, bufsiz2_r8, etc. moved from sdm_allocinit to sdm_numset
!! @li      2016-07-16 (S.Shima) [mod] Dimension of sd_itmp1,sd_itmp2,sd_itmp2 modified.
!! @li      2016-07-16 (S.Shima) [mod] Definition of prr_cres modified.
!! @li      2016-07-21 (S.Shima) [mod] Allocation of the working arrays for restart output moved from scale_atmos_phy_mp_sdm.F90
!! @li      2017-11-30 (S.Shima) [mod] sd_dtmp2
!! @li      2018-02-28 (S.Shima) [add] sd_dtmp3,4,5,6
!! @li      2018-06-30 (S.Shima) [add] rime mass and number of monomers as SD attributes
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_memmgr

  implicit none
  private
  public :: sdm_allocinit

contains
  !-----------------------------------------------------------------------------
  subroutine sdm_allocinit
    use scale_gridtrans, only: &
       I_XYZ, I_XYW,    &
       GTRANS_GSQRT, &
       GTRANS_J13G,  &
       GTRANS_J23G,  &
       GTRANS_J33G
    use scale_tracer_sdm, only: &
         QA_MP
    use scale_grid_index, only: &
         IE,IS,KE,KS,JE,JS,IA,KA,JA ! S:start, E: end of active grids. A: num of grids including HALO.
    use m_sdm_common
    integer :: n, s

    !--- Allocate arrays for SDM
    allocate(sdn_s2c(1:sdnum_s2c))
    allocate(sdri_s2c(1:sdnum_s2c))
    allocate(sdrj_s2c(1:sdnum_s2c))
    allocate(sdrk_s2c(1:sdnum_s2c))
    allocate(sdx_s2c(1:sdnum_s2c))
    allocate(sdy_s2c(1:sdnum_s2c))
    allocate(sdz_s2c(1:sdnum_s2c))
    allocate(sdr_s2c(1:sdnum_s2c))
    allocate(sdu_s2c(1:sdnum_s2c))
    allocate(sdv_s2c(1:sdnum_s2c))
    allocate(sdvz_s2c(1:sdnum_s2c))
    allocate(sdasl_s2c(1:sdnum_s2c,1:sdnumasl_s2c))
    allocate(sdn_fm(1:sdfmnum_s2c))
    allocate(sdri_fm(1:sdfmnum_s2c))
    allocate(sdrj_fm(1:sdfmnum_s2c))
    allocate(sdrk_fm(1:sdfmnum_s2c))
    allocate(sdx_fm(1:sdfmnum_s2c))
    allocate(sdy_fm(1:sdfmnum_s2c))
    allocate(sdz_fm(1:sdfmnum_s2c))
    allocate(sdr_fm(1:sdfmnum_s2c))
    allocate(sdvz_fm(1:sdfmnum_s2c))
    allocate(sdasl_fm(1:sdfmnum_s2c,1:sdnumasl_s2c))

    allocate(sdliqice_s2c(1:sdnum_s2c))
    allocate(sdliqice_fm(1:sdfmnum_s2c))

    if( sdm_cold ) then 
       allocate(sdice_s2c%re(1:sdnum_s2c))
       allocate(sdice_s2c%rp(1:sdnum_s2c))
       allocate(sdice_s2c%rho(1:sdnum_s2c))
!       allocate(sdice_s2c%t(1:sdnum_s2c))
       allocate(sdice_s2c%tf(1:sdnum_s2c))
       allocate(sdice_s2c%mrime(1:sdnum_s2c))
       allocate(sdice_s2c%nmono(1:sdnum_s2c))

       allocate(sdice_fm%re(1:sdfmnum_s2c))
       allocate(sdice_fm%rp(1:sdfmnum_s2c))
       allocate(sdice_fm%rho(1:sdfmnum_s2c))
!       allocate(sdice_fm%t(1:sdfmnum_s2c))
       allocate(sdice_fm%tf(1:sdfmnum_s2c))
       allocate(sdice_fm%mrime(1:sdfmnum_s2c))
       allocate(sdice_fm%nmono(1:sdfmnum_s2c))
    else 
       allocate(sdice_s2c%re(1:1))
       allocate(sdice_s2c%rp(1:1))
       allocate(sdice_s2c%rho(1:1))
!       allocate(sdice_s2c%t(1:1))
       allocate(sdice_s2c%tf(1:1))
       allocate(sdice_s2c%mrime(1:1))
       allocate(sdice_s2c%nmono(1:1))

       allocate(sdice_fm%re(1:1))
       allocate(sdice_fm%rp(1:1))
       allocate(sdice_fm%rho(1:1))
!       allocate(sdice_fm%t(1:1))
       allocate(sdice_fm%tf(1:1))
       allocate(sdice_fm%mrime(1:1))
       allocate(sdice_fm%nmono(1:1))
    end if

    allocate(sdrkl_s2c(IA,JA))
    allocate(sdrku_s2c(IA,JA))

    ! These variables are not needed for SCALE. They are for leap-frog scheme of CReSS.
    allocate(sdn_tmp(1:1))
    allocate(sdrk_tmp(1:1))
    allocate(sdx_tmp(1:1))
    allocate(sdy_tmp(1:1))
    allocate(sdz_tmp(1:1))
    allocate(sdr_tmp(1:1))
    allocate(sdu_tmp(1:1))
    allocate(sdv_tmp(1:1))
    allocate(sdvz_tmp(1:1))
    allocate(sdasl_tmp(1:1,1:1))

    allocate(rand_s2c(1:sdnum_s2c))
    allocate(sortid_s2c(1:sdnum_s2c))
    allocate(sortkey_s2c(1:sdnum_s2c))
    allocate(sortfreq_s2c(1:ni_s2c*nj_s2c*nk_s2c+1))
    allocate(sorttag_s2c(1:ni_s2c*nj_s2c*nk_s2c+2))

    allocate(rbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2))
    allocate(sbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2))
    allocate(rbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2))
    allocate(sbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2))
    allocate(rbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2))
    allocate(sbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2))
    if( sdm_cold ) then
       allocate(rbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2))
       allocate(sbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2))
    else
       allocate(rbuf_i4(1:1,1:1,1:2))
       allocate(sbuf_i4(1:1,1:1,1:2))
    end if

    allocate(sdm_itmp1(1:ni_s2c*nj_s2c*nk_s2c+2))
    allocate(sdm_itmp2(1:ni_s2c*nj_s2c*nk_s2c+2))
    allocate(sdm_itmp3(1:ni_s2c*nj_s2c*nk_s2c+2))

    allocate(sd_itmp1(1:sdnum_s2c))
    allocate(sd_itmp2(1:sdnum_s2c))
    allocate(sd_itmp3(1:sdnum_s2c))

    allocate(sd_dtmp1(1:sdnum_s2c))
    allocate(sd_dtmp2(1:sdnum_s2c))
    allocate(sd_dtmp3(1:sdnum_s2c))
    allocate(sd_dtmp4(1:sdnum_s2c))
    allocate(sd_dtmp5(1:sdnum_s2c))
    allocate(sd_dtmp6(1:sdnum_s2c))

    allocate( prr_crs(IA,JA,1:6) )
    allocate( dxiv_sdm(IA) )
    allocate( dyiv_sdm(JA) )
    allocate( dx_sdm(IA) )
    allocate( dy_sdm(JA) )

    allocate(QTRC_sdm(KA,IA,JA,QA_MP))

    ! Initialize allocated array
    QTRC_sdm(:,:,:,:) = 0.0_RP
    do n = 1, sdnum_s2c
       sortid_s2c(n) = 0
       sortkey_s2c(n) = 0
       rand_s2c(n) = 0.0_RP
       sdn_s2c(n) = int(0,kind=DP)
       sdx_s2c(n) = 0.0_RP
       sdy_s2c(n) = 0.0_RP
       sdz_s2c(n) = 0.0_RP
       sdr_s2c(n) = 0.0_RP
       sdu_s2c(n) = 0.0_RP
       sdv_s2c(n) = 0.0_RP
       sdvz_s2c(n) = 0.0_RP

       sd_itmp1(n) = 0
       sd_itmp2(n) = 0
       sd_itmp3(n) = 0

       sd_dtmp1(n) = 0.0_RP
       sd_dtmp2(n) = 0.0_RP
    enddo

    do s = 1, sdnumasl_s2c
       do n = 1, sdnum_s2c
          sdasl_s2c(n,s) = 0.0_RP
       enddo
    enddo

    if( abs(sdm_aslset) > 10 ) then
       do n = 1, sdfmnum_s2c
         sdn_fm(n) = int(0,kind=DP)
         sdx_fm(n) = 0.0_RP
         sdy_fm(n) = 0.0_RP
         sdz_fm(n) = 0.0_RP
         sdr_fm(n) = 0.0_RP
         sdvz_fm(n) = 0.0_RP
       enddo
       do s = 1, sdnumasl_s2c
            do n = 1, sdnum_s2c
              sdasl_fm(n,s) = 0.0_RP
            enddo
       enddo
    endif

    do n = 1, ni_s2c*nj_s2c*nk_s2c+1
       sortfreq_s2c(n) = 0
    enddo

    do n = 1, ni_s2c*nj_s2c*nk_s2c+2
       sorttag_s2c(n) = 0
       sdm_itmp1(n) = 0
       sdm_itmp2(n) = 0
       sdm_itmp3(n) = 0
    enddo

    allocate(sdn_s2c_restart(1:sdnum_s2c))
    allocate(sdrk_s2c_restart(1:sdnum_s2c))
    allocate(sdx_s2c_restart(1:sdnum_s2c))
    allocate(sdy_s2c_restart(1:sdnum_s2c))
    allocate(sdz_s2c_restart(1:sdnum_s2c))
    allocate(sdr_s2c_restart(1:sdnum_s2c))
    allocate(sdu_s2c_restart(1:sdnum_s2c))
    allocate(sdv_s2c_restart(1:sdnum_s2c))
    allocate(sdvz_s2c_restart(1:sdnum_s2c))
    allocate(sdasl_s2c_restart(1:sdnum_s2c,1:sdnumasl_s2c))
    allocate(sdliqice_s2c_restart(1:sdnum_s2c))
    if( sdm_cold ) then
       allocate(sdice_s2c_restart%re(1:sdnum_s2c))
       allocate(sdice_s2c_restart%rp(1:sdnum_s2c))
       allocate(sdice_s2c_restart%rho(1:sdnum_s2c))
       allocate(sdice_s2c_restart%tf(1:sdnum_s2c))
       allocate(sdice_s2c_restart%mrime(1:sdnum_s2c))
       allocate(sdice_s2c_restart%nmono(1:sdnum_s2c))
    else
       allocate(sdice_s2c_restart%re(1:1))
       allocate(sdice_s2c_restart%rp(1:1))
       allocate(sdice_s2c_restart%rho(1:1))
       allocate(sdice_s2c_restart%tf(1:1))
       allocate(sdice_s2c_restart%mrime(1:1))
       allocate(sdice_s2c_restart%nmono(1:1))
    end if

    return

   end subroutine sdm_allocinit
end module m_sdm_memmgr
