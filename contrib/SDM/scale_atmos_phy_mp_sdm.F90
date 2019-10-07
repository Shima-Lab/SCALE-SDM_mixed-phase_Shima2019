!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics by Super Droplet Method (SDM)
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
!! @li      2013-09-30 (Y.Sato)  [new] Implement from Original version of SDM
!! @li      2014-01-22 (Y.Sato)  [rev] Update for scale-0.0.0
!! @li      2014-05-04 (Y.Sato)  [rev] Update for scale-0.0.1
!! @li      2014-06-06 (S.Shima) [rev] Modify several bug 
!! @li      2014-06-07 (Y.Sato)  [rev] Remove dt=max(dt_sdm) and some other changes
!! @li      2014-06-09 (S.Shima) [rev] Check whether dt is the least common multiple of sdm_dtcmph(i) 
!!                                     that satisfies sdm_dtcmph(i)<= dt. Fixed a hidden bug: "/=" to "==".
!! @li      2014-06-13 (S.Shima) [rev] Common variables are separated into sdm_common.f90
!! @li      2014-06-14 (S.Shima) [rev] Check the initialization of the random number generator.
!! @li      2014-06-24 (S.Shima) [rev] Separated sdm_allocinit
!! @li      2014-06-25 (S.Shima) [rev] Bugfix and improvement of ATMOS_PHY_MP_sdm_restart_in and sdm_iniset
!! @li      2014-06-25 (S.Shima) [rev] Bugfix of sd position initialization, and many modification
!! @li      2014-06-25 (S.Shima) [rev] Bugfix: dx_sdm, dy_sdm, dxiv_sdn, dyiv_sdm restored
!! @li      2014-06-26 (S.Shima) [rev] sdm_getrklu and sdm_z2rk are separated
!! @li      2014-06-27 (S.Shima) [rev] sd data output functionality added
!! @li      2014-07-04 (S.Shima) [rev] Removed comment outputs for debugging
!! @li      2014-07-09 (S.Shima) [rev] Subroutines related to boundary conditions and MPI communications are totally revised
!! @li      2014-07-11 (S.Shima) [rev] Subroutines related to sdm_getvz are revised. Many bug fixes.
!! @li      2014-07-11 (S.Shima) [rev] Subroutines for conversion between fluid variables are separated into module m_sdm_fluidconv
!! @li      2014-07-11 (S.Shima) [rev] Subroutines to impose boundary conditions are separated into module m_sdm_boundary
!! @li      2014-07-11 (S.Shima) [rev] Motion (advection/sedimentation/precipitation) related subroutines are separated into the 
!!                                     module m_sdm_motion
!! @li      2014-07-12 (S.Shima) [rev] Add comments concerning when to diagnose QC and QR
!! @li      2014-07-12 (S.Shima) [rev] BUG of random number initialization removed
!! @li      2014-07-12 (S.Shima) [rev] sdm_sort repaired
!! @li      2014-07-12 (S.Shima) [rev] sdm_coales repaired
!! @li      2014-07-12 (S.Shima) [rev] sdm_coales tuned for FX (use compile option -Kprefetch_indirect)
!! @li      2014-07-12 (S.Shima) [rev] sdm_sort and sdm_getperm are separated into the module m_sdm_idutil
!! @li      2014-07-12 (S.Shima) [rev] sdm_coales is separated into the module m_sdm_coalescence
!! @li      2014-07-13 (S.Shima) [rev] sdm_condevp repaired
!! @li      2014-07-13 (S.Shima) [rev] sdm_sd2qcqr, sdm_sd2rhocr repaired
!! @li      2014-07-14 (S.Shima) [rev] sdm_condevp_updatefluid created modifying and repairing sdm_sd2qvptp
!! @li      2014-07-14 (S.Shima) [rev] Update HALO after call sdm_condevp_updatefluid
!! @li      2014-07-14 (S.Shima) [rev] sdm_condevp tuned for FX (use compile option -Kocl -Kprefetch_indirect)
!! @li      2014-07-14 (S.Shima) [rev] sdm_sd2rhow, sdm_sd2rhocr, sdm_sd2qcqr are separated into m_sdm_sd2fluid
!! @li      2014-07-14 (S.Shima) [rev] sdm_condevp,sdm_condevp_updatefluid are separated into m_sdm_condensation_water
!! @li      2014-07-14 (S.Shima) [rev] Removed unused subroutines, variables
!! @li      2014-07-18 (Y.Sato)  [mod] Modify the modules to calculate variable forradiation
!! @li      2014-07-18 (Y.Sato)  [add] Add history output for QTRC_sdm and input 0 to QTRC at the end of ATMOS_PHY_MP_sdm
!! @li      2014-07-18 (Y.Sato)  [mod] Modify a way to determine whether the grid stretch or not
!! @li      2014-07-18 (Y.Sato)  [add] Add ATMOS_PHY_MP_DENS, which is used in radiation code
!! @li      2014-07-22 (Y.Sato)  [mod] Modify the way to determine whether the grid stretch or not
!! @li      2014-07-22 (Y.Sato)  [mod] Modify the way to determine whether the lateral boundary is periodic or not
!! @li      2014-07-22 (Y.Sato)  [mod] Modify to write error information to normal output
!! @li      2014-07-22 (Y.Sato)  [mod] Modify bugs to calculate variable for radiation process
!! @li      2014-07-22 (Y.Sato)  [mod] Modify timing for assigning to QTRC_sdm and clear QTRC(2:QA)
!! @li      2014-07-22 (Y.Sato)  [mod] Modify the definition of kl and ku for calculating drate
!! @li      2014-07-22 (Y.Sato)  [mod] Modify the order of loop in a part
!! @li      2014-07-24 (Y.Sato)  [mod] Modify a bug for restart
!! @li      2014-07-25 (Y.Sato)  [rev] Move sdm_getrklu from sdm_iniset to ATMOS_PHY_MP_sdm_setup
!! @li      2014-07-25 (Y.Sato)  [rev] Add COMM_var, and COMM_wait for filling u_scale, v_scale, and w_scale
!! @li      2014-12-12 (Y.Sato)  [mod] Modify for using QTRC_sdm in sdm_sd2qcqr in sdm_iniset 
!! @li      2014-12-17 (Y.Sato)  [mod] Add initialization of prr_crs for Restart
!! @li      2015-06-22 (S.Shima) [add] Add section specification call for profiling (fipp and fapp)
!! @li      2015-06-27 (S.Shima) [add] Add more section specification call for profiling (fapp)
!! @li      2015-06-27 (S.Shima) [add] Store environmental variable PARALELL to num_threads 
!! @li      2015-07-30 (Y.Sato)  [add] Add "ifdef" for fapp and fipp module
!! @li      2015-09-08 (Y.Sato)  [mod] update for version SCALE 0.4.2
!! @li      2015-09-15 (Y.Sato)  [mod] update for version SCALE 0.4.3
!! @li      2016-07-11 (S.Shima) [mod] Initialization of sdliqice. Arguments of outascii.
!! @li      2016-07-13 (S.Shima) [mod] Modified for MPI communication of ice particles
!! @li      2016-07-14 (S.Shima) [add] Initialization of ice attributes
!! @li      2016-07-15 (S.Shima) [add] sdm_meltfreeze added
!! @li      2016-07-16 (S.Shima) [add] sdm_sd2qiqsqg added
!! @li      2016-07-16 (S.Shima) [mod] Dimension of sd_itmp1-3 modified
!! @li      2016-07-16 (S.Shima) [mod] sdm_sd2qcqr, sdm_sd2rhocr, sdm_sd2rhow modified
!! @li      2016-07-16 (S.Shima) [mod] sdm_sd2prec modified
!! @li      2016-07-17 (S.Shima) [add] Iteration of sdm_meltfreeze added
!! @li      2016-07-17 (S.Shima) [add] Heat exchange through melt/freeze added
!! @li      2016-07-18 (S.Shima) [mod] liqice argument is added to sdm_condevp because now it is only for liquid water phase particles
!! @li      2016-07-19 (S.Shima) [add] sdm_subldep.f90 added
!! @li      2016-07-20 (S.Shima) [mod] Temporal implementation of the prototype of ice particle advection
!! @li      2016-07-20 (S.Shima) [mod] Temporal implementation of the prototype of ice particle coalescence
!! @li      2016-07-21 (S.Shima) [fix] Small bug on when to copy SD data to restart array fixed
!! @li      2016-07-21 (S.Shima) [mod] Restart in/out are modified for supporting cold SDM
!! @li      2016-07-21 (S.Shima) [mod] Working array for restart output moved to sdm_common.f90 and sdm_memmgr.f90
!! @li      2017-02-08 (S.Shima) [mod] sdm_getvz is renamed to sdm_getvz_liq
!! @li      2017-02-08 (S.Shima) [add] sdm_getvz_ice
!! @li      2017-02-15 (S.Shima) [mod] for new sdm_subldep subroutine
!! @li      2017-02-15 (S.Shima) [mod] for the list vector of sdm_subldep subroutine
!! @li      2017-02-19 (S.Shima) [add] sdm_coales_cold
!! @li      2017-08-24 (S.Shima) [mod] argument of sdm_meltfreeze 
!! @li      2017-08-28 (S.Shima) [add] INAS density of Niemand et al. (2012) is implemented
!! @li      2017-09-27 (S.Shima) [mod] argument of sdm_coalescence_cold
!! @li      2017-10-01 (S.Shima) [mod] argument of sdm_subldep
!! @li      2017-11-30 (S.Shima) [mod] multiplicity of IN to sample rare events 
!! @li      2018-02-28 (S.Shima) [mod] predictor-corrector is used for SD motion eq
!! @li      2018-04-03 (S.Shima) [mod] argument of sdm_coales_cold
!! @li      2018-05-13 (S.Shima) [fix] latent heat release through riming
!! @li      2018-06-25 (S.Shima) [add] netcdf output
!! @li      2018-06-30 (S.Shima) [add] rime mass and number of monomers as SD attributes
!! @li      2019-01-09 (S.Shima) [add] sdm_fctr2multi
!! @li      2019-01-09 (S.Shima) [fix] argument declaration: sdm_calvar(3) -> sdm_calvar(5)
!! @li      2019-01-10 (S.Shima) [mod] no interpolation of scalar variables when evaluating terminal velocity
!! @li      2019-01-10 (S.Shima) [mod] the order of calculating cloud microphysical processes
!! @li      2019-01-12 (S.Shima) [add] momentum coupling (rhow)
!! @li      2019-01-12 (S.Shima) [mod] dry air density -> moist air density
!<  
!-------------------------------------------------------------------------------
#include "macro_thermodyn.h"
module scale_atmos_phy_mp_sdm
  !-----------------------------------------------------------------------------
  !
  !++ used modules ! For encapsulation, reduce the use of modules here as far as possible. 
  !   Modules should be called inside the subroutine here.
  !
  use mpi
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer_sdm
  use scale_process, only:  &
     mype => PRC_myrank,    &
     PRC_nprocs,            &
     PRC_IsMaster,          &
     PRC_MPIstop
  use scale_time, only:     &
     TIME_NOWSEC
  use gadg_algorithm, only: &
     gadg_count_sort
  use rng_uniform_mt, only: &
     c_rng_uniform_mt,      &
     rng_save_state,        &
     rng_load_state,        &
     gen_rand_array => rng_generate_array
  use m_sdm_common
  use m_sdm_numset
  use m_sdm_memmgr
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  public :: ATMOS_PHY_MP_sdm_setup
  public :: ATMOS_PHY_MP_sdm
  public :: ATMOS_PHY_MP_sdm_CloudFraction
  public :: ATMOS_PHY_MP_sdm_EffectiveRadius
  public :: ATMOS_PHY_MP_sdm_MixingRatio
  public :: ATMOS_PHY_MP_sdm_restart_out
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  real(RP), public, target :: ATMOS_PHY_MP_DENS(MP_QA) ! hydrometeor density [kg/m3]=[g/L]
  logical,  public, save   :: sd_rest_flg_out = .false. ! restart flg of Super Droplet

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_setup( MP_TYPE )
    use scale_les_process, only: &
       PRC_next,    &
       PRC_NUM_X,   &
       PRC_NUM_Y,   &
       PRC_W,       &
       PRC_E,       &
       PRC_S,       &
       PRC_N
    use scale_const, only: &
       PI     => CONST_PI,    &
       GRAV   => CONST_GRAV,  &
       dens_w => CONST_DWATR, &
       dens_i => CONST_DICE
!    use scale_specfunc, only: &
!       SF_gamma
!    use scale_comm, only: &
!       COMM_horizontal_mean
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_grid_real, only: & ! These should be all REAL_xx. Check later.
         REAL_FZ
    use scale_grid, only: &
         GRID_CDX,   &
         GRID_CDY,   &
         GRID_RCDX,  &
         GRID_RCDY,  &
         GRID_CDZ,   &
         GRID_FX,    &
         GRID_FY,    &
         DX,DY,DZ
!!$       CZ  => GRID_CZ,    &
!!$       FZ  => GRID_FZ,    &
!!$       FDX => GRID_FDX,   &
!!$       FDY => GRID_FDY,   &
!!$       FDZ => GRID_FDZ,   & 
!!       CBFZ => GRID_CBFZ, &
!!       CBFX => GRID_CBFX, &
!!       CBFY => GRID_CBFY, &
!!       ISG, IEG, JSG, JEG,&
    use scale_topography, only: &
       TOPO_exist
    use m_sdm_coordtrans, only: &
       sdm_getrklu
    implicit none
    character(len=H_SHORT), intent(in) :: MP_TYPE

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       docondensation, &
       doautoconversion, &
       doprecipitation

    NAMELIST / PARAM_ATMOS_PHY_MP_SDM / &
       RANDOM_IN_BASENAME,  &
       RANDOM_OUT_BASENAME, &
       SD_IN_BASENAME,      &
       SD_OUT_BASENAME,     &
       domovement,          &
       sdm_cold,            &
       domeltfreeze,   &
       dosublimation,   &
       donegative_fixer,    &
       sdm_dtcmph,          &
       sdm_rdnc,            &
       sdm_sdnmlvol,        &
       sdm_inisdnc,         &
       sdm_aslset,          &
       sdm_fctr2multi,      &
       sdm_aslmw,           &         
       sdm_zlower,          &
       sdm_zupper,          &
       sdm_extbuf,          &
       sdm_rqc2qr,          &
       sdm_rqi2qs,          &
       sdm_aslfmrate,       &
       sdm_aslfmdt,         &
       sdm_aslfmsdnc,       &
       sdm_colbrwn,         &
       sdm_colkrnl,         &
       sdm_mvexchg,         &
       sdm_nadjdt,          &
       sdm_nadjvar,         &
       sdm_dmpvar,          &
       sdm_dmpitva,         &
       sdm_dmpnskip,        & 
       sdm_dmpitvb,         & 
       sdm_dmpitvl,         &
       sdm_dmpsdsiz

!    real(RP) :: dtevl
!    real(RP) :: n0, dry_r
!    real(RP) :: delta1, delta2 !, sdn_tmp
    real(RP) :: buffact
    integer :: ierr
    integer :: i, j, ip, k, n, s
    integer :: bndsdmdim, bufsiz
    integer :: tmppe, ilcm, igcd
    character(len=17) :: fmt1="(A, '.', A, I*.*)"
    character(:), allocatable :: tmp_value
    integer len, status
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Super Droplet Method'

    if ( MP_TYPE /= 'SDM' ) then
       write(*,*) 'xxx ATMOS_TYPE_PHY_MP is not SDM. Check!'
       call PRC_MPIstop
    endif

    ! Get the number of threads of Auto Parallelization on K/FX10
    call get_environment_variable('PARALLEL', status = status, length = len)
    if (status /= 0) then
       num_threads=1
    else
       allocate (character(len) :: tmp_value)
       call get_environment_variable('PARALLEL', value = tmp_value)
       read(tmp_value,*) num_threads
       deallocate(tmp_value)
    end if

    ! [Bug] Stretched coordinate is not supported yet, but the check below is not working any more
    ! [Bug] Must be fixed in the future when introducing generalized coordinate 
    buffact = 0.0_RP
    do k = KS, KE
      buffact = max( buffact,GRID_CDZ(k)/DZ )
    enddo
    do j = JS, JE
      buffact = max( buffact,GRID_CDY(j)/DY )
    enddo
    do i = IS, IE
      buffact = max( buffact,GRID_CDX(i)/DX )
    enddo

    if( buffact < 1.0_RP .or. buffact > 1.0_RP ) then
       sthopt = 1
    else
       sthopt = 0
    endif

    if(sthopt==1) then
       write(*,*) 'ERROR: stretched coordinate is not yet supported!', buffact
       call PRC_MPIstop
    end if
    ! [Bug] Stretched coordinate is not supported yet, but the check above is not working any more

    if( TOPO_exist ) then
       trnopt = 2
    elseif( .not. TOPO_exist ) then
       trnopt = 0
    endif

    if(trnopt==2) then
       write(*,*) 'ERROR: terrain following coordinate is not yet supported!'
       call PRC_MPIstop
    end if

    !--- if the lateral boundary is not periodic, PRC_next has negative value
    if( minval( PRC_next,1 ) < 0 ) then
       write(*,*) 'ERROR: Only periodic B.C. is supported!'
       write(*,*) 'ERROR: Set PRC_PERIODIC_X=PRC_PERIODIC_Y=.true.'
       call PRC_MPIstop
    else
       wbc=1
       ebc=1
       sbc=1
       nbc=1
    end if

    nsub  = max( PRC_nprocs,1 )
    nisub = max( PRC_NUM_X,1 )
    njsub = max( PRC_NUM_Y,1 )

    !--- set process next to the process
    dstw_sub = PRC_next(PRC_W)
    dste_sub = PRC_next(PRC_E)
    srcw_sub = PRC_next(PRC_W)
    srce_sub = PRC_next(PRC_E)

    dsts_sub = PRC_next(PRC_S)
    dstn_sub = PRC_next(PRC_N)
    srcs_sub = PRC_next(PRC_S)
    srcn_sub = PRC_next(PRC_N)

    tag  = 0

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_SDM,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP_SDM. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP_SDM)

    if( sdm_cold == .false. ) then 
       if( IO_L ) write(IO_FID_LOG,*) 'Warm SDM is used'
    else if ( sdm_cold == .true. ) then
       if( IO_L ) write(IO_FID_LOG,*) 'Cold SDM is used'
    endif

    if( (sdm_cold == .false.) .and. (domeltfreeze == .true.) ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Warning: domeltfreeze option does not work for warm SDM'
    endif

    if( RANDOM_IN_BASENAME == '' ) then
       write(*,*) 'xxx Set Random number file! stop'
       write(*,*) 'To generate random number set, run scale_init with'
       write(*,*) 'flg_sdm = .true. in PARAM_MKINIT of init.conf'
       call PRC_MPIstop
    else
       write(fmt1(14:14),'(I1)') 6
       write(fmt1(16:16),'(I1)') 6
       write(RANDOM_IN_BASENAME,fmt1) trim(RANDOM_IN_BASENAME),'pe',mype
    endif
    fid_random_i = IO_get_available_fid()

    if( SD_IN_BASENAME == '' ) then
       sd_rest_flg_in = .false.
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found S.D. file, generate from random number.'
    else
       sd_rest_flg_in = .true.
       fid_sd_i = IO_get_available_fid()
       write(fmt1(14:14),'(I1)') 6
       write(fmt1(16:16),'(I1)') 6
       write(SD_IN_BASENAME,fmt1) trim(SD_IN_BASENAME),'pe',mype
    endif

    if( SD_OUT_BASENAME == '' ) then
       sd_rest_flg_out = .false.
    else
       sd_rest_flg_out = .true.
    endif

    if( docondensation )   sdm_calvar(1) = .true.
    if( doautoconversion ) sdm_calvar(2) = .true.
    if( domovement )       sdm_calvar(3) = .true.
    if( domeltfreeze )     sdm_calvar(4) = .true.
    if( dosublimation )    sdm_calvar(5) = .true.

! zlower+surface height is the lower boundary of SDs.
!!$     if( sdm_zlower < CZ(KS) ) then
!!$      if( mype == PRC_master )  write(*,*) "sdm_zlower was set to CZ(KS) because zlower < CZ(KS)"
!!$      sdm_zlower = CZ(KS)
!!$     endif
! 

!     if( sdm_zupper > CZ(KE) ) then
!      if( mype == PRC_master )  write(*,*) "sdm_zupper was set to CZ(KE) because zupper > CZ(KE)"
!      sdm_zupper = CZ(KE)
!     endif

    if( sdm_zupper > minval(REAL_FZ(KE-1,IS:IE,JS:JE)) ) then
!       if( mype == PRC_master )  write(*,*) "sdm_zupper was set to minval(REAL_FZ(KE-1)) because zupper > minval(REAL_FZ(KE-1))"
       if( PRC_IsMaster )  write(*,*) "sdm_zupper was set to minval(REAL_FZ(KE-1)) because zupper > minval(REAL_FZ(KE-1))"
       sdm_zupper = minval(REAL_FZ(KE-1,IS:IE,JS:JE))
    endif

    ! check whether sdm_dtcmph(1:3) > 0
     if(  ( (sdm_dtcmph(1) <= 0.0_RP) .and. docondensation   )   .or. &
         ( (sdm_dtcmph(2) <= 0.0_RP) .and. doautoconversion )   .or. &
         ( (sdm_dtcmph(3) <= 0.0_RP) .and. domovement       )   .or. &
         ( (sdm_dtcmph(4) <= 0.0_RP) .and. domeltfreeze     )   .or. &
         ( (sdm_dtcmph(5) <= 0.0_RP) .and. dosublimation    )        ) then
       write(*,*) 'ERROR: sdm_dtcmph(1:5) have to be positive'
       call PRC_MPIstop
     end if

    ! check whether dt (dt of mp) is larger than sdm_dtcmph(i).
    if( (dt < minval(sdm_dtcmph(:))) ) then
       write(*,*) 'ERROR: For now, sdm_dtcmph should be smaller than TIME_DTSEC_ATMOS_PHY_MP'
       call PRC_MPIstop
    end if

    ! aerosol nucleation and sd number adjustment functions are not supported yet
    if ( (abs(sdm_aslset) >= 10) .or. (sdm_nadjvar /= 0)) then
       write(*,*) 'ERROR: aerosol nucleation and sd number adjustment functions are not supported yet'
       write(*,*) 'ERROR: set sdm_aslset < 10 and sdm_nadjvar =0'
       call PRC_MPIstop
    end if

    ! rigorous momentum exchange function is not supported yet
    if ( sdm_mvexchg >= 2) then
       write(*,*) 'ERROR: Rigorous momentum exchange not yet supported. set sdm_mvexchg = 0 (none) or 1 (rhow only)'
       call PRC_MPIstop
    end if

    if( docondensation ) then
       nclstp(1)=10*int(1.E+5_RP*(dt+1.E-6_RP))            &
            /int(1.E+6_RP*(sdm_dtcmph(1)+1.E-7_RP))
       if(mod(10*int(1.E+5_RP*(dt+1.E-6_RP)),int(1.E+6_RP*(sdm_dtcmph(1)+1.E-7_RP))) /= 0) then
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_MPIstop
       end if
    else
       nclstp(1) = 1
    end if

    if( doautoconversion ) then
       nclstp(2)=10*int(1.E+5_RP*(dt+1.E-6_RP))            &
            /int(1.E+6_RP*(sdm_dtcmph(2)+1.E-7_RP))
       if(mod(10*int(1.E+5_RP*(dt+1.E-6_RP)),int(1.E+6_RP*(sdm_dtcmph(2)+1.E-7_RP))) /= 0) then
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_MPIstop
       end if
    else
       nclstp(2) = 1
    end if

    if( domovement ) then
       nclstp(3)=10*int(1.E+5_RP*(dt+1.E-6_RP))            &
            /int(1.E+6_RP*(sdm_dtcmph(3)+1.E-7_RP))
       if(mod(10*int(1.E+5_RP*(dt+1.E-6_RP)),int(1.E+6_RP*(sdm_dtcmph(3)+1.E-7_RP))) /= 0) then
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_MPIstop
       end if
    else
       nclstp(3) = 1
    end if

    if( domeltfreeze ) then
       nclstp(4)=10*int(1.E+5_RP*(dt+1.E-6_RP))            &
            /int(1.E+6_RP*(sdm_dtcmph(4)+1.E-7_RP))
       if(mod(10*int(1.E+5_RP*(dt+1.E-6_RP)),int(1.E+6_RP*(sdm_dtcmph(4)+1.E-7_RP))) /= 0) then
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_MPIstop
       end if
    else
       nclstp(4) = 1
    end if

    if( dosublimation ) then
       nclstp(5)=10*int(1.E+5_RP*(dt+1.E-6_RP))            &
            /int(1.E+6_RP*(sdm_dtcmph(5)+1.E-7_RP))
       if(mod(10*int(1.E+5_RP*(dt+1.E-6_RP)),int(1.E+6_RP*(sdm_dtcmph(5)+1.E-7_RP))) /= 0) then
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_MPIstop
       end if
    else
       nclstp(5) = 1
    end if

    nclstp(0)=min(nclstp(1),nclstp(2),nclstp(3),nclstp(4),nclstp(5))

    ilcm=0
    iterate_2: do
       ilcm=ilcm+1
       nclstp(0)=nclstp(0)*ilcm
       if( mod(nclstp(0),nclstp(1)) == 0 .and.           &
           mod(nclstp(0),nclstp(2)) == 0 .and.           &
           mod(nclstp(0),nclstp(3)) == 0 .and.           &
           mod(nclstp(0),nclstp(4)) == 0 .and.           &
           mod(nclstp(0),nclstp(5)) == 0   ) then
           exit iterate_2
       else
          nclstp(0)=nclstp(0)/ilcm
       end if
    end do iterate_2

    ! check whether dt is the least common multiple of sdm_dtcmph(i) that satisfies sdm_dtcmph(i)<= dt 
    !! find the smallest nclstp(1:3) that satisfies sdm_dtcmph(i)<= dt 
    igcd=maxval(nclstp(1:5))
    if((sdm_dtcmph(1)<= dt).and.docondensation)   igcd=min(nclstp(1),igcd)
    if((sdm_dtcmph(2)<= dt).and.doautoconversion) igcd=min(nclstp(2),igcd)
    if((sdm_dtcmph(3)<= dt).and.domovement)       igcd=min(nclstp(3),igcd)
    if((sdm_dtcmph(4)<= dt).and.domeltfreeze)     igcd=min(nclstp(4),igcd)
    if((sdm_dtcmph(5)<= dt).and.dosublimation)    igcd=min(nclstp(5),igcd)
    !! find the greatest common divisor of nclstp(1:5) that satisfies sdm_dtcmph(i)<= dt 
    do 
       if( igcd == 1) exit
       if( mod(nclstp(1),igcd) == 0 .and.           &
           mod(nclstp(2),igcd) == 0 .and.           &
           mod(nclstp(3),igcd) == 0 .and.           &
           mod(nclstp(4),igcd) == 0 .and.           &
           mod(nclstp(5),igcd) == 0   ) then
           exit
       end if
       igcd = igcd - 1
    end do
    !! if igcd>1, it meanst dt is not the least common multiple of sdm_dtcmph(i) that satisfies sdm_dtcmph(i)<= dt
    if( igcd > 1)then
       write(*,*) 'ERROR: TIME_DTSEC_ATMOS_PHY_MP should be the least comon multiple of sdm_dtcmph(1:5)', &
                  ' that are smaller than TIME_DTSEC_ATMOS_PHY_MP'
       call PRC_MPIstop
    end if

!!$    dx_sdm(1:IA) = CDX(1:IA)
!!$    dy_sdm(1:JA) = CDY(1:JA)
!!$    dz_sdm(1:KA) = CDZ(1:KA)
!!$    dx_sdm(1:IA) = FDX(1:IA)
!!$    dy_sdm(1:JA) = FDY(1:JA)
!!$    dz_sdm(1:KA) = FDZ(1:KA)
    xmax_sdm = GRID_FX(IE)-GRID_FX(IS-1)
    ymax_sdm = GRID_FY(JE)-GRID_FY(JS-1)

    allocate( zph_crs(KA,IA,JA) )

    !--- set number of super droplet etc...
    call sdm_numset(              &
      sdm_extbuf,                 &
      sdm_aslset, sdm_sdnmlvol,   &
      sdm_inisdnc,                &
      sdm_aslfmsdnc,              &
      sdm_zlower, sdm_zupper,     &
      minzph, sdininum_s2c,       &
      sdfmnum_s2c, sdnum_s2c,     &
      ni_s2c, nj_s2c, nk_s2c, zph_crs )

    call sdm_allocinit

    !### Get index[k/real] at "sdm_zlower" and "sdm_zupper"  ###!
    call sdm_getrklu(sdm_zlower,sdm_zupper,      &
                     sdrkl_s2c,sdrku_s2c)

    dx_sdm(1:IA) = GRID_CDX(1:IA)
    dy_sdm(1:JA) = GRID_CDY(1:JA)
    dxiv_sdm(1:IA) = GRID_RCDX(1:IA)
    dyiv_sdm(1:JA) = GRID_RCDY(1:JA)

    ATMOS_PHY_MP_DENS(I_mp_QC) = dens_w ! hydrometeor density [kg/m3]=[g/L]

    return
  end subroutine ATMOS_PHY_MP_sdm_setup
  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm( &
       DENS,      &
       MOMZ,      &
       MOMX,      &
       MOMY,      &
       RHOT,      &
       QTRC,      &
       CCN,       &
       EVAPORATE, &
       SFLX_rain, &
       SFLX_snow  )
    use scale_grid_index
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP !,&
    use scale_history, only: &
       HIST_in
    use scale_atmos_phy_mp_common, only: &
       MP_negative_fixer        => ATMOS_PHY_MP_negative_fixer,       &
       MP_precipitation         => ATMOS_PHY_MP_precipitation,        &
       MP_saturation_adjustment => ATMOS_PHY_MP_saturation_adjustment
    use scale_atmos_thermodyn, only: &
       CPw => AQ_CP, &
       THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,       &
       THERMODYN_rhot        => ATMOS_THERMODYN_rhot,       &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_gridtrans, only: &
       I_XYZ, I_XYW,    &
       GTRANS_GSQRT, &
       GTRANS_J13G,  &
       GTRANS_J23G,  &
       GTRANS_J33G
    use scale_const, only: &
       rw => CONST_DWATR, &
       CPdry => CONST_CPdry, &
       Rdry  => CONST_Rdry, &
       Rvap  => CONST_Rvap, &
       P00   => CONST_PRE00

    use m_sdm_io, only: &
       sdm_outasci,sdm_outnetcdf
    use m_sdm_coordtrans, only: &
       sdm_rk2z
    use m_sdm_fluidconv, only: &
       sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw
    use m_sdm_motion, only: &
       sdm_getvz_liq, sdm_getvz_ice
    use m_sdm_sd2fluid, only: &
         sdm_sd2qcqr,sdm_sd2qiqsqg
    implicit none
    real(RP), intent(inout) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)        !! Momentum [kg/s/m2]
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)    !! Ratio of mass of tracer to total mass[kg/kg]
    real(RP), intent(out)   :: EVAPORATE(KA,IA,JA)   !---- evaporated cloud number concentration [/m3]
    real(RP), intent(in)    :: CCN(KA,IA,JA)         !! cloud condensation nuclei calculated by aerosol module
    real(RP), intent(out)   :: SFLX_rain(IA,JA)      !! rain flux at surface (used in land and urban model)
    real(RP), intent(out)   :: SFLX_snow(IA,JA)      !! snow flux at surface (used in land and urban model)

    real(RP) :: rtmp4(KA,IA,JA)    ! Temporary array
    real(RP) :: rtmp5(KA,IA,JA)    ! Temporary array
    real(RP) :: rtmp6(KA,IA,JA)    ! Temporary array
    ! Work variables
    logical :: lsdmup       ! flag for updating water hydrometeor by SDM
    real(RP) :: sd_nc  ! averaged number concentration in a grid

    real(RP) :: flux_tot (KA,IA,JA)
    real(RP) :: flux_rain(KA,IA,JA)
    real(RP) :: flux_snow(KA,IA,JA)
    real(RP) :: crs_dtmp1(KA,IA,JA), crs_dtmp2(KA,IA,JA)
    real(RP) :: crs_dtmp3(KA,IA,JA), crs_dtmp4(KA,IA,JA)
    real(RP) :: crs_dtmp5(KA,IA,JA), crs_dtmp6(KA,IA,JA)
    integer  :: n, s, k, i, j, iq         ! index

    real(RP) :: pres_scale(KA,IA,JA) ! Pressure
    real(RP) :: t_scale(KA,IA,JA)    ! Temperature
    real(RP) :: rhoc_scale(KA,IA,JA) ! cloud water density
    real(RP) :: rhor_scale(KA,IA,JA) ! rain water density
    real(RP) :: rhoi_scale(KA,IA,JA) ! cloud ice density
    real(RP) :: rhos_scale(KA,IA,JA) ! snow density
    real(RP) :: rhog_scale(KA,IA,JA) ! graupel density
    real(RP) :: QHYD_sdm(KA,IA,JA)

    !---------------------------------------------------------------------------

#ifdef _FIPP_
    ! Section specification for fipp profiler
    call fipp_start()
#endif
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_all",0,0)
#endif
    
    ! QTRC except QV (and QDRY) are diagnosed from super-droplets
    ! To make this doubly sure, reset QTRC to zero
    do iq = QQS, QQE
       if(iq /= I_QV) then
          QTRC(:,:,:,iq) = 0.0_RP
       end if
    end do

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_initialization",0,0)
#endif
    if( sd_first ) then
      sd_first = .false.
      if( IO_L ) write(IO_FID_LOG,*) '*** S.D.: setup'
      if( sd_rest_flg_in ) then
         !---- read restart file
         call ATMOS_PHY_MP_sdm_restart_in
         prr_crs(1:IA,1:JA,1:6)=0.0_RP
     else
         !---- set initial condition
         call sdm_iniset(DENS, RHOT, QTRC,                   &
                      RANDOM_IN_BASENAME, fid_random_i,   &
                      xmax_sdm, ymax_sdm, sdm_dtcmph,     &
                      sdm_rdnc,sdm_sdnmlvol,sdm_aslset,   &
                      sdm_inisdnc,sdm_zlower,             &
                      sdm_zupper,sdm_calvar,              &
                      zph_crs,                            &
                      sdasl_s2c, sdx_s2c, sdy_s2c,        &
                      sdz_s2c, sdr_s2c,                   &
                      sdrk_s2c, sdvz_s2c,                 &
                      sdrkl_s2c, sdrku_s2c                )
         prr_crs(1:IA,1:JA,1:6)=0.0_RP
      end if
    endif

    !! [issue] We need to change here. SDs are copied every time step uselessly as below.
    if (sd_rest_flg_out == .true.) then
       ! For restart array
       sdn_s2c_restart(:) = sdn_s2c(:)
       sdrk_s2c_restart(:) = sdrk_s2c(:)
       sdx_s2c_restart(:) = sdx_s2c(:)
       sdy_s2c_restart(:) = sdy_s2c(:)
       sdz_s2c_restart(:) = sdz_s2c(:)
       sdr_s2c_restart(:) = sdr_s2c(:)
       sdu_s2c_restart(:) = sdu_s2c(:)
       sdv_s2c_restart(:) = sdv_s2c(:)
       sdvz_s2c_restart(:) = sdvz_s2c(:)
       sdasl_s2c_restart(:,:) = sdasl_s2c(:,:)
       rng_s2c_restart = rng_s2c
       sdliqice_s2c_restart(:) = sdliqice_s2c(:)
       if( sdm_cold ) then
          sdice_s2c_restart%re(:) = sdice_s2c%re(:)
          sdice_s2c_restart%rp(:) = sdice_s2c%rp(:)
          sdice_s2c_restart%rho(:) = sdice_s2c%rho(:)
          sdice_s2c_restart%tf(:) = sdice_s2c%tf(:)
          sdice_s2c_restart%mrime(:) = sdice_s2c%mrime(:)
          sdice_s2c_restart%nmono(:) = sdice_s2c%nmono(:)
       end if
    end if

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_initialization",0,0)
#endif
    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(Super Droplet Method)'

    ! -----
    if( .not. sdm_calvar(1) .and. &
        .not. sdm_calvar(2) .and. &
        .not. sdm_calvar(3) .and. &
        .not. sdm_calvar(4) .and. &
        .not. sdm_calvar(5)       ) return

    !--
    ! SD data output

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_out",0,0)
#endif
    if( (mod(sdm_dmpvar,10)==1) .and. sdm_dmpitva>0.0_RP .and. &
         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
             int(1.E+3_RP*(sdm_dmpitva+0.00010_RP))) == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) ' *** Output Super Droplet Data in ASCII'
       !! Evaluate diagnostic variables
       !!! z
       call sdm_rk2z(sdnum_s2c,sdx_s2c,sdy_s2c,sdrk_s2c,sdz_s2c,sdri_s2c,sdrj_s2c)
       !!! terminal velocity vz
       call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
       call sdm_getvz_liq(pres_scale,DENS,t_scale,            &
                           sdnum_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdr_s2c,sdvz_s2c,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation' )
       if( sdm_cold )then
          call sdm_getvz_ice(DENS,t_scale,            &
                           sdnum_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdice_s2c,sdvz_s2c,  &
                           sd_itmp1,'no_interpolation' )
       end if
       !! Output
       call sdm_outasci(TIME_NOWSEC,                               &
                        sdnum_s2c,sdnumasl_s2c,                    &
                        sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdr_s2c,sdasl_s2c,sdvz_s2c, &
                        sdice_s2c, &
                        sdm_dmpnskip)
    end if

    if( ((mod(sdm_dmpvar,100))/10==1) .and. sdm_dmpitvb>0.0_RP .and. &
         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
             int(1.E+3_RP*(sdm_dmpitvb+0.00010_RP))) == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) ' *** Output Super Droplet Data in NetCDF'
       !! Evaluate diagnostic variables
       !!! z
       call sdm_rk2z(sdnum_s2c,sdx_s2c,sdy_s2c,sdrk_s2c,sdz_s2c,sdri_s2c,sdrj_s2c)
       !!! terminal velocity vz
       call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
       call sdm_getvz_liq(pres_scale,DENS,t_scale,            &
                           sdnum_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdr_s2c,sdvz_s2c,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation' )
       if( sdm_cold )then
          call sdm_getvz_ice(DENS,t_scale,            &
                           sdnum_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdice_s2c,sdvz_s2c,  &
                           sd_itmp1,'no_interpolation' )
       end if
       !! Output
       call sdm_outnetcdf(TIME_NOWSEC,                               &
                        sdnum_s2c,sdnumasl_s2c,                    &
                        sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdr_s2c,sdasl_s2c,sdvz_s2c, &
                        sdice_s2c, &
                        sdm_dmpnskip)
    end if

    if((sdm_dmpvar/100)>1)then
       write(*,*) 'WARNING: sdm_dmpvar=100 not supported for now. Set only the 1st (ascii) and 2nd (netcdf) digits'
    end if

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_out",0,0)
#endif
    !== run SDM at future ==!
     call sdm_calc(MOMX,MOMY,MOMZ,DENS,RHOT,QTRC,                 & 
                   sdm_calvar,sdm_mvexchg,sdm_dtcmph, sdm_aslset,  &
                   prr_crs,zph_crs,                      &
                   lsdmup,ni_s2c,nj_s2c,nk_s2c,                   &
                   sdnum_s2c,sdnumasl_s2c,                        &
                   sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,    &
                   sdu_s2c,sdv_s2c,sdvz_s2c,sdr_s2c,sdasl_s2c,sdice_s2c,    &
                   sdrkl_s2c,sdrku_s2c,                           &
                   rng_s2c,rand_s2c,sortid_s2c,sortkey_s2c,       &
                   sortfreq_s2c,sorttag_s2c,                      &
                   bufsiz1,                                       &
                   bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4,   &
                   sdm_itmp1,sdm_itmp2,           &
                   sd_itmp1,sd_itmp2,sd_itmp3,sd_dtmp1,sd_dtmp2,sd_dtmp3,sd_dtmp4,     &
                   crs_dtmp1,crs_dtmp2,crs_dtmp3,crs_dtmp4,       &
                   crs_dtmp5,crs_dtmp6,                           &
                   rbuf_r8,sbuf_r8,rbuf_i8,sbuf_i8,rbuf_i2,sbuf_i2,rbuf_i4,sbuf_i4) 

     !== convert updated contravariant velocity of ==!
     !== super-droplets to {u,v,w} at future       ==!

!     if( sdm_mvexchg>0 ) then

!        call cnt2phy(idsthopt,idtrnopt,idmpopt,idmfcopt,          &
!                     ni,nj,nk,j31,j32,jcb8w,mf,                   &
!                     uf_crs,vf_crs,wcf_crs,wf_crs,                &
!                     rtmp4,rtmp5,rtmp6)
!     end if

     !== Diagnose QC and QR from super-droplets ==!
     !! note that when SDM is used QC := rhoc/(rhod+rhov), QR := rhor/(rhod+rhov)

#ifdef _FAPP_
     ! Section specification for fapp profiler
     call fapp_start("sdm_sd2qcqr",0,0)
#endif
     call sdm_sd2qcqr(DENS,QTRC_sdm(:,:,:,I_QC),QTRC_sdm(:,:,:,I_QR),        &
                      zph_crs,                                               &
                      sdnum_s2c,sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,                     &
                      sdr_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,&
                      rhoc_scale,rhor_scale,                                 &
                      sd_itmp1,sd_itmp2,crs_dtmp1,crs_dtmp2)
#ifdef _FAPP_
     ! Section specification for fapp profiler
     call fapp_stop("sdm_sd2qcqr",0,0)
#endif

     !== Diagnose QI, QS and QR from super-droplets ==!
     !! note that when SDM is used QI := rhoi/(rhod+rhov), QS := rhos/(rhod+rhov), QG := rhog/(rhod+rhov)
     if( sdm_cold ) then
#ifdef _FAPP_
        ! Section specification for fapp profiler
        call fapp_start("sdm_sd2qiqsqg",0,0)
#endif
        call sdm_sd2qiqsqg(DENS,QTRC_sdm(:,:,:,I_QI),QTRC_sdm(:,:,:,I_QS),QTRC_sdm(:,:,:,I_QG),        &
                      zph_crs,                                               &
                      sdnum_s2c,sdn_s2c,sdliqice_s2c, sdx_s2c,sdy_s2c,                     &
                      sdice_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,&
                      rhoi_scale,rhos_scale,rhog_scale,                                 &
                      sd_itmp1,sd_itmp2,sd_itmp3,crs_dtmp1,crs_dtmp2,crs_dtmp3)
#ifdef _FAPP_
        ! Section specification for fapp profiler
        call fapp_stop("sdm_sd2qiqsqg",0,0)
#endif
     end if

! not supported yet. S.Shima
!!$     ! Aerosol formation process of super-droplets
!!$
!!$     if( sdm_aslfmdt>0.0_RP .and. &
!!$         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
!!$             int(1.E+3_RP*(sdm_aslfmdt+0.00010_RP))) == 0 ) then
!!$
!!$        call sdm_aslform(DENS,RHOT,QTRC,                             &   
!!$                         sdm_calvar,sdm_aslset,                      &
!!$                         sdm_aslfmsdnc,sdm_sdnmlvol,                 &
!!$                         sdm_zupper,sdm_zlower,dtcl,                 &
!!$                         jcb,pbr_crs,ptbr_crs,ppf_crs,               &
!!$                         ptpf_crs,qvf_crs,zph_crs,rhod_crs,          &
!!$                         sdnum_s2c,sdnumasl_s2c,sdn_s2c,sdliqice_s2c,sdx_s2c,     &
!!$                         sdy_s2c,sdz_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdu_s2c,           &
!!$                         sdv_s2c,sdvz_s2c,sdr_s2c,sdasl_s2c,         &
!!$                         sdfmnum_s2c,sdn_fm,sdx_fm,sdy_fm,sdz_fm,    &
!!$                         sdri_fm,sdrj_fm,sdrk_fm,sdvz_fm,sdr_fm,sdasl_fm,            &
!!$                         ni_s2c,nj_s2c,nk_s2c,                       &
!!$                         sortid_s2c,sortkey_s2c,sortfreq_s2c,        &
!!$                         sorttag_s2c,rng_s2c,                        &
!!$!                         sorttag_s2c,                                &
!!$                         sdm_itmp1,sd_itmp1,sd_itmp2,sd_itmp3)
!!$
!!$     end if

! not supported yet. S.Shima
!!$     ! Adjust number of super-droplets
!!$
!!$     if( sdm_nadjdt>0.0_RP .and. &
!!$         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
!!$             int(1.E+3_RP*(sdm_nadjdt+0.00010_RP))) == 0 ) then
!!$
!!$       !== averaged number concentration in a grid ==!
!!$
!!$       sd_nc = sdininum_s2c/real(ni_s2c*nj_s2c*knum_sdm,kind=RP)
!!$
!!$       call sdm_adjsdnum(sdm_nadjvar,ni_s2c,nj_s2c,nk_s2c,          &
!!$                         sdnum_s2c,sdnumasl_s2c,sd_nc,              &
!!$                         sdn_s2c,sdx_s2c,sdy_s2c,sdr_s2c,           &
!!$                         sdasl_s2c,sdvz_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,               &
!!$                         sortid_s2c,sortkey_s2c,sortfreq_s2c,       &
!!$                         sorttag_s2c,rng_s2c,rand_s2c,                &
!!$!                         sorttag_s2c,rand_s2c,                      &
!!$                         sdm_itmp1,sdm_itmp2,sd_itmp1)
!!$
!!$      end if

    do i = IS, IE
    do j = JS, JE
       SFLX_rain(i,j) = prr_crs(i,j,1) * rw ! surface rain rate [kg/m2/s]
       SFLX_snow(i,j) = prr_crs(i,j,3) * rw ! surface snow rate [kg/m2/s]
       !! check the unit needed for SFLX
    enddo
    enddo

    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
        QHYD_sdm(k,i,j) = QTRC_sdm(k,i,j,I_QC)+QTRC_sdm(k,i,j,I_QR)
    enddo
    enddo
    enddo

    if( sdm_cold ) then
       do k = KS, KE
       do i = IS, IE
       do j = JS, JE
          QHYD_sdm(k,i,j) = QHYD_sdm(k,i,j) + QTRC_sdm(k,i,j,I_QI)+QTRC_sdm(k,i,j,I_QS)+QTRC_sdm(k,i,j,I_QG)
       enddo
       enddo
       enddo
    end if

!!$    call HIST_in( rhoa_sdm(:,:,:), 'RAERO', 'aerosol mass conc.', 'kg/m3', dt)
    call HIST_in( prr_crs(:,:,2),       'RAIN_sd',    'surface rain accumulation',     'm')
    call HIST_in( QTRC_sdm(:,:,:,I_QC), 'QC_sd',   'mixing ratio of cloud water in SDM',  'kg/kg')
    call HIST_in( QTRC_sdm(:,:,:,I_QR), 'QR_sd',   'mixing ratio of rain water in SDM',   'kg/kg')
    if( sdm_cold ) then
       call HIST_in( prr_crs(:,:,4),       'SNOW_sd',    'surface snow accumulation',    'm')
       call HIST_in( QTRC_sdm(:,:,:,I_QI), 'QI_sd',   'mixing ratio of cloud ice in SDM',  'kg/kg')
       call HIST_in( QTRC_sdm(:,:,:,I_QS), 'QS_sd',   'mixing ratio of snow in SDM',   'kg/kg')
       call HIST_in( QTRC_sdm(:,:,:,I_QG), 'QG_sd',   'mixing ratio of graupel in SDM',   'kg/kg')
    end if
    call HIST_in( prr_crs(:,:,6),       'PREC_sd',    'surface precipitation accumulation',    'm')
    call HIST_in( QHYD_sdm(:,:,:),      'QHYD_sd', 'mixing ratio of liquid water and ice in SDM', 'kg/kg')

#ifdef _FIPP_
    ! Section specification for fipp profiler
    call fipp_stop()
#endif
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_all",0,0)
#endif
    return
  end subroutine ATMOS_PHY_MP_sdm
  !-----------------------------------------------------------------------------
   subroutine sdm_iniset(DENS, RHOT, QTRC,                   &
                         RANDOM_IN_BASENAME, fid_random_i,   &
                         xmax_sdm, ymax_sdm, dtcmph,         &
                         sdm_rdnc,sdm_sdnmlvol,sdm_aslset,   &
                         sdm_inisdnc,sdm_zlower,             &
                         sdm_zupper,sdm_calvar,              &
!                         dx_sdm,dy_sdm,dz_sdm,               &
!                         nqw,jcb,                            &
!                         qwtr_crs,zph_crs,                   &
!                         jcb,                                &
                         zph_crs,                            &
                         sdasl_s2c, sdx_s2c, sdy_s2c,        &
                         sdz_s2c, sdr_s2c,                   &
                         sdrk_s2c, sdvz_s2c,                 &
                         sdrkl_s2c, sdrku_s2c                )
  !***********************************************************************
  ! Input variables
      use scale_const, only: &
        P00 => CONST_PRE00, &
        Rdry => CONST_Rdry, &
        Rvap => CONST_Rvap, &
        CPdry => CONST_CPdry
      use scale_process, only: &
        PRC_MPIstop, &
        mype => PRC_myrank
      use scale_atmos_thermodyn, only: &
        CPw => AQ_CP
      use scale_tracer, only: &
        QAD => QA
      use scale_grid, only: &
        GRID_FX,    &
        GRID_FY
      use m_sdm_coordtrans, only: &
        sdm_z2rk
      use m_sdm_fluidconv, only: &
        sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw
      use m_sdm_sd2fluid, only: &
           sdm_sd2qcqr,sdm_sd2qiqsqg
      use m_sdm_condensation_water, only: &
           sdm_condevp
      use m_sdm_meltfreeze, only: &
           sdm_meltfreeze

      real(RP), intent(in) :: DENS(KA,IA,JA) ! Density     [kg/m3]
      real(RP), intent(in) :: RHOT(KA,IA,JA) ! DENS * POTT [K*kg/m3]
      real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD) ! ratio of mass of tracer to total mass[kg/kg]
      character(len=H_LONG), intent(in) :: RANDOM_IN_BASENAME
      integer, intent(in) :: fid_random_i
      real(RP),intent(in) :: xmax_sdm, ymax_sdm
      real(RP),intent(in) :: dtcmph(5)  ! Time interval of cloud micro physics
      logical, intent(in) :: sdm_calvar(5)! Flag for cond./coll/move calculation
      real(RP),intent(in) :: sdm_rdnc     ! Number concentration of real droplets
      real(RP),intent(in) :: sdm_sdnmlvol ! Normal volume for number concentration of super droplets
      integer, intent(in) :: sdm_aslset   ! Option for aerosol species
      real(RP),intent(in) :: sdm_inisdnc  ! Initial number of super droplets per sdm_sdnmlvol
      real(RP),intent(in) :: sdm_zlower   ! Lower limitaion of initial SDs position
      real(RP),intent(in) :: sdm_zupper   ! Upper limitaion of initial SDs position
      real(RP),intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinates
      real(RP),intent(inout) :: sdasl_s2c(1:sdnum_s2c,1:sdnumasl_s2c)
      real(RP),intent(inout) :: sdx_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdy_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdz_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdr_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdrk_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdvz_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdrkl_s2c(IA,JA)
      real(RP),intent(inout) :: sdrku_s2c(IA,JA)
      ! Work variables
      real(RP) :: n0                            ! number of real droplets per unit volume and per aerosol radius
      real(RP) :: dry_r                         ! aerosol radius
      real(RP) :: delta1, delta2, sdn_tmp       ! temporary
      logical :: lsdmup                         ! flag for updating water hydrometeor by SDM
      integer :: iexced, sdnum_tmp1, sdnum_tmp2 ! temporary
      integer :: i, j, k, n, iq, np             ! index
      real(RP) :: crs_dtmp1(KA,IA,JA), crs_dtmp2(KA,IA,JA), crs_dtmp3(KA,IA,JA)
      integer :: sd_str, sd_end, sd_valid

      real(RP) :: pres_scale(KA,IA,JA)  ! Pressure
      real(RP) :: t_scale(KA,IA,JA)    ! Temperature
      real(RP) :: rhoc_scale(KA,IA,JA) ! cloud water density
      real(RP) :: rhor_scale(KA,IA,JA) ! rain water density
      real(RP) :: rhoi_scale(KA,IA,JA) ! cloud ice density
      real(RP) :: rhos_scale(KA,IA,JA) ! snow density
      real(RP) :: rhog_scale(KA,IA,JA) ! graupel density

      real(RP) :: sdm_dtevl  ! time step of {condensation/evaporation} process
      real(RP) :: sdm_dtcol  ! time step of {stochastic coalescence} process
      real(RP) :: sdm_dtadv  ! time step of {motion of super-droplets} process
      real(RP) :: sdm_dtmlt  ! time step of {melt/freeze of super-droplets} process
      real(RP) :: sdm_dtsbl  ! time step of {sublimation/deposition of super-droplets} process

     !
      real(RP) :: area, INAS_max, prob_INIA, INAS_tf, probdens_tf

      !---------------------------------------------------------------------

      ! Initialize and rename variables
      sdm_dtevl = real( sdm_dtcmph(1),kind=RP )  !! condensation/evaporation
      sdm_dtcol = real( sdm_dtcmph(2),kind=RP )  !! stochastic coalescence
      sdm_dtadv = real( sdm_dtcmph(3),kind=RP )  !! motion of super-droplets
      sdm_dtmlt = real( sdm_dtcmph(4),kind=RP )  !! melting/freezing
      sdm_dtsbl = real( sdm_dtcmph(5),kind=RP )  !! sublimation/depsition

      if( .not. sdm_calvar(1) .and. &
          .not. sdm_calvar(2) .and. &
          .not. sdm_calvar(3) .and. &
          .not. sdm_calvar(4) .and. &
          .not. sdm_calvar(5)        ) return

      !### Get random generator seed ###!
      !! Random number generator has already been initialized in scale-les/src/preprocess/mod_mkinit.f90
      !! Be careful. If unit (=fid_random_i) is specified, filename is ignored and the object is initialized by the unit.
      call rng_load_state( rng_s2c, trim(RANDOM_IN_BASENAME))
!      call rng_load_state( rng_s2c, trim(RANDOM_IN_BASENAME), fid_random_i )

      ! Initialized super-droplets.
      !### Get parameter for 3mode log-nomiral distribution ###!

      if( mod(sdm_aslset,10)==1 .or. mod(sdm_aslset,10)==3 ) then

         !! Derksen(2009)
         if( IO_L ) then
          write(IO_FID_LOG,*)  "    number concentration of (NH4)2SO4    : Derksen"
         endif
         n1_amsul   = n1_amsul_derksn
         n2_amsul   = n2_amsul_derksn
         rb1_amsul  = rb1_amsul_derksn
         rb2_amsul  = rb2_amsul_derksn
         sgm1_amsul = sgm1_amsul_derksn
         sgm2_amsul = sgm2_amsul_derksn
         rmax_amsul = rmax_amsul_derksn
         rmin_amsul = rmin_amsul_derksn

         if( mod(sdm_aslset,10)==1 ) then
            n3_nacl    = 1.E-10_RP     !! temporary initialize
            rb3_nacl   = 1.E-10_RP
            sgm3_nacl  = 1.E-10_RP
            rmax_nacl  = 1.E-10_RP
            rmin_nacl  = 1.E-10_RP
         end if

      else if( mod(sdm_aslset,10)==-1 .or. mod(sdm_aslset,10)==-3 ) then

         !! vanZanten(2010)
         if( IO_L ) then
          write(IO_FID_LOG,*) "    number concentration of (NH4)2SO4    : vanZanten"
         endif

         n1_amsul   = n1_amsul_zanten
         n2_amsul   = n2_amsul_zanten
         rb1_amsul  = rb1_amsul_zanten
         rb2_amsul  = rb2_amsul_zanten
         sgm1_amsul = sgm1_amsul_zanten
         sgm2_amsul = sgm2_amsul_zanten
         rmax_amsul = rmax_amsul_zanten
         rmin_amsul = rmin_amsul_zanten

         if( mod(sdm_aslset,10)==-1 ) then
            n3_nacl    = 1.E-10_RP     !! temporary initialize
            rb3_nacl   = 1.E-10_RP
            sgm3_nacl  = 1.E-10_RP
            rmax_nacl  = 1.E-10_RP
            rmin_nacl  = 1.E-10_RP
         end if

      end if

      if( mod(sdm_aslset,10)==2 .or. abs(mod(sdm_aslset,10))==3 ) then

         !! Derksen(2009)
         if( IO_L ) then
          write(IO_FID_LOG,*) "number concentration of NaCl aerosol : Derksen"
         endif
         if( mod(sdm_aslset,10)==2 ) then
            n1_amsul   = 1.E-10_RP     !! temporary initialize
            n2_amsul   = 1.E-10_RP
            rb1_amsul  = 1.E-10_RP
            rb2_amsul  = 1.E-10_RP
            sgm1_amsul = 1.E-10_RP
            sgm2_amsul = 1.E-10_RP
            rmax_amsul = 1.E-10_RP
            rmin_amsul = 1.E-10_RP
         end if
         n3_nacl   = n3_nacl_derksn
         rb3_nacl  = rb3_nacl_derksn
         sgm3_nacl = sgm3_nacl_derksn
         rmax_nacl = rmax_nacl_derksn
         rmin_nacl = rmin_nacl_derksn

      end if

      !### Get random number ###!
      call gen_rand_array( rng_s2c, sdx_s2c  )
      call gen_rand_array( rng_s2c, sdy_s2c  )
      call gen_rand_array( rng_s2c, sdz_s2c  )
      call gen_rand_array( rng_s2c, sdr_s2c  )
      call gen_rand_array( rng_s2c, sdvz_s2c )

      iexced = 1     !! check for memory size of int*8

      do k=1,sdnumasl_s2c
         call gen_rand_array( rng_s2c, sd_dtmp1 )
         do n=1,sdnum_s2c
            sdasl_s2c(n,k) = sd_dtmp1(n)
         end do
      end do

      if( sdm_cold ) then
         call gen_rand_array( rng_s2c, sdice_s2c%tf )
      end if

      ! Initialized all super-droplets as water droplet.
      !### status(liquid/ice) of super-droplets ###!
      sdliqice_s2c(1:sdnum_s2c) = STAT_LIQ

      !### Aerosol mass, muliplicity ###!
      do k=1,sdnumasl_s2c
         do n=1,sdnum_s2c
            i = mod(n-1,sdnumasl_s2c) + 1    !! select aerosol index
            if( abs(sdm_aslset)==12 ) then
               i = 2                         !! 1 : (NH4)2SO4, 2: NaCl
            end if
            if( k==i ) then
               !! match aerosol type
               if( abs(mod(sdm_aslset,10))==1 .or.                      &
                     ( k==1 .and. abs(mod(sdm_aslset,10))==3 ) ) then

                  !### (NH4)2SO4 [g] ###!
                  delta1 = log(rmax_amsul) - log(rmin_amsul)
                  dry_r  = exp( log(rmin_amsul)+delta1*sdasl_s2c(n,k) )
                  sdasl_s2c(n,k) = F_THRD * ONE_PI                 &
                                * (dry_r*dry_r*dry_r) * rho_amsul
                  !! n0(log(dry_r)) [m-3]
                  !! 2-mode log-noraml distribution for ammonium sulfate
                  delta1 = log(dry_r) - log(rb1_amsul)
                  delta1 = -(delta1*delta1)                             &
                               /(2.0_RP*log(sgm1_amsul)*log(sgm1_amsul))
                  delta2 = log(dry_r) - log(rb2_amsul)
                  delta2 = -(delta2*delta2)                             &
                               /(2.0_RP*log(sgm2_amsul)*log(sgm2_amsul))
                  n0 = (n1_amsul*exp(delta1))                           &
                                 /(sqrt(2.0_RP*ONE_PI)*log(sgm1_amsul))  &
                     + (n2_amsul*exp(delta2))                           &
                                 /(sqrt(2.0_RP*ONE_PI)*log(sgm2_amsul))
                  !! number per unit volume and per aerosol species
                  delta1 = real(sdm_inisdnc,kind=RP)                    &
                                    /real(sdm_sdnmlvol,kind=RP)
                  delta1 = delta1/real(sdnumasl_s2c,kind=RP)
                  !! continuous uniform distribution
                  delta2 = 1.0_RP/(log(rmax_amsul)-log(rmin_amsul))
                  !! muliplicity
                  sdn_tmp = n0/(delta1*delta2)
               else if( abs(mod(sdm_aslset,10))==2 .or.                 &
                     ( k==2 .and. abs(mod(sdm_aslset,10))==3 ) ) then
                  !### NaCl(seasalt,[g]) ###!
                  delta1 = log(rmax_nacl) - log(rmin_nacl)
                  dry_r  = exp( log(rmin_nacl) + delta1*sdasl_s2c(n,k) )
                  sdasl_s2c(n,k) = F_THRD * ONE_PI                 &
                                * (dry_r*dry_r*dry_r) * rho_nacl
                  !! n0(log(dry_r)) [m-3]
                  !! 1-mode log-noraml distribution for seasalt
                  delta1 = log(dry_r) - log(rb3_nacl)
                  delta1 = -(delta1*delta1)                             &
                               /(2.0_RP*log(sgm3_nacl)*log(sgm3_nacl))
                  n0 = (n3_nacl*exp(delta1))                            &
                                 /(sqrt(2.0_RP*ONE_PI)*log(sgm3_nacl))
                  !! number per unit volume and per aerosol species
                  delta1 = real(sdm_inisdnc,kind=RP)                    &
                                    /real(sdm_sdnmlvol,kind=RP)
                  delta1 = delta1/real(sdnumasl_s2c,kind=RP)
                  !! continuous uniform distribution
                  delta2 = 1.0_RP/(log(rmax_nacl)-log(rmin_nacl))
                  !! muliplicity
                  sdn_tmp = n0/(delta1*delta2)
               else
                  !### Other ###!
                  sdasl_s2c(n,k) = 1.0E-18_RP                           &
                                 + 1.0E-14_RP                           &
                                 * (log(1.0_RP/(1.0_RP-sdasl_s2c(n,k))))
                  sdn_tmp = real(sdm_rdnc,kind=RP)                      &
                          * real(sdm_sdnmlvol,kind=RP)                  &
                          / real(sdm_inisdnc,kind=RP)
               end if

	       !! Multiply the initial number density of soluble aerosol particles by sdm_fctr2multi
	       !! This only for soluble particles, not for insoluble particles
	       sdn_tmp = sdn_tmp * sdm_fctr2multi

               !! check muliplicity
               if( sdn_tmp<(2.0_RP**63.0_RP) ) then
                  sdn_s2c(n) = nint( sdn_tmp, kind=DP )
               else
                  iexced = -1
               end if
            else
               sdasl_s2c(n,k) = 0.0_RP
            end if
         end do
      end do

      ! Initialized SDs for ice phase (freezing temperature and multiplicity)
      if( sdm_cold ) then

         !###### soluble aerosol ######!
         do n=1,nint(sdininum_s2c*sdnumratio_soluble)
            sdice_s2c%tf(n) = -38.d0  ! homogeneous freezing limit [degC] 
            sdn_s2c(n) =  nint( real(sdn_s2c(n))/sdnumratio_soluble, kind=DP ) ! adjust multiplicity
         end do

         !###### insluble+soluble aerosol (internally mixed) ######!
         if(nint(sdininum_s2c*sdnumratio_soluble)+1 < nint(sdininum_s2c)) then

            !###### set parameters ######!
            ! set aerosol distribution parameters of soluble component
            n1_amsul   = n_mdust * (n1_amsul_zanten/(n1_amsul_zanten+n2_amsul_zanten))
            n2_amsul   = n_mdust * (n2_amsul_zanten/(n1_amsul_zanten+n2_amsul_zanten))
            rb1_amsul  = rb1_amsul_zanten
            rb2_amsul  = rb2_amsul_zanten
            sgm1_amsul = sgm1_amsul_zanten
            sgm2_amsul = sgm2_amsul_zanten
            rmax_amsul = rmax_amsul_zanten
            rmin_amsul = rmin_amsul_zanten

            !###### mass of soluble component ######!
            ! reset all chemical components
            do n=nint(sdininum_s2c*sdnumratio_soluble)+1,nint(sdininum_s2c)
               do k=1,sdnumasl_s2c
                  sdasl_s2c(n,k) = 0.0_RP
               end do
            end do
            !### (NH4)2SO4 [g] ###!
            k = 1
            ! prepare random numbers
            call gen_rand_array( rng_s2c, sd_dtmp1 )
            do n=nint(sdininum_s2c*sdnumratio_soluble)+1,nint(sdininum_s2c)
               delta1 = log(rmax_amsul) - log(rmin_amsul)
               dry_r  = exp( log(rmin_amsul)+delta1*sd_dtmp1(n) )
               sdasl_s2c(n,k) = F_THRD * ONE_PI                 &
                    * (dry_r*dry_r*dry_r) * rho_amsul
            end do

            !###### mass of insoluble component ######!
            !### Assume a monodisperse distribution of mineral dust with a diameter of mdust_dia ###!
            ! mdust_dia [m] is defined in sdm_common.f90

            !###### freezing temperature ######!
            ! prepare random numbers
            call gen_rand_array( rng_s2c, sd_dtmp1 )
            do n=nint(sdininum_s2c*sdnumratio_soluble)+1,nint(sdininum_s2c)
               delta1 = tfmax-tfmin
               sdice_s2c%tf(n) = tfmin+delta1*sd_dtmp1(n)
            end do

            !###### multiplicity ######!
            ! prepare random numbers
            call gen_rand_array( rng_s2c, sd_dtmp1 )
            call gen_rand_array( rng_s2c, sd_dtmp2 )
            k=1 ! internally mixed with (NH4)2SO4
            do n=nint(sdininum_s2c*sdnumratio_soluble)+1,nint(sdininum_s2c)
               !### contribution from soluble part ###!
               !! n0(log(dry_r)) [m-3]
               !! 2-mode log-noraml distribution for ammonium sulfate
               dry_r = (sdasl_s2c(n,k)/F_THRD/ONE_PI/rho_amsul)**O_THRD
               delta1 = log(dry_r) - log(rb1_amsul)
               delta1 = -(delta1*delta1)                             &
                    /(2.0_RP*log(sgm1_amsul)*log(sgm1_amsul))
               delta2 = log(dry_r) - log(rb2_amsul)
               delta2 = -(delta2*delta2)                             &
                    /(2.0_RP*log(sgm2_amsul)*log(sgm2_amsul))
               n0 = (n1_amsul*exp(delta1))                           &
                    /(sqrt(2.0_RP*ONE_PI)*log(sgm1_amsul))  &
                    + (n2_amsul*exp(delta2))                           &
                    /(sqrt(2.0_RP*ONE_PI)*log(sgm2_amsul))
               !! number of SDs per unit volume
               delta1 = (1.0_RP-sdnumratio_soluble)*real(sdm_inisdnc,kind=RP)   &
                    /real(sdm_sdnmlvol,kind=RP)
               !! continuous uniform distribution
               delta2 = 1.0_RP/(log(rmax_amsul)-log(rmin_amsul))
               !! muliplicity
               sdn_tmp = n0/(delta1*delta2)

               !### contribution from insluble part ###!
               area = ONE_PI*mdust_dia**2
               if(sd_dtmp1(n) < INIAsd_ratio) then
                  !! IN inactive mineral dust
                  INAS_max = a0_N12*exp(-a1_N12*tfmin+a2_N12) ! INAS of Niemand et al. (2012) 
                  prob_INIA = exp(-area*INAS_max)
                  sdn_tmp = sdn_tmp * prob_INIA / INIAsd_ratio
                  sdice_s2c%tf(n) = -38.0_RP  ! homogeneous freezing limit
               else
                  !! IN active mineral dust
                  INAS_tf = a0_N12*exp(-a1_N12*sdice_s2c%tf(n)+a2_N12)  ! INAS of Niemand et al. (2012) 
                  probdens_tf = area*a1_N12*INAS_tf*exp(-area*INAS_tf)
                  delta1 = tfmax-tfmin
                  sdn_tmp = sdn_tmp * probdens_tf * delta1 / (1.0_RP - INIAsd_ratio)
               end if

               !! check multiplicity
               if( sdn_tmp<(2.0_RP**63.0_RP) ) then
                  if(sdn_tmp>1.0_RP)then
                     sdn_s2c(n) = nint( sdn_tmp, kind=DP )
                  else ! to sample rare event
                     sdn_s2c(n) = 0 
                     if( sd_dtmp2(n) < sdn_tmp )then
                        sdn_s2c(n) = 1
                     end if
                  end if
               else
                  iexced = -1
               end if

            end do

         end if

      end if

      if( iexced<0 ) then
         write(*,*) "sdm_iniset, exceeded"
         call PRC_MPIstop
      end if

      !### position of super-droplets in horizontal ###!
      do n=1,sdnum_s2c
         sdx_s2c(n) = xmax_sdm * sdx_s2c(n)+GRID_FX(IS-1)
         sdy_s2c(n) = ymax_sdm * sdy_s2c(n)+GRID_FY(JS-1)
      end do

      !### position of super-droplets in vertical ###!
      !! valid super-droplets
      do n=1,nint(sdininum_s2c)
         if( sdn_s2c(n)>0 ) then
            sdz_s2c(n) = real(minzph+sdm_zlower,kind=RP)             &
                 + sdz_s2c(n)                                  &
                 * real(sdm_zupper-(minzph+sdm_zlower),kind=RP)
         else
            sdz_s2c(n) = INVALID     !!! check muliplicity
         end if
      end do

      !! invalid super-droplets
      do n=nint(sdininum_s2c)+1,sdnum_s2c
         sdz_s2c(n) = INVALID
      end do
      
!!$      sdnum_tmp1 = int( nint(sdininum_s2c)/nomp )
!!$      sdnum_tmp2 = mod( nint(sdininum_s2c),nomp )
!!$
!!$      do np=1,nomp
!!$
!!$         sd_str = int(sdnum_s2c/nomp)*(np-1) + 1
!!$         sd_end = int(sdnum_s2c/nomp)*np
!!$
!!$         if( np<=sdnum_tmp2 ) then
!!$            sd_valid = sd_str + sdnum_tmp1
!!$         else
!!$            sd_valid = sd_str + (sdnum_tmp1-1)
!!$         end if
!!$
!!$         !! valid super-droplets
!!$
!!$         do n=sd_str,sd_valid
!!$            if( sdn_s2c(n)>0 ) then
!!$               sdz_s2c(n) = real(minzph+sdm_zlower,kind=RP)             &
!!$                          + sdz_s2c(n)                                  &
!!$                          * real(sdm_zupper-(minzph+sdm_zlower),kind=RP)
!!$            else
!!$               sdz_s2c(n) = INVALID     !!! check muliplicity
!!$            end if
!!$         end do
!!$
!!$         !! invalid super-droplets
!!$         do n=sd_valid+1,sd_end
!!$            sdz_s2c(n) = INVALID
!!$         end do
!!$      end do

      !### index[k/real] of super-droplets               ###!
      !### modify position[z] of invalid super-droplets  ###!
      call sdm_z2rk(sdm_zlower,sdm_zupper,            &
                        sdnum_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c )

      !### initial equivalent radius                     ###!
      do n=1,sdnum_s2c
         sdr_s2c(n) = 1.0E-15_RP

!ORG     sdr_s2c(n) = 1.0e-5 * ( log(1.0/(1.0-sdr_s2c(n))) )**O_THRD
!ORG     sdr_s2c(n) = 1.0d-8

! temporary for test
!         sdr_s2c(n) = 3.0E-3_RP*sdr_s2c(n)
!         sdr_s2c(n) = exp((log(3.0E-3_RP)-log(1.0E-7_RP))*sdr_s2c(n)+log(1.0E-7_RP))
         
      end do

      !### ( at only condensation/evaporation process ) ###!
      if( sdm_calvar(1) ) then

         call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)

         call sdm_condevp(sdm_aslset,                                   &
                          sdm_aslmw,sdm_aslion,sdm_dtevl,               &
                          pres_scale,t_scale,QTRC(:,:,:,I_QV),          &
                          sdnum_s2c,sdnumasl_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,       &
                          sdr_s2c,sdasl_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c)

      end if

      !### ( melting/freezing process ) ###!
      if( sdm_cold .and. sdm_calvar(4) ) then

         call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)

         call sdm_meltfreeze(                              &
              t_scale,pres_scale,QTRC(:,:,:,I_QV),         &
              sdnum_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,      &
              sdr_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdice_s2c )

      end if

      !### diagnose terminal velocity (no need evaluate them here) ###!
!!$      call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
!!$      call sdm_getvz_liq(pres_scale,DENS,t_scale,                    &
!!$                     sdnum_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdr_s2c, &
!!$                     sdvz_s2c,sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation')

      !### Diagnose QC and QR from super-droplets ###!
      !! note that when SDM is used QC := rhoc/(rhod+rhov), QR := rhor/(rhod+rhov)
      call sdm_sd2qcqr(DENS,QTRC_sdm(:,:,:,I_QC),QTRC_sdm(:,:,:,I_QR),          &
                       zph_crs,                   &
                       sdnum_s2c,sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,        &
                       sdr_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,            &
                       rhoc_scale,rhor_scale,                      &
                       sd_itmp1,sd_itmp2,crs_dtmp1,crs_dtmp2            )

      !### Diagnose QI, QS, QG from super-droplets ###!
      if( sdm_cold ) then
         call sdm_sd2qiqsqg(DENS,QTRC_sdm(:,:,:,I_QI),QTRC_sdm(:,:,:,I_QS),QTRC_sdm(:,:,:,I_QG),        &
                      zph_crs,                                               &
                      sdnum_s2c,sdn_s2c,sdliqice_s2c, sdx_s2c,sdy_s2c,                     &
                      sdice_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,&
                      rhoi_scale,rhos_scale,rhog_scale,                                 &
                      sd_itmp1,sd_itmp2,sd_itmp3,crs_dtmp1,crs_dtmp2,crs_dtmp3)
      end if

      ! Output logfile about SDM
      if( mype==0 ) then

         if( IO_L ) then
          write(IO_FID_LOG,*)
          write(IO_FID_LOG,'(a)')"  ### [SDM] : Information of super-droplets ###"
          write(IO_FID_LOG,*)
         endif
!          write(moji,'(i25)')sdnum_s2c
!          write(*,'(a43,a25)')                                           &
!      &     "    allocate size for super-droplets     : ", adjustl(moji)

!          write(moji,'(i25)')nint(sdininum_s2c)
!          write(IO_FID_LOG,'(a43,a25)')                                 &
!      &     "    initial number of super-droplets     : ", adjustl(moji)

      end if

    return

  end subroutine sdm_iniset
  !-----------------------------------------------------------------------------
  subroutine sdm_calc(MOMX,MOMY,MOMZ,DENS,RHOT,QTRC,              & 
                      sdm_calvar,sdm_mvexchg,dtcl, sdm_aslset,    &
                      prec_crs,zph_crs,   &
                      lsdmup,ni_sdm,nj_sdm,nk_sdm,                &
                      sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,      &
                      sd_u,sd_v,sd_vz,sd_r,sd_asl,sdi,sd_rkl,sd_rku,  &
                      sd_rng,sd_rand,sort_id,sort_key,sort_freq,  &
                      sort_tag,                                   &
                      bufsiz1,                                    & 
                      bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4,      &
                      sdm_itmp1,sdm_itmp2,        &
                      sd_itmp1,sd_itmp2,sd_itmp3,sd_dtmp1,sd_dtmp2,sd_dtmp3,sd_dtmp4,   &
                      crs_val1p,crs_val1c,crs_val2p,crs_val2c,    &
                      crs_val3p,crs_val3c,                        &
                      rbuf_r8,sbuf_r8,rbuf_i8,sbuf_i8,rbuf_i2,sbuf_i2,rbuf_i4,sbuf_i4) 
   use scale_process, only: &
       PRC_MPIstop
   use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
   use scale_tracer, only: &
       QAD => QA
   use scale_comm, only: &
        COMM_vars8, &
        COMM_wait
   use scale_const, only: &
       GRAV   => CONST_GRAV
   use m_sdm_fluidconv, only: &
        sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw
   use m_sdm_sd2fluid, only: &
        sdm_sd2rhow, sdm_sd2prec, sdm_sd2rhosol
   use m_sdm_boundary, only: &
        sdm_jdginvdv, sdm_boundary
    use m_sdm_motion, only: &
        sdm_getvz_liq, sdm_getvz_ice, sdm_getvel, sdm_move
    use m_sdm_coalescence, only: &
        sdm_coales
    use m_sdm_coalescence_cold, only: &
        sdm_coales_cold
    use m_sdm_condensation_water, only: &
        sdm_condevp, sdm_condevp_updatefluid
    use m_sdm_meltfreeze, only: &
        sdm_meltfreeze, sdm_meltfreeze_updatefluid
    use m_sdm_subldep, only: &
        sdm_subldep, sdm_subldep_updatefluid

   real(RP), intent(inout) :: DENS(KA,IA,JA)        !! Density [kg/m3]
   real(RP), intent(inout) :: MOMZ(KA,IA,JA)        !! Momentum [kg/s/m2]
   real(RP), intent(inout) :: MOMX(KA,IA,JA)
   real(RP), intent(inout) :: MOMY(KA,IA,JA)
   real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
   real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)    !! Ratio of mass of tracer to total mass[kg/kg]
   ! Input variables
   logical,  intent(in) :: sdm_calvar(5)
   integer,  intent(in) :: sdm_aslset
   integer,  intent(in) :: sdm_mvexchg
   real(RP), intent(in) :: dtcl(1:5)     ! Time interval [sec] of condensation, coalescence, and droplet motion
                                         ! melting/freezing [coldsdm only], sublimation/deposition [coldsdm only]
   real(RP), intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinates
   integer, intent(in) :: ni_sdm                    ! SDM model dimension in x direction
   integer, intent(in) :: nj_sdm                    ! SDM model dimension in y direction
   integer, intent(in) :: nk_sdm                    ! SDM model dimension in z direction
   integer, intent(in) :: sd_num                    ! number of super-droplets
   integer, intent(in) :: sd_numasl     ! number of kind of chemical material contained as water-soluble aerosol in super droplets
   real(RP), intent(in) :: sd_rkl(IA,JA) ! index[k/real] at lower boundary in SDM
   real(RP), intent(in) :: sd_rku(IA,JA) ! index[k/real] at upper boundary in SDM
   integer, intent(in) :: bufsiz1       ! buffer size for MPI
   integer, intent(in) :: bufsiz2_r8 ! buffer size for MPI (real8)
   integer, intent(in) :: bufsiz2_i8 ! buffer size for MPI (int8)
   integer, intent(in) :: bufsiz2_i2 ! buffer size for MPI (int2)
   integer, intent(in) :: bufsiz2_i4 ! buffer size for MPI (int4)
   ! Input and output variables
   integer(DP), intent(inout) :: sd_n(1:sd_num)    ! multiplicity of super-droplets
   integer(i2), intent(inout) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
   real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
   real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
   real(RP), intent(inout) :: sd_ri(1:sd_num) ! face index-i(real) of super-droplets
   real(RP), intent(inout) :: sd_rj(1:sd_num) ! face index-j(real) of super-droplets
   real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets in vertical
   real(RP), intent(inout) :: sd_u(1:sd_num)  ! x-components velocity of super-droplets
   real(RP), intent(inout) :: sd_v(1:sd_num)  ! y-components velocity of super-droplets
   real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets /z velocity of super-droplets
   real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
   real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
   type(sdicedef), intent(inout) :: sdi   ! ice phase super-droplets
   type(c_rng_uniform_mt), intent(inout) :: sd_rng ! random number generator
   real(RP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
   integer, intent(inout) :: sort_id(1:sd_num)  ! id that super-droplets sorted by grids
   integer, intent(inout) :: sort_key(1:sd_num) ! sort key
   integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1) ! number of super-droplets in each grid
   integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)  ! accumulated number of super-droplets in each grid
   real(RP), intent(inout) :: prec_crs(IA,JA,1:6)! precipitation rate and accumlation [m]
                                               ! 1:rain rate, 2:rain accumulation
                                               ! 3:snow rate, 4:snow accumulation
                                               ! 5:precipitation rate, 6:precipitation accumulation!
   real(DP), intent(inout) ::                                   &
         &                             rbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Receiving buffer (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+4) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    real(DP), intent(inout) ::                                   &
         &                             sbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+4) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer(DP), intent(inout) ::                                &
         &                             rbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! reciving buffer for MPI (int8)
                       ! dim02 = 1 (multiplicity of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer(DP), intent(inout) ::                                &
         &                             sbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! sending buffer for MPI (int8)
                       ! dim02 = 1 (multiplicity of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer(i2), intent(inout) ::                                &
         &                             rbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! reciving buffer for MPI (int2)
                       ! dim02 = 1 (status(liq/ice) of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer(i2), intent(inout) ::                                &
         &                             sbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! sending buffer for MPI (int2)
                       ! dim02 = 1 (status(liq/ice) of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer, intent(inout) ::                                &
         &                             rbuf_i4(1:,1:,1:)
                       ! rbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or rbuf_i4(1:1,1:1,1:2) 
                       ! reciving buffer for MPI (int4)
                       ! dim02 = 1 (sdi_nmono of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    integer, intent(inout) ::                                &
         &                             sbuf_i4(1:,1:,1:)
                       ! sbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or sbuf_i4(1:1,1:1,1:2) 
                       ! sending buffer for MPI (int4)
                       ! dim02 = 1 (sdi_nmono of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
   ! Output variables
   logical, intent(out) :: lsdmup  ! flag for updating water hydrometeor by SDM
   integer, intent(out) :: sdm_itmp1(1:ni_sdm*nj_sdm*nk_sdm+2) ! temporary array of SDM dimension
   integer, intent(out) :: sdm_itmp2(1:ni_sdm*nj_sdm*nk_sdm+2) ! temporary array of SDM dimension
   integer, intent(out) :: sd_itmp1(1:sd_num) ! temporary array of the size of the number of super-droplets.
   integer, intent(out) :: sd_itmp2(1:sd_num) ! temporary array of the size of the number of super-droplets.
   integer, intent(out) :: sd_itmp3(1:sd_num) ! temporary array of the size of the number o fsuper-droplets.
   real(RP), intent(out) :: sd_dtmp1(1:sd_num) ! temporary array of the size of the number of super-droplets.
   real(RP), intent(out) :: sd_dtmp2(1:sd_num) ! temporary array of the size of the number of super-droplets.
   real(RP), intent(out) :: sd_dtmp3(1:sd_num) ! temporary array of the size of the number of super-droplets.
   real(RP), intent(out) :: sd_dtmp4(1:sd_num) ! temporary array of the size of the number of super-droplets.
   real(RP), intent(out) :: crs_val1p(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val1c(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val2p(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val2c(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val3p(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val3c(KA,IA,JA) ! temporary buffer of CReSS dimension
   ! Internal shared variables
   integer :: istep_sdm        ! step number of SDM
   integer :: istep_evl        ! step number of {condensation/evaporation} process
   integer :: istep_col        ! step number of {stochastic coalescence} process
   integer :: istep_adv        ! step number of {motion of super-droplets} process
   integer :: istep_mlt        ! step number of {melt/freeze of super-droplets} process
   integer :: istep_sbl        ! step number of {sublimation/deposition of super-droplets} process
   ! Work variables
   integer :: t, n     ! index
   integer :: k,i,j     ! index
   real(RP) :: u_scale(KA,IA,JA)   ! u components of velocity
   real(RP) :: v_scale(KA,IA,JA)   ! v components of velocity
   real(RP) :: w_scale(KA,IA,JA)   ! w components of velocity
   real(RP) :: pres_scale(KA,IA,JA)  ! Pressure
   real(RP) :: t_scale(KA,IA,JA)    ! Temperature 
   real(RP) :: sdm_dtevl  ! time step of {condensation/evaporation} process
   real(RP) :: sdm_dtcol  ! time step of {stochastic coalescence} process
   real(RP) :: sdm_dtadv  ! time step of {motion of super-droplets} process
   real(RP) :: sdm_dtmlt  ! time step of {melt/freeze of super-droplets} process
   real(RP) :: sdm_dtsbl  ! time step of {sublimation/deposition of super-droplets} process
   real(RP) :: tmp_mink
  !---------------------------------------------------------------------

      ! Initialize and rename variables
      sdm_dtevl = real( sdm_dtcmph(1),kind=RP )  !! condensation/evaporation
      sdm_dtcol = real( sdm_dtcmph(2),kind=RP )  !! stochastic coalescence
      sdm_dtadv = real( sdm_dtcmph(3),kind=RP )  !! motion of super-droplets
      sdm_dtmlt = real( sdm_dtcmph(4),kind=RP )  !! melting/freezing
      sdm_dtsbl = real( sdm_dtcmph(5),kind=RP )  !! sublimation/depsition

      istep_sdm = nclstp(0)                !! maximum step
      istep_evl = nclstp(1)                !! condensation/evaporation
      istep_col = nclstp(2)                !! stochastic coalescence
      istep_adv = nclstp(3)                !! motion of super-droplets
      istep_mlt = nclstp(4)                !! motion of super-droplets
      istep_sbl = nclstp(5)                !! motion of super-droplets

      lsdmup = .false.

      ! Calculate super-droplets process.
      !   1 : motion of super-droplets (advection, terminal velocity)
      !   2 : melting / freezing [cold-SDM only]
      !   3 : condensation / evaporation
      !   4 : sublimation / deposition [cold-SDM only]
      !   5 : stochastic coalescence
      do t=1,istep_sdm
         !### Run SDM  ###!

         !=== 1 : motion of super-droplets ===!
#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_start("sdm_motion",0,0)
#endif

         if( sdm_calvar(3) .and.                                 &
             mod(t,istep_sdm/istep_adv)==0 .and. sdm_dtadv>0.d0 ) then

            lsdmup = .true.

            ! get the terminal velocity of super-droplets
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
            !! evaluate the terminal velocity
            call sdm_getvz_liq(pres_scale,DENS,t_scale,            &
                           sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation' )
            if( sdm_cold )then
               call sdm_getvz_ice(DENS,t_scale,            &
                           sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sdi,sd_vz,  &
                           sd_itmp1,'no_interpolation' )
            end if

            ! for predictor-corrector
            sd_dtmp1(:) = sd_x(:)
            sd_dtmp2(:) = sd_y(:)
            sd_dtmp3(:) = sd_rk(:)
            sd_dtmp6(:) = sd_vz(:) 

            ! get the moving velocity of super-droplets
            !! diagnose necessary fluid variables
            call sdm_rho_mom2uvw(DENS,MOMX,MOMY,MOMZ,u_scale,v_scale,w_scale)
            !! update the HALO region of the fluid variables
            call COMM_vars8( u_scale(:,:,:), 1 )
            call COMM_vars8( v_scale(:,:,:), 2 )
            call COMM_vars8( w_scale(:,:,:), 3 )
            call COMM_wait ( u_scale(:,:,:), 1 )
            call COMM_wait ( v_scale(:,:,:), 2 )
            call COMM_wait ( w_scale(:,:,:), 3 )

            !! evaluate the velocity of each SD
            call sdm_getvel(u_scale,v_scale,w_scale,                   &
                            sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_u,sd_v,sd_vz)

            ! momentum coupling
            !! rhow coupling
            if( sdm_mvexchg>=1 ) then

               call sdm_sd2rhow(zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

               do j = JS, JE
                  do i = IS, IE
                     MOMZ(KS-1,i,j) = 0
                     do k = KS, KE-1
                        MOMZ(k,i,j) = MOMZ(k,i,j) - ( crs_val1p(k,i,j) + crs_val1p(k+1,i,j) )*0.5_RP*GRAV*sdm_dtadv
                     enddo
                     MOMZ(KE,i,j) = 0
                  enddo
               enddo

               call COMM_vars8( MOMZ(:,:,:), 1 )
               call COMM_wait ( MOMZ(:,:,:), 1 )

            end if

            !! calculate the current SD momentum field and update fluid
            !!(not supported yet)
            if( sdm_mvexchg>=2 ) then
               ! calculate the current SD momentum field (a)

               ! update fluid using (a) current and (b) previous SD momentum fields

            end if

            ! { motion } in SDM
            !! evaluate the motion eq.
!!$            !!!! explicit Euler scheme
!!$            call sdm_move(sdm_dtadv,                         &
!!$                          sd_num,sd_u,sd_v,sd_vz,sd_x,sd_y,sd_rk)
!!$
            !!!! predictor-corrector scheme
            !!!!!! predictor step
            call sdm_move(sdm_dtadv,                         &
                          sd_num,sd_u,sd_v,sd_vz,sd_dtmp1,sd_dtmp2,sd_dtmp3)
            !!!!!! corrector step
            tmp_mink = real(KS-1,kind=RP)
            sd_dtmp3(:) = max(tmp_mink,sd_dtmp3(:))
            call sdm_getvel(u_scale,v_scale,w_scale,                   &
                            sd_num,sd_dtmp1,sd_dtmp2,sd_ri,sd_rj,sd_dtmp3,sd_dtmp4,sd_dtmp5,sd_dtmp6)
            sd_dtmp4(:) = 0.5_RP*(sd_u(:) +sd_dtmp4(:))
            sd_dtmp5(:) = 0.5_RP*(sd_v(:) +sd_dtmp5(:))
            sd_dtmp6(:) = 0.5_RP*(sd_vz(:)+sd_dtmp6(:))
            call sdm_move(sdm_dtadv,                         &
                          sd_num,sd_dtmp4,sd_dtmp5,sd_dtmp6,sd_x,sd_y,sd_rk)

            ! lateral boundary routine in SDM
            !! judge super-droplets as invalid or valid in horizontal
            !! do MPI communication to send/receiv SDs
            call sdm_boundary(wbc,ebc,sbc,nbc,                           &
                             sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_rk,     &
                             sd_u,sd_v,sd_vz,sd_r,sd_asl,sdi,               &
                             bufsiz1,                                    &
                             bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4,      &
                             sd_itmp1,                              &
                             rbuf_r8,sbuf_r8,rbuf_i8,sbuf_i8,rbuf_i2,sbuf_i2,rbuf_i4,sbuf_i4) 

            ! judge super-droplets as invalid or valid in vartical
            call sdm_jdginvdv(sd_rkl,sd_rku,sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk)

            !! save the SDM momentum field at new position (b)
            !!(not supported yet)
            if( sdm_mvexchg>=2 ) then

            end if

         end if

#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_stop("sdm_motion",0,0)
#endif

         !=== 2 : melt / freeze ===!
#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_start("sdm_meltfreeze",0,0)
#endif
         if( sdm_cold .and. sdm_calvar(4) .and.                                 &
             mod(t,istep_sdm/istep_mlt)==0 .and. sdm_dtmlt>0.d0 ) then

            lsdmup = .true.

            ! get density of solid-water before melt/freeze
            !! cres_val1p is the rhosol before meltfreeze
            call sdm_sd2rhosol(zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

            ! { melting/freezing } in SDM
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
            !! update the phase of SDs
            call sdm_meltfreeze(                              &
              t_scale,pres_scale,QTRC(:,:,:,I_QV),         &
              sd_num,sd_liqice,sd_x,sd_y,      &
              sd_r,sd_ri,sd_rj,sd_rk,sdi )

            ! get density of solid-water after melt/freeze
            !! here cres_val1c is the rhosol after meltfreeze
            call sdm_sd2rhosol(zph_crs,crs_val1c,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2c,sd_itmp1)

            ! exchange the heat to fluid variables
            call sdm_meltfreeze_updatefluid(RHOT,QTRC,DENS,crs_val1p,crs_val1c)
            !! update the HALO region of the fluid variables
            call COMM_vars8( RHOT(:,:,:), 1 )
            call COMM_wait ( RHOT(:,:,:), 1 )

         end if

#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_stop("sdm_meltfreeze",0,0)
#endif

         !=== 3 : condensation / evaporation ===!
#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_start("sdm_condevp",0,0)
#endif
         if( sdm_calvar(1) .and.                                 &
             mod(t,istep_sdm/istep_evl)==0 .and. sdm_dtevl>0.d0 ) then

            lsdmup = .true.

            ! get density of liquid-water(qw) before process-1
            !! cres_val1p is the rhow before condevp
            call sdm_sd2rhow(zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

            ! { condensation/evaporation } in SDM
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
            !! update the equivalent radius of SDs
            call sdm_condevp(sdm_aslset,            &
                             sdm_aslmw,sdm_aslion,sdm_dtevl,      &
                             pres_scale,t_scale,QTRC(:,:,:,I_QV), &
                             sd_num,sd_numasl,sd_liqice,sd_x,sd_y,sd_r,sd_asl,&
                             sd_ri,sd_rj,sd_rk)

            ! get density of liquid-water(qw) after process-1
            !! here cres_val1c is the rhow after condevp
            call sdm_sd2rhow(zph_crs,crs_val1c,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2c,sd_itmp1)

            ! exchange the vapor and heat to fluid variables
            call sdm_condevp_updatefluid(RHOT,QTRC,DENS,crs_val1p,crs_val1c)
            !! update the HALO region of the fluid variables
            call COMM_vars8( RHOT(:,:,:), 1 )
            call COMM_vars8( QTRC(:,:,:,I_QV), 2 )
            call COMM_vars8( DENS(:,:,:), 3 )
            call COMM_wait ( RHOT(:,:,:), 1 )
            call COMM_wait ( QTRC(:,:,:,I_QV), 2 )
            call COMM_wait ( DENS(:,:,:), 3 )

         end if

#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_stop("sdm_condevp",0,0)
#endif

         !=== 4 : sublimation / deposition ===!
#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_start("sdm_subldep",0,0)
#endif
         if( sdm_cold .and. sdm_calvar(5) .and.                                 &
             mod(t,istep_sdm/istep_sbl)==0 .and. sdm_dtsbl>0.d0 ) then

            lsdmup = .true.

            ! get density of solid-water before sublimation/deposition
            !! cres_val1p is the rhosol before sublimation/deposition
            call sdm_sd2rhosol(zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

            ! { sublimation/deposition } in SDM
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
            !! evaluate the terminal velocity
            call sdm_getvz_ice(DENS,t_scale,            &
                           sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sdi,sd_vz,  &
                           sd_itmp1,'no_interpolation' )
            !! update the equivalent radius of SDs
            call sdm_subldep(            &
                             sdm_dtsbl,      &
                             pres_scale,t_scale,QTRC(:,:,:,I_QV),DENS,&
                             sd_num,sd_liqice,sd_x,sd_y,sdi,sd_vz,&
                             sd_ri,sd_rj,sd_rk,sd_itmp1)

            ! get density of solid-water after sublimation/deposition
            !! here cres_val1c is the rhosol after sublimation/deposition
            call sdm_sd2rhosol(zph_crs,crs_val1c,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2c,sd_itmp1)

            ! exchange the vapor and heat to fluid variables
            call sdm_subldep_updatefluid(RHOT,QTRC,DENS,crs_val1p,crs_val1c)
            !! update the HALO region of the fluid variables
            call COMM_vars8( RHOT(:,:,:), 1 )
            call COMM_vars8( QTRC(:,:,:,I_QV), 2 )
            call COMM_vars8( DENS(:,:,:), 3 )
            call COMM_wait ( RHOT(:,:,:), 1 )
            call COMM_wait ( QTRC(:,:,:,I_QV), 2 )
            call COMM_wait ( DENS(:,:,:), 3 )

         end if

#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_stop("sdm_subldep",0,0)
#endif

         !=== 5 : stochastic coalescence ===!
#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_start("sdm_coales",0,0)
#endif

         if( sdm_calvar(2) .and.                                 &
             mod(t,istep_sdm/istep_col)==0 .and. sdm_dtcol>0.d0 ) then

            lsdmup = .true.

            ! get the terminal velocity of super-droplets
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
            !! evaluate the terminal velocity
            call sdm_getvz_liq(pres_scale,DENS,t_scale,            &
                           sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation' )
            if( sdm_cold )then
               call sdm_getvz_ice(DENS,t_scale,            &
                           sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sdi,sd_vz,  &
                           sd_itmp1,'no_interpolation' )
            end if

            ! { coalescence } in SDM
            if( sdm_cold )then
               ! get density of solid-water before riming
               !! cres_val1p is the rhosol before riming
               call sdm_sd2rhosol(zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

               call sdm_coales_cold(sdm_colkrnl,sdm_colbrwn,sdm_aslset,         &
                            sdm_aslrho,sdm_dtcol,                       &
                            pres_scale,t_scale,QTRC(:,:,:,I_QV),DENS,&
                            zph_crs,                                    &
                            ni_sdm,nj_sdm,nk_sdm,sd_num,sd_numasl,      &
                            sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,     &
                            sdi,                                        & 
                            sort_id,sort_key,sort_freq,sort_tag,        &
                            sd_rng,sd_rand,                             &
                            sdm_itmp1,sdm_itmp2,                        &
                            sd_itmp1(1:sd_num),sd_itmp2(1:sd_num),  &
                            sd_dtmp1)

               ! get density of solid-water after riming
               !! here cres_val1c is the rhosol after riming
               call sdm_sd2rhosol(zph_crs,crs_val1c,                       &
                             sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2c,sd_itmp1)

               ! exchange the heat to fluid variables
               call sdm_meltfreeze_updatefluid(RHOT,QTRC,DENS,crs_val1p,crs_val1c)
               !! update the HALO region of the fluid variables
               call COMM_vars8( RHOT(:,:,:), 1 )
               call COMM_wait ( RHOT(:,:,:), 1 )

            else
               call sdm_coales(sdm_colkrnl,sdm_colbrwn,sdm_aslset,         &
                            sdm_aslrho,sdm_dtcol,                       &
                            pres_scale, t_scale,                        &
                            zph_crs,                                    &
                            ni_sdm,nj_sdm,nk_sdm,sd_num,sd_numasl,      &
                            sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,     &
                            sort_id,sort_key,sort_freq,sort_tag,        &
                            sd_rng,sd_rand,                             &
                            sdm_itmp1,sdm_itmp2,                        &
                            sd_itmp1(1:sd_num),sd_itmp2(1:sd_num),  &
                            sd_dtmp1)
            end if

         end if

#ifdef _FAPP_
         ! Section specification for fapp profiler
         call fapp_stop("sdm_coales",0,0)
#endif

      end do

    ! Convert super-droplets to precipitation

#ifdef _FAPP_
      ! Section specification for fapp profiler
      call fapp_start("sdm_sd2prec",0,0)
#endif
      call sdm_sd2prec(dt,                                &
                       prec_crs,                          &
                       sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sdi,sd_ri,sd_rj,sd_rk,  &
                       sd_itmp1,sd_itmp2,crs_val1c(1,1:IA,1:JA))

#ifdef _FAPP_
      ! Section specification for fapp profiler
      call fapp_stop("sdm_sd2prec",0,0)
#endif

    return
  end subroutine sdm_calc
!---------------------------------------------------------------------------------------
!  subroutine sdm_sd2momnt(ni,nj,nk,zph_crs,mu_sdm,mv_sdm,mw_sdm,  &
!                          sd_num,sd_n,sd_x,sd_y,sd_rk,            &
!                          sd_u,sd_v,sd_wc,sd_r,ilist)
!     use scale_const, only: &
!        rw => CONST_DWATR
!     use scale_grid, only: &
!        CX => GRID_CX, &
!        CY => GRID_CY
!     ! Input variables
!     integer, intent(in) :: ni  ! Model dimension in x direction
!     integer, intent(in) :: nj  ! Model dimension in y direction
!     integer, intent(in) :: nk  ! Model dimension in z direction
!     real(RP),intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinate
!     integer, intent(in) :: sd_num  ! number of super-droplets
!     integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
!     real(RP),intent(in) :: sd_x(1:sd_num)    ! x-coordinate of super-droplets
!     real(RP), intent(in) :: sd_y(1:sd_num)    ! y-coordinate of super-droplets
!     real(RP), intent(in) :: sd_rk(1:sd_num)   ! index[k/real] of super-droplets
!     real(RP), intent(in) :: sd_u(1:sd_num)    ! x-direction velocity of super-droplets
!     real(RP), intent(in) :: sd_v(1:sd_num)    ! y-direction velocity of super-droplets
!     real(RP), intent(in) :: sd_wc(1:sd_num)   ! zeta components of contravariant velocity of super-droplets
!     real(RP), intent(in) :: sd_r(1:sd_num)    ! equivalent radius of super-droplets
!     ! Output variables
!     real(RP), intent(out) :: mu_sdm(KA,IA,JA) ! a total of momentum of super-droplets in x-direction
!     real(RP), intent(out) :: mv_sdm(KA,IA,JA) ! a total of momentum of super-droplets in y-direction
!     real(RP), intent(out) :: mw_sdm(KA,IA,JA) ! a total of momentum of super-droplets in z-direction
!     integer,  intent(out) :: ilist(1:int(sd_num/nomp),1:nomp) ! buffer for list vectorization
!     ! Internal shared variables
!!     real(RP) :: dx_sdm              ! dx_sdm
!!     real(RP) :: dy_sdm              ! dy_sdm
!     ! Work variables for OpenMP
!     integer :: sd_str           ! index of divided loop by OpenMP
!     integer :: sd_end           ! index of divided loop by OpenMP
!     integer :: np               ! index for OpenMP
!     ! Work variables
!     real(RP) :: dcoef(KA,IA,JA)      ! coef.
!     real(RP) :: dist       ! temporary
!     real(RP) :: dtmp       ! temporary
!!     integer(DP) :: idx_sdm ! integer of 'dx_sdm'
!!     integer(DP) :: idy_sdm ! integer of 'dy_sdm'
!     integer :: nlist(1:nomp)    ! list number
!     integer :: cnt              ! counter
!     integer :: i, j, k, m, n    ! index
!     integer :: ix, jy           ! index
!    !---------------------------------------------------------------------
!
!      ! Initialize
!      do k = KS, KE
!      do i = IS, IE
!      do j = JS, JE
!       dcoef(k,i,j) = rw * F_THRD*ONE_PI / real(dx_sdm(i)*dy_sdm(i),kind=RP)
!!       idx_sdm = 10_i8 * floor( 1.e4*(dx_sdm+1.e-5), kind=i8 )
!!       idy_sdm = 10_i8 * floor( 1.e4*(dy_sdm+1.e-5), kind=i8 )
!      enddo
!      enddo
!      enddo
!
!      do np=1,nomp
!         nlist(np) = 0
!      end do
!
!      do k=1,nk
!      do j=0,nj+1
!      do i=0,ni+1
!         mu_sdm(k,i,j) = 0.0_RP
!         mv_sdm(k,i,j) = 0.0_RP
!         mw_sdm(k,i,j) = 0.0_RP
!      end do
!      end do
!      end do
!
!      ! Get index list for compressing buffer.
!      do np=1,nomp
!
!         sd_str = int(sd_num/nomp)*(np-1) + 1
!         sd_end = int(sd_num/nomp)*np
!         cnt    = 0
!         do n=sd_str,sd_end
!             if( sd_rk(n)>VALID2INVALID ) then
!               cnt = cnt + 1
!               ilist(cnt,np) = n
!            end if
!         end do
!         nlist(np) = cnt
!      end do
!
!      ! Get a total of momentum of super-droplets.
!
!      !### count momentum of super-droplets ###!
!      do np=1,nomp
!
!         if( nlist(np)>0 ) then
!            do m=1,nlist(np)
!               n = ilist(m,np)
!
!               dtmp = real(sd_n(n),kind=DP) * sd_r(n)*sd_r(n)*sd_r(n)
!
!               !=== momentum in x-direction ===!
!
!!               dist = sd_x(n) + 0.50_RP * real(dx_sdm,kind=RP)
!               dist = sd_x(n) + 0.50_RP * real(dx_sdm,kind=RP)
!               do ix = IS, IE
!!                if( sd_x(n) < CX(ix) ) then
!                if( dist < CX(ix) ) then
!                 i = ix
!                 exit
!                endif
!               enddo
!               do jy = JS, JE
!                if( sd_y(n) < CY(jy) ) then
!                 j = jy
!                 exit
!                endif
!               enddo
!!               i = int( floor( dist*1.d5   , kind=i8 )/idx_sdm ) + 2
!!               j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm ) + 2
!               k = floor( sd_rk(n) )
!
!               mu_sdm(i,j,k) = mu_sdm(i,j,k) + sd_u(n) * dtmp
!
!               !=== momentum in y-direction ===!
!
!               dist = sd_y(n) + 0.50_RP * real(dy_sdm,kind=RP)
!               do ix = IS, IE
!                if( sd_x(n) < CX(ix) ) then
!                 i = ix
!                 exit
!                endif
!               enddo
!               do jy = JS, JE
!!                if( sd_y(n) < CY(jy) ) then
!                if( dist < CY(jy) ) then
!                 j = jy
!                 exit
!                endif
!               enddo
!!               i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm ) + 2
!!               j = int( floor( dist*1.d5   , kind=i8 )/idy_sdm ) + 2
!               k = floor( sd_rk(n) )
!
!               mv_sdm(k,i,j) = mv_sdm(k,i,j) + sd_v(n) * dtmp
!
!               !=== momentum in zeta-direction ===!
!
!               do ix = IS, IE
!                if( sd_x(n) < CX(ix) ) then
!                 i = ix
!                 exit
!                endif
!               enddo
!               do jy = JS, JE
!                if( sd_y(n) < CY(jy) ) then
!                 j = jy
!                 exit
!                endif
!               enddo
!!               i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm ) + 2
!!               j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm ) + 2
!               k = floor( sd_rk(n) + 0.50_RP )
!               mw_sdm(k,i,j) = mw_sdm(k,i,j) + sd_wc(n) * dtmp
!
!            end do
!
!         end if
!
!      end do
!
!      ! Convert super-droplets to density of liquid-water.
!      do k=2,nk-1
!      do j=1,nj-1
!      do i=1,ni-1
!
!         dtmp = dcoef(k,i,j)/real(zph_crs(k+1,i,j)-zph_crs(k,i,j),kind=RP)
!         mu_sdm(k,i,j) = mu_sdm(k,i,j) * dtmp
!         mv_sdm(k,i,j) = mv_sdm(k,i,j) * dtmp
!         mw_sdm(k,i,j) = mw_sdm(k,i,j) * dtmp
!
!      end do
!      end do
!      end do
!
!    return
!  end subroutine sdm_sd2momnt
  !-----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  subroutine sdm_aslform(DENS,RHOT,QTRC,                          &   
                         sdm_calvar,sdm_aslset,                   &
                         sdm_aslfmsdnc,sdm_sdnmlvol,              &
                         sdm_zupper,sdm_zlower,sdm_dtcmph,        &
                         pbr_crs,ptbr_crs,pp_crs,         &
                         ptp_crs,qv_crs,zph_crs,rhod_crs,         &
                         sd_num,sd_numasl,sd_n,sd_x,sd_y,sd_z,    &
                         sd_ri,sd_rj,sd_rk,sd_u,sd_v,sd_vz,sd_r,sd_asl,       &
                         sd_fmnum,sd_fmn,sd_fmliqice,sd_fmx,sd_fmy,sd_fmz,    &
                         sd_fmri,sd_fmrj,sd_fmrk,sd_fmvz,sd_fmr,sd_fmasl,         &
                         ni_sdm,nj_sdm,nk_sdm,sort_id,sort_key,   &
                         sort_freq,sort_tag,sd_rng,               &
!                         sort_freq,sort_tag,                      &
                         sort_tag0,sd_itmp1,sd_itmp2,sd_itmp3)

    use scale_tracer, only: &
         QAD => QA
    use m_sdm_coordtrans, only: sdm_z2rk
    use m_sdm_fluidconv, only: &
         sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw
    use m_sdm_motion, only: &
         sdm_getvz_liq
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_idutil, only: &
         sdm_sort
    use m_sdm_condensation_water, only: &
        sdm_condevp
      ! Input variables
    real(RP), intent(in) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(in) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(in) :: QTRC(KA,IA,JA,QAD)    !! Ratio of mass of tracer to total mass[kg/kg]
      logical, intent(in) :: sdm_calvar(5) ! Control flag of calculation using SDM
      integer, intent(in) :: sdm_aslset    ! Option for aerosol species
      real(RP),intent(in) :: sdm_sdnmlvol  ! Normal volume for number concentration of super droplets
      real(RP),intent(in) :: sdm_aslfmsdnc ! Number of super droplets at aerosol formation per sdm_sdnmlvol
      real(RP),intent(in) :: sdm_zlower    ! Lower limitaion of initial droplet's position
      real(RP),intent(in) :: sdm_zupper    ! Upper limitaion of initial droplet's position
      real(RP),intent(in) :: sdm_dtcmph(1:5)     ! Time interval of cloud micro physics
!      real(RP), intent(in) :: jcb_crs(KA,IA,JA)   ! Jacobian at scalar points
      real(RP), intent(in) :: pbr_crs(KA,IA,JA)   ! Base state pressure
      real(RP), intent(in) :: ptbr_crs(KA,IA,JA)  ! Base state potential temperature
      real(RP), intent(in) :: pp_crs(KA,IA,JA)    ! Pressure perturbation
      real(RP), intent(in) :: ptp_crs(KA,IA,JA)   ! Potential temperature perturbation
      real(RP), intent(in) :: qv_crs(KA,IA,JA)    ! Water vapor mixing ratio
      real(RP), intent(in) :: zph_crs(KA,IA,JA)   ! Z-physical coordinates
      real(RP), intent(in) :: rhod_crs(KA,IA,jA)   ! dry air densiy
      integer, intent(in) :: ni_sdm   ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm   ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm   ! SDM model dimension in z direction
      integer, intent(in) :: sd_num   ! number of super-droplets
      integer, intent(in) :: sd_numasl! number of kind of chemical material contained as water-soluble aerosol in super droplets
      integer, intent(in) :: sd_fmnum ! number of super-droplets at aerosol formation
      ! Input and output variables
      integer(DP), intent(inout) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_x(1:sd_num)      ! x-coordinate of super-droplets
      real(RP), intent(inout) :: sd_y(1:sd_num)      ! y-coordinate of super-droplets
      real(RP), intent(inout) :: sd_z(1:sd_num)      ! z-coordinate of super-droplets
      real(RP), intent(inout) :: sd_ri(1:sd_num)     ! index[i/real] of super-droplets
      real(RP), intent(inout) :: sd_rj(1:sd_num)     ! index[j/real] of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num)     ! index[k/real] of super-droplets
      real(RP), intent(inout) :: sd_u(1:sd_num)      ! x-components velocity of super-droplets
      real(RP), intent(inout) :: sd_v(1:sd_num)      ! y-components velocity of super-droplets
      real(RP), intent(inout) :: sd_vz(1:sd_num)     ! terminal velocity of super-droplets
      real(RP), intent(inout) :: sd_r(1:sd_num)      ! equivalent radius of super-droplets
      real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl)  ! aerosol mass of super-droplets
      integer(DP), intent(inout) :: sd_fmn(1:sd_fmnum) ! multiplicity of super-droplets at aerosol formation
      integer(i2), intent(inout) :: sd_fmliqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
      real(RP), intent(inout) :: sd_fmx(1:sd_fmnum)  ! x-coordinate of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmy(1:sd_fmnum)  ! y-coordinate of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmz(1:sd_fmnum)  ! z-coordinate of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmri(1:sd_fmnum) ! index[i/real] of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmrj(1:sd_fmnum) ! index[j/real] of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmrk(1:sd_fmnum) ! index[k/real] of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmvz(1:sd_fmnum) ! terminal velocity of super-droplet at aerosol formation
      real(RP), intent(inout) :: sd_fmr(1:sd_fmnum)  ! equivalent radius of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmasl(1:sd_fmnum,1:sd_numasl)  ! aerosol mass of super-droplets at aerosol formation
      type(c_rng_uniform_mt), intent(inout) :: sd_rng     ! random number generator
      integer, intent(inout) :: sort_id(1:sd_num)         ! super-droplets sorted by SD-grids
      integer, intent(inout) :: sort_key(1:sd_num)        ! sort key
      integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1) ! number of super-droplets in each SD-grid
      integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)  ! accumulated number of super-droplets in each SD-grid
      ! Output variables
      integer, intent(out) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)   ! = sort_tag(n) - 1
      integer, intent(out) :: sd_itmp1(1:sd_num)  ! temporary array of the size of the number of super-droplets.
      integer, intent(out) :: sd_itmp2(1:sd_num)  ! temporary array of the size of the number of super-droplets.
      integer, intent(out) :: sd_itmp3(1:sd_num)  ! temporary array of the size of the number of super-droplets.
      ! Work variables
      real(RP) :: sd_fmnc      ! number concentration of super-droplets at aerosol formation
      integer :: max_key       ! max key to sort data
      integer :: n_invd        ! number of invalid droplets
      integer :: innum         ! temporary
      integer :: id_invd       ! index
      integer :: k, n

      real(RP) :: pres_scale(KA,IA,JA)  ! Pressure
      real(RP) :: t_scale(KA,IA,JA)    ! Temperature
      real(RP) :: sdm_dtevl  ! time step of {condensation/evaporation} process
      real(RP) :: sdm_dtcol  ! time step of {stochastic coalescence} process
      real(RP) :: sdm_dtadv  ! time step of {motion of super-droplets} process
      real(RP) :: sdm_dtmlt  ! time step of {melt/freeze of super-droplets} process
      real(RP) :: sdm_dtsbl  ! time step of {sublimation/deposition of super-droplets} process 
    !--------------------------------------------------------------------- !---------------------------------------------------------------------

      ! Initialize and rename variables
      sdm_dtevl = real( sdm_dtcmph(1),kind=RP )  !! condensation/evaporation
      sdm_dtcol = real( sdm_dtcmph(2),kind=RP )  !! stochastic coalescence
      sdm_dtadv = real( sdm_dtcmph(3),kind=RP )  !! motion of super-droplets
      sdm_dtmlt = real( sdm_dtcmph(4),kind=RP )  !! melting/freezing
      sdm_dtsbl = real( sdm_dtcmph(5),kind=RP )  !! sublimation/depsition


      call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
      call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

      ! Check active

      if( abs(sdm_aslset)<10 ) return
      ! Sorting valid super-droplets

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
                    sort_id,sort_key,sort_freq,sort_tag,'valid')

      max_key = ni_sdm * nj_sdm * knum_sdm + 1

      n_invd = sort_freq(max_key)    !! number of invalid droplets

      do n=1,max_key+1
         sort_tag0(n) = sort_tag(n) - 1  !! {1-xx} => {0-xx}
      end do

      ! Initialize


      ! Get random number

      call gen_rand_array( sd_rng, sd_fmx  )
      call gen_rand_array( sd_rng, sd_fmy  )
      call gen_rand_array( sd_rng, sd_fmz  )
      call gen_rand_array( sd_rng, sd_fmvz )

      do k=1,sd_numasl

         call gen_rand_array( sd_rng, sd_fmr )

         do n=1,sd_fmnum
            sd_fmasl(n,k) = sd_fmr(n)   !! sd_fmr : temporary
         end do

      end do

      ! Set aerosol mass and muliplicity of super-droplets

      !### Ammonium Sulfate [(NH4)2SO4] ###!

      sd_fmnc = real(sdm_aslfmsdnc,kind=RP)/real(sdm_sdnmlvol,kind=RP)
                                     !! number concentration of super-
                                     !! -droplets at aerosol formation

      call sdm_aslsulf(sdm_aslfmrate,sdm_aslfmdt,                &
                       sd_fmnum,sd_numasl,sd_fmn,sd_fmasl,sd_fmnc)

      ! Set position of super-droplets
      do n=1,sd_fmnum

         sd_fmx(n) = xmax_sdm * sd_fmx(n)
         sd_fmy(n) = ymax_sdm * sd_fmy(n)
         sd_fmz(n) = real(minzph+sdm_zlower,kind=RP) + sd_fmz(n)        &
                   * real(sdm_zupper-(minzph+sdm_zlower),kind=RP)

         sd_fmr(n) = 1.0E-15_RP

      end do

      !### index[k/real] of super-droplets ###!

      call sdm_z2rk(sdm_zlower,sdm_zupper,       &
                        sd_fmnum,sd_fmx,sd_fmy,sd_fmz,sd_fmri,sd_fmrj,sd_fmrk)

      ! Set radius of super-droplets
      if( sdm_calvar(1) ) then

         call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)

         call sdm_condevp(sdm_aslset,                              &
                          sdm_aslmw,sdm_aslion,sdm_dtevl,           &
                          pres_scale,t_scale,QTRC(:,:,:,I_QV),      &
                          sd_fmnum,sd_numasl,sd_fmliqice,sd_fmx,sd_fmy,sd_fmr,      &
                          sd_fmasl,sd_fmri,sd_fmrj,sd_fmrk)

      end if

     ! Set terminal velocity of super-droplets

      call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)

      call sdm_getvz_liq(pres_scale,DENS,t_scale,                         &
                     sd_fmnum,sd_fmliqice,sd_fmx,sd_fmy,sd_fmri,sd_fmrj,sd_fmrk,sd_fmr,sd_fmvz,   &
                     sd_itmp1(1:sd_fmnum),sd_itmp2(1:sd_fmnum),   &
                     sd_itmp3(1:sd_fmnum),'no_interpolation')

      ! Add new droplets as formed aerosol

      !### adjust number of super-droplets at aerosol formation ###!

      if( n_invd<sd_fmnum ) then

        if( IO_L ) then
         write(IO_FID_LOG,'(2a)')"  ### [SDM] : warning for buffer size limit",  &
          " @ stop to form super-droplets as aerosol ###"
        endif

         innum = nomp * int(n_invd/nomp)

      else

         innum = sd_fmnum

      end if

      do n=1,innum

         id_invd = sort_id( sort_tag0(max_key) + n )

         sd_n(id_invd)  = sd_fmn(n)

         sd_x(id_invd)  = sd_fmx(n)
         sd_y(id_invd)  = sd_fmy(n)
         sd_z(id_invd)  = sd_fmz(n)
         sd_rk(id_invd) = sd_fmrk(n)

         sd_u(id_invd)  = 0.0_RP
         sd_v(id_invd)  = 0.0_RP
         sd_vz(id_invd) = sd_fmvz(n)

         sd_r(id_invd)  = sd_fmr(n)

      end do

      do k=1,sd_numasl

         do n=1,innum

            id_invd = sort_id( sort_tag0(max_key) + n )

            sd_asl(id_invd,k) = sd_fmasl(n,k)

         end do

      end do


    return
  end subroutine sdm_aslform
  !----------------------------------------------------------------------------
  subroutine sdm_adjsdnum(sdm_nadjvar,ni_sdm,nj_sdm,nk_sdm,    &
                          sd_num,sd_numasl,sd_nc,                 &
                          sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk, &
                          sort_id,sort_key,sort_freq,sort_tag,    &
                          sd_rng,sd_rand,                         &
!                          sd_rand,                                &
                          sdm_itmp1,sdm_itmp2,sd_itmp1)
      use scale_process, only: &
           mype => PRC_myrank
      ! Input variables
      integer, intent(in) :: sdm_nadjvar
      integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
      integer, intent(in) :: sd_num  ! number of super-droplets
      integer, intent(in) :: sd_numasl  ! number of kind of chemical material contained as water-soluble aerosol in super droplets
      real(RP), intent(in) :: sd_nc ! averaged number concentration in a grid
      ! Input and output variables
      type(c_rng_uniform_mt), intent(inout) :: sd_rng  ! random number generator
      real(RP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
      integer, intent(inout) :: sort_id(1:sd_num)  ! super-droplets sorted by SD-grids
      integer, intent(inout) :: sort_key(1:sd_num) ! sort key
      integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1)  ! number of super-droplets in each SD-grid
      integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)   ! accumulated number of super-droplets in each SD-grid
      real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
      real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
      integer(DP), intent(inout) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
      real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
      real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
      real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
      real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
      ! Output variables
      integer, intent(out) :: sdm_itmp1(1:ni_sdm*nj_sdm*nk_sdm+2)  ! temporary array of SDM dimension
      integer, intent(out) :: sdm_itmp2(1:ni_sdm*nj_sdm*nk_sdm+2)  ! temporary array of SDM dimension
      integer, intent(out) :: sd_itmp1(1:sd_num)  ! temporary array of the size of the number of super-droplets.
      ! Parameter
      real(RP), parameter :: RATE4REMOVE = 1.40_RP ! upper limit rate to averaged number concentration in a grid for removing
      real(RP), parameter :: RATE4ADD  = 0.70_RP  ! lower limit rate to averaged number concentration in a grid for adding
      ! Work variables
      integer :: sdnum_upr ! upper limit number of super-droplets in a grid for removing
      integer :: sdnum_lwr ! lower limit number of super-droplets in a grid for adding
      !------------------------------------------------------------------7--

      if( sdm_nadjvar==0 ) return

      ! Remove super-droplets from grids with large number of super-droplets

      if( sdm_nadjvar==3 .or. sdm_nadjvar==2 ) then

         sdnum_upr = floor( RATE4REMOVE * sd_nc )
         call sdm_sdremove(ni_sdm,nj_sdm,nk_sdm,                        &
                           sdnum_upr,sd_num,sd_n,sd_x,sd_y,sd_ri,sd_rj,sd_rk,       &
                           sort_id,sort_key,sort_freq,sort_tag,         &
                           sd_rng,sd_rand,                              &
!                           sd_rand,                                     &
                           sdm_itmp1,sdm_itmp2,sd_itmp1)

      end if

      ! Add super-droplets to grids with small number of super-droplets

      if( sdm_nadjvar==3 .or. sdm_nadjvar==1 ) then

         sdnum_lwr = floor( RATE4ADD * sd_nc )

         call sdm_sdadd(ni_sdm,nj_sdm,nk_sdm,                           &
                        sdnum_lwr,sd_num,sd_numasl,                     &
                        sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,         &
                        sort_id,sort_key,sort_freq,sort_tag,            &
                        sd_rng,sd_rand,                                 &
!                        sd_rand,                                        &
                        sdm_itmp1,sdm_itmp2,sd_itmp1)

      end if

     ! Message for adjustment

      if( mype==0 .and. sdm_nadjvar==3 ) then
        if( IO_L ) then
         write(IO_FID_LOG,'(a,a,i4,a,i4,a)')                           &
                 "  ### [SDM] : adjust number of super-droplet ",      &
                       "( min",sdnum_lwr," -- max",sdnum_upr," ) ###"
        endif
      else if( mype==0 .and. sdm_nadjvar==1 ) then

        if( IO_L ) then
         write(IO_FID_LOG,'(a,a,i4,a)')                                &
                 "  ### [SDM] : adjust number of super-droplet ",      &
                           "( min",sdnum_lwr," -- max INFINITY ) ###"
        endif

      else if( mype==0 .and. sdm_nadjvar==2 ) then
        if( IO_L ) then
         write(IO_FID_LOG,'(a,a,i4,a)')                                &
                 "  ### [SDM] : adjust number of super-droplet ",      &
                                  "( min 0 -- max",sdnum_upr," ) ###"
        endif

      end if


    return
  end subroutine sdm_adjsdnum
  !----------------------------------------------------------------------------
  subroutine sdm_sdremove(ni_sdm,nj_sdm,nk_sdm,                   &
                          sdnum_upr,sd_num,sd_n,sd_x,sd_y,sd_ri,sd_rj,sd_rk,  &
                          sort_id,sort_key,sort_freq,sort_tag,    &
                          sd_rng,sd_rand,                         &
!                          sd_rand,                                &
                          sort_tag0,fsort_id,isd_perm)

!      use m_gadg_algorithm, only: &
!          gadg_count_sort
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_idutil, only: &
         sdm_sort,sdm_getperm
      ! Input variables
      integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
      integer, intent(in) :: sdnum_upr ! upper limit number of super-droplets in a grid  for removing
      integer, intent(in) :: sd_num    ! number of super-droplets
      ! Input and output variables
      type(c_rng_uniform_mt), intent(inout) :: sd_rng    ! random number generator
      real(RP), intent(inout) :: sd_rand(1:sd_num)  ! random numbers
      integer, intent(inout) :: sort_id(1:sd_num)        ! super-droplets sorted by SD-grids
      integer, intent(inout) :: sort_key(1:sd_num)       ! sort key
      integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1)  ! number of super-droplets in each SD-grid
      integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)   ! accumulated number of super-droplets in each SD-grid
      integer(DP), intent(inout) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_x(1:sd_num)     ! x-coordinate of super-droplets
      real(RP), intent(inout) :: sd_y(1:sd_num)     ! y-coordinate of super-droplets
      real(RP), intent(inout) :: sd_ri(1:sd_num)    ! index[i/real] of super-droplets
      real(RP), intent(inout) :: sd_rj(1:sd_num)    ! index[j/real] of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num)    ! index[k/real] of super-droplets
      ! Output variables
      integer, intent(out) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)  ! = sort_tag(n) - 1
      integer, intent(out) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
      integer, intent(out) :: isd_perm(1:sd_num,1:nomp)   ! random permutations
      ! Internal shared variables
      integer :: freq_max   ! get the maximum number of super-droplets in each grid
      ! Work variables
      real(RP) :: drate   ! temporary
      integer, allocatable :: fsort_tag(:)  ! buffer for sorting
      integer, allocatable :: fsort_freq(:) ! buffer for sorting
      integer :: gnum          ! grid number
      integer :: id_vd         ! index of valid super-droplets
      integer :: n_minus       ! number of reducing super-droplets
      integer :: sdnum_valid   ! number of valid super-droplets in each grid
      integer :: ip, m, n, s, t            ! index
     !---------------------------------------------------------------------

      call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
      call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

      ! Initialize

      gnum = ni_sdm * nj_sdm * knum_sdm

      freq_max = 1

      ! Sorting valid super-droplets

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
                    sort_id,sort_key,sort_freq,sort_tag,'valid')

      ! Initialize
      do n=1,gnum+2
         sort_tag0(n) = sort_tag(n) - 1  !! {1-xx} => {0-xx}
      end do

      ! Get the maximum number of super-droplets in each grid
      do m=1,gnum
         freq_max = max( freq_max, sort_freq(m) )
      end do


!write(*,*) sdnum_upr, freq_max, sort_freq(gnum+1)
      ! Sorting the grids by the number of super-droplets in each grid

      allocate( fsort_tag(0:freq_max+1) )
      allocate( fsort_freq(0:freq_max)  )

      call gadg_count_sort( sort_freq(1:gnum), 0, freq_max,             &
                            fsort_freq, fsort_tag, fsort_id )

      fsort_tag(freq_max+1) = fsort_tag(freq_max) + fsort_freq(freq_max)

      ! Get random number using random number generator

      call gen_rand_array( sd_rng, sd_rand )

      ! Get random permutation layout of super-droplets in each grid

      call sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,            &
                       sort_tag0,fsort_tag,fsort_id,sd_rand,isd_perm)

      ! Remove super-droplets from grids with large number of super-droplets

      do s=1,freq_max
         do t=fsort_tag(sdnum_upr+1),fsort_tag(freq_max+1)-1

            m = fsort_id(t)

            !### get number of removable droplets in each grid ###!
            sdnum_valid = sort_freq(m)   !! number of valid droplets

            n_minus = min( int(sdnum_valid/2),                          &
                           max(sdnum_valid-sdnum_upr,0) )

            if( s>sdnum_valid ) cycle    !! cycle at more than number of valid droplets

            !### get index of valid droplets ###!
            ip = isd_perm( sort_tag0(m)  + s ,1)   !! search forward
!SELECT     ip = isd_perm( sort_tag(m+1) - s )   !! search backward

            id_vd = sort_id( sort_tag0(m) + ip )

            !### Remove droplets ###!

            if( s<=n_minus ) then

               !== convert valid droplets to invalid droplet ==!

               sd_rk(id_vd) = INVALID

            else if( s>n_minus .and. s<=sdnum_valid ) then

               !== adjust multiplicity of valid droplets ==!

               drate = real(sdnum_valid,kind=RP)                        &
                                      /real(sdnum_valid-n_minus,kind=RP)

               sd_n(id_vd) = nint( real(sd_n(id_vd),kind=DP)*drate      &
                                                            , kind=DP )

            end if

         end do

      end do


     ! Deallocate

      deallocate( fsort_tag  )
      deallocate( fsort_freq )

    return
  end subroutine sdm_sdremove
  !----------------------------------------------------------------------------
  subroutine sdm_sdadd(ni_sdm,nj_sdm,nk_sdm,                      &
                       sdnum_lwr,sd_num,sd_numasl,                &
                       sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,    &
                       sort_id,sort_key,sort_freq,sort_tag,       &
                       sd_rng,sd_rand,                            &
!                       sd_rand,                                   &
                       sort_tag0,fsort_id,isd_perm)
!      use m_gadg_algorithm, only: &
!          gadg_count_sort
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_idutil, only: &
         sdm_sort,sdm_getperm
      ! Input variables
      integer, intent(in) :: ni_sdm ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm ! SDM model dimension in z direction
      integer, intent(in) :: sdnum_lwr  ! lower limit number of super-droplets in a grid for adding
      integer, intent(in) :: sd_num ! number of super-droplets
      integer, intent(in) :: sd_numasl  ! number of kind of chemical material contained as water-soluble aerosol in super droplets
      ! Input and output variables
      type(c_rng_uniform_mt), intent(inout) :: sd_rng   ! random number generator
      real(RP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
      integer, intent(inout) :: sort_id(1:sd_num)       ! super-droplets sorted by SD-grids
      integer, intent(inout) :: sort_key(1:sd_num)      ! sort key
      integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1)   ! number of super-droplets in each SD-grid
      integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)    ! accumulated number of super-droplets in each SD-grid
      real(RP), intent(inout) :: sd_x(1:sd_num)    ! x-coordinate of super-droplets
      real(RP), intent(inout) :: sd_y(1:sd_num)    ! y-coordinate of super-droplets
      integer(DP), intent(inout) :: sd_n(1:sd_num) ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_r(1:sd_num)    ! equivalent radius of super-droplets
      real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl)  ! aerosol mass of super-droplets
      real(RP), intent(inout) :: sd_vz(1:sd_num)   ! terminal velocity of super-droplets
      real(RP), intent(inout) :: sd_ri(1:sd_num)   ! index[i/real] of super-droplets
      real(RP), intent(inout) :: sd_rj(1:sd_num)   ! index[j/real] of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num)   ! index[k/real] of super-droplets
      ! Output variables
      integer, intent(out) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2) ! = sort_tag(n) - 1
      integer, intent(out) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
      integer, intent(out) :: isd_perm(1:sd_num,1:nomp)  ! random permutations
      ! Internal shared variables
      integer :: freq_max ! get the maximum number of super-droplets in each grid
      ! Work variables
      integer, allocatable :: fsort_tag(:)  ! buffer for sorting
      integer, allocatable :: fsort_freq(:) ! buffer for sorting
      integer, allocatable :: sort_freq_valid(:) ! number of valid super-droplets in each grid
      integer :: idx_nasl(1:20) ! index for vactorization
      integer :: gnum          ! grid number
      integer :: id_invd       ! index of invalid super-droplets
      integer :: id_vd         ! index of valid super-droplets
      integer :: iwarn         ! warning flag
      integer :: n_invd        ! number of invalid super-droplets
      integer :: n_plus        ! number of adding super-droplets
      integer :: sdnum_valid   ! number of valid super-droplets in each grid
      integer :: sdnum_split   ! number of splitable valid super-droplets in each grid
      integer :: ip            ! index
      integer :: m, n, q, s, t ! index
     !---------------------------------------------------------------------

      call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
      call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

      ! Initialize
      gnum = ni_sdm * nj_sdm * knum_sdm
      freq_max = 1

      ! Add super-droplets in the grid with small number of super-droplets

      ! Sorting super-droplets

      !### valid droplets ###!

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
                    sort_id,sort_key,sort_freq,sort_tag,'valid')

      !### valid droplets that multiplicity is over 2 ###!

      allocate( sort_freq_valid(1:ni_sdm*nj_sdm*nk_sdm+1) )

      do n=1,ni_sdm*nj_sdm*nk_sdm+1
         sort_freq_valid(n) = sort_freq(n)
      end do

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
                    sort_id,sort_key,sort_freq,sort_tag,'multi')

      ! Initialize
      do n=1,20

         if( n<=sd_numasl ) then
            idx_nasl(n) = n
         else
            idx_nasl(n) = sd_numasl
         end if

      end do

      do n=1,gnum+2
         sort_tag0(n) = sort_tag(n) - 1  !! {1-xx} => {0-xx}
      end do

      ! Get the maximum number of super-droplets in each grid

      do m=1,gnum
         freq_max = max( freq_max, sort_freq(m) )
      end do

      ! Sorting the grids by the number of super-droplets in each grid
      allocate( fsort_tag(0:freq_max+1) )
      allocate( fsort_freq(0:freq_max)  )

      call gadg_count_sort( sort_freq(1:gnum), 0, freq_max,             &
                            fsort_freq, fsort_tag, fsort_id )

      fsort_tag(freq_max+1) = fsort_tag(freq_max) + fsort_freq(freq_max)

      ! Get random number using random number generator

      call gen_rand_array( sd_rng, sd_rand )

      ! Get random permutation layout of super-droplets in each grid

      call sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,            &
                       sort_tag0,fsort_tag,fsort_id,sd_rand,isd_perm(1:sd_num,1))

      ! Search invalid super-droplets
      n_invd = sort_freq(gnum+1)   !! number of invalid droplet

      ! Add super-droplets to grids with small number of super-droplets

      q = 1        !! counter
      iwarn = -1   !! warning flg

      do s=1,int(sdnum_lwr/2)

         do t=fsort_tag(1),fsort_tag(sdnum_lwr)-1

            m = fsort_id(t)

            !### get number of addable droplets in each grid ###!

            sdnum_valid = sort_freq_valid(m)
                                        !! number of valid droplets
            sdnum_split = sort_freq(m)
                                        !! number of valid droplets that
                                        !! multiplicity is over 2
                                        !! ( splitable droplets )

            n_plus = min( sdnum_split, max(sdnum_lwr-sdnum_valid,0) )

            if( s>n_plus ) cycle        !! cycle at more than number
                                        !! of addable droplets

            !### get random permutation ###!
            ! 'sdm_getperm' target grid within two or more droplets
            ! ip = 1    : sort_freq(m)=1
            ! ip = perm : sort_freq(m)>1
            ip = min(sort_freq(m),2) - 1                 !! 0 or 1

            ip = (1-ip) + isd_perm(sort_tag0(m) +s,1)*ip   !! forward
!SELECT     ip = (1-ip) + isd_perm(sort_tag(m+1)+s)*ip   !! backward

            !### get index of valid and invalid droplets ###!

            id_vd   = sort_id( sort_tag0(m) + ip )
            id_invd = sort_id( sort_tag0(gnum+1) + min(q,n_invd) )

            !### Add droplets ###!

            !== split multiplicity of valid droplets ==!
            !== and copy other status to invalid one ==!

            sd_n(id_invd)  = sd_n(id_vd)/2
            sd_n(id_vd)    = sd_n(id_vd) - sd_n(id_invd)

            sd_x(id_invd)  = sd_x(id_vd)
            sd_y(id_invd)  = sd_y(id_vd)
            sd_r(id_invd)  = sd_r(id_vd)
            sd_vz(id_invd) = sd_vz(id_vd)
            sd_rk(id_invd) = sd_rk(id_vd)

            do n=1,20
              sd_asl(id_invd,idx_nasl(n)) = sd_asl(id_vd,idx_nasl(n))
            end do

            !### count and check number of invaid droplets ###!

            q = q + 1

            if( q>n_invd ) then
               iwarn = 1
            end if

         end do

      end do

      if( iwarn==1 ) then
        if( IO_L ) then
         write(IO_FID_LOG,'(2a)')"  ### [SDM] : warning for stopping to add",    &
               " super-droplets to avoid buffer limit exceeded ###"
        else
         write(*,'(2a)')"  ### [SDM] : warning for stopping to add",    &
               " super-droplets to avoid buffer limit exceeded ###"
        endif
      end if

      ! Deallocate

      deallocate( sort_freq_valid )

      deallocate( fsort_tag  )
      deallocate( fsort_freq )

    return
  end subroutine sdm_sdadd
  !----------------------------------------------------------------------------
  subroutine sdm_aslsulf(sd_aslfmrate,sd_aslfmdt,             &
                         sd_num,sd_numasl,sd_n,sd_asl,sd_fmnc)

      use scale_process, only: &
         PRC_MPIstop
      ! Input variables
      real(RP), intent(in) :: sd_aslfmrate   ! formation rate of aerosol
      real(RP), intent(in) :: sd_aslfmdt     ! time interval to form aerosol
      integer,  intent(in) :: sd_num         ! number of super-droplets
      integer,  intent(in) :: sd_numasl      ! number of kind of chemical material contained as 
                                             ! water-soluble aerosol in super droplets
      real(RP), intent(in) :: sd_fmnc        ! number concentration of super-droplets at aerosol formation
      ! Input and output variables
      integer(DP), intent(inout) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl)  ! aerosol mass of super-droplets
      ! Work variables
      real(DP) :: n_amsul
      real(RP) :: rb_amsul
      real(RP) :: rmax_amsul
      real(RP) :: rmin_amsul
      real(RP) :: sgm_amsul   ! parameter for 3mode log-nomiral distribution
      real(RP) :: n0          ! number of real droplets per unit volume and per aerosol radius
      real(RP) :: dry_r       ! aerosol radius
      real(RP) :: delta       ! temporary
      real(RP) :: sdn_tmp     ! temporary
      integer :: iexced       ! temporary
      integer :: n, t         ! index
     !-------------------------------------------------------------------

      ! Set ammonium sulfate aerosol mass and muliplicity of super-droplets

      iexced = 1

      rb_amsul   = rb_amsul_zanten_m
      sgm_amsul  = sgm_amsul_zanten_m
      rmax_amsul = rmax_amsul_zanten_m
      rmin_amsul = rmin_amsul_zanten_m

      n_amsul = real(sd_aslfmrate,RP) * real(sd_aslfmdt,RP)

      do n=1,sd_num

         !### 1 : (NH4)2SO4 in modified vanZanten(2010) ###!
         delta = log(rmax_amsul) - log(rmin_amsul)
         dry_r = exp( log(rmin_amsul) + delta*sd_asl(n,1) )
         sd_asl(n,1) = F_THRD * ONE_PI                             &
                              * (dry_r*dry_r*dry_r) * rho_amsul

         !! n0(log(dry_r)) [m-3]
         !! 1st mode log-noraml distribution for ammonium sulfate
         delta = log(dry_r) - log(rb_amsul)
         delta = -(delta*delta)/(2.0_RP*log(sgm_amsul)*log(sgm_amsul))

         n0 = (n_amsul*exp(delta))/(sqrt(2.0_RP*ONE_PI)*log(sgm_amsul))

         !! continuous uniform distribution

         delta = 1.0_RP/(log(rmax_amsul)-log(rmin_amsul))

         !! muliplicity

         sdn_tmp = n0/(sd_fmnc*delta)

         if( sdn_tmp<(2.0_RP**63_RP) ) then
            sd_n(n) = nint( sdn_tmp, DP )
         else
            iexced = -1
         end if

      end do

      if( iexced<0 ) then
         call PRC_MPIstop
      end if

     ! Set other aerosol mass

      do t=2,sd_numasl
      do n=1,sd_num
         sd_asl(n,t) = 0.0_RP
      end do
      end do

    return
  end subroutine sdm_aslsulf
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_sdm_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0, &
       TEMP0  )
    use scale_grid_index
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    use m_sdm_coordtrans, only: &
       sdm_x2ri, sdm_y2rj
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QAD) ! effective radius
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)   ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! density                   [kg/m3]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)       ! temperature [K]

    real(RP) :: sum3(KA,IA,JA), sum2(KA,IA,JA)
    real(RP) :: drate             ! temporary
    integer  :: k, i, j
    integer  :: tlist_l            ! total list number for cloud
    integer  :: lcnt               ! counter
    integer  :: kl                 ! index
    integer  :: ku                 ! index
    integer  :: m, n               ! index
    integer  :: ilist_l(1:sdnum_s2c)
    !---------------------------------------------------------------------------
    !-- reset
    Re(:,:,:,:) = 0.0_RP
    sum3(:,:,:) = 0.0_RP
    sum2(:,:,:) = 0.0_RP

    call sdm_x2ri(sdnum_s2c,sdx_s2c,sdri_s2c,sdrk_s2c)
    call sdm_y2rj(sdnum_s2c,sdy_s2c,sdrj_s2c,sdrk_s2c)

    ! Effective radius is defined by r3m/r2m=1.5/lambda.
    ! dcoef is cancelled by the division of sum3/sum2
!    do j  = JS, JE
!    do i  = IS, IE
!    do k  = KS, KE
!       dcoef(k,i,j) = dxiv_sdm(i) * dyiv_sdm(j) / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
!    enddo
!    enddo
!    enddo

    tlist_l = 0

    lcnt = 0

    do n=1,sdnum_s2c

       if( sdrk_s2c(n)<VALID2INVALID ) cycle

       lcnt = lcnt + 1
       ilist_l(lcnt) = n

    end do

    tlist_l = lcnt

    if( tlist_l>0 ) then

       do m=1,tlist_l
          n = ilist_l(m)

          i = floor(sdri_s2c(n))+1
          j = floor(sdrj_s2c(n))+1
          k = floor(sdrk_s2c(n))+1

          sum3(k,i,j) = sum3(k,i,j)            &
               + sdr_s2c(n) * sdr_s2c(n) * sdr_s2c(n)   &
               * real(sdn_s2c(n),kind=RP)
          sum2(k,i,j) = sum2(k,i,j)            &
               + sdr_s2c(n) * sdr_s2c(n)             &
               * real(sdn_s2c(n),kind=RP)
       end do

       !=== correction at the verical boundary ===!

       do j=JS, JE
       do i=IS, IE

          !! at lower boundary

          kl    = floor(sdrkl_s2c(i,j))+1
          drate = real(kl,kind=RP) - sdrkl_s2c(i,j)
          if( drate<0.50_RP ) then
             sum3(kl,i,j) = 0.0_RP           !! <50% in share
             sum2(kl,i,j) = 0.0_RP           !! <50% in share
          else
             sum3(kl,i,j) = sum3(kl,i,j)/drate
             sum2(kl,i,j) = sum2(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sdrku_s2c(i,j))+1
          drate = sdrku_s2c(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             sum3(ku,i,j) = 0.0_RP           !! <50% in share
             sum2(ku,i,j) = 0.0_RP           !! <50% in share
          else
             sum3(ku,i,j) = sum3(ku,i,j)/drate
             sum2(ku,i,j) = sum2(ku,i,j)/drate
          end if

       end do
       end do

       !=== convert super-droplets to density of cloud-water. ===!

       do k=KS,KE
       do i=IS,IE
       do j=JS,JE
          if( sum2(k,i,j) == 0.0_RP ) then
           Re(k,i,j,I_mp_QC) = 0.0_RP
          else
           Re(k,i,j,I_mp_QC) = sum3(k,i,j)/sum2(k,i,j)
          endif
       end do
       end do
       end do

    end if

    return
  end subroutine ATMOS_PHY_MP_sdm_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_sdm_CloudFraction( &
       cldfrac, &
       QTRC     )
    use scale_precision
    use scale_grid_index
    use scale_tracer, only: &
       QAD => QA
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QAD)

    real(RP) :: qhydro
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do k  = KS, KE
    do i  = IS, IE
    do j  = JS, JE
       qhydro = 0.0_RP
       do iq = QQS, QQE 
          qhydro = qhydro + QTRC_sdm(k,i,j,iq)
       enddo
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-EPS)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_sdm_CloudFraction
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_sdm_MixingRatio( &
       Qe,    &
       QTRC0  )
    use scale_precision
    use scale_grid_index
    use scale_tracer, only: &
       QAD => QA , &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Qe   (KA,IA,JA,MP_QAD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]

    integer  :: ihydro
    !---------------------------------------------------------------------------

!!    do ihydro = 1, 2
!!       Qe(:,:,:,I_mp_QC) = QTRC_sdm(:,:,:,1)+QTRC_sdm(:,:,:,2)
!!    enddo
    Qe(:,:,:,I_mp_QC) = QTRC_sdm(:,:,:,I_QC)+QTRC_sdm(:,:,:,I_QR)

    return
  end subroutine ATMOS_PHY_MP_sdm_MixingRatio
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_restart_in
    implicit none

    integer :: sdnum_dum, sdnumasl_dum, sdfmnum_dum
    integer :: dp_dum, rp_dum, n, m, ierr
    real(DP) :: otime

    !### Get random generator seed ###!
    !! Random number generator has already been initialized in scale-les/src/preprocess/mod_mkinit.f90
    !! Be careful. If unit (=fid_random_i) is specified, filename is ignored and the object is initialized by the unit.
    call rng_load_state( rng_s2c, trim(RANDOM_IN_BASENAME))
!    call rng_load_state( rng_s2c, trim(RANDOM_IN_BASENAME), fid_random_i )

    open (fid_sd_i, file = trim(SD_IN_BASENAME), &! action = "read", &
          access = "sequential", status = "old", form = "unformatted", &
          iostat = ierr)

    if( ierr /= 0 ) then
      write(*,*) "sdm_restart_in", "read error"
      call PRC_MPIstop
    endif 
    
    !--- read time and precision
    read(fid_sd_i) otime, rp_dum, dp_dum, sdnum_dum, sdnumasl_dum, sdfmnum_dum
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file of Super Droplet  '
    if( IO_L ) write(IO_FID_LOG,*) 'SD. restart now time =  ', real(otime,kind=DP)
    if( rp_dum /= RP .or. dp_dum /= DP .or. &
        sdnum_dum /= sdnum_s2c .or. sdnumasl_dum /= sdnumasl_s2c .or. &
        sdfmnum_dum /= sdfmnum_s2c ) then
       write(*,*) 'xxx RP, DP, sdnum_s2c, sdnumasl_s2c, sdfmnum_s2c in'
       write(*,*) 'is different from those in Param file!  stop'
       write(*,*) 'RP(in restart file) = ', rp_dum
       write(*,*) 'DP(in restart file) = ', dp_dum
       write(*,*) 'sdnum(in restart file) = ', sdnum_dum
       write(*,*) 'RP(in restart file) = ', sdnumasl_dum
       write(*,*) 'RP(in restart file) = ', sdfmnum_dum
       call PRC_MPIstop
    endif

    !--- read S.D.
    read(fid_sd_i) (sdn_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdrk_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdx_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdy_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdz_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdr_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdu_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdv_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdvz_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) ((sdasl_s2c(n,m),n=1,sdnum_s2c),m=1,sdnumasl_s2c)
    read(fid_sd_i) (sdliqice_s2c(n),n=1,sdnum_s2c)
    if( sdm_cold ) then
       read(fid_sd_i) (sdice_s2c%re(n),n=1,sdnum_s2c)
       read(fid_sd_i) (sdice_s2c%rp(n),n=1,sdnum_s2c)
       read(fid_sd_i) (sdice_s2c%rho(n),n=1,sdnum_s2c)
       read(fid_sd_i) (sdice_s2c%tf(n),n=1,sdnum_s2c)
       read(fid_sd_i) (sdice_s2c%mrime(n),n=1,sdnum_s2c)
       read(fid_sd_i) (sdice_s2c%nmono(n),n=1,sdnum_s2c)
    end if
    !--- read formation S.D.
    close(fid_sd_i)
    if( IO_L ) write(IO_FID_LOG,*) '*** Closed restart file of Super Droplet  '

    return
  end subroutine ATMOS_PHY_MP_sdm_restart_in
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_restart_out(otime)
    implicit none

    real(DP), intent(in) :: otime
    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=17) :: fmt3="(3A)"
    character(len=H_LONG) :: ftmp, ftmp2
    character(len=H_LONG) :: basename_sd_out
    character(len=H_LONG) :: basename_random
    character(len=H_LONG) :: basename_time
    integer :: n, m, ierr


    !--- output restart file of Super Droplet
    write(basename_time,'(F15.3)') otime
    do n = 1, 15
       if( basename_time(n:n) == ' ' ) basename_time(n:n) = '0'
    enddo
    fid_sd_o = IO_get_available_fid()

    if( SD_OUT_BASENAME == '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found S.D. restart file name. Default used..'
       write(fmt2(14:14),'(I1)') 6
       write(fmt2(16:16),'(I1)') 6
       write(ftmp,fmt3) 'SD_output', '_', trim(basename_time)
!       write(SD_OUT_BASENAME,fmt2) trim(ftmp), 'pe',mype ! Perhaps a bug?
       write(basename_sd_out,fmt2) trim(ftmp),'pe',mype
    else
       write(fmt2(14:14),'(I1)') 6
       write(fmt2(16:16),'(I1)') 6
       write(ftmp,fmt3) trim(SD_OUT_BASENAME), '_', trim(basename_time)
       write(basename_sd_out,fmt2) trim(ftmp),'pe',mype
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file of Super Droplet  ', trim(basename_sd_out)

    open (fid_sd_o, file = trim(basename_sd_out), & !action = "write", &
          access = "sequential", status = "replace", form = "unformatted", &
          iostat = ierr)

    if( ierr /= 0 ) then
      write(*,*) "sdm_restart_out", "Write error"
      call PRC_MPIstop
    endif 

    !--- write time and precision
    write(fid_sd_o) otime, RP, DP, sdnum_s2c, sdnumasl_s2c, sdfmnum_s2c
    !--- write S.D.
    write(fid_sd_o) (sdn_s2c_restart(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdrk_s2c_restart(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdx_s2c_restart(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdy_s2c_restart(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdz_s2c_restart(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdr_s2c_restart(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdu_s2c_restart(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdv_s2c_restart(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdvz_s2c_restart(n),n=1,sdnum_s2c)
    write(fid_sd_o) ((sdasl_s2c_restart(n,m),n=1,sdnum_s2c),m=1,sdnumasl_s2c)
    write(fid_sd_o) (sdliqice_s2c_restart(n),n=1,sdnum_s2c)
    if( sdm_cold ) then
       write(fid_sd_o) (sdice_s2c_restart%re(n),n=1,sdnum_s2c)
       write(fid_sd_o) (sdice_s2c_restart%rp(n),n=1,sdnum_s2c)
       write(fid_sd_o) (sdice_s2c_restart%rho(n),n=1,sdnum_s2c)
       write(fid_sd_o) (sdice_s2c_restart%tf(n),n=1,sdnum_s2c)
       write(fid_sd_o) (sdice_s2c_restart%mrime(n),n=1,sdnum_s2c)
       write(fid_sd_o) (sdice_s2c_restart%nmono(n),n=1,sdnum_s2c)
    end if

    close(fid_sd_o)
    if( IO_L ) write(IO_FID_LOG,*) '*** Closed restart file of Super Droplet  '

    !--- output restart file of Random number
    if( RANDOM_OUT_BASENAME == '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found random number output file name. Default used..'
       write(ftmp,fmt3) 'random_number_output', '_', trim(basename_time)
       write(basename_random,fmt2) trim(ftmp),'pe',mype
    else
       write(ftmp,fmt3) trim(RANDOM_OUT_BASENAME), '_', trim(basename_time)
       write(basename_random,fmt2) trim(ftmp),'pe',mype
    endif

    fid_random_o = IO_get_available_fid()
    if( IO_L ) write(IO_FID_LOG,*) '*** Output random number for SDM ***', trim(basename_random)
!    call rng_save_state( rng_s2c, trim(basename_random))
    call rng_save_state( rng_s2c_restart, trim(basename_random))
!    call rng_save_state( rng_s2c, trim(basename_random), fid_random_o )

    return
  end subroutine ATMOS_PHY_MP_sdm_restart_out
end module scale_atmos_phy_mp_sdm
!-------------------------------------------------------------------------------
