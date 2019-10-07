!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM common variables
!!
!! @par Description
!!          Common variables used in the Super Droplet Method (SDM)
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
!! @li      2014-06-13 (S.Shima) [new] Separated common variables from scale_atmos_phy_mp_sdm.F90
!! @li      2014-07-14 (S.Shima) [rev] Removed unused variables: KMIN,KPLS,...,
!! @li      2014-07-18 (Y.Sato)  [add] add QTRC_sdm
!! @li      2014-12-12 (Y.Sato)  [mod] modify LatHet and DNS_RL as those used in SCALE Library
!! @li      2014-12-12 (Y.Sato)  [mod] modify characteristics of aeorosl from ammonium sulfate to ammonium bisulfate
!! @li      2014-12-19 (Y.Sato)  [mod] modify the location for defining LatHet and DNS_RL, and modify some typo
!! @li      2014-12-25 (Y.Sato)  [mod] modify LatHet to LH0
!! @li      2015-06-27 (S.Shima) [add] Add num_threads to store the threads number of auto parallelization on K/FX10
!! @li      2016-07-11 (S.Shima) [add] A structure for ice particles added
!! @li      2016-07-13 (S.Shima) [mod] Modified for cold SDM
!! @li      2016-07-15 (S.Shima) [mod] sdm_meltingfreezgin -> sdm_meltfreeze
!! @li      2016-07-16 (S.Shima) [add] sdm_rqi2qs added
!! @li      2016-07-16 (S.Shima) [mod] Dimension of sd_itmp1-3 modified
!! @li      2016-07-16 (S.Shima) [mod] Definition of prr_cres modified.
!! @li      2016-07-17 (S.Shima) [add] dt for sdm_meltfreeze iteration added
!! @li      2016-07-19 (S.Shima) [add] Constants for sublimation/deposition added
!! @li      2016-07-21 (S.Shima) [add] Working array for restart output moved from scale_atmos_phy_mp_sdm.F90
!! @li      2017-05-04 (S.Shima) [add] Table of habit (inherent growth ratio)
!! @li      2017-08-28 (S.Shima) [add] INAS density of Niemand et al. (2012) is implemented
!! @li      2017-10-26 (S.Shima) [add] var_k_coef
!! @li      2017-11-30 (S.Shima) [add] sd_dtmp2
!! @li      2018-02-28 (S.Shima) [add] sd_dtmp3,4,5,6
!! @li      2018-06-30 (S.Shima) [add] rime mass and number of monomers as SD attributes
!! @li      2019-01-09 (S.Shima) [add] sdm_fctr2multi
!! @li      2019-01-11 (S.Shima) [mod] default value of sdm_mvexchg
!!
!< 
!-------------------------------------------------------------------------------
module m_sdm_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio

  use scale_const, only: &
     ONE_PI => CONST_PI, &
     rrst   => CONST_R, &       ! Gas constant [J/(K*mol)]
     mass_air => CONST_Mdry, &  ! Molecular mass of air [g/mol]
     GasV_C => CONST_Rvap, &    ! Gas Constant of vapor [J/K/kg]
     LHV0    => CONST_LHV0, &   ! 2.5008E6_RP !< latent heat of vaporizaion at 0C [J/kg] 
     LHS0    => CONST_LHS0, &   ! 2.8342E6_RP !< latent heat of sublimation at 0C [J/kg]
     CPvap  => CONST_CPvap, &
     CL     => CONST_CL, &
     TEM00  => CONST_TEM00, &
     DWATR  => CONST_DWATR

  use rng_uniform_mt, only: c_rng_uniform_mt

  !-----------------------------------------------------------------------------
  implicit none
  public
  !-----------------------------------------------------------------------------
  !
  !++ Parameters & variables
  !
  !-----------------------------------------------------------------------------
  character(len=H_LONG), save :: SD_IN_BASENAME = ''
  character(len=H_LONG), save :: SD_OUT_BASENAME = ''
  character(len=H_LONG), save :: RANDOM_IN_BASENAME = ''
  character(len=H_LONG), save :: RANDOM_OUT_BASENAME = ''
  type(c_rng_uniform_mt), save :: rng_s2c
  integer, save :: fid_sd_i, fid_sd_o
  integer, save :: fid_random_i, fid_random_o
  logical, save :: sd_rest_flg_in = .false. ! restart flg of Super Droplet
  logical, save :: sd_first = .true. 
  integer, save :: num_threads ! number of threads of auto parallelization on K/FX10
  integer, parameter :: i2=2 ! Kind of 2byte integer variable
  !
  !++ Basic variable for SDM
  !
  !------------------------------------------------------------------------------
  integer(DP), allocatable, save :: sdn_s2c(:)   ! multipilicity
  real(RP), allocatable, save :: sdri_s2c(:)     ! index-i(real) of s.d.
  real(RP), allocatable, save :: sdrj_s2c(:)     ! index-j(real) of s.d.
  real(RP), allocatable, save :: sdrk_s2c(:)     ! index-k(real) of s.d.
  real(RP), allocatable, save :: sdx_s2c(:)      ! x-cordinate of s.d.
  real(RP), allocatable, save :: sdy_s2c(:)      ! y-cordinate of s.d.
  real(RP), allocatable, save :: sdz_s2c(:)      ! z-cordinate of s.d.
  real(RP), allocatable, save :: sdr_s2c(:)      ! y-cordinate of s.d.
  real(RP), allocatable, save :: sdu_s2c(:)      ! x-components velocity of s.d.
  real(RP), allocatable, save :: sdv_s2c(:)      ! y-components velocity of s.d.
  real(RP), allocatable, save :: sdvz_s2c(:)     ! z-components velocity of s.d. and terminal veloicty of s.d.
  real(RP), allocatable, save :: sdasl_s2c(:,:)  ! aeosol mass of s.d.
  real(RP), allocatable, save :: sdrkl_s2c(:,:)  ! index-k(real) at 'sdm_zlower'
  real(RP), allocatable, save :: sdrku_s2c(:,:)  ! index-k(real) at 'sdm_upper'
  integer(DP), allocatable, save :: sdn_tmp(:)   ! multiplicity of super-droplets
  real(RP), allocatable, save :: sdrk_tmp(:)     ! index-k(real) of super-droplets
  real(RP), allocatable, save :: sdx_tmp(:)      ! x-coordinate of super-droplets
  real(RP), allocatable, save :: sdy_tmp(:)      ! y-coordinate of super-droplets
  real(RP), allocatable, save :: sdz_tmp(:)      ! z-coordinate of super-droplets
  real(RP), allocatable, save :: sdr_tmp(:)      ! equivalent radius of super-droplets
  real(RP), allocatable, save :: sdu_tmp(:)      ! x-components velocity of super-droplets
  real(RP), allocatable, save :: sdv_tmp(:)      ! y-components velocity of super-droplets
  real(RP), allocatable, save :: sdvz_tmp(:)     ! terminal velocity of super-droplets zeta components of contravariant velocity
  real(RP), allocatable, save :: sdasl_tmp(:,:)  ! aerosol mass of super-droplets
  ! SDM for aerosol formation
  integer(DP), allocatable, save :: sdn_fm(:)    ! multiplicity of super-droplets
  real(RP), allocatable, save :: sdri_fm(:)      ! index-i(real) of super-droplets
  real(RP), allocatable, save :: sdrj_fm(:)      ! index-j(real) of super-droplets
  real(RP), allocatable, save :: sdrk_fm(:)      ! index-k(real) of super-droplets
  real(RP), allocatable, save :: sdx_fm(:)       ! x-coordinate of super-droplets
  real(RP), allocatable, save :: sdy_fm(:)       ! y-coordinate of super-droplets
  real(RP), allocatable, save :: sdz_fm(:)       ! z-coordinate of super-droplets
  real(RP), allocatable, save :: sdr_fm(:)       ! equivalent radius of super-droplets
  real(RP), allocatable, save :: sdvz_fm(:)      ! terminal velocity of super-droplets zeta components of contravariant velocity
  real(RP), allocatable, save :: sdasl_fm(:,:)   ! aerosol mass of super-droplets formed by gas-to-particle conversion
!!$  real(RP), allocatable, save :: rhod_crs(:,:,:) ! dry air density
!!$  real(RP), allocatable, save :: rhoc_sdm(:,:,:) ! density of cloud water
!!$  real(RP), allocatable, save :: rhor_sdm(:,:,:) ! density of rain water
!!$  real(RP), allocatable, save :: rhoa_sdm(:,:,:) ! density of aerosol
  real(RP), allocatable, save :: rand_s2c(:)
  integer, allocatable, save :: sortid_s2c(:)
  integer, allocatable, save :: sortfreq_s2c(:)
  integer, allocatable, save :: sortkey_s2c(:)
  integer, allocatable, save :: sorttag_s2c(:)
  real(DP), allocatable :: rbuf_r8(:,:,:)
                       ! reciving buffer for MPI (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf,mrime
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
  real(DP), allocatable :: sbuf_r8(:,:,:)
                       ! sending buffer for MPI (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf,mrime
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
  integer(DP), allocatable :: rbuf_i8(:,:,:)
                       ! reciving buffer for MPI (int8)
                       ! dim02 = 1 (multiplicity of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
  integer(DP), allocatable :: sbuf_i8(:,:,:)
                       ! sending buffer for MPI (int8)
                       ! dim02 = 1 (multiplicity of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
  integer(i2), allocatable :: rbuf_i2(:,:,:)
                       ! reciving buffer for MPI (int2)
                       ! dim02 = 1 (status of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
  integer(i2), allocatable :: sbuf_i2(:,:,:)
                       ! sending buffer for MPI (int2)
                       ! dim02 = 1 (status of super-droplets)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
  integer, allocatable :: rbuf_i4(:,:,:)
                       ! reciving buffer for MPI (int4)
                       ! dim02 = 1 number of monomers (ice)
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
  integer, allocatable :: sbuf_i4(:,:,:)
                       ! sending buffer for MPI (int4)
                       ! dim02 = 1 number of monomers
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
  integer, allocatable, save :: sdm_itmp1(:)
  integer, allocatable, save :: sdm_itmp2(:)
  integer, allocatable, save :: sdm_itmp3(:)
  integer, allocatable, save :: sd_itmp1(:)
  integer, allocatable, save :: sd_itmp2(:)
  integer, allocatable, save :: sd_itmp3(:)
  real(RP), allocatable, save :: sd_dtmp1(:)
  real(RP), allocatable, save :: sd_dtmp2(:)
  real(RP), allocatable, save :: sd_dtmp3(:)
  real(RP), allocatable, save :: sd_dtmp4(:)
  real(RP), allocatable, save :: sd_dtmp5(:)
  real(RP), allocatable, save :: sd_dtmp6(:)
  real(RP), allocatable, save :: QTRC_sdm(:,:,:,:)
  ! Module variables (SDM for ice)
  integer(kind=i2), allocatable :: sdliqice_s2c(:)
!  integer(kind=i2), allocatable :: sdliqice_tmp(:)
  integer(kind=i2), allocatable :: sdliqice_fm(:)
                          ! phase(liquid/ice) of super-droplets
                          ! 01 = all liquid, 10 = all ice
                          ! 11 = mixture of ice and liquid
  type sdicedef
     real(RP), allocatable :: re(:)
     ! equatorial radius of ice crystal
     real(RP), allocatable :: rp(:)
     ! polar radius of ice crystal
     real(RP), allocatable :: rho(:)
     ! density of ice crystal
!!$     real(RP), allocatable :: t(:)
!!$     ! temperature of ice crystal
     real(RP), allocatable :: tf(:)
     ! freezing temperature of ice crystal
     real(RP), allocatable :: mrime(:)
     ! rime mass of ice crystal [kg]
     integer, allocatable :: nmono(:)
     ! number of monomers (primary ice crystals) in an ice particle []
  end type sdicedef
  type(sdicedef), save :: sdice_s2c
!  type(sdicedef), save :: sdice_tmp
  type(sdicedef), save :: sdice_fm ! for fomation of super-droplets

  !! These working arrays are needed due to the scale restart output timimg. Could be removed in the future scale version.  
  integer(DP), allocatable, save :: sdn_s2c_restart(:)      ! multipilicity
  real(RP), allocatable, save    :: sdrk_s2c_restart(:)     ! index-k(real) of s.d.
  real(RP), allocatable, save    :: sdx_s2c_restart(:)      ! x-cordinate of s.d.
  real(RP), allocatable, save    :: sdy_s2c_restart(:)      ! y-cordinate of s.d.
  real(RP), allocatable, save    :: sdz_s2c_restart(:)      ! z-cordinate of s.d.
  real(RP), allocatable, save    :: sdr_s2c_restart(:)      ! y-cordinate of s.d.
  real(RP), allocatable, save    :: sdu_s2c_restart(:)      ! x-components velocity of s.d.
  real(RP), allocatable, save    :: sdv_s2c_restart(:)      ! y-components velocity of s.d.
  real(RP), allocatable, save    :: sdvz_s2c_restart(:)     ! z-components velocity of s.d. and terminal veloicty of s.d.
  real(RP), allocatable, save    :: sdasl_s2c_restart(:,:)  ! aeosol mass of s.d.
  type(c_rng_uniform_mt), save   :: rng_s2c_restart
  integer(i2), allocatable, save :: sdliqice_s2c_restart(:) ! phase(liquid/ice) of s.d.
  type(sdicedef), save :: sdice_s2c_restart                 ! ice attributes of s.d.

  !------------------------------------------------------------------------------
  !
  !++ Basic variable for SDM
  !
  !------------------------------------------------------------------------------
  real(RP) :: xmax_sdm, ymax_sdm      ! Distance in X, Y-direction
  real(RP) :: rkumax, rkumin, rklmax  ! Maxinum and Minuimum real number of sd_rku and sd_rkl
  real(RP) :: minzph                  ! Minimum of zph(i,j,2)
  integer :: knum_sdm                ! Maximum number of SDM grid in vertical
!  integer :: ni=IE-IS+1, nj=JE-JS+1, nz=KE-KS+1     ! Number of grid for each dirrection with halo
  integer :: stat
  integer :: dstw_sub, dste_sub, dsts_sub, dstn_sub
  integer :: srcw_sub, srce_sub, srcs_sub, srcn_sub
  !------------------------------------------------------------------------------
  !
  !++ Basic parameter for SDM
  !
  !------------------------------------------------------------------------------
!!$  real(RP), parameter :: LatHet = 2.453E+6_RP   ! Latent heat of water at 293K [J/kg]
!!$  real(RP), parameter :: LatHet = LHV0 - ( CPvap - CL )*TEM00   ! Latent heat of water at 0K [J/kg]
  real(RP), parameter :: LatHet = LHV0           ! Latent heat of water at 273.15K [J/kg] 
  real(RP), parameter :: Heat_C = 2.550E-2_RP   ! thermal conductivity at 293K, 100kPa [J/m s K]
  !--- Y.Sato
!!$  real(RP), parameter :: DNS_RL = 998.203_RP ! Density of liquid water at 293K [kg m-3]
  real(RP), parameter :: DNS_RL = DWATR         ! Density of liquid water [kg m-3]
  !--- Y.Sato
  real(RP), parameter :: Diff_C = 2.52E-5_RP    ! Diffusion Constant
  !--- Y.Sato
  real(RP), parameter :: LatGas = LatHet / GasV_C
  real(RP), parameter :: L_RL_K = LatHet * DNS_RL / Heat_C
  !--- Y.Sato
  real(RP), parameter :: RLRv_D = DNS_RL * GasV_C / Diff_C
  real(RP), parameter :: LatHet_S = LHS0   ! 2.8342E6_RP !< latent heat of sublimation at 0C [J/kg]
  real(RP), parameter :: LatGas_S = LatHet_S / GasV_C
  real(RP), parameter :: L_RL_K_S = LatHet_S * DNS_RL / Heat_C
!  real(RP), parameter :: exp_K = GasD_C / cp_C    ! exponent k used in adiabatic process
  real(RP), parameter :: m2micro = 1.0E+6_RP    ! Convert unit [m] -> [micron]
  real(RP), parameter :: micro2m = 1.0E-6_RP    ! Convert unit [micron] -> [m]
  real(RP), parameter :: TWO_PI = 2.0_RP * 3.14159265358979_RP
  real(RP), parameter :: O_THRD = 1.0_RP / 3.0_RP
  real(RP), parameter :: F_THRD = 4.0_RP / 3.0_RP
  real(RP), parameter :: O_SIX  = 1.0_RP / 6.0_RP
  real(RP), parameter :: T_EIGH = 3.0_RP / 8.0_RP
  ! thereshold value between valid super-droplet and invalid super-droplets
  real(RP), parameter :: VALID2INVALID = -999.0_RP 
  integer(DP), parameter :: VALID2INVALID_i8 = -990_DP
                     ! thereshold value between valid super-droplets
                     ! and invalid super-droplets
  real(RP), parameter :: PRECIPI = -999.9_RP  ! value indicated as super-droplets is precipitation
  real(RP), parameter :: PREC2INVALID = -999.99_RP ! thereshold value between precipitation and invalid super-droplets
  real(RP), parameter :: INVALID = -999.999_RP ! value indicated as invalid super-droplets
  integer(DP), parameter :: INVALID_i8 = -999_DP ! value indicated as invalid super-droplets
  integer(i2), parameter :: INVALID_i2 = -999_i2 ! value indicated as invalid super-droplets

  !------------------------------------------------------------------------------
  !
  !++ Parameter relating to microphysics
  !
  !------------------------------------------------------------------------------
  real(RP), parameter :: boltz = 1.38066E-23_RP  ! Boltzmann constant
  ! Molecular mass of sea salt (NaCl) contained as water-solble aerosol in S.D. [g]
  real(RP), parameter :: mass_nacl = 58.44277_RP 
!!  ! Molecular mass of ammonium sulfate contained as water-solble aerosol in S.D. [g]
!!  real(RP), parameter :: mass_amsul = 132.14_RP  
  ! Molecular mass of ammonium bisulfate contained as water-solble aerosol in S.D. [g]
  real(RP), parameter :: mass_amsul = 115.11_RP  
  ! Degree of ion dissociation of ammonium sulfate and sea salt contained 
  ! as water-soluble aerosol in super droplets
  real(RP), parameter :: ion_amsul = 2.0_RP, ion_nacl = 2.0_RP 
!!  ! Density of ammonium sulfate and sea salt [g/m3] contained as water-soluble aerosol in super droplets
!!  real(RP), parameter :: rho_amsul = 1.769E+6_RP, rho_nacl = 2.165E+6_RP
  ! Density of ammonium bisulfate and sea salt [g/m3] contained as water-soluble aerosol in super droplets
  real(RP), parameter :: rho_amsul = 1.78E+6_RP, rho_nacl = 2.165E+6_RP
  ! parameter for 3mode log-normal distribution of aerosol
  real(RP)  :: n1_amsul                                  !! [1/m3]
  real(RP)  :: n2_amsul
  real(RP)  :: n3_nacl
  real(RP)  :: rb1_amsul                                 !! [m]
  real(RP)  :: rb2_amsul
  real(RP)  :: rb3_nacl
  real(RP)  :: sgm1_amsul                                !! [-]
  real(RP)  :: sgm2_amsul
  real(RP)  :: sgm3_nacl
  real(RP)  :: rmax_amsul                                !! [m]
  real(RP)  :: rmax_nacl
  real(RP)  :: rmin_amsul
  real(RP)  :: rmin_nacl
  ! parameter for 3mode log-normal distribution under RICO-observation : Derksen(2009)
  real(RP), parameter :: n1_amsul_derksn   = 118.E+6_RP   !! [1/m3]
  real(RP), parameter :: n2_amsul_derksn   = 11.E+6_RP
  real(RP), parameter :: n3_nacl_derksn    = 0.1E+6_RP
  real(RP), parameter :: rb1_amsul_derksn  = 1.9E-08_RP   !! [m]
  real(RP), parameter :: rb2_amsul_derksn  = 5.6E-08_RP
  real(RP), parameter :: rb3_nacl_derksn   = 8.0E-07_RP
  real(RP), parameter :: sgm1_amsul_derksn = 3.3_RP       !! [-]
  real(RP), parameter :: sgm2_amsul_derksn = 1.6_RP
  real(RP), parameter :: sgm3_nacl_derksn  = 2.2_RP
  real(RP), parameter :: rmax_amsul_derksn = 5.0E-06_RP   !! [m]
  real(RP), parameter :: rmax_nacl_derksn  = 5.0E-06_RP
  real(RP), parameter :: rmin_amsul_derksn = 2.0E-09_RP
  real(RP), parameter :: rmin_nacl_derksn  = 1.0E-07_RP
  ! parameter for 3mode log-nomiral distribution under RICO-observation : vanZanten(2010)
  real(RP), parameter :: n1_amsul_zanten   = 90.E+6_RP    !! [1/m3]
  real(RP), parameter :: n2_amsul_zanten   = 15.E+6_RP
  real(RP), parameter :: rb1_amsul_zanten  = 3.0E-08_RP   !! [m]
  real(RP), parameter :: rb2_amsul_zanten  = 1.4E-07_RP
  real(RP), parameter :: sgm1_amsul_zanten = 1.28_RP      !! [-]
  real(RP), parameter :: sgm2_amsul_zanten = 1.75_RP
  real(RP), parameter :: rmax_amsul_zanten = 5.0E-06_RP   !! [m]
  real(RP), parameter :: rmin_amsul_zanten = 1.0E-08_RP
  ! parameter for 3mode log-nomiral distribution under RICO-observation : modified vanZanten(2010)
  real(RP), parameter :: rb_amsul_zanten_m   = 3.0E-08_RP   !! [m]
  real(RP), parameter :: sgm_amsul_zanten_m  = 1.28_RP      !! [-]
  real(RP), parameter :: rmax_amsul_zanten_m = 1.0E-07_RP   !! [m]
  real(RP), parameter :: rmin_amsul_zanten_m = 1.0E-08_RP
  real(RP), parameter :: Es_T_A = 2.53E+11_RP  ! A[kg/m s2] : Es = A * exp(-B/T)
  real(RP), parameter :: Es_T_B = 5.42E+3_RP   ! B[K] : Es = A * exp(-B/T)
  real(RP), parameter :: CurveF = 3.3E-7_RP    !  Curvature term Const. a=CoruveF/T [m K]
  real(RP), parameter :: ASL_FF = 4.3E-6_RP    !  Cloud Condensation term Const. b=ASL_FF*(Ion*M)/ms [m3]
  ! aerosol distribution parameters of insoluble component (only for cold SDM)
  !! For now, assume a monodisperse distribution of mineral dust
  real(RP) :: n_mdust   =   1.0E+6_RP  ! number densitiy of mineral dust [/m3]
  real(RP) :: mdust_dia =   1.0e-6_RP  ! diameter of a mineral dust [m]
  ! the INAS density parameters of Niemand et al. (2012)
  real(RP) :: a0_N12 = 1.0_RP    ! [/m2]
  real(RP) :: a1_N12 = 0.517_RP  ! [/degC]
  real(RP) :: a2_N12 = 8.934_RP  ! []
  real(RP) :: tfmax = -12.0_RP ! [degC]
  real(RP) :: tfmin = -36.0_RP ! [degC]
  real(RP) :: sdnumratio_soluble = 0.5_RP ! ratio of SDs used for soluble aerosol particles
                                          ! The rest is used for insoluble+soluble internally mixed aerosol particles
  real(RP) :: INIAsd_ratio = 0.05_RP ! the ratio of (IN inactive SD)/(insluble+soluble aerosol SD)
                                     ! Note that this can be chosen arbitrarily from [0,1).
                                     ! The multiplicty will be adjusted later to reproduce the real particle distribution.
  ! parameter to determine the projected area of ice particles: var_k = exp( -var_k_coef*sdi%rp(n)/sdi%re(n) )
  real(RP), parameter :: var_k_coef = 1.0_RP

  !------------------------------------------------------------------------------
  !
  !++ Terminal velocity
  !
  !------------------------------------------------------------------------------
  real(RP) :: vz_b(7), vz_c(6)
  data vz_b / -0.318657E+1_RP,  0.9926960_RP,    -0.153193E-2_RP, -0.987059E-3_RP, &
              -0.578878e-3_RP,  0.855176E-4_RP,  -0.327815E-5_RP /
  data vz_c / -0.500015E+1_RP,  0.523778E+1_RP,  -0.204914E+1_RP,              &
               0.47529400_RP,  -0.542819E-1_RP,   0.238449E-2_RP /
  !------------------------------------------------------------------------------
  !
  !++ Collision/Coalescence
  !
  !------------------------------------------------------------------------------
  real(RP) :: ratcol(21), r0col(15), ecoll(15,21)
  data r0col /   6.0_RP,  8.0_RP, 10.0_RP, 15.0_RP,  20.0_RP,  25.0_RP,  30.0_RP,  &
                40.0_RP, 50.0_RP, 60.0_RP, 70.0_RP, 100.0_RP, 150.0_RP, 200.0_RP, &
                300.0_RP /
  data ratcol / 0.00_RP, 0.050_RP, 0.10_RP, 0.150_RP, 0.20_RP, 0.250_RP,        &
                0.30_RP, 0.350_RP, 0.40_RP, 0.450_RP, 0.50_RP, 0.550_RP,        &
                0.60_RP, 0.650_RP, 0.70_RP, 0.750_RP, 0.80_RP, 0.850_RP,        &
                0.90_RP, 0.950_RP, 1.00_RP /
  data ecoll /                                                      &
      0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP  &
     ,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0030_RP,0.0030_RP,0.0030_RP,0.0040_RP,0.0050_RP  &
     ,0.0050_RP,0.0050_RP,0.0100_RP,0.1000_RP,0.0500_RP,0.2000_RP,0.5000_RP,0.7700_RP,0.8700_RP,0.9700_RP  &
     ,0.0070_RP,0.0070_RP,0.0070_RP,0.0080_RP,0.0090_RP,0.0100_RP,0.0100_RP,0.0700_RP,0.4000_RP,0.4300_RP  &
     ,0.5800_RP,0.7900_RP,0.9300_RP,0.9600_RP,1.0000_RP,0.0090_RP,0.0090_RP,0.0090_RP,0.0120_RP,0.0150_RP  &
     ,0.0100_RP,0.0200_RP,0.2800_RP,0.6000_RP,0.6400_RP,0.7500_RP,0.9100_RP,0.9700_RP,0.9800_RP,1.0000_RP  &
     ,0.0140_RP,0.0140_RP,0.0140_RP,0.0150_RP,0.0160_RP,0.0300_RP,0.0600_RP,0.5000_RP,0.7000_RP,0.7700_RP  &
     ,0.8400_RP,0.9500_RP,0.9700_RP,1.0000_RP,1.0000_RP,0.0170_RP,0.0170_RP,0.0170_RP,0.0200_RP,0.0220_RP  &
     ,0.0600_RP,0.1000_RP,0.6200_RP,0.7800_RP,0.8400_RP,0.8800_RP,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0300_RP,0.0300_RP,0.0240_RP,0.0220_RP,0.0320_RP,0.0620_RP,0.2000_RP,0.6800_RP,0.8300_RP,0.8700_RP  &
     ,0.9000_RP,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0250_RP,0.0250_RP,0.0250_RP,0.0360_RP,0.0430_RP  &
     ,0.1300_RP,0.2700_RP,0.7400_RP,0.8600_RP,0.8900_RP,0.9200_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0270_RP,0.0270_RP,0.0270_RP,0.0400_RP,0.0520_RP,0.2000_RP,0.4000_RP,0.7800_RP,0.8800_RP,0.9000_RP  &
     ,0.9400_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0300_RP,0.0300_RP,0.0300_RP,0.0470_RP,0.0640_RP  &
     ,0.2500_RP,0.5000_RP,0.8000_RP,0.9000_RP,0.9100_RP,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0400_RP,0.0400_RP,0.0330_RP,0.0370_RP,0.0680_RP,0.2400_RP,0.5500_RP,0.8000_RP,0.9000_RP,0.9100_RP  &
     ,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0350_RP,0.0350_RP,0.0350_RP,0.0550_RP,0.0790_RP  &
     ,0.2900_RP,0.5800_RP,0.8000_RP,0.9000_RP,0.9100_RP,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0370_RP,0.0370_RP,0.0370_RP,0.0620_RP,0.0820_RP,0.2900_RP,0.5900_RP,0.7800_RP,0.9000_RP,0.9100_RP  &
     ,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0370_RP,0.0370_RP,0.0370_RP,0.0600_RP,0.0800_RP  &
     ,0.2900_RP,0.5800_RP,0.7700_RP,0.8900_RP,0.9100_RP,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0370_RP,0.0370_RP,0.0370_RP,0.0410_RP,0.0750_RP,0.2500_RP,0.5400_RP,0.7600_RP,0.8800_RP,0.9200_RP  &
     ,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0370_RP,0.0370_RP,0.0370_RP,0.0520_RP,0.0670_RP  &
     ,0.2500_RP,0.5100_RP,0.7700_RP,0.8800_RP,0.9300_RP,0.9700_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0370_RP,0.0370_RP,0.0370_RP,0.0470_RP,0.0570_RP,0.2500_RP,0.4900_RP,0.7700_RP,0.8900_RP,0.9500_RP  &
     ,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0360_RP,0.0360_RP,0.0360_RP,0.0420_RP,0.0480_RP  &
     ,0.2300_RP,0.4700_RP,0.7800_RP,0.9200_RP,1.0000_RP,1.0200_RP,1.0200_RP,1.0200_RP,1.0200_RP,1.0200_RP  &
     ,0.0400_RP,0.0400_RP,0.0350_RP,0.0330_RP,0.0400_RP,0.1120_RP,0.4500_RP,0.7900_RP,1.0100_RP,1.0300_RP  &
     ,1.0400_RP,1.0400_RP,1.0400_RP,1.0400_RP,1.0400_RP,0.0330_RP,0.0330_RP,0.0330_RP,0.0330_RP,0.0330_RP  &
     ,0.1190_RP,0.4700_RP,0.9500_RP,1.3000_RP,1.7000_RP,2.3000_RP,2.3000_RP,2.3000_RP,2.3000_RP,2.3000_RP  &
     ,0.0270_RP,0.0270_RP,0.0270_RP,0.0270_RP,0.0270_RP,0.1250_RP,0.5200_RP,1.4000_RP,2.3000_RP,3.0000_RP  &
     ,4.0000_RP,4.0000_RP,4.0000_RP,4.0000_RP,4.0000_RP /
  !------------------------------------------------------------------------------
  !
  !++ For Namelist
  !
  !------------------------------------------------------------------------------
  real(RP), save :: sdm_dtcmph(5)              ! Time interval [sec] of condensation, coalescence, and droplet motion
                                               ! melting/freezing [coldsdm only], sublimation/deposition [coldsdm only]
  logical, save :: sdm_calvar(5)              ! flag for calculation of condensation, coalecscence, and droplet motion
  real(RP), save :: sdm_rdnc     = 1000.0_RP   ! Number of real droplet
  real(RP), save :: sdm_sdnmlvol = 1000.0_RP   ! Normal volume for number concentration of S.D. [m3]
  real(RP), save :: sdm_inisdnc  = 1000.0_RP   ! Number of super droplets per sdm_nmlvol at initial [1/m^3]
  real(RP), save :: sdm_zlower   = 5.0_RP      ! Lowest limitation of super droplet [m]
  real(RP), save :: sdm_zupper   = 2.E+3_RP    ! Hightest limitation of super droplet [m]
  integer, save :: sdm_extbuf   = 30          ! Rate to buffer size of super droplets for extra [%]
  real(RP), save :: sdm_rqc2qr   = 1.E-5_RP    ! Threshold of radius between qc and qr [m]
  real(RP), save :: sdm_rqi2qs   = 40.E-6_RP   ! Threshold of radius between qi and qs [m]
  integer, save :: sdm_colkrnl  = 2           ! Kernel type for coalescence process (0:Golovin, 1:Long, 2:Halls,3:No col_effi hydrodynamic)
  integer, save :: sdm_colbrwn  = 0           ! Flag of Brownian Coagulation and Scavenge process
  integer, save :: sdm_mvexchg  = 1           ! flag of exchange momentum betweeen super-droplets and fluid 0:No, 1:rhow only
  integer, save :: sdm_aslset   = 1.0_RP      ! Conrol flag to set the species and way of chemical material as water-soluble aerosol
  real(RP), save :: sdm_fctr2multi = 1.0_RP   ! Factor to increase the initial number density of soluble aerosol particles. This will not be applied to insoluble aersol.
  real(RP), save :: sdm_aslfmdt  = 0.1_RP      ! time interval [sec] of aerosol nucleation
  real(RP), save :: sdm_aslfmsdnc = 1000.0_RP  ! Number of S.D. per sdnmvol as aeroosl(s)
  real(RP), save :: sdm_aslfmrate = 0.0_RP ! Formation rate of super droplets as aerosol [1/(m^3 s)]
  real(RP), save :: sdm_aslmw(1:20)      ! User specified molecular mass of chemical material contained as water-soluble aerosol in S.D.
  real(RP), save :: sdm_aslion(1:20)     ! User specified ion of chemical material contained as water-soluble aerosol in S.D.
  real(RP), save :: sdm_aslrho(1:20)     ! User specified density of chemical material contained as water-soluble aerosol in S.D.
  real(RP), save :: sdm_nadjdt = 0.1_RP  ! Time interval of adjust S.D. number [s]
  integer, save :: sdm_nadjvar = 0      ! Control flag of adjustment super-droplet number in each grid
                                            ! 0:No adjust number of droplet
                                            ! 1:adjust number of droplet by adding
                                            ! 2:adjust number of droplet by removal 
                                            ! 3:adjust number of droplet by adding and removal
  integer, save :: sdm_dmpvar  = 000   ! Control flag to output Super Droplets
                                       ! 1st digit (10^0) corresponds to text output
                                       ! 2nd digit (10^1) corresponds to binary output
                                       ! 3rd digit (10^2) corresponds to binary output of large droplets 
                                       ! 0: off, 1: on, 2: output with sort data (only for binary output)
                                       ! currently only 001 is supported   
  real(RP), save :: sdm_dmpitva  = 0.e0 ! Time interval of text output [s]
  integer, save :: sdm_dmpnskip = 1     ! Base skip to store super droplets in text format
  real(RP), save :: sdm_dmpitvb  = 0.e0  ! Time interval of binary output of all droplets [s]
  real(RP), save :: sdm_dmpitvl  = 0.e0  ! Time interval of binary output of large droplets [s]
  real(RP), save :: sdm_dmpsdsiz = 0.e0  ! Threshold radius to store large super droplets in binary format [m]

  data sdm_dtcmph / 0.1_RP,0.1_RP,0.1_RP,0.1_RP,0.1_RP /
  data sdm_calvar / .false.,.false.,.false.,.false.,.false. /
  data sdm_aslmw  / 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP  /
  data sdm_aslion / 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP  /
  data sdm_aslrho / 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP  /
  real(RP), allocatable, save  :: prr_crs(:,:,:) ! prr_crs(IA,JA,1:6) precipitation rate [m/s] and accumlation [m]
                                               ! 1:rain rate, 2:rain accumulation
                                               ! 3:snow rate, 4:snow accumulation
                                               ! 5:precipitation rate, 6:precipitation accumulation
  real(RP), allocatable, save  :: zph_crs(:,:,:), dxiv_sdm(:), dyiv_sdm(:)!, dziv_sdm(:)
  real(RP), allocatable, save  :: dx_sdm(:), dy_sdm(:)!, dz_sdm(:)   ! Dx, Dy, Dz for SDM (normally they are equal to those of SCALE)
!  integer, parameter :: nqw = QQA
  integer, save :: sdfmnum_s2c
  real(RP), save:: sdininum_s2c
  integer, save :: sdnumasl_s2c, sdnum_s2c
  integer, save :: bufsiz1      ! Buffer size for MPI
  integer, save :: bufsiz2_r8 ! buffer size for MPI (real8)
  integer, save :: bufsiz2_i8 ! buffer size for MPI (int8)
  integer, save :: bufsiz2_i2 ! buffer size for MPI (int2)
  integer, save :: bufsiz2_i4 ! buffer size for MPI (int4)
  integer, save :: ni_s2c, nj_s2c, nk_s2c
  integer, save :: tag
  integer, parameter :: nomp = 1
  integer, save :: wbc=1, ebc=1, sbc=1, nbc=1  ! only periodic boundary is applied
  integer, save :: nsub   ! Number of sub domain in group domain
  integer, save :: nisub, njsub
  integer, save :: sthopt=1, trnopt=0
!  real(RP), save :: crs_dtmp1(KA,IA,JA), crs_dtmp2(KA,IA,JA)
!  real(RP), save :: crs_dtmp3(KA,IA,JA), crs_dtmp4(KA,IA,JA)
!  real(RP), save :: crs_dtmp5(KA,IA,JA), crs_dtmp6(KA,IA,JA)
!  real(RP), allocatable, save :: qwtr_crs(:,:,:,:)
!  real(RP), allocatable, save :: j31(:,:,:) ! z-x components of Jacobian
!  real(RP), allocatable, save :: j32(:,:,:) ! z-y components of Jacobian
!  real(RP), allocatable, save :: jcb(:,:,:) ! Jacobian at scalar points
!  real(RP), allocatable, save :: jcb8w(:,:,:) ! Jacobian at w points
!  real(RP), allocatable, save :: mf(:,:)  ! Map scale factors
!!$  integer,  allocatable, save :: KMIN1(:), IMIN1(:), JMIN1(:)
!!$  integer,  allocatable, save :: KPLS1(:), IPLS1(:), JPLS1(:)
  integer, save :: nclstp(0:5)

  logical, save :: docondensation = .true.
  logical, save :: doautoconversion = .true.
  logical, save :: doprecipitation  = .true.
  logical, save :: domovement = .true.
  logical, save :: donegative_fixer  = .true.  ! apply negative fixer?

  logical, save :: sdm_cold  = .false.  ! switch to turn on cold SDM
  logical, save :: domeltfreeze = .false.
  logical, save :: dosublimation = .false.

! Module variables (status of super-droplet)
  integer(kind=i2), parameter :: STAT_LIQ = 01_i2  !! all liquid
  integer(kind=i2), parameter :: STAT_ICE = 10_i2  !! all ice
  integer(kind=i2), parameter :: STAT_MIX = 11_i2  !! mixture of ice
                                                   !! and liquid

! Table of the inherent growth ratio (Chen and Lamb 1994)
! T(n) = (1.375-n)/4. Here, T is the temperature in [degreeC], and n is the index.
  real(RP), save :: tb_habit(1:121)
  data tb_habit / &
       1.000000e+00_RP, 9.799412e-01_RP, 9.600636e-01_RP, 9.399397e-01_RP, 9.200258e-01_RP, &
       9.001191e-01_RP, 8.800351e-01_RP, 8.570378e-01_RP, 8.340652e-01_RP, 8.109611e-01_RP, &
       7.830689e-01_RP, 7.550922e-01_RP, 7.030723e-01_RP, 5.370318e-01_RP, 4.677351e-01_RP, &
       5.248075e-01_RP, 6.309573e-01_RP, 8.128305e-01_RP, 1.096478e+00_RP, 1.479108e+00_RP, &
       1.905461e+00_RP, 2.089296e+00_RP, 2.290868e+00_RP, 2.398833e+00_RP, 2.454709e+00_RP, &
       2.426610e+00_RP, 2.371374e+00_RP, 2.290868e+00_RP, 2.137962e+00_RP, 1.995262e+00_RP, &
       1.862087e+00_RP, 1.737801e+00_RP, 1.621810e+00_RP, 1.513561e+00_RP, 1.396368e+00_RP, &
       1.288250e+00_RP, 1.188502e+00_RP, 1.096478e+00_RP, 1.000000e+00_RP, 9.225714e-01_RP, &
       8.511380e-01_RP, 7.852356e-01_RP, 7.244360e-01_RP, 6.683439e-01_RP, 6.165950e-01_RP, &
       5.754399e-01_RP, 5.370318e-01_RP, 5.011872e-01_RP, 4.677351e-01_RP, 4.365158e-01_RP, &
       4.073803e-01_RP, 3.801894e-01_RP, 3.548134e-01_RP, 3.311311e-01_RP, 3.162278e-01_RP, &
       3.019952e-01_RP, 2.917427e-01_RP, 2.851018e-01_RP, 2.818383e-01_RP, 2.786121e-01_RP, &
       2.754229e-01_RP, 2.786121e-01_RP, 2.818383e-01_RP, 2.851018e-01_RP, 2.917427e-01_RP, &
       2.985383e-01_RP, 3.090295e-01_RP, 3.198895e-01_RP, 3.311311e-01_RP, 3.467369e-01_RP, &
       3.672823e-01_RP, 3.935501e-01_RP, 4.265795e-01_RP, 4.570882e-01_RP, 4.897788e-01_RP, &
       5.248075e-01_RP, 5.623413e-01_RP, 6.095369e-01_RP, 6.606934e-01_RP, 7.161434e-01_RP, &
       7.852356e-01_RP, 8.609938e-01_RP, 9.549926e-01_RP, 1.047129e+00_RP, 1.148154e+00_RP, &
       1.258925e+00_RP, 1.380384e+00_RP, 1.496236e+00_RP, 1.603245e+00_RP, 1.698244e+00_RP, &
       1.778279e+00_RP, 1.840772e+00_RP, 1.883649e+00_RP, 1.905461e+00_RP, 1.905461e+00_RP, &
       1.883649e+00_RP, 1.862087e+00_RP, 1.840772e+00_RP, 1.798871e+00_RP, 1.737801e+00_RP, &
       1.698244e+00_RP, 1.640590e+00_RP, 1.584893e+00_RP, 1.549173e+00_RP, 1.513910e+00_RP, &
       1.476046e+00_RP, 1.452112e+00_RP, 1.428894e+00_RP, 1.412863e+00_RP, 1.393157e+00_RP, &
       1.376892e+00_RP, 1.361131e+00_RP, 1.348963e+00_RP, 1.336903e+00_RP, 1.327089e+00_RP, &
       1.317953e+00_RP, 1.308881e+00_RP, 1.302867e+00_RP, 1.293898e+00_RP, 1.287953e+00_RP, &
       1.279087e+00_RP /

end module m_sdm_common
