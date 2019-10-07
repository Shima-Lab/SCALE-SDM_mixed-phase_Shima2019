!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Coalescence subroutines for the SDM
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
!! @li      2017-02-15 (S.Shima) [new] copied from sdm_coalescence.f90
!! @li      2017-02-20 (S.Shima) [mod] riming and aggregation implemented
!! @li      2017-09-13 (S.Shima) [mod] subroutine sdm_outcome_drolet_coalescence() introduced
!! @li      2017-09-13 (S.Shima) [mod] subroutine sdm_outcome_riming() introduced
!! @li      2017-09-13 (S.Shima) [mod] subroutine sdm_outcome_aggregation() introduced
!! @li      2017-09-14 (S.Shima) [mod] a new aggregation outcome scheme by S.Shima (2017) is implemented
!! @li      2017-09-14 (S.Shima) [mod] subcycle is introduced for calculating riming/aggregation outcome
!! @li      2017-09-24 (S.Shima) [mod] subroutine of the coalescence kernel of Golovin (1963) created
!! @li      2017-09-24 (S.Shima) [mod] subroutine of the coalescence efficiency of Long (1974) created
!! @li      2017-09-24 (S.Shima) [fix] bug of the coalescence efficiency of Long (1974): r_large -> r_small
!! @li      2017-09-24 (S.Shima) [mod] subroutine of the coalescence efficiency of Davis (1972), Jonas (1972), and Hall (1980) created
!! @li      2017-09-24 (S.Shima) [mod] loop structure of kernel evaluation was changed
!! @li      2017-09-24 (S.Shima) [add] aggregation efficiency of Field et al. (2006) implemented
!! @li      2017-09-27 (S.Shima) [add] riming efficiency of Beard and Grover (1974), Erfani and Mitchell (2017) implemented
!! @li      2017-09-27 (S.Shima) [mod] sdm_outcome_aggregation() -> sdm_outcome_aggregation_shima2()
!! @li      2017-09-27 (S.Shima) [add] sdm_outcome_aggregation_Chen_Lamb_1994_mod
!! @li      2017-10-02 (S.Shima) [mod] riming density chagend from 300kg/m^3 to 100kg/m^3 following Jensen and Harrington (2015)
!! @li      2017-10-02 (S.Shima) [mod] aggregation is suppreseed if min(rho1,rho2)<10kg/m^3 (tentative)
!! @li      2017-10-26 (S.Shima) [mod] Following Connoly et al. 2012, projected area of particles is used to evaluate the geometric cross section of aggregation kernel
!! @li      2017-11-06 (S.Shima) [mod] When using Beard and Grover (1974) for droplets collecting ice, use Stokes number instead of mixed Froude number 
!! @li      2017-11-06 (S.Shima) [mod] Following Beard and Grover (1974), fix K0=0.21 for Nre>400.
!! @li      2017-11-13 (S.Shima) [mod] use volume weighted density to determine the contact angle of aggregation
!! @li      2017-11-30 (S.Shima) [mod] artificial limiter to suppress aggregation of very low density ice (mimicking breakup) 
!! @li      2017-11-30 (S.Shima) [fix] wrong eq to calculate sd_m2_tmp in sdm_outcome_riming()
!! @li      2017-12-03 (S.Shima) [add] sdm_outcome_aggregation_shima4()
!! @li      2018-02-21 (S.Shima) [mod] sdm_outcome_riming() to allow multiple riming
!! @li      2018-02-21 (S.Shima) [add] sdm_outcome_aggregation_shima4_multi()
!! @li      2018-04-03 (S.Shima) [add] riming density formula of Heymsfield and Pflaum (1985)
!! @li      2018-06-21 (S.Shima) [add] sdm_outcome_aggregation_shima5_multi
!! @li      2018-06-30 (S.Shima) [add] calculation of rime mass and number of monomers
!! @li      2018-07-03 (S.Shima) [fix] calculation of rime mass and number of monomers
!! @li      2019-01-10 (S.Shima) [mod] mimick tumbling if the ice is quasi spherical (Jensen and Harrington, 2015)
!! @li      2019-01-10 (S.Shima) [mod] new geometric cross-sectional formula for collision-riming kernel (Shima et al., 2019)
!! @li      2019-01-12 (S.Shima) [mod] dry air density -> moist air density
!
!<
!-------------------------------------------------------------------------------
module m_sdm_coalescence_cold
  use scale_precision

  implicit none
  private
  public :: sdm_coales_cold

contains
  subroutine sdm_coales_cold(sdm_colkrnl,sdm_colbrwn,          &
                        sdm_aslset,sdm_aslrho,                 &
                        sdm_dtcol,            &
                        pres_scale, t_scale,qv_scale,rhom_scale,&
                        zph_crs,                &
                        ni_sdm,nj_sdm,nk_sdm,sd_num,sd_numasl, &
                        sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,&
                        sdi,                                   &
                        sort_id,sort_key,sort_freq,sort_tag,   &
                        sd_rng,sd_rand,                        &
                        sort_tag0,fsort_id,icp,sd_perm,c_rate  )
    use gadg_algorithm, only: &
         gadg_count_sort
    use rng_uniform_mt, only: &
         c_rng_uniform_mt, gen_rand_array
    use scale_grid_index, only: &
         IA,JA,KA
    use scale_const, only:  &
         grav_mks => CONST_GRAV, &
         t0 => CONST_TEM00, &    ! 0 degC in [K]
         rhoi_mks => CONST_DICE  ! density of ice [kg/m^3]
    use m_sdm_common, only: &
         VALID2INVALID,INVALID,knum_sdm, &
         sdicedef,STAT_LIQ,STAT_ICE,     &
         rho_amsul,rho_nacl,ONE_PI,dxiv_sdm,dyiv_sdm,F_THRD,O_THRD,rrst,boltz,mass_air,i2,var_k_coef
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_idutil, only: &
         sdm_sort, sdm_getperm
    !  Input variables
    integer, intent(in) :: sdm_colkrnl   ! Kernel type for coalescence process
    integer, intent(in) :: sdm_colbrwn   ! Control flag of Brownian Coagulation and Scavenging process
    integer, intent(in) :: sdm_aslset    ! Control flag to set species and way of chemical material as water-soluble aerosol
    real(RP),intent(in) :: sdm_aslrho(20)! User specified density of chemical material contained as water-soluble aerosol in super droplets
    real(RP), intent(in) :: sdm_dtcol   ! tims step of {stochastic coalescence} process
    real(RP), intent(in) :: pres_scale(KA,IA,JA)  ! Pressure
    real(RP), intent(in) :: t_scale(KA,IA,JA)     ! Temperature
    real(RP), intent(in) :: qv_scale(KA,IA,JA)    ! water-vapor mixing ratio[kg/kg]
    real(RP), intent(in) :: rhom_scale(KA,IA,JA)  ! moist air density [kg/m^3]
    real(RP), intent(in) :: zph_crs(KA,IA,JA) ! z physical coordinate
    integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
    integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
    integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
    integer, intent(in) :: sd_num  ! number of super-droplets
    integer, intent(in) :: sd_numasl ! Number of kind of chemical material contained as water-soluble aerosol in super droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    ! Input and output variables
    type(c_rng_uniform_mt), intent(inout) :: sd_rng ! random number generator
    real(RP),intent(inout) :: sd_rand(1:sd_num) ! random numbers
    integer, intent(inout) :: sort_id(1:sd_num) ! super-droplets sorted by SD-grids
    integer, intent(inout) :: sort_key(1:sd_num) ! sort key
    integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1) ! number of super-droplets in each SD-grid
    integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2) ! accumulated number of super-droplets in each SD-grid
    integer(DP), intent(inout) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(i2), intent(inout) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    type(sdicedef), intent(inout) :: sdi   ! ice phase super-droplets
    ! Output variables
    integer, intent(out) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)  ! = sort_tag(n) - 1
    integer, intent(out) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
    integer, intent(out) :: icp(1:sd_num) ! index of coalescence pair
    integer, intent(out) :: sd_perm(1:sd_num) ! random permutations
    real(RP), intent(out) :: c_rate(1:sd_num) ! coalescence probability
    ! Internal shared variables
    real(RP) :: sd_aslrho(1:22) ! Density of chemical material contained as water-soluble aerosol in super droplets
    integer  :: sd_ncol ! how many times coalescence occurs
    integer :: freq_max ! get the maximum number of super-droplets in each grid
    integer :: hfreq_max ! hfreq_max / 2
    integer :: ipremium  ! premium coef. for coalescence
    ! Work variables
    real(RP) :: dmask(1:22) ! mask for vactorization
    real(RP) :: sd_asl1(1:sd_numasl)
    real(RP) :: sd_asl2(1:sd_numasl) ! aerosol mass of super-droplets with large/small multiplicity
    real(RP) :: sd_cc1  ! slip correction of super-droplets
    real(RP) :: sd_cc2  ! with large/small multiplicity
    real(RP) :: sd_dia1 ! diameter of super-droplets
    real(RP) :: sd_dia2 ! with large/small multiplicity
    real(RP) :: sd_lmd1 ! mean free path of super-droplets
    real(RP) :: sd_lmd2 ! with large/small multiplicity
    real(RP) :: sd_m1   ! mass of super-droplets
    real(RP) :: sd_m2   ! with large/small multiplicity
    real(RP) :: sd_v1   ! volume of super-droplets
    real(RP) :: sd_v2   ! with large/small multiplicity
    real(RP) :: sd_r1   ! radius of super-droplets
    real(RP) :: sd_r2   ! with large/small multiplicity
    real(RP) :: sd_rk1  ! index[k/real] of super-droplets with large multiplicity
    real(RP) :: sd_rw1  ! radius of water parts in super-droplets
    real(RP) :: sd_rw2  ! with large/small multiplicity
    real(RP) :: sd_tmasl1
    real(RP) :: sd_tmasl2 ! total mass of aerosol part in super droplets with large/small multiplicity
    real(RP) :: sd_tvasl1 ! total volume of aerosol part in super
    real(RP) :: sd_tvasl2 ! droplets with large/small multiplicity
    real(RP) :: sd_c1   ! temporary
    real(RP) :: sd_c2
    real(RP) :: sd_d1
    real(RP) :: sd_d2
    real(RP) :: sd_g1
    real(RP) :: sd_g2
    real(RP) :: sd_re1
    real(RP) :: sd_re2
    real(RP) :: sd_rp1
    real(RP) :: sd_rp2
    real(RP) :: sd_rho1
    real(RP) :: sd_rho2
    real(RP) :: sd_tf1
    real(RP) :: sd_tf2
    integer  :: sd_nmono1
    integer  :: sd_nmono2
    real(RP) :: sd_mrime1
    real(RP) :: sd_mrime2
    real(RP) :: sd_vz1
    real(RP) :: sd_vz2
    real(RP) :: sd_eqr1  ! equivalent radius
    real(RP) :: sd_eqr2  ! equivalent radius
    real(RP) :: sd_maxD1 ! maximum dimension of the particle [m]
    real(RP) :: sd_maxD2 ! maximum dimension of the particle [m]
    real(RP) :: sd_phi1  ! sd_rp/sd_re: aspect ratio
    real(RP) :: sd_phi2  ! sd_rp/sd_re: aspect ratio
    real(RP) :: sd_area_ce1 ! circumscribed ellipse area of the particle projected to the flow direction
    real(RP) :: sd_area_ce2 ! circumscribed ellipse area of the particle projected to the flow direction
    real(RP) :: sd_area1 ! A: area of the particle projected to the flow direction [m^2]
    real(RP) :: sd_area2 ! A: area of the particle projected to the flow direction [m^2]
    real(RP) :: var_k1   ! k:=exp(-var_k_coef*c/a). Needed to evaluate A.
    real(RP) :: var_k2   ! k:=exp(-var_k_coef*c/a). Needed to evaluate A.
    real(RP) :: lmd_crs ! air mean free path
    real(RP) :: p_crs   ! pressure
    real(RP) :: t_crs   ! temperarure
    real(RP) :: vis_crs ! dynamic viscosity
    real(RP) :: sumdia  ! sum of variables of a pair of droplets
    real(RP) :: sumd
    real(RP) :: sumc
    real(RP) :: sumg
    real(RP) :: sumr
    real(RP) :: k12     ! brownian coagulation coefficient
    real(RP) :: dvz     ! difference in terminal velocity of a pair of super-droplets
    real(RP) :: dtmp    ! temporary
    real(RP) :: frac    ! fraction parts
    real(RP) :: ivvol(KA,IA,JA)   ! inverse of a grid volume
    real(RP) :: tdeg    ! temperature in degree
    real(RP) :: sd_rtmp ! temporary
    real(RP) :: size_ratio ! particle size ratio
    real(RP) :: nre        ! Reynolds number
    real(RP) :: n_mixFr    ! mixed Froude number
    real(RP) :: n_kn       ! Knudsen number
    real(RP) :: n_st       ! Stokes number
    real(RP) :: C_sc       ! Cunningham slip correction factor
    real(RP) :: sd_dvisc   ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
    real(RP) :: sd_t       ! temperature [K] at the location of super-droplets
    real(RP) :: sd_rd      ! moist air density [kg/m3] at the location of super-droplets
    real(RP) :: sd_p       ! pressure of the grid contained the SD [Pa]
    real(RP) :: sd_qv      ! water-vapor of the grid contained the SD [kg/kg]
    real(RP) :: ri,rj,rk
    real(RP) :: rim_eff_plate, rim_eff_column
    real(RP) :: nre_column

    integer(DP) :: sd_nmax  ! maximum multiplicity
    integer(DP) :: sd_n1    ! multiplicity of super-droplets with large multiplicity
    integer(DP) :: sd_n2    ! multiplicity of super-droplets  with small multiplicity
    integer(i2) :: sd_li1, sd_li2

    integer, allocatable :: fsort_tag(:) ! buffer for sorting
    integer, allocatable :: fsort_freq(:) ! buffer for sorting

    integer :: idx_nasl(1:22)  ! index for vactorization

    integer :: gnum          ! grid number

    integer :: in1, in2, in3, in4 ! index
    integer :: i, j, k, m, n, s   ! index
    integer :: tc, tp             ! index
    integer :: i_col

    integer :: sort_tag0m
    integer :: sort_freqm
    integer :: icptc, icptp

    !
    real(RP), parameter :: lambda_mfp = 6.62E-8 ! mean free path of air [m]
    real(RP), parameter :: para_Asc = 2.51E0    ! parameters of Cunningham slip correction factor
    real(RP), parameter :: para_Bsc = 0.800E0
    real(RP), parameter :: para_Csc = 0.550E0 

    integer :: skip_flag ! flag to save whether the aggregation is canceled or not

    !--------------------------------------------------------------------
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_coales_cold",1,1)
#endif
    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

    ! Initialize
    ! ni_sdm=IE-IS+1, nj_sdm=JE-JS+1, knum_sdm=floor(rkumax)+1)-KS+1
    gnum = ni_sdm * nj_sdm * knum_sdm

    freq_max = 1

    !### aerosol type ###!
    ! why we need sd_aslrho? Isn't it same to sdm_aslrho?? Check later
    do n=1,22
       sd_aslrho(n) = 1.0_RP
    end do

    if( abs(mod(sdm_aslset,10))==1 ) then

       !### numasl=1 @ init+rest : (NH4)2SO4 ###!

       sd_aslrho(1) = real(rho_amsul,kind=RP)

    else if( abs(mod(sdm_aslset,10))==2 ) then

       if( abs(sdm_aslset)==2 ) then

          !### numasl=1 @ init : NaCl ###!

          sd_aslrho(1) = real(rho_nacl,kind=RP)

       else if( abs(sdm_aslset)==12 ) then

          !### numasl=2 @ init : NaCl, rest : (NH4)2SO4 ###!

          sd_aslrho(1) = real(rho_amsul,kind=RP)
          sd_aslrho(2) = real(rho_nacl,kind=RP)

       end if
    else if( abs(mod(sdm_aslset,10))==3 ) then

       !### numasl>=2 @ init+rest : (NH4)2SO4, NaCl, ... ###!

       sd_aslrho(1) = real(rho_amsul,kind=RP)
       sd_aslrho(2) = real(rho_nacl,kind=RP)

       !         do n=1,20
       !            call getrname( id_sdm_aslrho + (n-1), sd_aslrho(n+2) )
       !         end do
       do n=1,20
          sd_aslrho(n+2) = sdm_aslrho(n)
       end do

    end if

    ! Sorting super-droplets.
    
    call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
         sort_id,sort_key,sort_freq,sort_tag,'valid')

    ! Initialize
    
    do n=1,22

       if( n<=sd_numasl ) then
          idx_nasl(n) = n
          dmask(n) = 1.0_RP
       else
          idx_nasl(n) = sd_numasl
          dmask(n) = 0.0_RP
       end if

    end do

    do n=1,gnum+2
       sort_tag0(n) = sort_tag(n) - 1  !! {1-xx} => {0-xx}
    end do

    do n=1,sd_num
       c_rate(n) = 0.0_RP
    end do

    ! Get the maximum number of super-droplets in each grid
    
    do m=1,gnum
       freq_max = max( freq_max, sort_freq(m) )
    end do

    hfreq_max = int(freq_max/2)

    ! Sorting the grids by the number of super-droplets in each grid

    allocate( fsort_tag(0:freq_max+1) )
    allocate( fsort_freq(0:freq_max)  )

    call gadg_count_sort( sort_freq(1:gnum), 0, freq_max, &
         fsort_freq, fsort_tag, fsort_id )

    fsort_tag(freq_max+1) = fsort_tag(freq_max) + fsort_freq(freq_max)

    ! Get random number using random number generator
    
    call gen_rand_array( sd_rng, sd_rand )

    ! Get random permutation layout of super-droplets in each grid

    call sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,            &
         &                 sort_tag0,fsort_tag,fsort_id,sd_rand,sd_perm)

    ! Get random number using random number generator

    call gen_rand_array( sd_rng, sd_rand )

    ! Select a pair of super-droples in random permutation layout

!!$      do n=1,hfreq_max
!!$
!!$         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)

!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle
       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2

          in1 = 2 * n
          in2 = in1 - 1

          !### select pair of super-droplets in each grid ###!

!!$            in3 = sd_perm( sort_tag0(m) + in1 )
!!$            in4 = sd_perm( sort_tag0(m) + in2 )
          in3 = sd_perm( sort_tag0m + in1 )
          in4 = sd_perm( sort_tag0m + in2 )

          !### set the random index ###!
!!$            tc = sort_tag0(m) + n
          tc = sort_tag0m + n
!!$            tp = tc + int(sort_freq(m)/2)
          tp = tc + int(sort_freqm/2)

!!$            icp(tc) = sort_id( sort_tag0(m) + in3 )
!!$            icp(tp) = sort_id( sort_tag0(m) + in4 )
          icp(tc) = sort_id( sort_tag0m + in3 )
          icp(tp) = sort_id( sort_tag0m + in4 )

       end do

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Get effective coalescence/riming/aggregation kernel for "Gravitational Settling" 
    
!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle

       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2

          tc = sort_tag0m + n
          tp = tc + sort_freqm/2

          icptc = icp(tc)
          icptp = icp(tp)

          sd_li1 = sd_liqice( icptc ) 
          sd_li2 = sd_liqice( icptp ) 

          if(     (sd_li1 .eq. STAT_LIQ) .and.(sd_li2 .eq. STAT_LIQ)) then
             ! coalescence of two droplets
             if( sdm_colkrnl==0 ) then
                !### Golovin's kernel (m^3/s) ###!

                sd_r1 = sd_r(icptc)
                sd_r2 = sd_r(icptp)
                call coalescence_kernel_golovin(c_rate(tc),sd_r1,sd_r2)

             else if( sdm_colkrnl==1 ) then
                !### Long's kernel (m^3/s)  ###!

                sd_r1 = sd_r(icptc)
                sd_r2 = sd_r(icptp)
                call coalescence_efficiency_long(c_rate(tc),sd_r1,sd_r2)

                sumr = sd_r1 + sd_r2
                dvz = abs( sd_vz(icptc) - sd_vz(icptp) )
                c_rate(tc) = c_rate(tc) * ONE_PI * (sumr*sumr) * dvz

             else if( sdm_colkrnl==2 ) then
                !### Davis,Jonas,and Hall's kernel (m^3/s) ###!

                sd_r1 = sd_r(icptc)
                sd_r2 = sd_r(icptp)
                call coalescence_efficiency_davis_jonas_hall(c_rate(tc),sd_r1,sd_r2)

                sumr = sd_r1 + sd_r2
                dvz = abs( sd_vz(icptc) - sd_vz(icptp) )
                c_rate(tc) = c_rate(tc) * ONE_PI * (sumr*sumr) * dvz

             else if( sdm_colkrnl==3 ) then
                !### no coalescence effeciency hydrodynamic kernel (m^3/s) ###!

                sd_r1 = sd_r(icptc)
                sd_r2 = sd_r(icptp)
                sumr = sd_r1 + sd_r2
                dvz = abs( sd_vz(icptc) - sd_vz(icptp) )
                c_rate(tc) = c_rate(tc) * ONE_PI * (sumr*sumr) * dvz

             end if
          else if((sd_li1 .eq. STAT_LIQ) .and.(sd_li2 .eq. STAT_ICE)) then
             ! riming of water droplet and ice particle

!!$             !### temporary (Davis,Jonas,Hall)
!!$             sd_r1 = sd_r(icptc)
!!$             sd_r2 = max(sdi%re(icptp),sdi%rp(icptp))
!!$             call coalescence_efficiency_davis_jonas_hall(c_rate(tc),sd_r1,sd_r2)
!!$             dvz = abs( sd_vz(icptc) - sd_vz(icptp) )
!!$             c_rate(tc) = c_rate(tc) * ONE_PI * (sd_r1+sdi%re(icptp)) &
!!$                                              * (sd_r1+sd_r2) &
!!$                                              * dvz

             sd_vz1 = sd_vz(icptc)
             sd_r1  = sd_r(icptc)

             sd_vz2 = sd_vz(icptp)
             sd_re2 = sdi%re(icptp)
             sd_rp2 = sdi%rp(icptp)
             sd_rho2 = sdi%rho(icptp)
             sd_eqr2= (sd_re2**2 * sd_rp2)**(O_THRD)
             sd_maxD2 = 2.0d0 * max(sd_re2,sd_rp2) ! maximum dimension of the particle [m]

             !! coversion to center grid index
             ri = sd_ri(icptc)
             rj = sd_rj(icptc)
             rk = sd_rk(icptc)
             i = floor(ri)+1
             j = floor(rj)+1
             k = floor(rk)+1

             ! moist air density [kg/m^3] at the location of super-droplets
             sd_rd = rhom_scale(k,i,j)

             ! temperature[K] at the location of super-droplets
             sd_t = t_scale(k,i,j)

             ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
             !== (Pruppacher & Klett,1997) ==!
             tdeg = sd_t - t0     !! [K] => [degC]
             if( tdeg>=0.0d0 ) then
                sd_dvisc = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-5_RP
             else
                sd_dvisc = ( 1.718_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-5_RP
             end if

             if(sd_vz1>sd_vz2) then ! collector is droplet

                !### Beard ang Grover (1974)
                size_ratio = sd_eqr2/sd_r1
                !###### Raynolds number
                nre = sd_rd * (2.0d0*sd_r1) * sd_vz1 / sd_dvisc
                !###### Knudsen number
                n_kn = lambda_mfp/(2.0*sd_eqr2)
                !###### Cunningham slip correction factor
                C_sc = 1.0d0 + n_kn * (para_Asc + para_Bsc * exp(-para_Csc/n_kn) )
                !###### Stokes number
                n_st = size_ratio**2 * sd_rho2 * nre * C_sc / (9.0*sd_rd)

                call riming_efficiency_Beard_Grover_1974(c_rate(tc),size_ratio,nre,n_st)

                dvz = abs( sd_vz1 - sd_vz2 )
                c_rate(tc) = c_rate(tc) * ONE_PI * (sd_r1+sd_re2) &
                                                 * (sd_r1+max(sd_re2,sd_rp2)) &
                                                 * dvz
             else ! collector is ice

                !### Beard ang Grover (1974)
                size_ratio = sd_r1/sd_eqr2
                !###### Raynolds number
                nre = sd_rd * sd_maxD2 * sd_vz2 / sd_dvisc
                !###### mixed Froude number
                n_mixFr = (sd_vz2 - sd_vz1) * sd_vz1 / ((sd_maxD2/2.0d0) * grav_mks)
                call riming_efficiency_Beard_Grover_1974(c_rate(tc),size_ratio,nre,n_mixFr)

                !### Erfani and Mitchell (2017)
                sd_phi2 = sd_rp2/sd_re2
                if(sd_phi2<1.0d0)then ! oblate
                   call riming_efficiency_Erfani_Mitchell_2017_plate(rim_eff_plate,nre,n_mixFr)
                   c_rate(tc) = sd_phi2*c_rate(tc) +(1.0d0-sd_phi2)*rim_eff_plate
                else ! prolate
                   !###### Raynolds number based on column width
                   nre_column = sd_rd * (2.0d0*sd_re2) * sd_vz2 / sd_dvisc
                   call riming_efficiency_Erfani_Mitchell_2017_column(rim_eff_column,nre_column,n_mixFr)
                   c_rate(tc) = (1.0d0/sd_phi2)*c_rate(tc) +(1.0d0-(1.0d0/sd_phi2))*rim_eff_column
                end if

                dvz = abs( sd_vz1 - sd_vz2 )

                ! circumscribed ellipse area of the particle projected to the flow direction
                sd_area_ce2 = ONE_PI * sd_re2 * max(sd_re2,sd_rp2)

                ! area of the particle projected to the flow direction
                var_k2 = exp( -var_k_coef*sd_phi2)
                sd_area2 = sd_area_ce2 * (sd_rho2/rhoi_mks)**var_k2

                ! collision-riming kernel: E_rime*A_g*|vj-vk|
                ! geometric cross-sectional area A_g is calculated by subtracting the indentation (Shima et al. 2019)
                c_rate(tc) = c_rate(tc) * (  (ONE_PI * (sd_r1+sd_re2) &
                                                     * (sd_r1+max(sd_re2,sd_rp2))) &
                                           - (sd_area_ce2 - sd_area2) &
                                          ) * dvz

             end if

          else if((sd_li1 .eq. STAT_ICE) .and.(sd_li2 .eq. STAT_LIQ)) then
             ! riming of ice particle and water droplet

!!$             !### temporary (Davis,Jonas,Hall)
!!$             sd_r1 = max(sdi%re(icptc),sdi%rp(icptc))
!!$             sd_r2 = sd_r(icptp)
!!$             call coalescence_efficiency_davis_jonas_hall(c_rate(tc),sd_r1,sd_r2)
!!$             dvz = abs( sd_vz(icptc) - sd_vz(icptp) )
!!$             c_rate(tc) = c_rate(tc) * ONE_PI * (sdi%re(icptc)+sd_r2) &
!!$                                              * (sd_r1+sd_r2) &
!!$                                              * dvz

             sd_vz2 = sd_vz(icptp)
             sd_r2  = sd_r(icptp)

             sd_vz1 = sd_vz(icptc)
             sd_re1 = sdi%re(icptc)
             sd_rp1 = sdi%rp(icptc)
             sd_rho1= sdi%rho(icptc)
             sd_eqr1= (sd_re1**2 * sd_rp1)**(O_THRD)
             sd_maxD1 = 2.0d0 * max(sd_re1,sd_rp1) ! maximum dimension of the particle [m]

             !! coversion to center grid index
             ri = sd_ri(icptc)
             rj = sd_rj(icptc)
             rk = sd_rk(icptc)
             i = floor(ri)+1
             j = floor(rj)+1
             k = floor(rk)+1

             ! moist air density [kg/m^3] at the location of super-droplets
             sd_rd = rhom_scale(k,i,j)

             ! temperature[K] at the location of super-droplets
             sd_t = t_scale(k,i,j)

             ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
             !== (Pruppacher & Klett,1997) ==!
             tdeg = sd_t - t0     !! [K] => [degC]
             if( tdeg>=0.0d0 ) then
                sd_dvisc = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-5_RP
             else
                sd_dvisc = ( 1.718_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-5_RP
             end if

             if(sd_vz2>sd_vz1) then ! collector is droplet

                !### Beard ang Grover (1974)
                size_ratio = sd_eqr1/sd_r2
                !###### Raynolds number
                nre = sd_rd * (2.0d0*sd_r2) * sd_vz2 / sd_dvisc
                !###### Knudsen number
                n_kn = lambda_mfp/(2.0*sd_eqr1)
                !###### Cunningham slip correction factor
                C_sc = 1.0d0 + n_kn * (para_Asc + para_Bsc * exp(-para_Csc/n_kn) )
                !###### Stokes number
                n_st = size_ratio**2 * sd_rho1 * nre * C_sc / (9.0*sd_rd)

                call riming_efficiency_Beard_Grover_1974(c_rate(tc),size_ratio,nre,n_st)

                dvz = abs( sd_vz2 - sd_vz1 )
                c_rate(tc) = c_rate(tc) * ONE_PI * (sd_r2+sd_re1) &
                                                 * (sd_r2+max(sd_re1,sd_rp1)) &
                                                 * dvz
             else ! collector is ice

                !### Beard ang Grover (1974)
                size_ratio = sd_r2/sd_eqr1
                !###### Raynolds number
                nre = sd_rd * sd_maxD1 * sd_vz1 / sd_dvisc
                !###### mixed Froude number
                n_mixFr = (sd_vz1 - sd_vz2) * sd_vz2 / ((sd_maxD1/2.0d0) * grav_mks)
                call riming_efficiency_Beard_Grover_1974(c_rate(tc),size_ratio,nre,n_mixFr)

                !### Erfani and Mitchell (2017)
                sd_phi1 = sd_rp1/sd_re1
                if(sd_phi1<1.0d0)then ! oblate
                   call riming_efficiency_Erfani_Mitchell_2017_plate(rim_eff_plate,nre,n_mixFr)
                   c_rate(tc) = sd_phi1*c_rate(tc) +(1.0d0-sd_phi1)*rim_eff_plate
                else ! prolate
                   !###### Raynolds number based on column width
                   nre_column = sd_rd * (2.0d0*sd_re1) * sd_vz1 / sd_dvisc
                   call riming_efficiency_Erfani_Mitchell_2017_column(rim_eff_column,nre_column,n_mixFr)
                   c_rate(tc) = (1.0d0/sd_phi1)*c_rate(tc) +(1.0d0-(1.0d0/sd_phi1))*rim_eff_column
                end if

                dvz = abs( sd_vz2 - sd_vz1 )

                ! circumscribed ellipse area of the particle projected to the flow direction
                sd_area_ce1 = ONE_PI * sd_re1 * max(sd_re1,sd_rp1)

                ! area of the particle projected to the flow direction
                var_k1 = exp( -var_k_coef*sd_phi1)
                sd_area1 = sd_area_ce1 * (sd_rho1/rhoi_mks)**var_k1

                ! collision-riming kernel: E_rime*A_g*|vj-vk|
                ! geometric cross-sectional area A_g is calculated by subtracting the indentation (Shima et al. 2019)
                c_rate(tc) = c_rate(tc) * (  (ONE_PI * (sd_r2+sd_re1) &
                                                     * (sd_r2+max(sd_re1,sd_rp1))) &
                                           - (sd_area_ce1 - sd_area1) &
                                          ) * dvz

             end if

          else if((sd_li1 .eq. STAT_ICE) .and.(sd_li2 .eq. STAT_ICE)) then
             ! aggregation of two ice particles

             sd_vz1 = sd_vz(icptc)
             sd_re1 = sdi%re(icptc)
             sd_rp1 = sdi%rp(icptc)
             sd_rho1 = sdi%rho(icptc)
             sd_maxD1 = 2.0d0 * max(sd_re1,sd_rp1) ! maximum dimension of the particle [m]

             sd_vz2 = sd_vz(icptp)
             sd_re2 = sdi%re(icptp)
             sd_rp2 = sdi%rp(icptp)
             sd_rho2 = sdi%rho(icptp)
             sd_maxD2 = 2.0d0 * max(sd_re2,sd_rp2) ! maximum dimension of the particle [m]

             !### Field et al, (2006), Morrison and Grabowski (2010)
             c_rate(tc) = 0.1_RP
!             c_rate(tc) = 0.0_RP

!!$             !! artificial limiter to suppress aggregation of very low density ice
!!$             if(min(sd_rho1,sd_rho2)<10.0d0)then
!!$                c_rate(tc) = 0.0_RP
!!$             end if

!             c_rate(tc) = c_rate(tc) * (1.0-tanh(-log(min(sdi%rho(icptc),sdi%rho(icptp))/100.0)))/2.0

             dvz = abs( sd_vz2 - sd_vz1 )

             ! A: area of the particle projected to the flow direction [m^2]
             var_k1 = exp( -var_k_coef*sd_rp1/sd_re1 )
             sd_area1 = ONE_PI * sd_re1 * (sd_maxD1/2.0d0) * (sd_rho1/rhoi_mks)**var_k1

             var_k2 = exp( -var_k_coef*sd_rp2/sd_re2 )
             sd_area2 = ONE_PI * sd_re2 * (sd_maxD2/2.0d0) * (sd_rho2/rhoi_mks)**var_k2

             c_rate(tc) = c_rate(tc) * (sqrt(sd_area1)+sqrt(sd_area2))**2 &
                                     * dvz

!!$             c_rate(tc) = c_rate(tc) * ONE_PI * (sd_re1+max(sd_re2,sd_rp2)) &
!!$                                                 * (sd_re2+max(sd_re1,sd_rp1)) &
!!$                                              * dvz

          end if

       end do

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Convert kernel to probability for "Gravitational Settling" 

!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle
       
       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2
             
          tc = sort_tag0m + n
          tp = tc + sort_freqm/2

          icptc = icp(tc)
          icptp = icp(tp)

          !! Get location of Super-Droplets
          ! index in center grid
          i = floor(sd_ri(icptc))+1
          j = floor(sd_rj(icptc))+1
          k = floor(sd_rk(icptc))+1

          ivvol(k,i,j) = 1.0_RP * dxiv_sdm(i) * dyiv_sdm(j) &
               / (zph_crs(k,i,j)-zph_crs(k-1,i,j))

          c_rate(tc) = c_rate(tc) * real(sdm_dtcol,kind=RP) * ivvol(k,i,j)

       end do
    end do

    ! Get effective collision for "Brownian Coagulation and Scavenging
    ! (Seinfeld & Pandis,2006)" 
    ! This is mechanisim is effective for droplets less than micrometer-size ( below 1um )
    if( sdm_colbrwn>0 ) then
!!$        do n=1,hfreq_max
!!$
!!$          do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)
!!$
!!$            tc = sort_tag0(m) + n
!!$            tp = tc + sort_freq(m)/2
!!$
!OCL NORECURRENCE
       do m=1,gnum
          if( sort_freq(m) <= 1 ) cycle

          sort_tag0m = sort_tag0(m)
          sort_freqm = sort_freq(m)
          
          do n=1,sort_freqm/2

             tc = sort_tag0m + n
             tp = tc + sort_freqm/2

             icptc = icp(tc)
             icptp = icp(tp)

             !### information of super-droplets ###

             !! radius of water parts in droplets [m]

             sd_li1 = sd_liqice( icptc )
             sd_li2 = sd_liqice( icptp )

             if(sd_li1 .eq. STAT_LIQ) then
                sd_rw1 = sd_r(icptc)
             else
                sd_rw1 = (sdi%re(icptc)**2 * sdi%rp(icptc))**(1.0d0/3.0d0) 
                ! equivalent radius
                ! just an idea. no rigorous justification to use this.
             end if

             if(sd_li2 .eq. STAT_LIQ) then
                sd_rw2 = sd_r(icptp)
             else
                sd_rw2 = (sdi%re(icptp)**2 * sdi%rp(icptp))**(1.0d0/3.0d0) 
                ! equivalent radius
                ! just an idea. no rigorous justification to use this.
             end if

             !! mass and volume of aerosol parts in droplets

             sd_tmasl1 = 0.0_RP
             sd_tmasl2 = 0.0_RP
             sd_tvasl1 = 0.0_RP
             sd_tvasl2 = 0.0_RP

             do k=1,22

                s = idx_nasl(k)

                sd_tmasl1 = sd_tmasl1 + sd_asl(icptc,s) * dmask(k)
                sd_tmasl2 = sd_tmasl2 + sd_asl(icptp,s) * dmask(k)

                sd_tvasl1 = sd_tvasl1                                    &
                     + sd_asl(icptc,s)/sd_aslrho(s) * dmask(k)
                sd_tvasl2 = sd_tvasl2                                    &
                     + sd_asl(icptp,s)/sd_aslrho(s) * dmask(k)

             end do

             sd_tmasl1 = sd_tmasl1  * 1.E-3_RP    !! [g]=>[kg]
             sd_tmasl2 = sd_tmasl2  * 1.E-3_RP

             !! diameter and mass and volume of droplets

             dtmp = ONE_PI * F_THRD

             sd_v1 = sd_tvasl1 + dtmp * (sd_rw1*sd_rw1*sd_rw1)
             sd_v2 = sd_tvasl2 + dtmp * (sd_rw2*sd_rw2*sd_rw2)

             sd_m1 = sd_tmasl1 + dtmp * (sd_rw1*sd_rw1*sd_rw1) * sdi%rho(icptc)
             sd_m2 = sd_tmasl2 + dtmp * (sd_rw2*sd_rw2*sd_rw2) * sdi%rho(icptp)

             sd_dia1 = (6.0_RP*sd_v1/ONE_PI)**O_THRD
             sd_dia2 = (6.0_RP*sd_v2/ONE_PI)**O_THRD

             !### location of super-droplets ###

             ! index in center grid
             i = floor(sd_ri(icptc))+1
             j = floor(sd_rj(icptc))+1
             k = floor(sd_rk(icptc))+1

             ivvol(k,i,j) = 1.0_RP * dxiv_sdm(i) * dyiv_sdm(j) &
                  / (zph_crs(k,i,j)-zph_crs(k-1,i,j))

             p_crs = pres_scale(k,i,j) !! [Pa]
             t_crs = t_scale(k,i,j)    !! [K]

             !### dynamic viscosity [Pa*s]  ###!
             !### (Pruppacher & Klett,1997) ###!
             tdeg = t_crs - t0     !! [K] => [degC]
             if( tdeg>=0.0_RP ) then
                vis_crs = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-5_RP
             else
                vis_crs = ( 1.7180_RP + 4.9E-3_RP*tdeg                        &
                     -1.2E-5_RP*tdeg*tdeg ) * 1.E-5_RP
             end if

             !### air mean free path [m] ###!
             dtmp = dsqrt(8.0_RP*mass_air*1.E-3_RP/(ONE_PI*rrst*t_crs))
             lmd_crs = (2.0_RP*vis_crs)/(p_crs*dtmp)

             !### slip correction of droplets [-]  ###!
             dtmp   = 1.2570_RP + 0.40_RP * exp(-0.550_RP*sd_dia1/lmd_crs)
             sd_cc1 = 1.0_RP + (2.0_RP*lmd_crs*dtmp)/(sd_dia1)
             dtmp   = 1.2570_RP + 0.40_RP * exp(-0.550_RP*sd_dia2/lmd_crs)
             sd_cc2 = 1.0_RP + (2.0_RP*lmd_crs*dtmp)/(sd_dia2)

             !### diffusion term [m*m/s] ###!
             dtmp = (boltz*t_crs)/(3.0_RP*ONE_PI*vis_crs)
             sd_d1 = dtmp * (sd_cc1/sd_dia1)
             sd_d2 = dtmp * (sd_cc2/sd_dia2)

             !### velocity term [m/s] ###!
             dtmp = (8.0_RP*boltz*t_crs)/ONE_PI
             sd_c1 = dsqrt(dtmp/sd_m1)
             sd_c2 = dsqrt(dtmp/sd_m2)

             !### mean free path of droplets [m] ###!
             dtmp = 8.0_RP/ONE_PI
             sd_lmd1 = dtmp * (sd_d1/sd_c1)
             sd_lmd2 = dtmp * (sd_d2/sd_c2)

             !### length term [m] ###!
             dtmp = (sd_dia1+sd_lmd1)*(sd_dia1+sd_lmd1)*(sd_dia1+sd_lmd1)&
                  - (sd_dia1*sd_dia1+sd_lmd1*sd_lmd1)**1.50_RP
             sd_g1 = dtmp/(3.0_RP*sd_dia1*sd_lmd1) - sd_dia1
             dtmp = (sd_dia2+sd_lmd2)*(sd_dia2+sd_lmd2)*(sd_dia2+sd_lmd2)&
                  - (sd_dia2*sd_dia2+sd_lmd2*sd_lmd2)**1.50_RP
             sd_g2 = dtmp/(3.0_RP*sd_dia2*sd_lmd2) - sd_dia2

             !### Brownian Coagulation Coefficient K12 [m3/s] ###!
             sumdia = sd_dia1 + sd_dia2
             sumd   = sd_d1   + sd_d2
             sumc   = dsqrt( sd_c1*sd_c1 + sd_c2*sd_c2 )
             sumg   = dsqrt( sd_g1*sd_g1 + sd_g2*sd_g2 )

             dtmp = sumdia/(sumdia+2.0_RP*sumg) + (8.0_RP*sumd)/(sumdia*sumc)
             k12 = 2.0_RP*ONE_PI * sumdia*sumd/dtmp

             !### add effective collision [-] ###!
             c_rate(tc) = c_rate(tc)                                     &
                  + k12 * real(sdm_dtcol,kind=RP) * ivvol(k,i,j)
          end do

       end do

    end if

    ! Get total effective collision of droplets
    
!!$      do n=1,hfreq_max
!!$         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)
!!$
!!$            tc = sort_tag0(m) + n
!!$            tp = tc + sort_freq(m)/2
!!$
!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle

       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2

          tc = sort_tag0m + n
          tp = tc + sort_freqm/2

          icptc = icp(tc)
          icptp = icp(tp)

          sd_nmax  = max( sd_n(icptc), sd_n(icptp) )
          ! maximum multiplicity
          ipremium = sort_freqm - 1 + iand(sort_freqm,1)
          ! IAND(sort_freq(i),1) => even:0, odd:1

          c_rate(tc) = c_rate(tc) * real( sd_nmax*ipremium, kind=RP )
       end do
    end do

    ! Stochastic coalescence process.
!!$      do n=1,hfreq_max
!!$         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
!!$
!!$            m = fsort_id(t)
!!$
!!$            tc = sort_tag0(m) + n
!!$            tp = tc + sort_freq(m)/2

!OCL NORECURRENCE
    do m=1,gnum
       if( sort_freq(m) <= 1 ) cycle

       sort_tag0m = sort_tag0(m)
       sort_freqm = sort_freq(m)

       do n=1,sort_freqm/2

          tc = sort_tag0m + n
          tp = tc + sort_freqm/2

          icptc = icp(tc)
          icptp = icp(tp)

          !### set coalescence count ###!
          sd_ncol = int( c_rate(tc), kind=RP )
          frac = c_rate(tc) - real( sd_ncol, kind=RP )

          ! judge coalescence by random number and fractional part

          if( sd_rand(tc) < frac ) then
             sd_ncol =  sd_ncol + 1
          end if

          if( sd_ncol<=0 ) cycle  !! no coalesecense

          !### coalescence procudure ###!

          if( sd_n(icptc) > sd_n(icptp) ) then

             ! SD with larger multiplicity
             sd_n1  = sd_n( icptc )
             sd_r1  = sd_r( icptc )
             sd_rk1 = sd_rk( icptc )
             sd_li1 = sd_liqice( icptc )
             sd_re1 = sdi%re( icptc )  
             sd_rp1 = sdi%rp( icptc )
             sd_rho1 = sdi%rho( icptc )
             sd_tf1 = sdi%tf( icptc )
             sd_mrime1 = sdi%mrime( icptc )
             sd_nmono1 = sdi%nmono( icptc )

             do k=1,22
                s = idx_nasl(k)
                sd_asl1(s) = sd_asl( icptc,s )
             end do

             sd_vz1 = sd_vz(icptc)

             ! SD with smaller multiplicity
             sd_n2  = sd_n( icptp )
             sd_r2  = sd_r( icptp )
             sd_li2 = sd_liqice( icptp )
             sd_re2 = sdi%re( icptp )  
             sd_rp2 = sdi%rp( icptp )
             sd_rho2 = sdi%rho( icptp )
             sd_tf2 = sdi%tf( icptp )
             sd_mrime2 = sdi%mrime( icptp )
             sd_nmono2 = sdi%nmono( icptp )

             do k=1,22
                s = idx_nasl(k)
                sd_asl2(s) = sd_asl( icptp,s )
             end do

             sd_vz2 = sd_vz(icptp)

          else

             ! SD with larger multiplicity
             sd_n1  = sd_n( icptp )
             sd_r1  = sd_r( icptp )
             sd_rk1 = sd_rk( icptp )
             sd_li1 = sd_liqice( icptp )
             sd_re1 = sdi%re( icptp )  
             sd_rp1 = sdi%rp( icptp )
             sd_rho1 = sdi%rho( icptp )
             sd_tf1 = sdi%tf( icptp )
             sd_mrime1 = sdi%mrime( icptp )
             sd_nmono1 = sdi%nmono( icptp )

             do k=1,22
                s = idx_nasl(k)
                sd_asl1(s) = sd_asl( icptp,s )
             end do

             sd_vz1 = sd_vz(icptp)

             ! SD with smaller multiplicity
             sd_n2  = sd_n( icptc )
             sd_r2  = sd_r( icptc )
             sd_li2 = sd_liqice( icptc )
             sd_re2 = sdi%re( icptc )  
             sd_rp2 = sdi%rp( icptc )
             sd_rho2 = sdi%rho( icptc )
             sd_tf2 = sdi%tf( icptc )
             sd_mrime2 = sdi%mrime( icptc )
             sd_nmono2 = sdi%nmono( icptc )

             do k=1,22
                s = idx_nasl(k)
                sd_asl2(s) = sd_asl( icptc,s )
             end do

             sd_vz2 = sd_vz(icptc)

          end if

          !! calculate the attributes after coalescence
          sd_ncol = min( sd_ncol, int(sd_n1/sd_n2,kind=RP) )

          if( (sd_li1 .eq. STAT_LIQ) .and. (sd_li2 .eq. STAT_LIQ) ) then ! droplet-droplet coalescence

             call sdm_outcome_drolet_coalescence(sd_ncol, sd_r1, sd_r2)

          else if( ((sd_li1 .eq. STAT_ICE) .and. (sd_li2 .eq. STAT_LIQ)) .or. &
                 & ((sd_li1 .eq. STAT_LIQ) .and. (sd_li2 .eq. STAT_ICE))        ) then ! riming

             !! coversion to center grid index
             ri = sd_ri(icptc)
             rj = sd_rj(icptc)
             rk = sd_rk(icptc)
             i = floor(ri)+1
             j = floor(rj)+1
             k = floor(rk)+1

             ! temperature[K] at the location of super-droplets
             sd_t = t_scale(k,i,j)

             ! moist air density [kg/m^3] at the location of super-droplets
             sd_rd = rhom_scale(k,i,j)

             ! pressure of the grid contained the SD [Pa]
             sd_p  = pres_scale(k,i,j)
             
             ! water-vapor of the grid contained the SD [kg/kg]
             sd_qv = qv_scale(k,i,j)

             call sdm_outcome_riming(sd_ncol,sd_r1,sd_r2,    & 
                                           & sd_li1,sd_li2,  &
                                           & sd_re1,sd_re2,  &
                                           & sd_rp1,sd_rp2,  &
                                           & sd_rho1,sd_rho2,&
                                           & sd_vz1,sd_vz2,&
                                           & sd_mrime1,sd_mrime2,&
                                           & sd_nmono1,sd_nmono2,&
                                           & sd_t,sd_rd,sd_p,sd_qv)

!!$             do i_col=1,sd_ncol
!!$                call sdm_outcome_riming(     sd_r1,sd_r2,    & 
!!$                                           & sd_li1,sd_li2,  &
!!$                                           & sd_re1,sd_re2,  &
!!$                                           & sd_rp1,sd_rp2,  &
!!$                                           & sd_rho1,sd_rho2)
!!$             end do
!!$
          else if( (sd_li1 .eq. STAT_ICE) .and. (sd_li2 .eq. STAT_ICE) ) then ! aggregation

             call sdm_outcome_aggregation_shima5_multi(sd_ncol,   &
                                                     & sd_re1,sd_re2,  &
                                                     & sd_rp1,sd_rp2,  &
                                                     & sd_rho1,sd_rho2)

             ! update the number of monomers
             sd_nmono2 = sd_nmono2 + sd_ncol*sd_nmono1
             ! update the rime mass
             sd_mrime2 = sd_mrime2 + real(sd_ncol,kind=RP)*sd_mrime1

!!$             skip_flag=0
!!$             do i_col=1,sd_ncol
!!$
!!$                call sdm_outcome_aggregation_shima4(                &
!!$                                           & sd_re1,sd_re2,  &
!!$                                           & sd_rp1,sd_rp2,  &
!!$                                           & sd_rho1,sd_rho2,skip_flag)

!!$                call sdm_outcome_aggregation_Chen_Lamb_1994_mod(                &
!!$                                           & sd_re1,sd_re2,  &
!!$                                           & sd_rp1,sd_rp2,  &
!!$                                           & sd_rho1,sd_rho2,skip_flag)

!!$                call sdm_outcome_aggregation_shima2(                & 
!!$                                           & sd_re1,sd_re2,  &
!!$                                           & sd_rp1,sd_rp2,  &
!!$                                           & sd_rho1,sd_rho2,skip_flag)

!!$                !! artificial limiter to suppress aggregation of very low density ice (mimicking breakup)
!!$                if(skip_flag .eq. 1) then
!!$                   sd_ncol = i_col-1
!!$                   exit
!!$                end if
!!$
!!$             end do

          end if

          ! update freezing temperature
          sd_tf2 = max(sd_tf1,sd_tf2)

          ! update aerosol components
          do k=1,22
             s = idx_nasl(k)
             dtmp = sd_asl1(s) * real(sd_ncol,kind=RP)
             sd_asl2(s) = sd_asl2(s) + dmask(k) * dtmp
          end do

          !! calculate the new multiplicity
          if( sd_n1 > sd_n2*sd_ncol ) then

             sd_n1 = sd_n1 - sd_n2*sd_ncol

          else
             !! coalescent SDs with same multiplicity

             sd_n1 = int( sd_n2/2, kind=RP )
             sd_n2 = sd_n2 - sd_n1

             if( (sd_li1 .eq. STAT_LIQ) .and. (sd_li2 .eq. STAT_LIQ) ) then
                sd_r1 = sd_r2
                sd_li1 = sd_li2
             else if( (sd_li1 .eq. STAT_ICE) .and. (sd_li2 .eq. STAT_LIQ) ) then
                sd_rho1= sd_rho2
                sd_re1 = sd_re2
                sd_rp1 = sd_rp2
                sd_li1 = sd_li2
                sd_mrime1 = sd_mrime2
                sd_nmono1 = sd_nmono2
             else if( (sd_li1 .eq. STAT_LIQ) .and. (sd_li2 .eq. STAT_ICE) ) then
                sd_rho1= sd_rho2
                sd_re1 = sd_re2
                sd_rp1 = sd_rp2
                sd_li1 = sd_li2
                sd_mrime1 = sd_mrime2
                sd_nmono1 = sd_nmono2
             else if( (sd_li1 .eq. STAT_ICE) .and. (sd_li2 .eq. STAT_ICE) ) then
                sd_rho1= sd_rho2
                sd_re1 = sd_re2
                sd_rp1 = sd_rp2
                sd_li1 = sd_li2
                sd_mrime1 = sd_mrime2
                sd_nmono1 = sd_nmono2
             end if

             sd_tf1 = sd_tf2

             do k=1,22
                s = idx_nasl(k)
                sd_asl1(s) = sd_asl2(s)
             end do

             !! invalid by collisions between SDs with
             !! sd_n1=sd_n2*sd_ncol and sd_n2=1

             if( sd_n1==0 ) then
                sd_rk1 = INVALID
             end if

          end if

          !! This never happens
!!$            !! check muliplicity
!!$
!!$            if( sd_n1>(2.0_RP**63._RP) .or. sd_n2>(2.0_RP**63._RP) ) then
!!$               iexced = -1
!!$               cycle
!!$            end if

          if( sd_n(icptc) > sd_n(icptp) ) then
             sd_n( icptc )  = sd_n1
             sd_r( icptc )  = sd_r1
             sd_rk( icptc ) = sd_rk1
             sd_liqice( icptc ) = sd_li1
             sdi%re( icptc )  = sd_re1
             sdi%rp( icptc )  = sd_rp1
             sdi%rho( icptc ) = sd_rho1
             sdi%tf( icptc )  = sd_tf1
             sdi%mrime( icptc )  = sd_mrime1
             sdi%nmono( icptc )  = sd_nmono1

             sd_n( icptp )  = sd_n2
             sd_r( icptp )  = sd_r2
             sd_liqice( icptp ) = sd_li2
             sdi%re( icptp )  = sd_re2
             sdi%rp( icptp )  = sd_rp2
             sdi%rho( icptp ) = sd_rho2
             sdi%tf( icptp )  = sd_tf2
             sdi%mrime( icptp )  = sd_mrime2
             sdi%nmono( icptp )  = sd_nmono2

             do k=1,22
                s = idx_nasl(k)
                sd_asl( icptc,s ) = sd_asl1(s)
                sd_asl( icptp,s ) = sd_asl2(s)
             end do


          else

             sd_n( icptp )  = sd_n1
             sd_r( icptp )  = sd_r1
             sd_rk( icptp ) = sd_rk1
             sd_liqice( icptp ) = sd_li1
             sdi%re( icptp )  = sd_re1
             sdi%rp( icptp )  = sd_rp1
             sdi%rho( icptp ) = sd_rho1
             sdi%tf( icptp )  = sd_tf1
             sdi%mrime( icptp )  = sd_mrime1
             sdi%nmono( icptp )  = sd_nmono1

             sd_n( icptc )  = sd_n2
             sd_r( icptc )  = sd_r2
             sd_liqice( icptc ) = sd_li2
             sdi%re( icptc )  = sd_re2
             sdi%rp( icptc )  = sd_rp2
             sdi%rho( icptc ) = sd_rho2
             sdi%tf( icptc )  = sd_tf2
             sdi%mrime( icptc )  = sd_mrime2
             sdi%nmono( icptc )  = sd_nmono2

             do k=1,22
                s = idx_nasl(k)
                sd_asl( icptp,s ) = sd_asl1(s)
                sd_asl( icptc,s ) = sd_asl2(s)
             end do

          end if

       end do

    end do

    ! Deallocate
    deallocate( fsort_tag  )
    deallocate( fsort_freq )

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_coales_cold",1,1)
#endif
    return
  end subroutine sdm_coales_cold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outcome_drolet_coalescence(sd_ncol,sd_r1,sd_r2)
    ! calculate the outcome of droplet-droplet coalescence
    use m_sdm_common, only: O_THRD

    integer, intent(in)     :: sd_ncol ! how many times coalescence occurs
    real(RP), intent(in)    :: sd_r1   ! radius of super-droplets with larger multiplicity
    real(RP), intent(inout) :: sd_r2   ! radius of super-droplets with smaller multiplicity

    ! update r
    sd_r2 = exp( O_THRD                                      &
         * log(sd_r1**3*real(sd_ncol,kind=RP)+sd_r2**3) ) 

    return
  end subroutine sdm_outcome_drolet_coalescence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outcome_riming(sd_ncol,sd_r1,sd_r2,    & 
                                      & sd_li1,sd_li2,  &
                                      & sd_re1,sd_re2,  &
                                      & sd_rp1,sd_rp2,  &
                                      & sd_rho1,sd_rho2, &
                                      & sd_vz1,sd_vz2, &
                                      & sd_mrime1,sd_mrime2, &
                                      & sd_nmono1,sd_nmono2, &
                                      & sd_t,sd_rd,sd_p,sd_qv)
    ! calculate the outcome of riming
    use scale_const, only:  &
         rhow_mks => CONST_DWATR, &  ! density of liquid water [kg/m^3]
         rhoi_mks => CONST_DICE      ! density of ice [kg/m^3]

    use m_sdm_common, only: O_THRD,F_THRD,ONE_PI,STAT_ICE,STAT_LIQ,i2

    integer, intent(in)     :: sd_ncol ! how many times coalescence occurs
    real(RP), intent(in)    :: sd_r1   ! radius of liquid droplet with larger multiplicity
    real(RP), intent(in)    :: sd_r2   ! radius of liquid droplet with smaller multiplicity
    real(RP), intent(in)    :: sd_re1  ! equatorial radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_re2  ! equatorial radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rp1  ! polar radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rp2  ! polar radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rho1 ! densitiy of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rho2 ! densitiy of ice particle with larger multiplicity   
    real(RP), intent(in)    :: sd_vz1  ! terminal velocity of super-droplet with larger multiplicity
    real(RP), intent(in)    :: sd_vz2  ! terminal velocity of super-droplet with smaller multiplicity
    real(RP), intent(in)    :: sd_mrime1  ! rime mass of super-droplet with larger multiplicity
    real(RP), intent(inout) :: sd_mrime2  ! rime mass of super-droplet with smaller multiplicity
    integer,  intent(in)    :: sd_nmono1  ! number of monomers of super-droplet with larger multiplicity
    integer,  intent(inout) :: sd_nmono2  ! number of monomers of super-droplet with smaller multiplicity
    integer(i2), intent(in) :: sd_li1 ! phase of super-droplet with larger multiplicity
    integer(i2), intent(inout) :: sd_li2 ! phase of super-droplet with smaller multiplicity
    real(RP), intent(in)  :: sd_t      ! temperature[K] at the location of super-droplets
    real(RP), intent(in)  :: sd_rd     ! moist air density [kg/m3] at the location of super-droplets
    real(RP), intent(in)  :: sd_p      ! pressure of the grid contained the SD [Pa]
    real(RP), intent(in)  :: sd_qv     ! water-vapor of the grid contained the SD [kg/kg]

    real(RP) :: sd_m2_tmp
    real(RP) :: sd_li2_tmp,sd_re2_tmp,sd_rp2_tmp,sd_rho2_tmp,sd_v2_tmp, sd_mrime2_tmp
    integer  :: sd_nmono2_tmp
    real(RP) :: rho_rime ! riming density [kg/m^3]
 
    if( (sd_li1 .eq. STAT_ICE) .and. (sd_li2 .eq. STAT_LIQ) ) then ! riming
       if( sd_r2 > max(sd_re1,sd_rp1) )then ! droplet is larger than ice
          sd_m2_tmp = F_THRD * ONE_PI * (sd_re1**2*sd_rp1*sd_rho1*real(sd_ncol,kind=RP) + sd_r2**3*rhow_mks)
          ! update rho
          sd_rho2_tmp = rhoi_mks
          ! update re
          sd_re2_tmp = (sd_m2_tmp/sd_rho2_tmp/(F_THRD*ONE_PI))**(1.0d0/3.0d0)
          ! update rp
          sd_rp2_tmp = sd_re2_tmp
          ! update phase
          sd_li2_tmp = STAT_ICE
          ! update rime mass (Large droplet freezes on contant with small ice is also considered as riming)
          sd_mrime2_tmp = sd_mrime1*real(sd_ncol,kind=RP) + F_THRD*ONE_PI*(sd_r2**3*rhow_mks)
          ! update number of monomers
          sd_nmono2_tmp = sd_nmono1*sd_ncol
       else ! ice is larger than droplet
          call sdm_rho_riming_Heymsfield_Pflaum_1985(rho_rime,sd_r2,sd_re1,sd_rp1,sd_vz2,sd_vz1,sd_t,sd_rd,sd_p,sd_qv)
!          if(sd_re1 > sd_rp1) then ! oblate
          ! mimick tumbling if the ice is quasi spherical (Jensen and Harrington, 2015)
          if((0.8d0 * sd_re1 > sd_rp1) .or. ((1.25d0 * sd_re1 > sd_rp1) .and. (sd_rp1 >= sd_re1))) then ! oblate and quasi spherical prolate
             sd_m2_tmp = F_THRD * ONE_PI * (sd_re1**2*sd_rp1*sd_rho1*real(sd_ncol,kind=RP) + sd_r2**3*rhow_mks)
             sd_v2_tmp = F_THRD * ONE_PI * (sd_re1**2*sd_rp1*real(sd_ncol,kind=RP) + sd_r2**3*(rhow_mks/rho_rime))
             ! update re
             sd_re2_tmp = sd_re1
             ! update rho
             sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
             ! update rp
             sd_rp2_tmp = sd_v2_tmp/F_THRD/ONE_PI/(sd_re2_tmp**2)
             ! update phase
             sd_li2_tmp = STAT_ICE
             ! update rime mass
             sd_mrime2_tmp = sd_mrime1*real(sd_ncol,kind=RP) + F_THRD*ONE_PI*(sd_r2**3*rhow_mks)
             ! update number of monomers
             sd_nmono2_tmp = sd_nmono1*sd_ncol
          else ! prolate and quasi spherical oblate
             sd_m2_tmp = F_THRD * ONE_PI * (sd_re1**2*sd_rp1*sd_rho1*real(sd_ncol,kind=RP) + sd_r2**3*rhow_mks)
             sd_v2_tmp = F_THRD * ONE_PI * (sd_re1**2*sd_rp1*real(sd_ncol,kind=RP) + sd_r2**3*(rhow_mks/rho_rime))
             ! update rp
             sd_rp2_tmp = sd_rp1
             ! update rho
             sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
             ! update re
             sd_re2_tmp = sqrt(sd_v2_tmp/F_THRD/ONE_PI/sd_rp2_tmp)
             ! update phase
             sd_li2_tmp = STAT_ICE
             ! update rime mass
             sd_mrime2_tmp = sd_mrime1*real(sd_ncol,kind=RP) + F_THRD*ONE_PI*(sd_r2**3*rhow_mks)
             ! update number of monomers
             sd_nmono2_tmp = sd_nmono1*sd_ncol
          end if
       end if
    else if( (sd_li1 .eq. STAT_LIQ) .and. (sd_li2 .eq. STAT_ICE) ) then ! riming
       if( sd_r1 > max(sd_re2,sd_rp2) )then ! droplet is larger than ice
          sd_m2_tmp = F_THRD * ONE_PI * (sd_r1**3*rhow_mks*real(sd_ncol,kind=RP) + sd_re2**2*sd_rp2*sd_rho2)
          ! update rho
          sd_rho2_tmp = rhoi_mks
          ! update re
          sd_re2_tmp = (sd_m2_tmp/sd_rho2_tmp/(F_THRD*ONE_PI))**(1.0d0/3.0d0)
          ! update rp
          sd_rp2_tmp = sd_re2_tmp
          ! update phase
          sd_li2_tmp = STAT_ICE
          ! update rime mass (Large droplet freezes on contant with small ice is also considered as riming)
          sd_mrime2_tmp = F_THRD*ONE_PI*(sd_r1**3*rhow_mks)*real(sd_ncol,kind=RP) + sd_mrime2
          ! update number of monomers
          sd_nmono2_tmp = sd_nmono2
       else ! ice is larger than droplet
          call sdm_rho_riming_Heymsfield_Pflaum_1985(rho_rime,sd_r1,sd_re2,sd_rp2,sd_vz1,sd_vz2,sd_t,sd_rd,sd_p,sd_qv)
!          if(sd_re2 > sd_rp2) then ! oblate
          ! mimick tumbling if the ice is quasi spherical (Jensen and Harrington, 2015)
          if((0.8d0 * sd_re2 > sd_rp2) .or. ((1.25d0 * sd_re2 > sd_rp2) .and. (sd_rp2 >= sd_re2))) then ! oblate and quasi spherical prolate
             sd_m2_tmp = F_THRD * ONE_PI * (sd_r1**3*rhow_mks*real(sd_ncol,kind=RP) + sd_re2**2*sd_rp2*sd_rho2)
             sd_v2_tmp = F_THRD * ONE_PI * (sd_r1**3*(rhow_mks/rho_rime)*real(sd_ncol,kind=RP) + sd_re2**2*sd_rp2)
             ! update re
             sd_re2_tmp = sd_re2
             ! update rho
             sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
             ! update rp
             sd_rp2_tmp = sd_v2_tmp/F_THRD/ONE_PI/(sd_re2_tmp**2)
             ! update phase
             sd_li2_tmp = STAT_ICE
             ! update rime mass
             sd_mrime2_tmp = F_THRD*ONE_PI*(sd_r1**3*rhow_mks)*real(sd_ncol,kind=RP) + sd_mrime2
             ! update number of monomers
             sd_nmono2_tmp = sd_nmono2
          else ! prolate and quasi spherical oblate
             sd_m2_tmp = F_THRD * ONE_PI * (sd_r1**3*rhow_mks*real(sd_ncol,kind=RP) + sd_re2**2*sd_rp2*sd_rho2)
             sd_v2_tmp = F_THRD * ONE_PI * (sd_r1**3*(rhow_mks/rho_rime)*real(sd_ncol,kind=RP) + sd_re2**2*sd_rp2)
             ! update rp
             sd_rp2_tmp = sd_rp2
             ! update rho
             sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
             ! update re
             sd_re2_tmp = sqrt(sd_v2_tmp/F_THRD/ONE_PI/sd_rp2_tmp)
             ! update phase
             sd_li2_tmp = STAT_ICE
             ! update rime mass
             sd_mrime2_tmp = F_THRD*ONE_PI*(sd_r1**3*rhow_mks)*real(sd_ncol,kind=RP) + sd_mrime2
             ! update number of monomers
             sd_nmono2_tmp = sd_nmono2
          end if
       end if
    end if

    sd_li2  = sd_li2_tmp
    sd_re2  = sd_re2_tmp
    sd_rp2  = sd_rp2_tmp
    sd_rho2 = sd_rho2_tmp 
    sd_mrime2= sd_mrime2_tmp
    sd_nmono2= sd_nmono2_tmp

    return
  end subroutine sdm_outcome_riming
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outcome_aggregation_shima2(                 & 
                                           & sd_re1,sd_re2,  &
                                           & sd_rp1,sd_rp2,  &
                                           & sd_rho1,sd_rho2)
    ! calculate the outcome of aggregation
    use m_sdm_common, only: F_THRD,ONE_PI

    real(RP), intent(in)    :: sd_re1  ! equatorial radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_re2  ! equatorial radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rp1  ! polar radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rp2  ! polar radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rho1 ! densitiy of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rho2 ! densitiy of ice particle with larger multiplicity   
 
    real(RP),parameter :: const_beta = 2.2d0     ! power low of the mass of snow aggregate. (Schmitt and Hemsfield 2010)
    real(RP),parameter :: const_lambda = 3.0d0   ! factor for oblate collector (Shima 2017) 
    real(RP),parameter :: const_eta    = 1.5d0   ! factor for prolate collector (Shima 2017)
    
    real(RP) :: sd_re2_tmp,sd_rp2_tmp,sd_rho2_tmp
    real(RP) :: sd_v2_tmp,sd_m2_tmp,sd_m1_tmp,sd_m2_org_tmp

    if( max(sd_re1,sd_rp1) > max(sd_re2,sd_rp2) )then ! ice1 is bigger than ice2 
       if(sd_re1 > sd_rp1) then ! collector is oblate
          sd_re2_tmp  = sd_re1
          sd_m2_tmp   = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1 + sd_re2**2*sd_rp2*sd_rho2)
          sd_m1_tmp   = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1)
          sd_rp2_tmp  = sd_rp1*(sd_m2_tmp/sd_m1_tmp)**(const_lambda/const_beta)
          sd_v2_tmp   = F_THRD*ONE_PI*(sd_re2_tmp**2*sd_rp2_tmp)
          sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
       else ! collector is prolate
          sd_rp2_tmp  = sd_rp1
          sd_m2_tmp   = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1 + sd_re2**2*sd_rp2*sd_rho2)
          sd_m1_tmp   = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1)
          sd_re2_tmp  = sd_re1*(sd_m2_tmp/sd_m1_tmp)**(const_eta/const_beta)
          sd_v2_tmp   = F_THRD*ONE_PI*(sd_re2_tmp**2*sd_rp2_tmp)
          sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
       end if
    else ! max(sd_re1,sd_rp1) <= max(sd_re2,sd_rp2) ! ice2 is bigger than ice1 
       if(sd_re2 > sd_rp2) then ! collector is oblate
          sd_re2_tmp  = sd_re2
          sd_m2_tmp   = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1 + sd_re2**2*sd_rp2*sd_rho2)
          sd_m2_org_tmp = F_THRD*ONE_PI*(sd_re2**2*sd_rp2*sd_rho2)
          sd_rp2_tmp  = sd_rp2*(sd_m2_tmp/sd_m2_org_tmp)**(const_lambda/const_beta)
          sd_v2_tmp   = F_THRD*ONE_PI*(sd_re2_tmp**2*sd_rp2_tmp)
          sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
       else ! collector is prolate
          sd_rp2_tmp  = sd_rp2
          sd_m2_tmp   = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1 + sd_re2**2*sd_rp2*sd_rho2)
          sd_m2_org_tmp = F_THRD*ONE_PI*(sd_re2**2*sd_rp2*sd_rho2)
          sd_re2_tmp  = sd_re2*(sd_m2_tmp/sd_m2_org_tmp)**(const_eta/const_beta)
          sd_v2_tmp   = F_THRD*ONE_PI*(sd_re2_tmp**2*sd_rp2_tmp)
          sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
       end if
    end if

    sd_re2  = sd_re2_tmp
    sd_rp2  = sd_rp2_tmp
    sd_rho2 = sd_rho2_tmp 

    return
  end subroutine sdm_outcome_aggregation_shima2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outcome_aggregation_Chen_Lamb_1994_mod(     & 
                                           & sd_re1,sd_re2,  &
                                           & sd_rp1,sd_rp2,  &
                                           & sd_rho1,sd_rho2,skip_flag)
    ! calculate the outcome of aggregation
    use scale_const, only:  &
         rhoi_mks => CONST_DICE  ! density of ice [kg/m^3]
    use m_sdm_common, only: F_THRD,ONE_PI

    real(RP), intent(in)    :: sd_re1  ! equatorial radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_re2  ! equatorial radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rp1  ! polar radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rp2  ! polar radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rho1 ! densitiy of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rho2 ! densitiy of ice particle with larger multiplicity   
    integer, intent(out) :: skip_flag  ! flag to save whether the aggregation is canceled or not
 
    real(RP),parameter :: const_s = 0.6d0  ! parameter to determine the distance of two particles (Chen and Lamb 1994)
    real(RP) :: const_theta_max            ! maximum contact angle of two particles (Chen and Lamb 1994)

    real(RP) :: theta         ! contact angle of two particles
    real(RP) :: height,width  ! height and width of the collected ice, rotated by an angle of theta 
    real(RP) :: distance_L    ! distance between the center of the two ice particles
    real(RP) :: excess_length ! the lenght of the collected particle that runs off the edge of the collector particle 
    real(RP) :: sd_rho_ave    ! volume weighted average density of the two ice particles

    real(RP) :: sd_re2_tmp,sd_rp2_tmp,sd_rho2_tmp
    real(RP) :: sd_v2_tmp,sd_m2_tmp

    const_theta_max = ONE_PI/4.0d0 ! maximum contact angle of two particles (Chen and Lamb 1994)

    sd_rho_ave = (sd_re1**2*sd_rp1*sd_rho1 + sd_re2**2*sd_rp2*sd_rho2)/(sd_re1**2*sd_rp1 + sd_re2**2*sd_rp2)

    if( max(sd_re1,sd_rp1) > max(sd_re2,sd_rp2) )then ! ice1 is bigger than ice2 
       theta = (const_theta_max)*(1.0d0-sd_rho_ave/rhoi_mks) ! contact angle
       height = sqrt(2.0d0 * ((sd_re2**2+sd_rp2**2) - abs(sd_re2**2-sd_rp2**2)*cos(2.0d0*theta)) ) 
       width  = sqrt(2.0d0 * ((sd_re2**2+sd_rp2**2) + abs(sd_re2**2-sd_rp2**2)*cos(2.0d0*theta)) ) 

       if(sd_re1 > sd_rp1) then ! collector is oblate
          distance_L = const_s * (sd_re1 + max(sd_re2,sd_rp2))
          excess_length = max( (distance_L+width/2.0d0)-sd_re1, 0.0d0 )
          sd_rp2_tmp = sd_rp1 + height/2.0d0
          sd_re2_tmp = sqrt( sd_re1*(sd_re1+excess_length/2.0d0) ) 
          sd_m2_tmp   = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1 + sd_re2**2*sd_rp2*sd_rho2)
          sd_v2_tmp   = F_THRD*ONE_PI*(sd_re2_tmp**2*sd_rp2_tmp)
          sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
       else ! collector is prolate
          distance_L = const_s * (sd_rp1 + max(sd_re2,sd_rp2))
          excess_length = max( (distance_L+width/2.0d0)-sd_rp1, 0.0d0 )
          sd_rp2_tmp = sd_rp1 + excess_length/2.0d0
          sd_re2_tmp = sqrt( sd_re1*(sd_re1+height/2.0d0) ) 
          sd_m2_tmp   = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1 + sd_re2**2*sd_rp2*sd_rho2)
          sd_v2_tmp   = F_THRD*ONE_PI*(sd_re2_tmp**2*sd_rp2_tmp)
          sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
       end if

    else ! max(sd_re1,sd_rp1) <= max(sd_re2,sd_rp2) ! ice2 is bigger than ice1 
       theta = (const_theta_max)*(1.0d0-sd_rho_ave/rhoi_mks) ! contact angle
       height = sqrt(2.0d0 * ((sd_re1**2+sd_rp1**2) - abs(sd_re1**2-sd_rp1**2)*cos(2.0d0*theta)) ) 
       width  = sqrt(2.0d0 * ((sd_re1**2+sd_rp1**2) + abs(sd_re1**2-sd_rp1**2)*cos(2.0d0*theta)) ) 

       if(sd_re2 > sd_rp2) then ! collector is oblate
          distance_L = const_s * (sd_re2 + max(sd_re1,sd_rp1))
          excess_length = max( (distance_L+width/2.0d0)-sd_re2, 0.0d0 )
          sd_rp2_tmp = sd_rp2 + height/2.0d0
          sd_re2_tmp = sqrt( sd_re2*(sd_re2+excess_length/2.0d0) ) 
          sd_m2_tmp   = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1 + sd_re2**2*sd_rp2*sd_rho2)
          sd_v2_tmp   = F_THRD*ONE_PI*(sd_re2_tmp**2*sd_rp2_tmp)
          sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
       else ! collector is prolate
          distance_L = const_s * (sd_rp2 + max(sd_re1,sd_rp1))
          excess_length = max( (distance_L+width/2.0d0)-sd_rp2, 0.0d0 )
          sd_rp2_tmp = sd_rp2 + excess_length/2.0d0
          sd_re2_tmp = sqrt( sd_re2*(sd_re2+height/2.0d0) ) 
          sd_m2_tmp   = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1 + sd_re2**2*sd_rp2*sd_rho2)
          sd_v2_tmp   = F_THRD*ONE_PI*(sd_re2_tmp**2*sd_rp2_tmp)
          sd_rho2_tmp = sd_m2_tmp/sd_v2_tmp
       end if
    end if

    !! artificial limiter to suppress aggregation of very low density ice (mimicking breakup)
    skip_flag=0
    if(sd_rho2_tmp<10.0d0) then
       skip_flag=1
       return
    end if

    sd_re2  = sd_re2_tmp
    sd_rp2  = sd_rp2_tmp
    sd_rho2 = sd_rho2_tmp 

    return
  end subroutine sdm_outcome_aggregation_Chen_Lamb_1994_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outcome_aggregation_shima4(     &
                                           & sd_re1,sd_re2,  &
                                           & sd_rp1,sd_rp2,  &
                                           & sd_rho1,sd_rho2,skip_flag)
    ! calculate the outcome of aggregation
    use scale_const, only:  &
         rhoi_mks => CONST_DICE  ! density of ice [kg/m^3]
    use m_sdm_common, only: F_THRD,ONE_PI

    real(RP), intent(in)    :: sd_re1  ! equatorial radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_re2  ! equatorial radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rp1  ! polar radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rp2  ! polar radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rho1 ! densitiy of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rho2 ! densitiy of ice particle with larger multiplicity
    integer, intent(out) :: skip_flag  ! flag to save whether the aggregation is canceled or not

    real(RP),parameter :: rhoi_min = 10.0_RP ! minimum apparent density of ice aggregate [kg/m^3]

    real(RP) :: sd_rho_ave    ! volume weighted average density of the two ice particles

    real(RP) :: sd_re_new,sd_rp_new,sd_rho_new
    real(RP) :: sd_v_new,sd_m_new
    real(RP) :: sd_v_new_min,sd_re_new_min,sd_re_new_max,sd_rp_new_min,sd_rp_new_max

    sd_m_new  = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1 + sd_re2**2*sd_rp2*sd_rho2)
    sd_v_new_min  = F_THRD*ONE_PI*(sd_re1**2*sd_rp1 + sd_re2**2*sd_rp2)
    sd_rho_ave = sd_m_new/sd_v_new_min

    if( max(sd_re1,sd_rp1) > max(sd_re2,sd_rp2) )then ! ice1 is bigger than ice2 
       if(sd_re1 > sd_rp1) then ! collector is oblate
          sd_re_new = sd_re1
          sd_rp_new_min = sd_v_new_min/(F_THRD*ONE_PI*sd_re_new**2)
          sd_rp_new_max = sd_rp1 + min(sd_re2,sd_rp2)
          sd_rp_new = ((rhoi_mks-sd_rho_ave)*sd_rp_new_min+(sd_rho_ave-rhoi_min)*sd_rp_new_max)/(rhoi_mks-rhoi_min)
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       else ! collector is prolate
          sd_rp_new = sd_rp1
          sd_re_new_min = sqrt(sd_v_new_min/(F_THRD*ONE_PI*sd_rp_new))
          sd_re_new_max = sqrt(max(sd_re1,sd_re2,sd_rp2)*(sd_re1 + min(sd_re2,sd_rp2)))
          sd_re_new = ((rhoi_mks-sd_rho_ave)*sd_re_new_min+(sd_rho_ave-rhoi_min)*sd_re_new_max)/(rhoi_mks-rhoi_min)
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       end if
    else ! max(sd_re1,sd_rp1) <= max(sd_re2,sd_rp2) ! ice2 is bigger than ice1
       if(sd_re2 > sd_rp2) then ! collector is oblate
          sd_re_new = sd_re2
          sd_rp_new_min = sd_v_new_min/(F_THRD*ONE_PI*sd_re_new**2)
          sd_rp_new_max = sd_rp2 + min(sd_re1,sd_rp1)
          sd_rp_new = ((rhoi_mks-sd_rho_ave)*sd_rp_new_min+(sd_rho_ave-rhoi_min)*sd_rp_new_max)/(rhoi_mks-rhoi_min)
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       else ! collector is prolate
          sd_rp_new = sd_rp2
          sd_re_new_min = sqrt(sd_v_new_min/(F_THRD*ONE_PI*sd_rp_new))
          sd_re_new_max = sqrt(max(sd_re2,sd_re1,sd_rp1)*(sd_re2 + min(sd_re1,sd_rp1)))
          sd_re_new = ((rhoi_mks-sd_rho_ave)*sd_re_new_min+(sd_rho_ave-rhoi_min)*sd_re_new_max)/(rhoi_mks-rhoi_min)
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       end if
    end if

    sd_re2  = sd_re_new
    sd_rp2  = sd_rp_new
    sd_rho2 = sd_rho_new

    skip_flag = 0

    return
  end subroutine sdm_outcome_aggregation_shima4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outcome_aggregation_shima4_multi(sd_ncol,   &
                                           & sd_re1,sd_re2,  &
                                           & sd_rp1,sd_rp2,  &
                                           & sd_rho1,sd_rho2)
    ! calculate the outcome of aggregation
    use scale_const, only:  &
         rhoi_mks => CONST_DICE  ! density of ice [kg/m^3]
    use m_sdm_common, only: F_THRD,ONE_PI

    integer, intent(in)     :: sd_ncol ! how many times coalescence occurs
    real(RP), intent(in)    :: sd_re1  ! equatorial radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_re2  ! equatorial radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rp1  ! polar radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rp2  ! polar radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rho1 ! densitiy of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rho2 ! densitiy of ice particle with larger multiplicity

    real(RP),parameter :: rhoi_min = 10.0_RP ! minimum apparent density of ice aggregate [kg/m^3]

    real(RP) :: sd_rho_ave    ! volume weighted average density of the two ice particles

    real(RP) :: sd_re_new,sd_rp_new,sd_rho_new
    real(RP) :: sd_v_new,sd_m_new
    real(RP) :: sd_v_new_min,sd_re_new_min,sd_re_new_max,sd_rp_new_min,sd_rp_new_max

    sd_m_new  = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1*real(sd_ncol,kind=RP) + sd_re2**2*sd_rp2*sd_rho2)
    sd_v_new_min  = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*real(sd_ncol,kind=RP) + sd_re2**2*sd_rp2)
    sd_rho_ave = sd_m_new/sd_v_new_min

    if( max(sd_re1,sd_rp1) > max(sd_re2,sd_rp2) )then ! ice1 is bigger than ice2 
       if(sd_re1 > sd_rp1) then ! collector is oblate
          sd_re_new = sd_re1
          sd_rp_new_min = sd_v_new_min/(F_THRD*ONE_PI*sd_re_new**2)
          sd_rp_new_max = sd_rp1*real(sd_ncol,kind=RP) + min(sd_re2,sd_rp2)
          sd_rp_new = ((rhoi_mks-sd_rho_ave)*sd_rp_new_min+(sd_rho_ave-rhoi_min)*sd_rp_new_max)/(rhoi_mks-rhoi_min)
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       else ! collector is prolate
          sd_rp_new = sd_rp1
          sd_re_new_min = sqrt(sd_v_new_min/(F_THRD*ONE_PI*sd_rp_new))
          sd_re_new_max = sqrt(max(sd_re1,sd_re2,sd_rp2)*(sd_re1*real(sd_ncol,kind=RP) + min(sd_re2,sd_rp2)))
          sd_re_new = ((rhoi_mks-sd_rho_ave)*sd_re_new_min+(sd_rho_ave-rhoi_min)*sd_re_new_max)/(rhoi_mks-rhoi_min)
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       end if
    else ! max(sd_re1,sd_rp1) <= max(sd_re2,sd_rp2) ! ice2 is bigger than ice1
       if(sd_re2 > sd_rp2) then ! collector is oblate
          sd_re_new = sd_re2
          sd_rp_new_min = sd_v_new_min/(F_THRD*ONE_PI*sd_re_new**2)
          sd_rp_new_max = sd_rp2 + min(sd_re1,sd_rp1)*real(sd_ncol,kind=RP)
          sd_rp_new = ((rhoi_mks-sd_rho_ave)*sd_rp_new_min+(sd_rho_ave-rhoi_min)*sd_rp_new_max)/(rhoi_mks-rhoi_min)
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       else ! collector is prolate
          sd_rp_new = sd_rp2
          sd_re_new_min = sqrt(sd_v_new_min/(F_THRD*ONE_PI*sd_rp_new))
          sd_re_new_max = sqrt(max(sd_re2,sd_re1,sd_rp1)*(sd_re2 + min(sd_re1,sd_rp1)*real(sd_ncol,kind=RP)))
          sd_re_new = ((rhoi_mks-sd_rho_ave)*sd_re_new_min+(sd_rho_ave-rhoi_min)*sd_re_new_max)/(rhoi_mks-rhoi_min)
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       end if
    end if

    sd_re2  = sd_re_new
    sd_rp2  = sd_rp_new
    sd_rho2 = sd_rho_new

    return
  end subroutine sdm_outcome_aggregation_shima4_multi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outcome_aggregation_shima5_multi(sd_ncol,   &
                                           & sd_re1,sd_re2,  &
                                           & sd_rp1,sd_rp2,  &
                                           & sd_rho1,sd_rho2)
    ! calculate the outcome of aggregation
    use scale_const, only:  &
         rhoi_mks => CONST_DICE  ! density of ice [kg/m^3]
    use m_sdm_common, only: F_THRD,ONE_PI

    integer, intent(in)     :: sd_ncol ! how many times coalescence occurs
    real(RP), intent(in)    :: sd_re1  ! equatorial radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_re2  ! equatorial radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rp1  ! polar radius of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rp2  ! polar radius of ice particle with smaller multiplicity
    real(RP), intent(in)    :: sd_rho1 ! densitiy of ice particle with larger multiplicity
    real(RP), intent(inout) :: sd_rho2 ! densitiy of ice particle with smaller multiplicity

    real(RP),parameter :: rhoi_min = 10.0_RP ! minimum apparent density of ice aggregate [kg/m^3]

    real(RP) :: sd_rho_ave    ! volume weighted average density of the two ice particles

    real(RP) :: sd_re_new,sd_rp_new,sd_rho_new
    real(RP) :: sd_v_new,sd_m_new
    real(RP) :: sd_v_new_min,sd_re_new_min,sd_re_new_max,sd_rp_new_min,sd_rp_new_max

    sd_m_new  = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*sd_rho1*real(sd_ncol,kind=RP) + sd_re2**2*sd_rp2*sd_rho2)
    sd_v_new_min  = F_THRD*ONE_PI*(sd_re1**2*sd_rp1*real(sd_ncol,kind=RP) + sd_re2**2*sd_rp2)
    sd_rho_ave = sd_m_new/sd_v_new_min

    if( max(sd_re1,sd_rp1) > max(sd_re2,sd_rp2) )then ! ice1 is bigger than ice2 
       if(sd_re1 > sd_rp1) then ! collector is oblate
          sd_re_new = sd_re1
          sd_rp_new_min = sd_v_new_min/(F_THRD*ONE_PI*sd_re_new**2)
          sd_rp_new_max = sd_rp1*real(sd_ncol,kind=RP) + min(sd_re2,sd_rp2)
          sd_rp_new = (rhoi_mks-rhoi_min) / &
               &      ((rhoi_mks-sd_rho_ave)/sd_rp_new_min + &
               &       (sd_rho_ave-rhoi_min)/sd_rp_new_max  &
               &      )
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       else ! collector is prolate
          sd_rp_new = sd_rp1
          sd_re_new_min = sqrt(sd_v_new_min/(F_THRD*ONE_PI*sd_rp_new))
          sd_re_new_max = sqrt(max(sd_re1,sd_re2,sd_rp2)*(sd_re1*real(sd_ncol,kind=RP) + min(sd_re2,sd_rp2)))
          sd_re_new = sqrt((rhoi_mks-rhoi_min) / &
               &           ((rhoi_mks-sd_rho_ave)/sd_re_new_min**2 + &
               &            (sd_rho_ave-rhoi_min)/sd_re_new_max**2   &
               &           )                                         &
               &          )
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       end if
    else ! max(sd_re1,sd_rp1) <= max(sd_re2,sd_rp2) ! ice2 is bigger than ice1
       if(sd_re2 > sd_rp2) then ! collector is oblate
          sd_re_new = sd_re2
          sd_rp_new_min = sd_v_new_min/(F_THRD*ONE_PI*sd_re_new**2)
          sd_rp_new_max = sd_rp2 + min(sd_re1,sd_rp1)*real(sd_ncol,kind=RP)
          sd_rp_new = (rhoi_mks-rhoi_min) / &
               &      ((rhoi_mks-sd_rho_ave)/sd_rp_new_min + &
               &       (sd_rho_ave-rhoi_min)/sd_rp_new_max  &
               &      )
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       else ! collector is prolate
          sd_rp_new = sd_rp2
          sd_re_new_min = sqrt(sd_v_new_min/(F_THRD*ONE_PI*sd_rp_new))
          sd_re_new_max = sqrt(max(sd_re2,sd_re1,sd_rp1)*(sd_re2 + min(sd_re1,sd_rp1)*real(sd_ncol,kind=RP)))
          sd_re_new = sqrt((rhoi_mks-rhoi_min) / &
               &           ((rhoi_mks-sd_rho_ave)/sd_re_new_min**2 + &
               &            (sd_rho_ave-rhoi_min)/sd_re_new_max**2   &
               &           )                                         &
               &          )
          sd_v_new   = F_THRD*ONE_PI*(sd_re_new**2*sd_rp_new)
          sd_rho_new = sd_m_new/sd_v_new
       end if
    end if

    sd_re2  = sd_re_new
    sd_rp2  = sd_rp_new
    sd_rho2 = sd_rho_new

    return
  end subroutine sdm_outcome_aggregation_shima5_multi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine coalescence_kernel_golovin(c_kernel,sd_r1,sd_r2)
    ! coalescence kernel of Golovin (1963)
    use m_sdm_common, only: F_THRD,ONE_PI

    real(RP), intent(out)   :: c_kernel  ! coalescence kernel
    real(RP), intent(in)    :: sd_r1  ! radius of particle 1
    real(RP), intent(in)    :: sd_r2  ! radius of particle 2

    c_kernel = 1500.0_RP * F_THRD * ONE_PI              &
         * ( sd_r1*sd_r1*sd_r1 + sd_r2*sd_r2*sd_r2 )

    return
  end subroutine coalescence_kernel_golovin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine coalescence_efficiency_long(c_eff,sd_r1,sd_r2)
    ! coalescence efficiency of Long (1974)
    ! note that the form of Long kernel is much simpler
    real(RP), intent(out)   :: c_eff  ! coalescence efficiency
    real(RP), intent(in)    :: sd_r1  ! radius of particle 1
    real(RP), intent(in)    :: sd_r2  ! radius of particle 2

    real(RP) :: r_large,r_small

    r_large = max( sd_r1, sd_r2 )  !! large
    r_small = min( sd_r1, sd_r2 )  !! small

    if( r_large <= 5.E-5_RP ) then
       c_eff = 4.5E+8_RP * ( r_large*r_large )                 &
            * ( 1.0_RP - 3.E-6_RP/(max(3.01E-6_RP,r_small)) )
    else
       c_eff = 1.d0
    end if

    return
  end subroutine coalescence_efficiency_long
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine coalescence_efficiency_davis_jonas_hall(c_eff,sd_r1_org,sd_r2_org)
    ! coalescence efficiency of Davis (1972), Jonas (1972), and Hall (1980)
    ! This is proposed in Seesselberg et al. (1996)
    use m_sdm_common, only: &
        m2micro,r0col,ratcol,ecoll

    real(RP), intent(out)   :: c_eff  ! coalescence efficiency
    real(RP), intent(in)    :: sd_r1_org  ! radius of particle 1
    real(RP), intent(in)    :: sd_r2_org  ! radius of particle 2

    real(RP) :: sd_r1, sd_r2
    real(RP) :: rq      ! radius ratio of a pair of super-droplets
    real(RP) :: ek      ! temporary
    real(RP) :: p       ! temporary
    real(RP) :: q       ! temporary
    integer :: irr, iqq           ! index of coef.

    sd_r1 = max( sd_r1_org, sd_r2_org )  !! large
    sd_r2   = min( sd_r1_org, sd_r2_org )  !! small

    rq    = sd_r2 / sd_r1
    sd_r1 = sd_r1 * m2micro    !! [m] => [micro-m]
    sd_r2 = sd_r2 * m2micro    !! [m] => [micro-m]

    !! Get index of the array {r0,rat}.

    if( sd_r1 <= r0col(1) ) then
       irr = 1
    else if( sd_r1 <= r0col(2) ) then
       irr = 2
    else if( sd_r1 <= r0col(3) ) then
       irr = 3
    else if( sd_r1 <= r0col(4) ) then
       irr = 4
    else if( sd_r1 <= r0col(5) ) then
       irr = 5
    else if( sd_r1 <= r0col(6) ) then
       irr = 6
    else if( sd_r1 <= r0col(7) ) then
       irr = 7
    else if( sd_r1 <= r0col(8) ) then
       irr = 8
    else if( sd_r1 <= r0col(9) ) then
       irr = 9
    else if( sd_r1 <= r0col(10) ) then
       irr = 10
    else if( sd_r1 <= r0col(11) ) then
       irr = 11
    else if( sd_r1 <= r0col(12) ) then
       irr = 12
    else if( sd_r1 <= r0col(13) ) then
       irr = 13
    else if( sd_r1 <= r0col(14) ) then
       irr = 14
    else if( sd_r1 <= r0col(15) ) then
       irr = 15
    else
       irr = 16
    end if
    
    if( rq <= ratcol(2) ) then
       iqq = 2
    else if( rq <= ratcol(3) ) then
       iqq = 3
    else if( rq <= ratcol(4) ) then
       iqq = 4
    else if( rq <= ratcol(5) ) then
       iqq = 5
    else if( rq <= ratcol(6) ) then
       iqq = 6
    else if( rq <= ratcol(7) ) then
       iqq = 7
    else if( rq <= ratcol(8) ) then
       iqq = 8
    else if( rq <= ratcol(9) ) then
       iqq = 9
    else if( rq <= ratcol(10) ) then
       iqq = 10
    else if( rq <= ratcol(11) ) then
       iqq = 11
    else if( rq <= ratcol(12) ) then
       iqq = 12
    else if( rq <= ratcol(13) ) then
       iqq = 13
    else if( rq <= ratcol(14) ) then
       iqq = 14
    else if( rq <= ratcol(15) ) then
       iqq = 15
    else if( rq <= ratcol(16) ) then
       iqq = 16
    else if( rq <= ratcol(17) ) then
       iqq = 17
    else if( rq <= ratcol(18) ) then
       iqq = 18
    else if( rq <= ratcol(19) ) then
       iqq = 19
    else if( rq <= ratcol(20) ) then
       iqq = 20
    else
       iqq = 21
    end if
    !! Get c_rate

    if( irr>=16 ) then

       q  = (rq-ratcol(iqq-1)) / (ratcol(iqq)-ratcol(iqq-1))
       ek = (1.0_RP-q)*ecoll(15,iqq-1) + q*ecoll(15,iqq)

       c_eff = min( ek, 1.0_RP )

    else if( irr>=2 .and. irr<16 ) then

       p = (sd_r1-r0col(irr-1))/(r0col(irr)-r0col(irr-1))
       q = (rq-ratcol(iqq-1))/(ratcol(iqq)-ratcol(iqq-1))

       c_eff = (1.0_RP-p)*(1.0_RP-q)*ecoll(irr-1,iqq-1)     &
            + p*(1.0_RP-q)*ecoll(irr,iqq-1)              &
            + q*(1.0_RP-p)*ecoll(irr-1,iqq)              &
            + p*q*ecoll(irr,iqq)

    else

       q = (rq-ratcol(iqq-1))/(ratcol(iqq)-ratcol(iqq-1))

       c_eff = (1.0_RP-q)*ecoll(1,iqq-1) + q*ecoll(1,iqq)

    end if

    return
  end subroutine coalescence_efficiency_davis_jonas_hall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine riming_efficiency_Beard_Grover_1974(E_rime,p,nre,n_mixFr)
    ! riming efficiency of Beard and Grover (1974)
    use m_sdm_common, only: ONE_PI
  
    real(DP),intent(out) :: E_rime

    real(DP),parameter :: A0 = -0.1007d0, &
                          A1 = -0.358d0,  &
                          A2 =  0.0261d0

    real(DP),parameter :: B0 =  0.1465d0, &
                          B1 =  1.302d0,  &
                          B2 = -0.607d0,  &
                          B3 =  0.293d0

    real(DP),intent(in) :: p         ! (radius of small particle)/(radius of large droplet)
    real(DP),intent(in) :: nre       ! Reynolds number
    real(DP),intent(in) :: n_mixFr   ! mixed Froude number

    real(DP) :: G,F,K0,H,yc0,z

    if( nre > 400)then
       K0 = 0.21d0
    else
       F = log(nre)
       G = A0 + A1*F + A2*F**2
       K0 = exp(G)
    end if
    z = log(n_mixFr/K0)
    H = max(B0 + B1*z + B2*z**2 + B3*z**3,0.0d0)
    yc0 = (2.0d0/ONE_PI)*atan(H)
    E_rime = (yc0+p)**2/(1.0d0+p)**2

  end subroutine riming_efficiency_Beard_Grover_1974
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine riming_efficiency_Erfani_Mitchell_2017_plate(E_rime,nre,n_mixFr)
  implicit none
  
  real(DP),intent(out) :: E_rime

  real(DP),intent(in) :: nre       ! Reynolds number
  real(DP),intent(in) :: n_mixFr   ! mixed Froude number

  real(DP) :: nre_lmtd ! min(nre,120.0d0)
  real(DP) :: n_mixFr_lmtd ! min(nre,35.0d0)
  real(DP) :: K_thres, K_crit
  real(DP) :: E_rime1,E_rime2
  
  if(nre<1.0d0) then
     E_rime = 0.0d0
     return
  end if
  
  nre_lmtd = min(nre,120.0d0)
  n_mixFr_lmtd = min(n_mixFr,35.0d0)
  
  K_thres = -5.07d-10*nre_lmtd**5 &
       &    +1.73d-7* nre_lmtd**4 &
       &    -2.17d-5* nre_lmtd**3 &
       &    +0.0013d0*nre_lmtd**2 &
       &    -0.037d0* nre_lmtd    &
       &    +0.8355

  if     (nre_lmtd <=  10.0d0)then
     K_crit = 1.250d0*nre_lmtd**(-0.350d0)
  else if(nre_lmtd <=  40.0d0)then
     K_crit = 1.072d0*nre_lmtd**(-0.301d0)
  else if(nre_lmtd <= 120.0d0)then
     K_crit = 0.356d0*nre_lmtd**(-0.003d0)
  end if

  if(n_mixFr_lmtd <= 0.35d0) then
     E_rime = 0.787d0*n_mixFr_lmtd**(0.988d0)*max(0.263d0*log(nre_lmtd)-0.264d0, 0.0d0)
  else if(n_mixFr_lmtd <= K_thres) then
     E_rime = (0.7475d0*log10(n_mixFr_lmtd)+0.620d0)*max(0.263d0*log(nre_lmtd)-0.264d0, 0.0d0)
!!$  else 
!!$     E_rime = sqrt( max(1.0d0-(log10(n_mixFr_lmtd/K_crit)-sqrt(5.0d0))**2/5.0d0, 0.0d0) )
!!$  end if
  else
     E_rime1 = (0.7475d0*log10(K_thres)+0.620d0)*max(0.263d0*log(nre_lmtd)-0.264d0, 0.0d0)
     E_rime2 = sqrt( max(1.0d0-(log10(n_mixFr_lmtd/K_crit)-sqrt(5.0d0))**2/5.0d0, 0.0d0) )
     E_rime = max(E_rime1,E_rime2)
  end if
  
end subroutine riming_efficiency_Erfani_Mitchell_2017_plate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine riming_efficiency_Erfani_Mitchell_2017_column(E_rime,nre,n_mixFr)
  implicit none
  
  real(DP),intent(out) :: E_rime

  real(DP),intent(in) :: nre       ! Reynolds number (based on column width, not length)
  real(DP),intent(in) :: n_mixFr   ! mixed Froude number

  real(DP) :: nre_lmtd ! min(nre,20.0d0)
  real(DP) :: n_mixFr_lmtd ! min(nre,20.0d0)
  real(DP) :: K_thres, K_crit
  real(DP) :: para_r
  
  if(nre<0.2d0) then
     E_rime = 0.0d0
     return
  end if
  
  nre_lmtd = min(nre,20.0d0)
  n_mixFr_lmtd = min(n_mixFr,20.0d0)

  if(nre_lmtd <= 2.0d0)then
     K_thres = 0.0251d0*nre_lmtd**2 &
          &   -0.0144d0*nre_lmtd    &
          &   +0.811d0
  else
     K_thres =-0.0003d0*nre_lmtd**3 &
          &   +0.0124d0*nre_lmtd**2 &
          &   -0.1634d0*nre_lmtd    &
          &   +1.0075d0
  end if

  if(nre_lmtd <= 1.7d0)then ! eq(19) of the paper is opposite
     para_r = 0.7422d0*nre_lmtd**(0.2111d0)
  else
     para_r = 0.8025d0*nre_lmtd**(0.0604d0)
  end if

  if(nre_lmtd <= 1.7d0)then
     K_crit = 0.7797d0*nre_lmtd**(-0.009d0)
  else
     K_crit = 1.0916d0*nre_lmtd**(-0.635d0)
  end if

  if(n_mixFr_lmtd <= K_thres) then
     if(nre_lmtd <= 3.0d0)then
        E_rime = 0.787d0*n_mixFr_lmtd**(0.988d0) &
             &  *(-0.0121d0*nre_lmtd**2 + 0.1297d0*nre_lmtd + 0.0598d0)
     else
        E_rime = 0.787d0*n_mixFr_lmtd**(0.988d0) &
             &  *(-0.0005d0*nre_lmtd**2 + 0.1028d0*nre_lmtd + 0.0359d0)
     end if
  else
     E_rime = para_r*sqrt( max(1.0d0-(log10(n_mixFr_lmtd/K_crit)-sqrt(3.5d0))**2/3.5d0, 0.0d0) )
  end if
  
end subroutine riming_efficiency_Erfani_Mitchell_2017_column
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sdm_rho_riming_Heymsfield_Pflaum_1985(rho_rime,sd_r,sd_re,sd_rp,sd_vz_drop,sd_vz_ice,sd_t,sd_rd,sd_p,sd_qv)
  use scale_const, only:  &
       rhow_mks => CONST_DWATR, &  ! density of liquid water [kg/m^3]
       t0 => CONST_TEM00           ! 0 degC in [K]
  use m_sdm_common, only: &
       O_THRD

  real(RP), intent(out) :: rho_rime   ! riming density [kg/m^3]
  real(RP), intent(in)  :: sd_r       ! droplet radius [m]
  real(RP), intent(in)  :: sd_re      ! equatorial radius of the ice [m]
  real(RP), intent(in)  :: sd_rp      ! polar radius of the ice [m]
  real(RP), intent(in)  :: sd_vz_drop ! terminal velocity of the droplet [m/s]
  real(RP), intent(in)  :: sd_vz_ice  ! terminal velocity of the ice [m/s]
  real(RP), intent(in)  :: sd_t       ! temperature[K] at the location of super-droplets
  real(RP), intent(in)  :: sd_rd      ! moist air density [kg/m3] at the location of super-droplets
  real(RP), intent(in)  :: sd_p       ! pressure of the grid contained the SD [Pa]
  real(RP), intent(in)  :: sd_qv      ! water-vapor of the grid contained the SD [kg/kg]

  real(RP) :: drop_radi_um   ! droplet radius [um]
  real(RP) :: impact_vel  ! impact velocity [m/s]
  real(RP) :: temp_icesfc ! surface temperature of ice [degC]
  real(RP) :: ratio_impact_vel ! V_impact/V_infinity

  real(RP) :: tdeg        ! temperaure in [degC]
  real(RP) :: sd_dvisc    ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
  real(RP) :: sd_maxD_ice ! maximum dimension of the ice [m]
  real(RP) :: sd_eqr_ice  ! equivalent radius of ice [m]
  real(RP) :: var_Y
  real(RP) :: n_re   ! Reynolds number
  real(RP) :: n_st   ! Stokes number

  !################## droplet radius in [um]
  drop_radi_um = sd_r * 1.0d6
  
  !################## impact velocity [m/s]
  ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
  !== (Pruppacher & Klett,1997) ==!
  tdeg = sd_t - t0     !! [K] => [degC]
  if( tdeg>=0.0d0 ) then
     sd_dvisc = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-5_RP
  else
     sd_dvisc = ( 1.718_RP + 4.9E-3_RP*tdeg              &
          - 1.2E-5_RP*tdeg*tdeg ) * 1.E-5_RP
  end if

  !###### Raynolds number
  sd_maxD_ice = 2.0d0 * max(sd_re,sd_rp) ! maximum dimension of the ice particle [m]
  n_re = sd_rd * sd_maxD_ice * sd_vz_ice / sd_dvisc

  !###### Stokes number
  sd_eqr_ice = (sd_re**2 * sd_rp)**(O_THRD) ! equivalent radius of ice [m]
  n_st = 2.0d0*sd_vz_ice*sd_r**2*rhow_mks/(9.0d0*sd_dvisc*sd_eqr_ice)

  !###### impact velocity [m/s]
  call sdm_impact_velocity_Rasummsen_Heymsfield_1985(ratio_impact_vel,n_re,n_st)
  impact_vel = abs( sd_vz_ice - sd_vz_drop ) * ratio_impact_vel

  !################## surface temperature of ice [degC]
  call sdm_ice_surface_temperature(temp_icesfc,sd_t,sd_p,sd_qv) ! surface temp in [K]
  temp_icesfc = temp_icesfc - t0     !! [K] => [degC]
  temp_icesfc = min(temp_icesfc,-0.01d0) ! limiter not to exceed 0 degC

  !################## evaluate the Y variable
  var_Y = -drop_radi_um*impact_vel/temp_icesfc

  !################## rho_rime in [g/cm^3]
  if( (temp_icesfc <= -5.0d0) .or. (var_Y < 1.6d0) )then
     rho_rime = 0.30d0*var_Y**0.44d0
  else
     rho_rime = exp( -0.03115d0 - 1.7030d0*var_Y + 0.9116d0*var_Y**2 - 0.1224*var_Y**3 )
  end if
  rho_rime = max(rho_rime,0.1d0)
  rho_rime = min(rho_rime,0.91d0)

  !################## rho_rime in [kg/m^3]
  rho_rime = rho_rime*1.0d3 

  return
end subroutine sdm_rho_riming_Heymsfield_Pflaum_1985
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sdm_impact_velocity_Rasummsen_Heymsfield_1985(ratio_impact_vel,n_re,n_st)
  real(RP), intent(out) :: ratio_impact_vel ! V_impact/V_infinity
  real(RP), intent(in) :: n_re   ! Reynolds number
  real(RP), intent(in) :: n_st   ! Stokes number

  real(RP) :: var_w ! log10(n_st)

  var_w = log10(n_st)

  if     (n_re<( 10.0d0+ 30.0d0)*0.5d0)then
     if  (n_st< 0.4d0)then
        ratio_impact_vel = 0.0d0
     else if(n_st<10.0d0)then
        ratio_impact_vel = 0.1701d0 + 0.7246d0*var_w + 0.2257d0*var_w**2 - 1.13d0*var_w**3 + 0.5756d0*var_w**4
     else
        ratio_impact_vel = 0.57d0
     end if
  else if(n_re<( 30.0d0+100.0d0)*0.5d0)then
     if  (n_st< 0.1d0)then
        ratio_impact_vel = 0.0d0
     else if(n_st<10.0d0)then
        ratio_impact_vel = 0.2927d0 + 0.5085d0*var_w - 0.03453d0*var_w**2 - 0.2184d0*var_w**3 + 0.03595d0*var_w**4
     else
        ratio_impact_vel = 0.59d0
     end if
  else if(n_re<(100.0d0+300.0d0)*0.5d0)then
     if  (n_st< 0.1d0)then
        ratio_impact_vel = 0.0d0
     else if(n_st<10.0d0)then
        ratio_impact_vel = 0.3272d0 + 0.4907d0*var_w - 0.09452d0*var_w**2 - 0.1906d0*var_w**3 + 0.07105d0*var_w**4
     else
        ratio_impact_vel = 0.61d0
     end if
  else
     if  (n_st< 0.1d0)then
        ratio_impact_vel = 0.0d0
     else if(n_st<10.0d0)then
        ratio_impact_vel = 0.356d0 + 0.4738d0*var_w - 0.1233d0*var_w**2 - 0.1618d0*var_w**3 + 0.08087d0*var_w**4
     else
        ratio_impact_vel = 0.63d0
     end if
  end if

  ratio_impact_vel = max(ratio_impact_vel,0.0d0)

  return
end subroutine sdm_impact_velocity_Rasummsen_Heymsfield_1985
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sdm_ice_surface_temperature(temp_icesfc,t_sd,p_sd,qv_sd)
  use scale_atmos_saturation, only: &
       ATMOS_SATURATION_pres2qsat_ice, &
       ATMOS_SATURATION_psat_ice
  use m_sdm_common, only: &
!!$       F_THRD, ONE_PI, &
!!$       VALID2INVALID, STAT_ICE, i2, sdicedef, &
       Diff_C, & ! diffusion constant of vapor [m^2/s]  
       GasV_C, & ! Gas Constant of vapor [J/K/kg]
       Heat_C, & ! thermal conductivity at 293K, 100kPa [J/m s K]
       LatHet_S  ! latent heat of sublimation at 0C [J/kg]

  real(RP), intent(out) :: temp_icesfc ! surface temperature of ice assuming the quasi-equilibrium of diffusion [K]
  real(RP), intent(in)  :: t_sd        ! anvient temperature [K]
  real(RP), intent(in)  :: p_sd        ! pressure of the grid contained the SD [Pa]
  real(RP), intent(in)  :: qv_sd       ! water-vapor of the grid contained the SD [kg/kg]

  real(RP) :: qvsi_sd   ! Saturation mixing ratio over ice of the grid contained the SD [kg/kg]  
  real(RP) :: esi_sd    ! Saturation vapor pressure over ice of the grid contained the SD [kg/kg]
  real(RP) :: satri_sd  ! Degree of super-saturation over ice of the grid contained the SD 
  real(RP) :: ssi       ! super saturation
  real(RP) :: Fac_ice_dd_invrhoi ! Fd/rhoi
  real(RP) :: Fac_ice_kk_invrhoi ! Fk/rhoi
  real(RP) :: ivt_sd    ! 1.d0 / t_sd
  real(RP) :: delta_rho   ! delta_rho = rho_vapor_ambient - rho_vapor_surface [kg/m^3] 
  real(RP) :: delta_temp  ! delta_temp = temp_ambient - temp_surface [K] 

  !!$          !! max of ssi ( = super saturation over ice at water saturation )
!!$          call ATMOS_SATURATION_pres2qsat_liq( qv_sd_max,t_sd,p_sd )
!!$          satri_sd_max  = qv_sd_max / qvsi_sd
!!$          ssi_max  = satri_sd_max - 1.0_RP
!!$          !! excess vapor density [g/m^3]: rhov_infty - rhov_surface, but limited by water saturation

  !###### ice saturation ratio ######!
  call ATMOS_SATURATION_pres2qsat_ice( qvsi_sd,t_sd,p_sd )
  satri_sd  = qv_sd / qvsi_sd
  ssi  = satri_sd - 1.0_RP

  !###### kinetic correction ######!
  ivt_sd = 1.0_RP / t_sd
  Fac_ice_kk_invrhoi = ( ivt_sd * LatHet_S / GasV_C - 1.0_RP ) * (LatHet_S/Heat_C) * ivt_sd

  !###### psychrometric correction ######!
  call ATMOS_SATURATION_psat_ice( esi_sd,t_sd )
  Fac_ice_dd_invrhoi = GasV_C / Diff_C * t_sd / esi_sd

  !###### delta_rho = rho_vapor_ambient - rho_vapor_surface [kg/m^3] 
  delta_rho  = ssi / ( Fac_ice_dd_invrhoi + Fac_ice_kk_invrhoi ) / Diff_C

  !###### delta_temp = temp_ambient - temp_surface [K] 
  delta_temp = (LatHet_S*Diff_C/Heat_C)*delta_rho

  !###### surface temperature [K]
  temp_icesfc = t_sd + delta_temp

  return
end subroutine sdm_ice_surface_temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_sdm_coalescence_cold
