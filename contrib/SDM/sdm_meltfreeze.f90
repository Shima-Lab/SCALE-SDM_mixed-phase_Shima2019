!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Melting and freezing of water particles
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
!! @li      2016-07-15 (S.Shima) [new] sdm_meltfreeze created.
!! @li      2016-07-15 (S.Shima) [add] sdm_meltfreeze_updatefluid added
!! @li      2017-02-09 (S.Shima) [mod] conversion between spheroid ice <-> sphere drop 
!! @li      2017-02-15 (S.Shima) [mod] to preserve mass of water when a spheroid ice melt to a droplet 
!! @li      2017-02-19 (S.Shima) [mod] to avoid confusion, reset sdi%rho(n) = 0.0d0 when it melts
!! @li      2017-08-24 (S.Shima) [mod] require Sw>=1 for immersion/deposition/homogeneous freezing 
!! @li      2018-06-30 (S.Shima) [add] calculation of rime mass and number of monomers
!<
!-------------------------------------------------------------------------------
module m_sdm_meltfreeze
  use scale_precision

  implicit none
  private
  public :: sdm_meltfreeze, sdm_meltfreeze_updatefluid

contains
  subroutine sdm_meltfreeze(           &
       t_scale,pres_scale,qv_scale,    &
       sd_num,sd_liqice,sd_x,sd_y,     &
       sd_r,sd_ri,sd_rj,sd_rk,sdi      )
    
    use scale_const, only: &
         t0   => CONST_TEM00, &      ! 0 degC in [K] 
         dens_i_mks => CONST_DICE, &   ! density of ice [kg/m^3]
         dens_w_mks => CONST_DWATR     ! density of liquid water [kg/m^3]
    use scale_grid_index, only: &
         IA,JA,KA
    use scale_atmos_saturation, only: &
        ATMOS_SATURATION_pres2qsat_liq
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         VALID2INVALID, DNS_RL, i2, sdicedef, STAT_LIQ, STAT_ICE

    ! Input variables
    real(RP), intent(in) :: t_scale(KA,IA,JA) ! Temperature [K]
    real(RP), intent(in) :: pres_scale(KA,IA,JA) ! Pressure [Pa]
    real(RP), intent(in) :: qv_scale(KA,IA,JA)   ! Water vapor mixing ratio [kg/kg]
    integer,  intent(in) :: sd_num      ! number of super-droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)   ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)   ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)   ! index[k/real] of super-droplets
    ! Input and output variables
    integer(i2), intent(inout) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(inout) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    type(sdicedef), intent(inout) :: sdi   ! ice phase super-droplets
    ! Work variables
    real(RP) :: p_sd      ! pressure of the grid contained the SD    [Pa]
    real(RP) :: t_sd      ! temperature of the grid contained the SD [K]
    real(RP) :: qv_sd     ! water-vapor of the grid contained the SD [kg/kg]
    real(RP) :: qvs_sd    ! Saturation mixing ratio of the grid contained the SD. [kg/kg]
    real(RP) :: rwori_onethird = (dens_w_mks/dens_i_mks)**(1.0d0/3.0d0)

    integer :: i, j, k, n    ! index
    !---------------------------------------------------------------------
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_meltfreeze",1,1)
#endif
    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

    

    ! Melting and freezing of water
!OCL INDEPENDENT
    do n=1,sd_num

       !### Skip invalid super-droplets ###!

       if( sd_rk(n)<VALID2INVALID ) cycle
       
       !### Get the location and variables of Super-Droplets ###!
       
       i = floor(sd_ri(n))+1
       j = floor(sd_rj(n))+1
       k = floor(sd_rk(n))+1

       !== temperarure [K] ==!

       t_sd  = t_scale(k,i,j)

       !== pressure [Pa] ==!

       p_sd  = pres_scale(k,i,j)

       !== mixing ratio [kg/kg] ==!

       qv_sd = qv_scale(k,i,j)

       !== saturation mixing ratio over water [kg/kg] ==!

       call ATMOS_SATURATION_pres2qsat_liq( qvs_sd,t_sd,p_sd )

       !### freezing ###!
       
       !###### immersion/condensation/homogeneous freezing ######!
       if( (sd_liqice(n)==STAT_LIQ) .and. (t_sd<(sdi%tf(n)+t0)) .and. (qv_sd>=qvs_sd) ) then

          sd_liqice(n) = STAT_ICE
          sdi%re(n)    = sd_r(n) * rwori_onethird
          sdi%rp(n)    = sd_r(n) * rwori_onethird
          sdi%rho(n)   = dens_i_mks
          sdi%mrime(n) = 0.0d0
          sdi%nmono(n) = 1
          sd_r(n)      = 0.0d0

       end if

       !### melting ###!
         
       if( (sd_liqice(n)==STAT_ICE) .and. (t_sd>t0) ) then

          sd_liqice(n) = STAT_LIQ
          sd_r(n)      = (sdi%re(n)**2 * sdi%rp(n) * sdi%rho(n) / dens_w_mks)**(1.0d0/3.0d0)
          sdi%rho(n)   = 0.0d0
          sdi%re(n)    = 0.0d0
          sdi%rp(n)    = 0.0d0
          sdi%mrime(n) = 0.0d0
          sdi%nmono(n) = 0
          
       end if

    end do

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_meltfreeze",1,1)
#endif
    return
  end subroutine sdm_meltfreeze
!-----------------------------------------------------------------------------
  subroutine sdm_meltfreeze_updatefluid(RHOT,QTRC,DENS,rhosolp_sdm,rhosolc_sdm)
    use scale_tracer, only: &
         I_QV,QAD=>QA
    use scale_atmos_thermodyn, only: &
        ATMOS_THERMODYN_templhf
    use scale_grid_index, only: &
         IA,JA,KA,IS,IE,JS,JE,KS,KE
    use m_sdm_fluidconv, only: &
         sdm_rhot_qtrc2cpexnr, &
         sdm_rhot_qtrc2p_t
    use m_sdm_common, only: &
         LatHet
    ! Input variables
    real(RP), intent(in) :: rhosolp_sdm(KA,IA,JA) ! density of solid water at preveous step
    real(RP), intent(in) :: rhosolc_sdm(KA,IA,JA) ! density of solid water at current step
    real(RP), intent(in) :: QTRC(KA,IA,JA,QAD)    !! Ratio of mass of tracer to total mass[kg/kg]
    real(RP), intent(in) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    ! Output variables
    real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    ! Work variables
    real(RP) :: cpexnr(KA,IA,JA)
    real(RP) :: delta_rhosol
    integer  :: i, j, k ! index
    real(RP) :: lhf(KA,IA,JA), pre(KA,IA,JA), temp(KA,IA,JA)
    !-------------------------------------------------------------------7--
     
    ! no exchange vapor

    ! exchange heat
    call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pre,temp)
    call ATMOS_THERMODYN_templhf(lhf,temp)

    call sdm_rhot_qtrc2cpexnr(RHOT,QTRC,DENS,cpexnr)

    do k=KS,KE
    do j=JS,JE
    do i=IS,IE
       delta_rhosol = rhosolc_sdm(k,i,j)-rhosolp_sdm(k,i,j)

       RHOT(k,i,j) = RHOT(k,i,j) + lhf(k,i,j)*delta_rhosol / cpexnr(k,i,j)  
    end do
    end do
    end do

    return
  end subroutine sdm_meltfreeze_updatefluid
end module m_sdm_meltfreeze
