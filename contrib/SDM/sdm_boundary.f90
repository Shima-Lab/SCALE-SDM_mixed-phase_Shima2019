!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Boundary condition management module
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
!! @li      2014-07-11 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2014-07-18 (Y.Sato)  [mod] Modify the output message in sdm_boundary
!! @li      2014-07-18 (Y.Sato)  [mod] Modify a bug in sdm_getbufsy
!! @li      2015-09-08 (Y.Sato)  [mod] update for version SCALE 0.2.4
!! @li      2016-07-13 (S.Shima) [mod] Modified for MPI communication of ice particles
!! @li      2018-05-13 (S.Shima) [fix] a bug in sdm_getbufsy
!! @li      2018-06-26 (S.Shima) [fix] a bug in sdm_boundary
!! @li      2018-06-30 (S.Shima) [add] rime mass and number of monomers as SD attributes
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_boundary
  use scale_precision

  implicit none
  private
  public :: sdm_jdginvdv, sdm_boundary

contains
  subroutine sdm_jdginvdv(sd_rkl,sd_rku,sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk)
    use scale_grid_index, only: &
         IA,JA
    use scale_grid, only: &
         FX => GRID_FX, &
         CZ => GRID_CX, &
         FY => GRID_FY
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         VALID2INVALID, INVALID, PRECIPI, rkumin, rklmax

    ! Input variables
    real(RP), intent(in) :: sd_rkl(IA,JA)  ! index[k/real] at 'sdm_zlower'
    real(RP), intent(in) :: sd_rku(IA,JA)  ! index[k/real] at 'sdm_zupper
    integer, intent(in) :: sd_num          ! number of super-droplets
    ! Input and output variables
    real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
    real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
    real(RP), intent(out) :: sd_ri(1:sd_num)   ! index-i(real) of super-droplets
    real(RP), intent(out) :: sd_rj(1:sd_num)   ! index-j(real) of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    ! Work variables
    
    real(RP) :: ri    ! real index [i] of super-droplets
    real(RP) :: rj    ! real index [j] of super-droplets
    real(RP) :: rkl   ! interpolated 'sd_rkl'
    real(RP) :: rku   ! interpolated 'sd_rku'
    
    real(RP) :: sXm   ! variable for inteporation
    real(RP) :: sXp   ! variable for inteporation
    real(RP) :: sYm   ! variable for inteporation
    real(RP) :: sYp   ! variable for inteporation
    
    integer :: iXm         ! index for inteporation
    integer :: iXp         ! index for inteporation
    integer :: iYm         ! index for inteporation
    integer :: iYp         ! index for inteporation
    
    integer :: n, i, j     ! index
    !---------------------------------------------------------------------
    
    !### get horizontal face index(real) of super-droplets ###!
    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)
    
    ! Judge super-droplets as invalid one or valid one in vertical
    do n=1,sd_num
       
       !### skip invalid super-droplets ###!
       
       if( sd_rk(n)<VALID2INVALID ) cycle
       
       !### convert invalid by outflowing from upper boundary ###!
       
       if( sd_rk(n)>=rkumin ) then
          
          ri = sd_ri(n)
          rj = sd_rj(n)
          
          !! Interpolation in certer grid
          iXm = floor(ri+0.5d0) !! real index in face grid to center grid
          iXp = iXm + 1
          sXm = (ri+0.5d0) - real(iXm,kind=RP)
          sXp = 1.0_RP - sXm

          iYm = floor(rj+0.5d0)
          iYp = iYm + 1
          sYm = (rj+0.5d0) - real(iYm,kind=RP)
          sYp = 1.0_RP - sYm

          rku =  sd_rku(iXm,iYm) * ( SXp * SYp )                       &
               + sd_rku(iXp,iYm) * ( SXm * SYp )                       &
               + sd_rku(iXm,iYp) * ( SXp * SYm )                       &
               + sd_rku(iXp,iYp) * ( SXm * SYm )

          if( sd_rk(n)>=rku ) then
             sd_rk(n) = INVALID
          end if

       end if

       !### convert precipitation by outflowing ###!
       !### from lower boundary                 ###!

       if( sd_rk(n)<=rklmax ) then
            
          ri = sd_ri(n)
          rj = sd_rj(n)

          !! Interpolation in certer grid
          iXm = floor(ri+0.5d0) !! real index in face grid to center grid
          iXp = iXm + 1
          sXm = (ri+0.5d0) - real(iXm,kind=RP)
          sXp = 1.0_RP - sXm

          iYm = floor(rj+0.5d0)
          iYp = iYm + 1
          sYm = (rj+0.5d0) - real(iYm,kind=RP)
          sYp = 1.0_RP - sYm
          
          rkl =  sd_rkl(iXm,iYm) * ( SXp * SYp )                       &
               + sd_rkl(iXp,iYm) * ( SXm * SYp )                       &
               + sd_rkl(iXm,iYp) * ( SXp * SYm )                       &
               + sd_rkl(iXp,iYp) * ( SXm * SYm )

          if( sd_rk(n)<rkl ) then
             sd_rk(n) = PRECIPI
          end if
          
       end if

    end do
    return
  end subroutine sdm_jdginvdv
  !----------------------------------------------------------------------------
  subroutine sdm_boundary(wbc,ebc,sbc,nbc,                         &
                         sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_rk,   &
                         sd_u,sd_v,sd_vz,sd_r,sd_asl,sdi,           &
                         bufsiz1,                               &
                         bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4,      &
                         sd_itmp1,                              &
                         rbuf_r8,sbuf_r8,rbuf_i8,sbuf_i8,rbuf_i2,sbuf_i2,rbuf_i4,sbuf_i4) 
    use scale_process, only: &
         PRC_MPIstop, &
         mype => PRC_myrank
    use scale_grid, only: &
         GRID_FX, &
         GRID_FY
    use scale_stdio, only: &
         IO_L, IO_FID_LOG
    use scale_grid_index, only: &
         IS,IE,JS,JE
    use m_sdm_common, only: &
         stat,nisub,njsub, &
         VALID2INVALID,i2,sdicedef

    ! Input variables
    integer, intent(in) :: wbc ! Option for west boundary conditions
    integer, intent(in) :: ebc ! Option for east boundary conditions
    integer, intent(in) :: sbc ! Option for south boundary conditions
    integer, intent(in) :: nbc ! Option for north boundary conditions
    integer, intent(in) :: sd_num ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of kind of chemical materia contained as water-soluble aerosol in super droplets
    integer, intent(in) :: bufsiz1  ! buffer size for MPI
    integer, intent(in) :: bufsiz2_r8 ! buffer size for MPI (real8)
    integer, intent(in) :: bufsiz2_i8 ! buffer size for MPI (int8)
    integer, intent(in) :: bufsiz2_i2 ! buffer size for MPI (int2)
    integer, intent(in) :: bufsiz2_i4 ! buffer size for MPI (int4)
    ! Input and output variables
    integer(DP), intent(inout) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
    integer(i2), intent(inout) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
    real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets in vertical
    real(RP), intent(inout) :: sd_u(1:sd_num)  ! x-components velocity of super-droplets
    real(RP), intent(inout) :: sd_v(1:sd_num)  ! y-components velocity of super-droplets
    real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    real(RP), intent(inout) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    type(sdicedef), intent(inout) :: sdi   ! ice phase super-droplets

    real(DP), intent(inout) ::                                   &
         &                             rbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Receiving buffer (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf,mrime
                       ! dim03 = 1:west, 2:east / 1:south, 2:north
    real(DP), intent(inout) ::                                   &
         &                             sbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf,mrime
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
    integer, intent(out) :: sd_itmp1(1:sd_num) ! temporary array of the size of the number of super-droplets.

    integer :: n     ! index
    !---------------------------------------------------------------------
    
    if( (wbc/=1).or.(ebc/=1).or.(sbc/=1).or.(nbc/=1))then
       if( IO_L ) write(IO_FID_LOG,*) 'ERROR: Only periodic B.C. is supported!'
       if( IO_L ) write(IO_FID_LOG,*) 'ERROR: Set PRC_PERIODIC_X=PRC_PERIODIC_Y=.true.'
       call PRC_MPIstop
    end if

    ! Initialize
    stat = 0
    
    ! Set the boundary condition of super-droplets in x direction
    if( nisub>=2 ) then
       !== Exchange the value horizontally in x direction ==!
       
       call sdm_putbufsx(wbc,ebc,sd_num,sd_numasl,              &
            sd_n,sd_liqice,sd_x,sd_y,sd_rk,sd_u,sd_v,sd_vz,    &
            sd_r,sd_asl,sdi,                               &
            bufsiz1,bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4, &
            stat,sd_itmp1,                            &
            sbuf_r8,sbuf_i8,sbuf_i2,sbuf_i4)

       ! In case of exsiting outflow super-droplets
       call sdm_shiftsx(wbc,ebc,bufsiz1, &
            bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4, &
            sbuf_r8,rbuf_r8,sbuf_i8,rbuf_i8,sbuf_i2,rbuf_i2,sbuf_i4,rbuf_i4)

       call sdm_getbufsx(wbc,ebc,sd_num,sd_numasl,              &
            sd_n,sd_liqice,sd_x,sd_y,sd_rk,sd_u,sd_v,sd_vz,    &
            sd_r,sd_asl,sdi,                               &
            bufsiz1,bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4, &
            stat,sd_itmp1,                            &
            rbuf_r8,rbuf_i8,rbuf_i2,rbuf_i4)

    else if((wbc==1).and.(ebc==1).and.(nisub==1))then
       !### Apply periodic boundary condition
       do n=1,sd_num
          if( sd_rk(n)>VALID2INVALID ) then
             if( sd_x(n)>=GRID_FX(IE) )then
                sd_x(n) = sd_x(n)-GRID_FX(IE)+GRID_FX(IS-1)
             else if( sd_x(n)<GRID_FX(IS-1) )then
                sd_x(n) = sd_x(n)-GRID_FX(IS-1)+GRID_FX(IE)
             end if
          end if
       end do
    end if

    ! Set the boundary condition of super-droplets in y direction
    if( njsub>=2 ) then
       !== Exchange the value horizontally in y direction ==!

       call sdm_putbufsy(sbc,nbc,sd_num,sd_numasl,              &
            sd_n,sd_liqice,sd_x,sd_y,sd_rk,sd_u,sd_v,sd_vz,    &
            sd_r,sd_asl,sdi,                               &
            bufsiz1,bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4, &
            stat,sd_itmp1,                            &
            sbuf_r8,sbuf_i8,sbuf_i2,sbuf_i4)

       ! In case of exsiting outflow super-droplets
       call sdm_shiftsy(sbc,nbc,bufsiz1,      &
            bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4, &
            sbuf_r8,rbuf_r8,sbuf_i8,rbuf_i8,sbuf_i2,rbuf_i2,sbuf_i4,rbuf_i4)
       
       call sdm_getbufsy(sbc,nbc,sd_num,sd_numasl,              &
            sd_n,sd_liqice,sd_x,sd_y,sd_rk,sd_u,sd_v,sd_vz,    &
            sd_r,sd_asl,sdi,                               &
            bufsiz1,bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4, &
            stat,sd_itmp1,                            &
            rbuf_r8,rbuf_i8,rbuf_i2,rbuf_i4)

    else if((nbc==1).and.(sbc==1).and.(njsub==1))then
       !### Apply periodic boundary condition
       do n=1,sd_num
          if( sd_rk(n)>VALID2INVALID ) then
             if( sd_y(n)>=GRID_FY(JE) )then
                sd_y(n) = sd_y(n)-GRID_FY(JE)+GRID_FY(JS-1)
             else if( sd_y(n)<GRID_FY(IS-1) )then
                sd_y(n) = sd_y(n)-GRID_FY(JS-1)+GRID_FY(JE)
             end if
          end if
       end do
    end if

    ! Check send/receive error of super-droplets
    if( stat /= 0 ) then
       ! Does this work?? Can mype>=1 write a message to the LOG file?
       write(*,*)"  ### [SDM] : send/receive error of super-droplets at MPI rank=", &
            &         mype," ###"
       call PRC_MPIstop
    end if

    return
  end subroutine sdm_boundary
  !----------------------------------------------------------------------------
  subroutine sdm_putbufsx(wbc,ebc,sd_num,sd_numasl,         &
       sd_n,sd_liqice,sd_x,sd_y,sd_rk,sd_u,sd_v,sd_vz, &
       sd_r,sd_asl,sdi,                           &
       bufsiz1,bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4, &
       stat,ilist,                   &
       sbuf_r8,sbuf_i8,sbuf_i2,sbuf_i4)
    use scale_grid, only: &
         GRID_FX
    use scale_grid_index, only: &
         IS,IE
    use m_sdm_common, only: &
         nisub,INVALID,INVALID_i8,INVALID_i2,VALID2INVALID,i2,sdm_cold,sdicedef

    ! Input variables
    integer, intent(in) :: wbc    ! Option for west boundary conditions
    integer, intent(in) :: ebc    ! Option for east boundaty conditions
    integer, intent(in) :: sd_num ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of kind of chemical material contained as water-soluble aerosol in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_y(1:sd_num)    ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_u(1:sd_num) ! x-components velocity of super-droplets
    real(RP), intent(inout) :: sd_v(1:sd_num) ! y-components velocity of super-droplets
    real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplet
    integer, intent(in) :: bufsiz1  ! buffer size for MPI
    integer, intent(in) :: bufsiz2_r8 ! buffer size for MPI (real8)
    integer, intent(in) :: bufsiz2_i8 ! buffer size for MPI (int8)
    integer, intent(in) :: bufsiz2_i2 ! buffer size for MPI (int2)
    integer, intent(in) :: bufsiz2_i4 ! buffer size for MPI (int4)
    ! Input and output variable
    integer, intent(inout) :: stat  ! Runtime status
    real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    ! Output variable
    real(DP), intent(out) :: sbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer in x direction (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf,mrime
                       ! dim03 = 1 : west, 2: east
    integer(DP), intent(out) ::                                  &
     &                              sbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! Sending buffer in x direction (int8)
                       ! dim02 = 1 : multiplicity 
                       ! dim03 = 1 : west, 2: east
    integer(i2), intent(out) ::                                  &
     &                              sbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! Sending buffer in x direction (int2)
                       ! dim02 = 1 : status (liq/ice)
                       ! dim03 = 1 : west, 2: east
    integer, intent(out) ::                                  &
     &                              sbuf_i4(1:,1:,1:)
                       ! sbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or sbuf_i4(1:1,1:1,1:2) 
                       ! Sending buffer in x direction (int4)
                       ! dim02 = 1 : nmono (ice)
                       ! dim03 = 1 : west, 2: east
    integer, intent(out) :: ilist(1:sd_num)
    ! buffer for list vectorization

    ! Work variables
    integer :: tesend           ! total list number
    integer :: twsend           ! total list number
    integer :: ecnt             ! counter
    integer :: wcnt             ! counter
    integer :: m, n    ! index
    integer :: nasl    ! index
    integer :: nwsdm            ! number of warm SDM variables
    !---------------------------------------------------------------------

    ! Initialize
    tesend = 0
    twsend = 0

    nwsdm = 7 + sd_numasl  !! warm SDM variables
      
    sbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2) = INVALID
    sbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2) = INVALID_i8
    sbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2) = INVALID_i2

    ! Put the sending buffer in x direction.

    if( nisub>=2 ) then

       !### Fill in the sending buffer with the value ###!
       !### in the west halo regions.                 ###!

       !== periodic boundary conditions ==!
       !! get number and index of outflowed super-droplet

       wcnt   = 0
       do n=1,sd_num
          if( sd_x(n)<GRID_FX(IS-1) .and. sd_rk(n)>VALID2INVALID ) then
             wcnt = wcnt + 1
             ilist(wcnt) = n
          end if
       end do
       twsend  = wcnt

       !! check send buffer
       if( stat==0 .and. (twsend>bufsiz1) ) then
          stat = -1
       end if

       !! send data
       if( twsend>0 ) then
          do m=1,twsend
             n = ilist(m)

             sbuf_r8(m,1,1) = sd_x(n)-GRID_FX(IS-1)
             sbuf_r8(m,2,1) = sd_y(n)
             sbuf_r8(m,3,1) = sd_rk(n)
             sbuf_r8(m,4,1) = sd_u(n)
             sbuf_r8(m,5,1) = sd_v(n)
             sbuf_r8(m,6,1) = sd_vz(n)
             sbuf_r8(m,7,1) = sd_r(n)

             sbuf_i8(m,1,1) = sd_n(n)
             sbuf_i2(m,1,1) = sd_liqice(n)
          end do

          do nasl=1,sd_numasl
             do m=1,twsend
                n = ilist(m)

                sbuf_r8(m,7+nasl,1) = sd_asl(n,nasl)

             end do
          end do

          !!! cold-SDM only 
          if( sdm_cold ) then

             do m=1,twsend
                n = ilist(m)

                sbuf_r8(m,nwsdm+1,1) = sdi%re(n)
                sbuf_r8(m,nwsdm+2,1) = sdi%rp(n)
                sbuf_r8(m,nwsdm+3,1) = sdi%rho(n)
                sbuf_r8(m,nwsdm+4,1) = sdi%tf(n)
                sbuf_r8(m,nwsdm+5,1) = sdi%mrime(n)
                sbuf_i4(m,1,1)       = sdi%nmono(n)

             end do

          end if

          !== convert to invalid ==!
          do m=1,twsend
             n = ilist(m)

             sd_rk(n) = INVALID
          end do
       end if


       !### Fill in the sending buffer with the value ###!
       !### in the east halo regions                  ###!
       
       !== periodic boundary conditions ==!
       !! get number and index of outflowed super-droplet
       
       ecnt   = 0
       do n=1,sd_num
          if( sd_x(n)>=GRID_FX(IE) .and. sd_rk(n)>VALID2INVALID ) then
             ecnt = ecnt + 1
             ilist(ecnt) = n
          end if
       end do
       tesend  = ecnt

       !! check send buffer
       if( stat==0 .and. (tesend>bufsiz1) ) then
          stat = -1
       end if

       !! send data
       if( tesend>0 ) then
          do m=1,tesend
             n = ilist(m)

             sbuf_r8(m,1,2) = sd_x(n)-GRID_FX(IE)
             sbuf_r8(m,2,2) = sd_y(n)
             sbuf_r8(m,3,2) = sd_rk(n)
             sbuf_r8(m,4,2) = sd_u(n)
             sbuf_r8(m,5,2) = sd_v(n)
             sbuf_r8(m,6,2) = sd_vz(n)
             sbuf_r8(m,7,2) = sd_r(n)

             sbuf_i8(m,1,2) = sd_n(n)
             sbuf_i2(m,1,2) = sd_liqice(n)
          end do

          do nasl=1,sd_numasl
             do m=1,tesend
                n = ilist(m)

                sbuf_r8(m,7+nasl,2) = sd_asl(n,nasl)

             end do
          end do

          !!! cold-SDM only 
          if( sdm_cold ) then

             do m=1,tesend
                n = ilist(m)

                sbuf_r8(m,nwsdm+1,2) = sdi%re(n)
                sbuf_r8(m,nwsdm+2,2) = sdi%rp(n)
                sbuf_r8(m,nwsdm+3,2) = sdi%rho(n)
                sbuf_r8(m,nwsdm+4,2) = sdi%tf(n)
                sbuf_r8(m,nwsdm+5,2) = sdi%mrime(n)
                sbuf_i4(m,1,2)       = sdi%nmono(n)

             end do

          end if

          !== convert to invalid ==!
          do m=1,tesend
             n = ilist(m)

             sd_rk(n) = INVALID
          end do
       end if

    end if

    return
  end subroutine sdm_putbufsx
  !----------------------------------------------------------------------------
  subroutine sdm_shiftsx(wbc,ebc,           &
       bufsiz1,bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4,&
       sbuf_r8,rbuf_r8,sbuf_i8,rbuf_i8,sbuf_i2,rbuf_i2,sbuf_i4,rbuf_i4)
    use mpi
    use scale_process, only: &
         PRC_LOCAL_COMM_WORLD
    use m_sdm_common, only: &
         nisub, &
         dstw_sub, dste_sub, srcw_sub, srce_sub, &
         INVALID,INVALID_i8,INVALID_i2, &
         tag,i2,sdm_cold
    ! Input variables
    integer, intent(in) :: wbc      ! Option for west boundary conditions
    integer, intent(in) :: ebc      ! Option for east boundaty conditions
    integer, intent(in) :: bufsiz1  ! buffer size for MPI
    integer, intent(in) :: bufsiz2_r8 ! buffer size for MPI (real8)
    integer, intent(in) :: bufsiz2_i8 ! buffer size for MPI (int8)
    integer, intent(in) :: bufsiz2_i2 ! buffer size for MPI (int2)
    integer, intent(in) :: bufsiz2_i4 ! buffer size for MPI (int4)
    real(DP), intent(in)                                         &
     &                      :: sbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer in x direction (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf,mrime
                       ! dim03 = 1 : west, 2: east
    integer(DP), intent(in)                                      &
     &                      :: sbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! Sending buffer in x direction (int8)
                       ! dim02 = 1 : multiplicity
                       ! dim03 = 1 : west, 2: east
    integer(i2), intent(in)                                      &
     &                      :: sbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! Sending buffer in x direction (int2)
                       ! dim02 = 1 : status (liq/ice)
                       ! dim03 = 1 : west, 2: east
    integer, intent(in)                                      &
     &                      :: sbuf_i4(1:,1:,1:)
                       ! sbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or sbuf_i4(1:1,1:1,1:2) 
                       ! Sending buffer in x direction (int4)
                       ! dim02 = 1 : nmono (ice)
                       ! dim03 = 1 : west, 2: east
    ! Output variable
    real(DP), intent(out)                                        &
     &                      :: rbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer in x direction (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,t,tf
                       ! dim03 = 1 : west, 2: east
    integer(DP), intent(out)                                     &
     &                      :: rbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! Sending buffer in x direction (int8)
                       ! dim02 = 1 : multiplicity
                       ! dim03 = 1 : west, 2: east
    integer(i2), intent(out)                                     &
     &                      :: rbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! Sending buffer in x direction (int2)
                       ! dim02 = 1 : status (liq/ice)
                       ! dim03 = 1 : west, 2: east
    integer, intent(out)                                     &
     &                      :: rbuf_i4(1:,1:,1:)
                       ! rbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or rbuf_i4(1:1,1:1,1:2) 
                       ! Sending buffer in x direction (int4)
                       ! dim02 = 1 : nmono (ice)
                       ! dim03 = 1 : west, 2: east

    integer :: dstw     ! West sending distnation
    integer :: dste     ! East sending distnation
    integer :: srcw     ! West receiving source
    integer :: srce     ! East receiving source
    integer :: statsw   ! West sending request status
    integer :: statse   ! East sending request status
    integer :: statrw   ! West receiving request status
    integer :: statre   ! East receiving request status
    integer :: siz_r8   ! Sending and receiving buffer size
    integer :: siz_i8   ! Sending and receiving buffer size
    integer :: siz_i2   ! Sending and receiving buffer size
    integer :: siz_i4   ! Sending and receiving buffer size
    integer :: ierr     ! Error descriptor
    integer :: stat(mpi_status_size) ! Runtime status table
    !------------------------------------------------------------------

    ! Initialize
    rbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2) = INVALID
    rbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2) = INVALID_i8
    rbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2) = INVALID_i2

    ! Exchange the value in x direction.

    if( nisub>=2 ) then

       !### Set the buffre size ###!

       siz_r8 = bufsiz1 * bufsiz2_r8
       siz_i8 = bufsiz1 * bufsiz2_i8
       siz_i2 = bufsiz1 * bufsiz2_i2
       siz_i4 = 0
       if( sdm_cold ) then
          siz_i4 = bufsiz1 * bufsiz2_i4
       end if

       !### Set the processor element number for sending ###!
       !### distination and receiving source             ###!

       dstw = dstw_sub
       dste = dste_sub
       srcw = srcw_sub
       srce = srce_sub

       !### Send/Receive real8 ###!

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function (towards west)###!

       call mpi_isend(sbuf_r8(1,1,1),siz_r8,MPI_DOUBLE_PRECISION,dstw,tag, &
            PRC_LOCAL_COMM_WORLD,statsw,ierr)

       call mpi_irecv(rbuf_r8(1,1,1),siz_r8,MPI_DOUBLE_PRECISION,srce,tag, &
            PRC_LOCAL_COMM_WORLD,statre,ierr)

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function (towards east)###!

       call mpi_isend(sbuf_r8(1,1,2),siz_r8,MPI_DOUBLE_PRECISION,dste,tag, &
            PRC_LOCAL_COMM_WORLD,statse,ierr)

       call mpi_irecv(rbuf_r8(1,1,2),siz_r8,MPI_DOUBLE_PRECISION,srcw,tag, &
            PRC_LOCAL_COMM_WORLD,statrw,ierr)

       !### Call the waiting MPI function ###!

       call mpi_wait(statsw,stat,ierr)
       call mpi_wait(statse,stat,ierr)
       call mpi_wait(statrw,stat,ierr)
       call mpi_wait(statre,stat,ierr)

       !### Send/Receive int8 ###!

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function (towards west)###!

       call mpi_isend(sbuf_i8(1,1,1),siz_i8,MPI_LONG_LONG_INT,dstw,tag, &
            PRC_LOCAL_COMM_WORLD,statsw,ierr)

       call mpi_irecv(rbuf_i8(1,1,1),siz_i8,MPI_LONG_LONG_INT,srce,tag, &
            PRC_LOCAL_COMM_WORLD,statre,ierr)

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function (towards east)###!

       call mpi_isend(sbuf_i8(1,1,2),siz_i8,MPI_LONG_LONG_INT,dste,tag, &
            PRC_LOCAL_COMM_WORLD,statse,ierr)

       call mpi_irecv(rbuf_i8(1,1,2),siz_i8,MPI_LONG_LONG_INT,srcw,tag, &
            PRC_LOCAL_COMM_WORLD,statrw,ierr)

       !### Call the waiting MPI function ###!

       call mpi_wait(statsw,stat,ierr)
       call mpi_wait(statse,stat,ierr)
       call mpi_wait(statrw,stat,ierr)
       call mpi_wait(statre,stat,ierr)

       !### Send/Receive int2 ###!

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function (towards west)###!

       call mpi_isend(sbuf_i2(1,1,1),siz_i2,MPI_INTEGER2,dstw,tag, &
            PRC_LOCAL_COMM_WORLD,statsw,ierr)

       call mpi_irecv(rbuf_i2(1,1,1),siz_i2,MPI_INTEGER2,srce,tag, &
            PRC_LOCAL_COMM_WORLD,statre,ierr)

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function (towards east)###!

       call mpi_isend(sbuf_i2(1,1,2),siz_i2,MPI_INTEGER2,dste,tag, &
            PRC_LOCAL_COMM_WORLD,statse,ierr)

       call mpi_irecv(rbuf_i2(1,1,2),siz_i2,MPI_INTEGER2,srcw,tag, &
            PRC_LOCAL_COMM_WORLD,statrw,ierr)

       !### Call the waiting MPI function ###!

       call mpi_wait(statsw,stat,ierr)
       call mpi_wait(statse,stat,ierr)
       call mpi_wait(statrw,stat,ierr)
       call mpi_wait(statre,stat,ierr)

       !### Send/Receive int4 ###!

       if( sdm_cold ) then

          !### Incliment the message tag ###!

          tag = tag + 1

          !### Call the sending and receiving MPI function (towards west)###!

          call mpi_isend(sbuf_i4(1,1,1),siz_i4,MPI_INTEGER,dstw,tag, &
               PRC_LOCAL_COMM_WORLD,statsw,ierr)

          call mpi_irecv(rbuf_i4(1,1,1),siz_i4,MPI_INTEGER,srce,tag, &
               PRC_LOCAL_COMM_WORLD,statre,ierr)

          !### Incliment the message tag ###!

          tag = tag + 1

          !### Call the sending and receiving MPI function (towards east)###!

          call mpi_isend(sbuf_i4(1,1,2),siz_i4,MPI_INTEGER,dste,tag, &
               PRC_LOCAL_COMM_WORLD,statse,ierr)

          call mpi_irecv(rbuf_i4(1,1,2),siz_i4,MPI_INTEGER,srcw,tag, &
               PRC_LOCAL_COMM_WORLD,statrw,ierr)

          !### Call the waiting MPI function ###!

          call mpi_wait(statsw,stat,ierr)
          call mpi_wait(statse,stat,ierr)
          call mpi_wait(statrw,stat,ierr)
          call mpi_wait(statre,stat,ierr)

       end if

    end if

    return
  end subroutine sdm_shiftsx
  !--------------------------------------------------------------------------
  subroutine sdm_getbufsx(wbc,ebc,sd_num,sd_numasl,         &
       sd_n,sd_liqice,sd_x,sd_y,sd_rk,sd_u,sd_v,sd_vz, &
       sd_r,sd_asl,sdi,                            &
       bufsiz1,bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4, &
       stat,ilist,                               &
       rbuf_r8,rbuf_i8,rbuf_i2,rbuf_i4)

    use scale_grid, only: &
         GRID_FX
    use scale_grid_index, only: &
         IS,IE
    use m_sdm_common, only: &
         nisub,VALID2INVALID,VALID2INVALID_i8,i2,sdm_cold,sdicedef
    ! Input variables
    integer, intent(in) :: wbc       ! Option for west boundary conditions
    integer, intent(in) :: ebc       ! Option for east boundaty conditions
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of kind of chemical material contained as water-soluble aerosol in super droplets
    integer, intent(in) :: bufsiz1   ! buffer size for MPI
    integer, intent(in) :: bufsiz2_r8 ! buffer size for MPI (real8)
    integer, intent(in) :: bufsiz2_i8 ! buffer size for MPI (int8)
    integer, intent(in) :: bufsiz2_i2 ! buffer size for MPI (int2)
    integer, intent(in) :: bufsiz2_i4 ! buffer size for MPI (int4)
    real(DP), intent(in) :: rbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer in x direction (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf,mrime
                       ! dim03 = 1 : west, 2: east
    integer(DP), intent(in) ::                                  &
     &                              rbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! Sending buffer in x direction (int8)
                       ! dim02 = 1 : multiplicity 
                       ! dim03 = 1 : west, 2: east
    integer(i2), intent(in) ::                                  &
     &                              rbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! Sending buffer in x direction (int2)
                       ! dim02 = 1 : status (liq/ice)
                       ! dim03 = 1 : west, 2: east
    integer, intent(in) ::                                  &
     &                              rbuf_i4(1:,1:,1:)
                       ! rbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or rbuf_i4(1:1,1:1,1:2) 
                       ! Sending buffer in x direction (int4)
                       ! dim02 = 1 : nmono (ice)
                       ! dim03 = 1 : west, 2: east
    integer, intent(inout) :: stat ! Runtime status
    integer(DP), intent(inout) :: sd_n(1:sd_num)    ! multiplicity of super-droplets
    integer(i2), intent(inout) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
    real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    real(RP), intent(inout) :: sd_u(1:sd_num)  ! x-components velocity of super-droplets
    real(RP), intent(inout) :: sd_v(1:sd_num)  ! y-components velocity of super-droplets
    real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    type(sdicedef), intent(inout) :: sdi   ! ice phase super-droplet
    ! Output variable
    integer, intent(out) :: ilist(1:sd_num) ! buffer for list vectorization
    ! Work variables
    integer :: terecv           ! total list number
    integer :: twrecv           ! total list number
    integer :: tesend           ! total list number
    integer :: twsend           ! total list number

    integer :: ecnt             ! counter
    integer :: wcnt             ! counter

    integer :: m, n, nasl       ! index
    integer :: nwsdm            ! number of warm SDM variables
    !-------------------------------------------------------------------

    ! Initialize
    
    terecv = 0
    twrecv = 0
    tesend = 0
    twsend = 0

    nwsdm = 7 + sd_numasl  !! warm SDM variables

    ! Get the receiving buffer in x direction.

    if( nisub>=2 ) then

       ! Get the receiving buffer sent toward the west

       !== periodic boundary conditions ==!
       !! check the num of sent sd
       do n=1,bufsiz1
          if( rbuf_i8(n,1,1)>VALID2INVALID_i8 ) then
             twsend = twsend + 1
          end if
       end do

       !! check the num of acceptable sd
       wcnt   = 0
       do n=1,sd_num
          if( sd_rk(n)<VALID2INVALID ) then
             wcnt = wcnt + 1
             ilist(wcnt) = n
          end if
       end do
       twrecv = wcnt

       !! check send/receive
       if( stat==0 .and.                                           &
            &         (twrecv<twsend)) then
          stat = -1
       end if

       !! send/receive
       if( twrecv>0 .and. twsend>0 ) then
          do m=1,twsend
             n = ilist(m)

             sd_x(n)    = rbuf_r8(m,1,1)+GRID_FX(IE)
             sd_y(n)    = rbuf_r8(m,2,1)
             sd_rk(n)   = rbuf_r8(m,3,1)
             sd_u(n)    = rbuf_r8(m,4,1)
             sd_v(n)    = rbuf_r8(m,5,1)
             sd_vz(n)   = rbuf_r8(m,6,1)
             sd_r(n)    = rbuf_r8(m,7,1)

             sd_n(n)      = rbuf_i8(m,1,1)
             sd_liqice(n) = rbuf_i2(m,1,1)
          end do

          do nasl=1,sd_numasl
             do m=1,twsend
                n = ilist(m)

                sd_asl(n,nasl) = rbuf_r8(m,7+nasl,1)

             end do
          end do

          !!! cold-SDM only 
          if( sdm_cold ) then
          
             do m=1,twsend
                n = ilist(m)

                sdi%re(n)  = rbuf_r8(m,nwsdm+1,1)
                sdi%rp(n)  = rbuf_r8(m,nwsdm+2,1)
                sdi%rho(n) = rbuf_r8(m,nwsdm+3,1)
                sdi%tf(n)  = rbuf_r8(m,nwsdm+4,1)
                sdi%mrime(n)  = rbuf_r8(m,nwsdm+5,1)
                sdi%nmono(n)  = rbuf_i4(m,1,1)

             end do

          end if

       end if


       ! Get the receiving buffer sent toward the east
       
       !== periodic boundary conditions ==!
       !! check the num of sent sd
       do n=1,bufsiz1
          if( rbuf_i8(n,1,2)>VALID2INVALID_i8 ) then
             tesend = tesend + 1
          end if
       end do

       !! check the num of acceptable sd
       ecnt   = 0
       do n=1,sd_num
          if( sd_rk(n)<VALID2INVALID ) then
             ecnt = ecnt + 1
             ilist(ecnt) = n
          end if
       end do
       terecv = ecnt

       !! check send/receive
       if( stat==0 .and.                                           &
            &         (terecv<tesend)) then
          stat = -1
       end if

       !! send/receive
       if( terecv>0 .and. tesend>0 ) then
          do m=1,tesend
             n = ilist(m)

             sd_x(n)    = rbuf_r8(m,1,2)+GRID_FX(IS-1)
             sd_y(n)    = rbuf_r8(m,2,2)
             sd_rk(n)   = rbuf_r8(m,3,2)
             sd_u(n)    = rbuf_r8(m,4,2)
             sd_v(n)    = rbuf_r8(m,5,2)
             sd_vz(n)   = rbuf_r8(m,6,2)
             sd_r(n)    = rbuf_r8(m,7,2)

             sd_n(n)      = rbuf_i8(m,1,2)
             sd_liqice(n) = rbuf_i2(m,1,2)

          end do

          do nasl=1,sd_numasl
             do m=1,tesend
                n = ilist(m)

                sd_asl(n,nasl) = rbuf_r8(m,7+nasl,2)

             end do
          end do

          !!! cold-SDM only 
          if( sdm_cold ) then
             
             do m=1,tesend
                n = ilist(m)

                sdi%re(n)  = rbuf_r8(m,nwsdm+1,2)
                sdi%rp(n)  = rbuf_r8(m,nwsdm+2,2)
                sdi%rho(n) = rbuf_r8(m,nwsdm+3,2)
                sdi%tf(n)  = rbuf_r8(m,nwsdm+4,2)
                sdi%mrime(n)  = rbuf_r8(m,nwsdm+5,2)
                sdi%nmono(n)  = rbuf_i4(m,1,2)

             end do

          end if

       end if

    end if
      
    return
  end subroutine sdm_getbufsx
  !----------------------------------------------------------------------------
  subroutine sdm_putbufsy(sbc,nbc,sd_num,sd_numasl,         &
       sd_n,sd_liqice,sd_x,sd_y,sd_rk,sd_u,sd_v,sd_vz, &
       sd_r,sd_asl,sdi,                           &
       bufsiz1,bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4, &
       stat,ilist,                   &
       sbuf_r8,sbuf_i8,sbuf_i2,sbuf_i4)

    use scale_grid, only: &
         GRID_FY
    use scale_grid_index, only: &
         JS,JE
    use m_sdm_common, only: &
         njsub,INVALID,INVALID_i8,INVALID_i2,VALID2INVALID,i2,sdm_cold,sdicedef
    ! Input variables
    integer, intent(in) :: sbc   ! Option for west boundary conditions
    integer, intent(in) :: nbc   ! Option for east boundaty conditions
    integer, intent(in) :: sd_num  ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of kind of chemical materia contained as water-soluble aerosol in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(inout) :: sd_u(1:sd_num) ! x-components velocity of super-droplets
    real(RP), intent(inout) :: sd_v(1:sd_num) ! y-components velocity of super-droplets
    real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl)  ! aerosol mass of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplet
    integer, intent(in) :: bufsiz1  ! buffer size for MPI
    integer, intent(in) :: bufsiz2_r8 ! buffer size for MPI (real8)
    integer, intent(in) :: bufsiz2_i8 ! buffer size for MPI (int8)
    integer, intent(in) :: bufsiz2_i2 ! buffer size for MPI (int2)
    integer, intent(in) :: bufsiz2_i4 ! buffer size for MPI (int4)
    ! Input and output variable
    integer, intent(inout) :: stat ! Runtime status
    real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    ! Output variable
    real(DP), intent(out) :: sbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer in x direction (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf,mrime
                       ! dim03 = 1 : south, 2: north
    integer(DP), intent(out) ::                                  &
     &                              sbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! Sending buffer in x direction (int8)
                       ! dim02 = 1 : multiplicity 
                       ! dim03 = 1 : south, 2: north
    integer(i2), intent(out) ::                                  &
     &                              sbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! Sending buffer in x direction (int2)
                       ! dim02 = 1 : status (liq/ice)
                       ! dim03 = 1 : south, 2: north
    integer, intent(out) ::                                  &
     &                              sbuf_i4(1:,1:,1:)
                       ! sbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or sbuf_i4(1:1,1:1,1:2) 
                       ! Sending buffer in x direction (int4)
                       ! dim02 = 1 : nmono (ice)
                       ! dim03 = 1 : south, 2: north
    ! dim03 = 1 : south, 2: north
    integer, intent(out) :: ilist(1:sd_num)  ! buffer for list vectorization
    ! Work variables
    integer :: tssend           ! total list number
    integer :: tnsend           ! total list number
    integer :: scnt             ! counter
    integer :: ncnt             ! counter
    integer :: m, n    ! index
    integer :: nasl    ! index
    integer :: nwsdm            ! number of warm SDM variables
    !---------------------------------------------------------------------

    ! Initialize
    tssend = 0
    tnsend = 0
    
    nwsdm = 7 + sd_numasl  !! warm SDM variables

    sbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2) = INVALID
    sbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2) = INVALID_i8
    sbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2) = INVALID_i2

    ! Put the sending buffer in y direction.
    
    if( njsub>=2 ) then

       !### Fill in the sending buffer with the value ###!
       !### in the south halo regions.                ###!

       !== periodic boundary conditions ==!
       !! get number and index of outflowed super-droplet

       scnt   = 0
       do n=1,sd_num
          if( sd_y(n)<GRID_FY(JS-1) .and. sd_rk(n)>VALID2INVALID ) then
             scnt = scnt + 1
             ilist(scnt) = n
          end if
       end do
       tssend  = scnt

       !! check send buffer
       if( stat==0 .and.                                           &
            &         (tssend>bufsiz1) ) then
          stat = -1
       end if

       !! send data
       if( tssend>0 ) then
          do m=1,tssend
             n = ilist(m)

             sbuf_r8(m,1,1) = sd_x(n)
             sbuf_r8(m,2,1) = sd_y(n)-GRID_FY(JS-1)
             sbuf_r8(m,3,1) = sd_rk(n)
             sbuf_r8(m,4,1) = sd_u(n)
             sbuf_r8(m,5,1) = sd_v(n)
             sbuf_r8(m,6,1) = sd_vz(n)
             sbuf_r8(m,7,1) = sd_r(n)

             sbuf_i8(m,1,1) = sd_n(n)
             sbuf_i2(m,1,1) = sd_liqice(n)
          end do

          do nasl=1,sd_numasl
             do m=1,tssend
                n = ilist(m)

                sbuf_r8(m,7+nasl,1) = sd_asl(n,nasl)

             end do
          end do

          !!! cold-SDM only 
          if( sdm_cold ) then

             do m=1,tssend
                n = ilist(m)

                sbuf_r8(m,nwsdm+1,1) = sdi%re(n)
                sbuf_r8(m,nwsdm+2,1) = sdi%rp(n)
                sbuf_r8(m,nwsdm+3,1) = sdi%rho(n)
                sbuf_r8(m,nwsdm+4,1) = sdi%tf(n)
                sbuf_r8(m,nwsdm+5,1) = sdi%mrime(n)
                sbuf_i4(m,1,1) = sdi%nmono(n)

             end do
          end if

          !== convert to invalid ==!
          do m=1,tssend
             n = ilist(m)

             sd_rk(n) = INVALID
          end do
       end if


       !### Fill in the sending buffer with the value ###!
       !### in the north halo regions                 ###!
       
       !== periodic boundary conditions ==!
       !! get number and index of outflowed super-droplet
       
       ncnt   = 0
       do n=1,sd_num
          if( sd_y(n)>=GRID_FY(JE) .and. sd_rk(n)>VALID2INVALID ) then
             ncnt = ncnt + 1
             ilist(ncnt) = n
          end if
       end do
       tnsend  = ncnt

       !! check send buffer
       if( stat==0 .and.                                           &
            &         (tnsend>bufsiz1) ) then
          stat = -1
       end if

       !! send data
       if( tnsend>0 ) then
          do m=1,tnsend
             n = ilist(m)

             sbuf_r8(m,1,2) = sd_x(n)
             sbuf_r8(m,2,2) = sd_y(n)-GRID_FY(JE)
             sbuf_r8(m,3,2) = sd_rk(n)
             sbuf_r8(m,4,2) = sd_u(n)
             sbuf_r8(m,5,2) = sd_v(n)
             sbuf_r8(m,6,2) = sd_vz(n)
             sbuf_r8(m,7,2) = sd_r(n)

             sbuf_i8(m,1,2) = sd_n(n)
             sbuf_i2(m,1,2) = sd_liqice(n)
          end do

          do nasl=1,sd_numasl
             do m=1,tnsend
                n = ilist(m)

                sbuf_r8(m,7+nasl,2) = sd_asl(n,nasl)

             end do
          end do

          !!! cold-SDM only 
          if( sdm_cold ) then

             do m=1,tnsend
                n = ilist(m)

                sbuf_r8(m,nwsdm+1,2) = sdi%re(n)
                sbuf_r8(m,nwsdm+2,2) = sdi%rp(n)
                sbuf_r8(m,nwsdm+3,2) = sdi%rho(n)
                sbuf_r8(m,nwsdm+4,2) = sdi%tf(n)
                sbuf_r8(m,nwsdm+5,2) = sdi%mrime(n)
                sbuf_i4(m,1,2) = sdi%nmono(n)

             end do
          end if

          !== convert to invalid ==!
          do m=1,tnsend
             n = ilist(m)

             sd_rk(n) = INVALID
          end do
       end if

    end if

    return
  end subroutine sdm_putbufsy
  !----------------------------------------------------------------------------
  subroutine sdm_shiftsy(sbc,nbc,               &
       bufsiz1,bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4,&
       sbuf_r8,rbuf_r8,sbuf_i8,rbuf_i8,sbuf_i2,rbuf_i2,sbuf_i4,rbuf_i4)
    use mpi
    use m_sdm_common, only: &
         njsub, &
         dsts_sub, dstn_sub, srcs_sub, srcn_sub, &
         INVALID,INVALID_i8,INVALID_i2, &
         tag,i2,sdm_cold

    ! Input variables
    integer, intent(in) :: sbc      ! Option for south boundary conditions
    integer, intent(in) :: nbc      ! Option for north boundaty conditions
    integer, intent(in) :: bufsiz1  ! buffer size for MPI
     integer, intent(in) :: bufsiz2_r8 ! buffer size for MPI (real8)
    integer, intent(in) :: bufsiz2_i8 ! buffer size for MPI (int8)
    integer, intent(in) :: bufsiz2_i2 ! buffer size for MPI (int2)
    integer, intent(in) :: bufsiz2_i4 ! buffer size for MPI (int4)
    real(DP), intent(in)                                         &
     &                      :: sbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer in x direction (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf,mrime
                       ! dim03 = 1 : south, 2: north
    integer(DP), intent(in)                                      &
     &                      :: sbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! Sending buffer in x direction (int8)
                       ! dim02 = 1 : multiplicity
                       ! dim03 = 1 : south, 2: north
    integer(i2), intent(in)                                      &
     &                      :: sbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! Sending buffer in x direction (int2)
                       ! dim02 = 1 : status (liq/ice)
                       ! dim03 = 1 : south, 2: north
    integer, intent(in)                                      &
     &                      :: sbuf_i4(1:,1:,1:)
                       ! sbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or sbuf_i4(1:1,1:1,1:2) 
                       ! Sending buffer in x direction (int4)
                       ! dim02 = 1 : nmono (ice)
                       ! dim03 = 1 : south, 2: north
    ! Output variable
    real(DP), intent(out)                                        &
     &                      :: rbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer in x direction (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,t,tf,mrime
                       ! dim03 = 1 : south, 2: north
    integer(DP), intent(out)                                     &
     &                      :: rbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! Sending buffer in x direction (int8)
                       ! dim02 = 1 : multiplicity
                       ! dim03 = 1 : south, 2: north
    integer(i2), intent(out)                                     &
     &                      :: rbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! Sending buffer in x direction (int2)
                       ! dim02 = 1 : status (liq/ice)
                       ! dim03 = 1 : south, 2: north
    integer, intent(out)                                     &
     &                      :: rbuf_i4(1:,1:,1:)
                       ! rbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or rbuf_i4(1:1,1:1,1:2) 
                       ! Sending buffer in x direction (int4)
                       ! dim02 = 1 : nmono (ice)
                       ! dim03 = 1 : south, 2: north

    ! Internal shared variables
    integer :: dsts     ! South sending distnation
    integer :: dstn     ! North sending distnation
    integer :: srcs     ! South receiving source
    integer :: srcn     ! North receiving source
    integer :: statss   ! South sending request status
    integer :: statsn   ! North sending request status
    integer :: statrs   ! South receiving request status
    integer :: statrn   ! Morth receiving request status
    integer :: siz_r8   ! Sending and receiving buffer size
    integer :: siz_i8   ! Sending and receiving buffer size
    integer :: siz_i2   ! Sending and receiving buffer size
    integer :: siz_i4   ! Sending and receiving buffer size
    integer :: ierr     ! Error descriptor

    integer :: stat(mpi_status_size) ! Runtime status table
    !--------------------------------------------------------------------

    ! Initialize
    rbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2) = INVALID
    rbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2) = INVALID_i8
    rbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2) = INVALID_i2

    ! Exchange the value in y direction.

    if( njsub>=2 ) then

       !### Set the buffre size ###!

       siz_r8 = bufsiz1 * bufsiz2_r8
       siz_i8 = bufsiz1 * bufsiz2_i8
       siz_i2 = bufsiz1 * bufsiz2_i2
       siz_i4 = 0
       if( sdm_cold ) then
          siz_i4 = bufsiz1 * bufsiz2_i4
       end if

       !### Set the processor element number for sending ###!
       !### distination and receiving source             ###!

       dsts = dsts_sub
       dstn = dstn_sub
       srcs = srcs_sub
       srcn = srcn_sub

       !### Send/Receive real8 ###!

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function ###!

       call mpi_isend(sbuf_r8(1,1,1),siz_r8,MPI_DOUBLE_PRECISION,dsts,tag, &
            MPI_COMM_WORLD,statss,ierr)

       call mpi_irecv(rbuf_r8(1,1,1),siz_r8,MPI_DOUBLE_PRECISION,srcn,tag, &
            MPI_COMM_WORLD,statrn,ierr)

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function ###!

       call mpi_isend(sbuf_r8(1,1,2),siz_r8,MPI_DOUBLE_PRECISION,dstn,tag, &
            MPI_COMM_WORLD,statsn,ierr)

       call mpi_irecv(rbuf_r8(1,1,2),siz_r8,MPI_DOUBLE_PRECISION,srcs,tag, &
            MPI_COMM_WORLD,statrs,ierr)

       !### Call the waiting MPI function ###!

       call mpi_wait(statss,stat,ierr)
       call mpi_wait(statsn,stat,ierr)
       call mpi_wait(statrs,stat,ierr)
       call mpi_wait(statrn,stat,ierr)

       !### Send/Receive int8 ###!

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function ###!

       call mpi_isend(sbuf_i8(1,1,1),siz_i8,MPI_LONG_LONG_INT,dsts,tag, &
            MPI_COMM_WORLD,statss,ierr)

       call mpi_irecv(rbuf_i8(1,1,1),siz_i8,MPI_LONG_LONG_INT,srcn,tag, &
            MPI_COMM_WORLD,statrn,ierr)

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function ###!

       call mpi_isend(sbuf_i8(1,1,2),siz_i8,MPI_LONG_LONG_INT,dstn,tag, &
            MPI_COMM_WORLD,statsn,ierr)

       call mpi_irecv(rbuf_i8(1,1,2),siz_i8,MPI_LONG_LONG_INT,srcs,tag, &
            MPI_COMM_WORLD,statrs,ierr)

       !### Call the waiting MPI function ###!

       call mpi_wait(statss,stat,ierr)
       call mpi_wait(statsn,stat,ierr)
       call mpi_wait(statrs,stat,ierr)
       call mpi_wait(statrn,stat,ierr)

       !### Send/Receive int2 ###!

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function ###!

       call mpi_isend(sbuf_i2(1,1,1),siz_i2,MPI_INTEGER2,dsts,tag, &
            MPI_COMM_WORLD,statss,ierr)

       call mpi_irecv(rbuf_i2(1,1,1),siz_i2,MPI_INTEGER2,srcn,tag, &
            MPI_COMM_WORLD,statrn,ierr)

       !### Incliment the message tag ###!

       tag = tag + 1

       !### Call the sending and receiving MPI function ###!

       call mpi_isend(sbuf_i2(1,1,2),siz_i2,MPI_INTEGER2,dstn,tag, &
            MPI_COMM_WORLD,statsn,ierr)

       call mpi_irecv(rbuf_i2(1,1,2),siz_i2,MPI_INTEGER2,srcs,tag, &
            MPI_COMM_WORLD,statrs,ierr)

       !### Call the waiting MPI function ###!

       call mpi_wait(statss,stat,ierr)
       call mpi_wait(statsn,stat,ierr)
       call mpi_wait(statrs,stat,ierr)
       call mpi_wait(statrn,stat,ierr)

       !### Send/Receive int4 ###!

       if( sdm_cold ) then

          !### Incliment the message tag ###!

          tag = tag + 1

          !### Call the sending and receiving MPI function ###!

          call mpi_isend(sbuf_i4(1,1,1),siz_i4,MPI_INTEGER,dsts,tag, &
               MPI_COMM_WORLD,statss,ierr)

          call mpi_irecv(rbuf_i4(1,1,1),siz_i4,MPI_INTEGER,srcn,tag, &
               MPI_COMM_WORLD,statrn,ierr)

          !### Incliment the message tag ###!
          
          tag = tag + 1

          !### Call the sending and receiving MPI function ###!

          call mpi_isend(sbuf_i4(1,1,2),siz_i4,MPI_INTEGER,dstn,tag, &
               MPI_COMM_WORLD,statsn,ierr)

          call mpi_irecv(rbuf_i4(1,1,2),siz_i4,MPI_INTEGER,srcs,tag, &
               MPI_COMM_WORLD,statrs,ierr)

          !### Call the waiting MPI function ###!

          call mpi_wait(statss,stat,ierr)
          call mpi_wait(statsn,stat,ierr)
          call mpi_wait(statrs,stat,ierr)
          call mpi_wait(statrn,stat,ierr)

       end if

    end if

    return
  end subroutine sdm_shiftsy
  !----------------------------------------------------------------------------
  subroutine sdm_getbufsy(sbc,nbc,sd_num,sd_numasl,         &
       sd_n,sd_liqice,sd_x,sd_y,sd_rk,sd_u,sd_v,sd_vz, &
       sd_r,sd_asl,sdi,                            &
       bufsiz1,bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4, &
       stat,ilist,                               &
       rbuf_r8,rbuf_i8,rbuf_i2,rbuf_i4)
    use scale_grid, only: &
         GRID_FY
    use scale_grid_index, only: &
         JS,JE
    use m_sdm_common, only: &
         njsub,VALID2INVALID,VALID2INVALID_i8,i2,sdm_cold,sdicedef
    ! Input variables
    integer, intent(in) :: sbc       ! Option for south boundary conditions
    integer, intent(in) :: nbc       ! Option for north boundaty conditions

    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of kind of chemical material contained as water-soluble aerosol in super droplets
    integer, intent(in) :: bufsiz1   ! buffer size for MPI
    integer, intent(in) :: bufsiz2_r8 ! buffer size for MPI (real8)
    integer, intent(in) :: bufsiz2_i8 ! buffer size for MPI (int8)
    integer, intent(in) :: bufsiz2_i2 ! buffer size for MPI (int2)
    integer, intent(in) :: bufsiz2_i4 ! buffer size for MPI (int4)
    real(DP), intent(in) :: rbuf_r8(1:bufsiz1,1:bufsiz2_r8,1:2)
                       ! Sending buffer in x direction (real8)
                       ! dim02 = 1 - ( 7+sd_numasl (+5) )
                       !   : [liq] x,y,rk,u,v,wc(vz),r,asl(1:sd_numasl)
                       !   : [ice] re,rp,rho,tf,mrime
                       ! dim03 = 1 : west, 2: east
    integer(DP), intent(in) ::                                  &
     &                              rbuf_i8(1:bufsiz1,1:bufsiz2_i8,1:2)
                       ! Sending buffer in x direction (int8)
                       ! dim02 = 1 : multiplicity 
                       ! dim03 = 1 : west, 2: east
    integer(i2), intent(in) ::                                  &
     &                              rbuf_i2(1:bufsiz1,1:bufsiz2_i2,1:2)
                       ! Sending buffer in x direction (int2)
                       ! dim02 = 1 : status (liq/ice)
                       ! dim03 = 1 : west, 2: east
    integer, intent(in) ::                                  &
     &                              rbuf_i4(1:,1:,1:)
                       ! rbuf_i4(1:bufsiz1,1:bufsiz2_i4,1:2) or rbuf_i4(1:1,1:1,1:2) 
                       ! Sending buffer in x direction (int4)
                       ! dim02 = 1 : nmono (ice)
                       ! dim03 = 1 : west, 2: east
    integer, intent(inout) :: stat  ! Runtime status
    integer(DP), intent(inout) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    integer(i2), intent(inout) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
    real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    real(RP), intent(inout) :: sd_u(1:sd_num)  ! x-components velocity of super-droplets
    real(RP), intent(inout) :: sd_v(1:sd_num)  ! y-components velocity of super-droplets
    real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl)  ! aerosol mass of super-droplets
    type(sdicedef), intent(inout) :: sdi   ! ice phase super-droplet
    ! Output variable
    integer, intent(out) :: ilist(1:sd_num)  ! buffer for list vectorization
    ! Work variables
    integer :: tsrecv           ! total list number
    integer :: tnrecv           ! total list number
    integer :: tssend           ! total list number
    integer :: tnsend           ! total list number

    integer :: scnt             ! index
    integer :: ncnt             ! index

    integer :: m, n, nasl       ! index
    integer :: nwsdm           ! number of warm SDM variables
    !---------------------------------------------------------------------

    ! Initialize
    
    tsrecv = 0
    tnrecv = 0
    tssend = 0
    tnsend = 0

    nwsdm = 7 + sd_numasl  !! warm SDM variables

    ! Get the receiving buffer in y direction.

    if( njsub>=2 ) then

       ! Get the receiving buffer sent toward the south

       !== periodic boundary conditions ==!
       !! check the num of sent sd
       do n=1,bufsiz1
          if( rbuf_i8(n,1,1)>VALID2INVALID_i8 ) then
             tssend = tssend + 1
          end if
       end do

       !! check the num of acceptable sd
       scnt   = 0
       do n=1,sd_num
          if( sd_rk(n)<VALID2INVALID ) then
             scnt = scnt + 1
             ilist(scnt) = n
          end if
       end do
       tsrecv = scnt

       !! check send/receive
       if( stat==0 .and.                                           &
            (tsrecv<tssend) ) then
          stat = -1
       end if

       !! send/receive
       if( tsrecv>0 .and. tssend>0 ) then
          do m=1,tssend
             n = ilist(m)

             sd_x(n)    = rbuf_r8(m,1,1)
             sd_y(n)    = rbuf_r8(m,2,1)+GRID_FY(JE)
             sd_rk(n)   = rbuf_r8(m,3,1)
             sd_u(n)    = rbuf_r8(m,4,1)
             sd_v(n)    = rbuf_r8(m,5,1)
             sd_vz(n)   = rbuf_r8(m,6,1)
             sd_r(n)    = rbuf_r8(m,7,1)

             sd_n(n)      = rbuf_i8(m,1,1)
             sd_liqice(n) = rbuf_i2(m,1,1)
          end do

          do nasl=1,sd_numasl
             do m=1,tssend
                n = ilist(m)

                sd_asl(n,nasl) = rbuf_r8(m,7+nasl,1)

             end do
          end do

          !!! cold-SDM only 
          if( sdm_cold ) then
          
             do m=1,tssend
                n = ilist(m)

                sdi%re(n)  = rbuf_r8(m,nwsdm+1,1)
                sdi%rp(n)  = rbuf_r8(m,nwsdm+2,1)
                sdi%rho(n) = rbuf_r8(m,nwsdm+3,1)
                sdi%tf(n)  = rbuf_r8(m,nwsdm+4,1)
                sdi%mrime(n)  = rbuf_r8(m,nwsdm+5,1)
                sdi%nmono(n)  = rbuf_i4(m,1,1)

             end do

          end if

       end if

       ! Get the receiving buffer sent toward the north
       
       !== periodic boundary conditions ==!
       !! check the num of sent sd
       do n=1,bufsiz1
          if( rbuf_i8(n,1,2)>VALID2INVALID_i8 ) then
             tnsend = tnsend + 1
          end if
       end do

       !! check the num of acceptable sd
       ncnt   = 0
       do n=1,sd_num
          if( sd_rk(n)<VALID2INVALID ) then
             ncnt = ncnt + 1
             ilist(ncnt) = n
          end if
       end do
       tnrecv = ncnt

       !! check send/receive
       if( stat==0 .and.                                           &
            &         (tnrecv<tnsend) ) then
          stat = -1
       end if

       !! send/receive
       if( tnrecv>0 .and. tnsend>0 ) then
          do m=1,tnsend
             n = ilist(m)

             sd_x(n)    = rbuf_r8(m,1,2)
             sd_y(n)    = rbuf_r8(m,2,2)+GRID_FY(JS-1)
             sd_rk(n)   = rbuf_r8(m,3,2)
             sd_u(n)    = rbuf_r8(m,4,2)
             sd_v(n)    = rbuf_r8(m,5,2)
             sd_vz(n)   = rbuf_r8(m,6,2)
             sd_r(n)    = rbuf_r8(m,7,2)

             sd_n(n)      = rbuf_i8(m,1,2)
             sd_liqice(n) = rbuf_i2(m,1,2)
          end do

          do nasl=1,sd_numasl
             do m=1,tnsend
                n = ilist(m)

                sd_asl(n,nasl) = rbuf_r8(m,7+nasl,2)

             end do
          end do

          !!! cold-SDM only 
          if( sdm_cold ) then
             
             do m=1,tnsend
                n = ilist(m)

                sdi%re(n)  = rbuf_r8(m,nwsdm+1,2)
                sdi%rp(n)  = rbuf_r8(m,nwsdm+2,2)
                sdi%rho(n) = rbuf_r8(m,nwsdm+3,2)
                sdi%tf(n)  = rbuf_r8(m,nwsdm+4,2)
                sdi%mrime(n)  = rbuf_r8(m,nwsdm+5,2)
                sdi%nmono(n)  = rbuf_i4(m,1,2)

             end do

          end if

       end if
    end if

    return
  end subroutine sdm_getbufsy
end module m_sdm_boundary
