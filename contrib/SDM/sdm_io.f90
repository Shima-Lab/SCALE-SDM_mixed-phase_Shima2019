!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Input and Output of the SDM variables
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
!! @li      2014-06-27 (S.Shima) [new] sdm_outasci is added
!! @li      2016-07-11 (S.Shima) [mod] modified to support sdice output
!! @li      2018-06-25 (S.Shima) [add] netcdf output
!! @li      2018-06-30 (S.Shima) [add] rime mass and number of monomers as SD attributes
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_io

  implicit none
  private
  public :: sdm_outasci,sdm_outnetcdf

contains
  subroutine sdm_outasci(otime,sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sdi,sdn_dmpnskip)
    use scale_precision
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank, &
         PRC_MPIstop
    use m_sdm_common, only: &
         i2, sdm_cold, STAT_LIQ, STAT_ICE, STAT_MIX, sdicedef

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of chemical species contained in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(kind=i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_z(1:sd_num) ! z-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    integer, intent(in) :: sdn_dmpnskip ! Base skip to store super droplets in text format

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=17) :: fmt3="(3A)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=H_LONG) :: basename_time
    integer :: fid_sdm_o
    integer :: n, m, ierr
    character(len=80) :: fmt     ! output formate
    character(len=5)  :: cstat   ! status character
    
    !--- output Super Droplets in ASCII format
    write(basename_time,'(F15.3)') otime
    do n = 1, 15
       if( basename_time(n:n) == ' ' ) basename_time(n:n) = '0'
    enddo
    fid_sdm_o = IO_get_available_fid()

    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    write(ftmp,fmt3) 'SD_output', '_ASCII_', trim(basename_time)
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    open (fid_sdm_o, file = trim(basename_sd_out), & !action = "write", &
          access = "sequential", status = "replace", form = "formatted", &
          iostat = ierr)

    if( ierr /= 0 ) then
      write(*,*) "sdm_ascii_out", "Write error"
      call PRC_MPIstop
    endif 

    if( .not.sdm_cold ) then

       write(fid_sdm_o,'(3a,i2.2,2a)') '# x[m],y[m],z[m],vz[m],',       &
            &                'radius(droplet)[m],',                        &
            &                'mass_of_aerosol_in_droplet(1:',sd_numasl,')[g],',&
            &                'multiplicity[-],status[-],index'

       write(fmt,'( "(", i2.2, "e16.8,i20,a5,i10)" )')(sd_numasl+5)

       do m=1,sd_num,sdn_dmpnskip

          cstat = '  LIQ'
          
          write(fid_sdm_o,trim(fmt)) sd_x(m),           &
               &                   sd_y(m),           &
               &                   sd_z(m), sd_vz(m), sd_r(m),        &
               &                   (sd_asl(m,n),n=1,sd_numasl),       &
               &                   sd_n(m), cstat, m
       end do

    else
!       write(fid_sdm_o,'(3a,i2.2,6a)') '# x[m],y[m],z[m],vz[m],',       &
       write(fid_sdm_o,'(3a,i2.2,5a)') '# x[m],y[m],z[m],vz[m],',       &
            &            'radius(droplet)[m],',                            &
            &            'mass_of_aerosol_in_droplet/ice(1:',sd_numasl,')[g],',&
            &            'radius_eq(ice)[m],radius_pol(ice)[m],',              &
            &            'density(droplet/ice)[kg/m3],',                       &
!            &            'temperature[K],',                                    &
            &            'freezing_temp.(ice)[deg],',                          &
            &            'multiplicity[-],status[-],index,rime_mass[kg],num_of_monomers[-]'

!       write(fmt,'( "(", i2.2, "e16.8,i20,a5,i10)" )')(sd_numasl+10)
       write(fmt,'( "(", i2.2, "e16.8,i20,a5,i10,e16.8,i10)" )')(sd_numasl+9)

       do m=1,sd_num,sdn_dmpnskip

          if( sd_liqice(m)==STAT_LIQ ) then
             cstat = '  LIQ'
          else if( sd_liqice(m)==STAT_ICE ) then
             cstat = '  ICE'
          else if( sd_liqice(m)==STAT_MIX ) then
             cstat = ' MELT'
          end if

          write(fid_sdm_o,trim(fmt)) sd_x(m),           &
               &                   sd_y(m),           &
               &                   sd_z(m), sd_vz(m), sd_r(m),        &
               &                   (sd_asl(m,n),n=1,sd_numasl),       &
               &                   sdi%re(m), sdi%rp(m), sdi%rho(m),  &
!               &                   sdi%t(m), sdi%tf(m),               &
               &                   sdi%tf(m),               &
               &                   sd_n(m), cstat, m, sdi%mrime(m), sdi%nmono(m)
       end do


    end if

    close(fid_sdm_o)
    if( IO_L ) write(IO_FID_LOG,*) '*** Closed output file (ASCII) of Super Droplet'

    return

  end subroutine sdm_outasci
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sdm_outnetcdf(otime,sd_num,sd_numasl,sd_n,sd_liqice,sd_x,sd_y,sd_z,sd_r,sd_asl,sd_vz,sdi,sdn_dmpnskip)
    use netcdf
    use scale_precision
    use scale_stdio
    use scale_process, only: &
         mype => PRC_myrank
    use m_sdm_common, only: &
         i2, sdm_cold, STAT_LIQ, STAT_ICE, STAT_MIX, sdicedef

    implicit none

    real(DP), intent(in) :: otime
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer, intent(in) :: sd_numasl ! number of chemical species contained in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(kind=i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_z(1:sd_num) ! z-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    integer, intent(in) :: sdn_dmpnskip ! Base skip to store super droplets in text format

    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=17) :: fmt3="(3A)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_sd_out
    character(len=H_LONG) :: basename_time
    integer :: fid_sdm_o
    integer :: n, m, ierr
    character(len=80) :: fmt     ! output formate
    character(len=5)  :: cstat   ! status character
    integer :: nf90_real_precision
    integer :: ncid, sd_num_id, sd_numasl_id
    integer :: sd_x_id, sd_y_id, sd_z_id, sd_vz_id, sd_r_id, sd_asl_id, sd_n_id, sd_liqice_id
    integer :: sdi_re_id, sdi_rp_id, sdi_rho_id, sdi_tf_id, sdi_mrime_id, sdi_nmono_id

    !--- output Super Droplets in NetCDF format
    write(basename_time,'(F15.3)') otime
    do n = 1, 15
       if( basename_time(n:n) == ' ' ) basename_time(n:n) = '0'
    enddo
    fid_sdm_o = IO_get_available_fid()

    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    write(ftmp,fmt3) 'SD_output', '_NetCDF_', trim(basename_time)
    write(basename_sd_out,fmt2) trim(ftmp), 'pe',mype

    ! Check the presision
    if(RP == SP)then
       nf90_real_precision = NF90_FLOAT
    else
       nf90_real_precision = NF90_DOUBLE
    end if

    ! Open the output file
    call check_netcdf( nf90_create(trim(basename_sd_out),NF90_NETCDF4, ncid) )

    ! Definition of dimensions and variables
    !!! sd_num
    call check_netcdf( nf90_def_dim(ncid, "sd_num", sd_num, sd_num_id) )
    !!! sd_numasl
    call check_netcdf( nf90_def_dim(ncid, "sd_numasl", sd_numasl, sd_numasl_id) )

    !!! sd_x
    call check_netcdf( nf90_def_var(ncid, "sd_x", nf90_real_precision, sd_num_id, sd_x_id) )
    call check_netcdf( nf90_put_att(ncid, sd_x_id, 'long_name', 'x-coordinate') )
    call check_netcdf( nf90_put_att(ncid, sd_x_id, 'units', 'm') )
    !!! sd_y
    call check_netcdf( nf90_def_var(ncid, "sd_y", nf90_real_precision, sd_num_id, sd_y_id) )
    call check_netcdf( nf90_put_att(ncid, sd_y_id, 'long_name', 'y-coordinate') )
    call check_netcdf( nf90_put_att(ncid, sd_y_id, 'units', 'm') )
    !!! sd_z
    call check_netcdf( nf90_def_var(ncid, "sd_z", nf90_real_precision, sd_num_id, sd_z_id) )
    call check_netcdf( nf90_put_att(ncid, sd_z_id, 'long_name', 'z-coordinate') )
    call check_netcdf( nf90_put_att(ncid, sd_z_id, 'units', 'm') )
    !!! sd_vz
    call check_netcdf( nf90_def_var(ncid, "sd_vz", nf90_real_precision, sd_num_id, sd_vz_id) )
    call check_netcdf( nf90_put_att(ncid, sd_vz_id, 'long_name', 'terminal velocity') )
    call check_netcdf( nf90_put_att(ncid, sd_vz_id, 'units', 'm/s') )
    !!! sd_r
    call check_netcdf( nf90_def_var(ncid, "sd_r", nf90_real_precision, sd_num_id, sd_r_id) )
    call check_netcdf( nf90_put_att(ncid, sd_r_id, 'long_name', 'equivalent radius of liquid droplets') )
    call check_netcdf( nf90_put_att(ncid, sd_r_id, 'units', 'm') )
    !!! sd_asl
    call check_netcdf( nf90_def_var(ncid, "sd_asl", nf90_real_precision, (/sd_num_id, sd_numasl_id/), sd_asl_id) )
    call check_netcdf( nf90_put_att(ncid, sd_asl_id, 'long_name', 'aerosol mass') )
    call check_netcdf( nf90_put_att(ncid, sd_asl_id, 'units', 'g') )

    !!! sd_n
    call check_netcdf( nf90_def_var(ncid, "sd_n", NF90_INT64, sd_num_id, sd_n_id) )
    call check_netcdf( nf90_put_att(ncid, sd_n_id, 'long_name', 'multiplicity') )
    call check_netcdf( nf90_put_att(ncid, sd_n_id, 'units', '') )
    !!! sd_liqice
    call check_netcdf( nf90_def_var(ncid, "sd_liqice", NF90_SHORT, sd_num_id, sd_liqice_id) )
    call check_netcdf( nf90_put_att(ncid, sd_liqice_id, 'long_name', 'status of droplets: 01=liquid, 10=ice, 11=mixture') )
    call check_netcdf( nf90_put_att(ncid, sd_liqice_id, 'units', '') )

    if( sdm_cold ) then
       !!! sdi%re
       call check_netcdf( nf90_def_var(ncid, "sdi_re", nf90_real_precision, sd_num_id, sdi_re_id) )
       call check_netcdf( nf90_put_att(ncid, sdi_re_id, 'long_name', 'equatorial radius of ice crystals') )
       call check_netcdf( nf90_put_att(ncid, sdi_re_id, 'units', 'm') )
       !!! sdi%rp
       call check_netcdf( nf90_def_var(ncid, "sdi_rp", nf90_real_precision, sd_num_id, sdi_rp_id) )
       call check_netcdf( nf90_put_att(ncid, sdi_rp_id, 'long_name', 'polar radius of ice crystals') )
       call check_netcdf( nf90_put_att(ncid, sdi_rp_id, 'units', 'm') )
       !!! sdi%rho
       call check_netcdf( nf90_def_var(ncid, "sdi_rho", nf90_real_precision, sd_num_id, sdi_rho_id) )
       call check_netcdf( nf90_put_att(ncid, sdi_rho_id, 'long_name', 'density of ice crystals') )
       call check_netcdf( nf90_put_att(ncid, sdi_rho_id, 'units', 'kg/m3') )
       !!! sdi%tf
       call check_netcdf( nf90_def_var(ncid, "sdi_tf", nf90_real_precision, sd_num_id, sdi_tf_id) )
       call check_netcdf( nf90_put_att(ncid, sdi_tf_id, 'long_name', 'freezing temperature of particles') )
       call check_netcdf( nf90_put_att(ncid, sdi_tf_id, 'units', 'degC') )
       !!! sdi%mrime
       call check_netcdf( nf90_def_var(ncid, "sdi_mrime", nf90_real_precision, sd_num_id, sdi_mrime_id) )
       call check_netcdf( nf90_put_att(ncid, sdi_mrime_id, 'long_name', 'rime mass') )
       call check_netcdf( nf90_put_att(ncid, sdi_mrime_id, 'units', 'kg') )
       !!! sdi%nmono
       call check_netcdf( nf90_def_var(ncid, "sdi_nmono", NF90_INT, sd_num_id, sdi_nmono_id) )
       call check_netcdf( nf90_put_att(ncid, sdi_nmono_id, 'long_name', 'number of monomers (primary ice crystals)') )
       call check_netcdf( nf90_put_att(ncid, sdi_nmono_id, 'units', '') )
    end if

    !!! End of definition
    call check_netcdf( nf90_enddef(ncid) )

    ! Save data
    !!! sd_x
    call check_netcdf( nf90_put_var(ncid, sd_x_id, sd_x) )
    !!! sd_y
    call check_netcdf( nf90_put_var(ncid, sd_y_id, sd_y) )
    !!! sd_z
    call check_netcdf( nf90_put_var(ncid, sd_z_id, sd_z) )
    !!! sd_vz
    call check_netcdf( nf90_put_var(ncid, sd_vz_id, sd_vz) )
    !!! sd_r
    call check_netcdf( nf90_put_var(ncid, sd_r_id, sd_r) )
    !!! sd_asl
    call check_netcdf( nf90_put_var(ncid, sd_asl_id, sd_asl) )

    !!! sd_n
    call check_netcdf( nf90_put_var(ncid, sd_n_id, sd_n) )
    !!! sd_liqice
    call check_netcdf( nf90_put_var(ncid, sd_liqice_id, sd_liqice) )

    if( sdm_cold ) then
       !!! sdi_re
       call check_netcdf( nf90_put_var(ncid, sdi_re_id, sdi%re) )
       !!! sdi_rp
       call check_netcdf( nf90_put_var(ncid, sdi_rp_id, sdi%rp) )
       !!! sdi_rho
       call check_netcdf( nf90_put_var(ncid, sdi_rho_id, sdi%rho) )
       !!! sdi_tf
       call check_netcdf( nf90_put_var(ncid, sdi_tf_id, sdi%tf) )
       !!! sdi_mrime
       call check_netcdf( nf90_put_var(ncid, sdi_mrime_id, sdi%mrime) )
       !!! sdi_nmono
       call check_netcdf( nf90_put_var(ncid, sdi_nmono_id, sdi%nmono) )
    end if

    ! Close the output file
    call check_netcdf( nf90_close(ncid) )

    if( IO_L ) write(IO_FID_LOG,*) '*** Closed output file (NetCDF) of Super Droplet'

    return

  end subroutine sdm_outnetcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_netcdf(status)
    use netcdf
    use scale_process, only: &
         PRC_MPIstop

    integer, intent (in) :: status
    
    if(status /= nf90_noerr) then 
       write(*,*) "sdm_netcdf_out: Write error ",status
       call PRC_MPIstop
    end if
  end subroutine check_netcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_sdm_io
