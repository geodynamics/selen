! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'config'
! reads parameters from the configuration file
! 
! DM Feb 25, 2021
! DM Sep 13, 2025
!
!
 module config_data
 implicit none
!
 integer :: cfg_l_max
 integer :: cfg_mp
!
 integer :: cfg_nn
 real(4) :: cfg_delta_t
!
 integer :: cfg_nv
!
 character(200) :: cfg_model_file
!
! 0: legacy format, 1: new format
! (when changing this parameter, delete also 'taboo.o', 'fluid-ln.o' and recompile)
 integer, parameter :: cfg_model_fmt = 0
!
 end module config_data
!
!
!
 subroutine config
 use config_data
 implicit none
!
 character(200) :: buffer   
!
 if (iargc().ne.1) then
    call getarg(0, buffer)
    write(*,*)
    !write(*,*) 'Usage: taboo.exe cfg_file'
    write(*,*) 'Usage: '//trim(adjustl(buffer))//' cfg_file'
    write(*,*)
    stop
 end if
!
 call getarg(1,buffer)
!
 open(90,file=trim(adjustl(buffer)),status='old')
!
!---------------------------------------------- Parsing configuration file
!
!call read_data_line(90,buffer)   ;   read(buffer,*) cfg_l_min
 call read_data_line(90,buffer)   ;   read(buffer,*) cfg_l_max
!
 call read_data_line(90,buffer)   ;   read(buffer,*) cfg_mp
!
 call read_data_line(90,buffer)   ;   read(buffer,*) cfg_nn
 call read_data_line(90,buffer)   ;   read(buffer,*) cfg_delta_t
!
 call read_data_line(90,buffer)   ;   read(buffer,*) cfg_nv
!
 call read_data_line(90,buffer)   ;   cfg_model_file=trim(adjustl(buffer))
!
 close(90)
!
!
 end subroutine config
!
