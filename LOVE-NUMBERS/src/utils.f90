! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! subroutine 'read_data_line'
! Reads from file opened at LUN the next line not starting with '!' or '#'
!
! Initial version DM February 24, 2020
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
    subroutine read_data_line(lun,buffer)
    implicit none
    character(*)   :: buffer
    character :: ch
    integer :: lun
    integer :: idx
!
910 read(lun,'(a)') buffer
    ch = buffer(1:1)
    if ((ch.eq.'#').or.(ch.eq.'!')) go to 910
!
    if( (index( buffer, '#' ).ne.0) .or. (index( buffer,'!' ).ne.0) ) then
! 
       if( index(buffer,'#').eq.0 ) then
           idx = index(buffer,'!')
       elseif( index(buffer,'!').eq.0 ) then
           idx = index(buffer,'#')
       else
           idx = min( index( buffer, '#' ), index( buffer, '!' ) )
       endif	   
!
       buffer = buffer( 1:(idx-1) )
!
    end if
!
    return
!
    end
!
!
!
    subroutine count_rows(lun,n)
    implicit none
!
! --- Returns as N the number of rows in text file opened at LUN
!     DM 28.04.2014
!
    integer lun,n
!
    rewind(lun)
!    
    n=0
901 read(lun,*,end=902)
    n=n+1
    go to 901
    902 continue
!
    return
!
    end subroutine count_rows
!    
!
!
!
!
    subroutine count_headers(lun,n)
    implicit none
!
! --- Returns as N the number of header lines in text file opened at LUN
!     Header lines begin with '#' or '!'
!     DM 01.05.2022
!
    integer lun,n
    character :: c
!
    rewind(lun)
!    
    n=0
901 read(lun,*,end=902) c
    if ( ( c.ne.'#' ) .and. ( c.ne.'!' ) ) go to 902
    n=n+1
    go to 901
    902 continue
!
    return
!
    end subroutine count_headers
!
!
!    
