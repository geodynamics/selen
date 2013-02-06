!
! this is program "PX-REBUILD.F90"   
!
! Created DM 09-10-2012
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Copyright (C) 2008 Giorgio Spada, Florence Colleoni, and Paolo Stocchi 
!
! This file is part of SELEN. 
!  
! SELEN is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or at your option) any later 
! version. 
!
! SELEN is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!		
!
! This script creates the main latitudes file "px-lat.dat" from
! file "px-table.dat".
!
!
!
 PROGRAM PXREBUILD
 IMPLICIT NONE 
 INCLUDE "data.inc"
 INTEGER I, J, K, NANCH
 INTEGER, ALLOCATABLE :: IANCH(:), IWET(:), IMIN(:), IMAX(:)
 REAL*8, ALLOCATABLE :: LON(:), LAT(:), LATA(:)
 CHARACTER CJUNK 
!
!
! --- Allocate dynamic vectors
!
 ALLOCATE( IANCH(NP), IWET(NP) )
 ALLOCATE( LON(NP), LAT(NP) )
 ALLOCATE( IMIN(NP), IMAX(NP) )
 ALLOCATE( LATA(NP) )
!
!
! --- Read pixel table
!
  open(2,file='px-table.dat',status='old') 
  do i=1, 4 
    read(2,*) 
  enddo
!
  do i=1,np
    read(2,*) lon(i), lat(i), ianch(i), k, iwet(i)
  end do
  close(2)
!
!
! --- Identify anchor pixels
!
 nanch=1
 imin(nanch)=1 
 imax(nanch)=1 
 lata(nanch)=lat(1)
!
 do i=2, np                
    if( ianch(i) .ne. ianch(i-1) ) then
       nanch=nanch+1
       imax(nanch-1)=i-1
       imin(nanch)  =i
       lata(nanch)  =lat(i)
    end if 
 end do
!
 imax(nanch)=np
!
!
! --- Write output files
!
 open(3,file='px-lat.dat',status='unknown')
!
 Write(3,*) "-------------------------------------------"  
 Write(3,*) " Main pixels latitudes and range of pixels "
 Write(3,*) "  pixel, lat (deg), min range, max range   "   
 Write(3,*) "-------------------------------------------" 
!
 do j=1, nanch
     write(3,*) j, lata(j), imin(j), imax(j)
 enddo
!
 close(3)
!
! 
!
 open(4,file='px.dat',status='unknown')
!
 do j=1, np
     write(4,*) lon(j), lat(j), j
 enddo
!
 close(4)
!
!
!
 open(5,file='drya.dat',status='unknown')
 open(6,file='weta.dat',status='unknown')
!
 do j=5,6
    Write(j,*) "----------------------------------------------"
    Write(j,*) " Pixels lon-lats and with increasing latitude "
    Write(j,*) "  pixel lon(deg), lat (deg), anchor, pixel #  "
    Write(j,*) "----------------------------------------------"
 end do
!
 do j=1, np
     write(5+iwet(j),*) lon(j), lat(j), ianch(j), j
 enddo
!
 close(5)
 close(6)
!
!
!
 open(8,file='pxa.dat',status='unknown')
!
 Write(8,*) "----------------------------------------------"
 Write(8,*) " Pixels lon-lats and with increasing latitude "
 Write(8,*) "  pixel lon(deg), lat (deg), anchor, pixel #  "
 Write(8,*) "----------------------------------------------"
!
 do j=1, np
     write(8,*) lon(j), lat(j), ianch(j), j
 enddo
!
 close(8)
!
! ---
!
 open(9,file='anchor.tmp',status='unknown')
 write(9,*) nanch
 close(9)
!
!
! --- Deallocate dynamic vectors
!
 DEALLOCATE( IANCH, IWET )
 DEALLOCATE( LON, LAT )
 DEALLOCATE( IMIN, IMAX )
 DEALLOCATE( LATA )
!
!
!
 end program pxrebuild
!
