!
! this is program "PX-REC.F90"   
!
! Last modified GS 02-10-2008
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Reviewed DM June 2011 - Dynamic memory allocation
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
! This script composes a pixels table "px-table" using 
! wet and dry pixels distribution obtained via "px.gmt". 
!
!
!
 PROGRAM PXREC
 IMPLICIT NONE 
 INCLUDE "data.inc"
 INTEGER I, J, K
 INTEGER NHL
 INTEGER, ALLOCATABLE :: IANCH(:), IWET(:)
 REAL*8, ALLOCATABLE :: LON(:), LAT(:)
 REAL*8 :: XLON, XLAT
 CHARACTER CJUNK 
!
!
! --- Allocate dynamic vectors
!
 ALLOCATE( IANCH(NP), IWET(NP) )
 ALLOCATE( LON(NP), LAT(NP) )
!
!
! --- Reading "wet" and "dry" pixels information 
!
	open(1,file='weta.dat',status='unknown')
	open(2,file='drya.dat',status='unknown') 
        call count_header_lines(1,nhl)
		do j=1, nhl
			read(1,'(a30)') cjunk 
			read(2,'(a30)') cjunk 	
		enddo
	do j=1, np 
		read(1,*,end=10) xlon, xlat, k, i
		lon(i) = xlon
                lat(i) = xlat
                ianch(i) = k
                iwet(i) = 1
	enddo
10 	continue
	do j=1, np 
		read(2,*,end=11) xlon, xlat, k, i
                lon(i) = xlon
                lat(i) = xlat
                ianch(i) = k
                iwet(i) = 0
	enddo 
11 	continue
	close(1)
	close(2) 
!
!
! --- Composing "px-table.dat" from "wet" and "dry" files
!	
!Write(*,*) " "
!Write(*,*) " ---> Composing a pixels table from wet & dry ta"
	Open(3,file='px-table.dat',status='unknown') 
	Write(3,*) "------------------------------------------------"
	Write(3,*) " * Table with all the info about pixelization * "
	Write(3,*) "  lon(deg), lat(deg), anchor, pixel#, wet(1/0)  "
	Write(3,*) "------------------------------------------------"
!
	do 500 i=1, np
!
	     write(3,*) lon(i), lat(i), ianch(i), i, iwet(i)
!
500	continue 
!
	close(3)
!
!
! --- Release memory
!
 DEALLOCATE( IANCH, IWET ) 
 DEALLOCATE( LON, LAT )
!
 end program pxrec
!
!
!
