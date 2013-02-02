!
! This is program "SH_OF.F90" 
!
! Last modified GS 04-11-2008  Intel port 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Reviewed DM June 2011 - Dynamic memory allocation
! *** Modified DM June 2011 - Parallel execution
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
! --------------------------------------------------------------------------- 
! Computes the harmonic coefficients of the ocean function up to degree lmax.
! Input files are the pixels file "px-table.dat" and the set of CSHs stored 
! in the binary file "sh.bin". The harmonic coefficients are reported on file 
! 'shof.dat' in the format j(l,m), real part, imaginary part. 
! ---------------------------------------------------------------------------
!
! Input files:
!	- px-table.dat 
!	- sh.bin
!	
! Output files:
!	- shof.dat
!
!
! INCLUDE "harmonics.f90"
 PROGRAM OF 
 IMPLICIT NONE 
 INCLUDE "data.inc"
 CHARACTER*12 JUNK
 INTEGER IJUNK
 INTEGER I, J, K, MJ, NWET, NDRY, NPIX, NA
 INTEGER, ALLOCATABLE :: MM(:), WET(:), ANC(:) 
 REAL*8, ALLOCATABLE :: LONP(:), LATP(:)
 REAL*8 LLON, LLAT
 REAL*8 OCR, OCI
 REAL*8, ALLOCATABLE :: ALF(:,:)
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:), OC(:)
!
!
! --- Determining the number of pixels
!
  open(1,file='anchor.tmp',status='old')
  read(1,*) na
  close(1)
  Write(*,*) "    - Found ", na, " anchor pixels in file px-lat.dat"
  npix=np
!
!
!
! --- Allocate memory space
!
 ALLOCATE( MM(JMAX) )
 ALLOCATE( OC(JMAX) )
 ALLOCATE( WET(NPIX), ANC(NPIX) )
 ALLOCATE( LONP(NPIX) )
 ALLOCATE( LATP(NPIX) )
 ALLOCATE( LONG_TABLE(0:LMAX,NPIX) )
 ALLOCATE( ALF(JMAX,NA) ) 
!
!
!
! --- Read pixels and ALFs table
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
    Open(1,file='px-table.dat',status='unknown') 
    Do i=1, 4 
       Read(1,'(a12)')junk
    Enddo
    nwet=0
    ndry=0
    do i=1, npix 
       read(1,*) lonp(i), latp(i), anc(i), k, wet(i) 
       if(wet(i)==1)nwet=nwet+1
       if(wet(i)==0)ndry=ndry+1	
    enddo
    Close(1) 
    Write(*,*) "    - The number of pixels in the px-table is ", npix 
    Write(*,*) "    - The number of WET pixels in the px-table is ", nwet 
    Write(*,*) "    - The number of DRY pixels in the px-table is ", ndry 
    npix = np
!
! --- Reading the ALFs table...  
!
    write(*,*) '    - Reading ALFs & TRIGs from file sh.bin'
    open(3,file='sh.bin',status='unknown',form='unformatted') 
    read(3)ALF
    read(3)LONG_TABLE
    close(3) 
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
!
!
! --- Pre-computing the degree 'm' corresponding to 'J'
!
 Write(*,*) '    - Pre-computing the harmonic order'
 do j=1, jmax 
        mm(j)=mj(j) 
 enddo	
!
!	
! --- Computing the OF SH coefficients- I first evaluate the "Continent 
!     function" (CF) SH coefficients since this involves a smaller number 
!     of pixels. Then I transform them into the OF coefficients taking 
!     into account that the following relationship holds: CF + OF = 1.   
!	
 write(*,*) '    - Building the ocean function coefficients...'
!
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) SHARED(OC,WET,ALF,ANC,LONG_TABLE,MM,NPIX)
 do j=1,jmax 
	OC(J)=0.
 	do i=1,npix 
	   if(wet(i)==0) OC(J) = OC(J) + ALF(J,ANC(I))*CONJG(LONG_TABLE(MM(J),I))
 	enddo
 enddo 
!$OMP END PARALLEL DO
!
!
!
 OC=OC/FLOAT(NP)
!
! --- Writing the result on "sh-of.dat"
!
 open(1,file='shof.dat',status='unknown')  
 do j=1, jmax   
   	if(j==1) write(1,*) j, 1.- real(oc(j)), -aimag(oc(j))
   	if(j/=1) write(1,*) j,   - real(oc(j)), -aimag(oc(j))
 enddo
 close(1) 
!
!
! --- Deallocate memory space
!
 DEALLOCATE( MM, WET, ANC )
 DEALLOCATE( LONP, LATP )
 DEALLOCATE( ALF )
 DEALLOCATE( LONG_TABLE, OC )
!
!
!
 END PROGRAM OF 
!
!
!
