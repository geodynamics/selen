!
! This is program "SH.F90" 
! 
! Last change: GS April 11, 2008 "INTEL port" 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!)
! *** Reviewed GS & FC November 2009 - Porting under gfortran  
! *** Reviewed DM June 2011 - Dynamic arrays  
! *** Modified DM June 2011 - Parallel execution
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
! WARRANTY; without even the implied warranty of MERCHANTABlILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! --------------------------------------------------------------------------- 
! This program computes the ALFs at the latitudes of the "anchor" pixels, and 
! the TRIGs at pixels longitudes for degrees from 0 to LMAX. ALFs and TRIGs 
! are written on the binary, unformatted file 'sh.bin'. ### GS 15-09-2007 ###
! --------------------------------------------------------------------------- 
!
! Input files:
!	- px-lat.dat 
!       - px-table.dat 
!
! Output files:
! 	- sh.bin 
!
!
! INCLUDE "harmonics.f90"
 PROGRAM SH
 IMPLICIT NONE 
 INCLUDE "data.inc"
 INTEGER, PARAMETER :: A_LARGE_INTEGER=10**6
 CHARACTER*12 JUNK 
 INTEGER I, J, L, NA, IJUNK, NPX
 REAL*8 LATA, LONG
 REAL*8, ALLOCATABLE :: ALF(:,:) !ALF(JMAX,NANCH)
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:) ! LONG_TABLE(0:LMAX,NP)
 REAL*8 COSDD, SINDD
 REAL*8, ALLOCATABLE :: LAT(:), LON(:)
!
!
!
!
! --- Output file for LAFs and Trigs 
  open(3,file='sh.bin',status='unknown',form='unformatted') 
!
  write(*,*)'    - Maximum degree is: ', lmax 
  write(*,*)'    - JMAX is:', (lmax+1)*(lmax+2)/2
  write(*,*)'    - Resolution is: ', res 
!
! --- Getting info about the number of anchor pixels 
!
  open(1,file='anchor.tmp',status='old')
  read(1,*) na
  close(1)
  Write(*,*) "    - Found ", na, " anchor pixels in file px-lat.dat"
  npx = np
!
!
! --- Allocate memory space
  allocate( alf(jmax,na) )
  allocate( long_table(0:lmax,npx) )
  allocate( lat(na) )
  allocate( lon(npx) )
! 
!
! --- Reading pixels
!
  open(2,file='px-lat.dat',status='unknown') 
  do i=1, 4 
  	read(2,*) JUNK 
  enddo
!  
  do i=1, na
    read(2,*) ijunk, lat(i)
  end do
  close(2)
!
  open(2,file='px-table.dat',status='unknown') 
  do i=1, 4 
  	read(2,*) JUNK 
  enddo
!
  do i=1,np
    read(2,*) lon(i)
  end do
  close(2)
!
!
!
!
! --- Computing the ALF at the latitudes of anchor pixels 
!
  Write(*,*) "    - Pre-computing ALFs at latitudes of anchor pixels"
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! CAUTION! PLMBAR_MOD pre-computes some expressions each time it is 
! called with a new lmax and stores them in a SAVEd array for reuse 
! in subsequent calls. This is thread-safe IF and ONLY IF (1) plmbar_mod
! is called BEFORE the parallel region and (2) all the calls in the
! parallel region use the same lmax.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
call PLMBAR_MOD( lmax, lat(1), ALF(:,1) )
!
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NA,ALF,LAT)
!$OMP DO SCHEDULE(GUIDED)
  do i=2, na   
!
!	If(nanch>=10000) then 
!        	      If(mod(i,50000)==0) Write(*,*) "    - ", i, " pixels of ", na	
!        elseif(nanch>=1000.and.nanch<10000) then 
!        	      If(mod(i,5000)==0)  Write(*,*) "    - ", i, " pixels of ", na	!
!	              Endif
!
!       read(2,*,end=9) ijunk, lata    
       lata = lat(i)
       call PLMBAR_MOD (lmax, lata, ALF(:,i))        
  enddo
!$OMP END PARALLEL
  write(3) ALF
!
!
! --- Computing the TRIGs at pixels  
!
  Write(*,*) "    - Pre-computing TRIG functions at pixels"
!
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(LONG_TABLE,LON,NPX) SCHEDULE(GUIDED)
  do i=1, npx  
    long = lon(i)
  	do l=0, lmax
	   long_table(l,i)=cmplx(cosdd(l*long),sindd(l*long))
  	enddo	
  enddo
!$OMP END PARALLEL DO
!  close(2) 
  write(3) long_table
!
  close(3)
!
   write(*,*) '    - ALFs & TRIGs are written on file sh.bin'
! 
!
! --- Deallocate memory space
  deallocate( alf )
  deallocate( long_table )
  deallocate( lon, lat )
!
!
!
  END PROGRAM SH
!
!
!
