!
! This is program "RSLC.F90"  
!
! Last modified GS 04-29-2008 "Intel port" & ALPS new paper ...
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
! ----------------------------------------------------------------
!  Computes the RSL curves at a number of <<virtual RSL sites>> 
! ----------------------------------------------------------------
! 
! Input files:
!	- "lonlat_rslc.dat"
! 	- shs.bin
!       - shu.bin 
! 	- shrslc.bin
!
!
! Output files:
!	- 'rslc_cont.dat' 
!
!
! INCLUDE "harmonics.f90"
 PROGRAM SL
 IMPLICIT NONE
 INCLUDE "data.inc"
 CHARACTER*20, PARAMETER :: FILEOUT1 = 'rslc-cont.dat'
 CHARACTER*20, PARAMETER :: FILEOUT2 = 'rslc.dat'
 CHARACTER*20 DATE, TIMC                
 CHARACTER*22 FILENAME 
 CHARACTER*42 JUNK
 COMPLEX*16, ALLOCATABLE :: YY(:,:), S(:,:), U(:,:), N(:,:) 
 REAL*8, ALLOCATABLE :: LONS(:), LATS(:), SLC(:,:)
 REAL*8 LON, LAT, RSL, TIME
 INTEGER I, J, K, DOM 
!    
!
!
!
!
! --- Allocate dynamic arrays
 ALLOCATE( YY(JMAX,NRSLC), S(JMAX,0:NN), U(JMAX,0:NN), N(JMAX,0:NN) )
 ALLOCATE( LONS(NRSLC), LATS(NRSLC), SLC(NRSLC,0:NN) )
!
! --- Reads the coordinates of RSL "sites"
        OPEN(1,FILE=RSLC_FILE,STATUS='unknown')
!Write(*,*) & 
!"    - rslc.f: reading the <<RSL sites>> co-ordinates from file ", rslc_file
        do i=1, NRSLC 
   		read(1,*) lons(i), lats(i) 		   
        enddo
	close(1) 
!
! --- rslc.f: reading the CSHs on file ', filename 
	filename='shrslc.bin'
!Write(*,*) & 
!"    - rslc.f: reading the SH at the <<RSL sites>> from file ", filename
	open(7,file=filename,status='unknown',form='unformatted') 
	read(7) yy
	close(7)
	do j=1,jmax 
		yy(j,:)=yy(j,:)*(2-dom(j))
	enddo
!
!
! --- rslc.f: reading the sealevel CSH coefficients ', filename 
	filename='shs.bin'
!Write(*,*) & 
!"    - rslc.f: reading the sealevel SH coefficients from file ", filename
	open(3,file=filename,status='unknown',form='unformatted') 
	read(3) S
	close(3)
!
!
! --- rslc.f: reading the vert. disp. CSH coefficients ', filename 
	filename='shu.bin'
!Write(*,*) & 
!"    - rslc.f: reading the vert. disp. SH coefficients from file ", filename
	open(3,file=filename,status='unknown',form='unformatted') 
	read(3) U
	close(3)
!
!
! --- rslc.f: reading the geoid undulation CSH coefficients ', filename 
	filename='shn.bin'
!Write(*,*) & 
!"    - rslc.f: reading the geoid SH coefficients from file ", filename
	open(3,file=filename,status='unknown',form='unformatted') 
	read(3) N
	close(3)
!
!
! --- rslc.f: sealevel change at the sites...'
  	Write(*,*) "    - Computing sealevel change at virtual RSL sites in file ", & 
	trim(adjustl(fileout1))
	do i=1, nrslc
		slc(i,:) = 0. 
	do j=1, jmax
		slc(i,:) = slc(i,:) + real(s(j,:)*yy(j,i)) 
	enddo
	enddo
!
!
! --- rslc.f: relative sea level change at the sites...' 
	Open(11,file=FILEOUT1,status='unknown')
!Write(*,*) & 
!"    - rslc.f: Computing RSL change at the <<RSL sites>> from file ", fileout1  
 	do i=1, nrslc
	Write(11,'(a16,i5,a1,i5,5x,a10,2(f10.4,1x))') & 	
	'RSL <<site>>    ', i, '/', nrslc, "  lon-lat:", lons(i), lats(i)
		do k=0, nn
		     rsl = -(slc(i,nn)-slc(i,nn-k))
		     Write(11,'(f10.4,1x,f10.4)') float(k), rsl
		enddo
	enddo
	Close(11) 
!
!
! --- rscl.f: filter the RSL values corresponding to the desired time BP 
	Open(3, file=FILEOUT2,status='unknown') 	
!
!     Header for FILEOUT2
	Write(3,*) '01 Relative- Sea Level at time ', time_bpc 
 	Write(3,*) '02 File ', trim(adjustl(FILEOUT1))	
	Write(3,*) '03 Run =', trim(adjustl(run))
	Write(3,*) '04 Lmax =', LMAX
	Write(3,*) '05 Res =', RES
	Write(3,*) '06 Ice model =', trim(adjustl(ICE_MODEL))
	Write(3,*) '07 NV & CODE =',NV, CDE 
	Write(3,*) '08 Iterations =', SMAX 
	Write(3,*) '09 Mode of solution: ', IMODE
!
	Open(11,file=FILEOUT1,status='unknown')
	do 10 i=1, nrslc 
	read(11,'(a42,2(f10.4,1x))') junk, lon, lat 
		do 20 j=0, nn 
		read(11,*) time, rsl 
		if(time==time_bpc) then 
				   Write(3,'(3(f10.4,1x))') lon, lat, rsl 
				   Endif		 
20		continue 
10      continue 
        close(11)
        close(3)	
!
 DEALLOCATE( YY, S, U, N )
 DEALLOCATE( LONS, LATS, SLC )
!
 END PROGRAM SL
!
!
!
