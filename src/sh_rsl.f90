!
! This is program "SH_RSL.F90" 
!
! Last modified GS 04-11-2008  Intel port
! Updated on August 2008 for v. 2.7 
! Re-touched on August & September 2008 for version 2.7 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! Reviewed GS August 2009 -   
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
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
!  Computes the 4pi-normalized, complex spherical harmonics at
!  Relative Sea Level sites listed in file RSL_DATABASE - 
!
! Input files:
!	- "sealevel.dat"
! Output files: 
!	- shrsl.bin
!
!
! INCLUDE "harmonics.f90"
 PROGRAM S
 IMPLICIT NONE 
 INCLUDE "data.inc"
 CHARACTER*22 NAME, CJUNK 
 INTEGER I, J, NDATA, CODE(NRSL) 
 REAL*8 AJUNK, LONS(NRSL), LATS(NRSL)
 COMPLEX*16 Y(JMAX,NRSL) 
 CHARACTER*10  LATSC10, LONSC10 
 CHARACTER*200 LINEP
 CHARACTER*100 SS(2)
 INTEGER NOUT
        
!
!
! ========================================================================
! New of August 13, 2008 (SELEN v. 2.7)
!
! This program is designed to read RSL databasdes of two possible formats:
!
! 1) The usual "sealevel.dat" archive of Peltier and Tushingham 
! 2) An italian Holocene database from Antonoili et al. (Quat. Int., 2008) 
! ========================================================================
!
!
!
!------ Opening the selected database 
!
! ****************************************
   If    (rsl_database_format=='0') then 
! ****************************************
!
   	OPEN(1,FILE=RSL_DATABASE,STATUS='unknown')
!
   	i=1
2  	READ (1,6,END=1) code(i), lats(i), lons(i), ndata, name 
   	IF(lons(i)<=0.) lons(i)=360.+lons(i)
   	do j=1, ndata
  	   	READ (1 ,*) ajunk
   	enddo
   	i=i+1
   	IF(i<=nrsl) GOTO 2
1  	CLOSE(1)
6  	FORMAT(1X,I3,1X,F5.1,1X,F6.1,1X,I2,1X,A22)
!
! ****************************************
   Elseif (rsl_database_format=='1') then 
! ****************************************
!
   	OPEN(1,FILE=RSL_DATABASE,STATUS='unknown')
!  
		Do i=1, nrsl 
			read(1,'(a22)')cjunk 
			read(1,'(a22)')cjunk 
			read(1,'(a22)')cjunk 
			read(1,*) lats(i), lons(i)  		
			read(1,'(a22)')cjunk
			read(1,'(a22)')cjunk 	
		Enddo
  	Close(1)
!
! ****************************************
   Elseif (rsl_database_format=='2') then 
! ****************************************
! 
   	OPEN(1,FILE=RSL_DATABASE,STATUS='unknown')
!  
		Do i=1, nrsl 
			read(1,'(a22)')cjunk 
			read(1,'(a22)')cjunk 
			read(1,'(a22)')cjunk 
!
			read(1,'(a200)') linep 	
!		
			call scan_string (linep, 2, ss, nout)				
!
			LATSC10=trim(adjustl(ss(1))) 
			LONSC10=trim(adjustl(ss(2)))			
!
 			call FROMDDMMSS_2_DEGREES (LATSC10, LATS(I)) 
 			call FROMDDMMSS_2_DEGREES (LONSC10, LONS(I)) 
!			
!read(1,*) lats(i), lons(i)  		
!
			read(1,'(a22)')cjunk
			read(1,'(a22)')cjunk 	
		Enddo
  	Close(1)
!
   Endif     
!
!
!
  write(*,*) '    - Computing the harmonics at the RSL sites'
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! CAUTION! PLMBAR_MOD pre-computes some expressions each time it is 
! called with a new lmax and stores them in a SAVEd array for reuse 
! in subsequent calls. This is thread-safe IF and ONLY IF (1) plmbar_mod
! is called BEFORE the parallel region and (2) all the calls in the
! parallel region use the same lmax. Sbr. HARMO calls PLMBAR_MOD.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
 call harmo(lmax, lons(1), lats(1), y(:,1))
!   
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(lons,lats,y) &
!$OMP    SCHEDULE(GUIDED)
 do i=2, nrsl    
!If(mod(i,25)==0) write(*,*) '    - sh_rsl.f:', i, 'RSL sites of', nrsl
       call harmo(lmax, lons(i), lats(i), y(:,i)) 
  enddo
!$OMP END PARALLEL DO
!
!
!
!write(*,*) "    - sh_rsl.f: Writing the harmonics on file 'shrsl.bin'"
  open(7,file='shrsl.bin',status='unknown',form='unformatted') 
  write(7) y; close(7) 
!
!
!
 END PROGRAM S
!
!
!
