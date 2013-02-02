!
! STOKES.F90 
!
! Last modified GS 04-11-2008 [Intel Port]
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! July 2010: The routine is modified significantly. It now provides the
!            Stokes coefficients in a fully-normalized form, following  the
!            GRACE conventions (fully normalized with NO CS phase)- This has
!            been done after a discussion with Luoise and Jens about the GIA
!            correction on the GRACE coefficients. No factorial is now needed,
!            and the limit on the maximum harmonic degree is removed. 
! ***  Revised GS August 2010 - g95 - Double precision implementation
! === Revised GS & FC May 21 2011 - DELTA parameter in ice history     
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
! Computes the present-day rate of variation of the Stokes coefficients
! in a range of degrees specified in "data.inc". Revised JULY 2010. 
!
! Input files:
!	- shn.bin
!
! Output files:
! 	- stokes.dat
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! INCLUDE "harmonics.f90"
 PROGRAM STOKES 
 IMPLICIT NONE 
 INCLUDE "data.inc"
 INTEGER J_INDEX, I, L, M, DOM             ! degree, order, and delta_{0m}
 REAL*4,  PARAMETER :: ERADIUS=6.371E6     ! Radius of the earth, m
 REAL*8 CDOT, SDOT                         ! dot(c_lm, s_lm) 
 REAL*8 PHILM
 COMPLEX*16 N(JMAX,0:NN)    	           ! 4-pi normalized, complex, geoid coefficients 
 COMPLEX*16 RATEN                           ! Complex SC
 CHARACTER*20 date, timc                   ! date and time
 CHARACTER*100 RIGA		           ! A string
!
 REAL*8, PARAMETER :: DDELTA=DELTA 
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!--- Output file for SC  
  	open(1,file='stokes.dat',status='unknown') 
	
! --- Header information about the SELEN settings 
        call DATE_AND_TIME (date,timc)      
        Write(1,*) date(1:4), '.', date(5:6), '.', date(7:8), & 
	          ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)
	Write(1,*) "Ice model: ", trim(adjustl(ice_model))
	Write(1,*) "Number of mantle layers: ", NV 
	Write(1,*) "Model code (see TABOO User guide): ", CDE
	Write(1,*) "Viscosity model: ", visco_model
	Write(1,*) "SLE iterations: ", SMAX 
	Write(1,*) "SLE mode of solution: ", IMODE
	Write(1,*) "Maximum harmonic degree: ", LMAX 
	Write(1,*) "Tegmark resolution: ", RES 
!	
	Write(1,*) ' ' 
	Write(1,*) ' ---------------------------------------------------------------- '
	Write(1,*) '  degree, degree, order,      Fully normalized coefficients       '
	Write(1,*) '     j       l      m         c_lm and s_lm (x 1E-11 yr^-1)       '   
	Write(1,*) ' ---------------------------------------------------------------- '
	Write(1,*) ' ' 
!
!
!--- Reading the 4pi-normalized geoid coefficients...
	open(3,file='shn.bin',status='unknown',form='unformatted') 
        read(3) N
	close(3)
!
!--- Loop on the degrees and orders... 
  		do l=stmin, stmax 
  		do m=0, l 
!
!--- rate of geoid change in  mm/yr (or m/kyr)  
  		RATEN = (N(j_index(l,m),NN)-N(j_index(l,m),NN-1))/DDELTA
!   
!--- now in (yrs)^(-1)
		RATEN = RATEN/(ERADIUS*1E3)  
!
!--- now in units of 10^-11               	       
   		RATEN = 1.E11*RATEN
!
!
! OLD piece of Fortran script 
! ===========================
!
! > --- normalization                           	       
! >	RATEN = RATEN* & 
! >	(2.-dom(j_index(l,m)))*sqrt(2.*float(l)+1.)*sqrt(fct(l-m))/sqrt(fct(l+m))       
! >	cdot=   real (RATEN)
! >	sdot= - aimag(RATEN)      
!
!
! NEW piece of Fortran script 
! ===========================
!
	        philm= ((-1)**m)*sqrt(2.-dom(j_index(l,m)))

	        cdot=    real (RATEN)* philm
	        sdot= - aimag (RATEN)* philm
!
  		write(1,'(1x,i6,4x,i4,3x,i4,8x,2(E12.6,6x))') & 
		j_index(l,m), l, m, cdot, sdot 
!
  	enddo
        enddo
!
 close(1) 
!
 end PROGRAM STOKES 
!
