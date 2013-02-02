!
! OF_DV.F90 
! 
! Last change: GS April 11, 2008 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
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
! --------------------------------------------------------------------------- 
! This program computes the degree variance (DV) of the ocean function (OF), 
! defined as: DV(OF)(at degree 'l') = \sum_{m=0,l}(2-\delta_{0m})|a_{lm}|^2 
! where \delta is Kronecker delta, and a_{lm} is the (4-pi) normalized degree 
! 'l' order 'm' SH coefficient of the OF expansion. Multiply DV by 4*pi/(2l+1) 
! in order to get the (squared) spectrum S^2(l) at degree 'l'- GS Feb 26 08-  
! --------------------------------------------------------------------------- 
!
! Input files:
!	- shof.dat 
!
! Output files:
! 	- ofdv.dat 
!       - dv-power.tmp
!
!
! INCLUDE "harmonics.f90"
 PROGRAM OFDV
 IMPLICIT NONE
 INCLUDE "data.inc"
!
 INTEGER, PARAMETER :: L_MIN=1, L_MAX=LMAX 
 INTEGER I, J, K, L, M, LJ, MJ, NDATA 
 REAL*8 REOC, IMOC, C2S2, DV(0:LMAX)
 REAL*8 SUMTT, SUMTY, SUMT, SUMY, B, P
 REAL*8 Y(L_MIN:L_MAX), T(L_MIN:L_MAX)
 REAL*8 C(0:LMAX,0:LMAX), S(0:LMAX,0:LMAX) 
!
!
!
! --- Reads the fully normalized coefficients
!
 open(1,file='shof.dat',status='unknown')  
 do j=1, jmax  
 	read(1,*) k, reoc, imoc 
	l=lj(j) 
	m=mj(j)	
	c(l,m)=reoc 
	s(l,m)=imoc
!write(*,*) l, m, c(l,m), s(l,m) 	
 enddo
 close(1)
!
!
! --- Computes the degree variance 
!
      do l=l_min,l_max
      dv(l)=0
      do m=0,l
           c2s2=c(l,m)**2+s(l,m)**2 
      	   if(m==0) dv(l)=dv(l)+    c2s2
      	   if(m/=0) dv(l)=dv(l)+ 2.*c2s2  
      enddo
!
	    t(l)=log10(float(l))
	    y(l)=log10(dble(dv(l)))
!	    	    
      enddo
      close(2) 
!
!
! --- Linear interpolation in the range [lmin:lmax] 
!
	sumt=0.
	sumy=0.
	sumty=0.
	sumtt=0.
	do l=l_min, l_max
		sumt=sumt+t(l)
		sumy=sumy+y(l)
		sumty=sumty+t(l)*y(l) 
		sumtt=sumtt+t(l)*t(l) 		
	enddo
!
	ndata=l_max-l_min+1	
!
!
! --- Best fit = (10^b)*degree**p 
!
! --- Intercept (b)
	b=(sumy*sumtt-sumt*sumty)/(float(ndata)*sumtt-sumt**2)
!
! --- Steepness (p) 		 	
	p= (float(ndata)*sumty-sumt*sumy)/(float(ndata)*sumtt-sumt**2)
!
! --- Writing data on <<ofdv.dat>>
      open (20,file='ofdv.dat',status='unknown')
      Write(20,*) '--------------------------------------'
      Write(20,*) '    deg      DV         Interp. DV '	
      Write(20,*) '--------------------------------------'
      do l=l_min,l_max
      		write(20,'(4x,i3,4(4x,E10.4))') l, dv(l), 10**b*float(l)**p  
      enddo
      close(20) 
!
! --- This string is later used by "ofdv.gmt"
	open (4,file='dv-power.tmp',status='unknown')  
	write(4,'(a26,f6.3)') "20 3e-2 16 0 2 BL Power= ", p
	close(4)  
!
!
 END PROGRAM OFDV 
!
!
!
