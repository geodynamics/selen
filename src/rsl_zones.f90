!
! This is file "RSL_ZONES.F90"
!
! Last modified GS 04-11-2008 "Intel port"
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Reviewed DM June 2011 - Dynamic memory allocation
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
! -------------------------------------- 
!  Determines user-defined  "RSL zones" 
! -------------------------------------- 
! 
! Input files:
!       - wet.dat 
!	- sh.bin 
! 	- shs.bin
!
!
! Outout files:
!  	- rslzones*.dat 
!  	- lonlat*.dat 
!       - rslzones*ave.dat 
!
!
! INCLUDE "harmonics.f90"
 PROGRAM CLARK
 IMPLICIT NONE
 INCLUDE "data.inc"
 INTEGER, PARAMETER :: NREGIONS = 4
 INTEGER :: NANCH
 INTEGER ICO, DOM, NWET, NPIX, NDATA, IUNIT  
 INTEGER I, J, K, N, IK, MJ, N10, N11, N20, N21
 INTEGER, ALLOCATABLE :: WET(:), ANC(:), MM(:) 
 REAL*8 RSL
 REAL*8, ALLOCATABLE :: LONP(:), LATP(:), SLC(:,:), & 
        PRSL(:,:), AVE(:), ALF(:,:) 
 REAL*8 X10, X11, X20, X21
 COMPLEX*16, ALLOCATABLE :: S(:,:), LONG_TABLE(:,:) 
 CHARACTER*13 CJUNK
!    
!
!
! --- Number of "anchor pixels"
!
  open(1,file='anchor.tmp',status='old')
  read(1,*) nanch
  close(1)
!
! --- Allocate memory space
!
 ALLOCATE( WET(NP), ANC(NP), MM(JMAX) )
 ALLOCATE( LONP(NP), LATP(NP), SLC(NP,0:NN) ) 
 ALLOCATE( PRSL(NP,0:NN), AVE(0:NN), ALF(JMAX,NANCH) )
 ALLOCATE( S(JMAX,0:NN), LONG_TABLE(0:LMAX,NP) )
!
! --- Examining the pixels table & extracting information
!
	Open(1,file='px-table.dat',status='unknown') 
	Do i=1, 4 
		Read(1,'(a13)')cjunk
	Enddo
	nwet=0
	do i=1, 2*np 
		read(1,*,end=10) lonp(i), latp(i), anc(i), k, wet(i) 
		if(wet(i)==1)nwet=nwet+1
	enddo
10      Close(1)
!	Write(*,*)'    - There are ', NWET, ' wet pixels in the pixels table'
        close(12)  	
!
!
! --- Pre-computing the degree 'm' corresponding to 'J'
!
! 	Write(*,*) '    - rsl_zones.f: Pre-computing the harmonic order'
 	do j=1, jmax 
        	mm(j)=mj(j)
 	enddo	
!
!
!------- Virtual RSL sites are used to determine the shape of RSL zones 
!	 are the pixels themselves, where the SH's are already known... 
!
! --- Reading the ALFs table...  
! 	write(*,*) '    - rsl_zones.f: reading ALFs & TRIGs from file <<sh.bin>>'
 	open(3,file='sh.bin',status='unknown',form='unformatted') 
		read(3)ALF
		read(3)LONG_TABLE
 	close(3) 
 	DO J=1, JMAX 
		ALF(J,:)=ALF(J,:)*(2-DOM(J))
 	ENDDO
!
!
!------- Reading the sealevel CSH coefficients at any virtual RSL site...', filename 
!
!	Write(*,*) '    - rsl_zones.f90: reading the sea level CSH coefficients'
	open(3,file='shs.bin',status='unknown',form='unformatted') 
	read(3) S
	close(3)
!
	Write(*,*) '    - Computing sea level change at WET pixels'
!
	do i=1, np 
!
	If(np>=10000) then 
        	      If(mod(i,50000)==0) Write(*,*) "    - ", i, " pixels of ", np	
        elseif(np>=1000.and.np<10000) then 
        	      If(mod(i,5000)==0)  Write(*,*) "    - ", i, " pixels of ", np	
	              Endif     	     	
!	
	slc(i,:) = 0.  
	if(wet(i)==1) then 
		do j=1, jmax
		slc(i,:) = slc(i,:) + ALF(J,ANC(I))*real(s(j,:)*long_table(mm(j),i)) 
		enddo
	endif
	enddo
	close(22)
!
!
! ----- RSL
! Write(*,*) '    - rsl_zones.f90: computing RSL at the pixels...'
 do i=1, NP 
 if(wet(i)==1)then 
! ---- k is kyrs BP
 	do k=0, NN   
		prsl(i,k)= -(slc(i,nn)-slc(i,nn-k))	
 	enddo
 endif
 enddo
!
!
!
!
!
! ------ Files for RSL histories ...
!
 open(10,file='rslzones-10.dat',status='unknown') 
 open(11,file='rslzones-11.dat',status='unknown') 
 open(20,file='rslzones-20.dat',status='unknown') 
 open(21,file='rslzones-21.dat',status='unknown') 
!
!
! ------ Files for RSL sites coordinates 
!
 open(50,file='lonlat-10.dat',status='unknown') 
 open(51,file='lonlat-11.dat',status='unknown') 
 open(60,file='lonlat-20.dat',status='unknown') 
 open(61,file='lonlat-21.dat',status='unknown') 
!
! ----- Counters 
! 
 n10=0; n11=0  
 n20=0; n21=0 
!
!
! ----- Loop on the pixels ("RSL sites") 
!
 DO 2000 I=1, NP 
!
 If(wet(i)==1) then 
!
 ico=0
 	do k=1, NN-1 
 		if(prsl(i,k)*prsl(i,k+1).le.0.)  ico=ico+1
 	enddo
!
!
 if(prsl(i,NN).le.0.) THEN      ! RSL was below present datum at the LGM 
!
!
 if(ico==0) then 
 	iunit=10
	n10=n10+1    ! No zero crossing 
 endif
 if(ico>=1) then 
 	iunit=11
 	n11=n11+1    ! One or more zero crossings 
 endif	
!
!
 ELSE      		       ! RSL was above present datum at the LGM 
!
!			  
 if(ico==0) then 
 	iunit=20
 	n20=n20+1    ! No zero crossing 
 endif
 if(ico>=1) then 
 	iunit=21
 	n21=n21+1    ! One or more zero crossing 
 Endif
!
!		   
ENDIF
!
 write(iunit+40,*) lonp(i), latp(i)
!
 write(iunit,'(a1,f10.4,f10.4,i4,i4)') ">", lonp(i), latp(i), iunit, ico 
 do k=0, NN
	write(iunit,*)    k, prsl(i,k) 
 enddo	
!
  Endif
!
  2000 CONTINUE 
!
!
!  Write(*,*) '    - rsl_zones.f90: Computed ', n10+n11+n20+n21,' RSL curves' 
  Write(*,*) '    - Number of wet pixels: ', NWET
  Write(*,*) '    - Statistics of the RSL zones'
  x10=float(n10)/float(nwet)*100.
  x11=float(n11)/float(nwet)*100.
  x20=float(n20)/float(nwet)*100.
  x21=float(n21)/float(nwet)*100.
  Write(*,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~'
  Write(*,*) '    - type, size(px), freq(%) '
  Write(*,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~'
  Write(*,3) '    - 10', n10, x10
  Write(*,3) '    - 11', n11, x11
  Write(*,3) '    - 20', n20, x20
  Write(*,3) '    - 21', n21, x21
  Write(*,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~'
  Write(*,3) '        ', n10+n11+n20+n21, x10+x11+x20+x21
  Write(*,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~'
3 Format (a11,i8,2x,f6.2)
!
  close(10); close(11) 
  close(20); close(21) 
  close(50); close(51) 
  close(60); close(61) 
!
!
! An average RSL curve for each region 
!
  do 111 i=1, nregions 
!
!  
  if(i==1)     then 
  			open(13,file='rslzones-10.dat',status='unknown') 
  			open(29,file='rslzones-10-ave.dat',status='unknown') 			
			ndata=n10 
  elseif(i==2) then 
  			open(13,file='rslzones-11.dat',status='unknown')  
  			open(29,file='rslzones-11-ave.dat',status='unknown') 			
			ndata=n11 
  elseif(i==3) then     
  			open(13,file='rslzones-20.dat',status='unknown') 
  			open(29,file='rslzones-20-ave.dat',status='unknown') 			
			ndata=n20 
  elseif(i==4) then 
  			open(13,file='rslzones-21.dat',status='unknown')
  			open(29,file='rslzones-21-ave.dat',status='unknown') 			
			ndata=n21
!
  endif 
!
  ave(:)=0.0
!
  do 222 n=1, ndata
	read(13,'(a13)') cjunk  
!	
	do k=0, nn 
        	read(13,*) ik, rsl
		ave(ik)=ave(ik)+rsl 
	enddo
!  
  222 continue 
  close(13) 
!
  do k=0, nn 
  		  If(ndata/=0) then 
                  	         write(29,*) k, ave(k)/float(ndata)  
		  	       else  
                  	         write(29,*) k, ' 0.'
			       endif
  enddo
!
  close(29) 
!  
!
111 continue 
!
! --- Deallocate memory space
!
 DEALLOCATE( WET, ANC, MM )
 DEALLOCATE( LONP, LATP, SLC ) 
 DEALLOCATE( PRSL, AVE, ALF )
 DEALLOCATE( S, LONG_TABLE )
!
  END PROGRAM CLARK
!
!
!
