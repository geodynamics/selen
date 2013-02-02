!
! This is file "WNW.F90" 
!
! Last modified GS 04-11-2008 "Intel port"
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
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
! --------------------------------------------
! Computes the window function to degree jmax 
! --------------------------------------------
!
!
! Input files: 
!	- sh.bin
!
! Output files: 
!	- wnw.dat
!
!
! INCLUDE "harmonics.f90"
 PROGRAM WINDOW 
 IMPLICIT NONE 
 INCLUDE "data.inc"
 INTEGER :: NANCH
 INTEGER I, J, MJ
 INTEGER, ALLOCATABLE :: MM(:), ANC(:) 
 INTEGER NB1, N1_5, N5_10, N10_100, NA100 
 REAL*8 XJUNK
 REAL*8, ALLOCATABLE :: ALF(:,:), ERR(:) 
 REAL*8 ERROR, PERR
 COMPLEX*16, ALLOCATABLE :: PROD(:), LONG_TABLE(:,:) 
 CHARACTER*30 JUNK 
!
!
!
! --- Number of anchor pixels
!
  open(1,file='anchor.tmp',status='old')
  read(1,*) nanch
  close(1)
!
! --- Allocate memory
!
  allocate( mm(jmax), anc(np) )
  allocate( ALF(JMAX,NANCH) )
  allocate( ERR(JMAX) )
  allocate( prod(jmax) )
  allocate( long_table(0:lmax,np) )
!
! --- Pre-computing the degree 'm' corresponding to 'J'
!
!Write(*,*) '    - wnw.f: Pre-computing the harmonic order'
 	do j=1, jmax 
 	       mm(j)=mj(j) 
 	enddo	
!
!
! --- Reading the ALFs table...  
!write(*,*) '    - Reading ALFs & TRIGs from file sh.bin'
 	open(3,file='sh.bin',status='unknown',form='unformatted') 
		read(3)ALF
		read(3)LONG_TABLE
 	close(3) 
!
!
! --- Reading the pixels table ...
 	open(2,file='px-table.dat',status='unknown')
 	do i=1, 4 
    		read(2,'(a30)') JUNK 	
 	enddo 
 	do i=1, np 
 		read(2,*) xjunk, xjunk, anc(i)
	enddo
	close(2) 
!
!
! --- Computing the "error-per-degree" 
	Write(*,*) "    - Computing the error-per-degree"
 	open(11,file='wnw.dat',status='unknown') 
!
 	do j=1, jmax 
!
	If(jmax>=10000) then 
        	        If(mod(jmax,10000)==0) Write(*,*) "    - wnw.f: ", j, " harmonics of ", jmax	
        elseif(jmax>=1000.and.jmax<10000) then 
        	        If(mod(jmax,1000)==0)  Write(*,*) "    - wnw.f: ", j, " harmonics of ", jmax 	
	                Endif
!		    
		err(j)=0. 
		prod(j)=(0.,0.) 
			do i=1, np 
			   prod(j)=prod(j) + & 
			   ALF(J,ANC(I))**2*CONJG(LONG_TABLE(MM(J),I))*LONG_TABLE(MM(J),I)
			enddo
		err(j)= abs(prod(j)-float(np))/float(np)*100. 
		write(11,*) j, err(j) 
		error = error + err(j) 
 	enddo 
 	close(11)
!
     nb1=0 
    n1_5=0 
   n5_10=0
 n10_100=0
   na100=0
!
 do j=1, jmax
 	perr=err(j)  
 	if(perr<1.) 			then 
			nb1=nb1+1
	elseif(perr>=1.and.perr<5.) 	then
		        n1_5=n1_5+1
	elseif(perr>=5.and.perr<10.) 	then
			n5_10=n5_10+1
	elseif(perr>=10.and.perr<100.) 	then
			n10_100=n10_100+1
	elseif(perr>=100.) 		then
			na100=na100+1
	endif			
 enddo			
!   
 Write(*,*) '    - There are ', JMAX, ' harmonics on file <<sh.bin>>'
 Write(*,*) '      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
 Write(*,*) '              Error table for surface integrals            '
 Write(*,*) '      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
 Write(*,*) '    - ------- Error < 1%   ---> ', nb1, '/', jmax  			
 Write(*,*) '    - 1%   <= Error < 5%   ---> ', n1_5, '/', jmax  			
 Write(*,*) '    - 5%   <= Error < 10%  ---> ', n5_10, '/', jmax  			
 Write(*,*) '    - 10%  <= Error < 100% ---> ', n10_100, '/', jmax  			
 Write(*,*) '    - 100% <= Error ----------> ', na100, '/', jmax 
 Write(*,*) '      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' 
 Write(*,'(A42,1x,F10.4,1X,A1)') '    - Average relative error per degree =', error/float(jmax), '%' 
 Write(*,*) '      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' 			
 if(error/float(jmax)>=1.) write(*,*) '    +++> WARNING: the average error exceeds 1%'
 if(error/float(jmax)>=5.) write(*,*) '    +++> WARNING: the average error exceeds 5%'
!
! --- Release memory
!
  deallocate( mm, anc )
  deallocate( ALF )
  deallocate( ERR )
  deallocate( prod )
  deallocate( long_table )
!  
end program WINDOW 
!
!
!
