!
! This is program "PX.F90"   
!
! Last modified GS 04-11-2008 [Intel port]
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Reviewed DM June 2011 - Dynamic memory allocation, QuickSort
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
!
!							 
! --------------------------------------------------------------------------- 
! Given resolution "R", this program generates  N_p = 40R(R-1)+12 pixels on 
! the unit sphere, using the algorithm by Max Tegmark ('An icosahedron-based 
! method for pixelizing the celestial sphere', ApJ Letters, 470, L81).   
! ---------------------------------------------------------------------------        					   !
!
! February 2008- The code has been modified to improve performance in pixels
!		 manipulation. Added 
!                                                                         
!   Input file:   - "data.inc"       	
!				    									  !
!   Output files: - "px.dat"
! 	  	  - "px-lat.dat"
!		  - "pxa.dat" 	
!                 - "anchor.tmp"   NEW number of anchor pixels
!                                                                          
!
MODULE COMMON
IMPLICIT NONE 
INTEGER, PARAMETER :: RES_MAX=5000                         ! max admissible resolution
INTEGER, PARAMETER :: N_MAX=2*RES_MAX*(RES_MAX-1)*20+12    ! max number of pixels
REAL*8  R(0:19,3,3), V(0:11,3)  			   ! rotation matrices & corner vectors 
END MODULE COMMON 
!
!
!
!

!INCLUDE "harmonics.f90"
PROGRAM PIX
!
  USE COMMON
  IMPLICIT NONE
  INCLUDE "data.inc"
  INTEGER PIXEL      			      ! pixel label		      
  INTEGER, PARAMETER :: N=2*RES*(RES-1)*20+12 ! current number of pixel        
  INTEGER I, J, NL
  INTEGER, ALLOCATABLE :: IMIN(:), IMAX(:)    ! indices   
  INTEGER, ALLOCATABLE :: IANCH(:)            ! index of anchor-point
  REAL*8 LONP, LATP
  REAL*8, ALLOCATABLE :: LON(:), LAT(:)   ! Pixels lon-lat 
  REAL*8 VECT(3)     			          ! unit vector of the pixel 
  REAL*8 X,Y,Z       			          ! pixel Cartesian coordinates 
  REAL*8, ALLOCATABLE :: LATEFF(:)        ! Anchor latitudes 
  REAL*8, PARAMETER :: EPS = 1.0e-8       ! Tolerance in anchor latitude identification
!
!
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! === Step 0: Checking some bounds... 
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
    IF(res > res_max) then
                      write(*, *) '    - The resolution exceeds the max allowed'
                      WRITE(*, *) '    - ************ JOB ABORTED *************'
          	      call Stop_config 	
		      Stop
    ENDIF
!
    IF(res < 10     ) then
                      write(*, *) '    - The resolution must be >= 10'
                      WRITE(*, *) '    - P******* JOB ABORTED ********'
          	      call Stop_config 	
		      Stop		      
    ENDIF
!
    write(*, *) '    - The resolution is:', res
    write(*, *) '    - The number of pixels is:', n
!
    allocate( IMIN(N), IMAX(N), IANCH(N) )   
    allocate( LON(N), LAT(N) )
    allocate( LATEFF(N) )
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! === Step 1: Pixels coordinates are computed using Max Tegmark's code. 
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! --- Rotation matrices and corners
!
	call compute_matrices (R)
	call compute_corners  (V)
!
        OPEN(1,FILE='px.dat',STATUS='unknown')
!
        write(*,*) '    - Computing the pixels coordinates...'
!
        do 10 pixel=0,n-1
	  	call pixel2vector (pixel,res,R,v,vect)
	        x=vect(1) ; y=vect(2) ; z=vect(3)
! 
! --- Polar pixel 
          	if(x==0..and.y==0.) then 
          		lonp=0.
                	IF(z>=0.) latp = +90.
                	IF(z<=0.) latp = -90.
! --- Ordinary pixel 
	  		else                         	  
	  		lonp = atan2(y,x)*180./pi
          		if (lonp < 0.) lonp=360.+ lonp
          		latp = 90.-acos(z)*180./pi 
          	endif
          WRITE(1,*) lonp, latp, pixel+1
10      Continue 
        CLOSE(1) 
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
! === Step 2: The pixels are sorted for increasing latitude
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
! --- Reading "px.dat" to get info about lon-lat of pixels...
      open(1,file='px.dat',status='unknown') 
      do i=1, n
		read(1,*,end=2) lon(i), lat(i)  
      enddo
 2    Close(1) 
!
! --- Sorting pixels by increasing latitude 
!      Write(*,*) '    - Bubble-sorting pixels by increasing latitude...'
!      call BUBBLE_SORT_MOD(n, lat, lon)
      Write(*,*) '    - Sorting pixels by increasing latitude...'
      call QSORT_MOD(n, lat, lon)
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
! === Step 3: Finding "main latitudes" & pixels having these latitudes 
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!
 open(3,file='px-lat.dat',status='unknown')
 Write(3,*) "-------------------------------------------"  
 Write(3,*) " Main pixels latitudes and range of pixels "
 Write(3,*) "  pixel, lat (deg), min range, max range   "   
 Write(3,*) "-------------------------------------------" 
 nl=1 
 imin(nl)=1 
 imax(nl)=1 
 ianch(1)=1
 lateff(nl)=lat(1)
 do 20 i=2, np                
	if(abs(lat(i)-lat(i-1)).lt.eps) then 
	    imax(nl)=i
	else 
	    nl=nl+1 
	    lateff(nl)=lat(i)
	    imin(nl)=i 
	endif
        ianch(i)=nl
20 continue 
imax(nl)=np 
!
write(*,*) '    - There are', nl, ' main latitudes' 
write(*,*) '    - Number(pixels)/Number(main pixels)= ', & 
                  float(np)/float(nl)
!
 do j=1, nl 
     write(3,*) j, lateff(j), imin(j), imax(j)
 enddo
!
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
! === Step 4: Associating each pixel to its "anchor pixel". 
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
      write(*,*) '    - Building file pxa.dat' 
      open(8,file='pxa.dat',status='unknown') 
      Write(8,*) "----------------------------------------------"  
      Write(8,*) " Pixels lon-lats and with increasing latitude "
      Write(8,*) "  pixel lon(deg), lat (deg), anchor, pixel #  "   
      Write(8,*) "----------------------------------------------"         
      do 5 i=1, np 	
!      do 5 j=1, nl 	
!	   if(imin(j)<=i.and.i<=imax(j)) & 
!	   write(8,*) lon(i), lateff(j), j, i 	 
          write(8,*) lon(i), lat(i), ianch(i), i
5     continue
      close(8) 
!
!
!
    deallocate( IMIN, IMAX, IANCH )   
    deallocate( LON, LAT )
    deallocate( LATEFF )
!
!
! ---- Write the number of anchors to file anchor.tmp
!
    open(9,file='anchor.tmp',status='unknown')
    write(9,*) nl
    close(9)
!
!
!
END program pix
!
!
!
!
!
!




