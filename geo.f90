!
! This is program "GEO.F90" 
!
! Created by GS September 2008 for version 2.7 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Revised GS July 2010 - g95 - Double precision implementation
! Revised GS December 2010- Longitude and latitude are now default 
! === Revised GS & FC May 21 2011 - DELTA parameter in ice history     

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
!  This program first computes the 4pi-normalized, complex spherical harmonics 
!  at geodetic sites listed in file GEODETIC_DATABASE. Subsequantly, it computes 
!  the 3D displacement (UP, NORTH, EAST), and the scalars S and N at the same 
!  sites. This set of five quantites specifies completely the "geodetic state" 
!  of the site. 
!
!  Input files:
!	- A list of geodetic sites (e. g., "gps-sites.dat")
!	- 'shs.bin' 
!	- 'shu.bin' 
!       - 'shn.bin'
!       - 'shv.bin'
!  
!  Output files: 
!	- 'geodetic-predictions.dat'	
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! INCLUDE "harmonics.f90"
 PROGRAM GEO
 IMPLICIT NONE 
 INCLUDE "data.inc"
 CHARACTER*11, PARAMETER :: UF='UNFORMATTED', F='FORMATTED', UN='UNKNOWN' 
 CHARACTER*22 CJUNK
 CHARACTER*22, ALLOCATABLE :: NAME(:) 
 CHARACTER*40 LITHOC, VVVV(NV)
 CHARACTER*20 DATE, TIMC            
 INTEGER I, J, DOM
 INTEGER, ALLOCATABLE :: CODE(:)  
 REAL*8 AJUNK
 REAL*8, ALLOCATABLE :: LONS(:), LATS(:)
 COMPLEX*16, ALLOCATABLE :: CS(:,:), CU(:,:), CN(:,:), CV(:,:)
 COMPLEX*16, ALLOCATABLE :: Y(:,:), GRD_THETA_Y(:,:), GRD_LAMBDA_Y(:,:)
 REAL*8, ALLOCATABLE :: RATE_S(:), RATE_N(:), RATE_UP(:), RATE_NO(:), RATE_EA(:)        
! 
 REAL*8, PARAMETER :: DDELTA=DELTA 
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
! -- Allocate dynamic arrays
  allocate( name(ngeod) )
  allocate( code(ngeod) )
  allocate( lons(ngeod), lats(ngeod) )
  allocate( CS(JMAX,0:NN), CU(JMAX,0:NN), CN(JMAX,0:NN), CV(JMAX,0:NN) )
  allocate( Y(JMAX,NGEOD), GRD_THETA_Y(JMAX,NGEOD), GRD_LAMBDA_Y(JMAX,NGEOD) )
  allocate( RATE_S(NGEOD), RATE_N(NGEOD), RATE_UP(NGEOD), RATE_NO(NGEOD), RATE_EA(NGEOD) )
!
  write(*,*) '    - Reading geodetic sites coordinates'
  OPEN(1,FILE=GEODETIC_DATABASE,STATUS='unknown')
  Do i=1, ngeod
       read(1,*)          code(i) 
       read(1,'(a22)')    name(i) 
!
! The first datum is LONGITUDE, the second is LATITUDE   
!
       read(1,*) lons(i), lats(i)  
       if(lons(i).le.0.) lons(i)=360.+lons(i)  
!		
       read(1,'(a22)')cjunk 	
!
  Enddo
  Close(1)
!
!
  write(*,*) '    - Computing the harmonics at the geodetic sites'
!
 call harmo            (lmax, lons(1), lats(1), y(:,1)) 
 call grad_theta_harmo (lmax, lons(1), lats(1), grd_theta_y(:,1))
 call grad_lambda_harmo(lmax, lons(1), lats(1), grd_lambda_y(:,1))
!
!!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(LONS,LATS,Y,GRD_THETA_Y,GRD_LAMBDA_Y) SCHEDULE(GUIDED)  
  do i=2, NGEOD    
     If(mod(i,15)==0) write(*,*) '    - geo.f:', i, ' sites of', ngeod
       call harmo            (lmax, lons(i), lats(i), y(:,i)) 
       call grad_theta_harmo (lmax, lons(i), lats(i), grd_theta_y(:,i))
       call grad_lambda_harmo(lmax, lons(i), lats(i), grd_lambda_y(:,i))
  enddo
!!$OMP END PARALLEL DO
!

! 
! Write(*,*) '    - geo.f: Reading the harmonics from binary files'
 open (101,file='shs.bin',status='UNKNOWN',form='UNFORMATTED') 
 open (102,file='shu.bin',status='UNKNOWN',form='UNFORMATTED') 
 open (103,file='shn.bin',status='UNKNOWN',form='UNFORMATTED') 
 open (104,file='shv.bin',status='UNKNOWN',form='UNFORMATTED') 
 read (101) CS
 read (102) CU
 read (103) CN
 read (104) CV
 close(101) ; close(102) ; close(103) ; close(104) 
!
 do 1 i=1, ngeod 
!
       rate_s(i) =0.  ! sea level 
       rate_n(i) =0.  ! geoid undulation 
       rate_up(i)=0.  ! velocity 'up' 
       rate_no(i)=0.  ! velocity 'north' 
       rate_ea(i)=0.  ! velocity 'east'  
! 
       do 2 j=1, jmax 
            AJUNK= (2.-dom(j))      
   	    rate_s (i) = rate_s (i) + AJUNK*real(((cs(j,nn)-cs(j,nn-1))/DDELTA)*y(j,i))  
   	    rate_n (i) = rate_n (i) + AJUNK*real(((cn(j,nn)-cn(j,nn-1))/DDELTA)*y(j,i))  
   	    rate_up(i) = rate_up(i) + AJUNK*real(((cu(j,nn)-cu(j,nn-1))/DDELTA)*y(j,i))  
   	    rate_no(i) = rate_no(i) + AJUNK*real(((cv(j,nn)-cv(j,nn-1))/DDELTA)*grd_theta_y(j,i))  
   	    rate_ea(i) = rate_ea(i) + AJUNK*real(((cv(j,nn)-cv(j,nn-1))/DDELTA)*grd_lambda_y(j,i)) 
 2    continue 
!
      		rate_no(i) = -rate_no(i)
!
 1    continue 
!
!
!
!------- File for scattered geodetic predictions 
 open(82,file='geodetic-predictions.dat',status='unknown')
!
! --- Header information about the SELEN settings 
!
 open(76, file=VISCO_MODEL,status='unknown') 
 read(76, '(a40)') lithoc 
 	do i=1, nv 
          read(76, '(a40)') vvvv(i)  	 
	enddo
 close(76) 
 call DATE_AND_TIME (date,timc) 
WRITE(82,*)'# SELEN 2.7'     
WRITE(82,*)'# ', date(1:4), '.', date(5:6), '.', date(7:8), & 
           ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)
WRITE(82,*)"# Ice model: ", trim(adjustl(ice_model))
WRITE(82,*)"# Geodetic database:   ", GEODETIC_DATABASE 
WRITE(82,*)"# Number of mantle layers: ", NV 
WRITE(82,*)"# Model code (see TABOO User guide): ", CDE
WRITE(82,*)"# Viscosity model: ", visco_model
WRITE(82,*)"# LITHO thickness (km): ", LITHOC 
WRITE(82,*)"# Mantle viscosity from bottom to top (* 1E21 Pa.s): "
DO I=1, NV 
   WRITE(82,*)"# ", vvvv(i) 
ENDDO
WRITE(82,*)"# SLE mode of solution: ", IMODE
WRITE(82,*)"# SLE iterations: ", SMAX 
WRITE(82,*)"# Maximum harmonic degree: ", LMAX 
WRITE(82,*)"# Tegmark resolution: ", RES 
WRITE(82,*)"#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
WRITE(82,*)"# 	             /// Geodetic predictions for all the sites in database ///                        "
WRITE(82,*)"#                                                                                                  "
WRITE(82,*)"#  Code   lon       lat      ------ all velocities are in units of mm/yr ------                    "
WRITE(82,*)"#         deg       deg       UP         NORTH      EAST       S-dot      N-dot    Site name       "  			     
WRITE(82,*)"#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
WRITE(82,*)"#	 "
do i=1, ngeod 	
    Write(82,'(3x,i4,2(1x,f8.2),5(f10.2,1x),4x,a22)') & 
		  code(i), lons(i), lats(i), rate_up(i), rate_no(i), rate_ea(i), rate_s(i), rate_n(i), name(i) 
enddo
 close(82)
!
! -- Dellocate dynamic arrays
  deallocate( name )
  deallocate( code )
  deallocate( lons, lats )
  deallocate( CS, CU, CN, CV )
  deallocate( Y, GRD_THETA_Y, GRD_LAMBDA_Y )
  deallocate( RATE_S, RATE_N, RATE_UP, RATE_NO, RATE_EA )

!
!
 END PROGRAM GEO
!
!
!
!
!
!
