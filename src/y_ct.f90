!
!
!
!
 SUBROUTINE HARMO_Y_CT (LMAX, LON, LAT, Y_CT)
 USE SHTOOLS 
 IMPLICIT NONE 
 INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
 INTEGER J, MJ, J_INDEX, LMAX  
 COMPLEX*16, PARAMETER :: I = CMPLX(0.,1.)
 COMPLEX*16 Y_CT(*) 
 REAL*8 LON, LAT
 REAL*8, ALLOCATABLE :: DPLM(:) 
 REAL*8 COSDD, SINDD
 INTEGER :: JMAX
!
!
! Given longitude LON and latitude LAT  - in degrees - this routine computes 
! the lambda-component of the gradient of *all* the 4-pi normalized complex 
! harmonics $\cal Y_{lm}(\theta,\lambda)$, with \theta=colatitude and \lambda=
! longitude, to degree LMAX given. Uses the (modified) SHTOOLS Legendre 
! functions by PlmBar (see Sbr PLMBAR_MOD_D1). - Created by GS August 8, 2008 - 
!
! - Modified by GS Apr 2020 -- dynamically allocate PLM
!
!*********************************************************************
!   Modified GS on MAy 1, 2020 
!   Modified DM on Mar 11, 2021 - Fixed "90-lat" that was "lat-90"
!
!   The routine computes "y_ct = cot(theta)*(d/dtheta)y" 
!*********************************************************************
!
!
! ---- Tests the longitude and latitude bounds 
  If(abs(lat)>90.) then 
 	        Write(* ,*) "Error in Sbr. GRAD_LAMBDA_HARMO: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1		
  elseif (lon<0.or.lon>360.) then
 	        Write(* ,*) "Error in Sbr. GRAD_LAMBDA_HARMO: Longitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2	
  Endif   
!  
!
! ---- Allocate temporary arrays
!
   JMAX=(LMAX+1)*(LMAX+2)/2
   ALLOCATE( DPLM(JMAX) )
!
! ---- Builds the SHTOOLS Legendre functions, with Condon-Shortley phase 
!
   CALL PLMBAR_MOD_D1 (LMAX, LAT, DPLM)  
!
! 
! ---- Computes the 4-pi normalized Spherical Harmonics up to degree LMAX
!
    		do j=1, j_index(lmax,lmax)
!
       			y_ct(j) = &
	                         dplm(j)* &
	                         cmplx(cosdd(float(mj(j))*lon),sindd(float(mj(j))*lon))*&
				 cosdd(90-lat)/sindd(90-lat)
!
    		enddo
!
!
! ---- Release memory
!
    DEALLOCATE( DPLM )
!
END SUBROUTINE  HARMO_Y_CT
!
!
!
