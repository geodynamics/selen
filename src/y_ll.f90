!
!
!
!
 SUBROUTINE HARMO_Y_LL (LMAX, LON, LAT, Y_LL)
 USE SHTOOLS 
 IMPLICIT NONE 
 INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
 INTEGER J, MJ, LJ, J_INDEX, LMAX  
 COMPLEX*16 Y_LL(*) 
 REAL*8 LON, LAT
 REAL*8, ALLOCATABLE :: PLM(:)  
 REAL*8 COSDD, SINDD
 INTEGER :: JMAX
!
!
! Given longitude LON and latitude LAT  - in degrees - this routine computes 
! the theta-component of the gradient of *all* the 4-pi normalized complex 
! harmonics $\cal Y_{lm}(\theta,\lambda)$, with \theta=colatitude and \lambda=
! longitude, to degree LMAX. This subroutine uses the (modified) SHTOOLS Legendre 
! functions by PlmBar (see Sbr PLMBAR_MOD). - Created by GS August 8, 2008 - 
!
! CREATED by GS on MAY 1, 2020. The program computes y_tt
!
! - Modified by GS Apr 2020 -- dynamically allocate DPLM
!
!**********************************************************
!   Modified GS on MAy 1, 2020 
!
!   The routine computes "y_ll = (1/sin theta)d/dlambda (1/sin theta)d/dlambda y"
!**********************************************************
!
!
! ---- Tests the longitude and latitude bounds 
  If(abs(lat)>90.) then 
 	        Write(* ,*) "Error in Sbr. GRAD_THETA_HARMO: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1		
  elseif (lon<0.or.lon>360.) then
 	        Write(* ,*) "Error in Sbr. GRAD_THETA_HARMO: Longitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2	
  Endif   
!  
!
!
! ---- Allocate temporary arrays
!
   JMAX=(LMAX+1)*(LMAX+2)/2
   ALLOCATE( PLM (JMAX) )   
!
! ---- Builds the SHTOOLS Legendre functions, with Condon-Shortley phase 
!
   CALL PLMBAR_MOD (LMAX, LAT,  PLM)  
!
! 
! ---- Computes the 4-pi normalized Spherical Harmonics up to degree LMAX
!
    		do j=1, j_index(lmax,lmax)
!
       			y_ll(j) = plm(j)* &
			(-float(mj(j)**2)/sindd(90.-lat)/sindd(90.-lat))* &
			cmplx(cosdd(float(mj(j))*lon),sindd(float(mj(j))*lon))
!
    		enddo
!
!
! ---- Release memory
!
    DEALLOCATE( PLM )
!
END SUBROUTINE HARMO_Y_LL
!
!
!
!
!
