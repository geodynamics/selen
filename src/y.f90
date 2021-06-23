!
!
!
 SUBROUTINE HARMO_Y (LMAX, LON, LAT, Y)
 USE SHTOOLS 
 IMPLICIT NONE 
!INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
 INTEGER J, MJ, J_INDEX, LMAX  
 COMPLEX*16 Y(*) 
 REAL*8 LON, LAT
 REAL*8, ALLOCATABLE :: PLM(:) 
 REAL*8 COSDD, SINDD
 INTEGER :: JMAX
!
!
! Given longitude LON and latitude LAT  - in degrees - this routine computes 
! *all* the 4-pi normalized complex harmonics $\cal Y_{lm}(\theta,\lambda)$, 
! with \theta=colatitude and \lambda=longitude, to degree LMAX given. Uses  
! the (modified) SHTOOLS Legendre functions by PlmBar (see Sbr PLMBAR_MOD). 
! - Last modified by GS 12-14-2007 - 
! - Modified by DM Apr 2019 -- dynamically allocate PLM
!
!**********************************************************
!   Modified GS on MAy 1, 2020 
!
!   The routine computes "y"
!**********************************************************
!
! ---- Tests the longitude and latitude bounds 
  If(abs(lat)>90.) then 
 	        Write(* ,*) "Error in Sbr. HARMO: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1		
  elseif (lon<0.or.lon>360.) then
 	        Write(* ,*) "Error in Sbr. HARMO: Longitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2	
  Endif   
!
!
! ---- Allocate temporary arrays
!
   JMAX=(LMAX+1)*(LMAX+2)/2
   ALLOCATE( PLM(JMAX) )
!  
!
! ---- Builds the SHTOOLS Legendre functions, with Condon-Shortley phase 
!
    CALL PLMBAR_MOD(LMAX, LAT, PLM)  
!
! 
! ---- Computes the 4-pi normalized Spherical Harmonics up to degree LMAX
!
    do j=1, j_index(lmax,lmax)
!
       	y(j) = plm(j)*cmplx(cosdd(float(mj(j))*lon),sindd(float(mj(j))*lon)) 
!
    enddo
!
!
! ---- Release memory
!
    DEALLOCATE( PLM )
!
END SUBROUTINE HARMO_Y
!
!
!
