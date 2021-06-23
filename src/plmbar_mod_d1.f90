!
!
!
!
!
!
!
  SUBROUTINE PLMBAR_MOD_D1 (LMAX, LAT, DPLM)
  USE SHTOOLS 
  IMPLICIT NONE 
  INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
  INTEGER J, J_INDEX, MJ, JMAX, LMAX
  REAL*8 LAT
  REAL*8 Z, PLM(JJMAX), DPLM(JJMAX)
  REAL*8 COSDD
!
! Given latitude LAT (in degrees), computes the detivatives wrt 'theta'
! of *all* the fully normalized associated Legendre functions \bar P_{lm}(x), 
! with x=cos(colatitude), to degree LMAX. Uses the SHTOOLS Legendre functions 
! by PlmBar_d1, rescales by sqrt(2-delta(0,m)), and transforms. 
!   
! Created by GS August 7 2008 for version 2.8 of SELEN  
!
!
! ---- Tests the latitude bounds 
!
  If(abs(lat).gt.90.) then 
 	        Write(88,*) "Error in Sbr. PLMBAR_MOD_D1: Latitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. PLMBAR_MOD_D1: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1
  elseif(lmax>llmax) then 
		Write(88,*) "Error in Sbr. PLMBAR_MOD_D1: The degree exceeds 256", lmax 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. PLMBAR_MOD_D1: The degree exceeds 256", lmax 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2		
  Endif  
!
!
! ---- "z" is cos(theta), theta is co-latitude 
!
	z=cosdd(90.-lat)   
!
! ---- Builds the SHTOOLS Legendre functions,
!      including the "Condon-Shortley phase"...  
!
  	call PlmBar_d1(plm, dplm, lmax, z, -1)  
!
!
! ---- Scales consistently to our convention for Legendre functions, 
!      dividing by factor sqrt(2-delta(0,m)), AND transofrms to obtain 
!      a derivative wrt 'theta' from the derivative wrt 'cos theta'
!  
  	do j=1, j_index(lmax,lmax)
     		if(mj(j)/=0) dplm(j)=dplm(j)/sqrt(2.)
	                     dplm(j)=-sqrt(1.-z**2)*dplm(j)		
  	enddo
!
  END SUBROUTINE PLMBAR_MOD_D1
!
!
!
!
!