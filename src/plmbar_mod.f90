!
!
!
  SUBROUTINE PLMBAR_MOD (LMAX, LAT, PLM)
  USE SHTOOLS 
  IMPLICIT NONE 
  INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
  INTEGER J, J_INDEX, MJ, JMAX, LMAX
  REAL*8 LAT
  REAL*8 Z, PLM(JJMAX)
  REAL*8 COSDD
!
! Given latitude LAT (in degrees), computes *all* the fully normalized 
! associated Legendre functions \bar P_{lm}(x), with x=cos(colatitude), 
! to degree LMAX. Uses the SHTOOLS Legendre functions by PlmBar, and 
! rescales by sqrt(2-delta(0,m)).   Last modified by GS 10-9-2007 - 
!
!
! ---- Tests the latitude bounds 
!
  If(abs(lat).gt.90.) then 
 	        Write(88,*) "Error in Sbr. PLMBAR_MOD: Latitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. PLMBAR_MOD: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1
  elseif(lmax>llmax) then 
		Write(88,*) "Error in Sbr. PLMBAR_MOD: The degree exceeds 256", lmax 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. PLMBAR_MOD: The degree exceeds 256", lmax 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2		
  Endif  
!
!
!
! ---- "z" is cos(theta), theta is co-latitude 
!
	z=cosdd(90.-lat) 
!
! ---- Builds the SHTOOLS Legendre functions including the Condon-Shortley phase 
!
  	call PlmBar(plm, lmax, z, -1)  
!
!
! ---- Scales to obtain our Legendre functions, dividing by sqrt(2-delta(0,m))
!  
  	do j=1, j_index(lmax,lmax)
     		if(mj(j)/=0) plm(j)=plm(j)/sqrt(2.)
  	enddo
!
  END SUBROUTINE PLMBAR_MOD 
!
!
!
