!
!
!
 FUNCTION SINDD(ALFA)
 implicit NONE
!
! --- Assuming that ALFA is given in DEGREES, it returns the SINE of ALFA
!     *********************** GS and FC 11-02-2009 **********************
!	
 REAL*8,  PARAMETER :: PI=3.14159265358979323840D0 
 REAL*8 SINDD, ALFA, ALFARAD		  
 ALFARAD=ALFA*PI/180D0	 
 SINDD=SIN(ALFARAD)	 
 END FUNCTION SINDD 
!
!
!
