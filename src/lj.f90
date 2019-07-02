!
!
!
   INTEGER FUNCTION LJ (J)
   IMPLICIT NONE 
   INTEGER J, L, M, K  
   REAL*4 Z
!
! Given the harmonic index j=l(l+1)/2+m+1, returns l (degree) 
! -GS 10.09.07 
!
 If(j<0) then 
	  Write(88,*) "LJ: The J-index is out of bounds" 
 	  Write(88,*) "The program will STOP ----------"
 	  Write(*, *) "LJ: The J-index is out of bounds"  
 	  Write(*, *) "The program will STOP ----------"	   
          call stop_config 
	  Stop
 endif
!
   m=0
   k=j 
10 z=(-1.0+sqrt(8.*float(k)-7.))/2.
   lj=aint(z) 
   if(z-aint(z)==0.) return 
   k=k-1
   m=m+1
   goto 10
   return 
   end function LJ 
!   
!
!  
