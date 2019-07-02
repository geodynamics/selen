!
!
!
   INTEGER FUNCTION DOM (J)
   IMPLICIT NONE 
   INTEGER J, MJ
!
! Given the index j=l(l+1)/2+m+1, returns 
! delta(0,m), with m=order - GS 14.09.07 
!
 If(j<0) then 
!Write(88,*) "DOM: The J-index is out of bounds" 
!Write(88,*) "The program will STOP -----------"
 	  Write(*, *) "DOM: The J-index is out of bounds"  
 	  Write(*, *) "The program will STOP -----------"	   
!call stop_config 
	  Stop
 endif
!
   DOM=0 
   If(mj(j)==0) DOM=1
!
   END FUNCTION DOM
!
!
!
