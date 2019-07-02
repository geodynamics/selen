!
!
FUNCTION J_INDEX(L,M) 
IMPLICIT NONE
INTEGER J_INDEX, L, M 
!
! given l (degree) and m(order), this function computes 
! the j-index "J", with J=l(l+1)/2+m+1 - GS 14-09-2007
!
 If(l<0.or.m<0.or.m>l) then 
	  Write(88,*) "J_INDEX: The degree and order are out of bounds" 
 	  Write(88,*) "The program will STOP -------------------------"
 	  Write(*, *) "J_INDEX: The degree and order are out of bounds"  
 	  Write(*, *) "The program will STOP -------------------------"	   
          call stop_config 
	  Stop
 endif
!
  J_INDEX = l*(l+1)/2+m+1
!
END FUNCTION J_INDEX 
!
!
!
