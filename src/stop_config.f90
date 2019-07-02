!
!
!
        Subroutine STOP_CONFIG 
	implicit NONE
!
! --- Creates a STOP file that stops the execution of makeseles.sh 
!     ******************* GS and FC 13-12-2007 *******************
!
	open (11,file='stop.dat',status='new')
	Write(11,*) "Errors have been detected - SELEN will stop"
	close(11)
!
	write(* ,*) 'STOP_CONFIG: SELEN will STOP!' 
	write(88,*) 'STOP_CONFIG: SELEN will STOP!' 
	STOP 
!
	End subroutine STOP_CONFIG 	
!
!
!
