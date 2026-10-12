!
!
!
!
!
!====================================
SUBROUTINE HQR_DRIVER (A, M, RTR, RTI)
USE TYPES
USE FMZM
IMPLICIT NONE  
!
! -----------------------------------------------
! THIS IS A DRIVER FOR THE ROOT FINDING PROCEDURE 
! GS JANUARY 7, 2008. =BASED ON NETLIB ROUTINES=
! -----------------------------------------------
!
!INTEGER, PARAMETER :: MAXM = 499   
INTEGER :: J, K, M, IGH, LOW, IERR  
!TYPE(FM) :: A (M + 1), RTR (M), RTI (M),  & 
!           H (MAXM, MAXM), ASCALE(MAXM), & 
!	   XR, XI   
TYPE(FM) :: A (M + 1), RTR (M), RTI (M)
TYPE(FM) :: XR, XI
TYPE(FM), ALLOCATABLE :: H (:, :), ASCALE(:) 
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! 
!If (m .gt. maxm .or. a(m+1) .eq. 0._qp) Stop
If (a(m+1) .eq. 0._qp) Stop
!
 allocate( h(m,m), ascale(m) )
!
! Building the Hessemberg matrix .....
 Do K = 1, M  
    H(1,K) = -A(M+1-K)/A(M+1)  
    DO J = 2, M  
        H(J,K) = 0._QP  
    enddo 
    IF (K.NE.M) H(K+1,K) = 1._QP  
 enddo  
!
! Computing the roots by the "Eigenvalues method" .....
!
 Call balanc_nlib (m, h, low, igh, ascale)  
 Call hqr_nlib (m, low, igh, h, rtr, rti, ierr)
!
! In normal conditions ierr=0
! Write(*,*) ierr
!
! Ascending sorting of real parts of  
! roots using a code by John Burkardt ..... 
!
 Call bubble_sort_qp(m, rtr)
!
 deallocate(h, ascale)
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! 
RETURN  
END
!
!
!

!
