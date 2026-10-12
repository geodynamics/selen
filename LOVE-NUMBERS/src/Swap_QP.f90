!
!
!
!
!===================================
  SUBROUTINE SWAP_QP(X,Y)
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The routine Swaps two real values 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Adapted to quadruple precision by GS on January 2008
! Retouched by GS on December 2019 
! 
!  Modified:
!    01 May 2000
!  Author:
!    John Burkardt
!  Parameters:
!    Input/output, real X, Y.  
!    On output, the values of X and
!    Y have been interchanged.
!---------------------------------------- 
 USE TYPES
 USE FMZM
 IMPLICIT NONE
 SAVE
!
 TYPE(FM) X, Y, Z 
!
 z = x
 x = y
 y = z
!
 RETURN
 END
!---------------------------------------- 
!
!
!
!
!
!
