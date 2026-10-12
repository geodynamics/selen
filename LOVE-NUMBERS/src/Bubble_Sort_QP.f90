
!=====================================================================
  SUBROUTINE BUBBLE_SORT_QP(N, A)
!
! The routine sorts a real vector using !bubble sorting" 
!  ----------------------------------------------------
!  Adapted to quadruple precision by GS on January 2008
!  ----------------------------------------------------
!  Discussion:
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!  Modified:
!    01 February 2001
!  Author:
!    John Burkardt
!  Parameters:
!    Input, integer N, the number of entries in the array.
!    Input/output, real A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!---------------------------------------- 
!
 USE TYPES
 USE FMZM
 IMPLICIT NONE
 SAVE
!
 INTEGER I, J, N
 TYPE(FM) A(N)
!
 do i = 1, n-1
   do j = i+1, n
     if ( a(i) > a(j) ) then
       call swap_qp (a(i), a(j))
     end if
   end do
 end do
!
 return
 end
!---------------------------------------- 
!
!
!
!
!
!
!
