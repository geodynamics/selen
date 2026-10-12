!
!
!====================================
  SUBROUTINE READ_DEG2 
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Reads degree 2 data from an external file 
!  produced by TABOO.f90 
!  GS October 3, 2018 
!  GS DEcember 16, 2019 in mortis A.R. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  USE DATA
  USE FMZM
  IMPLICIT NONE
!
  Integer, parameter :: iu=1
  Integer I 
!
  OPEN(iu,FILE=FILE_DEG2,STATUS='UNKNOWN') 
!
! >>> Reading the file header 
!
  Call rh(iu,6)
!
! >>> Reading k_2_e_load 
!
  Call rh(iu,2)
  Call FMREAD(iu,KL_E)
!
! >>> Reading k_2_f_load 
!
  Call rh(iu,2)
  Call FMREAD(iu,KL_F)
!
! >>> Reading k_2_i_load 
!
  Call rh(iu,2)
  DO I=1, M 
     Call FMREAD(iu,KL(I))
  ENDDO
!
! >>> Secular coefficients 
!
  Call rh(iu,2)
  DO I=1, M+1 
     Call FMREAD(iu,COEFF(I))
  ENDDO
  DO I=0, M 
     P(I)=COEFF(M+1-I) 
  ENDDO
!
! >>> Roots of the secular equation
!
  Call rh(iu,2)       
  DO I=1, M 
     Call FMREAD(iu,S(I))
  ENDDO   
!
! >>> Reading h_2_e_tide
!      
  Call rh(iu,2)       
  Call FMREAD(iu,HT_E)
!
! >>> Reading h_2_f_tide  
!   
  Call rh(iu,2)       
  Call FMREAD(iu,HT_F)
!
! >>> Reading h_2_i_tide 
!
  Call rh(iu,2)
  DO I=1, M 
     Call FMREAD(iu,HT(I))
  ENDDO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GS 301119
!
! >>> Reading l_2_e_tide
!      
  Call rh(iu,2)       
   Call FMREAD(iu,LT_E)
!
! >>> Reading l_2_f_tide
!     
  Call rh(iu,2)       
   Call FMREAD(iu,LT_F)
!
! >>> Reading l_2_i_tide 
!
  Call rh(iu,2)
  DO I=1, M 
     Call FMREAD(iu,LT(I))
  ENDDO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GS 301119
!
! >>> Reading k_2_e_tide
!     
  Call rh(iu,2)       
  Call FMREAD(iu,KT_E)
!
! >>> Reading k_2_f_tide
!     
  Call rh(iu,2)       
  Call FMREAD(iu,KT_F)
!
! >>> Reading k_2_i_tide
!
  Call rh(iu,2)
  DO I=1, M 
     Call FMREAD(iu,KT(I))
  ENDDO
!
 CLOSE(iu)  
!
  END SUBROUTINE READ_DEG2 
!
!
!
!
