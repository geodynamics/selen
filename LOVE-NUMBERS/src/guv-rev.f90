!
! /////////////////////////////////////////////////////////////////////////////
!
! This program computes the U, G and V rotational Green's functions for SELEN 4
! 
! The array dimensions are consistent with the <<REVISED rotational Theory>> 
!
! REVISED ROTATIONAL THEORY 
! REVISED ROTATIONAL THEORY 
! REVISED ROTATIONAL THEORY 
! REVISED ROTATIONAL THEORY 
! REVISED ROTATIONAL THEORY 
!
! /////////////////////////////////////////////////////////////////////////////
! Giorgio Spada October 2018. 
! Revised GS decembre 2019 
!
  MODULE DATAG
  IMPLICIT NONE 
!
! Parameters formerly in the input 'include' file - DM MAR 01, 2021
  INTEGER :: NV
  INTEGER :: NN
  REAL*4  :: DELTA_T
  CHARACTER*100 :: PATH_2_MODEL
  INTEGER :: MP
!
! INPUT file name with degree l=2 LNs, computed by "Love-Numbers-REV0.f90"
  Character*100, PARAMETER :: FILE_DEG2="l2-love-numbers.dat"
!
! INPUT file for the PMTF^ 
  Character*100, PARAMETER :: FILE_PMTF="PMTF-REV.dat"   
!
! OUTPUT file for the gamma functions......
  CHARACTER*72, PARAMETER :: ROT_FILE="gamma-REV.dat" 
!
! Various constants and variables  
  INTEGER :: M
!
  REAL*8, ALLOCATABLE :: A(:)
  REAL*8, ALLOCATABLE :: S(:)
!
  REAL*8 :: G_ELA_P, G_SEC_P 
  REAL*8, ALLOCATABLE :: G_RES_P(:)
  REAL*8, ALLOCATABLE :: G_RES_PP(:)
!
  REAL*8 :: U_ELA_P, U_SEC_P 
  REAL*8, ALLOCATABLE :: U_RES_P(:)
  REAL*8, ALLOCATABLE :: U_RES_PP(:)
!
! REV GS 011219 for "Horizontals"
!
  REAL*8 :: V_ELA_P, V_SEC_P 
  REAL*8, ALLOCATABLE :: V_RES_P(:)
  REAL*8, ALLOCATABLE :: V_RES_PP(:)
!
! REV DM 201025 for "GRACE geoid"
! 
  REAL*8 :: H_ELA_P, H_SEC_P
  REAL*8, ALLOCATABLE :: H_RES_P(:)
  REAL*8, ALLOCATABLE :: H_RES_PP(:)
!
  REAL*8, PARAMETER :: ERAD     =  6.371D6
  REAL*8, PARAMETER :: OMEGA    =  7.292115D-5
  REAL*8, PARAMETER :: SIGMA_21 =  OMEGA**2*ERAD**2/SQRT(30D0)
  REAL*8, PARAMETER :: GRA_REF  =  9.81D0 
!
  CHARACTER*20 DATE, TIMC
!
  END MODULE DATAG
!
!
!
!  Include "./Fortran-Utilities/RH.f90"
!
!
!
  PROGRAM GUROT 
  USE DATAG
  USE CONFIG_DATA
  IMPLICIT NONE 
!
  INTEGER I, J, K, N 
!
  REAL*8, ALLOCATABLE :: GAMMA_G(:)
  REAL*8, ALLOCATABLE :: GAMMA_U(:)
! REV GS 011219 for "Horizontals"
  REAL*8, ALLOCATABLE :: GAMMA_V(:)
! REV DM 201025 for "GRACE geoid"
  REAL*8, ALLOCATABLE :: GAMMA_H(:)
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
! Configuration
  call config
!
! --- Copy config parameters to local variables
!
  mp = cfg_mp
  nv = cfg_nv
  nn = cfg_nn
  delta_t      = cfg_delta_t
  path_2_model = trim(adjustl(cfg_model_file))
!
  M=4*NV
!
  ALLOCATE( A(M) )
  ALLOCATE( S(M) )
!
  ALLOCATE( G_RES_P(M)  )
  ALLOCATE( G_RES_PP(M) )
!
  ALLOCATE( U_RES_P(M)  )
  ALLOCATE( U_RES_PP(M) )
!
  ALLOCATE( V_RES_P(M)  )
  ALLOCATE( V_RES_PP(M) )  
! 
  ALLOCATE( H_RES_P(M)  )
  ALLOCATE( H_RES_PP(M) )
!
  ALLOCATE( GAMMA_G(0:NN+1) )
  ALLOCATE( GAMMA_U(0:NN+1) )
  ALLOCATE( GAMMA_V(0:NN+1) )
  ALLOCATE( GAMMA_H(0:NN+1) )
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
! Reading data from file FILE_DEG2 
  call READ_DATA 
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
! Computing Gamma and copying it on a file 
!
 OPEN(1,file=ROT_FILE,status="unknown") 
! 
 Call DATE_AND_TIME (date,timc)
 Write(1,*) "-----------------------------------------------------------------"
 Write(1,*) "File: ", trim(adjustl(ROT_FILE)), &
           " ", date(1:4), '.', date(5:6), '.', date(7:8), & 
           " ", timc(1:2), '.', timc(3:4), '.', timc(5:6)
 Write(1,"(A8,A100)")  " Model: ", PATH_2_MODEL 
 Write(1,"(A7,1x,I6)") " MP:    ", MP
 Write(1,*) "-----------------------------------------------------------------"	   
!
! 
 Write(*,*) "- Gamma is copied on file: ", trim(adjustl(ROT_FILE))
!  
  DO 1 n=0, nn+1
!
    GAMMA_G(n) = G_ELA_P + G_SEC_P*FLOAT(n)*DELTA_T 
    GAMMA_U(n) = U_ELA_P + U_SEC_P*FLOAT(n)*DELTA_T 
    GAMMA_V(n) = V_ELA_P + V_SEC_P*FLOAT(n)*DELTA_T 
    GAMMA_H(n) = H_ELA_P + H_SEC_P*FLOAT(n)*DELTA_T 
!
     DO i=1, M 
!
      GAMMA_G(n) = GAMMA_G(n) + & 
              (G_RES_P (i)/A(i))*(EXP(A(i)*FLOAT(n)*DELTA_T)-1D0)
!
      GAMMA_G(n) = GAMMA_G(n) + & 
              (G_RES_PP(i)/S(i))*(EXP(S(i)*FLOAT(n)*DELTA_T)-1D0)	  
!
      GAMMA_U(n) = GAMMA_U(n) + & 
              (U_RES_P (i)/A(i))*(EXP(A(i)*FLOAT(n)*DELTA_T)-1D0)
!
      GAMMA_U(n) = GAMMA_U(n) + & 
              (U_RES_PP(i)/S(i))*(EXP(S(i)*FLOAT(n)*DELTA_T)-1D0)	  
!
      GAMMA_V(n) = GAMMA_V(n) + & 
              (V_RES_P (i)/A(i))*(EXP(A(i)*FLOAT(n)*DELTA_T)-1D0)
!
      GAMMA_V(n) = GAMMA_V(n) + & 
              (V_RES_PP(i)/S(i))*(EXP(S(i)*FLOAT(n)*DELTA_T)-1D0)	  
!
      GAMMA_H(n) = GAMMA_H(n) + & 
              (H_RES_P (i)/A(i))*(EXP(A(i)*FLOAT(n)*DELTA_T)-1D0)
!
      GAMMA_H(n) = GAMMA_H(n) + & 
              (H_RES_PP(i)/S(i))*(EXP(S(i)*FLOAT(n)*DELTA_T)-1D0)	  
!
  Enddo
!
    GAMMA_G(n) = GAMMA_G(n) * SIGMA_21/GRA_REF
!
    GAMMA_U(n) = GAMMA_U(n) * SIGMA_21/GRA_REF
!
    GAMMA_V(n) = GAMMA_V(n) * SIGMA_21/GRA_REF
!
    GAMMA_H(n) = GAMMA_H(n) * SIGMA_21/GRA_REF
!
    write(1,*) n, GAMMA_G(n), GAMMA_U(n), GAMMA_V(n), GAMMA_H(n)
!
1 CONTINUE
!
  CLOSE(1)
!
! ********************
  END PROGRAM GUROT
! ********************
!
!
!
!
!
  SUBROUTINE READ_DATA 
  USE DATAG
  IMPLICIT NONE 
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Reads:
!  - B^ and B^^ terms from file "PMTF-REV-MP.dat" 
!    GS October 16, 2018 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  INTEGER, PARAMETER :: IU=1
  INTEGER I, J 
!
!  Write(*,*) "-Reading (in DP) data from: ", trim(adjustl(FILE_PMTF))
!
  OPEN(IU,FILE=FILE_PMTF,STATUS='UNKNOWN') 
!
  call rh(iu,5)
!
  call rh(iu,1)
  call rh(iu,3)
  call rh(iu,M)
!
  call rh(iu,1)
  call rh(iu,3)
  call rh(iu,M)
!  
  call rh(iu,1)
  call rh(iu,1)
  call rh(iu,M)
!
  call rh(iu,1)
  call rh(iu,1)
  READ(1,*)G_ELA_P 
  READ(1,*)G_SEC_P 
  DO I=1, M 
     READ(1,*) J, G_RES_P(I),  A(I)      
  ENDDO
!
  call rh(iu,1)
  call rh(iu,1)
  DO I=1, M 
     READ(1,*) J, G_RES_PP(I), S(I)      
  ENDDO
!
  call rh(iu,1)
  call rh(iu,1)
  READ(1,*)U_ELA_P 
  READ(1,*)U_SEC_P 
  DO I=1, M 
     READ(1,*) J, U_RES_P(I),  A(I)      
  ENDDO
!
  call rh(iu,1)
  call rh(iu,1)
  DO I=1, M 
     READ(1,*) J, U_RES_PP(I), S(I)      
  ENDDO
!
!--------------------------------- ! REV GS 251219 for "Horizontals"
!
  call rh(iu,1)
  call rh(iu,1)
  READ(1,*)V_ELA_P 
  READ(1,*)V_SEC_P 
  DO I=1, M
     READ(1,*) J, V_RES_P(I),  A(I)	 
  ENDDO

  call rh(iu,1)
  call rh(iu,1)
  DO I=1, M 
     READ(1,*) J, V_RES_PP(I), S(I)	 
  ENDDO
!
!--------------------------------- ! REV GS 251219 for "Horizontals"
!
!--------------------------------- ! REV GS 301025 for "GRACE Geoid"
!
  call rh(iu,1)
  call rh(iu,1)
  READ(1,*)H_ELA_P 
  READ(1,*)H_SEC_P 
  DO I=1, M
     READ(1,*) J, H_RES_P(I)
  ENDDO
!
  H_RES_PP = G_RES_PP
!
!--------------------------------- ! REV GS 301025 for "GRACE Geoid"
  CLOSE(IU) 
!
  END SUBROUTINE READ_DATA 
!
!
