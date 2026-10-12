!
! //////////////////////////////////////////////////////////////////
! This program computes the Polar Motion Excitation Function (PMTF) 
! and the PRIMED PMTF^ once the degree 2 Love numbers are given. The
! multi-precision routines "FM" are employed to compute efficiently
! the roots and the residues... 
!
! There are M-1 rotational modes A, and the secular residue is /= zero. 
!
! The program also computes the "B" coefficients in the expression 
! of the rotational Response function... 
! 
! The TRADITIONAL rotational theory is used: the fluid "k" Love number is:
!
! kf = ktf[LT,s=0]  
! //////////////////////////////////////////////////////////////////
! ******************************************************************
! Giorgio Spada October 2018. 
! REVISED December 2019 
!
  MODULE DATA
  USE FMZM 
  IMPLICIT NONE 
!
  INTEGER :: NV
  INTEGER :: MP
  CHARACTER*100 :: PATH_2_MODEL
! INCLUDE "TABOO.inc"
!
! INPUT file with degree l=2 LNs, obtained by "TABOO.f90"
  Character*100, PARAMETER :: FILE_DEG2="l2-data.dat"   
!
! OUTPUT file for the PMTF, PMTF' and other data  
  Character*100, PARAMETER :: FILE_PMTF="PMTF.dat"   
!
! Multi-precition (MP) digits
  INTEGER :: IP
! 
! Various constants and variables  
  INTEGER :: M
  TYPE(FM) :: KL_E, KL_F
  TYPE(FM) :: KT_E, KT_F  
  TYPE(FM) :: HT_E, HT_F
  TYPE(FM), ALLOCATABLE :: KL(:), KT(:), HT(:)
! GS 301119
! GS 161219
  TYPE(FM) :: LT_E, LT_F
  TYPE(FM), ALLOCATABLE :: LT(:)
  TYPE(FM), ALLOCATABLE :: S(:)
  TYPE(FM), ALLOCATABLE :: COEFF(:)
  TYPE(FM), ALLOCATABLE :: P(:)
  TYPE(FM), ALLOCATABLE :: H(:)      
  TYPE(FM), ALLOCATABLE :: B(:,:)    
  TYPE(FM), ALLOCATABLE :: C(:)
  TYPE(FM), ALLOCATABLE :: NUME(:)
  TYPE(FM), ALLOCATABLE :: DENO(:) 
  TYPE(FM), ALLOCATABLE :: COEF_REAL(:)
  TYPE(FM), ALLOCATABLE :: ROOT_REAL(:)
  TYPE(FM), ALLOCATABLE :: ROOT_IMAG(:)
  TYPE(FM), ALLOCATABLE :: NUME_P(:) 
  TYPE(FM), ALLOCATABLE :: D_DENO(:)
  TYPE(FM) :: NNN, DDD
  CHARACTER*20 DATE, TIMC
!
  TYPE(FM), ALLOCATABLE :: F(:)
  TYPE(FM), ALLOCATABLE :: G(:)
!
! PMTF components
  TYPE(FM) :: A_ELA, A_SEC 
  TYPE(FM), ALLOCATABLE :: A_RES(:)
  TYPE(FM), ALLOCATABLE :: A_ROT(:)
!
! PMTF^ components
  TYPE(FM) :: A_ELA_P, A_SEC_P 
  TYPE(FM), ALLOCATABLE :: A_RES_P(:)
  TYPE(FM), ALLOCATABLE :: A_RES_PP(:)
! 
! B^ coefficients for the G^rot Green's function 
  TYPE(FM) :: G_ELA_P, G_SEC_P 
  TYPE(FM), ALLOCATABLE :: G_RES_P(:)
  TYPE(FM), ALLOCATABLE :: G_RES_PP(:)
!
! B^ coefficients for the U^rot Green's function 
  TYPE(FM) :: U_ELA_P, U_SEC_P 
  TYPE(FM), ALLOCATABLE :: U_RES_P(:)
  TYPE(FM), ALLOCATABLE :: U_RES_PP(:)
!
! B^ coefficients for the V^rot Green's function 
! GS 301119
  TYPE(FM) :: V_ELA_P, V_SEC_P 
  TYPE(FM), ALLOCATABLE :: V_RES_P(:)
  TYPE(FM), ALLOCATABLE :: V_RES_PP(:)
!
! B^ coefficients for the G^rot Green function (w/o direct centrifugal potential)
  TYPE(FM) :: H_ELA_P, H_SEC_P
  TYPE(FM), ALLOCATABLE :: H_RES_P(:)
!
  END MODULE DATA
!
! ******************************************************************
!
!
!
!Include "./Fortran-Utilities/Balanc_Nlib.f90"
!Include "./Fortran-Utilities/Bubble_Sort_QP.f90"
!Include "./Fortran-Utilities/HQR_Nlib.f90"
!Include "./Fortran-Utilities/HQR_Driver.f90"
!Include "./Fortran-Utilities/Swap_QP.f90"
!Include "./Fortran-Utilities/RH.f90"
!Include "./Fortran-Utilities/Read_DEG2.f90"
!
!
! \\\\\\\\\\\\\\
  Program PMTF
! //////////////
!  
  USE FMZM
  USE DATA
  USE CONFIG_DATA
  IMPLICIT NONE
!
  INTEGER I, J, K, JP, I1, I2 
!
  character(200) :: buffer   
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! ///////////////////////////////////////////////////
!
!!!!!!!!!!!!!!!!!!!!!! CONFIGURATION
!
 call config 
!
! --- Copy config parameters to local variables
!
 mp = cfg_mp
 nv = cfg_nv
 path_2_model = trim(adjustl(cfg_model_file))
!
 IP = MP
!
 M = 4*NV
!
 ALLOCATE( KL(M), KT(M), HT(M) ) 
 ALLOCATE( LT(M) )
 ALLOCATE( S(M) )
 ALLOCATE( COEFF(M+1) )
 ALLOCATE( P(0:M) )
 ALLOCATE( H(M) )  
 ALLOCATE( B(0:M-1,M) )    
 ALLOCATE( C(0:M-1) )
 ALLOCATE( NUME(0:M) )
 ALLOCATE( DENO(0:M-1) ) 
 ALLOCATE( COEF_REAL(M) )
 ALLOCATE( ROOT_REAL(1:M-1) )
 ALLOCATE( ROOT_IMAG(1:M-1) )
 ALLOCATE( NUME_P(0:M-2) )
 ALLOCATE( D_DENO(0:M-2) )
!
 ALLOCATE( F(1:M-1) )
 ALLOCATE( G(1:M) )
!
 ALLOCATE( A_RES(M-1) )
 ALLOCATE( A_ROT(M-1) )
!
 ALLOCATE( A_RES_P(M-1) )
 ALLOCATE( A_RES_PP(M) )
!
 ALLOCATE( G_RES_P(M-1) )
 ALLOCATE( G_RES_PP(M) )
!
 ALLOCATE( U_RES_P(M-1) )
 ALLOCATE( U_RES_PP(M) )
!
 ALLOCATE( V_RES_P(M-1) )
 ALLOCATE( V_RES_PP(M) )
!
 ALLOCATE( H_RES_P(M-1) )
!
!
 !!!!!!!!!!!!!!!!!!!!! END OF CONFIGURATION
!
!
!!!!!!!!!!!!!!!!!!!!!!
  Call FM_SET(IP) 
!!!!!!!!!!!!!!!!!!!!!!
!
! Reading data from file FILE_DEG2 
  Call READ_DEG2 
!
! --------------------------------------------
!  Write(*,*) "= PMTF"
! --------------------------------------------!
!
!  Definition of the (real) h_i constants (i=1,...,M)
!  Write (*,*)"- Computing the 'h', 'b' and 'c' constants"  
  DO I=1, M 
     H(I) = -KT(I)/KT_F/S(I) 
  ENDDO
!
!
! Division of polynomial "p" by "s-s_i"; array 'b' is real 
! The other possible recurrence is unstable ...    
  DO I=1, M 
     B(M-1,I)=P(M) 
        DO K=2, M  
           B(M-K,I)=P(M-K+1)+B(M-K+1,I)*S(I)
       ENDDO
  ENDDO
!
!
! Definition of the (real) c_k constants (k=0,...,M-1)
  DO K=0, M-1
     C(K)=TO_FM(0) 
 	DO I=1, M 
 	   C(K)=C(K)+H(I)*B(K,I)
	ENDDO
  ENDDO
!
!
! We now turn the fraction upside down and we re-define the 
! polynomials on top and on bottom. Both have real coefficients.  
  NUME(:)=P(:)
  DENO(:)=C(:)
!
!  
! Prepares the (real) coefficients for the root finder   
  DO I=1, M
     COEF_REAL(I)= DENO(M-I)
  ENDDO 
!
! Elastic residue  
  A_ELA = NUME(M)/DENO(M-1) 
!
! Secular residue 
  A_SEC = NUME(0)/DENO(0) 
! 
  Write(*,*) "= Some particular values: "
  Write(*,*) "- Elastic and secular A  ", to_dp(A_ELA), to_dp(A_SEC) 
!
! Multi-precision root finder... (the same as in TABOO_HP) 
  CALL HQR_DRIVER (COEF_REAL, M-1, ROOT_REAL, ROOT_IMAG)  
!
!
! Rotational roots 
  DO I=1, M-1
     A_ROT(I)=TO_FM(1)/ROOT_REAL(I)
!    Write(*,*) i, to_dp(ROOT_REAL(i)), to_dp(ROOT_IMAG(i)) 
  ENDDO
!
!
! Defining a new polynomial at the numerator... 
  DO K=0, M-2
     NUME_P(K)=NUME(K+1) - A_ELA*DENO(K) - A_SEC*DENO(K+1) 
  ENDDO 
!
!
! Denominator is s*'deno'. The derivative is thus a degree M-2 polynomial
  DO K=0, M-2
     D_DENO(K)=DENO(K+1)*FLOAT(K+1) 
  ENDDO
!
!
! The rotational modes are computed using the traditional formula
  DO K=1, M-1  
     NNN=TO_FM(0) 
     DDD=TO_FM(0) 
       DO J=0, M-2 
           NNN = NNN + NUME_P(J)*A_ROT(K)**J 
           DDD = DDD + D_DENO(J)*A_ROT(K)**J 
       ENDDO	
  A_RES(K)=NNN/DDD
  ENDDO
!
! --------------------------------------------
!  Write(*,*) "= PMTF^"
! --------------------------------------------
!
! Primed elastic residue  
  A_ELA_P = A_ELA*(TO_FM(1)+KL_E)
!
! Primed seular residue  
  A_SEC_P = A_SEC*(TO_FM(1)+KL_F)
!
  Write(*,*) "- Elastic and secular A^ ", to_dp(A_ELA_P), to_dp(A_SEC_P) 
!
!
! Here, I follow the formulation of my revised theory of summer 2018... 
! Array F
  DO J=1, M-1 
     F(J)=TO_FM(0)
    DO JP=1, M 
      F(J)=F(J)+A_RES(J)*KL(JP)/(A_ROT(J)-S(JP))    
    ENDDO   
  ENDDO
!
!
! Here, I follow the formulation of my revised theory of summer 2018... 
! Array G
  DO J=1, M 
     G(J)=TO_FM(0)
    DO JP=1, M-1 
      G(J)=G(J)+A_RES(JP)*KL(J)/(A_ROT(JP)-S(J))    
    ENDDO   
  ENDDO
!
! Computing the primed rotational residues  
!  Write(*,*) "- A_i^ and A_i^^..."
  DO J=1, M-1 
     A_RES_P(J) = A_RES(J)*(TO_FM(1) + KL_E) + F(J)
  ENDDO
!
! Computing the doubly primed rotational residues  
!  Write(*,*) "- A_i^^"
  DO J=1, M  
     A_RES_PP(J) = KL(J)*(A_ELA + A_SEC/S(J)) - G(J)
  ENDDO
!
!
! ---------------------------------------
!  Write(*,*) "= G^rot"
! ---------------------------------------
!
! Elastic 
  G_ELA_P = (TO_FM(1) + KT_E) * A_ELA_P
  H_ELA_P = (TO_FM(0) + KT_E) * A_ELA_P
!  Write(*,*) "- Elastic G^ ", to_sp(G_ELA_P)
!
! Secular 
  G_SEC_P = (TO_FM(1) + KT_F) * A_SEC_P
  H_SEC_P = (TO_FM(0) + KT_F) * A_SEC_P
! Write(*,*) "- Secular G^ ", to_sp(G_SEC_P)
!
  Write(*,*) "- Elastic and secular G^ ", to_dp(G_ELA_P), to_dp(G_SEC_P)
!
! Array F
! Write(*,*) "- Auxiliary Array F"
  F(:)=TO_FM(0)
  DO J=1, M-1 
     F(J)=TO_FM(0)
    DO JP=1, M 
      F(J)=F(J)+A_RES_P(J)*KT(JP)/(A_ROT(J)-S(JP))    
    ENDDO 
  ENDDO
!
! Array G
! Write(*,*) "- Auxiliary Array G"
  G(:)=TO_FM(0)
  DO J=1, M 
     G(J)=TO_FM(0)
    DO JP=1, M-1 
      G(J)=G(J)+A_RES_P(JP)*KT(J)/(A_ROT(JP)-S(J))    
    ENDDO 
  ENDDO
!
! Array B^ for G^rot 
  G_RES_P(:)=TO_FM(0)
  DO J=1, M-1 
     G_RES_P(J) = (TO_FM(1) + KT_E)*A_RES_P(J) + F(J)
  ENDDO
!
! Array B^ for G^rot w/o direct centrifugal term
  H_RES_P(:)=TO_FM(0)
  DO J=1, M-1
     H_RES_P(J) = (TO_FM(0) + KT_E)*A_RES_P(J) + F(J)
  ENDDO
!
! Array B^^ for G^rot 
  G_RES_PP(:)=TO_FM(0)
  DO J=1, M 
     G_RES_PP(J) = KT(J) * (A_ELA_P + A_SEC_P/S(J)) - G(J)
  ENDDO
!
!
! ---------------------------------------
!  Write(*,*) "= U^rot"
! ---------------------------------------
!
! Elastic 
  U_ELA_P = (TO_FM(0) + HT_E) * A_ELA_P
!
! Secular 
  U_SEC_P = (TO_FM(0) + HT_F) * A_SEC_P
!
  Write(*,*) "- Elastic and secular U^ ", to_dp(U_ELA_P), to_dp(U_SEC_P)
!
! Array F
!  Write(*,*) "- Auxiliary Array F"
  F(:)=TO_FM(0)
  DO J=1, M-1 
     F(J)=TO_FM(0)
    DO JP=1, M 
      F(J)=F(J)+A_RES_P(J)*HT(JP)/(A_ROT(J)-S(JP))    
    ENDDO 
  ENDDO
!
! Array G
!  Write(*,*) "- Auxiliary Array G"
  G(:)=TO_FM(0)
  DO J=1, M 
     G(J)=TO_FM(0)
    DO JP=1, M-1 
      G(J)=G(J)+A_RES_P(JP)*HT(J)/(A_ROT(JP)-S(J))    
    ENDDO 
  ENDDO
!
! Array B^ for U^rot 
  U_RES_P(:)=TO_FM(0)
  DO J=1, M-1 
     U_RES_P(J) = (TO_FM(0) + HT_E)*A_RES_P(J) + F(J)
  ENDDO
!
! Array B^^ for U^rot 
  U_RES_PP(:)=TO_FM(0)
  DO J=1, M 
     U_RES_PP(J) = HT(J) * (A_ELA_P + A_SEC_P/S(J)) - G(J)
  ENDDO
!
!
! --------------------------------------- GS on December 2019 
!  Write(*,*) "= V^rot"
! --------------------------------------- GS on December 2019 
!
! Elastic 
  V_ELA_P = (TO_FM(0) + LT_E) * A_ELA_P
!
! Secular 
  V_SEC_P = (TO_FM(0) + LT_F) * A_SEC_P
!
  Write(*,*) "- Elastic and secular V^ ", to_dp(V_ELA_P), to_dp(V_SEC_P)
!
! Array F
!  Write(*,*) "- Auxiliary Array F"
  F(:)=TO_FM(0)
  DO J=1, M-1 
     F(J)=TO_FM(0)
    DO JP=1, M 
      F(J)=F(J)+A_RES_P(J)*LT(JP)/(A_ROT(J)-S(JP))    
    ENDDO 
  ENDDO
!
! Array G
!  Write(*,*) "- Auxiliary Array G"
  G(:)=TO_FM(0)
  DO J=1, M 
     G(J)=TO_FM(0)
    DO JP=1, M-1 
      G(J)=G(J)+A_RES_P(JP)*LT(J)/(A_ROT(JP)-S(J))    
    ENDDO 
  ENDDO
!
! Array B^ for V^rot 
  V_RES_P(:)=TO_FM(0)
  DO J=1, M-1 
     V_RES_P(J) = (TO_FM(0) + LT_E)*A_RES_P(J) + F(J)
  ENDDO
!
! Array B^^ for V^rot 
  V_RES_PP(:)=TO_FM(0)
  DO J=1, M 
     V_RES_PP(J) = LT(J) * (A_ELA_P + A_SEC_P/S(J)) - G(J)
  ENDDO
!
! --------------------------------------------! GS 301119
! --------------------------------------------! GS 301119
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! /////////////////////////////////////////////
!
 Write(*,*) "= Copying ALL terms on file: ", trim(adjustl(FILE_PMTF))
!
 open (1,file=FILE_PMTF,status='unknown')
!
 Call DATE_AND_TIME (date,timc)
 Write(1,*) "-----------------------------------------------------------------"
 Write(1,*) trim(adjustl(FILE_PMTF)), &
           " ", date(1:4), '.', date(5:6), '.', date(7:8), & 
           " ", timc(1:2), '.', timc(3:4), '.', timc(5:6)
 Write(1,"(A8,A100)")  " Model: ", PATH_2_MODEL 
 Write(1,"(A7,1x,I6)") " MP:    ", MP
 Write(1,*) "-----------------------------------------------------------------"	   
 Write(1,*) ""
!
! A 
!
 Write(1,*) "PMTF:"
 Write(1,*) TO_DP(A_ELA), "A_e (elastic)" 
 Write(1,*) TO_DP(A_SEC), "A_s (secular)"
 DO I=1, M-1 
     if(i==1) Write(1,*) i, TO_DP(A_RES(I)), TO_DP(A_ROT(I)), "i, A_i, a_i (i=1,...,M-1)" 
     if(i/=1) Write(1,*) i, TO_DP(A_RES(I)), TO_DP(A_ROT(I))  
 ENDDO 
!
! A^
!
 Write(1,*) ""
 Write(1,*) "PMTF^:"
 Write(1,*) TO_DP(A_ELA_P), "A_e^ (elastic, primed)" 
 Write(1,*) TO_DP(A_SEC_P), "A_s^ (secular, primed)"
 DO I=1, M-1 
     if(i==1)Write(1,*) i, TO_DP(A_RES_P(I)), TO_DP(A_ROT(I)), "i, A_i^, a_i (i=1,...,M-1)" 
     if(i/=1)Write(1,*) i, TO_DP(A_RES_P(I)), TO_DP(A_ROT(I)) 
 ENDDO 
!
! A^^
!
 Write(1,*) ""
 Write(1,*) "A^^ terms: "
 DO I=1, M 
     if(i==1)Write(1,*) i, TO_DP(A_RES_PP(I)), TO_DP(S(I)), "i, A_i^^, s_i (i=1,...,M)" 
     if(i/=1)Write(1,*) i, TO_DP(A_RES_PP(I)), TO_DP(S(I)) 
 ENDDO 
!
! G^
!
 Write(1,*) ""
 Write(1,*) "G^ terms:"
 Write(1,*) TO_DP(G_ELA_P), "G_e^ (elastic, primed)" 
 Write(1,*) TO_DP(G_SEC_P), "G_s^ (secular, primed)"
 DO I=1, M-1 
     if(i==1)Write(1,*) i, TO_DP(G_RES_P(I)), TO_DP(A_ROT(I)), "i, G_i^, a_i (i=1,...,M-1)" 
     if(i/=1)Write(1,*) i, TO_DP(G_RES_P(I)), TO_DP(A_ROT(I)) 
 ENDDO 
!
! G^^
!
 Write(1,*) ""
 Write(1,*) "G^^ terms:"
 DO I=1, M 
     if(i==1)Write(1,*) i, TO_DP(G_RES_PP(I)), TO_DP(S(I)), "i, G_i^^, s_i (i=1,...,M)" 
     if(i/=1)Write(1,*) i, TO_DP(G_RES_PP(I)), TO_DP(S(I)) 
 ENDDO 
!
! U^
!
 Write(1,*) ""
 Write(1,*) "U^ terms:"
 Write(1,*) TO_DP(U_ELA_P), "U_e^ (elastic, primed)" 
 Write(1,*) TO_DP(U_SEC_P), "U_s^ (secular, primed)"
 DO I=1, M-1 
     if(i==1)Write(1,*) i, TO_DP(U_RES_P(I)), TO_DP(A_ROT(I)), "i, U_i^, a_i (i=1,...,M-1)" 
     if(i/=1)Write(1,*) i, TO_DP(U_RES_P(I)), TO_DP(A_ROT(I)) 
 ENDDO 
!
! U^^
!
 Write(1,*) ""
 Write(1,*) "U^^ terms:"
 DO I=1, M 
     if(i==1)Write(1,*) i, TO_DP(U_RES_PP(I)), TO_DP(S(I)), "i, U_i^^, s_i (i=1,...,M)" 
     if(i/=1)Write(1,*) i, TO_DP(U_RES_PP(I)), TO_DP(S(I)) 
 ENDDO 
!
! V^
!
 Write(1,*) ""
 Write(1,*) "V^ terms:"
 Write(1,*) TO_DP(V_ELA_P), "V_e^ (elastic, primed)" 
 Write(1,*) TO_DP(V_SEC_P), "V_s^ (secular, primed)"
 DO I=1, M-1 
     if(i==1)Write(1,*) i, TO_DP(V_RES_P(I)), TO_DP(A_ROT(I)), "i, V_i^, a_i (i=1,...,M-1)" 
     if(i/=1)Write(1,*) i, TO_DP(V_RES_P(I)), TO_DP(A_ROT(I)) 
 ENDDO 
!
! V^^
!
 Write(1,*) ""
 Write(1,*) "V^^ terms:"
 DO I=1, M 
     if(i==1)Write(1,*) i, TO_DP(V_RES_PP(I)), TO_DP(S(I)), "i, V_i^^, s_i (i=1,...,M)" 
     if(i/=1)Write(1,*) i, TO_DP(V_RES_PP(I)), TO_DP(S(I)) 
 ENDDO 
!
! H^ (G^ without direct centrifugal potential)
!
 Write(1,*) ""
 Write(1,*) "H^ terms:"
 Write(1,*) TO_DP(H_ELA_P), "H_e^ (elastic, primed)"
 Write(1,*) TO_DP(H_SEC_P), "H_s^ (secular, primed)"
 DO I=1, M-1
     if(i==1)Write(1,*) i, TO_DP(H_RES_P(I)), TO_DP(A_ROT(I)), "i, H_i^, a_i (i=1,...,M)"
     if(i/=1)Write(1,*) i, TO_DP(H_RES_P(I)), TO_DP(A_ROT(I))
 ENDDO
!
 CLOSE(1) 
! 
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! //////////////////////////////////////////////////////////////////
!
!
! \\\\\\\\\\\\\\\\\\\\  
  end program PMTF  
! \\\\\\\\\\\\\\\\\\\\
!
!
