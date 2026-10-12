!
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Copyright (C) 2003 Giorgio Spada 
!
!  
! TABOO is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or at your option) any later 
! version. 
!
! TABOO is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with TABOO.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! ###########################################################################
!
!
!
!---------------------------------------------
 MODULE TYPES 
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Module for the types of variables in TABOO
! Originally written by Carlo Giunchi 20??
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
 IMPLICIT NONE 
     SAVE
     INTEGER, PARAMETER  :: I4B   = SELECTED_INT_KIND(9)
     INTEGER, PARAMETER  :: SP    = KIND(1.0)
     INTEGER, PARAMETER  :: DP    = KIND(1.0D0)
     INTEGER, PARAMETER  :: QP    = DP  
 END MODULE TYPES 
!---------------------------------------------
! 
!
!
!---------------------------------------------
 MODULE DATA_1
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Basic declarations 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
 USE TYPES 
 USE FMZM
 IMPLICIT NONE 
 SAVE
! 
! ====== Begin of arameters formerly in file "TABOO.inc"
!
 INTEGER :: L_MIN
 INTEGER :: L_MAX
!
! ------ MULTI-precision parameter 
 INTEGER :: MP
!
! ------ Number of time steps and their length 
 INTEGER :: NN 
 REAL*4 :: DELTA_T
!
! ------ Earth model   
 INTEGER :: NV
 CHARACTER*100 :: MODEL_NAME
 CHARACTER*100 :: PATH_2_MODEL
!
! ====== End of parameters formerly un file "TABOO.inc"
!
! 
 CHARACTER*100, PARAMETER :: LOG_FILE_NAME = "taboo.log"
 Character*100, PARAMETER :: DE2_FILE_NAME = "l2-data.dat"
!  
 REAL*8, PARAMETER :: Newton_const=6.6732D-11
 REAL*8, PARAMETER :: PI=3.1415926535897932384D0
 INTEGER :: NROOTS
! 
 TYPE(FM), ALLOCATABLE :: RAD(:)  
 TYPE(FM), ALLOCATABLE :: RHO(:)  
 TYPE(FM), ALLOCATABLE :: RIG(:)  
 TYPE(FM), ALLOCATABLE :: VIS(:)  
 TYPE(FM), ALLOCATABLE :: GRA(:)  
 TYPE(FM), ALLOCATABLE :: ACO(:)  
!
 TYPE(FM) :: GGG, PIG, EMASS
!
 TYPE(FM) :: AC(6,6), AD(6,6), BC(6,6), BD(6,6)
 TYPE(FM), ALLOCATABLE :: CC(:,:,:,:)    !CC(0:2,1:NV,6,6) 
 TYPE(FM), ALLOCATABLE :: COEFMAT(:,:,:) !COEFMAT(0:2*NV,6,6)
 TYPE(FM), ALLOCATABLE :: CO(:)          !CO(0:3*2*NV)
 TYPE(FM), ALLOCATABLE :: RR(:,:,:)      !RR(0:2*NV,3,3)
 TYPE(FM) :: K0(6,6), K1(6,6), K2(6,6), A(6,6), B(6,6), C(6,6)
 TYPE(FM) :: D(6,6), ELAMAT (6, 6), SINIST_2(3,6), SINIST_1(3,6)
 TYPE(FM), ALLOCATABLE :: RRRR(:,:,:)    !RRRR(0:2*NV,3,6)
 TYPE(FM), ALLOCATABLE :: QQQQ(:,:,:)    !QQQQ(0:2*NV,3,6)
 TYPE(FM) :: CMB(6,3), BCS(3)
 TYPE(FM), ALLOCATABLE :: QQ(:,:,:)       !QQ(0:2*NV,3,3), 
 TYPE(FM), ALLOCATABLE :: AA(:)           !AA(NROOTS+1) 
 TYPE(FM), ALLOCATABLE :: NORM_COEFF(:,:) !NORM_COEFF(0:l_max,NROOTS+1)
 TYPE(FM), ALLOCATABLE :: RT1(:), RT2(:)  !RT1(1:NROOTS), RT2(1:NROOTS)
 TYPE(FM), ALLOCATABLE :: S(:,:)          !S(l_min:l_max,0:NROOTS)  
 TYPE(FM), ALLOCATABLE :: R_R(:,:,:), Q_Q(:,:,:)      !R_R(3,3,0:NROOTS), Q_Q(3,3,0:NROOTS),  
 TYPE(FM), ALLOCATABLE :: RADJ(:,:,:), QR(:,:,:)      !RADJ(3,3,0:NROOTS), QR (3,3,0:NROOTS)   
 TYPE(FM) :: COINF,  X_EL(3)
 TYPE(FM), ALLOCATABLE :: CTMP(:)         !CTMP(0:3*2*NV-1)  
 TYPE(FM), ALLOCATABLE :: DERPO(:)        !DERPO(1:NROOTS)
 TYPE(FM) :: XX(3)
 TYPE(FM), ALLOCATABLE :: XR(:,:)         !XR(3,NROOTS) 
 TYPE(FM) :: DETER ; EXTERNAL DETER  
!
 TYPE(FM), ALLOCATABLE :: H_E(:), L_E(:), K_E(:)    !H_E (0:l_max), L_E(0:l_max), K_E(0:l_max)  
 TYPE(FM), ALLOCATABLE :: H_F(:), L_F(:), K_F(:)    !H_F (0:l_max), L_F(0:l_max), K_F(0:l_max) 
 TYPE(FM), ALLOCATABLE :: H_V(:,:)                  !H_V (0:l_max,1:NROOTS) 
 TYPE(FM), ALLOCATABLE :: L_V(:,:)                  !L_V (0:l_max,1:NROOTS) 
 TYPE(FM), ALLOCATABLE :: K_V(:,:)                  !K_V (0:l_max,1:NROOTS) 	  
!
 TYPE(FM), ALLOCATABLE :: H_E_LOAD(:), L_E_LOAD(:), K_E_LOAD(:)  !H_E_LOAD (0:l_max), L_E_LOAD(0:l_max), K_E_LOAD(0:l_max)  
 TYPE(FM), ALLOCATABLE :: H_F_LOAD(:), L_F_LOAD(:), K_F_LOAD(:)  !H_F_LOAD (0:l_max), L_F_LOAD(0:l_max), K_F_LOAD(0:l_max) 
 TYPE(FM), ALLOCATABLE :: H_V_LOAD(:,:)                          !H_V_LOAD (0:l_max,1:NROOTS) 
 TYPE(FM), ALLOCATABLE :: L_V_LOAD(:,:)                          !L_V_LOAD (0:l_max,1:NROOTS) 
 TYPE(FM), ALLOCATABLE :: K_V_LOAD(:,:)                          !K_V_LOAD (0:l_max,1:NROOTS) 
!
 TYPE(FM), ALLOCATABLE :: H_E_TIDE(:), L_E_TIDE(:), K_E_TIDE(:)  !H_E_TIDE (0:l_max), L_E_TIDE(0:l_max), K_E_TIDE(0:l_max)  
 TYPE(FM), ALLOCATABLE :: H_F_TIDE(:), L_F_TIDE(:), K_F_TIDE(:)  !H_F_TIDE (0:l_max), L_F_TIDE(0:l_max), K_F_TIDE(0:l_max) 
 TYPE(FM), ALLOCATABLE :: H_V_TIDE(:,:)                          !H_V_TIDE (0:l_max,1:NROOTS) 
 TYPE(FM), ALLOCATABLE :: L_V_TIDE(:,:)                          !L_V_TIDE (0:l_max,1:NROOTS) 
 TYPE(FM), ALLOCATABLE :: K_V_TIDE(:,:)                          !K_V_TIDE (0:l_max,1:NROOTS) 
!
 REAL*8, ALLOCATABLE :: R_K(:,:), R_H(:,:), R_L(:,:), TEMP(:,:)  !R_K (0:l_max,1:NROOTS), R_H (0:l_max,1:NROOTS), R_L (0:l_max,1:NROOTS), TEMP(0:l_max,1:NROOTS) 
!
 END MODULE DATA_1 
!---------------------------------------------
!
!
!
!
!***************************
!***************************
 PROGRAM E_and_BETA_Arrays 
!***************************
!***************************
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This routine computes the "E" and the "Beta" arrays for use in SELEN4
!
! This is accomplished in four steps:
! 1) Input of the Earth model parameters 
! 2) Noprmalization of the parameters 
! 3) Evaluation of the Love numbers spectra for loading and tidal forcing
! 4) Computation of the "E" and the "Beta" arrays
! ... In progress ... GS December 2019  
! ...
! ...
! ...
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
 use DATA_1	
 use CONFIG_DATA
 USE FMZM
 IMPLICIT NONE 
 SAVE
!
! >>> Default file name for the convolved load-deformation coefficients 
 Character*100, PARAMETER :: FILE_S="ebs.dat", FILE_U="ebu.dat"
 Character*100, PARAMETER :: FILE_G="ebg.dat", FILE_V="ebv.dat"
!
 Character*20 DATE, TIMC 
!
 Integer(i4b) i, j, k, l, m, i_love
! 
 REAL*4, ALLOCATABLE :: E_S(:), E_U(:), E_G(:), E_V(:) ! E_S(0:l_max), E_U(0:l_max), E_G(0:l_max), E_V(0:l_max) 
!       	 
 REAL*4, ALLOCATABLE :: BETA_S(:,:), BETA_U(:,:)       ! BETA_S (0:l_max,0:NN+1), BETA_U (0:l_max,0:NN+1)  
 REAL*4, ALLOCATABLE :: BETA_G(:,:), BETA_V(:,:)       ! BETA_G (0:l_max,0:NN+1), BETA_V (0:l_max,0:NN+1) 
!
 CHARACTER*100 FILEX 
!
!==================================
!
! --- Configuration (NEW)
!
 call config
!
! --- Copy config parameters to TABOO variables
!
 l_min = 1
 l_max = cfg_l_max
 mp    = cfg_mp
 nv    = cfg_nv
 nn    = cfg_nn
 delta_t      = cfg_delta_t
 path_2_model = trim(adjustl(cfg_model_file)) 
 model_name   = path_2_model
!
 NROOTS=4*NV
!
 Call FM_SET(MP) 
!
 ALLOCATE( RAD(0:NV+1) )
 ALLOCATE( RHO(0:NV+1) )
 ALLOCATE( RIG(0:NV+1) )
 ALLOCATE( VIS(0:NV+1) )
 ALLOCATE( GRA(0:NV+1) )
 ALLOCATE( ACO(0:NV+1) )
!
 ALLOCATE( CC(0:2,1:NV,6,6) )
 ALLOCATE( COEFMAT(0:2*NV,6,6) )
 ALLOCATE( CO(0:3*2*NV) )
 ALLOCATE( RR(0:2*NV,3,3) )
 ALLOCATE( RRRR(0:2*NV,3,6), QQQQ(0:2*NV,3,6) )
 ALLOCATE( QQ(0:2*NV,3,3) ) 
 ALLOCATE( AA(NROOTS+1) ) 
 ALLOCATE( NORM_COEFF(0:l_max,NROOTS+1) )
 ALLOCATE( RT1(1:NROOTS), RT2(1:NROOTS) )
 ALLOCATE( S(l_min:l_max,0:NROOTS) ) 
 ALLOCATE( R_R(3,3,0:NROOTS), Q_Q(3,3,0:NROOTS) )
 ALLOCATE( RADJ(3,3,0:NROOTS), QR (3,3,0:NROOTS) )
 ALLOCATE( CTMP(0:3*2*NV-1) )
 ALLOCATE( DERPO(1:NROOTS) )
 ALLOCATE( XR(3,NROOTS) ) 
!
 ALLOCATE( H_E (0:l_max), L_E(0:l_max), K_E(0:l_max) )
 ALLOCATE( H_F (0:l_max), L_F(0:l_max), K_F(0:l_max) )
 ALLOCATE( H_V (0:l_max,1:NROOTS) )
 ALLOCATE( L_V (0:l_max,1:NROOTS) )
 ALLOCATE( K_V (0:l_max,1:NROOTS) )	  
!
 ALLOCATE( H_E_LOAD (0:l_max), L_E_LOAD(0:l_max), K_E_LOAD(0:l_max) )
 ALLOCATE( H_F_LOAD (0:l_max), L_F_LOAD(0:l_max), K_F_LOAD(0:l_max) )
 ALLOCATE( H_V_LOAD (0:l_max,1:NROOTS) )
 ALLOCATE( L_V_LOAD (0:l_max,1:NROOTS) )
 ALLOCATE( K_V_LOAD (0:l_max,1:NROOTS) )
!
 ALLOCATE( H_E_TIDE (0:l_max), L_E_TIDE(0:l_max), K_E_TIDE(0:l_max) )
 ALLOCATE( H_F_TIDE (0:l_max), L_F_TIDE(0:l_max), K_F_TIDE(0:l_max) )
 ALLOCATE( H_V_TIDE (0:l_max,1:NROOTS) )
 ALLOCATE( L_V_TIDE (0:l_max,1:NROOTS) )
 ALLOCATE( K_V_TIDE (0:l_max,1:NROOTS) )
!
 ALLOCATE( R_K (0:l_max,1:NROOTS) )
 ALLOCATE( R_H (0:l_max,1:NROOTS) )
 ALLOCATE( R_L (0:l_max,1:NROOTS) )
 ALLOCATE( TEMP(0:l_max,1:NROOTS) )
!
 ALLOCATE( E_S(0:l_max), E_U(0:l_max), E_G(0:l_max), E_V(0:l_max) )
 ALLOCATE( BETA_S (0:l_max,0:NN+1), BETA_U (0:l_max,0:NN+1) )
 ALLOCATE( BETA_G (0:l_max,0:NN+1), BETA_V (0:l_max,0:NN+1) )
!
! --- End of configuration
!
!
!
 OPEN(10,file=LOG_FILE_NAME,status="unknown")
!
 Call DATE_AND_TIME (date,timc)
 Write(10,*)  '>', ' date=', date(1:4), '.', date(5:6), '.', date(7:8), & 
                   ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)
!
!
 Call Input_Earth_Model 
!
 Call Normalization
!
 Do i_love=1, 0, -1  
! 
     If(i_love==1) Write(*,*) "> Computing the loading LNs at harmonic degree..."  
     If(i_love==0) Write(*,*) "> Computing the tidal LNs at harmonic degree..."  
! 
     Call Love_Numbers_Spectra (i_love) 
!    
 enddo 
!
 Call Print_Degree2_Love_Numbers  
!
 Call Print_Love_Numbers 
!
!
!
! >>> Computing the E(l) constants 
!
 E_G(:)=0d0
 E_U(:)=0d0
 E_V(:)=0d0 
! 
 E_S(:)=0d0
!
 do l=l_min, l_max 
     E_G(l) = to_dp(TO_FM(1.0) + k_e_load(l))/ float(2*l+1)
     E_U(l) = to_dp(TO_FM(0.0) + h_e_load(l))/ float(2*l+1)
     E_V(l) = to_dp(TO_FM(0.0) + l_e_load(l))/ float(2*l+1)
!
     E_S(l) = E_G(l)-E_U(l) 
!
 enddo
!
!
! >>> Computing the Beta(l,k) constants 
!
 BETA_G(:,:) = 0d0 
 BETA_U(:,:) = 0d0 
 BETA_V(:,:) = 0d0 
!
 BETA_S(:,:) = 0d0
!
 Do k=0,nn+1 
!
   Do l=l_min, l_max 
!
     BETA_G (l,k) = 0d0 
     BETA_U (l,k) = 0d0 
     BETA_V (l,k) = 0d0 
!
     BETA_S (l,k) = 0d0
!
     do m=1, nroots
!
         BETA_G (l,k) = BETA_G (l,k) + ( r_k(l,m))*(1d0 - exp(-k*delta_t/temp(l,m)) ) / (2d0*float(l)+1d0) 
         BETA_U (l,k) = BETA_U (l,k) + ( r_h(l,m))*(1d0 - exp(-k*delta_t/temp(l,m)) ) / (2d0*float(l)+1d0)
         BETA_V (l,k) = BETA_V (l,k) + ( r_l(l,m))*(1d0 - exp(-k*delta_t/temp(l,m)) ) / (2d0*float(l)+1d0)
!				    
     enddo 
!
         BETA_S (l,k) = BETA_G(l,k) - BETA_U(l,k)   			
!
   enddo
 enddo
!
!
! >>> Writing the arrays on the output files,
!     for later use with SELEN4 or other apps...
!
 Do 11 i=1, 4 
!
     if(i==1) filex = FILE_G 
     if(i==2) filex = FILE_U 
     if(i==3) filex = FILE_V 
     if(i==4) filex = FILE_S
!
 open(1,file=filex,status='unknown')
!
     Write(1,*) "-----------------------------------------------------------------"
     Write(1,*) 'File ', trim(adjustl(filex)), ' created by TABOO on ', & 
                 'day: ', date(1:4), '.', date(5:6), '.', date(7:8),   & 
                 ' time:', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
     Write(1,"(A8,A100)")  " Model: ", PATH_2_MODEL 
     Write(1,"(A7,1x,I6)") " L_MIN: ", L_MIN
     Write(1,"(A7,1x,I6)") " L_MAX: ", L_MAX 
     Write(1,"(A7,1x,I6)") " MP:    ", MP
     Write(1,"(A7,1x,I6)") " NN:    ", NN
     Write(1,*) "-----------------------------------------------------------------" 
!  
 do l=0, l_max
     if(i==1) write(1,*) l, E_G(l), (BETA_G(l,k), k=0, nn+1)  
     if(i==2) write(1,*) l, E_U(l), (BETA_U(l,k), k=0, nn+1)   
     if(i==3) write(1,*) l, E_V(l), (BETA_V(l,k), k=0, nn+1)  
     if(i==4) write(1,*) l, E_S(l), (BETA_S(l,k), k=0, nn+1)   
 enddo
!
   close(1)
!
11 Continue
!
!
 close(10) 
!
 STOP
 END PROGRAM E_and_BETA_Arrays
!
!
!
!
!
!
!============================= 
 SUBROUTINE PRINT_LOVE_NUMBERS
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The routine prints on files the Love numbers spectra
! for all harmonic degrees in the chosen range  
!
! Edited by GS on December 2019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
 USE TYPES
 Use DATA_1
 Use FMZM 
 IMPLICIT NONE
 SAVE
!
 Character*20 DATE, TIMC 
!
 CHARACTER*100 FILEX 
!
 CHARACTER*100, PARAMETER :: file_h_load = "h-load.dat"
 CHARACTER*100, PARAMETER :: file_k_load = "k-load.dat"
 CHARACTER*100, PARAMETER :: file_l_load = "l-load.dat"
!
 CHARACTER*100, PARAMETER :: file_h_tide = "h-tide.dat"
 CHARACTER*100, PARAMETER :: file_k_tide = "k-tide.dat"
 CHARACTER*100, PARAMETER :: file_l_tide = "l-tide.dat"
!
 CHARACTER*100, PARAMETER :: file_spectrum = "spectrum.dat"
!
 CHARACTER*100, PARAMETER :: FMT1="(I4,1X,200(1X,E20.8))"
 CHARACTER*100, PARAMETER :: FMT2="(I4,1X,5(E15.7,1X))"
!
 Integer i, j, k, l 
!
!
! >>> Storing the loading and tidal Love numbers 
!
 Do 1 i=1, 6
!
     if(i==1) filex = file_h_load 
     if(i==2) filex = file_k_load 
     if(i==3) filex = file_l_load 
!
     if(i==4) filex = file_h_tide 
     if(i==5) filex = file_k_tide 
     if(i==6) filex = file_l_tide
!
     open(10,file=filex,status="unknown")  
! 
     call DATE_AND_TIME (date,timc)
!
     Write(10,*) "-----------------------------------------------------------------"
     Write(10,*) 'File ', trim(adjustl(filex)), ' created by TABOO on ', & 
                 'day: ', date(1:4), '.', date(5:6), '.', date(7:8),   & 
                 ' time:', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
     Write(10,"(A8,A100)")  " Model: ", PATH_2_MODEL 
     Write(10,"(A7,1x,I6)") " L_MIN: ", L_MIN
     Write(10,"(A7,1x,I6)") " L_MAX: ", L_MAX 
     Write(10,"(A7,1x,I6)") " MP:    ", MP
     Write(10,*) "-----------------------------------------------------------------"
!
     Do l=l_min, l_max 
!
        if(i==1) write (10,FMT1) l, to_dp(h_e_load(l)), to_dp(h_f_load(l)), & 
	                           (to_dp(h_v_load(l,k)), k=1,nroots)   
        if(i==2) write (10,FMT1) l, to_dp(k_e_load(l)), to_dp(k_f_load(l)), & 
	                           (to_dp(k_v_load(l,k)), k=1,nroots)   
        if(i==3) write (10,FMT1) l, to_dp(l_e_load(l)), to_dp(l_f_load(l)), & 
	                           (to_dp(l_v_load(l,k)), k=1,nroots)   
!
     enddo
!
     Do l=2, 2 
        if(i==4) write (10,FMT1) l, to_dp(h_e_tide(l)), to_dp(h_f_tide(l)), & 
	                           (to_dp(h_v_tide(l,k)), k=1,nroots)   
        if(i==5) write (10,FMT1) l, to_dp(k_e_tide(l)), to_dp(k_f_tide(l)), & 
	                           (to_dp(k_v_tide(l,k)), k=1,nroots)   
        if(i==6) write (10,FMT1) l, to_dp(l_e_tide(l)), to_dp(l_f_tide(l)), & 
	                           (to_dp(l_v_tide(l,k)), k=1,nroots)   
     enddo
! 
   Close(10)
!
 1 Continue 
!
!
! >>> Computing and storing the "spectrum of relaxation" 
!
 open(3,file=file_spectrum,status="unknown")  
     Write(3,*) "-----------------------------------------------------------------"
     Write(3,*) 'File spectrum.dat, created by TABOO on ', & 
                 date(1:4), '.', date(5:6), '.', date(7:8), & 
	        ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  

     Write(3,"(A8,A100)")  " Model: ", PATH_2_MODEL 
     Write(3,"(A7,1x,I6)") " L_MIN: ", L_MIN
     Write(3,"(A7,1x,I6)") " L_MAX: ", L_MAX 
     Write(3,"(A7,1x,I6)") " MP:    ", MP
     Write(3,*) "-----------------------------------------------------------------"
     Write(3,*) '1st column: l = Harmonic degree                '
     Write(3,*) '2nd   "   : LOG10 (l)                          '
     Write(3,*) '3th   "   : s, (kyrs**(-1))                    '
     Write(3,*) '4rd   "   : LOG10(-s), with s in kyrs**(-1)    '
     Write(3,*) '5th   "   : Relaxation time = -1000./s, (yrs)  '
     Write(3,*) '6th   "   : LOG10 (Relaxation time (yrs))      '
     Write(3,*) "-----------------------------------------------------------------"
!
 Temp(:,:)=0d0
!
 Do l=l_min, l_max      
     Write(3,'(a1)') '>'
     Do k = 1, nroots        					       
     write(3,FMT2)      l,			   & 
                        log10 (float (l)),	   &
                        to_dp(s(l,k)),  	   &
		        log10(to_dp(-s(l,k))),     &
		        (-1000.0/to_dp(s(l,k))),   &
		        log10(-1000.0/to_dp(s(l,k)))
!
! >>> Relaxation time in kyrs
! 
     Temp(l,k)=(-1d0/to_dp(s(l,k)))
!
     enddo 
 enddo
  Close(3) 
!
!
! >>> Computing the negative of the "Normalized residues" for LLNs
!
 r_h(:,:)=0d0
 r_l(:,:)=0d0
 r_k(:,:)=0d0
!
 Do l=l_min, l_max 
     Do k=1, nroots
            r_h(l,k) = to_dp(-h_v_load(l,k)/s(l,k))
            r_l(l,k) = to_dp(-l_v_load(l,k)/s(l,k))
            r_k(l,k) = to_dp(-k_v_load(l,k)/s(l,k))
     enddo
 enddo
!
!
 END SUBROUTINE PRINT_LOVE_NUMBERS
!
!
!
!
!
!
!=========================================== 
 SUBROUTINE LOVE_NUMBERS_SPECTRA (I_LOVE)
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The routine computes the Love numbers spectra
! for all harmonic degrees in the chosen range 
!
! I_love=1 -> Load Love numbers 
! I_love=2 -> Tide Love numbers 
!
! Edited by GS on December 2019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
 USE TYPES
 Use DATA_1
 Use FMZM 
 IMPLICIT NONE
 SAVE
!
 INTEGER(I4B) I_LOVE, L, INU 
! INTEGER(I4B), PARAMETER :: NL=NV+1
 INTEGER(I4B):: NL
 TYPE(FM) :: CCC(6,6), EEE(6,6)
 TYPE(FM) :: PJ1(3,6), PJ2(3,6)
!
 Integer i, j, k, ii, jj, kk, m, jmin, jmax 
!
 NL=NV+1 
!
!
! Initialization of the Love numbers 
!
 h_e(:)   = TO_FM(0.0)
 h_f(:)   = TO_FM(0.0)
 h_v(:,:) = TO_FM(0.0)
!
 k_e(:)   = TO_FM(0.0)
 k_f(:)   = TO_FM(0.0)
 k_v(:,:) = TO_FM(0.0)
!
 l_e(:)   = TO_FM(0.0)
 l_f(:)   = TO_FM(0.0)
 l_v(:,:) = TO_FM(0.0)
!
 if(I_love==1) then 
    jmin=l_min
    jmax=l_max
 endif
 if(I_love==0) then 
    jmin=2
    jmax=2
 endif
!
!++++++++++++++++++++ 
  DO 1 L = JMIN, JMAX
!++++++++++++++++++++ 
!
!
  if(l==jmin.or.l==jmax) write(*,"(4X,A3,2X,I5)") "l= ", l  
  if(mod(l,10)==0)       write(*,"(4X,A3,2X,I5)") "l= ", l  
!
!
! >>> Propagator 
!
 DO 2 k=1, nl-1  
!
 call MATPROD (l, rad(nl-k), rad(nl-k-1), rho(nl-k), aco(nl-k), gra(nl-k), gra(nl-k-1)) 
!
 do i=1, 6
   do j=1, 6 
     k0(i,j) = ac(i,j) + bd(i,j) + bc(i,j)*rig(nl-k) + ad(i,j)/rig(nl-k)
     k1(i,j) = ad(i,j)/vis(nl-k)
     k2(i,j) = bc(i,j)*(-rig(nl-k)**2/vis(nl-k))  
   enddo
 enddo
!
 do i=1, 6
   do j=1, 6 
     cc(0,k,i,j) = k1(i,j)*(rig(nl-k)/vis(nl-k))  
     cc(1,k,i,j) = k0(i,j)*(rig(nl-k)/vis(nl-k))+k1(i,j)+k2(i,j)
     cc(2,k,i,j) = k0(i,j)   
   enddo
 enddo
!
 2 CONTINUE  
!
!
! >>> Elastic product
!
 Call DIRETTA (l, rad(nl),   rho(nl), aco(nl), gra(nl),   a, b)  
 Call INVERSA (l, rad(nl-1), rho(nl), aco(nl), gra(nl-1), c, d)
!
!
 Do i=1, 6
     Do j=1, 6 
       ccc(i,j) = a(i,j)+b(i,j)*rig(nl) 
       eee(i,j) = c(i,j)+d(i,j)/rig(nl)
     enddo
 enddo
 elamat = MATMUL(ccc,eee)  
!
! >>> Projectors 
!
  pj1(:,:) = TO_FM(0.0)
  pj2(:,:) = TO_FM(0.0)
!
 If(l==1) then 
     pj1(1,1) = TO_FM(1.0)
     pj1(2,2) = TO_FM(1.0)
     pj1(3,5) = TO_FM(1.0)
     pj2(1,3) = TO_FM(1.0)
     pj2(2,4) = TO_FM(1.0)
     pj2(3,5) = TO_FM(1.0)
 else 
     pj1(1,1) = TO_FM(1.0)
     pj1(2,2) = TO_FM(1.0)
     pj1(3,5) = TO_FM(1.0)
     pj2(1,3) = TO_FM(1.0)
     pj2(2,4) = TO_FM(1.0)
     pj2(3,6) = TO_FM(1.0)
 Endif
! 
 sinist_2 = MATMUL(pj2, elamat)
 sinist_1 = MATMUL(pj1, elamat)  
! 
 Call Promat 
!
! >>> Left products
!
 RRRR(:,:,:)=TO_FM(0.0)
 QQQQ(:,:,:)=TO_FM(0.0)
 Do inu = 0, 2*nv  
     RRRR (inu, :, :) = MATMUL(sinist_2, coefmat (inu,:,:))  
     QQQQ (inu, :, :) = MATMUL(sinist_1, coefmat (inu,:,:))  
 Enddo  
!
! >>> Core mantle boundary conditions 
!
 Call COREBO (l, rad(0), rho(0), aco(0), cmb)  
!
! >>> Right products
!
 RR(:,:,:)=TO_FM(0.0)
 QQ(:,:,:)=TO_FM(0.0)
 Do inu = 0, 2*nv  
     RR(inu, :, :) = MATMUL (rrrr (inu, :, :), cmb)  
     QQ(inu, :, :) = MATMUL (qqqq (inu, :, :), cmb)  
 Enddo  
!
! >>> Computing the determinant
!
 call det_tu
!
! Computing the coefficients of the secular polynomial
!
 Do i=1, nroots+1  
     aa(i)=co(6*nv+1-i)  
 enddo  
!
 Do i=1, nroots+1  
     norm_coeff(l,i)=aa(i)/aa(1) 
 enddo
!
! >>> Finding the roots of the secular polynomial 
!
 Call hqr_driver (aa, nroots, rt1, rt2)  
!
! Computing the 's' roots                         
!
 s(l,0) = -TO_FM('1E20')      
 Do k = 1, nroots
     s(l,k) = TO_FM(1.0)/rt1(k)  
 enddo 
!
! >>> Arrays R and Q
!
 R_R(:,:,:)=TO_FM(0.0)
 Q_Q(:,:,:)=TO_FM(0.0)   
! 
 Do m = 0, nroots
! 
 Do ii=1, 3 
     Do jj=1, 3   
          r_r (ii, jj, m) = rr (2 * nv, ii, jj)  
          q_q (ii, jj, m) = qq (2 * nv, ii, jj)        	
     enddo
 enddo
!
 Do kk = 2*nv-1, 0, -1  
     do ii=1,3 
     do jj=1,3  
	   r_r(ii, jj, m) = r_r (ii, jj, m) * s (l,m) + rr (kk, ii, jj)  
	   q_q(ii, jj, m) = q_q (ii, jj, m) * s (l,m) + qq (kk, ii, jj)  
     enddo
     enddo
 enddo 
!	  	
 enddo   
!
! >>> Adjoint of R
!
 Do m=0,nroots  
     RADJ (1,1,m) =   DETER (r_r (2, 2, m), r_r (2, 3, m), r_r (3, 2, m), r_r (3, 3, m) )
      RADJ (2,1,m) = - DETER (r_r (2, 1, m), r_r (2, 3, m), r_r (3, 1, m), r_r (3, 3, m) )
       RADJ (3,1,m) =   DETER (r_r (2, 1, m), r_r (2, 2, m), r_r (3, 1, m), r_r (3, 2, m) )
         RADJ (1,2,m) = - DETER (r_r (1, 2, m), r_r (1, 3, m), r_r (3, 2, m), r_r (3, 3, m) )
          RADJ (2,2,m) =   DETER (r_r (1, 1, m), r_r (1, 3, m), r_r (3, 1, m), r_r (3, 3, m) )
         RADJ (3,2,m) = - DETER (r_r (1, 1, m), r_r (1, 2, m), r_r (3, 1, m), r_r (3, 2, m) )
       RADJ (1,3,m) =   DETER (r_r (1, 2, m), r_r (1, 3, m), r_r (2, 2, m), r_r (2, 3, m) )
      RADJ (2,3,m) = - DETER (r_r (1, 1, m), r_r (1, 3, m), r_r (2, 1, m), r_r (2, 3, m) )
     RADJ (3,3,m) =   DETER (r_r (1, 1, m), r_r (1, 2, m), r_r (2, 1, m), r_r (2, 2, m) )
 Enddo  
!

! >>> Product Q times R
!
 Do m = 0, nroots  
     qr (:,:,m) = MATMUL (q_q (:,:,m),radj (:,:,m))  
 enddo  
!
!
! >>> Surface boundary conditions (loading: 1;  tidal: 0) 
!
 Call Boundary_Condition (l, i_love)
! 
! >>> Secular polynomial in the limit for --> -oo
! 
 coinf = co(3*2*nv)    
 Do k = 3*2*nv-1, 0,-1  
     coinf = coinf * s(l,0) + co(k)  
 enddo  
!
! >>> Elastic part of the solution vector
!
 Do i=1, 3
     x_el(i) = TO_FM(0.0) 
         Do j=1, 3 	
	    x_el(i)=x_el(i) + qr(i,j,0)*bcs(j)	
     enddo
 enddo	
 do i=1, 3    
     x_el(i) = x_el(i)/coinf 
 enddo
!
! >>> Derivative of the secular polynomial
!
 Do ii = 0, 3*2*nv-1  
     Ctmp(ii) = co (ii+1)*float(ii+1)  
 enddo  
 Do m = 1, nroots  
     derpo(m) = ctmp(3*2*nv-1)  
     Do k = 3*2*nv-2, 0, -1  
         derpo(m) = derpo(m)*s(l,m) + ctmp(k)  
     Enddo  
 enddo  
!
! >>> Viscoelastic residues
!
 Do m = 1, nroots  
     Do i=1, 3 	
     xx(i) = TO_FM(0.0) 
        Do j=1, 3 
            xx(i) = xx(i) + qr(i,j,m)*bcs(j)  
        enddo 
     xr(i,m)=xx(i)/derpo(m) 
   enddo     
 enddo 
!
! -----------------------------  
! >>> Love numbers h, l, and k  
! -----------------------------  
!
! >>> Love number 'h' 
!
 h_e(l) = x_el(1)*EMASS/rad(nl)  
 Do k=1, nroots 
     h_v(l,k) = xr(1,k)*EMASS/rad(nl)  
 enddo 
 h_f(l) = h_e(l)  
 Do m = 1, nroots
     h_f(l) = h_f(l)-h_v(l,m)/s(l,m)  
 enddo
!
! >>> Love number 'l'   
!
 l_e(l) = x_el(2)*EMASS/rad(nl)  
 Do k=1, nroots 
     l_v(l,k) = xr(2,k)*EMASS/rad(nl)  
 enddo 
 l_f(l) = l_e(l)  
 Do m = 1, nroots
     l_f(l) = l_f(l)-l_v(l,m)/s(l,m)  
 enddo
!
! >>> Love number 'k'  
!
 k_e (l) = -1._qp-(EMASS/rad(nl)/gra(nl))*x_el(3)  
 Do k=1, nroots 
     k_v (l,k) = -(EMASS/rad(nl)/gra(nl))*xr(3,k)  
 enddo 
 k_f(l) = k_e(l) 
 Do m = 1, nroots
     k_f(l) = k_f(l)-k_v(l,m)/s(l,m)  
 enddo
!
!
! -----------------------------  
! >>> Love numbers h, l, and k  
! ----------------------------- 
!
! >>> Defines the *loading Love numbers* for all degrees 
!
 If(i_love==1) then 
!  
 h_e_load(l) = h_e(l)
 h_f_load(l) = h_f(l)
 l_e_load(l) = l_e(l)
 l_f_load(l) = l_f(l)
 k_e_load(l) = k_e(l)
 k_f_load(l) = k_f(l)
 Do k=1, nroots
     h_v_load(l,k) = h_v(l,k)
     l_v_load(l,k) = l_v(l,k)
     k_v_load(l,k) = k_v(l,k)
 enddo
!
! >>> Defines the *tidal Love numbers* for all degrees 
!
 elseif(i_love==0) then 
! 
 h_e_tide(l) = h_e(l)
 h_f_tide(l) = h_f(l)
 l_e_tide(l) = l_e(l)
 l_f_tide(l) = l_f(l)
 k_e_tide(l) = k_e(l)
 k_f_tide(l) = k_f(l)
 Do k=1, nroots
     h_v_tide(l,k) = h_v(l,k)
     l_v_tide(l,k) = l_v(l,k)
     k_v_tide(l,k) = k_v(l,k)
 enddo
! 
 endif
!
!
! >>> End of the loop on the harmonic degrees 
!
!++++++++++++++++++++ 
1 CONTINUE 
!++++++++++++++++++++ 
!
 END SUBROUTINE LOVE_NUMBERS_SPECTRA
!
!
!
!
!
!
!======================================
 Subroutine Print_Degree2_Love_Numbers 
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This routine fills the file that contains
! degree l=2 data for later use - GS 141219
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
 USE DATA_1
 IMPLICIT NONE 
 Integer i, j, k  
!
 Character*20 DATE, TIMC 
 Integer, parameter :: ll=2
!
!
 OPEN(56,file=DE2_FILE_NAME,status="unknown")
!
 Call DATE_AND_TIME (date,timc)
!
  Write(56,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  Write(56,*) 'This is file: ', trim(adjustl(DE2_FILE_NAME))
  Write(56,*) 'created by TABOO.f90 (in multi-precision) on: ' 
  Write(56,*)  date(1:4), '.', date(5:6), '.', date(7:8), ' time=',timc(1:2), '.', timc(3:4), '.', timc(5:6)
  WRITE(56,*) 'The file contains l=2 Love numbers for the computation of the PMTF '
  Write(56,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
!
! >>> Filling the degree 2 file ! Filling the degree 2 file 
!
! >>> k_2 loading 
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a10)')'>',  'k_2_e_load'
     CALL FMWRIT (56,k_e_load(ll))
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a10)')'>',  'k_2_f_load' 
     CALL FMWRIT (56,k_f_load(ll))
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a20)')'>',  'k_2_i_load (kyrs^-1)' 
     do k=1, nroots    
         CALL FMWRIT (56,k_v_load(ll,k))
     enddo     
!
! >>> Normalized Cìcoefficients  
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a32)')'>',  'Coefficients of secular equation'    
     Do k=1, nroots +1	
        CALL FMWRIT (56,norm_coeff(ll,k))
     enddo 
!
! >>> Roots   
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a13)')'>',  's_i (kyrs^-1)'
     Do k=1, nroots 	
       CALL FMWRIT (56,s(ll,k))
     enddo
!
! >>> h_2 tidal    
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a10)')'>',  'h_2_e_tide'
     CALL FMWRIT (56,h_e_tide(ll))
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a10)')'>',  'h_2_f_tide' 
     CALL FMWRIT (56,h_f_tide(ll))		   
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a20)')'>',  'h_2_i_tide (kyrs^-1)' 
     Do k=1, nroots    
         CALL FMWRIT (56,h_v_tide(ll,k))
     enddo     
!
! >>> l_2 tidal    
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a10)')'>',  'l_2_e_tide'
     CALL FMWRIT (56,l_e_tide(ll))
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a10)')'>',  'l_2_f_tide' 
     CALL FMWRIT (56,l_f_tide(ll))		   
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a20)')'>',  'l_2_i_tide (kyrs^-1)' 
     Do k=1, nroots    
         CALL FMWRIT (56,l_v_tide(ll,k))
     enddo     
!
! >>> k_2 tidal    
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a10)')'>',  'k_2_e_tide'
     CALL FMWRIT (56,k_e_tide(ll))
!
     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a10)')'>',  'k_2_f_tide' 
     CALL FMWRIT (56,k_f_tide(ll))		   

     Write(56,'(a1)')'>'
     Write(56,'(a1,1x,a20)')'>',  'k_2_i_tide (kyrs^-1)' 
     do k=1, nroots    
         CALL FMWRIT (56,k_v_tide(ll,k))
     enddo     
!
! >>> Done Filling the degree 2 file 
!
 Close(56) 
!
 END Subroutine Print_Degree2_Love_Numbers 
! 
!
!
!
!
!
!=========================================== 
 SUBROUTINE NORMALIZATION
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This routine normalizes the inout model parameters
! Very old stuff, re-edited by GS on December 2019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
 Use DATA_1
 Use FMZM 
 Implicit NONE 
 SAVE
! 
 TYPE(FM) :: aux1, rad_ref, rho_ref, rig_ref, t_ref  
!
 Integer i, j, k 
!
 GGG=TO_FM(Newton_const)
 PIG=TO_FM(Pi)
!
! -------------------------------------------
!
 EMASS = 4._qp*pig*rho(0)*rad(0)**3/3._qp  
 do k = 1, nv+1  
     EMASS = EMASS + (4._qp/3._qp)*pig*(rad(k)**3-rad(k-1)**3)*rho(k)
 enddo 
     Write(10,*) "> The Earth mass (kg) is: ", to_sp(emass)
!
 gra(nv+1) = ggg*EMASS/rad(nv+1)/rad(nv+1)  
     Write(10,*) "> The surface gravity (m/s^2) is: ", to_sp(gra(nv+1))
!
 gra(0) = ggg*(4._qp/3._qp)*pig*rho(0)*rad(0)    
     Write(10,*) "> The gravity at the CMB (m/s^2) is: ", to_sp(gra(0))
! 
 gra(0) = ggg*(4._qp/3._qp)*pig*rho(0)*rad(0)  
 do k = 1, nv+1
     aux1 = gra(0)*rad(0)*rad(0)/ggg  
     do j = 1, k  
          aux1 = aux1+ (4._qp/3._qp)*pig*(rad(j)**3-rad(j-1)**3)*rho(j)
     enddo  
     gra(k) = aux1*ggg/rad(k)/rad(k)      
 enddo  
!
 Write(10,*) "> Gravity at the interfaces (m/s^2) from bottom to top" 
 do k=0, nv+1 
     Write(10,*) k, to_sp(gra(k))  
 enddo
!
! Write(*,*) "NORM> Setting up the normalization scheme" 
!
! Reference radius is Earth's radius 
     rad_ref = rad(nv+1)   
!
! Reference rigidity is the deepmost mantle layer rigidity 
     rig_ref = rig(nv+1)    
!
! Reference density is the deepmost mantle layer density
     rho_ref = rho(1)     
!
! Number of seconds in 1 kyr
     t_ref = 1000._qp*365.25_qp*24._qp*3600._qp  
!
!
!Write(*,*) "> Normalizing the model parameters" 
 Write(10,*)"> Normalizing the model parameters" 
 do k = 0, nv+1  
     rad (k) = rad (k) / rad_ref
     rho (k) = rho (k) / rho_ref
     rig (k) = rig (k) / rig_ref
     vis (k) = vis (k) / rig_ref / t_ref
 enddo
!
! Write(*,*) "NORM> Normalizing the gravity constant"  
 ggg = ggg*rho_ref**2*rad_ref**2/rig_ref  
!    
!
 gra(0) = ggg*(4._qp/3._qp)*pig*rho(0)*rad(0)  
!
 do k = 1, nv+1
     aux1 = gra(0)*rad(0)*rad(0)/ggg  
     do j = 1, k  
          aux1=aux1+(4._qp/3._qp)*pig*(rad(j)**3-rad(j-1)**3)*rho(j)
     enddo  
     gra(k) = aux1*ggg/rad(k)/rad(k) 
 enddo   
!
 do k = 0, nv+1  
     aco(k) =(4._qp/3._qp)*pig*ggg*rho(k)  
 enddo  
!
!
 Write(10,*) "> Normalized rad, rho, rig, vis, gra, from bottom to top" 
 do k=0, nv+1 
     Write(10,*) k, to_sp(rad(k)), to_sp(rho(k)), & 
                    to_sp(rig(k)), to_sp(vis(k)), & 
		    to_sp(gra(k))   
 enddo
!
!
 EMASS = 4._qp*pig*rho(0)*rad(0)**3/3._qp  
 do k = 1, nv+1  
      EMASS = EMASS + (4._qp/3._qp)*pig*(rad(k)**3-rad(k-1)**3)*rho(k)
 enddo 
 Write(10,*) "> Normalized mass of the Earth: ", to_dp(EMASS) 
!
 END SUBROUTINE NORMALIZATION
!-------------------------------------------- 
!
!
!
!
!
!
!=========================================== 
 SUBROUTINE INPUT_EARTH_MODEL
!
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This routine reads the model parameters from an external file
! Conventions: 
! rho(0), rig(0), vis(0) = Core density, rigidity and viscosity  
! rho(nv+1),  rig(nv+1),  vis(nv+1) = Lito density, rigidity and viscosity  
! rad(0) is the CMB radius 
! rad(nv+1) is Earth's radius 
! Edited by GS on december 2019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
 use DATA_1
 use FMZM 
 use config_data
 implicit NONE 
 SAVE
!
 INTEGER :: NH
 Character*100 JUNK
 Integer i, j, k 
 character*100 riga
 REAL(DP) dp_jun, dp_rad, dp_rho, dp_rig, dp_vis
!
 character*100, parameter :: backup_model_file="run-model-parameters.dat"
!
!
 write(*,*) ""  
 write(*,*) "> Reading Earth model data from file ", trim(adjustl(PATH_2_MODEL))
 write(*,*)
!
 OPEN(1,file=PATH_2_MODEL,status="old")
 Open(2,file=backup_model_file,status="unknown")
!OPEN(33,file="fort.33",status="unknown")
!
 call count_headers(1, nh)
 rewind(1)
!
 do i=1, nh 
 	read (1, "(a100)") junk
! 	write(33,"(a100)") junk	
	write(2, "(a100)") junk
 enddo
 do i=nv+1, 0, -1 
     read (1, ("(A100)"),end=12) riga 
     write(2, ("(A100)")       ) riga      
!        do j=1, 100 
!	  if (riga(j:j)=="'") riga(j:j)=" "
!        enddo 
!    write(33,"(A100)") riga 
 enddo  
12 rewind(1)
   close(2)
!12 close(1)
!   close(2)  
!   close(33)
!
! 
!OPEN(33,file="fort.33",status="old")
 do i=1, nh 
! 	read(33,"(a100)") junk
        read(1,*)
 enddo   
 do i=nv+1, 0, -1 
     !read(33,*) j, dp_jun, dp_rad, dp_rho, dp_rig, dp_vis  
     if    ( cfg_model_fmt==0 ) then
        read(1,*) j, dp_jun, dp_rad, dp_rho, dp_rig, dp_vis  
     elseif( cfg_model_fmt==1 ) then
        read(1,*) j, dp_rad, dp_rho, dp_rig, dp_vis       ! New model file format, DM 16092025
     else
        stop 'unknown model file format'
     end if
     rad(i) = to_fm(dp_rad)
     rho(i) = to_fm(dp_rho)
     rig(i) = to_fm(dp_rig)
     vis(i) = to_fm(dp_vis)
 enddo
 rig(0)=to_fm(0.0)
!Close(33) 
 close(1)
!
 Do k = 1, nv+1  
     If (rho(k)-rho(k-1)>0.0_qp)  then
     Write (*,*) "INPUT> *** WARNING *** density inversions exist!"
 Endif
 Enddo
!
 Write(10,*) "> Rad, rho, rig, and vis from bottom to top" 
 do k=0, nv+1 
     Write(10,*) k, to_sp(rad(k)), to_sp(rho(k)), & 
                    to_sp(rig(k)), to_dp(vis(k))
 enddo
!
! 
 END SUBROUTINE INPUT_EARTH_MODEL
! -------------------------------------------
!
!
!
!
!
!
!========================================== 
SUBROUTINE Boundary_Condition (L,CHAR)  
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
! This routine computes the forcing vector 
! (surface boundary condition) -- 
! GS Manu many years ago -- 
! Updated In march 2010 by GS for degree 1
! Retouched by GS on December 2019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!
 USE DATA_1 
 IMPLICIT NONE
 SAVE
! 
 INTEGER (I4B) :: char, l  
 TYPE(FM) :: g_amma  
!
 PIG=TO_FM(Pi)
!
 g_amma=(2._qp*float(l)+1._qp)/4._qp/pig/rad(nv+1)/rad(nv+1)  
!
 if (char==1) then  
     bcs(1) = -gra(nv+1)*g_amma  
     bcs(2) = 0._qp  
      If(l.ne.1) bcs(3)= -4._qp*pig*ggg*g_amma  
      If(l.eq.1) bcs(3)= 0._qp  
 endif  
!
 if (char==0) then
     bcs(1) =  0._qp  
     bcs(2) =  0._qp  
     bcs(3) = -4._qp*pig*ggg*g_amma  
 endif  
!
END SUBROUTINE Boundary_Condition
!-------------------------------------------
!
!
!
!
!
!
!========================================== 
SUBROUTINE DET_TU
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
! Computes the 'secular determinant'
! Retouched by GS on December 2019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!
USE DATA_1
IMPLICIT NONE  
SAVE
!
INTEGER (I4B) :: JJJ, IS, K, L, M, MU, NU
!
!
DO JJJ = 0, 3 * 2 * NV  
   CO (JJJ) = 0._QP  
     DO K = 0, 2 * NV  
        DO L = 0, 2 * NV  
            DO M = 0, 2 * NV  
                IF (K + L + M == JJJ) THEN  
                   CO (JJJ) = CO (JJJ) + &
		   RR (K, 1, 1) * RR (L, 2, 2) * RR (M, 3, 3) - &
		    RR (K, 1, 1) * RR (L, 2, 3) * RR (M, 3, 2) - &
		   RR (K, 1, 2) * RR (L, 2, 1) * RR (M, 3, 3) + &
		    RR (K, 1, 2) * RR (L, 2, 3) * RR (M, 3, 1) + &
		   RR (K, 1, 3) * RR (L, 2, 1) * RR (M, 3, 2) - &
		    RR (K, 1, 3) * RR (L, 2, 2) * RR (M, 3, 1)
                 ENDIF  
            ENDDO  
        ENDDO  
     ENDDO  
ENDDO
END SUBROUTINE DET_TU
!-------------------------------------------
!
!
!
!
!
!
!========================================== 
SUBROUTINE COREBO (N, R, RO, A, G)
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
! Boundary conditions at the Core-mantle boundary 
! Originally written by GS on 19xx 
! Retouched by GS lastly in December 2019
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!
 USE TYPES 
 USE FMZM
 IMPLICIT NONE 
 SAVE
!
 INTEGER (I4B) :: N
 TYPE(FM) :: G (6, 3)  
 TYPE(FM) :: R, RO, A 
!
 g = TO_FM(0.0)
!  
 g(1,1)= -(r**(n-1))/a  
 g(1,3)= 1._qp  
 g(2,2)= 1._qp  
 g(3,3)= ro*a*r  
 g(5,1)= r**n  
 g(6,1)= 2._qp*(float(n)-1._qp)*r**(n-1)  
 g(6,3)= 3._qp*a  
!
END SUBROUTINE COREBO
!------------------------------------------- 
!
!
!
!
!
!
!====================================================
SUBROUTINE MATPROD (N, RB, RD, RHU, ZA, GRX, GRY) 
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This routine performs the roduct between an 
! inverse and a direct fundamental matrix
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! - Created by SPI and GS @INGV on ???
! - Revised on April 2010 by GS
! - Ported by GS on July 2010 to multi-precision 
! - Touched by GS on December 2019 
!
USE DATA_1 
USE FMZM
IMPLICIT NONE
SAVE
!
INTEGER (I4B) :: I, J, N  
TYPE(FM) :: GRX, GRY, RB, RD, RHU, ZA
TYPE(FM) :: V(24)
TYPE(FM) :: H1, H2, H3, H4, H5
TYPE(FM) :: RN
TYPE(FM) :: A1, A2, A3, A4, B1, B2, B3, B4 
TYPE(FM) :: C1, C2, C3, C4, D1, D2, D3, D4
TYPE(FM) :: X,Y
!
! ####################################################
!
 ac = TO_FM(0.0)  
 ad = TO_FM(0.0)  
 bc = TO_FM(0.0) 
 bd = TO_FM(0.0)
!
 rn = 1._qp*n    
!
 x=rb
 y=rd
!
 a1= rn / (2._qp * (2._qp * rn + 3._qp) )  
 a2= 1._qp  
 a3= (rn + 1._qp) / (2._qp * (2._qp * rn - 1._qp) )  
 a4= 1._qp  
!
 b1= (rn + 3._qp) / (2._qp * (2._qp * rn + 3._qp) * (1._qp * rn + 1._qp) )
 b2= 1._qp / rn  
 b3= ( - rn + 2._qp) / (2._qp * rn * (2._qp * rn - 1._qp) )  
 b4= - 1._qp / (rn + 1._qp)  
!
 c1= (rn * rn - rn - 3._qp) / (2._qp * rn + 3._qp)  
 c2= 2._qp * (rn - 1._qp)  
 c3= ( - rn * rn - 3._qp * rn + 1._qp) / (2._qp * rn - 1._qp)  
 c4= - 2._qp * (rn + 2._qp)  
!
 d1= rn * (rn + 2._qp) / ( (2._qp * rn + 3._qp) * (rn + 1._qp) )  
 d2= 2._qp * (rn - 1._qp) / rn  
 d3= (rn * rn - 1._qp) / (rn * (2._qp * rn - 1._qp) )  
 d4= 2._qp * (1._qp * rn + 2._qp) / (rn + 1._qp)  
!
 h1= rn + 1._qp  
 h2= 1._qp + 2._qp * rn  
 h3= 2._qp * rn - 1._qp  
 h4= 3._qp + 2._qp * rn  
 h5= rn - 1._qp  
!
 v(1)= h1 / h2  
 v(2)= - 2._qp * h1 * (rn + 2._qp) / h2  
 v(3)= rn * h1 / 2._qp / h3 / h2  
 v(4)= - (rn**2 + 3._qp * rn - 1._qp) * rn / h2 / h3  
 v(5)= rn / h2  
 v(6)= 2._qp * rn * h5 / h2  
 v(7)= rn * h1 / 2._qp / h2 / h4  
 v(8)= h1 * (rn**2 - rn - 3._qp) / h2 / h4  
 v(9) = - 2._qp * rn * h1 * (2._qp + rn) / h2  
 v(10)= - h5 * rn * h1 * h1 / h3 / h2  
 v(11)= - 2._qp * h5 * rn * h1 / h2  
 v(12)= - rn**2 * h1 * (2._qp + rn) / h2 / h4  
 v(13)= - v (1)  
 v(14)= - v (3)  
 v(15)= - v (5)  
 v(16)= - v (7)  
 v(17)= - rn * h1 / h2  
 v(18)= (2._qp - rn) * h1 * rn / 2._qp / h3 / h2  
 v(19)= rn * h1 / h2  
 v(20)= rn * h1 * (3._qp + rn) / 2._qp / h2 / h4  
 v(21)= - v (1)  
 v(22)= - v (3)  
 v(23)= - v (5)  
 v(24)= - v (7)  
!
 ac (1,1) =  (a1*v(2)*(x/y)**(n+1) -a2*v(4)*(x/y)**(n-1) + & 
              a3*v(6)*(y/x)**(n) -a4*v(8)*(y/x)**(n+2))	      
 ac (1,2) = -(a1*v(9)*(x/y)**(n+1)-a2*v(10)*(x/y)**(n-1) + & 
              a3*v(11)*(y/x)**(n)-a4*v(12)*(y/x)**(n+2))
 ac (2,1) =  (b1*v(2)*(x/y)**(n+1) -b2*v(4)*(x/y)**(n-1) + & 
              b3*v(6)*(y/x)**(n) -b4*v(8)*(y/x)**(n+2))
 ac (2,2) = -(b1*v(9)*(x/y)**(n+1)-b2*v(10)*(x/y)**(n-1) + & 
              b3*v(11)*(y/x)**(n)-b4*v(12)*(y/x)**(n+2))
!
 ac (3,1) =  rhu*grx*ac(1,1)-(3._qp*rhu*za*y/(2._qp*rn+1._qp))* & 
                             ((x/y)**(n)-(y/x)**(n+1))
 ac (3,2) =  rhu*grx*ac(1,2)
 ac (3,5) =  rhu*(y/x)**(n+1)
 ac (3,6) =  (rhu*y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))
!
 ac (5,1) = -(3._qp*za*y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))
 ac (5,5) =  (y/x)**(n+1)
 ac (5,6) =  (y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))
!
 ac (6,1) =  3._qp*za*(a1*v(2)*(x/y)**(n+1)-(1+a2*v(4))*(x/y)**(n-1)+ & 
             a3*v(6)*(y/x)**(n) -a4*v(8)*(y/x)**(n+2))
 ac (6,2) =  3._qp*za*ac(1,2)
 ac (6,6) =  (x/y)**(n-1)
!
 ad (1,1) =  rhu*gry*( a1*v(1)*x*(x/y)**(n) -a2*v(3)*x*(x/y)**(n-2) + & 
                       a3*v(5)*y*(y/x)**(n) -a4*v(7)*y*(y/x)**(n+2))
 ad (1,3) =      (a1*v(13)*x*(x/y)**(n)-a2*v(14)*x*(x/y)**(n-2) + & 
                  a3*v(15)*y*(y/x)**(n)-a4*v(16)*y*(y/x)**(n+2))
 ad (1,4) =     -(a1*v(17)*x*(x/y)**(n)-a2*v(18)*x*(x/y)**(n-2) + & 
                  a3*v(19)*y*(y/x)**(n)-a4*v(20)*y*(y/x)**(n+2))
 ad (1,5) = -rhu*(a1*v(21)*x*(x/y)**(n)-a2*v(22)*x*(x/y)**(n-2) + & 
                  a3*v(23)*y*(y/x)**(n)-a4*v(24)*y*(y/x)**(n+2))
!
 ad (2,1) =  rhu*gry*(b1*v(1)*x*(x/y)**(n) -b2*v(3)*x*(x/y)**(n-2)  + & 
 		     b3*v(5)*y*(y/x)**(n) -b4*v(7)*y*(y/x)**(n+2))
 ad (2,3) =          (b1*v(13)*x*(x/y)**(n)-b2*v(14)*x*(x/y)**(n-2) + & 
		     b3*v(15)*y*(y/x)**(n)-b4*v(16)*y*(y/x)**(n+2))
 ad (2,4) =         -(b1*v(17)*x*(x/y)**(n)-b2*v(18)*x*(x/y)**(n-2) + & 
		     b3*v(19)*y*(y/x)**(n)-b4*v(20)*y*(y/x)**(n+2))
 ad (2,5) =     -rhu*(b1*v(21)*x*(x/y)**(n)-b2*v(22)*x*(x/y)**(n-2) + & 
		     b3*v(23)*y*(y/x)**(n)-b4*v(24)*y*(y/x)**(n+2))
!
 ad (3,1) =  rhu*grx*ad(1,1)
 ad (3,3) =  rhu*grx*ad(1,3)
 ad (3,4) =  rhu*grx*ad(1,4)
 ad (3,5) =  rhu*grx*ad(1,5)
!
 ad (6,1) =  3._qp*za*ad(1,1)
 ad (6,3) =  3._qp*za*ad(1,3)
 ad (6,4) =  3._qp*za*ad(1,4)
 ad (6,5) =  3._qp*za*ad(1,5)
!
 bc (3,1) =  (c1*v(2)*(1._qp/y)*(x/y)**(n) -c2*v(4)*(1._qp/y)*(x/y)**(n-2) + &
	      c3*v(6)*(1._qp/x)*(y/x)**(n) -c4*v(8)*(1._qp/x)*(y/x)**(n+2))	     
 bc (3,2) = -(c1*v(9)*(1._qp/y)*(x/y)**(n)-c2*v(10)*(1._qp/y)*(x/y)**(n-2) + & 
              c3*v(11)*(1._qp/x)*(y/x)**(n) -c4*v(12)*(1._qp/x)*(y/x)**(n+2))
 bc (4,1) =  (d1*v(2)*(1._qp/y)*(x/y)**(n) -d2*v(4)*(1._qp/y)*(x/y)**(n-2) + & 
              d3*v(6)*(1._qp/x)*(y/x)**(n) -d4*v(8)*(1._qp/x)*(y/x)**(n+2))
 bc (4,2) = -(d1*v(9)*(1._qp/y)*(x/y)**(n)-d2*v(10)*(1._qp/y)*(x/y)**(n-2) + & 
              d3*v(11)*(1._qp/x)*(y/x)**(n)-d4*v(12)*(1._qp/x)*(y/x)**(n+2))
!
 bd (3,1) =  rhu*gry*(c1*v(1)*(x/y)**(n) -c2* v(3)*(x/y)**(n-2)+ & 
                      c3*v(5)*(y/x)**(n+1) -c4*v(7)*(y/x)**(n+3))
 bd (3,3) =          (c1*v(13)*(x/y)**(n)-c2*v(14)*(x/y)**(n-2)+ & 
                      c3*v(15)*(y/x)**(n+1)-c4*v(16)*(y/x)**(n+3))
 bd (3,4) =         -(c1*v(17)*(x/y)**(n)-c2*v(18)*(x/y)**(n-2)+ & 
                      c3*v(19)*(y/x)**(n+1)-c4*v(20)*(y/x)**(n+3))
 bd (3,5) = -rhu*    (c1*v(21)*(x/y)**(n)-c2*v(22)*(x/y)**(n-2)+ & 
                      c3*v(23)*(y/x)**(n+1)-c4*v(24)*(y/x)**(n+3))
!
 bd (4,1) = rhu*gry*(d1*v(1)*(x/y)**(n) -d2* v(3)*(x/y)**(n-2)+ & 
                     d3*v(5)*(y/x)**(n+1) -d4*v(7)*(y/x)**(n+3))
 bd (4,3) =         (d1*v(13)*(x/y)**(n)-d2*v(14)*(x/y)**(n-2)+ & 
                     d3*v(15)*(y/x)**(n+1)-d4*v(16)*(y/x)**(n+3))
 bd (4,4) =        -(d1*v(17)*(x/y)**(n)-d2*v(18)*(x/y)**(n-2)+ & 
                     d3*v(19)*(y/x)**(n+1)-d4*v(20)*(y/x)**(n+3))
 bd (4,5) = -rhu*   (d1*v(21)*(x/y)**(n)-d2*v(22)*(x/y)**(n-2)+ & 
                     d3*v(23)*(y/x)**(n+1)-d4*v(24)*(y/x)**(n+3))
!
END SUBROUTINE MATPROD
!-------------------------------------------- 
!
!
!
!
!
!
!======================================== 
FUNCTION DETER (A,B,C,D)
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes the determinant od a 2x2 table
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Created by GS ???
! Edited by GS on December 2019
!
 USE TYPES
 USE FMZM
 IMPLICIT NONE 
 SAVE
 TYPE(FM) :: A, B, C, D, DETER  
!
  DETER = A*D-B*C  
!
 END FUNCTION DETER
!---------------------------------------- 
!
!
!
!
!
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
!
!
!=====================================================================
SUBROUTINE HQR_DRIVER (A, M, RTR, RTI)
!   
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is a driver for the root finding procedure 
! GS january 7, 2008. =based on netlib routines=
! Retouched by GS on December 2019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
 USE TYPES
 USE FMZM
 IMPLICIT NONE 
 SAVE 
!
! Uses balanc_nlib 
! Uses hqr_nlib 
! Uses bubble_sort_qp
!
!INTEGER, PARAMETER :: MAXM = 499 
!TYPE(FM) :: H (MAXM, MAXM), ASCALE(MAXM)   
 TYPE(FM) :: A (M + 1), RTR (M), RTI (M) 
 INTEGER  :: J, K, M, IGH, LOW, IERR  
 TYPE(FM) :: XR, XI   
 TYPE(FM), ALLOCATABLE :: H(:,:), ASCALE(:)
!
!If (m .gt.maxm.or.a(m+1).eq.0._qp) then 
 If (a(m+1).eq.0._qp) then
 Write(*,*) "Bad arguments in hqr_driver "
   Stop 11011
 Endif
!
 allocate( h(m,m), ascale(m) )
!
! Builds the Hessemberg matrix .....
 Do k = 1, m  
    h(1,k) = -a(m+1-k)/a(m+1)  
    do j = 2, m  
        h(j,k) = 0._qp  
    end do  
    if (k.ne.m) h (k+1,k) = 1._qp  
 enddo  
!
! Computes the roots by the "Eigenvalues method" .....
 call balanc_nlib  (m, h, low, igh, ascale)  
 call hqr_nlib     (m, low, igh, h, rtr, rti, ierr)
!
! (in normal conditions ierr=0)
!
! Ascending sorting of real parts of roots
 call bubble_sort_qp(m, rtr)
!
 deallocate(h, ascale)
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
!=====================================================================
 SUBROUTINE BALANC_NLIB(N,A,LOW,IGH,ASCALE)
 USE TYPES 
 USE FMZM
 IMPLICIT NONE
 SAVE
!
 INTEGER I,J,K,L,M,N,IGH,LOW,IEXC
 TYPE(FM) A(N,N),ASCALE(N)
 TYPE(FM) C,F,G,R,S,B2,ARADIX
 LOGICAL NOCONV
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Tested and ported to quadruple precision by GS on January 6 2008.  
!    Original is from http://www.netlib.org/eispack/3090/double/balanc.f
!    Retouched by GS on December 1019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     this subroutine is a translation of the algol procedure balance,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!     this subroutine balances a real matrix and isolates
!     eigenvalues whenever possible.
!     on input
!	 nm must be set to the row dimension of two-dimensional
!	   array parameters as declared in the calling program
!	   dimension statement.
!	 n is the order of the matrix.
!	 a contains the input matrix to be balanced.
!     on output
!	 a contains the balanced matrix.
!	 low and igh are two integers such that a(i,j)
!	   is equal to zero if
!	    (1) i is greater than j and
!	    (2) j=1,...,low-1 or i=igh+1,...,n.
!	 scale contains information determining the
!	    permutations and scaling factors used.
!     suppose that the principal submatrix in rows low through igh
!     has been balanced, that p(j) denotes the index interchanged
!     with j during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by d(i,j).  then
!	 scale(j) = p(j),    for j = 1,...,low-1
!		  = d(j,j),	 j = low,...,igh
!		  = p(j)	 j = igh+1,...,n.
!     the order in which the interchanges are made is n to igh+1,
!     then 1 to low-1.
!     note that 1 is returned for igh if igh is zero formally.
!     the algol procedure exc contained in balance appears in
!     balanc  in line.  (note that the algol roles of identifiers
!     k,l have been reversed.)
!     Questions and comments should be directed to Alan K. Cline,
!     Pleasant Valley Software, 8603 Altus Cove, Austin, TX 78759.
!     Electronic mail to cline@cs.utexas.edu.
!     this version dated january 1989. (for the IBM 3090vf)
!     ------------------------------------------------------------------
!
     aradix = 16._qp
!
     b2 = aradix * aradix
     k = 1
     l = n
     go to 100
!     .......... in-line procedure for row and
!		 column exchange ..........
  20 ascale(m) = j
     if (j .eq. m) go to 50
!
     do 30 i = 1, l
 	f = a(i,j)
 	a(i,j) = a(i,m)
 	a(i,m) = f
  30 continue
!
     do 40 i = k, n
 	f = a(j,i)
 	a(j,i) = a(m,i)
 	a(m,i) = f
  40 continue
!
  50 go to (80,130), iexc
!     .......... search for rows isolating an eigenvalue
!		 and push them down ..........
  80 if (l .eq. 1) go to 280
     l = l - 1
!     .......... for j=l step -1 until 1 do -- ..........
 100 do 120 j = l, 1, -1
!
 	do 110 i = 1, l
 	   if (i .eq. j) go to 110
 	   if (a(j,i) .ne. 0._qp) go to 120
 110	continue
!
 	m = l
 	iexc = 1
 	go to 20
 120 continue
!
     go to 140
!     .......... search for columns isolating an eigenvalue
!		 and push them left ..........
 130 k = k + 1
!
 140 do 170 j = k, l
!
 	do 150 i = k, l
 	   if (i .eq. j) go to 150
 	   if (a(i,j) .ne. 0._qp) go to 170
 150	continue
!
 	m = k
 	iexc = 2
 	go to 20
 170 continue
!     .......... now balance the submatrix in rows k to l ..........
     do 180 i = k, l
 180 ascale(i) = 1._qp
!     .......... iterative loop for norm reduction ..........
 190 noconv = .false.
!
     do 270 i = k, l
 	c = 0._qp
 	r = 0._qp
!
 	do 200 j = k, l
 	   if (j .eq. i) go to 200
 	   c = c + abs(a(j,i))
 	   r = r + abs(a(i,j))
 200	continue
!     .......... guard against zero c or r due to underflow ..........
 	if (c .eq. 0._qp .or. r .eq. 0._qp) go to 270
 	g = r / aradix
 	f = 1._qp
 	s = c + r
 210	if (c .ge. g) go to 220
 	f = f * aradix
 	c = c * b2
 	go to 210
 220	g = r * aradix
 230	if (c .lt. g) go to 240
 	f = f / aradix
 	c = c / b2
 	go to 230
!     .......... now balance ..........
 240	if ((c + r) / f .ge. 0.95_qp * s) go to 270
 	g = 1._qp / f
 	ascale(i) = ascale(i) * f
 	noconv = .true.
!
 	do 250 j = k, n
 250	a(i,j) = a(i,j) * g
!
 	do 260 j = 1, l
 260	a(j,i) = a(j,i) * f
!
 270 continue
!
        if (noconv) go to 190
!
  280 low = k
      igh = l
      return
      end
!---------------------------------------- 
!
!
!
!
!
!=====================================================================
 SUBROUTINE HQR_NLIB(N,LOW,IGH,H,WR,WI,IERR)
 USE TYPES
 USE FMZM
 IMPLICIT NONE
 SAVE
!c
!c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!c Tested and ported to quadruple precision by GS on January 6 2008.  
!c Original is from http://www.netlib.org/eispack/3090vf/single/hqr.f
!c Retouched by GS on December 2019. The Fortran 77 flavour untouched
!c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!c
!C  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
!c
 integer i,j,k,l,m,n,en,ll,mm,na,igh,itn,its,low,mp2,enm2,ierr
 type(fm) h(n,n),wr(n),wi(n), p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
 logical notlas
!c
!c     this subroutine is a translation of the algol procedure hqr,
!c     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!c     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!c     this subroutine finds the eigenvalues of a real
!c     upper hessenberg matrix by the qr method.
!c     on input
!c        nm must be set to the row dimension of two-dimensional
!c          array parameters as declared in the calling program
!c          dimension statement.
!c        n is the order of the matrix.
!	 low and igh are integers determined by the balancing
!	   subroutine  balanc.  if  balanc  has not been used,
!	   set low=1, igh=n.
!	 h contains the upper hessenberg matrix.  information about
!	   the transformations used in the reduction to hessenberg
!	   form by  elmhes  or  orthes, if performed, is stored
!	   in the remaining triangle under the hessenberg matrix.
!     on output
!	 h has been destroyed.  therefore, it must be saved
!	   before calling  hqr  if subsequent calculation and
!	   back transformation of eigenvectors is to be performed.
!	 wr and wi contain the real and imaginary parts,
!	   respectively, of the eigenvalues.  the eigenvalues
!	   are unordered except that complex conjugate pairs
!	   of values appear consecutively with the eigenvalue
!	   having the positive imaginary part first.  if an
!	   error exit is made, the eigenvalues should be correct
!	   for indices ierr+1,...,n.
!	 ierr is set to
!	   zero       for normal return,
!	   j	      if the limit of 30*n iterations is exhausted
!		      while the j-th eigenvalue is being sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!     this version dated september 1989.
! ---------------------------------------------------------------------
!
     ierr = 0
     norm = 0._qp
     k = 1
!     .......... store roots isolated by balanc
!		 and compute matrix norm ..........
     do 50 i = 1, n
!
 	do 40 j = k, n
  40	norm = norm + abs(h(i,j))
!
 	k = i
 	if (i .ge. low .and. i .le. igh) go to 50
 	wr(i) = h(i,i)
 	wi(i) = 0._qp
  50 continue
!
     en = igh
     t = 0._qp
     itn = 30*n
!     .......... search for next eigenvalues ..........
  60 if (en .lt. low) go to 1001
     its = 0
     na = en - 1
     enm2 = na - 1
!     .......... look for single small sub-diagonal element
!		 for l=en step -1 until low do -- ..........
  70 do 80 ll = low, en
 	l = en + low - ll
 	if (l .eq. low) go to 100
 	s = abs(h(l-1,l-1)) + abs(h(l,l))
 	if (s .eq. 0._qp) s = norm
 	tst1 = s
 	tst2 = tst1 + abs(h(l,l-1))
 	if (tst2 .eq. tst1) go to 100
  80 continue
!     .......... form shift ..........
 100 x = h(en,en)
     if (l .eq. en) go to 270
     y = h(na,na)
     w = h(en,na) * h(na,en)
     if (l .eq. na) go to 280
     if (itn .eq. 0) go to 1000
     if (its .ne. 10 .and. its .ne. 20) go to 130
!     .......... form exceptional shift ..........
     t = t + x
!
     do 120 i = low, en
 120 h(i,i) = h(i,i) - x
!
     s = abs(h(en,na)) + abs(h(na,enm2))
     x = 0.75_qp * s
     y = x
     w = -0.4375_qp * s * s
 130 its = its + 1
     itn = itn - 1
!     .......... look for two consecutive small
!		 sub-diagonal elements.
!		 for m=en-2 step -1 until l do -- ..........
     do 140 mm = l, enm2
 	m = enm2 + l - mm
 	zz = h(m,m)
 	r = x - zz
 	s = y - zz
 	p = (r * s - w) / h(m+1,m) + h(m,m+1)
 	q = h(m+1,m+1) - zz - r - s
 	r = h(m+2,m+1)
 	s = abs(p) + abs(q) + abs(r)
 	p = p / s
 	q = q / s
 	r = r / s
 	if (m .eq. l) go to 150
 	tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
 	tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
 	if (tst2 .eq. tst1) go to 150
 140 continue
!
 150 mp2 = m + 2
!
     do 160 i = mp2, en
 	h(i,i-2) = 0._qp
 	if (i .eq. mp2) go to 160
 	h(i,i-3) = 0._qp
 160 continue
!     .......... double qr step involving rows l to en and
!		 columns m to en ..........
     do 260 k = m, na
 	notlas = k .ne. na
 	if (k .eq. m) go to 170
 	p = h(k,k-1)
 	q = h(k+1,k-1)
 	r = 0._qp
 	if (notlas) r = h(k+2,k-1)
 	x = abs(p) + abs(q) + abs(r)
 	if (x .eq. 0._qp) go to 260
 	p = p / x
 	q = q / x
 	r = r / x
 170	s = sign(sqrt(p*p+q*q+r*r),p)
 	if (k .eq. m) go to 180
 	h(k,k-1) = -s * x
 	go to 190
 180	if (l .ne. m) h(k,k-1) = -h(k,k-1)
 190	p = p + s
 	x = p / s
 	y = q / s
 	zz = r / s
 	q = q / p
 	r = r / p
 	if (notlas) go to 225
!     .......... row modification ..........
 	do 200 j = k, EN
 	   p = h(k,j) + q * h(k+1,j)
 	   h(k,j) = h(k,j) - p * x
 	   h(k+1,j) = h(k+1,j) - p * y
 200	continue
!
 	j = min0(en,k+3)
!     .......... column modification ..........
 	do 210 i = L, j
 	   p = x * h(i,k) + y * h(i,k+1)
 	   h(i,k) = h(i,k) - p
 	   h(i,k+1) = h(i,k+1) - p * q
 210	continue
 	go to 255
 225	continue
!     .......... row modification ..........
 	do 230 j = k, EN
 	   p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
 	   h(k,j) = h(k,j) - p * x
 	   h(k+1,j) = h(k+1,j) - p * y
 	   h(k+2,j) = h(k+2,j) - p * zz
 230	continue
!
 	j = min0(en,k+3)
!     .......... column modification ..........
 	do 240 i = L, j
 	   p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
 	   h(i,k) = h(i,k) - p
 	   h(i,k+1) = h(i,k+1) - p * q
 	   h(i,k+2) = h(i,k+2) - p * r
 240	continue
 255	continue
!
 260 continue
!
     go to 70
!     .......... one root found ..........
 270 wr(en) = x + t
     wi(en) = 0._qp
     en = na
     go to 60
!     .......... two roots found ..........
 280 p = (y - x) / 2._qp
     q = p * p + w
     zz = sqrt(abs(q))
     x = x + t
     if (q .lt. 0._qp) go to 320
!     .......... real pair ..........
     zz = p + sign(zz,p)
     wr(na) = x + zz
     wr(en) = wr(na)
     if (zz .ne. 0._qp) wr(en) = x - w / zz
     wi(na) = 0._qp
     wi(en) = 0._qp
     go to 330
!     .......... complex pair ..........
 320 wr(na) = x + p
     wr(en) = x + p
     wi(na) = zz
     wi(en) = -zz
 330 en = enm2
     go to 60
!     .......... set error -- all eigenvalues have not
!		 converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
!---------------------------------------- 
!
!
!
!
!
!=========================================================
SUBROUTINE INVERSA (N, R, RHO, ZA, GRA, A, B)  
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes the "direct" propagator - GS 19??
! Retouched by GS on December 2019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
 USE TYPES 
 USE FMZM
 IMPLICIT NONE
 SAVE
!
 INTEGER (I4B) :: I, J, N  
 TYPE(FM) :: GRA, R, RHO, ZA
 TYPE(FM) :: A (6, 6), B (6, 6)  
 TYPE(FM) :: H1, H2, H3, H4, H5, RN
 TYPE(FM) :: V(24)
!
!
 a(:,:) = TO_FM(0.0)  
 b(:,:) = TO_FM(0.0) 
!
 rn = float (n)  
!
 h1 = rn + 1._qp  
 h2 = 1._qp + 2._qp * rn  
 h3 = 2._qp * rn - 1._qp  
 h4 = 3._qp + 2._qp * rn  
 h5 = rn - 1._qp  
!
 v(1) = h1 / h2  
 v(2) = - 2._qp * h1 * (rn + 2._qp) / h2  
 v(3) = rn * h1 / 2._qp / h3 / h2  
 v(4) = - (rn**2 + 3._qp * rn - 1._qp) * rn / h2 / h3  
 v(5) = rn / h2  
 v(6) = 2._qp * rn * h5 / h2  
 v(7) = rn * h1 / 2._qp / h2 / h4  
 v(8) = h1 * (rn**2 - rn - 3._qp) / h2 / h4  
 v(9) = - 2._qp * rn * h1 * (2._qp + rn) / h2  
 v(10) = - h5 * rn * h1 * h1 / h3 / h2  
 v(11) = - 2._qp * h5 * rn * h1 / h2  
 v(12) = - rn**2 * h1 * (2._qp + rn) / h2 / h4  
 v(13) = -v(1)  
 v(14) = -v(3)  
 v(15) = -v(5)  
 v(16) = -v(7)  
 v(17) = - rn * h1 / h2  
 v(18) = (2._qp - rn) * h1 * rn / 2._qp / h3 / h2  
 v(19) = rn * h1 / h2  
 v(20) = rn * h1 * (3._qp + rn) / 2._qp / h2 / h4  
 v(21) = -v(1)  
 v(22) = -v(3)  
 v(23) = -v(5)  
 v(24) = -v(7)  
 !
 a(1,1) = r**( - n - 1) * v (2)  
 a(2,1) = - r**( - n + 1) * v (4)  
 a(3,1) = 3._qp * za * r**( - n + 1) / (2._qp * float (n) + 1._qp)
 a(4,1) = r**(n) * v (6)  
 a(5,1) = - r**(n + 2) * v (8)  
 a(6,1) = - 3._qp * za * r**(n + 2) / (2._qp * float (n) + 1._qp)  
 !
 a(1,2) = - r**( - n - 1) * v (9)  
 a(2,2) = r**( - n + 1) * v (10)  
 a(4,2) = - r**(n) * v (11)  
 a(5,2) = r**(n + 2) * v (12)  
 !
 a(6,5) = - r**(n + 1)  
 a(3,6) = - r**( - n + 1) / (2._qp * float (n) + 1._qp)  
 a(6,6) = r**(n + 2) / (2._qp * float (n) + 1._qp)  
 !
 b(1,1) = r**( - n) * rho * gra * v (1)  
 b(2,1) = - r**( - n + 2) * rho * gra * v (3)  
 b(4,1) = r**(n + 1) * rho * gra * v (5)  
 b(5,1) = - r**(n + 3) * rho * gra * v (7)  
 !
 b(1,3) = r**( - n) * v (13)  
 b(2,3) = - r**( - n + 2) * v (14)  
 b(4,3) = r**(n + 1) * v (15)  
 b(5,3) = - r**(n + 3) * v (16)  
 !
 b(1,4) = - r**( - n) * v (17)  
 b(2,4) = r**( - n + 2) * v (18)  
 b(4,4) = - r**(n + 1) * v (19)  
 b(5,4) = r**(n + 3) * v (20)  
 !
 b(1,5) = - rho * r**( - n) * v (21)  
 b(2,5) = rho * r**( - n + 2) * v (22)  
 b(4,5) = - rho * r**(n + 1) * v (23)  
 b(5,5) = rho * r**(n + 3) * v (24)  
!
 END SUBROUTINE INVERSA
!---------------------------------------- 
!
!
!
!
!
!
!=========================================================
 SUBROUTINE DIRETTA (N, R, RHO, ZA, GRA, A, B)  
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes the "direct" propagator - GS 19??
! Retouched by GS on December 2019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
 USE TYPES 
 USE FMZM
 IMPLICIT NONE
 SAVE
!
 INTEGER (I4B) :: I, J, N  
 TYPE(FM) :: A (6, 6), B (6, 6)
 TYPE(FM) :: GRA, R, RHO, ZA, RN
 TYPE(FM) :: A1, A2, A3, A4, B1, B2, B3, B4 
 TYPE(FM) :: C1, C2, C3, C4, D1, D2, D3, D4
!
 a(:,:) = TO_FM(0.0) 
 b(:,:) = TO_FM(0.0) 
!
 rn = float (n)  
!
 a1 = rn / (2._qp * (2._qp * rn + 3._qp) )  
 a2 = 1._qp  
 a3 = (rn + 1._qp) / (2._qp * (2._qp * rn - 1._qp) )  
 a4 = 1._qp  
!
 b1 = (rn + 3._qp) / (2._qp * (2._qp * rn + 3._qp) * (1._qp * rn + 1._qp) )
 b2 = 1._qp / rn  
 b3 = ( - rn + 2._qp) / (2._qp * rn * (2._qp * rn - 1._qp) )  
 b4 = - 1._qp / (rn + 1._qp)  
!
 c1 = (rn * rn - rn - 3._qp) / (2._qp * rn + 3._qp)  
 c2 = 2._qp * (rn - 1._qp)  
 c3 = ( - rn * rn - 3._qp * rn + 1._qp) / (2._qp * rn - 1._qp)  
 c4 = - 2._qp * (rn + 2._qp)  
!
 d1 = rn * (rn + 2._qp) / ( (2._qp * rn + 3._qp) * (rn + 1._qp) )  
 d2 = 2._qp * (rn - 1._qp) / rn  
 d3 = (rn * rn - 1._qp) / (rn * (2._qp * rn - 1._qp) )  
 d4 = 2._qp * (1._qp * rn + 2._qp) / (rn + 1._qp)  
!
 a(1,1) = a1 * r** (n + 1)  
 a(1,2) = a2 * r** (n - 1)  
 a(1,4) = a3 * r** ( - n)  
 a(1,5) = a4 * r** ( - n - 2)  
 a(2,1) = b1 * r** (n + 1)  
 a(2,2) = b2 * r** (n - 1)  
 a(2,4) = b3 * r** ( - n)  
 a(2,5) = b4 * r** ( - n - 2)  
 a(3,1) = a (1, 1) * rho * gra  
 a(3,2) = a (1, 2) * rho * gra  
 a(3,3) = - rho * r** (n)  
 a(3,4) = a (1, 4) * rho * gra  
 a(3,5) = a (1, 5) * rho * gra  
 a(3,6) = - rho * r** ( - n - 1)  
 a(5,3) = a (3, 3) / rho  
 a(5,6) = a (3, 6) / rho  
 a(6,1) = 3._qp * za * a (1, 1)  
 a(6,2) = 3._qp * za * a (1, 2)  
 a(6,3) = - (2._qp * float (n) + 1._qp) * r** (n - 1)  
 a(6,4) = 3._qp * za * a (1, 4)  
 a(6,5) = 3._qp * za * a (1, 5)  
!
 b(3,1) = c1 * r** (n)  
 b(3,2) = c2 * r** (n - 2)  
 b(3,4) = c3 * r** ( - n - 1)  
 b(3,5) = c4 * r** ( - n - 3)  
!
 b(4,1) = d1 * r** (n)  
 b(4,2) = d2 * r** (n - 2)  
 b(4,4) = d3 * r** ( - n - 1)  
 b(4,5) = d4 * r** ( - n - 3)  
!
END SUBROUTINE DIRETTA
!---------------------------------------- 
!
!
!
!
!
!=========================================================
 SUBROUTINE PROMAT  
 USE TYPES
 USE DATA_1
 IMPLICIT NONE
 SAVE
!
 INTEGER (I4B) :: I, J, K, L, N, II, JJ 
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This routine performs the product between propagators
! GS on 19??
! Retouched by GS on December 2019
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
 TYPE(FM) :: AI (0:2*NV,6,6),BI (0:2,6,6),CI(0:2*NV,6,6) 
 TYPE(FM) :: AAA(6,6), BBB(6,6), CCC(6,6), DDD(6,6) 
! 
 ai (0:2, :, :) = cc (:,1,:,:)  
 if (nv/=1) then  
    do i = 2, nv      
    n = 2 * i  
        do k = 0, n  
        ci (k, :, :) = to_fm(0.0)  
        bi (:, :, :) = cc (:, i, :, :)  
            do l = 0, n - 2  
                do j = 0, 2  
                    if (l + j == k) then 
		    ccc(:,:) = matmul(ai(l,:,:),bi(j,:,:)) 									
			do ii=1, 6 
			do jj=1, 6 
			     ci(k,ii,jj) = ci(k,ii,jj) + ccc(ii,jj) 
			enddo
			enddo
!
		    endif
!
                enddo  
            enddo  
        enddo  
    ai (0:n, :, :) = ci (0:n, :, :)  
    enddo  
!
 do ii=0, 2*nv     
     coefmat(ii,:,:) = ci(ii,:,:) 
 enddo 
!
 else  
!
     coefmat (0:2, :, :) = ai (0:2, :, :)  
!
 endif  
!
 END SUBROUTINE PROMAT
!---------------------------------------- 
!
!
!
!
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Copyright (C) 2003 Giorgio Spada 
!
!  
! TABOO is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or at your option) any later 
! version. 
!
! TABOO is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with TABOO.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! ###########################################################################
!
