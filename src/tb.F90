!
!
! This is program "TB.F90" 
!
! Last modified GS 04-11-2008 [INTEL port]
! Modified on August 2008, for version 2.7 (horizontal movements) 
! Also touched October 2008 for the PAULSON et al. profile
! ----->>> New as from August 2008: This program also computes arrays 
! Beta^v(l,k) and E^v(l), pertaining to the horizontal movements.
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Reviewed GS May 2010 - Porting under g95 
! *** Reviewed GS May 30 2010 - Degree 1 interface 
! *** Reviewed DM Dec 22 2010 - Adjusted REAL sizes for gfortran compatibility
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Copyright (C) 2008 Giorgio Spada, Florence Colleoni, and Paolo Stocchi 
!
! This file is part of SELEN. 
!  
! SELEN is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or at your option) any later 
! version. 
!
! SELEN is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! =============================================================================
!
! Computes the arrays Beta(l,k) and E(l), and Beta^u(l,k) and E^u(l), the only 
! rheology--dependent quantities involved in the solution of the SLE by the PS 
! method.  Here the program TABOO is used as a subroutine to compute the LDCs. 
!	
!
! Input file:
!	- task_1.dat 
!         [built here according to the input parameters of "data.inc"] 
!
! Output files:  
! 	- spectrum.dat
!	- ih.dat, il.dat, & ik.dat (ihh.dat, ill.dat, & ikk.dat
!	- h.dat,   l.dat, &  k.dat ( hh.dat,  ll.dat, &  kk.dat)
!	- ebs.dat    ! Sea level green function 
!	- ebu.dat    ! Vertical displacement Green function
!       - ebn.dat    ! Geoid elevation Green function
!       - ebv.dat    ! Hhorizontal displacement Green function
!
! =============================================================================
!
!

! INCLUDE "harmonics.f90"
 PROGRAM BB
 IMPLICIT NONE 
 INCLUDE "data.inc"
!
 INTEGER, PARAMETER :: L1=1, LMP1=L1+1, L2=LMAX
!
 INTEGER, PARAMETER :: NV_MAX=9             ! max. number of viscoelastic layers 
 INTEGER, PARAMETER :: NROOTS_MAX=4*NV_MAX  ! max. number of roots    
 INTEGER, PARAMETER :: NIMP=100
 INTEGER J, L, K, M, P, KK, NOUT, NROOTS  
 REAL*4 R_K (0:LMAX,NROOTS_MAX), & 
        R_H (0:LMAX,NROOTS_MAX), &
	R_L (0:LMAX,NROOTS_MAX), & 
        TEMP(0:LMAX,NROOTS_MAX), & 
        H_E (0:LMAX), &
	K_E (0:LMAX), &
	L_E (0:LMAX), & 
	VIS_LM, VIS_UM, VIS_TZ, VIS, JUNK       
 REAL*8 VSC(NV)
! 
! Imported from AM.F90:---------
 INTEGER ND
 REAL*4 ES(0:LMAX), & 
        EU(0:LMAX), &
        EN(0:LMAX), &  
        EV(0:LMAX) 
 REAL*4 BETAS(0:LMAX,0:NN), & 
        BETAU(0:LMAX,0:NN), & 
        BETAN(0:LMAX,0:NN), & 
        BETAV(0:LMAX,0:NN)  
 REAL*4  DEN
! ------------------------------
 REAL*8 LTHIC 
 CHARACTER*1 CHAR       
 CHARACTER*100 SS(NIMP)
 CHARACTER*200 ROW
!
! =============================================================================
!
!
! # Reading the lithospheric thickness & viscosity profile from "visco.dat"
!
 open(12,file=visco_model,status='unknown') 
 read(12,'(a200)') row 
      call scan_string (row, 1, ss, nout)
      call CHAR100_2_REAL(ss(1), LTHIC)	
 do j=1, nv 
     read(12,'(a200)') row
     call scan_string (row, 1, ss, nout)
     call CHAR100_2_REAL(ss(1),vsc(j))
 enddo 
 close(12)
!
!
! # Building the input file <<task_1.dat>>... 
!
 open(1,file='task_1.dat',status='unknown') 
!  
! General settings
 write(1,'(a6)' )'Active' 
 write(1,'(a16)')'Harmonic_Degrees'
 write(1,*) l1, l2   ! min and max degrees  
 write(1,*) '0'      ! verbose mode 
 write(1,*) '1'      ! i_loading 
!
! Earth model 
 write(1,'(a10)')'Make_Model' 
!
 NROOTS=4*NV   	     ! Number of Roots = 4*Number of v.e. layers  
 write(1,*)  NV      ! 
 write(1,*)  CDE     ! 
 write(1,*)  LTHIC   ! thickness of the lithosphere (km) 
 write(1,*) '0' 
! 
 do j=1, nv
     write(1,*) vsc(j)
 enddo
!
! Requesting normalized residues 
 write(1,'(a19)')'Normalized_Residues'  
 do k=1,3
     write(1,*) '1' 
 enddo
!
! Requesting LDCs 
 write(1,'(a15)')'El_Fluid_Viscel'  
 do k=1,3
     write(1,*) '1' 
 enddo
 close(1) 
!
!
!  
  write(*,*) '    - Calling TABOO' 
! ******************************* 
!
  call TABOO
!
! *******************************

!
!

!
!
! Initialization of the "BETA" arrays
  BETAS(:,:)=0. 
  BETAU(:,:)=0. 
  BETAN(:,:)=0.
  BETAV(:,:)=0.
!
!
! Reading the normalized residues and relaxation times
   open(1,file='ih.dat',status='old')
   open(2,file='ik.dat',status='old')
   open(3,file='il.dat',status='old')
   open(4,file='spectrum.dat',status='old')
!
   do k=1,2
           read(1,'(a1)') char
           read(2,'(a1)') char
           read(3,'(a1)') char	 
   enddo
   do k=1,7
	   read(4,'(a1)') char
   enddo
!
   r_h(:,:) =0.
   r_k(:,:) =0.
   r_l(:,:) =0. 
   temp(:,:)=0. 
!
        do l = L1, L2
!
	   do j=1, 4  
           	read(j, '(a1)') char
	   enddo
!     
	   do m=1, nroots
	      read  (1, '(i4,1x,e15.7)') p,  r_h (l,m)
	      read  (2, '(i4,1x,e15.7)') p,  r_k (l,m)
	      read  (3, '(i4,1x,e15.7)') p,  r_l (l,m)
	      read  (4, '(i4,1x,5(e15.7,1x))') kk, junk, junk, junk, temp(l,m)
!
              r_h(l,m)=-r_h(l,m)      		! **  I have (h/s) in memory  ** !
	      r_k(l,m)=-r_k(l,m)      		! **  I have (k/s) in memory  ** !
	      r_l(l,m)=-r_l(l,m)      		! **  I have (k/s) in memory  ** !
!	      
	      if(imode==2.or.imode==3.or.imode==4) then 
	      					r_h(l,m)=0.
	      					r_k(l,m)=0.
						r_l(l,m)=0.
			   			   endif 
	      if(imode==7)			   then
	      				        r_k(l,m)=0.
						   endif
!	      	      
	      temp(l,m) = temp(l,m)/1000.    	! ** relaxation times in kyrs ** !
!	      
	   enddo
        enddo
	close(1); close(2); close(3); close(4)  
!
!
!
!
!--- beta.f: reading the elastic LDCs... 
!
   OPEN(1,FILE='h.dat',STATUS='old')
   OPEN(2,FILE='k.dat',STATUS='old')
   OPEN(3,FILE='l.dat',STATUS='old')
!
   do J=1, 2
      READ(1,'(a1)') char
      READ(2,'(a1)') char
      READ(3,'(a1)') char
   enddo
!
   do l=L1, L2
      	READ (1,'((i4,1x,24(1x,e20.8)))') p, h_e(l)
      	READ (2,'((i4,1x,24(1x,e20.8)))') p, k_e(l)
      	READ (3,'((i4,1x,24(1x,e20.8)))') p, l_e(l)
   enddo
!
   CLOSE(1)
   CLOSE(2)
   CLOSE(3) 
!
!
!###################################################
! Building the "E" arrays (Elastic Green functions)
!###################################################
!
! Initialization of the "E" arrays 
   ES(:)=0.   
   EU(:)=0.
   EN(:)=0.
   EV(:)=0.
!
! --------------------------------------------------------------
! ====== Setting the ELASTIC Green Functions for Degree 1 ======
! --------------------------------------------------------------
!
 IF(DEG1==1) THEN 				
!
! ... if the degree 1 needs to be included in SLE modeling, the Green 
! functions at degree nd=1 are set by the corresponding Love numbers 
!
   nd=1
!
! ---- GSC or Elastic solution
	If(imode==1.or.imode==2.or.imode==5)  then 
	EU(nd) = h_e(nd) 
        EN(nd) = 0. 
        EV(nd) = l_e(nd)
	ES(nd) = EN(nd) - EU(nd)	
	Endif
!
! ---- Eustatic or solution 
	If(imode==3)  then 
		EU(nd) = 0.
		EN(nd) = 0. 
		EV(nd) = 0. 	
	        ES(nd) = 0.
	Endif 
!
! ---- Woodward solution 
	If(imode==4)  then 
		EU(nd) = 0.
		EN(nd) = 1. 
		EV(nd) = 0. 	
	        ES(nd) = 0.
	Endif 	
!
   den=2.*float(nd)+1. 	
!
   EU(nd) = EU(nd)/den 
   EN(nd) = EN(nd)/den 
   EV(nd) = EV(nd)/den 
   ES(nd) = ES(nd)/den 	
!
 ENDIF 
!
!
! ----------------------------------------------------------------
! ====== Setting the ELASTIC Green Functions for Degree > 1 ======
! ----------------------------------------------------------------
! 
 Do 10 l=LMP1, LMAX
!
	EN(l) =   1.+k_e(l)
	EU(l) =      h_e(l)    
	EV(l) =      l_e(l) 
      	ES(l) = EN(l)-EU(l)    
!       
! ---- Eustatic solution 
	if(imode==3) then 
			EN(l) = 0.
			EU(l) = 0.
			EV(l) = 0.
	                ES(l) = EN(l)-EU(l) 
		     endif 
!
! ---- Woodward solution			
	if(imode==4) then 
			EN(l) = 1.
			EU(l) = 0.
			EV(l) = 0.
	                ES(l) = EN(l)-EU(l) 
		     endif
!
        den=	(2.*float(l)+1.) 
!
        EU(l) = EU(l)/den 
	EN(l) = EN(l)/den 
	EV(l) = EV(l)/den 
	ES(l) = ES(l)/den 
!
10 Continue
!
!
!
!############################################################
! Building the "BETA" arrays (Visco-Elastic Green functions)
!############################################################
!
!
 BETAS(:,:)=0. 
 BETAU(:,:)=0. 
 BETAN(:,:)=0.
 BETAV(:,:)=0.
!
!
! --------------------------------------------------------------------
! ====== Setting the Visco-Elastic Green Functions for Degree 1 ======
! --------------------------------------------------------------------
!
! ... if the degree 1 needs to be included in SLE modeling, the Green 
! functions at degree nd=1 are set by the corresponding Love numbers 
!
 IF(DEG1==1) THEN 
!
	do k=0, NN 
!
	l=1 
!
                BETAN(l,k)=0. 
                BETAU(l,k)=0. 
                BETAV(l,k)=0. 
!
   		den = 2.*float(l)+1. 
!
        		do m=1, nroots				    
! 
			BETAU(l,k) = BETAU(l,k) + & 
			            (r_h(l,m))* & 
			            (1. - exp(-k*delta/temp(l,m)))/den
!
			BETAV(l,k) = BETAV(l,k) + & 
			            (r_l(l,m))* & 
			            (1. - exp(-k*delta/temp(l,m)))/den
!				    
			enddo
!			
	BETAS(l,k) = BETAN(l,k) - BETAU(l,k)		
!			
	enddo
!	
 ENDIF
!
!
! ----------------------------------------------------------------------
! ====== Setting the Visco-Elastic Green Functions for Degree > 1 ======
! ----------------------------------------------------------------------
!
	do k=0, NN 
!
        do l=lmp1, lmax 
!
                BETAU(l,k)=0.0 
                BETAN(l,k)=0.0 
                BETAV(l,k)=0.0 
!
   	        den = 2.*float(l)+1.  
!
        		do m=1, nroots
!
			BETAN(l,k) = BETAN(l,k) + & 
			            (r_k(l,m))* & 
			            (1.0 - exp(-k*delta/temp(l,m)))/ den 
! 
			BETAU(l,k) = BETAU(l,k) + & 
			            (r_h(l,m))* & 
			            (1.0 - exp(-k*delta/temp(l,m)))/ den 
!
			BETAV(l,k) = BETAV(l,k) + & 
			            (r_l(l,m))* & 
			            (1.0 - exp(-k*delta/temp(l,m)))/ den 
!				    
			enddo
!			
	      BETAS(l,k)=BETAN(l,k) - BETAU(l,k)   			
!			
	      enddo
	enddo


!
! ----------------------------------- 
! ====== Reporting the results ======
! ----------------------------------- 
!
!
 open(1,file='ebs.dat',status='unknown')
 open(2,file='ebu.dat',status='unknown')
 open(3,file='ebn.dat',status='unknown')
 open(4,file='ebv.dat',status='unknown')
!  
 do l=0, lmax
 	 write(1,*) l, ES(l), (betas(l,k), k=0, nn)
 	 write(2,*) l, EU(l), (betau(l,k), k=0, nn)
 	 write(3,*) l, EN(l), (betan(l,k), k=0, nn)		 
 	 write(4,*) l, EV(l), (betav(l,k), k=0, nn)
 enddo
!
 close(1)
 close(2)
 close(3)
 close(4) 
!
!
!
!================
  end program BB 
!================
!
!
!
!
!
!   +------------------------------------------------------+
!   |                                                      |
!   |                                                      |
!   |                +---++---++---++------+               |
!   |                | T || A || B ||  OO  |               |
!   |                +---++---++---++------+               |
!   |                                                      |
!   |        ( a posT glAcial reBOund calculatOr )         |
!   |							   |
!   |                by Giorgio Spada et al. 	           |
!   |                                                      |
!   |							   |	 
!   |              Release 2.0 -  January 2005             |
!   |                                                      |
!   |     						   |
!   |               http://samizdat.mines.edu 		   |
!   |					 	           |
!   |                                                      |
!   +------------------------------------------------------+
!                         
!
!
!
 MODULE STRATA 
!
! Module for the types of variables in TABOO
! Originally written by Carlo Giunchi  
! Revised GS May 2010 for g95 implementation 
! Reviewed DM Dec 22 2010 - Adjusted REAL sizes for gfortran compatibility
!
 IMPLICIT NONE 
 SAVE
 Integer, parameter  :: i4b   = selected_int_kind(9)
 Integer, parameter  :: sp    = kind(1.0)
 Integer, parameter  :: dp    = kind(1.0d0)
#ifdef SUPPORTS_QUAD_16
 Integer, parameter  :: qp    = kind(1.0q0)
#elif SUPPORTS_QUAD_10
 Integer, parameter  :: qp    = 10
#else
#error TABOO requires quad precision floating point.
#error Please use a compiler that supports this.
#endif
 Real(qp), parameter :: pi    = 3.141592653589793238462643383279502884197_qp
 END MODULE STRATA 
!
!
!
!
!
 MODULE COMMON
!
! This module declares variables common to all of the three TABOO Tasks 
! Revised GS May 2010 for g95 implementation 
!
 USE STRATA 
 IMPLICIT NONE 
 SAVE 
!
 Real(dp), parameter :: raggio   = 6.371D6     ! Earth radius (m)
 Real(dp), parameter :: emass	 = 5.97D24     ! Earth mass (kg)
 Real(dp), parameter :: rhoea	 = 5.51157D3   ! Average density (kg/m3)
 Real(dp), parameter :: rhoice   = 1.D3        ! Ice density (kg/m3)
!
 Integer(i4b), parameter :: nv_max     =  24   ! maximum allowed number of v.e. layers 
 Integer(i4b), parameter :: nroots_max =  96   ! maximum allowed number of modes  
 Integer(i4b), parameter :: llmin      =   0   ! minimum allowed degree 
 Integer(i4b), parameter :: llmax      =1024   ! maximum allowed degree 
 Integer(sp)  :: vec (llmin:llmax,nroots_max)	! Marker for fake modes 
! 
 Integer(i4b) :: IV	        ! Verbose (1) or Silent (0) mode
 Integer(i4b) :: lmin, lmax     ! Current values of lmin e lmax, taken from TASK_*.DAT	 
 Integer(i4b) :: NV	        ! NV is the number of viscoelastic layers  
 Integer(i4b) :: NROOTS         ! NROOTS is the number of viscoelastic modes 
 Integer(i4b) :: Only_Elastic   ! Switch for elastic fields 
!
 Real(qp) :: r  (0:nv_max+1)    ! Radii 
 Real(qp) :: rho(0:nv_max+1)    ! Density 
 Real(qp) :: rmu(0:nv_max+1)    ! Shear Moduli 
 Real(qp) :: vis(0:nv_max+1)    ! Viscosity 
 Real(qp) :: pon(0:nv_max+1)    ! Gravity  
!
 Real(qp) :: h_e(llmin:llmax), l_e(llmin:llmax), k_e(llmin:llmax)  ! Elastic modes
 Real(qp) :: h_f(llmin:llmax), l_f(llmin:llmax), k_f(llmin:llmax)  ! Fluid modes 
 Real(qp) :: h_v (llmin:llmax,1:nroots_max) ! Viscoelastic residues h
 Real(qp) :: l_v (llmin:llmax,1:nroots_max) !	    "	       "    l 
 Real(qp) :: k_v (llmin:llmax,1:nroots_max) !	    "	       "    k 
 Real(qp) :: s (llmin:llmax,0:nroots_max)   ! Roots of the secular equation, in kyr**(-1)  
!
 Real(sp) :: r_h (llmin:llmax,1:nroots_max)   ! Normalized residue for h
 Real(sp) :: r_l (llmin:llmax,1:nroots_max)   ! Normalized residue for l
 Real(sp) :: r_k (llmin:llmax,1:nroots_max)   ! Normalized residue for k
 Real(sp) :: tek (llmin:llmax,1:nroots_max)   ! Relaxation time in k-yrs
!
 End module COMMON 
!
!
!
!
!
!
 MODULE COMMON_FOR_SPECTRUM
!
! Variables in use by Sbr. 'spectrum' and in the routines depending on it  
! Revised GS May 2010 for g95 implementation 
!
 USE COMMON
 USE STRATA 
 IMPLICIT NONE 
!
 Real (qp) :: pivo	 ! Largest coeff. of the secular polynomial 
 Real (qp) :: rubens	 ! For the elastic part of the solution 
 Real (qp) :: ggg	 ! Newton constant
 Real (qp) :: xmass	 ! Earth mass 
!
 Real (qp) :: aco (0:nv_max+1) ! <<A>> coefficients (bottom to top)  
 Real (qp) :: g (0:nv_max+1)   ! gravity at the interfaces (bottom to top) 
!
 Real (qp) :: a (6, 6), b (6, 6), c (6, 6), d (6, 6)  ! 6*6 propagators 
 Real (qp) :: ac(6, 6), ad(6, 6), bc(6, 6), bd(6, 6)  ! Outputs of MATPROD 
 Real (qp) :: k0 (6, 6), k1 (6, 6), k2 (6, 6)	      ! 6*6 propagators 
 Real (qp) :: matela (6, 6)			      ! Elastic product 
 Real (qp) :: sinist_2 (3, 6), sinist_1 (3, 6)        ! Left products 
 Real (qp) :: rr (0:2 * nv_max, 3, 3)		      ! Matrix to be inverted 
 Real (qp) :: co (0:3 * 2 * nv_max), aa(nroots_max + 1), op(nroots_max+1) ! Polynomial coefficints 
!
 Real (qp) :: rad (llmin:llmax,nroots_max)  ! Roots in years 
 Real (qp) :: cc (0:2, 1:nv_max, 6, 6)      ! Propagator 
 Real (qp) :: coefmat (0:2 * nv_max, 6, 6)  ! Propagator 
 Real (qp) :: ctmp (0:3 * 2 * nv_max - 1)   ! Derivative of the secular poly. 
 Real (qp) :: rrrr (0:2 * nv_max, 3, 6), qqqq (0:2 * nv_max, 3, 6)  ! Left products 
 Real (qp) :: cmb(6,3), bcs (3) 	    ! CMB and surface Boundary conditions 
 Real (qp) :: qq (0:2 * nv_max, 3, 3)	    ! Q matrix 
 Real (qp) :: r_r (3, 3, 0:nroots_max), q_q (3, 3, 0:nroots_max)  ! Ausilium di R e Q  
 Real (qp) :: qr (3, 3, 0:nroots_max)	    ! Product Q*R 
 Real (qp) :: raggiu (3, 3, 0:nroots_max)   ! Adjoint 
 Real (qp) :: derpo (1:nroots_max)	    ! Derivative of the secular poly. 
!
 Real (qp) :: rt1 (1:nroots_max), rt2 (1:nroots_max)  ! Real and Imag. parts of the roots, in kyears  
 Real (qp) ::  x_el (3), xx (3), xr (3, nroots_max)   ! Solution in vector form 
 Real (qp), external :: ded			      ! Useful to compute the determinant
!
 END MODULE COMMON_FOR_SPECTRUM 
!
!
!
!
!
!
 SUBROUTINE TABOO  
! 
! =============================================================================
!
! Derived by program TABOO  
! Revised GS May 2010 for g95 implementation 
!
 USE COMMON
 IMPLICIT NONE 
!
 Integer (i4b) :: oo(3) 	
 Character*20 date, TIMC
 Character*6  ST
!
! =============================================================================
!
 Open(99,file='taboo.log',status='unknown')
!
 Call DATE_AND_TIME (date,timc)      
 Write(99,*) '# ', date(1:4), '.',  date(5:6), '.',  date(7:8),  & 
         ' time=', timc(1:2), '.',  timc(3:4), '.',  timc(5:6) 
!
 Write(99,*)'Detecting the Active INPUT file !'
!
 OO(:)=0  
!
 OPEN(1,file='task_1.dat',status='unknown') 
 READ(1,'(A6)') ST
 IF( ST == 'Active' ) then
         OO(1) = 1
         Write(99,*)'task_1.dat is active '
                      else
         Write(99,*)'task_1.dat is NOT active '
 Endif
 CLOSE(1)
!
 If (oo(1)==1) CALL TASK_1 
!
 call DATE_AND_TIME (date,timc) 
 Write(99,*) '# Closing this file (taboo.log) on', date(1:4), '.', date(5:6), '.', date(7:8), & 
	                                 ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
!
 Close(99) 
!	
 END SUBROUTINE TABOO
!
!
!
!
!
 SUBROUTINE TASK_1         
!	
 USE COMMON
 IMPLICIT NONE
!
! Revised GS May 2010
!
 INTEGER (i4b) :: i_loading    ! loading/tidal switch
 INTEGER (i4b) :: ideja        ! checks for errors in the model KWs
 INTEGER (i4b) :: inone        ! checks for errors in the model KWs
 INTEGER (i4b) :: iarmo        ! detects the Harmonic_Degrees KW
 INTEGER (i4b) :: iexte        ! detects the kw External_Model  
 integer (i4b) :: CODE         ! Code to identify the models
 integer (i4b) :: ILM          ! (0/1) Controls the lower mantle stratification
 integer (i4b) :: ih, il, ik   ! Inputs for the (h/s) analysis
 integer (i4b) :: J, l, kk, i, k, m, ITASK      ! do-loops control variables
 integer (i4b) :: NUMD                          ! number of (harmonic) degree
 integer (i4b) :: NPUN                          ! number of time points
 integer (i4b), parameter :: max_num_righe = 20000  ! Max. n. of rows of file task_1.dat
 integer (i4b), parameter :: npun_max = 2001        ! Max number of time points
 integer (i4b), parameter :: numd_max = 6           ! Max. number of Harmonics
 integer (i4b) :: D(numd_max)                       ! Vector of harmonic degrees
 REAL(DP) :: T_MIN, T_MAX, P_MIN, P_MAX, ALFA, BETA, TIME   ! For the Heaviside TH
 REAL(DP) :: RESPH, RESPL, RESPK                            ! For the HEaviside TH
 real(sp) :: ajunk 
 Real(QP) :: LT 	 ! Lithospheric Thickness
 CHARACTER*30 KEYWORD    ! A keyword
 CHARACTER*20 date, timc ! date and time
!
! =============================================================================
!
! # These integers are used to monitor the keywords sequence .... 
!
 iarmo=0
 ideja=0
 inone=0
 iexte=0
!
 Write(99,*) 'Opening file  task_1.dat ...  '
! 
 open(1,file='task_1.dat',status='unknown')
!
 Write(99,*) 'Reading file task_1.dat ...  '
! 
!
 DO 101 ITASK = 1, max_num_righe
!
!
 Read(1,'(a30)',End=102) KEYWORD
!
!
!
!
 If(KEYWORD(1:16) == 'Harmonic_Degrees') Then
!
!
!
 iarmo=1
!
 Read(1,*)  LMIN, LMAX
 READ(1,*)  IV
!
 IF(IV/=0.and.IV/=1) THEN
 WRITE(99,*) 'ERROR in sbr task_1.dat:     The VERBOSE switch '
 WRITE(99,*) 'can only assume values 0 and 1. ** JOB ABORTED* '
 STOP
 endif
!
 call DATE_AND_TIME (date,timc)
 IF    (IV==0) THEN
   Write(99,*) '# Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 ELSEIF(IV==1) then
   Write(99,*) '# Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 ENDIF
!
 IF(iv==1) then
     WRITE(*,*) '= = = = = = = = = = = = = = ='
     Write(*,*) '   Reading file task_1.dat   '
     WRITE(*,*) '= = = = = = = = = = = = = = ='
 ENDIF
!
 READ(1,*)  i_loading
!
 IF((i_loading/=0.and.i_loading/=1).and.IV==0) THEN
 WRITE(99,*)'ERROR in sbr task_1.dat:	  The i_loading switch '
 WRITE(99,*)'can only assume values 0 and 1. **** JOB ABORTED* '
 STOP
 ENDIF
 IF((i_loading/=0.and.i_loading/=1).and.IV==1) THEN
 WRITE(99,*)'ERROR in sbr task_1.dat:	  The i_loading switch '
 WRITE(99,*)'can only assume values 0 and 1. **** JOB ABORTED* '
 WRITE(*, *)'ERROR in sbr task_1.dat:	  The i_loading switch '
 WRITE(*, *)'can only assume values 0 and 1. **** JOB ABORTED* '
 STOP
 ENDIF
!
!
 WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
 IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
 IF(IV==1 .AND. i_loading==1) &
		 WRITE(* ,*) 'Loading analysis'
 IF(IV==0 .AND. i_loading==1) &
		 WRITE(99,*) 'Loading analysis'
 IF(IV==1 .AND. i_loading==0) &
		 WRITE(* ,*) 'Tidal analysis'
 IF(IV==0 .AND. i_loading==0) &
                 WRITE(99,*) 'Tidal analysis'
                 Write(99,*) 'Lmin and Lmax ', Lmin, Lmax
       IF(IV==1) Write (*,*) 'Lmin and Lmax ', Lmin, Lmax
!
!
! # Checking Lmin and Lmax
!
 IF(IV==0) then
  If(lmin<llmin .or. lmax>llmax .or. lmax<lmin .or. lmin<0 .or. lmax<0) Then
  Write(99,*)'ERROR IN SBR. TASK_1:     One of the following forbidden conditions'
  Write(99,*)'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
  Write(99,*)'Check the input file task_2.dat ****** JOB ABORTED ****************'
  Stop
  Endif
 Endif
 IF(IV==1) then
  If(lmin<llmin .or. lmax>llmax .or. lmax<lmin .or. lmin<0 .or. lmax<0) Then
  Write(99,*)'ERROR IN SBR. TASK_1:     One of the following forbidden conditions'
  Write(99,*)'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
  Write(99,*)'Check the input file task_2.dat ****** JOB ABORTED ****************'
  Write(*,*) 'ERROR IN SBR. TASK_1      One of the following forbidden conditions'
  Write(*,*) 'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
  Write(*,*) 'Check the input file task_2.dat ****** JOB ABORTED ****************'
  Stop
  Endif
 Endif
!
 Endif  ! on KW Harmonic_Degrees
!
!
!
 IF(KEYWORD(1:10) == 'Make_Model') THEN   
!
!
!
 ideja=1
 inone=1
!
 IF(iarmo==0 .and. IV==0) then
  WRITE(99,*) ' ERROR IN SBR. TASK_1: the KW Harmonic_Degrees '
  WRITE(99,*) ' is not active. **** JOB ABORTED ************* '
  stop
 endif
 IF(iarmo==0 .and. IV==1) then
  WRITE(99,*) ' ERROR IN SBR. TASK_1: the KW Harmonic_Degrees '
  WRITE(99,*) ' is not active. **** JOB ABORTED ************* '
  WRITE(*,*)  ' ERROR IN SBR. TASK_1: the KW Harmonic_Degrees '
  WRITE(*,*)  ' is not active. **** JOB ABORTED ************* '
  STOP
 endif
!
              WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
 IF(IV==1)    WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
              Write(99,*)    'Building the model '
 IF(IV==1)    Write (*,*)    'Building the model '
!
 READ(1,*) NV 
 READ(1,*) CODE
 READ(1,*) ajunk     ! Some CYGWIN gfortrans have trouble reading REAL(10)
 LT = ajunk          ! DM Aug 1, 2011
 READ(1,*) ILM 
!      
 Write(99,*) 'Number of VE layers      ', NV  
 Write(99,*) 'Model CODE is            ', CODE 
!
 IF(IV==1) then
   Write(*,*) 'Number of VE layers      ', NV
   Write(*,*) 'Model CODE is            ', CODE
 ENDIF
!
 Write(99,'(a42,1x,f12.4,1x,a3)') &
 'From input, the lithospheric thickness is', LT,  ' km'
 Write(99,'(a47,1x,i3)')       &
 'The ILM parameter (for NV=7 and 9) is set to  ', ILM
 IF(IV==1) THEN
 Write(*, '(a42,1x,f12.4,1x,a3)') &
 'From input, the lithospheric thickness is', LT,  ' km'
 Write(*, '(a47,1x,i3)')       &
 'The ILM parameter (for NV=7 and 9) is set to  ', ILM
 ENDIF
!
 Write(99,*) 'Mantle viscosity from BOTTOM to TOP (/1E21) '
 IF(IV==1) WRITE(*,*)  'Mantle viscosity from BOTTOM to TOP (/1E21) '
 Do k=1, nv
      Read(1,*) ajunk    ! Some CYGWIN gfortrans have trouble reading REAL(10)
      vis(k) = ajunk     ! DM Aug 1 2011
      Write(99, '(a19,i4,1x,a1,F10.4)') &
           'Viscosity of layer', k, '=', vis (k)
 IF(IV==1) Write(*,  '(a19,i4,1x,a1,F10.4)') &
           'Viscosity of layer', k, '=', vis (k)
 vis(k) = vis(k) * 1e21_qp    	  
 Enddo 
! 
 CALL SPEC (Code, Ilm, LT)       
!
 call SPECTRUM (i_loading)
!
            Write(99,*)'Writing the spectrum on file spectrum.dat'
 IF(IV==1)  Write(*,* )'Writing the spectrum on file spectrum.dat'
!
  open(3,file='spectrum.dat',status='unknown')
  open(4,file='ss.dat',     status='unknown')
!
  call DATE_AND_TIME (date,timc)
!
  Write(3,*) 'File spectrum.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  
  Write(3,*) '1st column: l = Harmonic degree                '
  Write(3,*) '2nd   "   : LOG10 (l)                          '
  Write(3,*) '3th   "   : s, (kyrs**(-1))                    '
  Write(3,*) '4rd   "   : LOG10(-s), with s in kyrs**(-1)    '
  Write(3,*) '5th   "   : Relaxation time = -1000./s, (yrs)  '
  Write(3,*) '6th   "   : LOG10 (Relaxation time (yrs))      '
  Write(4,*) 'File spectrum.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  
  Write(4,*) '1st column: l = Harmonic degree                '
  Write(4,*) '2nd   "   : LOG10 (l)                          '
  Write(4,*) '3th   "   : s, (kyrs**(-1))                    '
  Write(4,*) '4rd   "   : LOG10(-s), with s in kyrs**(-1)    '
  Write(4,*) '5th   "   : Relaxation time = -1000./s, (yrs)  '
  Write(4,*) '6th   "   : LOG10 (Relaxation time (yrs))      '
!
  DO l=lmin, lmax 
!      
  Write(3,'(a1)') '>'
!
  do kk = 1, nroots  
!      					       
  write(3,'(i4,1x,5(e15.7,1x))') l, log10(float (l)) , s(l,kk), log10(-s(l,kk)), & 
      				        (-1000.0/s(l,kk)), log10(-1000.0/s(l,kk))      					           
  enddo 
!
  if(l<=10.or.mod(l,5)==0.or.mod(l,50)==0) then  
  Write(4,'(a1)') '>'
  do kk = 1, nroots  
  write (4, '(i4,1x,5(e15.7,1x))') l,                  & 
                                        log10 (float (l)), &
                                        s(l,kk),            &
					log10(-s(l,kk)),    &
				        (-1000.0/s(l,kk)), &
				        log10(-1000.0/s(l,kk))
 enddo
 endif
!
 ENDDO
 close(3)
 close(4) 
!
 ENDIF  !  Endif on KW Make_Model
!
!
!
 If(KEYWORD(1:19)=='Normalized_Residues') Then
!
!
!
 IF(inone==0.and.IV==0) then
           WRITE(99,*)'ERROR IN SBR. TASK_1:  One of the Kws Make_Model     and '
           WRITE(99,*)'External_Model MUST be active before Normalized_Residues '
           WRITE(99,*)'can be executed. ******* JOB ABORTED ******************* ';STOP
 endif
 IF(inone==0.and.IV==1) then
           WRITE(99,*)'ERROR IN SBR. TASK_1:  One of the Kws Make_Model     and '
           WRITE(99,*)'External_Model MUST be active before Normalized_Residues '
           WRITE(99,*)'can be executed. ******* JOB ABORTED ******************* '
           WRITE(*,*) 'ERROR IN SBR. TASK_1:  One of the Kws Make_Model     and '
           WRITE(*,*) 'External_Model MUST be active before Normalized_Residues '
           WRITE(*,*) 'can be executed. ******* JOB ABORTED ******************* ';STOP
 endif
!
 WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
 IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
 Write(99,*) 'Computing the normalized residues '
 IF(IV==1)   WRITE(*,*)  'Computing the normalized residues '
!
 Read(1,*) ih
 Read(1,*) il
 Read(1,*) ik 
!
 if(ik==1.and.iv==0.and.iexte==1) then 
        write(99,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
        write(99,*) 'WARNING in sbr TASK1: When the kw External_Model is active, '
	write(99,*) 'the k ldcs are zero by definition. I set IK==0 and continue '
        write(99,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
	IK=0
 endif 	 
 if(ik==1.and.iv==1.and.iexte==1) then 
        write(* ,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '      
        write(* ,*) 'WARNING in sbr TASK1: When the kw External_Model is active, '
	write(* ,*) 'the k ldcs are zero by definition. I set IK==0 and continue '
        write(* ,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
	IK=0
 endif 
!
 If(ih==0.and.il==0.and.ik==0.and.IV==0) Then
	Write(99,*) 'WARNING IN SBR. TASK_1:       No normalized residues '
	Write(99,*) 'will be printed: ih=il=ik =0 in input file task_1.dat'
 Endif
!
 If(ih==0.and.il==0.and.ik==0.and.IV==1) Then
	Write(* ,*) 'WARNING IN SBR. TASK_1:       No normalized residues '
	Write(* ,*) 'will be printed: ih=il=ik =0 in input file task_1.dat'
 Endif 
!
 IF(IH==1) THEN 
          open(3,file='ih.dat', status='unknown')
          open(4,file='ihh.dat',status='unknown')	  
!
            WRITE(99,*) 'Writing on ih.dat the h normalized residues'
  IF(IV==1) WRITE(* ,*) 'Writing on ih.dat the h normalized residues'
!
  call DATE_AND_TIME (date,timc)
  Write(3,*) 'File ih.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  Write(4,*) 'File ihh.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 	     
  If(i_loading==1) write(3,*) 'These h normalized residues are of loading type'
  If(i_loading==0) write(3,*) 'These h normalized residues are of _tidal_ type'
  If(i_loading==1) write(4,*) 'These h normalized residues are of loading type'
  If(i_loading==0) write(4,*) 'These h normalized residues are of _tidal_ type'  	     
!
 DO l=lmin, lmax 
 Write(3,'(a1)') '>'
!
 do kk = 1, nroots  
       write (3, '(i4,1x,e15.7)') l, - h_v(l,kk)/s(l,kk)
 enddo
!
 if(l.le.10.or.mod(l,5)==0.or.mod(l,50)==0) then 
 Write(4,'(a1)') '>'
 do kk = 1, nroots  
	    write (4, '(i4,1x,e15.7)') l,  abs(h_v(l,kk)/s(l,kk))
 enddo
 endif
 ENDDO
 close(3)
 close(4)   
!
 ENDIF 
!
 IF(IL==1) THEN
 	  open(3,file='il.dat',status='unknown')
          open(4,file='ill.dat',status='unknown')	  
!
            WRITE(99,*) 'Writing on il.dat the l normalized residues'
 IF(IV==1) WRITE(* ,*) 'Writing on il.dat the l normalized residues'
!
 call DATE_AND_TIME (date,timc)
 Write(3,*) 'File il.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(4,*) 'File ill.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 	     
 If(i_loading==1) write(3,*) 'These l normalized residues are of loading type'
 If(i_loading==0) write(3,*) 'These l normalized residues are of _tidal_ type'  	     
 If(i_loading==1) write(4,*) 'These l normalized residues are of loading type'
 If(i_loading==0) write(4,*) 'These l normalized residues are of _tidal_ type'
!	  	  
 DO l=lmin, lmax 
 Write(3,'(a1)') '>'
 do kk = 1, nroots  
            write (3, '(i4,1x,e15.7)') l, - l_v(l,kk)/s(l,kk)
 enddo
!
 if(l.le.10.or.mod(l,5)==0.or.mod(l,50)==0) then 
 Write(4,'(a1)') '>'
 do kk = 1, nroots  
 write (4, '(i4,1x,e15.7)') l, abs(l_v(l,kk)/s(l,kk))
 enddo
 endif
 ENDDO
 close(3)
 close(4)   
 ENDIF 
!
 IF(IK==1) THEN 
!
 open(3,file='ik.dat',status='unknown') 
 open(4,file='ikk.dat',status='unknown')	  
!
 WRITE(99,*) 'Writing on ik.dat the k normalized residues'
 IF(IV==1) WRITE(* ,*) 'Writing on ik.dat the k normalized residues'
!
 call DATE_AND_TIME (date,timc)
 Write(3,*) 'File ik.dat, created by Task#1 of TABOO on ', & 
		      date(1:4), '.', date(5:6), '.', date(7:8), & 
            ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(4,*) 'File ikk.dat, created by Task#1 of TABOO on ', & 
		      date(1:4), '.', date(5:6), '.', date(7:8), & 
            ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 If(i_loading==1) write(3,*) 'These k normalized residues are of loading type'
 If(i_loading==0) write(3,*) 'These k normalized residues are of _tidal_ type'  	    
!	  
!
 DO l=lmin, lmax 
 Write(3,'(a1)') '>'
!
 do kk = 1, nroots  
       write (3, '(i4,1x,e15.7)') l, - k_v(l,kk)/s(l,kk)
 enddo
!
 if(l.le.10.or.mod(l,5)==0.or.mod(l,50)==0) then 
 Write(4,'(a1)') '>'
 do kk = 1, nroots  
 write (4, '(i4,1x,e15.7)') l,  abs(k_v(l,kk)/s(l,kk))
 enddo
 endif
 ENDDO
 close(3)
 close(4)   
!
 ENDIF 
!
 Endif ! Endif on KW Normalized_Residues
!
! 
! 
 If(KEYWORD(1:15)=='El_Fluid_Viscel') Then
!
 WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
 IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!      
 read(1,*) ih
 read(1,*) il
 read(1,*) ik
!
 If(ih==0.and.il==0.and.ik==0.and.IV==0) Then
 Write(99,*) 'WARNING IN SBR. TASK_1: No El., Fluid, Viscel amplitudes'
 Write(99,*) 'will be printed:    ih=il=ik =0 in input file task_1.dat'
 Endif
 If(ih==0.and.il==0.and.ik==0.and.IV==1) Then
 Write(* ,*) 'WARNING IN SBR. TASK_1: No El., Fluid, Viscel amplitudes'
 Write(* ,*) 'will be printed:    ih=il=ik =0 in input file task_1.dat'
 Endif 
!
 if(ik==1.and.iv==0.and.iexte==1) then 
        write(99,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
        write(99,*) 'WARNING in sbr TASK1: When the kw External_Model is active, '
	write(99,*) 'the k ldcs are zero by definition. I set IK==0 and continue '
        write(99,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
	IK=0
 endif 	 
 if(ik==1.and.iv==1.and.iexte==1) then 
        write(* ,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '      
        write(* ,*) 'WARNING in sbr TASK1: When the kw External_Model is active, '
	write(* ,*) 'the k ldcs are zero by definition. I set IK==0 and continue '
        write(* ,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
	IK=0
 endif 
!
 IF(IH==1) THEN 
!
 open(3,file='h.dat', status='unknown')
 open(4,file='hh.dat',status='unknown')	  
!
 WRITE(99,*) 'Writing on h.dat the elastic, fluid, and v-elastic h'
 IF(IV==1) WRITE(* ,*) 'Writing on h.dat the elastic, fluid, and v-elastic h'
!
 call DATE_AND_TIME (date,timc)
 Write(3,*) 'File h.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(4,*) 'File h.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 If(i_loading==1) write(3,*) 'These h numbers are of loading type'
 If(i_loading==0) write(3,*) 'These h numbers are of _tidal_ type' 
 If(i_loading==1) write(4,*) 'These h numbers are of loading type'
 If(i_loading==0) write(4,*) 'These h numbers are of _tidal_ type'  
!	  
 DO l=lmin, lmax 
 write (3, '((i4,1x,96(1x,e20.8)))') &
 l, h_e(l), h_f(l), (h_v (l,k), k = 1, nroots)    
!          
 if(l.le.10.or.mod(l,5)==0.or.mod(l,50)==0) &
 write (4, '((i4,1x,96(1x,e20.8)))') &
 l, h_e(l), h_f(l), (h_v (l,k), k = 1, nroots) 	  
!	  
 ENDDO 
 ENDIF
 close(3)
 close(4)  
!
 IF(IL==1) THEN
!
 WRITE(99,*) 'Writing on l.dat the elastic, fluid, and v-elastic l'
 IF(IV==1) WRITE(* ,*) 'Writing on l.dat the elastic, fluid, and v-elastic l'
!
 open(3,file='l.dat',status='unknown')
 open(4,file='ll.dat',status='unknown')	  
!
 call DATE_AND_TIME (date,timc)
 Write(3,*) 'File l.dat, created by Task#1 of TABOO on ', & 
		      date(1:4), '.', date(5:6), '.', date(7:8), & 
            ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(4,*) 'File l.dat, created by Task#1 of TABOO on ', & 
		      date(1:4), '.', date(5:6), '.', date(7:8), & 
            ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 If(i_loading==1) write(3,*) 'These l numbers are of loading type'
 If(i_loading==0) write(3,*) 'These l numbers are of _tidal_ type'  
 If(i_loading==1) write(4,*) 'These l numbers are of loading type'
 If(i_loading==0) write(4,*) 'These l numbers are of _tidal_ type'		    
!	  
 DO l=lmin, lmax 
 write (3, '((i4,1x,96(1x,e20.8)))') &
 l, l_e(l), l_f(l), (l_v (l,k), k = 1, nroots)      
!
 if(l.le.10.or.mod(l,5)==0.or.mod(l,50)==0) & 
 write (4, '((i4,1x,96(1x,e20.8)))') &
 l, l_e(l), l_f(l), (l_v (l,k), k = 1, nroots)
!
 ENDDO 
 ENDIF
 close(3)
 close(4)   
!
 IF(IK==1) THEN
!
 WRITE(99,*) 'Writing on k.dat the elastic, fluid, and v-elastic k'
 IF(IV==1) WRITE(* ,*) 'Writing on k.dat the elastic, fluid, and v-elastic k'
!
 open(3,file='k.dat',status='unknown') 
 open(4,file='kk.dat',status='unknown')	  
!
!
 call DATE_AND_TIME (date,timc)
 Write(3,*) 'File k.dat, created by Task#1 of TABOO on ', & 
		      date(1:4), '.', date(5:6), '.', date(7:8), & 
            ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(4,*) 'File k.dat, created by Task#1 of TABOO on ', & 
		      date(1:4), '.', date(5:6), '.', date(7:8), & 
            ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 If(i_loading==1) write(3,*) 'These k numbers are of loading type'
 If(i_loading==0) write(3,*) 'These k numbers are of _tidal_ type'
 If(i_loading==1) write(4,*) 'These k numbers are of loading type'
 If(i_loading==0) write(4,*) 'These k numbers are of _tidal_ type'  
!	  
 DO l=lmin, lmax 
 write (3, '((i4,1x,96(1x,e20.8)))') &
 l, k_e(l), k_f(l), (k_v (l,k), k = 1, nroots)      
!
 if(l.le.10.or.mod(l,5)==0.or.mod(l,50)==0) & 
 write (4, '((i4,1x,96(1x,e20.8)))') &
 l, k_e(l), k_f(l), (k_v (l,k), k = 1, nroots) 
!
 ENDDO 
 ENDIF
 close(3)
 close(4)   
!
 Endif  ! Endif on KW El_Fluid_Viscel
!
!
101 CONTINUE     ! End of the do-loop on the input file task_1.dat
!
102 CONTINUE
!
 CLOSE(1) 
!
!
 call DATE_AND_TIME (date,timc)
 IF(IV==0) THEN
 Write(99,*) '# Task#1 of TABOO closed on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 ELSEIF(IV==1) then
     Write(99,*) '# Task#1 of TABOO closed on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 ENDIF
!
 IF(IV==1) Write(*,*) 'For more details, see file taboo.log '	
!
!
 END SUBROUTINE TASK_1
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
 SUBROUTINE SPECTRUM (I_LOADING)
!
 use COMMON_FOR_SPECTRUM
 implicit NONE
!
! The sbr SPECTRUM computes the relaxation spectrum, as well as
! the ldcs (or the tLns) of the selected Earth model.
! - Developed by GS and Carlo Giunchi in ~2003 
! - Reviewed by GS May 2010 for g95 implementation 
!
 Integer (i4b) :: nl	    ! Number of radii is always equal to nv+1
 INTEGER (i4b) :: i_loading ! loading/tidal switch (0/1)
 Integer (i4b) :: indx, inu, i, ii, j, k, kk, l, ll, m  
 Real(qp) pj1(3,6), pj2(3,6)
!
! ###########################################################################
!
 nl = nv+1 
!
!# Normalization
 call DEFPA  
!
!# Loop on the harmonic degrees
 do 1 l = lmin, lmax
!
! # Propagator
 do 2 k = 1, nl - 1  
 call matprod (l, r (nl - k), r (nl - k - 1), & 
               rho(nl - k), aco (nl - k), g (nl - k), g (nl - k - 1)) 
 k0 = ac + bd + bc * rmu (nl - k) + ad / rmu (nl - k)
 k1 = ad / vis (nl - k)  
 k2 = bc * ( - rmu (nl - k) **2 / vis (nl - k) )  
 cc (0, k, :, :) = k1 * rmu (nl - k) / vis (nl - k)  
 cc (1, k, :, :) = k0 * rmu (nl - k) / vis (nl - k) + k1 + k2
 cc (2, k, :, :) = k0  
 2 CONTINUE   
!
! # Elastic product
!
 call diretta (l, r (nl),     rho (nl), aco (nl), g (nl),     a, b)  
 call inversa (l, r (nl - 1), rho (nl), aco (nl), g (nl - 1), c, d)
 matela = matmul (a + rmu (nl) * b, c + d / rmu (nl) )  
!
!#Projectors "1" and "2" 
 pj1(:,:)=0.0_qp 
 pj2(:,:)=0.0_qp
!
 If(l==1) then 
  pj1(1,1)=1.0_qp
  pj1(2,2)=1.0_qp
  pj1(3,5)=1.0_qp
  pj2(1,3)=1.0_qp
  pj2(2,4)=1.0_qp
  pj2(3,5)=1.0_qp
 Else 
  pj1(1,1)=1.0_qp
  pj1(2,2)=1.0_qp
  pj1(3,5)=1.0_qp
  pj2(1,3)=1.0_qp
  pj2(2,4)=1.0_qp
  pj2(3,6)=1.0_qp
 Endif
!
 sinist_2 = matmul (pj2, matela); sinist_1 = matmul (pj1, matela)  
! 
!# Returns coefmat 
 call promat  
!
!# Left products
 Do inu = 0, 2 * nv  
     rrrr (inu, :, :) = matmul (sinist_2, coefmat (inu, :, :) )  
     qqqq (inu, :, :) = matmul (sinist_1, coefmat (inu, :, :) )  
 Enddo  
!
!# Right products
 Call corebo (l, r (0), rho (0), aco (0), cmb)  
 Do inu = 0, 2 * nv  
     rr (inu, :, :) = matmul (rrrr (inu, :, :), cmb)  
     qq (inu, :, :) = matmul (qqqq (inu, :, :), cmb)  
 Enddo  
!
!# Computes the determinant
    call det_tu
!
    do i = 1, nroots + 1  
    indx = 6 * nv + 1  
!
!# Coefficients of the secular polynomial
    aa (i) = co (indx - i)      
    enddo  
    pivo = aa (1) 
!    	
!# Roots of the secular polynomial: rt1 = REAL part, rt2 = IMAGINARY part
!  Finds the roots of the secular polynomial by the eigenvalue method, 
!  using "NETLIB" routines 
 Call hqr_driver (aa, nroots, rt1, rt2)  
!
!# VEC==0 for non physical modes, VEC==1 for 'good' ones 
 Do i = 1, nroots 
    if(abs(rt2(i)) >= 1.d-30 .or. rt1(i) >= 0._qp) then
      vec(l,i)= 0._sp  
      else 
      vec(l,i)= 1._sp  
      endif 
 Enddo  
!
 Do kk=1,nroots
 If(abs(rt2(kk)) >= 1.d-30) then 
  write(99,*) 'WARNING from sbr spectrum: A root with abs(imaginary part) <='
  WRITE(99,*) '1E-30 kyrs**(-1)  has been found  in the relaxation spectrum '  
  WRITE(99,*) 'degree = ', l, 'root # = ',kk
  WRITE(99,*) 'R.P.= ', rt1(kk)
  WRITE(99,*) 'I.P.= ', rt2(kk) 
  if(iv==1) Write(*,*) 'There are''non physical modes''Check file taboo.log'
 Endif
 If(rt1(kk) >= 0._qp) then
  write(99,*) 'WARNING from sbr spectrum: A root with real part >= 0'
  WRITE(99,*) 'has been found in the relaxation spectrum ********** '
  WRITE(99,*) 'degree = ', l, 'root # = ',kk
  WRITE(99,*) 'R.P.= ', rt1(kk)
  WRITE(99,*) 'I.P.= ', rt2(kk) 
  if(iv==1) Write(*,*) 'There are''non physical modes''Check file taboo.log'
 Endif  
 Enddo
!    
!# Pole in "-oo", "s" is in kyrs^(-1), "rad" is in years                            
#ifdef XLF
 s(l,0) = - 1.e8_qp   
#else
 s(l,0) = - 1.e10_qp  
#endif
 Do k = 1, nroots
   s(l,k) = 1.0_qp/ rt1 (k)           
   rad(l,k) = -1000.0_qp/s(l,k)    
 Enddo  
!
!# Arrays R and Q
 Do m = 0, nroots  
      r_r (:, :, m) = rr (2 * nv, :, :)  
      q_q (:, :, m) = qq (2 * nv, :, :)  	
      Do kk = 2 * nv - 1, 0, - 1  
            r_r (:, :, m) = r_r (:, :, m) * s (l,m) + rr (kk, :, :)  
            q_q (:, :, m) = q_q (:, :, m) * s (l,m) + qq (kk, :, :)  
      Enddo  	  	
 Enddo  
!
!# Adjoint of R
 Do m=0,nroots  
  raggiu (1, 1, m) =   ded (r_r (2, 2, m), r_r (2, 3, m), r_r (3, 2, m), r_r (3, 3, m) )
  raggiu (2, 1, m) = - ded (r_r (2, 1, m), r_r (2, 3, m), r_r (3, 1, m), r_r (3, 3, m) )
  raggiu (3, 1, m) =   ded (r_r (2, 1, m), r_r (2, 2, m), r_r (3, 1, m), r_r (3, 2, m) )
  raggiu (1, 2, m) = - ded (r_r (1, 2, m), r_r (1, 3, m), r_r (3, 2, m), r_r (3, 3, m) )
  raggiu (2, 2, m) =   ded (r_r (1, 1, m), r_r (1, 3, m), r_r (3, 1, m), r_r (3, 3, m) )
  raggiu (3, 2, m) = - ded (r_r (1, 1, m), r_r (1, 2, m), r_r (3, 1, m), r_r (3, 2, m) )
  raggiu (1, 3, m) =   ded (r_r (1, 2, m), r_r (1, 3, m), r_r (2, 2, m), r_r (2, 3, m) )
  raggiu (2, 3, m) = - ded (r_r (1, 1, m), r_r (1, 3, m), r_r (2, 1, m), r_r (2, 3, m) )
  raggiu (3, 3, m) =   ded (r_r (1, 1, m), r_r (1, 2, m), r_r (2, 1, m), r_r (2, 2, m) )
 Enddo  
!
!# Product QR
 Do m = 0, nroots  
      qr (:, :, m) = matmul (q_q (:, :, m), raggiu (:, :, m) )  
 Enddo  
!
!# Surface boundary conditions
!   i_loading =1 --> Loading
!   i_loading =0 --> Tidal
 Call bf (l, i_loading)
!
! # Secular polynomial at s-> -oo
  rubens = co (3 * 2 * nv)  
 Do k = 3 * 2 * nv - 1, 0, - 1  
       rubens = rubens * s (l,0) + co (k)  
 Enddo  
!
! # Elastic part of the solution
 x_el = matmul (qr (:, :, 0), bcs)  
 x_el = x_el / rubens  
!
! # Derivative of the secular polynomial
 Do ii = 0, 3 * 2 * nv - 1  
       Ctmp (ii) = co (ii + 1) * float (ii + 1)  
 Enddo  
!
 Do m = 1, nroots  
      derpo (m) = ctmp (3*2*nv-1)  
      Do k = 3 * 2 * nv - 2, 0, -1  
           derpo (m) = derpo (m) * s (l,m) + ctmp (k)  
      Enddo  
 Enddo  
!
! # Viscoelastic residues
 Do m = 1, nroots  
      xx = matmul(qr (:, :, m), bcs)  
      xr(1,m) = xx (1) / derpo (m)  
      xr(2,m) = xx (2) / derpo (m)  
      xr(3,m) = xx (3) / derpo (m)  
 Enddo  
!
!# Love number 'l'
 l_e(l) = x_el(2)*xmass/r(nl)  
 Do k=1, nroots 
     l_v(l,k) = xr(2,k)*xmass/r(nl)  
 Enddo 
 l_f (l) = l_e(l)  
 Do m = 1, nroots
     l_f(l) = l_f(l) - l_v(l,m)/s(l,m)  
 Enddo
!
!# Love number 'h'
 h_e(l) = x_el(1)*xmass/r(nl)  
 Do k=1, nroots 
     h_v(l,k) = xr(1,k)*xmass/r(nl)  
 Enddo 
 h_f(l) = h_e(l)  
 Do m = 1, nroots
     h_f(l) = h_f(l) - h_v(l,m)/s(l,m)  
 Enddo
!
!# Love number 'k'
 k_e (l) = - 1._qp - (xmass/r(nl)/g(nl))*x_el(3)  
 Do k=1, nroots 
     k_v (l,k)  = - (xmass/r(nl)/g(nl))*xr(3,k)  
 Enddo 
 k_f(l) = k_e(l) 
 Do m = 1, nroots
      k_f(l) = k_f(l) - k_v(l,m)/s(l,m)  
 Enddo
!
!# End of the loop on the harmonic degrees
 1 CONTINUE 
!
 END SUBROUTINE SPECTRUM 
!
!
!
subroutine defpa  
!
! # Normalizes the variables and constants
!
use  COMMON_FOR_SPECTRUM; implicit none
integer (i4b) ::  j, k 
!
real (qp) :: massa_in (nv+1)
real (qp) :: ra0, rhor, rig, t0 
real (qp), parameter :: haskell = 1.e21_qp   ! HASKELL unit 
!
!
!
    ggg =0.667e-10_qp    ! Newton constant
!
!
!
!
! * Earth mass
!
xmass = 4._qp * pi * rho (0) * r (0) **3 / 3._qp  
do k = 1, nv+1  
xmass = xmass + (4._qp / 3._qp) * pi * (r (k) **3 - r (k - 1) **3) &
 * rho (k)
enddo  
!
!
!
! * Surface gravity
!
g (nv+1) = ggg * xmass / r (nv+1) / r (nv+1)  
!
!
! * Gravity at the interfaces
!
g (0) = ggg * (4._qp / 3._qp) * pi * rho (0) * r (0)  
do k = 1, nv+1
massa_in (k) = g (0) * r (0) * r (0) / ggg  
do j = 1, k  
massa_in (k) = massa_in (k) + (4._qp / 3._qp) * pi * (r (j) **3 - &
 r (j - 1) **3) * rho (j)
enddo  
g (k) = massa_in (k) * ggg / r (k) / r (k)  
enddo  
!	
!
! * Real gravity at the interfaces
!
    pon(:) = g(:)  
!
!
 Write(99,*) '-----------------------------------'
 Write(99,*) ' Gravity at the interfaces (m/s/s) '
 Write(99,*) '       (from bottom to top)        '
 Write(99,*) '-----------------------------------'
 Do k=0, nv+1 
       Write(99,'(i4,2x,E14.8)') k, pon(k) 
 Enddo 
!
!
!
 Write(99,'(A5,A40,E14.8)') '>>>>>', 'Average density of the Earth model is: ', & 
 		          xmass/((4.0_qp/3.0_qp)*pi*r(nv+1)**3)  
!
!	 
! ************************
! * Normalization scheme *
! ************************
!
ra0 = r (nv+1)    ! reference radius
!
rhor = rho (1)    ! reference density
!
rig = rmu (nv+1)  ! reference shear modulus
!
t0 = 1000._qp * 365.25_qp * 24._qp * 3600._qp  ! seconds in 1 kyr
!
!
! * Normalized parameters
!
do k = 0, nv+1  
    r (k) =   r (k) / ra0
  rho (k) = rho (k) / rhor
  rmu (k) = rmu (k) / rig
  vis (k) = vis (k) / rig / t0
enddo  
!
! * Normalized ravity constant
!
ggg = ggg * rhor**2 * ra0**2 / rig  
!
!
! * Normalized gravity at the interfaces
!
g (0) = ggg * (4._qp / 3._qp) * pi * rho (0) * r (0)  
do k = 1, nv+1  
massa_in (k) = g (0) * r (0) * r (0) / ggg  
do j = 1, k  
massa_in (k) = massa_in (k) + (4._qp / 3._qp) * pi * (r (j) **3 - &
 r (j - 1) **3) * rho (j)
enddo  
g (k) = massa_in (k) * ggg / r (k) / r (k)  
enddo  
!
! * 'A' constants
!
do k = 0, nv+1  
aco (k) = (4._qp * pi / 3._qp) * rho (k) * ggg  
enddo  
!
! * Normalized mass of the Earth
!
xmass = 4._qp * pi * rho (0) * r (0) **3 / 3._qp  
do k = 1, nv+1  
xmass = xmass + (4._qp / 3._qp) * pi * (r (k) **3 - r (k - 1) **3) &
 * rho (k)
enddo  
!
end subroutine DEFPA 
!
!
!
!
!
!
subroutine COREBO (n, r, ro, a, g)
!
! # Boundary conditions at the CMB
!
use strata 
implicit none 
integer (i4b) :: n
real (qp) :: g (6, 3)  
real (qp) :: r, ro, a 
!
g = 0._qp  
g (1, 1) = - (r** (n - 1) ) / a  
g (1, 3) = 1._qp  
g (2, 2) = 1._qp  
g (3, 3) = ro * a * r  
g (5, 1) = r**n  
g (6, 1) = 2._qp * (float (n) - 1._qp) * r** (n - 1)  
g (6, 3) = 3._qp * a  
!
end subroutine COREBO 
!
!
!
!
!
!
subroutine det_tu
!
! # Secular determinant
!
use  COMMON_FOR_SPECTRUM
implicit none  
integer (i4b) :: inex, is, k, l, m, mu, nu
!
!
do inex = 0, 3 * 2 * nv  
co (inex) = 0._qp  
    do k = 0, 2 * nv  
        do l = 0, 2 * nv  
            do m = 0, 2 * nv  
                if (k + l + m == inex) then  
                   co (inex) = co (inex) + &
		   rr (k, 1, 1) * rr (l, 2, 2) * rr (m, 3, 3) - &
		   rr (k, 1, 1) * rr (l, 2, 3) * rr (m, 3, 2) - &
		   rr (k, 1, 2) * rr (l, 2, 1) * rr (m, 3, 3) + &
		   rr (k, 1, 2) * rr (l, 2, 3) * rr (m, 3, 1) + &
		   rr (k, 1, 3) * rr (l, 2, 1) * rr (m, 3, 2) - &
		   rr (k, 1, 3) * rr (l, 2, 2) * rr (m, 3, 1)
                 endif  
            enddo  
        enddo  
    enddo  
enddo
end subroutine det_tu
!
!
!
!
!
!
function ded (a, b, c, d)
!
! # Determinant od a 2x2 table
!
use strata
implicit none  
real (qp) :: a, b, c, d, ded  
ded = a * d-b * c  
end function ded
!
!
!
!
!
!
subroutine bf (l,char)  
!
! # Surface boundary conditions
!
use COMMON_FOR_SPECTRUM 
implicit none  
integer (i4b) :: char, l  
!
!
real (qp) :: gamma  
gamma = (2._qp * float (l) + 1._qp) / 4._qp / pi / r (nv+1) / r (nv+1)  
!
!
! ********** For the implementation of "degree one" - March 2010 **********
!
if (char == 1) then  
    bcs (1) = - g (nv+1) * gamma  
    bcs (2) = 0._qp  
	If(l.ne.1) bcs (3) = - 4._qp * pi * ggg * gamma  
	If(l.eq.1) bcs (3) = 0._qp  
endif  
!
! ********** For the implementation of "degree one" - March 2010 **********
!
if (char == 0) then
    bcs (1) = 0._qp  
    bcs (2) = 0._qp  
    bcs (3) = -4._qp * pi * ggg * gamma  
endif  
end subroutine bf
!
!
!
!
!
!
subroutine inversa (n, r, rho, za, gra, a, b)  
!
! # Inverse of the fundamental matrix
!
use strata  
implicit none
integer (i4b) :: i, j, n  
real (qp) :: gra, r, rho, za
real (qp) :: a (6, 6), b (6, 6)  
real (qp) :: v(24)
real (qp) :: h1, h2, h3, h4, h5
real (qp) :: rn
!
!
a(:,:) = 0._qp 
b(:,:) = 0._qp 
!
!
rn = float (n)  
!
h1 = rn + 1._qp  
h2 = 1._qp + 2._qp * rn  
h3 = 2._qp * rn - 1._qp  
h4 = 3._qp + 2._qp * rn  
h5 = rn - 1._qp  
!
v (1) = h1 / h2  
v (2) = - 2._qp * h1 * (rn + 2._qp) / h2  
v (3) = rn * h1 / 2._qp / h3 / h2  
v (4) = - (rn**2 + 3._qp * rn - 1._qp) * rn / h2 / h3  
v (5) = rn / h2  
v (6) = 2._qp * rn * h5 / h2  
v (7) = rn * h1 / 2._qp / h2 / h4  
v (8) = h1 * (rn**2 - rn - 3._qp) / h2 / h4  
v (9) = - 2._qp * rn * h1 * (2._qp + rn) / h2  
v (10) = - h5 * rn * h1 * h1 / h3 / h2  
v (11) = - 2._qp * h5 * rn * h1 / h2  
v (12) = - rn**2 * h1 * (2._qp + rn) / h2 / h4  
v (13) = - v (1)  
v (14) = - v (3)  
v (15) = - v (5)  
v (16) = - v (7)  
v (17) = - rn * h1 / h2  
v (18) = (2._qp - rn) * h1 * rn / 2._qp / h3 / h2  
v (19) = rn * h1 / h2  
v (20) = rn * h1 * (3._qp + rn) / 2._qp / h2 / h4  
v (21) = - v (1)  
v (22) = - v (3)  
v (23) = - v (5)  
v (24) = - v (7)  
!
a (1, 1) = r**( - n - 1) * v (2)  
a (2, 1) = - r**( - n + 1) * v (4)  
a (3, 1) = 3._qp * za * r**( - n + 1) / (2._qp * float (n) + 1._qp)
a (4, 1) = r**(n) * v (6)  
a (5, 1) = - r**(n + 2) * v (8)  
a (6, 1) = - 3._qp * za * r**(n + 2) / (2._qp * float (n) + 1._qp)  
!
a (1, 2) = - r**( - n - 1) * v (9)  
a (2, 2) = r**( - n + 1) * v (10)  
a (4, 2) = - r**(n) * v (11)  
a (5, 2) = r**(n + 2) * v (12)  
!
a (6, 5) = - r**(n + 1)  
a (3, 6) = - r**( - n + 1) / (2._qp * float (n) + 1._qp)  
a (6, 6) = r**(n + 2) / (2._qp * float (n) + 1._qp)  
!
b (1, 1) = r**( - n) * rho * gra * v (1)  
b (2, 1) = - r**( - n + 2) * rho * gra * v (3)  
b (4, 1) = r**(n + 1) * rho * gra * v (5)  
b (5, 1) = - r**(n + 3) * rho * gra * v (7)  
!
b (1, 3) = r**( - n) * v (13)  
b (2, 3) = - r**( - n + 2) * v (14)  
b (4, 3) = r**(n + 1) * v (15)  
b (5, 3) = - r**(n + 3) * v (16)  
!
b (1, 4) = - r**( - n) * v (17)  
b (2, 4) = r**( - n + 2) * v (18)  
b (4, 4) = - r**(n + 1) * v (19)  
b (5, 4) = r**(n + 3) * v (20)  
!
b (1, 5) = - rho * r**( - n) * v (21)  
b (2, 5) = rho * r**( - n + 2) * v (22)  
b (4, 5) = - rho * r**(n + 1) * v (23)  
b (5, 5) = rho * r**(n + 3) * v (24)  
!
end subroutine inversa
!
!
!
!
!
!
subroutine diretta (n, r, rho, za, gra, a, b)  
!
! # Direct matrix
!
use strata 
implicit none
integer (i4b) :: i, j, n  
real (qp) :: a (6, 6), b (6, 6)
real (qp) :: gra, r, rho, za
real (qp) :: rn
real (qp) :: a1, a2, a3, a4, b1, b2, b3, b4 
real (qp) :: c1, c2, c3, c4, d1, d2, d3, d4
!
!
a(:,:) = 0._qp 
b(:,:) = 0._qp 
!
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
a (1, 1) = a1 * r** (n + 1)  
a (1, 2) = a2 * r** (n - 1)  
a (1, 4) = a3 * r** ( - n)  
a (1, 5) = a4 * r** ( - n - 2)  
a (2, 1) = b1 * r** (n + 1)  
a (2, 2) = b2 * r** (n - 1)  
a (2, 4) = b3 * r** ( - n)  
a (2, 5) = b4 * r** ( - n - 2)  
a (3, 1) = a (1, 1) * rho * gra  
a (3, 2) = a (1, 2) * rho * gra  
a (3, 3) = - rho * r** (n)  
a (3, 4) = a (1, 4) * rho * gra  
a (3, 5) = a (1, 5) * rho * gra  
a (3, 6) = - rho * r** ( - n - 1)  
a (5, 3) = a (3, 3) / rho  
a (5, 6) = a (3, 6) / rho  
a (6, 1) = 3._qp * za * a (1, 1)  
a (6, 2) = 3._qp * za * a (1, 2)  
a (6, 3) = - (2._qp * float (n) + 1._qp) * r** (n - 1)  
a (6, 4) = 3._qp * za * a (1, 4)  
a (6, 5) = 3._qp * za * a (1, 5)  
!
b (3, 1) = c1 * r** (n)  
b (3, 2) = c2 * r** (n - 2)  
b (3, 4) = c3 * r** ( - n - 1)  
b (3, 5) = c4 * r** ( - n - 3)  
!
b (4, 1) = d1 * r** (n)  
b (4, 2) = d2 * r** (n - 2)  
b (4, 4) = d3 * r** ( - n - 1)  
b (4, 5) = d4 * r** ( - n - 3)  
!
end subroutine diretta
!
!
!
!
!
!
subroutine promat  
use COMMON_FOR_SPECTRUM 
implicit none
integer (i4b) :: i, j, k, l, n
!
! # Product od propagators
!
real (qp) :: ai (0:2 * nv, 6, 6), bi (0:2, 6, 6), ci (0:2 * nv, 6, 6)  
ai (0:2, :, :) = cc (:, 1, :, :)  
if (nv/=1) then  
    do i = 2, nv  
    n = 2 * i  
        do k = 0, n  
        ci (k, :, :) = 0._qp  
        bi (:, :, :) = cc (:, i, :, :)  
            do l = 0, n - 2  
                do j = 0, 2  
                    if (l + j == k) then  
                    ci (k, :, :) = ci (k, :, :) + &
                          matmul (ai (l, :, :), bi (j, :, :) )
                    endif  
                enddo  
            enddo  
        enddo  
    ai (0:n, :, :) = ci (0:n, :, :)  
    enddo  
    coefmat = ci  
else  
    coefmat (0:2, :, :) = ai (0:2, :, :)  
endif  
end subroutine promat
!
!
!
!
!
!
      subroutine hqr_nlib(nm,n,low,igh,h,wr,wi,ierr)
      use strata
      implicit NONE
!c
!c ------------------------------------------------------------------------
!c    Tested and ported to quadruple precision by GS on January 6 2008.  
!c    Original is from http://www.netlib.org/eispack/3090vf/single/hqr.f
!c ------------------------------------------------------------------------
!c
!C  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
!c
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      real(qp) h(nm,n),wr(n),wi(n), p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
      logical notlas
!c
!c     this subroutine is a translation of the algol procedure hqr,
!c     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!c     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!c
!c     this subroutine finds the eigenvalues of a real
!c     upper hessenberg matrix by the qr method.
!c
!c     on input
!c
!c        nm must be set to the row dimension of two-dimensional
!c          array parameters as declared in the calling program
!c          dimension statement.
!c
!c        n is the order of the matrix.
!c
!c        low and igh are integers determined by the balancing
!c          subroutine  balanc.  if  balanc  has not been used,
!c          set low=1, igh=n.
!c
!c        h contains the upper hessenberg matrix.  information about
!c          the transformations used in the reduction to hessenberg
!c          form by  elmhes  or  orthes, if performed, is stored
!c          in the remaining triangle under the hessenberg matrix.
!c
!c     on output
!c
!c        h has been destroyed.  therefore, it must be saved
!c          before calling  hqr  if subsequent calculation and
!c          back transformation of eigenvectors is to be performed.
!c
!c        wr and wi contain the real and imaginary parts,
!c          respectively, of the eigenvalues.  the eigenvalues
!c          are unordered except that complex conjugate pairs
!c          of values appear consecutively with the eigenvalue
!c          having the positive imaginary part first.  if an
!c          error exit is made, the eigenvalues should be correct
!c          for indices ierr+1,...,n.
!c
!c        ierr is set to
!c          zero       for normal return,
!c          j          if the limit of 30*n iterations is exhausted
!c                     while the j-th eigenvalue is being sought.
!c
!c     questions and comments should be directed to burton s. garbow,
!c     mathematics and computer science div, argonne national laboratory
!c
!c     this version dated september 1989.
!c
!c     ------------------------------------------------------------------
!c
      ierr = 0
      norm = 0._qp
      k = 1
!c     .......... store roots isolated by balanc
!c                and compute matrix norm ..........
      do 50 i = 1, n
!c
         do 40 j = k, n
   40    norm = norm + abs(h(i,j))
!c
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0._qp
   50 continue
!c
      en = igh
      t = 0._qp
      itn = 30*n
!c     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 1001
      its = 0
      na = en - 1
      enm2 = na - 1
!c     .......... look for single small sub-diagonal element
!c                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. 0._qp) s = norm
         tst1 = s
         tst2 = tst1 + abs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
!c     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
!c     .......... form exceptional shift ..........
      t = t + x
!c
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
!c
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75_qp * s
      y = x
      w = -0.4375_qp * s * s
  130 its = its + 1
      itn = itn - 1
!c     .......... look for two consecutive small
!c                sub-diagonal elements.
!c                for m=en-2 step -1 until l do -- ..........
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
!c
  150 mp2 = m + 2
!c
      do 160 i = mp2, en
         h(i,i-2) = 0._qp
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0._qp
  160 continue
!c     .......... double qr step involving rows l to en and
!c                columns m to en ..........
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
  170    s = sign(sqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
!c     .......... row modification ..........
         do 200 j = k, EN
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
!c
         j = min0(en,k+3)
!c     .......... column modification ..........
         do 210 i = L, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
         go to 255
  225    continue
!c     .......... row modification ..........
         do 230 j = k, EN
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
!c
         j = min0(en,k+3)
!c     .......... column modification ..........
         do 240 i = L, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
  255    continue
!c
  260 continue
!c
      go to 70
!c     .......... one root found ..........
  270 wr(en) = x + t
      wi(en) = 0._qp
      en = na
      go to 60
!c     .......... two roots found ..........
  280 p = (y - x) / 2._qp
      q = p * p + w
      zz = sqrt(abs(q))
      x = x + t
      if (q .lt. 0._qp) go to 320
!c     .......... real pair ..........
      zz = p + sign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0._qp) wr(en) = x - w / zz
      wi(na) = 0._qp
      wi(en) = 0._qp
      go to 330
!c     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
!c     .......... set error -- all eigenvalues have not
!c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end



      subroutine balanc_nlib(nm,n,a,low,igh,scale)
      use strata 
      implicit NONE
      integer i,j,k,l,m,n,nm,igh,low,iexc
      real(qp) a(nm,n),scale(n)
      real(qp) c,f,g,r,s,b2,radix
      logical noconv
!c
!c ------------------------------------------------------------------------
!c    Tested and ported to quadruple precision by GS on January 6 2008.  
!c    Original is from http://www.netlib.org/eispack/3090/double/balanc.f
!c ------------------------------------------------------------------------
!c
!c     this subroutine is a translation of the algol procedure balance,
!c     num. math. 13, 293-304(1969) by parlett and reinsch.
!c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!c
!c     this subroutine balances a real matrix and isolates
!c     eigenvalues whenever possible.
!c
!c     on input
!c
!c        nm must be set to the row dimension of two-dimensional
!c          array parameters as declared in the calling program
!c          dimension statement.
!c
!c        n is the order of the matrix.
!c
!c        a contains the input matrix to be balanced.
!c
!c     on output
!c
!c        a contains the balanced matrix.
!c
!c        low and igh are two integers such that a(i,j)
!c          is equal to zero if
!c           (1) i is greater than j and
!c           (2) j=1,...,low-1 or i=igh+1,...,n.
!c
!c        scale contains information determining the
!c           permutations and scaling factors used.
!c
!c     suppose that the principal submatrix in rows low through igh
!c     has been balanced, that p(j) denotes the index interchanged
!c     with j during the permutation step, and that the elements
!c     of the diagonal matrix used are denoted by d(i,j).  then
!c        scale(j) = p(j),    for j = 1,...,low-1
!c                 = d(j,j),      j = low,...,igh
!c                 = p(j)         j = igh+1,...,n.
!c     the order in which the interchanges are made is n to igh+1,
!c     then 1 to low-1.
!c
!c     note that 1 is returned for igh if igh is zero formally.
!c
!c     the algol procedure exc contained in balance appears in
!c     balanc  in line.  (note that the algol roles of identifiers
!c     k,l have been reversed.)
!c
!c     Questions and comments should be directed to Alan K. Cline,
!c     Pleasant Valley Software, 8603 Altus Cove, Austin, TX 78759.
!c     Electronic mail to cline@cs.utexas.edu.
!c
!c     this version dated january 1989. (for the IBM 3090vf)
!c
!c     ------------------------------------------------------------------
!c
      radix = 16._qp
!c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
!c     .......... in-line procedure for row and
!c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
!c
      do 30 i = 1, l
         f = a(i,j)
         a(i,j) = a(i,m)
         a(i,m) = f
   30 continue
!c
      do 40 i = k, n
         f = a(j,i)
         a(j,i) = a(m,i)
         a(m,i) = f
   40 continue
!c
   50 go to (80,130), iexc
!c     .......... search for rows isolating an eigenvalue
!c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
!c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 j = l, 1, -1
!c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (a(j,i) .ne. 0._qp) go to 120
  110    continue
!c
         m = l
         iexc = 1
         go to 20
  120 continue
!c
      go to 140
!c     .......... search for columns isolating an eigenvalue
!c                and push them left ..........
  130 k = k + 1
!c
  140 do 170 j = k, l
!c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (a(i,j) .ne. 0._qp) go to 170
  150    continue
!c
         m = k
         iexc = 2
         go to 20
  170 continue
!c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1._qp
!c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
!c
      do 270 i = k, l
         c = 0._qp
         r = 0._qp
!c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(a(j,i))
            r = r + abs(a(i,j))
  200    continue
!c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0._qp .or. r .eq. 0._qp) go to 270
         g = r / radix
         f = 1._qp
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
!c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95_qp * s) go to 270
         g = 1._qp / f
         scale(i) = scale(i) * f
         noconv = .true.
!c
         do 250 j = k, n
  250    a(i,j) = a(i,j) * g
!c
         do 260 j = 1, l
  260    a(j,i) = a(j,i) * f
!c
  270 continue
!c
      if (noconv) go to 190
!c
  280 low = k
      igh = l
      return
      end



!
!
!
!
subroutine hqr_driver (a, m, rtr, rti)
use strata
implicit none  
!
! -----------------------------------------------
! This is a driver for the root finding procedure 
! GS January 7, 2008. =Based on NETLIB routines=
! -----------------------------------------------
!
INTEGER, PARAMETER :: MAXM = 499   
INTEGER :: J, K, M, IGH, LOW, IERR  
REAL (QP) :: A (M + 1), RTR (M), RTI (M), & 
             H (MAXM, MAXM), SCALE(MAXM), & 
	     XR, XI   
!
!
If (m .gt. maxm .or. a(m+1) .eq. 0._qp) then 
Write(*,*) "Bad arguments in hqr_driver "
Write(*,*) "The program will STOP ***** "; Stop
Endif
!
! ..... Builds the Hessemberg matrix .....
!
do k = 1, m  
    h(1,k) = -a(m+1-k)/a(m+1)  
    do j = 2, m  
        h(j,k) = 0._qp  
    end do  
    if (k.ne.m) h (k+1,k) = 1._qp  
end do  
!
! ..... Computes the roots by the "Eigenvalues method" .....
!
 call balanc_nlib  (maxm, m, h, low, igh, scale)  
 call hqr_nlib     (maxm, m, low, igh, h, rtr, rti, ierr)
! In normal conditions ierr=0
! Write(*,*) ierr
!
! ..... Ascending sorting of real parts of  
!       roots using a code by John Burkardt ..... 
!
 call bubble_sort_qp(m, rtr)
!
return  
end
!
!
!
!
!
!
!*******************************************************************************
!
!! RVEC_SORT_BUBBLE_A ascending sorts a real vector using bubble sort.
!
!  ----------------------------------------------------
!  Adapted to quadruple precision by GS on January 2008
!  ----------------------------------------------------
!
!  Discussion:
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
!
  SUBROUTINE BUBBLE_SORT_QP(N, A)
  USE STRATA
  implicit none
  integer i, j, n
  real(qp) a(n)
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
!
!
!
!
!
!
!
!*******************************************************************************
!
!! R_SWAP swaps two real values.
!
!  ----------------------------------------------------
!  Adapted to quadruple precision by GS on January 2008
!  ----------------------------------------------------
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
!
  SUBROUTINE SWAP_QP(X,Y)
  USE STRATA
  implicit none
!
  real(qp) x, y, z 
!
  z = x
  x = y
  y = z
!
  return
end
!
!
!
!
!
!
subroutine matprod (n, rb, rd, rhu, za, grx, gry) 
!
! # An inverse times a direct fundamental matrix
!
use common_for_spectrum 
implicit none
integer (i4b) :: i, j, n  
real (qp) :: grx, gry, rb, rd, rhu, za
real (qp) :: v(24)
real (qp) :: h1, h2, h3, h4, h5
real (qp) :: rn
real (qp) :: a1, a2, a3, a4, b1, b2, b3, b4 
real (qp) :: c1, c2, c3, c4, d1, d2, d3, d4
real (qp) :: x,y
ac = 0._qp  
ad = 0._qp  
bc = 0._qp  
bd = 0._qp  
!
rn = float (n)  
x=rb
y=rd
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
h1 = rn + 1._qp  
h2 = 1._qp + 2._qp * rn  
h3 = 2._qp * rn - 1._qp  
h4 = 3._qp + 2._qp * rn  
h5 = rn - 1._qp  
!
v (1) = h1 / h2  
v (2) = - 2._qp * h1 * (rn + 2._qp) / h2  
v (3) = rn * h1 / 2._qp / h3 / h2  
v (4) = - (rn**2 + 3._qp * rn - 1._qp) * rn / h2 / h3  
v (5) = rn / h2  
v (6) = 2._qp * rn * h5 / h2  
v (7) = rn * h1 / 2._qp / h2 / h4  
v (8) = h1 * (rn**2 - rn - 3._qp) / h2 / h4  
v (9) = - 2._qp * rn * h1 * (2._qp + rn) / h2  
v (10) = - h5 * rn * h1 * h1 / h3 / h2  
v (11) = - 2._qp * h5 * rn * h1 / h2  
v (12) = - rn**2 * h1 * (2._qp + rn) / h2 / h4  
v (13) = - v (1)  
v (14) = - v (3)  
v (15) = - v (5)  
v (16) = - v (7)  
v (17) = - rn * h1 / h2  
v (18) = (2._qp - rn) * h1 * rn / 2._qp / h3 / h2  
v (19) = rn * h1 / h2  
v (20) = rn * h1 * (3._qp + rn) / 2._qp / h2 / h4  
v (21) = - v (1)  
v (22) = - v (3)  
v (23) = - v (5)  
v (24) = - v (7)  
!
ac (1,1) =  (a1*v(2)*(x/y)**(n+1) -a2*v(4)*(x/y)**(n-1) +a3*v(6)*(y/x)**(n) -a4*v(8)*(y/x)**(n+2))
ac (1,2) = -(a1*v(9)*(x/y)**(n+1)-a2*v(10)*(x/y)**(n-1)+a3*v(11)*(y/x)**(n)-a4*v(12)*(y/x)**(n+2))

ac (2,1) =  (b1*v(2)*(x/y)**(n+1) -b2*v(4)*(x/y)**(n-1) +b3*v(6)*(y/x)**(n) -b4*v(8)*(y/x)**(n+2))
ac (2,2) = -(b1*v(9)*(x/y)**(n+1)-b2*v(10)*(x/y)**(n-1)+b3*v(11)*(y/x)**(n)-b4*v(12)*(y/x)**(n+2))

ac (3,1) =  rhu*grx*ac(1,1)-(3._qp*rhu*za*y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))
ac (3,2) =  rhu*grx*ac(1,2)
ac (3,5) =  rhu*(y/x)**(n+1)
ac (3,6) =  (rhu*y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))

ac (5,1) = -(3._qp*za*y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))
ac (5,5) =  (y/x)**(n+1)
ac (5,6) =  (y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))

ac (6,1) =  3._qp*za*(a1*v(2)*(x/y)**(n+1)-(1+a2*v(4))*(x/y)**(n-1)+a3*v(6)*(y/x)**(n) -a4*v(8)*(y/x)**(n+2))
ac (6,2) =  3._qp*za*ac(1,2)
ac (6,6) =  (x/y)**(n-1)

ad (1,1) =  rhu*gry*( a1*v(1)*x*(x/y)**(n) -a2*v(3)*x*(x/y)**(n-2) +a3*v(5)*y*(y/x)**(n) -a4*v(7)*y*(y/x)**(n+2))
ad (1,3) =          (a1*v(13)*x*(x/y)**(n)-a2*v(14)*x*(x/y)**(n-2)+a3*v(15)*y*(y/x)**(n)-a4*v(16)*y*(y/x)**(n+2))
ad (1,4) =         -(a1*v(17)*x*(x/y)**(n)-a2*v(18)*x*(x/y)**(n-2)+a3*v(19)*y*(y/x)**(n)-a4*v(20)*y*(y/x)**(n+2))
ad (1,5) = -rhu*    (a1*v(21)*x*(x/y)**(n)-a2*v(22)*x*(x/y)**(n-2)+a3*v(23)*y*(y/x)**(n)-a4*v(24)*y*(y/x)**(n+2))

ad (2,1) =  rhu*gry*( b1*v(1)*x*(x/y)**(n) -b2*v(3)*x*(x/y)**(n-2) +b3*v(5)*y*(y/x)**(n) -b4*v(7)*y*(y/x)**(n+2))
ad (2,3) =          (b1*v(13)*x*(x/y)**(n)-b2*v(14)*x*(x/y)**(n-2)+b3*v(15)*y*(y/x)**(n)-b4*v(16)*y*(y/x)**(n+2))
ad (2,4) =         -(b1*v(17)*x*(x/y)**(n)-b2*v(18)*x*(x/y)**(n-2)+b3*v(19)*y*(y/x)**(n)-b4*v(20)*y*(y/x)**(n+2))
ad (2,5) = -rhu*    (b1*v(21)*x*(x/y)**(n)-b2*v(22)*x*(x/y)**(n-2)+b3*v(23)*y*(y/x)**(n)-b4*v(24)*y*(y/x)**(n+2))

ad (3,1) =  rhu*grx*ad(1,1)
ad (3,3) =  rhu*grx*ad(1,3)
ad (3,4) =  rhu*grx*ad(1,4)
ad (3,5) =  rhu*grx*ad(1,5)

ad (6,1) =  3._qp*za*ad(1,1)
ad (6,3) =  3._qp*za*ad(1,3)
ad (6,4) =  3._qp*za*ad(1,4)
ad (6,5) =  3._qp*za*ad(1,5)

bc (3,1) =  (c1*v(2)*(1._qp/y)*(x/y)**(n) -c2*v(4)*(1._qp/y)*(x/y)**(n-2) +c3*v(6)*(1._qp/x)*(y/x)**(n) &
            -c4*v(8)*(1._qp/x)*(y/x)**(n+2))
bc (3,2) = -(c1*v(9)*(1._qp/y)*(x/y)**(n)-c2*v(10)*(1._qp/y)*(x/y)**(n-2)+c3*v(11)*(1._qp/x)*(y/x)**(n) &
           -c4*v(12)*(1._qp/x)*(y/x)**(n+2))

bc (4,1) =  (d1*v(2)*(1._qp/y)*(x/y)**(n) -d2*v(4)*(1._qp/y)*(x/y)**(n-2) +d3*v(6)*(1._qp/x)*(y/x)**(n) &
           -d4*v(8)*(1._qp/x)*(y/x)**(n+2))
bc (4,2) = -(d1*v(9)*(1._qp/y)*(x/y)**(n)-d2*v(10)*(1._qp/y)*(x/y)**(n-2)+d3*v(11)*(1._qp/x)*(y/x)**(n) &
           -d4*v(12)*(1._qp/x)*(y/x)**(n+2))


bd (3,1) =  rhu*gry*(c1*v(1)*(x/y)**(n) -c2* v(3)*(x/y)**(n-2) +c3*v(5)*(y/x)**(n+1) -c4*v(7)*(y/x)**(n+3))
bd (3,3) =          (c1*v(13)*(x/y)**(n)-c2*v(14)*(x/y)**(n-2)+c3*v(15)*(y/x)**(n+1)-c4*v(16)*(y/x)**(n+3))
bd (3,4) =         -(c1*v(17)*(x/y)**(n)-c2*v(18)*(x/y)**(n-2)+c3*v(19)*(y/x)**(n+1)-c4*v(20)*(y/x)**(n+3))
bd (3,5) = -rhu*    (c1*v(21)*(x/y)**(n)-c2*v(22)*(x/y)**(n-2)+c3*v(23)*(y/x)**(n+1)-c4*v(24)*(y/x)**(n+3))

bd (4,1) = rhu*gry*(d1*v(1)*(x/y)**(n) -d2* v(3)*(x/y)**(n-2) +d3*v(5)*(y/x)**(n+1) -d4*v(7)*(y/x)**(n+3))
bd (4,3) =         (d1*v(13)*(x/y)**(n)-d2*v(14)*(x/y)**(n-2)+d3*v(15)*(y/x)**(n+1)-d4*v(16)*(y/x)**(n+3))
bd (4,4) =        -(d1*v(17)*(x/y)**(n)-d2*v(18)*(x/y)**(n-2)+d3*v(19)*(y/x)**(n+1)-d4*v(20)*(y/x)**(n+3))
bd (4,5) = -rhu*   (d1*v(21)*(x/y)**(n)-d2*v(22)*(x/y)**(n-2)+d3*v(23)*(y/x)**(n+1)-d4*v(24)*(y/x)**(n+3))

end subroutine matprod
!
!
!
!
!
!
!
!
!
!
!
!
!
!
subroutine prem (radsup, radinf, rigi, dens)  
!
! # This program gives the volume-averaged mean values of the
!   density,  and the lame parameters  mu and  lambda for the
!   prem-models for the reference periods at  1 s  and  200 s
!
!
! Input: 
!
!   radsup = bottom radius (km)
!   radinf = top radius  (km)
!
! Output: 
!
!   rigi = shear modulus of the layer (in units of 10**11 Pa)
!   dens = layer density (kg/m**3)
!
use strata
implicit none  
integer (i4b) :: i, idepth, itop, k, mod
real (qp) :: radsup, radinf, rigi, dens  
real (qp) :: mu (100), kappa (100), rho (100), r (100)  
real (qp) :: mean1 (0:100), mean2 (0:100), mean3 (0:100), depth, top
!
open (19, file = './DATA/'//'prem1.dat'  , status = 'old')  
open (29, file = './DATA/'//'prem200.dat', status = 'old')  
!
top = radsup  
depth = radinf  
!
!
124 format (4(e14.5))  
!
do mod = 29, 29  ! Uses the PREM at 200 s
!
    do i = 1, 94  
    read (mod, * ) r (i), kappa (i), mu (i), rho (i)  
    enddo  
    if (depth < 0._qp .or. depth > 6371._qp) then  
        print * , 'sorry, outside earth'  
        print *, top, depth  
        goto 13  
    endif  
    if (depth > top) then  
        print * , 'sorry, top is smaller than depth'  
        print *, top, depth  
        goto 13  
    endif  
    do i = 1, 94  
        if (r (i) - top >= 0._qp) then  
            itop = i  
            goto 4  
        endif  
    end do  
4   do i = 1, 94  
        if (r (i) - depth >= 0._qp) then  
            idepth = i - 1  
                if ( (r (i)  == r (i + 1) .and.r (i)  == depth) .or.r &
                   (i)  == depth) then
                   idepth = i  
                 endif  
             goto 6  
         endif  
     end do  
6     continue  
    if (mod == 19) then  
!         print*,'reference period of 1 s:'
!       print*,'------------------------'
    else  
!         print*,'reference period of 200 s:'
!       print*,'--------------------------'
    endif  
    if (r (idepth)  == r (itop) ) then  
!         print*,'homogeneous layering'
        mean1 (itop - 1) = mu (itop)  
        mean2 (itop - 1) = kappa (itop) - 2._qp / 3._qp * mu (itop)  
        mean3 (itop - 1) = rho (itop)  
        goto 11  
    endif  
!       print*,'reference top layer =',r(itop),' km'
!       print*,'reference lower layer =',r(idepth),' km'
    mean1 (idepth - 1) = 0._qp  
    mean2 (idepth - 1) = 0._qp  
    mean3 (idepth - 1) = 0._qp  
    do k = idepth, itop - 1  
        mean1 (k) = mean1 (k - 1) + (r (k + 1) **3 - r (k) **3) &
         / (r (itop) **3 - r (idepth) **3) * mu (k + 1)
        mean2 (k) = mean2 (k - 1) + (r (k + 1) **3 - r (k) **3) &
         / (r (itop) **3 - r (idepth) **3) * (kappa (k) - 2._qp / &
         3._qp * mu (k + 1) )
        mean3 (k) = mean3 (k - 1) + (r (k + 1) **3 - r (k) **3) &
         / (r (itop) **3 - r (idepth) **3) * rho (k + 1)
    end do  
!  11   print*,'mean value of mu = ',mean1(itop-1)/1000._qp,'x 10^11 pa'
!       print*,'mean value of lambda = ',mean2(itop-1)/1000._qp,'x 10^11 pa'
!       print*,'mean value of rho = ',mean3(itop-1)*1000._qp,'kg per m^3'
!
!
!
   11     rigi = mean1 (itop - 1) / 1.e3_qp  
    dens = mean3 (itop - 1) * 1.e3_qp  
!
!
!
end do  
close (19)  
close (29)  
13 continue  
end subroutine prem
!
!
!
!
!
!*************************************
!
  subroutine SPEC (CODE, ILM, LT)
!
!*************************************
!
use COMMON_FOR_SPECTRUM; implicit NONE 
!
!
! # La sbr SPEC helps to select the Earth model.
!
!   Input:   NV = Number of viscoelastic layers (from TASK_*.DAT)
!
!   Output:  NROOTS = Number of relaxation modes expected
!
!            R   (0:NV+1) = Radii of the interfaces
!            RHO (0:NV+1) = Density of the layers
!            RMU (0:NV+1) = Rigidity of the layers
!
!
!          where:
!
!            R(0)      = radius of the CMB
!            R(NV+1)   = radius of the Earth surface
!            RHO(0)    = density of the core
!            RHO(NV+1) = density of the lithosphere
!            RMU(0)    = rigidity of the core = 0.
!            RMU(NV+1) = rigidity of the lithosphere
!
!
!
!
!
!
! Features of the PREM model:
!
! o) We use the PREM 200s values
! o) IN PREM, the core radius is 3480 km
! o)  "   " , the earth radius is 6371 km
! o)  "   " , the "670" discontinuity is at radius R=5701 km
! o)  "   " , at the "670" discontinuity, d(rho)/rho is 12%'
! o)  "   " , there is a discontinuity at the depth 24.4 km (MOHO)
! o)  "   " , the "420" discontinuity is at radius R=5971 (it is in fact a '400')
! o)  "   " , there is a discontinuity at the depth of = 220 km (R= 6151 km)
! 
!
!
! 
!
!
real     (qp) :: LT      ! Lithospheric thickness from TASK_*.DAT (km)
real     (qp) :: rcmb    ! Radius of the CMB (fixed) 
real     (qp) :: rade    ! Radius of the Earth (fixed) 
real     (qp) :: JUNK    ! A junk 
real     (qp) :: DEL     ! Increment for the Lm 
real     (qp) :: TLM     ! Thickness of the Lower MAntle
!
parameter (rcmb = 3480._qp, rade = 6371._qp)  ! Fixed parameters 
!
integer (i4b) :: k       ! Do-loop index 
integer (i4b) :: CODE    ! A Code to label the various models, given NV
integer (i4b) :: NLM     ! Number of Lm layers 
integer (i4b) :: ILM     ! (0/1) Controls the Lm profile for NV= 7 and NV=9
INTEGER (i4b) :: idens   ! density inversions switch (0/1)
!
!
!
!
If(NV/=1.and.nv/=2.and.nv/=3.and.nv/=4.and.nv/=7.and.nv/=9.and.iv==0)Then
Write(99,*)'ERROR IN SBR. SPEC: The required model is not available '
WRITE(99,*)'in the models library. **** JOB ABORTED *************** ';STOP
Endif
If(NV/=1.and.nv/=2.and.nv/=3.and.nv/=4.and.nv/=7.and.nv/=9.and.iv==1)Then
Write (*,*)'ERROR IN SBR. SPEC: The required model is not available '
WRITE (*,*)'in the models library. **** JOB ABORTED *************** ';STOP
Endif
!
!
!
!
!
! #--------------------#
! #       NV = 1       #
! #--------------------#
!
  IF (NV == 1) THEN 
!
!                                  
! NV=1  CODE=0 ---> Averaged PREM, 30 <= LT <= 300 km.
! NV=1  CODE=1 ---> Yuen Sabadini Boschi [1982], LT= 100 km.
! NV=1  CODE=2 ---> Ondrej Cadek, LT= 70 km.
! NV=1  CODE=3 ---> Ondrej Cadek, LT= 200 km.
! NV=1  CODE=4 ---> Giunchi & Spada [2000] BUT self-g., 30 <= LT <= 300 km
! 
!
If ((CODE<0 .or. CODE>4) .and. iv==0) Then
Write(99,*) 'ERROR in SPEC: The model CODE is not available '
Write(99,*) '**** JOB ABORTED ******************************';STOP
Endif
If ((CODE<0 .or. CODE>4) .and. iv==1) Then
Write(*,*)  'ERROR in SPEC: The model CODE is not available '
Write(*,*)  '**** JOB ABORTED ******************************';STOP
Endif
!
!
!
!
!
if(CODE == 0) then    !Averaged PREM, 30 <= LT <= 300 km.                                        
!
!
nroots=4*nv
!
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC:  The LT parameter is out of bounds'
Write(99,*) 'It  must  be  30 <= LT <= 300 km  for NV=1 and CODE=0'
Write(99,*) '******** JOB ABORTED ********************************';STOP
Endif
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==1) Then
Write (*,*) 'ERROR in Sbr SPEC:  The LT parameter is out of bounds'
Write (*,*) 'It  must  be  30 <= LT <= 300 km  for NV=1 and CODE=0'
Write (*,*) '******** JOB ABORTED ********************************';STOP
Endif  
!
r (0)  = rcmb; r (1) = rade-LT; r (2) = rade 
!
 call prem (r (0), 0._qp, JUNK,      rho(0)  )    ! Core     
 call prem (r (2), r (1), rmu (2),   rho(2)  )    ! Litho
 call prem (r (1), r (0), rmu (1),   rho(1)  )    ! Mantle
 rmu(0) = 0._qp 
!
endif    						    
!
!+++++++++++++++++++++++
if(CODE == 1) then   !Yuen Sabadini Boschi [1982], LT= 100 km. 				            
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10927._qp; rmu(0)=0._qp       
r (1)      = 6271._qp; rho(1)=4314._qp ; rmu(1)=1.45_qp     
r (2)      = 6371._qp; rho(2)=2689._qp ; rmu(2)=0.28_qp    
!
r (0)      = 3481._qp; rho(0)=10925._qp; rmu(0)=0._qp       
r (1)      = 6251._qp; rho(1)=4500._qp ; rmu(1)=1.45_qp     
r (2)      = 6371._qp; rho(2)=4500._qp ; rmu(2)=1.45_qp   
!
endif                                                       
!+++++++++++++++++++++++
!
!+++++++++++++++++++++++
!if(CODE == 2) then   !Ondrej Cadek, LT= 70 km.  				            
!nroots=4*nv
!
!r (0)      = 3480._qp; rho(0)=11000._qp; rmu(0)=0._qp       
!r (1)      = 6251._qp; rho(1)=4500._qp ; rmu(1)=1.45_qp     
!r (2)      = 6371._qp; rho(2)=4500._qp ; rmu(2)=1.45_qp    
!
!endif

if(CODE == 2) then    !TEST with P WU, January 2006.   				            
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=11000._qp; rmu(0)=0._qp       
r (1)      = 6251._qp; rho(1)=4500._qp ; rmu(1)=1.45_qp     
r (2)      = 6371._qp; rho(2)=4500._qp ; rmu(2)=1.45_qp    
!
endif  

     
!+++++++++++++++++++++++
!
!+++++++++++++++++++++++
if(CODE == 3) then  !Ondrej Cadek, LT= 200 km.  				            
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10927._qp; rmu(0)=0._qp       
r (1)      = 6171._qp; rho(1)=4500._qp ; rmu(1)=1.45_qp     
r (2)      = 6371._qp; rho(2)=4500._qp ; rmu(2)=1.45_qp    
!
endif   
!+++++++++++++++++++++++
!
!
!
if(CODE == 4) then  !Giunchi & Spada [2000] BUT self-g., 30 <= LT <= 300 km
nroots=4*nv
!
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds '
Write(99,*) 'It  must  be  30 <= LT <= 300 km  for NV=1 and CODE=4'
Write(99,*) '**** JOB ABORTED ************************************';STOP
Endif
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds '
Write(*,*)  'It  must  be  30 <= LT <= 300 km  for NV=1 and CODE=4'
Write(*,*)  '**** JOB ABORTED ************************************';STOP
Endif  
!
r (0)      = 3480._qp     ;  rho(0)=10977._qp;  rmu(0)=0._qp
r (1)      = 6371._qp - LT;  rho(1)=4518._qp ;  rmu(1)=1.45_qp
r (2)      = 6371._qp     ;  rho(2)=3300._qp ;  rmu(2)=0.28_qp
!
endif
!
ENDIF    ! Endif on NV=1
!
!
!
! #--------------------#
! #       NV = 2       #
! #--------------------#
!
  IF (NV==2) THEN 
!
!
! NV=2  CODE=0 ---> Averaged PREM, 30 <= LT <= 300 km.
! NV=2  CODE=1 ---> Averaged PREM, BUT d(rho)/rho at 670 from PREM, 30 <= LT <= 300 km.
! NV=2  CODE=2 ---> Yuen Sabadini Boschi [1982], LT= 100 km.
! NV=2  CODE=3 ---> Bills & James [1997], LT= 100 km.
! NV=2  CODE=4 ---> Lefftz Sabadini Legros  [1994], LT= 150 km.
! NV=2  CODE=5 ---> Piersanti Postseismic, LT= 80 km.
! NV=2  CODE=6 ---> Ricard Sabadini Spada [1992], Model h, LT= 100 km.
!
!
If((CODE<0 .or. CODE>6) .and. iv==0) Then
Write(99,*)'ERROR in Sbr SPEC: The model CODE is not available '
Write(99,*)'**** JOB ABORTED **********************************';STOP
Endif
If((CODE<0 .or. CODE>6) .and. iv==1) Then
Write(*,*) 'ERROR in Sbr SPEC: The model CODE is not available '
Write(*,*) '**** JOB ABORTED **********************************';STOP
Endif
!
!
!
If(CODE == 0) then  !Averaged PREM, 30 <= LT <= 300 km.
!
          Write(99,*) 'Averaged PREM, 30 <= LT <= 300 km'
IF(iv==1) Write (*,*) 'Averaged PREM, 30 <= LT <= 300 km'
!
!
nroots=4*nv
!
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  30 <= LT <= 300 km  for NV=2  and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
Endif
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==1) Then
Write(*,*) 'ERROR in Sbr SPEC: The LT parameter is out of  bounds'
Write(*,*) 'It must be  30 <= LT <= 300 km  for NV=2  and CODE=0'
Write(*,*) '**** JOB ABORTED ************************************';STOP
Endif  
!
r(3)    = rade            
r(2)    = rade-LT         
r(1)    = r(3)-670._qp   
r(0)    = rcmb           
!
call prem(r(0), 0._qp, rmu(0), rho(0) )  ! Core 
call prem(r(1), r(0),  rmu(1), rho(1) )  ! Mant. Inf.
call prem(r(2), r(1),  rmu(2), rho(2) )  ! Mant. Sup. 
call prem(r(3), r(2),  rmu(3), rho(3) )  ! Litho
rmu (0) = 0._qp  
!
              endif
!
!
If(CODE == 1) then  !Averaged PREM, BUT d(rho)/rho at 670 from PREM, 30 <= LT <= 300 km.
!
          write(99,*) &
'Averaged PREM, BUT d(rho)/rho at 670 is from PREM, 30<=LT<=300 km'
IF(iv==1) WRITE(*,*)  &
'Averaged PREM, BUT d(rho)/rho at 670 is from PREM, 30<=LT<=300 km'
!
nroots=4*nv
!
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  30 <= LT <= 300 km  for NV=2  and CODE=1'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
Endif
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==1) Then
Write(*,*) 'ERROR in Sbr SPEC:  The LT parameter is out of bounds'
Write(*,*) 'It must be  30 <= LT <= 300 km  for NV=2  and CODE=1'
Write(*,*) '**** JOB ABORTED ************************************';STOP
Endif  
!
!
r(3)    = rade           
r(2)    = rade-LT         
r(1)    = r(3)-670._qp    
r(0)    = rcmb           
!
call prem(r(0), 0._qp, rmu(0), rho(0) )  ! Core 
call prem(r(1), r(0), rmu(1), rho(1) )  ! Mant. Inf.
call prem(r(2), r(1), rmu(2), rho(2) )  ! Mant. Sup. 
call prem(r(3), r(2), rmu(3), rho(3) )  ! Lito 
rmu (0) = 0._qp  
!
! 
rho(1) = (4.41241_qp+4.38071_qp)/2._qp*1000._qp
rho(2) = (3.99214_qp+3.98399_qp)/2._qp*1000._qp
!
             Endif 
!
!
!
If(CODE == 2) then !Yuen Sabadini Boschi [1982], LT= 100 km.
!
                 Write(99,*) 'Yuen Sabadini Boschi [1982], LT= 100 km'
IF(IV==1)        Write(* ,*) 'Yuen Sabadini Boschi [1982], LT= 100 km'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10927._qp; rmu(0)=0._qp       ! Core 
r (1)      = 5701._qp; rho(1)=4919._qp ; rmu(1)=2.17_qp     ! Lower mantle
r (2)      = 6271._qp; rho(2)=4430._qp ; rmu(2)=0.837_qp    ! Upper mantle
r (3)      = 6371._qp; rho(3)=2689._qp ; rmu(3)=0.282_qp    ! Litosphere  
!
Endif 
!
!
!
if(CODE == 3) then  !Bills & James [1997], LT= 100 km.
!
           Write(99,*) 'Bills & James [1997], LT= 100 km'
IF(iv==1)  Write(* ,*) 'Bills & James [1997], LT= 100 km'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10925._qp; rmu(0)=0._qp       ! Core 
r (1)      = 5701._qp; rho(1)=4508._qp ; rmu(1)=1.99_qp     ! Lower mantle
r (2)      = 6271._qp; rho(2)=4120._qp ; rmu(2)=0.954_qp    ! Upper mantle
r (3)      = 6371._qp; rho(3)=2771._qp ; rmu(3)=0.315_qp    ! Litosphere  
!
endif 
!
!
!
if(CODE == 4) then  !Lefftz Sabadini Legros  [1994], LT= 150 km.
!
!
           write(99,*) 'Lefftz Sabadini Legros  [1994], LT= 150 km'
IF(iv==1)  write(* ,*) 'Lefftz Sabadini Legros  [1994], LT= 150 km'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10987._qp; rmu(0)=0._qp        ! Core 
r (1)      = 5701._qp; rho(1)=4904._qp ; rmu(1)=2.225_qp     ! Lower mantle
r (2)      = 6221._qp; rho(2)=3666._qp ; rmu(2)=0.9169_qp    ! Upper mantle
r (3)      = 6371._qp; rho(3)=3232._qp ; rmu(3)=0.6114_qp    ! Litosphere  
!
endif 
!
!
If(CODE == 5) then  !Piersanti Postseismic, LT= 80 km.
!
          Write(99,*)'Postseismic Piersanti model, LT= 80 km '
IF(iv==1) WRITE (*,*)'Postseismic Piersanti model, LT= 80 km '

nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10932._qp; rmu(0)=0._qp        ! Core 
r (1)      = 5701._qp; rho(1)=4878._qp ; rmu(1)=2.171_qp     ! Lower mantle
r (2)      = 6291._qp; rho(2)=3614._qp ; rmu(2)=0.8464_qp    ! Upper mantle
r (3)      = 6371._qp; rho(3)=3115._qp ; rmu(3)=0.5597_qp    ! Litosphere  
!
endif 
!
if(CODE == 6) then !Ricard Sabadini Spada [1992], Model h, LT= 100 km.
!
          WRITE(99,*)  'Ricard Sabadini Spada [1992], Model h, LT= 100 km'
IF(iv==1) WRITE(*, *)  'Ricard Sabadini Spada [1992], Model h, LT= 100 km'
!
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10926._qp; rmu(0)= 0._qp      ! Core 
r (1)      = 5701._qp; rho(1)=4508._qp ; rmu(1)= 1.51_qp    ! Lower mantle
r (2)      = 6271._qp; rho(2)=4120._qp ; rmu(2)= 1.38_qp    ! Upper mantle
r (3)      = 6371._qp; rho(3)=4120._qp ; rmu(3)= 1.38_qp    ! Litosphere  
!
endif 
!
!
ENDIF   ! Endif on NV=2
!
!
!
!
! #--------------------#
! #       NV = 3       #
! #--------------------#
!
!
  IF (NV == 3) THEN 
!
!
! NV=3  CODE=0 ---> Averaged PREM, 30 <= LT <= 300 km.
! NV=3  CODE=1 ---> Not completely Averaged Prem, 30 <= LT <= 300 km.
! NV=3  CODE=2 ---> Cianetti Giunchi Spada [2002], LT=  120 km.
! NV=3  CODE=3 ---> James & Morgan [1990] GRL (table #1), LT= 120 Km.
! NV=3  CODE=4 ---> Similar to W. R. Peltier [1985], LT=  120 km.
! NV=3  CODE=5 ---> Paul Johnston [benchmark, 1997], LT= 70 km.
! NV=3  CODE=6 ---> Paulson, Zhong, and Wahr (2008) 
! NV=3, CODE=7 ---> 3 -layers discretization of peltier's VM2
!
!
If( (CODE<0 .or. CODE>7) .and. iv==0) Then
Write(99,*)'ERROR in Sbr SPEC: The CODE is not available'
Write(99,*)'**** JOB ABORTED ***************************';stop
Endif
If( (CODE<0 .or. CODE>7) .and. iv==1) Then
Write(*,*) 'ERROR in Sbr SPEC: The CODE is not available'
Write(*,*) '**** JOB ABORTED ***************************';Stop
Endif
!
!
!
!
if(CODE == 0) then  ! Averaged PREM, 30 <= LT <= 300 km.
!
           Write(99,*) 'Averaged PREM, 30 <= LT <= 300 km'
IF(iv==1)  Write(* ,*) 'Averaged PREM, 30 <= LT <= 300 km'
!
nroots=4*nv
!
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  30 <= LT <= 300 km   for NV=3 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
Endif
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  30 <= LT <= 300 km   for NV=3 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
Endif  
!
r(4)    = rade
r(3)    = rade-LT
r(2)    = rade-400._qp   ! Discontinuity at 400 km depth
r(1)    = rade-670._qp   ! Discontinuity at 670 km depth
r(0)    = rcmb           ! CMB radius
!
call prem(r(0), 0._qp, rmu(0), rho(0))  ! Core 
call prem(r(1), r(0), rmu(1), rho(1))  ! Lower mantle
call prem(r(2), r(1), rmu(2), rho(2))  ! Transition Zone
call prem(r(3), r(2), rmu(3), rho(3))  ! Upper mantle
call prem(r(4), r(3), rmu(4), rho(4))  ! Litho
rmu (0) = 0._qp  
!
endif 
!
!
if(CODE == 7) then  ! Averaged PREM, 30 <= LT <= 300 km.
!
           Write(99,*) 'Averaged PREM, 30 <= LT <= 300 km'
IF(iv==1)  Write(* ,*) 'Averaged PREM, 30 <= LT <= 300 km'
!
nroots=4*nv
!
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  30 <= LT <= 300 km   for NV=3 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
Endif
If((LT < 30._qp .OR. LT > 300._qp) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  30 <= LT <= 300 km   for NV=3 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
Endif  
!
r(4)    = rade
r(3)    = rade-LT
r(2)    = 5700._qp       ! Discontinuity at 400 km depth
r(1)    = 5000._qp       !  
r(0)    = rcmb           ! CMB radius
!
call prem(r(0), 0._qp, rmu(0), rho(0))  ! Core 
call prem(r(1), r(0),  rmu(1), rho(1))  ! Lower mantle
call prem(r(2), r(1),  rmu(2), rho(2))  ! Transition Zone
call prem(r(3), r(2),  rmu(3), rho(3))  ! Upper mantle
call prem(r(4), r(3),  rmu(4), rho(4))  ! Litho
rmu (0) = 0._qp  
!
endif 




!
!
!
If(CODE == 1) then  ! Not completely Averaged Prem, 30 <= LT <= 300 km.
!
Write(99,*) '3-layers mantle model NOT PREM averaged          '
Write(99,*) 'Lower mantle density = PREM value at depth =670- '
Write(99,*) 'Upper mantle density = PREM value at depth =400+ '
Write(99,*) 'Transition zone density = average PREM value     '
Write(99,*) 'The shear modulus is PREM-averaged               '
!
IF(iv==1) THEN
Write(*,*)  '3-layers mantle model NOT PREM averaged          '
Write(*,*)  'Lower mantle density = PREM value at depth =670- '
Write(*,*)  'Upper mantle density = PREM value at depth =400+ '
Write(*,*)  'Transition zone density = average PREM value     '
Write(*,*)  'The shear modulus is PREM-averaged               '
endif
!
nroots=4*nv
!
r(4)    = rade
r(3)    = rade-LT
r(2)    = rade-400._qp   ! Discontinuity at 400 km depth
r(1)    = rade-670._qp   ! Discontinuity at 670 km depth
r(0)    = rcmb           ! CMB radius
!
call prem(r(0), 0._qp, rmu(0), rho(0))  ! Core 
call prem(r(1), r(0),  rmu(1), rho(1))  ! LM 
call prem(r(2), r(1),  rmu(2), rho(2))  ! TZ 
call prem(r(3), r(2),  rmu(3), rho(3))  ! UM 
call prem(r(4), r(3),  rmu(4), rho(4))  ! Litho
rmu (0) = 0._qp   
!
!
rho(2) = (4.41241_qp+4.38071_qp)/2._qp*1000._qp
rho(3) = (3.54325_qp+3.51639_qp)/2._qp*1000._qp
!
endif 
!
!                      
!
If(CODE == 2) then ! Cianetti Giunchi Spada [2002], LT=  90 km. <==========
!
          Write(99,*) '3-- layer mantle model by CIANETTI et al. [2002]'
IF(iv==1) Write (*,*) '3-- layer mantle model by CIANETTI et al. [2002]'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10925._qp; rmu(0)=0._qp        ! Core 
r (1)      = 5701._qp; rho(1)=4508._qp ; rmu(1)=2.000_qp     ! Lower mantle
r (2)      = 5951._qp; rho(2)=4220._qp ; rmu(2)=1.100_qp     ! TZ  
r (3)      = 6281._qp; rho(3)=4120._qp ; rmu(3)=0.950_qp     ! Upper mantle    6281.
r (4)      = 6371._qp; rho(4)=4120._qp ; rmu(4)=0.730_qp     ! Litosphere  
!
endif 
!
!
!
!
!
if(CODE == 3) then ! James & Morgan [1990] GRL (table #1), LT= 120 Km.
!
          write(99,*) '3-- layer mantle model by James and Morgan (1990) (Tab. #1)  '
IF(iv==1) WRITE(*, *) '3-- layer mantle model by James and Morgan (1990) (Tab. #1)  '
nroots=4*nv
!
r (0)      = 3486._qp; rho(0)=11110._qp; rmu(0)=0._qp        ! Core 
r (1)      = 5701._qp; rho(1)=4900._qp ; rmu(1)=2.300_qp     ! Lower mantle
r (2)      = 5951._qp; rho(2)=3800._qp ; rmu(2)=1.450_qp     ! TZ  
r (3)      = 6251._qp; rho(3)=3550._qp ; rmu(3)=0.710_qp     ! Upper mantle
r (4)      = 6371._qp; rho(4)=2900._qp ; rmu(4)=0.400_qp     ! Litosphere  
!
endif 
!
!
!
if(CODE == 4) then ! Similar to W. R. Peltier [1985], LT=  120 km.
!
write(99,*) 'Similar to the 3-- layer mantle model by Peltier (1985) (Tab. #2)'
IF(iv==1) &
write (*,*) 'Similar to the 3-- layer mantle model by Peltier (1985) (Tab. #2)'
!
!
!
!
nroots=4*nv
!
r(4)    = rade           ! Radius of the Earth  
r(0)    = rcmb           ! Radius of the CMB 
rmu(0)  = 0._qp           
!
 call prem(r(0), 0._qp, rmu(0), rho(0))                        ! Core 
r (1) = rade-670._qp; rho(1)=4372._qp ; rmu(1)=rho(1)*6117._qp**2/1.e11_qp ! Lm
r (2) = rade-420._qp; rho(2)=4100._qp ; rmu(2)=rho(2)*5475._qp**2/1.e11_qp ! Tz
r (3) = rade-120._qp; rho(3)=3959._qp ; rmu(3)=rho(3)*5219._qp**2/1.e11_qp ! Um
 call prem(r(4), r(3), rmu(4), rho(4))                        ! Litho
!
endif 
!
!
!
if(CODE == 5) then   ! Paul Johnston [benchmark, 1997], LT= 70 km.
!
write(99,*) '3-- layer mantle model by Paul Johnston (benchmark, 1997)'
IF(iv==1) &
write( *,*) '3-- layer mantle model by Paul Johnston (benchmark, 1997)'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10750._qp; rmu(0)=0._qp          ! Core 
r (1)      = 5701._qp; rho(1)=4978._qp ; rmu(1)=2.2834_qp      ! Lower mantle
r (2)      = 5951._qp; rho(2)=3871._qp ; rmu(2)=1.0549_qp      ! TZ  
r (3)      = 6301._qp; rho(3)=3438._qp ; rmu(3)=0.70363_qp     ! Upper mantle
r (4)      = 6371._qp; rho(4)=3037._qp ; rmu(4)=0.50605_qp     ! Litosphere
!
endif 
!
!
if(CODE == 6) then   ! Paulson, Zhong, and Wahr (GJI, 2008) - See their Table 1 
!
write(99,*) '3-- layer mantle model by  Paulson, Zhong, and Wahr (GJI, 2008) - See their Table 1'
IF(iv==1) &
write( *,*) '3-- layer mantle model by  Paulson, Zhong, and Wahr (GJI, 2008) - See their Table 1'
!
nroots=4*nv
!
r (0)      = 3503.5_qp; rho(0)=9900._qp ; rmu(0)=0._qp        ! Core 
r (1)      = 5700._qp; rho(1)=4617._qp ; rmu(1)=2.254_qp     ! Lower mantle
r (2)      = 5960._qp; rho(2)=4227._qp ; rmu(2)=1.227_qp     ! TZ  
r (3)      = 6250._qp; rho(3)=4047._qp ; rmu(3)=0.759_qp     ! Upper mantle
r (4)      = 6370._qp; rho(4)=4047._qp ; rmu(4)=0.759_qp     ! Litosphere
!
endif 
!
!
ENDIF  !  Endif on NV = 3.
!
!
!
!
!
! #--------------------#
! #       NV = 4       #
! #--------------------#
!
!
! 
  IF (NV == 4) THEN 
!
!                                          
! NV=4  CODE=0 ---> Averaged PREM, 40 <= LT <= 150 km.
! NV=4  CODE=1 ---> Averaged PREM, LM not evaraged in Rho, 40 <= LT <= 150 km.
! NV=4  CODE=2 ---> Averaged, LM not averaged in Rho nor Mu, 40 <= LT <= 150 km.     
!
!
if((CODE <0 .or. CODE >3) .and. iv==0) Then
Write(99,*)'ERROR in SPEC: The model CODE is not available'
Write(99,*)'**** JOB ABORTED *****************************';STOP
Endif
if((CODE <0 .or. CODE >3) .and. iv==1) Then
Write (*,*)'ERROR in SPEC: The model CODE is not available'
Write (*,*)'**** JOB ABORTED *****************************';STOP
Endif 
!
!
!
If(CODE == 3) then 
!
!
Write(99,*)'4--layers PREM- averaged mantle for VM2'
!
IF(iv==1) THEN
Write(99,*)'4--layers PREM- averaged mantle for VM2'
ENDIF
!
nroots=4*nv
!
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=4 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=4 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
r(5)    = rade
r(4)    = rade-LT
r(3)    = 5650._qp 
r(2)    = 5050._qp 
r(1)    = 4600._qp 
r(0)    = rcmb 
!
!
 call prem(r(0), 0._qp, rmu(0), rho(0))  ! Core 
 call prem(r(1), r(0), rmu(1), rho(1))  ! Lower mantle
 call prem(r(2), r(1), rmu(2), rho(2))  ! Transition zone
 call prem(r(3), r(2), rmu(3), rho(3))  ! Shallow upper mantle (1)
 call prem(r(4), r(3), rmu(4), rho(4))  ! Shallow upper mantle (2)
 call prem(r(5), r(4), rmu(5), rho(5))  ! Lithosphere
rmu (0) = 0._qp  
!
endif 




!
!
!
If(CODE == 0) then 
!
!
Write(99,*)'4--layers PREM- averaged mantle model with:  '
Write(99,*)'Elastic lithosphere 40 <= LT <= 150 km       '
Write(99,*)'Shallow upper mantle 1  (70 <= Thick <=180)  '
Write(99,*)'Shallow upper mantle 2  (Thick =  180 km)    '
Write(99,*)'Transition zone         (Thick =  270 km)    '
Write(99,*)'Lower mantle down to the CMB                 '
!
IF(iv==1) THEN
Write(*,*) '4--layers PREM- averaged mantle model with:  '
Write(*,*) 'Elastic lithosphere 40 <= LT <= 150 km       '
Write(*,*) 'Shallow upper mantle 1  (70 <= Thick <=180)  '
Write(*,*) 'Shallow upper mantle 2  (Thick =  180 km)    '
Write(*,*) 'Transition zone         (Thick =  270 km)    '
Write(*,*) 'Lower mantle down to the CMB                 '
ENDIF
!
nroots=4*nv
!
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=4 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=4 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
r(5)    = rade
r(4)    = rade-LT
r(3)    = rade-220._qp   ! 200 km depth discontinuity
r(2)    = rade-400._qp   ! 400       "        "
r(1)    = rade-670._qp   ! 670       "        "
r(0)    = rcmb           ! CMB radius
!
!
 call prem(r(0), 0._qp, rmu(0), rho(0))  ! Core 
 call prem(r(1), r(0), rmu(1), rho(1))  ! Lower mantle
 call prem(r(2), r(1), rmu(2), rho(2))  ! Transition zone
 call prem(r(3), r(2), rmu(3), rho(3))  ! Shallow upper mantle (1)
 call prem(r(4), r(3), rmu(4), rho(4))  ! Shallow upper mantle (2)
 call prem(r(5), r(4), rmu(5), rho(5))  ! Lithosphere
rmu (0) = 0._qp  
!
endif 
!
!
If(CODE == 1) Then
!
!
Write(99,*) '4-- layers mantle PREM- averaged mantle model:   '
Write(99,*) 'Elastic lithosphere with 40 <= LT <= 150 km      '
Write(99,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)   '
Write(99,*) 'Shallow upper mantle 2  (Thick   =  180 km)      '
Write(99,*) 'Transition zone         (Thick   =  270 km)      '
Write(99,*) 'Lower mantle density = that of PREM at 670-      '
Write(99,*) 'Lower mantle rigidity is PREM- averaged dal PREM '
!
IF(iv==1) THEN
Write(*,*) '4-- layers mantle PREM- averaged mantle model:   '
Write(*,*) 'Elastic lithosphere with 40 <= LT <= 150 km      '
Write(*,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)   '
Write(*,*) 'Shallow upper mantle 2  (Thick   =  180 km)      '
Write(*,*) 'Transition zone         (Thick   =  270 km)      '
WRITE(*,*) 'Lower mantle density = that of PREM at 670-      '
Write(*,*) 'Lower mantle rigidity is PREM- averaged dal PREM '
ENDIF
!
nroots=4*nv
!
!
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=4 and CODE=1'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=4 and CODE=1'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
!
r(5)    = rade
r(4)    = rade-LT
r(3)    = rade-220._qp   !  Discontinuita' at 220 km depth.
r(2)    = rade-400._qp   !        "           400 km "   "
r(1)    = rade-670._qp   !        "           670 km "   "
r(0)    = rcmb           !  Radius of the CMB 
!
 call prem(r(0), 0._qp,  rmu(0), rho(0))  ! Core 
 call prem(r(1), r(0),  rmu(1), rho(1))  ! LM 
 call prem(r(2), r(1),  rmu(2), rho(2))  ! TZ 
 call prem(r(3), r(2),  rmu(3), rho(3))  ! Sm(1)
 call prem(r(4), r(3),  rmu(4), rho(4))  ! Sm(2) 
 call prem(r(5), r(4),  rmu(5), rho(5))  ! Lito  
!
rmu (0) = 0._qp  
rho (1) = (4.41241+4.38071)/2._qp*1000._qp   ! Lower mantle density
!
endif 
!
!
if(CODE == 2) Then
!
Write(99,*) '4-- layers mantle PREM- averaged mantle model:   '
Write(99,*) 'Elastic lithosphere with 40 <= LT <= 150 km      '
Write(99,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)   '
Write(99,*) 'Shallow upper mantle 2  (Thick   =  180 km)      '
Write(99,*) 'Transition zone         (Thick   =  270 km)      '
Write(99,*) 'Lower mantle density = that of PREM at 670-      '
Write(99,*) 'Lower mantle rigidity = that of PREM at 670-     '
!
IF(iv==1) then
Write (*,*) '4-- layers mantle PREM- averaged mantle model:   '
Write (*,*) 'Elastic lithosphere with 40 <= LT <= 150 km      '
Write (*,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)   '
Write (*,*) 'Shallow upper mantle 2  (Thick   =  180 km)      '
Write (*,*) 'Transition zone         (Thick   =  270 km)      '
Write (*,*) 'Lower mantle density = that of PREM at 670-      '
Write (*,*) 'Lower mantle rigidity = that of PREM at 670-     '
ENDIF
!
!
nroots=4*nv
!
!
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=4 and CODE=2'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=4 and CODE=2'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
r(5)    = rade
r(4)    = rade-LT
r(3)    = rade-220._qp   !  Discontinuita' at 220 km depth.
r(2)    = rade-400._qp   !        "           400 km "   "
r(1)    = rade-670._qp   !        "           670 km "   "
r(0)    = rcmb           !  Radius of the CMB 
!
call prem(r(0), 0._qp,  rmu(0), rho(0))  ! Core 
call prem(r(1), r(0),  rmu(1), rho(1))  ! LM 
call prem(r(2), r(1),  rmu(2), rho(2))  ! TZ 
call prem(r(3), r(2),  rmu(3), rho(3))  ! Sm(1)
call prem(r(4), r(3),  rmu(4), rho(4))  ! Sm(2) 
call prem(r(5), r(4),  rmu(5), rho(5))  ! Litho
!
rmu (0) = 0._qp  
rmu (1) = (1.584_qp + 1.639_qp)/2._qp   ! Lower mantle density
rho (1) = (4.41241+4.38071)/2._qp*1000._qp ! Lower mantle rigidity
!
endif 
!
!
ENDIF  ! Endif on NV=4
!
!
!
!
!
! #--------------------#
! #       NV = 7       #
! #--------------------#
!
!
!
  IF (NV == 7) THEN 
! 
!                                      
! NV=7  CODE=0 ---> Averaged PREM, 40 <= LT <= 150 km. 3 SM layers, 4 LM layers.
!
!
if((CODE <0 .or. CODE >0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The CODE is not available'
Write(99,*) '**** JOB ABORTED ***************************';STOP
Endif 
if((CODE <0 .or. CODE >0) .and. iv==1) Then
Write (*,*) 'ERROR in Sbr SPEC: The CODE is not available'
Write (*,*) '**** JOB ABORTED ***************************';STOP
Endif 
!
!
!
If(CODE == 0) then 
!
Write(99,*) '7-layers PREM averaged mantle model:  '
Write(99,*) 'Elastic lithosphere with thickness 40 <= LT <= 150 km  '
Write(99,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)         '
Write(99,*) 'Shallow upper mantle 2  (Thick   =  180 km)            '
Write(99,*) 'Transition zone         (Thick   =  270 km)            '
WRITE(99,*) 'Lower mantle is PREM averaged on 4 layers              '
!
IF(iv==1) THEN
Write(*,*) '7-layers PREM averaged mantle model:  '
Write(*,*) 'Elastic lithosphere with thickness 40 <= LT <= 150 km  '
Write(*,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)         '
Write(*,*) 'Shallow upper mantle 2  (Thick   =  180 km)            '
Write(*,*) 'Transition zone         (Thick   =  270 km)            '
WRITE(*,*) 'Lower mantle is PREM averaged on 4 layers              '
ENDIF
!
nroots=4*nv
!
!
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=7 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=7 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
r(8)    = rade                   ! Radius of the Earth 
r(7)    = rade-LT                ! Litho with thickness 40 <= LT <= 150 km   
r(6)    = rade-220._qp           ! Discontinuity at 220 km depth
r(5)    = rade-400._qp           !       "          400 km "   " 
r(4)    = rade-670._qp           !       "          670 km "   " 
r(0)    = rcmb 
!
!
NLM = 4            !  Number of v.e. layers in the Lower Mantle  (4) 
TLM = r(4)-r(0)    !  Thickness of the Lower Mantle 
!
If (ILM == 1) DEL = TLM/Float(NLM)  ! For ILM=1, the lower mantle layers
                                    ! have all the same thickness
!
If (ILM == 0) DEL = 370._qp          ! For ILM=0, the lower mantle layers
                                    ! have a thickness increasing with depth
!
Do K = 1, NLM-1 
     R(k) = R(k-1) + DEL + & 
     Float(NLM-k)*(TLM - NLM*DEL)/float(NLM)/float(NLM-1)/0.5_qp 
Enddo      
!
!
Write(99,*) 'Thickness of the lower mantle layers from BOTTOM to TOP (km):'
Write(99,'(f10.4)')   r(1)-r(0)
Write(99,'(f10.4)')   r(2)-r(1)
Write(99,'(f10.4)')   r(3)-r(2)
Write(99,'(f10.4)')   r(4)-r(3)
IF(iv==1) THEN
Write(*,*)  'Thickness of the lower mantle layers from BOTTOM to TOP (km):'
Write(*,'(f10.4)')   r(1)-r(0)
Write(*,'(f10.4)')   r(2)-r(1)
Write(*,'(f10.4)')   r(3)-r(2)
Write(*,'(f10.4)')   r(4)-r(3)
ENDIF
!
!
call prem(r(0), 0._qp, rmu(0), rho(0))  ! Core 
call prem(r(1), r(0), rmu(1), rho(1))  ! LM(1) 
call prem(r(2), r(1), rmu(2), rho(2))  ! LM(2) 
call prem(r(3), r(2), rmu(3), rho(3))  ! LM(3) 
call prem(r(4), r(3), rmu(4), rho(4))  ! LM(4) 
call prem(r(5), r(4), rmu(5), rho(5))  ! TZ 
call prem(r(6), r(5), rmu(6), rho(6))  ! Sm(1)
call prem(r(7), r(6), rmu(7), rho(7))  ! Sm(2) 
call prem(r(8), r(7), rmu(8), rho(8))  ! Lito  
rmu (0) = 0._qp  
!
endif 
!
ENDIF  ! Endif on NV=7
!
!
!
!
!
!
! #--------------------#
! #       NV = 9       #
! #--------------------#
!
!
  IF (NV == 9) THEN 
! 
!                                      
! CODE =0 --> Prem mediato (sm=3, lm=6)                     
!
!
if((CODE <0 .or. CODE >0) .and. iv==0) Then
Write(99,*)'ERROR in Sbr SPEC: The  CODE is not available'
Write(99,*)'**** JOB ABORTED ****************************';STOP
Endif
if((CODE <0 .or. CODE >0) .and. iv==1) Then
Write(*,*) 'ERROR in Sbr SPEC: The  CODE is not available'
Write(*,*) '**** JOB ABORTED ****************************';STOP
Endif 
!
!
!
!
!
If(CODE == 0) then 
!
Write(99,*) '9-layers PREM averaged mantle model:  '
Write(99,*) 'Elastic lithosphere with thickness 40 <= LT <= 150 km  '
Write(99,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)         '
Write(99,*) 'Shallow upper mantle 2  (Thick   =  180 km)            '
Write(99,*) 'Transition zone         (Thick   =  270 km)            '
WRITE(99,*) 'Lower mantle is PREM averaged on 4 layers              '
!
IF(iv==1) THEN
Write(*,*) '9-layers PREM averaged mantle model:  '
Write(*,*) 'Elastic lithosphere with thickness 40 <= LT <= 150 km  '
Write(*,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)         '
Write(*,*) 'Shallow upper mantle 2  (Thick   =  180 km)            '
Write(*,*) 'Transition zone         (Thick   =  270 km)            '
WRITE(*,*) 'Lower mantle is PREM averaged on 6 layers              '
ENDIF
!
nroots=4*nv
!
!
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=9 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40._qp .OR. LT > 150._qp) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=9 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
!
!
r(10)   = rade                   ! Radius of the Earth   
r(9)    = rade-LT                ! Litho with thickness 40 <= LT <= 150 km     
r(8)    = rade-220._qp           ! Discontinuity at 220 km depth
r(7)    = rade-400._qp           !       "          400 km "   " 
r(6)    = rade-670._qp           !       "          670 km "   " 
r(0)    = rcmb                   ! Radius of the CMB 
!
!
!
NLM = 6            !  Number of v.e. layers in the Lower Mantle  (6) 
TLM = r(6)-r(0)    !  Thickness of the Lower Mantle 
!
If (ILM == 1) DEL = TLM/Float(NLM)  ! For ILM=1, the lower mantle layers
                                    ! have all the same thickness
!
If (ILM == 0) DEL = 300._qp          ! For ILM=0, the lower mantle layers
                                    ! have a thickness increasing with depth
!
!
Do K = 1, NLM-1 
     R(k) = R(k-1) + DEL + & 
     Float(NLM-k)*(TLM - NLM*DEL)/float(NLM)/float(NLM-1)/0.5_qp 
Enddo      
!
!

Write(99,*) 'Lower mantle layers thickness, from BOTTOM to TOP (km):'
Write(99,'(f10.4)')   r(1)-r(0)
Write(99,'(f10.4)')   r(2)-r(1)
Write(99,'(f10.4)')   r(3)-r(2)
Write(99,'(f10.4)')   r(4)-r(3)
Write(99,'(f10.4)')   r(5)-r(4)
Write(99,'(f10.4)')   r(6)-r(5)
!
IF(iv==1) THEN
Write(99,*) 'Lower mantle layers thickness, from BOTTOM to TOP (km):'
Write(99,'(f10.4)')   r(1)-r(0)
Write(99,'(f10.4)')   r(2)-r(1)
Write(99,'(f10.4)')   r(3)-r(2)
Write(99,'(f10.4)')   r(4)-r(3)
Write(99,'(f10.4)')   r(5)-r(4)
Write(99,'(f10.4)')   r(6)-r(5)
endif

!
!
call prem(r(0),  0._qp, rmu(0),  rho(0))   ! Core 
call prem(r(1),  r(0), rmu(1),  rho(1))   ! LM 
call prem(r(2),  r(1), rmu(2),  rho(2))   ! LM 
call prem(r(3),  r(2), rmu(3),  rho(3))   ! LM 
call prem(r(4),  r(3), rmu(4),  rho(4))   ! LM 
call prem(r(5),  r(4), rmu(5),  rho(5))   ! LM 
call prem(r(6),  r(5), rmu(6),  rho(6))   ! LM 
call prem(r(7),  r(6), rmu(7),  rho(7))   ! TZ 
call prem(r(8),  r(7), rmu(8),  rho(8))   ! Sm(1)
call prem(r(9),  r(8), rmu(9),  rho(9))   ! Sm(2) 
call prem(r(10), r(9), rmu(10), rho(10))  ! Lito  
rmu (0) = 0._qp  
endif 
!
!
ENDIF  !  Endif on  NV=9
!
!
!
!
!
!
!
!
!
!
!
!
!
!  +-----------------------------------------+
!  |  Conversion of R, RMU & RHO in SI units |
!  +-----------------------------------------+
!
    r  (:) =r  (:)*1.e3_qp     ! R in m
    rmu(:) =rmu(:)*1.e11_qp    ! RMU in Pa 
    rho(:) =rho(:)             ! Rho in Kg/m**3
!
!
!
Write(99,*) '---------------------------------------------------------------'
Write(99,*) ' Radii, densities, shear moduli & viscosity from bottom to top '
Write(99,*) '        SI units: m, kg/m**3, Pa, Pa.s, respectively           '
Write(99,*) '---------------------------------------------------------------'
!
Do k = 0, nv+1   
  Write(99,'(i4, 1x, 5(2x, E14.8))') k, r(k), rho(k), rmu(k), vis(k) 
Enddo  
!
!
!
if(lt .ne. (r(nv+1)-r(nv))/1000._qp .and. iv==0) then
!write(99,*) lt, (r(nv+1)-r(nv))/1000._qp  
write(99,*) ' + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  '
write(99,*) 'WARNING from sbr SPEC:   The LT value given in input is not ' 
write(99,*) 'consistent with the default for the current values of NV &  '
write(99,*) 'and CODE.  The litho thickness is set to the default value. '
write(99,*) ' + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  '
endif 
if(lt .ne. (r(nv+1)-r(nv))/1000._qp .and. iv==1) then 
!write(*,*) lt, (r(nv+1)-r(nv))/1000._qp  
write(* ,*) ' + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  '
write(* ,*) 'WARNING from sbr SPEC:   The LT value given in input is not ' 
write(* ,*) 'consistent with the default for the current values of NV &  '
write(* ,*) 'and CODE.  The litho thickness is set to the default value. '
write(* ,*) ' + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  '
endif 
!
          Write(99,'(a51,F12.4)') &
'Lithospheric thickness effectively employed (km) =', (r(nv+1)-r(nv))/1000._qp 
if(iv==1) Write(* ,'(a51,F12.4)') &
'Lithospheric thickness effectively employed (km) =', (r(nv+1)-r(nv))/1000._qp 
!
!
!
!  +-------------------------------+
!  | Test for density inversions   |
!  +-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=+
!
WRITE(99,*) 'Looking for density inversions...'
!
idens=0
Do k = 1, nv+1  
if (rho (k) - rho (k - 1)  >  0._qp)  then
idens=1
Write ( 99, * ) 'WARNING from Sbr. SPEC: a density inversion exists !'
endif
Enddo
!
IF(idens==0)             WRITE(99,*) 'No density inversions found'
IF(idens==0 .and. iv==1) WRITE(* ,*) 'No density inversions found'
!
END SUBROUTINE SPEC 
!
!
!
!
