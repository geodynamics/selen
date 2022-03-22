!
!-*-*-*-*-*-*-*-*-*-*-*-*-*  
!
  PROGRAM SELEN4 
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*  
!
! REVISION 22 - Mar 18, 2022
!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! This program upgrades previous version (8) by the 
! revised theory in preparation of the paper on GMD 
! 
! GS March 14, 2019: File created and first tested
! GS March 29, 2019: Compliance with the theory 
! GS March 29, 2019: Reproduces the version 8 results
! GS April 04, 2019: Slight changes ... 
! GS April 05, 2019: Slight changes ... 
! GS April 06, 2019: Slight changes and new outputs 
! DM April 14, 2019: Slight fixes ...
! DM May 1, 2019:    Bug fixes ...
! DM May 3, 2019:    Configuration from command-line cfg file
! DM Jun 19, 2019:   Removed the reference to the S Green Function
! DM Jun 27, 2019:   Changes in the post-processor
! DM Jul 11, 2019:   Changes in the post-processor
! DM Jul 26, 2019:   Changes in the post-processor
! GS Apr 19, 2020:   Implementation of "horizontals"
! DM Sep 21, 2020:   Future GIA
! DM Mar 07, 2021:   New format for SRFs and RRFs
! DM Mar 18, 2022:   Checks if the config file exists
! DM Mar 18, 2022:   Rearranged loops for optimal performance 
! DM Mar 19, 2022:   Further optimizations
!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!
 IMPLICIT NONE 
!
 INCLUDE "parameters.inc"	
!		    
!
!<<<<<<<<<<<<<<<<<<<  CONFIGURATION PARAMETERS  >>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<< (originally in 'data.sle') >>>>>>>>>>>>>>>>>>>
!
 CHARACTER*100 :: CFG_F
 CHARACTER*200 :: BUFFER
!
!
! >>>>>>> Basic parameters --------------------------------
!
 INTEGER :: RES    
 INTEGER :: LMAX
 INTEGER :: NN
 INTEGER :: NV
 REAL*8  :: DELTA
 INTEGER :: IROT
 INTEGER :: IEXT_MAX
 INTEGER :: IINT_MAX
! 
! >>>>>>> Input Files required for the execution ------------
! 
 CHARACTER*100 :: PT 
! 
 CHARACTER*100 :: PIX_F 
 CHARACTER*100 :: SHA_F
!CHARACTER*100 :: SGF_F
 CHARACTER*100 :: UGF_F
 CHARACTER*100 :: VGF_F 
 CHARACTER*100 :: GGF_F
 CHARACTER*100 :: ROT_F
 CHARACTER*100 :: ICE_F
 CHARACTER*100 :: TOP_F
 CHARACTER*100 :: OFP_F
 CHARACTER*100 :: TGS_F
 CHARACTER*100 :: GEO_F
 CHARACTER*100 :: RSL_F
 CHARACTER*100 :: PMT_F
!
! >>>>>>> Output files --------------------------------------
! 
 CHARACTER*100 :: F_OUT_TOP 
 CHARACTER*100 :: F_OUT_OFU 
 CHARACTER*100 :: F_OUT_SHS 
 CHARACTER*100 :: F_OUT_SHU
 CHARACTER*100 :: F_OUT_SHV 
 CHARACTER*100 :: F_OUT_SHG 
 CHARACTER*100 :: F_OUT_SHN 
 CHARACTER*100 :: F_OUT_TPW 
 CHARACTER*100 :: F_OUT_SAV 
!
!
!<<<<<<<<<<<<<<<<  END OF CONFIGURATION PARAMETERS  >>>>>>>>>>>>>>>>
!
!
!
!------------------------------------------------------------------------
! Declarations... Declarations... Declarations... Declarations... 
!------------------------------------------------------------------------
!
 INTEGER :: NP, JMAX, NMODES
!
 INTEGER III, DOM, IINT, IEXT, J_INDEX
 INTEGER I, J, K, L, M, N, P, NA, LJ, MJ
 INTEGER :: KREF, NN0
 LOGICAL :: IOFLAG
!
 INTEGER, ALLOCATABLE :: OF(:,:), OFP(:), MM(:), LL(:), DM(:)
!
 COMPLEX*16, ALLOCATABLE :: SH_RAP(:,:), SH_RBP(:,:), SH_RCP(:,:)
 COMPLEX*16, ALLOCATABLE :: SH_KA(:,:), SH_KB(:,:), SH_KC(:,:)
 COMPLEX*16, ALLOCATABLE :: SH_XA(:,:), SH_XB(:,:), SH_XC(:,:) 
 COMPLEX*16, ALLOCATABLE :: SIN_COS(:,:), SH_U(:,:), SH_W(:,:)      
 COMPLEX*16, ALLOCATABLE :: SH_X(:,:), SH_G(:,:), SH_OFU(:,:) 
 COMPLEX*16, ALLOCATABLE :: SH_Z(:,:,:), SH_S(:,:,:), SH_ZAV(:,:)
 COMPLEX*16, ALLOCATABLE :: SH_SAV(:,:), SH_LL21(:), PSI_RIG(:)         
 COMPLEX*16, ALLOCATABLE :: SH_GR(:,:), SH_UR(:,:), SH_RR(:,:) 
 COMPLEX*16, ALLOCATABLE :: SH_RRP(:,:), SH_KR(:,:), SH_N(:,:)
 COMPLEX*16, ALLOCATABLE :: SH_LA21(:), SH_LB21(:), SH_LC21(:)
 COMPLEX*16, ALLOCATABLE :: SH_VR(:,:), SH_V(:,:)
!
 REAL*4, ALLOCATABLE :: BETA_S(:,:), BETA_U(:,:), BETA_G(:,:), BETA_V(:,:)
 REAL*8, ALLOCATABLE :: RA_AVE(:), RB_AVE(:), RC_AVE(:)
 REAL*4, ALLOCATABLE :: ES(:), EU(:), EG(:), EV(:), RR_AVE(:)
 REAL*8, ALLOCATABLE :: S_EQU(:), S_OFU(:), S_AVE(:)   	 
 REAL*8, ALLOCATABLE :: TOPO(:,:), TOPOP(:), CHI(:,:)
 REAL*8, ALLOCATABLE :: GAMMA_G(:), GAMMA_U(:), GAMMA_V(:)
 REAL*8, ALLOCATABLE :: PLM(:,:) 
 REAL*8, ALLOCATABLE :: Q(:), H(:,:), C_CON(:)
 INTEGER, ALLOCATABLE :: ANC(:)
!
 REAL*8 XN, XD
!
!--- Number of header lines in GF files
 INTEGER, PARAMETER :: NH_GF=8
 INTEGER, PARAMETER :: NH_GA=5
! 
!------------------------------------------------------------------------
! End of Declarations... End of Declarations... End of Declarations... 
!------------------------------------------------------------------------
!
!
!
!
! ========================= BEGIN OF CONFIGURATION SECTION =========================
! ========================= BEGIN OF CONFIGURATION SECTION =========================
! ========================= BEGIN OF CONFIGURATION SECTION =========================
!
 write(*,*)
 write(*,*) ' ****'
 write(*,*) ' **** This is the SELEN4 SLE solver'
 write(*,*) ' ****'
 write(*,*)
!
 if( iargc().ne.1 ) then
    write(*,*) ''
    write(*,*) ' Usage: sle.exe CONFIG_FILE '
    write(*,*) ''
    stop
 end if
!
 call getarg(1,cfg_f)
!
 inquire(file=trim(cfg_f),exist=ioflag)
 if( .not.ioflag ) then
    write(*,*) ''
    write(*,*) ' ERROR: Cannot open configuration file "', trim(cfg_f) ,'"'
    write(*,*) ''
    stop
 end if
!
 open(99,file=trim(cfg_f),status='old')
!
 write(*,*) ''
 write(*,*) ' >>>> Using configuration file: ', trim(cfg_f)
 write(*,*) ''
!
 open(99,file=trim(cfg_f),status='old')
!
! >>>>>>>>>>> Basic parameters
!
 call read_data_line(99,buffer)   ;   read(buffer,*) res
 call read_data_line(99,buffer)   ;   read(buffer,*) lmax
 call read_data_line(99,buffer)   ;   read(buffer,*) nn
 call read_data_line(99,buffer)   ;   read(buffer,*) nv
 call read_data_line(99,buffer)   ;   read(buffer,*) delta
 call read_data_line(99,buffer)   ;   read(buffer,*) irot
 call read_data_line(99,buffer)   ;   read(buffer,*) kref 
 call read_data_line(99,buffer)   ;   read(buffer,*) iext_max
 call read_data_line(99,buffer)   ;   read(buffer,*) iint_max
!
! >>>>>>>>>>> Input files
!
 call read_data_line(99,buffer)   ;   PT=trim(adjustl(buffer))
!
 call read_data_line(99,buffer)   ;   PIX_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   SHA_F=trim(pt)//trim(adjustl(buffer))
!call read_data_line(99,buffer)   ;   SGF_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   UGF_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   GGF_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   VGF_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   ROT_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   ICE_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   TOP_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   OFP_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   TGS_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   GEO_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   RSL_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   PMT_F=trim(pt)//trim(adjustl(buffer))
!
! >>>>>>>>>>> Output files
!
 call read_data_line(99,buffer)   ;   F_OUT_TOP=trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   F_OUT_OFU=trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   F_OUT_SHS=trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   F_OUT_SHU=trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   F_OUT_SHG=trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   F_OUT_SHN=trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   F_OUT_SHV=trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   F_OUT_TPW=trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   F_OUT_SAV=trim(adjustl(buffer))
!
 close(99)
!
! >>>>>>>>>>> Parameters depending on RES, LMAX, NV
!
 NP=2*RES*(RES-1)*20+12
 JMAX=(LMAX+1)*(LMAX+2)/2
 NMODES=4*NV
!
! >>>>>>>>>>> Set time-step for reference topo
!
 IF( KREF<0 )     KREF=NN+1+KREF  
 NN0 = KREF
!
!
! ========================== END OF CONFIGURATION SECTION ==========================
! ========================== END OF CONFIGURATION SECTION ==========================
! ========================== END OF CONFIGURATION SECTION ==========================
!
!
!
! ************************************************************************
! ************************************************************************
! ************************************************************************
!                           Execution starts here!
! ************************************************************************
! ************************************************************************
! ************************************************************************
!
!
 write(*,*) ''
 write(*,*) ' ------------------'
 write(*,*) ' >>>> 1. Data input' 
 write(*,*) ' ------------------'
!
 write(*,*) 
 write(*,*) ' >>>> Using these input data: '
 write(*,*) ' ---- Pixelized grid:     ', trim(adjustl(PIX_F))
 write(*,*) ' ---- SHs at the pixels:  ', trim(adjustl(SHA_F))
 write(*,*) ' ---- Geoid displac. GF:  ', trim(adjustl(GGF_F))
 write(*,*) ' ---- Vert. displac. GF:  ', trim(adjustl(UGF_F))
 write(*,*) ' ---- Horz. displac. GF:  ', trim(adjustl(VGF_F))
 write(*,*) ' ---- Rotational GFs      ', trim(adjustl(ROT_F))
 write(*,*) ' ---- Ice sheets history: ', trim(adjustl(ICE_F))
 write(*,*) ' ---- Present topo:       ', trim(adjustl(TOP_F))
 write(*,*) ' ---- Present OF:         ', trim(adjustl(OFP_F))
!
 OPEN(33,file=f_out_top, status="unknown",form="unformatted",ACCESS="SEQUENTIAL")
 OPEN(34,file=f_out_ofu, status="unknown",form="unformatted")
 OPEN(35,file=f_out_shs, status="unknown",form="unformatted")
 OPEN(36,file=f_out_shu, status="unknown",form="unformatted")
 OPEN(37,file=f_out_shg, status="unknown",form="unformatted")
 OPEN(38,file=f_out_shn, status="unknown",form="unformatted")
 OPEN(32,file=f_out_shv, status="unknown",form="unformatted")
 OPEN(39,file=f_out_tpw, status="unknown",form="formatted")
 OPEN(40,file=f_out_sav, status="unknown",form="formatted")
!
 write(*,*) 
 IF(IROT==0) write(*,*) ' >>>> No Rotation'
 IF(IROT==1) write(*,*) ' >>>> Adopting the classic rotation theory'
 IF(IROT==2) write(*,*) ' >>>> Adopting the revised rotation theory'
!
! ************************************************************************
! ************************************************************************
! ************************************************************************
!
!
 write(*,*) 
 write(*,*) ' >>>> Configuration of this run: ' 
 write(*,"(A27,3(I7,1X))") " ---- Ext/int iterations: ", IEXT_MAX, IINT_MAX
 write(*,"(A29,3(I7,1X))") " ---- Degrees l_max, j_max: ", LMAX, JMAX
 write(*,"(A29,3(I7,1X))") " ---- Number of v/e layers: ", NV
 if(NV==0) write(*,*) " ---- SELEN is running in <<Elastic Mode>> "
!
!
!----------------------------------------------------: Data input                                  
 write(*,*) ' ---- INPUT: Tegmark grid pixels data'
 ALLOCATE ( ANC(NP) )
 open(1,file=PIX_F) 
 do p=1, np 
     read (1,*) xn, xn, anc(p) 
 enddo
 na=anc(np)
 write(*,"(A36,3(I7,1X))") " ---- Res, pixels and main pixels: ", res, np, na
 Close(1)
!
!
!----------------------------------------------------: Data input                                  
 write(*,*) ' ---- INPUT: Spherical harmonics Y_jp'
 ALLOCATE ( PLM(JMAX,NA), SIN_COS(0:LMAX,NP) ) 
 open(1,file=SHA_F,form='unformatted') 
     read(1)PLM 
    read(1)sin_cos
 Close(1) 
 ALLOCATE ( MM(JMAX),LL(JMAX),DM(JMAX) )
 do j=1, jmax
    mm(j)=mj(j) ; ll(j)=lj(j) ; dm(j)=2-dom(j)
 enddo
!
!
!----------------------------------------------------: Data input                                  
 write(*,*) ' ---- INPUT: Present day topography T_pN'
 ALLOCATE ( TOPOP(NP) )  
 open(1,file=TOP_F) 
 do p=1, np  
     read(1,*) xn, xn, topop(p)  
 enddo
 Close(1)  
!
!
!----------------------------------------------------: Data input                                  
 IF(NV==0) THEN 
 write(*,*) ' ---- INPUT: Present day OF O_pN'
 ALLOCATE ( OFP(NP) )  
 open(1,file=OFP_F) 
 do p=1, np  
     read(1,*) xn, xn, ofp(p)  
 enddo
 Close(1)
 ENDIF 
!
!
!----------------------------------------------------: Data input                                  
 write(*,*) ' ---- INPUT: Ice thickness history I_pn'
 ALLOCATE ( H(NP,0:NN+1) ) 
 open(1,file=ICE_F) 
  do p=1, np  
    read(1,*) j, xn, xn, xn, (h(p,n),n=0,nn+1)
  enddo
 Close(1) 
!
!
!----------------------------------------------------: Data input                                  
 write(*,*) ' ---- INPUT: Coefficients beta^g_ln'
 write(*,*) ' ---- INPUT: Coefficients beta^u_ln'
 write(*,*) ' ---- INPUT: Coefficients beta^v_ln'
 ALLOCATE ( beta_g(0:lmax,0:nn+1), eg(0:lmax) )   
 ALLOCATE ( beta_u(0:lmax,0:nn+1), eu(0:lmax) )   
 ALLOCATE ( beta_s(0:lmax,0:nn+1), es(0:lmax) )   
 ALLOCATE ( beta_v(0:lmax,0:nn+1), ev(0:lmax) )
 beta_s(:,:)=0d0
 beta_g(:,:)=0d0
 beta_u(:,:)=0d0 
 beta_v(:,:)=0d0 
 open(1,file=GGF_F) 
 open(2,file=UGF_F)  
 open(3,file=VGF_F)  
 do i=1, nh_gf
    read(1,*)
	read(2,*)
	read(3,*)
 end do
 do i=0, lmax 
    read(1,*) l, eg(l), (beta_g(l,k), k=0,nn+1) 
    read(2,*) l, eu(l), (beta_u(l,k), k=0,nn+1) 
    read(3,*) l, ev(l), (beta_v(l,k), k=0,nn+1) 
 enddo
 do l=0, lmax  
  es(l) = eg(l)-eu(l)
    do k=0, nn+1 
	beta_s(l,k) = beta_g(l,k)-beta_u(l,k)
    enddo
 enddo
 Close(1)
 Close(2)  
 Close(3)  
!  
!
!----------------------------------------------------: Data input                                  
 write(*,*) ' ---- INPUT: Coefficients gamma^g_n'
 write(*,*) ' ---- INPUT: Coefficients gamma^u_n'
 write(*,*) ' ---- INPUT: Coefficients gamma^v_n'
 ALLOCATE ( GAMMA_G(0:NN+1), GAMMA_U(0:NN+1), GAMMA_V(0:NN+1) )     
 IF(IROT==0) then 
   gamma_g(:)=0D0 ; gamma_u(:)=0D0 ; gamma_v(:)=0D0
 endif
 IF(IROT==1.or.IROT==2) then 
 open(1,file=ROT_F) 
 do i=1, nh_ga
    read(1,*)
 end do
 do n=0, nn+1
    read(1,*) j, gamma_g(n), gamma_u(n), gamma_v(n)
 enddo
 ENDIF
 Close(1)
!
 write(*,*) 
 write(*,*) ' >>>> INPUT: Done reading data from external files' 
!
!
!
!
 write(*,*) ''
 write(*,*) ' ----------------------- '
 write(*,*) ' >>>> 2. Initialisation  ' 
 write(*,*) ' ----------------------- '
 write(*,*) ''
!
!
 IF(NV.NE.0) THEN 
! .........................................................
! In the ELASTIC case, the present topography is not needed 
! .........................................................
 write(*,*) ' ---- INIT: T_pn at', np, 'pixels'  
!----------------------------------------------------: Init: T_pn                                   
 ALLOCATE ( TOPO(NP,0:NN+1) )
 do p=1, np 
    do n=0, nn+1 
       topo(p,n)=dble(topop(p)) 
    enddo
 enddo
!
 Write(*,*) " ==== OUT: T_pn"
 WRITE(33)TOPO
!
 ENDIF
!
 write(*,*) ' ---- INIT: O_pn'  
!----------------------------------------------------: Init: O_pn                                   
!
 ALLOCATE( OF(NP,0:NN+1) ) 
!
 IF(NV.NE.0) THEN 
! ...............................................................
! In the general case of "moving shorelines", the startup OF is
! defined by the startup TOPO and by the startup ice THICKNESS. 
! ...............................................................
 do p=1, np  
     do n=0, nn+1 
       if(topo(p,n)+rho_iw*dble(h(p,n))< 0) of(p,n)=1 
       if(topo(p,n)+rho_iw*dble(h(p,n))>=0) of(p,n)=0 
    enddo
 enddo
!
 Write(*,*) " ==== OUT: O_pn"
 WRITE(34)OF
!
 ELSE 
! .................................................................
! In the ELASTIC case, the shorelines DO NOT MOVE. The startup OF 
! is obtained from an external file, and kept constant at all time
! steps ---- GS 011118   
! .................................................................
 do n=0, nn+1 
    of(:,n)=ofp(:)
 enddo
!
 Write(*,*) " ==== OUT: O_pn"
 WRITE(34)OF
!
 ENDIF
!
!
 write(*,*) ' ---- INIT: S^equ_n and S^ofu_n'  
!----------------------------------------------------: Init: S^equ_n                                 
 ALLOCATE ( S_EQU(0:NN+1) )
 do n=0,nn+1
    s_equ(n)=0d0 
    xn=0d0 
    xd=0d0 
       do p=1,np 				
          xn=xn + dble(h(p,n))*dble(1-of(p,n))-&
	          dble(h(p,0))*dble(1-of(p,0)) 
	  xd=xd + dble(of(p,n))
 	enddo
    s_equ(n)=-rho_iw*(xn/xd)
 enddo
!
!----------------------------------------------------: Init: S^ofu_n                                   
 ALLOCATE ( S_OFU(0:NN+1) ) 
 do n=0,nn+1
    s_ofu(n)=0d0 
    xn=0d0 
    xd=0d0 
       do p=1,np 				
          xn=xn + topo(p,0)*(dble(of(p,n))-dble(of(p,0)))	                
          xd=xd + dble(of(p,n))
       enddo
    s_ofu(n)=(xn/xd)
 enddo
!
 ALLOCATE ( S_AVE(0:NN+1) )
!
 s_ave(:) = s_equ(:) + s_ofu(:)
! 
!DEALLOCATE ( H )
!
!
 write(*,*) ' ---- INIT: S^ave_jn'  
!----------------------------------------------------: Init: S^ave_jn 
 ALLOCATE ( SH_SAV(JMAX,0:NN+1) )
 SH_SAV = (0.d0,0.d0)
 SH_SAV(J_INDEX(0,0),:) = S_AVE(:) 
!
!
 write(*,*) ' ---- INIT: S^(0)_jn' 
!----------------------------------------------------: Init: S^(0)_jn                                 
 ALLOCATE ( SH_S(JMAX,0:NN+1,0:IEXT_MAX) )
 SH_S(:,:,0) = SH_SAV(:,:)
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Memory allocation
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ [to be optimised] 
!
 ALLOCATE ( SH_OFU (JMAX,0:NN+1) )
 ALLOCATE ( SH_ZAV (JMAX,0:NN+1) )
!
 ALLOCATE ( SH_Z (JMAX,0:NN+1,0:IINT_MAX) )
 ALLOCATE ( SH_W (JMAX,0:NN+1) )
 ALLOCATE ( SH_X (JMAX,0:NN+1) )
!
 ALLOCATE ( SH_LA21 (0:NN+1) )
 ALLOCATE ( SH_LB21 (0:NN+1) )
 ALLOCATE ( SH_LC21 (0:NN+1) )
 ALLOCATE ( SH_LL21 (0:NN+1) )
!
 ALLOCATE ( SH_RAP (JMAX,0:NN+1) )  
 ALLOCATE ( SH_RBP (JMAX,0:NN+1) )
 ALLOCATE ( SH_RCP (JMAX,0:NN+1) )
 ALLOCATE ( SH_RRP (JMAX,0:NN+1) )
!
 ALLOCATE ( SH_KA (JMAX,0:NN+1) )
 ALLOCATE ( SH_KB (JMAX,0:NN+1) )
 ALLOCATE ( SH_KC (JMAX,0:NN+1) )
!
 ALLOCATE ( RB_AVE (0:NN+1) )
 ALLOCATE ( RA_AVE (0:NN+1) )
 ALLOCATE ( RC_AVE (0:NN+1) )
 ALLOCATE ( RR_AVE (0:NN+1) )
!
 ALLOCATE ( sh_gr (JMAX,0:NN+1) )
 ALLOCATE ( sh_ur (JMAX,0:NN+1) )
 ALLOCATE ( sh_vr (JMAX,0:NN+1) )
 ALLOCATE ( sh_rr (JMAX,0:NN+1) )
 ALLOCATE ( sh_kr (JMAX,0:NN+1) )
! 
 ALLOCATE ( Q(NP), C_CON(0:NN+1) )
 ALLOCATE ( PSI_RIG(0:NN+1) )
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Memory allocation
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ [to be optimised] 
!
!
! ====================================  External iteration
! ====================================  External iteration
! ====================================  External iteration
!
   DO 20000 IEXT=1, IEXT_MAX 

! ====================================  External iteration
! ====================================  External iteration
! ====================================  External iteration
!
 write(*,*) ""
 write(*,*) " ++++++++++++++++++++++++++++++++++++++++++++++++"
 write(*,*) " Ext loop: ", IEXT, 'of', IEXT_MAX     
 write(*,*) " ++++++++++++++++++++++++++++++++++++++++++++++++"
 write(*,*) ""
!
!
!---------------------------------------------------------: Ext loop: Q_p  
 write(*,*) ' ---- EXT: Computing Q_p' 
!ALLOCATE ( H(NP,0:NN+1) ) 
!open(1,file=ICE_F) 
!do p=1, np  
!   read(1,*) j, xn, xn, xn, (h(p,n),n=0,nn+1)
!enddo
!Close(1) 
!
 Q(:)=-(rho_ir*dble(h(:,0))+rho_wr*dble(topo(:,0)))
!
!DEALLOCATE ( H )
!
!
!---------------------------------------------------------: Ext loop: O_jn  
 write(*,*) ' ---- EXT: Computing O_jn'
 sh_ofu=(0d0,0d0)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,P) &
!$OMP SHARED(SH_OFU,PLM,ANC,OF,SIN_COS,MM,NN,NP,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1 
 if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1,np
       if(of(p,n)==1) then 
	      do j=1, jmax 
              sh_ofu(j,n) = sh_ofu(j,n) + & 
                  PLM(j,anc(p))*conjg(sin_cos(mm(j),p)) 
          end do
       end if
    end do
 end do
!$OMP END PARALLEL DO
 sh_ofu(:,:)=sh_ofu(:,:)/float(np)   
!
!
!---------------------------------------------------------: Ext loop: Z^ave_n 
 write(*,*) ' ---- EXT: Computing Z^av_n'
 do n=0,nn+1 
    sh_zav(:,n)=s_ave(n)*sh_ofu(:,n)  
 enddo   
!
!
 write(*,*) ' ---- - - - - - - - - - - - - - - '
!
!
!---------------------------------------------------------: Ext loop: W_jn  
 write(*,*) ' ---- EXT: Computing W_jn'
!ALLOCATE ( H(NP,0:NN+1) ) 
!open(1,file=ICE_F) 
!do p=1, np  
!   read(1,*) j, xn, xn, xn, (h(p,n),n=0,nn+1)
!enddo
!Close(1) 
! 
 sh_w(:,:)=(0d0,0d0)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(N,J,P) &
!$OMP SHARED(SH_W,H,PLM,ANC,OF,SIN_COS,MM,NN,NP,JMAX) &        
!$OMP SCHEDULE(GUIDED)
 do n=0,nn+1
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1, np
       if((h(p,n)-h(p,0)/=0).and.(of(p,n)/=1)) then	   
          do j=1, jmax 
             sh_w(j,n) = sh_w(j,n) + & 
 		         dble(h(p,n)-h(p,0))*dble(1-of(p,n))* & 
 		         PLM(j,anc(p))*conjg(sin_cos(mm(j),p)) 
	      end do
       end if
    end do
 end do
 sh_w(:,:)=sh_w(:,:)/float(np)   
!DEALLOCATE ( H )
!
!
!---------------------------------------------------------: Ext loop: L^a_21n  
 write(*,*) ' ---- EXT: Computing L^a_21n'
 sh_la21(:) = rhoi*sh_w(j_index(2,1),:)
!
!
!---------------------------------------------------------: Ext loop: R^a_jn  
 write(*,*) ' ---- EXT: Computing R^a_jn'
!
! +++ We define it as SH_RAP, and re-define below +++
!
 sh_rap(:,:)=(0d0,0d0)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_RAP,SH_W,LL,beta_s,ES,NN,NV,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_rap(j,n)= dble(es(ll(j)))*sh_w(j,n) 
       if(n>=1)then  
	      if(nv.ne.0) then          
            do k=0, n-1 
               sh_rap(j,n) = sh_rap(j,n)+&
                 (sh_w(j,k+1)-sh_w(j,k))*dble(beta_s(ll(j),n-k))
	        enddo 
          endif
	   endif
    enddo
 enddo	
!$OMP END PARALLEL DO
 sh_rap(:,:)=sh_rap(:,:)*3d0*rho_ie
! 
! 
!---------------------------------------------------------: Ext loop: <R^a>_n   
 write(*,*) ' ---- EXT: Computing <R^a>_n'
 ra_ave(:)=0D0
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(N,J) &
!$OMP SHARED(SH_OFU,RA_AVE,SH_RAP,DM,NN,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0,nn+1 
     ra_ave(n)=0D0 
     do j=1, jmax 
          ra_ave(n)=ra_ave(n)+ &
	         dm(j)*real(conjg(sh_ofu(j,n))*sh_rap(j,n))&
	         /sh_ofu(j_index(0,0),n) 	
     end do 
 end do
!$OMP END PARALLEL DO
!
!
!---------------------------------------------------------: Ext loop: R^ap_jn  
 write(*,*) ' ---- EXT: Computing R^pa_jn'
 sh_rap(j_index(0,0),:)=sh_rap(j_index(0,0),:)-ra_ave(:)
!
!
!---------------------------------------------------------: Ext loop: chi^a_pn  
 write(*,*) ' ---- EXT: Computing chi^a_pn'
 ALLOCATE ( CHI(NP,0:NN+1) )
 chi = 0d0
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,N,J) &
!$OMP SHARED(CHI,SH_RAP,PLM,SIN_COS,DM,MM,ANC,NP,NN,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
   if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
   do p=1,np 
     do j=1,jmax 
          chi(p,n)=chi(p,n)+&
             PLM(j,anc(p))*dm(j)*&
	         real(sh_rap(j,n)*sin_cos(mm(j),p))	
     end do 
   end do
 end do
!$OMP END PARALLEL DO
 chi(:,:)=chi(:,:)/float(np)
!
!
!---------------------------------------------------------: Ext loop: K^a_jn  
 write(*,*) ' ---- EXT: Computing K^a_jn' 
 sh_ka(:,:)=(0D0,0D0) 
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,N,J) &
!$OMP SHARED(SH_KA,PLM,CHI,SIN_COS,MM,OF,ANC,JMAX,NP,NN) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1 
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1, np 	
       if(of(p,n)==1) then 
          do j=1,jmax
		      sh_ka(j,n) = sh_ka(j,n)+ & 
                  PLM(j,anc(p))*chi(p,n)*conjg(sin_cos(mm(j),p))	
          end do
	   end if
    end do
 end do 
 !$OMP END PARALLEL DO
 DEALLOCATE ( CHI )
!
!
 write(*,*) ' ---- - - - - - - - - - - - - - - '
!
!
!---------------------------------------------------------: Ext loop: X_jn  
 write(*,*) ' ---- EXT: Computing X_jn'
 sh_x=(0d0,0d0) 
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(N,J,P) &
!$OMP SHARED(SH_X,PLM,ANC,OF,Q,SIN_COS,MM,NN,JMAX,NP) &        
!$OMP SCHEDULE(GUIDED)
 do n=0,nn+1
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1,np
       if((of(p,n)-of(p,0))/=0) then
           do j=1, jmax 
	          sh_x(j,n)= sh_x(j,n) + & 
		            dble(of(p,n)-of(p,0))*Q(p)* & 
                    PLM(j,anc(p))*conjg(sin_cos(mm(j),p)) 
           end do
       end if 
    end do
 end do
!$OMP END PARALLEL DO
 sh_x(:,:)=sh_x(:,:)/float(np)   
!
!
!---------------------------------------------------------: Ext loop: L^c_21n  
 write(*,*) ' ---- EXT: Computing L^c_21n'
 sh_lc21(:) = rhor*sh_x(j_index(2,1),:)
!
!
!---------------------------------------------------------: Ext loop: R^c_jn  
 write(*,*) ' ---- EXT: Computing R^c_jn'
!
! +++ We define it as SH_RCP, and re-define below +++
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_RCP,SH_X,LL,beta_s,ES,JMAX,NN,NV) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_rcp(j,n)= dble(es(ll(j)))*sh_x(j,n) 
       if(n>=1)then  
	      if(nv.ne.0) then          
            do k=0, n-1 
                sh_rcp(j,n) = sh_rcp(j,n)+&
               (sh_x(j,k+1)-sh_x(j,k))*dble(beta_s(ll(j),n-k))
	        enddo 
          endif
	   endif
    enddo
 enddo	
!$OMP END PARALLEL DO
 sh_rcp(:,:)=sh_rcp(:,:)*(3d0*rho_re)
! 
! 
!---------------------------------------------------------: Ext loop: <R^c>_n   
 write(*,*) ' ---- EXT: Computing <R^c>_n'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(N,J) &
!$OMP SHARED(SH_OFU,RC_AVE,SH_RCP,DM,NN,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0,nn+1 
    rc_ave(n)=0D0 
       do j=1, jmax 
          rc_ave(n)=rc_ave(n)+ &
          dm(j)*real(conjg(sh_ofu(j,n))*sh_rcp(j,n))& 
	  /sh_ofu(j_index(0,0),n) 	
       enddo 
 enddo
!$OMP END PARALLEL DO
!
!
!---------------------------------------------------------: Ext loop: R^cp_jn  
 write(*,*) ' ---- EXT: Computing R^pc_jn'
 sh_rcp(j_index(0,0),:)=sh_rcp(j_index(0,0),:)-rc_ave(:)
!
!                                                                
!---------------------------------------------------------: Ext loop: chi^c_pn  
 write(*,*) ' ---- EXT: Computing chi^c_pn'
 ALLOCATE ( CHI(NP,0:NN+1) )
 chi = 0d0
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,N,J) &
!$OMP SHARED(CHI,SH_RCP,PLM,SIN_COS,DM,MM,ANC,NP,NN,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1 
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1,np 
       do j=1,jmax 
           chi(p,n)=chi(p,n)+&
              PLM(j,anc(p))*dm(j)*real(sh_rcp(j,n)*&
	          sin_cos(mm(j),p))	
       end do 
    end do
 end do
!$OMP END PARALLEL DO
 chi(:,:)=chi(:,:)/float(np)
!
!
!---------------------------------------------------------: Ext loop: K^c_jn  
 write(*,*) ' ---- EXT: Computing K^c_jn' 
 sh_kc = (0d0, 0d0)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,N,J) &
!$OMP SHARED(SH_KC,PLM,CHI,SIN_COS,MM,OF,ANC,JMAX,NN,NP) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1 
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1, np 	
       if(of(p,n)==1) then
	      do j=1,jmax 
             sh_kc(j,n) = sh_kc(j,n)+ & 
                PLM(j,anc(p))*chi(p,n)*conjg(sin_cos(mm(j),p))	
          end do
       end if
    end do
 end do
!$OMP END PARALLEL DO
 DEALLOCATE ( CHI )
!
!
!
 write(*,*) ' ---- - - - - - - - - - - - - - - '
!
!
!------------------------------------------: Initialising the internal loop   
!
!
 write(*,*) ' ---- Initializing the INTERNAL LOOP'
 write(*,*) ' ---- EXT: Computing Z^0_jn' 
 sh_z(:,:,:) = (0d0,0d0)
 sh_z(:,:,0) = sh_zav(:,:) 
!
!

! ============================  Internal iteration
! ============================  Internal iteration
! ============================  Internal iteration
!
   DO 10000 IINT=1, IINT_MAX 

! ============================  Internal iteration
! ============================  Internal iteration
! ============================  Internal iteration
!
  write(*,*) ""
  write(*,*) " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write(*,*) " Internal loop ", IINT, 'of', IINT_MAX     
  write(*,*) " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write(*,*) ""
!
!---------------------------------------------------------: INT loop: L^b_21n  
 write(*,*) ' ---- INT: Computing L^b_21n'
 sh_lb21(:)=rhow*sh_z(j_index(2,1),:,iint-1)
!
!
!---------------------------------------------------------: INT loop: R^b_jn  
 write(*,*) ' ---- INT: Computing R^b_jn'
!
! +++ We define it as SH_RBP, and re-define below +++
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_RBP,LL,beta_s,ES,SH_Z,IINT,JMAX,NN,NV) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_rbp(j,n)= dble(es(ll(j)))*sh_z(j,n,iint-1) 
         if(n>=1)then  
	    if(nv.ne.0) then          
            do k=0, n-1 
               sh_rbp(j,n) = sh_rbp(j,n)+&
              (sh_z(j,k+1,iint-1)-sh_z(j,k,iint-1))*dble(beta_s(ll(j),n-k))
            enddo 
            endif
         endif
    enddo
 enddo	
! $OMP END PARALLEL DO
 sh_rbp(:,:)=sh_rbp(:,:)*3d0*rho_we
! 
! 
!---------------------------------------------------------: INT loop: <R^b>_n   
 write(*,*) ' ---- INT: Computing <R^b>_n'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(N,J) &
!$OMP SHARED(SH_OFU,RB_AVE,SH_RBP,DM,NN,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0,nn+1 
    rb_ave(n)=0D0 
       do j=1, jmax 
          rb_ave(n)=rb_ave(n)+ &
	  dm(j)*real(conjg(sh_ofu(j,n))*&
	  sh_rbp(j,n))/sh_ofu(j_index(0,0),n) 	
       enddo 
 enddo
!$OMP END PARALLEL DO
!
!
!---------------------------------------------------------: INT loop: R^bp_jn  
 write(*,*) ' ---- INT: Computing R^pb_jn'
 sh_rbp(j_index(0,0),:)=sh_rbp(j_index(0,0),:)-rb_ave(:)
!
!                                                                                                                  
!---------------------------------------------------------: INT loop: chi^b_pn  
 write(*,*) ' ---- INT: Computing chi^b_pn'
 ALLOCATE ( CHI(NP,0:NN+1) )
 chi=0d0
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,N,J) &
!$OMP SHARED(CHI,SH_RBP,PLM,SIN_COS,DM,MM,ANC,NP,NN,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1 
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1,np 
        do j=1,jmax 
           chi(p,n)=chi(p,n)+&
              PLM(j,anc(p))*dm(j)*real(sh_rbp(j,n)*sin_cos(mm(j),p))	
        end do 
    end do
 end do
!$OMP END PARALLEL DO
 chi(:,:)=chi(:,:)/float(np)
!
!
!---------------------------------------------------------: INT loop: K^b_jn  
 write(*,*) ' ---- INT: Computing K^b_jn' 
 sh_kb=(0D0,0D0) 
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,N,J) &
!$OMP SHARED(SH_KB,PLM,CHI,SIN_COS,MM,OF,ANC,JMAX,NN,NP) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1 
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1, np 	
       if(of(p,n)==1) then
	      do j=1,jmax 
             sh_kb(j,n)=sh_kb(j,n)+ & 
                PLM(j,anc(p))*chi(p,n)*conjg(sin_cos(mm(j),p))	
          end do
       end if
    end do	   
 end do
!$OMP END PARALLEL DO
 DEALLOCATE ( CHI )
!
!---------------------------------------------------------: INT loop: L_jn
 write(*,*) ' ---- INT: Computing L_21n'
 sh_ll21(:) = sh_la21(:) + sh_lb21(:) + sh_lc21(:) 
!
!
!---------------------------------------------------------: INT loop: (G, U)^rot_jn
 IF(IROT==1.or.IROT==2) THEN
   write(*,*) ' ---- INT: Computing G^rot_jn, U^rot_jn and V^rot_jn'
   sh_gr(:,:)=(0d0,0d0)
   sh_ur(:,:)=(0d0,0d0)
   sh_vr(:,:)=(0d0,0d0)   
   do n=0, nn+1 
     sh_gr(j_index(2,1),n)=(0d0,0d0)
     sh_ur(j_index(2,1),n)=(0d0,0d0)
     sh_vr(j_index(2,1),n)=(0d0,0d0)
     if(n>=1) then 
        do k=0, n-1 
          sh_gr(j_index(2,1),n) = sh_gr(j_index(2,1),n) + &
	                   C21_PSI*(sh_ll21(k+1)-sh_ll21(k))*gamma_g(n-k)
          sh_ur(j_index(2,1),n) = sh_ur(j_index(2,1),n) + & 
	                   C21_PSI*(sh_ll21(k+1)-sh_ll21(k))*gamma_u(n-k)
          sh_vr(j_index(2,1),n) = sh_vr(j_index(2,1),n) + &
	                   C21_PSI*(sh_ll21(k+1)-sh_ll21(k))*gamma_v(n-k)
        enddo
     endif 
   enddo
!
!
!---------------------------------------------------------: INT loop: R^rot_jn   
   write(*,*) ' ---- INT: Computing R^rot_jn' 
   sh_rr(:,:)=(0d0,0d0)
   sh_rr(j_index(2,1),:)=sh_gr(j_index(2,1),:)-sh_ur(j_index(2,1),:) 
 ELSE
   sh_gr(:,:)=(0d0,0d0)
   sh_ur(:,:)=(0d0,0d0)
   sh_vr(:,:)=(0d0,0d0)
   sh_rr(:,:)=(0d0,0d0) 
 ENDIF
!
!
!---------------------------------------------------------: INT loop: <R^rot>_n 
 IF(IROT==1.or.IROT==2) THEN
 write(*,*) ' ---- INT: Computing <R^rot>'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(N,J) &
!$OMP SHARED(SH_OFU,RR_AVE,sh_rr,DM,NN,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0,nn+1 
   RR_AVE(n)=0D0
    do j=1, jmax 
       RR_AVE(n)=RR_AVE(n)+ &
       dm(j)*real(conjg(sh_ofu(j,n))*sh_rr(j,n))/&
       sh_ofu(j_index(0,0),n) 	
    enddo 
 enddo
!$OMP END PARALLEL DO
 ELSE
  RR_AVE(:)=0D0 
 ENDIF
!
!
!---------------------------------------------------------: INT loop: R^rotp_jn 
 SH_RRP(:,:)=(0d0,0d0) 
 IF(IROT==1.or.IROT==2) THEN
 write(*,*) ' ---- INT: Computing R^rotp_jn'
  SH_RRP(:,:)=sh_rr(:,:)
  SH_RRP(j_index(0,0),:)=sh_rr(j_index(0,0),:)-RR_AVE(:)
 ENDIF
!
!
!---------------------------------------------------------: INT loop: chi^rot_pn  
 write(*,*) ' ---- INT: Computing chi^rot_pn'
 ALLOCATE ( CHI(NP,0:NN+1) )
 CHI = 0d0
 IF(IROT==1.or.IROT==2)THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,N,J) &
!$OMP SHARED(CHI,SH_RRP,PLM,SIN_COS,DM,MM,ANC,NP,NN,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1 
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1,np 
        do j=1,jmax 
           chi(p,n)=chi(p,n)+&
              PLM(j,anc(p))*dm(j)*real(SH_RRP(j,n)*sin_cos(mm(j),p))	
        end do 
    end do
 end do
!$OMP END PARALLEL DO
 chi(:,:)=chi(:,:)/float(np)
 ELSE
     chi(:,:)=0d0
 ENDIF
!
!
!---------------------------------------------------------: INT loop: K^rot_jn  
 write(*,*) ' ---- INT: Computing K^rot_jn' 
 sh_kr(:,:)=(0D0,0D0) 
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,N,J) &
!$OMP SHARED(sh_kr,PLM,CHI,SIN_COS,MM,OF,ANC,JMAX,NP,NN) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1, np 	
        if(of(p,n)==1) then
		   do j=1,jmax 
              sh_kr(j,n) = sh_kr(j,n)+ & 
                 PLM(j,anc(p))*chi(p,n)*conjg(sin_cos(mm(j),p))	
	       end do
        end if 
    end do
 end do
!$OMP END PARALLEL DO
 DEALLOCATE ( CHI )
!
!
 write(*,*) ' ---- - - - - - - - - - - - - - - '
!
!
!------------------------------------------------------: INT loop: Psi^rig_n  
 write(*,*) ' ---- INT: Computing Psi^rig_n' 
 PSI_RIG(:)=(0d0,0d0)
 IF(IROT==1.or.IROT==2) PSI_RIG(:)=C21_PSI*conjg(sh_ll21(:))     
!
!
!------------------------------------------------------: INT loop: c_n  
 write(*,*) ' ---- INT: Computing c_n' 
 c_con(:) = s_ave(:) - ra_ave(:) - & 
                         rb_ave(:) - & 
                           rc_ave(:) - & 
                             rr_ave(:) 
!
!------------------------------------------------------: INT loop: Updating Z_jn 
 write(*,*) ' ---- INT: Computing Z^(IINT)_jn'
 SH_Z(:,:,IINT) =  sh_zav(:,:) +  & 
                     sh_kr(:,:) + &
                       sh_ka(:,:)  +  &
                         sh_kb(:,:)  +  &
                           sh_kc(:,:) 
!
 write(*,*) 
 write(*,*) " ====================" 
 write(*,*) " END of INTERNAL loop " 
 write(*,*) " ====================" 
! 
! End of INTERNAL loop ! End of INTERNAL loop ! End of INTERNAL loop 
! +++++++++++++++++++
   10000 CONTINUE 
! +++++++++++++++++++
!
!
!
 write(*,*) 
!-----------------------------------------------------: Ext loop: S^(IEXT)_jn
 write(*,*) ' ---- EXT: Computing S^(IEXT)_jn'
 SH_S(:,:,IEXT) = sh_sav(:,:) +   & 
                    sh_rrp(:,:) +   & 
		      sh_rap(:,:) +     &
		        sh_rbp(:,:) +      & 
		          sh_rcp(:,:)
!
!----------------------------------------------------: Ext loop: S^(IEXT)_pn 
 write(*,*) ' ---- EXT: Computing S^(IEXT)_pn' 
!
! * ** *** We use "CHI" to save memory *** ** * 
! 
 ALLOCATE ( CHI(NP,0:NN+1) )   
 chi = 0d0
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,N,J) &
!$OMP SHARED(CHI,SH_S,PLM,SIN_COS,DM,MM,ANC,IEXT,NP,NN,JMAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1 
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1, np 
 	   do j=1, jmax 
          chi(p,n)=chi(p,n) + & 
             PLM(j,anc(p))*dm(j)*&
	         real(SH_S(J,N,IEXT)*sin_cos(mm(j),p))	
	   end do
    end do
 end do
!$OMP END PARALLEL DO
!
!
!----------------------------------------------------: Ext loop: T_pn 
 write(*,*) ' ---- EXT: new topography T_pn'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,N) &
!$OMP SHARED(TOPO,TOPOP,CHI,NP,NN,NN0) &
!$OMP SCHEDULE(GUIDED)
!
! * ** *** We use "CHI" to save memory *** ** * 
! 
 do n=0, nn+1 
    if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
    do p=1, np 
 !if(p==1.or.p==np)write(*,"(A9,1X,I7)") "p=", p 
       topo(p,n) = dble(topop(p))-(chi(p,n)-chi(p,nn0)) 
    enddo
 enddo
!$OMP END PARALLEL DO
 DEALLOCATE ( CHI )
!
 Write(*,*) " ==== OUT: T_pn"
 WRITE(33)TOPO
!
!
!----------------------------------------------------: Ext loop: O_pn 
 write(*,*) ' ---- EXT: new Ocean Function O_pn' 
!ALLOCATE ( H(NP,0:NN+1) ) 
!open(1,file=ICE_F) 
!do p=1, np  
!   read(1,*) j, xn, xn, xn, (h(p,n),n=0,nn+1)
!enddo
!Close(1) 
! 
 do n=0, nn+1 
 if(n==0.or.n==nn)write(*,"(A9,1X,I3)") "n=", n
 !if(p==1.or.p==np)write(*,"(A9,1X,I7)") "p=", p 
   do p=1, np  
     of(p,n)=0
     if(topo(p,n)+rho_iw*dble(h(p,n))< 0) of(p,n)=1 
     if(topo(p,n)+rho_iw*dble(h(p,n))>=0) of(p,n)=0 
   enddo
 enddo
!
 Write(*,*) " ==== OUT: O_pn"
 WRITE(34)OF
!
!
!----------------------------------------------------: Ext loop: S^ave_n                                  
 write(*,*) ' ---- EXT: New S^ave_n'  
! ALLOCATE ( S_EQU(0:NN+1) )
 do n=0,nn+1
    s_equ(n)=0d0 
    xn=0d0 ; xd=0d0 
       do p=1,np 				
          xn=xn + dble(h(p,n))*dble(1-of(p,n))-&
                  dble(h(p,0))*dble(1-of(p,0)) 
	  xd=xd + dble(of(p,n))
 	enddo
    s_equ(n)=-rho_iw*(xn/xd)
 enddo
!----------------------------------------------------: Ext loop: S^ofu_n                                  
! ALLOCATE ( S_OFU(0:NN+1) )
 do n=0,nn+1
    s_ofu(n)=0d0 
    xn=0d0 ; xd=0d0 
       do p=1,np 				
          xn=xn + topo(p,0)*(dble(of(p,n))-dble(of(p,0)))	                
	  xd=xd + dble(of(p,n))
       enddo
    s_ofu(n)=(xn/xd)
 enddo
!----------------------------------------------------: Ext loop: S^ave_n                                  
!
 s_ave(:) = s_equ(:) + s_ofu(:)
!
!----------------------------------------------------: Ext loop: S^ave_jn 
 SH_SAV(:,:) = (0D0,0D0)
 SH_SAV(J_INDEX(0,0),:) = S_AVE(:) 
!
!DEALLOCATE ( H )
!
 write(*,*) 
 write(*,*) " o~o~o~o~o~o~o~o~o~o~o" 
 write(*,*) " END of Ext loop " 
 write(*,*) " o~o~o~o~o~o~o~o~o~o~+" 
! 
! End of Ext loop ! End of Ext loop ! End of Ext loop 
! +++++++++++++++++++
   20000 CONTINUE 
! +++++++++++++++++++
! End of Ext loop ! End of Ext loop ! End of Ext loop 
!
 write(*,*) 
 write(*,*) ' **** END of ALL ITERATIONS'
 write(*,*) 
!
!
!
 Write(*,*) " ==== OUT: S_jn"
 WRITE(35)SH_S 
! 
!
!===========!===========!===========!===========!======= : Closing 
!===========!===========!===========!===========!======= : Closing 
!===========!===========!===========!===========!======= : Closing 
!
!------------------------------------------------------- : Generic SH's 
 ALLOCATE ( SH_XA(JMAX,0:NN+1), & 
            SH_XB(JMAX,0:NN+1), & 
            SH_XC(JMAX,0:NN+1) )
!
!
!
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++ : Horizontal displacement 
!
!------------------------------------------------------- : Closing: V^a_jn
 write(*,*) ' >>>> Closing: V^a_jn'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_XA,SH_W,LL,beta_v,EV,JMAX,NN) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_xa(j,n)= dble(ev(ll(j)))*sh_w(j,n) 
       if(n>=1)then
          do k=0, n-1 
             sh_xa(j,n) = sh_xa(j,n)+&
            (sh_w(j,k+1)-sh_w(j,k))*dble(beta_v(ll(j),n-k))
         enddo 
       endif
    enddo
 enddo	
!$OMP END PARALLEL DO
 sh_xa(:,:)=sh_xa(:,:)*3d0*rho_ie
!
!
!-------------------------------------------------------: Closing: V^b_jn  
 write(*,*) ' >>>> Closing: V^b_jn'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_XB,SH_Z,LL,beta_v,EV,JMAX,NN,IINT_MAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_xb(j,n)= dble(ev(ll(j)))*sh_z(j,n,IINT_MAX) 
       if(n>=1)then
          do k=0, n-1 
             sh_xb(j,n) = sh_xb(j,n) + &
             (sh_z(j,k+1,IINT_MAX)-sh_z(j,k,IINT_MAX))*&
	     			dble(beta_v(ll(j),n-k))
          enddo 
       endif
    enddo
 enddo	
! $OMP END PARALLEL DO
 sh_xb(:,:)=sh_xb(:,:)*3d0*rho_we
!
!
!------------------------------------------------------- : Closing: V^c_jn
 write(*,*) ' >>>> Closing: V^c_jn'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_XC,SH_X,LL,beta_v,EV,JMAX,NN) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_xc(j,n)= dble(ev(ll(j)))*sh_x(j,n) 
       if(n>=1)then
          do k=0, n-1 
             sh_xc(j,n) = sh_xc(j,n)+&
             (sh_x(j,k+1)-sh_x(j,k))*dble(beta_v(ll(j),n-k))
         enddo 
	endif
    enddo
 enddo	
! $OMP END PARALLEL DO
 sh_xc(:,:)=sh_xc(:,:)*3d0*rho_re
!
!
 ALLOCATE ( SH_V(JMAX,0:NN+1) )
!
 sh_v(:,:)= sh_xa(:,:)+ sh_xb(:,:)+ sh_xc(:,:) + sh_vr(:,:)
!
 Write(*,*) " ==== OUT: V_jn"
 WRITE(32)SH_V 
!
!
!
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++ : Vertical displacement 			  
!
!------------------------------------------------------- : Closing: U^a_jn
 write(*,*) ' >>>> Closing: U^a_jn'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_XA,SH_W,LL,beta_u,EU,JMAX,NN) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_xa(j,n)= dble(eu(ll(j)))*sh_w(j,n) 
       if(n>=1)then
          do k=0, n-1 
             sh_xa(j,n) = sh_xa(j,n)+&
            (sh_w(j,k+1)-sh_w(j,k))*dble(beta_u(ll(j),n-k))
         enddo 
       endif
    enddo
 enddo	
!$OMP END PARALLEL DO
 sh_xa(:,:)=sh_xa(:,:)*3d0*rho_ie
!
!
!-------------------------------------------------------: Closing: U^b_jn  
 write(*,*) ' >>>> Closing: U^b_jn'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_XB,SH_Z,LL,beta_u,EU,JMAX,NN,IINT_MAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_xb(j,n)= dble(eu(ll(j)))*sh_z(j,n,IINT_MAX) 
       if(n>=1)then
          do k=0, n-1 
             sh_xb(j,n) = sh_xb(j,n) + &
             (sh_z(j,k+1,IINT_MAX)-sh_z(j,k,IINT_MAX))*&
	     dble(beta_u(ll(j),n-k))
          enddo 
       endif
    enddo
 enddo	
! $OMP END PARALLEL DO
 sh_xb(:,:)=sh_xb(:,:)*3d0*rho_we
!
!
!------------------------------------------------------- : Closing: U^c_jn
 write(*,*) ' >>>> Closing: U^c_jn'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_XC,SH_X,LL,beta_u,EU,JMAX,NN) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_xc(j,n)= dble(eu(ll(j)))*sh_x(j,n) 
       if(n>=1)then
          do k=0, n-1 
             sh_xc(j,n) = sh_xc(j,n)+&
             (sh_x(j,k+1)-sh_x(j,k))*dble(beta_u(ll(j),n-k))
         enddo 
	endif
    enddo
 enddo	
! $OMP END PARALLEL DO
 sh_xc(:,:)=sh_xc(:,:)*3d0*rho_re
!
!
 ALLOCATE ( SH_U(JMAX,0:NN+1) )
!
 sh_u(:,:)= sh_xa(:,:)+ sh_xb(:,:)+ sh_xc(:,:) + sh_ur(:,:)
!
 Write(*,*) " ==== OUT: U_jn"
 WRITE(36)SH_U
!
!
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++ : Geoid height  
!
!------------------------------------------------------- : Closing: G^a_jn
 write(*,*) ' >>>> Closing: G^a_jn'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_XA,SH_W,LL,beta_g,EG,JMAX,NN) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_xa(j,n)= dble(eg(ll(j)))*sh_w(j,n) 
       if(n>=1)then
          do k=0, n-1 
             sh_xa(j,n) = sh_xa(j,n)+&
            (sh_w(j,k+1)-sh_w(j,k))*dble(beta_g(ll(j),n-k))
          enddo 
       endif
    enddo
 enddo	
! $OMP END PARALLEL DO
 sh_xa(:,:)=sh_xa(:,:)*3d0*rho_ie
!
!
!-------------------------------------------------------: Closing: G^b_jn  
 write(*,*) ' >>>> Closing: G^b_jn'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_XB,SH_Z,LL,beta_g,EG,JMAX,NN,IINT_MAX) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_xb(j,n)= dble(eg(ll(j)))*sh_z(j,n,IINT_MAX) 
       if(n>=1)then
          do k=0, n-1 
             sh_xb(j,n) = sh_xb(j,n) + &
             (sh_z(j,k+1,IINT_MAX)-sh_z(j,k,IINT_MAX))*&
	     dble(beta_g(ll(j),n-k))
          enddo 
       endif
    enddo
 enddo	
! $OMP END PARALLEL DO
 sh_xb(:,:)=sh_xb(:,:)*3d0*rho_we
!
!
!------------------------------------------------------- : Closing: G^c_jn
 write(*,*) ' >>>> Closing: G^c_jn'
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,N,K) &
!$OMP SHARED(SH_XC,SH_X,LL,beta_g,EG,JMAX,NN) &
!$OMP SCHEDULE(GUIDED)
 do n=0, nn+1
    do j=1, jmax 
       sh_xc(j,n)= dble(eg(ll(j)))*sh_x(j,n) 
       if(n>=1)then
          do k=0, n-1 
             sh_xc(j,n) = sh_xc(j,n)+&
             (sh_x(j,k+1)-sh_x(j,k))*dble(beta_g(ll(j),n-k))
          enddo 
       endif
    enddo
 enddo	
!$OMP END PARALLEL DO
 sh_xc(:,:)=sh_xc(:,:)*3d0*rho_re
!
!
 ALLOCATE ( SH_G(JMAX,0:NN+1) )
!
 sh_g(:,:)= sh_xa(:,:)+ sh_xb(:,:)+ sh_xc(:,:) + sh_gr(:,:)
!
 Write(*,*) " ==== OUT: G_jn"
 WRITE(37)SH_G
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++ : Absolute SLC  
!
! ------------------------------------------------------ : Closing: N^_jn
 write(*,*) ' >>>> Closing: N_jn' 
!
 ALLOCATE( SH_N(JMAX,0:NN+1) )   
!
 sh_n(:,:)=sh_g(:,:)

 sh_n(j_index(0,0),:) = sh_g(j_index(0,0),:) + c_con(:) 
!
 Write(*,*) " ==== OUT: N_jn"
 WRITE(38)SH_N
!
! ------------------------------------------------------ : Closing: psi^rig_n
 Write(*,*) " ==== OUT: psi^rig_n"
 do k=0, nn+1 
     write(39,*) k,  psi_rig(k) 		  
 enddo
!
! ------------------------------------------------------ : Closing: s^ave_n
 Write(*,*) " ==== OUT: s^ave_n"
 do k=0, nn+1
     write(40,*) k, s_equ(k), s_ofu(k), s_ave(k)  		  
 enddo
!
!
 write(*,*) 
 write(*,*) " -------------------------" 
 write(*,*) "     END of EXECUTION      " 
 write(*,*) " -------------------------" 
 write(*,*) ""
 write(*,*) " good bye and good luck..." 
 write(*,*) ""
!
 Close (33) 
  Close (34) 
   Close (35) 
    Close (36) 
    Close (37) 
   Close (38) 
  Close (39) 
 Close (40) 
!
 STOP
!
! ++++++++++++++++++++++
 END Program SELEN4
! ++++++++++++++++++++++
!
!
