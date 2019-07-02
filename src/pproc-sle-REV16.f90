!
!-*-*-*-*-*-*-*-*
!
  PROGRAM PPROC 
!
!-*-*-*-*-*-*-*-*
!
! REVISION 16 - Jun 27, 2019
!
! Postprocessing the SLE 
! Author: GS July 2016 
! Revised GS Aug 2018... 
! Revised GS April 2019, for version 10 
! Revised DM May 2, 2019, to skip the PM analysis if IROT=0
! Revised DM May 3, 2019, to read the configuration from the CMD-line
! Revised DM Jun 19, 2019:   Removed the reference to the S Green Function
! Revised DM Jun 27, 2019:   Changed the TIME_STAMP subroutine to handle
!                            time steps different from 1.0kyr
!                            Increased format width for 'l' & 'm' in stokes.dat
!
!/////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////
!
!
 IMPLICIT NONE 
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
 CHARACTER*100 :: GGF_F
 CHARACTER*100 :: ROT_F
 CHARACTER*100 :: ICE_F
 CHARACTER*100 :: TOP_F
 CHARACTER*100 :: OFP_F
 CHARACTER*100 :: TGS_F
 CHARACTER*100 :: RSL_F
 CHARACTER*100 :: PMT_F
!
! >>>>>>> Output files --------------------------------------
!
 CHARACTER*100 :: F_OUT_TOP
 CHARACTER*100 :: F_OUT_OFU
 CHARACTER*100 :: F_OUT_SHS
 CHARACTER*100 :: F_OUT_SHU
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
 INTEGER :: NP
 INTEGER :: JMAX
 INTEGER :: NMODES
!
 REAL*8,  PARAMETER :: ERAD=6.371D6
 REAL*8,  PARAMETER :: RHOW=1000D0    
 REAL*8,  PARAMETER :: RHOI=931D0    
 REAL*8,  PARAMETER :: RHOR=900D0   
 REAL*8,  PARAMETER :: RHOE=5514D0
 REAL*8,  PARAMETER :: CCC=8.0394d37 
 REAL*8,  PARAMETER :: AAA=8.0131d37 
 REAL*8,  PARAMETER :: CMA=CCC-AAA
 REAL*8,  PARAMETER :: PI=3.14159265358979323840D0 
 REAL*8,  PARAMETER :: C0=(4D0*PI/3D0)*SQRT(6D0/5D0)
 INTEGER, PARAMETER :: NLAT=180, NLON=360  
 REAL*8,  PARAMETER :: ERADIUS=6.371D6   
!
 REAL*8 ALFA
 INTEGER :: ST_LMIN, ST_LMAX 
! 
 CHARACTER*32 FILE_TOPO, FILE_ICE, FILE_IFL, FILE_IGB 
 CHARACTER*32 FILE_IGA, FILE_OCE, FILE_CON
 CHARACTER*4 CJUNK 
 CHARACTER*4, ALLOCATABLE :: TBP(:)
 REAL*8 DDELTA, XJ, PSI_X, PSI_Y, C_STOK, S_STOK 
 REAL*8 AVERAGE_N, AVERAGE_U, AVERAGE_S, AVERAGE_G
 INTEGER I, N, J, K, L, M, P, IT, NA, NBP
 INTEGER J_INDEX, NPU, ILAT, ILON, NIN, MJ, LJ, DOM  
 COMPLEX*16 , ALLOCATABLE :: PSI(:) 
 COMPLEX*16 :: RATE_G
!
 REAL*8, ALLOCATABLE :: SDOT_PIX(:), UDOT_PIX(:), NDOT_PIX(:), GDOT_PIX(:)
 REAL*8, ALLOCATABLE :: SDOT(:,:),   UDOT(:,:),  NDOT(:,:), GDOT(:,:)
 REAL*8, ALLOCATABLE :: IN_TOPO(:,:), TOPO(:,:,:), T(:,:)
 REAL*8, ALLOCATABLE :: LONP(:), LATP(:),  PLM(:,:) 
 REAL*8, ALLOCATABLE :: LONG(:), LATG(:),  PLMG(:,:)
!
 COMPLEX*16, ALLOCATABLE :: SH_S(:,:,:), SH_U(:,:), SH_G(:,:), SH_N(:,:)
 COMPLEX*16, ALLOCATABLE :: SIN_COS_G(:,:), SIN_COS(:,:)
!
 INTEGER, ALLOCATABLE :: ANC(:), MM(:), LL(:), DM(:), OF(:,:), H(:,:)
!
!/////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////
!
! ~~~~~~~ Input files and data ~~~~~~~
!
! =====================================
!  Input files and data for TIDE GAUGES
! ===================================== 
 INTEGER, PARAMETER :: NH_TG = 10, LARGE_N=10**5
 INTEGER NTG, psmsl(large_n), code(large_n)
 CHARACTER*36 tg_name(large_n)
 REAL*8 LAT_TG(large_n), LON_TG(large_n)
 REAL*8 TG_SDOT, TG_UDOT, TG_NDOT, TG_GDOT, X
 COMPLEX*16, ALLOCATABLE :: YTG(:) 
 CHARACTER*72, PARAMETER :: SH_TG_F  = "sh_tg.bin"
!
! ===================================
!  Input files and data for RSL SITEs 
! =================================== 
 INTEGER, PARAMETER :: very_large_integer=10**5
 INTEGER, PARAMETER :: NALFA=52
 INTEGER, PARAMETER :: MAXD =100  
 INTEGER :: NRSL
 CHARACTER*200 ROW 
 CHARACTER*1 ALFAB(NALFA)
 DATA ALFAB/'A','B','C','D','E','F','G','H','I','J','K','L','M', & 
            'N','O','P','Q','R','S','T','U','V','Y','Y','Z','W', & 
            'a','b','c','d','e','f','g','h','i','j','k','l','m', &	     
            'n','o','p','q','r','s','t','u','v','y','y','z','w'/
 REAL*4 AJUNK
 INTEGER, ALLOCATABLE :: CODES(:), NDATA(:)
 REAL*8, ALLOCATABLE :: LATS(:), LONS(:), SLC(:,:), RSL(:,:)
 REAL*8, ALLOCATABLE :: EPOCH(:,:), D_EPOCH(:,:), RSL_OBS(:,:), DRSL_OBS(:,:)
 COMPLEX*16, ALLOCATABLE :: YRSL(:,:) 
 CHARACTER*3 STRING
 CHARACTER*32, ALLOCATABLE :: NAME(:)
 CHARACTER*72, PARAMETER :: SH_RSL_F  = "sh_rsl.bin"
 CHARACTER*72, PARAMETER :: ALL_RSL_CURVES  = "RSL.DAT"
!
! ========================================
!  Input files and data for POLAR MOTION   The dimension is nmodes, in any case. RGS 21 OCT 2018 
! ======================================== 
 INTEGER MP
 CHARACTER*4 JUNK
 REAL*8 AEP, ASP
 REAL*8 PRIMED_A, ROOT_A
 REAL*8, PARAMETER :: RAD2DEG=180D0/PI
 COMPLEX*16, ALLOCATABLE :: AP(:), A_ROOT(:) 
 COMPLEX*16, ALLOCATABLE :: PM(:), PM_DOT(:)
 REAL*8, ALLOCATABLE :: M1(:), M2(:) 
 REAL*8, ALLOCATABLE :: M1_DOT(:), M2_DOT(:)
 REAL*8 POLAR_DISP, POLAR_RATE
 REAL*8 TIME_BP, TIME_CC
!
! Output files
! ------------------------------------------------- 
! ******* Name of fingerprints files and fps stats
! -------------------------------------------------  
 CHARACTER*100, PARAMETER :: SDOT_F="sdot.dat"
 CHARACTER*100, PARAMETER :: UDOT_F="udot.dat"
 CHARACTER*100, PARAMETER :: NDOT_F="ndot.dat"
 CHARACTER*100, PARAMETER :: GDOT_F="gdot.dat"
!
 CHARACTER*100, PARAMETER :: SDOT_F_PIX="sdot.pix"
 CHARACTER*100, PARAMETER :: UDOT_F_PIX="udot.pix"
 CHARACTER*100, PARAMETER :: NDOT_F_PIX="ndot.pix"
 CHARACTER*100, PARAMETER :: GDOT_F_PIX="gdot.pix"

! -----------------------------------------
! ******* Name of file with TG predictions
! -----------------------------------------  
 CHARACTER*100, PARAMETER :: TG_OUT_F="tg.dat"

! ----------------------------------------------
! ******* Name of file with Stokes coefficients
! ----------------------------------------- ---- 
 CHARACTER*100, PARAMETER :: STOKES_F="stokes.dat"
!
! ----------------------------------------------
! ******* Name of file with polar motion and rate of polar motion 
! ----------------------------------------- ---- 
 CHARACTER*20, PARAMETER :: FILE_PM ="m.dat"
 CHARACTER*20, PARAMETER :: FILE_PMD="m.dot"    
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
 write(*,*) ' **** This is the SELEN4 post-processor'
 write(*,*) ' ****'
 write(*,*)
!
 call getarg(1,cfg_f)
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
 call read_data_line(99,buffer)   ;   ROT_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   ICE_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   TOP_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   OFP_F=trim(pt)//trim(adjustl(buffer))
 call read_data_line(99,buffer)   ;   TGS_F=trim(pt)//trim(adjustl(buffer))
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
 ALFA = (180D0/PI)*ACOS(1.-2./FLOAT(NP))
 ST_LMIN=0
 ST_LMAX=LMAX
!
! >>>>>>>>>>> Allocate variable whose size depends on input parameters
!
 ALLOCATE( M1(0:NN), M2(0:NN) )
 ALLOCATE( M1_DOT(0:NN), M2_DOT(0:NN) )
 ALLOCATE( PM(0:NN), PM_DOT(0:NN) )
 ALLOCATE( TBP(0:NN) )
 ALLOCATE( PSI(0:NN+1) )
 ALLOCATE( YTG(JMAX) )
!
!
! ========================== END OF CONFIGURATION SECTION ==========================
! ========================== END OF CONFIGURATION SECTION ==========================
! ========================== END OF CONFIGURATION SECTION ==========================
!
!
!
!
!
!
! -------------------- ! -------------------- ! --------------------
!		       !		      !
!        DATA	       !	DATA	      !        DATA
!		       !		      !
!--------------------- !--------------------- !---------------------
!
 write(*,*) 
 write(*,*) " **** WORKING on BASIC DATA: "
 write(*,*) ' ---- Loading the Tegmark grid data...'
!
 ALLOCATE ( LONP(NP),LATP(NP),ANC(NP) )
 open(1,file=PIX_F) 
   do p=1, np 
      read (1,*) lonp(p), latp(p), anc(p) 
   enddo
 na=anc(np)
 write(*,*) " ---- Tegmark resolution: ", RES
 write(*,*) " ---- number of pixels: ", np
 write(*,*) " ---- number of main pixels: ", na
 close(1)
!
 write(*,*) ' ---- Time stamps ' 
 do n=0, nn 
    CALL TIME_STAMP( nn, n, DELTA, cjunk )
    tbp(n)=cjunk  
 enddo
!
 write(*,*) ' ---- Degres and orders ' 
 ALLOCATE ( MM(JMAX),LL(JMAX),DM(JMAX) )
 do j=1, jmax
    mm(j)=mj(j)
    ll(j)=lj(j)
    dm(j)=2-dom(j)
 enddo
!
 write(*,*) ' ---- Reading pixelized harmonics from file: ', trim(adjustl(SHA_F))
 ALLOCATE ( PLM(JMAX,NA), SIN_COS(0:LMAX,NP) ) 
 open(1,file=SHA_F,form='unformatted') 
 read(1)PLM
 read(1)sin_cos
 close(1) 
!
 write(*,*) ' ---- Loading the ice thickness data'
 ALLOCATE ( H(NP,0:NN+1) ) 
 h(:,:)=0
 open(1,file=ICE_F) 
 do p=1, np  
    read(1,*) j, xj, xj, xj, (h(p,n),n=0,nn+1)
 enddo
 CLOSE(1) 
!
!!!!!!!!!!!!!!!!!!!!!!!!
! ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! --------------------  ! -------------------- ! -------------------- 
!			!		       !
!     TOPOGRAPHY	!     TOPOGRAPHY       !     TOPOGRAPHY
!			!		       !
!---------------------	!--------------------- !---------------------
!
!
 write(*,*) 
 write(*,*) " **** WORKING on TOPOGRAPHY: "
!
 OPEN(11,FILE=F_OUT_TOP,status="unknown",form="unformatted",ACCESS="SEQUENTIAL")
 write(*,*) ' ---- Reading TOPO data from file: ', trim(adjustl(F_OUT_TOP)) 
 ALLOCATE ( TOPO(NP,0:NN+1,0:IEXT_MAX) )
 do i=0, iext_max
    ALLOCATE (  IN_TOPO(NP,0:NN+1) )
    read(11) in_topo	
    topo(:,:,i)=in_topo(:,:)
    DEALLOCATE ( IN_TOPO )
 enddo
 CLOSE(11) 
!
 write(*,*) " ---- Individual files for topography since the LGM..." 
 DO N=0, NN
    FILE_TOPO="topo."//TBP(N)//".dat" 
    if(n==0.or.n==NN)WRITE(*,*) N, TBP(N), " ", FILE_TOPO
    OPEN(101, FILE=FILE_TOPO,STATUS="UNKNOWN") 
       DO P=1, NP 	
          WRITE(101,*) LONP(P), LATP(P), TOPO(P,N,IEXT_MAX)
       ENDDO 
    CLOSE(101) 
 ENDDO
!
!
! -------------------- ! -------------------- ! -------------------- 
!		       !		      !
!       ICE MODEL      !       ICE MODEL      !       ICE MODEL    
!		       !		      !
!--------------------- !--------------------- !---------------------
!
 write(*,*) 
 write(*,*) " **** WORKING on the ICE MODEL: "
!
 write(*,*) " ---- Creating individual files for ice since the LGM..." 
 DO N=0, NN
    FILE_ICE="ice."//TBP(N)//".dat" 
    if(n==0.or.n==NN)    WRITE(*,*) N, TBP(N), " ", FILE_ICE
    OPEN(102, FILE=FILE_ICE,STATUS="UNKNOWN") 
       DO P=1, NP 	
         WRITE(102,*) LONP(P), LATP(P), H(P,N)
       ENDDO 
    CLOSE(102) 
 ENDDO
!
!
! --------------------  ! -------------------- ! -------------------- 
!			!		       !
!    OCEAN FUNCTION	!    OCEAN FUNCTION    !    OCEAN FUNCTION
!			!		       !
!---------------------	!--------------------- !---------------------
!
 write(*,*) 
 write(*,*) " **** WORKING on the OCEAN FUNCTION: "
!
 ALLOCATE ( T(NP ,0:NN) )
 ALLOCATE ( OF(NP,0:NN) )
!
 write(*,*) " ---- Re-building the OF from the TOPO and ICE data" 
 DO N=0, NN
    FILE_ICE = "ice."//TBP(N)//".dat" 
    FILE_TOPO="topo."//TBP(N)//".dat" 
!
    if(n==0.or.n==nn)WRITE(*,*) N, TBP(N), " ", FILE_ICE, FILE_TOPO 
    OPEN(102, FILE=FILE_ICE ,STATUS="UNKNOWN") 
    OPEN(103, FILE=FILE_TOPO,STATUS="UNKNOWN") 
!
! Floating ice 
    FILE_IFL="ice_floating."//TBP(N)//".dat"
    OPEN(104, FILE=FILE_IFL,STATUS="UNKNOWN") 
!
! Grounded ice, below sea level 
    FILE_IGB="ice_grounded_below."//TBP(N)//".dat"
    OPEN(105, FILE=FILE_IGB,STATUS="UNKNOWN") 
!
! Grounded ice, above sea level 
    FILE_IGA="ice_grounded_above."//TBP(N)//".dat"
    OPEN(106, FILE=FILE_IGA,STATUS="UNKNOWN") 
!
! Water
    FILE_OCE="ocean."//TBP(N)//".dat"
    OPEN(107, FILE=FILE_OCE,STATUS="UNKNOWN") 
!
! Ground
    FILE_CON="continent."//TBP(N)//".dat"
    OPEN(108, FILE=FILE_CON,STATUS="UNKNOWN") 
!
!
       DO P=1, NP 	
         READ(102,*) LONP(P), LATP(P), H(P,N)
         READ(103,*) LONP(P), LATP(P), T(P,N)
!
         IF(DBLE(T(P,N))+(RHOI/RHOW)*DBLE(H(P,N))< 0) OF(P,N)=1 
         IF(DBLE(T(P,N))+(RHOI/RHOW)*DBLE(H(P,N))>=0) OF(P,N)=0 
!
! "floating ice":   
         if(of(p,n)==1.and.t(p,n)<0.and.h(p,n)/=0)  write(104,*) lonp(p), latp(p)
!
! "ice grounded below sea level":   
         if(of(p,n)==0.and.t(p,n)<0.and.h(p,n)/=0)  write(105,*) lonp(p), latp(p)
!
! "ice grounded above sea level": 
         if(of(p,n)==0.and.t(p,n)>=0.and.h(p,n)/=0) write(106,*) lonp(p), latp(p)
!
! "ocean":  
         if(of(p,n)==1.and.t(p,n)<0.and.h(p,n)==0)  write(107,*) lonp(p), latp(p)
!
! "continent"  
         if(of(p,n)==0.and.t(p,n)>=0.and.h(p,n)==0) write(108,*) lonp(p), latp(p)
!
       ENDDO 
!
    CLOSE(102) 
    CLOSE(103)
    CLOSE(104)
    CLOSE(105)
    CLOSE(106)
!
 ENDDO
!

! -------------------- ! -------------------- ! --------------------
!		       !		      !
!     FINGERPRINTS     !     FINGERPRINTS     !     FINGERPRINTS
!		       !		      !
!--------------------- !--------------------- !---------------------
!
!
 write(*,*) 
 write(*,*) " **** WORKING on the fingerprints "
!
 write(*,*) ' ---- Reading SHS data from file: ', trim(adjustl(F_OUT_SHS)) 
 write(*,*) ' ---- Reading SHU data from file: ', trim(adjustl(F_OUT_SHU)) 
 write(*,*) ' ---- Reading SHN data from file: ', trim(adjustl(F_OUT_SHN)) 
 write(*,*) ' ---- Reading SHG data from file: ', trim(adjustl(F_OUT_SHG)) 
!
 ALLOCATE ( SH_S(JMAX,0:NN+1,0:IEXT_MAX) )
 ALLOCATE ( SH_U(JMAX,0:NN+1) )
 ALLOCATE ( SH_N(JMAX,0:NN+1) )
 ALLOCATE ( SH_G(JMAX,0:NN+1) )
!
 OPEN (13,FILE=F_OUT_SHS, status="unknown", form="unformatted"); READ (13) SH_S; CLOSE(13)
 OPEN (13,FILE=F_OUT_SHU, status="unknown", form="unformatted"); READ (13) SH_U; CLOSE(13) 
 OPEN (13,FILE=F_OUT_SHN, status="unknown", form="unformatted"); READ (13) SH_N; CLOSE(13) 
 OPEN (13,FILE=F_OUT_SHG, status="unknown", form="unformatted"); READ (13) SH_G; CLOSE(13) 
!
 ALLOCATE( SDOT_PIX(NP) )
 ALLOCATE( UDOT_PIX(NP) )
 ALLOCATE( NDOT_PIX(NP) )
 ALLOCATE( GDOT_PIX(NP) )
!
 DDELTA=2D0*DELTA
!
 write(*,*) ' ---- Computing the fingeprints at the pixels...'  
!
 DO P=1, NP 
!
 if(mod(np,5000)==0) write(*,*) p, "of ", np
!
 sdot_pix(p)=0d0	
 udot_pix(p)=0d0	
 ndot_pix(p)=0d0	
 gdot_pix(p)=0d0	
!
 DO J=1, JMAX 
!
    sdot_pix(p) = sdot_pix(p) + dm(j)*PLM(j,anc(p))* &
         dble(((sh_s(j,nn+1,IEXT_MAX)-sh_s(j,nn-1,IEXT_MAX))/DDELTA)*sin_cos(mm(j),p))
!
    udot_pix(p) = udot_pix(p) + dm(j)*PLM(j,anc(p))* &
         dble(((sh_u(j,nn+1)-sh_u(j,nn-1))/DDELTA)*sin_cos(mm(j),p))
!
    ndot_pix(p) = ndot_pix(p) + dm(j)*PLM(j,anc(p))* &
         dble(((sh_n(j,nn+1)-sh_n(j,nn-1))/DDELTA)*sin_cos(mm(j),p))
!
    gdot_pix(p) = gdot_pix(p) + dm(j)*PLM(j,anc(p))* &
         dble(((sh_g(j,nn+1)-sh_g(j,nn-1))/DDELTA)*sin_cos(mm(j),p))
!
    ENDDO	
 ENDDO	
!
 open(7,file=SDOT_F_PIX,status='unknown')
 open(8,file=UDOT_F_PIX,status='unknown')
 open(9,file=NDOT_F_PIX,status='unknown')
 open(4,file=GDOT_F_PIX,status='unknown')
!
 do p=1,np
    write(7,*) lonp(p), latp(p), sdot_pix(p)       
    write(8,*) lonp(p), latp(p), udot_pix(p)       
    write(9,*) lonp(p), latp(p), ndot_pix(p)       
    write(4,*) lonp(p), latp(p), gdot_pix(p)       
 enddo
!
 close(7)
 close(8)
 close(9)
 close(4) 
!
!
 write(*,*) ''
 write(*,*) ' ---- Computing the present-day OCEAN averages'  
 OPEN(121,file="fps-stats.dat",status="unknown")
 average_s = 0d0 
 average_u = 0d0 
 average_n = 0d0 
 average_g = 0d0 
 nin=0
 do p=1, np 
    if(of(p,nn)==1) then 
       average_s = average_s + sdot_pix(p)
       average_u = average_u + udot_pix(p)
       average_n = average_n + ndot_pix(p)
       average_g = average_g + gdot_pix(p)
       nin=nin+1
    endif
 enddo
 average_s= average_s / float(nin)
 average_u= average_u / float(nin)
 average_n= average_n / float(nin)
 average_g= average_g / float(nin)
!
 write(*,*)   " ---- Average of S-DOT over the oceans (mm/yr): ", average_s
 write(*,*)   " ---- Average of U-DOT over the oceans (mm/yr): ", average_u
 write(*,*)   " ---- Average of N-DOT over the oceans (mm/yr): ", average_n
 write(*,*)   " ---- Average of G-DOT over the oceans (mm/yr): ", average_g
 write(*,  *)"" 
 write(121,*)"" 
 write(121,*) " ---- Average of S-DOT over the oceans (mm/yr): ", average_s
 write(121,*) " ---- Average of U-DOT over the oceans (mm/yr): ", average_u
 write(121,*) " ---- Average of N-DOT over the oceans (mm/yr): ", average_n
 write(121,*) " ---- Average of G-DOT over the oceans (mm/yr): ", average_g
!
 write(*,*) ' ---- Computing the WHOLE EARTH averages'  
 average_s=0d0 
 average_u=0d0 
 average_n=0d0 
 average_g=0d0 
 nin=0
 do p=1, np 
    average_s = average_s + sdot_pix(p)
    average_u = average_u + udot_pix(p)
    average_n = average_n + ndot_pix(p)
    average_g = average_g + gdot_pix(p)
    nin=nin+1
 enddo
 average_s = average_s / float(nin)
 average_u = average_u / float(nin)
 average_n = average_n / float(nin)
 average_g = average_g / float(nin)
 write(*,*)   " ---- Average of S-DOT over the earth (mm/yr): ", average_s
 write(*,*)   " ---- Average of U-DOT over the earth (mm/yr): ", average_u
 write(*,*)   " ---- Average of N-DOT over the earth (mm/yr): ", average_n
 write(*,*)   " ---- Average of G-DOT over the earth (mm/yr): ", average_g
 write(*,*)   "" 
 write(121,*) "" 
 write(121,*) " ---- Average of S-DOT over the earth (mm/yr): ", average_s
 write(121,*) " ---- Average of U-DOT over the earth (mm/yr): ", average_u
 write(121,*) " ---- Average of N-DOT over the earth (mm/yr): ", average_n
 write(121,*) " ---- Average of G-DOT over the earth (mm/yr): ", average_g
 close(121)
!
!
 DEALLOCATE ( SDOT_PIX, UDOT_PIX, NDOT_PIX, GDOT_PIX ) 
 DEALLOCATE ( PLM, SIN_COS, SH_U, SH_S, SH_N, SH_G )
!
!
!
!
!
! -------------------- ! -------------------- ! --------------------
!		       !		      !
!          RSL         !          RSL         !          RSL        
!		       !		      !
!--------------------- !--------------------- !---------------------
!
 write(*,*) 
 write(*,*) " **** WORKING on the Relative Sea level (RSL) sites "
!
 write(*,*) ' ---- Reading data from file: ', trim(adjustl(RSL_F)) 
!
 open(10,file=RSL_F,status='unknown') 
 nrsl=0
 do 33433 j=1, very_large_integer
    read(10,'(a200)',end=41492) row 
       do k=1, 200
          do n=1, NALFA
             if(row(k:k)==alfab(n)) then 
                nrsl=nrsl+1
                goto 33433
             endif 
          enddo
       enddo
33433 continue
41492 close(10)  
 write(*,*) ' ---- Number of sites in the file: ', NRSL
!
  write(*,*) ' ---- Reading code, lat, lon, and names of sites... ' 
  OPEN(1,FILE=RSL_F,STATUS='unknown')
  ALLOCATE ( CODES(NRSL), LATS(NRSL), LONS(NRSL), NAME(NRSL), NDATA(NRSL) )
  ALLOCATE (   EPOCH(NRSL,MAXD), & 
             D_EPOCH(NRSL,MAXD), & 
	     RSL_OBS(NRSL,MAXD), & 
	     DRSL_OBS(NRSL,MAXD)  )  
  i=1
!2 READ (1,6,END=1) codes(i), lats(i), lons(i), ndata(i), name(i) 
2 READ (1,*,END=1) codes(i), lats(i), lons(i), ndata(i), name(i) 
  IF(lons(i)<=0.) lons(i)=360.+lons(i)
  do k=1, ndata(i)
     READ (1 ,*) epoch(i,k), d_epoch(i,k), rsl_obs(i,k), drsl_obs(i,k)
  enddo
  i=i+1
  IF(i<=nrsl) GOTO 2
1 CLOSE(1)
6 FORMAT(1X,I3,1X,F5.1,1X,F6.1,1X,I2,1X,A22)
!
  ALLOCATE ( YRSL(JMAX,NRSL) )
  write(*,*) ' ---- Computing the SHs at the RSL sites coordinates'
  do i=1, nrsl    
     call HARMO ( lmax, lons(i), lats(i), yrsl(:,i) ) 
  enddo
!
 write(*,*) ' ---- Reading SHS data from file: ', trim(adjustl(F_OUT_SHS)) 
 ALLOCATE ( SH_S(JMAX,0:NN+1,0:IEXT_MAX) )
 OPEN (13,FILE=F_OUT_SHS, status="unknown", form="unformatted")
 READ (13) SH_S
 CLOSE(13)
!
  write(*,*) ' ---- Computing sea level change at the RSL sites coordinates'
  ALLOCATE ( SLC(NRSL,0:NN) )
  DO I=1, NRSL
     DO N=0, NN
        SLC(I,N) = 0D0 
           DO J=1, JMAX
              SLC(I,N) = SLC(I,N) + DM(J)*DBLE(SH_S(J,N,IEXT_MAX)*YRSL(J,I)) 
	   ENDDO
     ENDDO
  ENDDO
!
 write(*,*) ' ---- Computing RSL at the RSL sites coordinates...'
 write(*,*) ' ---- ... and dumping all curves in one file (rsl.dat)'
 ALLOCATE( RSL(NRSL, 0:NN) )
 OPEN(33,FILE=ALL_RSL_CURVES,STATUS="UNKNOWN")
 DO I=1, NRSL 
 write(33,*) codes(i), name(i)
    DO K=0, NN
       RSL(I,K)= -(SLC(I,NN)-SLC(I,NN-K))
       WRITE(33,*) FLOAT(K)*DELTA*1000.0, RSL(I,K)	 
    ENDDO
 ENDDO 
 CLOSE(33)	
!
 write(*,*) ' ---- Computing synthetic RSL curves at the sites'
 write(*,*) ' ---- (one RSL curve in one file)'
 DO I=1, NRSL
    !OPEN  (81,FILE='junk.dat',STATUS='unknown'); WRITE (81,'(i3)') codes(i); CLOSE(81)
    !OPEN  (81,FILE='junk.dat',STATUS='unknown'); READ  (81,'(a3)') string ;  CLOSE(81)
    WRITE(string,'(i3)') codes(i)
    OPEN  (91,FILE='rslp'//'-'//string//'.dat',STATUS='unknown')
    WRITE (91,'(a1,1x,a22,a5,a3,1x,a8,2(1x,f10.5))') &
    '#', NAME(i), 'code=', string, ' lon/lat:', lons(i), lats(i)	      
    WRITE (91,'(a2,a17)') '# ', 'years BP, RSL (m)'
       DO K=0, NN
          WRITE(91,"(F10.2,1X,F11.3)") 500.*k, RSL(I,K)	 
       ENDDO
 ENDDO
!
 write(*,*) ' ---- Files with RSL data at individual sites...'
 DO I=1, NRSL
    !OPEN (81,FILE='junk.dat',STATUS='unknown'); WRITE (81,'(i3)') codes(i); CLOSE (81)
    !OPEN (81,FILE='junk.dat',STATUS='unknown'); READ  (81,'(a3)') string;   CLOSE (81)
    WRITE(string,'(i3)') codes(i)    
    OPEN (91,FILE='rsld-'//string//'.dat',STATUS='unknown')
    WRITE(91,'(a1,1x,a22,1x,a3,a11,2(1x,f10.5))') &
         '#', NAME(i), string, '  lon/lat= ', lons(i), lats(i)
    WRITE(91,'(a1,a29)') '#','  kyrs BP   RSL (m)  DRSL (m)'
    DO K=1, NDATA(I)
       WRITE (91 ,'(F10.2,2(F10.5))') EPOCH(I,K), RSL_OBS(I,K), DRSL_OBS(I,K)
    ENDDO
    CLOSE(91)
ENDDO
!
 DEALLOCATE ( SH_S )
!
!
!
!
! -------------------- ! -------------------- ! --------------------
!		       !		      !
!     TIDE GAUGES      !     TIDE GAUGES      !     TIDE GAUGES
!		       !		      !
!--------------------- !--------------------- !---------------------
! GS July 10 and 11, 2016
! Rev August 2, same year 
!
 write(*,*) 
 write(*,*) " **** WORKING on TIDE GAUGES: "
!
 write(*,*) " ---- Counting the number of TGs from file: "
 write(*,*) "      ", trim(adjustl(TGS_F))
 open(1,file=TGS_F,status="unknown") 
 do i=1, nh_tg 
    read(1,"(A4)")cjunk
 enddo
 ntg=0
 do i=1, large_n 
    read(1,"(A4)",end=111)cjunk  
    ntg=ntg+1
 enddo
 111 continue 
 write(*,*) " ---- Number of TGs: ", ntg
 close(1) 
!
 write(*,*) " ---- Computing the SHs at the TGs locations..." 
 open(1,file=TGS_F,status="unknown") 
 do i=1, nh_tg 
    read(1,"(A4)")cjunk
 enddo
 open(17,file=SH_TG_F,form="unformatted",access="SEQUENTIAL")
 do i=1, ntg 
    read(1,*) lat_tg(i), lon_tg(i), psmsl(i), code(i), x, x, x, tg_name(i)  
    call harmo(lmax, lon_tg(i), lat_tg(i), ytg) 
    write(17) ytg   
 enddo 
 close(17)
 close(1)
!
!
 write(*,*) " ---- Computing  Sdot, Udot, Ndot and Gdot at TGs locations" 
!
 open(17,file=SH_TG_F,form="unformatted",access="SEQUENTIAL")
!
 write(*,*) ' ---- Reading SHS data from file: ', trim(adjustl(F_OUT_SHS)) 
 write(*,*) ' ---- Reading SHU data from file: ', trim(adjustl(F_OUT_SHU)) 
 write(*,*) ' ---- Reading SHN data from file: ', trim(adjustl(F_OUT_SHN)) 
 write(*,*) ' ---- Reading SHG data from file: ', trim(adjustl(F_OUT_SHG)) 
 ALLOCATE ( SH_S(JMAX,0:NN+1,0:IEXT_MAX) )
 ALLOCATE ( SH_U(JMAX,0:NN+1) )
 ALLOCATE ( SH_N(JMAX,0:NN+1) )
 ALLOCATE ( SH_G(JMAX,0:NN+1) )
 OPEN (13,FILE=F_OUT_SHS, status="unknown", form="unformatted"); READ (13) SH_S; CLOSE(13)
 OPEN (13,FILE=F_OUT_SHU, status="unknown", form="unformatted"); READ (13) SH_U; CLOSE(13)
 OPEN (13,FILE=F_OUT_SHN, status="unknown", form="unformatted"); READ (13) SH_N; CLOSE(13)
 OPEN (13,FILE=F_OUT_SHG, status="unknown", form="unformatted"); READ (13) SH_G; CLOSE(13)
!
 open(34,file=TG_OUT_F,status="unknown")
 write(34,*) "    Lon (deg)      Lat (deg)        S-dot          U-dot      &
         N-dot         G-dot        PSMSL ID and name" 
 DDELTA=2D0*DELTA
 do i=1, ntg 
    if(i==1)  write(*,*) " ---- ", "TG ", i, " of ", ntg
    if(i==ntg)write(*,*) " ---- ", "TG ", i, " of ", ntg
    read(17) ytg
    TG_SDOT = 0D0
    TG_UDOT = 0D0 
    TG_NDOT = 0D0 
    TG_GDOT = 0D0 
       do j=1, jmax
          TG_SDOT = TG_SDOT + &
	            DM(J)*REAL(((SH_S(J,NN+1,IEXT_MAX)-SH_S(J,NN-1,IEXT_MAX))/DDELTA)*YTG(J)) 
          TG_UDOT = TG_UDOT + & 
	            DM(J)*REAL(((SH_U(J,NN+1)-SH_U(J,NN-1))/DDELTA)*YTG(J)) 
          TG_NDOT = TG_NDOT + & 
	            DM(J)*REAL(((SH_N(J,NN+1)-SH_N(J,NN-1))/DDELTA)*YTG(J))    
          TG_GDOT = TG_GDOT + & 
	            DM(J)*REAL(((SH_G(J,NN+1)-SH_G(J,NN-1))/DDELTA)*YTG(J))    
       enddo
    write(34,"(6(F14.8,1X),I6,1X,A36)") & 
    lon_tg(i), lat_tg(i), TG_SDOT, TG_UDOT, TG_NDOT, TG_GDOT, psmsl(i), tg_name(i)  
 enddo
 close(17)
 close(34)
 write(*,*) ' ---- For the predictions at TGs, see file: ', trim(adjustl(TG_OUT_F)) 
!
 DEALLOCATE ( SH_S, SH_U, SH_N, SH_G ) 
!
!
!

!
! ----------------------- !! ----------------------- !! ----------------------- ! 
!		          !!			     !! 			!
!  STOKES's Coefficients  !!  STOKES's Coefficients  !!  STOKES's Coefficients  ! 
!                         !!			     !! 			!
!------------------------ !!------------------------ !!------------------------ !
! GS July 30, 2016
!
 write(*,*) 
 write(*,*) " **** WORKING on the STOKEs COEFFICIENTS: "
!
 write(*,*) ' ---- Reading SHG data from file: ', trim(adjustl(F_OUT_SHU)) 
 ALLOCATE ( SH_G(JMAX,0:NN+1) )
 OPEN (13,FILE=F_OUT_SHG, status="unknown", form="unformatted")
 READ (13) SH_G
 CLOSE(13)
!
 DDELTA=2D0*DELTA
!
!
 OPEN(17,file=STOKES_F,status="unknown")
!
 DO L= ST_LMIN, ST_LMAX  
    DO M=0, L 
!
       RATE_G = (SH_G(J_INDEX(L,M),NN+1)-SH_G(J_INDEX(L,M),NN-1))/DDELTA
!
! ---- now in (yrs)^(-1)
       RATE_G = RATE_G/(ERADIUS*1E3)  
!
! ---- now in units of 10^-11               	       
       RATE_G = 1.E11*RATE_G
!  
! ---- GRACE-compliant normalization
       C_STOK = +   REAL (RATE_G) * ((-1)**M)*SQRT(2.-DOM(J_INDEX(L,M)))
       S_STOK = - AIMAG (RATE_G) * ((-1)**M)*SQRT(2.-DOM(J_INDEX(L,M)))
!
       WRITE(17,'(4X,I9,6X,I4,5X,I4,8X,4(E12.6,6X))') & 
       J_INDEX(L,M), L, M, C_STOK, S_STOK, SQRT(C_STOK**2+S_STOK**2)
!
    ENDDO
 ENDDO   
!
 CLOSE(17)
!
 write(*,*) ' ---- For the STOKES COEFFICIENTS, see file: ', trim(adjustl(STOKES_F)) 
!
 DEALLOCATE ( SH_G )
!
!
!
!
! --------------------  ! -------------------- ! -------------------- 
!			!		       !
!     POLAR MOTION	!     POLAR MOTION     !    POLAR MOTION
!			!		       !
!---------------------	!--------------------- !---------------------
!
 IF(IROT/=0) THEN
!
 write(*,*)
!
 write(*,*) ' ---- Reading PSI data from file: ', trim(adjustl(F_OUT_TPW)) 
 open(17,file=F_OUT_TPW,status="unknown")
 do n=0, nn+1
    read(17,*) k, psi(n)
!    PSI(N) = DCMPLX(PSI_X,PSI_Y)
 enddo
 close(17) 
!
! ---- Traditional theory of Earth rotation
 IF(IROT==1) MP=NMODES-1 
 IF(IROT==1) write(*,*) ' **** ANALYSIS OF POLAR MOTION (TRADITIONAL THEORY) ' 
!
! ---- Revised theory of Earth rotation
 IF(IROT==2) MP=NMODES
 IF(IROT==2) write(*,*) ' **** ANALYSIS OF POLAR MOTION (REVISED THEORY) ' 
! 
 ALLOCATE( AP(MP), A_ROOT(MP) )
!
    write(*,*) ' ---- Number of rotational modes: ', MP 
    write(*,*) ' ---- Reading PMTF data from file: ', trim(adjustl(PMT_F)) 
    OPEN (1,file=PMT_F,status='unknown') 
    READ (1,"(A10)") JUNK 
    READ (1,"(A10)") JUNK 
    READ (1,"(A10)") JUNK 
    READ (1,"(A10)") JUNK 
    READ (1,"(A10)") JUNK 
    DO I=1, MP 
       READ (1,"(A10)") JUNK 
    ENDDO
    READ (1,"(A10)") JUNK 
    READ (1,"(A10)") JUNK 
    READ (1,*) AEP
    READ (1,*) ASP          
    DO I=1, MP 
       READ (1,*) J, PRIMED_A, ROOT_A 
       AP(I)     = DCMPLX(PRIMED_A , 0.0)
       A_ROOT(I) = DCMPLX(ROOT_A   , 0.0)   	
    ENDDO
    CLOSE(1)  
!
!
 write(*,*) ' ---- Computing the polar motion ' 
 do n=0, nn 
       pm(n)=(0d0,0d0) 
          if(n>=1) then 
             do k=0, n-1 
                pm(n) = pm(n) + (psi(k+1)-psi(k))*(aep + asp*float(n-k)*delta)   
                   do j=1, mp 
                      pm(n) = pm(n) + (psi(k+1)-psi(k))* &
		      (ap(j)/a_root(j))*(exp(a_root(j)*float(n-k)*delta)-1d0)  
                   enddo
             enddo
         endif
 enddo
 write(*,*) ' ---- Computing the rate of polar motion ' 
 do n=0, nn 
       pm_dot(n)=(0d0,0d0) 
          if(n>=1) then 
             do k=0, n-1 
                pm_dot(n) = pm_dot(n) + (psi(k+1)-psi(k))*asp    
                   do j=1, mp 
                      pm_dot(n) = pm_dot(n) + (psi(k+1)-psi(k))* &
		      ap(j)*exp(a_root(j)*float(n-k)*delta)  
                   enddo
             enddo
          endif 
 enddo	
!
!
 write(*,*) ' **** Polar displacement and its rate ' 
!
 write(*,*) ' ---- Output file for displacement: ', trim(adjustl(FILE_PM))  
 OPEN (55,file=FILE_PM, status='unknown')   
 CALL PRINT_HEADER(55)
!
 write(*,*) ' ---- Output file for the rate of displacement: ', trim(adjustl(FILE_PMD))  
 OPEN (57,file=FILE_PMD,status='unknown')   
 CALL PRINT_HEADER(57)
!
 do n=0, nn 
!
    m1(n) = real  (pm(n))
    m2(n) = aimag (pm(n))
!
    m1_dot(n) = real (pm_dot(n))
    m2_dot(n) = aimag(pm_dot(n)) 
!
! Polar displacement (deg)  
    polar_disp = sqrt(m1(n)**2+m2(n)**2)*rad2deg
!
! Rate of polar displacement (deg/Ma)    
    polar_rate = ((m1(n)*m1_dot(n) + m2(n)*m2_dot(n))/sqrt(m1(n)**2+m2(n)**2))*rad2deg*1000.
!
! time BP 
    time_bp=float(nn-n)*delta
    time_cc=float(n)*delta
!
    write(55,'(2(4x, f8.4), 10(1x,e14.5))') & 
    time_cc, & 
    time_bp,  & 
    m1(n)*rad2deg, & 
    m2(n)*rad2deg, & 
    atan2(m2(n),   m1(n)   )*rad2deg,  &
    polar_disp 
!        
    write(57,'(2(4x, f8.4), 10(1x,e14.5))') & 
    time_cc, & 
    time_bp,  & 	
    m1_dot(n)*rad2deg*1000., & 
    m2_dot(n)*rad2deg*1000., & 
    atan2(m2_dot(n),m1_dot(n))*rad2deg,  & 
    polar_rate  

enddo
!
 END IF         ! On IROT/=0
!
!!!!!!!!!!!!!!!!!!!!! 
 6969 CONTINUE
!!!!!!!!!!!!!!!!!!!!!
!
 write(*,*) 
 write(*,*) ' **** DONE '
 write(*,*) 
!
! ----------------
 STOP
 END program PPROC
! ----------------
!
!
!
!
!
!
!
!
 SUBROUTINE TIME_STAMP ( nn, n, delta, etichetta ) 
 implicit NONE 
 integer n, nn, nbp
 character*4 etichetta
 real(8)  tbp, delta
!  
     nbp=nn-n 
     tbp=dble(nbp)*delta
!
     if( tbp.lt.10. ) then
         write(etichetta,'(a1,f3.1)')  '0', tbp
     else
         write(etichetta,'(f4.1)') tbp
     endif
!
!      if(0.le.nbp.and.nbp.le.19) then 
!        !open(3,file="junk.dat",status="unknown")        
!        if(mod(nbp,2)==0) then 
!      	!write(3,"(A1,I1,A2)")"0", nbp/2,".0"
!        !close(3)
!	write(etichetta,"(A1,I1,A2)")"0", nbp/2,".0"        
!        endif
!        if(mod(nbp,2)/=0) then 
!     	!write(3,"(A1,I1,A2)")"0", nbp/2,".5"
!        !close(3)
!     	write(etichetta,"(A1,I1,A2)")"0", nbp/2,".5"
!        endif
!      elseif(nbp.ge.10) then 
!        !open(3,file="junk.dat",status="unknown") 
!        if(mod(nbp,2)==0) then 
!       	!write(3,"(I2,A2)")nbp/2,".0"
!        !close(3)
!        write(etichetta,"(I2,A2)")nbp/2,".0"
!	endif
!        if(mod(nbp,2)/=0) then 
!      	!write(3,"(I2,A2)")nbp/2,".5"
!        !close(3)
!	write(etichetta,"(I2,A2)")nbp/2,".5"
!        endif
!      endif       
!        !open(3,file="junk.dat",status="unknown") 
!	!read(3,"(A4)")etichetta 
!	!close(3)        
!
 End SUBROUTINE TIME_STAMP
!
!
!
	FUNCTION SINDD(ALFA)
	implicit NONE
!
! --- Assuming that ALFA is given in DEGREES, it returns the SINE of ALFA
!     *********************** GS and FC 11-02-2009 **********************
!	
	REAL*8,  PARAMETER :: PI=3.14159265358979323840D0 
	REAL*8 SINDD, ALFA, ALFARAD              
	ALFARAD=ALFA*PI/180.	
	SINDD=SIN(ALFARAD)	
	END FUNCTION SINDD 
!
!
!
	FUNCTION COSDD(ALFA)
	implicit NONE
!
! --- Assuming that ALFA is given in DEGREES, it returns the COSINE of ALFA
!     ************************ GS and FC 11-02-2009 ***********************
!	
	REAL*8,  PARAMETER :: PI=3.14159265358979323840D0 
	REAL*8 COSDD, ALFA, ALFARAD               
	ALFARAD=ALFA*PI/180.	
	COSDD=COS(ALFARAD)	
	END FUNCTION COSDD
!
!
!
!####################
 subroutine RJ(i,j)
!####################
!
! Reads "i" junk lines from unit "j" 
! GS 2011
!
 implicit NONE 
 integer i, j, l  
 character*50 junk
! 
 do l=1, i 
   read(j,'(a50)')junk 
 enddo
!
 end subroutine RJ
!
!
!
!
!
!
     SUBROUTINE PRINT_HEADER(II)
     IMPLICIT NONE
     INTEGER II
     CHARACTER*20 DATE, TIMC
!
     call DATE_AND_TIME (date,timc)
!
     IF(II==55) THEN 
!       
! ---Header for m.dat
 Write(55,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
 Write(55,*) '# File <<m.dat>>, created by program LIOUVILLE.F90 on ', & 
              date(1:4), '.', date(5:6), '.', date(7:8), & 
	      '      time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(55,*) "# Note: the Chandler Wobble is NOT included"
 Write(55,*) "# " 
 Write(55,*) "#     time      time BP        m1             m2           Arg (m)        Mod (m)"	
 Write(55,*) "#      ka         ka          deg            deg             deg            deg  "	  
 Write(55,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
!     
     ENDIF
!
     IF(II==57) THEN 

! ---Header for m.dot
 Write(57,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
 Write(57,*) '# File <<m.dot>>, created by program LIOUVILLE.F90 on ', & 
              date(1:4), '.', date(5:6), '.', date(7:8), & 
              '      time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(57,*) "# Note: the Chandler Wobble is NOT included"
 Write(57,*) "# " 
 Write(57,*) "#     time      time BP       m1_dot         m2_dot      Arg (m_dot)    Mod (m_dot)"   
 Write(57,*) "#      ka         ka          deg/Ma         deg/Ma          deg          deg/Ma   "      
 Write(57,*) "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " 
!      
     ENDIF
!   
     END SUBROUTINE PRINT_HEADER 
!
!
!
!

