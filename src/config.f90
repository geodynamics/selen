! 
!*****************************
! This is program "config.f90" 
!*****************************
!
! Last modified GS & FC July 8, 2008  = INTEL PORT = V 2.6 
! Also touched a number of times during July, 2008-
! Re-touched on February 6 2009 for version 2.7 (on the way back from Liverpool) 
! Reviewed GS & FC July 2009 - v3.1-  "Varying coastlines" (reprise!) 
! Reviewed GS & FC July 2009 - v3.1-  "ALMA coupling"
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! Also revised April 2010 by GS for completing the g95 version in 3.1 
! Updated on April 13, 2010 for degree 1 implementation
! Updated on May 16, 2010 (eb* files move to ALMA love number repository)
! Originally designed for v 3.1 - Now ported to 2.9 
! Modified GS April 2011 for the GMT switch 
! Modified DM April 2011 for case sensitiveness + gfortran
! Modified DM June 2011 for parallel execution
! Modified DM & GS August 2011 - 3D velocity at sites and variable time step (29.8)
! Modified DM August 2011 - Stokes coefficients (29.8-1)
! Modified DM August 2011 - Some optimizations in the compilation process (29.8-2)
! Modified DM Jan 2012 - Alaska ice model
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
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! --------------------------------------------------------------------------- 
!
!   Program "config.f90" is the interface between the configuration file 
!   "config.dat" and the bash script "selen.sh" which drives the execution 
!   of SELEN. During execution, "config.f90" writes a summary of the SELEN 
!   settings on monitor. The outputs of "config.f90" are: 
!
!   - "data.inc", a Fortran script with general declarations, which serves 
!                 as input for most of the SELEN Fortran programs, 
! 
!   - "selen.sh", a bash shell script that compiles and executes the SELEN 
!                 programs and organizes their outputs according to the 
!		  settings in "config.dat", 
!
!   - several GMT scripts based on the options in "config.dat". 
!
!   Last change GS & FC February 6 2009 
!
! --------------------------------------------------------------------------- 
!
!
!
PROGRAM CONFIG
IMPLICIT NONE 
!
!
! =======================================================
!                 --- Declarations ---
! =======================================================
!
! -Some bounds  
      INTEGER, PARAMETER :: NALFA=52, MAX_LAYERS=100, NIMP=100, & 
  	                    VERY_LARGE_INTEGER=10000, LARGE_INTEGER=200, & 
			    RES_MIN=12, NRG=10, NMAX_3D_REGIONS=10, RES_MAX=48, & 
			    NP_MAX=2*RES_MAX*(RES_MAX-1)*20+12
			     	
!			   
! -Do-loop indices & other integers 
      INTEGER I, J, K, N, II, NP		
!
! -Some counters 
      INTEGER HEADER_LINES, ICE_FLAG, IMUL, NLINES, NNV, & 
      	      NDEGREE, NRESOLUTION, NOUT, NRSL, NRSLC, NICE, & 
	      NTIDEGAUGES, NRSL_DATA, LEN_ICE, LEN_RSL, & 
	      LEN_VISCO, LEN_TGAUGES, LEN_RSLC, LEN_GEOD, NGEOD, & 
	      LEN_3D_REGIONS, N_3D_REGIONS, NP_3D_REGIONS(NMAX_3D_REGIONS), & 
	      LEN_TRED_REGIONS_NAME(NMAX_3D_REGIONS), & 
	      LEN_FILE_TOPO, LEN_FILE_PXTOPO, ITER_CI, NNINC, NNTHREAD, &
	      NNTASK, LEN_TMP
!	    
! -A few real variables 	  
      REAL(8) :: DELTA = 1.0	     ! Default ice step is 1kyr
      REAL*8 LATS, LONS, RJUNK, RSL_DATUM, D_RSLDATUM, & 
             TIME_RSL, TIME_RSL1, TIME_RSL2, TIME_BPC
!
! -A string useful for compilation commands
      CHARACTER*65 :: Compile
!
! -Sequential and parallel compilation commands
      CHARACTER*65 :: CompileSeq
      CHARACTER*65 :: CompileSmp, CompileMpi
!
! -MPI execution prefix
      CHARACTER*100 :: RunCmd
!
!      CHARACTER*55, PARAMETER :: &
!      SeqCompile = "gfortran -m64 -w    "
!      CHARACTER*55, PARAMETER :: &
!      ParCompile = "gfortran -m64 -w -fopenmp -Wl,-stack_size,0x10000000   "      
!
!      CHARACTER*12, PARAMETER :: & 
!      SeqCompile = "ifort       "
!      CHARACTER*22, PARAMETER :: &
!      ParCompile = "ifort        -openmp  "
!
	 
	 
	 
!
! -The alphabet 
      CHARACTER*1 ALFA(NALFA)
      DATA ALFA/'A','B','C','D','E','F','G','H','I','J','K','L','M', & 
                'N','O','P','Q','R','S','T','U','V','Y','Y','Z','W', & 
	        'a','b','c','d','e','f','g','h','i','j','k','l','m', & 		   	     
		'n','o','p','q','r','s','t','u','v','y','y','z','w'/
!
! -Strings for basic SELEN settings variables
      CHARACTER*10  RESOLUTION, DEGREE, TITLICE, DEGREE_ST_MIN, DEGREE_ST_MAX  
      CHARACTER*10  TIME_BPCC, MIN_RSLC, MAX_RSLC, RSL_INT
      CHARACTER*10  LONMINC, LONMAXC, LATMINC, LATMAXC
      CHARACTER*10  RADIUS_ZOF
      CHARACTER*10  ICE_DELTA 
      CHARACTER*1   ITER, MODE
      CHARACTER*10  ITER_C  
      CHARACTER*3   NINC, NV, CODE
      CHARACTER*3   MINIMUM_PERIOD  
      CHARACTER*50  NAME_OF_REGION
      CHARACTER*100 WDIR
      CHARACTER*10  LATSC10, LONSC10 
      CHARACTER*3   NTHREAD, NTASK
!
! -File and GMT scripts names 
      CHARACTER*20 & 
      SH_FILE, 	      ICE_FILE,            ICE_FILENAME,      SHAPE_FILE,        &
      SHOF_FILE,      FILE_GMT,                      & 
      FILE1_GMT,      FILE2_GMT,           FILE_CAP, & 
      RSLC_LONLAT_FILE, &       
      FILE_TOPO,        FILE_PXTOPO
      CHARACTER*30 RSL_FILE,     RSL_DATABASE,       & 
      		   FILE_REGION,  FILE_REGION_LONLAT, & 
                   TGAUGES_FILE, TGAUGES_DATABASE,   & 
		   FILE_3D,      GEODETIC_DATABASE,  & 
		   VISCO_FILE, 			     &
      FILE_3D_REGIONS,  TRED_REGIONS_DATABASE,       & 
      NAME_3D_REGIONS(NMAX_3D_REGIONS), PIX_3D_FILENAMES(NMAX_3D_REGIONS), &
      TRED_REGIONS_NAME(NMAX_3D_REGIONS)             	 
      CHARACTER*100 VISCO_FILENAME, ALMA_FILE_MODE_2
      CHARACTER*100 ALMA_FILE_PWA_H, ALMA_FILE_PWA_L, ALMA_FILE_PWA_K
      CHARACTER*100 ALMA_LOGFILE_PWA 
      CHARACTER*100 SHORT_VISCO_FILENAME
      CHARACTER*100 TMP_PREFIX

!                   	
!
! -PAth to ICE-MODELS 
      CHARACTER*100 XFILE
!
! -Rheological parameters 
      CHARACTER*100 VSTRING, LITHO, VSC(MAX_LAYERS)
!
! -Input text strings 
      CHARACTER*100 SS(NIMP)
      CHARACTER*200 LINE, LINEP, ROW, JUNK
      CHARACTER*10 CJUNK, AJUNK, SS2(2) 
      CHARACTER*1 RSL_DATABASE_FORMAT 
!
! -Options   
      CHARACTER*1 & 
      OPTION_SH,       OPTION_SF,     OPTION_WI,      OPTION_PX,           OPTION_OF,    & 
      OPTION_RI,       OPTION_LN,     OPTION_GM,      OPTION_TGPLOT,       OPTION_ST,    &
      OPTION_CL,       OPTION_LB,     OPTION_US,      OPTION_OH,           OPTION_OR,    &
      OPTION_RSLDB,    OPTION_RSLP,   OPTION_RSL,     OPTION_RSLZ,         OPTION_RSLMF, & 
      OPTION_RSLSCA,   OPTION_TG,     OPTION_TGSCA,   OPTION_BIN,          OPTION_ESL,   & 
      OPTION_RSLC,     OPTION_TGA,    OPTION_RSLA,    OPTION_RSLTAB,       OPTION_OFDV,  & 
      OPTION_RM(0:NRG),OPTION_ROF,    OPTION_3D,      OPTION_3D_REGIONS,   OPTION_TOPO,  & 
      OPTION_PXTOPO,   OPTION_NM,     OPTION_PW,      OPTION_LOVE_NUMBERS, OPTION_PTMAP, & 
      OPTION_PWA,      OPTION_DEG1,   OPTION_GMT,     OPTION_SYS,   &
      OPTION_MPI,      OPTION_OMP,    OPTION_PURGE,   OPTION_NLS
!
! -More options   
      CHARACTER*2      OPTION_RFRAME            
!
! -Repository name 
      CHARACTER*12 depot 
!
! -Other character constants 
      CHARACTER*40 for_argument
! 
! -Repository label       
      CHARACTER*4 RUN
      CHARACTER*6 RUNP
!
! -Date and time...
      CHARACTER*20 date, timc   
!
! -Some logical switches 
      LOGICAL LEX, LEXX(4) 
!      
! -Log file
      CHARACTER*30, parameter :: logfile="selen.log"
!
! -Water and ice density 
!  REAL*8,  PARAMETER :: RHO_EARTH=5511.57
      REAL*8,  PARAMETER :: RHO_WATER=1000.00 
      REAL*8,  PARAMETER :: RHO_ICE  = 931.00
!
!
! ************************************************************************
! ************************************************************************
! ************************************************************************
! ************************************************************************
! ************************************************************************
!                           Execution starts here!
! ************************************************************************
! ************************************************************************
! ************************************************************************
! ************************************************************************
! ************************************************************************
!
!
      Open(88,file=logfile,status='unknown') 
      call DATE_AND_TIME (date,timc)      
      Write(88,*) " "
      Write(88,*) '# This is file "selen.log", created by "config.f90" on ', & 
             date(1:4), '.', date(5:6), '.', date(7:8), ' ', & 
	     timc(1:2), '.', timc(3:4), '.', timc(5:6) 
      Write(88,*) "  "
!
!
! ====================================
!
! 	Part #1 : 
! 	1/1: General settings 
! 	2/1: Output settings
! 	3/1: Organizing the outputs
!
! ====================================
!
!
! 	# Part 1/1: Reading the general settings "config.dat"   
!
!
  open(1,file='config.dat',status='unknown')
!
!
  Write(88,*) ""
  Write(88,*) "+----------------------------------------------+"
  Write(88,*) "|              Settings of SELEN               |"
  Write(88,*) "+----------------------------------------------+"
!
!
  DO 100 I=1,LARGE_INTEGER
!
  read(1,'(a200)',end=444) line 
!
!
!
! ###### WORKING DIRECTORY ######  
!
  IF(line(1:3)=='100') THEN 
	call scan_string (line, 1, ss, nout)
	wdir=ss(1) 
	Write(88,*) "The working directory is ", wdir 
	Open (55,file='working-directory.txt',status='unknown') 
	Write(55,'(a100)') wdir 
	close(55) 
  ENDIF
!
! ###### PURGING OPTION ######  
!
  IF(line(1:3)=='110') THEN 
	call scan_string (line, 1, ss, nout)
	option_purge=ss(1) 
	if(option_purge=='y')Write(88,*) "The wdir will be purged before and after execution..." 
	if(option_purge=='n')Write(88,*) "The wdir will NOT be purged before and after execution..." 
  ENDIF  
!
!
! ###### MPI / OpenMP setup ######  
!  
  IF(line(1:3)=='120') THEN 
	call scan_string (line, 2, ss, nout)
	option_omp=ss(1) 
	nthread   =ss(2)
        Call CHAR3_2_INT(nthread, nnthread)
	if(option_omp=='y') Write(88,*) "OpenMP is enabled with ",nthread," threads/task"
	If(option_omp/='y'.and.option_omp/='n') then 
		Write(* ,*) "For the OpenMP switch, only y/n are valid options"
		Write(88,*) "For the OpenMP switch, only y/n are valid options"
		Call stop_config 
		Stop	
	Endif
  ENDIF  
!
!
! ###### Compiler setup ###### 
!
if(line(1:3)=="121") THEN
    option_mpi='n'
    call scan_string (line, 1, ss, nout)
    option_sys=ss(1)
    ! --------------------------------------------------------------------------
    If( option_sys == "1" ) then
       Write(88,*) "Configuring SELEN for gfortran on Mac OS X"
       CompileSeq = "gfortran -w -m64 -DGNU "                    
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
         CompileSmp = "gfortran -m64 -w -fopenmp -DGNU -Wl,-stack_size,0x10000000 "   
         CompileMpi = CompileSmp
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
         CompileMpi = "mpif90 -m64 -w -fopenmp -DGNU -DMPI "   
         CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
         CompileMpi = "mpif90 -m64 -w -fopenmp -DGNU -DMPI -Wl,-stack_size,0x10000000 " 
         CompileSmp = "gfortran -m64 -w -fopenmp -DGNU -Wl,-stack_size,0x10000000 "          
       endif       
    ElseIf( option_sys == "2" ) then
       Write(88,*) "Configuring SELEN for Intel ifort on Mac OS X"
       CompileSeq = "ifort "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
         CompileSmp = "ifort -openmp -Wl,-stack_size,0x10000000 "
         CompileMpi = CompileSmp
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
         CompileMpi = "mpif90 -DMPI "    
         CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
         CompileMpi = "mpif90 -DMPI -openmp -Wl,-stack_size,0x10000000 "   
         CompileSmp = "ifort -openmp -Wl,-stack_size,0x10000000 "
       endif       
    ElseIf( option_sys == "3" ) then
       Write(88,*) "Configuring SELEN for gfortran on Linux"
       CompileSeq = "gfortran -w -DGNU "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
         CompileSmp = "gfortran -w -fopenmp -DGNU "   
	 CompileMpi = CompileSmp 
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
         CompileMpi = "mpif90 -w -DMPI -DGNU "  
	 CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
         CompileMpi = "mpif90 -w -DMPI -DGNU -fopenmp "     
         CompileSmp = "gfortran -w -fopenmp -DGNU "   	 
       endif       
    ElseIf( option_sys == "4" ) then
       Write(88,*) "Configuring SELEN for Intel ifort on Linux"
       CompileSeq = "ifort "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
         CompileSmp = "ifort -openmp "    	   
	 CompileMpi = CompileSmp
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
         CompileMpi = "mpif90 -DMPI "
	 CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
         CompileMpi = "mpif90 -openmp -DMPI "    
         CompileSmp = "ifort -openmp "    	   	   
       endif   
    ElseIf( option_sys == "5" ) then
       Write(88,*) "Configuring SELEN for Intel ifort on Linux"
       CompileSeq = "ifort "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
         CompileSmp = "ifort -openmp "   
	 CompileMpi = CompileSmp  	   
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
         CompileMpi = "mpiifort -DMPI "
	 CompileSmp = CompileSeq      
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
         CompileMpi = "mpiifort -openmp -DMPI "
         CompileSmp = "ifort -openmp "   	       
       endif   
    ElseIf( option_sys == "6" ) then
       Write(88,*) "Configuring SELEN for IBM XLF on SP6"
       CompileSeq = "xlf90_r -WF,-DXLF "      
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
          CompileSmp = "xlf90_r -qsmp=omp -qtune=pwr6 -qarch=pwr6 -WF,-DXLF "
          CompileMpi = CompileSmp       
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
          CompileMpi = "mpxlf90_r -qtune=pwr6 -qarch=pwr6 -WF,-DXLF,-DMPI "
          CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
          CompileMpi = "mpxlf90_r -qsmp=omp -qtune=pwr6 -qarch=pwr6 -WF,-DXLF,-DMPI "
          CompileSmp = "xlf90_r -qsmp=omp -qtune=pwr6 -qarch=pwr6 -WF,-DXLF "          
       endif  
    ElseIf( option_sys == "7" ) then
       Write(88,*) "Configuring SELEN for gfortran on CYGWIN"
       CompileSeq = "gfortran -w -DGNU "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
         CompileSmp = "gfortran -w -fopenmp -DGNU "   
	 CompileMpi = CompileSmp 
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
         CompileMpi = "mpif90 -w -DMPI -DGNU "  
	 CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
         CompileMpi = "mpif90 -w -DMPI -DGNU -fopenmp "     
         CompileSmp = "gfortran -w -fopenmp -DGNU "   	 
       endif              
    ElseIf( option_sys == "8" ) then
       Write(88,*) "Configuring SELEN for g95"
       CompileSeq = "g95 -w "
       if( ( option_mpi=='n' ) .and. ( option_omp=='y' ) ) then
 	 Write(* ,*) "Error: OpenMP is not supported by the g95 compiler !"
	 Write(88,*) "Error: OpenMP is not supported by the g95 compiler !"
	 Call stop_config 
	 Stop	       
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='n' ) ) then
         CompileMpi = "mpif90 -w -DMPI "  
	 CompileSmp = CompileSeq
       endif
       if( ( option_mpi=='y' ) .and. ( option_omp=='y' ) ) then
 	 Write(* ,*) "Error: OpenMP is not supported by the g95 compiler !"
	 Write(88,*) "Error: OpenMP is not supported by the g95 compiler !"
	 Call stop_config 
	 Stop	       
       endif       
    Else   
	Write(* ,*) "Error: unknown platform"
	Write(88,*) "Error: unknown platform"
	Call stop_config 
	Stop	       
    EndIf
    if( option_mpi == "y" ) then
       if( option_sys == "6" ) then
!           RunCmd = "poe "
            RunCmd = ""
       elseif( option_sys == "5" ) then
           if( option_omp == "n" ) then
               RunCmd = "mpirun -np "//ntask
           else
               RunCmd = "mpirun -ppn 1 -np "//ntask//" -env OMP_NUM_THREADS "//trim(nthread)//" "
           endif
       elseif( option_sys == "7" ) then
           RunCmd = "mpirun -np "//ntask
       else
           if( option_omp == "n" ) then
               RunCmd = "mpirun -np "//ntask
           else
               RunCmd = "mpirun -bynode -np "//ntask//" -x OMP_NUM_THREADS="//trim(nthread)//" "
	       endif
       endif
    else
       RunCmd=""
    endif
    if( (option_mpi == "n") .and. (option_omp == "n" ) ) then
          CompileMpi = CompileSeq
	  CompileSmp = CompileSeq    
    endif
endif
!
!
!
!
! ###### Node-local storage ######
!
IF(line(1:3)=="980") THEN 
	call scan_string (line, 2, ss, nout)
        option_nls = ss(1) 
        tmp_prefix = ss(2) 
        if((option_mpi=='n') .and. (option_nls=='y')) then
		     Write(88,*) "Node-local storage is available only for MPI jobs"
		     Write(88,*) "Option will be ignored"
		     option_nls='n'
		Endif   
	    if(option_nls=='y') then 
		     Write(88,*) "SELEN will store PX and SH files under: ", trim(adjustl(tmp_prefix))
      	 	 tmp_prefix=trim(adjustl(tmp_prefix))//'/'
        	 len_tmp=len(trim(adjustl(tmp_prefix))) 
      	 	 tmp_prefix="'"//trim(adjustl(tmp_prefix))//"'"
        Endif
ENDIF
!
!
! ###### Solution of the SLE ######  
!
! >>>>>>  Iterations 
!
IF(line(1:3)=="130") THEN 
	call scan_string (line, 2, ss, nout)
	iter   = ss(1)
	mode   = ss(2)	
!
	if    (mode=='1') then 
		Write(88,*) "The SLE will be solved in mode 1 (fully self-gravitating)"
	elseif(mode=='2') then 
		Write(88,*) "The SLE will be solved in mode 2 (elastic approximation)"
	elseif(mode=='3') then 
		Write(88,*) "The SLE will be solved in mode 3 (eustatic approximation)"	
	elseif(mode=='4') then 
		Write(88,*) "The SLE will be solved in mode 4 (Woodward approximation)"		
	elseif(mode=='5') then 
		Write(88,*) "The SLE will be solved in mode 5 (ice load neglected)"		
	endif					
!	
	If(iter/='0') then 
	Write(88,*) "The SLE will be solved with ", iter, " iteration(s)"
	else
	Write(88,*) "The SLE will *not* be solved"
        endif   
ENDIF
!
!
!
! ###### Maximum harmonic degree ###### 
!
IF(line(1:3)=="140") THEN 
	call scan_string (line, 1, ss, nout)
	degree=ss(1) 
 	call CHAR10_2_INT(degree, ndegree)
 	write(88,*) 'Maximum degree is: ', ndegree
 	write(88,*) 'Total number of harmonics:', (ndegree+1)*(ndegree+2)/2
ENDIF 
!
!
! ###### Dealing with harmonic degree "1" ###### 
!
IF(line(1:3)=="145") THEN 
!
	call scan_string (line, 2, ss, nout)
	OPTION_DEG1   = ss(1) 	
	OPTION_RFRAME = ss(2) 
!
	If    (option_deg1=='y') then 
	write(88,*) 'Harmonic degree 1 IS included Yahouh!'	
!
 	If    (option_rframe=='CM') & 
	write(88,*) 'Reference frame of the MASS CENTER of the WHOLE EARTH'	
!
 	If    (option_rframe=='CE') & 
	write(88,*) 'Reference frame of the MASS CENTER of the SOLID EARTH'	
!
 	If    (option_rframe/='CE'.and.option_rframe/='CM') THEN 
		Write(88,*) "Please select a valid option (CE/CM) for the reference frame" 
		Write(*, *) "Please select a valid option (CE/CM) for the reference frame" 
        	CALL STOP_CONFIG 
        Endif	
!	
	elseif(option_deg1=='n') then 
	write(88,*) 'Harmonic degree 1 is NOT included'
!
	elseif(option_deg1/='n'.and.option_deg1/='y') then 
	Write(88,*) "Please select a valid option (y/n) for harmonic degree 1" 
	Write(*, *) "Please select a valid option (y/n) for harmonic degree 1" 
        CALL STOP_CONFIG 
	Endif 
!
ENDIF 

!
!
! ###### Ice model ######
!
! >>>>>> Ice file name 
!
IF(line(1:3)=="170") THEN 
	call scan_string (line, 1, ss, nout)
	ice_file=ss(1)
! 
	Write(88,*) 'The ice file is: ', ice_file	
!
ice_flag=0
!
If    (ice_file=='ice3g.dat'    .or.ice_file=='ice3g_and.dat'.or.& 
       ice_file=='ice3g_ant.dat'.or.ice_file=='ice3g_bal.dat'.or.&
       ice_file=='ice3g_bar.dat'.or.ice_file=='ice3g_bri.dat'.or.& 
       ice_file=='ice3g_gro.dat'.or.ice_file=='ice3g_ice.dat'.or.&  
       ice_file=='ice3g_nam.dat'.or.ice_file=='ice3g_sib.dat')  then 
!
       ice_flag=1 ; ninc='18' ; header_lines=20 ; IMUL=2 ; titlice='ICE3G'
!
ELSEIF(ice_file=='disk.dat'.or.ice_file=='disk_on.dat'.or.&
                               ice_file=='disk_off.dat')        then
!
       ice_flag=1 ; ninc='18' ; header_lines=20 ; IMUL=2 ; titlice='DISK'
!
ELSEIF(ice_file=='icap.dat'.or.ice_file=='icap_on.dat'.or.&
                               ice_file=='icap_off.dat')        then
!
       ice_flag=1 ; ninc='18' ; header_lines=20 ; IMUL=2 ; titlice='ICAP'

!      
ELSEIF(ice_file=='ice1.dat'    .or.ice_file=='ice1_eup.dat'.or. & 
       ice_file=='ice1_gro.dat'.or.ice_file=='ice1_nam.dat')    then 
!
       ice_flag=1 ; ninc='18' ; header_lines=20 ; IMUL=1 ; titlice='ICE1'
!
ELSEIF(ice_file=='ice5g.dat'    .or.ice_file=='ice5g_and.dat'.or.&
       ice_file=='ice5g_ant.dat'.or.ice_file=='ice5g_fen.dat'.or.& 
       ice_file=='ice5g_gre.dat'.or.ice_file=='ice5g_icl.dat'.or.&
       ice_file=='ice5g_lau.dat'.or.ice_file=='ice5g_nwz.dat')  then 
!
       ice_flag=1 ; ninc='21' ; header_lines=24 ; IMUL=1 ; titlice='ICE5G'		
!
ELSEIF(ice_file=='ice5g26.dat')  then 
!
       ice_flag=1 ; ninc='26' ; header_lines=29 ; IMUL=1 ; titlice='ICE5G26'
!
ELSEIF (ice_file=='alpst.dat'   .or.  ice_file=='alpsf.dat'  .or. & 
	ice_file=='alpsh.dat'   .or.  ice_file=='alpsc.dat')    then 	 
! 
       ice_flag=1 ; ninc='21' ; header_lines=4  ; IMUL=1 ; titlice='ALPS'
!
ELSEIF (ice_file(1:6)=='alaska' ) then
!
       ice_flag=1 ; ninc='40' ; header_lines=4  ; IMUL=1 ; titlice='ALASKA'
!       
ELSEIF (ice_file=='ij05mod.dat') then 
!
       ice_flag=1 ; ninc='18' ; header_lines=18 ; IMUL=2 ; titlice='IJ05MOD'	
!
ELSEIF (ice_file(1:5)=='anu05') then 
!
       ice_flag=1 ; ninc='30' ; header_lines=26 ; IMUL=1 ;  titlice='ANU05'
!
ENDIF
!
	  if(ice_flag==0) then
 	  Write(88,*) "The ice file apparently does not exist" 
 	  Write(*, *) "The ice file apparently does not exist" 
          call stop_config 
	  Stop
	  Endif
!
! ---- Counts the ice elements listed in the ice file 
	  xfile='../ICE-MODELS/'//trim(adjustl(ice_file))
          call ice_count(xfile, imul, header_lines, nice)
   	  Write(88,*) 'There are ', nice, 'elements in the ice file ', ice_file 
!
	  len_ice=len(trim(adjustl(ice_file))) 
!
	  ice_filename="'"//trim(adjustl(ice_file))//"'"
!
! ----  Completes the configuration of "alma.inc" 
	If(option_pw=='y') then 
		OPEN(120,file=trim(adjustl(wdir))//"/ALMA/alma.inc",status='old')		
		Call CHAR3_2_INT(NINC, NNINC)
		Write(120,*)"integer, parameter :: p= ",    NNINC, " ! Number of time points minus one" 
		Write(120,*)"real, parameter :: m1=0, m2=", NINC, " ! tmin and tmax"
		Write(120,*)"! "		 
		Write(120,*)"! ********* File closed by config.f90 *********"
		Write(120,*)"! "		 		 
		CLOSE(120)	
	Endif	  
!
ENDIF
!
!
!
IF(line(1:3)=="172") THEN                       ! In progress ... In progress ... In progress ... 
	call scan_string (line, 1, ss, nout)
!
	ice_delta=ss(1)   
!
	Call CHAR10_2_REAL(ice_delta, delta)	
! 
! ....  Any warning here ???????
!
	Write(88,*) 'The ice time step (kyr) is :', delta	
!       Write(* ,*) 'The ice time step (kyr) is :', delta	
!
ENDIF                                           ! In progress ... In progress ... In progress ...
!
! >>>>>>  "Shape factors" file 
!
IF(line(1:3)=="171") THEN 
	call scan_string (line, 2, ss, nout)
	option_sf  = ss(1)
  	shape_file = ss(2) 
!
	if(option_sf=='y') then 
!		Write(88,*) "*New* shape factors will be computed"
!
        Write(88,*) "A *new* ice SH dechomposition will be determined ..."
        Write(88,*) " - this accounts for the ice shape factors AND for the history of ice thickness "
        Write(88,*) " - This computation is time-consuming. Saving the data for future computations  "
        Write(88,*) "   based on the SAME ice and on themaximum harmonic DEGREE is higly rechommended."
!
        else
!		Write(88,*) "Found a pre-built shape factors file ", & 
!		            trim(adjustl(shape_file)) 
!
        Write(88,*) "A pre-built file is possibly available for ice SHs: ", & 
                     trim(adjustl(shape_file)) 
!
	Endif     
ENDIF
!
!
!
!
! ###### Earth model parameters ######  
!
!
!
! >>>>>>  Number of mantle Maxwell VE layers and code 
!         Mantle viscosity profile & lithospheric thickness
!
        IF(line(1:3)=="160") THEN 
!
	option_nm  = 'n'
!
	call scan_string (line, 3, ss, nout)
!       option_nm  = ss(1)
	nv         = ss(1)
	code       = ss(2) 
	visco_file = '../VSC/'//ss(3)	
!
        option_nm='y'
!
        if(option_nm=='y') then     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
        Write(88,*) "Love numbers & spectra will be computed by Normal Modes (TABOO)"		
	Write(88,*) "... number of mantle v/e layers (NV): ", NV 
        Write(88,*) "... Earth model CODE (see TABOO User guide): ", CODE
!
!	 
! ---   Number of lines in the viscosity file 
!
	IF(NV/='0') THEN 
!
        INQUIRE(FILE=visco_file,EXIST=lex)
!
        If(lex) then 
	Write(88,*) "The viscosity file is: ", trim(adjustl(visco_file))
	        Else 
	Write(88,*) "The file ", trim(adjustl(visco_file)), " has NOT been found in ./VSC"
        Write(* ,*) "The file ", trim(adjustl(visco_file)), " has NOT been found in ./VSC"
        call Stop_Config
        	Endif
!
	open(12,file=visco_file,status='unknown') 
!
	LEN_VISCO=len(trim(adjustl(visco_file))) 
	visco_filename="'"//trim(adjustl(visco_file))//"'"
	NLINES=0
		do j=1, large_integer
		read(12,*,end=88099) junk
!		write(*,'(a80)')JUNK
		NLINES=NLINES+1
		enddo
	88099   close(12) 	
	call CHAR3_2_INT(nv, nnv)		 
	if(nnv==NLINES-1) then 
		Write(88,*) "The number of lines in file ", trim(adjustl(visco_file)), & 
			   " is consistent with NV=", NNV
		           else
		Write(88,*) "The number of lines in file ", trim(adjustl(visco_file)), & 
		 	    "is NOT consistent with NV=", NNV
		Write(88,*) "Please check config.dat and ", trim(adjustl(visco_file))
		Write(*,*)  "The number of lines in file ", trim(adjustl(visco_file)), & 
		 	    " is NOT consistent with NV=", NNV
		Write(*,*)  "Please check config.dat and ", trim(adjustl(visco_file))			  
	       call Stop_config 			    
	endif
!	
! ---   Reading litho thickness and viscosity values from "visco_file"
!
	open(12,file=visco_file,status='unknown') 
	read(12,'(a200)') row 
	call scan_string (row, 1, ss, nout)
	litho=ss(1) 
!	
		do j=1, nnv 
			read(12,'(a200)') row
			call scan_string (row, 1, ss, nout)
			vsc(j)=ss(1)			
		enddo
	close(12)
!
! ---   Viscosity string for output on monitor & GMT plots
!
	open(44,file='visco.tmp',status='unknown') 
	if(nnv==1) then 
		   write(44,*) "/",trim(vsc(1)),"/" 
        	   else 
		   Write(44,*) "/",trim(vsc(1))," ",& 
		   		  (trim(vsc(j))," ",j=2,nnv-1),&
		   		   trim(vsc(nnv)),"/"		       	       
	endif
	close(44) 
!	
	open(44,file='visco.tmp',status='unknown') 
	read(44,'(a100)') vstring 
	close(44) 	     	
!
	ENDIF	
!
	ENDIF                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
ENDIF
!
!
!
!
!
!
! ###### Tegmark resolution ######
!
IF(line(1:3)=="150") THEN 
	call scan_string (line, 1, ss, nout)
        resolution=ss(1)
!
 	call CHAR10_2_INT(resolution, nresolution)
!
 	np=2*nresolution*(nresolution-1)*20+12
!
 	Write(88,*) "Tegmark's resolution, r: ", nresolution 
 	Write(88,*) "Number of pixels: ", np
!
!
! ---- A lower bound for Resolution 
 	If(nresolution.lt.res_min) then 
 	Write(*, *) "------ Resolution it TOO LOW -------"
 	Write(*, *) "  Please use a resolution RES>=12   "
 	Write(88,*) "------ Resolution it TOO LOW -------"
 	Write(88,*) "  Please use a resolution RES>=12   "
        Write(88,*) 
        call Stop_Config
	stop 
 	endif 

!
! ---- Testing Tegmark's condition 
 	If(np>=ndegree**2/3.and.np/=0) then 
 	Write(88,*) "------ 'Tegmark condition' is met ------"
 	Write(88,*) "(enough pixels for this harmonic degree)"
				       else 
 	Write(88,*) " ------- Tegmark condition is NOT met ------- "
 	Write(88,*) " (not enough pixels for this harmonic degree) "
 	endif 
ENDIF 
!
!
!
!
! ###### Spherical harmonics (SH) file ######
!
IF(line(1:3)=="180") THEN 
	call scan_string (line, 2, ss, nout)
        option_sh = ss(1) 
        sh_file   = ss(2) 
	if(option_sh=='y') then 
		              Write(88,*) "SELEN will prepare the new SH file: ", sh_file
		           else
			      INQUIRE(FILE=sh_file,EXIST=lex)
			      If(lex) then 
			           Write(88,*) "SELEN will use the pre-built SH file: ", trim(adjustl(sh_file))
			      Else 
			           Write(88,*) "The SH file ", trim(adjustl(sh_file)), " has NOT been found"
			           Write(* ,*) "The SH file ", trim(adjustl(sh_file)), " has NOT been found"
			           call Stop_Config
			      Endif		   
       Endif 
ENDIF
!
!
!
!
!
! ###### OF SH dechomposition ######
!
IF(line(1:3)=="190") THEN 
	call scan_string (line, 2, ss, nout)
        option_oh = ss(1) 
        shof_file = ss(2) 
!
        option_rof ='r'
!
	if(option_oh=='y') then 
		 Write(88,*) "SELEN will prepare the new SH OF file: ", shof_file
		           else
			      INQUIRE(FILE=shof_file,EXIST=lex)			   
			      If(lex) then 
			           Write(88,*) "SELEN will use the pre-built SH OF file: ", trim(adjustl(shof_file))
			      Else 
			           Write(88,*) "The SH OF file ", trim(adjustl(shof_file)), " has NOT been found"
			           Write(* ,*) "The SH OF file ", trim(adjustl(shof_file)), " has NOT been found"
			           call Stop_Config
			      Endif		   
        Endif 
ENDIF
!
!
!
!
! ###### Reading the depot Label ######  		
!
IF(line(1:3)=="195") THEN
	call scan_string (line, 1, ss, nout)
		run=ss(1)
		write(88,*) 'The depot label is: ', run	
		depot="./depot-"//trim(adjustl(run))
		runp="'"//trim(adjustl(run))//"'"
ENDIF  
!
!
100 CONTINUE
!
444 close(1)
!
!
!
! # Part 2/1 Reading the output settings from "config.dat"   
!
!
  open(1,file='config.dat',status='unknown')
!
!
  Write(88,*) ""
  Write(88,*) "+----------------------------------------------+"
  Write(88,*) "|               Outputs of SELEN               |"
  Write(88,*) "+----------------------------------------------+"
!
!
  DO 200 I=1,LARGE_INTEGER
!
!
  read(1,'(a200)',end=555) line 
!
!
! ###### How to deal with GMT scripts ###### 
!
IF(line(1:3)=="200") THEN 
	call scan_string (line, 1, ss, nout)
	option_GMT = ss(1) 
	If(option_gmt=='n') & 
	Write(88,*) "The execution of GMT scripts is switched off"
	If(option_gmt/='y'.and.option_gmt/='n') then 
		Write(* ,*) "For the GMT switch, only y/n are valid options"
		Write(88,*) "For the GMT switch, only y/n are valid options"
		Call stop_config 
		Stop	
	Endif
	
	
	
ENDIF
!
!
! ###### Pixelization ###### 
!
! >>>>>> Producing pixelization maps
!
IF(line(1:3)=="205") THEN 
	call scan_string (line, 1, ss, nout)
	option_px = ss(1) 
	If(option_px=='y') & 
	Write(88,*) "Maps of the pixelizations will be drawn"
ENDIF
!
! >>>>>> Computing and drawing the function 
!
IF(line(1:3)=="206") THEN 
	call scan_string (line, 1, ss, nout)
	option_wi = ss(1) 
	If(option_wi=='y') &
	Write(88,*) "The window function will be computed & plotted"
ENDIF
!
!
!
!
! ###### Ocean function rechonstruction and mapping ###### 
!
IF(line(1:3)=="210") THEN 
	call scan_string (line, 1, ss, nout)
	option_of = ss(1) 
	If(option_of=='y') & 
	Write(88,*) "The OF will be rechonstructed & mapped"
ENDIF
!
! ###### Plot of Ocean Function "degree variance" ###### 
!
IF(line(1:3)=="215") THEN 
	call scan_string (line, 1, ss, nout)
	option_ofdv = ss(1) 
	If(option_ofdv=='y') & 
	Write(88,*) "The OF DV will be computed & plotted"
ENDIF
!
!
!
!
!
! ###### Ice sheets graphical outputs ###### 
!
! >>>>>> Maps of original ice sheets  
!
IF(line(1:3)=="220") THEN 
	call scan_string (line, 1, ss, nout)
	option_or = ss(1)
	If(option_or=='y') & 
	Write(88,*) "Maps of original ice sheets will be drawn"
ENDIF
!
! >>>>>> Plot of Equivalent Sea Level (ESL) vs time   
!
IF(line(1:3)=="221") THEN 
	call scan_string (line, 1, ss, nout)
	option_esl = ss(1)
	If(option_esl=='y') & 
	Write(88,*) "A plot of Equivalent Sea level will be drawn"
ENDIF
!
! >>>>>> Ice sheets rechonstruction from SH coefficients & mapping
!
IF(line(1:3)=="222") THEN 
	call scan_string (line, 1, ss, nout)
	option_ri = ss(1)
	If(option_ri=='y') & 
	Write(88,*) "The ice sheets will be rechonstructed from SH coefficients and mapped"
ENDIF
!
!
!
!
! ###### Earth model Spectral properties ###### 
!
! >>>>>> Diagrams of LDCs, relaxation spectrum & residues
!	
IF(line(1:3)=="230") THEN 
	call scan_string (line, 1, ss, nout)
	option_ln = ss(1) 
	If(option_ln=='y') & 
	Write(88,*) "The Love numbers will be plotted"
ENDIF
!
!
!
!
! ###### Relative Sea level (RSL) analysis ###### 
!
! >>>>>>  RSL database
!
IF(line(1:3)=="240") THEN 
!
	call scan_string (line, 3, ss, nout)
!
	OPTION_RSLA = SS(1) 
!	
        If(option_rsla=='y') then 
!
		RSL_FILE = '../DATA/'//SS(2)
!		
		RSL_DATABASE_FORMAT = SS(3)		
!
        	INQUIRE(FILE=RSL_FILE,EXIST=lex)
!
        	if(lex)then 
!
		len_rsl=len(trim(adjustl(rsl_file))) 
!
		rsl_database="'"//trim(adjustl(rsl_file))//"'"
!
! ---   Counts the number of RSL sites in database  
!
	IF(RSL_DATABASE_FORMAT=='0')    then 
!
		open(10,file=rsl_file,status='old') 
!
		nrsl=0
		do 33433 j=1, very_large_integer
		read(10,'(a200)',end=414) row 
			do k=1, 200
			do n=1, NALFA
				if(row(k:k)==alfa(n)) then 
				nrsl=nrsl+1
				goto 33433
				endif 
			enddo
			enddo
		33433 continue
414	close(10) 
!
	Elseif(RSL_DATABASE_FORMAT=='1'.or.RSL_DATABASE_FORMAT=='2') then 
!	
		open(10,file=rsl_file,status='old') 
!
		nrsl=0
		do 33533 j=1, very_large_integer
		read(10,'(a200)',end=484) row 
!		
!               write(*,'(a100)') row 
!		
			do k=1, 200
			do n=1, NALFA
				if(row(k:k)==alfa(n)) then 
				nrsl=nrsl+1
				goto 33533
				endif 
			enddo
			enddo
			
		33533 continue
484	close(10)		
	If(RSL_DATABASE_FORMAT=='1')nrsl=nrsl/2 
	If(RSL_DATABASE_FORMAT=='2')nrsl=nrsl/3 
!			
	Endif
!
        Write(88,*) 'The RSL database is: ', trim(rsl_file), ' ---', nrsl, ' RSL sites'
!
	If(nrsl==0) then 
		Write(88,*) "The RSL database has no elements"
		Write(*, *) "The RSL database has no elements"
		Call stop_config 
		Stop	
	Endif
!
!
! --- File for scattered RSl data 
 		open(82,file='scatter-data.dat',status='unknown') 
!
! --- Reading the coordinates of the RSL sites
!
	        If    (RSL_DATABASE_FORMAT=='0') then 	
!
  		OPEN(44,FILE=RSL_FILE,STATUS='old')
!
     		ii=1 
299  		READ (44, 699, END=199) j, lats, lons, nrsl_data     
   		IF(lons<=0.) lons=360.+lons		    
   		do k=1, nrsl_data 
           		READ (44 ,*) time_rsl, rjunk, rsl_datum, d_rsldatum
	   	Write(82,*) time_rsl/1000., rsl_datum, d_rsldatum 
   		enddo
   		ii=ii+1
   		IF(ii <= nrsl) GOTO 299
199  CLOSE(44)
699  FORMAT(1X,I3,1X,F5.1,1X,F6.1,1X,I2,1X,A22)
!
!
		Elseif(RSL_DATABASE_FORMAT=='1') then
!
  		OPEN(44,FILE=RSL_FILE,STATUS='old')
!
		Do ii=1, nrsl 
			READ(44,'(a10)')CJUNK 
			READ(44,'(a10)')AJUNK 
			READ(44,'(a10)')CJUNK 
			READ(44,*) LATS, LONS 
!									
			IF     (AJUNK(1:1)=="-") THEN 			
				READ(44,*) TIME_RSL, RJUNK, RSL_DATUM, D_RSLDATUM 
							
			ELSEIF (AJUNK(1:1)=="*") THEN 	
					
                                READ(44,*) TIME_RSL1, TIME_RSL2, RSL_DATUM, D_RSLDATUM 
				
			        TIME_RSL =  MIN(TIME_RSL1,TIME_RSL2)+ &
				           (MAX(TIME_RSL1,TIME_RSL2)-MIN(TIME_RSL1,TIME_RSL2))/2.
					   			
			ELSEIF (AJUNK(1:1)/="*".AND.AJUNK(1:1)/="-") THEN 
!						
				Write(88,*) "The RSL database is badly configured"
				Write(*, *) "The RSL database is badly configured"
				Call stop_config 
				Stop					
			ENDIF			 		
			READ(44,'(a10)')CJUNK 	
 			WRITE(82,*) TIME_RSL/1000., RSL_DATUM, D_RSLDATUM		
!			
		Enddo
		Close(44) 
!
		Elseif(RSL_DATABASE_FORMAT=='2') then	
! 
  		OPEN(44,FILE=RSL_FILE,STATUS='old')
!
		Do ii=1, nrsl 
			READ(44,'(a10)')CJUNK 
			READ(44,'(a10)')AJUNK 
			READ(44,'(a10)')CJUNK 
!
			READ(44,'(a200)') linep 
!			
			call scan_string (linep, 2, ss, nout)	
!			
			LATSC10=trim(adjustl(ss(1))) 
			LONSC10=trim(adjustl(ss(2)))
!							
			call FROMDDMMSS_2_DEGREES (LATSC10, LATS) 
			call FROMDDMMSS_2_DEGREES (LONSC10, LONS) 
!									
			IF     (AJUNK(1:1)=="-") THEN 			
				READ(44,*) TIME_RSL, RJUNK, RSL_DATUM, D_RSLDATUM 
							
			ELSEIF (AJUNK(1:1)=="*") THEN 	
					
                                READ(44,*) TIME_RSL1, TIME_RSL2, RSL_DATUM, D_RSLDATUM 
				
			        TIME_RSL =  MIN(TIME_RSL1,TIME_RSL2)+ &
				           (MAX(TIME_RSL1,TIME_RSL2)-MIN(TIME_RSL1,TIME_RSL2))/2.
					   			
			ELSEIF (AJUNK(1:1)/="*".AND.AJUNK(1:1)/="-") THEN 
			 
			
				Write(88,*) "The RSL database is badly configured"
				Write(*, *) "The RSL database is badly configured"
				Call stop_config 
				Stop					
			ENDIF			 		
			READ(44,'(a10)')CJUNK 	
 			WRITE(82,*) TIME_RSL/1000., RSL_DATUM, D_RSLDATUM		
!			
		Enddo
		Close(44) 
!
	ENDIF 
! 
     CLOSE(82) 	     
!
    	 Else
		Write(88,*) "The RSL database does not exists"
		Write(*, *) "The RSL database does not exists"
		Call stop_config 
		Stop	
!     
     	 Endif 	
!
     Else 
		Write(88,*) "NO RSL analysis will be performed "
     Endif     
!
ENDIF
!
!
!
!
IF(option_rsla=='y') THEN 
!
!
!
! >>>>>>  Plotting the distribution of RSL sites in database 
!
IF(line(1:3)=="241") THEN 
	call scan_string (line, 1, ss, nout)
	option_rsldb   = ss(1)
	If(option_rsldb=='y') Write(88,*) "The distribution of RSL stations will be plotted"
ENDIF
!
!
! >>>>>> Computing RSL predictions and individual plots of predictions vs. observations 
!
IF(line(1:3)=="242") THEN 
	call scan_string (line, 2, ss, nout)
	option_rsl   = ss(1) 
	option_rslp  = ss(2) 
	If(option_rsl =='y') Write(88,*) "RSL predictions will be obtained"
	If(option_rslp=='y') Write(88,*) "RSL predictions will be plotted against data"
	If(option_rslp=='y'.and.option_rsl=='n') then 
			Write(88,*) "Plotting RSl predictions vs. data requires computing ..." 
			Write(88,*) "... RSL predictions -- please modify file config.dat' --"
                        Write(* ,*) "Plotting RSl predictions vs. data requires computing ..." 
			Write(* ,*) "... RSL predictions -- please modify file config.dat' --"
			call stop_config 			    
			stop 
			endif
ENDIF
!
!
! >>>>>> Scatterplot of RSL observations vs. predictions 
!
IF(line(1:3)=="243") THEN 
	call scan_string (line, 1, ss, nout)
	option_rslsca   = ss(1) 
	If(option_rslsca=='y') Write(88,*) "A RSL scatterplot (global data vs predictions) will be drawn"
	If(option_rslsca=='y'.and.option_rsl=='n') then 
			Write(88,*) "Plotting a RSL 'scatterplot' showing data & predictions requires ..."
			Write(88,*) "... computing RSL predictions first. Please modify file 'config.dat'"
			Write(*,*)  "Plotting a RSL 'scatterplot' showing data & predictions requires ..."
			Write(*,*)  "... computing RSL predictions first. Please modify file 'config.dat'"
			call stop_config 			    			
			stop
			endif				
ENDIF
!
!
! >>>>>> Misfit analysis for RSL   
!
IF(line(1:3)=="244") THEN 
	call scan_string (line, 1, ss, nout)
	option_rslmf   = ss(1) 
	If(option_rslmf=='y') Write(88,*) "A misfit analysis for RSL will be performed"
	If(option_rslmf=='y'.and.option_rsl=='n') then 
			Write(88,*) "Plotting the RSL <<misfit histogram>> requires computing "
			Write(88,*) "... RSL predictions first. Please modify file config.dat'"
                        Write(*, *) "Plotting the RSL <<misfit histogram>> requires computing "
			Write(*, *) "... RSL predictions first. Please modify file config.dat'"
			Call Stop_Config 
			stop
			endif		
ENDIF
!
!
!
! >>>>>> RSL table with all observations & predictions   
!
IF(line(1:3)=="245") THEN 
	call scan_string (line, 1, ss, nout)
	option_rsltab   = ss(1) 
	If(option_rsltab=='y') Write(88,*) "A table with all data and predictions will be printed" 
	If(option_rsltab=='y'.and.option_rsl=='n') then 
			Write(88,*) "Printing a RSL table with all predictions requires computing ..."
			Write(88,*) "... RSL predictions first. ** Please modify file config.dat ** '"
			Write(* ,*) "Printing a RSL table with all predictions requires computing ..."
			Write(* ,*) "... RSL predictions first. ** Please modify file config.dat ** '"
			Call Stop_Config 
			stop
			endif	
ENDIF
!
ENDIF    ! on "option_rsla"
!
!
!
! >>>>>> RSL "regions"    
!
! >>> Global RSL zones
!
IF(line(1:3)=="250") THEN 
	call scan_string (line, 1, ss, nout)
	option_rslz   = ss(1) 
	If(option_rslz=='y') Write(88,*) "<<RSL zones>> will be determined"
ENDIF
!
! >>> Regional RSL contour lines 
!
IF(line(1:3)=="251") THEN 
	call scan_string (line, 2, ss, nout)
	option_rslc     = ss(1) 
!
	file_region     = '../DATA/'//ss(2) 
!
	if(option_rslc=='y') then 	
!
        INQUIRE(FILE=file_region,EXIST=lex)
!
	If(lex) then
!
                Write(88,*) "<<RSL contour lines>> will be drawn" 
		Write(88,*) "The RSL region is described in file ", file_region
!
! --- Counting the number of points ("virtual RSL sites") in rectangle
!
	        file_region_lonlat='lonlat_rslc.dat'
!
! --- Sbr. Scan_region scans the user-supplied input file "file_region" in order to:
!		- Check lon-lat bounds 
!		- Check the contour interval 
!		- Count the number of "Virtual RSL" sites
! 	        - Report lon-lat coordinates on file "file_region_lonlat"
!		---> NOTICE that the rectangle is enlarged by 2 
!		     degrees NSEW to avoid border effects with GMT 
!
  	        Call scan_region (file_region, file_region_lonlat, nrslc, time_bpc, time_bpcc, & 
				  LONMINC, LONMAXC, LATMINC, LATMAXC, & 
				  min_rslc, max_rslc, RSL_INT, & 
				  name_of_region)				 	
!
	        LEN_RSLC = len(trim(adjustl(file_region_lonlat)))
	        RSLC_LONLAT_FILE     = "'"//trim(adjustl(file_region_lonlat))//"'"	
!
                Write(88,*) "There are ", nrslc, " points in file ", trim(adjustl(file_region)) 
                Write(88,*) "Time of analysis is ", time_bpc, " ka" 
		Write(88,*) "Coordinates of RSL <<sites>> are reported on file ", & 
		             trim(adjustl(file_region_lonlat))
!		
   	        else 
                	Write(*, *) "File ", trim(adjustl(file_region)), " apparently"
			Write(*, *) "does not exist -  Please check file <config.dat>"
                	Write(88,*) "File ", trim(adjustl(file_region)), " apparently"
			Write(88,*) "does not exist -  Please check file <config.dat>"
			Call Stop_Config 
			stop		 
		endif 
	endif 
ENDIF
!
!
!
!
! ###### Sea level change at tide-gauge stations ###### 
!
!
! >>>>>> Tide-gauges database 
!
IF(line(1:3)=="260") THEN 
!
	call scan_string (line, 2, ss, nout)
!	
 	OPTION_TGA    = SS(1)
!
	TGAUGES_FILE  = '../DATA/'//SS(2) 
!
        if(option_tga=='y') then 
!
! ---   Counts the number of tide-gauge sites in database  
	
        INQUIRE(FILE=TGAUGES_FILE,EXIST=lex)
!
	If(lex)then 
!
		Write(88,*) 'The tide gauges database is: ', TGAUGES_FILE
!
		LEN_TGAUGES=len(trim(adjustl(TGAUGES_FILE))) 
		TGAUGES_DATABASE="'"//trim(adjustl(TGAUGES_FILE))//"'"
!
		open(10,file=TGAUGES_FILE,status='old') 
		NTIDEGAUGES=0
		do 33833 j=1, very_large_integer
			read(10,'(a200)',end=464) row 
			if(row(1:1)/='#') NTIDEGAUGES=NTIDEGAUGES+1
		33833 continue
		464 continue 
        	Write(88,*) 'There are ', NTIDEGAUGES, 'tide-gauges rechords in file ', TGAUGES_FILE
 		close(10)
!
		else
!	
		Write(88,*) 'No tide-gauge database has been found...'	
		Write(*, *) 'No tide-gauge database has been found...'	
                Call Stop_Config 
                stop		 		
!			
		endif	
!
	endif
!
ENDIF
!
!
IF(option_tga=='y') THEN 
!
! >>>>>> Plotting tide-gauge stations distribution 
!
IF(line(1:3)=="261") THEN 
	call scan_string (line, 1, ss, nout)
	option_tgplot = ss(1) 
	If(option_tgplot=='y')&
	Write(88,*) 'The tide gauge sites distribution will be mapped'
ENDIF
!
! >>>>>> Tide gauge data scatterplot & averages
!
IF(line(1:3)=="262") THEN 
	call scan_string (line, 1, ss, nout)
	option_tgsca = ss(1) 
	If(option_tgsca=='y')&
	Write(88,*) 'Tide gauge data scatterplots will be drawn'
ENDIF
!
!
! >>>>>> Predictions at tide-gauges 
!
IF(line(1:3)=="262") THEN 
	call scan_string (line, 1, ss, nout)
	option_tg = ss(1) 
	If(option_tg=='y')& 
	Write(88,*) 'Predictions of SLC at tide-gauges will be obtained'
ENDIF
!
ENDIF
!
!
!
!
! ###### Present-day rates ###### 
!
! >>>>>> *Global* maps of dot S, U, and N  
!
IF(line(1:3)=="270") THEN 
	call scan_string (line, 1, ss, nout)
	option_gm = ss(1) 
	If(option_gm=='y') then  
	Write(88,*) "\dot S, U, and N at present-time will be computed..."
	Write(88,*) "... and GLOBAL maps will be drawn."
	Endif	
ENDIF
!
!
! >>>>>> *Regional* maps of dot S, U, and N  
!
IF(line(1:3)=="280") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(0) = ss(1) 
!
	If(option_rm(0)=='y'.and.option_rof=='r') then 
	Write(88,*) "\dot S, U, and N at present-time will be computed..."
	Write(88,*) "... and some REGIONAL maps will be drawn."
	Endif
!
	If(option_rm(0)=='y'.and.option_rof=='z') then 
	Write(88,*) "NO REGIONAL maps for a zonal OF"
	option_rm(0)='n'	
	Endif
ENDIF
!
IF(line(1:3)=="281") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(1) = ss(1) 
	If(option_rm(1)=='y'.and.option_rm(0)=='y') & 
	Write(88,*) "Regional Map of \dot S, U, and N for ITALY"
ENDIF
IF(line(1:3)=="282") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(2) = ss(1) 
	If(option_rm(2)=='y'.and.option_rm(0)=='y') & 
        Write(88,*) "Regional Map of \dot S, U, and N for the MEDITERRANEAN"
ENDIF
IF(line(1:3)=="283") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(3) = ss(1) 
	If(option_rm(3)=='y'.and.option_rm(0)=='y') & 
	Write(88,*) "Regional Map of \dot S, U, and N for EUROPE"
ENDIF
IF(line(1:3)=="284") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(4) = ss(1) 
	If(option_rm(4)=='y'.and.option_rm(0)=='y') & 	
	Write(88,*) "Regional Map of \dot S, U, and N for FENNOSCANDIA"
ENDIF
IF(line(1:3)=="285") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(5) = ss(1) 
        If(option_rm(5)=='y'.and.option_rm(0)=='y') & 
	Write(88,*) "Regional Map of \dot S, U, and N for GREENLAND"
ENDIF
IF(line(1:3)=="286") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(6) = ss(1) 
	If(option_rm(6)=='y'.and.option_rm(0)=='y') & 
        Write(88,*) "Regional Map of \dot S, U, and N for NORTH AMERICA"
ENDIF
IF(line(1:3)=="287") THEN 
	call scan_string (line, 1, ss, nout)
	option_rm(7) = ss(1) 
	If(option_rm(7)=='y'.and.option_rm(0)=='y') & 
	Write(88,*) "Regional Map of \dot S, U, and N for ANTARCTICA"
ENDIF
!
!
! ###### 3D VELOCITY ###### 
!
! >>>>>> 3D velocity and dot S & N at user-supplied geodetic sites...  
!
IF(line(1:3)=="275") THEN 
!
	call scan_string (line, 2, ss, nout)
!
	option_3d = 'n'
!
	option_3d = ss(1) 
!	
        If(option_3d=='y') THEN 
!
!	If(ice_type/='h') then   
!		Write(88,*) "3D velocity at sites is only av available for Holocene ice"	
!		Write(*, *) "3D velocity at sites is only av available for Holocene ice"	
!		Call Stop_Config 
!	Endif 
!	
	file_3d = '../DATA/'//ss(2) 
!
        Write(88,*) "Point estimates for Glacial Isostatic Adjustment (GIA)"	   
	Write(88,*) "3D velocities and \dot S and N at sites on file ", & 
		     trim(adjustl(file_3d))
!
! ---   Counts the number of sites in database  
!
        INQUIRE(FILE=FILE_3D,EXIST=lex)
!
        If(lex)then 
!
		LEN_GEOD=len(trim(adjustl(FILE_3D))) 
		GEODETIC_DATABASE="'"//trim(adjustl(FILE_3D))//"'"
!
        	open(10,file=file_3d,status='old') 
!
        	ngeod=0
		do 33733 j=1, very_large_integer
		read(10,'(a200)',end=474) row 
			do k=1, 200
			do n=1, NALFA
				if(row(k:k)==alfa(n)) then 
				ngeod=ngeod+1
				goto 33733
				endif 
			enddo
			enddo
		33733 continue
474	close(10) 
!
        Write(88,*) 'There are ', NGEOD, 'geodetic points in file ', FILE_3D
!
        else
!	
        Write(88,*) 'No geodetic database has been found...'	
        Write(*, *) 'No geodetic database has been found...'	
        Call Stop_Config ; stop	
!
	ENDIF   ! On the existence of the database ... 
!
	ENDIF   ! On the option ... 
!
 ENDIF 
!
!
!
! ###### Stokes coefficients (SC) ###### 
!
! >>>>>> Rate of change of SC, and range of degrees  
!
IF(line(1:3)=="290") THEN 
	call scan_string (line, 3, ss, nout)
	option_st     = ss(1)
	degree_st_min = ss(2)
	degree_st_max = ss(3)
	If(option_st=='y') Write(88,*) 'Range of degrees for the Stokes coefficients: ', & 
	    	  	   trim(degree_st_min), '-', trim(degree_st_max)
ENDIF
!
!
!
200 CONTINUE
!
555 close(1)
!
!
!
!
!
! ====================================
!
! 	Part #2 : building "selen.sh"
!
! 	1/2: Compilation commands 
! 	2/2:  Execution commands
!
! ====================================
!
!
!
  Write(88,*) ""
  Write(88,*) "Compilation and execution commands are written on selen.sh"
!
  open(2,file='selen.sh',status='unknown',recl=256)
!
! 
! >>>>>> A time stamp on "selen.sh"
!
  call DATE_AND_TIME (date,timc)      
  Write(2,*) " "
  Write(2,*) '# This is file "selen.sh", created by "config.f90" on ', & 
             date(1:4), '.', date(5:6), '.', date(7:8), ' ', & 
	     timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  Write(2,*) "  "
!
!
! >>>>>>  # Part 1/2: compilation 
!
Write(2,*) "echo ''"
Write(2,*) "echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'"
Write(2,*) "echo '     SELEN, a Sea levEL EquatioN solver, Version 2.9     '"
Write(2,*) "echo '                  gfortran VERSION                       '" 
Write(2,*) "echo '      Web page: http://fcolleoni.free.fr/SELEN.html      '"	
Write(2,*) "echo '   http://www.fis.uniurb.it/spada/SELEN_minipage.html    '"   
Write(2,*) "echo '   Send comments, requests of help and suggestions to:   '" 
Write(2,*) "echo '                <giorgio.spada@gmail.com>                '"
Write(2,*) "echo '                            -                            '"
Write(2,*) "echo '                    Copyright(C) 2008                    '"    
Write(2,*) "echo '     Giorgio Spada, Florence Colleoni & Paolo Stocchi    '"
Write(2,*) "echo '                          * * *                          '"
Write(2,*) "echo '     This programs comes with  ABSOLUTELY NO WARRANTY    '"
Write(2,*) "echo ' This is free software and you are welcome to distribute '"
Write(2,*) "echo '              it under certain conditions.               '"
Write(2,*) "echo '    For details, visit  <http://www.gnu.org/licenses/>   '"
Write(2,*) "echo '                  or edit file COPYING                   '"
Write(2,*) "echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'"
Write(2,*) "echo ''"
!
Write(2,*) " "
Write(2,*) " "
Write(2,*) "echo '---------------------------------------'"
Write(2,*) "echo ' >>> 0. Compiling the SELEN programs...'"
Write(2,*) "echo '---------------------------------------'"
!
!
!
!
!
! --- Compilation: shtools module
! Write(2,*) "echo '--- shtools.f90'"
! Write(2,*) trim(CompileSeq), " -O3 -c shtools.f90" 
!
!
! --- Compilation: harmonics module
! Write(2,*) "echo '--- harmonics.f90'"
! Write(2,*) trim(CompileSeq), " -O3 -c harmonics.f90" 
!
!
! --- Compilation: SLE
! Write(2,*) "echo '--- sle.f90'"
! Write(2,*) trim(CompileMpi), " sle.f90 harmonics.o -o sle.exe -O3" 
!

!
! --- Compilation: Ice sheets ----
!
! --- Computation of ice Shape factors and SH dechomposition 
If(option_sf=='y') then    
!
! Write(2,*) "echo '--- shape_factors.f90'" 
! Write(2,*) trim(CompileSmp), " shape_factors.f90 harmonics.o -o shapefactors.exe -O3" 
!
! Write(2,*) "echo '--- shice.f90'"
! Write(2,*) trim(CompileSmp), " shice.f90 harmonics.o -o shice.exe"
!
Endif
!
! --- Equivalent sea level 
If     (option_esl=='y') then 
! 	Write(2,*) "echo '--- esl.f90'"
! 	Write(2,*) trim(CompileSeq), " esl.f90 harmonics.o -o esl.exe"
endif 
!
! --- Ice sheets contours 
If(option_or=='y') then 
! Write(2,*) "echo '--- ms.f90'"
! Write(2,*) trim(CompileSeq), " ms.f90 harmonics.o -o ms.exe"
endif
!
!
! --- Compilation of Love numbers tools 
If    (option_nm=='y') then 
!
! --- TABOO
! Write(2,*) "echo '--- tb.f90'"
! Write(2,*) trim(CompileSeq), " tb.F90 harmonics.o -o tb.exe"  
!
elseif(option_pw=='y') then 
 continue 
!
elseif(option_pwa=='y') then 
 continue
!
endif 
!
! --- Compilation: Pixelization (i) 
! Write(2,*) "echo '--- px.f90'"      
! Write(2,*) trim(CompileSeq), " px.f90 harmonics.o -o px.exe" 
!
! --- Compilation: Pixelization (ii) 
! Write(2,*) "echo '--- px_rec.f90'"      
! Write(2,*) trim(CompileSeq), " px_rec.f90 harmonics.o  -o pxrec.exe" 
!
! --- Compilation: Pixelization partitioning
if (option_mpi=='y') then
!   Write(2,*) "echo '--- px_part.f90'"      
!   Write(2,*) trim(CompileSeq),  " px_part.f90 harmonics.o -o pxpart.exe" 
endif
!
! --- Compilation: Copy to local storage
if (option_nls=='y') then
!   Write(2,*) "echo '--- px_copy.f90'"      
!   Write(2,*) trim(CompileMpi),  " px_copy.f90 harmonics.o -o pxcopy.exe" 
endif
!
! --- Compilation: Parallel wet/dry pixel separation
if (option_mpi=='y') then
!   Write(2,*) "echo '--- px_select.f90'"      
!   Write(2,*) trim(CompileMpi),  " px_select.f90 harmonics.o -o pxselect.exe" 
endif
!
! --- Compilation: Spherical harmonics 
If(option_sh=='y') then 
! Write(2,*) "echo '--- sh.f90'"      
! Write(2,*) trim(CompileMpi), " sh.f90 harmonics.o -o sh.exe" 
Endif
!
! --- Compilation: Window function
If(option_wi=='y') then 
! Write(2,*) "echo '--- wnw.f90'"      
! Write(2,*) trim(CompileSeq), " wnw.f90 harmonics.o -O3 -o wnw.exe" 
Endif
!     
! --- Compilation: Ocean function (OF) harmonics 
If(option_oh=='y') then 
! Write(2,*) "echo '--- sh_of.f90'"   
! Write(2,*) trim(CompileMpi), " sh_of.f90 harmonics.o -o shof.exe"   
Endif 
!
! --- Compilation: OF DV computation  
If(option_ofdv=='y') then 
! Write(2,*) "echo '--- of_dv.f90'"  
! Write(2,*) trim(CompileSeq), " of_dv.f90 harmonics.o -o ofdv.exe"
Endif
!
! --- Compilation: OF Rechonstruction 
If(option_of=='y') then 
! Write(2,*) "echo '--- rec_of.f90'"  
! Write(2,*) trim(CompileMpi), " rec_of.f90 harmonics.o -o recof.exe"
Endif
!
! --- Compilation: Ice sheets Rechonstruction 
If(option_ri=='y') then 
! Write(2,*) "echo '--- rec_ice.f90'" 
! Write(2,*) trim(CompileMpi), " rec_ice.f90 harmonics.o -o recice.exe"
Endif 
!
! --- Compilation: Global maps 
If(option_gm=='y') then 
! Write(2,*) "echo '--- gmaps.f90'"
! Write(2,*) trim(CompileMpi), " gmaps.f90 harmonics.o -o gmaps.exe"
Endif
!
! --- Compilation: Regional maps 
If(option_rm(0)=='y') then 
! Write(2,*) "echo '--- rmaps.f90'"
! Write(2,*) trim(CompileMpi), " rmaps.f90 harmonics.o -o rmaps.exe"
Endif
!
! --- Compilation: 3D veolcity and \dot S, N, U at sites
If(option_3d=='y') then 
! Write(2,*) "echo '--- geo.f90'"
! Write(2,*) trim(CompileSmp), " geo.f90 harmonics.o -o geo.exe"
Endif
!
! --- Compilation: 3D veolcity and \dot U at pixels on maps
If(option_3d_regions=='y') then 
! Write(2,*) "echo '--- geo_maps.f90'"
! Write(2,*) trim(CompileSeq), " geo_maps.f90 harmonics.o -o geomaps.exe"
Endif
!
! --- Compilation: SH at Relative Sea Level (RSL) sites 
If(option_rsl=='y') then
! Write(2,*) "echo '--- sh_rsl.f90'"
! Write(2,*) trim(CompileSeq), " sh_rsl.f90 harmonics.o -o shrsl.exe" 
!
! --- Compilation: RSL at RSL sites 
! Write(2,*) "echo '--- rsl.f90'" 
! Write(2,*) trim(CompileSeq), " rsl.f90 harmonics.o -o rsl.exe"
Endif
!
! --- Compilation: RSL Zones 
If(option_rslz=='y') then
! Write(2,*) "echo '--- rsl_zones.f90'"
! Write(2,*) trim(CompileSeq), " rsl_zones.f90 harmonics.o -o rslzones.exe" 
Endif
!
! --- Compilation: SH at virtual sites for RSL contours
If(option_rslc=='y') then
! Write(2,*) "echo '--- sh_rslc.f90'"
! Write(2,*) trim(CompileSeq), " sh_rslc.f90 harmonics.o -o shrslc.exe" 
!
! --- Compilation: RSL at virtual sites for RSL contours 
! Write(2,*) "echo '--- rslc.f90'" 
! Write(2,*) trim(CompileSeq), " rslc.f90 harmonics.o -o rslc.exe"
Endif
!
! --- Compilation: SH at tide gauges
If(option_tg=='y') then
! Write(2,*) "echo '--- sh_tgauges.f90'" 
! Write(2,*) trim(CompileSeq), " sh_tgauges.f90 harmonics.o -o shtgauges.exe" 
!
! --- Compilation: SL change at tide gauges
! Write(2,*) "echo '--- tgauges.f90'" 
! Write(2,*) trim(CompileSeq), " tgauges.f90 harmonics.o -o tgauges.exe"
endif
!
! --- Stokes coefficients
If(option_st=='y') then
! Write(2,*) "echo '--- stokes.f90'"
! Write(2,*) trim(CompileSeq), " stokes.f90 harmonics.o -o stokes.exe" 
Endif
!	      
! --- End of compilation
!
!
!
Write(2,*) " "
Write(2,*) " "
Write(2,*) "echo ''"
Write(2,*) "echo '----------------------------'"
Write(2,*) "echo ' >>> 1. Executing SELEN  ...'"
Write(2,*) "echo '----------------------------'"
!
!
!
! ===================================================
! --- Creating folders into the ./depot directory ---
! ===================================================
!
! ================================
! EXE 00 --- Working directory --- 
! ================================
Write(2,*) " "
Write(2,*) "#echo --------------------------------------------------------"
Write(2,*) "echo                                                          "
Write(2,*) " echo '---> Working directory: '", trim(adjustl(wdir)) 
Write(2,*) "#echo --------------------------------------------------------"


Write(2,*) " "
Write(2,*) "#echo ------------------------------------------------------"
Write(2,*) "echo"
Write(2,*) " echo '---> Output data will be stored into directory '", trim(adjustl(depot))
Write(2,*) "#echo ------------------------------------------------------"
Write(2,*) "if [ ! -d ", depot, " ]"
Write(2,*) "then"
Write(2,*) "mkdir ",depot//"/" 
Write(2,*) "mkdir ",depot//"/bin/" 
Write(2,*) "mkdir ",depot//"/px/" 
Write(2,*) "mkdir ",depot//"/wnw/"
Write(2,*) "mkdir ",depot//"/log/"	
Write(2,*) "mkdir ",depot//"/of/"
Write(2,*) "mkdir ",depot//"/of/degree_variance/" 
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))//"/original/"
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))//"/esl/"
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))//"/rechonstructed/"	
Write(2,*) "mkdir ",depot//"/"//trim(adjustl(titlice))//"/sh/"
Write(2,*) "mkdir ",depot//"/TABOO/"
Write(2,*) "mkdir ",depot//"/Love-Numbers-by-TABOO/"
Write(2,*) "mkdir ",depot//"/rsl/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-sites/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-curves/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-curves/ps/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-curves/pdf/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-scplot/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-misfit/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-table/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-zones/"
Write(2,*) "mkdir ",depot//"/rsl/rsl-contours/"
Write(2,*) "mkdir ",depot//"/tgauges/"
Write(2,*) "mkdir ",depot//"/tgauges/tgauges-sites/"
Write(2,*) "mkdir ",depot//"/tgauges/tgauges-scplots/"
Write(2,*) "mkdir ",depot//"/tgauges/tgauges-predictions/"
Write(2,*) "mkdir ",depot//"/gmaps/"
Write(2,*) "mkdir ",depot//"/geod/"
Write(2,*) "mkdir ",depot//"/geod/sites/"
Write(2,*) "mkdir ",depot//"/geod/3dmaps/"
Write(2,*) "mkdir ",depot//"/geod/3dmaps/northamerica/"
Write(2,*) "mkdir ",depot//"/geod/3dmaps/fennoscandia/"
Write(2,*) "mkdir ",depot//"/geod/3dmaps/mediterranean/"
Write(2,*) "mkdir ",depot//"/rmaps/"
Write(2,*) "mkdir ",depot//"/rmaps/sun_data"
Write(2,*) "mkdir ",depot//"/rmaps/Italy"
Write(2,*) "mkdir ",depot//"/rmaps/Mediterranean"
Write(2,*) "mkdir ",depot//"/rmaps/Europe"
Write(2,*) "mkdir ",depot//"/rmaps/Fennoscandia"
Write(2,*) "mkdir ",depot//"/rmaps/Greenland"
Write(2,*) "mkdir ",depot//"/rmaps/North_America"
Write(2,*) "mkdir ",depot//"/rmaps/Antarctica"
Write(2,*) "mkdir ",depot//"/stokes/"
Write(2,*) "else"
Write(2,*) "echo"
Write(2,*) "echo '+++>  WARNING/ a repository already exists with name: '", depot 
Write(2,*) "    "
Write(2,*) "echo '+++>  [existing data will be overwritten]'"
Write(2,*) "    "
Write(2,*) "fi"
!
!
! ============================== 
! --- Purging, if requested  ---
! ==============================
if(option_purge=='y') then
  Write(2,*) " "
  Write(2,*) "#echo --------------------------------------"
  Write(2,*) " echo"
  Write(2,*) " echo '---> Purging the working directory: '"  
  Write(2,*) "#echo --------------------------------------"
  Write(2,*) " echo"
  Write(2,*) "/bin/rm *.grd"
  Write(2,*) "/bin/rm *.jun"  
  Write(2,*) "/bin/rm junk*"  
  Write(2,*) "/bin/rm *.mod" 
  Write(2,*) "/bin/rm *.obj" 
  Write(2,*) "/bin/rm *.old"
  Write(2,*) "/bin/rm *.bak"
Endif
!
! ==========================================================
! --- Setting the number of threads for OpenMP execution ---
! ==========================================================
!
if( (option_omp=='y') ) then
  Write(2,*) " "
  Write(2,*) "#echo --------------------------------------------------------"
  Write(2,*) "echo                                                          "
  Write(2,*) " echo '---> Number of threads for OpenMP execution: '", nthread 
  Write(2,*) "#echo --------------------------------------------------------"
  Write(2,*) " echo"
  Write(2,*) "export OMP_NUM_THREADS=",trim(nthread)
Endif
!
! ==========================================================
! --- Adjusting the unformatted record length on IBM XLF ---
! ==========================================================
!
if(option_sys=='6') then
  Write(2,*) " "
  Write(2,*) "#echo --------------------------------------------------------"
  Write(2,*) "echo                                                          "
  Write(2,*) " echo '---> Setting XLFRTEOPTS '"
  Write(2,*) "#echo --------------------------------------------------------"
  Write(2,*) " echo"
  Write(2,*) 'export XLFRTEOPTS="uwidth=64"'
Endif
!
! ==========================================================
! --- Setting the MP_PROCS variable for MPI jobs on IBM XLF ---
! ==========================================================
!
if( (option_sys=='6') .and. (option_mpi=='y') )then
  Write(2,*) " "
  Write(2,*) "#echo --------------------------------------------------------"
  Write(2,*) "echo                                                          "
  Write(2,*) " echo '---> Setting MP_PROCS to '", ntask
  Write(2,*) "#echo --------------------------------------------------------"
  Write(2,*) " echo"
  Write(2,*) 'export MP_PROCS=',trim(ntask)
Endif
!
!
!
! ====================================
! EXE 01 --- Pixelizing the sphere --- 
! ====================================
!
Write(2,*) " "
Write(2,*) "#echo --------------------------------------------------------"
Write(2,*) "echo                                                          "
Write(2,*) " echo '---> PX.F90: Hicosahedral pixelization of the sphere  '"
Write(2,*) "#echo --------------------------------------------------------"
Write(2,*) "px.exe"
Write(2,*) "cp px.dat ", depot//"/px"
Write(2,*) "cp pxa.dat ", depot//"/px"
Write(2,*) "cp px-lat.dat ", depot//"/px"
!
!
! ========================================================
! EXE 02 --- Pixelizing the PRESENT-DAY ocean function ---
! ========================================================
!
 Write(2,*) " "
 if( option_mpi=='y' ) then
    Write(2,*) "#echo ---------------------------------------"
    Write(2,*) "echo" 
    Write(2,*) " echo '---> pxselect.exe: Separating wet from dry pixels'"
    Write(2,*) "#echo ---------------------------------------"
 else
    Write(2,*) "#echo ---------------------------------------"
    Write(2,*) "echo" 
    Write(2,*) " echo '---> px.gmt: Separating wet from dry pixels'"
    Write(2,*) "#echo ---------------------------------------"
 endif 
!
!--- A realistic Ocean Function
 if(option_rof=='r') then 
	     Write(2,*) "#echo -------------------------------------"
             Write(2,*) " echo '     - Realistic ocean function'"
	     Write(2,*) "#echo -------------------------------------"
             Call OF_pixelization_real
 Endif
! 
!--- A zonal Ocean Function
 if(option_rof=='z') then 
 CONTINUE
 Endif
!
 if( option_mpi=='y' ) then
    Write(2,*) trim(RunCmd)//" pxselect.exe"
 else
    Write(2,*) "sh px.gmt"
 endif
 Write(2,*) "cp weta.dat ", depot//"/px" 
 Write(2,*) "cp drya.dat ", depot//"/px" 
! 
 Write(2,*) " "
 Write(2,*) "#echo ---------------------------------------"
 Write(2,*) "echo" 
 Write(2,*) " echo '---> PX_REC.F90: Merging the wet & dry pixels tables'"
 Write(2,*) "#echo ---------------------------------------" 
 Write(2,*) "pxrec.exe"
 Write(2,*) "cp px-table.dat ", depot//"/px" 
! 
 If(option_mpi=='y') then
    Write(2,*) " "
    Write(2,*) "#echo ---------------------------------------"
    Write(2,*) "echo" 
    Write(2,*) " echo '---> PX_PART.F90: Partitioning the pixels for parallel execution'"
    Write(2,*) "#echo ---------------------------------------" 
    Write(2,*) "pxpart.exe"
    Write(2,*) "cp px-partition.dat ", depot//"/px" 
 EndIf
!
 If(option_nls=='y') then
    Write(2,*) " "
    Write(2,*) "#echo ---------------------------------------"
    Write(2,*) "echo" 
    Write(2,*) " echo '---> Copying pixel files to node local storage'"
    Write(2,*) "#echo ---------------------------------------" 
    Write(2,*) trim(RunCmd)," pxcopy.exe"
 Endif
!
 Write(2,*) "mv px.gmt ", depot//"/px"
!
!
	If(option_px=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo -----------------------------------"
 	Write(2,*) "echo" 	
        Write(2,*) " echo '---> pxmap.gmt: Producing pixelization maps'"
	Write(2,*) "#echo -----------------------------------"
		file_gmt="pxmap.gmt"
		Call make_pxmap (nresolution, file_gmt)
		if(option_gmt=='y') Write(2,*) "sh pxmap.gmt"
			if(option_gmt=='y') Write(2,*) "ps2pdf px.ps"
			if(option_gmt=='y') Write(2,*) "ps2pdf px-sphere.ps"			
			if(option_gmt=='y') Write(2,*) "mv px.pdf ", depot//"/px"
			if(option_gmt=='y') Write(2,*) "mv px.ps ", depot//"/px"
			if(option_gmt=='y') Write(2,*) "ps2pdf px-sphere.ps"			
			if(option_gmt=='y') Write(2,*) "mv px-sphere.pdf ", depot//"/px"
			if(option_gmt=='y') Write(2,*) "mv px-sphere.ps ", depot//"/px"
			Write(2,*) "mv ", trim(adjustl(file_gmt)), " ", depot//"/px"
	Endif	
!
!
! ==================================================
! EXE 03 --- Computing the spherical harmonics  --- 
! ==================================================
! 
If(option_sh=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo ----------------------------------------"
	Write(2,*) "echo" 
	Write(2,*) " echo '---> SH.F90: Building the spherical harmonics'"
	Write(2,*) "#echo ----------------------------------------"
	Write(2,*) trim(RunCmd)," sh.exe"
	Write(2,*) "cp sh.bin ", trim(adjustl(sh_file)) 
!
	           else
	Write(2,*) " "		   
	Write(2,*) "#echo ------------------------------------------------------"
	Write(2,*) "echo" 	
	Write(2,*) " echo '+++> A SH file already exists with name:' ", trim(adjustl(sh_file)) 
	Write(2,*) "#echo ------------------------------------------------------"	
	Write(2,*) "cp ", trim(adjustl(sh_file)), " sh.bin"
!
endif	    
!
!	 
! =========================================================
! EXE 04 ---  Evaluating & plotting the WINDOW function ---
! =========================================================
!
If(option_wi=='y') then
	Write(2,*) " " 
	Write(2,*) "#echo -----------------------------------------------------"
	Write(2,*) "echo" 
	Write(2,*) " echo '--->  WNW.F90: SH orthonormality evaluating the window function'"
	Write(2,*) "#echo -----------------------------------------------------"
	Write(2,*) "wnw.exe"
	file_gmt="wnw.gmt"
	call make_wnw (resolution, degree, file_gmt)
!
	If(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))		
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/wnw"
		Write(2,*) "mv wnw.dat ", depot//"/wnw"
		If(option_gmt=='y') Write(2,*) "ps2pdf wnw.ps "
		If(option_gmt=='y') Write(2,*) "mv wnw.pdf ", depot//"/wnw"
		If(option_gmt=='y') Write(2,*) "mv wnw.ps ", depot//"/wnw"
		 Write(2,*) "mv tmpp4.dat ", depot//"/wnw"
		 Write(2,*) "mv tmpp5.dat ", depot//"/wnw"
		 Write(2,*) "mv tmpp6.dat ", depot//"/wnw"
		 Write(2,*) "mv tmpp7.dat ", depot//"/wnw"
		 Write(2,*) "/bin/rm tmpw0.dat"
		 Write(2,*) "/bin/rm tmpw1.dat"
		 Write(2,*) "/bin/rm tmpw2.dat"
		 Write(2,*) "/bin/rm tmpw3.dat"
!
Endif
!
!
! ===================================================================
! EXE 05 --- Evaluating the SH coefficients of the PRESENT DAY OF ---
! ===================================================================
!
If(option_oh=='y') then 
!
! A new SH dechomposition of the present-day OF is performed 
	Write(2,*) " "
	Write(2,*) "#echo ------------------------------------------"
	Write(2,*) "echo                                         " 
	Write(2,*) " echo '---> SH_OF.F90: SH expansion of the present-day ocean function'"
	Write(2,*) "#echo ------------------------------------------"
	Write(2,*) trim(RunCmd)," shof.exe"
	Write(2,*) "cp shof.dat ", trim(adjustl(shof_file)) 
	Write(2,*) "cp ", trim(adjustl(shof_file)), " ", depot//"/of"
		   else
!		   
! An existing SH dechomposition is employed  
!
	Write(2,*) " "
	Write(2,*) "#echo -----------------------------------------------------------"
	Write(2,*) "echo" 	
	Write(2,*) " echo '+++> An OF SH file already exists with name:' ", trim(adjustl(shof_file))
	Write(2,*) "#echo -----------------------------------------------------------" 
	Write(2,*) "cp ", trim(adjustl(shof_file)), " shof.dat"
	Write(2,*) "cp ", trim(adjustl(shof_file)), " ", depot//"/of"
Endif
!
!
! ==================================================
! EXE 06 --- Present-day OF Degree Variance (DV) ---
! ==================================================
!
If(option_ofdv=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo ---------------------------"
	Write(2,*) "echo"
	Write(2,*) " echo '---> Computing the present-day OF DV'"
	Write(2,*) "#echo ---------------------------"
!
! Computing the DV 
	Write(2,*) "ofdv.exe"  
!
! Producing a script for plotting the degree variance 
	file_gmt="ofdv.gmt"
 	Call make_ofdvmap (degree, file_gmt) 
!
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))
		if(option_gmt=='y') Write(2,*) "ps2pdf ofdv.ps"
		if(option_gmt=='y') Write(2,*) "mv ofdv.ps ",      depot//"/of/degree_variance"
		if(option_gmt=='y') Write(2,*) "mv ofdv.pdf ",     depot//"/of/degree_variance"						
		Write(2,*) "mv dv-power.tmp ", depot//"/of/degree_variance"
		Write(2,*) "mv ofdv.dat ",     depot//"/of/degree_variance"
		Write(2,*) "mv ofdv.gmt ",     depot//"/of/degree_variance"
Endif
!
!
! ============================================================
! EXE 07 --- Reconstruction and mapping the present-day OF ---
! ============================================================
!
If(option_of=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo -----------------------------------------------------"
	Write(2,*) "echo" 
	Write(2,*) " echo '---> REC_OF.F90: Reconstructing and mapping the ocean function'"
	Write(2,*) "#echo -----------------------------------------------------"
!
! Rechonstruction 
	Write(2,*) trim(RunCmd)," recof.exe"  
!		
! Mapping 	
	file_gmt="of.gmt"
 	Call make_ofmap (degree, option_rof, radius_zof, file_cap, file_gmt)    		
!
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
				    Write(2,*) "mv recof.dat ", depot//"/of"
                if(option_gmt=='y') Write(2,*) "ps2pdf of.ps"
!		
		If(option_gmt=='y'.and.option_rof=='z') Write(2,*) "cp ", trim(adjustl(file_cap)), " "//depot//"/of" 

		if(option_gmt=='y') Write(2,*) "mv of.ps ", depot//"/of"
		if(option_gmt=='y') Write(2,*) "mv of.pdf ", depot//"/of"
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/of"
!
		If(option_gmt=='y'.and.option_rof=='r') Write(2,*) "mv pale.cpt ", " "//depot//"/of"		  
Endif 
!
!!
! ===============================================
! --- Copying the ice file from the ice store --- 
! ===============================================
	Write(2,*) " "
        Write(2,*) "#echo --------------------------------------------------------------------"
	Write(2,*) "echo"   
        Write(2,*) " echo '---> Importing '", trim(adjustl(ice_file))," from ICE-MODELS/"   
        Write(2,*) "#echo --------------------------------------------------------------------"  	
	Write(2,*) "cp ../ICE-MODELS/"//trim(adjustl(ice_file)), " ", "./"//trim(adjustl(ice_file))
!
!
! ======================================
! EXE 09 --- Computing shape factors --- 
! ======================================
!
IF(option_sf=='y') THEN 
!
  Write(2,*) " "
  Write(2,*) "#echo ------------------------------------------------------------"  
  Write(2,*) "echo                                                                   " 
  Write(2,*) " echo '---> SHAPE_FACTORS.F90: Computing the shape factors for model: '", trim(adjustl(ice_file))   
  Write(2,*) "#echo ------------------------------------------------------------"  
  Write(2,*) "shapefactors.exe"
!
  Write(2,*) " "
  Write(2,*) "#echo -------------------------------------------"
  Write(2,*) "echo                                                                "   
  Write(2,*) " echo '---> SHICE.F90: Computing SH coefficients for the ice model'"  
  Write(2,*) "#echo -------------------------------------------"  
  Write(2,*) "shice.exe" 
!
  Write(2,*) "cp shice.dat ", depot//"/"//trim(adjustl(titlice))//"/sh/"
!
  Write(2,*) "cp shice.dat ", trim(adjustl(shape_file))
!  
  If(option_ri=='n') Write(2,*) "mv sht*.dat ", depot//"/"//trim(adjustl(titlice))//"/sh/"
!
  		   ELSE
!
   Write(2,*) " "
   Write(2,*) "#echo ------------------------------------------------------------" 
   Write(2,*) "echo                                                                    "  
   Write(2,*) " echo '---> Ice harmonics are pre-computed in file: '", shape_file
   Write(2,*) "#echo ------------------------------------------------------------"  
!
   Write(2,*) "cp ", trim(adjustl(shape_file)), " shice.dat"
!
  If(option_ri=='y') then 
   Write(2,*) "#echo ------------------------------------------------------------" 
   Write(2,*) "echo                                                              "  
   Write(2,*) " echo '---> Copying from '", depot, " the ice harmonics for rechonstruction"
   Write(2,*) "#echo ------------------------------------------------------------"  
!  
  Write(2,*) "cp  ", depot//"/"//trim(adjustl(titlice))//"/sh/"//"sht*.dat  ."
!
  Endif
!
  ENDIF 
!
!
! ================================================== 
! EXE 12 --- Mapping of the original ice sheets ---
! ==================================================
!
If(option_or=='y') then 
!
Write(2,*) " "
Write(2,*) "#echo ---------------------------------------------" 
Write(2,*) "echo"    
Write(2,*) "echo '---> MS.F90: Creating multi-segment files for ice sheets maps'"
Write(2,*) "#echo ---------------------------------------------"  
Write(2,*) "ms.exe" 
!
!
Write(2,*) " "
Write(2,*) "#echo ------------------------------------------------"  
Write(2,*) "echo"   
Write(2,*) " echo '---> mapice.gmt: Creating ps images of original ice sheets'"
Write(2,*) "#echo ------------------------------------------------"  	   	
!	   	
!--- Creates a "rainbow palette" suitable for plotting the 
!    rechonstructed ice distribution, named "ice-pal.cpt": the 
!    same palette is also used for the rechonstructed ice 
!    distribution...(see below)
!
!
   file_gmt="mapice.gmt"
   Call make_icemap (ninc, titlice, option_rof, FILE_CAP, file_gmt)
!	
   if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))
   		       Write(2,*) "mv msg* ", depot//"/"//trim(adjustl(titlice))//"/original/"
   		       Write(2,*) "mv rtmp*.dat ", depot//"/"//trim(adjustl(titlice))//"/original/"   
   if(option_gmt=='y') Write(2,*) "mv mapice*.ps ", depot//"/"//trim(adjustl(titlice))//"/original/"
   if(option_gmt=='y') Write(2,*) "mv mapice*.pdf ", depot//"/"//trim(adjustl(titlice))//"/original/"
   if(option_gmt=='y') Write(2,*) "mv pice*.cpt ", depot//"/"//trim(adjustl(titlice))//"/original/"
   Write(2,*) "mv ",  trim(adjustl(file_gmt)), & 
   	       " "//depot//"/"//trim(adjustl(titlice))//"/original/"
   If(option_rof=='z'.and.option_gmt=='y')& 
      Write(2,*)"cp ", trim(adjustl(file_cap)), " ", & 
      depot//"/"//trim(adjustl(titlice))//"/original/"
!
ENDIF
!
!
! ============================================================
! EXE 13---  Rechonstruction and mapping of the ice sheets ---
! ============================================================
!
If(option_ri=='y') then 
!
Write(2,*) " "
Write(2,*) "#echo -------------------------------------------------"  
Write(2,*) "echo"   
Write(2,*) " echo '---> REC_ICE.F90: Reconstruction of ice sheets distribution'"
Write(2,*) "#echo -------------------------------------------------"  
!
! --- Rechonstruction ...
!
   Write(2,*) trim(RunCmd)," recice.exe"
!	 
! --- ... and mapping of rechonstruction 
! 
Write(2,*) " "
Write(2,*) "#echo ---------------------------------------------------------"  
Write(2,*) "echo"   
Write(2,*) " echo '---> recice.gmt: Creating ps images of reconstructed ice sheets'"
Write(2,*) "#echo ---------------------------------------------------------"  
!

!	   	
!--- Creates a "rainbow palette" suitable for plotting the 
!    rechonstructed ice distribution, named "ice-pal.cpt": the 
!    same palette is also used for the original ice sheets
!    distribution...(see above). 
!
!
   File_gmt="recice.gmt"	
   Call make_recicemap (ninc, degree, titlice, option_rof, file_cap, file_gmt)
!
   if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))
   Write(2,*) "mv rect*.dat ", depot//"/"//trim(adjustl(titlice))//"/rechonstructed/"
   Write(2,*) "mv tmp0*.dat ", depot//"/"//trim(adjustl(titlice))//"/rechonstructed/"
   if(option_gmt=='y') Write(2,*) "mv recice*.ps ", depot//"/"//trim(adjustl(titlice))//"/rechonstructed/"
   if(option_gmt=='y') Write(2,*) "mv recice*.pdf ", depot//"/"//trim(adjustl(titlice))//"/rechonstructed/"
   if(option_gmt=='y') Write(2,*) "mv pice*.cpt ", depot//"/"//trim(adjustl(titlice))//"/rechonstructed/"
   Write(2,*) "mv ", trim(adjustl(file_gmt)), & 
              " "//depot//"/"//trim(adjustl(titlice))//"/rechonstructed/" 
!
   Write(2,*) "mv sht*.dat ", depot//"/"//trim(adjustl(titlice))//"/sh/"            
!
ENDIF
!
!
!
!
! =========================================
! EXE 14 --- Computing the Love numbers ---
! =========================================
! 
!************************
IF(OPTION_NM=='y') THEN 
!************************
!
Write(2,*) "#echo  ---------------------------------------------------------------------" 
Write(2,*) "echo " 
Write(2,*) " echo  '---> TB.F90: Load-deformation coefficients by TABOO (Normal Modes)'"
Write(2,*) "#echo  ---------------------------------------------------------------------"  
Write(2,*) "tb.exe"         
Write(2,*) "echo " 
Write(2,*) "echo '+++> WARNING/ in SELEN 2.9, the average density of the Earth is computed by the'"
Write(2,*) "echo '            / density structure of the input model, NOT using an a-priori'"
Write(2,*) "echo '            / value as done in previous versions - GS & FC July 27 2009 - '"
!
!
Write(2,*) "/bin/rm visco.tmp"
!
Write(2,*) "cp ", trim(adjustl(visco_file)), " ", depot//"/TABOO/"
Write(2,*) "cp ", trim(adjustl(visco_file)), " ", depot//"/Love-Numbers-by-TABOO/"  
!
Write(2,*) "cp task_1.dat ", depot//"/log/"
Write(2,*) "mv task_1.dat ", depot//"/TABOO/"
Write(2,*) "mv spectrum.dat ", depot//"/Love-Numbers-by-TABOO/"
Write(2,*) "mv ih.dat ", depot//"/Love-Numbers-by-TABOO/"
Write(2,*) "mv ik.dat ", depot//"/Love-Numbers-by-TABOO/"
Write(2,*) "mv il.dat ", depot//"/Love-Numbers-by-TABOO/"
Write(2,*) "mv h.dat ", depot//"/Love-Numbers-by-TABOO"
Write(2,*) "mv l.dat ", depot//"/Love-Numbers-by-TABOO"
Write(2,*) "mv k.dat ", depot//"/Love-Numbers-by-TABOO"
!
If(option_ln=='y') then 
Write(2,*) " "
Write(2,*) "#echo -----------------------------------------------------------"  
Write(2,*) "echo "
Write(2,*) " echo '---> ldcs.gmt: Plotting Love numbers and other spectral quantities'"
Write(2,*) "#echo -----------------------------------------------------------"  
!
	File_gmt="ldcs.gmt"
	CALL make_plot_ldc (nv, code, degree, vstring, file_gmt)
!
        if(option_gmt=='y') Write(2,*) "sh ", file_gmt
	Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ss.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ihh.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ikk.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ill.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv hh.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv kk.dat ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv ll.dat ", depot//"/Love-Numbers-by-TABOO/"
	if(option_gmt=='y') Write(2,*) "mv h.tmp ", depot//"/Love-Numbers-by-TABOO/"
	if(option_gmt=='y') Write(2,*) "mv k.tmp ", depot//"/Love-Numbers-by-TABOO/"
	if(option_gmt=='y') Write(2,*) "mv l.tmp ", depot//"/Love-Numbers-by-TABOO/"
	if(option_gmt=='y') Write(2,*) "mv spe.tmp ", depot//"/Love-Numbers-by-TABOO/"
	Write(2,*) "mv tmpg*.dat ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y')Write(2,*) "ps2pdf ela-flu.ps"
	if(option_gmt=='y')Write(2,*) "ps2pdf spectrum.ps"
	if(option_gmt=='y')Write(2,*) "ps2pdf n-residues.ps"
	if(option_gmt=='y')Write(2,*) "mv ela-flu.ps ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y')Write(2,*) "mv spectrum.ps ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y')Write(2,*) "mv n-residues.ps ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y')Write(2,*) "mv ela-flu.pdf ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y')Write(2,*) "mv spectrum.pdf ", depot//"/Love-Numbers-by-TABOO"
	if(option_gmt=='y')Write(2,*) "mv n-residues.pdf ", depot//"/Love-Numbers-by-TABOO"
!
Endif
!
!****************************
ELSEIF(OPTION_PW=='y') THEN 
!****************************
 continue
!****************************
ELSEIF(OPTION_PWA=='y') THEN 
!****************************
 continue
!******
ENDIF 
!******
!
!
!
!
!
! ================================================================================
! --- Solving the Sea Level Equation (SLE) in the USUAL way (FIXED coastlines) ---
! ================================================================================
!
Write(2,*) " "
Write(2,*) "#echo  ------------------------------------------------------------------- " 
Write(2,*) "echo " 
Write(2,*) " echo  '---> SLE.F90: Solving the Sea Level Equation - SLE - for FIXED coastlines'"
Write(2,*) "#echo  ------------------------------------------------------------------- "  
Write(2,*) trim(RunCmd)," sle.exe" 
!




!
! ========================================================= 
! EXE 11 --- Plotting Equivalent Sea Level (ESL) curves ---
! =========================================================
!
If(option_esl=='y') then 
!

    Write(2,*) " "
    Write(2,*) "#echo -------------------------------"  
    Write(2,*) "echo "   
    Write(2,*) " echo '---> ESL.F90: Plotting the ESL curve'"
    Write(2,*) "#echo -------------------------------"  
    Write(2,*) "esl.exe" 
!
   file_gmt='eslplot.gmt'
   Call make_eslplot (ninc, titlice, shof_file, file_gmt)
!
   if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt))
   if(option_gmt=='y') Write(2,*) "ps2pdf esl.ps"
   Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   if(option_gmt=='y') Write(2,*) "mv esl.ps ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   if(option_gmt=='y') Write(2,*) "mv esl.pdf ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   Write(2,*) "mv esl.dat ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   if(option_gmt=='y') Write(2,*) "mv esl.tmp ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   Write(2,*) "mv esl-thin.dat ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
   Write(2,*) "mv esl-tot.dat ", " "//depot//"/"//trim(adjustl(titlice))//"/esl/"
!
Endif 
!
!
! ========================================
! EXE 12 ---  Relative Sea Level (RSL) ---
! ========================================
!
! RSL database 
If(trim(adjustl(option_rsla))=='y') then
	Write(2,*) " "
	Write(2,*) "#echo ------------------------------------------"
	Write(2,*) "echo"
	Write(2,*) " echo '---> The RSL database is: '", trim(adjustl(rsl_file)) 
	Write(2,*) "#echo ------------------------------------------"
ENDIF
!
! Mapping the distribution of RSL sites ... 
	If(option_rsldb=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo ------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> rsl-s.gmt: Map of RSL sites'"
		Write(2,*) "#echo ------------------------"	
		file_gmt='rsl-s.gmt'
		call MAKE_RSLDB(RSL_FILE, RSL_DATABASE_FORMAT, NRSL, FILE_GMT)
		if(option_gmt=='y') write(2,*) "sh ", file_gmt 
		Write(2,*) "cp ", trim(adjustl(rsl_file)), " "//depot//"/rsl/rsl-sites"		
		Write(2,*) "mv rsl-s.gmt ", " "//depot//"/rsl/rsl-sites"
		if(option_gmt=='y') Write(2,*) "ps2pdf maprsl.ps"
		if(option_gmt=='y') Write(2,*) "mv maprsl.ps ", " "//depot//"/rsl/rsl-sites"
		if(option_gmt=='y') Write(2,*) "mv maprsl.pdf ", " "//depot//"/rsl/rsl-sites"
		Write(2,*) "mv lon-lat-rsl.dat ", " "//depot//"/rsl/rsl-sites"
		Write(2,*) "mv tmptitle ", " "//depot//"/rsl/rsl-sites"	
	ENDIF
!
! --- Computing synthetic RSL curves ... 
	If(option_rsl=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo  ---------------------------------------------------------"
		Write(2,*) "echo"	
		Write(2,*) " echo '---> RSL.F90: Predicting RSL at the sites of database: '", trim(rsl_file)
		Write(2,*) "#echo  ---------------------------------------------------------"	
	        Write(2,*) "shrsl.exe"
		Write(2,*) "rsl.exe" 		
!Write(2,*) "mv shrsl.bin ", depot//"/bin"
!
! --- ... and drawing RSL figures for each site		
		If(option_rslp=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo  ------------------------------------------"
		Write(2,*) "echo"		
		Write(2,*) " echo '---> rsl-curves.gmt: Postscript images of RSL data vs predictions'"
		Write(2,*) "#echo  ------------------------------------------"		
!		
		file_gmt="rsl-curves.gmt"
		CALL MAKE_RSL (NV, CODE, RSL_FILE, RSL_DATABASE_FORMAT, RUN, NINC, & 
		               NRSL, TITLICE, RESOLUTION, ITER, MODE, DEGREE, VSTRING, FILE_GMT, SHORT_VISCO_FILENAME)
		if(option_gmt=='y') Write(2,*) 'sh ', trim(adjustl(file_gmt))
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/rsl/rsl-curves/"
		Write(2,*) "mv tmpb* ", depot//"/rsl/rsl-curves/"
	        if(option_gmt=='y') Write(2,*) "mv rslp*.ps ", depot//"/rsl/rsl-curves/ps"
	        if(option_gmt=='y') Write(2,*) "mv rslp*.pdf ", depot//"/rsl/rsl-curves/pdf"
!		Write(2,*) "/bin/rm junky*"
		Endif
!
	        Write(2,*) "mv rsld*.dat ", depot//"/rsl/rsl-curves/"
	        Write(2,*) "mv rslp*.dat ", depot//"/rsl/rsl-curves/"
!
	Write(2,*) "cp ", trim(adjustl(rsl_file)), " "//depot//"/rsl/rsl-curves/"
!
	If(option_rslsca=='n') then 
		Write(2,*) "mv scatter-data.dat ",  " "//depot//"/rsl/rsl-scplot/"
		Write(2,*) "mv scatter-pred.dat ",  " "//depot//"/rsl/rsl-scplot/"
	endif
!
	ENDIF 
!
! --- A scatterplot of all RSL data vs all observations...  
	If(option_rslsca=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo  ---------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> Drawing a RSL scatterplot'"
		Write(2,*) "#echo  ---------------------------------"		
!		
		file_gmt="rsl-sca.gmt"
		CALL MAKE_RSLSCA (NV, CODE, RSL_FILE, RUN, NINC, NRSL, TITLICE, RESOLUTION, & 
			          ITER, MODE, DEGREE, VSTRING, FILE_GMT, SHORT_VISCO_FILENAME)
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
		if(option_gmt=='y') Write(2,*) "ps2pdf scatter-plot.ps"
		if(option_gmt=='y') Write(2,*) "mv scatter-plot.ps ",  depot//"/rsl/rsl-scplot/"
		if(option_gmt=='y') Write(2,*) "mv scatter-plot.pdf ", depot//"/rsl/rsl-scplot/"
		Write(2,*) "mv scatter-data.dat ", depot//"/rsl/rsl-scplot/"
		Write(2,*) "mv scatter-pred.dat ", depot//"/rsl/rsl-scplot/"
		Write(2,*) "mv title.tmp ", 	   depot//"/rsl/rsl-scplot/"
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " ", depot//"/rsl/rsl-scplot/"
	ENDIF 
!
! --- Misfit study for RSL 
	If(option_rslmf=='y') then
		Write(2,*) " "
		Write(2,*) "#echo  -------------------------------"
		Write(2,*) "echo"		
		Write(2,*) " echo '---> Misfit analysis for RSL'"	
		Write(2,*) "#echo  -------------------------------"		
		file_gmt='rsl-mis.gmt'
        	Call MAKE_RSLMIS (NV, CODE, RUN, NINC, NRSL, TITLICE, RESOLUTION, ITER, & 
			          MODE, DEGREE, VSTRING, FILE_GMT, SHORT_VISCO_FILENAME)
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
		if(option_gmt=='y') Write(2,*) "ps2pdf rsl-misfit.ps"
		if(option_gmt=='y') Write(2,*) "mv rsl-misfit.ps ",  depot//"/rsl/rsl-misfit/"
		if(option_gmt=='y') Write(2,*) "mv rsl-misfit.pdf ", depot//"/rsl/rsl-misfit/"
		Write(2,*) "mv rsl-mis.gmt ",    depot//"/rsl/rsl-misfit/"
		Write(2,*) "mv tmph.dat ",       depot//"/rsl/rsl-misfit/"
		Write(2,*) "mv mis.dat ",        depot//"/rsl/rsl-misfit/"
		Write(2,*) "mv gmis.dat ",       depot//"/rsl/rsl-misfit/"
	ENDIF 
!
! --- RSL table with all observations & predictions   
	If(option_rsltab=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo  ---------------------------------------------"
		Write(2,*) "echo"		
		Write(2,*) " echo '---> RSL data and predictions table'"	
		Write(2,*) "#echo  ---------------------------------------------"
		Write(2,*) "mv rsl-obs-vs-predictions.dat ",  depot//"/rsl/rsl-table/"
	ENDIF
!
! --- <<Global RSL zones>> 
	If(option_rslz=='y') then
		Write(2,*) " "
		Write(2,*) "#echo  -------------------------------------"
		Write(2,*) "echo"					            
		Write(2,*) " echo '---> rsl-zones.gmt: Global RSL zones'"
		Write(2,*) "#echo  -------------------------------------"				    
		file1_gmt='rsl-zones.gmt'
		file2_gmt='rsl-allzones.gmt'
        	Call MAKE_RSLZONES (NV, CODE, RUN, NINC, NRSL, TITLICE, RESOLUTION, ITER, & 
			            MODE, DEGREE, VSTRING, FILE1_GMT, FILE2_GMT, SHORT_VISCO_FILENAME)
		Write(2,*) "rslzones.exe"
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file1_gmt)) 
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file2_gmt)) 
		Write(2,*) "mv ", trim(adjustl(file1_gmt)), " "//depot//"/rsl/rsl-zones/" 
		Write(2,*) "mv ", trim(adjustl(file2_gmt)), " "//depot//"/rsl/rsl-zones/" 
		Write(2,*) "mv tmpz*", " "//depot//"/rsl/rsl-zones/" 				
		Write(2,*) "mv lonlat-*.dat ", " "//depot//"/rsl/rsl-zones/" 
		Write(2,*) "mv rslzones-*.dat ", " "//depot//"/rsl/rsl-zones/" 
		if(option_gmt=='y') Write(2,*) "mv rslzones-*.pdf ", " "//depot//"/rsl/rsl-zones/" 		
		if(option_gmt=='y') Write(2,*) "mv rslzones-*.ps ", " "//depot//"/rsl/rsl-zones/" 	
	ENDIF 
!
! --- Regional RSL contour lines... 
	If(option_rslc=='y') then
		Write(2,*) " "
		Write(2,*) "#echo  ----------------------------------"
		Write(2,*) "echo"					            			    
		Write(2,*) " echo '---> RSLC.F90: Regional RSL contour lines'"
		Write(2,*) "#echo  ----------------------------------"
	        Write(2,*) "shrslc.exe"
		Write(2,*) "rslc.exe" 
		file_gmt='rslc.gmt'
		Call MAKE_RSLC (TIME_BPCC, MIN_RSLC, MAX_RSLC, RSL_INT, & 
				LONMINC, LONMAXC, LATMINC, LATMAXC, & 
			        NV, CODE, RUN, NINC, NRSLC, TITLICE, RESOLUTION, & 
				ITER, MODE, DEGREE, VSTRING, NAME_OF_REGION, &
				FILE_GMT, SHORT_VISCO_FILENAME) 
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 	
 		Write(2,*) "cp ", trim(adjustl(file_region)), " "//depot//"/rsl/rsl-contours/"
		Write(2,*) "mv lonlat_rslc.dat", " "//depot//"/rsl/rsl-contours/"	
		Write(2,*) "mv rslc-cont.dat", " "//depot//"/rsl/rsl-contours/"
		Write(2,*) "mv rslc.dat", " "//depot//"/rsl/rsl-contours/"
		Write(2,*) "/bin/rm shrslc.bin"
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/rsl/rsl-contours/"	
		if(option_gmt=='y') Write(2,*) "mv rslc-map.ps", " "//depot//"/rsl/rsl-contours/"
		Write(2,*) "mv rslc*.tmp", " "//depot//"/rsl/rsl-contours/"							  
		if(option_gmt=='y') Write(2,*) "mv pal_rslc.cpt", " "//depot//"/rsl/rsl-contours/"	
	ENDIF
!
!
! ==============================================
! EXE 13 --- Sea level change at tide-gauges ---
! ==============================================
!
! Tide gauge database 
If(option_tga=='y') then 
	Write(2,*) " "
	Write(2,*) "#echo -------------------------------------------------"
	Write(2,*) "echo"
	Write(2,*) " echo '---> The tide-gauge database is: '", trim(adjustl(TGAUGES_FILE)) 
	Write(2,*) "#echo -------------------------------------------------"
ENDIF
!
!
! Plotting map showing the distribution of tide-gauges...
	If(option_tgplot=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo ------------------------------------------------"
		Write(2,*) "echo"		
		Write(2,*) " echo '---> tgauges.gmt: Plotting the distribution of tide gauges'"
		Write(2,*) "#echo ------------------------------------------------"	
		file_gmt="tgauges.gmt"
		Call MAKE_TGAUGES (TGAUGES_FILE, FILE_GMT)
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 	
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/tgauges/tgauges-sites/"	
		Write(2,*) "cp ", trim(adjustl(TGAUGES_FILE)), " "//depot//"/tgauges/tgauges-sites/"
		Write(2,*) "mv lon-lat-tgauges-*.dat ",  depot//"/tgauges/tgauges-sites/"
		Write(2,*) "mv titlege*.tmp ", depot//"/tgauges/tgauges-sites/"
		Write(2,*) "mv titleall.tmp ", depot//"/tgauges/tgauges-sites/"
		Write(2,*) "mv titletgauges.tmp ", depot//"/tgauges/tgauges-sites/"
		if(option_gmt=='y') Write(2,*) "ps2pdf map-tgauges.ps "
		if(option_gmt=='y') Write(2,*) "mv map-tgauges.ps ", depot//"/tgauges/tgauges-sites/"
		if(option_gmt=='y') Write(2,*) "mv map-tgauges.pdf ", depot//"/tgauges/tgauges-sites/"
	Endif
!
!
! --- Plotting a scatterplot of observations vs years of observations 
	If(option_tgsca=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo -----------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> tgauges-scpl.gmt: Scatterplot of tide gauges trends'"
		Write(2,*) "#echo -----------------------------------------"	
		file1_gmt="tgauges-scpl.gmt"
		Call MAKE_TGAUGESSCA (TGAUGES_FILE, FILE1_GMT)
!		
! --- Scatterplots and stats 
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file1_gmt)) 
		if(option_gmt=='y') Write(2,*) "ps2pdf tgauges-scpl.ps"
		Write(2,*) "mv ", trim(adjustl(file1_gmt)), " "//depot//"/tgauges/tgauges-scplots/"	
		Write(2,*) "cp ", trim(adjustl(TGAUGES_FILE)), " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "mv tgauges-scpl.dat", " "//depot//"/tgauges/tgauges-scplots/"	
		if(option_gmt=='y') Write(2,*) "mv tgauges-scpl.ps", " "//depot//"/tgauges/tgauges-scplots/"
		if(option_gmt=='y') Write(2,*) "mv tgauges-scpl.pdf", " "//depot//"/tgauges/tgauges-scplots/"	
		if(option_gmt=='y') Write(2,*) "mv tgauges-stat.dat", " "//depot//"/tgauges/tgauges-scplots/"	
		if(option_gmt=='y') Write(2,*) "mv tgauges-scpl.tmp", " "//depot//"/tgauges/tgauges-scplots/"	
		Write(2,*) "mv tmpctitle1", " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "mv tmpctitle2", " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "mv tmpctitle3", " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "mv tmpctitle4", " "//depot//"/tgauges/tgauges-scplots/"
		Write(2,*) "mv tmpctitle5", " "//depot//"/tgauges/tgauges-scplots/"
	Endif	
!
! --- S-dot predictions at TG stations... 
	If(option_tg=='y') then 
		Write(2,*) " "
		Write(2,*) "#echo --------------------------------------------------"
		Write(2,*) "echo"		
		Write(2,*) " echo '---> TGAUGES.F90: Dot-S, U, and N predictions at tide gauges'"
		Write(2,*) "#echo --------------------------------------------------"
		Write(2,*) "shtgauges.exe"
		Write(2,*) "tgauges.exe"					
!Write(2,*) "mv shtidegauges.bin ", depot//"/bin/"
!Write(2,*) "cp shs.bin ", depot//"/bin/"	
		Write(2,*) "mv ptidegauges.dat ", " "//depot//"/tgauges/tgauges-predictions/"			
	Endif	
!
!
! ===============================================================
! ! EXE 14 ---  S, U, and N-dot at present time (Global maps) ---
! ===============================================================
!
	if(option_gm=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo -----------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> GMAPS.F90: Global maps of dot S, U and N at present time'"
		Write(2,*) "#echo -----------------------------------------------------"		
		file_gmt='gmaps.gmt'
		Write(2,*) trim(RunCmd)," gmaps.exe"  
        	Call MAKE_GMAPS (TITLICE, RESOLUTION, NV, CODE, ITER, MODE, & 
				 DEGREE, VSTRING, OPTION_ROF, FILE_CAP, FILE_GMT, SHORT_VISCO_FILENAME)
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/gmaps" 
		Write(2,*) "mv sdotmap.dat ", depot//"/gmaps/"
		Write(2,*) "mv udotmap.dat ", depot//"/gmaps/"
		Write(2,*) "mv ndotmap.dat ", depot//"/gmaps/"
		if(option_gmt=='y') Write(2,*) "ps2pdf sdotmap.ps"
		if(option_gmt=='y') Write(2,*) "ps2pdf udotmap.ps"
		if(option_gmt=='y') Write(2,*) "ps2pdf ndotmap.ps"		
		if(option_gmt=='y') Write(2,*) "mv sdotmap.ps ", depot//"/gmaps/"
		if(option_gmt=='y') Write(2,*) "mv udotmap.ps ", depot//"/gmaps/"
		if(option_gmt=='y') Write(2,*) "mv ndotmap.ps ", depot//"/gmaps/"
		if(option_gmt=='y') Write(2,*) "mv sdotmap.pdf ", depot//"/gmaps/"
		if(option_gmt=='y') Write(2,*) "mv udotmap.pdf ", depot//"/gmaps/"
		if(option_gmt=='y') Write(2,*) "mv ndotmap.pdf ", depot//"/gmaps/"
		Write(2,*) "mv tmpf* ", depot//"/gmaps/"
                if(option_gmt=='y') Write(2,*) "mv pale.cpt ", depot//"/gmaps/"			
	Endif
!
!
!
! ==============================================================
! EXE 15 --- S, U, and N-dot at present time (Regional maps) ---
! ==============================================================
!
	if(option_rm(0)=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo -------------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> RMAPS.F90: Regional maps of dot S, U and N at present time'"
		Write(2,*) "#echo -------------------------------------------------------"
		Write(2,*) trim(RunCmd)," rmaps.exe"  
		Write(2,*) " cp sdotmap.dat ", depot//"/rmaps/sun_data"
		Write(2,*) " cp udotmap.dat ", depot//"/rmaps/sun_data"
		Write(2,*) " cp ndotmap.dat ", depot//"/rmaps/sun_data"
        	Call MAKE_RMAPS (TITLICE, RESOLUTION, NV, CODE, ITER, MODE, DEGREE, VSTRING, & 
		                 OPTION_RM, SHORT_VISCO_FILENAME)
			if(option_rm(1)=='y') then 
			Write(2,*) " echo '     - Italy'"
			 if(option_gmt=='y') Write(2,*) "sh italy.gmt"  
			 if(option_gmt=='y') Write(2,*) "ps2pdf sdot-italy.ps"	
			 if(option_gmt=='y') Write(2,*) "ps2pdf udot-italy.ps"	
			 if(option_gmt=='y') Write(2,*) "ps2pdf ndot-italy.ps"
				Write(2,*) " mv italy.gmt ", depot//"/rmaps/Italy"
				Write(2,*) "mv *dot*ital* ", depot//"/rmaps/Italy"
				Write(2,*) "mv *pale*ital* ",depot//"/rmaps/Italy"		
			Endif			
			if(option_rm(2)=='y') then 
			Write(2,*) " echo '     - Mediterranean'"
			 if(option_gmt=='y')  Write(2,*) "sh mediterranean.gmt"       
			 if(option_gmt=='y')  Write(2,*) "ps2pdf sdot-mediterranean.ps"       
			 if(option_gmt=='y')  Write(2,*) "ps2pdf udot-mediterranean.ps"       
			 if(option_gmt=='y')  Write(2,*) "ps2pdf ndot-mediterranean.ps"
				Write(2,*) " mv mediterranean.gmt ", depot//"/rmaps/Mediterranean"
				Write(2,*) "mv *pale*medi* ",        depot//"/rmaps/Mediterranean"
				Write(2,*) "mv *dot*medi* ",         depot//"/rmaps/Mediterranean"			
			Endif			
			if(option_rm(3)=='y') then 
			Write(2,*) " echo '     - Europe'"
			 if(option_gmt=='y')  Write(2,*) "sh europe.gmt"  
			 if(option_gmt=='y')  Write(2,*) "ps2pdf sdot-europe.ps"      
			 if(option_gmt=='y')  Write(2,*) "ps2pdf udot-europe.ps"      
			 if(option_gmt=='y')  Write(2,*) "ps2pdf ndot-europe.ps"
				Write(2,*) " mv europe.gmt ", depot//"/rmaps/Europe"
				Write(2,*) "mv *pale*euro* ", depot//"/rmaps/Europe"
				Write(2,*) "mv *dot*euro* ",  depot//"/rmaps/Europe"
			Endif							
			if(option_rm(4)=='y') then 
			Write(2,*) " echo '     - Fennoscandia'"
			 if(option_gmt=='y')  Write(2,*) "sh fennoscandia.gmt"        
			 if(option_gmt=='y')  Write(2,*) "ps2pdf sdot-fennoscandia.ps"        
			 if(option_gmt=='y')  Write(2,*) "ps2pdf udot-fennoscandia.ps"        
			 if(option_gmt=='y')  Write(2,*) "ps2pdf ndot-fennoscandia.ps"
				Write(2,*) " mv fennoscandia.gmt ", depot//"/rmaps/Fennoscandia"
				Write(2,*) "mv *pale*fenn* ",       depot//"/rmaps/Fennoscandia"
				Write(2,*) "mv *dot*fenn* ",        depot//"/rmaps/Fennoscandia"
			Endif					
			if(option_rm(5)=='y') then 
			Write(2,*) " echo '     - Greenland'"
			 if(option_gmt=='y')  Write(2,*) "sh greenland.gmt"   
			 if(option_gmt=='y')  Write(2,*) "ps2pdf sdot-greenland.ps"   
			 if(option_gmt=='y')  Write(2,*) "ps2pdf udot-greenland.ps"   
			 if(option_gmt=='y')  Write(2,*) "ps2pdf ndot-greenland.ps"
				Write(2,*) " mv greenland.gmt ", depot//"/rmaps/Greenland"
				Write(2,*) "mv *pale*gree* ",    depot//"/rmaps/Greenland"
				Write(2,*) "mv *dot*gree* ",     depot//"/rmaps/Greenland"
			Endif					
			if(option_rm(6)=='y') then 
			Write(2,*) " echo '     - North America'"
			 if(option_gmt=='y')  Write(2,*) "sh north-america.gmt"       
			 if(option_gmt=='y')  Write(2,*) "ps2pdf sdot-namerica.ps"    
			 if(option_gmt=='y')  Write(2,*) "ps2pdf udot-namerica.ps"    
			 if(option_gmt=='y')  Write(2,*) "ps2pdf ndot-namerica.ps"
				Write(2,*) " mv north-america.gmt ", depot//"/rmaps/North_America"
			        Write(2,*) "mv *pale*name* ",        depot//"/rmaps/North_America"
				Write(2,*) "mv *dot*name* ",         depot//"/rmaps/North_America"
			Endif
			if(option_rm(7)=='y') then 
			Write(2,*) " echo '     - Antarctica'"
			 if(option_gmt=='y')  Write(2,*) "sh antarctica.gmt"	      
			 if(option_gmt=='y')  Write(2,*) "ps2pdf sdot-antarctica.ps"  
			 if(option_gmt=='y')  Write(2,*) "ps2pdf udot-antarctica.ps"  
			 if(option_gmt=='y')  Write(2,*) "ps2pdf ndot-antarctica.ps"
				Write(2,*) " mv antarctica.gmt ", depot//"/rmaps/Antarctica"
				Write(2,*) "mv *pale*anta* ",     depot//"/rmaps/Antarctica"
				Write(2,*) "mv *dot*anta* ",      depot//"/rmaps/Antarctica"
			Endif	
	Endif
!
!
!
! ======================================================================
! EXE 16 --- 3D velocity & S and N-dot at present time at specific sites 
! ======================================================================
!
	if(option_3d=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo -----------------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> GEO.F90: 3D velocity and S and N-dot today at specific sites'"
		Write(2,*) "#echo -----------------------------------------------------------"
		Write(2,*) "geo.exe"  
		Write(2,*) "cp ", file_3d, depot//"/geod/sites" 		
		Write(2,*) "mv geodetic-predictions.dat ", depot//"/geod/sites" 		
	Endif 
!
!
!
! ========================================================
! EXE 17 --- 3D velocity at present time at points on maps  
! ========================================================
!
	if(option_3d_regions=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo -----------------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> GEO_MAPS.F90: 3D velocity and U-dot today at points on maps'"
		Write(2,*) "#echo -----------------------------------------------------------"
		Write(2,*) "geomaps.exe"  
	Endif 
!
!
!
!
! =========================================================================
! EXE 19 ---  Rate of variation of Stokes coefficients at present time  ---
! =========================================================================
!
	if(option_st=='y') then 
		Write(2,*) " " 
		Write(2,*) "#echo ----------------------------------------------------------------"
		Write(2,*) "echo"
		Write(2,*) " echo '---> STOKES.F90: Present-time rate of change of Stokes coefficients'"
		Write(2,*) "#echo ----------------------------------------------------------------"		
		file_gmt='stokes.gmt'
		Write(2,*) "stokes.exe"  					   
                Call MAKE_STOKES (RUN, NV, CODE, TITLICE, RESOLUTION, ITER, & 
			          MODE, DEGREE, VSTRING, SHORT_VISCO_FILENAME, FILE_GMT)
		if(option_gmt=='y') Write(2,*) "sh ", trim(adjustl(file_gmt)) 
		Write(2,*) "mv ", trim(adjustl(file_gmt)), " "//depot//"/stokes/"  
		if(option_gmt=='y') Write(2,*) "ps2pdf stokes.ps"
		if(option_gmt=='y') Write(2,*) "mv stokes.ps ", depot//"/stokes/"
		if(option_gmt=='y') Write(2,*) "mv stokes.pdf ", depot//"/stokes/"
		Write(2,*) "mv stokes.dat ", depot//"/stokes/"
		if(option_gmt=='y') Write(2,*) "mv cosine.tmp ", depot//"/stokes/"
		if(option_gmt=='y') Write(2,*) "mv sine.tmp ", depot//"/stokes/"
		Write(2,*) "mv stokes*.tmp ", depot//"/stokes/"	
		Write(2,*) "mv title_stokes.tmp ", depot//"/stokes/"		
	Endif
!
!
!
!
!
Write(2,*) " "
Write(2,*) " "
Write(2,*) "echo ''"
Write(2,*) "echo '------------------------------------'"
Write(2,*) "echo ' >>> 2. Cleaning up the directory...'"
Write(2,*) "echo '------------------------------------'"
!
If(option_rm(0)=='y')  then
		Write(2,*) " mv sdotmap.dat ", depot//"/rmaps/sun_data"
		Write(2,*) " mv udotmap.dat ", depot//"/rmaps/sun_data"
		Write(2,*) " mv ndotmap.dat ", depot//"/rmaps/sun_data"
Endif
!
Write(2,*) "if [ -f  ./gmtdefaults4  ]" 
Write(2,*) "then"
Write(2,*) "/bin/rm -v ./gmtdefaults4"
Write(2,*) "fi"
!
Write(2,*) "/bin/rm -v *brok*.dat"
!
If(option_nm=='y') Write(2,*) "mv taboo.log ",    depot//"/log/"
Write(2,*) "mv selen.log ",    depot//"/log/"
Write(2,*) "cp config.dat ",   depot//"/log/"
Write(2,*) "cp config.f90 ",   depot//"/log/"
Write(2,*) "cp data.inc ",     depot//"/log/"
Write(2,*) "cp selen.sh ",     depot//"/log/"
!
Write(2,*) "mv pxa.dat ",      depot//"/px"
Write(2,*) "mv weta.dat ",     depot//"/px/"
Write(2,*) "mv drya.dat ",     depot//"/px/"
Write(2,*) "mv px-lat.dat ",   depot//"/px/"
Write(2,*) "mv px.dat ",       depot//"/px/"
Write(2,*) "mv px-table.dat ", depot//"/px/"
!
If    (option_nm=='y')then
! 
	Write(2,*) "mv ebu.dat ", depot//"/Love-Numbers-by-TABOO"
	Write(2,*) "mv ebv.dat ", depot//"/Love-Numbers-by-TABOO"
	Write(2,*) "mv ebs.dat ", depot//"/Love-Numbers-by-TABOO"
	Write(2,*) "mv ebn.dat ", depot//"/Love-Numbers-by-TABOO"
elseif(option_pw=='y')then  
        continue
Endif
!
Write(2,*) "mv ",  trim(adjustl(ice_file)), &
           "   "//depot//"/"//trim(adjustl(titlice))//"/original/"	   
!
If    (option_topo=='y')then
        continue	
Endif
!

! ============================== 
! --- Purging, if requested  ---
! ==============================
if(option_purge=='y') then
  Write(2,*) " "
  Write(2,*) "/bin/rm *.exe"
  Write(2,*) "/bin/rm *.grd"
  Write(2,*) "/bin/rm *.jun"  
  Write(2,*) "/bin/rm *.log"
  Write(2,*) "/bin/rm *.tmp" 
  Write(2,*) "/bin/rm *.cpt"
  Write(2,*) "/bin/rm *tmp*"
  Write(2,*) "/bin/rm junk*"  
  Write(2,*) "/bin/rm *.mod" 
  Write(2,*) "/bin/rm *.obj" 
  Write(2,*) "/bin/rm *.o"
  Write(2,*) "/bin/rm *.old"
  Write(2,*) "/bin/rm *.bak"
  Write(2,*) "/bin/rm fort*"
Endif
!
Write(2,*) "echo ''"
Write(2,*) "echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'"
Write(2,*) "echo '     SELEN, a Sea levEL EquatioN solver, Version 2.9     '"
Write(2,*) "echo '                      g95  VERSION                       '" 
Write(2,*) "echo '      Web page: http://fcolleoni.free.fr/SELEN.html      '"	
Write(2,*) "echo '   http://www.fis.uniurb.it/spada/SELEN_minipage.html    '"       
Write(2,*) "echo '   Send comments, requests of help and suggestions to:   '" 
Write(2,*) "echo '                <giorgio.spada@gmail.com>                '"
Write(2,*) "echo '                            -                            '"
Write(2,*) "echo '                    Copyright(C) 2008                    '"    
Write(2,*) "echo '     Giorgio Spada, Florence Colleoni & Paolo Stocchi    '"
Write(2,*) "echo '                          * * *                          '"
Write(2,*) "echo '     This programs comes with  ABSOLUTELY NO WARRANTY    '"
Write(2,*) "echo ' This is free software and you are welcome to distribute '"
Write(2,*) "echo '              it under certain conditions.               '"
Write(2,*) "echo '    For details, visit  <http://www.gnu.org/licenses/>   '"
Write(2,*) "echo '                  or edit file COPYING                   '"
Write(2,*) "echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'"
Write(2,*) "echo ''"
!
Write(2,*) "echo ''"
Write(2,*) "echo ' >>> Outputs for this run are available in directory: '", trim(adjustl(depot))
!
! --- Closing "selen.sh"
   close(2)
!
! --- Closing "config.dat"
   close(1)
!
!
!
!
! =========================================================
!
! 	Part #3 : preparing the include file 'data.inc'
!
! =========================================================
!
! 
! >>>>>> A time stamp on 'data.inc'
!
  call DATE_AND_TIME (date,timc)      
!
open(3,file='data.inc',status='unknown',recl=256)
!
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!", " File <<data.inc>>, created by <<config.f90>> on ", & 
             date(1:4), '.', date(5:6), '.', date(7:8), ' ', & 
	     timc(1:2), '.', timc(3:4), '.', timc(5:6) 
Write(3,*)"!                  				                                "
Write(3,*)"! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
Write(3,*)"! Copyright (C) 2008 Giorgio Spada, Florence Colleoni, and Paolo Stocchi     "
Write(3,*)"! 								                "								
Write(3,*)"! This file is part of SELEN.                                                " 
Write(3,*)"!       							                " 
Write(3,*)"! SELEN is free software: you can redistribute it and/or modify it under the "
Write(3,*)"! terms of the GNU General Public License as published by the Free Software  " 
Write(3,*)"! Foundation, either version 3 of the License, or at your option) any later  " 
Write(3,*)"! version. 									"
Write(3,*)"! 									        "
Write(3,*)"! SELEN is distributed in the /hope that it will be useful, but WITHOUT ANY  "
Write(3,*)"! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  "
Write(3,*)"! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      "
Write(3,*)"! details. 									"
Write(3,*)"!   										"
Write(3,*)"! You should have received a copy of the GNU General Public License along    "
Write(3,*)"! with SELEN.  If not, see <http://www.gnu.org/licenses/>.                   "
Write(3,*)"! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
Write(3,*) "!"
Write(3,*) "!------------------- 1- General settings --------------------------"  
Write(3,*) "!"
Write(3,*) "! --- Repository name"
Write(3,*) "CHARACTER*4, PARAMETER :: RUN=",trim(adjustl(runp))
if( option_mpi == 'y' ) then 
  Write(3,*) "! --- Number of CPUs"
  Write(3,*) "INTEGER, PARAMETER :: NCPU=", ntask
endif
Write(3,*) "! --- Prefix for PX and SH files"
if( option_nls == 'y' ) then
  Write(3,*) "CHARACTER*",len_tmp,", PARAMETER :: TMP_PREFIX=", tmp_prefix
else
  Write(3,*) "CHARACTER*2, PARAMETER :: TMP_PREFIX='./'"
endif
Write(3,*) "! --- Tegmark resolution"
Write(3,*) "INTEGER, PARAMETER :: RES=", resolution
Write(3,*) "! --- Number of pixels"
Write(3,*) "INTEGER, PARAMETER :: NP=2*RES*(RES-1)*20+12"  
Write(3,*) "! --- Maximum harmonic degree"
Write(3,*) "INTEGER, PARAMETER :: LMAX=", degree   
!
Write(3,*) "! --- Dealing with the degree 1 Love numbers"
If(option_deg1=='y') Write(3,*) "INTEGER, PARAMETER :: DEG1= 1"
If(option_deg1=='n') Write(3,*) "INTEGER, PARAMETER :: DEG1= 0"
!
If (option_rframe=='CM') Write(3,*) "CHARACTER*2, PARAMETER :: RFRAME= 'CM'" 
If (option_rframe=='CE') Write(3,*) "CHARACTER*2, PARAMETER :: RFRAME= 'CE'" 
!
Write(3,*) "! --- Jmax index"
Write(3,*) "INTEGER, PARAMETER :: JMAX=(LMAX+1)*(LMAX+2)/2" 
Write(3,*) "! --- Pi"
Write(3,*) "REAL*8,  PARAMETER :: PI=3.14159265358979323840" 
Write(3,*) "! --- Earth, water, and ice densities (kg/m^3)"
!
Write(3,*) "REAL*8,  PARAMETER :: RHOW= ", RHO_WATER
Write(3,*) "REAL*8,  PARAMETER :: RHOI= ", RHO_ICE 
!
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 2- Ice models settings -----------------------"  
Write(3,*) "!"
Write(3,*) "! --- Ice sheet model file"
Write(3,*) "CHARACTER*",len_ice,", PARAMETER :: ICE_MODEL=", ice_filename
Write(3,*) "! --- Number of time steps"
Write(3,*) "INTEGER, PARAMETER :: NN=", NINC		
Write(3,*) "! --- Number of ice elements"
write(3,*) "INTEGER, PARAMETER :: NEL=", nice
Write(3,*) "! --- Time increment, ka"
Write(3,*) "REAL*8, PARAMETER :: DELTA=", DELTA   
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
If    (option_nm=='y') then
Write(3,*) "!------------------- 3- Normal modes (TABOO) settings -------------"  
Write(3,*) "!"
Write(3,*) "! --- Earth model (see TABOO user guide)"
Write(3,*) "INTEGER, PARAMETER :: NV=", nv, ", CDE=", code
Write(3,*) "! --- Viscosity & Litho model file"
write(3,*) "CHARACTER*",LEN_VISCO,", PARAMETER :: VISCO_MODEL= &"
write(3,*)  trim(adjustl(visco_filename))
elseif (option_pw=='y') then 
 continue
elseif (option_pwa=='y') then 
 continue
Endif 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 4- Sea level Equation settings ---------------"  
Write(3,*) "!"
Write(3,*) "! --- Number of iterations"
Write(3,*) "INTEGER, PARAMETER :: SMAX=", iter
Write(3,*) "! --- Mode of solution"
Write(3,*) "INTEGER, PARAMETER :: IMODE=", mode
If(option_topo=='y') then
 CONTINUE
Endif 
!
!
If(option_rsl=='y'.or.option_rslp=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 5- RSL analysis at site of a database --------"
Write(3,*) "!"
Write(3,*) "! --- RSL database"
Write(3,*) "CHARACTER*", len_rsl, ", PARAMETER :: RSL_DATABASE=", rsl_database
Write(3,*) "! --- Format of the RSL database"
Write(3,*) "CHARACTER* 1, PARAMETER :: RSL_DATABASE_FORMAT = ", & 
		"'"//trim(adjustl(RSL_DATABASE_FORMAT))//"'"
Write(3,*) "! --- Number of RSL sites in database"
Write(3,*) "INTEGER, PARAMETER :: NRSL=", nrsl
Endif
!
If(option_rslc=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 6- Regional RSL analysis ---------------------"
Write(3,*) "!"
Write(3,*) "! --- Virtual RSL sites lon-lat file"
Write(3,*) "CHARACTER*", LEN_RSLC, ", PARAMETER :: RSLC_FILE=", RSLC_LONLAT_FILE
Write(3,*) "! --- Number of virtual RSL sites"
Write(3,*) "INTEGER, PARAMETER :: NRSLC=", nrslc
Write(3,*) "! --- Time of analysis (ka)"
Write(3,*) "REAL, PARAMETER :: TIME_BPC=", time_bpc
Endif
!
If(option_tg=='y'.or.option_tgplot=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 7- Predictions at tide-gauges ----------------"
Write(3,*) "!"
Write(3,*) "! --- Tide-gauge database"
Write(3,*) "CHARACTER*", LEN_TGAUGES, ", PARAMETER :: TGAUGES_DATABASE=", & 
					              TGAUGES_DATABASE
Write(3,*) "! --- Number of tide gauges"
Write(3,*) "INTEGER, PARAMETER :: NTIDEGAUGES=", NTIDEGAUGES
Endif 
!
If(option_3d=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 8- Predictions at geodetic sites ----------------"
Write(3,*) "!"
Write(3,*) "! --- Database of geodetic points"
Write(3,*) "CHARACTER*", LEN_GEOD, ", PARAMETER :: GEODETIC_DATABASE=", & 
						   GEODETIC_DATABASE
Write(3,*) "! --- Number of tide gauges"
Write(3,*) "INTEGER, PARAMETER :: NGEOD=", NGEOD
Endif 
!
If(option_3d_regions=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------- 9- Regional 3D geodetic predictions ----------------"
Write(3,*) "!"
Write(3,*) "! --- File with description of regions"
Write(3,*) "CHARACTER*", LEN_3D_REGIONS, ", PARAMETER :: TRED_REGIONS_DATABASE=", &
							 TRED_REGIONS_DATABASE
Write(3,*) "! --- Number of regions"
Write(3,*) "INTEGER, PARAMETER :: N_3D_REGIONS=", N_3D_REGIONS
Write(3,*) "! --- Names of regions"
Write(3,*) "CHARACTER*20 TRED_REGIONS_NAME(N_3D_REGIONS)"
do j=1, N_3D_REGIONS
        Write(3,*) "DATA TRED_REGIONS_NAME(", j, ")", "/"//trim(adjustl(TRED_REGIONS_NAME(j)))//"/" 
enddo
Endif
!
!
If(option_st=='y')then 
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!------------------ 12- Stokes coefficients (SC) settings ---------"
Write(3,*) "!"
Write(3,*) "! --- Min. degree for the SC"
Write(3,*) "INTEGER, PARAMETER :: STMIN=", degree_st_min
Write(3,*) "! --- Max. degree for the SC"
Write(3,*) "INTEGER, PARAMETER :: STMAX=", degree_st_max
endif 
!
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
Write(3,*) "!"
!
! --- Closing "data.inc"
   close(3)
!
! --- Closing "selen.log"
   close(88)
!
STOP
END
!
!
!
!
!
!
 SUBROUTINE MAKE_OFDVMAP (DEGREE, FILE_GMT) 
 IMPLICIT NONE
 INTEGER NRESOLUTION 
 CHARACTER*3  DEGREE
 CHARACTER*20 FILE_GMT, PSFILENAME  
 CHARACTER*80 R_OPTION, B_OPTION, G_OPTION, & 
              W_OPTION, H_OPTION, S_OPTION, & 
	      J_OPTION, XY_OPTION
!
!
! # Revised GS June 2008 for upgrade to "Selen 2.6" 
!
!
  open(9,file=file_gmt,status='unknown')
!
  Write(9,*) "${GMT} gmtset PAPER_MEDIA A4+"
  Write(9,*) "${GMT} gmtset LABEL_FONT_SIZE 26p"
  Write(9,*) "${GMT} gmtset ANOT_FONT_SIZE 16p"
  Write(9,*) " "
!
!
! ======= Settings for various options ======
!
! --- Projection ---
	   J_OPTION =  "-JX16l/16l"
!
! --- Range option ---
           R_OPTION = "-R0.5/300/1e-4/1e-1"
!
! --- Basemap option ---
	   B_OPTION = "-Bf3a2g3:'Harmonic degree':/f3a1g3p:'OF degree variance':WSne"
!
! --- Size of symbol --- 
	   S_OPTION = "-Ss0.45"
!
! --- Color of symbol --- 
	   G_OPTION = "-G0/0/255"
!
! --- Style of interpolating line --- 
	   W_OPTION = "-W8/255/0/0ta" 
!
! --- X-Y offset 
	   XY_OPTION = "-X4 -Y4" 
!
! --- Postscript output filename 
 	   psfilename = "ofdv.ps"
!
!
! ======= Writing the GMT script ======
!

! - Psbasemap 	
  Write(9,*) "${GMT} psbasemap -U/-2/-2/'SELEN 2.9' -K", " ", & 
  	     trim(adjustl(XY_OPTION)), " ", & 
	     trim(adjustl(R_OPTION)),  " ", &    
	     trim(adjustl(B_OPTION)),  " ", &    
	     trim(adjustl(J_OPTION)),  " ", &
	     "> ", trim(adjustl(psfilename))  
!	   
! - A filter 
  Write(9,*) "awk '{print $1, $3}' ofdv.dat > predof.dat"
!
!
! - Plotting the interpolating line  
  Write(9,*) "${GMT} psxy predof.dat -O -K -H3 -B -R -JX ", & 
  	     trim(adjustl(W_OPTION)), " ", & 
	     ">> ", trim(adjustl(psfilename))  
!
! - Plotting symbols 
  Write(9,*) "${GMT} psxy ofdv.dat -O -K -H3 -B -R -JX ", & 
  	     trim(adjustl(S_OPTION)), " ", & 
  	     trim(adjustl(G_OPTION)), " ", & 
	     ">> ", trim(adjustl(psfilename))  
!
! - Plotting a string with harmonic degree	     
  Write(9,*) "echo '20 5e-2 18 0 0 BL MAX DEGREE=", trim(adjustl(DEGREE)), "'", & 
  	     " | ${GMT} pstext -N ", trim(adjustl(J_OPTION)), " ", & 
	                      trim(adjustl(R_OPTION)), " ", & 
			      " -G0 -W255 -O -K ", & 
	     		      ">> ", trim(adjustl(psfilename))  
!
! - Plotting a string with "power" of the interpolating line
  Write(9,*) "${GMT} pstext -O dv-power.tmp -N -JX -R -G0 -W255 >> ", trim(adjustl(psfilename)) 
!
!
! - Cleaning 
  Write(9,*) "/bin/rm predof.dat"
!
 close(9)
!
 End Subroutine Make_ofdvmap
!
!
!
!
!
!
 SUBROUTINE MAKE_RSLC   & 
 	   (TIMEBP,     & 
	    MINRSL,     &
	    MAXRSL,     &
	    RSLINT,     &
	    LONMINC,    &
	    LONMAXC,    &
	    LATMINC,    &
	    LATMAXC,    & 
	    NV,         &
	    CODE,       &
	    RUN,        &
	    NINC,       &
	    NRSLC,      & 
	    TITLICE,    &
	    RESOLUTION, & 
	    ITER,       &
	    MODE,       &
	    DEGREE,     & 
	    VSTRING,    &
	    NAME_REGION,&
	    FILE_GMT,   &
	    SHORT)      
!
!     ----------------------------------------------------------------
! --- A GMT script for plotting RSL CONTOUR LINES - GS January 29 2008
!     ** Inspired by the work of Paolo Stocchi in the Mediterranean ** 
!     ---------------------------------------------------------------- 
!     Revided July 2006 for version 2.6 
!     Revised on April 2010 by GS for version 3.1 ALMA & g95 
!
 IMPLICIT NONE
 REAL*8 LAT_TIME, LON_TIME, LAT_STRING, MID_LON 
 REAL*8 FLAT_MIN, FLAT_MAX, FLON_MIN, FLON_MAX 
 REAL*8 RSL 
 INTEGER I, J, NN, NRSLC 
 CHARACTER*50 NAME_REGION
 CHARACTER*100 TITLE
 CHARACTER*10 TIMEBP, MINRSL, MAXRSL, RSLINT, LONMINC, LONMAXC, LATMINC, LATMAXC  
 CHARACTER*80 T_OPTION, R_OPTION, B_OPTION, G_OPTION, & 
              W_OPTION, C_OPTION, S_OPTION, H_OPTION, & 
	      D_OPTION, J_OPTION, W1_OPTION, W2_OPTION, ANNOT 
 CHARACTER*10 DEGREE, RESOLUTION, TITLICE
 CHARACTER*10 PROJ, COPT, GSPAC 
 CHARACTER*2  LABEL 
 CHARACTER*3  NINC, NV, CODE
 CHARACTER*1  ITER, MODE 
 CHARACTER*20 DATE, TIMC
 CHARACTER*20 FILE_GMT
 CHARACTER*30 VSTRING
 CHARACTER*4  RUN
 CHARACTER*100 SHORT
!	
  open(9,file=file_gmt,status='unknown')
  Write(9,*) "${GMT} gmtset PAPER_MEDIA A4+"
  Write(9,*) "${GMT} gmtset LABEL_FONT_SIZE 24p"
  Write(9,*) "${GMT} gmtset ANOT_FONT_SIZE 10p"
  Write(9,*) "${GMT} gmtset FRAME_WIDTH 0.15c"
  Write(9,*) " "
!
!
!
! ======= Settings for various options ======
!
! --- Range option ---
           R_OPTION = "-R"//trim(adjustl(lonminc))& 
                     //"/"//trim(adjustl(lonmaxc))& 
		     //"/"//trim(adjustl(latminc))&
		     //"/"//trim(adjustl(latmaxc))  
!
! --- Cpt option ---
           T_OPTION = "-T"//trim(adjustl(minrsl))&
                     //"/"//trim(adjustl(maxrsl))&
		     //"/"//trim(adjustl(RSLINT)) 
!
! --- Resolution of coastlines ---
           D_OPTION = "-Di"		     
!
! --- Shades for oceans and lands ---
	   S_OPTION = "-S0/60/255" 
           G_OPTION = "-G150/220/150" 
!
! --- Number of header lines in "rslc.dat" ---
           H_OPTION = "-H9"	   	
! 
! --- C option: a thick sechond level contour
	   C_OPTION="-C10"	
!
! --- Grid spacing for "Surface" ---
	   GSPAC="-I0.1"
!
! --- Projection ---
	   J_OPTION="-JM12"	 
!
! --- Annotation interval ---
           ANNOT='4'
!
! --- B Option ---
           title='"'//trim(adjustl(name_region))//'"'
	   B_OPTION="-Ba"//trim(adjustl(ANNOT))//"/a"&
	                 //trim(adjustl(ANNOT))//"WSEn"&
			 //":."//trim(adjustl(title))//":"
!
! --- Pen attributes for first level of contour lines ---	
           W1_OPTION = "-W3/255/0/0"
!
! --- Pen attributes for sechond level of contour lines ---	
           W2_OPTION = "-W10/255"
!
!
! ======= Writing the GMT script ======
!	   
! --- Surface ---  
    	   Write(9,*) "${GMT} surface ", & 
	               trim(adjustl(GSPAC)), " ", & 
		       trim(adjustl(H_OPTION)), & 
		       " rslc.dat ", & 
		       trim(adjustl(R_OPTION)), & 
		       " -Gg.grd"
!
! --- A no-green palette ---
           Write(9,*) "${GMT} makecpt -Cno_green ", & 
	   	       trim(adjustl(GSPAC)), " ", & 
		       trim(adjustl(T_OPTION)), & 
		       " > pal_rslc.cpt"
!
! --- Psbasemap ---
	   Write(9,*) "${GMT} psbasemap -Y4 -X4 -R  ", & 
	   	       trim(adjustl(J_OPTION)), " ", & 
		       trim(adjustl(B_OPTION)), & 
		       "  -K > rslc-map.ps"
!
! --- Coastlines ---
	   Write(9,*) "${GMT} pscoast ", &
	   	       trim(adjustl(R_OPTION)), " ", &  
		       trim(adjustl(J_OPTION)), " ", & 
		       trim(adjustl(D_OPTION)), " ", & 
       	               " -B ", & 
		       trim(adjustl(S_OPTION)), " ", & 
		       trim(adjustl(G_OPTION)), " ", & 
		       " -O -K >> rslc-map.ps"		 
!
! --- First contour ---
	   Write(9,*) "${GMT} grdcontour -Cpal_rslc.cpt ", & 
	               trim(adjustl(W1_OPTION)), & 
		       " -R  -G4/10 g.grd -JM -B -O -K >> rslc-map.ps"		       
!
! --- Sechond contour ---
	   Write(9,*) "${GMT} grdcontour -U/0.5/0.5/'SELEN 2.9' ", & 
	               trim(adjustl(C_OPTION)), " ", & 
	               trim(adjustl(W2_OPTION)), & 
		       " -R -G4/4 -A1f10 g.grd -JM -B -O -K >> rslc-map.ps"
!
! --- Placing some text strings ---
           call CHAR10_2_FLOAT(lonminc, flon_min)
           call CHAR10_2_FLOAT(lonmaxc, flon_max)
 	   call CHAR10_2_FLOAT(latminc, flat_min)
 	   call CHAR10_2_FLOAT(latmaxc, flat_max)	  
	   mid_lon=flon_min+(flon_max-flon_min)/2. 	  
	   lat_string=flat_max+(flat_max-flat_min)/10.
	   lat_string=48.80 
	   lon_time=flon_max-(flon_max-flon_min)/30.
           lat_time=flat_min+(flat_max-flat_min)/10.  	   	   		       		      
!
	   If(code=='-1') then 
	   CONTINUE
	   Else 
 	   open(4,file='rslc1.tmp',status='unknown') 
 	   Write(4,*) mid_lon, lat_string, "10 0 2 BC -Ice model: ", trim(adjustl(TITLICE)), & 
 	   " -Viscosity: ", trim(adjustl(VSTRING)),&
 	   " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
 	   " -NV=",   trim(adjustl(NV)),      " -CODE=", trim(adjustl(CODE)),& 
 	   " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER)) 	   
 	   close(4) 	   
	   Endif
!
 	   open(4,file='rslc2.tmp',status='unknown') 
	   Write(4,*) lon_time, lat_time, "10 0 0 BR ", trim(adjustl(timebp)), " ka"
	   close(4)	   	   
!
 	   Write(9,*) "${GMT} pstext rslc1.tmp -N -JM -R -G0 -O -K >> rslc-map.ps" 
 	   Write(9,*) "${GMT} pstext rslc2.tmp -N -JM -R -G0 -W255 -O >> rslc-map.ps" 

!
 End Subroutine MAKE_RSLC
!
!
!
!
!
!
	SUBROUTINE SCAN_REGION (FILEIN, FILEOUT, N, TIME, TIMEC, & 
		                LONMINC, LONMAXC, LATMINC, LATMAXC, & 
				MIN_RSLC, MAX_RSLC, CINT, & 
				NAME_OF_REGION)
	IMPLICIT NONE
!	
! # Scans "FILEIN" to retrieve information upon the number of points and 
!   time for the rsl contour analysis     === GS January 28, 2007 ===
!
	CHARACTER*10 LONMINC, LONMAXC, LATMINC, LATMAXC, TIMEC, MIN_RSLC, MAX_RSLC, CINT 
	CHARACTER*30 FILEIN, FILEOUT
	CHARACTER*50 NAME_OF_REGION
        CHARACTER*100 SS(10) 
        CHARACTER*200 LINE
	INTEGER, PARAMETER :: LARGE_INTEGER = 100
	INTEGER, PARAMETER :: MAXP = 19999
        INTEGER J, N, IERR, NOUT
	INTEGER ILAT, ILON, I_LON_MAX, I_LAT_MAX
	REAL*8 TIME, LON_LAT_INC, LON, LAT, LON_MIN, LON_MAX, LAT_MIN, LAT_MAX
	REAL*8, PARAMETER :: FRAME_WIDTH = 2. 
!
!
!
! ==> 1) Reading bounds from the input file
!
	IERR=0
!
	open(70,file=filein,status='old')
!
		do j=1, large_integer 
!
	        read(70,'(a200)',end=1) line 
!
		if(line(1:2)=='10') then 
				      CALL scan_string(line, 1, ss, nout)
				      name_of_region=trim(adjustl(ss(1)))
				 endif	
		if(line(1:2)=='20') then 
				      CALL scan_string(line, 1, ss, nout)
				      timec=trim(adjustl(ss(1)))
				      CALL CHAR100_2_REAL(ss(1), time)
				 endif
		if(line(1:2)=='30') then 
				      CALL scan_string(line, 2, ss, nout)		
			              CALL CHAR100_2_REAL(ss(1), lon_min)
			              CALL CHAR100_2_REAL(ss(2), lon_max)	
				      If(lon_max<lon_min) IERR=1  	
				      LONMINC=trim(adjustl(ss(1))) 
				      LONMAXC=trim(adjustl(ss(2))) 	
            			 endif				    
		if(line(1:2)=='40') then 
				      CALL scan_string(line, 2, ss, nout)
				      CALL CHAR100_2_REAL(ss(1), lat_min)
			              CALL CHAR100_2_REAL(ss(2), lat_max)
				      If(lat_max<lat_min) IERR=1  
				      LATMINC=trim(adjustl(ss(1))) 
				      LATMAXC=trim(adjustl(ss(2))) 				 			      
				 endif
		if(line(1:2)=='50') then 
				      CALL scan_string(line, 1, ss, nout)
				      CALL CHAR100_2_REAL(ss(1), lon_lat_inc)
				      If(lon_lat_inc<=0.)IERR=1			
				 endif 
		if(line(1:2)=='60') then 
				      CALL scan_string(line, 3, ss, nout)
				      MIN_RSLC=trim(adjustl(ss(1)))
				      MAX_RSLC=trim(adjustl(ss(2)))
	 			      CINT=trim(adjustl(ss(3)))
				 endif 								 				 		  
		enddo
1       continue
        close(70)
!
	If(IERR==1) then 
        	Write(*, *) "File ", trim(adjustl(filein)), " is badly configured"
        	Write(88,*) "File ", trim(adjustl(filein)), " is badly configured"
        	Call Stop_Config 
	Endif	
	
!
! ==> 2) Computing lon-lat of points and writing on the outpout file... 
!
! --- Counts the number of points in the User-supplied file "filein" 

	
! NEW code as of November 2009 
	i_lon_max=aint((lon_max-lon_min)/lon_lat_inc) + 1
	i_lat_max=aint((lat_max-lat_min)/lon_lat_inc) + 1
	n=0 
	do 2 ilon=1, i_lon_max 
		do 2 ilat=1, i_lat_max 
		n=n+1	
2	continue 
        Write(88,*) "The number of points in file ", trim(adjustl(filein)), " is", n
!	
!	n=0
!	do 2 lon = lon_min, lon_max, lon_lat_inc 
!		do 2 lat = lat_min, lat_max, lon_lat_inc 
!		n=n+1	
!2       continue
!        Write(88,*) "The number of points in file ", trim(adjustl(filein)), " is", n
!	
! --- Adds a frame to avoid border effects using GMT  
	lon_min=lon_min - frame_width 
	lon_max=lon_max + frame_width 
	lat_min=lat_min - frame_width
	lat_max=lat_max + frame_width
!
! NEW code as of November 2009 
	i_lon_max=aint((lon_max-lon_min)/lon_lat_inc) + 1
	i_lat_max=aint((lat_max-lat_min)/lon_lat_inc) + 1
	n=0 
	do 3 ilon=1, i_lon_max 
		do 3 ilat=1, i_lat_max 
		n=n+1	
3	continue 
        Write(88,*) "After framing (2 degrees NSEW), the number of points is ", n	
!
!	n=0
!	do 3 lon = lon_min, lon_max, lon_lat_inc 
!		do 3 lat = lat_min, lat_max, lon_lat_inc 
!		n=n+1	
!3       continue
!        Write(88,*) "After framing (2 degrees NSEW), the number of points is ", n	
!
	If(n>=maxp) then 
        	Write(*, *) "The total number of points in file ", trim(adjustl(filein)), " is", n
		Write(*, *) "this exceeds the maximum allowed (i. e.,", maxp, ")" 
        	Write(88,*) "The total number of points in file ", trim(adjustl(filein)), " is", n 
		Write(88,*) "this exceeds the maximum allowed (i. e.,", maxp, ")" 
        	Call Stop_Config 
	Endif	
!
! --- Composes the final set of points and reports on "fileout"
!     Negative lon values are transformed in the range [0:360] 
!
        open(71,file=fileout,status='unknown') 	
	i_lon_max=aint((lon_max-lon_min)/lon_lat_inc) + 1
	i_lat_max=aint((lat_max-lat_min)/lon_lat_inc) + 1
!
	do 4 ilon=1, i_lon_max 
		lon=lon_min+(ilon-1)*lon_lat_inc		
!
		do 4 ilat=1, i_lat_max	
		lat=lat_min+(ilat-1)*lon_lat_inc							
!
			if(lon>=0.) then 
				write(71,*) lon, lat 
				else 
				write(71,*) lon+360., lat 		
			Endif		
4       continue
	close(71) 
!
!        open(71,file=fileout,status='unknown') 	
!	 do 4 lon = lon_min, lon_max, lon_lat_inc 
!		do 4 lat = lat_min, lat_max, lon_lat_inc 				
!			if(lon>=0.) then 
!				write(71,*) lon, lat 
!				else 
!				write(71,*) lon+360., lat 		
!			Endif		
!4       continue
!	 close(71) 
! 	
	End subroutine scan_region
!
!
!
!
!
!
   	Subroutine make_eslplot (ninc, titlice, shof_file, file_gmt)
        IMPLICIT NONE
!
!--- Creates a GMT script for plotting the Equivalent Sea Level as a function
!    of time, for the ice model with name "titlice".  == GS Urbino 8/12/07 ==
!
!    Revised July 4 2008 for v. 2.6 
!    Also revided July 8
!
!
        CHARACTER*100 B_OPTION
	CHARACTER*40 Y_LABEL, R_OPTION, ESL_STRING
	CHARACTER*20 FILE_GMT, SHOF_FILE
	CHARACTER*40 AWKSTRING
	CHARACTER*10 TITLICE
	CHARACTER*3  NINC
	REAL*8 VMAX, VMIN, ESL
	INTEGER K
!
!
        open(9,file=file_gmt,status='unknown')
        Write(9,*) "${GMT} gmtset PAPER_MEDIA A4+"
        Write(9,*) "${GMT} gmtset LABEL_FONT_SIZE 24p"
        Write(9,*) "${GMT} gmtset ANOT_FONT_SIZE 12p"
        Write(9,*) " "
!
!
! --- R- option 
        if(titlice=='ICE1')     R_OPTION="-R-2/23/0/220"
        if(titlice=='ICE3G')    R_OPTION="-R-2/23/0/220"
        if(titlice=='DISK')     R_OPTION="-R-2/23/0/220"
        if(titlice=='ICE5G')    R_OPTION="-R-2/23/0/220" 
        if(titlice=='ICE5G26')  R_OPTION="-R-2/28/0/220" 
        if(titlice=='IJ05MOD')  R_OPTION="-R-2/23/0/220" 
        if(titlice=='ANU05')    R_OPTION="-R-2/32/0/220" 
        if(titlice=='ICAP')     R_OPTION="-R-2/23/0/220"
        if(titlice=='ALPS')     R_OPTION="-R-2/23/0/2" 
        if(titlice=='ALASKA')   R_OPTION="-R-0.1/2/0/0.05" 
!
!
! --- y label 
	if(titlice=='ICE1')     y_label='ICE1 ESL (m)'
	if(titlice=='ICE3G')    y_label='ICE3G ESL (m)'
	if(titlice=='DISK')     y_label='DISK ESL (m)'
	if(titlice=='IJ05MOD')  y_label='IJ05MOD ESL (m)'
	if(titlice=='ICE5G')    y_label='ICE5G ESL (m)'
	if(titlice=='ICE5G26')  y_label='ICE5G ESL (m)'
        if(titlice=='ANU05')    y_label='ANU ESL (m)'
	if(titlice=='ALPS')     y_label='ALPS ESL (m)'
        if(titlice=='ALASKA')   y_label='ALASKA ESL (m)'
        if(titlice=='ICAP')     y_label='ICAP ESL (m)'

!
! --- B option 
       if( titlice=='ALASKA' ) then
          b_option=" -Ba0.5f0.1:'time (ka)':/a0.02f0.01WSen:"//"'"//trim(adjustl(y_label))//"':"
       else
          b_option=" -Ba2f1:'time (ka)':/a50f10WSen:"//"'"//trim(adjustl(y_label))//"':"
       endif
!	
! --- A psbasemap 
        Write(9,*) "${GMT} psbasemap -U'SELEN 2.9' -X6 -Y10 ", " ", & 
	            trim(adjustl(b_option)), " ", & 
		    trim(adjustl(R_OPTION)), " -JX14/9 -K > plot.ps"
 
! --- Draw thick and thin lines 
	if(ninc=='30') awkstring="'{print 30.-$1, $2}'"
	if(ninc=='26') awkstring="'{print 26.-$1, $2}'"
	if(ninc=='21') awkstring="'{print 21.-$1, $2}'"
	if(ninc=='18') awkstring="'{print 18.-$1, $2}'"
        if(ninc=='40') awkstring="'{print (40.-$1)*0.050125, $2}'"
!	
	Write(9,*) "awk "//trim(adjustl(awkstring))//" esl-thin.dat > esl.tmp"	
	Write(9,*) "${GMT} psxy esl.tmp -M -B -R -JX -W2/0 -O -K >> plot.ps"
!
	Write(9,*) "awk "//trim(adjustl(awkstring))//" esl.dat > esl.tmp"	
	Write(9,*) "${GMT} psxy esl.tmp -M -B -R -JX -W6/0 -O -K >> plot.ps"
!
! --- Inset with total Equivalent Sea level 
	 Write(9,*) "${GMT} pstext -N esl-tot.dat -JX -R -O >> plot.ps"
!
! --- The end 
         Write(9,*) "mv plot.ps esl.ps"
!
         close(9) 
!
   	 End Subroutine make_eslplot 

!
!
!
!
!
!
	SUBROUTINE OF_PIXELIZATION_REAL
        IMPLICIT NONE
	CHARACTER*20, PARAMETER :: FILENAME='px.gmt'
!				        
!--- REALISTIC OCEAN FUNCTION (according to GMT topography database)   
!
!--- Creates a GMT script named "px.gmt", which is used to separate wet from dry pixels. 
!    The wet/dry distribution of pixels *is* sensitive to the GMT options -D (coastlines 
!    resolution) and -A (spatial filter), so different choices may  give (slightly) 
!    different distributions, affecting all the computations in SELEN. 
!				        
!    [*** GS & FC Bologna 8/12/07 ***]
!    Last change: GS February 21, 2008
!    Name changed on July 19 for v. 2.6 
!
        open(13,file=filename,status='unknown')
!Write(13,*)"echo '     - ", trim(adjustl(filename))//":", " gmt-selecting pixels in <<pxa.dat>>'"
	Write(13,*)"#"
	Write(13,*)"# GMT script for sepating wet from dry pixels: px.gmt"
	Write(13,*)"#               REALISTIC OCEAN FUNCTION             " 
	Write(13,*)"# ====== Made by GS and FC on December 8, 2007 ======"	
	Write(13,*)"#" 
	Write(13,*)"${GMT} gmtselect pxa.dat -H4 -Df -R0/360/-90/90  -A0 -JQ180/200 -Nk/s/k/s/k > weta.dat"
	Write(13,*)"${GMT} gmtselect pxa.dat -H4 -Df -R              -A0 -JQ        -Ns/k/s/k/s > drya.dat"
	close(13) 		
	End Subroutine OF_PIXELIZATION_REAL
!
!
!
!   	
!
!
        SUBROUTINE MAKE_RAINBOW_PALETTE 
        IMPLICIT NONE
!
!--- Creates a "rainbow palette", named "ice-pal.cpt" suitable
!    for plotting both the original and the rechonstructed ice 
!    distribution...     *** GS and FC Bologna 7/12/07 ***
!
        open(13,file='ice-pal.cpt',status='unknown') 
	Write(13, '(a1)') "#"
	Write(13,'(a74)') "# A rainbow palette by Florence and Giorgio as of 7 Dec. 2007: ice-pal.cpt"
	Write(13, '(a1)') "#"	
	Write(13,'(a58)') "-500   255     0       255     0       255     0       255"
	Write(13,'(a58)') "0      164     0       255     500     164     0       255"	
	Write(13,'(a58)') "500    73      0       255     1000    73      0       255"
	Write(13,'(a58)') "1000   0      18       255     1500    0      18       255"
	Write(13,'(a58)') "1500   0     109       255     2000    0     109       255"
	Write(13,'(a58)') "2000   0     201       255     2500    0     201       255"
	Write(13,'(a58)') "2500   0     255       219     3000    0     255       219"
	Write(13,'(a58)') "3000   0     255       127     3500    0     255       127"
	Write(13,'(a57)') "3500   0     255        36     4000    0     255       36"
	Write(13,'(a56)') "4000   54    255         0     4500    54    255       0"
	Write(13,'(a56)') "4500   146   255         0     5000    146   255       0"
	Write(13,'(a56)') "5000   237   255         0     5500    237   255       0"
	Write(13,'(a56)') "5500   255   182         0     6000    255   182       0"
	Write(13,'(a56)') "6000   255   91          0     6500    255   91        0"
	Write(13,'(a17)') "B  255    0   255"
	Write(13,'(a17)') "F  255    0     0"
	Write(13,'(a17)') "N  128  128   128"	    
        close(13) 	
!
 	END SUBROUTINE MAKE_RAINBOW_PALETTE 
!
!
!
!
!
!
!
 SUBROUTINE MAKE_TGAUGESSCA (DATAFILE, FILE1_GMT)
 IMPLICIT NONE
!
! --- This routine perfoms simple statistical analyses on tide-gauge data 
!     and creates a two-frames scatterplot of tide gauges trends vs years 
!     of data 
!
!     *** Last modified GS December 8, 2007 *** 
!
!     === Re touched on July 2008 for SLEN 2.6 ===
! 
!
 CHARACTER*30 DATAFILE
 CHARACTER*20 FILE1_GMT, FILE2_GMT 
 INTEGER I, J, K, N, P, NH, IMIN, IMAX, PMIN, PMAX, IEARS_MAX, IEARS_MIN
 INTEGER, PARAMETER :: ILARGE=10000 ! A large number 
 INTEGER IEARS(ILARGE) 
 CHARACTER*100 RIGA		   ! A row 
 CHARACTER*8 PCODE                 ! Tide gauge code  
 CHARACTER*3 GLOSS_CODE            ! GLOSS code 
 CHARACTER*3 YEARS                 ! Number of years 
 CHARACTER*11 RANGE_YEARS	   ! Range of years 
 CHARACTER*9 DATUM		   ! Trend 
 CHARACTER*3 PLUS_MINUS 	   ! +/- 
 CHARACTER*5 ERROR		   ! Error 
 CHARACTER*7 STDV		   ! Std. deviation of residues 						       
 CHARACTER*2 LAT1, LAT2 	   ! Latitide degrees and minutes      
 CHARACTER*3 LON1		   ! Longitude degrees 
 CHARACTER*2 LON2		   ! Longitude minutes 
 CHARACTER*30 NAME(ILARGE)	   ! Name of the station 
 CHARACTER*1 WLAT, WLON 	   ! Where is the station (e, w, s, n) 
 CHARACTER*10 YMINC, YMAXC, DELTAC, IMINC, IMAXC, PMINC, PMAXC 
 REAL*8 LAT(ILARGE), NLAT, DLAT 	   ! Latitude  
 REAL*8 LON(ILARGE), NLON, DLON 	   ! Longitude 
 REAL*8 TREND(ILARGE), D_TREND(ILARGE)
 REAL*8 AVE(ILARGE), DAVE(ILARGE), SUMP(ILARGE)	
 REAL*8 AVEW, DAVEW, SUMPW	
 REAL*8 RANGE, DELTA
 REAL*8 AVE_RANGE, AVE_MAX, AVE_MIN, AVE_DELTA
 REAL*8 TREND_MAX, TREND_MIN 
 CHARACTER*100  B_OPTION, R_OPTION  
 CHARACTER*80 STRING1, STRING2
	  
!
!
! ===============================================
! ====== Analysis of tide gauges database =======
! ===============================================
!
! --- Open and reads the tide-gauges data file for getting info...
!
! --- Counting the header lines (beginning with '#'), and the data lines
       open(10,file=DATAFILE,status='unknown')
       nh=0 ; n=0
       do i=1, ilarge 
       		read(10,'(a100)',end=1) riga
       		if(riga(1:1)=='#') nh=nh+1
       		if(riga(1:1)/='#') n=n+1
       enddo 
 1     close(10)  
!
! --- Opening the tide-gauges data file... 
       open(10,file=DATAFILE,status='unknown')
!
! --- Reading again the header lines  
       do i=1, nh ; read(10,'(a100)',end=1) riga ; enddo   
!
! --- Loop on the tide-gauge stations 
       do 5 i=1, n
! 
! --- Reading one line   
       	     read  (10,200) pcode, gloss_code, years, range_years, & 
   		     	    datum , plus_minus, error, stdv, & 
			    lat1, lat2, wlat, lon1, lon2, wlon, name(i)  
! --- Reading format 
200    format(a8,1x,a3,1x,a3,2x,a11,1x,a8,1x, & 
              a3,1x,a5,1x,a7,3x,a2,1x,a2,1x,a1,1x,& 
              a3,1x,a2,1x,a1,3x,a30)			      
!
! --- Converts the range in floating-point format 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; write(50,'(a3)') years ; close(50) 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; read (50,*)      iears(i) ; close(50)
!
! --- Converts the trend in floating-point format 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; write(50,'(a9)') datum ; close(50) 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; read (50,*)      trend(i) ; close(50)	     
!
! --- Converts the error on trend in floating-point format 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; write(50,'(a5)') error ; close(50) 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; read (50,*)    d_trend(i) ; close(50)	     	     	  
!
! --- Extracts lon and lat of the station from the texst strings... 
             call  find_lon_lat (lat1, lat2, wlat, lon1, lon2, wlon, lon(i), lat(i)) 
!
! --- End of Loop
 5 continue 	
!
      close(10)		
!
!
! ==========================================
! --- Filter of data for stats and plots ---
! ==========================================
!
! --- Target file for stats
!
        open(73,file='tgauges-stat.dat',status='unknown')
	Write(73,*) " "    
!
! >>> #1 Creating a scatterplot file  
!
  	Open(28,file='tgauges-scpl.dat',status='unknown') 
!
	do i=1, n 
		write(28,*) iears(i), trend(i), d_trend(i), lon(i), lat(i), name(i)  
	enddo 
	Close(28) 
!
! >>> #2 Weighted average of the whole set of data  
!
	avew  = 0. 
	davew = 0. 
	sumpw = 0.
	do i=1, n      
		avew  = avew  + trend(i)   *float(iears(i))
		davew = davew + d_trend(i) *float(iears(i)) 
		sumpw = sumpw +             float(iears(i))
	enddo
	avew  = avew/sumpw 
	davew = davew/sumpw
	Write(73,*) "Average weighted rate of sea level change: "
	Write(73,'(f5.2,a1,a3,a1,f5.2,a6)') & 
		    avew, " ", plus_minus, " ", davew, ' mm/yr'
	Write(73,*) " "
!
! >>> #3 Average over subsets of increasing number of years 
!
	iears_max=-99999
	iears_min=+99999
	do i=1, n 
		if(iears(i)>=iears_max)iears_max=iears(i) 
		if(iears(i)<=iears_min)iears_min=iears(i) 		
	enddo
	Write(73,*) 'Max number of years of data = ', iears_max 
	Write(73,*) 'Min number of years of data = ', iears_min 
	Write(73,*) " "	
!
!
! -------------------------------
! --- GMT scripts for Figures ---
! -------------------------------
!	
! === Figure 1: A two-frame figure with scattered tide-gauges data
!		with and w/o error bars, with the weighted average 
!		based on all data. Useful for presentations. 
!
! --- Determines min and max trend of tide-gauge data
	trend_max=-99999
	trend_min=+99999
	do i=1, n 
		if(trend(i)>=trend_max)then
					trend_max=trend(i) 
					imax=i 
		Endif
		if(trend(i)<=trend_min)then 
					trend_min=trend(i) 
					imin=i
	        Endif			
	enddo
	Write(73,*) 'Max trend (mm/yr) = ', trend_max 
	Write(73,*) 'Min trend (mm/yr) = ', trend_min
	Write(73,*) " "	
	close(73) 
!
	trend_max=trend_max+d_trend(imax)
	trend_min=trend_min-d_trend(imin)	
!
	range = trend_max - trend_min 
!	
	if(range>=0  .and.range<50)  delta=10  
	if(range>=50 .and.range<100) delta=10	
	if(range>=100.and.range<200) delta=20
	if(range>200) 		  delta=50 
!
   	OPEN  (10,FILE='junk.dat',STATUS='unknown')
	WRITE (10,'(i10)') int(trend_min) ; WRITE (10,'(i10)') int(trend_max)
	WRITE (10,'(i10)') int(delta) 
	WRITE (10,'(i10)') int(iears_min) ; WRITE (10,'(i10)') int(iears_max)	
	CLOSE(10)

   	OPEN  (10,FILE='junk.dat',STATUS='unknown')
	READ (10,'(a10)') yminc ; READ (10,'(a10)') ymaxc
	READ (10,'(a10)') deltac
	READ (10,'(a10)') iminc ; READ (10,'(a10)') imaxc 	
	CLOSE(10)
	iminc='0'
!
! --- Target fiel for GMT script 
 	 open (19,file=file1_gmt,status='unknown')
	 Write(19,*) "${GMT} gmtset PAPER_MEDIA A4+" 
	 Write(19,*) "${GMT} gmtset HEADER_FONT_SIZE 24p"
	 Write(19,*) "${GMT} gmtset FRAME_WIDTH 0.1c"
	 Write(19,*) " "
!
! --- -R and -B options 
	 r_option = "-R"//trim(adjustl(iminc))//"/"//trim(adjustl(imaxc))//"/"//trim(adjustl(yminc))//"/"//trim(adjustl(ymaxc))
	 b_option = "-Ba30f10:'Lenght of rechord (yr)':/a"//trim(adjustl(deltac))//"f"//trim(adjustl(deltac))//"WSen:'Rate of SLC, mm/yr':"
!
! --- Basemap for the top frame 
	 Write(19,*) "${GMT} psbasemap -X4 -Y17 ", " "//trim(adjustl(b_option)), & 	 		                    
	 	 			            " "//trim(adjustl(r_option)), " -P -JX14/8  -K >  plot.ps" 
!
! --- A title for top frame
	 open (8,file='tmpctitle1',status='unknown')
	 Write(8,*) iears_max/2, trend_max+trend_max/6, " 16 0 1 BC Scatterplot from file: ", trim(adjustl(datafile)) ; close(8)  
	 Write(19,*) "${GMT} pstext -N tmpctitle1", " -JX -R -O -K >> plot.ps"
!
! --- Filters the useful data	 
         Write(19,*) "awk '{print $1, $2, $3}' tgauges-scpl.dat > tgauges-scpl.tmp"
!	 
! --- Top frame: data with errorbars  
	 Write(19,*) "${GMT} psxy tgauges-scpl.tmp ", "-Ey0.25/2 -B -R -JX -O -K >> plot.ps" 
!
! --- Basemap for the bottom frame 
	 Write(19,*) "${GMT} psbasemap -U'SELEN 2.9' -X0 -Y-12 ", " "//trim(adjustl(b_option)), & 
	 		                        " "//trim(adjustl(r_option)), " -P -JX14/8  -O -K >>  plot.ps" 	
! --- A title for bottom frame
	 open (8,file='tmpctitle2',status='unknown')
	 Write(8,*) iears_max/2, trend_max+trend_max/6, " 16 0 1 BC Scatterplot from file: ", trim(adjustl(datafile)) ; close(8)  
	 Write(19,*) "${GMT} pstext -N tmpctitle2", " -JX -R -O -K >> plot.ps"
!	 
! --- Bottom frame: data w/o errorbars, only dots   
	 Write(19,*) "${GMT} psxy tgauges-scpl.tmp ", "-Sc0.1 -G0 -B -R -JX -O -K >> plot.ps" 
!	 
! --- A note about the average rate, in the bottom frame
	open (8,file='tmpctitle3',status='unknown')
	Write(8,*) "Average (weighted) rate:"
	Write(8,'(f5.2,a1,a3,a1,f5.2,a6)') avew, " ", plus_minus, " ", davew, ' mm/yr'
	close(8) 
	open (8,file='tmpctitle3',status='unknown')
	read(8,'(a80)') string1 
	read(8,'(a80)') string2		
	close(8) 
	open (8,file='tmpctitle4',status='unknown')
	Write(8,*) iears_max-iears_max/4, trend_max-1*delta/1,   " 12 0 1 BC ", trim(adjustl(string1)) ; close(8)
	Write(19,*) "${GMT} pstext -N tmpctitle4", " -JX -R -O -K >> plot.ps" 
	open (8,file='tmpctitle5',status='unknown')
	Write(8,*) iears_max-iears_max/4, trend_max-3*delta/2, " 14 0 2 BC ", trim(adjustl(string2)) ; close(8)	
	Write(19,*) "${GMT} pstext -N tmpctitle5", " -G0 -JX -R -O >> plot.ps" 	
	Write(19,*) "mv plot.ps tgauges-scpl.ps"
!
! --- Closing 
	close(19) 
!
 END SUBROUTINE MAKE_TGAUGESSCA
!
!
!
!
!
!
 SUBROUTINE MAKE_TGAUGES (DATAFILE, FILE_GMT)
 IMPLICIT NONE
 CHARACTER*30 DATAFILE
 CHARACTER*20 FILE_GMT
 INTEGER I, J, K, N, NH, NALL, IEARS, NGE30, NGE60, NGE90 
!
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +   
! --- This routine creates a GMT script for plotting the spatial distribution 
!     of the tide gauges in database DATAFILE. Last modified GS Nov. 07 2007 
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +   
! Retouched on April 2010 by GS
!
 INTEGER, PARAMETER :: MAXN=100000 ! A large number  
 CHARACTER*100 RIGA		   ! A row 
 CHARACTER*8 PCODE                 ! Tide gauge code  
 CHARACTER*3 GLOSS_CODE            ! GLOSS code 
 CHARACTER*3 YEARS                 ! Number of years 
 CHARACTER*11 RANGE_YEARS	   ! Range of years 
 CHARACTER*9 DATUM		   ! Trend 
 CHARACTER*3 PLUS_MINUS 	   ! +/- 
 CHARACTER*5 ERROR		   ! Error 
 CHARACTER*7 STDV		   ! Std. deviation of residues 						       
 CHARACTER*2 LAT1, LAT2 	   ! Latitide degrees and minutes      
 CHARACTER*3 LON1		   ! Longitude degrees 
 CHARACTER*2 LON2		   ! Longitude minutes 
 CHARACTER*30 NAME		   ! Name of the station 
 CHARACTER*1 WLAT, WLON 	   ! Where is the station (e, w, s, n)  
 REAL*8 LAT, NLAT, DLAT 	   ! Latitude  
 REAL*8 LON, NLON, DLON 	   ! Longitude  		  
!
!
!
! ====== Part 1: Analysis of tide gauges database...
!
! --- Target files for lon-lat of tide gauge sites: 
!     	Unit 39: all stations, 
!     	Unit 40: stations with >= 30 years of data, 
!     	Unit 41: stations with >= 60 years of data, 
!     	Unit 42: stations with >= 90 years of data...
! 
       open(39,file='lon-lat-tgauges-all.dat', status='unknown')
       open(40,file='lon-lat-tgauges-ge30.dat',status='unknown')
       open(41,file='lon-lat-tgauges-ge60.dat',status='unknown')
       open(42,file='lon-lat-tgauges-ge90.dat',status='unknown')
!
! --- An header 
	Write(39,*) " lon-lat (deg) & name of tide gauges sites - all data from <<rlr-trends.txt>> "
	Write(40,*) " lon-lat (deg) & name of tide gauges sites - only stations with >= 30 yrs of data"
	Write(41,*) " lon-lat (deg) & name of tide gauges sites - only stations with >= 60 yrs of data"
	Write(42,*) " lon-lat (deg) & name of tide gauges sites - only stations with >= 90 yrs of data"	
!
! --- Open and reads the tide-gauges data file for getting info...
!
! --- Counting the header lines (beginning with '#'), and the data lines
       open(10,file=DATAFILE,status='unknown')
       nh=0 ; n=0
       do i=1, maxn 
       		read(10,'(a100)',end=1) riga
       		if(riga(1:1)=='#') nh=nh+1
       		if(riga(1:1)/='#') n=n+1
       enddo 
 1     close(10)  
!
! --- Counters
	nall=n ; nge30=0 ; nge60=0 ; nge90=0
!
! --- Opening the tide gauge data file... 
       open(10,file=DATAFILE,status='unknown')
!
! --- Reading again the header lines  
       do i=1, nh ; read(10,'(a100)',end=1) riga ; enddo  
!
! --- Loop on the tide gauges stations 
       do 5 i=1, n
! 
! --- Reading one line   
       	     read  (10,200) pcode, gloss_code, years, range_years, & 
   		     	    datum , plus_minus, error, stdv, & 
			    lat1, lat2, wlat, lon1, lon2, wlon, name      
!
! --- Converts the range in floating-point format 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; write(50,'(a3)') years ; close(50) 
	     open(50,file='junk-tgauges.tmp',status='unknown') ; read (50,*)      iears ; close(50)	     
!
! --- Extracts lon and lat of the station from the texst strings... 
             call  find_lon_lat (lat1, lat2, wlat, lon1, lon2, wlon, lon, lat) 
!
! --- Reports lon-lat on lon-lat files... 
             write(39,*) lon, lat, name 
	     if(iears>=30)  then 
	     			nge30 = nge30 + 1 
				write(40,*) lon, lat, name 
				endif			    
	     if(iears>=60)  then
	     			nge60 = nge60 + 1 
				write(41,*) lon, lat, name 
				endif 		    
	     if(iears>=90)  then 					
	     			nge90 = nge90 + 1 
				write(42,*) lon, lat, name 
	     endif				     			
!
! --- End of Loop
 5 continue 	
!
      close(10) ; do i=39, 42 ; close(i) ; enddo
!
! --- Reading format 
200    format(a8,1x,a3,1x,a3,2x,a11,1x,a8,1x, & 
              a3,1x,a5,1x,a7,3x,a2,1x,a2,1x,a1,1x,& 
              a3,1x,a2,1x,a1,3x,a30)
!
!
! ====== Part 2: Creates a GMT script for plotting a map of tide gauges sites  
!
! ------ Target GMT file for a map of RSL sites 
!
 	open(29,file=file_gmt,status='unknown')
!
	Write(29,*) "${GMT} gmtset PAPER_MEDIA A4+" 
	Write(29,*) "${GMT} gmtset HEADER_FONT_SIZE 24p"
	Write(29,*) "${GMT} gmtset FRAME_WIDTH 0.1c"
	Write(29,*) "${GMT} gmtset ANOT_FONT_SIZE 12p"
!
! ------ Four frames (Mercator projection)
!
! --- Frame with all data 
	Write(29,*) "${GMT} psbasemap -X5 -Y11 -Ba180/a80f80Wsen -R0/360/-80/80  -JM9 -K > map-tgauges.ps"
 	Write(29,*) "${GMT} pscoast -G0/120/0 -S0/0/220 -B -R -O -K  -JM -Dc  -A10000 >> map-tgauges.ps"
 	Write(29,*) "${GMT} psxy -H1 lon-lat-tgauges-all.dat -B -R -JM -Sc0.08 -G220 -O -K >> map-tgauges.ps"
	open(18,file='titleall.tmp',status='unknown') 
	Write(18,*) " 180 84.4 16 0 1 BC all tide gauge stations"
	Write(18,*) " 180 82.2 12 0 2 BC N=", nall ; close(18) 	
	Write(29,*) "${GMT} pstext -N titleall.tmp -G0 -JM -R -O -K>> map-tgauges.ps"
	Write(29,*) ""
!
! --- Title for all frames 
	open(18,file='titletgauges.tmp',status='unknown') 
	Write(18,*) "0 87.4 24 0 3 BL Distribution of tide gauge stations as of 1/22/07" ; close(18) 		
	Write(29,*) "${GMT} pstext -N titletgauges.tmp -G0 -JM -R -O -K >> map-tgauges.ps"
	Write(29,*) ""
!
! --- Frame with data with age >= 30 yrs  
	Write(29,*) "${GMT} psbasemap -X10 -Y0 -Ba180/a80f80wsEn -R -JM -O -K >> map-tgauges.ps"
 	Write(29,*) "${GMT} pscoast -G120 -S0/0/220 -B -R -O -K  -JM -Dc  -A10000 >> map-tgauges.ps"
 	Write(29,*) "${GMT} psxy -H1 lon-lat-tgauges-ge30.dat -B -R -JM -Sc0.08 -G220 -O -K >> map-tgauges.ps"
	open(18,file='titlege30.tmp',status='unknown') 
 	Write(18,*) "180 84.4 16 0 1 BC more than 30 years of data"
 	Write(18,*) "180 82.2 12 0 2 BC N=", nge30  ; close(18) 
 	Write(29,*) "${GMT} pstext -N titlege30.tmp -G0 -JM -R -O -K >> map-tgauges.ps"
	Write(29,*) ""
!
! --- Frame with data with age >= 60 yrs  
 	Write(29,*) "${GMT} psbasemap -X-10 -Y-9 -Ba180/a80f80WSen -R -JM -O -K -U"//"'SELEN 2.9'", " >> map-tgauges.ps"
 	Write(29,*) "${GMT} pscoast -G120 -S0/0/220 -B -R -O -K  -JM -Dc  -A10000 >> map-tgauges.ps"
 	Write(29,*) "${GMT} psxy -H1 lon-lat-tgauges-ge60.dat -B -R -JM -Sc0.08 -G220 -O -K >> map-tgauges.ps"
	open(18,file='titlege60.tmp',status='unknown') 
 	Write(18,*) "180 84.4 16 0 1 BC more than 60 years of data"
 	Write(18,*) "180 82.2 12 0 2 BC N=", nge60 ; close(18)  
 	Write(29,*) "${GMT} pstext -N titlege60.tmp -G0 -JM -R  -O -K >> map-tgauges.ps"
	Write(29,*) ""
!
! --- Frame with data with age >= 90 yrs 
	Write(29,*) "${GMT} psbasemap -X10 -Y0 -Ba180/a80f80wSEn -R -JM -O -K >> map-tgauges.ps"
 	Write(29,*) "${GMT} pscoast -G120 -S0/0/220 -B -R -O -K  -JM -Dc  -A10000 >> map-tgauges.ps"
 	Write(29,*) "${GMT} psxy -H1 lon-lat-tgauges-ge90.dat -B -R -JM -Sc0.08 -G220 -O -K >> map-tgauges.ps"
	open(18,file='titlege90.tmp',status='unknown')  	
 	Write(18,*) "180 84.4 16 0 1 BC more than 90 years of data"
	Write(18,*) "180 82.2 12 0 2 BC N=", nge90 ; close(18)  !
	Write(29,*) "${GMT} pstext -N titlege90.tmp -G0 -JM -R -O >> map-tgauges.ps"
!
	close(29) 	
!	
 END SUBROUTINE MAKE_TGAUGES
!
!
!
!
!
!
 SUBROUTINE MAKE_RMAPS (TITLICE, &
 			RESOLUTION, & 
			NV, & 
			CODE, & 
			ITER, & 
			MODE, & 
			DEGREE, & 
			VSTRING, & 
			OPT, & 
			SHORT_VISCO)
!			
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
! Prepares a GMT script for *** REGIONAL maps *** of dot-S, U & N at present time 
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
! --- Last change GS March 30 08 --- 
! Also revised April 2010 for ALMA coupling 
!
!
 IMPLICIT NONE
 INTEGER, PARAMETER :: NRG=10 
 CHARACTER*1  OPT(0:NRG)
 CHARACTER*80 TITRE
 CHARACTER*30 NAMEIN, NAMEOUT, NAMEF   
 CHARACTER*10 TITLICE, RESOLUTION, DEGREE
 CHARACTER*30 VSTRING
 CHARACTER*20 FILE1_GMT
 CHARACTER*20 R_OPTION, T_OPTION
 CHARACTER*3  NV, CODE
 CHARACTER*1  ITER, MODE 
 CHARACTER*33 TABLE_NAME, TABLE_NAME_SU, TABLE_NAME_N, SCRIPT_NAME  
 CHARACTER*33 TMP_FILE, DATA_FILE, PS_FILE 
 CHARACTER*99 TMP_TITLE
 CHARACTER*100 SHORT_VISCO
 INTEGER K, U, US 
!
!
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This subroutine prepares a GMT script for plotting dot S, N & U 
! across specific regions. Regions available (as of March 2008): 
!  
! -1) Italy  
! -2) Mediterranean 
! -3) Europe 
! -4) Fennoscandia 
! -5) Greenland 
! -6) North America  
! -7) Antarctica 
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!
!
  If(opt(1)=='y') then 
!
! ################      ################     ################     ################
! Region 1 - ITALY      Region 1 - ITALY     Region 1 - ITALY	  Region 1 - ITALY
! ################      ################     ################	  ################
!
  table_name = "pale-italy.cpt" ; script_name = "italy.gmt"
!   
  Call color_tables (table_name)  
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"${GMT} gmtset PAPER_MEDIA A4+"
  Write(us,*)"${GMT} gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"${GMT} gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 100 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-italy.tmp" ; tmp_title= "Rate of sea level change today - ITALY"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-italy.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-italy.tmp" ; tmp_title= "Rate of vertical uplift today - ITALY"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-italy.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-italy.tmp" ; tmp_title= "Rate of sea surface variation today - ITALY"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-italy.ps"
  Endif 
!
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"13 +49.0 24 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
	   CONTINUE
  Else
!
  Write(u,*)"13 +48.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"13 +30.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
!
!
  Write(us,*)"${GMT} psbasemap -P -X4 -Y10 -Ba2f2WSEn -R5/21/35/47 -JJ13/14 -K > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*) "${GMT} pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), & 
              " -C"//trim(adjustl(table_name)), " >> ",  trim(adjustl(ps_file))	    
  Write(us,*)"${GMT} pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} psscale -E -U/0.5/0.5/'SELEN 2.9' ", "-C"//trim(adjustl(table_name)), &
              " -Bf0.1a0.5/:mm/yr: -D7/-2/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)" "
!
  100 Continue 
!
  close(us)   
!
  ENDIF
!
! ----------------      ----------------     ----------------     ----------------
! Region 1 - ITALY      Region 1 - ITALY     Region 1 - ITALY	  Region 1 - ITALY
! ----------------      ----------------     ----------------	  ----------------
!
  If(opt(2)=='y') then 
!
! ########################     ########################    ########################
! Region 2 - MEDITERRANEAN     Region 2 - MEDITERRANEAN    Region 2 - MEDITERRANEAN
! ########################     ########################    ########################
!
  table_name = "pale-mediterranean.cpt" ; script_name = "mediterranean.gmt"
!   
  Call color_tables (table_name)  
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"${GMT} gmtset PAPER_MEDIA A4+"
  Write(us,*)"${GMT} gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"${GMT} gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 200 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-mediterranean.tmp" ; tmp_title= "Rate of sea level change today - MEDITERRANEAN"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-mediterranean.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-mediterranean.tmp" ; tmp_title= "Rate of vertical uplift today - MEDITERRANEAN"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-mediterranean.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-mediterranean.tmp" ; tmp_title= "Rate of sea surface variation today - MEDITERRANEAN"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-mediterranean.ps"
  Endif 
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"20 +53.0   22 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
	   CONTINUE
  Else
!
  Write(u,*)"20 +51.5 16 0 0 BC", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"20 +21.0 16 0 0 BC", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
!
!  
  Write(us,*)"${GMT} psbasemap -X3 -Y6 -Ba8f4WSEn -R-8/48/30/50 -JJ18/24 -K > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*)"${GMT} pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))	     	      
  Write(us,*)"${GMT} pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} psscale -E -U/0.5/0.5/'SELEN 2.9' ", "-C"//trim(adjustl(table_name)), &
			  " -Bf0.1a0.5/:mm/yr: -D12/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)" "
!
  200 Continue 
!
  close(us)  
!
  ENDIF  
!
! ------------------------     ------------------------    ------------------------
! Region 2 - MEDITERRANEAN     Region 2 - MEDITERRANEAN    Region 2 - MEDITERRANEAN
! ------------------------     ------------------------    ------------------------
!
  If(opt(3)=='y') then 
!
! #################    #################    #################    #################
! Region 3 - EUROPE    Region 3 - EUROPE    Region 3 - EUROPE	 Region 3 - EUROPE
! #################    #################    #################	 #################
!
  table_name = "pale-europe.cpt" ; script_name = "europe.gmt"
!   
  Call color_tables (table_name)  
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"${GMT} gmtset PAPER_MEDIA A4+"
  Write(us,*)"${GMT} gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"${GMT} gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 300 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-europe.tmp" ; tmp_title= "Rate of sea level change today - EUROPE"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-europe.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-europe.tmp" ; tmp_title= "Rate of vertical uplift today - EUROPE"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-europe.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-europe.tmp" ; tmp_title= "Rate of sea surface variation today - EUROPE"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-europe.ps"
  Endif 
!
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"20 +80.0 26 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
	   CONTINUE
  Else
!
  Write(u,*)"20 +77.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"20 -2.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
!
!
!  
  Write(us,*)"${GMT} psbasemap -X3 -Y5 -Ba20f10/a10f5WSEn -R-40/80/25/75 -JJ20/18 -K > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*)"${GMT} pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} psscale -E -U/0.5/0.5/'SELEN 2.9' ", "-C"//trim(adjustl(table_name)), &
			  " -Bf1a4/:mm/yr: -D9/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)" "
!
  300 Continue 
!
  close(us)  
!
  ENDIF
!
! -----------------    -----------------    -----------------    -----------------
! Region 3 - EUROPE    Region 3 - EUROPE    Region 3 - EUROPE	 Region 3 - EUROPE
! -----------------    -----------------    -----------------	 -----------------
!
  If(opt(4)=='y') then 
!
! #######################  #######################  #######################
! Region 4 - FENNOSCANDIA  Region 4 - FENNOSCANDIA  Region 4 - FENNOSCANDIA
! #######################  #######################  #######################
!
  table_name_su = "pale-fennoscandia-su.cpt" 
  Call color_tables (table_name_su)  
!
  table_name_n = "pale-fennoscandia-n.cpt" 
  Call color_tables (table_name_n) 
!
  script_name = "fennoscandia.gmt"
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"${GMT} gmtset PAPER_MEDIA A4+"
  Write(us,*)"${GMT} gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"${GMT} gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 400 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-fennoscandia.tmp" ; tmp_title= "Rate of sea level change today - FENNOSCANDIA"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-fennoscandia.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-fennoscandia.tmp" ; tmp_title= "Rate of vertical uplift today - FENNOSCANDIA"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-fennoscandia.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-fennoscandia.tmp" ; tmp_title= "Rate of sea surface variation today - FENNOSCANDIA"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-fennoscandia.ps"
  Endif 
!
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"30 +77.0 24 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
	   CONTINUE
  Else
!
  Write(u,*)"30 +76.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"30 +40.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................

!
!  
  if(u==12.or.u==13) table_name=table_name_su
  if(u==14) 	     table_name=table_name_n
!
  Write(us,*)"${GMT} psbasemap -X3 -Y5 -Ba10f5WSEn -R0/60/50/75 -JJ26/20 -K  > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*)"${GMT} pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  if(u==12.or.u==13) &   
  Write(us,*)"${GMT} psscale -E -U/0.5/0.5/'SELEN 2.9' ", "-C"//trim(adjustl(table_name)), &
              " -Bf1a5/:mm/yr: -D10/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  if(u==14) &   
  Write(us,*)"${GMT} psscale -E -U/0.5/0.5/'SELEN 2.9' ", "-C"//trim(adjustl(table_name)), &
              " -Bf2a1/:mm/yr: -D10/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))	      	      	      	      
  Write(us,*)" "
!
  400 Continue 
!
  close(us) 
!
  ENDIF
!
! -----------------------   -----------------------   -----------------------
! Region 4 - Fennoscandia   Region 4 - Fennoscandia   Region 4 - Fennoscandia
! -----------------------   -----------------------   -----------------------
!
  If(opt(5)=='y') then 
!
! ####################       ####################       ####################
! Region 5 - GREENLAND       Region 5 - GREENLAND       Region 5 - GREENLAND
! ####################       ####################       ####################
!
  table_name_su = "pale-greenland-su.cpt" 
  Call color_tables (table_name_su)  
!
  table_name_n = "pale-greenland-n.cpt" 
  Call color_tables (table_name_n) 
!
  script_name = "greenland.gmt"
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"${GMT} gmtset PAPER_MEDIA A4+"
  Write(us,*)"${GMT} gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"${GMT} gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 500 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-greenland.tmp" ; tmp_title= "Rate of sea level change today - GREENLAND"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-greenland.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-greenland.tmp" ; tmp_title= "Rate of vertical uplift today - GREENLAND"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-greenland.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-greenland.tmp" ; tmp_title= "Rate of sea surface variation today - GREENLAND"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-greenland.ps"
  Endif 
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"320  91   22 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
	   CONTINUE
  Else
!
  Write(u,*)"320  89.5 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"320  30   16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
! 
! 
  if(u==12.or.u==13) table_name=table_name_su
  if(u==14) 	     table_name=table_name_n
!
  Write(us,*)"${GMT} psbasemap -P -X3 -Y9 -Ba10f10/f10a10WSEn -R-80/0/50/87.5 -JJ-40/16 -K  > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*)"${GMT} pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  if(u==12.or.u==13) &   
  Write(us,*)"${GMT} psscale -U/0.5/13/'SELEN 2.9' -E ", "-C"//trim(adjustl(table_name)), &
              " -Bf1a5/:mm/yr: -D8/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  if(u==14) &   
  Write(us,*)"${GMT} psscale -U/0.5/13/'SELEN 2.9' -E ", "-C"//trim(adjustl(table_name)), &
              " -Bf.1a.5/:mm/yr: -D8/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))	      	      	      	      
  Write(us,*)" "
!
  500 Continue 
!
  close(us) 
!
  ENDIF
!
! --------------------       --------------------       --------------------
! Region 5 - GREENLAND       Region 5 - GREENLAND       Region 5 - GREENLAND
! --------------------       --------------------       --------------------
!
  If(opt(6)=='y') then   
!
! ########################    ########################    ######################## 
! Region 6 - North AMERICA    Region 6 - North AMERICA    Region 6 - North AMERICA 
! ########################    ########################    ######################## 
!
  table_name_su = "pale-namerica-su.cpt" 
  Call color_tables (table_name_su)  
!
  table_name_n = "pale-namerica-n.cpt" 
  Call color_tables (table_name_n) 
!
  script_name = "north-america.gmt"
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"${GMT} gmtset PAPER_MEDIA A4+"
  Write(us,*)"${GMT} gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"${GMT} gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 600 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-namerica.tmp" ; tmp_title= "Rate of sea level change today - N. AMERICA"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-namerica.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-namerica.tmp" ; tmp_title= "Rate of vertical uplift today - N. AMERICA"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-namerica.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-namerica.tmp" ; tmp_title= "Rate of sea surface variation today - N. AMERICA"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-namerica.ps"
  Endif 
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"270  92 20 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
	   CONTINUE
  Else
!
  Write(u,*)"270  90 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"270 -30 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
!
!  
  if(u==12.or.u==13) table_name=table_name_su
  if(u==14) 	     table_name=table_name_n
!
  Write(us,*)"${GMT} psbasemap -P -X3 -Y9 -Ba20f10/f10a10WSEn -R200/340/10/87.5 -JJ-90/16 -K  > ", & 
  	      trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscontour -I -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file)) 
  Write(us,*)"${GMT} pscontour -G16 -W1/255 -A -JJ -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscoast -R -JJ -Di -B -W2/0 -A1000 -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pstext -N -R -JJ -B -G0 -O -K ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))
  if(u==12.or.u==13) &   
  Write(us,*)"${GMT} psscale -U/0.5/0.5/'SELEN 2.9' -E ", "-C"//trim(adjustl(table_name)), &
              " -Bf1a4/:mm/yr: -D8/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))
  if(u==14) &   
  Write(us,*)"${GMT} psscale -U/0.5/0.5/'SELEN 2.9' -E ", "-C"//trim(adjustl(table_name)), &
              " -Bf1a2/:mm/yr: -D8/-1.5/10/1h -O >> ", & 
              trim(adjustl(ps_file))	      	      	      	      
  Write(us,*)" "
!
  600 Continue 
!
  close(us) 
!
  ENDIF
!
! ------------------------    ------------------------    ------------------------ 
! Region 6 - North AMERICA    Region 6 - North AMERICA    Region 6 - North AMERICA 
! ------------------------    ------------------------    ------------------------ 
!
  If(opt(7)=='y') then   
!
! #####################    #####################    #####################
! Region 7 - Antarctica    Region 7 - Antarctica    Region 7 - Antarctica
! #####################    #####################    #####################
!
  table_name_su = "pale-antarctica-su.cpt" 
  Call color_tables (table_name_su)  
!
  table_name_n = "pale-antarctica-n.cpt" 
  Call color_tables (table_name_n) 
!
  script_name = "antarctica.gmt"
!
  us=11 
  Open (us,file=script_name,status='unknown') 
  Write(us,*)"${GMT} gmtset PAPER_MEDIA A4+"
  Write(us,*)"${GMT} gmtset LABEL_FONT_SIZE 24p"
  Write(us,*)"${GMT} gmtset ANOT_FONT_SIZE 16p" 
  Write(us,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
  Write(us,*)" " 
!
  do 700 u=12, 14 
! 
  if    (u==12) then 
  	    tmp_file = "sdot-antarctica.tmp" ; tmp_title= "Rate of sea level change today - ANTARCTICA"
	    data_file= "sdotmap.dat"    ; ps_file  = "sdot-antarctica.ps"  
  elseif(u==13) then
  	    tmp_file = "udot-antarctica.tmp" ; tmp_title= "Rate of vertical uplift today - ANTARCTICA"
	    data_file= "udotmap.dat"    ; ps_file  = "udot-antarctica.ps"
  elseif(u==14) then
  	    tmp_file = "ndot-antarctica.tmp" ; tmp_title= "Rate of sea surface variation today - ANTARCTICA"
	    data_file= "ndotmap.dat"    ; ps_file  = "ndot-antarctica.ps"
  Endif 
!
! ..............................................................................
  Open (u,file=tmp_file,status='unknown') 
!
  Write(u,*)"0 +18.0 22 0 3 BC ", trim(adjustl(tmp_title))
!
  If(CODE=='-1') then 
	   CONTINUE
  Else
!
  Write(u,*)"0 +10.0 16 0 0 BC ", & 
  		               " -Ice model: ", trim(adjustl(titlice)), & 
  			       " -Viscosity profile: ", trim(adjustl(vstring))
  Write(u,*)"0 -115.0 16 0 0 BC ", & 
  	                       " -LMAX=", trim(adjustl(DEGREE)), & 
			       " -RES=",  trim(adjustl(RESOLUTION)), & 
			       " -NV=",   trim(adjustl(NV)), & 
			       " -CODE=", trim(adjustl(CODE)), & 
			       " -MODE=", trim(adjustl(MODE)), & 
			       " -ITER=", trim(adjustl(ITER))	
  Endif
!
  Close(u)			       		       
! ..............................................................................
!
!
  if(u==12.or.u==13) table_name=table_name_su
  if(u==14) 	     table_name=table_name_n
!
  Write(us,*)"${GMT} pscontour -P -X4 -Y8 -I -JE0/-90/12 -R0/360/-90/-52 -K ", trim(adjustl(data_file)), & 
  	     " -C"//trim(adjustl(table_name)), " > ", trim(adjustl(ps_file)) 
  Write(us,*)"${GMT} pscontour -G15 -W1/255 -A -JE -R -O -K ", trim(adjustl(data_file)), " -C"//trim(adjustl(table_name)), " >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pscoast -R -JE -Di -Ba30f30g30 -W2/255/0/0 -A0 -O -K >> ", & 
              trim(adjustl(ps_file))
  if(u==12.or.u==13) &   
  Write(us,*)"${GMT} psscale -E -X-6 -Y-7 ", "-C"//trim(adjustl(table_name)), &
              " -Bf1a2/:mm/yr: -D12/5/8/1h -O -K >> ", & 
              trim(adjustl(ps_file))
  if(u==14) &   
  Write(us,*)"${GMT} psscale -E -X-6 -Y-7 ", "-C"//trim(adjustl(table_name)), &
              " -Bf0.5a1/:mm/yr: -D12/5/8/1h -O -K >> ", & 
              trim(adjustl(ps_file))
  Write(us,*)"${GMT} pstext -X6 -N -R -JE -G0 -O -U/4/20.5/'SELEN 2.9' ", trim(adjustl(tmp_file)), " >> ", & 
              trim(adjustl(ps_file))	      	      	      	      
  Write(us,*)" "
!
  700 Continue 
!
  close(us) 
!
  ENDIF
!
! ---------------------    ---------------------    ---------------------
! Region 7 - Antarctica    Region 7 - Antarctica    Region 7 - Antarctica
! ---------------------    ---------------------    ---------------------
!
 END SUBROUTINE MAKE_RMAPS
!
!
!
!
!
!
 SUBROUTINE COLOR_TABLES (TABLE_NAME)
 IMPLICIT NONE
 CHARACTER*33 TABLE_NAME
 CHARACTER*5 FMT1, FMT2
!
! # Draws color tables for regional maps- 
!
 if(table_name=='pale-italy.cpt'          .or.& 
    table_name=='pale-mediterranean.cpt')     & 
    then 
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#    cpt created by: makecpt -Cno_green -T-1.5/1.5/0.1     '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-1.5    32      96      255     -1.4    32      96      255'
 	Write(9,fmt1) '-1.4    32      96      255     -1.3    32      96      255'
 	Write(9,fmt1) '-1.3    32      159     255     -1.2    32      159     255'
 	Write(9,fmt1) '-1.2    32      159     255     -1.1    32      159     255'
 	Write(9,fmt1) '-1.1    32      191     255     -1      32      191     255'
 	Write(9,fmt1) '-1      32      191     255     -0.9    32      191     255'
 	Write(9,fmt1) '-0.9    0       207     255     -0.8    0       207     255'
 	Write(9,fmt1) '-0.8    0       207     255     -0.7    0       207     255'
 	Write(9,fmt1) '-0.7    42      255     255     -0.6    42      255     255'
 	Write(9,fmt1) '-0.6    42      255     255     -0.5    42      255     255'
 	Write(9,fmt1) '-0.5    85      255     255     -0.4    85      255     255'
 	Write(9,fmt1) '-0.4    85      255     255     -0.3    85      255     255'
 	Write(9,fmt1) '-0.3    127     255     255     -0.2    127     255     255'
 	Write(9,fmt1) '-0.2    127     255     255     -0.1    127     255     255'
 	Write(9,fmt1) '-0.1    170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      0.1     255     255     84 '
 	Write(9,fmt1) '0.1     255     255     84      0.2     255     255     84 '
 	Write(9,fmt1) '0.2     255     240     0       0.3     255     240     0  '
 	Write(9,fmt1) '0.3     255     240     0       0.4     255     240     0  '
 	Write(9,fmt1) '0.4     255     191     0       0.5     255     191     0  '
	Write(9,fmt1) '0.5     255     191     0       0.6     255     191     0  '
	Write(9,fmt1) '0.6     255     168     0       0.7     255     168     0  '
	Write(9,fmt1) '0.7     255     168     0       0.8     255     168     0  '
	Write(9,fmt1) '0.8     255     138     0       0.9     255     138     0  '
	Write(9,fmt1) '0.9     255     138     0       1       255     138     0  '
	Write(9,fmt1) '1       255     112     0       1.1     255     112     0  '
	Write(9,fmt1) '1.1     255     112     0       1.2     255     112     0  '
	Write(9,fmt1) '1.2     255     77      0       1.3     255     77      0  '
	Write(9,fmt1) '1.3     255     77      0       1.4     255     77      0  '
	Write(9,fmt1) '1.4     255     0       0       1.5     255     0       0  '
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9)
!
 	Elseif(table_name=='pale-europe.cpt')  then
! 
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-8/8/1         '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-8      32      96      255     -7      32      96      255'
 	Write(9,fmt1) '-7      32      159     255     -6      32      159     255'
 	Write(9,fmt1) '-6      32      191     255     -5      32      191     255'
 	Write(9,fmt1) '-5      0       207     255     -4      0       207     255'
 	Write(9,fmt1) '-4      42      255     255     -3      42      255     255'
 	Write(9,fmt1) '-3      85      255     255     -2      85      255     255'
 	Write(9,fmt1) '-2      127     255     255     -1      127     255     255'
 	Write(9,fmt1) '-1      170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      1       255     255     84 '
 	Write(9,fmt1) '1       255     240     0       2       255     240     0  '
 	Write(9,fmt1) '2       255     191     0       3       255     191     0  '
 	Write(9,fmt1) '3       255     168     0       4       255     168     0  '
 	Write(9,fmt1) '4       255     138     0       5       255     138     0  '
 	Write(9,fmt1) '5       255     112     0       6       255     112     0  '
 	Write(9,fmt1) '6       255     77      0       7       255     77      0  '
 	Write(9,fmt1) '7       255     0       0       8       255     0       0  '
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
! 
 	Elseif(table_name=='pale-fennoscandia-su.cpt'.or.    & 
	       table_name=='pale-greenland-su.cpt'.or. &  
	       table_name=='pale-antarctica-su.cpt')       then
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-10/10/1       '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-10     32      96      255     -9      32      96      255'
 	Write(9,fmt1) '-9      32      96      255     -8      32      96      255'
 	Write(9,fmt1) '-8      32      159     255     -7      32      159     255'
 	Write(9,fmt1) '-7      32      191     255     -6      32      191     255'
 	Write(9,fmt1) '-6      0       207     255     -5      0       207     255'
 	Write(9,fmt1) '-5      42      255     255     -4      42      255     255'
 	Write(9,fmt1) '-4      42      255     255     -3      42      255     255'
 	Write(9,fmt1) '-3      85      255     255     -2      85      255     255'
 	Write(9,fmt1) '-2      127     255     255     -1      127     255     255'
 	Write(9,fmt1) '-1      170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      1       255     255     84 '
 	Write(9,fmt1) '1       255     255     84      2       255     255     84 '
 	Write(9,fmt1) '2       255     240     0       3       255     240     0  '
 	Write(9,fmt1) '3       255     191     0       4       255     191     0  '
 	Write(9,fmt1) '4       255     168     0       5       255     168     0  '
 	Write(9,fmt1) '5       255     138     0       6       255     138     0  '
 	Write(9,fmt1) '6       255     138     0       7       255     138     0  '
 	Write(9,fmt1) '7       255     112     0       8       255     112     0  '
 	Write(9,fmt1) '8       255     77      0       9       255     77      0  '
 	Write(9,fmt1) '9       255     0       0       10      255     0       0  '
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
!
 	Elseif(table_name=='pale-fennoscandia-n.cpt'.or.   & 
	       table_name=='pale-antarctica-n.cpt')    then
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-2/2/0.2       '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-2      32      96      255     -1.8    32      96      255'
 	Write(9,fmt1) '-1.8    32      96      255     -1.6    32      96      255'
 	Write(9,fmt1) '-1.6    32      159     255     -1.4    32      159     255'
 	Write(9,fmt1) '-1.4    32      191     255     -1.2    32      191     255'
 	Write(9,fmt1) '-1.2    0       207     255     -1      0       207     255'
 	Write(9,fmt1) '-1      42      255     255     -0.8    42      255     255'
 	Write(9,fmt1) '-0.8    42      255     255     -0.6    42      255     255'
 	Write(9,fmt1) '-0.6    85      255     255     -0.4    85      255     255'
 	Write(9,fmt1) '-0.4    127     255     255     -0.2    127     255     255'
 	Write(9,fmt1) '-0.2    170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      0.2     255     255     84 '
 	Write(9,fmt1) '0.2     255     255     84      0.4     255     255     84 '
 	Write(9,fmt1) '0.4     255     240     0       0.6     255     240     0  '
 	Write(9,fmt1) '0.6     255     191     0       0.8     255     191     0  '
 	Write(9,fmt1) '0.8     255     168     0       1       255     168     0  '
 	Write(9,fmt1) '1       255     138     0       1.2     255     138     0  '
 	Write(9,fmt1) '1.2     255     138     0       1.4     255     138     0  '
 	Write(9,fmt1) '1.4     255     112     0       1.6     255     112     0  '
 	Write(9,fmt1) '1.6     255     77      0       1.8     255     77      0  '
 	Write(9,fmt1) '1.8     255     0       0       2       255     0       0  '
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
!
 	Elseif(table_name=='pale-greenland-n.cpt')  then
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-1/1/0.2       '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-1      32      96      255     -0.8    32      96      255'
 	Write(9,fmt1) '-0.8    32      159     255     -0.6    32      159     255'
 	Write(9,fmt1) '-0.6    0       207     255     -0.4    0       207     255'
 	Write(9,fmt1) '-0.4    42      255     255     -0.2    42      255     255'
 	Write(9,fmt1) '-0.2    127     255     255     0       127     255     255'
 	Write(9,fmt1) '0       255     255     84      0.2     255     255     84 '
 	Write(9,fmt1) '0.2     255     240     0       0.4     255     240     0  '
 	Write(9,fmt1) '0.4     255     168     0       0.6     255     168     0  '
 	Write(9,fmt1) '0.6     255     138     0       0.8     255     138     0  '
 	Write(9,fmt1) '0.8     255     77      0       1       255     77      0  '
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
!
 	Elseif(table_name=='pale-namerica-su.cpt')  then
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-24/24/2       '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-24     32      96      255     -22     32      96      255'
 	Write(9,fmt1) '-22     32      96      255     -20     32      96      255'
 	Write(9,fmt1) '-20     32      159     255     -18     32      159     255'
 	Write(9,fmt1) '-18     32      191     255     -16     32      191     255'
 	Write(9,fmt1) '-16     32      191     255     -14     32      191     255'
 	Write(9,fmt1) '-14     0       207     255     -12     0       207     255'
 	Write(9,fmt1) '-12     42      255     255     -10     42      255     255'
 	Write(9,fmt1) '-10     42      255     255     -8      42      255     255'
 	Write(9,fmt1) '-8      85      255     255     -6      85      255     255'
 	Write(9,fmt1) '-6      127     255     255     -4      127     255     255'
 	Write(9,fmt1) '-4      127     255     255     -2      127     255     255'
 	Write(9,fmt1) '-2      170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      2       255     255     84 '
 	Write(9,fmt1) '2       255     255     84      4       255     255     84 '
 	Write(9,fmt1) '4       255     240     0       6       255     240     0  '
 	Write(9,fmt1) '6       255     191     0       8       255     191     0  '
 	Write(9,fmt1) '8       255     191     0       10      255     191     0  '
 	Write(9,fmt1) '10      255     168     0       12      255     168     0  '
 	Write(9,fmt1) '12      255     138     0       14      255     138     0  '
 	Write(9,fmt1) '14      255     138     0       16      255     138     0  '
	Write(9,fmt1) '16      255     112     0       18      255     112     0  '
	Write(9,fmt1) '18      255     77      0       20      255     77      0  '
	Write(9,fmt1) '20      255     77      0       22      255     77      0  '
	Write(9,fmt1) '22      255     0       0       24      255     0       0  ' 
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
!
 	Elseif(table_name=='pale-namerica-n.cpt')  then
!
 	open(9,file=table_name,status='unknown') 
	fmt1='(a59)'
	fmt2='(a27)'
	Write(9,fmt1) '#----------------------------------------------------------'
 	Write(9,fmt1) '#      cpt created by: makecpt -Cno_green -T-10/10/1       '
	Write(9,fmt1) '#----------------------------------------------------------'	
 	Write(9,fmt1) '-10     32      96      255     -9      32      96      255'
 	Write(9,fmt1) '-9      32      96      255     -8      32      96      255'
 	Write(9,fmt1) '-8      32      159     255     -7      32      159     255'
 	Write(9,fmt1) '-7      32      191     255     -6      32      191     255'
 	Write(9,fmt1) '-6      0       207     255     -5      0       207     255'
 	Write(9,fmt1) '-5      42      255     255     -4      42      255     255'
 	Write(9,fmt1) '-4      42      255     255     -3      42      255     255'
 	Write(9,fmt1) '-3      85      255     255     -2      85      255     255'
 	Write(9,fmt1) '-2      127     255     255     -1      127     255     255'
 	Write(9,fmt1) '-1      170     255     255     0       170     255     255'
 	Write(9,fmt1) '0       255     255     84      1       255     255     84 '
 	Write(9,fmt1) '1       255     255     84      2       255     255     84 '
 	Write(9,fmt1) '2       255     240     0       3       255     240     0  '
 	Write(9,fmt1) '3       255     191     0       4       255     191     0  '
 	Write(9,fmt1) '4       255     168     0       5       255     168     0  '
 	Write(9,fmt1) '5       255     138     0       6       255     138     0  '
 	Write(9,fmt1) '6       255     138     0       7       255     138     0  '
 	Write(9,fmt1) '7       255     112     0       8       255     112     0  '
 	Write(9,fmt1) '8       255     77      0       9       255     77      0  '
 	Write(9,fmt1) '9       255     0       0       10      255     0       0  ' 
 	Write(9,fmt2) 'B       32      96      255'
 	Write(9,fmt2) 'F       255     0         0'
 	Write(9,fmt2) 'N       128     128     128'
 	close(9) 
!
 Endif
!
 END SUBROUTINE COLOR_TABLES
!
!
!
!
!
!
 SUBROUTINE MAKE_GMAPS (TITLICE, & 
 			RESOLUTION, &
			NV, & 
			CODE, & 
			ITER, & 
 	                MODE, & 
			DEGREE, & 
			VSTRING, & 
			OPTION_ROF, & 
			FILECAP, & 
			FILE1_GMT, &
			SHORT_VISCO)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Prepares a GMT script for *** GLOBAL MAPS *** maps of dot-S, U & N at present time 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Last change GS March 24 08 
! Modified on June 24, 2008 for SELEN 2.6 GS 
! Modified on July 26, 2008 for SELEN 2.6 GS 
! Upgraded April 2010 by GS for ALMA on g95 
!
!
 IMPLICIT NONE
 CHARACTER*100 SHORT_VISCO 
 CHARACTER*80  TITRE
 CHARACTER*30  NAMEIN, NAMEOUT, NAMEF   
 CHARACTER*10  TITLICE, RESOLUTION, DEGREE
 CHARACTER*30  VSTRING
 CHARACTER*20  FILECAP, FILE1_GMT
 CHARACTER*20  R_OPTION, T_OPTION
 CHARACTER*3   NV, CODE
 CHARACTER*1   OPTION_ROF
 CHARACTER*1   ITER, MODE 
 INTEGER K
!
!
!
 open(9,file=file1_gmt,status='unknown')
 Write(9,*)"${GMT} gmtset PAPER_MEDIA A4+"
 Write(9,*)"${GMT} gmtset LABEL_FONT_SIZE 24p"
 Write(9,*)"${GMT} gmtset ANOT_FONT_SIZE 16p"
 Write(9,*)"${GMT} gmtset FRAME_WIDTH 0.1c"

 Write(9,*)" "
 Write(9,*)"# ---- Global maps of dot S, U, and N at present time  ----" 		      
!
 R_OPTION="-R0/360/-80/80"      ! Range of the *Mercator* map common to all maps 
!
do 1 k=1, 3
 	if    (k==1) then 
	       Write(9,*) ""
	       namein ="sdotmap.dat"; nameout="sdotmap.ps"; namef="tmpfs.dat"
!
	       TITRE="Rate of sea level change today"
        elseif(k==2) then 
	       Write(9,*) ""
	       namein ="udotmap.dat"; nameout="udotmap.ps"; namef="tmpfu.dat"
!
	       TITRE="Rate of vertical displacement today"
        elseif(k==3) then 
	       Write(9,*) ""			
	       namein ="ndotmap.dat"; nameout="ndotmap.ps"; namef="tmpfn.dat"
!
	       TITRE="Rate of sea surface variation today"	
	endif
!
T_OPTION="-T-1.0/1.0/0.2"         ! Range of the palette 
!
Write(9,*) " "
Write(9,*) "${GMT} makecpt -Cno_green ", trim(adjustl(T_OPTION)), " > pale.cpt" 
Write(9,*) "${GMT} psbasemap -X3 -Y5 -Ba180/a40WSEn -Jm0.018i ", & 
            trim(adjustl(R_OPTION)), " -K > ", " ", trim(adjustl(nameout))  
Write(9,*) "${GMT} pscontour -I -Jm -O -K ", trim(adjustl(R_OPTION)), & 
           " ", trim(adjustl(namein)), " -Cpale.cpt  >> ", trim(adjustl(nameout)) 
If(option_rof=='r')then 
  Write(9,*) "${GMT} pscoast -Jm -Dc -B -W2/0       -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
  else
  Write(9,*) "${GMT} pscoast -Jm -Dc -B -W1/240 -A1000 -O -K ", trim(adjustl(R_OPTION)), " >> ", trim(adjustl(nameout))  
Endif
If(option_rof=='z')& 
Write(9,*) "${GMT} psxy -R -Jm ", trim(adjustl(filecap)), " -M -W4/0 -B -A -O -K >> ", trim(adjustl(nameout))  
!
open (4,file=NAMEF,status='unknown')  

If(CODE=='-1')    THEN 
	   CONTINUE
ELSEIF(CODE/='-1')THEN 
Write(4,*) "180 +84.2 26 0 3 BC ", trim(adjustl(TITRE)) 
Write(4,*) "180 +82.1 16 0 0 BC  -Ice model: ", trim(adjustl(TITLICE)),& 
 	   "  -Viscosity profile: ", trim(adjustl(VSTRING))  
Write(4,*) "180 -87.5   13  0 0 BC ",& 
 	    " -LMAX=", trim(adjustl(DEGREE)),& 
 	    " -RES=", trim(adjustl(RESOLUTION)),& 
	    " -NV=", trim(adjustl(NV)),& 
	    " -CODE=", trim(adjustl(CODE)),& 
	    " -MODE=", trim(adjustl(MODE)),&   
	    " -ITER=", trim(adjustl(ITER));  close(4) 
Write(9,*) "${GMT} pstext -N -R -Jm -B ", trim(adjustl(NAMEF)), " -G0 -O -K >> ", trim(adjustl(nameout)) 
Write(9,*) "${GMT} psscale -U/0.5/0.5/'SELEN 2.9' -E -Cpale.cpt -B1f0.25a0.5/:mm/yr: ", &
           "-D8.25/-1/10/1h -O >> ", trim(adjustl(nameout))
ENDIF
!
1 CONTINUE 
!
 close(9) 		      
!
 END SUBROUTINE MAKE_GMAPS		      
!
!
!
!
!
!
 SUBROUTINE MAKE_WNW (RESOLUTION, DEGREE, FILE1_GMT)
 IMPLICIT NONE
 CHARACTER*3 RESOLUTION, DEGREE
 CHARACTER*20 R_OPTION, FILE1_GMT 
 CHARACTER*6 NPC
 CHARACTER*5 JMAXC
 INTEGER LMAX, NP, RES, JMAX, JMAXP  
!
! ---- Prepares a GMT script for plotting the window function
!         Last modification - Giorgio Spada October 2007 
!
 open(9,file=file1_gmt,status='unknown')
 Write(9,*) "${GMT} gmtset PAPER_MEDIA A4+"
 Write(9,*) "${GMT} gmtset ANOT_FONT_SIZE 16p"
 Write(9,*) "${GMT} gmtset LABEL_FONT_SIZE 26p"
 Write(9,*) "#"
 Write(9,*) "# ---- Window function ----" 
 open(3,file='tmpw0.dat',status='unknown'); Write(3,'(a3)')degree; close(3) 
 open(3,file='tmpw0.dat',status='unknown'); Read(3,*)lmax; close(3) 
 open(3,file='tmpw1.dat',status='unknown'); Write(3,'(a3)')resolution; close(3) 
 open(3,file='tmpw1.dat',status='unknown'); Read(3,*)res; close(3) 
 NP=2*RES*(RES-1)*20+12 
 JMAX=(LMAX+1)*(LMAX+2)/2
 JMAXP=JMAX+30
 open(3,file='tmpw2.dat',status='unknown'); Write(3,'(i5)') jmaxp; close(3) 
 open(3,file='tmpw2.dat',status='unknown'); Read (3,'(a5)') jmaxc; close(3) 
 open(3,file='tmpw3.dat',status='unknown'); Write(3,'(i5)') NP;    close(3) 
 open(3,file='tmpw3.dat',status='unknown'); Read (3,'(a5)') NPC;   close(3) 
!
 R_OPTION="-R3.5/"//trim(adjustl(JMAXC))//"/1e-4/1e2"
!
 Write(9,*) "${GMT} psbasemap -X6 -Y4 -U/0.5/12.75/'SELEN 2.9' ", R_OPTION, & 
            "  -Ba1f1p:'degree*(degree+1)/2+order+1':/f1a10:@~e@~,%::,%::.'Window function':WSne -JX-18l/14l -K > wnw.ps"
 open (4,file='tmpp4.dat',status='unknown') 
 Write(4,*) " 4 3e1 16 0 2 BR -RES=", trim(adjustl(RESOLUTION)), & 
            " -NP=", trim(adjustl(NPC)), " -LMAX= ", trim(adjustl(DEGREE)), "  " ; close(4)
 Write(9,*) "${GMT} pstext tmpp4.dat -N -JX -R -G0 -O -K >> wnw.ps"	     
 Write(9,*) "${GMT} psxy wnw.dat -B -R -JX -Sc0.1 -G0 -K -O >> wnw.ps"
 open (4,file='tmpp5.dat',status='unknown') 
 Write(4,*) JMAXP, " 1"
 Write(4,*)  " 1 1 " ; CLOSE(4) 
 Write(9,*) "${GMT} psxy tmpp5.dat -B -R -JX -W4ta -G0 -K -O >> wnw.ps"
 open (4,file='tmpp6.dat',status='unknown') 
 Write(4,*) JMAXP, " 5"
 Write(4,*)  " 1 5 " ; CLOSE(4) 
 Write(9,*) "${GMT} psxy tmpp6.dat -B -R -JX -W4to -G0 -K -O >> wnw.ps"
 open (4,file='tmpp7.dat',status='unknown') 
 Write(4,*) JMAXP, " 10"
 Write(4,*)  " 1 10 " ; CLOSE(4) 
 Write(9,*) "${GMT} psxy tmpp7.dat -B -R -JX -W2to -G0 -O >> wnw.ps"
!
 close(9) 
!
 END SUBROUTINE MAKE_WNW
!
!
!
!
!
!
 SUBROUTINE MAKE_PLOT_LDC (NV, CODE, DEGREE, VSTRING, FILE1_GMT)
 IMPLICIT NONE
 CHARACTER*3 NV, CODE
 CHARACTER*3 DEGREE
 CHARACTER*20 FILE1_GMT
 CHARACTER*22 R_OPTION1, R_OPTION2, R_OPTION3, R_OPTION4 
 CHARACTER*30 VSTRING
 INTEGER LMAX
!
! ---- Prepares the GMT script for plotting the LDCs and the spectrum 
!             Last modification Giorgio Spada November 2007 
!	      Revised  June 15, 2007 (V 2.6 INTEL) 
!             Revised May 2010 g95
!
!
! ---- Target GMT file for plotting the spectral properties of the chosen Earth model 
!
open(9,file=file1_gmt,status='unknown')
!	
Write(9,*) "${GMT} gmtset PAPER_MEDIA A4+"
Write(9,*) "${GMT} gmtset ANOT_FONT_SIZE 16p"
Write(9,*) "${GMT} gmtset LABEL_FONT_SIZE 18p"
Write(9,*) "#"
!
! ---- The x-range varies according to the maximum degree but the y-range is fixed ... 
open(3,file='tmpg0.dat',status='unknown'); Write(3,'(a3)')degree; close(3) 
open(3,file='tmpg0.dat',status='unknown'); Read(3,*)LMAX; close(3) 
! for the elastic LDC  	
if(0<=LMAX.and.LMAX<=10)    R_OPTION1=" -R0.8/10/-0.40/0.05"
if(0<=LMAX.and.LMAX<=100)   R_OPTION1="-R0.8/100/-0.40/0.05"
if(100<=LMAX.and.LMAX<=150) R_OPTION1="-R0.8/150/-0.40/0.05"
if(150<=LMAX.and.LMAX<=200) R_OPTION1="-R0.8/200/-0.40/0.05"
if(200<=LMAX.and.LMAX<=300) R_OPTION1="-R0.8/300/-0.40/0.05"
!
! for the fluid LDC
if(0<=LMAX.and.LMAX<=10)    R_OPTION2=" -R0.8/10/-1.25/0.25"
if(0<=LMAX.and.LMAX<=100)   R_OPTION2="-R0.8/100/-1.25/0.25"
if(100<=LMAX.and.LMAX<=150) R_OPTION2="-R0.8/150/-1.25/0.25"
if(150<=LMAX.and.LMAX<=200) R_OPTION2="-R0.8/200/-1.25/0.25"
if(200<=LMAX.and.LMAX<=300) R_OPTION2="-R0.8/300/-1.25/0.25"
!
! for the Spectrum 
if(  0<=LMAX.and.LMAX<=10)  R_OPTION3="-R0.8/10/1e1/1e9"
if(  0<=LMAX.and.LMAX<=100) R_OPTION3="-R0.8/100/1e1/1e9"
if(100<=LMAX.and.LMAX<=150) R_OPTION3="-R0.8/150/1e1/1e9"
if(150<=LMAX.and.LMAX<=200) R_OPTION3="-R0.8/200/1e1/1e9"
if(200<=LMAX.and.LMAX<=300) R_OPTION3="-R0.8/300/1e1/1e9"

! for the Normalized residues  
if(0<=LMAX.and.LMAX<=10)    R_OPTION4=" -R0.8/10/1e-7/1e1"
if(0<=LMAX.and.LMAX<=100)   R_OPTION4="-R0.8/100/1e-7/1e1"
if(100<=LMAX.and.LMAX<=150) R_OPTION4="-R0.8/150/1e-7/1e1"
if(150<=LMAX.and.LMAX<=200) R_OPTION4="-R0.8/200/1e-7/1e1"
if(200<=LMAX.and.LMAX<=300) R_OPTION4="-R0.8/300/1e-7/1e1"
!
Write(9,*) "#"
Write(9,*) "# ---------- Elastic LDCs vs degree (left frame of 'ela-flu.ps')"
Write(9,*) "${GMT} psbasemap -X5 -Y5 -U/8/-3/'SELEN 2.9' ", R_OPTION1, & 
	   " -Ba1f1p:'Harmonic degree, n':/f0.05a0.05:'Elastic LDC':WSen -JX9l/12 -K > ela-flu.ps" 
open (4,file='tmpg1.dat',status='unknown') 
Write(4,*) "1 0.09 18 0 2 BL Earth model: NV=", trim(adjustl(NV)), " CODE=", trim(adjustl(CODE))
Write(4,*) "1 0.07 16 0 2 BL Viscosity profile: ", trim(adjustl(VSTRING)), "*1E21 Pa.s"; close(4) 
Write(9,*) "${GMT} pstext tmpg1.dat -N -JX -R -G0 -CM -O -K >> ela-flu.ps"
!
Write(9,*) "awk '{print $1, $2/(2*$1+1)}' hh.dat > h.tmp" 
Write(9,*) "awk '{print $1, $2         }' ll.dat > l.tmp"  
Write(9,*) "awk '{print $1, $2         }' kk.dat > k.tmp" 
Write(9,*) "${GMT} psxy -H2 h.tmp -B -R -JX -Ss0.4  -G0  -K -O >> ela-flu.ps" 
Write(9,*) "${GMT} psxy -H2 l.tmp -B -R -JX -Sc0.35 -G0  -K -O >> ela-flu.ps" 
Write(9,*) "${GMT} psxy -H2 k.tmp -B -R -JX -Si0.4  -G0  -K -O >> ela-flu.ps" 
!
open (4,file='tmpg2.dat',status='unknown') 
Write(4,*) "20 -0.20 18 0 0 BL h/(2*n+1)"
Write(4,*) "20 -0.23 18 0 0 BL l" 
Write(4,*) "20 -0.26 18 0 0 BL k"; close(4) 
Write(9,*) "${GMT} pstext tmpg2.dat -N -JX -R -G0 -CM -O -K >> ela-flu.ps"
!
open (4,file='tmpg3.dat',status='unknown') 
Write(4,*) "15 -0.20"; close(4) 
Write(9,*) "${GMT} psxy -JX -R -G0 -Ss0.4  -O -K tmpg3.dat >> ela-flu.ps"
open (4,file='tmpg4.dat',status='unknown') 
Write(4,*) "15 -0.23"; close(4) 
Write(9,*) "${GMT} psxy -JX -R -G0 -Sc0.35 -O -K tmpg4.dat >> ela-flu.ps"
open (4,file='tmpg5.dat',status='unknown') 
Write(4,*) "15 -0.26"; close(4)
Write(9,*) "${GMT} psxy -JX -R -G0 -Si0.4  -O -K tmpg5.dat >> ela-flu.ps"
! 
Write(9,*) "#"
Write(9,*) "# ---------- Fluid LDCs vs degree (right frame of 'ela-flu.ps')"
Write(9,*) "${GMT} psbasemap -X11 ", R_OPTION2, & 
	    " -Ba1f1p:'Harmonic degree, n':/f0.25a0.25:'Fluid LDC':wSEn -JX -K -O >> ela-flu.ps" 
Write(9,*) "awk '{print $1, $3/(2*$1+1)}' hh.dat > h.tmp" 
Write(9,*) "awk '{print $1, $3         }' ll.dat > l.tmp"  
Write(9,*) "awk '{print $1, $3         }' kk.dat > k.tmp" 
Write(9,*) "${GMT} psxy -H2 h.tmp -B -R -JX -Ss0.4  -G0  -K -O >> ela-flu.ps" 
Write(9,*) "${GMT} psxy -H2 l.tmp -B -R -JX -Sc0.35 -G0  -K -O >> ela-flu.ps"
Write(9,*) "${GMT} psxy -H2 k.tmp -B -R -JX -Si0.4  -G0     -O >> ela-flu.ps"
!
Write(9,*) "#"
Write(9,*) "# ----  Relaxation spectrum (i. e., relaxation times vs. degree), 'spectrum.ps'"
Write(9,*) "${GMT} psbasemap -X5 -Y5 -U/0.5/11/'SELEN 2.9' ", R_OPTION3, & 
	   " -Ba1f1p:'Harmonic degree, n':/f1a1p:'Relaxation time (yrs)':WSen -JX9l/12l  -K >  spectrum.ps" 
open (4,file='tmpg6.dat',status='unknown') 
Write(4,*) "1 0.8e10 14 0 2 BL Earth model: NV=", trim(adjustl(NV)), " CODE=", trim(adjustl(CODE))
Write(4,*) "1 0.3e10 12 0 2 BL Viscosity profile: ", trim(adjustl(VSTRING)), "*1E21 Pa.s"; close(4) 
Write(9,*) "${GMT} pstext -N -JX -R -G0 -CM -O -K tmpg6.dat >> spectrum.ps"
Write(9,*) "awk '{print $1, $5}' ss.dat > spe.tmp" 
Write(9,*) "${GMT} psxy -M -H7 spe.tmp -B -R -JX -Sc0.15 -G0 -O >> spectrum.ps" 
!
Write(9,*) "#"
Write(9,*) "# ---------- Normalized residues vs degree, all three in figure 'n-residues.ps'"
Write(9,*) "# ---------- h"
Write(9,*) "${GMT} psbasemap  -X3 -Y4 -U/0.25/0.25/'SELEN 2.9' ", R_OPTION4, & 
	   " -Ba1f1p/f1a1p:'Normalized residue':WSen -JX7l/11l  -K > n-residues.ps" 
Write(9,*) "${GMT} psxy -M -H2 ihh.dat -B -R -JX -Ss0.25 -W4 -G0   -O -K >> n-residues.ps" 
open (4,file='tmpg7.dat',status='unknown') 
Write(4,*) "1.3 2e0 20 0 0 BL (h)"; close(4) 
Write(9,*) "${GMT} pstext -N tmpg7.dat -JX -R -G0 -CM -O -K >> n-residues.ps"
Write(9,*) "#"
open (4,file='tmpg8.dat',status='unknown') 
Write(4,*) "1 0.8e2 18 0 2 BL Earth model: NV=", trim(adjustl(NV)), " CODE=", trim(adjustl(CODE))
Write(4,*) "1 0.3e2 16 0 2 BL Viscosity profile: ", trim(adjustl(VSTRING)), "*1E21 Pa.s"; close(4) 
Write(9,*) "${GMT} pstext -N tmpg8.dat -JX -R -G0 -CM -O -K >> n-residues.ps"
Write(9,*) "# ---------- l"
Write(9,*) "${GMT} psbasemap  -X8 -R -Ba1f1p:'Harmonic degree, n':/f1a1p:'Normalized residue':wSen -JX -K -O >> n-residues.ps"
Write(9,*) "${GMT} psxy -M -H2 ill.dat -B -R -JX -Sc0.2 -W4 -G0 -O -K >> n-residues.ps"
open (4,file='tmpg9.dat',status='unknown') 
Write(4,*) "1.3 2e0 20 0 0 BL (l)"; close(4)
Write(9,*) "${GMT} pstext -N tmpg9.dat -JX -R -G0 -CM -O -K >> n-residues.ps"
Write(9,*) "#"
Write(9,*) "# ---------- k"
Write(9,*) "${GMT} psbasemap -X8 -R -Ba1f1p/f1a1p:'Normalized residue':wSen -JX -O -K >> n-residues.ps" 
Write(9,*) "${GMT} psxy -M -H2 ikk.dat -B -R -JX -Si0.25 -W4 -G0   -O -K >> n-residues.ps" 
open (4,file='tmpg10.dat',status='unknown') 
Write(4,*) "1.3 2e0 20 0 0 BL (k)"; close(4) 
Write(9,*) "${GMT} pstext -N tmpg10.dat -JX -R -G0 -CM -O >> n-residues.ps"
!
 close(9)
!
 END SUBROUTINE make_plot_ldc
!
!
!
!
!
!
 SUBROUTINE MAKE_RSLZONES (NV, CD, RUN, NINC, NRSL, TITLICE, RESOLUTION, ITER, & 
 		          MODE, DEGREE, VSTRING, FILE1_GMT, FILE2_GMT, SHORT_VISCO_FILENAME)
!
! --- A GMT script for plotting the RSL zones - GS November 12 2007 
!     Revised JULY 04 2008 for version 2.6 of SELEN
!     Revised by GS April 2010 version 3.1 ALMA & g95  
!
 IMPLICIT NONE
 INTEGER, PARAMETER :: NREGIONS=4 
 CHARACTER*20 R_OPTION(NREGIONS), G_OPTION(NREGIONS), W_OPTION(NREGIONS) 
 CHARACTER*20 FILE1_GMT, FILE2_GMT
 CHARACTER*30 VSTRING
 CHARACTER*10 DEGREE, RESOLUTION, TITLICE   
 CHARACTER*3 STRING, NV, CD
 CHARACTER*4 RUN    
 CHARACTER*2 LABEL(NREGIONS)
 CHARACTER*3 NINC
 CHARACTER*1 ITER, MODE
 CHARACTER*100 SHORT_VISCO_FILENAME 
 INTEGER I, NRSL 
!
!
! ---- A range for each  region 
 DATA LABEL  /'10','11','20','21'/
!  
! ---- A colour for each region 
 DATA G_OPTION /" -G0/0/255", " -G107/199/231", " -G226/0/122", " -G105/62/142"/ 
!
! ---- A colour for the "average RSL curve" 
 DATA W_OPTION /" -W8ta/0/0/255", " -W8ta/107/199/231 ", " -W8ta/226/0/122 ", " -W8ta/105/62/142 "/ 
! 
! ---- An R option for each region 
        R_OPTION(1) = " -R0/"//trim(adjustl(ninc))//"/-200/50"
	R_OPTION(2) = " -R0/"//trim(adjustl(ninc))//"/-200/50"
	R_OPTION(3) = " -R0/"//trim(adjustl(ninc))//"/-50/500"
	R_OPTION(4) = " -R0/"//trim(adjustl(ninc))//"/-50/150"
!
!
!
! ======== PART 1 - a single frame for each zone 
!
! --- Target GMT file for plotting <<RSL zones>> 
!
 	open(19,file=file1_gmt,status='unknown')

	Write(19,*) " ${GMT} gmtset PAPER_MEDIA A4+" 
	Write(19,*) " ${GMT} gmtset FRAME_WIDTH 0.1c"
	Write(19,*) " " 
!Write(19,*) "echo '     - ", trim(adjustl(file1_gmt))//":", " Creating ps images of RSL zones - one image per frame....'"
Write(19,*) " " 
!
! ---- Here we consider 4 zones (10, 11, and 20, 21) Some of these may be empty... this depends 
!      on ice chronology, mantle viscosity and the range of harmonic degrees. Default input files 
!      are rslzones-**.dat and lonlat-**.dat. 
!
	DO I=1, NREGIONS 
!
! --- xy plot of RSL 
!
        Write(19,*) " ${GMT} psbasemap -U'SELEN 2.9' -X3 -Y10.5 -Ba2f1:'time (ka)':/a50f50WSen:'RSL (m)': ", & 
	              r_option(i), "-JX14/9  -K >  plot.ps" 
 	Write(19,*) " ${GMT} psxy rslzones-"//label(i)//".dat", " -M -H -B -R -JX -W0.05/0 -O -K >> plot.ps"
	Write(19,*) " ${GMT} psxy rslzones-"//label(i)//"-ave.dat", " -M -H -B -R -JX ", W_OPTION(I), " -O -K >> plot.ps"	
!
! --- first sphere 
 	Write (19,*) " ${GMT} psbasemap -Y0 -X15 -Ba90/a85f90WSEN -R0/360/-90/90  -JG-20/70/9 -O -K  >> plot.ps"
 	Write (19,*) " ${GMT} psxy lonlat-"//label(i)//".dat", " -B -R -JG -Sh0.32 ", g_option(i), " -O -K >> plot.ps"
 	Write (19,*) " ${GMT} pscoast -G255 -B -R -JG -Dc -W2 -A10000 -O -K >> plot.ps"
!
! --- sechond sphere 
 	Write (19,*) " ${GMT} psbasemap -Y-10 -X0 -B -R -JG-20/-70/9 -O -K  >> plot.ps"
 	Write (19,*) " ${GMT} psxy lonlat-"//label(i)//".dat", " -B -R -JG -Sh0.32 ", g_option(i), " -O -K >> plot.ps"
 	Write (19,*) " ${GMT} pscoast -G255 -B -R -JG -Dc -W2 -A10000 -O -K >> plot.ps"
!
! --- a new basemap for titles... 
	Write (19,*) " ${GMT} psbasemap -X-15 -Y10 -Bf1000wesn   -R0/10/0/10     -JX14/9  -O -K >>  plot.ps"
!
	If(cd=='-1')then 
	   CONTINUE
	Else 
 	open (8,file='tmpz'//"."//label(i),status='unknown')	
		Write(8,*) ".5 -5  26 0 0 BL  RSL curves for zone ", label(i)
		Write(8,*) ".5 -6 16  0 0 BL  Ice model: ", trim(adjustl(titlice))
		Write(8,*) ".5 -7 16  0 0 BL  Repository: ", "./depot-/", trim(adjustl(run))
		Write(8,*) ".5 -8 16  0 0 BL  Viscosity profile: ", trim(adjustl(VSTRING))
		Write(8,*) ".5 -9 16 0 0 BL  -LMAX =", trim(adjustl(DEGREE)), & 
					       "-RES =", trim(adjustl(RESOLUTION)), & 
					       " -ITER =", trim(adjustl(ITER)), & 
					       " -Mode =", trim(adjustl(MODE)), & 
					       " -NV =", trim(adjustl(NV)), & 
					       " -CODE =", trim(adjustl(CD))			  	         
 	close(8) 	
	Endif 		
!		
	Write (19,*) " ${GMT} pstext -N ", 'tmpz'//"."//label(i), " -JX -R -O >> plot.ps"
!
	Write (19,*) " ps2pdf plot.ps"	
!
	Write (19,*) " mv plot.pdf rslzones-"//label(i)//".pdf"	
!
	Write (19,*) " mv plot.ps rslzones-"//label(i)//".ps"	
!
	Write (19,*) " "	
	
	ENDDO 
!
!
! ======== PART 2 - a frame with all zones
!
! --- Target GMT file for plotting <<RSL zones>> 
!
 	open(19,file=file2_gmt,status='unknown')

	Write(19,*) " ${GMT} gmtset PAPER_MEDIA A4+" 
	Write(19,*) " ${GMT} gmtset FRAME_WIDTH 0.1c"
	Write(19,*) " ${GMT} gmtset ANOT_FONT_SIZE 12p"
	Write(19,*) " ${GMT} gmtset LABEL_FONT_SIZE 18p"
	Write(19,*) " " 
!Write(19,*) & 
!"echo '     - ", trim(adjustl(file2_gmt))//":", " Creating ps images of RSL zones - one images in one frame....'"
Write(19,*) " " 
!
! --- xy plots of RSL zones  
!
        Write(19,*) " ${GMT} psbasemap -X3 -Y14 -Ba2f1:'time (ka)':/a50f50Wsen:'RSL (m)':  ", r_option(1), "-JX6.5/5  -K >  plot.ps" 
 	Write(19,*) " ${GMT} psxy rslzones-"//label(1)//".dat", " -M -H -B -R -JX -Sc0.01 -G0/0/255 -O -K >> plot.ps"	
	Write(19,*) " ${GMT} psxy rslzones-"//label(1)//"-ave.dat", " -M -H -B -R -JX ", W_OPTION(1), " -O -K >> plot.ps"
 	open (8,file='tmpzu1.dat',status='unknown')	
		Write(8,*) "1 -180 18 0 2 BL  zone 10"  	         
 	close(8) 
	Write (19,*) " ${GMT} pstext -N tmpzu1.dat -JX -R -O -K >> plot.ps"
!
        Write(19,*) " ${GMT} psbasemap -X8.25 -Y0 -Ba2f1:'time (ka)':/a50f50Wsen  ", r_option(2), "-JX6.5/5  -O -K >>  plot.ps" 
 	Write(19,*) " ${GMT} psxy rslzones-"//label(2)//".dat", " -M -H -B -R -JX -Sc0.01 -G107/199/231 -O -K >> plot.ps"
	Write(19,*) " ${GMT} psxy rslzones-"//label(2)//"-ave.dat", " -M -H -B -R -JX ", W_OPTION(2), " -O -K >> plot.ps"
 	open (8,file='tmpzu2.dat',status='unknown')	
		Write(8,*) "1 -180 18 0 2 BL  11"  	         
 	close(8) 
	Write (19,*) " ${GMT} pstext -N tmpzu2.dat -JX -R -O -K >> plot.ps"	
!
        Write(19,*) " ${GMT} psbasemap -U'SELEN 2.9' -X-8.25 -Y-6  -Ba2f1:'time (ka)':/a100f100WSen:'RSL (m)':  ", & 
	              r_option(3), "-JX6.5/5  -O -K >>  plot.ps" 
 	Write(19,*) " ${GMT} psxy rslzones-"//label(3)//".dat", " -M -H -B -R -JX -Sc0.01 -G226/0/122 -O -K >> plot.ps"
	Write(19,*) " ${GMT} psxy rslzones-"//label(3)//"-ave.dat", " -M -H -B -R -JX ", W_OPTION(3), " -O -K >> plot.ps"
 	open (8,file='tmpzu3.dat',status='unknown')	
		Write(8,*) "1 400 18 0 2 BL  20"  	         
 	close(8) 
	Write (19,*) " ${GMT} pstext -N tmpzu3.dat -JX -R -O -K >> plot.ps"		
!
        Write(19,*) " ${GMT} psbasemap -X8.25 -Y0    -Ba2f1:'time (ka)':/a50f50WSen  ", r_option(4), "-JX6.5/5  -O -K >>  plot.ps" 
 	Write(19,*) " ${GMT} psxy rslzones-"//label(4)//".dat", " -M -H -B -R -JX -Sc0.01 -G105/62/142  -O -K  >> plot.ps"
	Write(19,*) " ${GMT} psxy rslzones-"//label(4)//"-ave.dat", " -M -H -B -R -JX ", W_OPTION(4), " -O -K  >> plot.ps"
 	open (8,file='tmpzu4.dat',status='unknown')	
		Write(8,*) "1 115 18 0 2 BL  21"  	         
 	close(8) 
	Write (19,*) " ${GMT} pstext -N tmpzu4.dat -JX -R -O -K >> plot.ps"		
!
! --- first sphere 
 	Write (19,*) " ${GMT} psbasemap -X7.5 -Y3 -Ba90/a85f90WSEN -R0/360/-90/90  -JG-20/70/9 -O -K  >> plot.ps"
 	Write (19,*) " ${GMT} psxy lonlat-"//label(1)//".dat", " -B -R -JG -Sh0.32 ", g_option(1), " -O -K >> plot.ps"
 	Write (19,*) " ${GMT} psxy lonlat-"//label(2)//".dat", " -B -R -JG -Sh0.32 ", g_option(2), " -O -K >> plot.ps"
 	Write (19,*) " ${GMT} psxy lonlat-"//label(3)//".dat", " -B -R -JG -Sh0.32 ", g_option(3), " -O -K >> plot.ps"
 	Write (19,*) " ${GMT} psxy lonlat-"//label(4)//".dat", " -B -R -JG -Sh0.32 ", g_option(4), " -O -K >> plot.ps"
 	Write (19,*) " ${GMT} pscoast -G255 -B -R -JG -Dc -W2 -A10000 -O -K >> plot.ps"
!
! --- sechond sphere 
 	Write (19,*) " ${GMT} psbasemap -Y-10 -X0 -B -R -JG-20/-70/9 -O -K  >> plot.ps"
 	Write (19,*) " ${GMT} psxy lonlat-"//label(1)//".dat", " -B -R -JG -Sh0.32 ", g_option(1), " -O -K >> plot.ps"
 	Write (19,*) " ${GMT} psxy lonlat-"//label(2)//".dat", " -B -R -JG -Sh0.32 ", g_option(2), " -O -K >> plot.ps"
 	Write (19,*) " ${GMT} psxy lonlat-"//label(3)//".dat", " -B -R -JG -Sh0.32 ", g_option(3), " -O -K >> plot.ps"
 	Write (19,*) " ${GMT} psxy lonlat-"//label(4)//".dat", " -B -R -JG -Sh0.32 ", g_option(4), " -O -K >> plot.ps"
 	Write (19,*) " ${GMT} pscoast -G255 -B -R -JG -Dc -W2 -A10000 -O -K >> plot.ps"
!
! --- a new basemap for titles... 
	If(cd=='-1')then 
	   CONTINUE
	Else 
	Write (19,*) " ${GMT} psbasemap -X-7.5 -Y7 -Bf1000wesn   -R0/10/0/10     -JX6.5/5  -O -K >>  plot.ps"
 	open (8,file='tmpzu.dat',status='unknown')	

		Write(8,*) "-9 -7  18 0  0 BL  RSL zones "
		Write(8,*) "-9 -9  14  0 0 BL  Ice model: ", trim(adjustl(titlice))
		Write(8,*) "-9 -10.5 14  0 0 BL  Repository: ./depot-", trim(adjustl(run)) 
		Write(8,*) "-9 -12.0 14  0 0 BL  Viscosity profile: ", trim(adjustl(VSTRING))
		Write(8,*) "-9 -13.5 14  0 0 BL  -LMAX =", trim(adjustl(DEGREE)), & 
					       " -RES =", trim(adjustl(RESOLUTION)), & 
					       " -ITER =", trim(adjustl(ITER)), & 
					       " -Mode =", trim(adjustl(MODE)), & 
					       " -NV =", trim(adjustl(NV)), & 
					       " -CODE =", trim(adjustl(CD))			  	         
 	close(8) 	
	Endif 		
!
!
	Write (19,*) " ${GMT} pstext -N tmpzu.dat -JX -R -O >> plot.ps"
!	
	Write (19,*) " ps2pdf plot.ps"	
!
	Write (19,*) " mv plot.ps rslzones-all.ps"	
!
	Write (19,*) " mv plot.pdf rslzones-all.pdf"	
!
	close(19) 
!	
 END SUBROUTINE MAKE_RSLZONES
!
!
!
!
!
!
 SUBROUTINE MAKE_RSLMIS (NV, CD, RUN, NINC, NRSL, TITLICE, RESOLUTION, ITER, & 
 		      MODE, DEGREE, VSTRING, FILE_GMT, SHORT_VISCO_FILENAME)
!
! --- Prepares a GMT script for plotting a simple histogram
!     showing the RSL misfit distribution - GS January 23 2008 
!     Revised July 2008 for version 2.6 of SELEN
!     Revised April 2010 GS for version g95 
!
 IMPLICIT NONE
 CHARACTER*30 VSTRING
 CHARACTER*20 FILE_GMT
 CHARACTER*10 DEGREE, RESOLUTION, TITLICE   
 CHARACTER*4 RUN  
 CHARACTER*3 NV, CD 
 CHARACTER*3 NINC
 CHARACTER*1 ITER, MODE
 CHARACTER*100 SHORT_VISCO_FILENAME
 INTEGER NRSL 
!
!
! --- Target GMT file for misfit histogram   
!
 	open(19,file=file_gmt,status='unknown')

	Write(19,*) "${GMT} gmtset PAPER_MEDIA A4+" 
	Write(19,*) "${GMT} gmtset FRAME_WIDTH 0.1c"
!	Write(19,*) "echo       - Creating a postscript image of the misfit histogram..."	
!
  	Write(19,*) "${GMT} psbasemap -X3 -Y3 -U'SELEN 2.9' -JX16/10  -G255 -R0/100/0/50 -K > rsl-misfit.ps" 
	Write(19,*) "awk '{print $1, $2}' mis.dat > mis.tmp"
   	Write(19,*) "pshistogram mis.tmp -JX -R -B100a20:misfit:/a10:frequency:WSen -G140 -W2 -O -K >> rsl-misfit.ps" 
	Write(19,*) "/bin/rm mis.tmp"
!
 	open (8,file='tmph.dat',status='unknown')		  	         
	If(cd=='-1')then 
	   CONTINUE
	Else 
!
	Write(8,*) "45 40 12 0 0 BL  Ice model: ", trim(adjustl(titlice))
	Write(8,*) "45 37 12 0 0 BL  Repository: ./depot-", trim(adjustl(run))
	Write(8,*) "45 34 12 0 0 BL  Viscosity profile: ", trim(adjustl(VSTRING))
	Write(8,*) "45 31 10 0 0 BL  -LMAX =", trim(adjustl(DEGREE)), & 
				   " -RES =",  trim(adjustl(RESOLUTION)), & 
				   " -ITER =", trim(adjustl(ITER)), & 
				   " -Mode =", trim(adjustl(MODE))					  
	Write(8,*) "45 28 10 0 0 BL    -NV =", trim(adjustl(NV)),     & 
		       	           " -CODE =", trim(adjustl(CD))

	Write(8,*) "45 25 10 0 2 BL  Number of sites =", NRSL	
 	close(8) 
!	
	Endif
!
	Write(19,*) "${GMT} pstext tmph.dat -O -JX -R >> rsl-misfit.ps"	
	close(19) 
!	
 END SUBROUTINE MAKE_RSLMIS
!
!
!
!
!
!
 SUBROUTINE MAKE_RSLDB (DATAFILE, RSL_FORMAT, NRSL, FILE_GMT)
 IMPLICIT NONE
!
! Writes a GMT script for plotting a map of RSL sites distribution 
! The output filename is "File_GMT"   --- GS October 27 2007 ---
!
! Also modified in August 2009 for the implementation of SELEN 3 
!
 CHARACTER*1   CJUNK, RSL_FORMAT 
 CHARACTER*100 SS(100)
 CHARACTER*10  LATSC10, LONSC10 
 CHARACTER*50  LINE
 CHARACTER*20  FILE_GMT
 CHARACTER*30  DATAFILE 
 INTEGER I, J, NRSL, NOUT, CODE, NDATA
 REAL*8 LONS, LATS
 CHARACTER*200 LINEP 
!
!
! --- Open and reads the RSL database for getting info about lon/lat of data... 
!
! --- Target file for lon-lat of RSL sites 
 	open(37,file='lon-lat-rsl.dat',status='unknown')
!
! --- An header 
	Write(37,*) "Longitudes and latitudes of RSL sites (degrees)"
!
	If(RSL_FORMAT=='0') then 
!
        OPEN(1,FILE=DATAFILE,STATUS='unknown')
   	do i=1, nrsl 	
		READ (1, 6) code, lats, lons, ndata
		if(lons<=0.) lons=lons+360. 			 			 		 			 
		write(37,'(F6.1, 1X, F6.1)') lons, lats
		do j=1, ndata
			read(1,*) LINE			
		enddo
	enddo
	CLOSE(1) 
	CLOSE(37) 
6       FORMAT(1X,I3,1X,F5.1,1X,F6.1,1X,I2,1X,A22)
!
	ELSEIF(RSL_FORMAT=='1') then
!
        OPEN(1,FILE=DATAFILE,STATUS='unknown')
!	
   	DO I=1, NRSL 
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 
			        READ(1,*) LATS, LONS  
                        	IF(LONS<=0.) LONS=LONS+360. 	
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 
 			WRITE(37,'(F9.4, 1X, F9.4)') LONS, LATS
	ENDDO 
	CLOSE(1) 
	CLOSE(37) 
!
	ELSEIF(RSL_FORMAT=='2') then
!
        OPEN(1,FILE=DATAFILE,STATUS='unknown')
!	
   	DO I=1, NRSL 
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 			
			read(1,'(a200)') linep 	
!		
			call scan_string (linep, 2, ss, nout)				
!
			LATSC10=trim(adjustl(ss(1))) 
			LONSC10=trim(adjustl(ss(2)))			
!
 			call FROMDDMMSS_2_DEGREES (LATSC10, LATS) 
 			call FROMDDMMSS_2_DEGREES (LONSC10, LONS) 
!
                        IF(LONS<=0.) LONS=LONS+360. 	
			READ(1,'(A1)')CJUNK 
			READ(1,'(A1)')CJUNK 
 			WRITE(37,'(F9.4, 1X, F9.4)') LONS, LATS
	ENDDO 
	CLOSE(1) 
	CLOSE(37) 



        Endif 
!
!
! -------------------------------------------------------------
! --- Creates a GMT script for plotting a map of RSL sites  --- 
! -------------------------------------------------------------
!
! --- Target GMT file for a map of RSL sites 
!
 	open(29,file=file_gmt,status='unknown')
!
	Write(29,*) "${GMT} gmtset PAPER_MEDIA A4+" 
	Write(29,*) "${GMT} gmtset HEADER_FONT_SIZE 24p"
	Write(29,*) "${GMT} gmtset FRAME_WIDTH 0.1c"
        Write(29,*) "${GMT} gmtset ANOT_FONT_SIZE 12p"
	Write(29,*) ""
	Write(29,*) ""

	open(22,file='tmptitle',status='unknown') 
	Write(22,*) "180 87.4 22 0 1 BC Relative Sea Level (RSL) sites"
	Write(22,'(a29,a12,a7,i4)') "180 86.2 16 0 2 BC from file ", trim(adjustl(datafile)), " -NRSL=", NRSL 			
	close(22) 
!
! ------ Main frame (Mercator projection) 
!
	Write(29,*) "${GMT} psbasemap -X4 -Y5 -Ba180/a85f90WSEn -R0/360/-85/85  -JM12 -K > maprsl.ps" 	
	Write(29,*) "${GMT} pscoast -G120 -S0/0/220 -B -R -O -K  -JM -Dc  -A10000 >> maprsl.ps" 	
	Write(29,*) "${GMT} psxy -H1 lon-lat-rsl.dat -B -R -JM -Sc0.1 -G220 -O -K >> maprsl.ps"  
	Write(29,*) "${GMT} pstext -N tmptitle", " -G0 -O -K -JM -R >> maprsl.ps"
	Write(29,*) ""
!
! ------ Northern emisphere 
!	     	
 	Write(29,*) "${GMT} psbasemap -Y6 -X14 -Ba90/a85f90WSEN -R0/360/0/90  -JG-45/90/8 -O -K  >> maprsl.ps"  
	Write(29,*) "${GMT} pscoast -G120 -S0/0/220 -B -R -O -K -JG -Dc  -A10000 >> maprsl.ps" 	
	Write(29,*) "${GMT} psxy -H1 lon-lat-rsl.dat -B -R -JG -Sc0.1 -G220 -O -K >> maprsl.ps"  
	Write(29,*) ""
!
! ------Southern emisphere 
!	     	
 	Write(29,*) "${GMT} psbasemap -Y-9 -X0 -Ba90/a85f90WSEN -R0/360/-90/0  -JG-90/-90/8 -U/-10/1/'SELEN 2.9' -O -K  >> maprsl.ps"  
	Write(29,*) "${GMT} pscoast -G120 -S0/0/220 -B -R -O -K -JG -Dc  -A10000 >> maprsl.ps" 	
	Write(29,*) "${GMT} psxy -H1 lon-lat-rsl.dat -B -R -JG -Sc0.1 -G220 -O >> maprsl.ps"  
	Write(29,*) ""
!
	close(29) 	
!	
 END SUBROUTINE MAKE_RSLDB
!
!
!
!
!
!
 SUBROUTINE MAKE_RSL (NV, CD, DATAFILE, RSL_FORMAT, RUN, NINC, NRSL, TITLICE, & 
 	              RESOLUTION, ITER, MODE, DEGREE, VSTRING, FILE_GMT,      & 
		      SHORT_VISCO_FILENAME)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Writes a GMT script for plotting the observed RSL and predictions 
! --- GS February 21, 2008 ---
! Revised July 2008 for v 2.6 
! Revised August 2008 for v 2.7 
! Also revised August 2009 for v 3
! Updated on April 2010 by GS for v. 3.1 (g95)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
 IMPLICIT NONE
 CHARACTER*20 R_OPTION
 CHARACTER*50 B_OPTION
 CHARACTER*30 VSTRING
 CHARACTER*20 FILE_GMT, YMINC, YMAXC, DELTAC
 CHARACTER*30 DATAFILE
 CHARACTER*10 RESOLUTION, DEGREE, TITLICE 
 CHARACTER*4 RUN, STRING, CJUNK     
 CHARACTER*3 NV, CD
 CHARACTER*3 NINC
 CHARACTER*1 ITER, MODE, RSL_FORMAT 
 INTEGER I, J, K, NRSL 
 INTEGER, PARAMETER :: MAXSITES=1000, MAXDATA=1000
 CHARACTER*50 TITRE(MAXSITES), TITREA(MAXSITES), TITREB(MAXSITES)
 INTEGER CODE(MAXSITES), NDATA(MAXSITES)
 REAL*8    TIME(MAXSITES,0:MAXDATA)    	! times bp for which data are available at site#k
 REAL*8   DTIME(MAXSITES,0:MAXDATA)    	! uncertainties on the above times...
 REAL*8     RSL(MAXSITES,0:MAXDATA)    	! rsl datum for site k ad the times above
 REAL*8    DRSL(MAXSITES,0:MAXDATA)    	! uncertainties on the rsl datum
 REAL*8 AJUNK, LONS(MAXSITES), LATS(MAXSITES)
 REAL*8 DMIN, DMAX, RANGE, DELTA 
 CHARACTER*3 CODEC(MAXSITES)
 CHARACTER*10 LATSC10, LONSC10 
 CHARACTER*100 SS(2)
 CHARACTER*100 SHORT_VISCO_FILENAME  
 CHARACTER*200 LINEP 
 INTEGER NOUT
!
!
!
! ------ Open and reads the sealevel database for getting info about data... 
!

	If    (RSL_FORMAT=='0') then 
!  
        do k=1, 2 
 		OPEN(1, FILE=DATAFILE,STATUS='unknown')
   	do i=1, nrsl 	
 		if(k==1) then 
			 READ (1, 6) code(i), lats(i), lons(i), ndata(i)
6  			 FORMAT(1X,I3,1X,F5.1,1X,F6.1,1X,I2,1X,A22)
			 if(lons(i)<=0.) lons(i)=lons(i)+360.
!			 
			 open(33,file='junk.dat',status='unknown'); write(33,'(i3)')  code(i); close(33) 
			 open(33,file='junk.dat',status='unknown'); read (33,'(a3)') codec(i); close(33) 			 			 
!			 			 
			 else
			 READ (1, '(A50)') TITRE(I)
			 endif 
!			 
		do j=1, ndata(i)
			read(1,*) time(i,j), dtime(i,j), rsl(i,j), drsl(i,j)
		enddo
	enddo
!
	CLOSE(1) 
	enddo
!
	Elseif(RSL_FORMAT=='1') then
!
 		OPEN(1, FILE=DATAFILE,STATUS='unknown')
!
   		DO I = 1, NRSL 
		    NDATA(I)=1
		    READ(1,*)CODE(I)
			 open(33,file='junk.dat',status='unknown'); write(33,'(i3)')  code(i); close(33) 
			 open(33,file='junk.dat',status='unknown'); read (33,'(a3)') codec(i); close(33)			
!
	            READ(1,'(A50)')TITREA(I) 
		    READ(1,'(A50)')TITREB(I)
		         TITRE(I)=TRIM(ADJUSTL(TITREA(I)))//" - "//TRIM(ADJUSTL(TITREB(I))) 
	            READ(1,*) LATS(I), LONS(I)  
                         IF(LONS(I)<=0.) LONS(I)=LONS(I)+360. 					
		    READ(1,*) TIME(I,1), DTIME(I,1), RSL(I,1), DRSL(I,1)
		    READ(1,'(A1)')CJUNK 
		ENDDO 
		CLOSE(1) 
!
	Elseif(RSL_FORMAT=='2') then
!
 		OPEN(1, FILE=DATAFILE,STATUS='unknown')
!
   		DO I = 1, NRSL 
		    NDATA(I)=1
		    READ(1,*)CODE(I)
			 open(33,file='junk.dat',status='unknown'); write(33,'(i3)')  code(i); close(33) 
			 open(33,file='junk.dat',status='unknown'); read (33,'(a3)') codec(i); close(33)			
!
	            READ(1,'(A50)')TITREA(I) 
		    READ(1,'(A50)')TITREB(I)
		         TITRE(I)=TRIM(ADJUSTL(TITREA(I)))//" - "//TRIM(ADJUSTL(TITREB(I))) 
!
			READ(1,'(A200)') LINEP 			
			CALL SCAN_STRING (LINEP, 2, SS, NOUT)
!	
			LATSC10=trim(adjustl(SS(1))) 
			LONSC10=trim(adjustl(SS(2)))				
! 
 			call FROMDDMMSS_2_DEGREES (LATSC10, LATS(I)) 
 			call FROMDDMMSS_2_DEGREES (LONSC10, LONS(I)) 
!
                         IF(LONS(I)<=0.) LONS(I)=LONS(I)+360. 					
		    READ(1,*) TIME(I,1), DTIME(I,1), RSL(I,1), DRSL(I,1)
		    READ(1,'(A1)')CJUNK 
		ENDDO 
		CLOSE(1) 
!
	ENDIF
!
!
!
! ------ Part 1: Creates a GMT script for RSL curves  
!
! ------ Target GMT file for RSL curves  
!
 	open(19,file=file_gmt,status='unknown')

	Write(19,*) "${GMT} gmtset PAPER_MEDIA A4+" 
	Write(19,*) "${GMT} gmtset FRAME_WIDTH 0.1c"
!
!Write(19,*) & 
!"echo '--- ", trim(adjustl(file_gmt))//":", " Creating ps images of RSL data and predictions'"
!
	DO 11111 I=1, NRSL
!
! ------- Based on the sealevel database, finds a range for plotting
!
	dmax=-99999.
	dmin=+99999.
	do j=1, ndata(i) 
		if(rsl(i,j)<=dmin) dmin=rsl(i,j) 
		if(rsl(i,j)>=dmax) dmax=rsl(i,j)  	
	enddo
	if(dmin<=0.and.dmax<=0) dmax=0 
	if(dmin>=0.and.dmax>=0) dmin=0 
!	 	
	range=dmax-dmin 
!	
	if(range>=0  .and.range<=10)  delta=1
	if(range>=10 .and.range<=50)  delta=5	  
	if(range>=50 .and.range<100) delta=10	
	if(range>=100.and.range<200) delta=20
	if(range>=200.and.range<300) delta=50
	if(range>=300.and.range<400) delta=100
	if(range>=400) 		  delta=150 
!
   	OPEN  (10,FILE='junky'//trim(adjustl(codec(i))),STATUS='unknown')
	WRITE (10,'(i10)') int(dmin-delta)
	WRITE (10,'(i10)') int(dmax+delta)
	WRITE (10,'(i10)') int(delta) 	
	CLOSE(10)
   	OPEN  (10,FILE='junky'//trim(adjustl(codec(i))),STATUS='unknown')
	READ (10,'(a10)') yminc
	READ (10,'(a10)') ymaxc
	READ (10,'(a10)') deltac 
	close(10)		
!
! ------ GMT options 
!
	r_option = & 
	"-R0/"//trim(adjustl(ninc))//"/"//trim(adjustl(yminc))//"/"//trim(adjustl(ymaxc))
	b_option = & 
	"-Ba2f1:'time (ka)':/a"//trim(adjustl(deltac))//"f"//trim(adjustl(deltac))//"WSen:'RSL (m)':"
!	
!If(mod(i,25)==0)Write(19,*) "echo  - RSL Code ", CODEC(i)
!	
	Write(19,*) "${GMT} psbasemap -U'SELEN 2.9' -X4 -Y14 ", b_option, r_option," -P -JX12.8/8  -K >  plot.ps" 
	Write(19,*) "${GMT} psxy -H2 ", "rslp"//"-"//trim(adjustl(CODEC(i)))//".dat", " -B -R -JX -W4 -O -K >> plot.ps" 
	Write(19,*) "${GMT} psxy -H2 ", "rsld"//"-"//trim(adjustl(CODEC(i)))//".dat", " -B -R -JX -Ss0.25 -G0 -O -K >> plot.ps" 
	Write(19,*) "${GMT} psxy -H2 ", "rsld"//"-"//trim(adjustl(CODEC(i)))//".dat", " -Ey0.25/4 -B -R -JX -O -K >> plot.ps" 
!
! --- A new basemap for titles... 
	Write(19,*) "${GMT} psbasemap -X0 -Y0 -Ba1000wesn -R0/18/0/10 -P -JX12.8/8 -O -K >>  plot.ps" 
!
	If(cd=='-1') then 
	   CONTINUE	
	Else 
 	open (8,file='tmpb'//"."//trim(adjustl(CODEC(i))),status='unknown')
	Write(8,*) "8 9.0 9 0 0 BL  ", titre(i)			  	         
	Write(8,*) "8 8.5 9 0 0 BL  Ice model: ", trim(adjustl(titlice))
	Write(8,*) "8 8.0 9 0 0 BL  Repository label: ", trim(adjustl(run))
	Write(8,*) "8 7.5 9 0 0 BL  Viscosity profile: ", trim(adjustl(VSTRING))
	Write(8,*) "8 7.0 8 0 0 BL     -LMAX =", trim(adjustl(DEGREE)),     & 
		       	               "-RES =", trim(adjustl(RESOLUTION)), & 
				     " -ITER =", trim(adjustl(ITER)),       & 
				     " -Mode =", trim(adjustl(MODE))	  
	Write(8,*) "8 6.5 8 0 0 BL       -NV =", trim(adjustl(NV)),     & 
		       	             " -CODE =", trim(adjustl(CD))				     	         
 	close(8) 		
	Endif
!		
	Write(19,*) "${GMT} pstext ", 'tmpb'//"."//trim(adjustl(CODEC(i))), " -O -JX -R >> plot.ps"	
!
	Write(19,*) "ps2pdf plot.ps"
!
	Write(19,*) "mv plot.pdf ", "rslp-"//trim(adjustl(CODEC(i)))//".pdf"
!
	Write(19,*) "mv plot.ps ", "rslp-"//trim(adjustl(CODEC(i)))//".ps"
!
	close(11) 
!
11111  CONTINUE
!
!Write(19,*) "echo       - RSL curves reported on files rslp-yyy.ps and rslp-yyy.pdf"		
	close(19)	
!	
 END SUBROUTINE MAKE_RSL
!
!
!
!
!
!
 SUBROUTINE MAKE_ICEMAP(NINC, TITLICE, OPTION_ROF, FILECAP, FILE)
 IMPLICIT NONE
 CHARACTER*1  OPTION_ROF
 CHARACTER*2  LABCHAR, CHARJUNK
 CHARACTER*3 NINC 
 CHARACTER*10 TITLICE  
 CHARACTER*20 FILE, FILECAP, FILENAME
 CHARACTER*20 C_TABLE, C_OPTION, B_OPTION 
 INTEGER I, N, NN
!
! --- Creates a GMT script for plotting the original ice sheets
!     distributions using various projections - 
!     					    GS November 14 2007
!					    Revised February 21, 2008- 
!					    Revised June 14, 2008. 
!					    Revised (IJ05) July 5, 2008 
!					    Also revised July 14, 2008 
!					    Revised July 26, 2008 
!
!
 open(17,file=file,status='unknown')
!
!
! # Some settings ... 
!
 Write(17,*) "${GMT} gmtset PAPER_MEDIA A4+"
 Write(17,*) "${GMT} gmtset HEADER_FONT_SIZE 24p" 
 Write(17,*) "${GMT} gmtset FRAME_WIDTH 0.1c" 
!
!
! ------ 
 C_TABLE="pice.cpt"
 C_OPTION=" -C"//trim(adjustl(C_TABLE))
 If(titlice=="IJ05MOD")then 
   Write(17,*) "${GMT} makecpt -Crainbow  -T-100/1000/100 -D > ", trim(adjustl(C_TABLE)) 
   B_OPTION = " -Bf1000a100g100/:: "			      
 endif 
 If(titlice/="IJ05MOD")then 
   Write(17,*) "${GMT} makecpt -Crainbow  -T-500/4500/500 -D > ", trim(adjustl(C_TABLE)) 
   B_OPTION = " -Bf4500a1000g500/:: "
 endif
!
 Write(17,*) ""
!
 open(10,file='junk.dat',status='unknown'); Write(10,'(a3)') NINC; close(10) 
 open(10,file='junk.dat',status='unknown'); Read(10,*) NN; close(10) 
!
! Write(17,*) "echo '--- ", trim(adjustl(file))//":", " Creating ps images of original ice sheets distributions'"
!
 do i=0, nn+1 
! 
 Write(17,*) ""
 If (i==0)  then 
 	Write(17,*) "echo '     - until'", nn, "ka "
!elseif (1<=i.and.i<=nn) then
!Write(17,*) "echo '     - Between '", nn-i+1, " and ", nn-i, "ka "	    
	    elseif (i.eq.nn+1) 	        then
 	Write(17,*) "echo '     - today'"
 Endif 
!
	open(3,file='junk.dat',status='unknown') 
	if(i<=9) write(3,'(a1,i1)') '0',i  
	if(i> 9) write(3,'(i2)')        i
 	close(3)
	open(3,file='junk.dat',status='unknown'); read(3,'(a2)') labchar; close(3)
!
 	IF(titlice=="ICE5G26")  filename='msg5-'//labchar//'.dat'
 	IF(titlice=="ICE5G")    filename='msg5-'//labchar//'.dat'
 	IF(titlice=="ICE3G")    filename='msg3-'//labchar//'.dat'
 	IF(titlice=="DISK")     filename='msgD-'//labchar//'.dat'
 	IF(titlice=="ICAP")     filename='msgC-'//labchar//'.dat'
 	IF(titlice=="IJ05MOD")  filename='msgI-'//labchar//'.dat'
 	IF(titlice=="ANU05")    filename='msgU-'//labchar//'.dat'
 	IF(titlice=="ICE1")     filename='msg1-'//labchar//'.dat'
 	IF(titlice=="ALPS")     filename='msgA-'//labchar//'.dat'
        IF(titlice=="ALASKA")   filename='msgK-'//labchar//'.dat'
!
!
! ------ Main frame (Mercator projection) 
!
	Write(17,*) "${GMT} psbasemap -X4 -Y6 -Ba180/a85f90WSEn -R0/360/-85/85  -JM12 -U/4/13.5/'SELEN 2.9' -K > map.ps" 	
	Write(17,*) "${GMT} psxy ", trim(adjustl(C_OPTION)), " -A -L -M ", trim(adjustl(filename)), " -O -K -JM -R >> map.ps" 
	If(option_rof=='r')then
		           Write(17,*) "${GMT} pscoast -R -JM -Dc  -B -W1/0   -A10000 -O -K >> map.ps" 
		           else
		           Write(17,*) "${GMT} pscoast -R -JM -Dc  -B -W1/240 -A10000 -O -K >> map.ps" 
	Endif			   
	Write(17,*) "${GMT} psscale -E ", trim(adjustl(C_OPTION)), " ", trim(adjustl(B_OPTION)), " -D6/-1/11/1h -O  -K >> map.ps"
!
! ------ Prepares a small title 
!
 	open(8,file='rtmp0'//labchar//'.dat',status='unknown')
	if(i==0) then 
	Write(8,'(a18,1x,a5,a18,i2,1x,a8)')    "180 86.5 18 0 1 BC ", & 
	      				     trim(adjustl(titlice)), " thickness until ", nn, " ka "
 	      	 elseif (1<=i.and.i<=nn)     then	
	Write(8,'(a18,1x,a5,a11,i2,a1,i2,a8)') "180 86.5 18 0 1 BC ", & 
		   	                     trim(adjustl(titlice)), " thickness ", nn-i+1,"-",nn-i, " ka "
	         elseif (i==nn+1) 	     then
	Write(8,'(a18,1x,a5,a16)')             "180 86.5 18 0 1 BC ", & 
				             trim(adjustl(titlice)), " thickness today"
	endif
	Write(8,*) "180 -89.2 14 0 2 BC ", adjustl(trim(titlice)), " thickness (m)"		
	close(8) 
	Write(17,*) "${GMT} pstext -N ", "rtmp0"//labchar//".dat", " -O -K -JM -R -G0 >> map.ps"

        if(option_rof=='z')& 
	Write(17,*) "${GMT} psxy -R -JM ", trim(adjustl(filecap)), " -M -W1/0 -B -A -O >> map.ps"
!
!
	if(option_rof=='r') then 
!
! ------ Northern emisphere 
!	     	
 	Write(17,*) "${GMT} psbasemap -Y6 -X14.5 -Ba90/a85f90WSEN -R0/360/0/90 -JG-45/90/8 -O -K >> map.ps"  
	Write(17,*) "${GMT} psxy ", trim(adjustl(C_OPTION)), " -A -L -M ", trim(adjustl(filename)), " -O  -K -JG -R >> map.ps"
	If(option_rof=='r')then
		           Write(17,*) "${GMT} pscoast -R -JG -Dc  -B -W1/0   -A10000 -O -K >> map.ps" 
		           else
		           Write(17,*) "${GMT} pscoast -R -JG -Dc  -B -W1/240 -A10000 -O -K >> map.ps" 
	Endif	
!
! ------ Southern emisphere 
!
 	Write(17,*) "${GMT} psbasemap -Y-9 -X0 -Ba90/a85f90WSEN -R0/360/-90/0  -JG-90/-90/8 -O -K  >> map.ps"  
	Write(17,*) "${GMT} psxy ", trim(adjustl(C_OPTION)), " -A -L -M ", trim(adjustl(filename)), " -O  -K -JG -R >> map.ps"
	If(option_rof=='r')then
		           Write(17,*) "${GMT} pscoast -R -JG -Dc  -B -W1/0   -A10000 -O >> map.ps" 
		           else
		           Write(17,*) "${GMT} pscoast -R -JG -Dc  -B -W1/240 -A10000 -O >> map.ps" 
	Endif	
!
	Endif
!	
! ------ A pdf image and a name for each frame  
!
	Write(17,*) "ps2pdf map.ps"
	Write(17,*) "mv map.pdf ", "mapice"//labchar//".pdf"
	Write(17,*) "mv map.ps ", "mapice"//labchar//".ps"
!
enddo
!
 close(17)
!
 END SUBROUTINE MAKE_ICEMAP
!
!
!
!
!
!
 SUBROUTINE MAKE_RECICEMAP(NINC, DEGREE, TITLICE, OPTION_ROF, FILECAP, FILE)
 IMPLICIT NONE
 CHARACTER*1 OPTION_ROF
 CHARACTER*2 LABCHAR, CHARJUNK
 CHARACTER*3 NINC  
 CHARACTER*10 DEGREE, TITLICE   
 CHARACTER*20 FILE, FILECAP, B_OPTION, C_TABLE, C_OPTION 
 CHARACTER*12 FILENAME
 INTEGER I, N, NN
!
!
! --- Creates a GMT script for plotting the rechonstructed ice sheets
!     distributions using various projections - GS November 07
! 
!     Revised July 05, 2008 - For INTEL version (2.6)  
!     Also revised July 14, 2008 
!     Revised July 26, 2008 
!
!
 open(1,file=file,status='unknown')
!
!
! # Some settings ... 
!
 Write(1,*) "${GMT} gmtset PAPER_MEDIA A4+"
 Write(1,*) "${GMT} gmtset HEADER_FONT_SIZE 24p" 
 Write(1,*) "${GMT} gmtset FRAME_WIDTH 0.1c" 
!
 Write(1,*) ""
!
!
 open(10,file='junk.dat',status='unknown'); Write(10,'(a3)') NINC; close(10) 
 open(10,file='junk.dat',status='unknown'); Read(10,*) NN; close(10) 
!
! Write(1,*) "echo '--- ", trim(adjustl(file))//":", " Creating ps images of rechonstructed ice sheets'"
!
 do i=0, nn+1 
! 
 Write(1,*) ""
 If (i==0)  then 
 	Write(1,*) "echo '     - until'", nn, "ka "
!elseif (1<=i.and.i<=nn) then
!Write(1,*) "echo '     - Between '", nn-i+1, " and ", nn-i, "ka "	    
	    elseif (i.eq.nn+1) 	        then
 	Write(1,*) "echo '     - today'"
 Endif 
!
	open(3,file='junk.dat',status='unknown') 
	if(i<=9) write(3,'(a1,i1)') '0',i  
	if(i> 9) write(3,'(i2)')        i
 	close(3)
	open(3,file='junk.dat',status='unknown'); read(3,'(a2)') labchar; close(3)
!
	filename='rect'//labchar//'.dat'
!
!
! ------ [NEW] Main frame (Mercator projection + PScontour) 
 C_TABLE="pice.cpt"
 C_OPTION=" -C"//trim(adjustl(C_TABLE))
 If(titlice=="IJ05MOD")then 
   Write(1,*) "${GMT} makecpt -Crainbow  -T-100/1000/100 -D > ", trim(adjustl(C_TABLE)) 
   B_OPTION = " -Bf1000a100g100/:: "			      
 endif 
 If(titlice/="IJ05MOD")then 
   Write(1,*) "${GMT} makecpt -Crainbow  -T-500/4500/500 -D > ", trim(adjustl(C_TABLE))  
   B_OPTION = " -Bf4500a1000g500/:: "
 endif
!
        Write(1,*) "${GMT} pscontour -X3 -Y5  -I -JM12 -Ba180/a85f90WSEn -R0/360/-85/85 ", filename, " ", & 
	            trim(adjustl(C_OPTION)), " -K > map.ps"
!
	If(option_rof=='r')then
		   Write(1,*) "${GMT} pscoast -U/4/13.5/'SELEN 2.9' -B -R -O -K -W2/255 -JM -Dc -A10000 >> map.ps" 
		           else
		   Write(1,*) "${GMT} pscoast -U/4/13.5/'SELEN 2.9' -B -R -O -K -W2/0/0/255 -JM -Dc -A10000 >> map.ps" 		   
	Endif
!
        if(option_rof=='z')& 
	Write(1,*) "${GMT} psxy -R -JM ", trim(adjustl(filecap)), " -M -W1/0 -B -A -O -K >> map.ps"
!        
        Write(1,*) "${GMT} psscale -E ", trim(adjustl(C_OPTION)), " ", trim(adjustl(B_OPTION)), & 
	           " -D6/-1/11/1h -O -K >> map.ps"
        Write(1,*) "${GMT} pstext -N ", " tmp0"//labchar//".dat",  " -JM -R -G0 -O -K >> map.ps"
!
! ------ Prepares a small title 
!
 	open(8,file='tmp0'//labchar//'.dat',status='unknown')
	if(i==0) then 
	Write(8,'(a18,1x,a7,a21,a4,a3,a7,i2,a9)')    "180 86.5 18 0 1 BC  ", trim(adjustl(titlice)), & 
		                                  " thickness to degree ", trim(adjustl(degree)), & 
					          " / ", " until ", nn, " ka  "
 	      	 elseif (1<=i.and.i<=nn) then	
	Write(8,'(a18,1x,a7,a21,a4,a3,i2,a1,i2,a9)') "180 86.5 18 0 1 BC  ", trim(adjustl(titlice)), & 
	        			          " thickness to degree ", trim(adjustl(degree)), & 
					          " / ", nn-i+1,"-",nn-i, " ka  "
	         elseif (i.eq.nn+1) then
	Write(8,'(a20,1x,a7,a21,a4,a3,a6)')          "180 86.5 18 0 1 BC  ", trim(adjustl(titlice)), & 
	        				  " thickness to degree ", trim(adjustl(degree)), & 
						  " / ", " today"
	endif	
	Write(8,*) "180 -89.2 14 0 2 BC ", trim(adjustl(titlice)), " thickness (m)"		
	close(8) 
!
!--- Southern emisphere 
!
 	Write(1,*) "${GMT} pscontour -I -Y1 -X14.5 -JA0/-90/10 -R0/360/-90/-50 ", filename, " ", & 
	            trim(adjustl(C_OPTION)), " -O -K >> map.ps"
! 	
	If(option_rof=='r')then
		   Write(1,*) "${GMT} pscoast -R -JA -Di -Ba45f45 -W2/255     -A0 -O >> map.ps"
		           else
		   Write(1,*) "${GMT} pscoast -R -JA -Di -Ba45f45 -W2/0/0/255 -A0 -O >> map.ps"		   
	Endif
!
!
! ------ A pdf image and a name for each frame  
!
	Write(1,*) "ps2pdf map.ps"
	Write(1,*) "mv map.pdf ", "recice"//labchar//".pdf"
	Write(1,*) "mv map.ps ",  "recice"//labchar//".ps"
!
enddo
!
 close(1)
!
 END SUBROUTINE MAKE_RECICEMAP
!
!
!
!
!
!
 SUBROUTINE MAKE_OFMAP(DEG, OPT, AMP, FILECAP, FILE)
 IMPLICIT NONE
 CHARACTER*1 OPT
 CHARACTER*10 DEG
 CHARACTER*10 AMP 
 CHARACTER*20 FILE, FILECAP 
!
! Prepares a GMT script for plotting the rechonstructed and 
! the original Ocean Function - 
! ---> Last changes GS Nov 14 2007
! ---> Revised June 25, 2008 for implementation of "Selen 2.6"
! ---> Revised July XX, 2008 for implementation of the "Zonal" OF
!
 open(1,file=file,status='unknown')
!
 If(opt=='r')then 
!
! --- Realistic OF 
!
   Write(1,*)"${GMT} gmtset PAPER_MEDIA A4+"
   Write(1,*)"${GMT} gmtset HEADER_FONT_SIZE 24p"
   Write(1,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
   Write(1,*)"${GMT} gmtset ANOT_FONT_SIZE 12p"
   Write(1,*)" "
!
   Write(1,*)"# Map of the original ocean function"
   Write(1,*)"${GMT} psbasemap -X3 -Y18 -P -Bf180a90/f90a60WSEN -R0/360/-85/85 -JQ180/16 -K > of.ps"
   Write(1,*)"${GMT} pscoast -R -JQ -Di -B -W1/255/255/255 -A1000 -S255/0/0 -G0/0/255 -O -K >> of.ps"
   Write(1,*)"echo '180 120 18 0 2 CM OCEAN FUNCTION 0 ka' | ${GMT} pstext -U/0.5/0.5/'SELEN 2.9' -N -R -JQ -O -K >> of.ps"
   Write(1,*)" "
!
   Write(1,*)"# Map of rechonstructed ocean function"
   Write(1,*)"${GMT} makecpt  -Cno_green -T-0.1/1.1/0.1 > pale.cpt"
   Write(1,*)"${GMT} psbasemap -X0 -Y-11 -P -B -R -JQ -K -O >> of.ps"
   Write(1,*)"${GMT} pscontour -I -JQ -R recof.dat -Cpale.cpt -O -K >> of.ps"
   Write(1,*)"${GMT} pscoast -R -JQ -Di -W1/255/255/255 -B -A1000 -O -K >> of.ps"
   Write(1,*)"echo '180 120 18 0 2 CM OCEAN FUNCTION 0 ka to degree ", trim(adjustl(deg)),"'", & 
   " | ${GMT} pstext -U/0.5/0.5/'SELEN 2.9' -N -R -JQ -O -K >> of.ps"
   Write(1,*)" "
!
   Write(1,*)"# A color table"
   Write(1,*)"${GMT} psscale -E -Cpale.cpt -B1f0.1a0.5g0.5/:O.F.: -D8/-2/8/1h -O >> of.ps"
!
 Else
!
! --- "Zonal" OF 
!
   Write(1,*)"${GMT} gmtset PAPER_MEDIA A4+"
   Write(1,*)"${GMT} gmtset HEADER_FONT_SIZE 24p"
   Write(1,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
   Write(1,*)"${GMT} gmtset ANOT_FONT_SIZE 12p"
   Write(1,*)" "
!
   Write(1,*)"# Map of the original ocean function"
   Write(1,*)"${GMT} psbasemap -X3 -Y18 -P -Bf180a90/f90a60WSEN -Gred -R0/360/-85/85 -JQ180/16 -K > of.ps"
   Write(1,*)"${GMT} psxy -R -JQ ", trim(adjustl(FILECAP)), " -M -W1/255 -Gblue -L -B -A -O -K >> of.ps"
   Write(1,*)"echo '180 120 18 0 2 CM OCEAN FUNCTION' | ${GMT} pstext -U/0.5/0.5/'SELEN 2.9' -N -R -JQ -O -K >> of.ps"
   Write(1,*)" "
!
   Write(1,*)"# Map of rechonstructed ocean function"
   Write(1,*)"${GMT} makecpt  -Cno_green -T0/1/0.1 > pale.cpt"
   Write(1,*)"${GMT} psbasemap -X0 -Y-11 -P -B -R -JQ -K -O >> of.ps"
   Write(1,*)"${GMT} pscontour -I -JQ -R recof.dat -Cpale.cpt -O -K >> of.ps" 
   Write(1,*)"${GMT} psxy -R -JQ ", trim(adjustl(FILECAP)), " -M -W1/255 -B -A -O -K >> of.ps"
   Write(1,*)"echo '180 120 18 0 2 CM OCEAN FUNCTION to degree ", trim(adjustl(deg)),"'", & 
   " | ${GMT} pstext -U/0.5/0.5/'SELEN 2.9' -N -R -JQ -O -K >> of.ps"
   Write(1,*)" "
!
   Write(1,*)"# A color table"
   Write(1,*)"${GMT} psscale -E -Cpale.cpt -B1f0.1a0.5g0.5/:O.F.: -D8/-2/8/1h -O >> of.ps"
!
 Endif
!
 Close(1) 
!
 END SUBROUTINE MAKE_OFMAP
!
!
!
!
!
!
 SUBROUTINE MAKE_PXMAP(NRES, FILE)
 IMPLICIT NONE
 INTEGER NP, NRES
 CHARACTER*20 FILE
 CHARACTER*20, PARAMETER :: SIZE_OF_PIXELS="-Sc0.04" 
! 
!
! Prepares a GMT script for plotting  wet, dry, and global 
! pixels distribution  ---- Last changes GS Feb 21 2008 ----
!
! Last reviewed July 1, 2008 - GS for port to SELEN 2.6 
! Now the script generated TWO plots: 1) wet/dry distributions
! and 2) a global "spherical" view of pixels - 
! 
 NP=2*nres*(nres-1)*20+12
!
Open(11,file=file,status='unknown') 

Write(11,*)"${GMT} gmtset PAPER_MEDIA A4+" 
Write(11,*)"${GMT} gmtset HEADER_FONT_SIZE 24p"
Write(11,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
Write(11,*)"${GMT} gmtset ANOT_FONT_SIZE 12p"
!
Write(11,*)"echo '     - wet pixels'"
Write(11,*)""
Write(11,*)"# Map of wet pixels distribution"	
Write(11,*)""	
Write(11,*)"${GMT} psbasemap -Y18 -P -Bf180a90/f90a60WSEN -R0/360/-85/85 -JQ180/16 -K > px.ps"
Write(11,*)"${GMT} pscoast -B -R -O -K -W1/255/0/0 -JQ -Df -A1000 >> px.ps"
Write(11,*)"${GMT} psxy -G255/0/0 -H4 weta.dat -O ", trim(adjustl(size_of_pixels)), " -JQ -R -K >> px.ps"
Write(11,*)"echo '180 120 18 0 0 CM WET pixels' | ${GMT} pstext -N -R -JQ -O -K >> px.ps"
Write(11,'(a27,i3,a4,i6,a1,a59)') "echo '10 -70 14 0 2 LM RES=", nres, ", N=", NP, "'", & 
" | ${GMT} pstext -N -R -JQ -O -K -G255/0/0 -W255 >> px.ps"
!
Write(11,*)"echo '     - dry pixels'"
Write(11,*)""
Write(11,*)"# Map of dry pixels distribution"	
Write(11,*)""	
Write(11,*)"${GMT} psbasemap -X0 -Y-11 -P -Bf180a90/f90a60WSEN -R0/360/-85/85 -JQ180/16 -K -O >>  px.ps"
Write(11,*)"${GMT} pscoast -B -R -O -K -W1/0/0/255 -JQ -Df -A1000 >> px.ps"
Write(11,*)"${GMT} psxy -G0/0/255 -H4 drya.dat -O ", trim(adjustl(size_of_pixels)), " -JQ -R -K >> px.ps"
Write(11,*)"${GMT} psbasemap -U/6/-2/'SELEN 2.9' -B -R -P -JQ -K -O >> px.ps"
Write(11,*)"echo '180 120 18 0 0 CM DRY pixels' | ${GMT} pstext -N -R -JQ -O -K >> px.ps"
Write(11,'(a27,i3,a4,i6,a1,a59)') "echo '10 -70 14 0 2 LM RES=", nres, ", N=", NP, "'", & 
" | ${GMT} pstext -N -R -JQ -O -G0/0/255 -W255 >> px.ps"
!
Write(11,*)"echo '     - Spherical map of wet and dry pixels'"
Write(11,*)""
Write(11,*)"# Spherical map of dry & wet pixels distribution"	
Write(11,*)""
Write(11,*)"${GMT} psbasemap -X4 -Y12 -U/0/-3/'SELEN 2.9' -Bg60f60/g60f60 -R0/360/-80/80 -P -JG-20/-44/14 -K > px-sphere.ps"
Write(11,*)"${GMT} pscoast -B -R -O -K -W1/220/220/220 -JG -Df -S0/60/255 -G100/250/100 -A1000 >> px-sphere.ps"
Write(11,*)"${GMT} psxy  -G255 -H4 pxa.dat -O -K ", trim(adjustl(size_of_pixels)),  " -JG -R >> px-sphere.ps"
Write(11,'(a27,i3,a4,i6,a1,a56)') "echo '0 -110 14 0 2 LM RES=", nres, ", N=", NP, "'", & 
" | ${GMT} pstext -N -R -JQ -O -G0 -W255 >> px-sphere.ps"
!
 close(11)
! 
  end subroutine make_pxmap
!
!
!
!
!
!
 SUBROUTINE MAKE_RSLSCA (NV, CODE, RSL_FILE, RUN, NINC, NRSL, TITLICE, RESOLUTION, & 
			          ITER, MODE, DEGREE, VSTRING, FILE_GMT, VISCO_FILENAME)
!
! Prepares a simple GMT script for drawing a scatterplot of 
! global RSL predictions vs. observations  - GS 31.10-.2007 
!
! Deeply Revised on July 2008 for the implementation of SELEN 2.6
! Also revised on April 2010 GS 
!
 IMPLICIT NONE
 INTEGER, PARAMETER :: MAXDATA=100000
 REAL*8 TIME, RSL, DMIN, DMAX, TMIN, TMAX, RANGE, DELTA 
 CHARACTER*100 VSTRING, VISCO_FILENAME
 CHARACTER*20 FILE_GMT, YMINC, YMAXC, TMINC, TMAXC, DELTAC, R_OPTION
 CHARACTER*30 RSL_FILE 
 CHARACTER*10 RESOLUTION, DEGREE, TITLICE, STRING 
 CHARACTER*50  B_OPTION
 CHARACTER*3 RUN
 CHARACTER*3 NV, CODE 
 CHARACTER*3 NINC
 CHARACTER*1 ITER, MODE
 INTEGER I, J, K, NRSL 
!
!
! ==== Creates a GMT script for a scatterplot of RSL data 
! 	 
! --- Determines min and max of cumulative RSL data 
open(11,file='scatter-data.dat',status='unknown') 
dmax=-99999.
dmin=+99999.
tmax=-99999.
tmin=+99999.
do j=1, maxdata
       read(11,*,end=88) time, rsl 
       if(rsl<=dmin)  dmin=rsl 
       if(rsl>=dmax)  dmax=rsl 
       if(time<=tmin) tmin=time
       if(time>=tmax) tmax=time        
enddo
88 continue
   close(11)  
!
	range=dmax-dmin  
!	
if(range>=0  .and.range<=50) delta=10  
if(range>=50 .and.range<100) delta=10	
if(range>=100.and.range<200) delta=20
if(range>=200.and.range<300) delta=50
if(range>=300.and.range<400) delta=100
if(range>400) 		     delta=150 
!
OPEN  (10,FILE='junk.dat',STATUS='unknown')
WRITE (10,'(i10)') int(dmin-delta/2.)
WRITE (10,'(i10)') int(dmax+delta/2.)
WRITE (10,'(i10)') int(delta) 
WRITE (10,'(i10)') int(tmin)-1		
WRITE (10,'(i10)') int(tmax)+1	
 CLOSE(10)
OPEN  (10,FILE='junk.dat',STATUS='unknown')
READ (10,'(a10)') yminc
READ (10,'(a10)') ymaxc
READ (10,'(a10)') deltac
READ (10,'(a10)') tminc 
READ (10,'(a10)') tmaxc 	
 CLOSE(10)
!
! --- Target GMT file for scatterplot of RSL data 
!   
 open (19,file=file_gmt,status='unknown')
 Write(19,*)"${GMT} gmtset PAPER_MEDIA A4+"
 Write(19,*)"${GMT} gmtset HEADER_FONT_SIZE 24p"
 Write(19,*)"${GMT} gmtset FRAME_WIDTH 0.1c"
 Write(19,*)"${GMT} gmtset ANOT_FONT_SIZE 16p"
 Write(19,*)"${GMT} gmtset LABEL_FONT_SIZE 16p"
!
 r_option = & 
 "-R"//trim(adjustl(tminc))//"/"//trim(adjustl(tmaxc))//"/"//trim(adjustl(yminc))//"/"//trim(adjustl(ymaxc))
!
 b_option = & 
 "-Ba2f1/a"//trim(adjustl(deltac))//"f"//trim(adjustl(deltac))//"WSen:'RSL (m)':"
!
 Write(19,*)"${GMT} psbasemap -X4.5 -Y18 ", b_option, r_option," -P -JX12/7 -K > plot.ps" 
 Write(19,*)"${GMT} psxy scatter-data.dat ", " -B -R -JX -Ss0.1 -G0 -O -K >> plot.ps" 
 Write(19,*)"${GMT} psxy scatter-data.dat ", "-Ey0.25/2 -B -R -JX -O -K >> plot.ps" 
 Write(19,*)"echo '", int(tmin)+(int(tmax)-int(tmin))/2., (11./10.)*int(dmax+delta/2.), " 22 0 0 BC RSL data from file: ", & 
 trim(adjustl(RSL_FILE)), "'", " | ${GMT} pstext -JX -R -O -K -N >> plot.ps "  
!
 b_option = & 
 "-Ba2f1:'time (ka)':/a"//trim(adjustl(deltac))//"f"//trim(adjustl(deltac))//"WSen:'RSL (m)':"
!
 Write(19,*)"${GMT} psbasemap -U/-2/-3/'SELEN 2.9' -Y-11 ", b_option, r_option, " -P -JX -O -K >> plot.ps"
 Write(19,*)"${GMT} psxy scatter-pred.dat -B -R -JX -Sc0.12 -G0 -O -K >> plot.ps"
 Write(19,*)"echo '", int(tmin)+(int(tmax)-int(tmin))/2., (13./10.)*int(dmax+delta/2.), " 22 0 2 BC RSL predictions for: ", & 
 trim(adjustl(TITLICE)), "'", " | ${GMT} pstext -JX -R -O -K -N >> plot.ps "   
!
 open (8,file='title.tmp',status='unknown')  
 If(code=='-1')then 
	   CONTINUE
 Else 
 Write(8,*) int(tmin)+(int(tmax)-int(tmin))/2., (11./10.)*int(dmax+delta/2.), " 13 0 0 BC ",& 
    " -LMAX=", trim(adjustl(DEGREE)),& 
    " -RES=", trim(adjustl(RESOLUTION)),& 
    " -NV=", trim(adjustl(NV)),& 
    " -CODE=", trim(adjustl(CODE)),& 
    " -MODE=", trim(adjustl(MODE)),&   
    " -ITER=", trim(adjustl(ITER)),& 
    "  -Visco: ", trim(adjustl(VSTRING))  
 close(8) 
 Endif 
 Write(19,*)"${GMT} pstext title.tmp -N -JX -R -O >> plot.ps"
 Write(19,*)"mv plot.ps scatter-plot.ps"
!
 Close(19) 
!
 END SUBROUTINE MAKE_RSLSCA 
!
!
!
!
!
!
 SUBROUTINE MAKE_STOKES (RUN, NV, CODE, TITLICE, RESOLUTION, & 
 			ITER, MODE, DEGREE, VSTRING, SHORT, FILE_GMT)
!
! --- Creates a GMT script for plotting the rates - of - variation
!     of cosine and sine stokes coeffciients -GS November 12 2007
!
!     Revised JULY 2008 for version 2.6 of SELEN - 
!     Revised Luky 2010 - Fully normalized coefficients---
!
 IMPLICIT NONE
 CHARACTER*10 TITLICE,  DEGREE, RESOLUTION     
 CHARACTER*20 R_OPTION, G_OPTION, W_OPTION, H_OPTION, FILE_GMT 
 CHARACTER*150 B_OPTION
 CHARACTER*30 VSTRING
 CHARACTER*60 TITRE 
 CHARACTER*3 RUN, STRING
 CHARACTER*3 NV, CODE
 CHARACTER*1 ITER, MODE 
 CHARACTER*100 SHORT
 INTEGER L, M, J_INDEX 
!
 open(9,file=file_gmt,status='unknown')
!
 Write(9,*) "${GMT} gmtset PAPER_MEDIA A4+"
 Write(9,*) "${GMT} gmtset ANOT_FONT_SIZE 16p"
 Write(9,*) "${GMT} gmtset LABEL_FONT_SIZE 20p"
 Write(9,*) " "
 Write(9,*) "# ------ A plot of dot(c_lm,s_lm) vs degree (at present time)  ----" 		      
 Write(9,*) " "
!
! --- Recommended y-range for plot 
 R_OPTION="-R0/48/-2/2"      
!
! --- Header lines in file "stokes.dat" 
 H_OPTION="-H15" 
!
! --- Plot title  
 TITRE="'Rate of change of fully-normalized Stokes coefficients'"
!
! --- B-option 
 B_OPTION="-Ba4f2:'j(l,m)=l(l+1)/2+m+1':/f0.1a0.5:'d/dt (c@-lm@-, s@-lm@-) x 10@+-11@+/yr'"//":WSne"
! 
!
!
 Write(9,*) "" ; Write(9,*) "#--- Base of plot"
 Write(9,*) "${GMT} psbasemap -X6 -Y4 -U'SELEN 3.2' "//trim(adjustl(R_OPTION))//" "//B_OPTION//"-JX18/14 -K > stokes.ps"
!
! A title (revised) 
 open(4,file='title_stokes.tmp',status='unknown') 
   write(4,*) "24 2.4 22 0 1 CB Rate of change of the fully-normalized Stokes coefficients"
 close(4) 
 Write(9,*) "${GMT} pstext title_stokes.tmp -N -JX -R -G0 -O -K >> stokes.ps"
!
 Write(9,*) "" ; Write(9,*) "#--- Extracts cosine and sine coefficients"
 Write(9,*) "awk '{print $1, $4}' stokes.dat > cosine.tmp"
 Write(9,*) "awk '{print $1, $5}' stokes.dat >  sine.tmp"
!
 Write(9,*) "" ; Write(9,*) "#--- Plots labels for zonal degrees..."
!
!
  open(4,file='stokes0.tmp',status='unknown') 
  do 2 l=2, 48
  do 2 m=0, l
 	if(m==0) then 
	Write(4,'(i4,1x,a16,i2,a1,i1)') j_index(l,m), " -1.8 12 0 1 CM ", l, "-", m 
	Write(4,'(i4,1x,a16,i2,a1,i1)') j_index(l,m), " +1.8 12 0 1 CM ", l, "-", m 
	Endif
2 continue 
  close(4) 
 Write(9,*) "${GMT} pstext stokes0.tmp -JX -R -G0 -O -W220 -K >> stokes.ps" 
!
 Write(9,*) "" ; Write(9,*) "#--- Plots open squares for sine components"
 Write(9,*) "${GMT} psxy  sine.tmp ", trim(adjustl(H_OPTION)), " -B -R -JX  -W4 -K -O >> stokes.ps"
 Write(9,*) "${GMT} psxy  sine.tmp ", trim(adjustl(H_OPTION)), " -B -R -JX -Ss0.35 -G0 -K -O >> stokes.ps"
 Write(9,*) "${GMT} psxy  sine.tmp ", trim(adjustl(H_OPTION)), " -B -R -JX -Ss0.25 -G255 -K -O >> stokes.ps"
!
 Write(9,*) "" ; Write(9,*) "#--- Plots filled squares for cosine components"
 Write(9,*) "${GMT} psxy cosine.tmp ", trim(adjustl(H_OPTION)), " -B -R -JX -Ss0.35 -G0 -K -O >> stokes.ps"
 Write(9,*) "${GMT} psxy cosine.tmp ", trim(adjustl(H_OPTION)), " -B -R -JX  -W3 -K -O >> stokes.ps" 
!

 If(code=='-1') then 
 open(4,file='stokes1.tmp',status='unknown') 
 Write(4,*) "24 2.2 14 0 2 BC -Ice model: ", trim(adjustl(TITLICE)), & 
 " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
 " -ALMA rheology:", trim(adjustl(SHORT)), &  
 " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER)) 
 close(4) 
 else
 open(4,file='stokes1.tmp',status='unknown') 
 Write(4,*) "24 2.2 14 0 2 BC -Ice model: ", trim(adjustl(TITLICE)), & 
 " -Viscosity profile: ", trim(adjustl(VSTRING)),&
 " -LMAX=", trim(adjustl(DEGREE)),  " -RES=",  trim(adjustl(RESOLUTION)),& 	  
 " -NV=",   trim(adjustl(NV)),      " -CODE=", trim(adjustl(CODE)),& 
 " -MODE=", trim(adjustl(MODE)),    " -ITER=", trim(adjustl(ITER))    
 close(4) 
 endif
!
!
 Write(9,*) "" ; Write(9,*) "#--- A subtitle with some parameters"
 Write(9,*) "${GMT} pstext stokes1.tmp -N -JX -R -G0 -O -K >> stokes.ps" 
!
 open(4,file='stokes2.tmp',status='unknown') ; write(4,*) "32 -1.00" ; close(4) 
 open(4,file='stokes3.tmp',status='unknown') ; write(4,*) "32 -1.25" ; close(4)  
 open(4,file='stokes4.tmp',status='unknown') 
 write(4,*) "33 -1.00  12 0 1 ML cosine" 
 write(4,*) "33 -1.25  12 0 1 ML  sine" 
 close(4) 
! 
 Write(9,*) "" ; Write(9,*) "#--- Plots a legend"
 Write(9,*) "${GMT} psxy stokes2.tmp -N -JX -R -G0 -Ss0.35 -K -O >> stokes.ps" 
 Write(9,*) "${GMT} psxy stokes3.tmp -N -JX -R -G0 -Ss0.35 -K -O >> stokes.ps" 
 Write(9,*) "${GMT} psxy stokes3.tmp -N -JX -R -G255 -Ss0.25 -K -O >> stokes.ps" 
 Write(9,*) "${GMT} pstext stokes4.tmp -N -JX -R -G0 -O  >> stokes.ps"
!
 CLOSE(9) 
!
 END SUBROUTINE MAKE_STOKES 
!
!
!
!
!
