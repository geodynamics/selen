!
! RSL.F90 
!
! Last modified GS 04-11-2008 "Intel port"
! Also modified by GS on July 2008 for v. 2.6
! Re-touched on August 2008 for version 2.7 
! Re-touched on February 6 2009 for version 2.7 (on the way back from Liverpool) 
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! Reviewed GS August 2009 -  
! *** Reviewed GS & FC November 2009 - Porting under gfortran  
! *** Reviewed GS April 2010 - Alma 
! *** Reviewed DM March 2013 - Dynamic memory allocation
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
! ----------------------------------------------------------------
!  Computes the RSL curves at the sites of the Peltier database... 
! ----------------------------------------------------------------
! 
! Input files:
!	- sealevel.dat
! 	- shs.bin
! 	- shrsl.bin
!
!
! Outout files:
!	- rslp-*.dat     ! RSL predictions 
! 	- rsld-*.dat     ! RSL data 
!	- mis.dat
! 	- gmis.dat
!	- scatter-pred.dat
!	- scatter-data.dat
!
!
! INCLUDE "harmonics.f90"
 PROGRAM SL
 IMPLICIT NONE
 INCLUDE "data.inc"
 CHARACTER*22 FILENAME, INPUT_FILENAME, OUTPUT_FILENAME
 CHARACTER*50 AJUNK
 CHARACTER(50), ALLOCATABLE :: TITREA(:), TITREB(:)
 CHARACTER(100), ALLOCATABLE :: TITRE(:)
 CHARACTER(22), ALLOCATABLE :: NAME(:)
 CHARACTER*20 DATE, TIMC       
 CHARACTER*3 STRING, CJUNK  
 CHARACTER*40 LITHOC
 CHARACTER(40), ALLOCATABLE :: VVVV(:)     
 INTEGER, PARAMETER :: MAXD =100  
 INTEGER I, J, K, K1, K2, DOM
 INTEGER, ALLOCATABLE :: CODE(:), NDATA(:)
 REAL*8, ALLOCATABLE ::    TIME(:,:)    	! times bp for which data are available at site#k
 REAL*8, ALLOCATABLE ::   DTIME(:,:)    	! uncertainties on the above times...
 REAL*8, ALLOCATABLE ::     RSL(:,:)    	! rsl datum for site k ad the times above
 REAL*8, ALLOCATABLE ::    DRSL(:,:)    	! uncertainties on the rsl datum
 REAL*8, ALLOCATABLE :: LONS(:), LATS(:), MISFIT(:)
 REAL*8 GLOBAL_MISFIT, INTERP  
 REAL*8, ALLOCATABLE :: SLC(:,:), PREDICTED_RSL(:,:) 
 COMPLEX*16, ALLOCATABLE :: YY(:,:), S(:,:) 
 REAL*8 TIME_RSL1, TIME_RSL2 
 CHARACTER*10  LATSC10, LONSC10 
 CHARACTER*200 LINEP
 CHARACTER*100 SS(2)
 INTEGER NOUT
!
!
!-- Allocate memory
!
 ALLOCATE( TITREA(NRSL), TITREB(NRSL) )
 ALLOCATE( TITRE(NRSL) )
 ALLOCATE( NAME(NRSL) )
 ALLOCATE( VVVV(NV) )
 ALLOCATE( CODE(NRSL), NDATA(NRSL) )
 ALLOCATE(  TIME(NRSL,0:MAXD) )
 ALLOCATE( DTIME(NRSL,0:MAXD) )
 ALLOCATE(   RSL(NRSL,0:MAXD) )
 ALLOCATE(  DRSL(NRSL,0:MAXD) )
 ALLOCATE( LONS(NRSL), LATS(NRSL), MISFIT(NRSL) )
 ALLOCATE( SLC(NRSL,0:NN), PREDICTED_RSL(NRSL,0:NN) )
 ALLOCATE( YY(JMAX,NRSL), S(JMAX,0:NN) )
!    
!
!
!
!
!--- rsl.f: reading the coordinates of the RSL sites
!
! *********************************************
	If    (RSL_DATABASE_FORMAT=='0') then 
! *********************************************
!
   	OPEN(1,FILE=RSL_DATABASE,STATUS='unknown')
!
   	i=1 
2  	READ (1, 6,END=1) code(i), lats(i), lons(i), ndata(i), NAME(i)     
   	IF(lons(i)<=0.) lons(i)=360.+lons(i)		    
   	do k=1, ndata(i) 
        	   READ (1 ,*) time(i,k), dtime(i,k), rsl(i,k), drsl(i,k)
   	enddo
   	i=i+1
   	IF(i <= nrsl) GOTO 2
1  	CLOSE(1)
6  	FORMAT(1X,I3,1X,F5.1,1X,F6.1,1X,I2,1X,A22)
!
! *********************************************
	Elseif(RSL_DATABASE_FORMAT=='1') then
! *********************************************
!
 		OPEN(1, FILE=RSL_DATABASE,STATUS='unknown')
!
   		do i = 1, nrsl 
		    ndata(i)=1
		    read(1,*)code(i)
	            read(1,'(a50)')titrea(i)
		    AJUNK=TITREA(I) 
		    read(1,'(a50)')titreb(i)
		    titre(i)=TRIM(ADJUSTL(TITREA(I)))//" / "//TRIM(ADJUSTL(TITREB(I)))
	            read(1,*) lats(i), lons(i)  
                         if(lons(i)<=0.) lons(i)=lons(i)+360. 	
!		    		    
			IF     (AJUNK(1:1)=="-") THEN 			
				READ(1,*) time(i,1), dtime(i,1), rsl(i,1), drsl(i,1) 
!							
			ELSEIF (AJUNK(1:1)=="*") THEN 
!			
                                READ(1,*)  TIME_RSL1, TIME_RSL2, rsl(i,1), drsl(i,1)  
!
				time(i,1)  =   MIN(TIME_RSL1,TIME_RSL2)+ &
				              (MAX(TIME_RSL1,TIME_RSL2)-MIN(TIME_RSL1,TIME_RSL2))/2.
!					     
			        dtime(i,1) =  (MAX(TIME_RSL1,TIME_RSL2)-MIN(TIME_RSL1,TIME_RSL2))/2.			
!	
			ELSEIF (AJUNK(1:1)/="*".AND.AJUNK(1:1)/="-") THEN 
				Write(88,*) "The RSL database is badly configured"
				Write(88,*) "The program will STOP --------------"
				Write(*, *) "The RSL database is badly configured"
				Write(*, *) "The program will STOP --------------"
				Call stop_config 
				Stop	
				
			ENDIF	
!
		    read(1,'(a3)')cjunk 
		enddo 
		close(1) 
!		
! *********************************************
	Elseif(RSL_DATABASE_FORMAT=='2') then
! *********************************************
!
 		OPEN(1, FILE=RSL_DATABASE,STATUS='unknown')
!
   		do i = 1, nrsl 
		    ndata(i)=1
		    read(1,*)code(i)
	            read(1,'(a50)')titrea(i)
		    AJUNK=TITREA(I) 
		    read(1,'(a50)')titreb(i)
		    titre(i)=TRIM(ADJUSTL(TITREA(I)))//" / "//TRIM(ADJUSTL(TITREB(I)))
!
			read(1,'(a200)') linep 	
!		
			call scan_string (linep, 2, ss, nout)				
!
			LATSC10=trim(adjustl(ss(1))) 
			LONSC10=trim(adjustl(ss(2)))			
!
 			call FROMDDMMSS_2_DEGREES (LATSC10, LATS(I)) 
 			call FROMDDMMSS_2_DEGREES (LONSC10, LONS(I)) 
!
!read(1,*) lats(i), lons(i)  
!
                         if(lons(i)<=0.) lons(i)=lons(i)+360. 	
!		    		    
			IF     (AJUNK(1:1)=="-") THEN 			
				READ(1,*) time(i,1), dtime(i,1), rsl(i,1), drsl(i,1) 
!							
			ELSEIF (AJUNK(1:1)=="*") THEN 
!			
                                READ(1,*)  TIME_RSL1, TIME_RSL2, rsl(i,1), drsl(i,1)  
!
				time(i,1)  =   MIN(TIME_RSL1,TIME_RSL2)+ &
				              (MAX(TIME_RSL1,TIME_RSL2)-MIN(TIME_RSL1,TIME_RSL2))/2.
!					     
			        dtime(i,1) =  (MAX(TIME_RSL1,TIME_RSL2)-MIN(TIME_RSL1,TIME_RSL2))/2.			
!	
			ELSEIF (AJUNK(1:1)/="*".AND.AJUNK(1:1)/="-") THEN 
				Write(88,*) "The RSL database is badly configured"
				Write(88,*) "The program will STOP --------------"
				Write(*, *) "The RSL database is badly configured"
				Write(*, *) "The program will STOP --------------"
				Call stop_config 
				Stop	
				
			ENDIF	
!
		    read(1,'(a3)')cjunk 
		enddo 
		close(1) 

	Endif
!
!
!
!
! --- rsl.f: reading the CSHs on file ', filename 
	open(7,file='shrsl.bin',status='unknown',form='unformatted') 
	read(7) yy
	close(7)
!
! --- rsl.f: Rescaling the CSHs  
	do j=1,jmax 
		yy(j,:) = yy(j,:)*(2-dom(j)) 
	enddo
!
!
! --- rsl.f: reading the sealevel CSH coefficients ', filename 
	filename='shs.bin'
	open(3,file=filename,status='unknown',form='unformatted') 
	read(3) S
	close(3)
!
!
! --- rsl.f: sealevel change at the sites...'
	do i=1, nrsl
		slc(i,:) = 0. 
	do j=1, jmax
		slc(i,:) = slc(i,:) + real(s(j,:)*yy(j,i)) 
	enddo
	enddo
!
 write(*,*) '    - Computing sinthetic RSL curves'
!
 do 11 i=1, nrsl
!
  	OPEN  (10,FILE='junk.dat',STATUS='unknown'); WRITE (10,'(i3)') code(i); CLOSE(10)
  	OPEN  (10,FILE='junk.dat',STATUS='unknown'); READ  (10,'(a3)') string;  CLOSE(10)
!
  	OPEN  (10,FILE='rslp'//'-'//string//'.dat',STATUS='unknown')
!
	IF          (RSL_DATABASE_FORMAT=='0') then	 
  	WRITE (10,'(a1,1x,a22,a5,a3,1x,a8,2(1x,f10.5))') &
              '#', NAME(i), 'code=', string, 'lon/lat:', lons(i), lats(i)	      

	      elseif(RSL_DATABASE_FORMAT=='1'.or.RSL_DATABASE_FORMAT=='2') then 	 
        WRITE (10,'(a1,1x,a100,a5,a3,1x,a8,2(1x,f10.5))') &
              '#', TITRE(i), 'code=', string, 'lon/lat:', lons(i), lats(i)
	      	      
	ENDIF
!	      
  	WRITE(10,'(a1,a18)') '#', '  kyrs BP, RSL (m)'
!
	do k=0, nn
	predicted_rsl(i,k)= -(slc(i,nn)-slc(i,nn-k))
		WRITE(10,*) float(k), predicted_rsl(i,k)	
	enddo
!
 close(10)
 11 CONTINUE 
!
!write(*,*) '    - Misfit between data and predictions' 
!
!------- File for individual misfits 
  open(66,file='mis.dat',status='unknown')
!
!------- File for global misfits 
  open(67,file='gmis.dat',status='unknown')
  global_misfit=0. 
!
!------- File for scattered RSL predictions 
 open(81,file='scatter-pred.dat',status='unknown') 
!
!------- File for scattered RSL predictions 
 open(82,file='rsl-obs-vs-predictions.dat',status='unknown')
!
! --- Header information about the SELEN settings 
!
 open(76, file=VISCO_MODEL,status='unknown') 
 read(76, '(a40)') lithoc 
 	do i=1, nv 
          read(76, '(a40)') vvvv(i)  	 
	enddo
 close(76) 
 call DATE_AND_TIME (date,timc) 
WRITE(82,*)'# SELEN 3.1 - g95 version'     
WRITE(82,*)'# ', date(1:4), '.', date(5:6), '.', date(7:8), & 
           ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)
WRITE(82,*)"# Ice model: ", trim(adjustl(ice_model))
WRITE(82,*)"# RSL database ", RSL_DATABASE 
!
	If(NV.NE.1.AND.CDE.NE.-1) THEN
		Write(82,*)"# Computations by TABOO " 
		WRITE(82,*)"# Number of mantle layers: ", NV 
		WRITE(82,*)"# Model code (see TABOO User guide): ", CDE
		WRITE(82,*)"# Viscosity from bottom to top (* 1E21 Pa.s): "
		WRITE(82,*)"# LITHO thickness (km): ", LITHOC 
		WRITE(82,*)"# Viscosity model: ", visco_model
		DO I=1, NV 
   			WRITE(82,*)"# ", vvvv(i) 
		ENDDO
	ELSE 
		Write(82,*)"# Computations by ALMA - See log files for model details" 
		WRITE(82,*)"# Rheological profile: ", visco_model
	ENDIF 
!
WRITE(82,*)"# SLE mode of solution: ", IMODE
WRITE(82,*)"# SLE iterations: ", SMAX 
WRITE(82,*)"# Maximum harmonic degree: ", LMAX 
WRITE(82,*)"# Tegmark resolution: ", RES 
WRITE(82,*)"#-------------------------------------------------------------------------------------------"
WRITE(82,*)"#  - - - - - - - - RSL observations vs. predictions for all the sites in database - - - - - "
WRITE(82,*)"#										                "
WRITE(82,*)"#  Code   lon     lat     time BP  time ERR  Obs RSL   RSL Err  RSL Pred     RSL Site name  "
WRITE(82,*)"#         deg     deg      kyrs      kyrs       m         m         m                       "  			     
WRITE(82,*)"#-------------------------------------------------------------------------------------------"
WRITE(82,*)"#	 "
!
  do j=1, nrsl  
  	misfit(j)=0.0    		
  	do i=1, ndata(j)       
  	k1 = int(time(j,i)/1000.0)   
  	k2 = k1+1                   
!
!------linearly interpolated value 
  	interp = predicted_rsl(j,k1)+(time(j,i)/1000.0-float(k1)  )*&
             (predicted_rsl(j,k2)-predicted_rsl(j,k1))     	       
  	misfit(j)=misfit(j)+ (interp-rsl(j,i))**2/(drsl(j,i)**2)
!	
	Write(81,*) time(j,i)/1000., interp 				
!
	if(i==1) then 
		if(rsl_database_format=='0')     then 
		  Write(82,'(2x,i4,2(f8.3),5(f10.3),6x,a22)') & 
		  code(j), & 
		  lons(j), & 
		  lats(j), & 
		  time(j,i)/1e3, & 
		  dtime(j,i)/1e3, & 
		  rsl(j,i), & 
		  drsl(j,i), & 
		  interp, & 
		  name(j) 
		elseif(rsl_database_format=='1'.or.rsl_database_format=='2') then
		  Write(82,'(2x,i4,2(f8.3),5(f10.3),6x,a100)') & 
		  code(j), & 
		  lons(j), & 
		  lats(j), & 
		  time(j,i)/1e3, & 
		  dtime(j,i)/1e3, & 
		  rsl(j,i), & 
		  drsl(j,i), & 
		  interp, & 
		  titre(j)
		endif				
	else 
	Write(82,'(2x,a4,1x,a7,a7,1x,5(f10.3),6x,a10)') & 
		"   ", & 
		"        ", & 
		"        ", & 
		time(j,i)/1e3, & 
                dtime(j,i)/1e3, & 
		rsl(j,i), & 
		drsl(j,i), & 
		interp, & 
		"        " 
	endif	
!
        enddo
!
!       Write(82,*) " "				
!
  	misfit(j)=misfit(j)/float(ndata(j))    
!	    
  	global_misfit=global_misfit+misfit(j)
!
  	write(66,'(f14.8,i4,1x,a20)') misfit(j), code(j), name(j) 
!
  enddo
!
  close(81) 
!
WRITE(82,*)" "
WRITE(82,*)"#-------------------------------------------------------------------------------------------"
WRITE(82,*)"#  Code   lon     lat     time BP  time ERR  Obs RSL   RSL Err  RSL Pred     RSL Site name  "
WRITE(82,*)"#         deg     deg      kyrs      kyrs       m         m         m                       "  			     
WRITE(82,*)"#-------------------------------------------------------------------------------------------"
!
  close(82) 
!
  global_misfit=global_misfit/float(nrsl) 
!
  write(67,*) global_misfit 
!
  close(66)
  close(67) 
!                                                                            
!
!
!write(*,*) '    - rsl.f: Preparing the rsld files...' !
do 12 i=1, nrsl
  	OPEN  (10,FILE='junk.dat',STATUS='unknown'); WRITE (10,'(i3)') code(i); CLOSE (10)
  	OPEN  (10,FILE='junk.dat',STATUS='unknown'); READ  (10,'(a3)') string;  CLOSE (10)
  	OPEN (10,FILE='rsld-'//string//'.dat',STATUS='unknown')
!
  	WRITE(10,'(a1,1x,a22,1x,a3,a11,2(1x,f10.5))') &
              '#', NAME(i), string, '  lon/lat= ', lons(i), lats(i)
  	WRITE(10,'(a1,a29)') '#','  kyrs BP   RSL (m)  DRSL (m)'
!
  	do j=1, ndata(i)
  		write (10 ,'(3(f10.5))') time(i,j)/1.e3, rsl(i,j), drsl(i,j)
  	enddo
  	CLOSE(10)
12 CONTINUE 
!
!
 DEALLOCATE( TITREA, TITREB )
 DEALLOCATE( TITRE )
 DEALLOCATE( NAME )
 DEALLOCATE( VVVV )
 DEALLOCATE( CODE, NDATA )
 DEALLOCATE(  TIME )
 DEALLOCATE( DTIME )
 DEALLOCATE(   RSL )
 DEALLOCATE(  DRSL )
 DEALLOCATE( LONS, LATS, MISFIT )
 DEALLOCATE( SLC, PREDICTED_RSL )
 DEALLOCATE( YY, S )
!
!
 END PROGRAM SL
!
