!
! This is file "TGAUGES.F90" 
! 
! Last modified GS 04-11-2008 "Intel port"
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! Reviewed again by GS on April 2010 Porting under g95 
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
! -------------------------------------------------------------------------
!  Computes the present-day rate of sealevel change at the PSMSL sites ... 
! -------------------------------------------------------------------------
!
! Input files:
!	- rlr-trends.txt  
! 	- shs.bin
! 	- shn.bin
! 	- shu.bin
! 	- shtidegauges.bin
!
! Output files:
!       - psmsl*.dat
!
!
!  INCLUDE "harmonics.f90"
  PROGRAM PSMSL
  IMPLICIT NONE
  INCLUDE "data.inc"
  CHARACTER*1, PARAMETER :: NORTH='N', SOUTH='S', WEST='W', EAST='E'
  INTEGER, PARAMETER :: MAXN=10000  ! A large number   
  CHARACTER*100 RIGA		    ! A row 
  CHARACTER*8 PCODE                 ! PSMSL code  
  CHARACTER*3 GLOSS_CODE            ! GLOSS code 
  CHARACTER*3 YEARS                 ! Number of years 
  CHARACTER*11 RANGE_YEARS          ! Range of years 
  CHARACTER*9 DATUM 	            ! Trend 
  CHARACTER*3 PLUS_MINUS 	    ! +/- 
  CHARACTER*5 ERROR	            ! Error 
  CHARACTER*7 STDV 	            ! Std. deviation of residues 							
  CHARACTER*2 LAT1, LAT2            ! Latitide degrees and minutes 	
  CHARACTER*3 LON1	            ! Longitude degrees 
  CHARACTER*2 LON2 		    ! Longitude minutes 
  CHARACTER*30 NAME                 ! Name of the station 
  CHARACTER*1 WLAT, WLON            ! Where is the station (e, w, s, n)  
  CHARACTER*20 DATE, TIMC           ! DATE AND TIME
  REAL*8 LAT, NLAT, DLAT            ! Latitude  
  REAL*8 LON, NLON, DLON            ! Longitude 
  REAL*8 PSDOT, PNDOT, PUDOT        ! Predicted rate 
  INTEGER I, J, N, NH, DOM
  COMPLEX*16 Y(JMAX), SA(JMAX,0:NN), & 
  		     NA(JMAX,0:NN), & 
		     UA(JMAX,0:NN)  
  REAL*8, PARAMETER :: DDELTA=DELTA 
!
!
!
!
! --- SH file for PSMSL sites (input) 
        open(17,file='shtidegauges.bin',status='unknown',form='unformatted')  
!
! --- Writes an header for the output file 	
	OPEN (99,FILE='ptidegauges.dat',STATUS='unknown')
!
! --- Header information about the SELEN settings 
        call DATE_AND_TIME (date,timc)      
        Write(99,*) date(1:4), '.', date(5:6), '.', date(7:8), & 
	          ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)
	Write(99,*) "Ice model: ", trim(adjustl(ice_model))
	If(NV.NE.1.AND.CDE.NE.-1) THEN
		Write(99,*) "Number of mantle layers: ", NV 
		Write(99,*) "Model code (see TABOO User guide): ", CDE
	ELSE
		Write(99,*) "Computations by ALMA - See log files for model details" 
	ENDIF 
	Write(99,*) "Viscosity model: ", visco_model
!
	open(32,file=visco_model,status='unknown') 
	do i=1,NV+1 
	read(32,'(a100)')riga 
	If(NV.NE.1.AND.CDE.NE.-1) THEN
	if(i==1)  then 
		Write(99,*) "Thickness of the lithosphere (km): ", trim(adjustl(riga))
		  else
		Write(99,*) "Viscosity (bottom-to-top, Haskell units) ", trim(adjustl(riga))
	          Endif	
	ENDIF	
	enddo ; close(32) 
	Write (99,*) "SLE iterations: ", SMAX 
	Write (99,*) "SLE mode of solution: ", IMODE
	Write (99,*) "Maximum harmonic degree: ", LMAX 
	Write (99,*) "Tegmark resolution: ", RES 
	WRITE (99,*) ''
	WRITE (99,*) &
	' PSMSL    yrs    range       trend  error     lon       lat       s-dot     n-dot     u-dot    PSMSL station name' 
	WRITE (99,*) &
	' code             yrs            mm/yr        deg       deg       mm/yr     mm/yr     mm/yr    =================='
	WRITE (99,*) ''
!
! --- tgauges.f: reading sealevel CSH coefficients from file', filename 
       open(3,file='shs.bin',status='unknown',form='unformatted') ; read(3) SA ; close(3)
       open(3,file='shn.bin',status='unknown',form='unformatted') ; read(3) NA ; close(3)
       open(3,file='shu.bin',status='unknown',form='unformatted') ; read(3) UA ; close(3)
!
! --- Counting the header lines (, beginning with '#'), and the data lines
       open(10,file=TGAUGES_DATABASE,status='unknown')
       nh=0 
       n=0
       do i=1, maxn 
       		read(10,'(a100)',end=1) riga
       		if(riga(1:1)=='#') nh=nh+1
       		if(riga(1:1)/='#') n=n+1
       enddo 
 1     close(10)  
!write(*,*) &
!'    - tgauges.f: file ', trim(adjustl(TGAUGES_DATABASE)), ' contains', nh, 'header lines'  
!write(*,*) &
!'    - tgauges.f: there are data for', n, 'stations in file ',  trim(adjustl(TGAUGES_DATABASE))
!
!
! --- Opening the PSMSL data file... 
       open(10,file=TGAUGES_DATABASE,status='unknown')
!
! --- Reading again the header lines  
       do i=1, nh ; read(10,'(a100)',end=1) riga ; enddo  
!
! === Reading all lines. GS = October 5 2007 - 
!
       write(*,*) '    - Rate of sealevel change at site'
!
       do 5 i=1, n 
!
       		read  (10,200) pcode, gloss_code, years, range_years, & 
   		      	       datum , plus_minus, error, stdv, & 
			       lat1, lat2, wlat, lon1, lon2, wlon, name        
!
! --- Extracts lon and lat of the station from the texst strings... 
                call  find_lon_lat (lat1, lat2, wlat, lon1, lon2, wlon, lon, lat) 
!
! --- SH st the current PSMSL site...
      	        read(17) y
!    
! --- Predicted rate at the site 'i' 
		If(n>=50) then 
      	        If(mod(i,50)==0.or.i==1) & 
		write(*,*) '    - ', i, 'of', n, name 
		elseif(10<=n.and.n<=50)then
      	        If(mod(i,10)==0.or.i==1) & 
		write(*,*) '    - ', i, 'of', n, name 		
		elseif(n<10)then
		write(*,*) '    - ', i, 'of', n, name 			       
		endif 
!		
		psdot = 0.
		pndot = 0. 
		pudot = 0. 
		do j=1, jmax
			psdot = psdot + real((sa(j,nn)-sa(j,nn-1))*y(j)) / ddelta
   		    pndot = pndot + real((na(j,nn)-na(j,nn-1))*y(j)) / ddelta
			pudot = pudot + real((ua(j,nn)-ua(j,nn-1))*y(j)) / ddelta
		enddo
!
  		WRITE (99,201)  pcode, & 
				years, & 
				range_years, & 
				datum, & 
				error, & 
				lon, lat, & 
				psdot, pndot, pudot, & 
				name 
!
 5    continue 
!
!
! --- Writing format  
201 format(a8,2x,a3,2x,a11,1x,a8,2x,a5,2x,f8.3,2x,f8.3,2x,f8.3,2x,f8.3,2x,f8.3,4x,a30) 
!
! --- Reading format 
200    format(a8,1x,a3,1x,a3,2x,a11,1x,a8,1x, & 
              a3,1x,a5,1x,a7,3x,a2,1x,a2,1x,a1,1x,& 
              a3,1x,a2,1x,a1,3x,a30)
!
 close(17)
 close(10)
 close(99)  
!
 END PROGRAM PSMSL
!
!
!
