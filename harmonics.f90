!
! This is the include file "HARMONICS.F90"
!
! HARMONICS.F90 is a set of Fortran utilities for Spherical Harmonics (SH) 
! analysis, and some other useful stuff, listed below. 
!
! *** Modified by GS & FC 06-13-2008  = INTEL PORT = V 2.6 
! *** Romberg integration program added July 09, 2008. 
! *** On July 19, 2008 I have added to_real10 --- GS 
! *** Updated with new routines on August 2008, for v. 2.7
! *** Some typos corrected on September 29, 2008. 
! *** On october 2, I have moved here the pixelization routines 
! *** Reviewed GS & FC July 2009 -  "Varying coastlines" & ALMA coupling
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Revised GS May 2010 - g95 - Ice breaker routine included 
!
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
! This file contains several units *mainly* related to spherical harmonics:  
! 
! - Subroutine PLegendre (*)          Legendre polynomials 
! - Subroutine PlmBar (*) 	      Normalized Legendre functions
! - Subroutine PlmBar_d1 (*) 	      Derivatives wrt colatitude of Legendre functions
! - Subroutine PLMBAR_MOD             (Modified) normalized Legendre functions 	
! - Subroutine PLMBAR_MOD_D1          (Modified) derivatives of Legendre functions 	
! - Subroutine Pleg                   Legendre polynomials              
! - Function FUNCLM                   Integrand for rectangular shape factors
! - Function J_INDEX 	              J=l*(l+1)/2+m+1 
! - Function DOM                      Delta(0,m) 
! - Function MJ			      Retrives m (order) from J
! - Function LJ 		      Retrives l (degree) from J
! - Subroutine HARMO 	              Complex spherical harmonics (CSH)
! - Subroutine grad_theta_HARMO       Colatitude component of gradient of CSH  
! - Subroutine grad_lambda_HARMO      Longitude component of gradient of CSH  
! - Subroutine Ice_count              Counts ice elements 
! - Subroutine Scan_string            Scans for valid strings
! - Subroutine find_lon_lat           Lon-lat infor from PSMSL trend file  
! - Function FCT 		      Factorial of small integers
! - Subroutine STOP_CONFIG            Stop configuration if errors are detected
! - Subroutine CHAR100_2_REAL         Converts from CHARACTER*100 to FLOATING POINT 
! - Subroutine CHAR10_2_REAL          Converts from CHARACTER*10 to FLOATING POINT 
! - Subroutine CHAR10_2_INT           Converts from CHARACTER*10 to INTEGER
! - Subroutine CHAR3_2_INT            Converts from CHARACTER*3 to INTEGER
! - Subroutine CHAR2_2_INT            Converts from CHARACTER*2 to INTEGER
! - Subroutine INT_2_CHAR3            Converts from INTEGER to CHARACTER*3
! - Subroutine ROMINT (**)            "Romberg integration"
! - Subroutine BUBBLE_SORT_MOD        A bubble sorting routine 
! - Subroutine BUBBLE_SORT_PAR        A parallel bubble sorting routine 
! - Subroutine QSORT_MOD              A QuickSort sorting routine 
! - Subroutine SWAP_MOD               A swapping routine 
! - Subroutine FINDPX                 Finding the pixels given the resolution 
! - Subroutine AVERAGE_EARTH_DENSITY  Computes the Earth's average density 
! - Subroutine FromDDMMSS_2_degrees   Transforms DD:MM:SS into decimal degrees
! - ALL the Tegmark routines          The full pixelization code by M. Tegmark 
! - Function SINDD                    Sine with argument in degrees 
! - Function COSDD                    Cosine with argument in degrees 
! - Subroutine ICE_BREAKER            Breaks the ice array in small parts 
! - Subroutine INTEGER_2_CHAR2        Int to CHAR*2, with left zerp padding

!
!   (*)  adapted from SHTOOLS, Copyright (c) 2005, Mark A. Wieczorek
!   (**) implemented by the ACM TOMS algorithm #TOMS351. The "early" algorithm 
!        used here is by John Burkardt (see http://people.scs.fsu.edu/~burkardt). 
!        No modifications have been made on the code available from http://
!        people.scs.fsu.edu/~burkardt/f77_src/toms351/toms351.html.  
! ----------------------------------------------------------------------------
! 
!
!
!
!
!
!
!
	FUNCTION SINDD(ALFA)
	implicit NONE
!
! --- Assuming that ALFA is given in DEGREES, it returns the SINE of ALFA
!     *********************** GS and FC 11-02-2009 **********************
!	
	REAL*8,  PARAMETER :: PI=3.14159265358979323840 
	REAL*8 SINDD, ALFA, ALFARAD              
	ALFARAD=ALFA*PI/180.	
	SINDD=SIN(ALFARAD)
	
	END FUNCTION SINDD 
!
!
!
!
!
!
	FUNCTION COSDD(ALFA)
	implicit NONE
!
! --- Assuming that ALFA is given in DEGREES, it returns the COSINE of ALFA
!     ************************ GS and FC 11-02-2009 ***********************
!	
	REAL*8,  PARAMETER :: PI=3.14159265358979323840 
	REAL*8 COSDD, ALFA, ALFARAD               
	ALFARAD=ALFA*PI/180.	
	COSDD=COS(ALFARAD)	
!
	END FUNCTION COSDD
!
!
!
!
!
!
        Subroutine STOP_CONFIG 
	implicit NONE
!
! --- Creates a STOP file that stops the execution of makeseles.sh 
!     ******************* GS and FC 13-12-2007 *******************
!
	open (11,file='stop.dat',status='new')
	Write(11,*) "Errors have been detected - SELEN will stop"
	close(11)
!
	write(* ,*) 'STOP_CONFIG: SELEN will STOP!' 
	write(88,*) 'STOP_CONFIG: SELEN will STOP!' 
	STOP 
!
	End subroutine STOP_CONFIG 	
!
!
!
!
!
!
 SUBROUTINE ICE_COUNT(FILENAME, IMUL, HEADER_LINES, N)
 IMPLICIT NONE
 CHARACTER*100 FILENAME, RIGA  
 INTEGER IMUL, HEADER_LINES, IUNIT, J, N
!
!
!
! Given the a valid ice file ''filename'' this routine counts the number of 
! ice elements as N (output). Imul (input) is multiplicity (i.e., number of 
! lines relative to one ice element, e. g., 1 or 2). Header_lines (input) is 
! the  number  of header lines in file Filename   ----  GS Oct 19 2007 ---- 
!
! write(*,'(a80)') trim(adjust(filename)) 
! read(*,*) 

open(18,file=filename,status='old')
n=0
do j=1, 100000
	read(18,'(a80)',end=1) riga
	n =n +1
enddo
1 n =(n - header_lines)/imul ; close(18) 
!
 END SUBROUTINE ICE_COUNT
!
!
!
!
!
!
 SUBROUTINE CHAR10_2_FLOAT(STRING, FLOAT_VAL)
 IMPLICIT NONE
 INTEGER, PARAMETER :: IU=30
 CHARACTER*10 STRING 
 REAL*8 FLOAT_VAL  
!
! Given the input CHARACTER*10 string STRING, this routine converts it in a floating 
! point value by copying it on a file and reading it back. ==== GS Feb 02 2007 ==== 
! Revised GS July 27, 2009. 
!
 open(iu,file='junk.dat',status='unknown') ; write(iu,*)   string ; close(iu) 
 open(iu,file='junk.dat',status='unknown') ; read(iu,*) float_val ; close(iu)  
!
 end SUBROUTINE CHAR10_2_FLOAT
!
!
!
!
!
!
 SUBROUTINE CHAR10_2_INT(STRING, INT_VAL)
 IMPLICIT NONE
 INTEGER, PARAMETER :: IU=30
 CHARACTER*10 STRING 
 INTEGER INT_VAL  
!
! Given the input CHARACTER*10 string STRING, this routine converts it in an 
! integer value by copying it on a file and reading it back. GS Oct 17 2007. 
! Revised GS July 27, 2009. 
!
 open(iu,file='junk.dat',status='unknown') ; write(iu,*) string ; close(iu) 
 open(iu,file='junk.dat',status='unknown') ; read(iu,*) int_val ; close(iu)  
!
 end SUBROUTINE CHAR10_2_INT
!
!
!
!
!
!
 SUBROUTINE CHAR2_2_INT(STRING, INT_VAL)
 IMPLICIT NONE
 INTEGER, PARAMETER :: IU=30
 CHARACTER*2 STRING 
 INTEGER INT_VAL  
!
! Given the input CHARACTER*2 string STRING, this routine converts it in an 
! integer value by copying it on a file and reading it back. GS Oct 17 2007. 
! Revised GS July 27, 2009. 
!
 open(iu,file='junk.dat',status='unknown') ; write(iu,*) string ; close(iu) 
 open(iu,file='junk.dat',status='unknown') ; read(iu,*) int_val ; close(iu)  
!
 end SUBROUTINE CHAR2_2_INT
!
!
!
!
!
!
 SUBROUTINE CHAR3_2_INT(STRING, INT_VAL)
 IMPLICIT NONE
 INTEGER, PARAMETER :: IU=30
 CHARACTER*3 STRING 
 INTEGER INT_VAL  
!
! Given the input CHARACTER*3 string STRING, this routine converts it in an 
! integer value by copying it on a file and reading it back. GS Oct 17 2007. 
! Revised GS July 27, 2009. 
!
 open(iu,file='junk.dat',status='unknown') ; write(iu,'(a3)') string ; close(iu) 
 open(iu,file='junk.dat',status='unknown') ; read (iu,*)     int_val ; close(iu)  
!
 end SUBROUTINE CHAR3_2_INT
!
!
!
!
!
!
 SUBROUTINE INT_2_CHAR3(INT_VAL,STRING)
 IMPLICIT NONE
 INTEGER, PARAMETER :: IU=30
 CHARACTER*3 STRING 
 INTEGER INT_VAL  
!
! Given the input CHARACTER*3 string STRING, this routine converts it in an 
! integer value by copying it on a file and reading it back. GS Oct 17 2007. 
! Revised GS July 27, 2009. 
!
 open(iu,file='junk.dat',status='unknown') ; write(iu,'(i3)')int_val ; close(iu)  
 open(iu,file='junk.dat',status='unknown') ; read (iu,'(a3)')string  ; close(iu) 
!
 end SUBROUTINE INT_2_CHAR3
!
!
!
!
!
!
 SUBROUTINE INTEGER_2_CHAR2(INT_VAL,STRING)
 IMPLICIT NONE
 INTEGER, PARAMETER :: IU=30
 CHARACTER*2 STRING 
 INTEGER INT_VAL  
!
! Given the input integer INT_VAL <= 99, this routine converts it 
! into a CHARACTER*2 string, with zero padding on left if INT_VAL <=9.  
! GS May 24 2010. 
!

 open (iu,file='junk.dat',status='unknown') 
       if(int_val<=9) write(iu,'(a1,i1)') '0', int_val  
       if(int_val> 9) write(iu,'(i2)')         int_val 
 close(iu)
 open (iu,file='junk.dat',status='unknown')
       read(iu,'(a2)') STRING
 close(iu) 
!
 end SUBROUTINE INTEGER_2_CHAR2
!
!
!
!
!
!
!
 SUBROUTINE CHAR100_2_REAL(STRING, REAL_VAL)
 IMPLICIT NONE
 INTEGER, PARAMETER :: IU=30
 CHARACTER*100 STRING 
 REAL*8 REAL_VAL  
!
! Given the input CHARACTER*100 string STRING, this routine converts it in Floating 
! point value by copying it on a file and reading it back. *** GS Oct 19 2007 *** 
! Revised GS July 27, 2009. 
!
 open(iu,file='junk.dat',status='unknown') ; write(iu,*)  string ; close(iu) 
 open(iu,file='junk.dat',status='unknown') ; read(iu,*) real_val;  close(iu)  
!
 end SUBROUTINE CHAR100_2_REAL
!
!
!
!
!
!
 SUBROUTINE CHAR10_2_REAL(STRING, REAL_VAL)
 IMPLICIT NONE
 INTEGER, PARAMETER :: IU=30
 CHARACTER*10 STRING 
 REAL*8 REAL_VAL  
!
! Given the input CHARACTER*10 string STRING, this routine converts it in Floating 
! point value by copying it on a file and reading it back. *** GS July 19 2008 *** 
! Revised GS July 27, 2009. 
!
 open(iu,file='junk.dat',status='unknown') ; write(iu,*)  string ; close(iu) 
 open(iu,file='junk.dat',status='unknown') ; read(iu,*) real_val;  close(iu)  
!
 end SUBROUTINE CHAR10_2_REAL
!
!
!
!
!
!
 SUBROUTINE SCAN_STRING (S, N, SS, M)
 IMPLICIT NONE 
 INTEGER I, J, N, M, APO(2*N) 
 CHARACTER*200 S ; CHARACTER*100 SS(N)
 CHARACTER*1, PARAMETER :: PRIME="'" 
!
! Given the input string s, this routine determines the sub-strings ss of s 
! which are delimited by consecutive pairs of primes ('). n (input) is the 
! expected number of substrings, and m (output) is the effective number 
! found. Apo is a pointer to the primes. The execution is stopped whenever 
! i) n=0, ii) m=0, iii) m/=n, or iv) if m is odd. The maximum lenght of the 
! substrings is 20, that of the input string is 200.  == GS Nov 17 2007 == 
!
! ----- Exits for n=0
 if(n==0) then 
 	  Write(88,*) "SCAN_STRING has nothing to do" 
 	  Write(88,*) "The program will STOP ----------------"
 	  Write(*, *) "SCAN_STRING has nothing to do" 
 	  Write(*, *) "The program will STOP ----------------"	   
          call stop_config 
	  Stop
 endif
!
! ----- Looking for primes (', not ") within the input string  
 i=0
 do j=1, 200
	if(s(j:j)==prime) then 
	 	i=i+1 ; apo(i)=j
	endif
 enddo 
 m=i/2 
!
! ----- Exits if no primes are found 
 if(i==0) then 
 	  Write(*,*) "SCAN_STRING has found NO primes" 
 	  Write(*,*) "The program will STOP ----------------"
 	  Write(*, *) "SCAN_STRING has found NO primes" 
	  Write(*,'(a100)') s 
 	  Write(*, *) "The program will STOP ----------------"	   
          call stop_config 
	  Stop
 endif 
!
! ----- Exits if an odd number of primes is found
 if(mod(i,2)/=0) then 
 	  Write(*,*) "SCAN_STRING has found an even number of primes:", i
 	  Write(*,*) "The program will STOP ----------------"
 	  Write(*, *) "SCAN_STRING has found an even number of primes:", i
 	  Write(*, *) "The program will STOP ----------------"	   
          call stop_config 
	  Stop
 endif
!
! ----- Exits if m/=n, otherwise determines the substrings
 if(m/=n) then 
 	  Write(*,*) "SCAN_STRING has found ", m, " substrings"
          Write(*,*) "SCAN_STRING  expected ", n, " substrings"
 	  Write(*,*) "The program will STOP ----------------"	
 	  Write(*, *) "SCAN_STRING has found ", m, " substrings"
          Write(*, *) "SCAN_STRING  expected ", n, " substrings"
 	  Write(*, *) "The program will STOP ----------------"
          call stop_config 
	  Stop
 else
!
!        Write(*,*) "SCAN_STRING has found ", m, " substrings"
!        Write(*,*) "SCAN_STRING  expected ", n, " substrings"	 
!
         do i=1, n 
	        ss(i)=s(apo(2*i-1)+1:apo(2*i)-1)
         enddo
 endif	
! 
 end subroutine scan_string 
!
!
!
!
!
!
 SUBROUTINE SCAN_STRING_SOFT (S, N, SS, M)
 IMPLICIT NONE 
 INTEGER I, J, N, M, APO(2*N) 
 CHARACTER*200 S ; CHARACTER*100 SS(N)
 CHARACTER*1, PARAMETER :: PRIME="'" 
!
! Given the input string s, this routine determines the sub-strings ss of s 
! which are delimited by consecutive pairs of primes ('). n (input) is the 
! expected number of substrings, and m (output) is the effective number 
! found. Apo is a pointer to the primes. The execution is stopped whenever 
! i) n=0, ii) m=0, iii) m/=n, or iv) if m is odd. The maximum lenght of the 
! substrings is 20, that of the input string is 200.  == GS Nov 17 2007 == 
!
! ----- Exits for n=0
 if(n==0) then 
 	  Write(88,*) "SCAN_STRING has nothing to do" 
 	  Write(88,*) "The program will STOP ----------------"
 	  Write(*, *) "SCAN_STRING has nothing to do" 
 	  Write(*, *) "The program will STOP ----------------"	   
          call stop_config 
	  Stop
 endif
!
! ----- Looking for primes (', not ") within the input string  
 i=0
 do j=1, 200
	if(s(j:j)==prime) then 
	 	i=i+1 ; apo(i)=j
	endif
 enddo 
 m=i/2 
!
!
         do i=1, n 
	        ss(i)=s(apo(2*i-1)+1:apo(2*i)-1)
         enddo
! 
 end subroutine scan_string_soft
!
!
!
!
!
!
  SUBROUTINE FIND_LON_LAT (LAT1, LAT2, WLAT, LON1, LON2, WLON, LON, LAT)
!
! This routine converts latitudes and longitudes written in the DM format in
! a floating point format. This is useful for using the data contained in file
! 'rlr-trends.txt', which comes form the PSMSL web site. The input character 
! variables are LAT1, LAT2, WLAT (s/n), LON1, LON2, WLON (e/w). The output is
! the couple (lon,lat) for a given site, in REAL*8 format. GS 06/10/2007 - 
!
! Fixed a significant BUG here on May 2010. Tere was a mistake in the conversion. 
!
! 
  IMPLICIT NONE 
  CHARACTER*1, PARAMETER :: NORTH='N', SOUTH='S', WEST='W', EAST='E'
  CHARACTER*1 WLAT, WLON         ! where is the station (e, w, s, n)  
  CHARACTER*2 LAT1, LAT2         ! latitide degrees and minutes 
  CHARACTER*3 LON1	         ! longitude degrees 
  CHARACTER*2 LON2 		 ! longitude minutes 
  REAL*8 DLON, DLAT, NLON, NLAT  ! auxiliary variables  	
  REAL*8 LON, LAT 		 ! Longitude and latitude (output) 
!
! --- Converts latitudes in floating point values. 
!     According to our conventions, latitude spans 
!     from -90 to 90 degrees 
!
! -- First part (degrees) 
   if(lat1(1:1)/='0') then 
   		      open(21,file='junk.dat',status='unknown') ; write(21,*) lat1(1:2) ; close(21) 
   		      open(21,file='junk.dat',status='unknown') ; read(21,*) nlat ; close(21) 
   elseif(lat1(2:2)/='0') then 
   		      open(21,file='junk.dat',status='unknown') ; write(21,*) lat1(2:2) ; close(21) 
   		      open(21,file='junk.dat',status='unknown') ; read(21,*) nlat ; close(21) 		      	
   endif 
!
! -- Second part (minutes) 
   if(lat2(1:1)/='0') then 
      		      open(21,file='junk.dat',status='unknown') ; write(21,*) lat2(1:2) ; close(21) 
   		      open(21,file='junk.dat',status='unknown') ; read(21,*) dlat ; close(21) 
   elseif(lat2(2:2)/='0') then 
		      open(21,file='junk.dat',status='unknown') ; write(21,*) lat2(2:2) ; close(21) 
   		      open(21,file='junk.dat',status='unknown') ; read(21,*) dlat ; close(21) 		      	
   endif   
   dlat=dlat/60. 
!
! -- Computing latitude... 
   	    if(wlat=='N') lat = +(nlat+dlat)  
   	    if(wlat=='S') lat = -(nlat+dlat)
!   		
!
! --- Converts longitude in floating point values. 
!     According to our conventions, longitude spans 
!     from 0 to 360 degrees  
!
! -- First part (degrees) 
   if(lon1(1:1)/='0') then 
   	              open(21,file='junk.dat',status='unknown') ; write(21,*) lon1(1:3) ; close(21) 
   		      open(21,file='junk.dat',status='unknown') ; read(21,*) nlon ; close(21) 
   elseif(lon1(2:2)/='0') then 
		      open(21,file='junk.dat',status='unknown') ; write(21,*) lon1(2:3) ; close(21) 
   		      open(21,file='junk.dat',status='unknown') ; read(21,*) nlon ; close(21) 		      	
   elseif(lon1(3:3)/='0') then 
		      open(21,file='junk.dat',status='unknown') ; write(21,*) lon1(3:3) ; close(21) 
   		      open(21,file='junk.dat',status='unknown') ; read(21,*) nlon ; close(21) 
   endif  
!
! -- Second part (minutes) 
   if(lon2(1:1)/='0') then 
   	              open(21,file='junk.dat',status='unknown') ; write(21,*) lon2(1:2) ; close(21) 
   		      open(21,file='junk.dat',status='unknown') ; read(21,*) dlon ; close(21) 
   elseif(lon2(2:2)/='0') then 
		      open(21,file='junk.dat',status='unknown') ; write(21,*) lon2(2:2) ; close(21) 
   		      open(21,file='junk.dat',status='unknown') ; read(21,*) dlon ; close(21) 		      	
   endif  
!
   if(lon1(1:1)==' '.and.lon1(2:2)=='0'.and.lon1(3:3)=='0'.and. &
                         lon2(1:1)=='0'.and.lon2(2:2)=='0')        dlon=0.
!			    
   dlon=dlon/60.
!
!
! -- Computing longitude... 
   	    if(wlon=='E') lon =    +(nlon+dlon)  
   	    if(wlon=='W') lon= 360.-(nlon+dlon)
!
   end subroutine find_lon_lat 
!
!
!
!
!
!
subroutine PLegendre(p, lmax, z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine evalutates all of the unnormalized Legendre polynomials 
!	up to degree lmax. 
!
!	Calling Parameters:
!		Out
!			p:	A vector of all unnormalized Legendgre polynomials evaluated at 
!				z up to lmax. The lenght must by greater or equal to (lmax+1).
!		IN
!			lmax:	Maximum degree to compute.
!			z:	Value within [-1, 1], cos(colatitude) or sin(latitude).
!
!	Notes:
!	
!	1.	The integral of Pl**2 over (-1,1) is 2/(2l+1).
!	2.	Values are calculated accoring to the following recursion scheme:
!			P_0(z) = 1.0, P_1(z) = z, and 
!			P_l(z) = (2l-1) * z * P_{l-1}(z) / l - (l-1) * P_{l-2}(z) / l
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek June 2004
!
! ----> Modified to SINGLE PRECISION by Giorgio Spada 2007 
! ----> Also modified for the management of ERROR conditions 
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	implicit none
	integer, intent(in) ::	lmax
	REAL*8, intent(out) ::	p(:)
       	REAL*8, intent(in) ::	z
       	REAL*8 ::	pm2, pm1, pl
      	integer ::	l

	if(size(p) < lmax+1) then
 	        Write(88,*) "Error --- PlegendreL" 
		Write(88,*) "P must be dimensioned as (LMAX+1) where LMAX is ", lmax
		Write(88,*) "Input array is dimensioned ", size(p)
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error --- PlegendreL" 
		Write(* ,*) "P must be dimensioned as (LMAX+1) where LMAX is ", lmax
		Write(* ,*) "Input array is dimensioned ", size(p)
		Write(* ,*) "The program will STOP ----------------"	   
                call stop_config 
	        Stop
     	elseif (lmax < 0) then 
	        Write(88,*) "Error --- PlegendreL" 
		Write(88,*) "LMAX must be greater than or equal to 0."
		Write(88,*) "Input value is ", lmax
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error --- PlegendreL" 
		Write(* ,*) "LMAX must be greater than or equal to 0."
		Write(* ,*) "Input value is ", lmax
		Write(* ,*) "The program will STOP ----------------"	   
                call stop_config 
	        Stop
     	elseif(abs(z) > 1.) then
	        Write(88,*) "Error --- PlegendreL" 
		Write(88,*) "ABS(Z) must be less than or equal to 1."
		Write(88,*) "Input value is ", z
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error --- PlegendreL" 
		Write(* ,*) "ABS(Z) must be less than or equal to 1."
		Write(* ,*) "Input value is ", z
		Write(* ,*) "The program will STOP ----------------"	   
                call stop_config 
	        Stop
     	endif
      	
   	pm2  = 1.
      	p(1) = 1.
      	
      	pm1  = z
      	p(2) = pm1
      	
      	do l = 2, lmax
         	pl = (float(2*l-1) * z * pm1 - float(l-1) * pm2)/float(l)
         	p(l+1) = pl
         	pm2  = pm1
         	pm1  = pl
      	enddo
!
end subroutine PLegendre
!
!
!
!
!
!
  SUBROUTINE PLEG (LMAX, LAT, LEGP)
  USE SHTOOLS 
  IMPLICIT NONE 
  INTEGER, PARAMETER :: LLMAX=512
  INTEGER L, LMAX 
  REAL*8 LAT, LEGP(0:LLMAX), P(1:LLMAX+1)
  REAL*8 COSDD
!
! Given latitude LAT (in degrees), computes *all* the Legendre polynomials P_{l}(x), 
! with x=cos(colatitude), to degree LMAX. Uses the SHTOOLS subroutine PLegendre in
! REAL*8 precision - Note that here the degree 0 is at place '0' in the array PLEG. 
! Last modified by GS 10-9-2007 -  
!
  If(abs(lat).gt.90.) then 
 	        Write(88,*) "Error in Sbr. PLEG: Latitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. PLEG: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1
  elseif(lmax>llmax) then 
		Write(88,*) "Error in Sbr. PLEG: The degree exceeds 512", lmax 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. PLEG: The degree exceeds 512", lmax 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2		
  Endif  
!
!
! --- Computes all the the P_l's to degree LMAX using Mark Wieczorek's code 
!     Note that p(1)=Leg(0), p(2)=Leg(1), ..., p(LMAX+1)=Leg(LMAX) 
!
  	call PLegendre(p, lmax, cosdd(90.-lat))
!  
!
! --- Shifts according to our conventions (degree 0 is in place 0 in the array)
!
	do l=0, lmax
   		legp(l)=p(l+1)
 	enddo 
! 
  END SUBROUTINE PLEG
!
!
!
!
!
!
  SUBROUTINE PLMBAR_MOD (LMAX, LAT, PLM)
  USE SHTOOLS 
  IMPLICIT NONE 
  INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
  INTEGER J, J_INDEX, MJ, JMAX, LMAX
  REAL*8 Z, LAT, PLM(JJMAX)
  REAL*8 COSDD
!
! Given latitude LAT (in degrees), computes *all* the fully normalized 
! associated Legendre functions \bar P_{lm}(x), with x=cos(colatitude), 
! to degree LMAX. Uses the SHTOOLS Legendre functions by PlmBar, and 
! rescales by sqrt(2-delta(0,m)).   Last modified by GS 10-9-2007 - 
!
!
! ---- Tests the latitude bounds 
!
  If(abs(lat).gt.90.) then 
 	        Write(88,*) "Error in Sbr. PLMBAR_MOD: Latitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. PLMBAR_MOD: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1
  elseif(lmax>llmax) then 
		Write(88,*) "Error in Sbr. PLMBAR_MOD: The degree exceeds 512", lmax 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. PLMBAR_MOD: The degree exceeds 512", lmax 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2		
  Endif  
!
!
!
! ---- "z" is cos(theta), theta is co-latitude 
!
	z=cosdd(90.-lat) 
!
! ---- Builds the SHTOOLS Legendre functions including the Condon-Shortley phase 
!
  	call PlmBar(plm, lmax, z, -1)  
!
!
! ---- Scales to obtain our Legendre functions, dividing by sqrt(2-delta(0,m))
!  
  	do j=1, j_index(lmax,lmax)
     		if(mj(j)/=0) plm(j)=plm(j)/sqrt(2.)
  	enddo
!
  END SUBROUTINE PLMBAR_MOD 
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
  SUBROUTINE PLMBAR_MOD_D1 (LMAX, LAT, DPLM)
  USE SHTOOLS 
  IMPLICIT NONE 
  INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
  INTEGER J, J_INDEX, MJ, JMAX, LMAX
  REAL*8 Z, LAT, PLM(JJMAX), DPLM(JJMAX)
  REAL*8 COSDD
!
! Given latitude LAT (in degrees), computes the detivatives wrt 'theta'
! of *all* the fully normalized associated Legendre functions \bar P_{lm}(x), 
! with x=cos(colatitude), to degree LMAX. Uses the SHTOOLS Legendre functions 
! by PlmBar_d1, rescales by sqrt(2-delta(0,m)), and transforms. 
!   
! Created by GS August 7 2008 for version 2.8 of SELEN  
!
!
! ---- Tests the latitude bounds 
!
  If(abs(lat).gt.90.) then 
 	        Write(88,*) "Error in Sbr. PLMBAR_MOD_D1: Latitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. PLMBAR_MOD_D1: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1
  elseif(lmax>llmax) then 
		Write(88,*) "Error in Sbr. PLMBAR_MOD_D1: The degree exceeds 512", lmax 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. PLMBAR_MOD_D1: The degree exceeds 512", lmax 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2		
  Endif  
!
!
! ---- "z" is cos(theta), theta is co-latitude 
!
	z=cosdd(90.-lat)   
!
! ---- Builds the SHTOOLS Legendre functions,
!      including the "Condon-Shortley phase"...  
!
  	call PlmBar_d1(plm, dplm, lmax, z, -1)  
!
!
! ---- Scales consistently to our convention for Legendre functions, 
!      dividing by factor sqrt(2-delta(0,m)), AND transofrms to obtain 
!      a derivative wrt 'theta' from the derivative wrt 'cos theta'
!  
  	do j=1, j_index(lmax,lmax)
     		if(mj(j)/=0) dplm(j)=dplm(j)/sqrt(2.)
	                     dplm(j)=-sqrt(1.-z**2)*dplm(j)		
  	enddo
!
  END SUBROUTINE PLMBAR_MOD_D1
!
!
!
!
!
!
RECURSIVE FUNCTION FCT (N) 
IMPLICIT NONE 
INTEGER J, N
REAL*8 FCT
!
! A naive recursive function for computing the
! factorial of small numbers (n<=24) GS 14-9-07
!
if(n<0.or.n>=24) then 
 	  Write(88,*) "<<N>> is out of bounds in sbr. factorial" 
 	  Write(88,*) "The program will STOP ----------------"
 	  Write(*, *) "<<N>> is out of bounds in sbr. factorial"  
 	  Write(*, *) "The program will STOP ----------------"	   
          call stop_config 
	  Stop
endif
!
if(n==0.or.n==1) then
	fct=1.
else 
	fct=float(n)  
	do j=n-1,1 
		fct=fct*float(j)
	enddo
endif 
!

	     
!if(n==0.or.n==1) then
!	fct=1.  
!		 else 	     
!	fct=float(n)*fct(n-1)
!endif  
!
END FUNCTION FCT 
!
!
!
!
!
!
FUNCTION J_INDEX(L,M) 
IMPLICIT NONE
INTEGER J_INDEX, L, M 
!
! given l (degree) and m(order), this function computes 
! the j-index "J", with J=l(l+1)/2+m+1 - GS 14-09-2007
!
 If(l<0.or.m<0.or.m>l) then 
	  Write(88,*) "J_INDEX: The degree and order are out of bounds" 
 	  Write(88,*) "The program will STOP -------------------------"
 	  Write(*, *) "J_INDEX: The degree and order are out of bounds"  
 	  Write(*, *) "The program will STOP -------------------------"	   
          call stop_config 
	  Stop
 endif
!
  J_INDEX = l*(l+1)/2+m+1
!
END FUNCTION J_INDEX 
!
!
!
!
!
!
   INTEGER FUNCTION DOM (J)
   IMPLICIT NONE 
   INTEGER J, MJ
!
! Given the index j=l(l+1)/2+m+1, returns 
! delta(0,m), with m=order - GS 14.09.07 
!
 If(j<0) then 
	  Write(88,*) "DOM: The J-index is out of bounds" 
 	  Write(88,*) "The program will STOP -----------"
 	  Write(*, *) "DOM: The J-index is out of bounds"  
 	  Write(*, *) "The program will STOP -----------"	   
          call stop_config 
	  Stop
 endif
!
   DOM=0 
   If(mj(j)==0) DOM=1
!
   END FUNCTION DOM
!
!
!
!
!
!
   INTEGER FUNCTION LJ (J)
   IMPLICIT NONE 
   INTEGER J, L, M, K  
   REAL*8 Z
!
! Given the harmonic index j=l(l+1)/2+m+1, returns l (degree) 
! -GS 10.09.07 
!
 If(j<0) then 
	  Write(88,*) "LJ: The J-index is out of bounds" 
 	  Write(88,*) "The program will STOP ----------"
 	  Write(*, *) "LJ: The J-index is out of bounds"  
 	  Write(*, *) "The program will STOP ----------"	   
          call stop_config 
	  Stop
 endif
!
   m=0
   k=j 
10 z=(-1.0+sqrt(8.*float(k)-7.))/2.
   lj=aint(z) 
   if(z-aint(z)==0.) return 
   k=k-1
   m=m+1
   goto 10
   return 
   end function LJ 
!   
!
!   
!
!
!   
   INTEGER FUNCTION MJ (J)
   IMPLICIT NONE 
   INTEGER J, L, M, K 
   REAL*8 Z
!
! Given the harmonic index j=l(l+1)/2+m+1, returns m (order) 
! -GS 10.09.07 
!
 If(j<0) then 
	  Write(88,*) "LJ: The J-index is out of bounds" 
 	  Write(88,*) "The program will STOP ----------"
 	  Write(*, *) "LJ: The J-index is out of bounds"  
 	  Write(*, *) "The program will STOP ----------"	   
          call stop_config 
	  Stop
 endif
!
   mj=0
   k=j 
10 z=(-1.0+sqrt(8.*float(k)-7.))/2.
   l=aint(z) 
   if(z-aint(z)==0.) return 
   k=k-1
   mj=mj+1
   goto 10
   return 
   end FUNCTION MJ 
!
!
!
!
!
!
 SUBROUTINE HARMO(LMAX, LON, LAT, ARMOY)
 USE SHTOOLS 
 IMPLICIT NONE 
 INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
 INTEGER J, MJ, J_INDEX, LMAX  
 COMPLEX*16 ARMOY(JJMAX) 
 REAL*8 LON, LAT, PLM(JJMAX) 
 REAL*8 COSDD, SINDD
!
!
! Given longitude LON and latitude LAT  - in degrees - this routine computes 
! *all* the 4-pi normalized complex harmonics $\cal Y_{lm}(\theta,\lambda)$, 
! with \theta=colatitude and \lambda=longitude, to degree LMAX given. Uses  
! the (modified) SHTOOLS Legendre functions by PlmBar (see Sbr PLMBAR_MOD). 
! - Last modified by GS 12-14-2007 - 
!
!
! ---- Tests the longitude and latitude bounds 
!
  If(abs(lat)>90.) then 
 	        Write(88,*) "Error in Sbr. HARMO: Latitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. HARMO: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1		
  elseif (lon<0.or.lon>360.) then
 	        Write(88,*) "Error in Sbr. HARMO: Longitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. HARMO: Longitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2	
  elseif(lmax>512) then 
		Write(88,*) "Error in Sbr. HARMO: The degree exceeds 512", lmax 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. HARMO: The degree exceeds 512", lmax 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
  Endif   
!  
!
! ---- Builds the SHTOOLS Legendre functions, with Condon-Shortley phase 
!
    CALL PLMBAR_MOD(LMAX, LAT, PLM)  
!
! 
! ---- Computes the 4-pi normalized Spherical Harmonics up to degree LMAX
!
    do j=1, j_index(lmax,lmax)
       	armoy(j) = plm(j)*cmplx(cosdd(float(mj(j))*lon),sindd(float(mj(j))*lon)) 
    enddo
!
END SUBROUTINE HARMO
!
!
!
!
!
!
 SUBROUTINE GRAD_THETA_HARMO(LMAX, LON, LAT, GRD_THETA_ARMOY)
 USE SHTOOLS 
 IMPLICIT NONE 
 INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
 INTEGER J, MJ, J_INDEX, LMAX  
 COMPLEX*16 GRD_THETA_ARMOY(JJMAX) 
 REAL*8 LON, LAT, DPLM(JJMAX) 
 REAL*8 COSDD, SINDD
!
!
! Given longitude LON and latitude LAT  - in degrees - this routine computes 
! the theta-component of the gradient of *all* the 4-pi normalized complex 
! harmonics $\cal Y_{lm}(\theta,\lambda)$, with \theta=colatitude and \lambda=
! longitude, to degree LMAX. This subroutine uses the (modified) SHTOOLS Legendre 
! functions by PlmBar (see Sbr PLMBAR_MOD). - Created by GS August 8, 2008 - 
!
!
! ---- Tests the longitude and latitude bounds 
!
  If(abs(lat)>90.) then 
 	        Write(88,*) "Error in Sbr. GRAD_THETA_HARMO: Latitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. GRAD_THETA_HARMO: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1		
  elseif (lon<0.or.lon>360.) then
 	        Write(88,*) "Error in Sbr. GRAD_THETA_HARMO: Longitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. GRAD_THETA_HARMO: Longitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2	
  elseif(lmax>512) then 
		Write(88,*) "Error in Sbr. GRAD_THETA_HARMO: The degree exceeds 512", lmax 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. GRAD_THETA_HARMO: The degree exceeds 512", lmax 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
  Endif   
!  
!
! ---- Builds the SHTOOLS Legendre functions, with Condon-Shortley phase 
!
    		CALL PLMBAR_MOD_D1(LMAX, LAT, DPLM)  
!
! 
! ---- Computes the 4-pi normalized Spherical Harmonics up to degree LMAX
!
    		do j=1, j_index(lmax,lmax)
       			grd_theta_armoy(j) = &
	        	dplm(j)*cmplx(cosdd(float(mj(j))*lon),sindd(float(mj(j))*lon)) 
    		enddo
!
END SUBROUTINE GRAD_THETA_HARMO
!
!
!
!
!
!
 SUBROUTINE GRAD_LAMBDA_HARMO(LMAX, LON, LAT, GRD_LAMBDA_ARMOY)
 USE SHTOOLS 
 IMPLICIT NONE 
 INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
 INTEGER J, MJ, J_INDEX, LMAX  
 COMPLEX*16, PARAMETER :: I = CMPLX(0.,1.)
 COMPLEX*16 GRD_LAMBDA_ARMOY(JJMAX) 
 REAL*8 LON, LAT, PLM(JJMAX) 
 REAL*8 COSDD, SINDD

!
!
! Given longitude LON and latitude LAT  - in degrees - this routine computes 
! the lambda-component of the gradient of *all* the 4-pi normalized complex 
! harmonics $\cal Y_{lm}(\theta,\lambda)$, with \theta=colatitude and \lambda=
! longitude, to degree LMAX given. Uses the (modified) SHTOOLS Legendre 
! functions by PlmBar (see Sbr PLMBAR_MOD_D1). - Created by GS August 8, 2008 - 
!
!
! ---- Tests the longitude and latitude bounds 
!
  If(abs(lat)>90.) then 
 	        Write(88,*) "Error in Sbr. GRAD_LAMBDA_HARMO: Latitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. GRAD_LAMBDA_HARMO: Latitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 1		
  elseif (lon<0.or.lon>360.) then
 	        Write(88,*) "Error in Sbr. GRAD_LAMBDA_HARMO: Longitude is out of range" 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. GRAD_LAMBDA_HARMO: Longitude is out of range" 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 2	
  elseif(lmax>512) then 
		Write(88,*) "Error in Sbr. GRAD_LAMBDA_HARMO: The degree exceeds 512", lmax 
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error in Sbr. GRAD_LAMBDA_HARMO: The degree exceeds 512", lmax 
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
  Endif   
!  
!
! ---- Builds the SHTOOLS Legendre functions, with Condon-Shortley phase 
!
    		CALL PLMBAR_MOD(LMAX, LAT, PLM)  
!
! 
! ---- Computes the 4-pi normalized Spherical Harmonics up to degree LMAX
!
    		do j=1, j_index(lmax,lmax)
       			grd_lambda_armoy(j) = &
	        			        plm(j)*I*mj(j)* &
			   			cmplx(cosdd(float(mj(j))*lon),sindd(float(mj(j))*lon))/ &
						sindd(90.-lat) 
    		enddo
!
END SUBROUTINE GRAD_LAMBDA_HARMO
!
!
!
!
!
!
SUBROUTINE PLMBAR(P, LMAX, Z, CSPHASE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function evalutates all of the normalized associated Legendre
!	functions up to degree lmax. The functions are initially scaled by 
!	10^280 sin^m in order to minimize the effects of underflow at large m 
!	near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299). 
!	On a Mac OSX system with a maximum allowable double precision value of 
!	2.225073858507203E-308 the scaled portion of the algorithm will not overflow 
!	for degrees less than or equal to 2800.
!
!	For each value of m, the rescaling factor is computed as rescalem=rescalem*sin(theta), 
!	with the intial value of rescalem being equal to 1/scalef (which is here equal 
!	to 10^280). This will gradually reduce this huge number to a tiny number, and will 
!	ultimately underflow. In order to prevent this underflow, when rescalem becomes less than
!	10^(-280), the subsequent rescaling factors of sin(theta) will be directly applied to Plm, and then this
!	number will be multipled by the old value of rescalem.
!
!	Temporary variables in saved in an allocated array. In order to explicitly deallocate this
!	memory, call this routine with a spherical harmonic degree of -1.
!
!	Calling Parameters:
!		OUT
!			p:		A vector of all associated Legendgre polynomials evaluated at 
!					z up to lmax. The length must by greater or equal to (lmax+1)*(lmax+2)/2.
!		OPTIONAL (IN)
!			csphase:	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!		IN
!			lmax:		Maximum spherical harmonic degree to compute.
!			z:		cos(colatitude) or sin(latitude).
!
!	Notes:
!	
!	1.	The employed normalization is the "geophysical convention." The integral of
!		(plm*cos(m theta))**2 or (plm*sin (m theta))**2 over all space is 4 pi.
!	2.	The integral of plm**2 over (-1,1) is 2 * (2 - delta(0,m))
!	3.	The index of the array p corresponds to l*(l+1)/2 + m + 1. As such
!		the array p should be dimensioned as (lmax+1)*(lmax+2)/2 in the 
!		calling routine.
!	4. 	The default is to exclude the Condon-Shortley phase of (-1)^m.
!
!
!	Dependencies:	CSPHASE_DEFAULT
!
!	Written by Mark Wieczorek September 25, 2005.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!
! ----> Modified by GS December 14, 2007 - error messages 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!
!	use SHTOOLS, only: CSPHASE_DEFAULT
!
!
	implicit none
	integer, intent(in) ::	lmax
	REAL*8, intent(out) ::	p(:)
       	REAL*8, intent(in) ::	z
       	integer, intent(in), optional :: csphase
       	REAL*8 ::	pmm, rescalem, phase, u, scalef
      	REAL*8, save, allocatable ::	f1(:), f2(:), sqr(:)
!!$OMP THREADPRIVATE(f1,f2,sqr)
      	integer ::	k, kstart, m, l, astat(3)
      	integer, save ::	lmax_old  = 0
!!$OMP THREADPRIVATE(lmax_old)
!	
        integer CSPHASE_DEFAULT
!
	if (lmax == -1) then
		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		lmax_old = 0
		return
	endif

	if (size(p) < (lmax+1)*(lmax+2)/2) then 
 	        Write(88,*) "Error --- PlmBar"
     		Write(88,*) "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 where LMAX is ", lmax
     		Write(88,*) "Input array is dimensioned ", size(p)
	 	Write(88,*) "The program will STOP ----------------"
 	        Write(* ,*) "Error --- PlmBar"
     		Write(* ,*) "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 where LMAX is ", lmax
     		Write(* ,*) "Input array is dimensioned ", size(p)
	 	Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 		
	elseif (lmax < 0) then 
     		Write(88,*) "Error --- PlmBar"
     		Write(88,*) "LMAX must be greater than or equal to 0."
     		Write(88,*) "Input value is ", lmax
		Write(88,*) "The program will STOP ----------------"
	     	Write(* ,*) "Error --- PlmBar"
     		Write(* ,*) "LMAX must be greater than or equal to 0."
     		Write(* ,*) "Input value is ", lmax
		Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 	
	elseif(abs(z) > 1.00) then
     		Write(88,*) "Error --- PlmBar"
     		Write(88,*) "ABS(Z) must be less than or equal to 1."
     		Write(88,*) "Input value is ", z
		Write(88,*) "The program will STOP ----------------"
	     	Write(* ,*) "Error --- PlmBar"
     		Write(* ,*) "ABS(Z) must be less than or equal to 1."
     		Write(* ,*) "Input value is ", z
		Write(*, *) "The program will STOP ----------------"
                call stop_config 
	        Stop 	
	endif     	
!     	
     	if (present(csphase)) then
     		if (csphase == -1) then
     			phase = -1.00
     		elseif (csphase == 1) then
     			phase = 1.00
     		else
     		Write(88,*) "PlmBar --- Error"
     		Write(88,*) "CSPHASE must be 1 (exclude) or -1 (include)."
     		Write(88,*) "Input value is ", csphase
		Write(88,*) "The program will STOP ----------------"
     		Write(* ,*) "PlmBar --- Error"
     		Write(* ,*) "CSPHASE must be 1 (exclude) or -1 (include)."
     		Write(* ,*) "Input value is ", csphase
		Write(* ,*) "The program will STOP ----------------"
                call stop_config 
	        Stop 
     		endif
     	else
     		phase = float(CSPHASE_DEFAULT)   !dble(CSPHASE_DEFAULT)
     	endif
     		
	scalef = 1.0e-20
	
	
	if (lmax /= lmax_old) then
		
		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		
		allocate(sqr(2*lmax+1), stat=astat(1))
		allocate(f1((lmax+1)*(lmax+2)/2), stat=astat(2))
		allocate(f2((lmax+1)*(lmax+2)/2), stat=astat(3))
		
		if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
     		Write(* ,*) "PlmBar --- Error"
     		Write(* ,*) "Problem allocating arrays SQR, F1 and F2", astat(1), astat(2), astat(3)
                call stop_config 
	        Stop 
		endif
			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		!	Precompute square roots of integers that are used several times.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		do l=1, 2*lmax+1
			sqr(l) = sqrt(float(l))   !sqrt(dble(l))
		enddo

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Precompute multiplicative factors used in recursion relationships
		! 		Plmbar(l,m) = x*f1(l,m)*Plmbar(l-1,m) - Plmbar(l-2,m)*f2(l,m)
		!		k = l*(l+1)/2 + m + 1
		!	Note that prefactors are not used for the case when m=l and m=l-1,
		!	as a different recursion is used for these two values.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		k = 3
	
		do l=2, lmax, 1
			k = k + 1
			f1(k) = sqr(2*l-1) * sqr(2*l+1) / float(l)   !dble(l)
			f2(k) = float(l-1) * sqr(2*l+1) / sqr(2*l-3) / float(l)   !dble(l)
			do m=1, l-2
				k = k+1
				f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                		f2(k) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                  			 / sqr(2*l-3) / sqr(l+m) / sqr(l-m) 
			enddo
			k = k + 2
		enddo
	
		lmax_old = lmax
	
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!	
	!	Calculate P(l,0). These are not scaled.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	u = sqrt((1.00-z)*(1.00+z)) ! sin(theta)
	
      	p(1) = 1.00
      	
      	if (lmax == 0) return
      	
      	p(2)  = sqr(3)*z
      	
      	k = 2

      	do l = 2, lmax, 1
         	k = k+l
         	p(k) = f1(k)*z*p(k-l)-f2(k)*p(k-2*l+1)
      	enddo

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate P(m,m), P(m+1,m), and P(l,m)
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
      	pmm = sqr(2)*scalef
      	rescalem = 1./scalef
      	kstart = 1

      	do m = 1, lmax - 1, 1
      	
      		rescalem = rescalem * u
      	
		! Calculate P(m,m)
        	kstart = kstart+m+1
         
         	pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
        	p(kstart) = pmm

		! Calculate P(m+1,m)
		k = kstart+m+1
	   	p(k) = z * sqr(2*m+3) * pmm

		! Calculate P(l,m)
               	do l = m+2, lmax, 1
               		k = k+l
                  	p(k) = z*f1(k)*p(k-l)-f2(k)*p(k-2*l+1)
                  	p(k-2*l+1) = p(k-2*l+1) * rescalem
               	enddo
               	
               	p(k) = p(k) * rescalem
               	p(k-lmax) = p(k-lmax) * rescalem
               	              
      	enddo
      	
      	! Calculate P(lmax,lmax)
      	
      	rescalem = rescalem * u
            	
        kstart = kstart+m+1
        p(kstart) = phase * pmm * sqr(2*lmax+1) / sqr(2*lmax) * rescalem
      		
end subroutine PlmBar
!
!
!
!
!
!
SUBROUTINE PLMBAR_D1(P, DP, LMAX, Z, CSPHASE)
!
! --- Implemented into SELEN by GS on August 7, 2008, with very slight 
!     modifications with respect to the original by M. Wieczoreck. I 
!     use REAL*8 precision, here, as for PLMBAR.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function evalutates all of the normalized associated Legendre
!	functions up to degree lmax. The functions are initially scaled by 
!	10^280 sin^m in order to minimize the effects of underflow at large m 
!	near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299). 
!	On a Mac OSX system with a maximum allowable double precision value of 
!	2.225073858507203E-308 the scaled portion of the algorithm will not overflow 
!	for degrees less than or equal to 2800.
!
!	For each value of m, the rescaling factor is computed as rescalem=rescalem*sin(theta), 
!	with the intial value of rescalem being equal to 1/scalef (which is here equal 
!	to 10^280). This will gradually reduce this huge number to a tiny number, and will 
!	ultimately underflow. In order to prevent this underflow, when rescalem becomes less than
!	10^(-280), the subsequent rescaling factors of sin(theta) will be directly applied to Plm, and then this
!	number will be multipled by the old value of rescalem.
!
!	Temporary variables in saved in an allocated array. In order to explicitly deallocate this
!	memory, call this routine with a spherical harmonic degree of -1.
!
!	Calling Parameters:
!		OUT
!			p:		A vector of all associated Legendgre polynomials evaluated at 
!					z up to lmax. The lenght must by greater or equal to (lmax+1)*(lmax+2)/2.
!			dp:		A vector of all first derivatives of the normalized Legendgre polynomials evaluated at 
!					z up to lmax with dimension (lmax+1).
!		IN
!			lmax:		Maximum spherical harmonic degree to compute.
!			z:		cos(colatitude) or sin(latitude).
!		OPTIONAL (IN)
!			csphase:	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!
!	Notes:
!	
!	1.	The employed normalization is the "geophysical convention." The integral of
!		(plm*cos(m theta))**2 or (plm*sin (m theta))**2 over all space is 4 pi.
!	2.	The integral of plm**2 over (-1,1) is 2 * (2 - delta(0,m))
!	3.	The index of the array p corresponds to l*(l+1)/2 + m + 1. As such
!		the array p should be dimensions as (lmax+1)*(lmax+2)/2 in the 
!		calling routine.
!	4. 	Derivatives are calculated using the unnormalized identities
!			P'l,m = ( (l+m) Pl-1,m - l z Plm ) / (1-z**2)	(for l>m), and
!			P'll = - l z Pll / (1-z**2)	(for l=m).
!	5.	The derivative is not defined at z=+-1 for all m>0, and is therefore not
!		calculated here.
!	6.	The default is to exlude the Condon-Shortley phase of (-1)^m.
!
!
!	Dependencies:	CSPHASE_DEFAULT
!
!	Written by Mark Wieczorek September 25, 2005.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!	
!	use SHTOOLS, only: CSPHASE_DEFAULT
!
	implicit none
	integer, intent(in) ::	lmax
	REAL*8, intent(out) ::	p(:), dp(:)
       	REAL*8, intent(in) ::	z
       	integer, intent(in), optional :: csphase
       	REAL*8 ::	pm2, pm1, pmm, plm, rescalem, phase, u, scalef
      	REAL*8, allocatable, save ::	f1(:), f2(:), sqr(:)
      	integer ::	k, kstart, m, l, sdim, astat(3)
      	integer, save ::	lmax_old = 0
!
        integer CSPHASE_DEFAULT
!
	if (lmax == -1) then
		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		lmax_old = 0
		return
	endif

      	sdim = (lmax+1)*(lmax+2)/2
      	
	if (size(p) < sdim) then 
		print*, "Error --- PlmBar_d1"
     		print*, "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 where LMAX is ", lmax
     		print*, "Input array is dimensioned ", size(p)
     		stop
     	elseif (size(dp) < sdim) then 
		print*, "Error --- PlmBar_d1"
     		print*, "DP must be dimensioned as (LMAX+1)*(LMAX+2)/2 where LMAX is ", lmax
     		print*, "Input array is dimensioned ", size(dp)
     		stop
     	elseif (lmax < 0) then 
     		print*, "Error --- PlmBar_d1"
     		print*, "LMAX must be greater than or equal to 0."
     		print*, "Input value is ", lmax
     		stop
     	elseif(abs(z) > 1.) then
     		print*, "Error --- PlmBar_d1"
     		print*, "ABS(Z) must be less than or equal to 1."
     		print*, "Input value is ", z
     		stop
     	elseif (abs(z) == 1.) then
     		print*, "Error --- PlmBar_d1"
     		print*, "Derivative can not be calculated at Z = 1 or -1."
     		print*, "Input value is ", z
     		stop
     	endif
     	
     	
     	if (present(csphase)) then
     		if (csphase == -1) then
     			phase = -1.
     		elseif (csphase == 1) then
     			phase = 1.
     		else
     			print*, "PlmBar_d1 --- Error"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)."
     			stop
     		endif
     	else
     		phase = dble(CSPHASE_DEFAULT)
     	endif
     		
	scalef = 1E-20


	if (lmax /= lmax_old) then
		
		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		
		allocate(sqr(2*lmax+1), stat=astat(1))
		allocate(f1((lmax+1)*(lmax+2)/2), stat=astat(2))
		allocate(f2((lmax+1)*(lmax+2)/2), stat=astat(3))
		
		if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
			print*, "PlmBar_d1 --- Error"
			print*, "Problem allocating arrays SQR, F1 and F2", astat(1), astat(2), astat(3)
			stop
		endif

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		!	Precompute square roots of integers that are used several times.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		do l=1, 2*lmax+1
			sqr(l) = sqrt(float(l))
		enddo

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Precompute multiplicative factors used in recursion relationships
		! 		Plmbar(l,m) = x*f1(l,m)*Plmbar(l-1,m) - Plmbar(l-2,m)*f2(l,m)
		!		k = l*(l+1)/2 + m + 1
		!	Note that prefactors are not used for the case when m=l and m=l-1,
		!	as a different recursion is used for these two values.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		k = 3
	
		do l=2, lmax, 1
			k = k + 1
			f1(k) = sqr(2*l-1) * sqr(2*l+1) / float(l)
			f2(k) = float(l-1) * sqr(2*l+1) / sqr(2*l-3) / float(l)
			do m=1, l-2
				k = k+1
				f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                		f2(k) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                  			 / sqr(2*l-3) / sqr(l+m) / sqr(l-m) 
			enddo
			k = k + 2
		enddo
		
		lmax_old = lmax
		
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!	
	!	Calculate P(l,0). These are not scaled.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	u = sqrt((1.-z)*(1.+z)) ! sin(theta)

      	pm2  = 1.	
      	p(1) = 1.
      	dp(1)= 0.
      	
      	if (lmax == 0) return
      	
      	pm1  = sqr(3)*z	
      	p(2) = pm1
      	dp(2) = sqr(3)
      		
	k = 2

      	do l = 2, lmax, 1
         	k = k+l
         	plm = f1(k)*z*pm1-f2(k)*pm2
         	p(k) = plm
         	dp(k) = float(l) * ( sqr(2*l+1) / sqr(2*l-1)  * &
      			pm1 - z * plm ) / u**2
         	pm2  = pm1
         	pm1  = plm
      	enddo

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate P(m,m), P(m+1,m), and P(l,m)
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
      	pmm  = sqr(2)*scalef
      	rescalem = 1.0/scalef
      	kstart = 1

      	do m = 1, lmax - 1, 1
      		
      		rescalem = rescalem*u

		! Calculate P(m,m)
        	kstart = kstart+m+1
         
         	pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
        	p(kstart) = pmm*rescalem
        	dp(kstart) = -float(m) * z * p(kstart) / u**2
        	pm2 = pmm

		! Calculate P(m+1,m)
		k = kstart+m+1
	   	pm1 = z * sqr(2*m+3) * pmm
	    	p(k) = pm1*rescalem
	    	dp(k) =  ( sqr(2*m+3) * p(k-m-1) - z * float(m+1) * p(k)) / u**2

		! Calculate P(l,m)
               	do l = m+2, lmax, 1
               		k = k+l
                  	plm  = z*f1(k)*pm1-f2(k)*pm2
                  	p(k) = plm*rescalem
                  	dp(k) = ( sqr(2*l+1) * sqr(l-m) * sqr(l+m) / sqr(2*l-1) * &
      				p(k-l) - z * float(l) * p(k)) / u**2
                  	pm2  = pm1
                  	pm1  = plm
               	enddo
              
      	enddo
      	
      	! Calculate P(lmax,lmax)
      
	rescalem = rescalem*u	
   	
        kstart = kstart+m+1
        pmm = phase * pmm * sqr(2*lmax+1) / sqr(2*lmax)
        p(kstart) = pmm*rescalem
        dp(kstart) = -float(lmax) * z * p(kstart) / u**2
      		
end subroutine PlmBar_d1
!
!
!
!
!
!
  FUNCTION FUNCLM (X)
!
! Computes the integrand needed to determine the shape factors for "quadrilateral" 
! disc elements. Uses the (modified) normalized associated Legendre functions at 
! degree 'l' and order 'm', which enter through a labeled common block. Here 'x' 
! is colatitude, in radians.  - Last modified GS 17/9/07
!
  IMPLICIT NONE
  INTEGER L, M 
  COMMON /AREA_DEG/L, M
  INTEGER, PARAMETER :: LLMAX=512, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
  INTEGER J_INDEX 
  REAL*8, PARAMETER :: PI=3.14159265358979323840     
  REAL*8 X, FUNCLM, PLM(JJMAX)			           
! 
! ---- Builds the Associated legendre functions up to harminic degree 'l' 
!      using SHTOOLS - The argument of the call is latitude (in degrees) 
!
  call plmbar_mod (l, 90.0-180.*x/pi, plm)
!
! ---- Extracts the harmonic of degree 'l' and order 'm' from Plm and scales
!      for sin(x), see Spada and Stocchi, Computers & Geosciences, (2007). 
!
  FUNCLM = sin(x)*plm(j_index(l,m)) 
!
  END FUNCTION FUNCLM 
!
!
!
!
!
!
      SUBROUTINE ROMINT (VAL, ERR, EPS, A, B, N, MAXE, F)
      IMPLICIT NONE
!
! ------------------------------------------------------------------------
! downloaded from:
! http://people.scs.fsu.edu/~burkardt/f77_src/toms351/toms351.html
! ------------------------------------------------------------------------
!
      REAL*8 A
      REAL*8 B
      REAL*8 BB
      REAL*8 EPS
      REAL*8 ERR
      EXTERNAL F
      REAL*8 F
      REAL*8 H
      INTEGER J
      INTEGER K
      INTEGER K0
      INTEGER K1
      INTEGER K2
      INTEGER KK
      INTEGER KKK
      INTEGER L
      INTEGER MAXE
      INTEGER N
      INTEGER N0
      INTEGER N1
      REAL*8 R
      REAL*8 RM(16)
      REAL*8 S
      REAL*8 S0
      REAL*8 S1
      REAL*8 T
      REAL*8 VAL
!
!  INITIAL TRAPEZOID RULE.
!
      T = ( B - A ) * ( F(A) + F(B) ) * 0.5
!
!  INITIAL RECTANGLE VALUE.
!
      RM(1) = ( B - A ) * F ( ( A + B ) * 0.5 )
!
      N = 2
      R = 4.0
!
      DO 11 K = 1, MAXE

      BB = ( R * 0.5 - 1.0 ) / ( R - 1.0 )
!
!  IMPROVED TRAPEZOID VALUE.
!
       T = RM(1) + BB * ( T - RM(1) )
!
!  DOUBLE NUMBER OF SUBDIVISIONS OF (A,B).
!
       N = 2 * N
       S = 0.0 
       H = ( B - A ) / FLOAT ( N )
!
!  CALCULATE RECTANGLE VALUE.
!
      IF ( N - 32 ) 1, 1, 2

1       N0 = N
        GO TO 4
2       N0 = 32
3       IF ( N - 512 ) 4, 4, 5
4       N1 = N
        GO TO 6	
5       N1 = 512
6       DO 9 K2 = 1, N, 512
          S1 = 0.0 
          KK = K2 + N1 - 1
          DO 8 K1 = K2, KK, 32
            S0 = 0.0 
            KKK = K1 + N0 - 1
            DO 7 K0 = K1, KKK, 2
              S0 = S0 + F ( A + FLOAT ( K0 ) * H )
7           CONTINUE
            S1 = S1 + S0
8         CONTINUE
          S = S + S1
9       CONTINUE
        RM(K+1) = 2.0 * H * S
!
!  END CALCULATION OF RECTANGLE VALUE.
!
       R = 4.0
!
!  FORM ROMBERG TABLE FROM RECTANGLE VALUES.
!
        DO 10 J = 1, K
          L = K + 1 - J
          RM(L) = RM(L+1) + ( RM(L+1) - RM(L) ) / ( R - 1.0  )
          R = 4.0  * R
10      CONTINUE

        ERR = ABS ( T - RM(1) ) * 0.5 
!
!  CONVERGENCE TEST.
!
        IF ( ERR - EPS ) 12, 12, 11

11    CONTINUE
!
      K = 0
!
12    VAL = ( T + RM(1) ) * 0.5 
      N = N + 1
      MAXE = K
!
      RETURN
      END
!
!
!

!
!
!
!
!
!								      
!				   ~~~ 				      
!				 --- --- 		              
!	 		-+-+-+-+-+-+-+-+-+-+-+-+-		      
! ----===========================================================---- 
!         The code which follows is entierely by Mark Tegmark ...     
! ----===========================================================---- 
! 			-+-+-+-+-+-+-+-+-+-+-+-+-		      
!				 --- --- 			      
!				   ~~~				      



!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     THESE SUBROUTINES ARE ALL YOU NEED TO CALL FROM OUTSIDE    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!	These subroutines convert between unit vectors and         !!!
!!!     pixel numbers, and are the only ones that the user of this !!!
!!!     package calls repeatedly:				   !!!
!!!	  subroutine vector2pixel(vector,resolution,R,v,pixel)     !!!
!!!	  subroutine pixel2vector(pixel,resolution,R,v,vector)     !!!
!!!                                                                !!!
!!!	These subroutines are called only once, in the beginning,  !!!
!!!     and compute the necessary rotation matrices and corner     !!!
!!!     vectors once and for all:                                  !!!
!!!	  subroutine compute_matrices(R)                           !!!
!!!	  subroutine compute_corners(v)                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	subroutine vector2pixel(vector,resolution,R,v,pixel)
	implicit none
	real*8 vector(3), R(0:19,3,3), v(0:11,3)
	real*8 A(3,3), vec(3), x, y
	integer resolution, pixel
	integer pix, face, pixperface, ifail
	if (resolution.lt.1) then
		write(*,*) "ERROR in the Tegmark modules: Resolution must exceed 0"
		STOP
	Endif
	pixperface = 2*resolution*(resolution-1)
	call find_face(vector,R,face)
	call getmatrix(face,R,A)
	call vecmatmul2(A,vector,vec)	
	x	= vec(1)/vec(3)
	y	= vec(2)/vec(3)
	call adjust(x,y)
	call tangentplanepixel(resolution,x,y,pix,ifail)
	if (ifail.gt.0) then 
	  ! Try the runner-up face:
	  call find_another_face(vector,R,face)
	  call getmatrix(face,R,A)
	  call vecmatmul2(A,vector,vec)	
	  x	= vec(1)/vec(3)
	  y	= vec(2)/vec(3)
	  call adjust(x,y)
	  call tangentplanepixel(resolution,x,y,pix,ifail)
	end if	
	pixel = face*pixperface + pix
	if (ifail.gt.0) then
	  ! The pixel wasn't in any of those two faces,
	  ! so it must be a corner pixel.
	  call find_corner(vector,v,pix)
	  pixel = 20*pixperface + pix
	end if
	return
	end

	subroutine pixel2vector(pixel,resolution,R,v,vector)
	! Returns a unit vector pointing towards pixel.
	! Resolution must be an even, positive integer.
	implicit none
	real*8 vector(3), R(0:19,3,3), v(0:11,3)
	real*8 A(3,3), x, y, norm
	integer resolution, pixel
	integer pix, face, pixperface
	if (resolution.lt.1) then 
		Write(*,*) "ERROR in the Tegmark modules: Resolution must exceed 0"
		STOP
	ENDIF
	pixperface = 2*resolution*(resolution-1)
	if (pixel.lt.0) then 
		Write(*,*)"ERROR in the Tegmark modules: negative pixel number"	
		STOP
	ENDIF
	if (pixel.ge.20*pixperface+12) then 
		Write(*,*)"ERROR in the Tegmark modules: pixel number too large"
		STOP 
	Endif
	if (pixperface.gt.0) then 
	  face = pixel/pixperface
	  if (face.gt.20) face = 20
	else ! There are no pixels at all on the faces - it's all just corners.
	  face = 20
	end if	
	pix = pixel - face*pixperface
	if (face.lt.20) then
	  ! The pixel is on one of the 20 faces:
	  call tangentplanevector(pix,resolution,x,y)
	  call unadjust(x,y)
	  norm 	= sqrt(x*x+y*y+1)
	  vector(1)	= x/norm
	  vector(2)	= y/norm
	  vector(3)	= 1./norm
	  call getmatrix(face,R,A)
	  call vecmatmul1(A,vector,vector)
	else
	  ! This is a corner pixel:
	  if (pix.gt.11) then 
	  	Write(*,*)"ERROR in the Tegmark modules: pixel number too big"
	  	STOP 
	  endif
	  vector(1) = v(pix,1)
	  vector(2) = v(pix,2)
	  vector(3) = v(pix,3)
	end if
	return
	end

	subroutine compute_matrices(R)
	! On exit, R will contain the 20 rotation matrices
	! that rotate the 20 icosahedron faces 
	! into the tangent plane
	! (the horizontal plane with z=1). 
	! Only called once, so speed is irrelevant.
	implicit none
	real*8 R(0:19,3,3)
	real*8 A(3,3), B(3,3), C(3,3), D(3,3), E(3,3)
	real*8 pi, sn, cs, ct, x
	integer i,j,n
	do i=1,3
	  do j=1,3
	    A(i,j) = 0.
	    B(i,j) = 0.
	    C(i,j) = 0.
	    D(i,j) = 0.
	  end do
	end do
	pi 	= 4.*atan(1.)
	x 	= 2.*pi/5.
	cs 	= cos(x)
	sn	= sin(x)
	A(1,1)	= cs
	A(1,2)	= -sn
	A(2,1)	= sn
	A(2,2)	= cs
	A(3,3) 	= 1.
	! A rotates by 72 degrees around the z-axis.
	x 		= pi/5.
	ct 		= cos(x)/sin(x)
	cs		= ct/sqrt(3.)
	sn		= sqrt(1-ct*ct/3.)
	C(1,1)	= 1
	C(2,2)	= cs
	C(2,3)	= -sn
	C(3,2)	= sn
	C(3,3) 	= cs
	! C rotates around the x-axis so that the north pole 	
	! ends up at the center of face 1.
	cs		= -0.5
	sn		= sqrt(3.)/2
	D(1,1)	= cs
	D(1,2)	= -sn
	D(2,1)	= sn
	D(2,2)	= cs
	D(3,3)	= 1.
	! D rotates by 120 degrees around z-axis.
	call matmul1(C,D,E)
	call matmul2(E,C,B)	! B = CDC^t
	! B rotates face 1 by 120 degrees. 
	do i=1,3
	  do j=1,3
	    E(i,j) = 0.
	  end do
	  E(i,i) = 1.
	end do	! Now E is the identity matrix.
	call putmatrix(0,R,E)	
	call matmul1(B,A,E)
	call matmul1(B,E,E)
	call putmatrix(5,R,E)	
	call matmul1(E,A,E)
	call putmatrix(10,R,E)
	call matmul1(E,B,E)
	call matmul1(E,B,E)
	call matmul1(E,A,E)
	call putmatrix(15,R,E)
	do n=0,15,5	
	  call getmatrix(n,R,E)
	  do i=1,4
	    call matmul1(A,E,E)
	    call putmatrix(n+i,R,E)
	  end do
	end do
	! Now the nth matrix in R will rotate 
	! face 1 into face n. 
	! Multiply by C so that they will rotate
	! the tangent plane into face n instead:
	do n=0,19
	  call getmatrix(n,R,E)
	  call matmul1(E,C,E)
	  call putmatrix(n,R,E)
	end do
	return
	end

	subroutine compute_corners(v)
	! On exit, v will contain unit vectors pointing toward 
	! the 12 icoshedron corners. 
	implicit none
	real*8 v(0:11,3)
	real*8 pi, z, rho, dphi
	integer i
	pi = 4.*atan(1.)
	dphi = 2.*pi/5.
	! First corner is at north pole:
	v(0,1) = 0.
	v(0,2) = 0.
	v(0,3) = 1.
	! The next five lie on a circle, with one on the y-axis:
	z = 0.447213595	! This is 1/(2 sin^2(pi/5)) - 1
	rho = sqrt(1.-z*z)
	do i=0,4
	  v(1+i,1) = -rho*sin(i*dphi)
	  v(1+i,2) =  rho*cos(i*dphi)
	  v(1+i,3) = z 
	end do
	! The 2nd half are simply opposite the first half:
	do i=0,5
	  v(6+i,1) = -v(i,1)
	  v(6+i,2) = -v(i,2)
	  v(6+i,3) = -v(i,3)
	end do
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     THE SUBROUTINES BELOW ARE SUBORDINATE TO THOSE ABOVE, AND  !!!
!!!     CAN BE SAFELY IGNORED BY THE GENERAL USER OF THE PACKAGE.  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!	These subroutines perform some standard linear algebra:    !!!
!!!	  subroutine matmul1(A,B,C)                                !!!
!!!	  subroutine matmul2(A,B,C)                                !!!	
!!!	  subroutine matmul3(A,B,C)                                !!!	
!!!	  subroutine vecmatmul1(A,b,c)	                           !!!	
!!!	  subroutine vecmatmul2(A,b,c)	                           !!!	
!!!                                                                !!!
!!!	These subroutines copy matrices in and out of storage:     !!!
!!!	  subroutine getmatrix(n,R,A)                              !!!
!!!	  subroutine putmatrix(n,R,A)                              !!!
!!!                                                                !!!
!!!     These subroutines help vector2pixel reduce the 3D sphere   !!!
!!!     problem to a problem on an equilateral triangle in the     !!!
!!!     z=1 tangent plane (an icosahedron face):                   !!!		
!!!	  subroutine find_face(vector,R,face)                      !!!
!!!	  subroutine find_another_face(vector,R,face)              !!!
!!!	  subroutine find_corner(vector,v,corner)                  !!!
!!!                                                                !!!
!!!     These subroutines pixelize this triangle with a regular    !!!
!!!     triangular grid:                                           !!!
!!!	  subroutine find_mn(pixel,resolution,m,n)                 !!!
!!!	  subroutine tangentplanepixel(resolution,x,y,pix,ifail)   !!!
!!!	  subroutine tangentplanevector(pix,resolution,x,y)        !!!
!!!                                                                !!!
!!!     These subroutines reduce the area equalization problem to  !!!
!!!     one on the right triangle in the lower right corner:       !!!
!!!	  subroutine find_sixth(x,y,rot,flip)                      !!!
!!!	  subroutine rotate_and_flip(rot,flip,x,y)	           !!!
!!!	  subroutine adjust(x,y)                                   !!!
!!!	  subroutine unadjust(x,y)                                 !!!
!!!                                                                !!!
!!!     These subroutines perform the area equalization mappings   !!!
!!!     on this right triangle:                                    !!!
!!!	  subroutine adjust_sixth(x,y)                             !!!
!!!	  subroutine unadjust_sixth(x,y)                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	subroutine matmul1(A,B,C)	
	! Matrix multiplication C = AB.
	! A, B and C are allowed to be physiclly the same.
	implicit none
	real*8 A(3,3), B(3,3), C(3,3), D(3,3), sum
	integer i,j,k
	sum = 0.
	do i=1,3
	  do j=1,3
	    sum = 0.
	    do k=1,3
	      sum = sum + A(i,k)*B(k,j)
	    end do
	    D(i,j) = sum
	  end do
	end do
	call copymatrix(D,C)
	return
	end

	subroutine matmul2(A,B,C)	
	! Matrix multiplication C = AB^t
	! A, B and C are allowed to be physically the same.
	implicit none
	real*8 A(3,3), B(3,3), C(3,3), D(3,3), sum
	integer i,j,k
	sum = 0.
	do i=1,3
	  do j=1,3
	    sum = 0.
	    do k=1,3
	      sum = sum + A(i,k)*B(j,k)
	    end do
	    D(i,j) = sum
	  end do
	end do
	call copymatrix(D,C)
	return
	end

	subroutine matmul3(A,B,C)	
	! Matrix multiplication C = A^t B
	! A, B and C are allowed to be physically the same.
	implicit none
	real*8 A(3,3), B(3,3), C(3,3), D(3,3), sum
	integer i,j,k
	sum = 0.
	do i=1,3
	  do j=1,3
	    sum = 0.
	    do k=1,3
	      sum = sum + A(k,i)*B(k,j)
	    end do
	    D(i,j) = sum
	  end do
	end do
	call copymatrix(D,C)
	return
	end

	subroutine vecmatmul1(A,b,c)	
	! Matrix multiplication c = Ab
	! b and c are allowed to be physically the same.
	implicit none
	real*8 A(3,3), b(3), c(3), d(3), sum
	integer i,j
	sum = 0.
	do i=1,3
	  sum = 0.
	  do j=1,3
	    sum = sum + A(i,j)*b(j)
	  end do
	  d(i) = sum
	end do
	call copyvector(d,c)
	return
	end

	subroutine vecmatmul2(A,b,c)	
	! Matrix multiplication c = A^tb
	! b and c are allowed to be physiclly the same.
	implicit none
	real*8 A(3,3), b(3), c(3), d(3), sum
	integer i,j
	sum = 0.
	do i=1,3
	  sum = 0.
	  do j=1,3
	    sum = sum + A(j,i)*b(j)
	  end do
	  d(i) = sum
	end do
	call copyvector(d,c)
	return
	end

	subroutine copymatrix(A,B)	
	! B = A
	implicit none
	real*8 A(3,3), B(3,3)
	integer i,j
	do i=1,3
	  do j=1,3
	    B(i,j) = A(i,j)
	  end do
	end do
	return
	end

	subroutine copyvector(a,b)	
	! b = a
	implicit none
	real*8 a(3), b(3)
	integer i
	do i=1,3
	  b(i) = a(i)
	end do
	return
	end

	subroutine getmatrix(n,R,A)	
	! A = the nth matrix in R
	implicit none
	real*8 R(0:19,3,3), A(3,3)
	integer i,j,n
	do i=1,3
	  do j=1,3
	    A(i,j) = R(n,i,j)
	  end do
	end do
	return
	end

	subroutine putmatrix(n,R,A)	
	! the nth matrix in R = A
	implicit none
	real*8 R(0:19,3,3), A(3,3)
	integer i,j,n
	do i=1,3
	  do j=1,3
	    R(n,i,j) = A(i,j)
	  end do
	end do
	return
	end

	subroutine find_face(vector,R,face)
	! Locates the face to which vector points.
	! Computes the dot product with the vectors
	! pointing to the center of each face and picks the
	! largest one. 
	! This simple routine can be substantially accelerated
	! by adding a bunch of if-statements, to avoid looping
	! over more than a few faces.
	implicit none
	real*8 vector(3), R(0:19,3,3), dot, max
	integer n,face,i
	max = -17.
	do n=0,19
	  dot = 0.
	  do i=1,3
	    dot = dot + R(n,i,3)*vector(i)
	  end do
	  if (dot.gt.max) then
	    face = n
	    max = dot
	  end if
	end do
	return
	end

	subroutine find_another_face(vector,R,face)
	! Computes the dot product with the vectors
	! pointing to the center of each face and picks the
	! largest one other than face.
	! This simple routine can be substantially accelerated
	! by adding a bunch of if-statements, to avoid looping
	! over more than a few faces.
	implicit none
	real*8 vector(3), R(0:19,3,3), dot, max
	integer n,face,facetoavoid,i
	facetoavoid = face
	max = -17.
	do n=0,19
	  if (n.ne.facetoavoid) then 
	    dot = 0.
	    do i=1,3
	      dot = dot + R(n,i,3)*vector(i)
	    end do
	    if (dot.gt.max) then
	      face = n
	      max = dot
	    end if
	  end if
	end do
	return
	end

	subroutine find_corner(vector,v,corner)
	! Locates the corner to which vector points.
	! Computes the dot product with the vectors
	! pointing to each corner and picks the
	! largest one. 
	! This simple routine can be substantially accelerated
	! by adding a bunch of if-statements, but that's pretty
	! pointless since it gets called so rarely.
	implicit none
	real*8 vector(3), v(0:11,3), dot, max
	integer corner,n,i
	max = -17.
	do n=0,11
	  dot = 0.
	  do i=1,3
	    dot = dot + v(n,i)*vector(i)
	  end do
	  if (dot.gt.max) then
	    corner = n
	    max = dot
	  end if
	end do
	return
	end

	subroutine find_mn(pixel,resolution,m,n)
	! Computes the integer coordinates (m,n) of the pixel 
	! numbered pix on the basic triangle.
	implicit none
	integer pixel, pix, resolution, m, n
	integer interiorpix , pixperedge
	pix 	    = pixel
	interiorpix = (2*resolution-3)*(resolution-1)
	pixperedge  = (resolution)-1
	if (pix.lt.interiorpix) then 
	  ! The pixel lies in the interior of the triangle.
	  m = (sqrt(1.+8.*pix)-1.)/2. + 0.5/resolution
	  ! 0.5/resolution was added to avoid problems with
	  ! rounding errors for the case when n=0.
	  ! As long as you don't add more than 2/m, you're OK.
	  n = pix - m*(m+1)/2
	  m = m + 2
	  n = n + 1
	  goto 555
	end if
	pix = pix - interiorpix 
	if (pix.lt.pixperedge) then
	  ! The pixel lies on the bottom edge.
	  m = 2*resolution-1
	  n = pix+1
	  goto 555
	end if
	pix = pix - pixperedge
	if (pix.lt.pixperedge) then
	  ! The pixel lies on the right edge.
	  m = 2*resolution-(pix+2)
	  n = m
	  goto 555
	end if
	pix = pix - pixperedge
	! The pixel lies on the left edge.
	m = pix+1
	n = 0
555	return
	end

	subroutine tangentplanepixel(resolution,x,y,pix,ifail)
	! Finds the hexagon in which the point (x,y) lies
	! and computes the corresponding pixel number pix.
	! Returns ifail=0 if (x,y) lies on the face, 
	! otherwise returns ifail=1.
	implicit none
	real*8 x, y, a, b, c, d, edgelength
	parameter (c=0.866025404)	! sqrt(3)/2
	parameter(edgelength=1.3231690765)
	! The edge length of the icosahedron is 
	! sqrt(9 tan^2(pi/5) - 3) when scaled so that
	! it circumscribes the unit sphere. 
	integer resolution, pix, ifail, i, j, k, m, n, r2
	r2	= 2*resolution
	a 	= 0.5*x
	b 	= c*y
	d 	= 0.5*edgelength/r2
	i 	= x/d 	  + r2
	j 	= (a+b)/d + r2
	k 	= (a-b)/d + r2
	m 	= (r2+r2-j+k-1)/3
	n 	= (i+k+1-r2)/3
	pix 	= (m-2)*(m-1)/2 + (n-1)
	ifail = 0
	if (m.eq.r2-1) then		! On bottom row
	  if ((n.le.0).or.(n.ge.resolution)) then
	    ifail=1
	  end if
	  goto 666				! Pix already correct
	end if
	if (n.eq.m) then			! On right edge
	  k = (r2-1) - m
	  if ((k.le.0).or.(k.ge.resolution)) then
	    ifail = 1
	  else
	    pix = (r2-2)*(resolution-1) + k - 1
 	  end if
	  goto 666
	end if
	if (n.eq.0) then			! On left edge
	  if ((m.le.0).or.(m.ge.resolution)) then
	    ifail = 1
	  else
	    pix = (r2-1)*(resolution-1) + m - 1
	  end if
	end if
666	return
	end

	subroutine tangentplanevector(pix,resolution,x,y)
	! Computes the coordinates (x,y) of the pixel 
	! numbered pix on the basic triangle.
	implicit none
	real*8 x, y, c1, c2, edgelength
	parameter(c1=0.577350269)	! 1/sqrt(3)
	parameter(c2=0.866025404)	! sqrt(3)/2
	parameter(edgelength=1.3231690765)
	! The edge length of the icosahedron is 
	! sqrt(9 tan^2(pi/5) - 3) when scaled so that
	! it circumscribes the unit sphere. 
	integer pix, resolution, m, n
	call find_mn(pix,resolution,m,n)
	x	= edgelength*(n-0.5*m)/(2*resolution-1)
	y 	= edgelength*(c1-(c2/(2*resolution-1))*m)
	return
	end

	subroutine find_sixth(x,y,rot,flip)
	! Find out in which sixth of the basic triangle
	! the point (x,y) lies, identified by the
	! two integers rot (=0, 1 or 2) and flip = 0 or 1).
	! rot and flip are defined such that the sixth is 
	! mapped onto the one at the bottom right by
	! these two steps:
	! 1. Rotate by 120 degrees anti-clockwise, rot times.
	! 2. Flip the sign of x if flip = 1, not if flip=0.
	! The if-statements below go through the six cases 
	! anti-clockwise, starting at the bottom right.
	implicit none
	real*8 x, y, c, d
	parameter(c=1.73205081)		! sqrt(3)
	integer rot, flip
	d = c*y
	if (x.ge.0) then
	  if (x.le.-d) then
	    rot  = 0
	    flip = 0
	  else
	    if (x.ge.d) then
	      rot  = 2
	      flip = 1
	    else 
	      rot  = 2
	      flip = 0
	    end if
	  end if
	else
	  if (x.ge.-d) then
	    rot  = 1
	    flip = 1
	  else
	    if (x.le.d) then
	      rot  = 1
	      flip = 0
	    else 
	      rot  = 0
	      flip = 1
	    end if
	  end if
	end if
	return
	end

	subroutine rotate_and_flip(rot,flip,x,y)
	implicit none
	real*8 x, y, x1, cs, sn, c
	integer rot, flip
	parameter(cs=-0.5)
	parameter(c=0.866025404)	! sqrt(3)/2
	if (rot.gt.0) then
	  if (rot.eq.1) then
	    sn = c	! Rotate 120 degrees anti-clockwise 
	  else 
	    sn = -c	! Rotate 120 degrees anti-clockwise 
	  end if
	  x1 = x
	  x  = cs*x1 - sn*y
	  y  = sn*x1 + cs*y
	end if
	if (flip.gt.0) x = -x
	return
	end

	subroutine adjust(x,y)
	! Maps the basic triangle onto itself in such a way
	! that pixels will have equal area when mapped onto
	! the sphere. 
	implicit none
	real*8 x, y
	integer rot, flip
	call find_sixth(x,y,rot,flip)
	call rotate_and_flip(rot,flip,x,y)
	call adjust_sixth(x,y)
	! Now rotate & flip the sixth back into its
	! original position:
	if ((flip.eq.0).and.(rot.gt.0)) then  
	  call rotate_and_flip(3-rot,flip,x,y)
	else
	  call rotate_and_flip(rot,flip,x,y)
	end if
	return
	end
	
	subroutine unadjust(x,y)
	! Performs the inverse of what adjust does. 
	implicit none
	real*8 x, y
	integer rot, flip
	call find_sixth(x,y,rot,flip)
	call rotate_and_flip(rot,flip,x,y)
	call unadjust_sixth(x,y)
	! Now rotate & flip the sixth back into its
	! original position:
	if ((flip.eq.0).and.(rot.gt.0)) then  
	  call rotate_and_flip(3-rot,flip,x,y)
	else
	  call rotate_and_flip(rot,flip,x,y)
	end if
	return
	end
	
	subroutine adjust_sixth(x,y)
	! Maps the basic right triangle (the sixth of the face that
	! is in the lower right corner) onto itself in such a way
	! that pixels will have equal area when mapped onto the sphere. 
	implicit none
	real*8 x, y, u, v, g, v2, root, trig, eps, scale
	parameter(eps=1.e-14, scale=1.09844)
	parameter(g=1.7320508075689)		! sqrt(3)
	u 		= x  + eps
	v		= -y + eps
	v2		= v*v
	root		= sqrt(1.+4.*v2)
	trig		= atan((g*root-g)/(root+3.))
	y		= sqrt(trig*2./g)
	x		= sqrt((1.+4.*v2)/(1.+u*u+v2))*u*y/v
	x		= scale*x
	y 		= -scale*y	
	return
	end
	
	subroutine unadjust_sixth(x,y)
	! Performs the inverse of what adjust_sixth does.
	implicit none
	real*8 x, y, u, v, g, v2, y2, tmp, trig, eps, scale
	parameter(eps=1.e-14, scale=1.09844)
	parameter(g=1.7320508075689)		! sqrt(3)
	u 		=  x/scale + eps
	v		= -y/scale + eps
	v2		= v*v
	trig		= tan(g*v2/2.)
	tmp		= (g+3.*trig)/(g-trig)
	y2		= (tmp*tmp-1.)/4.
	y		= sqrt(y2)
	tmp		= v2*(1.+4.*y2) - u*u*y2
	x	 	= u*y*sqrt((1.+y2)/tmp) 
	y 		= -y	
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     END OF THE ICOSAHEDRON PACKAGE FOR PIXELIZING THE SPHERE   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!
!
!
!
!
!
  SUBROUTINE QSORT_MOD(N, ARR, BRR)
  implicit none
!
! ### A QuickSort routine 
!
      INTEGER n,M,NSTACK
      real*8 arr(n), brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      real*8 a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          temp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK) stop 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 1>.
!
!
!
!
!
!
!
  SUBROUTINE BUBBLE_SORT_PAR(N, A, B)
  implicit none
!
! ### A (parallel) bubble-sorting routine 
!
  integer i, j, n
  real*8 a(n), b(n), z 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,z)
  do i=1,n  
      if( mod(i,2)==1 ) then
!$OMP DO
         do j=1,n/2
	    if( a(2*j-1) > a(2*j) ) then
	       z = a(2*j-1)
               a(2*j-1) = a(2*j)
               a(2*j) = z
	       z = b(2*j-1)
               b(2*j-1) = b(2*j)
               b(2*j) = z
            end if
	 end do
      else
!$OMP DO
         do j=1,n/2 + mod(n,2) - 1
	    if( a(2*j) > a(2*j+1) ) then
	       z = a(2*j)
               a(2*j) = a(2*j+1)
               a(2*j+1) = z
	       z = b(2*j)
               b(2*j) = b(2*j+1)
               b(2*j+1) = z	       
	    end if
	 end do
      end if
  end do
!$OMP END PARALLEL
  
  return
  end
!
!
!
!
!
!
  SUBROUTINE BUBBLE_SORT_MOD(N, A, B)
  implicit none
!
! ### A bubble-sorting routine 
!
  integer i, j, n
  real*8 a(n), b(n)
  do i = 1, n-1
    do j = i+1, n
      if ( a(i) > a(j) ) then
        call swap_mod(a(i), a(j))
	 call swap_mod(b(i), b(j))
      end if
    end do
  end do
  return
  end
!
!
!
!
!
!
  SUBROUTINE SWAP_MOD(X,Y)
  implicit none
!
! ### A swapping routine useful for sorting 
!
  real*8 x, y, z 
  z = x
  x = y
  y = z
  return
  end
!
!
!
!
!
!
	Subroutine FINDPX(RES,NP,LONP,LATP)
	implicit NONE 
!
! # Given a Resolution RES, this routine returns the number of pixels NP,
!   and the longitude and latitude arrays LONP, LATP (1:NP) with the pixels
!   coordinates. Written on October 2, 2008 for the implementation of the 
!   3D velocity maps for SELEN 2.7 - but perhaps useful also for other 
!   purposes... 
!
	INTEGER J, NP, RES
        REAL*8 R(0:19,3,3), V(0:11,3), VECT(3), X,Y,Z  
	REAL*8 LONP(1:NP), LATP(1:NP), LONPX, LATPX  
	REAL*8,  PARAMETER :: PI=3.14159265358979323840 
!
        call compute_matrices (r)
        call compute_corners  (v)
!
	np=40*res*(res-1)+12 
!
        do 1 j=0,np-1
	  	call pixel2vector (j,res,r,v,vect)
	        x=vect(1)
		y=vect(2)
		z=vect(3)
! 
! --- polar pixel 
          	if(x==0..and.y==0.) then 
          		lonpx=0.
                	if(z>=0.) latpx = +90.
                	if(z<=0.) latpx = -90.
! --- ordinary pixel 
	  		else                         	  
	  		lonpx = atan2(y,x)*180./pi 
          		if (lonpx < 0.) lonpx=360.+ lonpx
          		latpx = 90.-acos(z)*180./pi 
          	endif
	 lonp(j+1)=lonpx 
	 latp(j+1)=latpx 	
!
1       continue 
!
	end subroutine FINDPX 
!
!
!
!
!
!
 SUBROUTINE AVERAGE_EARTH_DENSITY(IND, RHO_EARTH) 
 IMPLICIT NONE 
! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Filters the Average Earth density from the TABOO (1) or ALMA (2) Log files - !
! GS July 27, 2009. 			                                       !
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
 INTEGER, PARAMETER :: A_LARGE_INTEGER=10000
 CHARACTER*80 LINE 
 CHARACTER*5  JUNK5
 CHARACTER*40 JUNK40
 INTEGER I, NL, IND  
 REAL*8 RHO_EARTH  
 CHARACTER*30 LOGNAME 
!
!
  If(ind==1)     then 
  	LOGNAME='taboo.log'
  elseif(ind==2) then
  	LOGNAME='alma-logfile.dat' 
  endif
  
  Open(103,file=LOGNAME,status='unknown') 
!  
  NL=0
  do 20 I=1, A_LARGE_INTEGER 
	read(103,'(a80)') line 
	if(line(1:5)==">>>>>")then 
		NL=I  	
		goto 30 
	endif	
20 continue  
30 close(103)
!
  Open(103,file=LOGNAME,status='unknown') 
	do 40 i=1, nl-1 
	read(103,'(a80)') line 
40 continue 
  READ(103,'(A5,A40,E14.8)') JUNK5, JUNK40, RHO_EARTH   
  close(103)    
!    
 END SUBROUTINE Average_Earth_density   
!
!
!
!
!
!
! 
 Subroutine FROMDDMMSS_2_DEGREES (SS, FLOAT_DEGREES)   
 IMPLICIT NONE 
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Converts angular data with the format degrees-minutes-seconds 
! into a decimal format. The input is a CHARACTER*10 string with
! format DDDsMMsSSY, where:
!
! - DDD gives degrees (e.g., 056 or 112) 
! - MM  gives minutes (e.g., 11) 
! - SS  gives seconds (e.g., 59)   
! - Y stands for one of E, W, N, S 
! - s  is a separator (":" is a good choice)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ************************************
! Warning: BOUNDS are NOT checked!!!
! ************************************
!
  CHARACTER*10 SS
  CHARACTER*1  LAB 
  REAL*8 FLOAT_DEGREES, X_MINUTES, X_SECONDS
!
! write(*,*) "********", i, ss
!
  open (30,file='junk-deg.dat',status='unknown') 
  if    (ss(1:1)=='0')  then 
  	write(30,'(a2)') ss(2:3) 
  elseif(ss(1:2)=='00') then 	
  	write(30,'(a1)') ss(3:3)		
  elseif(ss(1:1)/='0') then 
  	write(30,'(a3)') ss(1:3)     
  endif
  close(30) 
  open (30,file='junk-deg.dat',status='unknown') 
  read (30,*) float_degrees  
  close(30) 
!
  open (30,file='junk-deg.dat',status='unknown') 
  	if(ss(5:5)=='0') write(30,'(a1)') ss(6:6) 
  	if(ss(5:5)/='0') write(30,'(a2)') ss(5:6) 
  close(30) 
  open (30,file='junk-deg.dat',status='unknown') 
  read (30,*) x_minutes   
  close(30)
  float_degrees = float_degrees + x_minutes/60. 
!
  open (30,file='junk-deg.dat',status='unknown') 
  	if(ss(8:8)=='0') write(30,'(a1)') ss(9:9) 
  	if(ss(8:8)/='0') write(30,'(a2)') ss(8:9) 
  close(30) 
  open (30,file='junk-deg.dat',status='unknown') 
  read (30,*) x_seconds   
  close(30)
  float_degrees = float_degrees + x_seconds/60./60. 
!
  lab=ss(10:10) 
!
  if    (lab=='N'.or.lab=='n') then
 	float_degrees =      + float_degrees
  elseif(lab=='S'.or.lab=='s') then 
        float_degrees =      - float_degrees  
  elseif(lab=='E'.or.lab=='e') then   
	float_degrees =      + float_degrees 
  elseif(lab=='W'.or.lab=='w') then
        float_degrees = 360. - float_degrees   
  Endif
! 
 End Subroutine FromDDMMSS_2_degrees
!
!
!
!
!
!
   SUBROUTINE ICE_BREAKER (NDATA, SLOTS, LOW, HIG)
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
!  - Determines the number of SLOTS of length SIZE=1000 in the range 1:NDATA.  
!
!  - If NDATA is the number of ice elements in some ice input file, SLOTS is 
!    then the number of files in which the ice model is "broken". Each piece
!    contains SIZE=1000 elements but the last one may contain less data...  
!
!  - Written by GS May 21 2010 for G95 implementation of SELEN. 
!           This is a patch for the case of large harmonic degree to avoid
!           large memory requests at run time. 
!         
! ---------------------------------------------------------------------------
!
   IMPLICIT NONE 
   INTEGER, PARAMETER :: SIZE  = 1000
   INTEGER, PARAMETER :: NMAX  = 3E4
!
   INTEGER S, NDATA, SLOTS, SMIN, SMAX, LOW(1), HIG(1)
!         
! ---------------------------------------------------------------------------
!
   If(NDATA.GT.NMAX)  then 
   	Write(*,*) "ICE_BREAKER.F90: The ICE data file is too large - Maximum size:", NMAX 
	Write(*,*) "ICE_BREAKER.F90: The program will STOP"
	STOP
   Endif
!
   slots=-999
!
   if (ndata.le.size) slots= 1 
   if (ndata.gt.size) slots= ndata/size+1
!
   smin = 1
   smax = slots
!
   if(ndata.le.size) then 
  	   low (smin)=1 
  	   hig (smin)=NDATA		   
   endif
!
   if(ndata.gt.size) then  
!
   do 10 s=1, smax-1
  	   low  (s)=size*(s-1)+1
  	   hig  (s)=size*(s)
10 continue   
!
   low(smax)=size*(s-1)+1
   hig(smax)=NDATA
!
   endif
!
   If(SLOTS==-999)  then 
   	Write(*,*) "ICE_BREAKER.F90: Some unknown error occurred"
	Write(*,*) "ICE_BREAKER.F90: The program will STOP"
	STOP
   Endif
!
   END SUBROUTINE ICE_BREAKER
!
!
!
!
!
!
