 !
 !
 !
 ! File <<data.inc>>, created by <<config.f90>> on 2012.10.30 14.41.46
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
 ! SELEN is distributed in the /hope that it will be useful, but WITHOUT ANY  
 ! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  
 ! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      
 ! details. 									
 !   										
 ! You should have received a copy of the GNU General Public License along    
 ! with SELEN.  If not, see <http://www.gnu.org/licenses/>.                   
 ! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 !
 !------------------- 1- General settings --------------------------
 !
 ! --- Repository name
 CHARACTER*4, PARAMETER :: RUN='TEST'
 ! --- Prefix for PX and SH files
 CHARACTER*2, PARAMETER :: TMP_PREFIX='./'
 ! --- Tegmark resolution
 INTEGER, PARAMETER :: RES=44        
 ! --- Number of pixels
 INTEGER, PARAMETER :: NP=2*RES*(RES-1)*20+12
 ! --- Maximum harmonic degree
 INTEGER, PARAMETER :: LMAX=128       
 ! --- Dealing with the degree 1 Love numbers
 INTEGER, PARAMETER :: DEG1= 1
 ! --- Jmax index
 INTEGER, PARAMETER :: JMAX=(LMAX+1)*(LMAX+2)/2
 ! --- Pi
 REAL*8,  PARAMETER :: PI=3.14159265358979323840
 ! --- Earth, water, and ice densities (kg/m^3)
 REAL*8,  PARAMETER :: RHOW=    1000.0000000000000     
 REAL*8,  PARAMETER :: RHOI=    931.00000000000000     
 !
 !
 !------------------- 2- Ice models settings -----------------------
 !
 ! --- Ice sheet model file
 CHARACTER*           9 , PARAMETER :: ICE_MODEL='ice5g.dat'         
 ! --- Number of time steps
 INTEGER, PARAMETER :: NN=21 
 ! --- Number of ice elements
 INTEGER, PARAMETER :: NEL=       11388
 ! --- Time increment, ka
 REAL*8, PARAMETER :: DELTA=   1.0000000000000000     
 !
 !
 !
 !------------------- 3- Normal modes (TABOO) settings -------------
 !
 ! --- Earth model (see TABOO user guide)
 INTEGER, PARAMETER :: NV=3  , CDE=2  
 ! --- Viscosity & Litho model file
 CHARACTER*          18 , PARAMETER :: VISCO_MODEL= &
 './VSC/vsc_VM2r.dat'
 !
 !
 !------------------- 4- Sea level Equation settings ---------------
 !
 ! --- Number of iterations
 INTEGER, PARAMETER :: SMAX=3
 ! --- Mode of solution
 INTEGER, PARAMETER :: IMODE=1
 !
 !
 !------------------- 5- RSL analysis at site of a database --------
 !
 ! --- RSL database
 CHARACTER*          19 , PARAMETER :: RSL_DATABASE='./DATA/sealevel.dat'         
 ! --- Format of the RSL database
 CHARACTER* 1, PARAMETER :: RSL_DATABASE_FORMAT = '0'
 ! --- Number of RSL sites in database
 INTEGER, PARAMETER :: NRSL=         392
 !
 !
 !------------------- 7- Predictions at tide-gauges ----------------
 !
 ! --- Tide-gauge database
 CHARACTER*          21 , PARAMETER :: TGAUGES_DATABASE='./DATA/rlr-trends.txt'       
 ! --- Number of tide gauges
 INTEGER, PARAMETER :: NTIDEGAUGES=        1123
 !
 !
 !------------------ 12- Stokes coefficients (SC) settings ---------
 !
 ! --- Min. degree for the SC
 INTEGER, PARAMETER :: STMIN=2         
 ! --- Max. degree for the SC
 INTEGER, PARAMETER :: STMAX=20        
 !
 !
 !
 !
 !
 !
