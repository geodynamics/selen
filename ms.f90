!
! This is program "MS.F90" 
!
! Created by GS July 12 2008  -Intel version 2.6- 
! Modified a number of times during July 2008 GS
! Also modified on July 24, 2008 for "DISK" ice file
! Updated on August, 2008 for "ANU" ice file (2.7)
! Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Reviewed DM June 2011 - Dynamic memory allocation
! *** Reviewed DM Jan 2012 - ALASKA ice model
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
! ---------------------------------------------------------------------------
!  MS.F90 prepares multi-segmet files of the ice model for plotting purposes. 
!  This routine substitutes all the SH*_P.F90 program units that were written 
!  "ad-hoc" for individual ice models until v. 2.5 
! ---------------------------------------------------------------------------
!
!
!
!       INCLUDE "harmonics.f90"
       PROGRAM MULTI_SEGMENT 
       IMPLICIT NONE        
       INCLUDE "data.inc"
!
!#---- General declarations 
       INTEGER, PARAMETER :: NH1=20, NH3=20, NH5=24 
       INTEGER, PARAMETER :: NHA=4,  NHI=18, NHD=20
       INTEGER, PARAMETER :: NHU=26, NHC=20
       REAL*8, PARAMETER :: AMP5 = 1.0, AMP1 = 5.     !AMP5 = 0.70313
       REAL*8, PARAMETER :: RAD=6.371E6   
       REAL*8, PARAMETER :: HMIN=1.0    
       INTEGER, PARAMETER :: ILARGE=1000        
       INTEGER I, J, K, L, N, IJ, NH, I1, I2
       INTEGER, ALLOCATABLE :: ICR(:)
       REAL*8, ALLOCATABLE :: TETAC(:), LONGC(:), LATIC(:), ALFA(:) 
       REAL*8 AJ, AMP
       REAL*8, ALLOCATABLE :: H(:,:), CR(:)       
       CHARACTER*2 CJ
       REAL*8 COSDD, SINDD

!
!#---- Declarations for "ICE3G"
       INTEGER, PARAMETER :: DSTEP=10  
       REAL*8 XXX, OMEGA, COSTETA, SINTETA, TETA, LON, LAT

!
!
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!      - Input files:  
!        Ice_model (from 'data.inc') 


!
!      - Output files: 
         CHARACTER*20 FILEOUT 
!	  ...    
!#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!
!
!#---- Allocate memory space
!	    	
    allocate( ICR(0:NN) )
    allocate( TETAC(NEL), LONGC(NEL), LATIC(NEL), ALFA(NEL) )
    allocate( H(NEL,0:NN+1), CR(0:NN) ) 
!
!#---- Reading the header lines of the ice table
!	    	
  OPEN(10,FILE=ICE_MODEL,STATUS='UNKNOWN') 
!
   if(ice_model(1:4)=='ice1')nh=nh1 
   if(ice_model(1:4)=='ice3')nh=nh3
   if(ice_model(1:4)=='ice5')nh=nh5
   if(ice_model(1:4)=='alps')nh=nha
   if(ice_model(1:4)=='alas')nh=nha
   if(ice_model(1:4)=='ij05')nh=nhi
   if(ice_model(1:4)=='disk')nh=nhd
   if(ice_model(1:4)=='anu0')nh=nhu
   if(ice_model(1:4)=='icap')nh=nhc
!		  
   do j=1, nh 
       read(10,'(a20)') cj
   enddo	     
!
!	    
!#---- Reading the time-dependent thickness of the ice elements 
!
!	    	
 DO 10 I=1, NEL   
!
!	    
!#---- "ICE1" ----
!             	    
 IF    (ICE_MODEL(1:4)=='ice1') THEN 
  Read (10,*) i1, i2, (icr(k),k=nn,0,-1) 
  cr(:)=icr(:)	
  tetac(i)=float(i1) 
  longc(i)=float(i2) 		
!
!
!#---- "ICE3" ----   
!
 ELSEIF(ICE_MODEL(1:4)=='ice3') THEN 
  Read (10,111) ij, tetac(i), longc(i), alfa(i), aj, (cr(k),k=nn,11,-1)
  Read (10,112)                                      (cr(k),k=10, 0,-1)		   		
  latic(i)=90.-tetac(i) 
!
!
!#---- "DISK" ----   
!
 ELSEIF(ICE_MODEL(1:4)=='disk') THEN 
  Read (10,*) ij, tetac(i), longc(i), alfa(i), (cr(k),k=nn,11,-1)
  Read (10,*)                                  (cr(k),k=10, 0,-1)		   		
  latic(i)=90.-tetac(i) 
!
!
!#---- "ICAP" ----   
!
 ELSEIF(ICE_MODEL(1:4)=='icap') THEN 
  Read (10,*) ij, tetac(i), longc(i), alfa(i), (cr(k),k=nn,11,-1)
  Read (10,*)                                  (cr(k),k=10, 0,-1)		   		
  latic(i)=90.-tetac(i) 
!
!
!#---- "ICE5" ----
!
 ELSEIF(ICE_MODEL(1:4)=='ice5') THEN 
  If (ice_model/='ice5g26.dat') &  
    Read (10,113) ij, longc(i), ij, latic(i), (cr(k),k=nn,0,-1) 			
  If (ice_model=='ice5g26.dat') &  
    Read (10,114) ij, longc(i), ij, latic(i), (cr(k),k=nn,0,-1) 
  tetac(i)=90.-latic(i)
!
!		        		
!#---- "ALPS" ----
!
 ELSEIF(ICE_MODEL(1:4)=='alps'.or.ICE_MODEL(1:4)=='alas') THEN 
    Read (10,*) ij, longc(i), latic(i), alfa(i), (cr(k),k=nn,0,-1) 
    if( ice_model(1:4)=='alps' ) alfa(i)=0.09
    tetac(i)=90.-latic(i)
!
!
!#---- "IJ05" ----
!
 ELSEIF(ICE_MODEL(1:4)=='ij05') THEN 
    Read (10,*) ij, longc(i), latic(i), alfa(i), aj, (cr(k),k=nn,9,-1) 
    Read (10,*) 	                             (cr(k),k= 8,0,-1)  
    tetac(i)=90.-latic(i) 
!
!
!#---- "ANU05" ----
!  		    	                      
 ELSEIF(ICE_MODEL(1:4)=='anu0') THEN 

    Read (10,115) ij, longc(i), latic(i), alfa(i), (h(i,k),k=0,nn+1) 
    tetac(i)=90.-latic(i) 
!
!#---- "End of available ice models..." 
!
 ENDIF
!
!
!#---- Conversion to our default for ice thickness
!
 If(ICE_MODEL(1:4)/='anu0') then 
 	h(i,0)=cr(nn)              
      		do k=1,nn                   
         		h(i,k)=cr(nn-k)   
      		enddo 
 	h(i,nn+1) = h(i,nn) 
 Endif 
 
!
!
10 CONTINUE
!
   close(10) 
!
!
   Write(*,*)'    - Read ',nel, '  elements from ', ice_model
!
!
!
!
!#---- Multi-segmentation 
!
	DO 20 K=0,NN+1
!
	    open(3,file='junk.dat',status='unknown') 
	    if(k<=9) write(3,'(a1,i1)') '0',k  
	    if(k> 9) write(3,'(i2)')        k
 	    close(3)
!	
	    open(3,file='junk.dat',status='unknown')
	    read(3,'(a2)') cj
	    close(3)
!
            if(k==0.or.k==nn+1) & 
	    WRITE(*,*) '    - Time ', cj, ' of', nn+1  	     
!
!	    
!#---- "ICE1" or "ICE5 "----
!
           IF    (ICE_MODEL(1:4)=='ice1'.or.ICE_MODEL(1:4)=='ice5') THEN 
!
	   if(ice_model(1:4)=='ice1') then 
	   			      fileout='msg1-'//cj//'.dat'
				      amp=amp1 
	   endif
	   if(ice_model(1:4)=='ice5') then 
	   			      fileout='msg5-'//cj//'.dat'
				      amp=amp5 
	   endif	   	   				    
!
	   open(7,file=fileout,status='unknown') 
       	   do 30 i=1, nel
           if(h(i,k).ge.hmin) then        
       			write(7,'(a4,1x,F10.4)') '> -Z', h(i,k)
       			write(7,*) longc(i)-amp/2., 90.-(tetac(i)-amp/2.)   
       			write(7,*) longc(i)-amp/2., 90.-(tetac(i)+amp/2.)   
       			write(7,*) longc(i)+amp/2., 90.-(tetac(i)+amp/2.)
       			write(7,*) longc(i)+amp/2., 90.-(tetac(i)-amp/2.)
	   endif 
30	   continue
           close(7) 
!
!
!#---- "ICE3" or "ALPS" or "ALASKA" or "IJ05" or "DISK" or "ANU05" ----   
!
           ELSEIF(ICE_MODEL(1:4)=='ice3'.or.&
	          ICE_MODEL(1:4)=='alps'.or.&
	          ICE_MODEL(1:4)=='alas'.or.&
		  ICE_MODEL(1:4)=='ij05'.or.&
		  ICE_MODEL(1:4)=='anu0'.or.&
		  ICE_MODEL(1:4)=='disk'.or.&
		  ICE_MODEL(1:4)=='icap')      THEN 
!	   
	   if(ice_model(1:4)=='ice3') fileout='msg3-'//cj//'.dat'
	   if(ice_model(1:4)=='alps') fileout='msgA-'//cj//'.dat'
	   if(ice_model(1:4)=='alas') fileout='msgK-'//cj//'.dat'
	   if(ice_model(1:4)=='ij05') fileout='msgI-'//cj//'.dat'
	   if(ice_model(1:4)=='disk') fileout='msgD-'//cj//'.dat'
	   if(ice_model(1:4)=='anu0') fileout='msgU-'//cj//'.dat'
	   if(ice_model(1:4)=='icap') fileout='msgC-'//cj//'.dat'		   	   
!
	   open(7,file=fileout,status='unknown') 
!
       	   do 40 i=1, nel  
           if(h(i,k).ge.hmin) then        
!
!--- A general element		
 	   if(tetac(i).ne.180.and.tetac(i).ne.0.) then 
       		write(7,'(a4,1x,F10.4)') '> -Z', h(i,k) 
                do j=1, 360, dstep  
       	             omega=float(j)       
		     costeta=cosdd(tetac(i))*cosdd(alfa(i)) + & 
		             sindd(tetac(i))*sindd(alfa(i))*cosdd(omega)       
       		     xxx=sindd(alfa(i))*sindd(omega)*sindd(tetac(i))/&
		           (cosdd(alfa(i))-costeta*cosdd(tetac(i))) 		              
       		     lon=atan(xxx)*180./pi + longc(i)       
       		     sinteta=sindd(alfa(i))*sindd(omega)/sindd(lon-longc(i))       
       		     teta=atan2(sinteta,costeta)*180./pi 		      
		     lat=90.-teta 		   
		     write(7,*) lon, lat, h(i,k)        
       		enddo
	   endif	
!		
!--- A south pole element		
 	   if(tetac(i).eq.180.) then 
	        write(7,'(a4,1x,F10.4)') '> -Z', h(i,k)  
   		write(7,'(a7)')      "0   -90"
		write(7,'(a7,f14.8)')"0      ", -90.+alfa(i)		 
      		do j=1, 360, dstep  
       		     omega=float(j)     			
     		     lon=omega 
      		     teta=180.-alfa(i)
		     lat=90.-teta 
       		     write(7,*) lon, lat
 		enddo
		write(7,'(a7,f14.8)')"360    ", -90.+alfa(i) 
   		write(7,'(a7)')      "360 -90" 
	    endif
!
!--- A north pole element		
 	   if(tetac(i).eq.0.) then 
	        write(7,'(a4,1x,F10.4)') '> -Z', h(i,k)  
   		write(7,'(a7)')      "0   90"
		write(7,'(a7,f14.8)')"0     ", 90.-alfa(i)		 
      		do j=1, 360, dstep  
       		     omega=float(j)     			
     		     lon=omega 
      		     teta=alfa(i)
		     lat=90.-teta 
       		     write(7,*) lon, lat
 		enddo
		write(7,'(a7,f14.8)')"360   ", 90.-alfa(i) 
   		write(7,'(a7)')      "360 90" 
	    endif
	    endif
!	
40          continue
!
!
!#---- "End of available ice models... 
!
           ENDIF
!

!#---- "End of multi-segmentation time steps
!
 20    CONTINUE
!
!

!
!#---- Formats for ICE3G
111    format(i4,2f9.3,2f6.3,1x,8f5.0)
112    format(15f5.0)
!
!#---- Formats for ICE5G
113    format(i4,f10.4,1x,i4,f10.4,1x,22(f10.4))
114    format(i4,f10.4,1x,i4,f10.4,1x,27(f10.4)) 	
!
!#---- Formats for ANU05
115    format(i3,1x,3(f10.4,1x),32(f10.4,1x))
!
!
!#---- Dellocate memory space
!	    	
    deallocate( ICR )
    deallocate( TETAC, LONGC, LATIC, ALFA )
    deallocate( H, CR ) 
!
!
!#----- 
!
 END PROGRAM MULTI_SEGMENT
!
!
!
