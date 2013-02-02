!
! This is program "SLE.F90" 
!
!       GS 04-11-2008 "Intel port..."
! 		Modified by GS september 2008 for horizontal movements 
! 		Reviewed GS & FC July 2009 -  "Varying coastlines" (reprise!) 
! 		*** Reviewed GS & FC November 2009 - Porting under gfortran 
! 		Also reviewed in April 2010 for the implementation of degree 1 
! 		(remember the mode coupling issue...) 
! 		Revised May 26 GS for new ice SH  
!       Revised Jun 16, 2011 DM for dynamic memory allocation
!       Revised Jun 2011 DM for parallel execution
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
! SELEN is distributed in the /hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! ============================================
! Solves the SLE by the pseudo-spectral method 
! ============================================
!
! Input files:
!       - px-table.dat 
!	- shof.dat 
!	- shice.dat
!	- eb.dat
!	- ebu.dat
!       - ebv.dat  <<< new... 
! 	- sh.bin
!
! Output files:
!	- shu.bin 
!	- shn.bin 
!	- shs.bin 
!       - shz.bin 
!       - shv.bin <<< new... 
!
!
! INCLUDE "harmonics.f90"
 PROGRAM SLE 
 IMPLICIT NONE 
 INCLUDE "data.inc"
 INTEGER ijunk
 INTEGER :: npix, na
 CHARACTER*22, PARAMETER :: FILES='shs.bin', FILEU='shu.bin', & 
 			    FILEN='shn.bin', FILEZ='shz.bin', & 
			    FILEV='shv.bin'
 CHARACTER*12 HEADER 
 INTEGER I, J, K, L, P, IS, LJ, MJ, LI, DOM, IND
! INTEGER LL(JMAX), MM(JMAX), DM(JMAX), ANC(NP), WET(NP) 
! REAL*4 ALF(JMAX,NANCH), LONP(NP), LATP(NP), X(NP,0:NN)  
 INTEGER, ALLOCATABLE :: LL(:), MM(:), DM(:), ANC(:), WET(:)
 REAL*8, ALLOCATABLE :: ALF(:,:), LONP(:), LATP(:), X(:,:)  
! REAL*4 BETAS(0:LMAX,0:NN), ES(0:LMAX)
! REAL*4 BETAU(0:LMAX,0:NN), EU(0:LMAX)
! REAL*4 BETAN(0:LMAX,0:NN), EN(0:LMAX)
! REAL*4 BETAV(0:LMAX,0:NN), EV(0:LMAX)      
! REAL*4 RESH, IMSH, AAVV(0:NN), BBVV(0:NN)
 REAL*8, ALLOCATABLE :: BETAS(:,:), ES(:)
 REAL*8, ALLOCATABLE :: BETAU(:,:), EU(:)
 REAL*8, ALLOCATABLE :: BETAN(:,:), EN(:)
 REAL*8, ALLOCATABLE ::  BETAV(:,:), EV(:)      
 REAL*8 RESH, IMSH
 REAL*8, ALLOCATABLE :: AAVV(:), BBVV(:) 
! COMPLEX*16 LONG_TABLE(0:LMAX,NP), OC(JMAX) 
! COMPLEX*16 IIII(JMAX,0:NN), ZE(JMAX,0:NN), SE(JMAX,0:NN)
! COMPLEX*16 AAAA(JMAX,0:NN), AAAA_MOD(JMAX,0:NN)
! COMPLEX*16 BBBB(JMAX,0:NN), BBBB_MOD(JMAX,0:NN) 
! COMPLEX*16 HHHH(JMAX,0:NN), KKKK(JMAX,0:NN) 
! COMPLEX*16 S(JMAX,0:NN), N(JMAX,0:NN), U(JMAX,0:NN), V(JMAX,0:NN) 
! COMPLEX*16 Z(JMAX,0:NN,0:SMAX)
 COMPLEX*16, ALLOCATABLE :: LONG_TABLE(:,:), OC(:)
 COMPLEX*16, ALLOCATABLE :: IIII(:,:), ZE(:,:), SE(:,:)
 COMPLEX*16, ALLOCATABLE :: AAAA(:,:), AAAA_MOD(:,:)
 COMPLEX*16, ALLOCATABLE :: BBBB(:,:), BBBB_MOD(:,:) 
 COMPLEX*16, ALLOCATABLE :: HHHH(:,:), KKKK(:,:) 
 COMPLEX*16, ALLOCATABLE :: S(:,:), N(:,:), U(:,:), V(:,:)
 COMPLEX*16, ALLOCATABLE :: Z(:,:,:)
 REAL*8 RHOE, RHOI_O_RHOE_X3, RHOW_O_RHOE_X3, RHOI_O_RHOW 
!
!
! --- Determining the number of pixels
!
  open(1,file='anchor.tmp',status='old')
  read(1,*) na
  close(1)
  Write(*,*) "    - Found ", na, " anchor pixels in file px-lat.dat"
  npix=np
!
!
! ========================================================
! ========================================================
! ========================================================
!
! --- Allocate memory space
!
 ALLOCATE( LL(JMAX), MM(JMAX), DM(JMAX), ANC(NPIX), WET(NPIX) )
 ALLOCATE( ALF(JMAX,NA), LONP(NPIX), LATP(NPIX), X(NPIX,0:NN) )
 ALLOCATE( LONG_TABLE(0:LMAX,NPIX), OC(JMAX) )
 ALLOCATE( IIII(JMAX,0:NN), ZE(JMAX,0:NN), SE(JMAX,0:NN) )
 ALLOCATE( AAAA(JMAX,0:NN), AAAA_MOD(JMAX,0:NN) )
 ALLOCATE( BBBB(JMAX,0:NN), BBBB_MOD(JMAX,0:NN) )
 ALLOCATE( HHHH(JMAX,0:NN), KKKK(JMAX,0:NN) )
 ALLOCATE( Z(JMAX,0:NN,0:SMAX) )
 ALLOCATE( S(JMAX,0:NN), N(JMAX,0:NN), U(JMAX,0:NN), V(JMAX,0:NN) )
 ALLOCATE( BETAS(0:LMAX,0:NN), ES(0:LMAX) )
 ALLOCATE( BETAU(0:LMAX,0:NN), EU(0:LMAX) )
 ALLOCATE( BETAN(0:LMAX,0:NN), EN(0:LMAX) )
 ALLOCATE( BETAV(0:LMAX,0:NN), EV(0:LMAX) )     
 ALLOCATE( AAVV(0:NN), BBVV(0:NN) ) 
!
!
! ********************************************************
! ********************************************************
! ********************************************************
!
! --- Extracting the average Earth density from 
!     the TABOO or ALMA log files... GS July 09
!
 IF(CDE.GE. 0) IND=1 
 IF(CDE.EQ.-1) IND=2 
 CALL AVERAGE_EARTH_DENSITY(IND, RHOE)
 RHOI_O_RHOE_X3 = 3.*RHOI/RHOE 
 RHOW_O_RHOE_X3 = 3.*RHOW/RHOE 
 RHOI_O_RHOW    =    RHOI/RHOW   
!
!
! --- Pre-computing 'l' and 'm' corresponding to degree 'J'
!
	do j=1, jmax 
		mm(j)=mj(j) 
		ll(j)=lj(j)
		dm(j)=2-dom(j)
	enddo	
!
! --- Reading the ALFs table from <<sh.bin>>
! 
 	Write(*,*) '    - Reading the ALFs from file sh.bin'
 	open(3,file='sh.bin',status='unknown',form='unformatted') 
		read(3)ALF
		read(3)LONG_TABLE
 	Close(3) 
!
!
! --- Examining the pixels table & extracting information
!
	Open(1,file='px-table.dat',status='unknown') 
	Do i=1, 4 
		Read(1,'(a12)')header
	Enddo
	Do i=1, np 
		Read (1,*) lonp(i), latp(i), anc(i), k, wet(i) 	
	Enddo
	close(1)
!
!
! --- Reading the SH OF coefficients from shof.dat 
!
 	Write(*,*) '    - Reading the SH OF coeff. from shof.dat'
 	open(3,file='shof.dat',status='unknown')
 	do j=1, jmax   
		read(3,*) k, resh, imsh 
                oc(j)=cmplx(resh, imsh)	
 	enddo
 	close(3)
!
!
!
! --- Reading the <<I>> array 
	IIII(:,:)=0. 
	If(imode/=5) then 
 		Write(*,*) "    - Reading array 'I'"
		open(1,file='shice.dat',status='unknown') 
		read(1,*) IIII
		close(1)
    Endif
!
!
!
! --- Reading the E and BETA arrays 
 	Write(*,*) "    - Reading arrays 'E' and 'Beta'"
	open(1,file='ebs.dat', status='unknown')    
	open(2,file='ebu.dat',status='unknown')    
	open(3,file='ebn.dat',status='unknown')    
	open(4,file='ebv.dat',status='unknown')    

	do li=0, lmax
		read(1,*) l, Es(l), (betas(l,k), k=0,nn) 
		read(2,*) l, Eu(l), (betau(l,k), k=0,nn) 
		read(3,*) l, En(l), (betan(l,k), k=0,nn) 
		read(4,*) l, Ev(l), (betav(l,k), k=0,nn) 
	enddo
	close(4); close(3) ; close(2) ; close(1) 
!
!
! --- Computing the eustatic Z array...
	ze(:,:) = 0.  		
	do k=0,nn	
		ze(:,k) = - rhoi_o_rhow*(iiii(1,k)/oc(1))*oc(:)
	enddo
!
!
! --- Computing the eustatic S array...
	se(:,:) = 0.
	se(1,:) = - rhoi_o_rhow*(iiii(1,:)/oc(1)) 
!
!
! --- Computing the A array...
	aaaa(:,:)=0.
	do j=1, jmax 
 	    do k=0,NN
            aaaa(j,k) = ES(ll(j))*IIII(j,k)        
 	      do p=0, k
              if(p==0) aaaa(j,k) = aaaa(j,k)-(IIII(j,p)            )*BETAS(ll(j),k-p)
    	      if(p/=0) aaaa(j,k) = aaaa(j,k)-(IIII(j,p)-IIII(j,p-1))*BETAS(ll(j),k-p)  
 	      enddo 
	      aaaa(j,k)= RHOI_o_RHOE_x3*aaaa(j,k)
	      enddo
	enddo
!
!
! --- Computing the ocean average of A 
	do k=0, NN 
		aavv(k)=0.
		do j=1, jmax 
		aavv(k) = aavv(k) + & 
		          dm(j)*real(oc(j)*conjg(aaaa(j,k)))/oc(1)  
		enddo
	enddo 
!
!
! --- Computing the modified ocean average of A 
	aaaa_mod(:,:) = aaaa(:,:)
	aaaa_mod(1,:) = aaaa(1,:)-aavv(:) 
!
! print *,'Task ',my_rank,' starting loop 1'
!
! --- Computing the R-array...
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP     SHARED(X,ALF,ANC,DM,AAAA_MOD,LONG_TABLE,MM,NPIX) &
!$OMP     SCHEDULE(GUIDED)
	do i=1, npix
	  do k=0, NN
	    x(i,k)=0.			
	    do j=1, jmax  
           x(i,k) = x(i,k) + & 
			 ALF(j,anc(i))*dm(j)*real(aaaa_mod(j,k)*long_table(mm(j),i)) 
		enddo   
	  enddo
	enddo 
!$OMP END PARALLEL DO
!
! print *,'Task ',my_rank,' starting loop 2'

        hhhh = 0.
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(WET,HHHH,X,ALF,ANC,LONG_TABLE,MM,NPIX) &
!$OMP       SCHEDULE(GUIDED)
      do j=1,jmax
 	  do i=1,npix
	    if(wet(i)==1) then
	      do k=0, NN
                 hhhh(j,k) = hhhh(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i)) 	          
	      end do
	    end if
	  end do
	end do   
!$OMP END PARALLEL DO
!
! It was:
!	do k=0, NN
!	do j=1, jmax 
!	    hhhh(j,k)=0.  
!	    do i=1, np 
!	    if(wet(i)==1) hhhh(j,k) = & 
!			  hhhh(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
!	    enddo
!	enddo
!	enddo 
!
!
!
!
	hhhh(:,:)=hhhh(:,:)/float(np) + ze(:,:)
!
!
! --- Initializing the Z and S arrays 
	Z(:,:,0) = ZE(:,:)
	S(:,:)   = SE(:,:) 
!
!
! --- No recursion for the "Explicit approach"
	if(imode==6.or.imode==7.or.imode==3.or.SMAX==0) goto 2000     
!
!
!
! -----------------------
! ---    Recursion    ---
! -----------------------
!
	Write(*,*) "    - Starting the recursion"
!
 	do 1000 is = 1, SMAX     
!
!
        write(*,'(a12,i2,a3,i2)') '     - step ', is, ' of', SMAX
!
!
! --- Computing the 'B' array...
        bbbb(:,:)=0.
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,K,P) &
!$OMP     SHARED(BBBB,ES,LL,Z,IS,BETAS,rhow_o_rhoe_x3) &
!$OMP     SCHEDULE(GUIDED)
        do j=1, jmax 
	    do k=0,NN
	    bbbb(j,k) = ES(ll(j))*Z(j,k,is-1)        
	        do p=0, k
	        if(p==0) bbbb(j,k) = bbbb(j,k)-(Z(j,p,is-1)-0.           )*BETAS(ll(j),k-p)
	        if(p/=0) bbbb(j,k) = bbbb(j,k)-(Z(j,p,is-1)-Z(j,p-1,is-1))*BETAS(ll(j),k-p)		 
	        enddo 
            bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        enddo
	enddo
!$OMP END PARALLEL DO
!	
!
! --- Computing the ocean-average of array B array
	bbvv(:)=0.
	do k=0, NN 
		bbvv(k)=0.
		do j=1, jmax 
		bbvv(k) = bbvv(k) + & 
		          dm(j)*real(oc(j)*conjg(bbbb(j,k)))/oc(1)
		enddo
	enddo 
!
!
! --- Computing modified 'B' array
	bbbb_mod(:,:)=bbbb(:,:)
	bbbb_mod(1,:)=bbbb(1,:)-bbvv(:) 
!
!
! --- Computing array K...
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,K,J) &
!$OMP    SHARED(X,ALF,ANC,DM,BBBB_MOD,LONG_TABLE,MM,NPIX) SCHEDULE(GUIDED)
	do i=1, npix 
!	
	do k=0, NN
	    x(i,k)=0.			
	    do j=1, jmax  
            x(i,k) = x(i,k) + & 
			    ALF(j,anc(i))*dm(j)*real(bbbb_mod(j,k)*long_table(mm(j),i)) 
		enddo   
	enddo

	enddo 
!$OMP END PARALLEL DO
!
!
!
      kkkk = 0.
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J,K) &
!$OMP       SHARED(WET,KKKK,X,ALF,ANC,LONG_TABLE,MM,NPIX) &
!$OMP       SCHEDULE(GUIDED)
       do j=1,jmax
 	  do i=1,npix
	    if(wet(i)==1) then
	      do k=0, NN
                 kkkk(j,k) = kkkk(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i)) 	          
	      end do
	    end if
	  end do
	end do   
!$OMP END PARALLEL DO
!
!--IT WAS:---------------------------------------------------------------------------
!	do k=0, NN
!	do j=1, jmax 
!	    kkkk(j,k)=0.  
!	    do i=1, np 
!	    if(wet(i)==1) kkkk(j,k) = & 
!			  kkkk(j,k) + x(i,k)*alf(j,anc(i))*conjg(long_table(mm(j),i))  
!	    enddo
!	enddo
!	enddo 
!------------------------------------------------------------------------------------
!
!
	kkkk(:,:)=kkkk(:,:)/float(np) 
!
!
! --- Solving for arrays 'Z' and 'S' 
	Z(:,:,is) = HHHH(:,:) + KKKK(:,:) 
	S(:,:)    = AAAA_MOD(:,:) + SE(:,:) + BBBB_MOD(:,:)

!
! ------------------------------
! ---    End of recursion    ---
! ------------------------------
!
1000 CONTINUE      
!
!
!
2000 CONTINUE       
!
!
! --- Eustatic solution: U=0, S=N, V=0  
        if(imode==3.or.smax==0) then 
				U(:,:) = 0.
				S(:,:) = SE(:,:)
				N(:,:) = S(:,:)
				V(:,:) = 0. 
				goto 3000
				endif
!
! --- Array "B" for vertical displacement 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = EU(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAU(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAU(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for vertical displacement 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
 	 	aaaa(j,k) = EU(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAU(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAU(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*aaaa(j,k)
 		enddo
	enddo
!
!
! --- Vertical displacement
	U(:,:) = aaaa(:,:) + bbbb(:,:)

!
! --- Array "B" for Geoid heigth   
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = EN(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAN(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAN(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for Geoid heigth  
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
 	 	aaaa(j,k) = EN(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAN(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAN(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*aaaa(j,k)
 		enddo
	enddo
!
! --- Geoid undulations 
	N(:,:) = aaaa(:,:) + bbbb(:,:)
!
! --- Adding a constant to geoid undulations 
	N(1,:) = N(1,:) +  SE(1,:) - AAVV(:) - BBVV(:)  
!
! --- Geoid undulations (previous formulation) 
!       N(:,:) = S(:,:) + U(:,:)
!
!
! --- Array "B" for horizontal displacement 
	bbbb(:,:)=0.
 	do j=1, jmax 
		do k=0,NN
		bbbb(j,k) = EV(ll(j))*Z(j,k,SMAX)       
		do p=0, k
		if(p==0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-0.           )*BETAV(ll(j),k-p)
		if(p/=0)  bbbb(j,k) = bbbb(j,k) - (Z(j,p,SMAX)-Z(j,p-1,SMAX))*BETAV(ll(j),k-p)
		enddo 
        	bbbb(j,k)= RHOW_o_RHOE_X3*bbbb(j,k)
        	enddo
	enddo
!
! --- Array "A" for horizontal displacement 
 	aaaa(:,:)=0.
 	do j=1, jmax 
	 	do k=0,NN       
 	 	aaaa(j,k) = EV(ll(j))*IIII(j,k) 
 	 	do p=0, k
 	 	if(p==0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-0.         )*BETAV(ll(j),k-p) 	 
 	 	if(p/=0) aaaa(j,k) = aaaa(j,k) - (IIII(j,p)-IIII(j,p-1))*BETAV(ll(j),k-p)
 		enddo 
 		aaaa(j,k)= RHOI_o_RHOE_X3*aaaa(j,k)
 		enddo
	enddo
!
!
! --- Horizontal displacement
	V(:,:) = aaaa(:,:) + bbbb(:,:)
!
!
3000 CONTINUE 
!
!
!
! --- Writing SH coefficients of sea level variations, S 
 	open(3,file=files,status='unknown',form='unformatted') 
        write(3) S ; close(3)
!
! --- Writing SH coefficients of vertical displacement, U  
 	open(3,file=fileu,status='unknown',form='unformatted') 
        write(3) U ; close(3) 
!
! --- Writing SH coefficients of geoid undulations, N    
 	open(3,file=filen,status='unknown',form='unformatted') 
        write(3) N ; close(3) 
!
! --- Writing SH coefficients of 'reduced' sea level change, Z=OS 
 	open(3,file=filez,status='unknown',form='unformatted') 
        write(3) Z ; close(3)
!
! --- Writing SH coefficients of horizontal displacement, V 
 	open(3,file=filev,status='unknown',form='unformatted') 
        write(3) V ; close(3)
!
	Write(*,*) "    +++ SLE solved <3 +++" 
!
!
!
! ========================================================
! ========================================================
! ========================================================
!
! --- Free up memory space
!
 DEALLOCATE( LL, MM, DM, ANC, WET )
 DEALLOCATE( ALF, LONP, LATP, X )
 DEALLOCATE( LONG_TABLE, OC )
 DEALLOCATE( IIII, ZE, SE )
 DEALLOCATE( AAAA, AAAA_MOD )
 DEALLOCATE( BBBB, BBBB_MOD )
 DEALLOCATE( HHHH, KKKK )
 DEALLOCATE( Z )
 DEALLOCATE( S, N, U, V )
 DEALLOCATE( BETAS, ES )
 DEALLOCATE( BETAU, EU )
 DEALLOCATE( BETAN, EN )
 DEALLOCATE( BETAV, EV )     
 DEALLOCATE( AAVV, BBVV ) 
!
!
!
	End program SLE
!
!
!
!
