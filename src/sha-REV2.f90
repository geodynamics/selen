!
!
  PROGRAM SH
!+++++++++++++++++++++++++++++++++++++++
!
! Builds a table of spherical harmonics 
! DM April 18, 2019
! based on 'SHA-REV1.F90', GS June 15, 2016
!
!+++++++++++++++++++++++++++++++++++++++
  implicit NONE
!
  INTEGER :: RES
  INTEGER :: NP
  INTEGER :: LMAX
  INTEGER :: JMAX
  CHARACTER*72 :: PXFILE
  CHARACTER*72 :: PLFILE
  CHARACTER*72 :: SHFILE
  CHARACTER*200 :: BUFFER
  CHARACTER*20   :: CHR, CHL
!
  INTEGER, PARAMETER :: A_LARGE_INTEGER=10**6
!
  INTEGER I, J, K, L, NA
  REAL*8 COSDD, SINDD, LATA, LONG   
  REAL*8,     DIMENSION (:,:), ALLOCATABLE:: ALF
  COMPLEX*16, DIMENSION (:,:), ALLOCATABLE:: LONG_TABLE
  REAL*8, ALLOCATABLE :: LON(:), LAT(:)
  REAL*8  JUNK
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! --- Examining command-line arguments
!
  if( (iargc().ne.2) .and. (iargc().ne.3) .and. (iargc().ne.4) .and. (iargc().ne.5) ) then
     write (*,*) 
     write (*,*) "USAGE: sha.exe RES LMAX <px-file> <pxlat-file> <sh-file>"
     write (*,*)  
     write (*,*) "  RES                  Tegmark resolution"
     write (*,*) "  LMAX                 Maximum harmonic degree"
     write (*,*) "  px-file, pxlat-file  Files with pixels and main latitudes."
     write (*,*) "                         (if not specified, the default is "
     write (*,*) "                          px-Rnn.dat and px-lat-Rnn.dat )"
     write (*,*) "  sh-file              Output file with ALFs and TRIGs."
     write (*,*) "                         (if not specified, the default is "
     write (*,*) "                          sh-RnnLnn.dat)"
     write (*,*) 
     stop
  end if  
!  
  call getarg(1,buffer)   ; read(buffer,*) res
  call getarg(2,buffer)   ; read(buffer,*) lmax
!
  write(chr,*) res  ; chr=adjustl(chr)
  write(chl,*) lmax ; chl=adjustl(chl) 
!
  if( iargc().eq.2 ) then         ! No file is specified
      pxfile = 'px-R'//trim(chr)//'.dat'
      plfile = 'px-lat-R'//trim(chr)//'.dat'
      shfile='sh-R'//trim(chr)//'L'//trim(chl)//'.bin'
  elseif( iargc().eq.3 ) then     ! Only the output file is specified
      call getarg(3,shfile)
      pxfile = 'px-R'//trim(chr)//'.dat'
      plfile = 'px-lat-R'//trim(chr)//'.dat'
  elseif( iargc().eq.4 ) then     ! Only the pixel files are specified
      call getarg(3,pxfile)
      call getarg(4,plfile)
      shfile='sh-R'//trim(chr)//'L'//trim(chl)//'.bin'
  elseif( iargc().eq.5 ) then     ! All files are specified
      call getarg(3,pxfile)
      call getarg(4,plfile)
      call getarg(5,shfile)
  else
      write(*,*) ' Wrong number of arguments ! '
      stop
  end if
!
  NP=2*RES*(RES-1)*20+12
  JMAX=(LMAX+1)*(LMAX+2)/2
!
!
!
!
  write(*,*)" ---- Maximum degree: ", lmax 
  write(*,*)" ---- Resolution: ", res 
  write(*,*)" ---- Number of pixels: ", np 
!
!
  write(*,*) " ---- Reading pixels data from file: ", trim(adjustl(PXFILE))
  open(10,file=PXFILE,status="old")
  do i=1,np
     read(10,*) junk, junk, na
  end do
  close(10)
  Write(*,*) " ---- Number of MAIN pixels: ", NA
!
  ALLOCATE(LAT(NA))
!
  write(*,*) " ---- Reading pixels data from file: ", trim(adjustl(PLFILE))
  open(10,file=PLFILE,status="old")
  do i=1,na
     read(10,*) junk, lat(i) 
  end do
  close(10)
!
  ALLOCATE(LON(NP))
  write(*,*) " ---- Reading longitude data" 
  open(10,file=PXFILE,status="old")
  do i=1,np
     read(10,*) lon(i) 
  end do
  close(10)
!
!
!
  write(*,*) " ---- Opening file: ", trim(adjustl(SHFILE))
  open(20,file=SHFILE,status='unknown',form='unformatted') 
!
!
  Write(*,*) " ---- Pre-computing ALFs at latitudes of anchor pixels"
  ALLOCATE(ALF(JMAX,NA))
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! CAUTION! PLMBAR_MOD pre-computes some expressions each time it is 
! called with a new lmax and stores them in a SAVEd array for reuse 
! in subsequent calls. This is thread-safe IF and ONLY IF (1) plmbar_mod
! is called BEFORE the parallel region and (2) all the calls in the
! parallel region use the same lmax.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!  
  Write(*,*) " ---- First call to PLMBAR_MOD"
  call PLMBAR_MOD( lmax, lat(1), ALF(:,1) )
  Write(*,*) " ---- Other calls"
!
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NA,ALF,LAT,LMAX)
!$OMP DO SCHEDULE(GUIDED)
  do i=2, na  
     if(mod(i,5000)==0)write(*,*) i, " of ", na
     lata = lat(i)
     call PLMBAR_MOD (lmax, lata, ALF(:,i))        
  enddo
!$OMP END PARALLEL
!
  write(*,*) " ---- Writing the Legendre functions on file: ", trim(adjustl(SHFILE))
  WRITE(20) ALF
  DEALLOCATE(ALF) 
!
!
  Write(*,*) " ---- Pre-computing TRIG functions at pixels"
  ALLOCATE(LONG_TABLE(0:LMAX,NP))
!
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(LONG_TABLE,LON,NP,LMAX) SCHEDULE(GUIDED)
  do i=1, np  
     long = lon(i)
     do l=0, lmax
	  long_table(l,i)=cmplx(cosdd(l*long),sindd(l*long))
     enddo	
  enddo
!$OMP END PARALLEL DO
!
  write(*,*) " ---- Writing the sines and cosines on file: ", trim(adjustl(SHFILE))
  WRITE(20) LONG_TABLE
  DEALLOCATE(LONG_TABLE)
!
  CLOSE(20) 
!
  DEALLOCATE(LON)
  DEALLOCATE(LAT)
!
!
 write(*,*) " ***************************************************** "
 write(*,*) " - Copy or better move the following file into ../DATA " 
 write(*,*) trim(adjustl(SHFILE))
 write(*,*) " ***************************************************** "
 write(*,*) " Thank you!" 
!  
!
  end program SH
!
!
!
