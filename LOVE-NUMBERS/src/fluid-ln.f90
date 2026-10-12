! This program computes the fluid 'k' LN for the model, assuming 
! a v/e lithosphere.
! DM Sep 11, 2025
!
!
 module constants
 use fmzm
 implicit none
!
 type(fm) :: pi
 type(fm) :: Gnwt
!
 end module constants
!
!
!
 module fluid_ln_data
 use fmzm
 implicit none
!
 integer :: nla
 integer :: mp
!
 character(200) :: path_2_model
!
 type(fm), allocatable :: r(:)
 type(fm), allocatable :: rho(:)
!
 type(fm), allocatable :: gra(:)
 type(fm), allocatable :: mass(:)
!
 type(fm) :: emass
!
 type(fm) :: r0, rho0
!
 end module fluid_ln_data
!
!
!
 program fluid_ln
 use fmzm
 use fluid_ln_data
 use constants
 use config_data
 implicit none
!
 type(fm) :: yf(2)
 integer :: n
 integer :: i
!
 type(fm) :: Ydir(2,2), Yinv(2,2), lambda(2,2)
 type(fm) :: c, qtidal
 type(fm) :: x(3)
 type(fm), allocatable :: hh(:), ll(:), kk(:)
!
 integer :: nmin, nmax
!
! --- Configuration
!
 call config
!
! --- Copy config parameters to local variables
!
 mp  = cfg_mp
 nla = cfg_nv + 2
 path_2_model = trim(adjustl(cfg_model_file))
!
 call fm_set(mp)
!
 call fm_pi(pi)
 Gnwt = to_fm('6.674e-11')
!
 allocate( r(0:nla)   )
 allocate( gra(0:nla) )
 allocate( rho(nla)   )
 allocate( mass(nla)  )
!
 call define_model
!
 nmin = 2
 nmax = 2
!
 allocate( hh(nmin:nmax) )
 allocate( ll(nmin:nmax) )
 allocate( kk(nmin:nmax) )
! 
 do n=nmin, nmax
!
! --- Set the initial conditions at the CMB
!     (we factor out a term r**(n-1)
!   yf(1) = r(1)**n
!   yf(2) = 2*(n-1)*r(1)**(n-1)
    yf(1) = r(1)
    yf(2) = to_fm( 2*(n-1) )
!
! --- Propagate to the surface
    do i=2,nla 
!
       call direct_fluid_matrix (n,r(i  ),rho(i),gra(i  ),Ydir)
       call inverse_fluid_matrix(n,r(i-1),rho(i),gra(i-1),Yinv)
!
       lambda = matmul(Ydir,Yinv)
       yf     = matmul(lambda,yf)
!
    end do
!
! ---- Apply the boundary conditions at the surface
!
    qtidal = -4 * pi * Gnwt * ( 2*n+1 ) / (4*pi*r(nla)**2)
!
    c = qtidal/yf(2)
!
    x(1)=-c*yf(1)/gra(nla)
    x(2)= to_fm('0')
    x(3)= c*yf(1)
!
    hh(n) = x(1) * emass / r(nla)
    ll(n) = x(2) * emass / r(nla)
    kk(n) = ( -to_fm('1') - x(3) * emass / ( r(nla) * gra(nla) ) )
!
 end do
!
 open(10, file='fluid-ln.dat', status='unknown')
!
 write(10,'(a)') '# Fluid LNs for model '//trim(path_2_model)
 write(10,'(a)') '# n, h_n, l_n, k_n' 
 do n=nmin,nmax
    write(10,*) n, to_dp(hh(n)), to_dp(ll(n)), to_dp(kk(n))
 end do
!
 close(10)
!
 end program fluid_ln
!
!
!
!
!
 subroutine define_model
 use fmzm
 use fluid_ln_data
 use constants
 use config_data
 implicit none
!
 integer :: i,k
 integer :: nh
 character(50) :: vc(7)
!
 open(1,file=path_2_model,status='old')
!
 call count_headers(1,nh)
 rewind(1)
!
 do i=1,nh
     read(1,*) 
 end do
!
 do k=nla,1,-1
    if     (cfg_model_fmt == 0) then
       read(1,*) ( vc(i),i=1,7 )
       r(k)   = to_fm(vc(3))
       rho(k) = to_fm(vc(4))
    elseif (cfg_model_fmt == 1) then
       read(1,*)  ( vc(i),i=1,5 )      ! New file format, DM 16092025
       r(k)    = to_fm(vc(2))
       rho(k)  = to_fm(vc(3))
    else
       stop 'unknown model file format'
    end if
 end do
!
 r(0) = to_fm(0)
!
! Reference scales
 r0   = r(nla)              ! Reference radius
 rho0 = maxval(rho)         ! Reference density
!
! NORMALIZED radii  
 r   = r / r0   
!
! NORMALIZED densities
 rho = rho / rho0
!
! NORMALIZED Newton constant
!Gnwt = Gnwt * rho0**2 * r0**2 / mu0 
! If all the model layers are perfectly fluid, the solution only depends 
! upon the g/G ratio (there is no rigidity scale).
! Hence we can set G=1
 Gnwt = to_fm(1)
!
! Mass of the layers
 do i=1, nla
     mass(i) = rho(i) * (r(i)**3 - r(i-1)**3)
 end do 
 mass = mass * to_fm(4) / to_fm(3) * pi
!
! Mass of the Earth
 emass = sum(mass)
!
! Gravity at the boundaries 
 gra(0)=to_fm(0)
 do i=1, nla
     gra(i) = sum(mass(1:i)) / r(i)**2
 end do
 gra = gra * Gnwt
!
!
 return
!
 end subroutine define_model
!
!
!
 subroutine direct_fluid_matrix(n,r,rho,gra,Y)
 use fmzm
 use constants
 implicit none
!
 type(fm) :: Y(2,2)
 integer  :: n
 type(fm) :: r
 type(fm) :: rho
 type(fm) :: gra
!
!
 Y(1,1) = r**n
 Y(1,2) = r**(-n-1)
!
 Y(2,1) = (2*n+1) * r**(n-1) - 4*pi*Gnwt*rho/gra * r**n
 Y(2,2) =  -4*pi*Gnwt*rho/gra * r**(-n-1)
! 
!
end subroutine direct_fluid_matrix
!
!
!
 subroutine inverse_fluid_matrix(n,r,rho,gra,Yinv)
 use fmzm
 use constants
 implicit none
!
 type(fm) :: Yinv(2,2)
 integer  :: n
 type(fm) :: r
 type(fm) :: rho
 type(fm) :: gra
!
!
 Yinv(1,1) = 4*pi*Gnwt*rho/gra * r**(-n+1)
 Yinv(1,2) = r**(-n+1)
!
 Yinv(2,1) = ((2*n+1) - 4*pi*Gnwt*rho/gra * r) * r**(n+1)
 Yinv(2,2) = -r**(n+2)
!
 Yinv = Yinv / (2*n+1);
!
!
 end subroutine inverse_fluid_matrix
!
!
!
 subroutine inv2(A,Ainv)
 use fmzm
 implicit none
!
 type(fm) :: A(2,2)
 type(fm) :: Ainv(2,2)
 type(fm), save :: det
!
 det = A(1,1) * A(2,2) - A(1,2) * A(2,1)
!
 Ainv(1,1) =  A(2,2) / det
 Ainv(1,2) = -A(1,2) / det
 Ainv(2,1) = -A(2,1) / det
 Ainv(2,2) = -A(1,1) / det
!
 return
!
 end subroutine inv2
!
