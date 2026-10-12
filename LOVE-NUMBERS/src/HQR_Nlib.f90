!
!
!
!
     SUBROUTINE HQR_NLIB(N,LOW,IGH,H,WR,WI,IERR)
     USE TYPES
     USE FMZM
     IMPLICIT NONE
!
! ------------------------------------------------------------------------
!    Tested and ported to quadruple precision by GS on January 6 2008.  
!    Original is from http://www.netlib.org/eispack/3090vf/single/hqr.f
! ------------------------------------------------------------------------
!
!    RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
!
     INTEGER I,J,K,L,M,N,EN,LL,MM,NA,NM,IGH,ITN,ITS,LOW,MP2,ENM2,IERR
     TYPE(FM) H(N,N),WR(N),WI(N), P,Q,R,S,T,W,X,Y,ZZ,NORM,TST1,TST2
     LOGICAL NOTLAS
!
!     this subroutine is a translation of the algol procedure hqr,
!     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!
!     this subroutine finds the eigenvalues of a real
!     upper hessenberg matrix by the qr method.
!
!     on input
!
!	 nm must be set to the row dimension of two-dimensional
!	   array parameters as declared in the Calling program
!	   dimension statement.
!
!	 n is the order of the matrix.
!
!	 low and igh are integers determined by the balancing
!	   subroutine  balanc.  if  balanc  has not been used,
!	   set low=1, igh=n.
!
!	 h contains the upper hessenberg matrix.  information about
!	   the transformations used in the reduction to hessenberg
!	   form by  elmhes  or  orthes, if performed, is stored
!	   in the remaining triangle under the hessenberg matrix.
!
!     on output
!
!	 h has been destroyed.  therefore, it must be saved
!	   before Calling  hqr  if subsequent calculation and
!	   back transformation of eigenvectors is to be performed.
!
!	 wr and wi contain the real and imaginary parts,
!	   respectively, of the eigenvalues.  the eigenvalues
!	   are unordered except that complex conjugate pairs
!	   of values appear consecutively with the eigenvalue
!	   having the positive imaginary part first.  if an
!	   error exit is made, the eigenvalues should be correct
!	   for indices ierr+1,...,n.
!
!	 ierr is set to
!	   zero       for normal return,
!	   j	      if the limit of 30*n iterations is exhausted
!		      while the j-th eigenvalue is being sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated september 1989.
!
!     ------------------------------------------------------------------
!
     ierr = 0
     norm = 0._qp
     k = 1
!     .......... store roots isolated by balanc
!		 and compute matrix norm ..........
     do 50 i = 1, n
!
 	do 40 j = k, n
  40	norm = norm + abs(h(i,j))
!
 	k = i
 	if (i .ge. low .and. i .le. igh) go to 50
 	wr(i) = h(i,i)
 	wi(i) = 0._qp
  50 continue
!
     en = igh
     t = 0._qp
     itn = 30*n
!     .......... search for next eigenvalues ..........
  60 if (en .lt. low) go to 1001
     its = 0
     na = en - 1
     enm2 = na - 1
!     .......... look for single small sub-diagonal element
!		 for l=en step -1 until low do -- ..........
  70 do 80 ll = low, en
 	l = en + low - ll
 	if (l .eq. low) go to 100
 	s = abs(h(l-1,l-1)) + abs(h(l,l))
 	if (s .eq. 0._qp) s = norm
 	tst1 = s
 	tst2 = tst1 + abs(h(l,l-1))
 	if (tst2 .eq. tst1) go to 100
  80 continue
!     .......... form shift ..........
 100 x = h(en,en)
     if (l .eq. en) go to 270
     y = h(na,na)
     w = h(en,na) * h(na,en)
     if (l .eq. na) go to 280
     if (itn .eq. 0) go to 1000
     if (its .ne. 10 .and. its .ne. 20) go to 130
!     .......... form exceptional shift ..........
     t = t + x
!
     do 120 i = low, en
 120 h(i,i) = h(i,i) - x
!
     s = abs(h(en,na)) + abs(h(na,enm2))
     x = 0.75_qp * s
     y = x
     w = -0.4375_qp * s * s
 130 its = its + 1
     itn = itn - 1
!     .......... look for two consecutive small
!		 sub-diagonal elements.
!		 for m=en-2 step -1 until l do -- ..........
     do 140 mm = l, enm2
 	m = enm2 + l - mm
 	zz = h(m,m)
 	r = x - zz
 	s = y - zz
 	p = (r * s - w) / h(m+1,m) + h(m,m+1)
 	q = h(m+1,m+1) - zz - r - s
 	r = h(m+2,m+1)
 	s = abs(p) + abs(q) + abs(r)
 	p = p / s
 	q = q / s
 	r = r / s
 	if (m .eq. l) go to 150
 	tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
 	tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
 	if (tst2 .eq. tst1) go to 150
 140 continue
!
 150 mp2 = m + 2
!
     do 160 i = mp2, en
 	h(i,i-2) = 0._qp
 	if (i .eq. mp2) go to 160
 	h(i,i-3) = 0._qp
 160 continue
!     .......... double qr step involving rows l to en and
!		 columns m to en ..........
     do 260 k = m, na
 	notlas = k .ne. na
 	if (k .eq. m) go to 170
 	p = h(k,k-1)
 	q = h(k+1,k-1)
 	r = 0._qp
 	if (notlas) r = h(k+2,k-1)
 	x = abs(p) + abs(q) + abs(r)
 	if (x .eq. 0._qp) go to 260
 	p = p / x
 	q = q / x
 	r = r / x
 170	s = sign(sqrt(p*p+q*q+r*r),p)
 	if (k .eq. m) go to 180
 	h(k,k-1) = -s * x
 	go to 190
 180	if (l .ne. m) h(k,k-1) = -h(k,k-1)
 190	p = p + s
 	x = p / s
 	y = q / s
 	zz = r / s
 	q = q / p
 	r = r / p
 	if (notlas) go to 225
!     .......... row modification ..........
 	do 200 j = k, EN
 	   p = h(k,j) + q * h(k+1,j)
 	   h(k,j) = h(k,j) - p * x
 	   h(k+1,j) = h(k+1,j) - p * y
 200	continue
!
 	j = min0(en,k+3)
!     .......... column modification ..........
 	do 210 i = L, j
 	   p = x * h(i,k) + y * h(i,k+1)
 	   h(i,k) = h(i,k) - p
 	   h(i,k+1) = h(i,k+1) - p * q
 210	continue
 	go to 255
 225	continue
!     .......... row modification ..........
 	do 230 j = k, EN
 	   p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
 	   h(k,j) = h(k,j) - p * x
 	   h(k+1,j) = h(k+1,j) - p * y
 	   h(k+2,j) = h(k+2,j) - p * zz
 230	continue
!
 	j = min0(en,k+3)
!     .......... column modification ..........
 	do 240 i = L, j
 	   p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
 	   h(i,k) = h(i,k) - p
 	   h(i,k+1) = h(i,k+1) - p * q
 	   h(i,k+2) = h(i,k+2) - p * r
 240	continue
 255	continue
!
 260 continue
!
     go to 70
!     .......... one root found ..........
 270 wr(en) = x + t
     wi(en) = 0._qp
     en = na
     go to 60
!     .......... two roots found ..........
 280 p = (y - x) / 2._qp
     q = p * p + w
     zz = sqrt(abs(q))
     x = x + t
     if (q .lt. 0._qp) go to 320
!     .......... real pair ..........
     zz = p + sign(zz,p)
     wr(na) = x + zz
     wr(en) = wr(na)
     if (zz .ne. 0._qp) wr(en) = x - w / zz
     wi(na) = 0._qp
     wi(en) = 0._qp
     go to 330
!     .......... complex pair ..........
 320 wr(na) = x + p
     wr(en) = x + p
     wi(na) = zz
     wi(en) = -zz
 330 en = enm2
     go to 60
!     .......... set error -- all eigenvalues have not
!		 converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
!
!
!
!
!
