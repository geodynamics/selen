!
!
!
!=====================================================================
 SUBROUTINE BALANC_NLIB(N,A,LOW,IGH,ASCALE)
 USE TYPES 
 USE FMZM
 IMPLICIT NONE
 SAVE
!
 INTEGER I,J,K,L,M,N,NM,IGH,LOW,IEXC
 TYPE(FM) A(N,N),ASCALE(N)
 TYPE(FM) C,F,G,R,S,B2,ARADIX
 LOGICAL NOCONV
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Tested and ported to quadruple precision by GS on January 6 2008.  
!    Original is from http://www.netlib.org/eispack/3090/double/balanc.f
!    Retouched by GS on December 1019 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     this subroutine is a translation of the algol procedure balance,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!     this subroutine balances a real matrix and isolates
!     eigenvalues whenever possible.
!     on input
!	 nm must be set to the row dimension of two-dimensional
!	   array parameters as declared in the calling program
!	   dimension statement.
!	 n is the order of the matrix.
!	 a contains the input matrix to be balanced.
!     on output
!	 a contains the balanced matrix.
!	 low and igh are two integers such that a(i,j)
!	   is equal to zero if
!	    (1) i is greater than j and
!	    (2) j=1,...,low-1 or i=igh+1,...,n.
!	 scale contains information determining the
!	    permutations and scaling factors used.
!     suppose that the principal submatrix in rows low through igh
!     has been balanced, that p(j) denotes the index interchanged
!     with j during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by d(i,j).  then
!	 scale(j) = p(j),    for j = 1,...,low-1
!		  = d(j,j),	 j = low,...,igh
!		  = p(j)	 j = igh+1,...,n.
!     the order in which the interchanges are made is n to igh+1,
!     then 1 to low-1.
!     note that 1 is returned for igh if igh is zero formally.
!     the algol procedure exc contained in balance appears in
!     balanc  in line.  (note that the algol roles of identifiers
!     k,l have been reversed.)
!     Questions and comments should be directed to Alan K. Cline,
!     Pleasant Valley Software, 8603 Altus Cove, Austin, TX 78759.
!     Electronic mail to cline@cs.utexas.edu.
!     this version dated january 1989. (for the IBM 3090vf)
!     ------------------------------------------------------------------
!
     aradix = to_FM(16.0) 
!
     b2 = aradix * aradix
     k = 1
     l = n
     go to 100
!     .......... in-line procedure for row and
!		 column exchange ..........
  20 ascale(m) = j
     if (j .eq. m) go to 50
!
     do 30 i = 1, l
 	f = a(i,j)
 	a(i,j) = a(i,m)
 	a(i,m) = f
  30 continue
!
     do 40 i = k, n
 	f = a(j,i)
 	a(j,i) = a(m,i)
 	a(m,i) = f
  40 continue
!
  50 go to (80,130), iexc
!     .......... search for rows isolating an eigenvalue
!		 and push them down ..........
  80 if (l .eq. 1) go to 280
     l = l - 1
!     .......... for j=l step -1 until 1 do -- ..........
 100 do 120 j = l, 1, -1
!
 	do 110 i = 1, l
 	   if (i .eq. j) go to 110
 	   if (a(j,i) .ne. 0._qp) go to 120
 110	continue
!
 	m = l
 	iexc = 1
 	go to 20
 120 continue
!
     go to 140
!     .......... search for columns isolating an eigenvalue
!		 and push them left ..........
 130 k = k + 1
!
 140 do 170 j = k, l
!
 	do 150 i = k, l
 	   if (i .eq. j) go to 150
 	   if (a(i,j) .ne. 0._qp) go to 170
 150	continue
!
 	m = k
 	iexc = 2
 	go to 20
 170 continue
!     .......... now balance the submatrix in rows k to l ..........
     do 180 i = k, l
 180 ascale(i) = 1._qp
!     .......... iterative loop for norm reduction ..........
 190 noconv = .false.
!
     do 270 i = k, l
 	c = 0._qp
 	r = 0._qp
!
 	do 200 j = k, l
 	   if (j .eq. i) go to 200
 	   c = c + abs(a(j,i))
 	   r = r + abs(a(i,j))
 200	continue
!     .......... guard against zero c or r due to underflow ..........
 	if (c .eq. 0._qp .or. r .eq. 0._qp) go to 270
 	g = r / aradix
 	f = 1._qp
 	s = c + r
 210	if (c .ge. g) go to 220
 	f = f * aradix
 	c = c * b2
 	go to 210
 220	g = r * aradix
 230	if (c .lt. g) go to 240
 	f = f / aradix
 	c = c / b2
 	go to 230
!     .......... now balance ..........
 240	if ((c + r) / f .ge. 0.95_qp * s) go to 270
 	g = 1._qp / f
 	ascale(i) = ascale(i) * f
 	noconv = .true.
!
 	do 250 j = k, n
 250	a(i,j) = a(i,j) * g
!
 	do 260 j = 1, l
 260	a(j,i) = a(j,i) * f
!
 270 continue
!
        if (noconv) go to 190
!
  280 low = k
      igh = l
      return
      end
!---------------------------------------- 
!
!
!
