
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
* ======================================

      DOUBLE PRECISION FUNCTION dasum(N,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
* Purpose:
* =============
*
*    DASUM takes the sum of the absolute values.
*
* Arguments:
* ==========
*
*    N is INTEGER number of elements in input vector(s)
*
*    DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*
*    INCX is INTEGER storage spacing between elements of DX
*
* Further Details:
* =====================
*
*      jack dongarra, linpack, 3/11/78.
*      modified 3/93 to return if incx .le. 0.
*      modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dabs,mod
*     ..
      dasum = 0.0d0
      dtemp = 0.0d0
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
*        code for increment equal to 1
*
*
*        clean-up loop
*
         m = mod(n,6)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dtemp + dabs(dx(i))
            END DO
            IF (n.LT.6) THEN
               dasum = dtemp
               RETURN
            END IF
         END IF
         mp1 = m + 1
         DO i = mp1,n,6
            dtemp = dtemp + dabs(dx(i)) + dabs(dx(i+1)) +
     $              dabs(dx(i+2)) + dabs(dx(i+3)) +
     $              dabs(dx(i+4)) + dabs(dx(i+5))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            dtemp = dtemp + dabs(dx(i))
         END DO
      END IF
      dasum = dtemp
      RETURN
      END

      SUBROUTINE daxpy(N,DA,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
* Purpose:
* =============
*
*    DAXPY constant times a vector plus a vector.
*    uses unrolled loops for increments equal to one.
*
* Arguments:
* ==========
*
*    N is INTEGER number of elements in input vector(s)
*
*    DA is DOUBLE PRECISION. On entry, DA specifies the scalar alpha.
*
*    DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*
*    INCX is INTEGER storage spacing between elements of DX
*
*    DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*
*    INCY is INTEGER storage spacing between elements of DY
*
* Further Details:
* =====================
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      IF (n.LE.0) RETURN
      IF (da.EQ.0.0d0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         m = mod(n,4)
         IF (m.NE.0) THEN
            DO i = 1,m
               dy(i) = dy(i) + da*dx(i)
            END DO
         END IF
         IF (n.LT.4) RETURN
         mp1 = m + 1
         DO i = mp1,n,4
            dy(i) = dy(i) + da*dx(i)
            dy(i+1) = dy(i+1) + da*dx(i+1)
            dy(i+2) = dy(i+2) + da*dx(i+2)
            dy(i+3) = dy(i+3) + da*dx(i+3)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
          dy(iy) = dy(iy) + da*dx(ix)
          ix = ix + incx
          iy = iy + incy
         END DO
      END IF
      RETURN
      END

      SUBROUTINE dcopy(N,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
* Purpose:
* =============
*
*    DCOPY copies a vector, x, to a vector, y.
*    uses unrolled loops for increments equal to 1.
*
* Arguments:
* ==========
*
*    N is INTEGER number of elements in input vector(s)
*
*    DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*
*    INCX is INTEGER storage spacing between elements of DX
*
*    DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*
*    INCY is INTEGER storage spacing between elements of DY
*
* Further Details:
* =====================
*
*    jack dongarra, linpack, 3/11/78.
*    modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         m = mod(n,7)
         IF (m.NE.0) THEN
            DO i = 1,m
               dy(i) = dx(i)
            END DO
            IF (n.LT.7) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,7
            dy(i) = dx(i)
            dy(i+1) = dx(i+1)
            dy(i+2) = dx(i+2)
            dy(i+3) = dx(i+3)
            dy(i+4) = dx(i+4)
            dy(i+5) = dx(i+5)
            dy(i+6) = dx(i+6)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dy(iy) = dx(ix)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
      END

      DOUBLE PRECISION FUNCTION ddot(N,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
* Purpose:
* =============
*
*    DDOT forms the dot product of two vectors.
*    uses unrolled loops for increments equal to one.
*
* Arguments:
* ==========
*
*    N is INTEGER number of elements in input vector(s)
*
*    DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*
*    INCX is INTEGER storage spacing between elements of DX
*
*    DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*
*    INCY is INTEGER storage spacing between elements of DY
*
* Further Details:
* =====================
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      ddot = 0.0d0
      dtemp = 0.0d0
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         m = mod(n,5)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dtemp + dx(i)*dy(i)
            END DO
            IF (n.LT.5) THEN
               ddot=dtemp
            RETURN
            END IF
         END IF
         mp1 = m + 1
         DO i = mp1,n,5
          dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) +
     $            dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      ddot = dtemp
      RETURN
      END

      SUBROUTINE dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of
*  Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
*  Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(lda,*),B(ldb,*),C(ldc,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL NOTA,NOTB
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are
*     not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of
*     rows
*     and  columns of  A  and the  number of  rows  of  B
*     respectively.
*
      nota = lsame(transa,'N')
      notb = lsame(transb,'N')
      IF (nota) THEN
          nrowa = m
          ncola = k
      ELSE
          nrowa = k
          ncola = m
      END IF
      IF (notb) THEN
          nrowb = k
      ELSE
          nrowb = n
      END IF
*
*     Test the input parameters.
*
      info = 0
      IF ((.NOT.nota) .AND. (.NOT.lsame(transa,'C')) .AND.
     +    (.NOT.lsame(transa,'T'))) THEN
          info = 1
      ELSE IF ((.NOT.notb) .AND. (.NOT.lsame(transb,'C')) .AND.
     +         (.NOT.lsame(transb,'T'))) THEN
          info = 2
      ELSE IF (m.LT.0) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (k.LT.0) THEN
          info = 5
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 8
      ELSE IF (ldb.LT.max(1,nrowb)) THEN
          info = 10
      ELSE IF (ldc.LT.max(1,m)) THEN
          info = 13
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DGEMM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
     +    (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
*
*     And if  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          IF (beta.EQ.zero) THEN
              DO 20 j = 1,n
                  DO 10 i = 1,m
                      c(i,j) = zero
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 j = 1,n
                  DO 30 i = 1,m
                      c(i,j) = beta*c(i,j)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (notb) THEN
          IF (nota) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 50 i = 1,m
                          c(i,j) = zero
   50                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 60 i = 1,m
                          c(i,j) = beta*c(i,j)
   60                 CONTINUE
                  END IF
                  DO 80 l = 1,k
                      temp = alpha*b(l,j)
                      DO 70 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B + beta*C
*
              DO 120 j = 1,n
                  DO 110 i = 1,m
                      temp = zero
                      DO 100 l = 1,k
                          temp = temp + a(l,i)*b(l,j)
  100                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (nota) THEN
*
*           Form  C := alpha*A*B**T + beta*C
*
              DO 170 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 130 i = 1,m
                          c(i,j) = zero
  130                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 140 i = 1,m
                          c(i,j) = beta*c(i,j)
  140                 CONTINUE
                  END IF
                  DO 160 l = 1,k
                      temp = alpha*b(j,l)
                      DO 150 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B**T + beta*C
*
              DO 200 j = 1,n
                  DO 190 i = 1,m
                      temp = zero
                      DO 180 l = 1,k
                          temp = temp + a(l,i)*b(j,l)
  180                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END

      SUBROUTINE dgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(lda,*),X(*),Y(*)
*     ..
*
* Purpose:
* =============
*
* DGEMV  performs one of the matrix-vector operations
*
*    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
*
* where alpha and beta are scalars, x and y are vectors and A is an
* m by n matrix.
*
* Arguments:
* ==========
*
*    TRANS is CHARACTER*1
*     On entry, TRANS specifies the operation to be performed as
*     follows:
*
*        TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*        TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
*
*        TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
*    M is INTEGER
*     On entry, M specifies the number of rows of the matrix A.
*     M must be at least zero.
*
*    N is INTEGER
*     On entry, N specifies the number of columns of the matrix A.
*     N must be at least zero.
*
*    ALPHA is DOUBLE PRECISION.
*     On entry, ALPHA specifies the scalar alpha.
*
*    A is DOUBLE PRECISION array, dimension ( LDA, N )
*     Before entry, the leading m by n part of the array A must
*     contain the matrix of coefficients.
*
*    LDA is INTEGER
*     On entry, LDA specifies the first dimension of A as declared
*     in the calling (sub) program. LDA must be at least
*     max( 1, m ).
*
*    X is DOUBLE PRECISION array, dimension at least
*     ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*     and at least
*     ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*     Before entry, the incremented array X must contain the
*     vector x.
*
*    INCX is INTEGER
*     On entry, INCX specifies the increment for the elements of
*     X. INCX must not be zero.
*
*    BETA is DOUBLE PRECISION.
*     On entry, BETA specifies the scalar beta. When BETA is
*     supplied as zero then Y need not be set on input.
*
*    Y is DOUBLE PRECISION array, dimension at least
*     ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*     and at least
*     ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*     Before entry with BETA non-zero, the incremented array Y
*     must contain the vector y. On exit, Y is overwritten by the
*     updated vector y.
*
*    INCY is INTEGER
*     On entry, INCY specifies the increment for the elements of
*     Y. INCY must not be zero.
*
* Further Details:
* =====================
*
*  Level 2 Blas routine.
*  The vector and matrix arguments are not referenced when N = 0, or M = 0
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND.
     +    .NOT.lsame(trans,'C')) THEN
          info = 1
      ELSE IF (m.LT.0) THEN
          info = 2
      ELSE IF (n.LT.0) THEN
          info = 3
      ELSE IF (lda.LT.max(1,m)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      ELSE IF (incy.EQ.0) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DGEMV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
     +    ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (lsame(trans,'N')) THEN
          lenx = n
          leny = m
      ELSE
          lenx = m
          leny = n
      END IF
      IF (incx.GT.0) THEN
          kx = 1
      ELSE
          kx = 1 - (lenx-1)*incx
      END IF
      IF (incy.GT.0) THEN
          ky = 1
      ELSE
          ky = 1 - (leny-1)*incy
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (beta.NE.one) THEN
          IF (incy.EQ.1) THEN
              IF (beta.EQ.zero) THEN
                  DO 10 i = 1,leny
                      y(i) = zero
   10             CONTINUE
              ELSE
                  DO 20 i = 1,leny
                      y(i) = beta*y(i)
   20             CONTINUE
              END IF
          ELSE
              iy = ky
              IF (beta.EQ.zero) THEN
                  DO 30 i = 1,leny
                      y(iy) = zero
                      iy = iy + incy
   30             CONTINUE
              ELSE
                  DO 40 i = 1,leny
                      y(iy) = beta*y(iy)
                      iy = iy + incy
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (alpha.EQ.zero) RETURN
      IF (lsame(trans,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          jx = kx
          IF (incy.EQ.1) THEN
              DO 60 j = 1,n
                  temp = alpha*x(jx)
                  DO 50 i = 1,m
                      y(i) = y(i) + temp*a(i,j)
   50             CONTINUE
                  jx = jx + incx
   60         CONTINUE
          ELSE
              DO 80 j = 1,n
                  temp = alpha*x(jx)
                  iy = ky
                  DO 70 i = 1,m
                      y(iy) = y(iy) + temp*a(i,j)
                      iy = iy + incy
   70             CONTINUE
                  jx = jx + incx
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A**T*x + y.
*
          jy = ky
          IF (incx.EQ.1) THEN
              DO 100 j = 1,n
                  temp = zero
                  DO 90 i = 1,m
                      temp = temp + a(i,j)*x(i)
   90             CONTINUE
                  y(jy) = y(jy) + alpha*temp
                  jy = jy + incy
  100         CONTINUE
          ELSE
              DO 120 j = 1,n
                  temp = zero
                  ix = kx
                  DO 110 i = 1,m
                      temp = temp + a(i,j)*x(ix)
                      ix = ix + incx
  110             CONTINUE
                  y(jy) = y(jy) + alpha*temp
                  jy = jy + incy
  120         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END

      SUBROUTINE dger(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of
*  Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
*  Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(lda,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      parameter(zero=0.0d+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (m.LT.0) THEN
          info = 1
      ELSE IF (n.LT.0) THEN
          info = 2
      ELSE IF (incx.EQ.0) THEN
          info = 5
      ELSE IF (incy.EQ.0) THEN
          info = 7
      ELSE IF (lda.LT.max(1,m)) THEN
          info = 9
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DGER  ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (incy.GT.0) THEN
          jy = 1
      ELSE
          jy = 1 - (n-1)*incy
      END IF
      IF (incx.EQ.1) THEN
          DO 20 j = 1,n
              IF (y(jy).NE.zero) THEN
                  temp = alpha*y(jy)
                  DO 10 i = 1,m
                      a(i,j) = a(i,j) + x(i)*temp
   10             CONTINUE
              END IF
              jy = jy + incy
   20     CONTINUE
      ELSE
          IF (incx.GT.0) THEN
              kx = 1
          ELSE
              kx = 1 - (m-1)*incx
          END IF
          DO 40 j = 1,n
              IF (y(jy).NE.zero) THEN
                  temp = alpha*y(jy)
                  ix = kx
                  DO 30 i = 1,m
                      a(i,j) = a(i,j) + x(ix)*temp
                      ix = ix + incx
   30             CONTINUE
              END IF
              jy = jy + incy
   40     CONTINUE
      END IF
*
      RETURN
*
*     End of DGER  .
*
      END

      DOUBLE PRECISION FUNCTION dnrm2(N,X,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*
* Purpose:
* =============
*
* DNRM2 returns the euclidean norm of a vector via the function
* name, so that
*
*    DNRM2 := sqrt( x'*x )
*
* Arguments:
* ==========
*
*    N is INTEGER number of elements in input vector(s)
*
*    X is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*
*    INCX is INTEGER storage spacing between elements of DX
*
* Further Details:
* =====================
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC abs,sqrt
*     ..
      IF (n.LT.1 .OR. incx.LT.1) THEN
          norm = zero
      ELSE IF (n.EQ.1) THEN
          norm = abs(x(1))
      ELSE
          scale = zero
          ssq = one
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 ix = 1,1 + (n-1)*incx,incx
              IF (x(ix).NE.zero) THEN
                  absxi = abs(x(ix))
                  IF (scale.LT.absxi) THEN
                      ssq = one + ssq* (scale/absxi)**2
                      scale = absxi
                  ELSE
                      ssq = ssq + (absxi/scale)**2
                  END IF
              END IF
   10     CONTINUE
          norm = scale*sqrt(ssq)
      END IF
*
      dnrm2 = norm
      RETURN
*
*     End of DNRM2.
*
      END

      SUBROUTINE drot(N,DX,INCX,DY,INCY,C,S)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of
*  Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
*  Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*       code for both increments equal to 1
*
         DO i = 1,n
            dtemp = c*dx(i) + s*dy(i)
            dy(i) = c*dy(i) - s*dx(i)
            dx(i) = dtemp
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = c*dx(ix) + s*dy(iy)
            dy(iy) = c*dy(iy) - s*dx(ix)
            dx(ix) = dtemp
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
      END

      SUBROUTINE dscal(N,DA,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
* Purpose:
* =============
*
*    DSCAL scales a vector by a constant.
*    uses unrolled loops for increment equal to 1.
*
* Arguments:
* ==========
*
*    N is INTEGER number of elements in input vector(s)
*
*    DA is DOUBLE PRECISION On entry, DA specifies the scalar alpha.
*
*    DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*
*    INCX is INTEGER storage spacing between elements of DX
*
* Further Details:
* =====================
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
         m = mod(n,5)
         IF (m.NE.0) THEN
            DO i = 1,m
               dx(i) = da*dx(i)
            END DO
            IF (n.LT.5) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,5
            dx(i) = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
            dx(i+4) = da*dx(i+4)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            dx(i) = da*dx(i)
         END DO
      END IF
      RETURN
      END

      SUBROUTINE dswap(N,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
* Purpose:
* =============
*
*    DSWAP interchanges two vectors.
*    uses unrolled loops for increments equal to 1.
*
* Arguments:
* ==========
*
*    N is INTEGER number of elements in input vector(s)
*
*    DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*
*    INCX is INTEGER storage spacing between elements of DX
*
*    DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*
*    INCY is INTEGER storage spacing between elements of DY
*
* Further Details:
* =====================
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*       code for both increments equal to 1
*
*
*       clean-up loop
*
         m = mod(n,3)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dx(i)
               dx(i) = dy(i)
               dy(i) = dtemp
            END DO
            IF (n.LT.3) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,3
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
            dtemp = dx(i+1)
            dx(i+1) = dy(i+1)
            dy(i+1) = dtemp
            dtemp = dx(i+2)
            dx(i+2) = dy(i+2)
            dy(i+2) = dtemp
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = dx(ix)
            dx(ix) = dy(iy)
            dy(iy) = dtemp
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
      END

      SUBROUTINE dtrmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of
*  Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
*  Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(lda,*),B(ldb,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Test the input parameters.
*
      lside = lsame(side,'L')
      IF (lside) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      nounit = lsame(diag,'N')
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF ((.NOT.lsame(transa,'N')) .AND.
     +         (.NOT.lsame(transa,'T')) .AND.
     +         (.NOT.lsame(transa,'C'))) THEN
          info = 3
      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
          info = 4
      ELSE IF (m.LT.0) THEN
          info = 5
      ELSE IF (n.LT.0) THEN
          info = 6
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DTRMM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (m.EQ.0 .OR. n.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          DO 20 j = 1,n
              DO 10 i = 1,m
                  b(i,j) = zero
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lside) THEN
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*A*B.
*
              IF (upper) THEN
                  DO 50 j = 1,n
                      DO 40 k = 1,m
                          IF (b(k,j).NE.zero) THEN
                              temp = alpha*b(k,j)
                              DO 30 i = 1,k - 1
                                  b(i,j) = b(i,j) + temp*a(i,k)
   30                         CONTINUE
                              IF (nounit) temp = temp*a(k,k)
                              b(k,j) = temp
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 j = 1,n
                      DO 70 k = m,1,-1
                          IF (b(k,j).NE.zero) THEN
                              temp = alpha*b(k,j)
                              b(k,j) = temp
                              IF (nounit) b(k,j) = b(k,j)*a(k,k)
                              DO 60 i = k + 1,m
                                  b(i,j) = b(i,j) + temp*a(i,k)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*A**T*B.
*
              IF (upper) THEN
                  DO 110 j = 1,n
                      DO 100 i = m,1,-1
                          temp = b(i,j)
                          IF (nounit) temp = temp*a(i,i)
                          DO 90 k = 1,i - 1
                              temp = temp + a(k,i)*b(k,j)
   90                     CONTINUE
                          b(i,j) = alpha*temp
  100                 CONTINUE
  110             CONTINUE
              ELSE
                  DO 140 j = 1,n
                      DO 130 i = 1,m
                          temp = b(i,j)
                          IF (nounit) temp = temp*a(i,i)
                          DO 120 k = i + 1,m
                              temp = temp + a(k,i)*b(k,j)
  120                     CONTINUE
                          b(i,j) = alpha*temp
  130                 CONTINUE
  140             CONTINUE
              END IF
          END IF
      ELSE
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*B*A.
*
              IF (upper) THEN
                  DO 180 j = n,1,-1
                      temp = alpha
                      IF (nounit) temp = temp*a(j,j)
                      DO 150 i = 1,m
                          b(i,j) = temp*b(i,j)
  150                 CONTINUE
                      DO 170 k = 1,j - 1
                          IF (a(k,j).NE.zero) THEN
                              temp = alpha*a(k,j)
                              DO 160 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  160                         CONTINUE
                          END IF
  170                 CONTINUE
  180             CONTINUE
              ELSE
                  DO 220 j = 1,n
                      temp = alpha
                      IF (nounit) temp = temp*a(j,j)
                      DO 190 i = 1,m
                          b(i,j) = temp*b(i,j)
  190                 CONTINUE
                      DO 210 k = j + 1,n
                          IF (a(k,j).NE.zero) THEN
                              temp = alpha*a(k,j)
                              DO 200 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  200                         CONTINUE
                          END IF
  210                 CONTINUE
  220             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*A**T.
*
              IF (upper) THEN
                  DO 260 k = 1,n
                      DO 240 j = 1,k - 1
                          IF (a(j,k).NE.zero) THEN
                              temp = alpha*a(j,k)
                              DO 230 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      temp = alpha
                      IF (nounit) temp = temp*a(k,k)
                      IF (temp.NE.one) THEN
                          DO 250 i = 1,m
                              b(i,k) = temp*b(i,k)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              ELSE
                  DO 300 k = n,1,-1
                      DO 280 j = k + 1,n
                          IF (a(j,k).NE.zero) THEN
                              temp = alpha*a(j,k)
                              DO 270 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  270                         CONTINUE
                          END IF
  280                 CONTINUE
                      temp = alpha
                      IF (nounit) temp = temp*a(k,k)
                      IF (temp.NE.one) THEN
                          DO 290 i = 1,m
                              b(i,k) = temp*b(i,k)
  290                     CONTINUE
                      END IF
  300             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of DTRMM .
*
      END

      SUBROUTINE dtrmv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of
*  Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
*  Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(lda,*),X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      parameter(zero=0.0d+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
          info = 1
      ELSE IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND.
     +         .NOT.lsame(trans,'C')) THEN
          info = 2
      ELSE IF (.NOT.lsame(diag,'U') .AND. .NOT.lsame(diag,'N')) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,n)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DTRMV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (n.EQ.0) RETURN
*
      nounit = lsame(diag,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (incx.LE.0) THEN
          kx = 1 - (n-1)*incx
      ELSE IF (incx.NE.1) THEN
          kx = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (lsame(trans,'N')) THEN
*
*        Form  x := A*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 20 j = 1,n
                      IF (x(j).NE.zero) THEN
                          temp = x(j)
                          DO 10 i = 1,j - 1
                              x(i) = x(i) + temp*a(i,j)
   10                     CONTINUE
                          IF (nounit) x(j) = x(j)*a(j,j)
                      END IF
   20             CONTINUE
              ELSE
                  jx = kx
                  DO 40 j = 1,n
                      IF (x(jx).NE.zero) THEN
                          temp = x(jx)
                          ix = kx
                          DO 30 i = 1,j - 1
                              x(ix) = x(ix) + temp*a(i,j)
                              ix = ix + incx
   30                     CONTINUE
                          IF (nounit) x(jx) = x(jx)*a(j,j)
                      END IF
                      jx = jx + incx
   40             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 60 j = n,1,-1
                      IF (x(j).NE.zero) THEN
                          temp = x(j)
                          DO 50 i = n,j + 1,-1
                              x(i) = x(i) + temp*a(i,j)
   50                     CONTINUE
                          IF (nounit) x(j) = x(j)*a(j,j)
                      END IF
   60             CONTINUE
              ELSE
                  kx = kx + (n-1)*incx
                  jx = kx
                  DO 80 j = n,1,-1
                      IF (x(jx).NE.zero) THEN
                          temp = x(jx)
                          ix = kx
                          DO 70 i = n,j + 1,-1
                              x(ix) = x(ix) + temp*a(i,j)
                              ix = ix - incx
   70                     CONTINUE
                          IF (nounit) x(jx) = x(jx)*a(j,j)
                      END IF
                      jx = jx - incx
   80             CONTINUE
              END IF
          END IF
      ELSE
*
*        Form  x := A**T*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 100 j = n,1,-1
                      temp = x(j)
                      IF (nounit) temp = temp*a(j,j)
                      DO 90 i = j - 1,1,-1
                          temp = temp + a(i,j)*x(i)
   90                 CONTINUE
                      x(j) = temp
  100             CONTINUE
              ELSE
                  jx = kx + (n-1)*incx
                  DO 120 j = n,1,-1
                      temp = x(jx)
                      ix = jx
                      IF (nounit) temp = temp*a(j,j)
                      DO 110 i = j - 1,1,-1
                          ix = ix - incx
                          temp = temp + a(i,j)*x(ix)
  110                 CONTINUE
                      x(jx) = temp
                      jx = jx - incx
  120             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 140 j = 1,n
                      temp = x(j)
                      IF (nounit) temp = temp*a(j,j)
                      DO 130 i = j + 1,n
                          temp = temp + a(i,j)*x(i)
  130                 CONTINUE
                      x(j) = temp
  140             CONTINUE
              ELSE
                  jx = kx
                  DO 160 j = 1,n
                      temp = x(jx)
                      ix = jx
                      IF (nounit) temp = temp*a(j,j)
                      DO 150 i = j + 1,n
                          ix = ix + incx
                          temp = temp + a(i,j)*x(ix)
  150                 CONTINUE
                      x(jx) = temp
                      jx = jx + incx
  160             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of DTRMV .
*
      END

      SUBROUTINE dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(lda,*),B(ldb,*)
*     ..
*
* Purpose:
* =============
*
* DTRSM  solves one of the matrix equations
*
*    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
* where alpha is a scalar, X and B are m by n matrices, A is a unit, or
* non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*    op( A ) = A   or   op( A ) = A**T.
*
* The matrix X is overwritten on B.
*
* Arguments:
* ==========
*
*    SIDE is CHARACTER*1
*     On entry, SIDE specifies whether op( A ) appears on the left
*     or right of X as follows:
*
*        SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*        SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*    UPLO is CHARACTER*1
*     On entry, UPLO specifies whether the matrix A is an upper or
*     lower triangular matrix as follows:
*
*        UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*        UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*    TRANSA is CHARACTER*1
*     On entry, TRANSA specifies the form of op( A ) to be used in
*     the matrix multiplication as follows:
*
*        TRANSA = 'N' or 'n'   op( A ) = A.
*
*        TRANSA = 'T' or 't'   op( A ) = A**T.
*
*        TRANSA = 'C' or 'c'   op( A ) = A**T.
*
*    DIAG is CHARACTER*1
*     On entry, DIAG specifies whether or not A is unit triangular
*     as follows:
*
*        DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*        DIAG = 'N' or 'n'   A is not assumed to be unit
*                            triangular.
*
*    M is INTEGER
*     On entry, M specifies the number of rows of B. M must be at
*     least zero.
*
*    N is INTEGER
*     On entry, N specifies the number of columns of B.  N must be
*     at least zero.
*
*    ALPHA is DOUBLE PRECISION.
*     On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*     zero then  A is not referenced and  B need not be set before
*     entry.
*
*    A is DOUBLE PRECISION array, dimension ( LDA, k ),
*     where k is m when SIDE = 'L' or 'l'
*       and k is n when SIDE = 'R' or 'r'.
*     Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*     upper triangular part of the array  A must contain the upper
*     triangular matrix  and the strictly lower triangular part of
*     A is not referenced.
*     Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*     lower triangular part of the array  A must contain the lower
*     triangular matrix  and the strictly upper triangular part of
*     A is not referenced.
*     Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*     A  are not referenced either,  but are assumed to be  unity.
*
*    LDA is INTEGER
*     On entry, LDA specifies the first dimension of A as declared
*     in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*     LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*     then LDA must be at least max( 1, n ).
*
*    B is DOUBLE PRECISION array, dimension ( LDB, N )
*     Before entry,  the leading  m by n part of the array  B must
*     contain  the  right-hand  side  matrix  B,  and  on exit  is
*     overwritten by the solution matrix  X.
*
*    LDB is INTEGER
*     On entry, LDB specifies the first dimension of B as declared
*     in  the  calling  (sub)  program.   LDB  must  be  at  least
*     max( 1, m ).
*
* Further Details:
* =====================
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC max
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*
*     Test the input parameters.
*
      lside = lsame(side,'L')
      IF (lside) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      nounit = lsame(diag,'N')
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF ((.NOT.lsame(transa,'N')) .AND.
     +         (.NOT.lsame(transa,'T')) .AND.
     +         (.NOT.lsame(transa,'C'))) THEN
          info = 3
      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
          info = 4
      ELSE IF (m.LT.0) THEN
          info = 5
      ELSE IF (n.LT.0) THEN
          info = 6
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DTRSM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (m.EQ.0 .OR. n.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          DO 20 j = 1,n
              DO 10 i = 1,m
                  b(i,j) = zero
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lside) THEN
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
              IF (upper) THEN
                  DO 60 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 30 i = 1,m
                              b(i,j) = alpha*b(i,j)
   30                     CONTINUE
                      END IF
                      DO 50 k = m,1,-1
                          IF (b(k,j).NE.zero) THEN
                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
                              DO 40 i = 1,k - 1
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 70 i = 1,m
                              b(i,j) = alpha*b(i,j)
   70                     CONTINUE
                      END IF
                      DO 90 k = 1,m
                          IF (b(k,j).NE.zero) THEN
                              IF (nounit) b(k,j) = b(k,j)/a(k,k)
                              DO 80 i = k + 1,m
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*inv( A**T )*B.
*
              IF (upper) THEN
                  DO 130 j = 1,n
                      DO 120 i = 1,m
                          temp = alpha*b(i,j)
                          DO 110 k = 1,i - 1
                              temp = temp - a(k,i)*b(k,j)
  110                     CONTINUE
                          IF (nounit) temp = temp/a(i,i)
                          b(i,j) = temp
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 j = 1,n
                      DO 150 i = m,1,-1
                          temp = alpha*b(i,j)
                          DO 140 k = i + 1,m
                              temp = temp - a(k,i)*b(k,j)
  140                     CONTINUE
                          IF (nounit) temp = temp/a(i,i)
                          b(i,j) = temp
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
              IF (upper) THEN
                  DO 210 j = 1,n
                      IF (alpha.NE.one) THEN
                          DO 170 i = 1,m
                              b(i,j) = alpha*b(i,j)
  170                     CONTINUE
                      END IF
                      DO 190 k = 1,j - 1
                          IF (a(k,j).NE.zero) THEN
                              DO 180 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (nounit) THEN
                          temp = one/a(j,j)
                          DO 200 i = 1,m
                              b(i,j) = temp*b(i,j)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 j = n,1,-1
                      IF (alpha.NE.one) THEN
                          DO 220 i = 1,m
                              b(i,j) = alpha*b(i,j)
  220                     CONTINUE
                      END IF
                      DO 240 k = j + 1,n
                          IF (a(k,j).NE.zero) THEN
                              DO 230 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (nounit) THEN
                          temp = one/a(j,j)
                          DO 250 i = 1,m
                              b(i,j) = temp*b(i,j)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*inv( A**T ).
*
              IF (upper) THEN
                  DO 310 k = n,1,-1
                      IF (nounit) THEN
                          temp = one/a(k,k)
                          DO 270 i = 1,m
                              b(i,k) = temp*b(i,k)
  270                     CONTINUE
                      END IF
                      DO 290 j = 1,k - 1
                          IF (a(j,k).NE.zero) THEN
                              temp = a(j,k)
                              DO 280 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (alpha.NE.one) THEN
                          DO 300 i = 1,m
                              b(i,k) = alpha*b(i,k)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 k = 1,n
                      IF (nounit) THEN
                          temp = one/a(k,k)
                          DO 320 i = 1,m
                              b(i,k) = temp*b(i,k)
  320                     CONTINUE
                      END IF
                      DO 340 j = k + 1,n
                          IF (a(j,k).NE.zero) THEN
                              temp = a(j,k)
                              DO 330 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (alpha.NE.one) THEN
                          DO 350 i = 1,m
                              b(i,k) = alpha*b(i,k)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END

      INTEGER FUNCTION idamax(N,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
* Purpose:
* =============
*
*    IDAMAX finds the index of the first element having maximum absolute value.
*
* Arguments:
* ==========
*
*    N is INTEGER number of elements in input vector(s)
*
*    DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*
*    INCX is INTEGER storage spacing between elements of SX
*
* Further Details:
* =====================
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dabs
*     ..
      idamax = 0
      IF (n.LT.1 .OR. incx.LE.0) RETURN
      idamax = 1
      IF (n.EQ.1) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         dmax = dabs(dx(1))
         DO i = 2,n
            IF (dabs(dx(i)).GT.dmax) THEN
               idamax = i
               dmax = dabs(dx(i))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         ix = 1
         dmax = dabs(dx(1))
         ix = ix + incx
         DO i = 2,n
            IF (dabs(dx(ix)).GT.dmax) THEN
               idamax = i
               dmax = dabs(dx(ix))
            END IF
            ix = ix + incx
         END DO
      END IF
      RETURN
      END

      LOGICAL FUNCTION lsame(CA,CB)
*
*  -- Reference BLAS level1 routine (version 3.1) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER CA,CB
*     ..
*
* Purpose:
* =============
*
* LSAME returns .TRUE. if CA is the same letter as CB regardless of
* case.
*
* Arguments:
* ==========
*
*    CA is CHARACTER*1
*    CB is CHARACTER*1
*    CA and CB specify the single characters to be compared.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ichar
*     ..
*     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
*     ..
*
*     Test if the characters are equal
*
      lsame = ca .EQ. cb
      IF (lsame) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      zcode = ichar('Z')
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      inta = ichar(ca)
      intb = ichar(cb)
*
      IF (zcode.EQ.90 .OR. zcode.EQ.122) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
          IF (inta.GE.97 .AND. inta.LE.122) inta = inta - 32
          IF (intb.GE.97 .AND. intb.LE.122) intb = intb - 32
*
      ELSE IF (zcode.EQ.233 .OR. zcode.EQ.169) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
          IF (inta.GE.129 .AND. inta.LE.137 .OR.
     +        inta.GE.145 .AND. inta.LE.153 .OR.
     +        inta.GE.162 .AND. inta.LE.169) inta = inta + 64
          IF (intb.GE.129 .AND. intb.LE.137 .OR.
     +        intb.GE.145 .AND. intb.LE.153 .OR.
     +        intb.GE.162 .AND. intb.LE.169) intb = intb + 64
*
      ELSE IF (zcode.EQ.218 .OR. zcode.EQ.250) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
          IF (inta.GE.225 .AND. inta.LE.250) inta = inta - 32
          IF (intb.GE.225 .AND. intb.LE.250) intb = intb - 32
      END IF
      lsame = inta .EQ. intb
*
*     RETURN
*
*     End of LSAME
*
      END

      SUBROUTINE xerbla( SRNAME, INFO )
*
*  -- Reference BLAS level1 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
*     ..
*
* Purpose:
* =============
*
* XERBLA  is an error handler for the LAPACK routines.
* It is called by an LAPACK routine if an input parameter has an
* invalid value.  A message is printed and execution stops.
*
* Installers may consider modifying the STOP statement in order to
* call system-specific exception-handling facilities.
*
* Arguments:
* ==========
*
*    SRNAME is CHARACTER*(*)
*    The name of the routine which called XERBLA.
*
*    INFO is INTEGER
*    The position of the invalid parameter in the parameter list
*    of the calling routine.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          len_trim
*     ..
*     .. Executable Statements ..
*
      WRITE( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
*
      stop
*
 9999 FORMAT( ' ** On entry to ', a, ' parameter number ', i2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END

