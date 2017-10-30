
SUBROUTINE Extrapolate(d, n, pts, p, proj, rnorm, err)
USE REAL_PRECISION
IMPLICIT NONE
! Extrapolate outside the convex hull of a point set in d-space, by     
! projecting onto the convex hull of the point set.
! INTEGER (IN) d : Dimension of space.                                  
! INTEGER (IN) n : Number of points in set.                             
! DOUBLE (IN) pts(d,n) : Set of n points in d-space.                    
! DOUBLE (IN) p(d) : A single point in d-space, outside the convex hull.
! DOUBLE (OUT) proj(d) : The projection of p onto the convex hull of
!                        pts.
! DOUBLE (OUT) rnorm : The 2-norm of the residual vector: ||p-proj||_2.
! INTEGER (OUT) err : Error flag. SUCCESS = 0.

! inputs
INTEGER, INTENT(IN) :: d, n
REAL(KIND=R8), INTENT(IN) :: pts(d,n), p(d)
! outputs
REAL(KIND=R8) :: proj(d), rnorm
INTEGER, INTENT(OUT) :: err
! local vars
REAL(KIND=R8), ALLOCATABLE :: A(:,:), prgopt(:), x(:), work(:)
INTEGER, ALLOCATABLE :: iwork(:)
INTEGER :: i, j
! blas subroutine
EXTERNAL :: DGEMV
! external subroutine from ACM: TOMS-587
EXTERNAL :: WNNLS
! allocate
ALLOCATE(A(d+1,n+1), prgopt(1), x(n), work(d+1+n*5),&
        iwork(d+1+n), STAT=err)
IF (err /= 0) RETURN
! initialize work arrays
prgopt(1) = REAL(1.,R8)
iwork(1) = d+1+5*n
iwork(2) = d+1+n
! set convex condition
A(1, :) = REAL(1.,R8)
! copy in points
A(2:d+1,1:n) = pts(:,:)
! create interpolation point
A(2:d+1,n+1) = p(:)
! get solution
CALL WNNLS(A,d+1,1,d,n,0,prgopt,x,rnorm,err,iwork,work)
IF (err /= 0) RETURN
! compute proj
CALL DGEMV('N',d,n,REAL(1.,R8),pts,d,x,1,REAL(0.,R8),proj,1)
! free heap
DEALLOCATE(A, prgopt, X, work, iwork)
! done
RETURN
END SUBROUTINE Extrapolate

