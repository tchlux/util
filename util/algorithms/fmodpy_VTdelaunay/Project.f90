
SUBROUTINE Project(d, n, pts, p, proj, rnorm, err)
!
! Project a point outside the convex hull of a point set onto the convex
! hull by solving a constrained least squares problem. Uses WNNLS, as
! proposed/coded in ACM TOMS Algorithm #587.
! INTEGER (IN) d : Dimension of space.                                  
! INTEGER (IN) n : Number of points in set.                             
! DOUBLE (IN) pts(d,n) : Set of n points in d-space.                    
! DOUBLE (IN) p(d) : A single point in d-space, outside the convex hull.
! DOUBLE (OUT) proj(d) : The projection of p onto the convex hull of
!                        pts.
! DOUBLE (OUT) rnorm : The 2-norm of the residual vector: ||p-proj||_2.
! INTEGER (OUT) err : Error flag. SUCCESS = 0.
!
! AUTHOR: Tyler Chang                                   
! Last Update: September, 2017                               
!
USE REAL_PRECISION
IMPLICIT NONE

! inputs
INTEGER, INTENT(IN) :: d, n
REAL(KIND=R8), INTENT(IN) :: pts(d,n), p(d)
! outputs
REAL(KIND=R8) :: proj(d), rnorm
INTEGER, INTENT(OUT) :: err
! local vars
REAL(KIND=R8) :: A(d+1,n+1), prgopt(1), x(n), work(d+1+n*5)
INTEGER :: iwork(d+1+n), i, j
! blas subroutine
EXTERNAL :: DGEMV
! external subroutine from ACM: TOMS-587
EXTERNAL :: WNNLS

! initialize work arrays
prgopt(1) = 1.0_R8
iwork(1) = d+1+5*n
iwork(2) = d+1+n
! set convex condition
A(1, :) = 1.0_R8
! copy in points
A(2:d+1,1:n) = pts(:,:)
! create interpolation point
A(2:d+1,n+1) = p(:)
! get solution
CALL WNNLS(A,d+1,1,d,n,0,prgopt,x,rnorm,err,iwork,work)
IF (err /= 0) THEN
        err = 30
        RETURN
END IF
! compute proj
CALL DGEMV('N',d,n,1.0_R8,pts,d,x,1,0.0_R8,proj,1)

! done
RETURN
END SUBROUTINE Project

