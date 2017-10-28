MODULE VTgeometry
! FORTRAN 2003 module for numerically stable geometry   
! computations in arbitrary dimension d                 
! Uses BLAS and LAPACK routines for linear equations    
! Compile w/ -llapack -lblas                            
! Tested w/ gfortran 4.8.5 on x86_64-redhat-linux       
! AUTHOR: Tyler Chang                                   
! Last Update: July, 2017                               

USE real_precision
IMPLICIT NONE

CONTAINS
! Subroutines:

SUBROUTINE CircleRad(d, n, list, pts, rad, center, err)
! Get the circumcircle for d+1 points in d dimensions.
! integer (in) d: number of dimensions                  
! integer (in) n: number of points                      
! integer (in) list(d+1): indices of d+1 points to use  
! double (in) pts(d,n): list of points                  
! integer (out) rad: radius of circumcircle             
! double (out) center(d): center of circumcircle        
! integer (out) error: flag, success = 0                

! inputs
INTEGER, INTENT(IN) :: d, n, list(d+1)
REAL(KIND=R8), INTENT(IN) :: pts(d,n)
! outputs
REAL(KIND=R8), INTENT(OUT) :: rad, center(d)
INTEGER, INTENT(OUT) :: err
! local vars
REAL(KIND=R8), ALLOCATABLE :: A(:,:)
INTEGER :: i
INTEGER, ALLOCATABLE :: ipiv(:)
! blas functions & lapack subroutines
REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
EXTERNAL :: DGESV
! allocate data
ALLOCATE(A(d,d), ipiv(d), STAT=err)
IF (err /= 0) RETURN
! build linear system Ac = b; where c = center,
! the rows of A are normal vectors v to the bisecting
! planes, and b is dot(v, x)
DO i = 2, d+1
        ! compute normal vectors and store in A
        A(i-1,:) = pts(:,list(i)) - pts(:,list(1))
        ! compute right hand side and store in center
        center(i-1) = DDOT(d, A(i-1,:), 1, A(i-1,:), 1) / REAL(-2.,R8)
END DO
! solve linear system to get center
CALL DGESV(d, 1, A, d, ipiv, center, d, err)
IF (err /= 0) RETURN
! get radius and center
rad = DNRM2(d, center, 1)
center = center + pts(:,list(1))
! free heap
DEALLOCATE(A, ipiv)
RETURN
END SUBROUTINE CircleRad

SUBROUTINE MakePlane(d, n, list, pts, plane, err)
! Given d points in d space, find the equation for the  
! intersecting plane.                                   
! integer (in) d : number of dimensions                 
! integer (in) n : total number of points               
! integer (in) list(d) : indices for (d) points         
!                        which define a plane           
! double (in) pts(d,n) : list of points                 
! double (out) plane(d+1) : coefficients for plane,     
!                           intercept in [d+1]          
! integer (out) err : error flag (0=success)            

! inputs
INTEGER, INTENT(IN) :: d, n, list(d)
REAL(KIND=R8), INTENT(IN) :: pts(d,n)
! outputs
REAL(KIND=R8), INTENT(OUT) :: plane(d+1)
INTEGER, INTENT(OUT) :: err
! data
REAL(KIND=R8), ALLOCATABLE :: A(:,:), work(:)
INTEGER :: i
INTEGER, ALLOCATABLE :: jpiv(:)
! blas functions & lapack subroutines
REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
EXTERNAL :: DGEQP3, DTRSM
! allocate data
ALLOCATE(A(d-1,d), work(3*d+1), jpiv(d), STAT=err)
IF(err /= 0) RETURN
! build a matrix of vectors in the plane
DO i = 2, d
        A(i-1,:) = pts(:,list(i)) - pts(:,list(1))
END DO
! get QR
jpiv=0
CALL DGEQP3(d-1,d,A,d-1,jpiv,plane,work,3*d+1,err)
IF(err /= 0) RETURN
! do triangular solve to get plane
! store in work array
work(1:d-1) = A(:,d)
CALL DTRSM('L','U','N','N',d-1,1,REAL(-1.,R8),A,d-1,work,d-1)
work(d) = REAL(1.,R8)
! undo pivots
DO i = 1,d
        plane(jpiv(i)) = work(i)
END DO
! normalize
plane(1:d) = plane(1:d) / DNRM2(d,plane(1:d),1)
! calculate intercept
plane(d+1) = DDOT(d,plane(1:d),1,pts(:,list(1)),1)
! free heap
DEALLOCATE(A, work, jpiv)
RETURN
END SUBROUTINE MakePlane

FUNCTION PtInSimp(d, n, simp, pts, p, weights, err) RESULT(tf)
! Determine whether a simplex contains p.               
! integer (in) d : number of dimensions                 
! integer (in) n : total number of points               
! integer (in) simp(d+1) : indices for (d+1) points     
!                        which define a simplex         
! double (in) pts(d,n) : list of points                 
! double (in) p(d) : point to check for                 
! double (out) weights(d+1) : interpolation weights       
!                           associated with vertices.
!                           When tf = TRUE,     
!                           weights(k) > 0 for all k    
! integer (out) err : error flag (0=success)            
! logical (result) tf : is p in s?                      

! inputs
INTEGER, INTENT(IN) :: d, n, simp(d+1)
REAL(KIND=R8), INTENT(IN) :: pts(d,n), p(d)
! outputs
REAL(KIND=R8), INTENT(OUT) :: weights(d+1)
INTEGER, INTENT(OUT) :: err
! return value
LOGICAL :: tf
! local vars
REAL(KIND=R8), ALLOCATABLE :: A(:,:), b(:)
INTEGER :: i
INTEGER, ALLOCATABLE :: ipiv(:)
! lapack subroutines
EXTERNAL :: DGESV
! init result
tf = .FALSE.
! allocate data
ALLOCATE(A(d,d), ipiv(d), STAT=err)
IF (err /= 0) RETURN
! build linear system to get weights
DO i = 1, d
        A(1:d, i) = pts(:,simp(i+1)) - pts(:,simp(1))
END DO
weights(2:d+1) = p(:) - pts(:,simp(1))
! solve linear system to get weights
CALL DGESV(d, 1, A, d, ipiv, weights(2:d+1), d, err)
IF (err /= 0) RETURN
weights(1) = REAL(1.,R8) - SUM(weights(2:d+1))
! check if weights define a convex combination
IF (ALL(weights >= REAL(0.,R8))) tf = .TRUE.
! free heap
DEALLOCATE(A, ipiv)
RETURN
END FUNCTION PtInSimp

SUBROUTINE Project(d, p, q, v, proj, dist)
! Project a point onto a plane defined by a normal
! vector and point in the plane.
! integer (in) d : number of dimensions                 
! double (in) p(d) : point to project onto plane
! double (in) q(d) : point in plane
! double (in) v(d) : unit length vector normal to plane
! double (out) dist : euclidean distance between p and
!                     its projection on the plane
! double (out) proj : projection of p on the plane

! inputs
INTEGER, INTENT(IN) :: d
REAL(KIND=R8), INTENT(IN) :: p(d), q(d), v(d)
! outputs
REAL(KIND=R8), INTENT(OUT) :: proj(d), dist
! blas functions
REAL (KIND=R8), EXTERNAL :: DDOT, DNRM2
! calculate
proj = p - v * DDOT(d, p - q, 1, v, 1)
dist = DNRM2(d, p - proj, 1)
! done
RETURN
END SUBROUTINE Project

END MODULE VTgeometry
