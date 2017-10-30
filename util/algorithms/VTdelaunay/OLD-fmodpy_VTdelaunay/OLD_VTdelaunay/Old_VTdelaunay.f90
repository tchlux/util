
MODULE VTdelaunay
! FORTRAN 2003 module for computing a delaunay simplex
! containing a given point.
! Compile w/ -llapack -lblas                            
! Tested w/ gfortran 4.8.5 on x86_64-redhat-linux       
! AUTHOR: Tyler Chang                                   
! Last Update: August, 2017                               
USE real_precision

IMPLICIT NONE
PRIVATE
PUBLIC :: delaunayp

CONTAINS

SUBROUTINE delaunayp(d, n, pts, p, eps, simp, weights, err)
! Compute the simplex in the delaunay triangulation containing the point
! p. Inspired by DeWall algorithm. Compute a delaunay simplex "close to"
! p by solving a series of LS problems. Proceed to connect faces of the
! triangulation until a simplex is found containing p.                  
! INTEGER (IN) d : Dimension of space.                                  
! INTEGER (IN) n : Number of points in set.                             
! DOUBLE (IN) pts(d,n) : Set of n points in d-space.                    
! DOUBLE (IN) p(d) : A single point in d-space to interpolate using     
!                    the delaunay triangulation.                        
! DOUBLE (IN) eps : A small number defining the tolerance for accepting 
!                   "nearly" containing simplices. 0 < eps <<< 1.
!                   Recommended value: EPSILON ** (1.0 / DOUBLE(d))
! INTEGER (OUT) simp(d+1) : The indices of the points which make up the 
!                           vertices of the simplex containing p.       
! DOUBLE (OUT) weights(d+1) : The weights for the vertices of simp for  
!                             interpolating at p.                       
! INTEGER (OUT) err : Error flag. SUCCESS = 0.

! inputs
INTEGER, INTENT(IN) :: d, n
REAL(KIND=R8), INTENT(IN) :: pts(d,n), p(d), eps
! outputs
INTEGER, INTENT(OUT) :: simp(d+1)
REAL(KIND=R8), INTENT(OUT) :: weights(d+1)
INTEGER, INTENT(OUT) :: err
! local vars
REAL(KIND=R8) :: currRad, minRad, side1, side2, rnorm
REAL(KIND=R8), ALLOCATABLE :: A(:,:), b(:), currCenter(:), LQ(:,:), &
        minCenter(:), plane(:), proj(:), work(:)
INTEGER :: budget, i, j, k, tmp
INTEGER, ALLOCATABLE :: piv(:)
! blas functions
REAL(KIND=R8), EXTERNAL :: DDOT
! external subroutine
EXTERNAL :: Project
! allocate data
ALLOCATE(A(d,d), b(d), currCenter(d), LQ(d,d), minCenter(d), piv(d), &
        plane(d+1), proj(d), work(5*d), STAT=err)
IF (err /= 0) RETURN
budget = 1000
! initialize the projection
proj = p
rnorm = REAL(0., R8)
! construct a delaunay simplex "near" p
CALL MakeFirstSimp()
IF (err /= 0 .OR. simp(d+1) == 0) RETURN
! begin search for simplex containing p
DO k = 1, budget
        ! check if contains p
        IF (PtInSimp()) EXIT
        ! swap out least weighted vertex
        i = MINLOC(weights(1:d+1), DIM=1)
        tmp = simp(i)
        simp(i) = simp(d+1)
        ! if i /= 1, simply drop a row from linear system
        IF(i /= 1) THEN
                A(i-1,:) = A(d,:)
                b(i-1) = b(d)
        ! if i == 1, must reconstruct A & b
        ELSE
                DO j=1,d
                        A(j,:) = pts(:,simp(1)) - pts(:,simp(j+1))
                        b(j) = DDOT(d, A(j,:), 1, A(j,:), 1) / REAL(-2.,R8) 
                END DO
        END IF
        ! compute next simplex
        CALL MakeSimplex()
        IF (err /= 0) RETURN
        ! check for extrapolation condition
        IF (simp(d+1) == 0) THEN
                ! project p onto convex hull
                CALL Project(d, n, pts, p, proj, rnorm, err)
                IF (err /= 0) RETURN
                ! remake previous simplex
                simp(d+1) = tmp
                A(d,:) = pts(:,simp(1)) - pts(:,tmp)
                b(d) = DDOT(d, A(d,:), 1, A(d,:), 1) / REAL(-2.,R8)
        END IF
END DO
! set error flag
IF (rnorm > eps) err = 1
IF (k > budget) err = 2
! free heap
DEALLOCATE(A, b, currCenter, LQ, minCenter, piv, plane, proj, work)
! done
RETURN
! end subroutine

! internal subroutines & functions
CONTAINS

SUBROUTINE MakeFirstSimp()
! Build first face by solving sequence of LS problems ranging from 2Xd  
! up to [d-1]Xd. At each step the center is given by the minimum norm   
! solution to Ac = b. Then the radius is given by r = ||c||_2.          
! A^(m) = [ x2 - x1 | x3 - x1 | ... | xm - x1 ] where x1 - xm-1 are the 
! current points in the face. The xm which results in the minimum r is  
! added to A^(m+1) until an entire face (d points) has been constructed.

! blas functions & lapack subroutines
REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
EXTERNAL :: DGELSS
! find first point
simp = 0
minRad = HUGE(minRad)
DO i = 1, n
        ! check distance to p
        currRad = DNRM2(d, pts(:,i) - proj(:), 1)
        IF (currRad < minRad) THEN
                minRad = currRad
                simp(1) = i
        END IF
END DO
! check that a point was found
IF (simp(1) == 0) RETURN
! find second point
minRad = HUGE(minRad)
DO i = 1, n
        ! avoid repeats
        IF (i == simp(1)) CYCLE
        ! check radius
        currRad = DNRM2(d, pts(:,i)-pts(:,simp(1)), 1)
        IF (currRad < minRad) THEN
                minRad = currRad
                simp(2) = i
        END IF
END DO
! check that a point was found
IF (simp(2) == 0) RETURN
! set up first row of LS system
A(1,:) = pts(:,simp(1)) - pts(:,simp(2))
b(1) = DDOT(d, A(1,:), 1, A(1,:), 1) / REAL(-2.,R8)
! loop over remaining d-1 points in first simplex
DO i = 2, d
        ! re-init radius for each iteration
        minRad = HUGE(minRad)
        ! find next point to add
        DO j = 1, n
                ! check point is not already in face
                IF (ANY(simp == j)) CYCLE
                ! add curr point to LS system
                A(i,:) = pts(:,simp(1)) - pts(:,j)
                b(i) = DDOT(d, A(i,:), 1, A(i,:), 1) / REAL(-2.,R8)
                ! solve LS problem for min norm solution using SVD
                LQ(1:i,:) = A(1:i,:)
                currCenter(1:i) = b(1:i)
                CALL DGELSS(i, d, 1, LQ, d, currCenter, d, plane, eps, &
                        tmp, work, 5*d, err)
                IF (err /= 0) RETURN
                ! if rank deficient, cycle
                IF (tmp < i) CYCLE
                ! calculate radius
                currRad = DNRM2(d, currCenter, 1)
                ! check for new minimum
                IF (currRad < minRad) THEN
                        minRad = currRad
                        simp(i+1) = j
                END IF
        END DO
        ! check that a point was found
        IF (simp(i+1) == 0) RETURN
        ! add next point to LS system
        A(i,:) = pts(:,simp(1)) - pts(:,simp(i+1))
        b(i) = DDOT(d, A(i,:), 1, A(i,:), 1) / REAL(-2.,R8) 
END DO
! done
err = 0
RETURN
END SUBROUTINE MakeFirstSimp

SUBROUTINE MakeSimplex()
! Given a face of a delaunay simplex, complete the simplex by adding a  
! point on the same side of the face as p. The new point must minimize  
! the radius of the circumcircle. Determine the radius by solving the   
! linear system Ac = b, then taking r = ||c||_2. Note that              
! A = [ x2 - x1 | x3 - x1 | ... | x(d+1) - x1 ] where x1 ... xd are the 
! points in the face and x(d+1) is a potential new point (completing the
! simplex). b_i = dot(A_i,A_i) / -2.                                    

! blas functions & lapack subroutines
REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
EXTERNAL :: DGESV
! calculate halfspace
CALL MakePlane()
IF(err /= 0) RETURN
! calculate side of the viable points
side1 =  DDOT(d,plane(1:d),1,proj(:),1) - plane(d+1)
! normalize magnitude of p
side1 = SIGN(REAL(1.,R8),side1)
! initialize center, radius, and simplex
simp(d+1) = 0
minCenter = REAL(0.,R8)
minRad = HUGE(minRad)
! loop through all points
DO i = 1, n
        ! check point(i) is inside current minimum
        IF (DNRM2(d, pts(:,i) - minCenter(:), 1) > minRad) CYCLE
        ! check point(i) on viable side of halfspace
        side2 = DDOT(d,plane(1:d),1,pts(:,i),1) - plane(d+1)
        IF (side1 * side2 < eps .OR. ANY(simp == i)) CYCLE
        ! add point i to linear system
        A(d,:) = pts(:,simp(1)) - pts(:,i)
        b(d) = DDOT(d, A(d,:), 1, A(d,:), 1) / REAL(-2.,R8)
        LQ = A
        currCenter = b
        ! solve linear system to get center
        CALL DGESV(d, 1, LQ, d, piv, currCenter, d, err)
        IF (err < 0) RETURN
        IF (err > 0) CYCLE
        currRad = DNRM2(d, currCenter, 1)
        ! check if new minimum
        IF (currRad < minRad) THEN
                minRad = currRad
                minCenter = currCenter + pts(:,simp(1))
                simp(d+1) = i
        END IF
END DO
! reset error flag
err = 0
! check for extrapolation condition
IF(simp(d+1) == 0) RETURN
! otherwise add point to linear system
A(d,:) = pts(:,simp(1)) - pts(:,simp(d+1))
b(d) = DDOT(d, A(d,:), 1, A(d,:), 1) / REAL(-2.,R8)
! done
RETURN
END SUBROUTINE MakeSimplex

SUBROUTINE MakePlane()
! Construct a plane intersecting the current d points. The coefficients 
! are determined by the normal vector: v and an intercept, stored in    
! plane(1:d) and plane(d+1) respectively. To find the normal vector, the
! current d points are transformed into d-1 vectors in the plane.       
! Together, these d-1 d-vectors have nullity rank-1. The normal vector  
! is any nontrivial vector in the nullspace.                            

! blas functions & lapack subroutines
REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
EXTERNAL :: DGEQP3, DTRSM
! check that d-1 > 0
IF (d > 1) THEN
        ! get QR
        piv=0
        LQ = A
        CALL DGEQP3(d-1,d,LQ,d,piv,plane,work,5*d,err)
        IF(err /= 0) RETURN
        ! do triangular solve to get plane
        ! store in work array
        work(1:d-1) = LQ(1:d-1,d)
        CALL DTRSM('L','U','N','N',d-1,1,REAL(-1.,R8),LQ,d,work,d)
        work(d) = REAL(1.,R8)
        ! undo pivots
        DO i = 1,d
                plane(piv(i)) = work(i)
        END DO
        ! normalize
        plane(1:d) = plane(1:d) / DNRM2(d,plane(1:d),1)
        ! calculate intercept
        plane(d+1) = DDOT(d,plane(1:d),1,pts(:,simp(1)),1)
ELSE
        ! compute plane for d=1
        plane(1) = REAL(1.,R8)
        plane(2) = pts(1,simp(1))
END IF
! done 
RETURN
END SUBROUTINE MakePlane

FUNCTION PtInSimp() RESULT(tf)
! Determine if a point is in the current simplex. The simplex can also  
! be thought of as a cone of vectors s(i) - s(1) where each s(k) is a   
! vertex for the simplex. The d vectors defined thusly span the entire  
! d-dimensional space. If the linear combination which creates p - s(1) 
! is convex (all values are positive), then the cone contains p. Also,  
! the weights for s(2:d+1) for interpolating p will be given by the d   
! linear coefficients.                                                  

! return value
LOGICAL :: tf
! lapack subroutines
EXTERNAL :: DGESV
! init result
tf = .FALSE.
! solve linear system to get weights
LQ = -TRANSPOSE(A)
weights(2:d+1) = proj(:) - pts(:,simp(1))
CALL DGESV(d, 1, LQ, d, piv, weights(2:d+1), d, err)
IF (err /= 0) RETURN
! get last interpolation weight
weights(1) = REAL(1.,R8) - SUM(weights(2:d+1))
! check if weights define a convex combination
IF (ALL(weights >= -eps)) tf = .TRUE.
! done
RETURN
END FUNCTION PtInSimp

END SUBROUTINE delaunayp

END MODULE VTdelaunay


