
MODULE VTdelaunay
!
! Module containing DelaunayP() subroutine for computing a delaunay
! simplex containing a given point and R8 data type for ~64 bit
! arithmetic on most known machines.
!
! Written to FORTRAN 2003 standard. 
! Requires BLAS and LAPACK
!
! Tested w/ gfortran 4.8.5 on x86_64-redhat-linux       
! AUTHOR: Tyler Chang                                   
! Last Update: October, 2017                               
!
USE real_precision

IMPLICIT NONE

PUBLIC
PRIVATE :: GetDiameter

CONTAINS

SUBROUTINE DelaunayP(d, pts, p, work, simp, weights, err, &
        interp_in_opt, interp_out_opt,                    &
        eps_opt, rnorm_opt, budget_opt, extrap_opt)
!
! Compute a simplex in the delaunay triangulation containing the point(s)
! in p. Compute a delaunay simplex "close to" each point in p by solving
! a series of LS problems. Proceed to connect faces of the triangulation
! until a simplex is found containing each point in p. 
!
! INTEGER (IN) d : Dimension of space.                                  
!
! DOUBLE (IN) pts(d,n) : Set of n points in R^d to triangulate.
!       Note: n is assumed based on SIZE(pts).
!
! DOUBLE (IN) p(d,m) : Set of m points in R^d to interpolate using the
!       delaunay triangulation. 
!       Note: m is assumed based on SIZE(p)
!
! DOUBLE (INOUT) work(lwork) : Work array. The work array must have size
!       of AT LEAST 5*d with optimal performance at size m*d. 
!
! INTEGER (OUT) simp(d+1) : The indices of the points which make up the 
!       vertices of the simplex containing p.       
!
! DOUBLE (OUT) weights(d+1) : The weights for the vertices of simp for  
!       interpolating at p.                       
!
! INTEGER (OUT) err(m) : Vector of error flags corresponding to each
!       point in p:
!       SUCCESSFUL INTERPOLATION = 0
!       SUCCESSFUL EXTRAPOLATION = 1
!       OVER EXTRAPOLATION, NO SOLN COMPUTED = 2
!       See other codes in ErrorCodes.txt file.
!       Note: m must match SIZE(p,2)
!
! -------------------------OPTIONAL INPUTS/OUTPUTS-------------------------
!
! OPTIONAL, DOUBLE (IN) interp_in_opt(z,n) : A vector of (z) response values
!       for each of the n triangulation points in pts. Note: z is assumed
!       based on SIZE(interp_in_opt).
!
! OPTIONAL, DOUBLE (out) interp_out_opt(z,m) : A vector of (z) interpolated
!       values for each of the m interpolation points in p. Note: z and m 
!       must match the dimensions of interp_in_opt and p respectively.
!
! OPTIONAL, DOUBLE (IN) eps_opt : A small number defining the tolerance for
!       accepting "nearly" containing simplices. 0 < eps <<< scale of data.
!       Default value: eps = DIAMETER(PTS) * SQRT(EPSILON(KIND=R8)).
! OPTIONAL, DOUBLE (IN) extrap_opt : If extrap_opt = 0 (to working
!       precision), then no extrapolation is performed. Otherwise,
!       extrapolation will be done when p is outside the convex hull but
!       within extrap_opt * diameter(pts) distance of the convex hull.
!       By default, extrap_opt = 0.1 (10% of diameter extrapolation).
! OPTIONAL, DOUBLE (OUT) rnorm_opt(m) : The 2-norm of the residual vector
!       when projecting p onto the convex hull of pts. This value is always
!       calculated when extrapolation is turned on, but only returned when
!       this parameter is present.
! OPTIONAL, INTEGER (IN) budget_opt : The maximum number of simplices that
!       the user is willing to construct when searching for the simplex
!       containing p. In practice, solutions should be found in no more than
!       (d)log(n) iterations, with a typical case of 1~10. The default value
!       for budget is: budget = 1000. If failure to converge in this many
!       iterations, it is most likely the case that eps is too small for the
!       scale of this problem.
!

! inputs
INTEGER, INTENT(IN) :: d
REAL(KIND=R8), INTENT(IN) :: pts(:,:), p(:,:)
! work array
REAL(KIND=R8), INTENT(INOUT) :: work(:)
! outputs
INTEGER, INTENT(OUT) :: simp(:,:)
REAL(KIND=R8), INTENT(OUT) :: weights(:,:)
INTEGER, INTENT(OUT) :: err(:)
! optional arguments
REAL(KIND=R8), INTENT(IN), OPTIONAL:: interp_in_opt(:,:), eps_opt, extrap_opt 
REAL(KIND=R8), INTENT(OUT), OPTIONAL :: interp_out_opt(:,:), rnorm_opt(:)
INTEGER, INTENT(IN), OPTIONAL :: budget_opt
! local vars
REAL(KIND=R8) :: currRad, diam, eps, extrap, minRad, rnorm, side1, side2
INTEGER :: budget, i, j, k, lwork, m, mi, n, tmp
! local arrays (automatic) : O(d^2) extra memory
REAL(KIND=R8) :: A(d,d), b(d), center(d), LQ(d,d), plane(d+1), proj(d)
INTEGER :: piv(d)
! blas functions
REAL(KIND=R8), EXTERNAL :: DDOT
! external subroutine
EXTERNAL :: Project

! get problem dimensions
n = SIZE(pts,2)
m = SIZE(p,2)
! check for input size errors
IF ( d < 1 .OR. n < 1 .OR. m < 1        &
   .OR. SIZE(pts,1) /= d                &
   .OR. SIZE(p,1) /= d                  &
   .OR. SIZE(simp,1) /= d+1             &
   .OR. SIZE(weights,1) /= d+1          &
   .OR. SIZE(simp,2) /= m               &
   .OR. SIZE(weights,2) /= m            &
   .OR. SIZE(err,1) /= m) THEN
        err = 10
        RETURN
ELSE IF (n < d+1) THEN
        err = 11
        RETURN
END IF

! check work array size
lwork = SIZE(work,1)
IF (lwork < 5*d) THEN
        err = 14
        RETURN
END IF

! compute diameter and min-distance between points
CALL GetDiameter(d, n, pts, minRad, diam)

! check for optional inputs arguments
IF (PRESENT(interp_in_opt) .NEQV. &
    PRESENT(interp_out_opt) ) THEN
        err = 15
        RETURN
END IF
IF (PRESENT(interp_in_opt)) THEN
        ! if present, sizes must agree
        IF(SIZE(interp_in_opt,1) /=          &
           SIZE(interp_out_opt,1) .OR.       &
           SIZE(interp_in_opt,2) /= n .OR.   &
           SIZE(interp_out_opt,2) /= m) THEN
                err = 10
                RETURN
        END IF
        interp_out_opt = 0
END IF
IF (PRESENT(eps_opt)) THEN
        eps = eps_opt
        IF (eps <= 0) THEN
                err = 10
                RETURN
        END IF
ELSE
        eps = diam*SQRT(EPSILON(eps))
END IF
IF (PRESENT(budget_opt)) THEN
        budget = budget_opt
        IF (budget < 1) THEN
                err = 10
                RETURN
        END IF
ELSE
        budget = 1000
END IF
IF (PRESENT(extrap_opt)) THEN
        extrap = extrap_opt*diam
        IF (extrap < 0) THEN
                err = 10
        END IF
ELSE
        extrap = 0.1_R8*diam
END IF
IF (PRESENT(rnorm_opt)) THEN
        ! if present, the size of rnorm must match
        IF (SIZE(rnorm_opt,1) /= m) THEN
                err = 10
                RETURN
        END IF
        ! initialize to 0 values
        rnorm_opt = 0.0_R8
END IF
! check for degeneracies in point spacing
IF (minRad < eps) THEN
        err = 13
        RETURN
END IF

! initialize error codes to "TBD" values
err = 40

! outer loop over all points in p(:,:)
outer : DO mi = 1, m

! check if this value was already found
IF (err(mi) == 0) CYCLE outer

! initialize the projection
proj(:) = p(:,mi)
! reset the residual
rnorm = 0.0_R8

! check for previous good simplices
DO k = mi-1, 1, -1
        IF(err(k) <= 1) THEN
                ! use old simplex
                simp(:,mi) = simp(:,k)
                ! rebuild linear system
                DO j=1,d
                        A(j,:) = pts(:,simp(j+1,mi)) - pts(:,simp(1,mi))
                        b(j) = DDOT(d, A(j,:), 1, A(j,:), 1) / 2.0_R8 
                END DO
                ! continue with this value
                EXIT
        END IF
END DO
! none found, build a simplex "near" p(:,mi)
IF (k == 0) THEN
        CALL MakeFirstSimp()
        IF(err(mi) /= 0) CYCLE outer
END IF

! inner loop searching for simplex containing p(mi)
inner : DO k = 1, budget

        ! check if contains p(:,mi)
        IF (PtInSimp()) EXIT inner
        IF (err(mi) /= 0) CYCLE outer

        ! swap out least weighted vertex
        i = MINLOC(weights(1:d+1,mi), DIM=1)
        tmp = simp(i,mi)
        simp(i,mi) = simp(d+1,mi)

        ! if i /= 1, simply drop a row from linear system
        IF(i /= 1) THEN
                A(i-1,:) = A(d,:)
                b(i-1) = b(d)
        ! if i == 1, must reconstruct A & b
        ELSE
                DO j=1,d
                        A(j,:) = pts(:,simp(j+1,mi)) - pts(:,simp(1,mi))
                        b(j) = DDOT(d, A(j,:), 1, A(j,:), 1) / 2.0_R8 
                END DO
        END IF

        ! compute next simplex
        CALL MakeSimplex()
        IF (err(mi) /= 0) CYCLE outer

        ! check for extrapolation condition
        IF (simp(d+1,mi) == 0) THEN
                ! if extrapolation not allowed, do not proceed
                IF (extrap < eps) THEN
                        ! zero values
                        simp(:,mi) = 0
                        weights(:,mi) = 0
                        ! set error flag
                        err(mi) = 2
                        CYCLE outer
                END IF

                ! project p(:,mi) onto convex hull
                CALL Project(d, n, pts(:,:), p(:,mi), proj(:), &
                        rnorm, err(mi))
                IF (err(mi) /= 0) CYCLE outer

                ! check for over-extrapolation
                IF (rnorm > extrap) THEN
                        ! zero values
                        simp(:,mi) = 0
                        weights(:,mi) = 0
                        ! set error flag
                        err(mi) = 2
                        CYCLE outer
                END IF

                ! remake previous simplex
                simp(d+1,mi) = tmp
                A(d,:) = pts(:,tmp) - pts(:,simp(1,mi))
                b(d) = DDOT(d, A(d,:), 1, A(d,:), 1) / 2.0_R8
        END IF

! end of inner loop for finding one point
END DO inner

! check for budget violation
IF (k > budget) THEN 
        ! zero values
        simp(:,mi) = 0
        weights(:,mi) = 0
        ! set error flag
        err(mi) = 20
        CYCLE outer
END IF

! set extrapolation error flag
IF (rnorm > eps) THEN
        err(mi) = 1
END IF

! optional outputs
IF (PRESENT(rnorm_opt)) THEN
        rnorm_opt(mi) = rnorm
END IF

! end of outer loop over all points
END DO outer

! compute interpolation
IF (PRESENT(interp_in_opt)) THEN
        ! loop over all interpolation points
        DO mi = 1, m
                ! check for errors
                IF (err(mi) <= 1) THEN
                        ! compute weighted sum of vertices
                        DO k = 1, d+1
                            interp_out_opt(:,mi) = interp_out_opt(:,mi) &
                                    + interp_in_opt(:,simp(k,mi))       &
                                    * weights(k,mi)
                        END DO
                END IF
        END DO
END IF

! done
RETURN
! end subroutine

! internal subroutines & functions
CONTAINS

SUBROUTINE MakeFirstSimp()
!
! Build first simplex by solving sequence of LS problems ranging from
! 2Xd up to dXd. At each step the center is given by the minimum norm   
! solution to Ac = b. Then the radius is given by r = ||c||_2.          
! A^(m) = [ x2 - x1 | x3 - x1 | ... | xm - x1 ] where x1 - xm-1 are the 
! current points in the face. The xm which results in the minimum r is  
! added to A^(m+1) until an entire simplex (d+1 points) has been
! constructed.
!
! Dummy variables U and VT for storing U and V^T when computing SVD.
! Values are never initialized or used.
REAL(KIND=R8), ALLOCATABLE :: U(:), VT(:)
! blas functions & lapack subroutines
REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
EXTERNAL :: DGELS, DGESVD
! find first point
simp(:,mi) = 0
minRad = HUGE(minRad)
DO i = 1, n
        ! check distance to p
        currRad = DNRM2(d, pts(:,i) - proj(:), 1)
        IF (currRad < minRad) THEN
                minRad = currRad
                simp(1,mi) = i
        END IF
END DO
! check that a point was found
IF (simp(1,mi) == 0) THEN
        err(mi) = 91
        RETURN
END IF
! find second point
minRad = HUGE(minRad)
DO i = 1, n
        ! avoid repeats
        IF (i == simp(1,mi)) CYCLE
        ! check radius
        currRad = DNRM2(d, pts(:,i)-pts(:,simp(1,mi)), 1)
        IF (currRad < minRad) THEN
                minRad = currRad
                simp(2,mi) = i
        END IF
END DO
! check that a point was found
IF (simp(2,mi) == 0) THEN
        err(mi) = 91
        RETURN
END IF
! set up first row of LS system
A(1,:) = pts(:,simp(2,mi)) - pts(:,simp(1,mi))
b(1) = DDOT(d, A(1,:), 1, A(1,:), 1) / 2.0_R8
! loop over remaining d-1 points in first simplex
DO i = 2, d
        ! re-init radius for each iteration
        minRad = HUGE(minRad)
        ! find next point to add
        DO j = 1, n
                ! check point is not already in face
                IF (ANY(simp(:,mi) == j)) CYCLE
                ! add curr point to LS system
                A(i,:) = pts(:,j) - pts(:,simp(1,mi))
                b(i) = DDOT(d, A(i,:), 1, A(i,:), 1) / 2.0_R8
                ! solve LS problem for min norm solution
                LQ(1:i,:) = A(1:i,:)
                center(1:i) = b(1:i)
                ! solve LS problem using QR/LQ factorization
                CALL DGELS('N', i, d, 1, LQ, d, center, d, &
                        work, lwork, err(mi))
                ! check for errors
                IF (err(mi) < 0) THEN
                        err(mi) = 90
                        RETURN
                ! indicates rank-deficiency detected
                ELSE IF (err(mi) > 0) THEN
                        CYCLE
                END IF
                ! calculate radius
                currRad = DNRM2(d, center, 1)
                ! check for new minimum
                IF (currRad < minRad) THEN
                        minRad = currRad
                        simp(i+1,mi) = j
                END IF
        END DO
        ! check that a point was found
        IF (simp(i+1,mi) == 0) THEN
                EXIT
        END IF
        ! add next point to LS system
        A(i,:) = pts(:,simp(i+1,mi)) - pts(:,simp(1,mi))
        b(i) = DDOT(d, A(i,:), 1, A(i,:), 1) / 2.0_R8 
END DO
! If no rank deficiency was detected, double-check
IF (i > d) THEN
        ! compute the SVD
        LQ = A
        CALL DGESVD('N','N',d,d,LQ,d,plane(1:d),U,1,VT,1,work,lwork,err(mi))
        ! check for errors
        IF (err(mi) < 0) THEN
                err(mi) = 90
                RETURN
        ELSE IF (err(mi) > 0) THEN
                err(mi) = 92
                RETURN
        END IF
        ! Check for rank deficiency up to working precision.
        ! If found, re-compute first face using SVD.
        IF (plane(d)/plane(1) < eps) THEN
                CALL MakeFirstSimp_RankDef()
                IF (err(mi) /= 0) RETURN
        END IF
! Otherwise rank deficiency already detected
ELSE
        CALL MakeFirstSimp_RankDef()
        IF (err(mi) /= 0) RETURN
        
END IF
! done
err(mi) = 0
RETURN
END SUBROUTINE MakeFirstSimp

SUBROUTINE MakeFirstSimp_RankDef()
!
! Build first simplex by solving sequence of LS problems ranging from
! 2Xd up to dXd. At each step the center is given by the minimum norm   
! solution to Ac = b. Then the radius is given by r = ||c||_2.          
! A^(m) = [ x2 - x1 | x3 - x1 | ... | xm - x1 ] where x1 - xm-1 are the 
! current points in the face. The xm which results in the minimum r is  
! added to A^(m+1) until an entire simplex (d+1 points) has been
! constructed.
! Use the SVD of A to check for rank defficiency at each step.
!
! blas functions & lapack subroutines
REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
EXTERNAL :: DGELSS
! We already found the first 2 points, their indices are in simp(1)
! and simp(2) and the LS system components A(1,:) and b(1) are
! constructed.
simp(3:d+1,mi) = 0
! Loop over remaining d-1 points in first simplex
DO i = 2, d
        ! re-init radius for each iteration
        minRad = HUGE(minRad)
        ! find next point to add
        DO j = 1, n
                ! check point is not already in face
                IF (ANY(simp(:,mi) == j)) CYCLE
                ! add curr point to LS system
                A(i,:) = pts(:,j) - pts(:,simp(1,mi))
                b(i) = DDOT(d, A(i,:), 1, A(i,:), 1) / 2.0_R8
                ! solve LS problem for min norm solution
                LQ(1:i,:) = A(1:i,:)
                center(1:i) = b(1:i)
                ! Solve using SVD. Store effective rank in tmp.
                CALL DGELSS(i, d, 1, LQ, d, center, d, plane, eps, &
                        tmp, work, lwork, err(mi))
                ! check for errors
                IF (err(mi) < 0) THEN
                        err(mi) = 90
                        RETURN
                ELSE IF (err(mi) > 0) THEN
                        err(mi) = 92
                        RETURN
                END IF
                ! if rank deficient, cycle
                IF (tmp < i) CYCLE
                ! calculate radius
                currRad = DNRM2(d, center, 1)
                ! check for new minimum
                IF (currRad < minRad) THEN
                        minRad = currRad
                        simp(i+1,mi) = j
                END IF
        END DO
        ! check that a point was found
        IF (simp(i+1,mi) == 0) THEN
                err(mi) = 12
                RETURN
        END IF
        ! add next point to LS system
        A(i,:) = pts(:,simp(i+1,mi)) - pts(:,simp(1,mi))
        b(i) = DDOT(d, A(i,:), 1, A(i,:), 1) / 2.0_R8 
END DO
! done
err(mi) = 0
RETURN
END SUBROUTINE MakeFirstSimp_RankDef

SUBROUTINE MakeSimplex()
!
! Given a face of a delaunay simplex, complete the simplex by adding a  
! point on the same side of the face as p. The new point must minimize  
! the radius of the circumcircle. Determine the radius by solving the   
! linear system Ac = b, then taking r = ||c||_2. Note that              
! A = [ x2 - x1 | x3 - x1 | ... | x(d+1) - x1 ] where x1 ... xd are the 
! points in the face and x(d+1) is a potential new point (completing the
! simplex). b_i = dot(A_i,A_i) / -2.                                    
!
! blas functions & lapack subroutines
REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
EXTERNAL :: DGESV
! calculate halfspace
CALL MakePlane()
IF(err(mi) /= 0) RETURN
! calculate side of the viable points
side1 =  DDOT(d,plane(1:d),1,proj(:),1) - plane(d+1)
! normalize magnitude of p
side1 = SIGN(1.0_R8,side1)
! initialize center, radius, and simplex
simp(d+1,mi) = 0
center = 0_R8
minRad = HUGE(minRad)
! loop through all points
DO i = 1, n
        ! check point(i) is inside current ball
        IF (DNRM2(d, pts(:,i) - center(:), 1) > minRad) CYCLE
        ! check point(i) on viable side of halfspace
        side2 = DDOT(d,plane(1:d),1,pts(:,i),1) - plane(d+1)
        IF (side1 * side2 < eps .OR. ANY(simp(:,mi) == i)) CYCLE
        ! add point i to linear system
        A(d,:) = pts(:,i) - pts(:,simp(1,mi))
        b(d) = DDOT(d, A(d,:), 1, A(d,:), 1) / 2.0_R8
        LQ = A
        center = b
        ! solve linear system to get center
        CALL DGESV(d, 1, LQ, d, piv, center, d, err(mi))
        IF (err(mi) < 0) THEN
                err(mi) = 90
                RETURN
        ELSE IF (err(mi) > 0) THEN
                err(mi) = 21
                RETURN
        END IF
        ! update radius, center, and simplex
        currRad = DNRM2(d, center, 1)
        minRad = currRad
        center = center + pts(:,simp(1,mi))
        simp(d+1,mi) = i
END DO
! reset error flag
err(mi) = 0
! check for extrapolation condition
IF(simp(d+1,mi) == 0) RETURN
! otherwise add point to linear system
A(d,:) = pts(:,simp(d+1,mi)) - pts(:,simp(1,mi))
b(d) = DDOT(d, A(d,:), 1, A(d,:), 1) / 2.0_R8
! done
RETURN
END SUBROUTINE MakeSimplex

SUBROUTINE MakePlane()
!
! Construct a plane intersecting the current d points. The coefficients 
! are determined by the normal vector: v and an intercept, stored in    
! plane(1:d) and plane(d+1) respectively. To find the normal vector, the
! current d points are transformed into d-1 vectors in the plane.       
! Together, these d-1 d-vectors have nullity rank-1. The normal vector  
! is any nontrivial vector in the nullspace.                            
!
! blas functions & lapack subroutines
REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
EXTERNAL :: DGEQP3, DTRSM
! check that d-1 > 0
IF (d > 1) THEN
        ! get QR
        piv=0
        LQ = A
        CALL DGEQP3(d-1,d,LQ,d,piv,plane,work,lwork,err(mi))
        IF(err(mi) < 0) THEN
                err(mi) = 90
                RETURN
        ELSE IF (err(mi) > 0) THEN
                err(mi) = 21
                RETURN
        END IF
        ! do triangular solve to get plane
        ! store in work array
        work(1:d-1) = LQ(1:d-1,d)
        CALL DTRSM('L','U','N','N',d-1,1,-1.0_R8,LQ,d,work,d)
        work(d) = 1.0_R8
        ! undo pivots
        DO i = 1,d
                plane(piv(i)) = work(i)
        END DO
        ! normalize
        plane(1:d) = plane(1:d) / DNRM2(d,plane(1:d),1)
        ! calculate intercept
        plane(d+1) = DDOT(d,plane(1:d),1,pts(:,simp(1,mi)),1)
ELSE
        ! compute plane for d=1
        plane(1) = 1.0_R8
        plane(2) = pts(1,simp(1,mi))
END IF
! done 
RETURN
END SUBROUTINE MakePlane

FUNCTION PtInSimp() RESULT(tf)
!
! Determine if point is in the current simplex. The simplex can also  
! be thought of as a cone of vectors s(i) - s(1) where each s(k) is a   
! vertex for the simplex. The d vectors defined thusly span the entire  
! d-dimensional space. If the linear combination which creates p - s(1) 
! is convex (all values are positive), then the cone contains p. Also,  
! the weights for s(2:d+1) for interpolating p will be given by the d   
! linear coefficients.                                                  
! While doing so, also check if any future points are in the simplex.
! If so, resolve their solutions.
!
! return value
LOGICAL :: tf
! local data
INTEGER :: nrhs
! lapack subroutines
EXTERNAL :: DGESV
! calculate number of points we can check in this iteration
nrhs = MIN(lwork/d, m+1-mi)
! init result
tf = .FALSE.
! build matrix for LU factorization
LQ = TRANSPOSE(A)
! fill (nrhs) RHS equations (one for each point)
work(1:d) = proj(:) - pts(:,simp(1,mi))
DO i = mi+1, mi+nrhs-1
        work((i-mi)*d+1:(i-mi+1)*d) = p(:,i) &
                - pts(:,simp(1,mi))
END DO
! solve linear system to get ALL weights
CALL DGESV(d, nrhs, LQ, d, piv, work(1:nrhs*d), d, err(mi))
IF (err(mi) < 0) THEN
        err(mi) = 90
        RETURN
ELSE IF (err(mi) > 0) THEN
        err(mi) = 21
        RETURN
END IF
! get interpolation weights for p(:,mi)
weights(2:d+1,mi) = work(1:d)
weights(1,mi) = 1.0_R8 - SUM(weights(2:d+1,mi))
! check if weights define a convex combination
IF (ALL(weights(:,mi) >= -eps)) tf = .TRUE.
! get the rest of the interpolation weights
DO i = 1, nrhs-1
        ! check if a solution was already found
        IF (err(i+mi) == 40) THEN
                ! copy weights from work array
                weights(2:d+1,i+mi) = work(d*i+1:d*(i+1))
                weights(1,i+mi) = 1.0_R8 - SUM(weights(2:d+1,i+mi))
                ! copy simplex
                simp(:, i+mi) = simp(:,mi)
                ! check if weights define a convex combination
                IF (ALL(weights(:,mi+i) >= -eps)) err(i+mi) = 0
        END IF
END DO
! done
RETURN
END FUNCTION PtInSimp

END SUBROUTINE DelaunayP

! --- Private auxillary subroutines ---

SUBROUTINE GetDiameter(d, n, pts, minD, maxD)
!
! Compute the diameter of a point set by brute force.
! O(n^2) complexity.
! INTEGER (IN) d : Dimension of space.                                  
! INTEGER (IN) n : Number of points in set.                             
! DOUBLE (IN) pts(d,n) : Set of n points in d-space.                    
! DOUBLE (OUT) minD : Minimum distance between 2 points
! DOUBLE (OUT) maxD : Maximum distance between 2 points
! INTEGER(OUT) error : Error flag
!
! inputs
INTEGER, INTENT(IN) :: d, n
REAL(KIND=R8), INTENT(IN) :: pts(d,n)
! outputs
REAL(KIND=R8), INTENT(OUT) :: minD, maxD
! local vars
REAL(KIND=R8) :: distance
INTEGER :: i, j
! blas function
REAL(KIND=R8), EXTERNAL :: DNRM2
! initialize values
minD = HUGE(minD)
maxD = 0
! compute min and max distance in dataset
DO i = 1, n-1
        DO j = i+1, n
                distance = DNRM2(d, pts(1:d,i) - pts(1:d,j), 1)
                IF (distance > maxD) THEN
                        maxD = distance
                END IF
                IF (distance < minD) THEN
                        minD = distance
                END IF
        END DO
END DO
! done
RETURN
END SUBROUTINE GetDiameter

END MODULE VTdelaunay
