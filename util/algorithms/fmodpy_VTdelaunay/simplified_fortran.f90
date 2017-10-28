

MODULE VTDELAUNAY
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
USE REAL_PRECISION


PUBLIC
PRIVATE : : GETDIAMETER


SUBROUTINE DELAUNAYP ( D , PTS , P , WORK , SIMP , WEIGHTS , ERR , INTERP_IN_OPT , INTERP_OUT_OPT , EPS_OPT , RNORM_OPT , BUDGET_OPT , EXTRAP_OPT )
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
INTEGER , INTENT ( IN ) : : D
REAL ( KIND = R8 ) , INTENT ( IN ) : : PTS ( : , : ) , P ( : , : )
! work array
REAL ( KIND = R8 ) , INTENT ( INOUT ) : : WORK ( : )
! outputs
INTEGER , INTENT ( OUT ) : : SIMP ( : , : )
REAL ( KIND = R8 ) , INTENT ( OUT ) : : WEIGHTS ( : , : )
INTEGER , INTENT ( OUT ) : : ERR ( : )
! optional arguments
REAL ( KIND = R8 ) , INTENT ( IN ) , OPTIONAL : : INTERP_IN_OPT ( : , : ) , EPS_OPT , EXTRAP_OPT
REAL ( KIND = R8 ) , INTENT ( OUT ) , OPTIONAL : : INTERP_OUT_OPT ( : , : ) , RNORM_OPT ( : )
INTEGER , INTENT ( IN ) , OPTIONAL : : BUDGET_OPT
! local vars
REAL ( KIND = R8 ) : : CURRRAD , DIAM , EPS , EXTRAP , MINRAD , RNORM , SIDE1 , SIDE2
INTEGER : : BUDGET , I , J , K , LWORK , M , MI , N , TMP
! local arrays (automatic) : O(d^2) extra memory
REAL ( KIND = R8 ) : : A ( D , D ) , B ( D ) , CENTER ( D ) , LQ ( D , D ) , PLANE ( D + 1 ) , PROJ ( D )
INTEGER : : PIV ( D )
! blas functions
REAL ( KIND = R8 ) , EXTERNAL : : DDOT
! external subroutine
EXTERNAL : : PROJECT

! get problem dimensions
! check for input size errors

! check work array size

! compute diameter and min-distance between points

! check for optional inputs arguments
! if present, sizes must agree
! if present, the size of rnorm must match
! initialize to 0 values
! check for degeneracies in point spacing

! initialize error codes to "TBD" values

! outer loop over all points in p(:,:)

! check if this value was already found

! initialize the projection
! reset the residual

! check for previous good simplices
! use old simplex
! rebuild linear system
! continue with this value
! none found, build a simplex "near" p(:,mi)

! inner loop searching for simplex containing p(mi)

! check if contains p(:,mi)

! swap out least weighted vertex

! if i /= 1, simply drop a row from linear system
! if i == 1, must reconstruct A & b

! compute next simplex

! check for extrapolation condition
! if extrapolation not allowed, do not proceed
! zero values
! set error flag

! project p(:,mi) onto convex hull

! check for over-extrapolation
! zero values
! set error flag

! remake previous simplex

! end of inner loop for finding one point

! check for budget violation
! zero values
! set error flag

! set extrapolation error flag

! optional outputs

! end of outer loop over all points

! compute interpolation
! loop over all interpolation points
! check for errors
! compute weighted sum of vertices

! done
! end subroutine

! internal subroutines & functions

SUBROUTINE MAKEFIRSTSIMP ( )
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
REAL ( KIND = R8 ) , ALLOCATABLE : : U ( : ) , VT ( : )
! blas functions & lapack subroutines
REAL ( KIND = R8 ) , EXTERNAL : : DDOT , DNRM2
EXTERNAL : : DGELS , DGESVD
! find first point
! check distance to p
! check that a point was found
! find second point
! avoid repeats
! check radius
! check that a point was found
! set up first row of LS system
! loop over remaining d-1 points in first simplex
! re-init radius for each iteration
! find next point to add
! check point is not already in face
! add curr point to LS system
! solve LS problem for min norm solution
! solve LS problem using QR/LQ factorization
! check for errors
! indicates rank-deficiency detected
! calculate radius
! check for new minimum
! check that a point was found
! add next point to LS system
! If no rank deficiency was detected, double-check
! compute the SVD
! check for errors
! Check for rank deficiency up to working precision.
! If found, re-compute first face using SVD.
! Otherwise rank deficiency already detected

! done
END SUBROUTINE MAKEFIRSTSIMP

SUBROUTINE MAKEFIRSTSIMP_RANKDEF ( )
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
REAL ( KIND = R8 ) , EXTERNAL : : DDOT , DNRM2
EXTERNAL : : DGELSS
! We already found the first 2 points, their indices are in simp(1)
! and simp(2) and the LS system components A(1,:) and b(1) are
! constructed.
! Loop over remaining d-1 points in first simplex
! re-init radius for each iteration
! find next point to add
! check point is not already in face
! add curr point to LS system
! solve LS problem for min norm solution
! Solve using SVD. Store effective rank in tmp.
! check for errors
! if rank deficient, cycle
! calculate radius
! check for new minimum
! check that a point was found
! add next point to LS system
! done
END SUBROUTINE MAKEFIRSTSIMP_RANKDEF

SUBROUTINE MAKESIMPLEX ( )
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
REAL ( KIND = R8 ) , EXTERNAL : : DDOT , DNRM2
EXTERNAL : : DGESV
! calculate halfspace
! calculate side of the viable points
! normalize magnitude of p
! initialize center, radius, and simplex
! loop through all points
! check point(i) is inside current ball
! check point(i) on viable side of halfspace
! add point i to linear system
! solve linear system to get center
! update radius, center, and simplex
! reset error flag
! check for extrapolation condition
! otherwise add point to linear system
! done
END SUBROUTINE MAKESIMPLEX

SUBROUTINE MAKEPLANE ( )
!
! Construct a plane intersecting the current d points. The coefficients
! are determined by the normal vector: v and an intercept, stored in
! plane(1:d) and plane(d+1) respectively. To find the normal vector, the
! current d points are transformed into d-1 vectors in the plane.
! Together, these d-1 d-vectors have nullity rank-1. The normal vector
! is any nontrivial vector in the nullspace.
!
! blas functions & lapack subroutines
REAL ( KIND = R8 ) , EXTERNAL : : DDOT , DNRM2
EXTERNAL : : DGEQP3 , DTRSM
! check that d-1 > 0
! get QR
! do triangular solve to get plane
! store in work array
! undo pivots
! normalize
! calculate intercept
! compute plane for d=1
! done
END SUBROUTINE MAKEPLANE

FUNCTION PTINSIMP ( ) RESULT ( TF )
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
LOGICAL : : TF
! local data
INTEGER : : NRHS
! lapack subroutines
EXTERNAL : : DGESV
! calculate number of points we can check in this iteration
! init result
! build matrix for LU factorization
! fill (nrhs) RHS equations (one for each point)
! solve linear system to get ALL weights
! get interpolation weights for p(:,mi)
! check if weights define a convex combination
! get the rest of the interpolation weights
! check if a solution was already found
! copy weights from work array
! copy simplex
! check if weights define a convex combination
! done
END FUNCTION PTINSIMP

END SUBROUTINE DELAUNAYP

! --- Private auxillary subroutines ---

SUBROUTINE GETDIAMETER ( D , N , PTS , MIND , MAXD )
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
INTEGER , INTENT ( IN ) : : D , N
REAL ( KIND = R8 ) , INTENT ( IN ) : : PTS ( D , N )
! outputs
REAL ( KIND = R8 ) , INTENT ( OUT ) : : MIND , MAXD
! local vars
REAL ( KIND = R8 ) : : DISTANCE
INTEGER : : I , J
! blas function
REAL ( KIND = R8 ) , EXTERNAL : : DNRM2
! initialize values
! compute min and max distance in dataset
! done
END SUBROUTINE GETDIAMETER

END MODULE VTDELAUNAY

