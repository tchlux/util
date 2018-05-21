!
! This module (VTDELAUNAY_MOD) contains the DELAUNAYINTERP subroutine for
! computing the Delaunay simplices containing interpolation points Q in R^D
! given data points PTS. Other auxiliary routines are included for scaling
! points to the unit ball (SCALETOBALL) and projecting a point in R^D onto the
! convex hull of a point set PTS (PROJECT). Also, the REAL_PRECISION type R8 is
! included for ~64-bit computation on many architectures.
!
! The following BLAS subroutines and functions are directly referenced. All
! directly and indirectly called BLAS soubroutines are contained in the file
! blas.f:
!    - DDOT
!    - DGEMV
!    - DNRM2
!    - DTRSM
!
! The following LAPACK subroutines are directly referenced. All directly and
! indirectly called LAPACK subroutines are contained in the file lapack.f:
!    - DGELS
!    - DGELSS
!    - DGEQP3
!    - DGESV
!    - DGESVD
!
! The following SLATEC subroutine is directly referenced. DWNNLS and all its
! SLATEC dependencies are contained in the file slatec.f. Note, these files have
! been slightly edited, with all print statements and references to stderr being
! commented to avoid legacy issues. For a reference to DWNNLS, see ACM TOMS
! Algorithm 587 (Hanson and Haskell):
!    - DWNNLS
!
! The following code was written for the Fortran 2003 standard.
!
! File Name: VTDelaunay_v16.f90
! Primary Author: Tyler H. Chang
! Last Update: May, 2018
!
MODULE VTDELAUNAY_MOD
USE REAL_PRECISION
IMPLICIT NONE

PUBLIC

CONTAINS
SUBROUTINE DELAUNAYINTERP( D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
 INTERP_IN_OPT, INTERP_OUT_OPT, EPS_OPT, RNORM_OPT, IBUDGET_OPT, &
 EXTRAP_OPT, SEED_OPT )
! This is a serial implementation of an algorithm for efficiently performing
! interpolation in R^D via the Delaunay triangulation. The algorithm is fully
! described and analyzed in:
!
! T. H. Chang, L. T. Watson,  T. C.H. Lux, B. Li, L. Xu, A. R. Butt, K. W.
! Cameron, and Y. Hong. 2018. A polynomial time algorithm for multivariate
! interpolation in arbitrary dimension via the Delaunay triangulation. In
! Proceedings of the ACMSE 2018 Conference (ACMSE '18). ACM, New York, NY, USA.
! Article 12, 8 pages.
!
!
! On input:
!
! D is the dimension of the space for PTS and Q.
!
! N is the number of data points in PTS.
!
! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
! coordinates of a single data point in R^D.
!
! M is the number of interpolation points in Q.
!
! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
! coordinates of a single interpolation point in R^D.
!
!
! On output:
!
! PTS and Q have been rescaled and shifted. All the data points in PTS are now
! contained in the unit hyperball in R^D, and the points in Q have been shifted
! and scaled accordingly in relation to PTS.
!
! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns in
! PTS) for the D+1 vertices of the Delaunay simplex containing each
! interpolation point in Q.
!
! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
! point in Q as a convex combination of the D+1 corresponding vertices in SIMPS.
!
! IERR(1:M) contains integer valued error flags associated with the computation
! of each of the M interpolation points in Q. The error codes are listed below:
!
! 00 : Succesful interpolation.
! 01 : Succesful extrapolation (up to the allowed extrapolation distance).
! 02 : This point was outside the allowed extrapolation distance, the
!      corresponding entries in SIMPS and WEIGHTS contain zero values.
!
! 10 : The given dimension D was not valid.
! 11 : Too few data points to construct a triangulation (I.e., N < D+1).
! 12 : No interpolation points given (I.e., M < 1).
! 13 : The first dimension of PTS does not agree with the dimension D.
! 14 : The second dimension of PTS does not agree with the number of points N.
! 15 : The first dimension of Q does not agree with the dimension D.
! 16 : The second dimension of Q does not agree with the number of interpolation
!      points M.
! 17 : The first dimension of the output array SIMPS does not match the number
!      of vertices needed for a D-simplex (D+1).
! 18 : The second dimension of the output array SIMPS does not match the number
!      of interpolation points M.
! 19 : The first dimension of the output array WEIGHTS does not match the number
!      of vertices for a a D-simplex (D+1).
! 20 : The second dimension of the output array WEIGHTS does not match the
!      number of interpolation points M.
! 21 : The size of the error array IERR does not match the number of
!      interpolation points M.
! 22 : INTERP_IN_OPT cannot be present without INTERP_OUT_OPT or vice versa.
! 23 : The first dimension of INTERP_IN_OPT does not match the first dimension
!      of INTERP_OUT_OPT.
! 24 : The second dimension of INTERP_IN_OPT does not match the number of
!      data points PTS.
! 25 : The second dimension of INTERP_OUT_OPT does not match the number of
!      interpolation points M.
! 26 : The budget supplied in IBUDGET_OPT does not contain a valid positive
!      integer.
! 27 : The extrapolation distance supplied in EXTRAP_OPT cannot be negative.
! 28 : The size of the RNORM_OPT output array does not match the number of
!      interpolation points M.
! 29 : The seed simplex supplied in SEED_OPT is not a valid Delaunay simplex in
!      R^D.
!
! 30 : Two or more points in the data set PTS are too close together with
!      respect to the working precision (EPS_OPT).
! 31 : All the data points in PTS lie in some lower dimensional linear manifold
!      (up to the working precision), and no valid triangulation exists.
!
! 40 : An error caused VTDelaunay to terminate before this value could be
!      computed. Note: The corresponding entries in SIMPS and WEIGHTS may
!      contain garbage values.
!
! 50 : A memory allocation error occurred while allocating the work array.
!
! 60 : The budget was exceeded before the algorithm converged on this value. If
!      the dimension is high, try increasing IBUDGET_OPT. This error can also be
!      caused by a working precision EPS_OPT that is too fine for the
!      conditioning of the problem.
! 61 : A value that was judged appropriate later caused LAPACK to encounter a
!      singularity. Try increasing the value of EPS_OPT.
!
! 70 : The SLATEC subroutine DWNNLS failed to converge during the projection of
!      an extrapolation point onto the convex hull.
! 71 : The SLATEC subroutine DWNNLS has reported a usage error. Possible causes:
!      An automatic array failed to allocate when additional memory was
!      allocated for the projection.
!
! 80 : The LAPACK subroutine DGELS has reported an illegal value.
! 81 : The LAPACK subroutine DGESVD has reported an illegal value.
! 82 : The LAPACK subroutine DGELSS has reported an illegal value.
! 83 : The LAPACK subroutine DGESV has reported an illegal value.
! 84 : The LAPACK subroutine DGEQP3 has reported an illegal value.
! 85 : The LAPACK subroutine DGESVD failed to converge while performing a
!      singular value decomposition.
! 86 : The LAPACK subroutine DGELSS failed to converge while performing a
!      singular value decomposition.
!
!
! Optional arguments:
!
! INTERP_IN_OPT(1:IR,1:N) contains real valued response vectors for each of the
! data points in PTS on input. The first dimension of INTERP_IN_OPT is inferred
! to be the dimension of this response vector (and is inferred), and the second
! dimension must match N. If present, the response values will be computed for
! each interpolation point in Q, and stored in INTERP_OUT_OPT. Therefore, if
! INTERP_IN_OPT is present, then INTERP_OUT_OPT must also be present. If either
! is omitted, the response values will not be computed, and only the containing
! simplices and convex weights are returned.
!
! INTERP_OUT_OPT(1:IR,1:M) contains real valued response vectors for each
! interpolation point in Q on output. The first dimension of INTERP_OUT_OPT must
! match the first dimension of INTERP_IN_OPT, and the second dimension must
! match M. If present, the response values are computed as convex combinations
! of the response values at each data point (as supplied in INTERP_IN_OPT).
! Therefore, if INTERP_OUT_OPT is present, then INTERP_IN_OPT must also be
! present. If either is omitted, the response values will not be computed, and
! only the simplices and ocnvex weights are returned.
!
! EPS_OPT contains the working precision for the problem on input. By default,
! EPS_OPT is assigned \sqrt{\mu} where \mu denotes the unit roundoff for the
! machine. In general, any values that differ by less than EPS_OPT are judged as
! equal, and any weights that are greater than -EPS_OPT are judged as positive.
! EPS_OPT cannot take a value less than the default value of \sqrt{\mu}. If any
! value less than \sqrt{\mu} is supplied, the default value will be used instead
! automatically.
!
! EXTRAP_OPT contains the maximum extrapolation distance (relative to the
! diameter of PTS) on input. Any interpolation point in Q that is more than
! EXTRAP_OPT * DIAMETER(PTS) units outside the convex hull of PTS will not be
! computed and an error code of 2 will be returned. Note, that computing the
! projection can be expensive. Setting EXTRAP_OPT=0 will cause all extrapolation
! points to be ignored without ever computing a projection. By default,
! EXTRAP_OPT=0.1. (Extrapolate by up to 10% of the diameter of PTS).
!
! RNORM_OPT(1:M) contains the residual distances from any projection
! computations on output. If not present, these distances are still computed for
! each extrapolation point, but are never returned.
!
! IBUDGET_OPT contains the integer valued budget for performing flips while
! iterating toward the simplex containing each interpolation point in Q. This
! prevents VTDelaunay from falling into an infinite loop when an inappropriate
! value of EPS_OPT is given with respect to the problem conditioning. By
! default, IBUDGET_OPT=50000. However, for extremely high-dimensional problems
! and pathological inputs, the default value may be insufficient.
!
! SEED_OPT(1:D+1) contains the integer indices (i.e., columns of PTS) for an
! initial ``seed'' Delaunay simplex. If present, SEED_OPT will be used as the
! initial simplex, significantly reducing computational cost (especially for
! small values of M). However, if SEED_OPT is not supplied with a valid Delaunay
! simplex, it is an error. If omitted, the initial simplex will be generated.

! Input arguments.
INTEGER, INTENT(IN) :: D, N, M
REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:), Q(:,:) ! Rescaled on output.
! Output arguments.
INTEGER, INTENT(OUT) :: SIMPS(:,:)
REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
INTEGER, INTENT(OUT) :: IERR(:)

! Optional arguments.
REAL(KIND=R8), INTENT(IN), OPTIONAL:: INTERP_IN_OPT(:,:), EPS_OPT, EXTRAP_OPT 
REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT_OPT(:,:), RNORM_OPT(:)
INTEGER, INTENT(IN), OPTIONAL :: IBUDGET_OPT, SEED_OPT(:)
! Local copy of optional input arguments.
REAL(KIND=R8) :: EPS, EXTRAP
INTEGER :: IBUDGET, SEED(D+1)

! Local variables.
INTEGER :: I, J, K ! Loop iteration variables.
INTEGER :: ITMP ! Temporary variable for swapping indices.
INTEGER :: LWORK ! Size of the work array.
INTEGER :: MI ! Index of current interpolation point.

REAL(KIND=R8) :: CURRRAD ! Radius of the current circumsphere.
REAL(KIND=R8) :: MINRAD ! Minimum radius observed. 
REAL(KIND=R8) :: RNORM ! Euclidean norm of the projection residual.
REAL(KIND=R8) :: SIDE1, SIDE2 ! Signs (+/-1) denoting sides of a facet.

! Local arrays, requiring O(d^2) additional memory.
INTEGER :: IPIV(D) ! Pivot indices.

REAL(KIND=R8) :: A(D,D) ! The LHS of a linear system.
REAL(KIND=R8) :: B(D) ! The RHS of a linear system. 
REAL(KIND=R8) :: CENTER(D) ! The circumcenter of a simplex. 
REAL(KIND=R8) :: LQ(D,D) ! The LU or QR factorization of A.
REAL(KIND=R8) :: PLANE(D+1) ! The hyperplane containing a facet.
REAL(KIND=R8) :: PROJ(D) ! The projection of the current iterate.

! Work arrays.
REAL(KIND=R8), ALLOCATABLE :: WORK(:) ! Allocated with size LWORK.
REAL(KIND=R8) :: U(D,D), VT(D,D) ! Dummy arrays that will not be initialized.

! External functions and subroutines.
REAL(KIND=R8), EXTERNAL :: DDOT ! Inner product (BLAS).
REAL(KIND=R8), EXTERNAL :: DNRM2 ! Euclidean norm (BLAS).
EXTERNAL :: DGELS ! Least squares solve full-rank system (LAPACK).
EXTERNAL :: DGELSS ! Least squares solve rank-deficient system (LAPACK).
EXTERNAL :: DGEQP3 ! Perform a QR factorization (LAPACK). 
EXTERNAL :: DGESV ! Linear solve for exactly determined system (LAPACK).
EXTERNAL :: DGESVD ! Perform a singular value decomposition (LAPACK). 
EXTERNAL :: DTRSM ! Perform a triangular solve (BLAS).

! Check for input size and dimension errors.
IF (D < 1) THEN ! The dimension must satisfy: D > 0.
   IERR(:) = 10
   RETURN
END IF
IF (N < D+1) THEN ! Must have at least D+1 data points.
   IERR(:) = 11
   RETURN
END IF
IF (M < 1) THEN ! Must have at least one interpolation point.
   IERR(:) = 12
   RETURN
END IF
IF (SIZE(PTS,1) .NE. D) THEN ! Dimension of PTS array should match.
   IERR(:) = 13
   RETURN
END IF
IF (SIZE(PTS,2) .NE. N) THEN ! Number of data points should match.
   IERR(:) = 14
   RETURN
END IF
IF (SIZE(Q,1) .NE. D) THEN ! Dimension of Q should match.
   IERR(:) = 15
   RETURN
END IF
IF (SIZE(Q,2) .NE. M) THEN ! Number of interpolation points should match.
   IERR(:) = 16
   RETURN
END IF
IF (SIZE(SIMPS,1) .NE. D+1) THEN ! Need space for D+1 vertices per simplex.
   IERR(:) = 17
   RETURN
END IF
IF (SIZE(SIMPS,2) .NE. M) THEN  ! There will be M output simplices.
   IERR(:) = 18
   RETURN
END IF
IF (SIZE(WEIGHTS,1) .NE. D+1) THEN ! There will be D+1 weights per simplex.
   IERR(:) = 19
   RETURN
END IF
IF (SIZE(WEIGHTS,2) .NE. M) THEN ! One set of weights per simplex.
   IERR(:) = 20
   RETURN
END IF
IF (SIZE(IERR,1) .NE. M) THEN ! An error flag for each interpolation point.
   IERR(:) = 21
   RETURN
END IF

! Check for optional inputs arguments.
IF (PRESENT(INTERP_IN_OPT) .NEQV. PRESENT(INTERP_OUT_OPT)) THEN
   IERR(:) = 22
   RETURN
END IF
IF (PRESENT(INTERP_IN_OPT)) THEN ! Sizes must agree.
   IF (SIZE(INTERP_IN_OPT,1) .NE. SIZE(INTERP_OUT_OPT,1)) THEN
      IERR(:) = 23 
      RETURN
   END IF
   IF(SIZE(INTERP_IN_OPT,2) .NE. N) THEN
      IERR(:) = 24
      RETURN
   END IF
   IF (SIZE(INTERP_OUT_OPT,2) .NE. M) THEN
      IERR(:) = 25
      RETURN
   END IF
   INTERP_OUT_OPT(:,:) = 0.0_R8 ! Initialize output to zeros.
END IF
EPS = SQRT(EPSILON(EPS)) ! Get the machine unit roundoff constant.
IF (PRESENT(EPS_OPT)) THEN
   IF (EPS < EPS_OPT) THEN ! If the given precision is too small, ignore it.
      EPS = EPS_OPT
   END IF
END IF
IF (PRESENT(IBUDGET_OPT)) THEN
   IBUDGET = IBUDGET_OPT ! Use the given budget if present.
   IF (IBUDGET < 1) THEN
      IERR(:) = 26
      RETURN
   END IF
ELSE
   IBUDGET = 50000 ! Default value for budget.
END IF
IF (PRESENT(EXTRAP_OPT)) THEN
   EXTRAP = EXTRAP_OPT 
   IF (EXTRAP < 0) THEN ! Check that the extrapolation distance is legal.
      IERR(:) = 27
   END IF
ELSE
   EXTRAP = 0.1_R8 ! Default extrapolation distance (for normalized points).
END IF
IF (PRESENT(RNORM_OPT)) THEN
   IF ( SIZE(RNORM_OPT,1) .NE. M ) THEN ! The length of the array must match.
      IERR(:) = 28
      RETURN
   END IF
   RNORM_OPT(:) = 0.0_R8 ! Initialize output to zeros.
END IF

! Scale and center the data points in the unit ball, and transform the
! interpolation points similarly.
CALL SCALETOBALL(D, N, PTS, M, Q, MINRAD)
IF (MINRAD < EPS) THEN ! Check for degeneracies in points spacing.
   IERR(:) = 30
   RETURN
END IF

! Now check for the last optional input: a seed simplex.
IF (PRESENT(SEED_OPT)) THEN
   IF (SIZE(SEED_OPT,1) .NE. D+1) THEN ! Seed must have legal size.
      IERR(:) = 29
      RETURN
   END IF
   IF (ANY(SEED_OPT < 1) .OR. ANY(SEED_OPT > N)) THEN ! Must have legal values.
      IERR(:) = 29
      RETURN
   END IF
   SEED = SEED_OPT ! Get the seed.
   DO J=1,D ! Build linear system to check if the seed is empty.
      LQ(J,:) = PTS(:,SEED(J+1)) - PTS(:,SEED(1))
      CENTER(J) = DDOT(D, LQ(J,:), 1, LQ(J,:), 1) / 2.0_R8 
   END DO
   CALL DGESV(D, 1, LQ, D, IPIV, CENTER, D, I) ! Compute the (shifted) center.
   IF(I .NE. 0) THEN ! If an error occurs, the seed was not a legal simplex.
      IERR(:) = 29 
      RETURN
   END IF
   MINRAD = DNRM2(D, CENTER, 1) ! Compute the radius of the circumsphere.
   CENTER(:) = CENTER(:) + PTS(:,SEED(1)) ! Compute the circumcenter.
   ! Check that the circumball is empty (Delaunay property).
   DO J=1,N 
      IF (DNRM2(D, PTS(:,J) - CENTER(:), 1) < MINRAD - EPS) THEN
         IERR(:) = 29 ! A point was detected inside the circumball.
         RETURN
      END IF
   END DO
ELSE
   SEED(:) = 0 ! If no seed was given, initialize SEED to zeros.
END IF

! Query DGESVD for optimal work array size (LWORK).
LWORK = -1
CALL DGESVD('N','N',D,D,LQ,D,PLANE(1:D),U,1,VT,1,B,LWORK,IERR(1))
LWORK = MAX(M*D, 5*D, INT(B(1))) ! Compute the optimal work array size.
ALLOCATE(WORK(LWORK), STAT=I) ! Allocate WORK to size LWORK.
IF (I .NE. 0) THEN ! Check for memory allocation errors.
   IERR(:) = 50
   RETURN
END IF

! Initialize all error codes to "TBD" values.
IERR(:) = 40

! Outer loop over all interpolation points (in Q).
OUTER : DO MI = 1, M

   ! Check if this interpolation point was already found.
   IF (IERR(MI) == 0) CYCLE OUTER

   ! Initialize the projection and reset the residual.
   PROJ(:) = Q(:,MI)
   RNORM = 0.0_R8

   ! If there is no useable seed, then make a new simplex.
   IF(SEED(1) == 0) THEN
      CALL MAKEFIRSTSIMP()
      IF(IERR(MI) .NE. 0) CYCLE OUTER
   ! Otherwise, use the seed.
   ELSE
      ! Copy the seed to the current simplex.
      SIMPS(:,MI) = SEED(:)
      ! Rebuild the linear system.
      DO J=1,D
         A(J,:) = PTS(:,SIMPS(J+1,MI)) - PTS(:,SIMPS(1,MI))
         B(J) = DDOT(D, A(J,:), 1, A(J,:), 1) / 2.0_R8 
      END DO
   END IF

   ! Inner loop searching for a simplex containing the point Q(:,MI).
   INNER : DO K = 1, IBUDGET

      ! Check if the current simplex contains Q(:,MI).
      IF (PTINSIMP()) THEN
         SEED(:) = SIMPS(:,MI) ! Set the next seed to this simplex.
         EXIT INNER
      END IF
      IF (IERR(MI) .NE. 0) CYCLE OUTER ! Check for an error flag.

      ! Save each good simplex as the next seed.
      SEED(:) = SIMPS(:,MI)

      ! Swap out the least weighted vertex.
      I = MINLOC(WEIGHTS(1:D+1,MI), DIM=1)
      ITMP = SIMPS(I,MI)
      SIMPS(I,MI) = SIMPS(D+1,MI)

      ! If the least weighted vertex (I) is not the first vertex, then just
      ! drop row I from the linear system.
      IF(I .NE. 1) THEN
         A(I-1,:) = A(D,:)
         B(I-1) = B(D)
      ! However, if I == 1, then both A and B must be reconstructed.
      ELSE
         DO J=1,D
            A(J,:) = PTS(:,SIMPS(J+1,MI)) - PTS(:,SIMPS(1,MI))
            B(J) = DDOT(D, A(J,:), 1, A(J,:), 1) / 2.0_R8 
         END DO
      END IF

      ! Compute the next simplex (do one flip).
      CALL MAKESIMPLEX()
      IF (IERR(MI) .NE. 0) CYCLE OUTER

      ! If no vertex was found, then this is an extrapolation point.
      IF (SIMPS(D+1,MI) == 0) THEN
         ! If extrapolation is not allowed (EXTRAP_OPT=0), do not proceed.
         IF (EXTRAP < EPS) THEN
            ! Zero all output values.
            SIMPS(:,MI) = 0
            WEIGHTS(:,MI) = 0
            ! Set the error flag and skip this point.
            IERR(MI) = 2
            CYCLE OUTER
         END IF

         ! Otherwise, project the extrapolation point onto the convex hull.
         CALL PROJECT(D, N, PTS, PROJ, RNORM, IERR(MI))
         IF (IERR(MI) .NE. 0) CYCLE OUTER

         ! Check the value of RNORM for over-extrapolation.
         IF (RNORM > EXTRAP) THEN
            ! Zero all output values.
            SIMPS(:,MI) = 0
            WEIGHTS(:,MI) = 0
            ! Set the error flag and skip this point.
            IERR(MI) = 2
            CYCLE OUTER
         END IF

         ! Otherwise, remake the previous simplex and continue with the 
         ! projected value.
         SIMPS(D+1,MI) = ITMP
         A(D,:) = PTS(:,ITMP) - PTS(:,SIMPS(1,MI))
         B(D) = DDOT(D, A(D,:), 1, A(D,:), 1) / 2.0_R8
      END IF

   ! End of inner loop for finding each interpolation point.
   END DO INNER

   ! Check for budget violation conditions.
   IF (K > IBUDGET) THEN 
      ! Zero all output values.
      SIMPS(:,MI) = 0
      WEIGHTS(:,MI) = 0
      ! Set the error flag and skip this point.
      IERR(MI) = 60
      CYCLE OUTER
   END IF

   ! If the residual is nonzero, set the extrapolation flag.
   IF (RNORM > EPS) THEN
      IERR(MI) = 1
   END IF

   ! If present, record the RNORM_OPT output.
   IF (PRESENT(RNORM_OPT)) THEN
      RNORM_OPT(MI) = RNORM
   END IF

! End of outer loop over all interpolation points.
END DO OUTER

! Compute interpolation point response values.
IF (PRESENT(INTERP_IN_OPT)) THEN
   ! Loop over all interpolation points.
   DO MI = 1, M
      ! Check for errors.
      IF (IERR(MI) .LE. 1) THEN
         ! Compute the weighted sum of vertex response values.
         DO K = 1, D+1
            INTERP_OUT_OPT(:,MI) = INTERP_OUT_OPT(:,MI) &
             + INTERP_IN_OPT(:,SIMPS(K,MI)) * WEIGHTS(K,MI)
         END DO
      END IF
   END DO
END IF

RETURN

! Internal subroutines and functions.
CONTAINS

SUBROUTINE MAKEFIRSTSIMP()
! Iteratively construct the first simplex by choosing points that minimize the
! radius of the smallest circumball. Let ( P_1, P_2, ..., P_K ) denote the
! current list of vertices for the simplex. Let P* denote the candidate vertex
! to be added to the simplex. Let CENTER denote the circumcenter of the simplex. 
! Then
!
! X = CENTER - P_1
!
! is given by the minimum norm solution to the underdetermined linear system 
!
! AX = B, where
!
! A = [ P_2 - P_1; P_3 - P_1; ... ; P* - P_1 ] and
! B = [ < A_1, A_1 > / 2; < A_2, A_2 > / 2; ... ; < A_D, A_D > / 2 ]^T.
!
! Then the radius of the circumsphere can be computed from CURRRAD = || X ||.
! Then the next vertex is given by P_{K+1} = argmin_{P*} CURRRAD, where P* can
! be any point in PTS that is not already a vertex of the simplex.
!
! On output, this subroutine fully populates the matrix A and vector B, and
! fills SIMPS(:,MI) with the indices of a valid Delaunay simplex. This
! subroutine depends on a variation for handling degenerate or near-degenerate 
! input cases.

! Find the first point, i.e., the closest point to Q(:,MI).
SIMPS(:,MI) = 0
MINRAD = HUGE(MINRAD)
DO I = 1, N
   ! Check the distance to Q(:,MI)
   CURRRAD = DNRM2(D, PTS(:,I) - PROJ(:), 1)
   IF (CURRRAD < MINRAD) THEN
      MINRAD = CURRRAD
      SIMPS(1,MI) = I
   END IF
END DO
! Find the second point, i.e., the closest point to PTS(:,SIMP(1,MI)).
MINRAD = HUGE(MINRAD)
DO I = 1, N
   ! Skip repeated vertices.
   IF (I == SIMPS(1,MI)) CYCLE
   ! Check the diameter of the resulting circumsphere.
   CURRRAD = DNRM2(D, PTS(:,I)-PTS(:,SIMPS(1,MI)), 1)
   IF (CURRRAD < MINRAD) THEN
      MINRAD = CURRRAD
      SIMPS(2,MI) = I
   END IF
END DO
! Set up the first row of the LS system AX=B.
A(1,:) = PTS(:,SIMPS(2,MI)) - PTS(:,SIMPS(1,MI))
B(1) = DDOT(D, A(1,:), 1, A(1,:), 1) / 2.0_R8
! Loop to collect the remaining D+1 vertices for the first simplex.
DO I = 2, D
   MINRAD = HUGE(MINRAD) ! Re-initialize the radius for each iteration.
   ! Check each point P* in PTS.
   DO J = 1, N
      ! Check that this point is not already in the simplex.
      IF (ANY(SIMPS(:,MI) == J)) CYCLE
      ! Add P* to LS system, and compute the minimum norm solution.
      A(I,:) = PTS(:,J) - PTS(:,SIMPS(1,MI))
      B(I) = DDOT(D, A(I,:), 1, A(I,:), 1) / 2.0_R8
      LQ(1:I,:) = A(1:I,:)
      CENTER(1:I) = B(1:I)
      CALL DGELS('N', I, D, 1, LQ, D, CENTER, D, WORK, LWORK, IERR(MI))
      IF (IERR(MI) < 0) THEN ! LAPACK illegal input errors.
         IERR(MI) = 80
         RETURN
      ELSE IF (IERR(MI) > 0) THEN ! Errors caused by rank-deficiency.
         CYCLE ! If rank-deficient, skip this point.
      END IF
      ! Calculate the radius and compare it to the current minimum.
      CURRRAD = DNRM2(D, CENTER, 1)
      IF (CURRRAD < MINRAD) THEN
         MINRAD = CURRRAD
         SIMPS(I+1,MI) = J
      END IF
   END DO
   ! Check that a point was found. If not, skip to rank-deficiency code.
   IF (SIMPS(I+1,MI) == 0) THEN
      EXIT
   END IF
   ! If all operations were successful, add the best P* to the LS system.
   A(I,:) = PTS(:,SIMPS(I+1,MI)) - PTS(:,SIMPS(1,MI))
   B(I) = DDOT(D, A(I,:), 1, A(I,:), 1) / 2.0_R8 
END DO
! Perform a check for rank-deficiency using a singular-value decomposition.
IF (I > D) THEN
   ! Compute the SVD
   LQ = A
   CALL DGESVD('N','N',D,D,LQ,D,PLANE(1:D),U,1,VT,1,WORK,LWORK,IERR(MI))
   IF (IERR(MI) < 0) THEN ! Usage errors.
      IERR(MI) = 81
      RETURN
   ELSE IF (IERR(MI) > 0) THEN ! Failure to converge.
      IERR(MI) = 85
      RETURN
   END IF
   ! Check for rank-deficiency up to working precision.
   IF (PLANE(D)/PLANE(1) < EPS) THEN
      CALL MAKEFIRSTSIMP_RANKDEF() ! Redo entire construction using SVDs.
      IF (IERR(MI) .NE. 0) RETURN
   END IF
ELSE ! If I .LEQ. D, then rank-deficiency was already detected.
   CALL MAKEFIRSTSIMP_RANKDEF() ! Redo entire construction using SVDs.
   IF (IERR(MI) .NE. 0) RETURN
END IF
! Set error flag to 'success' for a normal return.
IERR(MI) = 0
RETURN
END SUBROUTINE MAKEFIRSTSIMP

SUBROUTINE MAKEFIRSTSIMP_RANKDEF()
! Iteratively construct the first simplex by choosing points that minimize the
! radius of the smallest circumball for degenerate and near degenerate inputs. 
! This subroutine is called by MAKEFIRSTSIMP(). Let ( P_1, P_2, ..., P_K ), P*, 
! A, B, CENTER, and X be defined as in MAKEFIRSTSIMP(). Assume that SIMPS(1,MI)
! and SIMPS(2,MI) have already been selected, and A(:,1) and B(1) have been
! filled appropriately.
!
! On output, this subroutine fully populates the matrix A and vector B, and
! fills SIMPS(:,MI) with the indices of a valid Delaunay simplex. This
! subroutine uses the LAPACK routine DGELSS to perform a singular-value
! decomposition during the evaluation of each P*. Using the ratio between the
! smallest and largest singular values, the condition number is judged for each
! P*. Any P* that results in a nearly singular simplex is skipped, making this
! version robust against (nearly) degenerate data.

SIMPS(3:D+1,MI) = 0 ! Zero vertices 3 through D+1 from MAKEFIRSTSIMP().
! Loop over the remaining D-1 points in the first simplex.
DO I = 2, D
   MINRAD = HUGE(MINRAD) ! Re-initialize the radius for each iteration.
   ! Check each point P* in PTS.
   DO J = 1, N
      ! Check that this point is not already in the simplex.
      IF (ANY(SIMPS(:,MI) == J)) CYCLE
      ! Add P* to LS system, and compute the minimum norm solution using SVD.
      A(I,:) = PTS(:,J) - PTS(:,SIMPS(1,MI))
      B(I) = DDOT(D, A(I,:), 1, A(I,:), 1) / 2.0_R8
      LQ(1:I,:) = A(1:I,:)
      CENTER(1:I) = B(1:I)
      CALL DGELSS(I, D, 1, LQ, D, CENTER, D, PLANE, EPS, ITMP, WORK, &
       LWORK, IERR(MI))
      IF (IERR(MI) < 0) THEN ! Usage errors.
         IERR(MI) = 82
         RETURN
      ELSE IF (IERR(MI) > 0) THEN ! Failure to converge.
         IERR(MI) = 86
         RETURN
      END IF
      IF (ITMP < I) CYCLE ! If A was not full-rank, then skip P*.
      ! Calculate the radius, and compare to the current minimum.
      CURRRAD = DNRM2(D, CENTER, 1)
      IF (CURRRAD < MINRAD) THEN
         MINRAD = CURRRAD
         SIMPS(I+1,MI) = J
      END IF
   END DO
   ! Check that a point was found. If not, then all the points must lie in a
   ! lower dimensional linear manifold (error case).
   IF (SIMPS(I+1,MI) == 0) THEN
      IERR(MI) = 31
      RETURN
   END IF
   ! If all operations successful, add the best P* to the LS system.
   A(I,:) = PTS(:,SIMPS(I+1,MI)) - PTS(:,SIMPS(1,MI))
   B(I) = DDOT(D, A(I,:), 1, A(I,:), 1) / 2.0_R8 
END DO
! Set error flag to 'success' for a normal return.
IERR(MI) = 0
RETURN
END SUBROUTINE MAKEFIRSTSIMP_RANKDEF

SUBROUTINE MAKESIMPLEX()
! Given a Dalaunay facet F, complete the simplex by adding a point from the same
! `side' of F as Q(:,MI). Assume SIMPS(1:D,MI) contains the vertices of F (P_1,
! P_2, ..., P_D), and assume the matrix A(1:D-1,:) and vector B(1:D-1) are
! filled appropriately (similarly as in MAKEFIRSTSIMP()). Then for any P* in
! PTS, let CENTER denote the circumcenter of the simplex with vertices (P_1,
! P_2, ..., P_D, P*). Then
!
! X = CENTER - P_1
!
! is given by the solution to the exactly determined linear system
!
! AX = B where
!
! A = [ P_2 - P_1; P_3 - P_1; ... ; P* - P_1 ] and
! B = [ < A_1, A_1 > / 2; < A_2, A_2 > / 2; ... ; < A_D, A_D > / 2 ]^T.
!
! Then CENTER = X + P_1 and RADIUS = || X ||. P_{D+1} will be given by the 
! candidate P* which satisfies both of the following:
!
! 1) Let PLANE denote the hyperplane containing F. Then P_{D+1} and Q(:,MI) 
! must be on the same side of PLANE.
!
! 2) The circumball about P must not contain any other points in PTS in its
! interior (Delaunay property).
! 
! The above are necessarry and sufficient conditions for flipping the Delaunay
! simplex, given that F is indeed a Delaunay facet.
!
! Upon output, SIMPS(:,MI) will contain a Delaunay simplex closer to Q(:,MI).
! Also, the matrix A and vector B will be updated accordingly. If
! SIMPS(D+1,MI)=0, then there were no points in PTS on the appropriate side of
! F, meaning that Q(:,MI) is an extrapolation point.

! Compute the hyperplane PLANE.
CALL MAKEPLANE()
IF(IERR(MI) .NE. 0) RETURN ! Check for errors.
! Compute the sign for the side of PLANE containing Q(:,MI).
SIDE1 =  DDOT(D,PLANE(1:D),1,PROJ(:),1) - PLANE(D+1)
SIDE1 = SIGN(1.0_R8,SIDE1)
! Initialize the center, radius, and simplex.
SIMPS(D+1,MI) = 0
CENTER(:) = 0_R8
MINRAD = HUGE(MINRAD)
! Loop through all points P* in PTS.
DO I = 1, N
   ! Check that P* is inside the current ball.
   IF (DNRM2(D, PTS(:,I) - CENTER(:), 1) > MINRAD) CYCLE ! If not, skip.
   ! Check that P* is on the appropriate halfspace.
   SIDE2 = DDOT(D,PLANE(1:D),1,PTS(:,I),1) - PLANE(D+1)
   IF (SIDE1 * SIDE2 < EPS .OR. ANY(SIMPS(:,MI) == I)) CYCLE ! If not, skip.
   ! Add P* to the linear system, and solve to get shifted CENTER.
   A(D,:) = PTS(:,I) - PTS(:,SIMPS(1,MI))
   B(D) = DDOT(D, A(D,:), 1, A(D,:), 1) / 2.0_R8
   LQ = A
   CENTER = B
   CALL DGESV(D, 1, LQ, D, IPIV, CENTER, D, IERR(MI))
   IF (IERR(MI) < 0) THEN ! LAPACK illegal input error.
      IERR(MI) = 83
      RETURN
   ELSE IF (IERR(MI) > 0) THEN ! Rank-deficiency detected.
      IERR(MI) = 61
      RETURN
   END IF
   ! Update the new radius, center, and simplex.
   CURRRAD = DNRM2(D, CENTER, 1)
   MINRAD = CURRRAD
   CENTER(:) = CENTER(:) + PTS(:,SIMPS(1,MI))
   SIMPS(D+1,MI) = I
END DO
IERR(MI) = 0 ! Reset the error flag to 'success' code.
! Check for extrapolation condition.
IF(SIMPS(D+1,MI) == 0) RETURN
! Add new point to the linear system.
A(D,:) = PTS(:,SIMPS(D+1,MI)) - PTS(:,SIMPS(1,MI))
B(D) = DDOT(D, A(D,:), 1, A(D,:), 1) / 2.0_R8
RETURN
END SUBROUTINE MAKESIMPLEX

SUBROUTINE MAKEPLANE()
! Construct a plane intersecting the first D vertices in SIMPS(:,MI). The plane
! is determined by its normal vector and an intercept. Let (P_1, P_2, ..., P_D)
! be the vertices of SIMPS(:,MI). To find the normal vector, search the 
! nullspace of the matrix
!
! A = [ P_2 - P_1; P_3 - P_1; ... ; P_D - P_1 ].
!
! Any nontrivial vector in ker(A) is a valid normal vector. To find, such a
! vector, solve AX = 0, subject to the contstraint that X_D = 1.
!
! Upon output, PLANE(1:D) contains the coefficients of the normal vector and
! PLANE(D+1) contains the intercept of the plane.

IF (D > 1) THEN ! Check that D-1 > 0, otherwise the plane is trivial.
   ! Compute the QR factorization.
   IPIV=0
   LQ = A
   CALL DGEQP3(D-1,D,LQ,D,IPIV,PLANE,WORK,LWORK,IERR(MI))
   IF(IERR(MI) < 0) THEN ! LAPACK illegal input error.
      IERR(MI) = 84
      RETURN
   ELSE IF (IERR(MI) > 0) THEN ! Rank-deficiency detected.
      IERR(MI) = 61
      RETURN
   END IF
   ! Do a triangular solve to get the normal vector, and store in a work array.
   WORK(1:D-1) = LQ(1:D-1,D)
   CALL DTRSM('L','U','N','N',D-1,1,-1.0_R8,LQ,D,WORK,D)
   WORK(D) = 1.0_R8
   ! Undo the pivots.
   DO I = 1,D
      PLANE(IPIV(I)) = WORK(I)
   END DO
   ! Normalize the normal vector to length 1.
   PLANE(1:D) = PLANE(1:D) / DNRM2(D,PLANE(1:D),1)
   ! Calculate the intercept of the plane.
   PLANE(D+1) = DDOT(D,PLANE(1:D),1,PTS(:,SIMPS(1,MI)),1)
ELSE ! Special case where D=1.
   PLANE(1) = 1.0_R8
   PLANE(2) = PTS(1,SIMPS(1,MI))
END IF
RETURN
END SUBROUTINE MAKEPLANE

FUNCTION PTINSIMP() RESULT(TF)
! Determine if any interpolation points are in the current simplex. The vertices
! of the simplex (P_1, P_2, ..., P_{D+1}) determine a cone of vectors given by
! V_I = P_{I+1} - P_1. Each interpolation point in Q can be expressed as a
! linear combination of V_I. If all the linear weights are positive for some Q*
! in Q, then Q* is in the convex cone. Furthermore, if all the weights sum to
! less than 1.0, then Q* is in the simplex with vertices \{P_I\}_{I=1}^{D+1}.
!
! If any interpolation points in Q are contained in SIMPS(:,MI), then those
! points are marked as solved and the values of SIMPS and WEIGHTS are updated
! appropriately. On output, WEIGHTS(:,MI) contains the affine weights for
! producing Q(:,MI) as an affine combintation of SIMPS(:,MI). If these weights
! are convex, then PTINSIMP() returns TRUE.

! Initialize the return value and local variables.
LOGICAL :: TF ! True/False value.
INTEGER :: NRHS ! Number of RHS for linear solve.

! Calculate the number of RHS for this system, and initialize NRHS.
NRHS = M + 1 - MI
TF = .FALSE.
! Prepare the appropriate linear system, using A^T and NRHS right-hand sides.
LQ = TRANSPOSE(A)
WORK(1:D) = PROJ(:) - PTS(:,SIMPS(1,MI))
DO I = MI+1, M
   WORK((I-MI)*D+1:(I-MI+1)*D) = Q(:,I) &
    - PTS(:,SIMPS(1,MI))
END DO
! Solve a single linear system to get the weights for all unsolved interpolation
! points.
CALL DGESV(D, NRHS, LQ, D, IPIV, WORK(1:NRHS*D), D, IERR(MI))
IF (IERR(MI) < 0) THEN ! LAPACK illegal input.
   IERR(MI) = 83
   RETURN
ELSE IF (IERR(MI) > 0) THEN ! Rank-deficiency detected.
   IERR(MI) = 61
   RETURN
END IF
! Get the affine weights for Q(:,MI).
WEIGHTS(2:D+1,MI) = WORK(1:D)
WEIGHTS(1,MI) = 1.0_R8 - SUM(WEIGHTS(2:D+1,MI))
! Check if the weights for Q(:,MI) are convex.
IF (ALL(WEIGHTS(:,MI) >= -EPS)) TF = .TRUE.
! Compute the affine weights for the rest of the interpolation points.
DO I = 1, NRHS-1
   IF (IERR(I+MI) == 40) THEN ! Check that no solution has already been found.
      ! Get the affine weights and simplices.
      WEIGHTS(2:D+1,I+MI) = WORK(D*I+1:D*(I+1))
      WEIGHTS(1,I+MI) = 1.0_R8 - SUM(WEIGHTS(2:D+1,I+MI))
      SIMPS(:, I+MI) = SIMPS(:,MI)
      ! Check if the weights define a convex combination.
      IF (ALL(WEIGHTS(:,MI+I) >= -EPS)) IERR(I+MI) = 0
   END IF
END DO
RETURN
END FUNCTION PTINSIMP

END SUBROUTINE DELAUNAYINTERP

!
! Auxiliary subroutines.
!
SUBROUTINE SCALETOBALL(D, N, PTS, M, Q, MINDIST)
! Rescale and transform data to be centered at the origin and contained in the
! unit ball. This subroutine has O(n^2) complexity.
!
!
! On input:
!
! D is the dimension of the space for PTS and Q.
!
! N is the number of data points in PTS.
!
! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
! coordinates of a single data point in R^D.
!
! M is the number of interpolation points in Q.
!
! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
! coordinates of a single interpolation point in R^D.
!
!
! On output:
!
! PTS and Q have been rescaled and shifted. All the data points in PTS are now
! contained in the unit hyperball in R^D, and the points in Q have been shifted
! and scaled accordingly in relation to PTS.
!
! MINDIST contains the minimum distance between any two data points in PTS.

! Input arguments.
INTEGER, INTENT(IN) :: D, N, M
! In/Out arguments.
REAL(KIND=R8), INTENT(INOUT) :: PTS(D,N), Q(D,M)
! Output arguments.
REAL(KIND=R8), INTENT(OUT) :: MINDIST

! Local variables.
INTEGER :: I, J ! Loop iteration variables.
REAL(KIND=R8) :: CENTER(D) ! The center of the data points PTS.
REAL(KIND=R8) :: DIAMETER ! The diameter of the data points PTS.
REAL(KIND=R8) :: DISTANCE ! The current distance.

! External functions.
REAL(KIND=R8), EXTERNAL :: DNRM2 ! Compute the Euclidean norm (BLAS).

! Initialize local values.
MINDIST = HUGE(MINDIST)
DIAMETER = 0.0_R8
CENTER(:) = 0.0_R8

! Compute the minimum and maximum distances and the center of the PTS.
DO I = 1, N ! Cycle through all pairs of points.
   DO J = I+1, N
      DISTANCE = DNRM2(D, PTS(:,I) - PTS(:,J), 1) ! Compute the distance.
      IF (DISTANCE > DIAMETER) THEN ! Compare to the current diameter.
         DIAMETER = DISTANCE
      END IF
      IF (DISTANCE < MINDIST) THEN ! Compare to the current minimum distance.
         MINDIST = DISTANCE
      END IF
   END DO
   CENTER(:) = CENTER(:) + PTS(:,I) ! Aggregate all values in CENTER.
END DO
CENTER = CENTER / REAL(N, KIND=R8) ! Compute the center.
! Now, center and scale PTS.
DO I = 1, N
   PTS(:,I) = 2.0_R8 * (PTS(:,I) - CENTER(:)) / DIAMETER + 1.0_R8
END DO
! Also, transform Q appropriately.
DO I = 1, M
   Q(:,I) = 2.0_R8 * (Q(:,I) - CENTER(:)) / DIAMETER + 1.0_R8
END DO
! Also, scale the minimum distance.
MINDIST = MINDIST / DIAMETER
RETURN
END SUBROUTINE SCALETOBALL

SUBROUTINE PROJECT(D, N, PTS, PROJ, RNORM, IERR)
! Project a point outside the convex hull of the point set onto the convex hull
! by solving an inequality constrained least squares problem. The solution to
! the least squares problem gives the projection as a convex combination of the
! data points. The projection can then be computed by performing a matrix
! vector multiplication.
!
!
! On input:
!
! D is the dimension of the space for PTS and Q.
!
! N is the number of data points in PTS.
!
! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
! coordinates of a single data point in R^D.
!
! PROJ(1:D) is a real valued vector containing the coordinates of a single point
! in R^D, to be projected onto the convex hull of PTS.
!
!
! On output:
!
! PROJ(1:D) has been projected onto the convex hull of PTS and now contains the
! coordinates of the projected value.
!
! RNORM contains the Euclidean-norm of the residual vector from the projection.
!
! IERR contains an integer valued error flag, consistent with the error codes
! for VTDelaunay. The error codes are listed below:
!
! 00 : Succesful projection.
!
! 10 : The given dimension D was not valid.
! 11 : Too few data points to construct a triangulation (I.e., N < D+1).
! 13 : The first dimension of PTS does not agree with the dimension D.
! 14 : The second dimension of PTS does not agree with the number of points N.
! 15 : The dimension of PROJ does not agree with the dimension D.
!
! 70 : The SLATEC subroutine DWNNLS failed to converge.
! 71 : The SLATEC subroutine DWNNLS has reported a usage error. Possible causes:
!      An automatic array failed to allocate.

! Input arguments.
INTEGER, INTENT(IN) :: D, N
REAL(KIND=R8), INTENT(IN) :: PTS(:,:)
! Output arguments.
REAL(KIND=R8), INTENT(INOUT) :: PROJ(:)
REAL(KIND=R8), INTENT(OUT) :: RNORM
INTEGER, INTENT(OUT) :: IERR

! Local variables.
INTEGER :: I, J ! Loop iteration variables.
INTEGER :: IWORK(D+1+N) ! Integer work array for DWNNLS.
REAL(KIND=R8) :: A(D+1,N+1) ! Input array for DWNNLS.
REAL(KIND=R8) :: PRGOPT(1) ! Program options for DWNNLS. 
REAL(KIND=R8) :: X(N) ! Output array for DWNNLS.
REAL(KIND=R8) :: WORK(D+1+N*5) ! Real valued work array for DWNNLS.

! External subroutines.
EXTERNAL :: DGEMV ! General matrix vector multiply (BLAS)
! DWNNLS subroutine from SLATEC library. Solves an inequality constrained least
! squares problem. For a reference, see ACM TOMS Algorithm 587.
EXTERNAL :: DWNNLS 

! Check for legal inputs. When called by VTDelaunay, this is redundant. However,
! this subroutine can also be called separately.
IF ( D < 1 ) THEN 
   IERR = 10
   RETURN
END IF
IF ( N < D+1 ) THEN
   IERR = 11
   RETURN
END IF
IF ( SIZE(PTS,1) .NE. D ) THEN
   IERR = 13
   RETURN
END IF
IF ( SIZE(PTS,2) .NE. N ) THEN
   IERR = 14
   RETURN
END IF
IF ( SIZE(PROJ,1) .NE. D ) THEN
   IERR = 15
   RETURN
END IF

! Initialize work arrays.
PRGOPT(1) = 1.0_R8
IWORK(1) = D+1+5*N
IWORK(2) = D+1+N
! Set convexity (equality) constraint.
A(1, :) = 1.0_R8
! Copy data points.
A(2:D+1,1:N) = PTS(:,:)
! Copy extrapolation point.
A(2:D+1,N+1) = PROJ(:)
! Compute the solution to the inequality constrained least squares problem to
! get the projection coefficients.
CALL DWNNLS(A,D+1,1,D,N,0,PRGOPT,X,RNORM,IERR,IWORK,WORK)
IF (IERR .EQ. 1) THEN ! Failure to converge.
   IERR = 70
   RETURN
ELSE IF (IERR .EQ. 2) THEN ! Illegal input detected.
   IERR = 71
   RETURN
END IF
! Compute the actual projection via matrix vector multiplication.
CALL DGEMV('N',D,N,1.0_R8,PTS,D,X,1,0.0_R8,PROJ,1)
RETURN
END SUBROUTINE PROJECT

END MODULE VTDELAUNAY_MOD
