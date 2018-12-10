MODULE REAL_PRECISION  ! HOMPACK90 module for 64-bit arithmetic.
INTEGER, PARAMETER:: R8=SELECTED_REAL_KIND(13)
END MODULE REAL_PRECISION

MODULE DELSPARSE_MOD
! This module contains the REAL_PRECISION R8 data type for 64-bit arithmetic
! and interface blocks for the DELAUNAYSPARSES and DELAUNAYSPARSEP 
! subroutines for computing the Delaunay simplices containing interpolation
! points Q in R^D given data points PTS.
USE REAL_PRECISION
PUBLIC

INTERFACE
   ! Interface for serial subroutine DELAUNAYSPARSES.
   SUBROUTINE DELAUNAYSPARSES( D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR,     &
                               INTERP_IN, INTERP_OUT, EPS, EXTRAP, RNORM, &
                               IBUDGET, CHAIN                             )
      USE REAL_PRECISION, ONLY : R8
      INTEGER, INTENT(IN) :: D, N
      REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
      INTEGER, INTENT(IN) :: M
      REAL(KIND=R8), INTENT(INOUT) :: Q(:,:)
      INTEGER, INTENT(OUT) :: SIMPS(:,:)
      REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
      INTEGER, INTENT(OUT) :: IERR(:)
      REAL(KIND=R8), INTENT(IN), OPTIONAL:: INTERP_IN(:,:)
      REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
      REAL(KIND=R8), INTENT(IN), OPTIONAL:: EPS, EXTRAP 
      REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
      INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
      LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
   END SUBROUTINE DELAUNAYSPARSES

   ! Interface for parallel subroutine DELAUNAYSPARSEP.
   SUBROUTINE DELAUNAYSPARSEP( D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR,     &
                               INTERP_IN, INTERP_OUT, EPS, EXTRAP, RNORM, &
                               IBUDGET, CHAIN, PMODE, CHUNKSIZE           )
      USE REAL_PRECISION, ONLY : R8
      INTEGER, INTENT(IN) :: D, N
      REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:)
      INTEGER, INTENT(IN) :: M
      REAL(KIND=R8), INTENT(INOUT) :: Q(:,:)
      INTEGER, INTENT(OUT) :: SIMPS(:,:)
      REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
      INTEGER, INTENT(OUT) :: IERR(:)
      REAL(KIND=R8), INTENT(IN), OPTIONAL:: INTERP_IN(:,:)
      REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
      REAL(KIND=R8), INTENT(IN), OPTIONAL:: EPS, EXTRAP 
      REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
      INTEGER, INTENT(IN), OPTIONAL :: IBUDGET, PMODE, CHUNKSIZE
      LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
   END SUBROUTINE DELAUNAYSPARSEP

   ! Interface for SLATEC subroutine DWNNLS.
   SUBROUTINE DWNNLS( W, MDW, ME, MA, N, L, PRGOPT, X, RNORM,  &
                               MODE, IWORK, WORK )
      USE REAL_PRECISION, ONLY : R8
      INTEGER :: IWORK(*), L, MA, MDW, ME, MODE, N
      REAL(KIND=R8) :: PRGOPT(*), RNORM, W(MDW,*), WORK(*), X(*)
   END SUBROUTINE DWNNLS

END INTERFACE

END MODULE DELSPARSE_MOD

SUBROUTINE DELAUNAYSPARSES( D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, &
  INTERP_IN, INTERP_OUT, EPS, EXTRAP, RNORM, IBUDGET, CHAIN )
! This is a serial implementation of an algorithm for efficiently performing
! interpolation in R^D via the Delaunay triangulation. The algorithm is fully
! described and analyzed in
!
! T. H. Chang, L. T. Watson,  T. C.H. Lux, B. Li, L. Xu, A. R. Butt, K. W.
! Cameron, and Y. Hong. 2018. A polynomial time algorithm for multivariate
! interpolation in arbitrary dimension via the Delaunay triangulation. In
! Proceedings of the ACMSE 2018 Conference (ACMSE '18). ACM, New York, NY,
! USA. Article 12, 8 pages.
!
!
! On input:
!
! D is the dimension of the space for PTS and Q.
!
! N is the number of data points in PTS.
!
! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
!    coordinates of a single data point in R^D.
!
! M is the number of interpolation points in Q.
!
! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
!    coordinates of a single interpolation point in R^D.
!
!
! On output:
!
! PTS and Q have been rescaled and shifted. All the data points in PTS
!    are now contained in the unit hyperball in R^D, and the points in Q
!    have been shifted and scaled accordingly in relation to PTS.
!
! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
!    in PTS) for the D+1 vertices of the Delaunay simplex containing each
!    interpolation point in Q.
!
! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
!    point in Q as a convex combination of the D+1 corresponding vertices
!    in SIMPS.
!
! IERR(1:M) contains integer valued error flags associated with the
!    computation of each of the M interpolation points in Q. The error
!    codes are:
!
! 00 : Succesful interpolation.
! 01 : Succesful extrapolation (up to the allowed extrapolation distance).
! 02 : This point was outside the allowed extrapolation distance; the
!      corresponding entries in SIMPS and WEIGHTS contain zero values.
!
! 10 : The dimension D must be positive.
! 11 : Too few data points to construct a triangulation (i.e., N < D+1).
! 12 : No interpolation points given (i.e., M < 1).
! 13 : The first dimension of PTS does not agree with the dimension D.
! 14 : The second dimension of PTS does not agree with the number of points N.
! 15 : The first dimension of Q does not agree with the dimension D.
! 16 : The second dimension of Q does not agree with the number of
!      interpolation points M.
! 17 : The first dimension of the output array SIMPS does not match the number
!      of vertices needed for a D-simplex (D+1).
! 18 : The second dimension of the output array SIMPS does not match the
!      number of interpolation points M.
! 19 : The first dimension of the output array WEIGHTS does not match the
!      number of vertices for a a D-simplex (D+1).
! 20 : The second dimension of the output array WEIGHTS does not match the
!      number of interpolation points M.
! 21 : The size of the error array IERR does not match the number of
!      interpolation points M.
! 22 : INTERP_IN cannot be present without INTERP_OUT or vice versa.
! 23 : The first dimension of INTERP_IN does not match the first
!      dimension of INTERP_OUT.
! 24 : The second dimension of INTERP_IN does not match the number of
!      data points PTS.
! 25 : The second dimension of INTERP_OUT does not match the number of
!      interpolation points M.
! 26 : The budget supplied in IBUDGET does not contain a positive
!      integer.
! 27 : The extrapolation distance supplied in EXTRAP cannot be negative.
! 28 : The size of the RNORM output array does not match the number of
!      interpolation points M.
!
! 30 : Two or more points in the data set PTS are too close together with
!      respect to the working precision (EPS), which would result in a
!      numerically degenerate simplex.
! 31 : All the data points in PTS lie in some lower dimensional linear
!      manifold (up to the working precision), and no valid triangulation
!      exists.
! 40 : An error caused DELAUNAYSPARSES to terminate before this value could
!      be computed. Note: The corresponding entries in SIMPS and WEIGHTS may
!      contain garbage values.
!
! 50 : A memory allocation error occurred while allocating the work array
!      WORK.
!
! 60 : The budget was exceeded before the algorithm converged on this
!      value. If the dimension is high, try increasing IBUDGET. This
!      error can also be caused by a working precision EPS that is too
!      small for the conditioning of the problem.
!
! 61 : A value that was judged appropriate later caused LAPACK to encounter a
!      singularity. Try increasing the value of EPS.
!
! 70 : Allocation error for the extrapolation work arrays.
! 71 : The SLATEC subroutine DWNNLS failed to converge during the projection
!      of an extrapolation point onto the convex hull.
! 72 : The SLATEC subroutine DWNNLS has reported a usage error.
!
!      The errors 72, 80--83 should never occur, and likely indicate a
!      compiler bug or hardware failure.
! 80 : The LAPACK subroutine DGEQP3 has reported an illegal value.
! 81 : The LAPACK subroutine DGETRF has reported an illegal value.
! 82 : The LAPACK subroutine DGETRS has reported an illegal value.
! 83 : The LAPACK subroutine DORMQR has reported an illegal value.
!
! Optional arguments:
!
! INTERP_IN(1:IR,1:N) contains real valued response vectors for each of
!    the data points in PTS on input. The first dimension of INTERP_IN is
!    inferred to be the dimension of these response vectors, and the
!    second dimension must match N. If present, the response values will
!    be computed for each interpolation point in Q, and stored in INTERP_OUT,
!    which therefore must also be present. If both INTERP_IN and INTERP_OUT
!    are omitted, only the containing simplices and convex combination
!    weights are returned.
! 
! INTERP_OUT(1:IR,1:M) contains real valued response vectors for each
!    interpolation point in Q on output. The first dimension of INTERP_OUT
!    must match the first dimension of INTERP_IN, and the second dimension
!    must match M. If present, the response values at each interpolation
!    point are computed as a convex combination of the response values
!    (supplied in INTERP_IN) at the vertices of a Delaunay simplex containing
!    that interpolation point.  Therefore, if INTERP_OUT is present, then
!    INTERP_IN must also be present.  If both are omitted, only the
!    simplices and convex combination weights are returned.
! 
! EPS contains the real working precision for the problem on input. By default,
!    EPS is assigned \sqrt{\mu} where \mu denotes the unit roundoff for the
!    machine. In general, any values that differ by less than EPS are judged
!    as equal, and any weights that are greater than -EPS are judged as
!    nonnegative.  EPS cannot take a value less than the default value of
!    \sqrt{\mu}. If any value less than \sqrt{\mu} is supplied, the default
!    value will be used instead automatically. 
! 
! EXTRAP contains the real maximum extrapolation distance (relative to the
!    diameter of PTS) on input. Interpolation at a point outside the convex
!    hull of PTS is done by projecting that point onto the convex hull, and
!    then doing normal Delaunay interpolation at that projection.
!    Interpolation at any point in Q that is more than EXTRAP * DIAMETER(PTS)
!    units outside the convex hull of PTS will not be done and an error code
!    of 2 will be returned. Note that computing the projection can be
!    expensive. Setting EXTRAP=0 will cause all extrapolation points to be
!    ignored without ever computing a projection. By default, EXTRAP=0.1
!    (extrapolate by up to 10% of the diameter of PTS). 
! 
! RNORM(1:M) contains the real unscaled projection (2-norm) distances from
!    any projection computations on output. If not present, these distances
!    are still computed for each extrapolation point, but are never returned.
! 
! IBUDGET on input contains the integer budget for performing flips while
!    iterating toward the simplex containing each interpolation point in
!    Q. This prevents DELAUNAYSPARSES from falling into an infinite loop when
!    an inappropriate value of EPS is given with respect to the problem
!    conditioning.  By default, IBUDGET=50000. However, for extremely
!    high-dimensional problems and pathological inputs, the default value
!    may be insufficient. 
!
! CHAIN is a logical input argument that determines whether a new first
!    simplex should be constructed for each interpolation point
!    (CHAIN=.FALSE.), or whether the simplex walks should be "daisy-chained."
!    By default, CHAIN=.FALSE. Setting CHAIN=.TRUE. is generally not
!    recommended, unless the size of the triangulation is relatively small
!    or the interpolation points are known to be tightly clustered.
!
!
! Subroutines and functions directly referenced from BLAS are
!      DDOT, DGEMV, DNRM2, DTRSM,
! and from LAPACK are
!      DGEQP3, DGETRF, DGETRS, DORMQR.
! The SLATEC subroutine DWNNLS is directly referenced. DWNNLS and all its
! SLATEC dependencies have been slightly edited to comply with the Fortran
! 2008 standard, with all print statements and references to stderr being
! commented out. For a reference to DWNNLS, see ACM TOMS Algorithm 587
! (Hanson and Haskell).  The module REAL_PRECISION from HOMPACK90 (ACM TOMS
! Algorithm 777) is used for the real data type. The REAL_PRECISION module,
! DELAUNAYSPARSES, and DWNNLS and its dependencies comply with the Fortran
! 2008 standard.  
! 
! Primary Author: Tyler H. Chang
! Last Update: September, 2018
! 
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE

! Input arguments.
INTEGER, INTENT(IN) :: D, N
REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:) ! Rescaled on output.
INTEGER, INTENT(IN) :: M
REAL(KIND=R8), INTENT(INOUT) :: Q(:,:) ! Rescaled on output.
! Output arguments.
INTEGER, INTENT(OUT) :: SIMPS(:,:)
REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
INTEGER, INTENT(OUT) :: IERR(:)
! Optional arguments.
REAL(KIND=R8), INTENT(IN), OPTIONAL:: INTERP_IN(:,:)
REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
REAL(KIND=R8), INTENT(IN), OPTIONAL:: EPS, EXTRAP 
REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
LOGICAL, INTENT(IN), OPTIONAL :: CHAIN 

! Local copies of optional input arguments.
REAL(KIND=R8) :: EPSL, EXTRAPL
INTEGER :: IBUDGETL
LOGICAL :: CHAINL

! Local variables.
INTEGER :: I, J, K ! Loop iteration variables.
INTEGER :: ITMP ! Temporary variable for swapping, looping, etc.
INTEGER :: LWORK ! Size of the work array.
INTEGER :: MI ! Index of current interpolation point.
REAL(KIND=R8) :: CURRRAD ! Radius of the current circumsphere.
REAL(KIND=R8) :: MINRAD ! Minimum circumsphere radius observed.
REAL(KIND=R8) :: PTS_SCALE ! Data scaling factor.
REAL(KIND=R8) :: PTS_DIAM ! Scaled diameter of data set.
REAL(KIND=R8) :: RNORML ! Euclidean norm of the projection residual.
REAL(KIND=R8) :: SIDE1, SIDE2 ! Signs (+/-1) denoting sides of a facet.

! Local arrays, requiring O(d^2) additional memory.
INTEGER :: IPIV(D) ! Pivot indices.
INTEGER :: SEED(D+1) ! Copy of the SEED simplex. Only used if CHAIN = .TRUE.
REAL(KIND=R8) :: A(D,D) ! The LHS of a linear system.
REAL(KIND=R8) :: B(D) ! The RHS of a linear system. 
REAL(KIND=R8) :: CENTER(D) ! The circumcenter of a simplex. 
REAL(KIND=R8) :: LQ(D,D) ! The LU or QR factorization of A.
REAL(KIND=R8) :: PLANE(D+1) ! The hyperplane containing a facet.
REAL(KIND=R8) :: PRGOPT_DWNNLS(1) ! Options array for DWNNLS.
REAL(KIND=R8) :: PROJ(D) ! The projection of the current iterate.
REAL(KIND=R8) :: TAU(D) ! Householder reflector constants.

! Extrapolation work arrays are only allocated if DWNNLS is called.
INTEGER, ALLOCATABLE :: IWORK_DWNNLS(:) ! Only for DWNNLS.
REAL(KIND=R8), ALLOCATABLE :: W_DWNNLS(:,:) ! Only for DWNNLS.
REAL(KIND=R8), ALLOCATABLE :: WORK(:) ! Allocated with size LWORK.
REAL(KIND=R8), ALLOCATABLE :: WORK_DWNNLS(:) ! Only for DWNNLS.
REAL(KIND=R8), ALLOCATABLE :: X_DWNNLS(:) ! Only for DWNNLS.

! External functions and subroutines.
REAL(KIND=R8), EXTERNAL :: DDOT  ! Inner product (BLAS).
REAL(KIND=R8), EXTERNAL :: DNRM2 ! Euclidean norm (BLAS).
EXTERNAL :: DGEMV  ! General matrix vector multiply (BLAS)
EXTERNAL :: DGEQP3 ! Perform a QR factorization with column pivoting (LAPACK).
EXTERNAL :: DGETRF ! Perform a LU factorization with partial pivoting (LAPACK).
EXTERNAL :: DGETRS ! Use the output of DGETRF to solve a linear system (LAPACK).
EXTERNAL :: DORMQR ! Apply householder reflectors to a matrix (LAPACK).
EXTERNAL :: DTRSM  ! Perform a triangular solve (BLAS).
EXTERNAL :: DWNNLS ! Solve an inequality constrained least squares problem
                   ! (SLATEC).

! Check for input size and dimension errors.
IF (D < 1) THEN ! The dimension must satisfy D > 0.
   IERR(:) = 10; RETURN; END IF
IF (N < D+1) THEN ! Must have at least D+1 data points.
   IERR(:) = 11; RETURN; END IF
IF (M < 1) THEN ! Must have at least one interpolation point.
   IERR(:) = 12; RETURN; END IF
IF (SIZE(PTS,1) .NE. D) THEN ! Dimension of PTS array should match.
   IERR(:) = 13; RETURN; END IF
IF (SIZE(PTS,2) .NE. N) THEN ! Number of data points should match.
   IERR(:) = 14; RETURN; END IF
IF (SIZE(Q,1) .NE. D) THEN ! Dimension of Q should match.
   IERR(:) = 15; RETURN; END IF
IF (SIZE(Q,2) .NE. M) THEN ! Number of interpolation points should match.
   IERR(:) = 16; RETURN; END IF
IF (SIZE(SIMPS,1) .NE. D+1) THEN ! Need space for D+1 vertices per simplex.
   IERR(:) = 17; RETURN; END IF
IF (SIZE(SIMPS,2) .NE. M) THEN  ! There will be M output simplices.
   IERR(:) = 18; RETURN; END IF
IF (SIZE(WEIGHTS,1) .NE. D+1) THEN ! There will be D+1 weights per simplex.
   IERR(:) = 19; RETURN; END IF
IF (SIZE(WEIGHTS,2) .NE. M) THEN ! One vector of weights per simplex.
   IERR(:) = 20; RETURN; END IF
IF (SIZE(IERR) .NE. M) THEN ! An error flag for each interpolation point.
   IERR(:) = 21; RETURN; END IF

! Check for optional arguments.
IF (PRESENT(INTERP_IN) .NEQV. PRESENT(INTERP_OUT)) THEN
   IERR(:) = 22; RETURN; END IF
IF (PRESENT(INTERP_IN)) THEN ! Sizes must agree.
   IF (SIZE(INTERP_IN,1) .NE. SIZE(INTERP_OUT,1)) THEN
      IERR(:) = 23 ; RETURN; END IF
   IF(SIZE(INTERP_IN,2) .NE. N) THEN
      IERR(:) = 24; RETURN; END IF
   IF (SIZE(INTERP_OUT,2) .NE. M) THEN
      IERR(:) = 25; RETURN; END IF
   INTERP_OUT(:,:) = 0.0_R8 ! Initialize output to zeros.
END IF
EPSL = SQRT(EPSILON(0.0_R8)) ! Get the machine unit roundoff constant.
IF (PRESENT(EPS)) THEN
   IF (EPSL < EPS) THEN ! If the given precision is too small, ignore it.
      EPSL = EPS
   END IF
END IF
IF (PRESENT(IBUDGET)) THEN
   IBUDGETL = IBUDGET ! Use the given budget if present.
   IF (IBUDGETL < 1) THEN
      IERR(:) = 26; RETURN; END IF
ELSE
   IBUDGETL = 50000 ! Default value for budget.
END IF
IF (PRESENT(EXTRAP)) THEN
   EXTRAPL = EXTRAP 
   IF (EXTRAPL < 0) THEN ! Check that the extrapolation distance is legal.
      IERR(:) = 27; RETURN; END IF
ELSE
   EXTRAPL = 0.1_R8 ! Default extrapolation distance (for normalized points).
END IF
IF (PRESENT(RNORM)) THEN
   IF (SIZE(RNORM,1) .NE. M) THEN ! The length of the array must match.
      IERR(:) = 28; RETURN; END IF
   RNORM(:) = 0.0_R8 ! Initialize output to zeros.
END IF
IF (PRESENT(CHAIN)) THEN
   CHAINL = CHAIN ! Turn chaining on, if necessarry.
   SEED(:) = 0 ! Initialize SEED in case it is needed.
ELSE
   CHAINL = .FALSE.
END IF

! Scale and center the data points and interpolation points.
CALL RESCALE(MINRAD, PTS_DIAM, PTS_SCALE)
IF (MINRAD < EPSL) THEN ! Check for degeneracies in points spacing.
   IERR(:) = 30; RETURN; END IF

! Query DGEQP3 for optimal work array size (LWORK).
LWORK = -1
CALL DGEQP3(D-1,D,LQ,D,IPIV,TAU,B,LWORK,IERR(1))
LWORK = INT(B(1)) ! Compute the optimal work array size.
ALLOCATE(WORK(LWORK), STAT=I) ! Allocate WORK to size LWORK.
IF (I .NE. 0) THEN ! Check for memory allocation errors.
   IERR(:) = 50; RETURN; END IF

! Initialize all error codes to "TBD" values.
IERR(:) = 40

! Outer loop over all interpolation points (in Q).
OUTER : DO MI = 1, M

   ! Check if this interpolation point was already found.
   IF (IERR(MI) .EQ. 0) CYCLE OUTER

   ! Initialize the projection and reset the residual.
   PROJ(:) = Q(:,MI)
   RNORML = 0.0_R8

   ! If there is no useable seed or if chaining is turned off, then make a new
   ! simplex.
   IF( (.NOT. CHAINL) .OR. SEED(1) .EQ. 0) THEN
      CALL MAKEFIRSTSIMP()
      IF(IERR(MI) .NE. 0) CYCLE OUTER
   ! Otherwise, use the seed.
   ELSE
      ! Copy the seed to the current simplex.
      SIMPS(:,MI) = SEED(:)
      ! Rebuild the linear system.
      DO J=1,D
         A(:,J) = PTS(:,SIMPS(J+1,MI)) - PTS(:,SIMPS(1,MI))
         B(J) = DDOT(D, A(:,J), 1, A(:,J), 1) / 2.0_R8 
      END DO
   END IF

   ! Inner loop searching for a simplex containing the point Q(:,MI).
   INNER : DO K = 1, IBUDGETL

      ! Check if the current simplex contains Q(:,MI).
      IF (PTINSIMP()) THEN
         ! If chaining is on, set the next seed to this simplex.
         IF (CHAINL) SEED(:) = SIMPS(:,MI)
         EXIT INNER
      END IF
      IF (IERR(MI) .NE. 0) CYCLE OUTER ! Check for an error flag.

      ! If chaining is on, save each good simplex as the next seed.
      IF (CHAINL) SEED(:) = SIMPS(:,MI)

      ! Swap out the least weighted vertex, but save its value in case it
      ! needs to be restored later.
      I = MINLOC(WEIGHTS(1:D+1,MI), DIM=1)
      ITMP = SIMPS(I,MI)
      SIMPS(I,MI) = SIMPS(D+1,MI)

      ! If the least weighted vertex (I) is not the first vertex, then just
      ! drop row I from the linear system.
      IF(I .NE. 1) THEN
         A(:,I-1) = A(:,D); B(I-1) = B(D)
      ! However, if I = 1, then both A and B must be reconstructed.
      ELSE
         DO J=1,D
            A(:,J) = PTS(:,SIMPS(J+1,MI)) - PTS(:,SIMPS(1,MI))
            B(J) = DDOT(D, A(:,J), 1, A(:,J), 1) / 2.0_R8 
         END DO
      END IF

      ! Compute the next simplex (do one flip).
      CALL MAKESIMPLEX()
      IF (IERR(MI) .NE. 0) CYCLE OUTER

      ! If no vertex was found, then this is an extrapolation point.
      IF (SIMPS(D+1,MI) .EQ. 0) THEN
         ! If extrapolation is not allowed (EXTRAP=0), do not proceed.
         IF (EXTRAPL < EPSL) THEN
            SIMPS(:,MI) = 0; WEIGHTS(:,MI) = 0  ! Zero all output values.
            ! Set the error flag and skip this point.
            IERR(MI) = 2; CYCLE OUTER
         END IF

         ! Otherwise, project the extrapolation point onto the convex hull.
         CALL PROJECT()
         IF (IERR(MI) .NE. 0) CYCLE OUTER

         ! Check the value of RNORML for over-extrapolation.
         IF (RNORML > EXTRAPL * PTS_DIAM) THEN
            SIMPS(:,MI) = 0; WEIGHTS(:,MI) = 0  ! Zero all output values.
            ! If present, record the unscaled RNORM output.
            IF (PRESENT(RNORM)) RNORM(MI) = RNORML*PTS_SCALE
            ! Set the error flag and skip this point.
            IERR(MI) = 2; CYCLE OUTER
         END IF

         ! Otherwise, restore the previous simplex and continue with the 
         ! projected value.
         SIMPS(D+1,MI) = ITMP
         A(:,D) = PTS(:,ITMP) - PTS(:,SIMPS(1,MI))
         B(D) = DDOT(D, A(:,D), 1, A(:,D), 1) / 2.0_R8
      END IF

   ! End of inner loop for finding each interpolation point.
   END DO INNER

   ! Check for budget violation conditions.
   IF (K > IBUDGETL) THEN 
      SIMPS(:,MI) = 0; WEIGHTS(:,MI) = 0  ! Zero all output values.
      ! Set the error flag and skip this point.
      IERR(MI) = 60; CYCLE OUTER
   END IF

   ! If the residual is nonzero, set the extrapolation flag.
   IF (RNORML > EPSL) IERR(MI) = 1

   ! If present, record the RNORM output.
   IF (PRESENT(RNORM)) RNORM(MI) = RNORML*PTS_SCALE

END DO OUTER  ! End of outer loop over all interpolation points.

! If INTERP_IN and INTERP_OUT are present, compute all values f(q).
IF (PRESENT(INTERP_IN)) THEN
   ! Loop over all interpolation points.
   DO MI = 1, M
      ! Check for errors.
      IF (IERR(MI) .LE. 1) THEN
         ! Compute the weighted sum of vertex response values.
         DO K = 1, D+1
            INTERP_OUT(:,MI) = INTERP_OUT(:,MI) &
             + INTERP_IN(:,SIMPS(K,MI)) * WEIGHTS(K,MI)
         END DO
      END IF
   END DO
END IF

! Free dynamic work arrays.
DEALLOCATE(WORK)
IF (ALLOCATED(IWORK_DWNNLS)) DEALLOCATE(IWORK_DWNNLS)
IF (ALLOCATED(WORK_DWNNLS))  DEALLOCATE(WORK_DWNNLS)
IF (ALLOCATED(W_DWNNLS))     DEALLOCATE(W_DWNNLS)
IF (ALLOCATED(X_DWNNLS))     DEALLOCATE(X_DWNNLS)

RETURN

CONTAINS    ! Internal subroutines and functions.

SUBROUTINE MAKEFIRSTSIMP()
! Iteratively construct the first simplex by choosing points that
! minimize the radius of the smallest circumball. Let (P_1, P_2, ..., P_K)
! denote the current list of vertices for the simplex. Let P* denote the
! candidate vertex to be added to the simplex. Let CENTER denote the
! circumcenter of the simplex.  Then
!
! X = CENTER - P_1
!
! is given by the minimum norm solution to the underdetermined linear system 
!
! A X = B, where
!
! A = [ P_2 - P_1, P_3 - P_1, ..., P_K - P_1, P* - P_1 ]^T and
! B = [ <A_{1.},A_{1.}>/2, <A_{2.},A_{2.}>/2, ..., <A_{K.},A_{K.}>/2 ]^T.
!
! Then the radius of the smallest circumsphere is CURRRAD = \| X \|,
! and the next vertex is given by P_{K+1} = argmin_{P*} CURRRAD, where P*
! ranges over points in PTS that are not already a vertex of the simplex.
!
! On output, this subroutine fully populates the matrix A^T and vector B,
! and fills SIMPS(:,MI) with the indices of a valid Delaunay simplex.

! Find the first point, i.e., the closest point to Q(:,MI).
SIMPS(:,MI) = 0
MINRAD = HUGE(0.0_R8)
DO I = 1, N
   ! Check the distance to Q(:,MI)
   CURRRAD = DNRM2(D, PTS(:,I) - PROJ(:), 1)
   IF (CURRRAD < MINRAD) THEN; MINRAD = CURRRAD; SIMPS(1,MI) = I; END IF
END DO
! Find the second point, i.e., the closest point to PTS(:,SIMPS(1,MI)).
MINRAD = HUGE(0.0_R8)
DO I = 1, N
   ! Skip repeated vertices.
   IF (I .EQ. SIMPS(1,MI)) CYCLE
   ! Check the diameter of the resulting circumsphere.
   CURRRAD = DNRM2(D, PTS(:,I)-PTS(:,SIMPS(1,MI)), 1)
   IF (CURRRAD < MINRAD) THEN; MINRAD = CURRRAD; SIMPS(2,MI) = I; END IF
END DO
! Set up the first row of the system A X = B.
A(:,1) = PTS(:,SIMPS(2,MI)) - PTS(:,SIMPS(1,MI))
B(1) = DDOT(D, A(:,1), 1, A(:,1), 1) / 2.0_R8
! Loop to collect the remaining D+1 vertices for the first simplex.
DO I = 2, D
   MINRAD = HUGE(0.0_R8) ! Re-initialize the radius for each iteration.
   ! Check each point P* in PTS.
   DO J = 1, N
      ! Check that this point is not already in the simplex.
      IF (ANY(SIMPS(:,MI) .EQ. J)) CYCLE
      ! Add P* to LS system, and compute the minimum norm solution.
      A(:,I) = PTS(:,J) - PTS(:,SIMPS(1,MI))
      B(I) = DDOT(D, A(:,I), 1, A(:,I), 1) / 2.0_R8
      LQ(:,1:I) = A(:,1:I)
      ! Compute A^T P = Q R.
      CALL DGEQP3(D, I, LQ, D, IPIV, TAU, WORK, LWORK, IERR(MI))
      IF(IERR(MI) < 0) THEN ! LAPACK illegal input error.
         IERR(MI) = 80; RETURN
      ELSE IF (ABS(LQ(I,I)) < EPSL) THEN ! A is rank-deficient.
         CYCLE ! If rank-deficient, skip this point.
      END IF
      ! Set CENTER = P^T B.
      FORALL (ITMP = 1:I) CENTER(ITMP) = B(IPIV(ITMP))
      ! Get the radius using R^T Q^T X = P^T B.
      CALL DTRSM('L', 'U', 'T', 'N', I, 1, 1.0_R8, LQ, D, CENTER, D)
      CENTER(I+1:D) = 0
      CALL DORMQR('L', 'N', D, 1, I, LQ, D, TAU, CENTER, D, WORK, LWORK, &
         IERR(MI))
      IF(IERR(MI) < 0) THEN ! LAPACK illegal input error.
         IERR(MI) = 83; RETURN
      END IF
      ! Calculate the radius and compare it to the current minimum.
      CURRRAD = DNRM2(D, CENTER, 1)
      IF (CURRRAD < MINRAD) THEN; MINRAD = CURRRAD; SIMPS(I+1,MI) = J; END IF
   END DO
   ! Check that a point was found. If not, then all the points must lie in a
   ! lower dimensional linear manifold (error case).
   IF (SIMPS(I+1,MI) .EQ. 0) THEN; IERR(MI) = 31; RETURN; END IF
   ! If all operations were successful, add the best P* to the LS system.
   A(:,I) = PTS(:,SIMPS(I+1,MI)) - PTS(:,SIMPS(1,MI))
   B(I) = DDOT(D, A(:,I), 1, A(:,I), 1) / 2.0_R8 
END DO
IERR(MI) = 0 ! Set error flag to 'success' for a normal return.
RETURN
END SUBROUTINE MAKEFIRSTSIMP

SUBROUTINE MAKESIMPLEX()
! Given a Delaunay facet F whose containing hyperplane does not contain
! Q(:,MI), complete the simplex by adding a point from PTS on the same `side'
! of F as Q(:,MI). Assume SIMPS(1:D,MI) contains the vertex indices of F
! (corresponding to data points P_1, P_2, ..., P_D in PTS), and assume the
! matrix A(1:D-1,:)^T and vector B(1:D-1) are filled appropriately (similarly
! as in MAKEFIRSTSIMP()). Then for any P* (not in the hyperplane containing
! F) in PTS, let CENTER denote the circumcenter of the simplex with vertices
! P_1, P_2, ..., P_D, P*. Then
!
! X = CENTER - P_1
!
! is given by the solution to the nonsingular linear system
!
! AX = B where
!
! A = [ P_2 - P_1, P_3 - P_1, ..., P_D - P_1, P* - P_1 ]^T and
! B = [ <A_{1.},A_{1.}>/2, <A_{2.},A_{2.}>/2, ..., <A_{D.},A_{D.}>/2 ]^T.
!
! Then CENTER = X + P_1 and RADIUS = \| X \|.  P_{D+1} will be given by the 
! candidate P* that satisfies both of the following:
!
! 1) Let PLANE denote the hyperplane containing F. Then P_{D+1} and Q(:,MI) 
! must be on the same side of PLANE.
!
! 2) The circumball about CENTER must not contain any points in PTS in its
! interior (Delaunay property).
! 
! The above are necessary and sufficient conditions for flipping the
! Delaunay simplex, given that F is indeed a Delaunay facet.
!
! On input, SIMPS(1:D,MI) should contain the vertex indices (column indices
! from PTS) of the facet F.  Upon output, SIMPS(:,MI) will contain the vertex
! indices of a Delaunay simplex closer to Q(:,MI).  Also, the matrix A^T and
! vector B will be updated accordingly. If SIMPS(D+1,MI)=0, then there were
! no points in PTS on the appropriate side of F, meaning that Q(:,MI) is an
! extrapolation point (not a convex combination of points in PTS).

! Compute the hyperplane PLANE.
CALL MAKEPLANE()
IF(IERR(MI) .NE. 0) RETURN ! Check for errors.
! Compute the sign for the side of PLANE containing Q(:,MI).
SIDE1 =  DDOT(D,PLANE(1:D),1,PROJ(:),1) - PLANE(D+1)
SIDE1 = SIGN(1.0_R8,SIDE1)
! Initialize the center, radius, and simplex.
SIMPS(D+1,MI) = 0
CENTER(:) = 0.0_R8
MINRAD = HUGE(0.0_R8)
! Loop through all points P* in PTS.
DO I = 1, N
   ! Check that P* is inside the current ball.
   IF (DNRM2(D, PTS(:,I) - CENTER(:), 1) > MINRAD) CYCLE ! If not, skip.
   ! Check that P* is on the appropriate halfspace.
   SIDE2 = DDOT(D,PLANE(1:D),1,PTS(:,I),1) - PLANE(D+1)
   IF (SIDE1*SIDE2 < EPSL .OR. ANY(SIMPS(:,MI) .EQ. I)) CYCLE ! If not, skip.
   ! Add P* to the linear system, and solve to get shifted CENTER.
   A(:,D) = PTS(:,I) - PTS(:,SIMPS(1,MI))
   B(D) = DDOT(D, A(:,D), 1, A(:,D), 1) / 2.0_R8
   LQ = A
   CENTER = B
   ! Compute A^T = LU
   CALL DGETRF(D, D, LQ, D, IPIV, IERR(MI))
   IF (IERR(MI) < 0) THEN ! LAPACK illegal input.
      IERR(MI) = 81; RETURN
   ELSE IF (IERR(MI) > 0) THEN ! Rank-deficiency detected.
      IERR(MI) = 61; RETURN
   END IF
   ! Use A^T = LU to solve A X = B, where X = CENTER - P_1.
   CALL DGETRS('T', D, 1, LQ, D, IPIV, CENTER, D, IERR(MI))
   IF (IERR(MI) < 0) THEN ! LAPACK illegal input.
      IERR(MI) = 82; RETURN
   END IF
   ! Update the new radius, center, and simplex.
   MINRAD = DNRM2(D, CENTER, 1)
   CENTER(:) = CENTER(:) + PTS(:,SIMPS(1,MI))
   SIMPS(D+1,MI) = I
END DO
IERR(MI) = 0 ! Reset the error flag to 'success' code.
! Check for extrapolation condition.
IF(SIMPS(D+1,MI) .EQ. 0) RETURN
! Add new point to the linear system.
A(:,D) = PTS(:,SIMPS(D+1,MI)) - PTS(:,SIMPS(1,MI))
B(D) = DDOT(D, A(:,D), 1, A(:,D), 1) / 2.0_R8
RETURN
END SUBROUTINE MAKESIMPLEX

SUBROUTINE MAKEPLANE()
! Construct a plane c^T x = \alpha containing the first D vertices indexed
! in SIMPS(:,MI). The plane is determined by its normal vector c and \alpha.
! Let (P_1, P_2, ..., P_D) be the vertices indexed in SIMPS(1:D,MI). A normal
! vector is any nonzero vector in ker(A^T), where the matrix
! 
! A^T = [ P_2 - P_1, P_3 - P_1, ..., P_D - P_1 ].
! 
! Since rank A^T = D-1, dim ker(A^T) = 1, and ker(A^T) can be found from a QR
! factorization of A^T:  A^T P = QR, where P permutes the columns of A^T. 
! Then the last column of Q is orthonormal to the range of A^T.
! 
! Upon output, PLANE(1:D) contains the normal vector c and PLANE(D+1)
! contains \alpha defining the plane.

IF (D > 1) THEN ! Check that D-1 > 0, otherwise the plane is trivial.
   ! Compute the QR factorization.
   IPIV=0
   LQ = A
   CALL DGEQP3(D, D-1, LQ, D, IPIV, TAU, WORK, LWORK, IERR(MI))
   IF(IERR(MI) < 0) THEN ! LAPACK illegal input error.
      IERR(MI) = 80; RETURN
   END IF
   ! The nullspace is given by the last column of Q.
   PLANE(1:D-1) = 0.0_R8
   PLANE(D) = 1.0_R8
   CALL DORMQR('L', 'N', D, 1, D-1, LQ, D, TAU, PLANE, D, WORK, &
    LWORK, IERR(MI))
   IF(IERR(MI) < 0) THEN ! LAPACK illegal input error.
      IERR(MI) = 83; RETURN
   END IF
   ! Calculate the constant \alpha defining the plane.
   PLANE(D+1) = DDOT(D,PLANE(1:D),1,PTS(:,SIMPS(1,MI)),1)
ELSE ! Special case where D=1.
   PLANE(1) = 1.0_R8
   PLANE(2) = PTS(1,SIMPS(1,MI))
END IF
RETURN
END SUBROUTINE MAKEPLANE

FUNCTION PTINSIMP() RESULT(TF)
! Determine if any interpolation points are in the current simplex, whose
! vertices (P_1, P_2, ..., P_{D+1}) are indexed by SIMPS(:,MI). These
! vertices determine a positive cone with generators V_I = P_{I+1} - P_1,
! I = 1, ..., D. For each interpolation point Q* in Q, Q* - P_1 can be
! expressed as a unique linear combination of the V_I.  If all these linear
! weights are nonnegative and sum to less than or equal to 1.0, then Q* is
! in the simplex with vertices {P_I}_{I=1}^{D+1}.
!
! If any interpolation points in Q are contained in the simplex whose
! vertices are indexed by SIMPS(:,MI), then those points are marked as solved
! and the values of SIMPS and WEIGHTS are updated appropriately. On output,
! WEIGHTS(:,MI) contains the affine weights for producing Q(:,MI) as an
! affine combination of the points in PTS indexed by SIMPS(:,MI). If these
! weights are nonnegative, then PTINSIMP() returns TRUE.

! Initialize the return value and local variables.
LOGICAL :: TF ! True/False value.
TF = .FALSE.

! Compute the LU factorization of the matrix A whose rows are P_{I+1} - P_1.
LQ = A
CALL DGETRF(D, D, LQ, D, IPIV, IERR(MI))
IF (IERR(MI) < 0) THEN ! LAPACK illegal input.
   IERR(MI) = 81; RETURN
ELSE IF (IERR(MI) > 0) THEN ! Rank-deficiency detected.
   IERR(MI) = 61; RETURN
END IF
! Solve A^T w = WORK to get the affine weights for Q(:,MI) or its projection.
WORK(1:D) = PROJ(:) - PTS(:,SIMPS(1,MI))
CALL DGETRS('N', D, 1, LQ, D, IPIV, WORK(1:D), D, IERR(MI))
IF (IERR(MI) < 0) THEN ! LAPACK illegal input.
   IERR(MI) = 82; RETURN
END IF
WEIGHTS(2:D+1,MI) = WORK(1:D)
WEIGHTS(1,MI) = 1.0_R8 - SUM(WEIGHTS(2:D+1,MI))
! Check if the weights for Q(:,MI) are nonnegative.
IF (ALL(WEIGHTS(:,MI) .GE. -EPSL)) TF = .TRUE.

! Compute the affine weights for the rest of the interpolation points.
DO I = MI+1, M
   ! Check that no solution has already been found.
   IF (IERR(I) .NE. 40) CYCLE
   ! Solve A^T w = WORK to get the affine weights for Q(:,I).
   WORK(2:D+1) = Q(:,I) - PTS(:,SIMPS(1,MI))
   CALL DGETRS('N', D, 1, LQ, D, IPIV, WORK(2:D+1), D, ITMP)
   IF (ITMP < 0) CYCLE ! Illegal input error that should never occurr.
   ! Check if the weights define a convex combination.
   WORK(1) = 1.0_R8 - SUM(WORK(2:D+1))
   IF (ALL(WORK(1:D+1) .GE. -EPSL)) THEN
      ! Copy the simplex indices and weights then flag as complete.
      SIMPS(:,I) = SIMPS(:,MI)
      WEIGHTS(:,I) = WORK(1:D+1)
      IERR(I) = 0
   END IF
END DO
RETURN
END FUNCTION PTINSIMP

SUBROUTINE PROJECT()
! Project a point outside the convex hull of the point set onto the convex hull
! by solving an inequality constrained least squares problem. The solution to
! the least squares problem gives the projection as a convex combination of the
! data points. The projection can then be computed by performing a matrix
! vector multiplication.

! Allocate work arrays.
IF (.NOT. ALLOCATED(IWORK_DWNNLS)) THEN
   ALLOCATE(IWORK_DWNNLS(D+1+N), STAT=IERR(MI))
   IF(IERR(MI) .NE. 0) THEN; IERR(MI) = 70; RETURN; END IF
END IF
IF (.NOT. ALLOCATED(WORK_DWNNLS)) THEN
   ALLOCATE(WORK_DWNNLS(D+1+N*5), STAT=IERR(MI))
   IF(IERR(MI) .NE. 0) THEN; IERR(MI) = 70; RETURN; END IF
END IF
IF (.NOT. ALLOCATED(W_DWNNLS)) THEN
   ALLOCATE(W_DWNNLS(D+1,N+1), STAT=IERR(MI))
   IF(IERR(MI) .NE. 0) THEN; IERR(MI) = 70; RETURN; END IF
END IF
IF (.NOT. ALLOCATED(X_DWNNLS)) THEN
   ALLOCATE(X_DWNNLS(N), STAT=IERR(MI))
   IF(IERR(MI) .NE. 0) THEN; IERR(MI) = 70; RETURN; END IF
END IF

! Initialize work array and settings values.
PRGOPT_DWNNLS(1) = 1.0_R8
IWORK_DWNNLS(1) = D+1+5*N
IWORK_DWNNLS(2) = D+1+N
W_DWNNLS(1, :) = 1.0_R8          ! Set convexity (equality) constraint.
W_DWNNLS(2:D+1,1:N) = PTS(:,:)   ! Copy data points.
W_DWNNLS(2:D+1,N+1) = PROJ(:)    ! Copy extrapolation point.
! Compute the solution to the inequality constrained least squares problem to
! get the projection coefficients.
CALL DWNNLS(W_DWNNLS, D+1, 1, D, N, 0, PRGOPT_DWNNLS, X_DWNNLS, RNORML, &
   IERR(MI), IWORK_DWNNLS, WORK_DWNNLS)
IF (IERR(MI) .EQ. 1) THEN ! Failure to converge.
   IERR(MI) = 71; RETURN
ELSE IF (IERR(MI) .EQ. 2) THEN ! Illegal input detected.
   IERR(MI) = 72; RETURN
END IF
! Compute the actual projection via matrix vector multiplication.
CALL DGEMV('N', D, N, 1.0_R8, PTS, D, X_DWNNLS, 1, 0.0_R8, PROJ, 1)
RETURN
END SUBROUTINE PROJECT

SUBROUTINE RESCALE(MINDIST, DIAMETER, SCALE)
! Rescale and transform data to be centered at the origin with unit
! radius. This subroutine has O(n^2) complexity.
!
! On output, PTS and Q have been rescaled and shifted. All the data
! points in PTS are centered with unit radius, and the points in Q
! have been shifted and scaled in relation to PTS.
!
! MINDIST is a real number containing the (scaled) minimum distance
!    between any two data points in PTS.
!
! DIAMETER is a real number containing the (scaled) diameter of the
!    data set PTS.
!
! SCALE contains the real factor used to transform the data and
!    interpolation points: scaled value = (original value -
!    barycenter of data points)/SCALE.

! Output arguments.
REAL(KIND=R8), INTENT(OUT) :: MINDIST, DIAMETER, SCALE

! Local variables.
REAL(KIND=R8) :: PTS_CENTER(D) ! The center of the data points PTS.
REAL(KIND=R8) :: DISTANCE ! The current distance.

! Initialize local values.
MINDIST = HUGE(0.0_R8)
DIAMETER = 0.0_R8
SCALE = 0.0_R8

! Compute barycenter of all data points.
PTS_CENTER(:) = SUM(PTS(:,:), DIM=2)/REAL(N, KIND=R8)
! Compute the minimum and maximum distances.
DO I = 1, N ! Cycle through all pairs of points.
   DO J = I + 1, N
      DISTANCE = DNRM2(D, PTS(:,I) - PTS(:,J), 1) ! Compute the distance.
      IF (DISTANCE > DIAMETER) THEN ! Compare to the current diameter.
         DIAMETER = DISTANCE
      END IF
      IF (DISTANCE < MINDIST) THEN ! Compare to the current minimum distance.
         MINDIST = DISTANCE
      END IF
   END DO
END DO
! Center the points.
FORALL (I = 1:N) PTS(:,I) = PTS(:,I) - PTS_CENTER(:)
! Compute the scale factor (for unit radius).
DO I = 1, N ! Cycle through all points again.
   DISTANCE = DNRM2(D, PTS(:,I), 1) ! Compute the distance from the center.
   IF (DISTANCE > SCALE) THEN ! Compare to the current radius.
      SCALE = DISTANCE
   END IF
END DO
! Scale the points to unit radius.
PTS = PTS / SCALE
! Also transform Q similarly.
FORALL (I = 1:M) Q(:,I) = (Q(:,I) - PTS_CENTER(:)) / SCALE
! Also scale the minimum distance and diameter.
MINDIST = MINDIST / SCALE
DIAMETER = DIAMETER / SCALE
RETURN
END SUBROUTINE RESCALE

END SUBROUTINE DELAUNAYSPARSES


SUBROUTINE DELAUNAYSPARSEP( D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, & 
  INTERP_IN, INTERP_OUT, EPS, EXTRAP, RNORM, IBUDGET, CHAIN,       &
  PMODE, CHUNKSIZE )
! This is a parallel implementation of an algorithm for efficiently performing
! interpolation in R^D via the Delaunay triangulation. The algorithm is fully
! described and analyzed in
!
! T. H. Chang, L. T. Watson,  T. C.H. Lux, B. Li, L. Xu, A. R. Butt, K. W.
! Cameron, and Y. Hong. 2018. A polynomial time algorithm for multivariate
! interpolation in arbitrary dimension via the Delaunay triangulation. In
! Proceedings of the ACMSE 2018 Conference (ACMSE '18). ACM, New York, NY,
! USA. Article 12, 8 pages.
!
!
! On input:
!
! D is the dimension of the space for PTS and Q.
!
! N is the number of data points in PTS.
!
! PTS(1:D,1:N) is a real valued matrix with N columns, each containing the
!    coordinates of a single data point in R^D.
!
! M is the number of interpolation points in Q.
!
! Q(1:D,1:M) is a real valued matrix with M columns, each containing the
!    coordinates of a single interpolation point in R^D.
!
!
! On output:
!
! PTS and Q have been rescaled and shifted. All the data points in PTS
!    are now contained in the unit hyperball in R^D, and the points in Q
!    have been shifted and scaled accordingly in relation to PTS.
!
! SIMPS(1:D+1,1:M) contains the D+1 integer indices (corresponding to columns
!    in PTS) for the D+1 vertices of the Delaunay simplex containing each
!    interpolation point in Q.
!
! WEIGHTS(1:D+1,1:M) contains the D+1 real valued weights for expressing each
!    point in Q as a convex combination of the D+1 corresponding vertices
!    in SIMPS.
!
! IERR(1:M) contains integer valued error flags associated with the
!    computation of each of the M interpolation points in Q. The error
!    codes are:
!
! 00 : Succesful interpolation.
! 01 : Succesful extrapolation (up to the allowed extrapolation distance).
! 02 : This point was outside the allowed extrapolation distance; the
!      corresponding entries in SIMPS and WEIGHTS contain zero values.
!
! 10 : The dimension D must be positive.
! 11 : Too few data points to construct a triangulation (i.e., N < D+1).
! 12 : No interpolation points given (i.e., M < 1).
! 13 : The first dimension of PTS does not agree with the dimension D.
! 14 : The second dimension of PTS does not agree with the number of points N.
! 15 : The first dimension of Q does not agree with the dimension D.
! 16 : The second dimension of Q does not agree with the number of
!      interpolation points M.
! 17 : The first dimension of the output array SIMPS does not match the number
!      of vertices needed for a D-simplex (D+1).
! 18 : The second dimension of the output array SIMPS does not match the
!      number of interpolation points M.
! 19 : The first dimension of the output array WEIGHTS does not match the
!      number of vertices for a a D-simplex (D+1).
! 20 : The second dimension of the output array WEIGHTS does not match the
!      number of interpolation points M.
! 21 : The size of the error array IERR does not match the number of
!      interpolation points M.
! 22 : INTERP_IN cannot be present without INTERP_OUT or vice versa.
! 23 : The first dimension of INTERP_IN does not match the first
!      dimension of INTERP_OUT.
! 24 : The second dimension of INTERP_IN does not match the number of
!      data points PTS.
! 25 : The second dimension of INTERP_OUT does not match the number of
!      interpolation points M.
! 26 : The budget supplied in IBUDGET does not contain a positive
!      integer.
! 27 : The extrapolation distance supplied in EXTRAP cannot be negative.
! 28 : The size of the RNORM output array does not match the number of
!      interpolation points M.
!
! 30 : Two or more points in the data set PTS are too close together with
!      respect to the working precision (EPS), which would result in a
!      numerically degenerate simplex.
! 31 : All the data points in PTS lie in some lower dimensional linear
!      manifold (up to the working precision), and no valid triangulation
!      exists.
! 40 : An error caused DELAUNAYSPARSEP to terminate before this value could
!      be computed. Note: The corresponding entries in SIMPS and WEIGHTS may
!      contain garbage values.
!
! 50 : A memory allocation error occurred while allocating the work array
!      WORK.
!
! 60 : The budget was exceeded before the algorithm converged on this
!      value. If the dimension is high, try increasing IBUDGET. This
!      error can also be caused by a working precision EPS that is too
!      small for the conditioning of the problem.
!
! 61 : A value that was judged appropriate later caused LAPACK to encounter a
!      singularity. Try increasing the value of EPS.
!
! 70 : Allocation error for the extrapolation work arrays.
! 71 : The SLATEC subroutine DWNNLS failed to converge during the projection
!      of an extrapolation point onto the convex hull.
! 72 : The SLATEC subroutine DWNNLS has reported a usage error.
!
!      The errors 72, 80--83 should never occur, and likely indicate a
!      compiler bug or hardware failure.
! 80 : The LAPACK subroutine DGEQP3 has reported an illegal value.
! 81 : The LAPACK subroutine DGETRF has reported an illegal value.
! 82 : The LAPACK subroutine DGETRS has reported an illegal value.
! 83 : The LAPACK subroutine DORMQR has reported an illegal value.
!
! 90 : The value of PMODE is not valid.
! 91 : The value of CHUNKSIZE is not valid.
!
! Optional arguments:
!
! INTERP_IN(1:IR,1:N) contains real valued response vectors for each of
!    the data points in PTS on input. The first dimension of INTERP_IN is
!    inferred to be the dimension of these response vectors, and the
!    second dimension must match N. If present, the response values will
!    be computed for each interpolation point in Q, and stored in INTERP_OUT,
!    which therefore must also be present. If both INTERP_IN and INTERP_OUT
!    are omitted, only the containing simplices and convex combination
!    weights are returned.
! 
! INTERP_OUT(1:IR,1:M) contains real valued response vectors for each
!    interpolation point in Q on output. The first dimension of INTERP_OU
!    must match the first dimension of INTERP_IN, and the second dimension
!    must match M. If present, the response values at each interpolation
!    point are computed as a convex combination of the response values
!    (supplied in INTERP_IN) at the vertices of a Delaunay simplex containing
!    that interpolation point.  Therefore, if INTERP_OUT is present, then
!    INTERP_IN must also be present.  If both are omitted, only the
!    simplices and convex combination weights are returned.
! 
! EPS contains the real working precision for the problem on input. By
!    default, EPS is assigned \sqrt{\mu} where \mu denotes the unit roundoff
!    for the machine. In general, any values that differ by less than EPS
!    are judged as equal, and any weights that are greater than -EPS are
!    judged as nonnegative.  EPS cannot take a value less than the default
!    value of \sqrt{\mu}. If any value less than \sqrt{\mu} is supplied,
!    the default value will be used instead automatically. 
! 
! EXTRAP contains the real maximum extrapolation distance (relative to the
!    diameter of PTS) on input. Interpolation at a point outside the convex
!    hull of PTS is done by projecting that point onto the convex hull, and
!    then doing normal Delaunay interpolation at that projection.
!    Interpolation at any point in Q that is more than EXTRAP * DIAMETER(PTS)
!    units outside the convex hull of PTS will not be done and an error code
!    of 2 will be returned. Note that computing the projection can be
!    expensive. Setting EXTRAP=0 will cause all extrapolation points to be
!    ignored without ever computing a projection. By default, EXTRAP=0.1
!    (extrapolate by up to 10% of the diameter of PTS). 
! 
! RNORM(1:M) contains the real unscaled projection (2-norm) distances from
!    any projection computations on output. If not present, these distances
!    are still computed for each extrapolation point, but are never returned.
! 
! IBUDGET on input contains the integer budget for performing flips while
!    iterating toward the simplex containing each interpolation point in Q.
!    This prevents DELAUNAYSPARSEP from falling into an infinite loop when
!    an inappropriate value of EPS is given with respect to the problem
!    conditioning.  By default, IBUDGET=50000. However, for extremely
!    high-dimensional problems and pathological inputs, the default value
!    may be insufficient. 
!
! CHAIN is a logical input argument that determines whether a new first
!    simplex should be constructed for each interpolation point
!    (CHAIN=.FALSE.), or whether the simplex walks should be "daisy-chained."
!    By default, CHAIN=.FALSE. Setting CHAIN=.TRUE. is generally not
!    recommended, unless the size of the triangulation is relatively small
!    or the interpolation points are known to be tightly clustered.
!
! CHUNKSIZE contains the integer valued size of dynamically scheduled chunks
!    during Level 1 parallelism. By default, the CHUNKSIZE is 1, for fully
!    dynamic scheduling. Changing CHUNKSIZE is only recommended if the
!    dimension is low or contiguous points in Q are known to be relatively
!    nearby.
!
! PMODE is an integer specifying the level of parallelism to be exploited.
!    If PMODE = 1, then parallelism is exploited at the level of the loop
!    over all interpolation points (Level 1 parallelism).
!    If PMODE = 2, then parallelism is exploited at the level of the loops
!    over data points when constructing/flipping simplices (Level 2
!    parallelism).
!    If PMODE = 3, then parallelism is exploited at both levels. Note: this
!    implies that the total number of threads active at any time could be up
!    to OMP_NUM_THREADS^2.
!    By default, PMODE is set to 1 if there is more than 1 interpolation
!    point and 2 otherwise.
!
! Subroutines and functions directly referenced from BLAS are
!      DDOT, DGEMV, DNRM2, DTRSM,
! and from LAPACK are
!      DGEQP3, DGETRF, DGETRS, DORMQR.
! The SLATEC subroutine DWNNLS is directly referenced. DWNNLS and all its
! SLATEC dependencies have been slightly edited to comply with the Fortran
! 2008 standard, with all print statements and references to stderr being
! commented out. For a reference to DWNNLS, see ACM TOMS Algorithm 587
! (Hanson and Haskell).  The module REAL_PRECISION from HOMPACK90 (ACM TOMS
! Algorithm 777) is used for the real data type. The REAL_PRECISION module,
! DELAUNAYSPARSEP, and DWNNLS and its dependencies comply with the Fortran
! 2008 standard.  
! 
! Primary Author: Tyler H. Chang
! Last Update: September, 2018
! 
USE REAL_PRECISION, ONLY : R8
IMPLICIT NONE

! Input arguments.
INTEGER, INTENT(IN) :: D, N
REAL(KIND=R8), INTENT(INOUT) :: PTS(:,:) ! Rescaled on output.
INTEGER, INTENT(IN) :: M
REAL(KIND=R8), INTENT(INOUT) :: Q(:,:) ! Rescaled on output.
! Output arguments.
INTEGER, INTENT(OUT) :: SIMPS(:,:)
REAL(KIND=R8), INTENT(OUT) :: WEIGHTS(:,:)
INTEGER, INTENT(OUT) :: IERR(:)
! Optional arguments.
REAL(KIND=R8), INTENT(IN), OPTIONAL:: INTERP_IN(:,:)
REAL(KIND=R8), INTENT(OUT), OPTIONAL :: INTERP_OUT(:,:)
REAL(KIND=R8), INTENT(IN), OPTIONAL:: EPS, EXTRAP 
REAL(KIND=R8), INTENT(OUT), OPTIONAL :: RNORM(:)
INTEGER, INTENT(IN), OPTIONAL :: IBUDGET, CHUNKSIZE, PMODE
LOGICAL, INTENT(IN), OPTIONAL :: CHAIN 

! Local copies of optional input arguments.
REAL(KIND=R8) :: EPSL, EXTRAPL
INTEGER :: IBUDGETL, CHUNKSIZEL
LOGICAL :: CHAINL, PLVL1, PLVL2

! Local variables.
LOGICAL :: PTINSIMP ! Tells if Q(:,MI) is in SIMPS(:,MI).
INTEGER :: I, J, K ! Loop iteration variables.
INTEGER :: IERR_PRIV ! Private copy of the error flag.
INTEGER :: ITMP ! Temporary variable for swapping, looping, etc.
INTEGER :: LAST ! Last index of the current chunk.
INTEGER :: LWORK ! Size of the work array.
INTEGER :: MI ! Index of current interpolation point.
INTEGER :: VERTEX_PRIV ! Private copy of next vertex to add.
REAL(KIND=R8) :: CURRRAD ! Radius of the current circumsphere.
REAL(KIND=R8) :: MINRAD ! Minimum circumsphere radius observed.
REAL(KIND=R8) :: MINRAD_PRIV ! Private copy of MINRAD.
REAL(KIND=R8) :: PTS_SCALE ! Data scaling factor.
REAL(KIND=R8) :: PTS_DIAM ! Scaled diameter of data set.
REAL(KIND=R8) :: RNORML ! Euclidean norm of the projection residual.
REAL(KIND=R8) :: SIDE1, SIDE2 ! Signs (+/-1) denoting sides of a facet.

! Local arrays, requiring O(d^2) additional memory.
INTEGER :: IPIV(D) ! Pivot indices.
INTEGER :: SEED(D+1) ! Copy of the SEED simplex. Only used if CHAIN = .TRUE.
REAL(KIND=R8) :: A(D,D) ! The LHS of a linear system.
REAL(KIND=R8) :: B(D) ! The RHS of a linear system. 
REAL(KIND=R8) :: CENTER(D) ! The circumcenter of a simplex. 
REAL(KIND=R8) :: LQ(D,D) ! The LU or QR factorization of A.
REAL(KIND=R8) :: PLANE(D+1) ! The hyperplane containing a facet.
REAL(KIND=R8) :: PRGOPT_DWNNLS(1) ! Options array for DWNNLS.
REAL(KIND=R8) :: PROJ(D) ! The projection of the current iterate.
REAL(KIND=R8) :: TAU(D) ! Householder reflector constants.

! Extrapolation work arrays are only allocated if DWNNLS is called.
INTEGER, ALLOCATABLE :: IWORK_DWNNLS(:) ! Only for DWNNLS.
REAL(KIND=R8), ALLOCATABLE :: W_DWNNLS(:,:) ! Only for DWNNLS.
REAL(KIND=R8), ALLOCATABLE :: WORK(:) ! Allocated with size LWORK.
REAL(KIND=R8), ALLOCATABLE :: WORK_DWNNLS(:) ! Only for DWNNLS.
REAL(KIND=R8), ALLOCATABLE :: X_DWNNLS(:) ! Only for DWNNLS.

! External functions and subroutines.
REAL(KIND=R8), EXTERNAL :: DDOT  ! Inner product (BLAS).
REAL(KIND=R8), EXTERNAL :: DNRM2 ! Euclidean norm (BLAS).
EXTERNAL :: DGEMV  ! General matrix vector multiply (BLAS)
EXTERNAL :: DGEQP3 ! Perform a QR factorization with column pivoting (LAPACK).
EXTERNAL :: DGETRF ! Perform a LU factorization with partial pivoting (LAPACK).
EXTERNAL :: DGETRS ! Use the output of DGETRF to solve a linear system (LAPACK).
EXTERNAL :: DORMQR ! Apply householder reflectors to a matrix (LAPACK).
EXTERNAL :: DTRSM  ! Perform a triangular solve (BLAS).
EXTERNAL :: DWNNLS ! Solve an inequality constrained least squares problem
                   ! (SLATEC).

! Check for input size and dimension errors.
IF (D < 1) THEN ! The dimension must satisfy D > 0.
   IERR(:) = 10; RETURN; END IF
IF (N < D+1) THEN ! Must have at least D+1 data points.
   IERR(:) = 11; RETURN; END IF
IF (M < 1) THEN ! Must have at least one interpolation point.
   IERR(:) = 12; RETURN; END IF
IF (SIZE(PTS,1) .NE. D) THEN ! Dimension of PTS array should match.
   IERR(:) = 13; RETURN; END IF
IF (SIZE(PTS,2) .NE. N) THEN ! Number of data points should match.
   IERR(:) = 14; RETURN; END IF
IF (SIZE(Q,1) .NE. D) THEN ! Dimension of Q should match.
   IERR(:) = 15; RETURN; END IF
IF (SIZE(Q,2) .NE. M) THEN ! Number of interpolation points should match.
   IERR(:) = 16; RETURN; END IF
IF (SIZE(SIMPS,1) .NE. D+1) THEN ! Need space for D+1 vertices per simplex.
   IERR(:) = 17; RETURN; END IF
IF (SIZE(SIMPS,2) .NE. M) THEN  ! There will be M output simplices.
   IERR(:) = 18; RETURN; END IF
IF (SIZE(WEIGHTS,1) .NE. D+1) THEN ! There will be D+1 weights per simplex.
   IERR(:) = 19; RETURN; END IF
IF (SIZE(WEIGHTS,2) .NE. M) THEN ! One vector of weights per simplex.
   IERR(:) = 20; RETURN; END IF
IF (SIZE(IERR) .NE. M) THEN ! An error flag for each interpolation point.
   IERR(:) = 21; RETURN; END IF

! Check for optional arguments.
IF (PRESENT(INTERP_IN) .NEQV. PRESENT(INTERP_OUT)) THEN
   IERR(:) = 22; RETURN; END IF
IF (PRESENT(INTERP_IN)) THEN ! Sizes must agree.
   IF (SIZE(INTERP_IN,1) .NE. SIZE(INTERP_OUT,1)) THEN
      IERR(:) = 23 ; RETURN; END IF
   IF(SIZE(INTERP_IN,2) .NE. N) THEN
      IERR(:) = 24; RETURN; END IF
   IF (SIZE(INTERP_OUT,2) .NE. M) THEN
      IERR(:) = 25; RETURN; END IF
   INTERP_OUT(:,:) = 0.0_R8 ! Initialize output to zeros.
END IF
EPSL = SQRT(EPSILON(0.0_R8)) ! Get the machine unit roundoff constant.
IF (PRESENT(EPS)) THEN
   IF (EPSL < EPS) THEN ! If the given precision is too small, ignore it.
      EPSL = EPS
   END IF
END IF
IF (PRESENT(IBUDGET)) THEN
   IBUDGETL = IBUDGET ! Use the given budget if present.
   IF (IBUDGETL < 1) THEN
      IERR(:) = 26; RETURN; END IF
ELSE
   IBUDGETL = 50000 ! Default value for budget.
END IF
IF (PRESENT(EXTRAP)) THEN
   EXTRAPL = EXTRAP 
   IF (EXTRAPL < 0) THEN ! Check that the extrapolation distance is legal.
      IERR(:) = 27; RETURN; END IF
ELSE
   EXTRAPL = 0.1_R8 ! Default extrapolation distance (for normalized points).
END IF
IF (PRESENT(RNORM)) THEN
   IF (SIZE(RNORM,1) .NE. M) THEN ! The length of the array must match.
      IERR(:) = 28; RETURN; END IF
   RNORM(:) = 0.0_R8 ! Initialize output to zeros.
END IF
IF (PRESENT(CHAIN)) THEN
   CHAINL = CHAIN ! Turn chaining on, if necessarry.
   SEED(:) = 0 ! Initialize SEED in case it is needed.
ELSE
   CHAINL = .FALSE.
END IF
! Set the PMODE.
PLVL1 = .FALSE.
PLVL2 = .FALSE.
IF (PRESENT(PMODE)) THEN ! Check PMODE for legal values.
   IF (PMODE .EQ. 1) THEN
      PLVL1 = .TRUE.
   ELSE IF (PMODE .EQ. 2) THEN
      PLVL2 = .TRUE.
   ELSE IF (PMODE .EQ. 3) THEN
      PLVL1 = .TRUE.; PLVL2 = .TRUE.
   ELSE
      IERR(:) = 90; RETURN
   END IF
ELSE ! The default setting for PMODE is level 1 parallelism if M > 1.
   IF (M > 1) THEN
      PLVL1 = .TRUE.
   ELSE
      PLVL2 = .TRUE.
   END IF
END IF
IF (PRESENT(CHUNKSIZE)) THEN
   IF (CHUNKSIZE .LE. 0) THEN ! Illegal chunk size.
      IERR(:) = 91; RETURN; END IF
   CHUNKSIZEL = CHUNKSIZE
ELSE
   CHUNKSIZEL = 1 ! Default to fully dynamic scheduling.
END IF

! Scale and center the data points and interpolation points.
CALL RESCALE(MINRAD, PTS_DIAM, PTS_SCALE)
IF (MINRAD < EPSL) THEN ! Check for degeneracies in points spacing.
   IERR(:) = 30; RETURN; END IF

! Query DGEQP3 for optimal work array size (LWORK).
LWORK = -1
CALL DGEQP3(D-1,D,LQ,D,IPIV,TAU,B,LWORK,IERR(1))
LWORK = INT(B(1)) ! Compute the optimal work array size.
ALLOCATE(WORK(LWORK), STAT=I) ! Allocate WORK to size LWORK.
IF (I .NE. 0) THEN ! Check for memory allocation errors.
   IERR(:) = 50; RETURN; END IF

! Initialize PRGOPT_DWNNLS in case of extrapolation.
PRGOPT_DWNNLS(1) = 1.0_R8

! Initialize all error codes to "TBD" values.
IERR(:) = 40

! Level 1 parallel region (over all interpolation points in Q).
!$OMP PARALLEL FIRSTPRIVATE(SEED), &
!$OMP& PRIVATE(I, J, K, ITMP, CURRRAD, MINRAD, RNORML, SIDE1, SIDE2), &
!$OMP& PRIVATE(IERR_PRIV, LAST, VERTEX_PRIV, MINRAD_PRIV, PTINSIMP), &
!$OMP& PRIVATE(IPIV, A, B, CENTER, LQ, PLANE, PROJ, TAU), &
!$OMP& PRIVATE(IWORK_DWNNLS, W_DWNNLS, WORK, WORK_DWNNLS, X_DWNNLS), &
!$OMP& IF(PLVL1)
!$OMP DO SCHEDULE(DYNAMIC, CHUNKSIZEL)
OUTER : DO MI = 1, M
   ! Check if this interpolation point was already found.
   IF (IERR(MI) .EQ. 0) CYCLE OUTER

   ! Initialize the projection and reset the residual.
   PROJ(:) = Q(:,MI)
   RNORML = 0.0_R8

   ! If there is no useable seed or if chaining is turned off, then make a new
   ! simplex.
   IF( (.NOT. CHAINL) .OR. SEED(1) .EQ. 0) THEN
!      CALL MAKEFIRSTSIMP(); IF(IERR(MI) .NE. 0) CYCLE OUTER
!
! Due to OpenMP's handling of variable scope, the parallel implementation of
! the subroutine MAKEFIRSTSIMP() has been in-lined here.
!
! SUBROUTINE MAKEFIRSTSIMP()
!
! Iteratively construct the first simplex by choosing points that
! minimize the radius of the smallest circumball. Let (P_1, P_2, ..., P_K)
! denote the current list of vertices for the simplex. Let P* denote the
! candidate vertex to be added to the simplex. Let CENTER denote the
! circumcenter of the simplex.  Then
!
! X = CENTER - P_1
!
! is given by the minimum norm solution to the underdetermined linear system 
!
! A X = B, where
!
! A = [ P_2 - P_1, P_3 - P_1, ..., P_K - P_1, P* - P_1 ]^T and
! B = [ <A_{1.},A_{1.}>/2, <A_{2.},A_{2.}>/2, ..., <A_{K.},A_{K.}>/2 ]^T.
!
! Then the radius of the smallest circumsphere is CURRRAD = \| X \|,
! and the next vertex is given by P_{K+1} = argmin_{P*} CURRRAD, where P*
! ranges over points in PTS that are not already a vertex of the simplex.
!
! On output, this subroutine fully populates the matrix A^T and vector B,
! and fills SIMPS(:,MI) with the indices of a valid Delaunay simplex.

! Initialize simplex and shared variables.
SIMPS(:,MI) = 0
MINRAD_PRIV = HUGE(0.0_R8)
MINRAD = HUGE(0.0_R8)

! Below is a Level 2 parallel region over N data points to find the
! first and second vertices SIMPS(1,MI) and SIMPS(2,MI).
!$OMP PARALLEL FIRSTPRIVATE(MINRAD_PRIV), PRIVATE(CURRRAD, VERTEX_PRIV), &
!$OMP& IF(PLVL2)
! Find the first point, i.e., the closest point to Q(:,MI).
!$OMP DO SCHEDULE(STATIC)
DO I = 1, N
   ! Check the distance to Q(:,MI)
   CURRRAD = DNRM2(D, PTS(:,I) - PROJ(:), 1)
   IF (CURRRAD < MINRAD_PRIV) THEN
      MINRAD_PRIV = CURRRAD; VERTEX_PRIV = I;
   END IF
END DO
!$OMP CRITICAL
IF (MINRAD_PRIV < MINRAD) THEN
   MINRAD = MINRAD_PRIV; SIMPS(1,MI) = VERTEX_PRIV;
END IF
! Find the second point, i.e., the closest point to PTS(:,SIMPS(1,MI)).
MINRAD_PRIV = HUGE(0.0_R8)
!$OMP END CRITICAL
!$OMP BARRIER
!$OMP SINGLE
MINRAD = HUGE(0.0_R8)
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
DO I = 1, N
   ! Skip repeated vertices.
   IF (I .EQ. SIMPS(1,MI)) CYCLE
   ! Check the diameter of the resulting circumsphere.
   CURRRAD = DNRM2(D, PTS(:,I)-PTS(:,SIMPS(1,MI)), 1)
   IF (CURRRAD < MINRAD_PRIV) THEN
      MINRAD_PRIV = CURRRAD; VERTEX_PRIV = I
   END IF
END DO
!$OMP CRITICAL
IF (MINRAD_PRIV < MINRAD) THEN
   MINRAD = MINRAD_PRIV; SIMPS(2,MI) = VERTEX_PRIV
END IF
!$OMP END CRITICAL
!$OMP END PARALLEL
! This is the end of the Level 2 parallel block.

! Set up the first row of the system A X = B.
A(:,1) = PTS(:,SIMPS(2,MI)) - PTS(:,SIMPS(1,MI))
B(1) = DDOT(D, A(:,1), 1, A(:,1), 1) / 2.0_R8
! Initialize the error flag.
IERR(MI) = 0

! Loop to collect the remaining D+1 vertices for the first simplex.
DO I = 2, D
   ! Re-initialize the radius for each iteration.
   MINRAD = HUGE(0.0_R8)
   MINRAD_PRIV = HUGE(0.0_R8)
   VERTEX_PRIV = 0

   ! This is another Level 2 parallel block over all P* in PTS.
   !$OMP PARALLEL FIRSTPRIVATE(MINRAD_PRIV, VERTEX_PRIV, B), &
   !$OMP& PRIVATE(CURRRAD, LQ, CENTER, WORK, TAU, IPIV, ITMP, IERR_PRIV), &
   !$OMP& IF(PLVL2)
   !$OMP DO SCHEDULE(STATIC)
   DO J = 1, N
      ! Check that this point is not already in the simplex.
      IF (ANY(SIMPS(:,MI) .EQ. J)) CYCLE
      ! Add P* to LS system, and compute the minimum norm solution.
      LQ(:,1:I-1) = A(:,1:I-1)
      LQ(:,I) = PTS(:,J) - PTS(:,SIMPS(1,MI))
      B(I) = DDOT(D, LQ(:,I), 1, LQ(:,I), 1) / 2.0_R8
      ! Compute A^T P = Q R.
      CALL DGEQP3(D, I, LQ, D, IPIV, TAU, WORK, LWORK, IERR_PRIV)
      IF(IERR_PRIV < 0) THEN ! LAPACK illegal input error.
         !$OMP ATOMIC WRITE
         IERR(MI) = 80
      ELSE IF (ABS(LQ(I,I)) < EPSL) THEN ! A is rank-deficient.
         CYCLE ! If rank-deficient, skip this point.
      END IF
      ! Set CENTER = P^T B.
      FORALL (ITMP = 1:I) CENTER(ITMP) = B(IPIV(ITMP))
      ! Get the radius using R^T Q^T X = P^T B.
      CALL DTRSM('L', 'U', 'T', 'N', I, 1, 1.0_R8, LQ, D, CENTER, D)
      CENTER(I+1:D) = 0
      CALL DORMQR('L', 'N', D, 1, I, LQ, D, TAU, CENTER, D, WORK, LWORK, &
         IERR_PRIV)
      IF(IERR_PRIV < 0) THEN ! LAPACK illegal input error.
         !$OMP ATOMIC WRITE
         IERR(MI) = 83
      END IF
      ! Calculate the radius and compare it to the current minimum.
      CURRRAD = DNRM2(D, CENTER, 1)
      IF (CURRRAD < MINRAD_PRIV) THEN
         MINRAD_PRIV = CURRRAD; VERTEX_PRIV = J
      END IF
   END DO
   !$OMP CRITICAL
   IF (MINRAD_PRIV < MINRAD) THEN
      MINRAD = MINRAD_PRIV; SIMPS(I+1,MI) = VERTEX_PRIV
   END IF
   !$OMP END CRITICAL
   !$OMP END PARALLEL
   ! End of Level 2 parallel block.

   ! Check the final error flag.
   IF (IERR(MI) .NE. 0) CYCLE OUTER
   ! Check that a point was found. If not, then all the points must lie in a
   ! lower dimensional linear manifold (error case).
   IF (SIMPS(I+1,MI) .EQ. 0) THEN; IERR(MI) = 31; CYCLE OUTER; END IF
   ! If all operations were successful, add the best P* to the LS system.
   A(:,I) = PTS(:,SIMPS(I+1,MI)) - PTS(:,SIMPS(1,MI))
   B(I) = DDOT(D, A(:,I), 1, A(:,I), 1) / 2.0_R8 
END DO
! RETURN
! END SUBROUTINE MAKEFIRSTSIMP
! This marks the end of the in-lined MAKEFIRSTSIMP() subroutine call.

   ! Otherwise, use the seed.
   ELSE
      ! Copy the seed to the current simplex.
      SIMPS(:,MI) = SEED(:)
      ! Rebuild the linear system.
      DO J=1,D
         A(:,J) = PTS(:,SIMPS(J+1,MI)) - PTS(:,SIMPS(1,MI))
         B(J) = DDOT(D, A(:,J), 1, A(:,J), 1) / 2.0_R8 
      END DO
   END IF

! Compute the last index of this chunk.
IF (.NOT. PLVL1) THEN
! If no Level 1 parallelism is exploited, then it is safe to check all future
! interpolation points.
   LAST = M 
ELSE
! If PMODE = 1 or 3, Level 1 parallelism is exploited, and it is only safe
! to check interpolation points up to the end of this chunk.
   LAST = MIN(MI - MOD(MI-1, CHUNKSIZEL) + CHUNKSIZEL - 1, M)
END IF

   ! Inner loop searching for a simplex containing the point Q(:,MI).
   INNER : DO K = 1, IBUDGETL

      ! Check if the current simplex contains Q(:,MI).
!
! Due to OpenMP's handling of variable scope, the parallel implementation of
! the subroutine PTINSIMP() has been in-lined here.
!
! FUNCTION PTINSIMP() RESULT(TF)
! Determine if any interpolation points are in the current simplex, whose
! vertices (P_1, P_2, ..., P_{D+1}) are indexed by SIMPS(:,MI). These
! vertices determine a positive cone with generators V_I = P_{I+1} - P_1,
! I = 1, ..., D. For each interpolation point Q* in Q, Q* - P_1 can be
! expressed as a unique linear combination of the V_I.  If all these linear
! weights are nonnegative and sum to less than or equal to 1.0, then Q* is
! in the simplex with vertices {P_I}_{I=1}^{D+1}.
!
! If any interpolation points in Q are contained in the simplex whose
! vertices are indexed by SIMPS(:,MI), then those points are marked as solved
! and the values of SIMPS and WEIGHTS are updated appropriately. On output,
! WEIGHTS(:,MI) contains the affine weights for producing Q(:,MI) as an
! affine combination of the points in PTS indexed by SIMPS(:,MI). If these
! weights are nonnegative, then PTINSIMP() returns TRUE.

! Initialize the return value and local variables.
PTINSIMP = .FALSE.

! Compute the LU factorization of the matrix A whose rows are P_{I+1} - P_1.
LQ = A
CALL DGETRF(D, D, LQ, D, IPIV, IERR(MI))
IF (IERR(MI) < 0) THEN ! LAPACK illegal input.
   IERR(MI) = 81; CYCLE OUTER
ELSE IF (IERR(MI) > 0) THEN ! Rank-deficiency detected.
        IERR(MI) = 61; CYCLE OUTER
END IF
! Solve A^T w = WORK to get the affine weights for Q(:,MI) or its projection.
WORK(1:D) = PROJ(:) - PTS(:,SIMPS(1,MI))
CALL DGETRS('N', D, 1, LQ, D, IPIV, WORK(1:D), D, IERR(MI))
IF (IERR(MI) < 0) THEN ! LAPACK illegal input.
   IERR(MI) = 82; CYCLE OUTER
END IF
WEIGHTS(2:D+1,MI) = WORK(1:D)
WEIGHTS(1,MI) = 1.0_R8 - SUM(WEIGHTS(2:D+1,MI))
! Check if the weights for Q(:,MI) are nonnegative.
IF (ALL(WEIGHTS(:,MI) .GE. -EPSL)) PTINSIMP = .TRUE.

! Level 2 parallel block over current chunk of interoplation points.
! Compute the affine weights for the interpolation points in this chunk.
! If PMODE=2, then LAST=M and this loop is carried out over all remaining
! interpolation points. Uses PLANE(:) as a work array.
!$OMP PARALLEL DO PRIVATE(PLANE, ITMP), SCHEDULE(STATIC), IF(PLVL2)
DO I = MI+1, LAST
   ! Check that no solution has already been found.
   IF (IERR(I) .NE. 40) CYCLE
   ! Solve A^T w = PLANE to get the affine weights for Q(:,I).
   PLANE(2:D+1) = Q(:,I) - PTS(:,SIMPS(1,MI))
   CALL DGETRS('N', D, 1, LQ, D, IPIV, PLANE(2:D+1), D, ITMP)
   IF (ITMP < 0) CYCLE ! Illegal input error that should never occurr.
   ! Check if the weights define a convex combination.
   PLANE(1) = 1.0_R8 - SUM(PLANE(2:D+1))
   IF (ALL(PLANE(1:D+1) .GE. -EPSL)) THEN
      ! Copy the simplex indices and weights then flag as complete.
      SIMPS(:,I) = SIMPS(:,MI)
      WEIGHTS(:,I) = PLANE(1:D+1)
      IERR(I) = 0
   END IF
END DO
!$OMP END PARALLEL DO
! End of Level 2 parallel block.
! RETURN
! END FUNCTION PTINSIMP
! This marks the end of the in-lined MAKEFIRSTSIMP() subroutine call.

      IF (PTINSIMP) THEN
         ! If chaining is on, set the next seed to this simplex.
         IF (CHAINL) SEED(:) = SIMPS(:,MI)
         EXIT INNER
      END IF

      ! If chaining is on, save each good simplex as the next seed.
      IF (CHAINL) SEED(:) = SIMPS(:,MI)

      ! Swap out the least weighted vertex, but save its value in case it
      ! needs to be restored later.
      I = MINLOC(WEIGHTS(1:D+1,MI), DIM=1)
      ITMP = SIMPS(I,MI)
      SIMPS(I,MI) = SIMPS(D+1,MI)

      ! If the least weighted vertex (I) is not the first vertex, then just
      ! drop row I from the linear system.
      IF(I .NE. 1) THEN
         A(:,I-1) = A(:,D); B(I-1) = B(D)
      ! However, if I = 1, then both A and B must be reconstructed.
      ELSE
         DO J=1,D
            A(:,J) = PTS(:,SIMPS(J+1,MI)) - PTS(:,SIMPS(1,MI))
            B(J) = DDOT(D, A(:,J), 1, A(:,J), 1) / 2.0_R8 
         END DO
      END IF

      ! Compute the next simplex (do one flip).
!      CALL MAKESIMPLEX(); IF (IERR(MI) .NE. 0) CYCLE OUTER
!
! Due to OpenMP's handling of variable scope, the parallel implementation of
! the subroutine MAKESIMPLEX() has been in-lined here.
!
! SUBROUTINE MAKESIMPLEX()
! Given a Delaunay facet F whose containing hyperplane does not contain
! Q(:,MI), complete the simplex by adding a point from PTS on the same `side'
! of F as Q(:,MI). Assume SIMPS(1:D,MI) contains the vertex indices of F
! (corresponding to data points P_1, P_2, ..., P_D in PTS), and assume the
! matrix A(1:D-1,:)^T and vector B(1:D-1) are filled appropriately (similarly
! as in MAKEFIRSTSIMP()). Then for any P* (not in the hyperplane containing
! F) in PTS, let CENTER denote the circumcenter of the simplex with vertices
! P_1, P_2, ..., P_D, P*. Then
!
! X = CENTER - P_1
!
! is given by the solution to the nonsingular linear system
!
! AX = B where
!
! A = [ P_2 - P_1, P_3 - P_1, ..., P_D - P_1, P* - P_1 ]^T and
! B = [ <A_{1.},A_{1.}>/2, <A_{2.},A_{2.}>/2, ..., <A_{D.},A_{D.}>/2 ]^T.
!
! Then CENTER = X + P_1 and RADIUS = \| X \|.  P_{D+1} will be given by the 
! candidate P* that satisfies both of the following:
!
! 1) Let PLANE denote the hyperplane containing F. Then P_{D+1} and Q(:,MI) 
! must be on the same side of PLANE.
!
! 2) The circumball about CENTER must not contain any points in PTS in its
! interior (Delaunay property).
! 
! The above are necessary and sufficient conditions for flipping the
! Delaunay simplex, given that F is indeed a Delaunay facet.
!
! On input, SIMPS(1:D,MI) should contain the vertex indices (column indices
! from PTS) of the facet F.  Upon output, SIMPS(:,MI) will contain the vertex
! indices of a Delaunay simplex closer to Q(:,MI).  Also, the matrix A^T and
! vector B will be updated accordingly. If SIMPS(D+1,MI)=0, then there were
! no points in PTS on the appropriate side of F, meaning that Q(:,MI) is an
! extrapolation point (not a convex combination of points in PTS).

! Compute a hyperplane c^T x = \alpha containing the first D vertices indexed
! in SIMPS(:,MI). The plane is determined by its normal vector c and \alpha.
! Let (P_1, P_2, ..., P_D) be the vertices indexed in SIMPS(1:D,MI). A normal
! vector is any nonzero vector in ker(A^T), where the matrix
! A^T = [ P_2 - P_1, P_3 - P_1, ..., P_D - P_1 ].
! Since rank A^T = D-1, dim ker(A^T) = 1, and ker(A^T) can be found from a QR
! factorization of A^T:  A^T P = QR, where P permutes the columns of A^T.  Then
! the last column of Q is orthonormal to the range of A^T.
IF (D > 1) THEN ! Check that D-1 > 0, otherwise the plane is trivial.
   ! Compute the QR factorization.
   IPIV=0
   LQ = A
   CALL DGEQP3(D, D-1, LQ, D, IPIV, TAU, WORK, LWORK, IERR(MI))
   IF(IERR(MI) < 0) THEN ! LAPACK illegal input error.
      IERR(MI) = 80; CYCLE OUTER
   END IF
   ! The nullspace is given by the last column of Q.
   PLANE(1:D-1) = 0.0_R8
   PLANE(D) = 1.0_R8
   CALL DORMQR('L', 'N', D, 1, D-1, LQ, D, TAU, PLANE, D, WORK, &
    LWORK, IERR(MI))
   IF(IERR(MI) < 0) THEN ! LAPACK illegal input error.
      IERR(MI) = 83; CYCLE OUTER
   END IF
   ! Calculate the constant \alpha defining the plane.
   PLANE(D+1) = DDOT(D,PLANE(1:D),1,PTS(:,SIMPS(1,MI)),1)
ELSE ! Special case where D=1.
   PLANE(1) = 1.0_R8
   PLANE(2) = PTS(1,SIMPS(1,MI))
END IF
! Compute the sign for the side of PLANE containing Q(:,MI).
SIDE1 =  DDOT(D,PLANE(1:D),1,PROJ(:),1) - PLANE(D+1)
SIDE1 = SIGN(1.0_R8,SIDE1)
! Initialize the center, radius, simplex, error flag, and OpenMP variabls.
SIMPS(D+1,MI) = 0
CENTER(:) = 0.0_R8
TAU(:) = 0.0_R8
IERR(MI) = 0
MINRAD = HUGE(0.0_R8)
MINRAD_PRIV = HUGE(0.0_R8)
VERTEX_PRIV = 0

! Begin Level 2 parallel loop over all points P* in PTS.
!$OMP PARALLEL FIRSTPRIVATE(MINRAD_PRIV, VERTEX_PRIV), &
!$OMP& PRIVATE(SIDE2, LQ, CENTER, IPIV, IERR_PRIV), &
!$OMP& IF(PLVL2)
!$OMP DO SCHEDULE(STATIC)
DO I = 1, N
   ! Check that P* is inside the current ball.
   IF (DNRM2(D, PTS(:,I) - CENTER(:), 1) > MINRAD_PRIV) CYCLE ! If not, skip.
   ! Check that P* is on the appropriate halfspace.
   SIDE2 = DDOT(D,PLANE(1:D),1,PTS(:,I),1) - PLANE(D+1)
   IF (SIDE1*SIDE2 < EPSL .OR. ANY(SIMPS(:,MI) .EQ. I)) CYCLE ! If not, skip.
   ! Add P* to the linear system, and solve to get shifted CENTER.
   LQ(:,1:D-1) = A(:,1:D-1)
   CENTER(1:D-1) = B(1:D-1)
   LQ(:,D) = PTS(:,I) - PTS(:,SIMPS(1,MI))
   CENTER(D) = DDOT(D, LQ(:,D), 1, LQ(:,D), 1) / 2.0_R8
   ! Compute A^T = LU.
   CALL DGETRF(D, D, LQ, D, IPIV, IERR_PRIV)
   IF (IERR_PRIV < 0) THEN ! LAPACK illegal input.
      !$OMP ATOMIC WRITE
      IERR(MI) = 81
   ELSE IF (IERR_PRIV > 0) THEN ! Rank-deficiency detected.
      !$OMP ATOMIC WRITE
      IERR(MI) = 61
   END IF
   ! Use A^T = LU to solve A X = B, where X = CENTER - P_1.
   CALL DGETRS('T', D, 1, LQ, D, IPIV, CENTER, D, IERR_PRIV)
   IF (IERR_PRIV < 0) THEN ! LAPACK illegal input.
      !$OMP ATOMIC WRITE
      IERR(MI) = 82
   END IF
   ! Update the new radius, center, and simplex.
   MINRAD_PRIV = DNRM2(D, CENTER, 1)
   CENTER(:) = CENTER(:) + PTS(:,SIMPS(1,MI))
   VERTEX_PRIV = I
END DO
!$OMP CRITICAL
! Check if PTS(:,VERTEX_PRIV) is inside the circumball.
IF (DNRM2(D, PTS(:,VERTEX_PRIV) - TAU(:), 1) < MINRAD) THEN
   MINRAD = MINRAD_PRIV
   TAU(:) = CENTER(:)
   SIMPS(D+1,MI) = VERTEX_PRIV
END IF
!$OMP END CRITICAL
!$OMP END PARALLEL
! End level 2 parallel region.

! Check for error flags.
IF(IERR(MI) .NE. 0) CYCLE OUTER
! Check for extrapolation condition.
IF(SIMPS(D+1,MI) .NE. 0) THEN
   ! Add new point to the linear system.
   A(:,D) = PTS(:,SIMPS(D+1,MI)) - PTS(:,SIMPS(1,MI))
   B(D) = DDOT(D, A(:,D), 1, A(:,D), 1) / 2.0_R8
END IF
! RETURN
! END SUBROUTINE MAKESIMPLEX
! End of in-lined code for MAKESIMPLEX().

      ! If no vertex was found, then this is an extrapolation point.
      IF (SIMPS(D+1,MI) .EQ. 0) THEN
         ! If extrapolation is not allowed (EXTRAP=0), do not proceed.
         IF (EXTRAPL < EPSL) THEN
            SIMPS(:,MI) = 0; WEIGHTS(:,MI) = 0  ! Zero all output values.
            ! Set the error flag and skip this point.
            IERR(MI) = 2; CYCLE OUTER
         END IF

         ! Otherwise, project the extrapolation point onto the convex hull.
!         CALL PROJECT(); IF (IERR(MI) .NE. 0) CYCLE OUTER
!
! Due to OpenMP's handling of variable scope, the parallel (identical to serial)
! implementation of the subroutine PROJECT() has been in-lined here.
!
! SUBROUTINE PROJECT()
! Project a point outside the convex hull of the point set onto the convex hull
! by solving an inequality constrained least squares problem. The solution to
! the least squares problem gives the projection as a convex combination of the
! data points. The projection can then be computed by performing a matrix
! vector multiplication.

! Allocate work arrays.
IF (.NOT. ALLOCATED(IWORK_DWNNLS)) THEN
   ALLOCATE(IWORK_DWNNLS(D+1+N), STAT=IERR(MI))
   IF(IERR(MI) .NE. 0) THEN; IERR(MI) = 70; CYCLE OUTER; END IF
END IF
IF (.NOT. ALLOCATED(WORK_DWNNLS)) THEN
   ALLOCATE(WORK_DWNNLS(D+1+N*5), STAT=IERR(MI))
   IF(IERR(MI) .NE. 0) THEN; IERR(MI) = 70; CYCLE OUTER; END IF
END IF
IF (.NOT. ALLOCATED(W_DWNNLS)) THEN
   ALLOCATE(W_DWNNLS(D+1,N+1), STAT=IERR(MI))
   IF(IERR(MI) .NE. 0) THEN; IERR(MI) = 70; CYCLE OUTER; END IF
END IF
IF (.NOT. ALLOCATED(X_DWNNLS)) THEN
   ALLOCATE(X_DWNNLS(N), STAT=IERR(MI))
   IF(IERR(MI) .NE. 0) THEN; IERR(MI) = 70; CYCLE OUTER; END IF
END IF

! Initialize work array and settings values.
IWORK_DWNNLS(1) = D+1+5*N
IWORK_DWNNLS(2) = D+1+N
W_DWNNLS(1, :) = 1.0_R8          ! Set convexity (equality) constraint.
W_DWNNLS(2:D+1,1:N) = PTS(:,:)   ! Copy data points.
W_DWNNLS(2:D+1,N+1) = PROJ(:)    ! Copy extrapolation point.
! Compute the solution to the inequality constrained least squares problem to
! get the projection coefficients.
CALL DWNNLS(W_DWNNLS, D+1, 1, D, N, 0, PRGOPT_DWNNLS, X_DWNNLS, RNORML, &
   IERR(MI), IWORK_DWNNLS, WORK_DWNNLS)
IF (IERR(MI) .EQ. 1) THEN ! Failure to converge.
   IERR(MI) = 71; CYCLE OUTER
ELSE IF (IERR(MI) .EQ. 2) THEN ! Illegal input detected.
   IERR(MI) = 72; CYCLE OUTER
END IF
! Compute the actual projection via matrix vector multiplication.
CALL DGEMV('N', D, N, 1.0_R8, PTS, D, X_DWNNLS, 1, 0.0_R8, PROJ, 1)
! RETURN
! END SUBROUTINE PROJECT
! End of in-lined code for PROJECT().

         ! Check the value of RNORML for over-extrapolation.
         IF (RNORML > EXTRAPL * PTS_DIAM) THEN
            SIMPS(:,MI) = 0; WEIGHTS(:,MI) = 0  ! Zero all output values.
            ! If present, record the unscaled RNORM output.
            IF (PRESENT(RNORM)) RNORM(MI) = RNORML*PTS_SCALE
            ! Set the error flag and skip this point.
            IERR(MI) = 2; CYCLE OUTER
         END IF

         ! Otherwise, restore the previous simplex and continue with the 
         ! projected value.
         SIMPS(D+1,MI) = ITMP
         A(:,D) = PTS(:,ITMP) - PTS(:,SIMPS(1,MI))
         B(D) = DDOT(D, A(:,D), 1, A(:,D), 1) / 2.0_R8
      END IF

   ! End of inner loop for finding each interpolation point.
   END DO INNER

   ! Check for budget violation conditions.
   IF (K > IBUDGETL) THEN 
      SIMPS(:,MI) = 0; WEIGHTS(:,MI) = 0  ! Zero all output values.
      ! Set the error flag and skip this point.
      IERR(MI) = 60; CYCLE OUTER
   END IF

   ! If the residual is nonzero, set the extrapolation flag.
   IF (RNORML > EPSL) IERR(MI) = 1

   ! If present, record the RNORM output.
   IF (PRESENT(RNORM)) RNORM(MI) = RNORML*PTS_SCALE

   ! Compute interpolation point response values.
   IF (PRESENT(INTERP_IN)) THEN
      ! Check for errors.
      IF (IERR(MI) .LE. 1) THEN
         ! Compute the weighted sum of vertex response values.
         DO K = 1, D+1
            INTERP_OUT(:,MI) = INTERP_OUT(:,MI) &
             + INTERP_IN(:,SIMPS(K,MI)) * WEIGHTS(K,MI)
         END DO
      END IF
   END IF

END DO OUTER  ! End of outer loop over all interpolation points.
!$OMP END DO

! Free optional work arrays.
IF (ALLOCATED(IWORK_DWNNLS)) DEALLOCATE(IWORK_DWNNLS)
IF (ALLOCATED(WORK_DWNNLS))  DEALLOCATE(WORK_DWNNLS)
IF (ALLOCATED(W_DWNNLS))     DEALLOCATE(W_DWNNLS)
IF (ALLOCATED(X_DWNNLS))     DEALLOCATE(X_DWNNLS)
!$OMP END PARALLEL 
! End of Level 1 parallel region.

! Free dynamic work arrays.
DEALLOCATE(WORK)

RETURN

CONTAINS    ! Internal subroutines and functions.

SUBROUTINE RESCALE(MINDIST, DIAMETER, SCALE)
! Rescale and transform data to be centered at the origin with unit
! radius. This subroutine has O(n^2) complexity.
!
! On output, PTS and Q have been rescaled and shifted. All the data
! points in PTS are centered with unit radius, and the points in Q
! have been shifted and scaled in relation to PTS.
!
! MINDIST is a real number containing the (scaled) minimum distance
!    between any two data points in PTS.
!
! DIAMETER is a real number containing the (scaled) diameter of the
!    data set PTS.
!
! SCALE contains the real factor used to transform the data and
!    interpolation points: scaled value = (original value -
!    barycenter of data points)/SCALE.

! Output arguments.
REAL(KIND=R8), INTENT(OUT) :: MINDIST, DIAMETER, SCALE

! Local variables.
REAL(KIND=R8) :: PTS_CENTER(D) ! The center of the data points PTS.
REAL(KIND=R8) :: DISTANCE ! The current distance.

! Initialize local values.
MINDIST = HUGE(0.0_R8)
DIAMETER = 0.0_R8
SCALE = 0.0_R8

! Compute barycenter of all data points.
PTS_CENTER(:) = SUM(PTS(:,:), DIM=2)/REAL(N, KIND=R8)
! Compute the minimum and maximum distances.
!$OMP PARALLEL DO PRIVATE(DISTANCE), REDUCTION(MAX:DIAMETER), &
!$OMP& REDUCTION(MIN:MINDIST), SCHEDULE(STATIC, 100)
DO I = 1, N ! Cycle through all pairs of points.
   DO J = I + 1, N
      DISTANCE = DNRM2(D, PTS(:,I) - PTS(:,J), 1) ! Compute the distance.
      IF (DISTANCE > DIAMETER) THEN ! Compare to the current diameter.
         DIAMETER = DISTANCE
      END IF
      IF (DISTANCE < MINDIST) THEN ! Compare to the current minimum distance.
         MINDIST = DISTANCE
      END IF
   END DO
END DO
!$OMP END PARALLEL DO
! Center the points.
FORALL (I = 1:N) PTS(:,I) = PTS(:,I) - PTS_CENTER(:)
! Compute the scale factor (for unit radius).
!$OMP PARALLEL DO PRIVATE(DISTANCE), REDUCTION(MAX:SCALE), SCHEDULE(STATIC)
DO I = 1, N ! Cycle through all points again.
   DISTANCE = DNRM2(D, PTS(:,I), 1) ! Compute the distance from the center.
   IF (DISTANCE > SCALE) THEN ! Compare to the current radius.
      SCALE = DISTANCE
   END IF
END DO
!$OMP END PARALLEL DO
! Scale the points to unit radius.
PTS = PTS / SCALE
! Also transform Q similarly.
FORALL (I = 1:M) Q(:,I) = (Q(:,I) - PTS_CENTER(:)) / SCALE
! Also scale the minimum distance and diameter.
MINDIST = MINDIST / SCALE
DIAMETER = DIAMETER / SCALE
RETURN
END SUBROUTINE RESCALE

END SUBROUTINE DELAUNAYSPARSEP
