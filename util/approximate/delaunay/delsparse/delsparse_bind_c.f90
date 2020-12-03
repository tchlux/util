! This automatically generated Fortran wrapper file allows codes
! written in Fortran to be called directly from C and translates all
! C-style arguments into expected Fortran-style arguments (with
! assumed size, local type declarations, etc.).


SUBROUTINE C_DELAUNAYSPARSES(D, N, PTS_DIM_1, PTS_DIM_2, PTS, M, Q_DIM_1, Q_DIM_2, Q, SIMPS_DIM_1, SIMPS_DIM_2, SIMPS, WEIGHTS_DIM_&
&1, WEIGHTS_DIM_2, WEIGHTS, IERR_DIM_1, IERR, INTERP_IN_PRESENT, INTERP_IN_DIM_1, INTERP_IN_DIM_2, INTERP_IN, INTERP_OUT_PRESENT, I&
&NTERP_OUT_DIM_1, INTERP_OUT_DIM_2, INTERP_OUT, EPS_PRESENT, EPS, EXTRAP_PRESENT, EXTRAP, RNORM_PRESENT, RNORM_DIM_1, RNORM, IBUDGE&
&T_PRESENT, IBUDGET, CHAIN_PRESENT, CHAIN) BIND(C)
  USE REAL_PRECISION , ONLY : R8
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: D
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN) :: PTS_DIM_1
  INTEGER, INTENT(IN) :: PTS_DIM_2
  REAL(KIND=R8), INTENT(INOUT), DIMENSION(PTS_DIM_1,PTS_DIM_2) :: PTS
  INTEGER, INTENT(IN) :: M
  INTEGER, INTENT(IN) :: Q_DIM_1
  INTEGER, INTENT(IN) :: Q_DIM_2
  REAL(KIND=R8), INTENT(INOUT), DIMENSION(Q_DIM_1,Q_DIM_2) :: Q
  INTEGER, INTENT(IN) :: SIMPS_DIM_1
  INTEGER, INTENT(IN) :: SIMPS_DIM_2
  INTEGER, INTENT(OUT), DIMENSION(SIMPS_DIM_1,SIMPS_DIM_2) :: SIMPS
  INTEGER, INTENT(IN) :: WEIGHTS_DIM_1
  INTEGER, INTENT(IN) :: WEIGHTS_DIM_2
  REAL(KIND=R8), INTENT(OUT), DIMENSION(WEIGHTS_DIM_1,WEIGHTS_DIM_2) :: WEIGHTS
  INTEGER, INTENT(IN) :: IERR_DIM_1
  INTEGER, INTENT(OUT), DIMENSION(IERR_DIM_1) :: IERR
  LOGICAL, INTENT(IN) :: INTERP_IN_PRESENT
  INTEGER, INTENT(IN) :: INTERP_IN_DIM_1
  INTEGER, INTENT(IN) :: INTERP_IN_DIM_2
  REAL(KIND=R8), INTENT(IN), DIMENSION(INTERP_IN_DIM_1,INTERP_IN_DIM_2) :: INTERP_IN
  LOGICAL, INTENT(IN) :: INTERP_OUT_PRESENT
  INTEGER, INTENT(IN) :: INTERP_OUT_DIM_1
  INTEGER, INTENT(IN) :: INTERP_OUT_DIM_2
  REAL(KIND=R8), INTENT(OUT), DIMENSION(INTERP_OUT_DIM_1,INTERP_OUT_DIM_2) :: INTERP_OUT
  LOGICAL, INTENT(IN) :: EPS_PRESENT
  REAL(KIND=R8), INTENT(IN) :: EPS
  LOGICAL, INTENT(IN) :: EXTRAP_PRESENT
  REAL(KIND=R8), INTENT(IN) :: EXTRAP
  LOGICAL, INTENT(IN) :: RNORM_PRESENT
  INTEGER, INTENT(IN) :: RNORM_DIM_1
  REAL(KIND=R8), INTENT(OUT), DIMENSION(RNORM_DIM_1) :: RNORM
  LOGICAL, INTENT(IN) :: IBUDGET_PRESENT
  INTEGER, INTENT(IN) :: IBUDGET
  LOGICAL, INTENT(IN) :: CHAIN_PRESENT
  LOGICAL, INTENT(IN) :: CHAIN

  INTERFACE
    SUBROUTINE DELAUNAYSPARSES(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, INTERP_IN, INTERP_OUT, EPS, EXTRAP, RNORM, IBUDGET, CHAIN)
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
      !    EPS is assigned \sqrt{{\mu}} where \mu denotes the unit roundoff for the
      !    machine. In general, any values that differ by less than EPS are judged
      !    as equal, and any weights that are greater than -EPS are judged as
      !    nonnegative.  EPS cannot take a value less than the default value of
      !    \sqrt{{\mu}}. If any value less than \sqrt{{\mu}} is supplied, the default
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
      ! Last Update: February, 2019
      !
      USE REAL_PRECISION , ONLY : R8
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: D
      INTEGER, INTENT(IN) :: N
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(:,:) :: PTS
      INTEGER, INTENT(IN) :: M
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(:,:) :: Q
      INTEGER, INTENT(OUT), DIMENSION(:,:) :: SIMPS
      REAL(KIND=R8), INTENT(OUT), DIMENSION(:,:) :: WEIGHTS
      INTEGER, INTENT(OUT), DIMENSION(:) :: IERR
      REAL(KIND=R8), INTENT(IN), OPTIONAL, DIMENSION(:,:) :: INTERP_IN
      REAL(KIND=R8), INTENT(OUT), OPTIONAL, DIMENSION(:,:) :: INTERP_OUT
      REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
      REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
      REAL(KIND=R8), INTENT(OUT), OPTIONAL, DIMENSION(:) :: RNORM
      INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
      LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
    END SUBROUTINE DELAUNAYSPARSES
  END INTERFACE

  IF (INTERP_IN_PRESENT) THEN
    IF (INTERP_OUT_PRESENT) THEN
      IF (EPS_PRESENT) THEN
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP)
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EPS=EPS)
              END IF
            END IF
          END IF
        END IF
      ELSE
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, EXTRAP=EXTRAP)
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, INT&
&ERP_OUT=INTERP_OUT)
              END IF
            END IF
          END IF
        END IF
      END IF
    ELSE
      IF (EPS_PRESENT) THEN
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, EXTRAP=EXTRAP, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, EXTRAP=EXTRAP, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, EXTRAP=EXTRAP)
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EPS&
&=EPS)
              END IF
            END IF
          END IF
        END IF
      ELSE
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EXT&
&RAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EXT&
&RAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EXT&
&RAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EXT&
&RAP=EXTRAP, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EXT&
&RAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EXT&
&RAP=EXTRAP, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EXT&
&RAP=EXTRAP, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, EXT&
&RAP=EXTRAP)
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, RNO&
&RM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, RNO&
&RM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, RNO&
&RM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, RNO&
&RM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, IBU&
&DGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, IBU&
&DGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, CHA&
&IN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN)
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF
  ELSE
    IF (INTERP_OUT_PRESENT) THEN
      IF (EPS_PRESENT) THEN
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, EXTRAP=EXTRAP, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, EXTRAP=EXTRAP)
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&PS=EPS)
              END IF
            END IF
          END IF
        END IF
      ELSE
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&XTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&XTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&XTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&XTRAP=EXTRAP, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&XTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&XTRAP=EXTRAP, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&XTRAP=EXTRAP, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, E&
&XTRAP=EXTRAP)
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, R&
&NORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, R&
&NORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, R&
&NORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, R&
&NORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, I&
&BUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, I&
&BUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT, C&
&HAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT)
              END IF
            END IF
          END IF
        END IF
      END IF
    ELSE
      IF (EPS_PRESENT) THEN
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP, &
&RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP, &
&RNORM=RNORM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP, &
&RNORM=RNORM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP, &
&RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP, &
&IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP, &
&IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP, &
&CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP)
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM, IB&
&UDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM, IB&
&UDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM, CH&
&AIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, IBUDGET=IBUDGET&
&, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, IBUDGET=IBUDGET&
&)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS)
              END IF
            END IF
          END IF
        END IF
      ELSE
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=RNO&
&RM, IBUDGET=IBUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=RNO&
&RM, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=RNO&
&RM, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=RNO&
&RM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, IBUDGET=I&
&BUDGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, IBUDGET=I&
&BUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, CHAIN=CHA&
&IN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP)
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM, IBUDGET=IBU&
&DGET, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM, IBUDGET=IBU&
&DGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM, CHAIN=CHAIN&
&)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM)
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, IBUDGET=IBUDGET, CHAIN=C&
&HAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, IBUDGET=IBUDGET)
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, CHAIN=CHAIN)
              ELSE
                CALL DELAUNAYSPARSES(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR)
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF
  END IF
END SUBROUTINE C_DELAUNAYSPARSES


SUBROUTINE C_DELAUNAYSPARSEP(D, N, PTS_DIM_1, PTS_DIM_2, PTS, M, Q_DIM_1, Q_DIM_2, Q, SIMPS_DIM_1, SIMPS_DIM_2, SIMPS, WEIGHTS_DIM_&
&1, WEIGHTS_DIM_2, WEIGHTS, IERR_DIM_1, IERR, INTERP_IN_PRESENT, INTERP_IN_DIM_1, INTERP_IN_DIM_2, INTERP_IN, INTERP_OUT_PRESENT, I&
&NTERP_OUT_DIM_1, INTERP_OUT_DIM_2, INTERP_OUT, EPS_PRESENT, EPS, EXTRAP_PRESENT, EXTRAP, RNORM_PRESENT, RNORM_DIM_1, RNORM, IBUDGE&
&T_PRESENT, IBUDGET, CHAIN_PRESENT, CHAIN, PMODE_PRESENT, PMODE) BIND(C)
  USE REAL_PRECISION , ONLY : R8
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: D
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN) :: PTS_DIM_1
  INTEGER, INTENT(IN) :: PTS_DIM_2
  REAL(KIND=R8), INTENT(INOUT), DIMENSION(PTS_DIM_1,PTS_DIM_2) :: PTS
  INTEGER, INTENT(IN) :: M
  INTEGER, INTENT(IN) :: Q_DIM_1
  INTEGER, INTENT(IN) :: Q_DIM_2
  REAL(KIND=R8), INTENT(INOUT), DIMENSION(Q_DIM_1,Q_DIM_2) :: Q
  INTEGER, INTENT(IN) :: SIMPS_DIM_1
  INTEGER, INTENT(IN) :: SIMPS_DIM_2
  INTEGER, INTENT(OUT), DIMENSION(SIMPS_DIM_1,SIMPS_DIM_2) :: SIMPS
  INTEGER, INTENT(IN) :: WEIGHTS_DIM_1
  INTEGER, INTENT(IN) :: WEIGHTS_DIM_2
  REAL(KIND=R8), INTENT(OUT), DIMENSION(WEIGHTS_DIM_1,WEIGHTS_DIM_2) :: WEIGHTS
  INTEGER, INTENT(IN) :: IERR_DIM_1
  INTEGER, INTENT(OUT), DIMENSION(IERR_DIM_1) :: IERR
  LOGICAL, INTENT(IN) :: INTERP_IN_PRESENT
  INTEGER, INTENT(IN) :: INTERP_IN_DIM_1
  INTEGER, INTENT(IN) :: INTERP_IN_DIM_2
  REAL(KIND=R8), INTENT(IN), DIMENSION(INTERP_IN_DIM_1,INTERP_IN_DIM_2) :: INTERP_IN
  LOGICAL, INTENT(IN) :: INTERP_OUT_PRESENT
  INTEGER, INTENT(IN) :: INTERP_OUT_DIM_1
  INTEGER, INTENT(IN) :: INTERP_OUT_DIM_2
  REAL(KIND=R8), INTENT(OUT), DIMENSION(INTERP_OUT_DIM_1,INTERP_OUT_DIM_2) :: INTERP_OUT
  LOGICAL, INTENT(IN) :: EPS_PRESENT
  REAL(KIND=R8), INTENT(IN) :: EPS
  LOGICAL, INTENT(IN) :: EXTRAP_PRESENT
  REAL(KIND=R8), INTENT(IN) :: EXTRAP
  LOGICAL, INTENT(IN) :: RNORM_PRESENT
  INTEGER, INTENT(IN) :: RNORM_DIM_1
  REAL(KIND=R8), INTENT(OUT), DIMENSION(RNORM_DIM_1) :: RNORM
  LOGICAL, INTENT(IN) :: IBUDGET_PRESENT
  INTEGER, INTENT(IN) :: IBUDGET
  LOGICAL, INTENT(IN) :: CHAIN_PRESENT
  LOGICAL, INTENT(IN) :: CHAIN
  LOGICAL, INTENT(IN) :: PMODE_PRESENT
  INTEGER, INTENT(IN) :: PMODE

  INTERFACE
    SUBROUTINE DELAUNAYSPARSEP(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, INTERP_IN, INTERP_OUT, EPS, EXTRAP, RNORM, IBUDGET, CHAIN, PM&
&ODE)
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
      !    default, EPS is assigned \sqrt{{\mu}} where \mu denotes the unit roundoff
      !    for the machine. In general, any values that differ by less than EPS
      !    are judged as equal, and any weights that are greater than -EPS are
      !    judged as nonnegative.  EPS cannot take a value less than the default
      !    value of \sqrt{{\mu}}. If any value less than \sqrt{{\mu}} is supplied,
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
      ! Last Update: February, 2019
      !
      USE REAL_PRECISION , ONLY : R8
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: D
      INTEGER, INTENT(IN) :: N
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(:,:) :: PTS
      INTEGER, INTENT(IN) :: M
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(:,:) :: Q
      INTEGER, INTENT(OUT), DIMENSION(:,:) :: SIMPS
      REAL(KIND=R8), INTENT(OUT), DIMENSION(:,:) :: WEIGHTS
      INTEGER, INTENT(OUT), DIMENSION(:) :: IERR
      REAL(KIND=R8), INTENT(IN), OPTIONAL, DIMENSION(:,:) :: INTERP_IN
      REAL(KIND=R8), INTENT(OUT), OPTIONAL, DIMENSION(:,:) :: INTERP_OUT
      REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
      REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
      REAL(KIND=R8), INTENT(OUT), OPTIONAL, DIMENSION(:) :: RNORM
      INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
      LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
      INTEGER, INTENT(IN), OPTIONAL :: PMODE
    END SUBROUTINE DELAUNAYSPARSEP
  END INTERFACE

  IF (INTERP_IN_PRESENT) THEN
    IF (INTERP_OUT_PRESENT) THEN
      IF (EPS_PRESENT) THEN
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, EXTRAP=EXTRAP)
                END IF
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EPS=EPS)
                END IF
              END IF
            END IF
          END IF
        END IF
      ELSE
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, EXTRAP=EXTRAP)
                END IF
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&NTERP_OUT=INTERP_OUT)
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    ELSE
      IF (EPS_PRESENT) THEN
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, EXTRAP=EXTRAP)
                END IF
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&PS=EPS)
                END IF
              END IF
            END IF
          END IF
        END IF
      ELSE
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, E&
&XTRAP=EXTRAP)
                END IF
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, R&
&NORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, R&
&NORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, R&
&NORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, R&
&NORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, R&
&NORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, R&
&NORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, R&
&NORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, R&
&NORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&BUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&BUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&BUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, I&
&BUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, C&
&HAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, C&
&HAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN, P&
&MODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_IN=INTERP_IN)
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF
  ELSE
    IF (INTERP_OUT_PRESENT) THEN
      IF (EPS_PRESENT) THEN
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, EXTRAP=EXTRAP)
                END IF
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EPS=EPS)
                END IF
              END IF
            END IF
          END IF
        END IF
      ELSE
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& EXTRAP=EXTRAP)
                END IF
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT,&
& PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, INTERP_OUT=INTERP_OUT)
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    ELSE
      IF (EPS_PRESENT) THEN
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, RNORM=RNORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, RNORM=RNORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, RNORM=RNORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, RNORM=RNORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, RNORM=RNORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, RNORM=RNORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, EXTRAP=EXTRAP&
&)
                END IF
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM, &
&IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM, &
&IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM, &
&IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM, &
&IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM, &
&CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM, &
&CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM, &
&PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, IBUDGET=IBUDG&
&ET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, IBUDGET=IBUDG&
&ET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, IBUDGET=IBUDG&
&ET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, IBUDGET=IBUDG&
&ET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, CHAIN=CHAIN, &
&PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EPS=EPS)
                END IF
              END IF
            END IF
          END IF
        END IF
      ELSE
        IF (EXTRAP_PRESENT) THEN
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=R&
&NORM, IBUDGET=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=R&
&NORM, IBUDGET=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=R&
&NORM, IBUDGET=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=R&
&NORM, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=R&
&NORM, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=R&
&NORM, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=R&
&NORM, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, RNORM=R&
&NORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, IBUDGET&
&=IBUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, IBUDGET&
&=IBUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, IBUDGET&
&=IBUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, IBUDGET&
&=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, CHAIN=C&
&HAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, CHAIN=C&
&HAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP, PMODE=P&
&MODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, EXTRAP=EXTRAP)
                END IF
              END IF
            END IF
          END IF
        ELSE
          IF (RNORM_PRESENT) THEN
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM, IBUDGET=I&
&BUDGET, CHAIN=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM, IBUDGET=I&
&BUDGET, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM, IBUDGET=I&
&BUDGET, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM, IBUDGET=I&
&BUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM, CHAIN=CHA&
&IN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM, CHAIN=CHA&
&IN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM, PMODE=PMO&
&DE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, RNORM=RNORM)
                END IF
              END IF
            END IF
          ELSE
            IF (IBUDGET_PRESENT) THEN
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, IBUDGET=IBUDGET, CHAIN&
&=CHAIN, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, IBUDGET=IBUDGET, CHAIN&
&=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, IBUDGET=IBUDGET, PMODE&
&=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, IBUDGET=IBUDGET)
                END IF
              END IF
            ELSE
              IF (CHAIN_PRESENT) THEN
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, CHAIN=CHAIN, PMODE=PMO&
&DE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, CHAIN=CHAIN)
                END IF
              ELSE
                IF (PMODE_PRESENT) THEN
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR, PMODE=PMODE)
                ELSE
                  CALL DELAUNAYSPARSEP(D=D, N=N, PTS=PTS, M=M, Q=Q, SIMPS=SIMPS, WEIGHTS=WEIGHTS, IERR=IERR)
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF
  END IF
END SUBROUTINE C_DELAUNAYSPARSEP


MODULE C_REAL_PRECISION
  IMPLICIT NONE


CONTAINS


  ! Getter and setter for R8.
  SUBROUTINE REAL_PRECISION_GET_R8(R8_LOCAL) BIND(C)
    USE REAL_PRECISION, ONLY: R8
    INTEGER :: R8_LOCAL
    R8_LOCAL = R8
  END SUBROUTINE REAL_PRECISION_GET_R8
END MODULE C_REAL_PRECISION


MODULE C_DELSPARSE_MOD
USE REAL_PRECISION
  IMPLICIT NONE


  INTERFACE
    SUBROUTINE DELAUNAYSPARSES(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, INTERP_IN, INTERP_OUT, EPS, EXTRAP, RNORM, IBUDGET, CHAIN)
      USE REAL_PRECISION , ONLY : R8
      INTEGER, INTENT(IN) :: D
      INTEGER, INTENT(IN) :: N
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(:,:) :: PTS
      INTEGER, INTENT(IN) :: M
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(:,:) :: Q
      INTEGER, INTENT(OUT), DIMENSION(:,:) :: SIMPS
      REAL(KIND=R8), INTENT(OUT), DIMENSION(:,:) :: WEIGHTS
      INTEGER, INTENT(OUT), DIMENSION(:) :: IERR
      REAL(KIND=R8), INTENT(IN), OPTIONAL, DIMENSION(:,:) :: INTERP_IN
      REAL(KIND=R8), INTENT(OUT), OPTIONAL, DIMENSION(:,:) :: INTERP_OUT
      REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
      REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
      REAL(KIND=R8), INTENT(OUT), OPTIONAL, DIMENSION(:) :: RNORM
      INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
      LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
    END SUBROUTINE DELAUNAYSPARSES
    SUBROUTINE DELAUNAYSPARSEP(D, N, PTS, M, Q, SIMPS, WEIGHTS, IERR, INTERP_IN, INTERP_OUT, EPS, EXTRAP, RNORM, IBUDGET, CHAIN, PM&
&ODE)
      ! Interface for parallel subroutine DELAUNAYSPARSEP.
      USE REAL_PRECISION , ONLY : R8
      INTEGER, INTENT(IN) :: D
      INTEGER, INTENT(IN) :: N
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(:,:) :: PTS
      INTEGER, INTENT(IN) :: M
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(:,:) :: Q
      INTEGER, INTENT(OUT), DIMENSION(:,:) :: SIMPS
      REAL(KIND=R8), INTENT(OUT), DIMENSION(:,:) :: WEIGHTS
      INTEGER, INTENT(OUT), DIMENSION(:) :: IERR
      REAL(KIND=R8), INTENT(IN), OPTIONAL, DIMENSION(:,:) :: INTERP_IN
      REAL(KIND=R8), INTENT(OUT), OPTIONAL, DIMENSION(:,:) :: INTERP_OUT
      REAL(KIND=R8), INTENT(IN), OPTIONAL :: EPS
      REAL(KIND=R8), INTENT(IN), OPTIONAL :: EXTRAP
      REAL(KIND=R8), INTENT(OUT), OPTIONAL, DIMENSION(:) :: RNORM
      INTEGER, INTENT(IN), OPTIONAL :: IBUDGET
      LOGICAL, INTENT(IN), OPTIONAL :: CHAIN
      INTEGER, INTENT(IN), OPTIONAL :: PMODE
    END SUBROUTINE DELAUNAYSPARSEP
    SUBROUTINE DWNNLS(W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE, IWORK, WORK)
      ! Interface for SLATEC subroutine DWNNLS.
      USE REAL_PRECISION , ONLY : R8
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(MDW,*) :: W
      INTEGER, INTENT(INOUT) :: MDW
      INTEGER, INTENT(INOUT) :: ME
      INTEGER, INTENT(INOUT) :: MA
      INTEGER, INTENT(INOUT) :: N
      INTEGER, INTENT(INOUT) :: L
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(*) :: PRGOPT
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(*) :: X
      REAL(KIND=R8), INTENT(INOUT) :: RNORM
      INTEGER, INTENT(INOUT) :: MODE
      INTEGER, INTENT(INOUT), DIMENSION(*) :: IWORK
      REAL(KIND=R8), INTENT(INOUT), DIMENSION(*) :: WORK
    END SUBROUTINE DWNNLS
  
  END INTERFACE

CONTAINS

END MODULE C_DELSPARSE_MOD

