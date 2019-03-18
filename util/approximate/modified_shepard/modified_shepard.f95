MODULE MODIFIED_SHEPARD
  USE ISO_FORTRAN_ENV, ONLY : REAL64
  IMPLICIT NONE

  PRIVATE
  PUBLIC SHEPMOD, SHEPMODVAL

CONTAINS

  SUBROUTINE SHEPMOD ( M, N, X, RW, IER )
    ! This subroutine computes a set of parameters defining a function that 
    ! interpolates N interpolation nodes at scattered X(i) in M dimensions. The 
    ! interpolant may be evaluated at an arbitrary point by the function SHEPMODVAL.
    ! The interpolation scheme is a modified Shepard method, and can
    ! use either the Shepard weighted fit or standarad shepard method.
    !
    ! Input parameters:  
    !   M is the dimension of the data.
    !   N is the number of nodes.
    !   X(M, N) contains the coordinates of the nodes.
    ! Output parameters:
    !   RW(N) is an array containing the radius of influence for each point.
    !   IER 
    !   = 0, if no errors were encountered.
    !   = 1, if N is too small relative to M.
    INTEGER, INTENT(IN) :: M, N
    REAL(KIND=REAL64), DIMENSION(M,N), INTENT(IN) :: X
    REAL(KIND=REAL64), DIMENSION(N), INTENT(OUT) :: RW
    INTEGER, INTENT(OUT) :: IER

    ! Local variables.         
    INTEGER :: I, J, K  ! loop control variables.
    INTEGER :: NP    ! number of nodes used for the local fit.
    REAL(KIND=REAL64), DIMENSION(N-1) :: DIST ! array that stores all the distances. 

    IER = 0
    ! Check the validity of M and N.    
    IF ( ( M < 1 ) .OR. ( N <= (M + 1) ) ) THEN
       IER = 1
       RETURN
    END IF
    ! Set the number of neighbors to define radius.
    NP = M + 1
    ! Calculate the radius of influence "RW" for all interpolation nodes.
    DO K = 1, N
       J = 0
       DO I = 1, N
          IF (I /= K) THEN
             ! Compute the squared distance to each *other* node.
             J = J + 1
             DIST(J) = DOT_PRODUCT( X(:, I) - X(:, K), X(:, I) - X(:, K) )
          END IF
       END DO
       ! Sort the squared distances to all other nodes.
       CALL SORT( DIST, NP, N - 1 )
       ! Only bother taking the square root of the choice distance.
       RW(K) = SQRT(DIST(NP))
    END DO
    RETURN
  END SUBROUTINE SHEPMOD

  SUBROUTINE SORT( DIST, NUM, LENGTH )
    ! The subroutine SORT sorts the real array DIST of length LENGTH in ascending 
    ! order for the smallest NUM elements.
    ! Local variables.
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:)::DIST
    INTEGER :: I, ITEMP, J, LENGTH, NUM
    REAL(KIND=REAL64) :: TEMP

    DO I = 2, LENGTH
       DO J = 1, MIN( I-1, NUM )
          IF ( DIST(I) < DIST(J) ) THEN
             TEMP = DIST(I)
             DIST(J+1:I) = DIST(J:I-1) 
             DIST(J) = TEMP
             EXIT 
          END IF
       END DO
    END DO
  END SUBROUTINE SORT


  SUBROUTINE SHEPMODVAL( XP, M, N, X, RW, WTS, IER )
    ! SHEPMODVAL returns the modified Shepard approximation at the
    ! point XP, using the points computed by SHEPMOD.
    !
    ! Input parameters:
    !  XP is the point at which the modified Shepard interpolant 
    !     function approximation function is to be evaluated.
    !  M is the dimension of the data.
    !  N is the number of interpolation points.
    !  X contains the interpolation nodes, by column.
    !  RW contains the radius of influence about each interpolation 
    !    node returned by SHEPMOD.
    ! Output parameter:
    !  WTS contains the weight associated with each interpolation node.
    !  IER 
    !   = 0, normal returns.
    !   = 1, if the point XP is outside the radius of influence RW(i) for all 
    !        nodes, in which case SHEPMODVAL is computed using the original Shepard 
    !        algorithm with the M+1 closest points.
    INTEGER, INTENT(IN)  :: M, N
    REAL(KIND=REAL64), DIMENSION(M), INTENT(IN) :: XP 
    REAL(KIND=REAL64), DIMENSION(M, N), INTENT(IN) :: X
    REAL(KIND=REAL64), DIMENSION(N), INTENT(IN) :: RW
    REAL(KIND=REAL64), DIMENSION(N), INTENT(OUT) :: WTS ! weight values. 
    INTEGER, INTENT(OUT):: IER

    ! Local variables.
    INTEGER :: I   ! number of nodal functions used in Shepard approximation.
    INTEGER :: J   ! temporary integer variable.
    INTEGER :: K   ! loop control variable.
    REAL(KIND=REAL64) :: DIST ! distance between XP and the interpolation nodes.
    REAL(KIND=REAL64) :: TEMP ! temporary real variable.
    REAL(KIND=REAL64) :: W    ! weight value, depends on D.    
    REAL(KIND=REAL64), DIMENSION(M+1) :: D ! weights in the original Shepard scheme.
                                           ! local approximations with positive weights.  
    INTEGER, DIMENSION(M+1) :: D_IDX ! Indices of points in D

    IER = 0
    I = 0
    J = 1
    ! Initialize all weights to be zero.
    WTS(:) = 0.0_REAL64
    ! Be prepared to interpolate using the nearest (M+1) points if
    ! there are no points whose radius of influence contain XP.
    D_IDX(:) = 0
    D(1:M+1) = 0.0_REAL64
    TEMP = 0.0_REAL64
    ! Cycle through interpolation nodes, checking distances to XP.
    DO K = 1, N
       DIST = SQRT( DOT_PRODUCT( XP(:) - X(:, K), XP(:) - X(:, K) ) )
       ! If XP is within the radius of influence of X(K)
       IF ( DIST < RW(K) ) THEN
          ! Set the weight for this point to 1 if it is exact.
          IF ( RW(K) - DIST == RW(K) ) THEN
             WTS(:) = 0
             WTS(K) = 1.0_REAL64
             RETURN
          END IF
          ! Otherwise store the weight based on standard scheme.
          ! WEIGHT = ((RADIUS - DIST_TO_DATA) / (RADIUS * DIST_TO_DATA))^2
          WTS(K) = ((RW(K) - DIST) / (RW(K) * DIST))**2
          ! Track the number of points inside the radius of influence.
          I = I + 1
       ! Otherwise, XP is *not* inside the radius of influence of X(K).
       ELSE
          ! Compute the weight as the inverse distance (will square later).
          W = 1.0_REAL64 / DIST
          ! If this point is closer than one of the elements of D, replace.
          IF ( W > TEMP ) THEN
             D(J) = W
             D_IDX(J) = K
             ! Get the next smallest influencer currently in D to be removed.
             J = MINLOC( D(1:M+1), DIM = 1 )
             TEMP = D(J)
          END IF
       END IF
    END DO
    ! I = 0 iff the point XP is not within the radius RW(K) of X(:,K) for all K; 
    ! IER = 1. Revert to using the (M+1) points with inverse distance
    ! squared weighting (original shepard method).
    IF ( I == 0 ) THEN
       D = D**2
       ! If all (M+1) points have been used, increment J to include all.
       IF (D(MINLOC(D(1:M+1), DIM = 1 )) .GT. 0.0_REAL64) J = M+1
       ! Assume weights D(1:J) are all nonzero, use Shepard method.
       WTS(D_IDX(1:J)) = D(:) / SUM( D(:) )
       IER = 1 
       RETURN
    ELSE
       WTS(1:N) = WTS(1:N) / SUM( WTS(1:N) )
    END IF

    RETURN
  END SUBROUTINE SHEPMODVAL

END MODULE MODIFIED_SHEPARD
