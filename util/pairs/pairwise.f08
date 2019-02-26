SUBROUTINE FUll(PTS, DISTS)
  ! Compute the pairwise distance between points in "PTS(N,D)" and
  ! return in array "DISTS" of length N(N-1)//2
  USE ISO_FORTRAN_ENV, ONLY : REAL64
  IMPLICIT NONE
  ! Arguments
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: PTS
  REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(PTS,1)*(SIZE(PTS,1) - 1) / 2) :: DISTS
  ! Local variables
  INTEGER :: I, P1, P2
  ! Compute the pairwise distances.
  !$omp parallel do private(P1,P2)
  DO I = 0, SIZE(DISTS)-1
     P1 = FLOOR(SQRT(0.25_REAL64 + 2._REAL64 * REAL(I,REAL64)) + 0.5_REAL64)
     P2 = FLOOR(REAL(I,REAL64) - REAL(P1,REAL64) * (REAL(P1,REAL64) - 1._REAL64) / 2._REAL64)
     DISTS(I+1) = NORM2(PTS(P1+1,:) - PTS(P2+1,:))
  END DO
END SUBROUTINE FULL


SUBROUTINE HALF(FROM_PTS, TO_PTS, DISTS)
  ! Compute the pairwise distance between points in "PTS(N,D)" and
  ! return in array "DISTS" of length N(N-1)//2
  USE ISO_FORTRAN_ENV, ONLY : REAL64
  IMPLICIT NONE
  ! Arguments
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: FROM_PTS, TO_PTS
  REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(FROM_PTS,2),SIZE(TO_PTS,2)) :: DISTS
  ! Local variables
  REAL(KIND=REAL64), DIMENSION(SIZE(FROM_PTS,2)) :: FROM_SQ
  REAL(KIND=REAL64), DIMENSION(SIZE(TO_PTS,2))   :: TO_SQ
  REAL(KIND=REAL64), DIMENSION(SIZE(FROM_PTS,2),SIZE(TO_PTS,2)) :: DOTS
  INTEGER :: P1, P2
  ! Compute the squared values.
  FROM_SQ(:) = SUM(FROM_PTS(:,:)**2, DIM=1)
  TO_SQ(:) = SUM(TO_PTS(:,:)**2, DIM=1)
  DOTS = 2 * MATMUL(TRANSPOSE(FROM_PTS), TO_PTS)
  ! Compute the pairwise sums of the squared values.
  FORALL (P2=1:SIZE(TO_PTS,2))
     DISTS(:,P2) = FROM_SQ(:) + TO_SQ(P2)
  END FORALL
  ! Compute the pairwise distances.     
  DISTS(:,:) = SQRT(DISTS(:,:) - DOTS(:,:))

END SUBROUTINE HALF
