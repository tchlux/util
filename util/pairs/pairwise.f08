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
  REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(FROM_PTS,1)*SIZE(TO_PTS,1)) :: DISTS
  ! Local variables
  INTEGER :: I, P1, P2
  ! Compute the pairwise distances.
  !$omp parallel do private(P1,P2)
  DO I = 0, SIZE(DISTS)-1
     P1 = I / SIZE(TO_PTS,1)
     P2 = I - P1*SIZE(TO_PTS,1)
     DISTS(I+1) = NORM2(FROM_PTS(P1+1,:) - TO_PTS(P2+1,:))
  END DO
END SUBROUTINE HALF
