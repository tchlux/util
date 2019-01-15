
SUBROUTINE WEIGHT(PTS, PT, WTS)
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  ! Inputs and outputs
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: PTS
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(SIZE(PTS,1)) :: PT
  REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(PTS,2)) :: WTS
  ! Local variables
  REAL(KIND=REAL64), DIMENSION(SIZE(WTS)) :: DISTS
  REAL(KIND=REAL64) :: MAX_DIST, MIN_DIST, EPS
  INTEGER :: I
  !$omp parallel do
  DO I = 1, SIZE(PTS, 2)
     DISTS(I) = NORM2(PTS(:,I) - PT(:))
  END DO
  ! Compute the epsilon (for checking "zero" values).
  EPS = SQRT(EPSILON(WTS(1)))
  ! Compute the maximum and minimum distance to points.
  MIN_DIST = MINVAL(DISTS)
  MAX_DIST = MAXVAL(DISTS)
  ! Make sure that the maximum distance is positive.
  IF (MAX_DIST .LT. EPS) THEN
     ! Set all weights to be equal.
     WTS(:) = 1._REAL64
  ELSE IF (MIN_DIST .LT. EPS) THEN
     ! Share the weights equally across those close values.
     WTS(:) = 0._REAL64
     WHERE (DISTS .LT. EPS) WTS = 1._REAL64
  ELSE
     ! Compute the weights based on inverse distances.
     WTS(:) = 1._REAL64 / DISTS(:)**2
  END IF
  ! Convexify the weights.
  WTS = WTS / SUM(WTS)
END SUBROUTINE WEIGHT
