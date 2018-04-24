SUBROUTINE EVAL_MESH(int_prods, app_int_prods, mesh)
  USE ISO_FORTRAN_ENV, ONLY : REAL64
  IMPLICIT NONE
  ! Subroutine for evaluating a standard voronoi mesh at a set of
  ! points given the products have already been calculated.
  REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: int_prods
  REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: app_int_prods
  REAL(KIND=REAL64), DIMENSION(SIZE(app_int_prods,1),SIZE(int_prods,1)), &
       INTENT(OUT) :: mesh
  ! Local variables
  REAL(KIND=REAL64) :: max_ratio, ratio, numerator, denominator
  INTEGER :: y, xi, c
  ! Loop for computing the mesh
  approx_points : DO y = 1, SIZE(app_int_prods,1)
     !$omp parallel do private(xi,c,max_ratio,ratio,numerator,denominator)
     interp_points : DO xi = 1, SIZE(int_prods,1)
        max_ratio = 0
        cell_centers : DO c = 1, SIZE(int_prods,1)
           ! Cycle when we are considering a cell's center and itself
           IF (c .EQ. xi) CYCLE
           ! Calculate the ratio 
           numerator = app_int_prods(y,c) - app_int_prods(y,xi) - &
                int_prods(xi,c) + int_prods(xi,xi)
           denominator = int_prods(c,c) - int_prods(c,xi) - &
                int_prods(xi,c) + int_prods(xi,xi)
           valid_denom : IF (ABS(denominator) .GT. EPSILON(ratio)) THEN
              ratio = numerator / denominator
              track_max : IF (ratio > max_ratio) THEN
                 max_ratio = ratio
              END IF track_max
           END IF valid_denom
        END DO cell_centers
        ! Store the ratio for this cell in the mesh
        mesh(y,xi) = MAX(REAL(1,KIND=REAL64) - max_ratio, REAL(0,KIND=REAL64))
     END DO interp_points
     !$omp end parallel do
     ! Normalize the mesh across the interpolation points (to be convex combinations)
     mesh(y,:) = mesh(y,:) / SUM(mesh(y,:))
  END DO approx_points
END SUBROUTINE EVAL_MESH
