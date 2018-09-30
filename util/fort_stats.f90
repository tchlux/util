MODULE FORT_STATS
CONTAINS
  ! The serial version of the basis vector update step.
  SUBROUTINE BASIS_UPDATE(count, basis, vec, vec_length)
    USE ISO_FORTRAN_ENV, ONLY: REAL64
    IMPLICIT NONE
    REAL(KIND=REAL64), INTENT(IN)                  :: count
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:) :: basis, vec
    REAL(KIND=REAL64), INTENT(OUT)                 :: vec_length
    ! Local variables
    INTEGER :: i
    REAL(KIND=REAL64) :: bv_dot, denominator

    ! Compute the length of the update vector, if it is small enough, break.
    vec_length = NORM2(vec)
    IF (vec_length .LE. 0) RETURN

    ! Flip the vector if it is not facing the same way as the basis.
    bv_dot = DOT_PRODUCT(basis, vec)
    IF (bv_dot .LT. 0) THEN
       vec(:) = -vec(:)
    END IF

    basis(:) = basis(:) + (vec(:) - basis(:)) / count

    ! Update the vector by removing the component in the direction of
    ! the basis vector. Skip if the basis vector is too short.
    denominator = NORM2(basis)**2
    IF (denominator .LE. SQRT(SQRT(EPSILON(bv_dot)))) RETURN
    vec = vec - DOT_PRODUCT(vec, basis) * (basis / denominator)
  END SUBROUTINE BASIS_UPDATE
END MODULE FORT_STATS
