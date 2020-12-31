
! This function identifies the set of row vectors in a Vandermonde
! matrix that should be kept. It outputs the row indices, however
! it expects the Vandermonde matrix to be provided in its transpose.
! Inspiration for greedy selection technique from:
! 
!    Bos, Len, et al. "Computing multivariate Fekete and Leja points
!    by numerical linear algebra." SIAM Journal on Numerical Analysis
!    48.5 (2010): 1984-1999. DOI 10.1137/090779024
! 
! USES LAPACK:
!    DGEQP3  --  (double precision QR decomposition with column pivoting)
! 
! INPUT:
!    VM_T  --  The transposed Vandermonde matrix. Given functions 
!              {f_1, ..., f_k} and points {p_1, ..., p_n}, we have
!                   VM_T_{i,j} = f_i(p_j)
! 
! OUTPUT:
!    INDS  --  The sequence of indices (of the points {p_i}) that
!              are most useful in determining an interpolant. E.g.
!              to make a Vandermonde matrix square that was
!              previously long (more points than functions), use the
!              points at indices JPVT(1:K) as the Fekete points.
!              
SUBROUTINE FEKETE_INDICES(VM_T, INDS, INFO)
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  EXTERNAL :: DGEQP3
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: VM_T
  INTEGER, INTENT(OUT), DIMENSION(SIZE(VM_T,2)) :: INDS
  INTEGER, INTENT(OUT)                          :: INFO
  ! Local variables.
  REAL(KIND=REAL64), DIMENSION(SIZE(VM_T,1))   :: TAU
  REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: WORK
  ! Set all columns to be "free" (movable).
  INDS(:) = 0
  ! Do a work array size query (store size in TAU).
  CALL DGEQP3(SIZE(VM_T,1), SIZE(VM_T,2), VM_T, SIZE(VM_T,1), &
       INDS, TAU, TAU, -1, INFO)
  ! Allocate the work array.
  ALLOCATE( WORK(INT(TAU(1))) )
  ! Call the QR (with column pivoting) function.
  CALL DGEQP3(SIZE(VM_T,1), SIZE(VM_T,2), VM_T, SIZE(VM_T,1), &
       INDS, TAU, WORK, SIZE(WORK), INFO)
  DEALLOCATE( WORK )
END SUBROUTINE FEKETE_INDICES


! Evaluate polynomials at points to construct a Vandermonde matrix,
! the lengths are the number of terms in each polynomial, the indices
! are the active index terms in each polynomial.
! 
! This should be super fast, aka even a matrix with 100M elements should
! take less than one second. For some reason it appears to slow down when
! (SIZE(POINTS, 1) == 1). That is the reason for including schedule(dynamic),
! which does not have a huge negative effect otherwise.
SUBROUTINE EVALUATE_VANDERMONDE(POINTS, DEGREES, INDICES, VANDERMONDE)
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64
  REAL(KIND=REAL64),   INTENT(IN),  DIMENSION(:,:) :: POINTS
  INTEGER(KIND=INT64), INTENT(IN),  DIMENSION(:)   :: DEGREES
  INTEGER(KIND=INT64), INTENT(IN),  DIMENSION(:,:) :: INDICES
  REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(SIZE(DEGREES),SIZE(POINTS,2)) :: VANDERMONDE
  ! Local variables.
  INTEGER :: I, J
  ! Compute the Vandermonde matrix values.
  !$omp parallel do schedule(dynamic) private(i,j)
  pts : DO I = 1,SIZE(POINTS,2)
     funcs : DO J = 1,SIZE(DEGREES)
        ! Assign the value of the Vandermonde matrix at this point.
        IF (DEGREES(J) .EQ. 0_INT64) THEN
           VANDERMONDE(J,I) = 1.0_REAL64
        ELSE
           VANDERMONDE(J,I) = PRODUCT(POINTS(INDICES(1:DEGREES(J),J), I))
        END IF
     END DO funcs
  END DO pts
  !$omp end parallel do
END SUBROUTINE EVALUATE_VANDERMONDE
