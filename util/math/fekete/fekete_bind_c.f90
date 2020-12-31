! This automatically generated Fortran wrapper file allows codes
! written in Fortran to be called directly from C and translates all
! C-style arguments into expected Fortran-style arguments (with
! assumed size, local type declarations, etc.).


SUBROUTINE C_FEKETE_INDICES(VM_T_DIM_1, VM_T_DIM_2, VM_T, INDS_DIM_1, INDS, INFO) BIND(C)
  USE ISO_FORTRAN_ENV , ONLY : REAL64
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: VM_T_DIM_1
  INTEGER, INTENT(IN) :: VM_T_DIM_2
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(VM_T_DIM_1,VM_T_DIM_2) :: VM_T
  INTEGER, INTENT(IN) :: INDS_DIM_1
  INTEGER, INTENT(OUT), DIMENSION(INDS_DIM_1) :: INDS
  INTEGER, INTENT(OUT) :: INFO

  INTERFACE
    SUBROUTINE FEKETE_INDICES(VM_T, INDS, INFO)
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
      !              {{f_1, ..., f_k}} and points {{p_1, ..., p_n}}, we have
      !                   VM_T_{{i,j}} = f_i(p_j)
      !
      ! OUTPUT:
      !    INDS  --  The sequence of indices (of the points {{p_i}}) that
      !              are most useful in determining an interpolant. E.g.
      !              to make a Vandermonde matrix square that was
      !              previously long (more points than functions), use the
      !              points at indices JPVT(1:K) as the Fekete points.
      !
      USE ISO_FORTRAN_ENV , ONLY : REAL64
      REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: VM_T
      INTEGER, INTENT(OUT), DIMENSION(SIZE(VM_T,2)) :: INDS
      INTEGER, INTENT(OUT) :: INFO
    END SUBROUTINE FEKETE_INDICES
  END INTERFACE

  CALL FEKETE_INDICES(VM_T, INDS, INFO)
END SUBROUTINE C_FEKETE_INDICES


SUBROUTINE C_EVALUATE_VANDERMONDE(POINTS_DIM_1, POINTS_DIM_2, POINTS, DEGREES_DIM_1, DEGREES, INDICES_DIM_1, INDICES_DIM_2, INDICES&
&, VANDERMONDE_DIM_1, VANDERMONDE_DIM_2, VANDERMONDE) BIND(C)
  USE ISO_FORTRAN_ENV , ONLY : REAL64 , INT64
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: POINTS_DIM_1
  INTEGER, INTENT(IN) :: POINTS_DIM_2
  REAL(KIND=REAL64), INTENT(IN), DIMENSION(POINTS_DIM_1,POINTS_DIM_2) :: POINTS
  INTEGER, INTENT(IN) :: DEGREES_DIM_1
  INTEGER(KIND=INT64), INTENT(IN), DIMENSION(DEGREES_DIM_1) :: DEGREES
  INTEGER, INTENT(IN) :: INDICES_DIM_1
  INTEGER, INTENT(IN) :: INDICES_DIM_2
  INTEGER(KIND=INT64), INTENT(IN), DIMENSION(INDICES_DIM_1,INDICES_DIM_2) :: INDICES
  INTEGER, INTENT(IN) :: VANDERMONDE_DIM_1
  INTEGER, INTENT(IN) :: VANDERMONDE_DIM_2
  REAL(KIND=REAL64), INTENT(OUT), DIMENSION(VANDERMONDE_DIM_1,VANDERMONDE_DIM_2) :: VANDERMONDE

  INTERFACE
    SUBROUTINE EVALUATE_VANDERMONDE(POINTS, DEGREES, INDICES, VANDERMONDE)
      ! Evaluate polynomials at points to construct a Vandermonde matrix,
      ! the lengths are the number of terms in each polynomial, the indices
      ! are the active index terms in each polynomial.
      !
      ! This should be super fast, aka even a matrix with 100M elements should
      ! take less than one second. For some reason it appears to slow down when
      ! (SIZE(POINTS, 1) == 1). That is the reason for including schedule(dynamic),
      ! which does not have a huge negative effect otherwise.
      USE ISO_FORTRAN_ENV , ONLY : REAL64 , INT64
      REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: POINTS
      INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: DEGREES
      INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:,:) :: INDICES
      REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(DEGREES),SIZE(POINTS,2)) :: VANDERMONDE
    END SUBROUTINE EVALUATE_VANDERMONDE
  END INTERFACE

  CALL EVALUATE_VANDERMONDE(POINTS, DEGREES, INDICES, VANDERMONDE)
END SUBROUTINE C_EVALUATE_VANDERMONDE

