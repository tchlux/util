PROGRAM GET_SIZE
USE ISO_FORTRAN_ENV , ONLY : REAL64 , INT64
  INTEGER(KIND=INT64) :: DEGREES
  REAL(KIND=REAL64) :: POINTS
  WRITE (*,*) SIZEOF(DEGREES)
  WRITE (*,*) SIZEOF(POINTS)
END PROGRAM GET_SIZE