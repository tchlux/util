PROGRAM GET_SIZE
USE ISO_FORTRAN_ENV , ONLY : REAL64
  REAL(KIND=REAL64) :: PTS
  WRITE (*,*) SIZEOF(PTS)
END PROGRAM GET_SIZE
