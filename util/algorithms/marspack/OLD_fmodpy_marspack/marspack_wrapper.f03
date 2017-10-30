MODULE MARSPACK_WRAPPER
  USE ISO_C_BINDING, ONLY: C_LONG, C_DOUBLE
  
  IMPLICIT NONE

CONTAINS


  SUBROUTINE c_MARS( N, P, X, Y, W, NK, MI, LX , FM_0, FM, IM_0, IM ) BIND(c)
    !        &, SP_0, SP, DP_0, DP, MM_0, MM ) BIND(c)
    INTEGER(KIND=C_LONG) :: N
    INTEGER(KIND=C_LONG) :: P
    REAL(KIND=C_DOUBLE), DIMENSION(N,P) :: X
    REAL(KIND=C_DOUBLE), DIMENSION(N) :: Y
    REAL(KIND=C_DOUBLE), DIMENSION(N) :: W
    INTEGER(KIND=C_LONG) :: NK
    INTEGER(KIND=C_LONG) :: MI
    INTEGER(KIND=C_LONG), DIMENSION(P) :: LX
    INTEGER(KIND=C_LONG), INTENT(IN) :: FM_0
    REAL(KIND=C_DOUBLE), INTENT(OUT), DIMENSION(FM_0) :: FM
    INTEGER(KIND=C_LONG), INTENT(IN) :: IM_0
    INTEGER(KIND=C_LONG), INTENT(OUT), DIMENSION(IM_0) :: IM
    ! INTEGER(KIND=C_LONG), INTENT(IN) :: SP_0
    ! REAL(KIND=C_DOUBLE), DIMENSION(SP_0) :: SP
    ! INTEGER(KIND=C_LONG), INTENT(IN) :: DP_0
    ! REAL(KIND=C_DOUBLE), DIMENSION(DP_0) :: DP
    ! INTEGER(KIND=C_LONG), INTENT(IN) :: MM_0
    ! INTEGER(KIND=C_LONG), DIMENSION(MM_0) :: MM

    ! Local variables that need to be made for data type matching
    INTEGER :: TEMP_N
    INTEGER :: TEMP_P
    REAL, DIMENSION(N,P) :: TEMP_X
    REAL, DIMENSION(N) :: TEMP_Y
    REAL, DIMENSION(N) :: TEMP_W
    INTEGER :: TEMP_NK
    INTEGER :: TEMP_MI
    INTEGER, DIMENSION(P) :: TEMP_LX
    REAL, DIMENSION(3+NK*(5*MI+0+6)+2*P+0) :: TEMP_FM
    INTEGER, DIMENSION(21+NK*(3*MI+8)) :: TEMP_IM
    DOUBLE PRECISION, DIMENSION(&
         N*(MAX(NK+1,2)+3)+MAX(3*N+5*NK+P,2*P,4*N)+2*P+4*NK ) :: TEMP_SP
    DOUBLE PRECISION, DIMENSION(&
         MAX(N*NK,(NK+1)*(NK+1))+MAX((NK+2)*(0+3),4*NK)     ) :: TEMP_DP
    INTEGER, DIMENSION(N*P+2*MAX(MI,0)) :: TEMP_MM
    ! WARNING: Automatically generated the following copies to handle a
    ! type-mismatch. NumPy uses 64-bit (KIND=8) INTEGERs and REALs.
    ! Consider converting 'MARS' to (KIND=8) for better performance.
    TEMP_N = N
    TEMP_P = P
    TEMP_X = X
    TEMP_Y = Y
    TEMP_W = W
    TEMP_NK = NK
    TEMP_MI = MI
    TEMP_LX = LX
    TEMP_FM = FM
    TEMP_IM = IM
    ! TEMP_SP = SP
    ! TEMP_DP = DP
    ! TEMP_MM = MM
    CALL MARS(TEMP_N, TEMP_P, TEMP_X, TEMP_Y, TEMP_W, TEMP_NK, TEMP_MI&
         &, TEMP_LX, TEMP_FM, TEMP_IM, TEMP_SP, TEMP_DP, TEMP_MM)
    ! Copy back into the output containers
    FM = TEMP_FM
    IM = TEMP_IM
  END SUBROUTINE c_MARS


  SUBROUTINE c_FMOD( M, N, X_1, X, FM_0, FM, IM_0, IM, F_0, F, SP_0, S&
       &P ) BIND(c)
    INTEGER(KIND=C_LONG) :: M
    INTEGER(KIND=C_LONG) :: N
    INTEGER(KIND=C_LONG), INTENT(IN) :: X_1
    REAL(KIND=C_DOUBLE), DIMENSION(N,X_1) :: X
    INTEGER(KIND=C_LONG), INTENT(IN) :: FM_0
    REAL(KIND=C_DOUBLE), DIMENSION(FM_0) :: FM
    INTEGER(KIND=C_LONG), INTENT(IN) :: IM_0
    INTEGER(KIND=C_LONG), DIMENSION(IM_0) :: IM
    INTEGER(KIND=C_LONG), INTENT(IN) :: F_0
    REAL(KIND=C_DOUBLE), INTENT(OUT), DIMENSION(F_0) :: F
    INTEGER(KIND=C_LONG), INTENT(IN) :: SP_0
    REAL(KIND=C_DOUBLE), DIMENSION(SP_0,2) :: SP
    ! Local variables that need to be made for data type matching
    INTEGER :: TEMP_M
    INTEGER :: TEMP_N
    REAL, DIMENSION(N,X_1) :: TEMP_X
    REAL, DIMENSION(FM_0) :: TEMP_FM
    INTEGER, DIMENSION(IM_0) :: TEMP_IM
    REAL, DIMENSION(F_0) :: TEMP_F
    REAL, DIMENSION(SP_0,2) :: TEMP_SP
    ! WARNING: Automatically generated the following copies to handle a
    ! type-mismatch. NumPy uses 64-bit (KIND=8) INTEGERs and REALs.
    ! Consider converting 'FMOD' to (KIND=8) for better performance.
    TEMP_M = M
    TEMP_N = N
    TEMP_X = X
    TEMP_FM = FM
    TEMP_IM = IM
    TEMP_F = F
    TEMP_SP = SP
    CALL FMOD(TEMP_M, TEMP_N, TEMP_X, TEMP_FM, TEMP_IM, TEMP_F, TEMP_S&
         &P)
    F = TEMP_F
  END SUBROUTINE c_FMOD
END MODULE MARSPACK_WRAPPER
