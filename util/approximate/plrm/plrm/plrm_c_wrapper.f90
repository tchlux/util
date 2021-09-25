! This automatically generated Fortran wrapper file allows codes
! written in Fortran to be called directly from C and translates all
! C-style arguments into expected Fortran-style arguments (with
! assumed size, local type declarations, etc.).


MODULE C_MATRIX_MULTIPLICATION
USE ISO_FORTRAN_ENV , ONLY : RT => REAL32
  IMPLICIT NONE


CONTAINS


  
  SUBROUTINE C_GEMM(OP_A, OP_B, OUT_ROWS, OUT_COLS, INNER_DIM, AB_MULT, A_DIM_1, A_DIM_2, A, A_ROWS, B_DIM_1, B_DIM_2, B, B_ROWS, C&
&_MULT, C_DIM_1, C_DIM_2, C, C_ROWS) BIND(C)
    USE MATRIX_MULTIPLICATION, ONLY: GEMM
    IMPLICIT NONE
    CHARACTER, INTENT(IN) :: OP_A
    CHARACTER, INTENT(IN) :: OP_B
    INTEGER, INTENT(IN) :: OUT_ROWS
    INTEGER, INTENT(IN) :: OUT_COLS
    INTEGER, INTENT(IN) :: INNER_DIM
    REAL(KIND=RT), INTENT(IN) :: AB_MULT
    INTEGER, INTENT(IN) :: A_DIM_1
    INTEGER, INTENT(IN) :: A_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(A_DIM_1,A_DIM_2) :: A
    INTEGER, INTENT(IN) :: A_ROWS
    INTEGER, INTENT(IN) :: B_DIM_1
    INTEGER, INTENT(IN) :: B_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(B_DIM_1,B_DIM_2) :: B
    INTEGER, INTENT(IN) :: B_ROWS
    REAL(KIND=RT), INTENT(IN) :: C_MULT
    INTEGER, INTENT(IN) :: C_DIM_1
    INTEGER, INTENT(IN) :: C_DIM_2
    REAL(KIND=RT), INTENT(OUT), DIMENSION(C_DIM_1,C_DIM_2) :: C
    INTEGER, INTENT(IN) :: C_ROWS
  
    CALL GEMM(OP_A, OP_B, OUT_ROWS, OUT_COLS, INNER_DIM, AB_MULT, A, A_ROWS, B, B_ROWS, C_MULT, C, C_ROWS)
  END SUBROUTINE C_GEMM
  
END MODULE C_MATRIX_MULTIPLICATION


MODULE C_PLRM
USE ISO_FORTRAN_ENV , ONLY : RT => REAL32
USE MATRIX_MULTIPLICATION , ONLY : GEMM
  USE ISO_FORTRAN_ENV, ONLY: INT64
  IMPLICIT NONE


CONTAINS


  ! Getter and setter for DISCONTINUITY.
  SUBROUTINE PLRM_GET_DISCONTINUITY(DISCONTINUITY_LOCAL) BIND(C)
    USE PLRM, ONLY: DISCONTINUITY
    REAL(KIND=RT) :: DISCONTINUITY_LOCAL
    DISCONTINUITY_LOCAL = DISCONTINUITY
  END SUBROUTINE PLRM_GET_DISCONTINUITY

  ! Getter and setter for INITIAL_SHIFT_RANGE.
  SUBROUTINE PLRM_GET_INITIAL_SHIFT_RANGE(INITIAL_SHIFT_RANGE_LOCAL) BIND(C)
    USE PLRM, ONLY: INITIAL_SHIFT_RANGE
    REAL(KIND=RT) :: INITIAL_SHIFT_RANGE_LOCAL
    INITIAL_SHIFT_RANGE_LOCAL = INITIAL_SHIFT_RANGE
  END SUBROUTINE PLRM_GET_INITIAL_SHIFT_RANGE

  ! Getter and setter for INITIAL_FLEX.
  SUBROUTINE PLRM_GET_INITIAL_FLEX(INITIAL_FLEX_LOCAL) BIND(C)
    USE PLRM, ONLY: INITIAL_FLEX
    REAL(KIND=RT) :: INITIAL_FLEX_LOCAL
    INITIAL_FLEX_LOCAL = INITIAL_FLEX
  END SUBROUTINE PLRM_GET_INITIAL_FLEX

  ! Getter and setter for INITIAL_OUTPUT_SCALE.
  SUBROUTINE PLRM_GET_INITIAL_OUTPUT_SCALE(INITIAL_OUTPUT_SCALE_LOCAL) BIND(C)
    USE PLRM, ONLY: INITIAL_OUTPUT_SCALE
    REAL(KIND=RT) :: INITIAL_OUTPUT_SCALE_LOCAL
    INITIAL_OUTPUT_SCALE_LOCAL = INITIAL_OUTPUT_SCALE
  END SUBROUTINE PLRM_GET_INITIAL_OUTPUT_SCALE

  ! Getter and setter for INITIAL_STEP.
  SUBROUTINE PLRM_GET_INITIAL_STEP(INITIAL_STEP_LOCAL) BIND(C)
    USE PLRM, ONLY: INITIAL_STEP
    REAL(KIND=RT) :: INITIAL_STEP_LOCAL
    INITIAL_STEP_LOCAL = INITIAL_STEP
  END SUBROUTINE PLRM_GET_INITIAL_STEP

  ! Getter and setter for INITIAL_STEP_MEAN_CHANGE.
  SUBROUTINE PLRM_GET_INITIAL_STEP_MEAN_CHANGE(INITIAL_STEP_MEAN_CHANGE_LOCAL) BIND(C)
    USE PLRM, ONLY: INITIAL_STEP_MEAN_CHANGE
    REAL(KIND=RT) :: INITIAL_STEP_MEAN_CHANGE_LOCAL
    INITIAL_STEP_MEAN_CHANGE_LOCAL = INITIAL_STEP_MEAN_CHANGE
  END SUBROUTINE PLRM_GET_INITIAL_STEP_MEAN_CHANGE

  ! Getter and setter for INITIAL_STEP_CURV_CHANGE.
  SUBROUTINE PLRM_GET_INITIAL_STEP_CURV_CHANGE(INITIAL_STEP_CURV_CHANGE_LOCAL) BIND(C)
    USE PLRM, ONLY: INITIAL_STEP_CURV_CHANGE
    REAL(KIND=RT) :: INITIAL_STEP_CURV_CHANGE_LOCAL
    INITIAL_STEP_CURV_CHANGE_LOCAL = INITIAL_STEP_CURV_CHANGE
  END SUBROUTINE PLRM_GET_INITIAL_STEP_CURV_CHANGE

  ! Getter and setter for STEP_GROWTH_RATE.
  SUBROUTINE PLRM_GET_STEP_GROWTH_RATE(STEP_GROWTH_RATE_LOCAL) BIND(C)
    USE PLRM, ONLY: STEP_GROWTH_RATE
    REAL(KIND=RT) :: STEP_GROWTH_RATE_LOCAL
    STEP_GROWTH_RATE_LOCAL = STEP_GROWTH_RATE
  END SUBROUTINE PLRM_GET_STEP_GROWTH_RATE

  ! Getter and setter for STEP_SHRINK_RATE.
  SUBROUTINE PLRM_GET_STEP_SHRINK_RATE(STEP_SHRINK_RATE_LOCAL) BIND(C)
    USE PLRM, ONLY: STEP_SHRINK_RATE
    REAL(KIND=RT) :: STEP_SHRINK_RATE_LOCAL
    STEP_SHRINK_RATE_LOCAL = STEP_SHRINK_RATE
  END SUBROUTINE PLRM_GET_STEP_SHRINK_RATE

  ! Getter and setter for MIN_STEPS_TO_STABILITY.
  SUBROUTINE PLRM_GET_MIN_STEPS_TO_STABILITY(MIN_STEPS_TO_STABILITY_LOCAL) BIND(C)
    USE PLRM, ONLY: MIN_STEPS_TO_STABILITY
    INTEGER :: MIN_STEPS_TO_STABILITY_LOCAL
    MIN_STEPS_TO_STABILITY_LOCAL = MIN_STEPS_TO_STABILITY
  END SUBROUTINE PLRM_GET_MIN_STEPS_TO_STABILITY

  ! Getter and setter for MDI.
  SUBROUTINE PLRM_GET_MDI(MDI_LOCAL) BIND(C)
    USE PLRM, ONLY: MDI
    INTEGER :: MDI_LOCAL
    MDI_LOCAL = MDI
  END SUBROUTINE PLRM_GET_MDI
  SUBROUTINE PLRM_SET_MDI(MDI_LOCAL) BIND(C)
    USE PLRM, ONLY: MDI
    INTEGER :: MDI_LOCAL
    MDI = MDI_LOCAL
  END SUBROUTINE PLRM_SET_MDI

  ! Getter and setter for MDS.
  SUBROUTINE PLRM_GET_MDS(MDS_LOCAL) BIND(C)
    USE PLRM, ONLY: MDS
    INTEGER :: MDS_LOCAL
    MDS_LOCAL = MDS
  END SUBROUTINE PLRM_GET_MDS
  SUBROUTINE PLRM_SET_MDS(MDS_LOCAL) BIND(C)
    USE PLRM, ONLY: MDS
    INTEGER :: MDS_LOCAL
    MDS = MDS_LOCAL
  END SUBROUTINE PLRM_SET_MDS

  ! Getter and setter for MNS.
  SUBROUTINE PLRM_GET_MNS(MNS_LOCAL) BIND(C)
    USE PLRM, ONLY: MNS
    INTEGER :: MNS_LOCAL
    MNS_LOCAL = MNS
  END SUBROUTINE PLRM_GET_MNS
  SUBROUTINE PLRM_SET_MNS(MNS_LOCAL) BIND(C)
    USE PLRM, ONLY: MNS
    INTEGER :: MNS_LOCAL
    MNS = MNS_LOCAL
  END SUBROUTINE PLRM_SET_MNS

  ! Getter and setter for MDO.
  SUBROUTINE PLRM_GET_MDO(MDO_LOCAL) BIND(C)
    USE PLRM, ONLY: MDO
    INTEGER :: MDO_LOCAL
    MDO_LOCAL = MDO
  END SUBROUTINE PLRM_GET_MDO
  SUBROUTINE PLRM_SET_MDO(MDO_LOCAL) BIND(C)
    USE PLRM, ONLY: MDO
    INTEGER :: MDO_LOCAL
    MDO = MDO_LOCAL
  END SUBROUTINE PLRM_SET_MDO

  ! Getter and setter for INPUT_VECS.
  SUBROUTINE PLRM_GET_INPUT_VECS(INPUT_VECS_ALLOCATED, INPUT_VECS_DIM_1, INPUT_VECS_DIM_2, INPUT_VECS_LOCAL) BIND(C)
    USE ISO_FORTRAN_ENV, ONLY: INT64
    USE PLRM, ONLY: INPUT_VECS
    LOGICAL, INTENT(OUT) :: INPUT_VECS_ALLOCATED
    INTEGER, INTENT(OUT) :: INPUT_VECS_DIM_1
    INTEGER, INTENT(OUT) :: INPUT_VECS_DIM_2
    INTEGER(KIND=INT64) :: INPUT_VECS_LOCAL
    INPUT_VECS_ALLOCATED = ALLOCATED(INPUT_VECS)
    IF (.NOT. INPUT_VECS_ALLOCATED) RETURN
    INPUT_VECS_DIM_1 = SIZE(INPUT_VECS, 1)
    INPUT_VECS_DIM_2 = SIZE(INPUT_VECS, 2)
    INPUT_VECS_LOCAL = LOC(INPUT_VECS(1,1))
  END SUBROUTINE PLRM_GET_INPUT_VECS
  SUBROUTINE PLRM_SET_INPUT_VECS(INPUT_VECS_DIM_1, INPUT_VECS_DIM_2, INPUT_VECS_LOCAL) BIND(C)
    USE PLRM, ONLY: INPUT_VECS
    INTEGER, INTENT(IN) :: INPUT_VECS_DIM_1
    INTEGER, INTENT(IN) :: INPUT_VECS_DIM_2
    REAL(KIND=RT), DIMENSION(INPUT_VECS_DIM_1,INPUT_VECS_DIM_2) :: INPUT_VECS_LOCAL
    IF (ALLOCATED(INPUT_VECS)) THEN
      DEALLOCATE(INPUT_VECS)
    END IF
    ALLOCATE(INPUT_VECS(1:INPUT_VECS_DIM_1,1:INPUT_VECS_DIM_2))
    INPUT_VECS = INPUT_VECS_LOCAL
  END SUBROUTINE PLRM_SET_INPUT_VECS

  ! Getter and setter for INPUT_SHIFT.
  SUBROUTINE PLRM_GET_INPUT_SHIFT(INPUT_SHIFT_ALLOCATED, INPUT_SHIFT_DIM_1, INPUT_SHIFT_LOCAL) BIND(C)
    USE ISO_FORTRAN_ENV, ONLY: INT64
    USE PLRM, ONLY: INPUT_SHIFT
    LOGICAL, INTENT(OUT) :: INPUT_SHIFT_ALLOCATED
    INTEGER, INTENT(OUT) :: INPUT_SHIFT_DIM_1
    INTEGER(KIND=INT64) :: INPUT_SHIFT_LOCAL
    INPUT_SHIFT_ALLOCATED = ALLOCATED(INPUT_SHIFT)
    IF (.NOT. INPUT_SHIFT_ALLOCATED) RETURN
    INPUT_SHIFT_DIM_1 = SIZE(INPUT_SHIFT, 1)
    INPUT_SHIFT_LOCAL = LOC(INPUT_SHIFT(1))
  END SUBROUTINE PLRM_GET_INPUT_SHIFT
  SUBROUTINE PLRM_SET_INPUT_SHIFT(INPUT_SHIFT_DIM_1, INPUT_SHIFT_LOCAL) BIND(C)
    USE PLRM, ONLY: INPUT_SHIFT
    INTEGER, INTENT(IN) :: INPUT_SHIFT_DIM_1
    REAL(KIND=RT), DIMENSION(INPUT_SHIFT_DIM_1) :: INPUT_SHIFT_LOCAL
    IF (ALLOCATED(INPUT_SHIFT)) THEN
      DEALLOCATE(INPUT_SHIFT)
    END IF
    ALLOCATE(INPUT_SHIFT(1:INPUT_SHIFT_DIM_1))
    INPUT_SHIFT = INPUT_SHIFT_LOCAL
  END SUBROUTINE PLRM_SET_INPUT_SHIFT

  ! Getter and setter for INPUT_FLEX.
  SUBROUTINE PLRM_GET_INPUT_FLEX(INPUT_FLEX_ALLOCATED, INPUT_FLEX_DIM_1, INPUT_FLEX_LOCAL) BIND(C)
    USE ISO_FORTRAN_ENV, ONLY: INT64
    USE PLRM, ONLY: INPUT_FLEX
    LOGICAL, INTENT(OUT) :: INPUT_FLEX_ALLOCATED
    INTEGER, INTENT(OUT) :: INPUT_FLEX_DIM_1
    INTEGER(KIND=INT64) :: INPUT_FLEX_LOCAL
    INPUT_FLEX_ALLOCATED = ALLOCATED(INPUT_FLEX)
    IF (.NOT. INPUT_FLEX_ALLOCATED) RETURN
    INPUT_FLEX_DIM_1 = SIZE(INPUT_FLEX, 1)
    INPUT_FLEX_LOCAL = LOC(INPUT_FLEX(1))
  END SUBROUTINE PLRM_GET_INPUT_FLEX
  SUBROUTINE PLRM_SET_INPUT_FLEX(INPUT_FLEX_DIM_1, INPUT_FLEX_LOCAL) BIND(C)
    USE PLRM, ONLY: INPUT_FLEX
    INTEGER, INTENT(IN) :: INPUT_FLEX_DIM_1
    REAL(KIND=RT), DIMENSION(INPUT_FLEX_DIM_1) :: INPUT_FLEX_LOCAL
    IF (ALLOCATED(INPUT_FLEX)) THEN
      DEALLOCATE(INPUT_FLEX)
    END IF
    ALLOCATE(INPUT_FLEX(1:INPUT_FLEX_DIM_1))
    INPUT_FLEX = INPUT_FLEX_LOCAL
  END SUBROUTINE PLRM_SET_INPUT_FLEX

  ! Getter and setter for INTERNAL_VECS.
  SUBROUTINE PLRM_GET_INTERNAL_VECS(INTERNAL_VECS_ALLOCATED, INTERNAL_VECS_DIM_1, INTERNAL_VECS_DIM_2, INTERNAL_VECS_DIM_3, INTERNA&
&L_VECS_LOCAL) BIND(C)
    USE ISO_FORTRAN_ENV, ONLY: INT64
    USE PLRM, ONLY: INTERNAL_VECS
    LOGICAL, INTENT(OUT) :: INTERNAL_VECS_ALLOCATED
    INTEGER, INTENT(OUT) :: INTERNAL_VECS_DIM_1
    INTEGER, INTENT(OUT) :: INTERNAL_VECS_DIM_2
    INTEGER, INTENT(OUT) :: INTERNAL_VECS_DIM_3
    INTEGER(KIND=INT64) :: INTERNAL_VECS_LOCAL
    INTERNAL_VECS_ALLOCATED = ALLOCATED(INTERNAL_VECS)
    IF (.NOT. INTERNAL_VECS_ALLOCATED) RETURN
    INTERNAL_VECS_DIM_1 = SIZE(INTERNAL_VECS, 1)
    INTERNAL_VECS_DIM_2 = SIZE(INTERNAL_VECS, 2)
    INTERNAL_VECS_DIM_3 = SIZE(INTERNAL_VECS, 3)
    INTERNAL_VECS_LOCAL = LOC(INTERNAL_VECS(1,1,1))
  END SUBROUTINE PLRM_GET_INTERNAL_VECS
  SUBROUTINE PLRM_SET_INTERNAL_VECS(INTERNAL_VECS_DIM_1, INTERNAL_VECS_DIM_2, INTERNAL_VECS_DIM_3, INTERNAL_VECS_LOCAL) BIND(C)
    USE PLRM, ONLY: INTERNAL_VECS
    INTEGER, INTENT(IN) :: INTERNAL_VECS_DIM_1
    INTEGER, INTENT(IN) :: INTERNAL_VECS_DIM_2
    INTEGER, INTENT(IN) :: INTERNAL_VECS_DIM_3
    REAL(KIND=RT), DIMENSION(INTERNAL_VECS_DIM_1,INTERNAL_VECS_DIM_2,INTERNAL_VECS_DIM_3) :: INTERNAL_VECS_LOCAL
    IF (ALLOCATED(INTERNAL_VECS)) THEN
      DEALLOCATE(INTERNAL_VECS)
    END IF
    ALLOCATE(INTERNAL_VECS(1:INTERNAL_VECS_DIM_1,1:INTERNAL_VECS_DIM_2,1:INTERNAL_VECS_DIM_3))
    INTERNAL_VECS = INTERNAL_VECS_LOCAL
  END SUBROUTINE PLRM_SET_INTERNAL_VECS

  ! Getter and setter for INTERNAL_SHIFT.
  SUBROUTINE PLRM_GET_INTERNAL_SHIFT(INTERNAL_SHIFT_ALLOCATED, INTERNAL_SHIFT_DIM_1, INTERNAL_SHIFT_DIM_2, INTERNAL_SHIFT_LOCAL) BI&
&ND(C)
    USE ISO_FORTRAN_ENV, ONLY: INT64
    USE PLRM, ONLY: INTERNAL_SHIFT
    LOGICAL, INTENT(OUT) :: INTERNAL_SHIFT_ALLOCATED
    INTEGER, INTENT(OUT) :: INTERNAL_SHIFT_DIM_1
    INTEGER, INTENT(OUT) :: INTERNAL_SHIFT_DIM_2
    INTEGER(KIND=INT64) :: INTERNAL_SHIFT_LOCAL
    INTERNAL_SHIFT_ALLOCATED = ALLOCATED(INTERNAL_SHIFT)
    IF (.NOT. INTERNAL_SHIFT_ALLOCATED) RETURN
    INTERNAL_SHIFT_DIM_1 = SIZE(INTERNAL_SHIFT, 1)
    INTERNAL_SHIFT_DIM_2 = SIZE(INTERNAL_SHIFT, 2)
    INTERNAL_SHIFT_LOCAL = LOC(INTERNAL_SHIFT(1,1))
  END SUBROUTINE PLRM_GET_INTERNAL_SHIFT
  SUBROUTINE PLRM_SET_INTERNAL_SHIFT(INTERNAL_SHIFT_DIM_1, INTERNAL_SHIFT_DIM_2, INTERNAL_SHIFT_LOCAL) BIND(C)
    USE PLRM, ONLY: INTERNAL_SHIFT
    INTEGER, INTENT(IN) :: INTERNAL_SHIFT_DIM_1
    INTEGER, INTENT(IN) :: INTERNAL_SHIFT_DIM_2
    REAL(KIND=RT), DIMENSION(INTERNAL_SHIFT_DIM_1,INTERNAL_SHIFT_DIM_2) :: INTERNAL_SHIFT_LOCAL
    IF (ALLOCATED(INTERNAL_SHIFT)) THEN
      DEALLOCATE(INTERNAL_SHIFT)
    END IF
    ALLOCATE(INTERNAL_SHIFT(1:INTERNAL_SHIFT_DIM_1,1:INTERNAL_SHIFT_DIM_2))
    INTERNAL_SHIFT = INTERNAL_SHIFT_LOCAL
  END SUBROUTINE PLRM_SET_INTERNAL_SHIFT

  ! Getter and setter for INTERNAL_FLEX.
  SUBROUTINE PLRM_GET_INTERNAL_FLEX(INTERNAL_FLEX_ALLOCATED, INTERNAL_FLEX_DIM_1, INTERNAL_FLEX_DIM_2, INTERNAL_FLEX_LOCAL) BIND(C)
    USE ISO_FORTRAN_ENV, ONLY: INT64
    USE PLRM, ONLY: INTERNAL_FLEX
    LOGICAL, INTENT(OUT) :: INTERNAL_FLEX_ALLOCATED
    INTEGER, INTENT(OUT) :: INTERNAL_FLEX_DIM_1
    INTEGER, INTENT(OUT) :: INTERNAL_FLEX_DIM_2
    INTEGER(KIND=INT64) :: INTERNAL_FLEX_LOCAL
    INTERNAL_FLEX_ALLOCATED = ALLOCATED(INTERNAL_FLEX)
    IF (.NOT. INTERNAL_FLEX_ALLOCATED) RETURN
    INTERNAL_FLEX_DIM_1 = SIZE(INTERNAL_FLEX, 1)
    INTERNAL_FLEX_DIM_2 = SIZE(INTERNAL_FLEX, 2)
    INTERNAL_FLEX_LOCAL = LOC(INTERNAL_FLEX(1,1))
  END SUBROUTINE PLRM_GET_INTERNAL_FLEX
  SUBROUTINE PLRM_SET_INTERNAL_FLEX(INTERNAL_FLEX_DIM_1, INTERNAL_FLEX_DIM_2, INTERNAL_FLEX_LOCAL) BIND(C)
    USE PLRM, ONLY: INTERNAL_FLEX
    INTEGER, INTENT(IN) :: INTERNAL_FLEX_DIM_1
    INTEGER, INTENT(IN) :: INTERNAL_FLEX_DIM_2
    REAL(KIND=RT), DIMENSION(INTERNAL_FLEX_DIM_1,INTERNAL_FLEX_DIM_2) :: INTERNAL_FLEX_LOCAL
    IF (ALLOCATED(INTERNAL_FLEX)) THEN
      DEALLOCATE(INTERNAL_FLEX)
    END IF
    ALLOCATE(INTERNAL_FLEX(1:INTERNAL_FLEX_DIM_1,1:INTERNAL_FLEX_DIM_2))
    INTERNAL_FLEX = INTERNAL_FLEX_LOCAL
  END SUBROUTINE PLRM_SET_INTERNAL_FLEX

  ! Getter and setter for OUTPUT_VECS.
  SUBROUTINE PLRM_GET_OUTPUT_VECS(OUTPUT_VECS_ALLOCATED, OUTPUT_VECS_DIM_1, OUTPUT_VECS_DIM_2, OUTPUT_VECS_LOCAL) BIND(C)
    USE ISO_FORTRAN_ENV, ONLY: INT64
    USE PLRM, ONLY: OUTPUT_VECS
    LOGICAL, INTENT(OUT) :: OUTPUT_VECS_ALLOCATED
    INTEGER, INTENT(OUT) :: OUTPUT_VECS_DIM_1
    INTEGER, INTENT(OUT) :: OUTPUT_VECS_DIM_2
    INTEGER(KIND=INT64) :: OUTPUT_VECS_LOCAL
    OUTPUT_VECS_ALLOCATED = ALLOCATED(OUTPUT_VECS)
    IF (.NOT. OUTPUT_VECS_ALLOCATED) RETURN
    OUTPUT_VECS_DIM_1 = SIZE(OUTPUT_VECS, 1)
    OUTPUT_VECS_DIM_2 = SIZE(OUTPUT_VECS, 2)
    OUTPUT_VECS_LOCAL = LOC(OUTPUT_VECS(1,1))
  END SUBROUTINE PLRM_GET_OUTPUT_VECS
  SUBROUTINE PLRM_SET_OUTPUT_VECS(OUTPUT_VECS_DIM_1, OUTPUT_VECS_DIM_2, OUTPUT_VECS_LOCAL) BIND(C)
    USE PLRM, ONLY: OUTPUT_VECS
    INTEGER, INTENT(IN) :: OUTPUT_VECS_DIM_1
    INTEGER, INTENT(IN) :: OUTPUT_VECS_DIM_2
    REAL(KIND=RT), DIMENSION(OUTPUT_VECS_DIM_1,OUTPUT_VECS_DIM_2) :: OUTPUT_VECS_LOCAL
    IF (ALLOCATED(OUTPUT_VECS)) THEN
      DEALLOCATE(OUTPUT_VECS)
    END IF
    ALLOCATE(OUTPUT_VECS(1:OUTPUT_VECS_DIM_1,1:OUTPUT_VECS_DIM_2))
    OUTPUT_VECS = OUTPUT_VECS_LOCAL
  END SUBROUTINE PLRM_SET_OUTPUT_VECS

  ! Getter and setter for OUTPUT_SHIFT.
  SUBROUTINE PLRM_GET_OUTPUT_SHIFT(OUTPUT_SHIFT_ALLOCATED, OUTPUT_SHIFT_DIM_1, OUTPUT_SHIFT_LOCAL) BIND(C)
    USE ISO_FORTRAN_ENV, ONLY: INT64
    USE PLRM, ONLY: OUTPUT_SHIFT
    LOGICAL, INTENT(OUT) :: OUTPUT_SHIFT_ALLOCATED
    INTEGER, INTENT(OUT) :: OUTPUT_SHIFT_DIM_1
    INTEGER(KIND=INT64) :: OUTPUT_SHIFT_LOCAL
    OUTPUT_SHIFT_ALLOCATED = ALLOCATED(OUTPUT_SHIFT)
    IF (.NOT. OUTPUT_SHIFT_ALLOCATED) RETURN
    OUTPUT_SHIFT_DIM_1 = SIZE(OUTPUT_SHIFT, 1)
    OUTPUT_SHIFT_LOCAL = LOC(OUTPUT_SHIFT(1))
  END SUBROUTINE PLRM_GET_OUTPUT_SHIFT
  SUBROUTINE PLRM_SET_OUTPUT_SHIFT(OUTPUT_SHIFT_DIM_1, OUTPUT_SHIFT_LOCAL) BIND(C)
    USE PLRM, ONLY: OUTPUT_SHIFT
    INTEGER, INTENT(IN) :: OUTPUT_SHIFT_DIM_1
    REAL(KIND=RT), DIMENSION(OUTPUT_SHIFT_DIM_1) :: OUTPUT_SHIFT_LOCAL
    IF (ALLOCATED(OUTPUT_SHIFT)) THEN
      DEALLOCATE(OUTPUT_SHIFT)
    END IF
    ALLOCATE(OUTPUT_SHIFT(1:OUTPUT_SHIFT_DIM_1))
    OUTPUT_SHIFT = OUTPUT_SHIFT_LOCAL
  END SUBROUTINE PLRM_SET_OUTPUT_SHIFT

  
  SUBROUTINE C_NEW_MODEL(DI, DS, NS, DO) BIND(C)
    USE PLRM, ONLY: NEW_MODEL
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: DI
    INTEGER, INTENT(IN) :: DS
    INTEGER, INTENT(IN) :: NS
    INTEGER, INTENT(IN) :: DO
  
    CALL NEW_MODEL(DI, DS, NS, DO)
  END SUBROUTINE C_NEW_MODEL
  

  
  SUBROUTINE C_INIT_MODEL(SEED_PRESENT, SEED) BIND(C)
    USE PLRM, ONLY: INIT_MODEL
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: SEED_PRESENT
    INTEGER, INTENT(IN) :: SEED
  
    IF (SEED_PRESENT) THEN
      CALL INIT_MODEL(SEED=SEED)
    ELSE
      CALL INIT_MODEL()
    END IF
  END SUBROUTINE C_INIT_MODEL
  

  
  SUBROUTINE C_RANDOM_UNIT_VECTORS(COLUMN_VECTORS_DIM_1, COLUMN_VECTORS_DIM_2, COLUMN_VECTORS) BIND(C)
    USE PLRM, ONLY: RANDOM_UNIT_VECTORS
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: COLUMN_VECTORS_DIM_1
    INTEGER, INTENT(IN) :: COLUMN_VECTORS_DIM_2
    REAL(KIND=RT), INTENT(OUT), DIMENSION(COLUMN_VECTORS_DIM_1,COLUMN_VECTORS_DIM_2) :: COLUMN_VECTORS
  
    CALL RANDOM_UNIT_VECTORS(COLUMN_VECTORS)
  END SUBROUTINE C_RANDOM_UNIT_VECTORS
  

  
  SUBROUTINE C_EVALUATE(INPUTS_DIM_1, INPUTS_DIM_2, INPUTS, OUTPUTS_DIM_1, OUTPUTS_DIM_2, OUTPUTS) BIND(C)
    USE PLRM, ONLY: EVALUATE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: INPUTS_DIM_1
    INTEGER, INTENT(IN) :: INPUTS_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(INPUTS_DIM_1,INPUTS_DIM_2) :: INPUTS
    INTEGER, INTENT(IN) :: OUTPUTS_DIM_1
    INTEGER, INTENT(IN) :: OUTPUTS_DIM_2
    REAL(KIND=RT), INTENT(OUT), DIMENSION(OUTPUTS_DIM_1,OUTPUTS_DIM_2) :: OUTPUTS
  
    CALL EVALUATE(INPUTS, OUTPUTS)
  END SUBROUTINE C_EVALUATE
  

  
  SUBROUTINE C_MINIMIZE_MSE(X_DIM_1, X_DIM_2, X, Y_DIM_1, Y_DIM_2, Y, STEPS, BATCH_SIZE_PRESENT, BATCH_SIZE, NUM_THREADS_PRESENT, N&
&UM_THREADS, STEP_SIZE_PRESENT, STEP_SIZE, KEEP_BEST_PRESENT, KEEP_BEST, MEAN_SQUARED_ERROR, RECORD_PRESENT, RECORD_DIM_1, RECORD, &
&LOGS_PRESENT, LOGS_DIM_1, LOGS_DIM_2, LOGS_DIM_3, LOGS_DIM_4, LOGS) BIND(C)
    USE VIS_PLRM , ONLY : LOG_STATUS
    USE PLRM, ONLY: MINIMIZE_MSE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: X_DIM_1
    INTEGER, INTENT(IN) :: X_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(X_DIM_1,X_DIM_2) :: X
    INTEGER, INTENT(IN) :: Y_DIM_1
    INTEGER, INTENT(IN) :: Y_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(Y_DIM_1,Y_DIM_2) :: Y
    INTEGER, INTENT(IN) :: STEPS
    LOGICAL, INTENT(IN) :: BATCH_SIZE_PRESENT
    INTEGER, INTENT(IN) :: BATCH_SIZE
    LOGICAL, INTENT(IN) :: NUM_THREADS_PRESENT
    INTEGER, INTENT(IN) :: NUM_THREADS
    LOGICAL, INTENT(IN) :: STEP_SIZE_PRESENT
    REAL(KIND=RT), INTENT(IN) :: STEP_SIZE
    LOGICAL, INTENT(IN) :: KEEP_BEST_PRESENT
    LOGICAL, INTENT(IN) :: KEEP_BEST
    REAL(KIND=RT), INTENT(OUT) :: MEAN_SQUARED_ERROR
    LOGICAL, INTENT(IN) :: RECORD_PRESENT
    INTEGER, INTENT(IN) :: RECORD_DIM_1
    REAL(KIND=RT), INTENT(OUT), DIMENSION(RECORD_DIM_1) :: RECORD
    LOGICAL, INTENT(IN) :: LOGS_PRESENT
    INTEGER, INTENT(IN) :: LOGS_DIM_1
    INTEGER, INTENT(IN) :: LOGS_DIM_2
    INTEGER, INTENT(IN) :: LOGS_DIM_3
    INTEGER, INTENT(IN) :: LOGS_DIM_4
    REAL(KIND=RT), INTENT(OUT), DIMENSION(LOGS_DIM_1,LOGS_DIM_2,LOGS_DIM_3,LOGS_DIM_4) :: LOGS
  
    IF (BATCH_SIZE_PRESENT) THEN
      IF (NUM_THREADS_PRESENT) THEN
        IF (STEP_SIZE_PRESENT) THEN
          IF (KEEP_BEST_PRESENT) THEN
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, STEP_SIZE=STEP_SIZE, KEEP_BEST=KEEP_BEST, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, STEP_SIZE=STEP_SIZE, KEEP_BEST=KEEP_BEST, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, STEP_SIZE=STEP_SIZE, KEEP_BEST=KEEP_BEST, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, STEP_SIZE=STEP_SIZE, KEEP_BEST=KEEP_BEST)
              END IF
            END IF
          ELSE
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, STEP_SIZE=STEP_SIZE, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, STEP_SIZE=STEP_SIZE, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, STEP_SIZE=STEP_SIZE, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, STEP_SIZE=STEP_SIZE)
              END IF
            END IF
          END IF
        ELSE
          IF (KEEP_BEST_PRESENT) THEN
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, KEEP_BEST=KEEP_BEST, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, KEEP_BEST=KEEP_BEST, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, KEEP_BEST=KEEP_BEST, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, KEEP_BEST=KEEP_BEST)
              END IF
            END IF
          ELSE
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, NUM_THREADS=&
&NUM_THREADS)
              END IF
            END IF
          END IF
        END IF
      ELSE
        IF (STEP_SIZE_PRESENT) THEN
          IF (KEEP_BEST_PRESENT) THEN
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, STEP_SIZE=ST&
&EP_SIZE, KEEP_BEST=KEEP_BEST, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, STEP_SIZE=ST&
&EP_SIZE, KEEP_BEST=KEEP_BEST, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, STEP_SIZE=ST&
&EP_SIZE, KEEP_BEST=KEEP_BEST, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, STEP_SIZE=ST&
&EP_SIZE, KEEP_BEST=KEEP_BEST)
              END IF
            END IF
          ELSE
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, STEP_SIZE=ST&
&EP_SIZE, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, STEP_SIZE=ST&
&EP_SIZE, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, STEP_SIZE=ST&
&EP_SIZE, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, STEP_SIZE=ST&
&EP_SIZE)
              END IF
            END IF
          END IF
        ELSE
          IF (KEEP_BEST_PRESENT) THEN
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, KEEP_BEST=KE&
&EP_BEST, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, KEEP_BEST=KE&
&EP_BEST, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, KEEP_BEST=KE&
&EP_BEST, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, KEEP_BEST=KE&
&EP_BEST)
              END IF
            END IF
          ELSE
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, RECORD=RECOR&
&D, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, RECORD=RECOR&
&D)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, BATCH_SIZE=BATCH_SIZE)
              END IF
            END IF
          END IF
        END IF
      END IF
    ELSE
      IF (NUM_THREADS_PRESENT) THEN
        IF (STEP_SIZE_PRESENT) THEN
          IF (KEEP_BEST_PRESENT) THEN
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, STEP_SIZE=&
&STEP_SIZE, KEEP_BEST=KEEP_BEST, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, STEP_SIZE=&
&STEP_SIZE, KEEP_BEST=KEEP_BEST, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, STEP_SIZE=&
&STEP_SIZE, KEEP_BEST=KEEP_BEST, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, STEP_SIZE=&
&STEP_SIZE, KEEP_BEST=KEEP_BEST)
              END IF
            END IF
          ELSE
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, STEP_SIZE=&
&STEP_SIZE, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, STEP_SIZE=&
&STEP_SIZE, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, STEP_SIZE=&
&STEP_SIZE, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, STEP_SIZE=&
&STEP_SIZE)
              END IF
            END IF
          END IF
        ELSE
          IF (KEEP_BEST_PRESENT) THEN
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, KEEP_BEST=&
&KEEP_BEST, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, KEEP_BEST=&
&KEEP_BEST, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, KEEP_BEST=&
&KEEP_BEST, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, KEEP_BEST=&
&KEEP_BEST)
              END IF
            END IF
          ELSE
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, RECORD=REC&
&ORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, RECORD=REC&
&ORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, NUM_THREADS=NUM_THREADS)
              END IF
            END IF
          END IF
        END IF
      ELSE
        IF (STEP_SIZE_PRESENT) THEN
          IF (KEEP_BEST_PRESENT) THEN
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, STEP_SIZE=STEP_SIZE, KEEP_BEST=KEEP&
&_BEST, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, STEP_SIZE=STEP_SIZE, KEEP_BEST=KEEP&
&_BEST, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, STEP_SIZE=STEP_SIZE, KEEP_BEST=KEEP&
&_BEST, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, STEP_SIZE=STEP_SIZE, KEEP_BEST=KEEP&
&_BEST)
              END IF
            END IF
          ELSE
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, STEP_SIZE=STEP_SIZE, RECORD=RECORD,&
& LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, STEP_SIZE=STEP_SIZE, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, STEP_SIZE=STEP_SIZE, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, STEP_SIZE=STEP_SIZE)
              END IF
            END IF
          END IF
        ELSE
          IF (KEEP_BEST_PRESENT) THEN
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, KEEP_BEST=KEEP_BEST, RECORD=RECORD,&
& LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, KEEP_BEST=KEEP_BEST, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, KEEP_BEST=KEEP_BEST, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, KEEP_BEST=KEEP_BEST)
              END IF
            END IF
          ELSE
            IF (RECORD_PRESENT) THEN
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, RECORD=RECORD, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, RECORD=RECORD)
              END IF
            ELSE
              IF (LOGS_PRESENT) THEN
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR, LOGS=LOGS)
              ELSE
                CALL MINIMIZE_MSE(X=X, Y=Y, STEPS=STEPS, MEAN_SQUARED_ERROR=MEAN_SQUARED_ERROR)
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF
  END SUBROUTINE C_MINIMIZE_MSE
  
END MODULE C_PLRM


MODULE C_VIS_PLRM
USE ISO_FORTRAN_ENV , ONLY : RT => REAL32
USE MATRIX_MULTIPLICATION , ONLY : GEMM
USE PLRM , ONLY : MDI , MDS , MDO , MNS , DISCONTINUITY , INPUT_VECS , INPUT_SHIFT , INPUT_FLEX , INTERNAL_VECS , INTERNAL_SHIFT , &
&INTERNAL_FLEX , OUTPUT_VECS , OUTPUT_SHIFT , EVALUATE
  IMPLICIT NONE


CONTAINS


  
  SUBROUTINE C_LOG_STATUS(X_DIM_1, X_DIM_2, X, Y_DIM_1, Y_DIM_2, Y, LOG_VALS_DIM_1, LOG_VALS_DIM_2, LOG_VALS_DIM_3, LOG_VALS) BIND(&
&C)
    USE VIS_PLRM, ONLY: LOG_STATUS
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: X_DIM_1
    INTEGER, INTENT(IN) :: X_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(X_DIM_1,X_DIM_2) :: X
    INTEGER, INTENT(IN) :: Y_DIM_1
    INTEGER, INTENT(IN) :: Y_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(Y_DIM_1,Y_DIM_2) :: Y
    INTEGER, INTENT(IN) :: LOG_VALS_DIM_1
    INTEGER, INTENT(IN) :: LOG_VALS_DIM_2
    INTEGER, INTENT(IN) :: LOG_VALS_DIM_3
    REAL(KIND=RT), INTENT(OUT), DIMENSION(LOG_VALS_DIM_1,LOG_VALS_DIM_2,LOG_VALS_DIM_3) :: LOG_VALS
  
    CALL LOG_STATUS(X, Y, LOG_VALS)
  END SUBROUTINE C_LOG_STATUS
  

  
  SUBROUTINE C_DISABLE_AND_COMPUTE_MSE(LAYER, POSITION, INPUTS_DIM_1, INPUTS_DIM_2, INPUTS, OUTPUTS_DIM_1, OUTPUTS_DIM_2, OUTPUTS, &
&MSE) BIND(C)
    USE VIS_PLRM, ONLY: DISABLE_AND_COMPUTE_MSE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LAYER
    INTEGER, INTENT(IN) :: POSITION
    INTEGER, INTENT(IN) :: INPUTS_DIM_1
    INTEGER, INTENT(IN) :: INPUTS_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(INPUTS_DIM_1,INPUTS_DIM_2) :: INPUTS
    INTEGER, INTENT(IN) :: OUTPUTS_DIM_1
    INTEGER, INTENT(IN) :: OUTPUTS_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(OUTPUTS_DIM_1,OUTPUTS_DIM_2) :: OUTPUTS
    REAL(KIND=RT), INTENT(OUT) :: MSE
  
    CALL DISABLE_AND_COMPUTE_MSE(LAYER, POSITION, INPUTS, OUTPUTS, MSE)
  END SUBROUTINE C_DISABLE_AND_COMPUTE_MSE
  

  
  SUBROUTINE C_COMPUTE_MSE(INPUTS_DIM_1, INPUTS_DIM_2, INPUTS, OUTPUTS_DIM_1, OUTPUTS_DIM_2, OUTPUTS, MSE) BIND(C)
    USE VIS_PLRM, ONLY: COMPUTE_MSE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: INPUTS_DIM_1
    INTEGER, INTENT(IN) :: INPUTS_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(INPUTS_DIM_1,INPUTS_DIM_2) :: INPUTS
    INTEGER, INTENT(IN) :: OUTPUTS_DIM_1
    INTEGER, INTENT(IN) :: OUTPUTS_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(OUTPUTS_DIM_1,OUTPUTS_DIM_2) :: OUTPUTS
    REAL(KIND=RT), INTENT(OUT) :: MSE
  
    CALL COMPUTE_MSE(INPUTS, OUTPUTS, MSE)
  END SUBROUTINE C_COMPUTE_MSE
  

  
  SUBROUTINE C_COMPUTE_VALUES(LAYER, POSITION, INPUTS_DIM_1, INPUTS_DIM_2, INPUTS, OUTPUTS_DIM_1, OUTPUTS_DIM_2, OUTPUTS, VALS_DIM_&
&1, VALS, GRADS_DIM_1, GRADS) BIND(C)
    USE VIS_PLRM, ONLY: COMPUTE_VALUES
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LAYER
    INTEGER, INTENT(IN) :: POSITION
    INTEGER, INTENT(IN) :: INPUTS_DIM_1
    INTEGER, INTENT(IN) :: INPUTS_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(INPUTS_DIM_1,INPUTS_DIM_2) :: INPUTS
    INTEGER, INTENT(IN) :: OUTPUTS_DIM_1
    INTEGER, INTENT(IN) :: OUTPUTS_DIM_2
    REAL(KIND=RT), INTENT(IN), DIMENSION(OUTPUTS_DIM_1,OUTPUTS_DIM_2) :: OUTPUTS
    INTEGER, INTENT(IN) :: VALS_DIM_1
    REAL(KIND=RT), INTENT(OUT), DIMENSION(VALS_DIM_1) :: VALS
    INTEGER, INTENT(IN) :: GRADS_DIM_1
    REAL(KIND=RT), INTENT(OUT), DIMENSION(GRADS_DIM_1) :: GRADS
  
    CALL COMPUTE_VALUES(LAYER, POSITION, INPUTS, OUTPUTS, VALS, GRADS)
  END SUBROUTINE C_COMPUTE_VALUES
  
END MODULE C_VIS_PLRM
