
! Module for matrix multiplication (absolutely crucial for PLRM speed).
MODULE MATRIX_MULTIPLICATION
  USE ISO_FORTRAN_ENV, ONLY: RT => REAL32
  IMPLICIT NONE

CONTAINS

  ! Convenience wrapper routine for calling matrix multiply.
  SUBROUTINE GEMM(OP_A, OP_B, OUT_ROWS, OUT_COLS, INNER_DIM, &
       AB_MULT, A, A_ROWS, B, B_ROWS, C_MULT, C, C_ROWS)
    CHARACTER, INTENT(IN) :: OP_A, OP_B
    INTEGER, INTENT(IN) :: OUT_ROWS, OUT_COLS, INNER_DIM, A_ROWS, B_ROWS, C_ROWS
    REAL(KIND=RT), INTENT(IN) :: AB_MULT, C_MULT
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: A, B
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: C
    ! Call external single-precision matrix-matrix multiplication
    !  (should be provided by hardware manufacturer, if not write custom).
    EXTERNAL :: SGEMM 
    CALL SGEMM(OP_A, OP_B, OUT_ROWS, OUT_COLS, INNER_DIM, &
       AB_MULT, A, A_ROWS, B, B_ROWS, C_MULT, C, C_ROWS)
  END SUBROUTINE GEMM

END MODULE MATRIX_MULTIPLICATION


! A piecewise linear regression model.
MODULE PLRM
  USE ISO_FORTRAN_ENV, ONLY: RT => REAL32
  USE MATRIX_MULTIPLICATION, ONLY: GEMM
  IMPLICIT NONE
  ! Parameters that define model properties.
  REAL(KIND=RT), PARAMETER :: DISCONTINUITY = 0.0_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_SHIFT_RANGE = 1.0_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_OUTPUT_SCALE = 0.001_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_STEP = 0.001_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_STEP_MEAN_CHANGE = 0.1_RT  
  REAL(KIND=RT), PARAMETER :: INITIAL_STEP_CURV_CHANGE = 0.01_RT  
  REAL(KIND=RT), PARAMETER :: STEP_GROWTH_RATE = 1.001_RT
  REAL(KIND=RT), PARAMETER :: STEP_SHRINK_RATE = 0.999_RT
  INTEGER,       PARAMETER :: MIN_STEPS_TO_STABILITY = 1

  ! MDI = Model dimension -- input
  ! MDS = Model dimension -- states (internal)
  ! MNS = Model number    -- states (internal)
  ! MDO = Model dimension -- output
  INTEGER :: MDI, MDS, MNS, MDO
  ! Model variables, for evaluation.
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INPUT_VECS     ! MDI, MDS
  REAL(KIND=RT), DIMENSION(:),     ALLOCATABLE :: INPUT_SHIFT    ! MDS
  REAL(KIND=RT), DIMENSION(:,:,:), ALLOCATABLE :: INTERNAL_VECS  ! MDS, MDS, MNS-1
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_SHIFT ! MDS, MNS-1
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: OUTPUT_VECS    ! MDS, MDO
  REAL(KIND=RT), DIMENSION(:),     ALLOCATABLE :: OUTPUT_SHIFT   ! MDO

  ! Total allocated memory is equal to:
  !   MDI MDS + 2 MDS + MDS (MDS(MNS-1) + 2(MNS-1) + MDO) + MDO
  !   MDS (MDI + 2) + MDS ((MNS-1)(MDS + 2) + MDO) + MDO

  PUBLIC :: NEW_MODEL, INIT_MODEL, RANDOM_UNIT_VECTORS, MINIMIZE_MSE
  PRIVATE :: SSE_GRADIENT

CONTAINS

  ! Allocate storage for a new model given paramters:
  !   DI -- Dimension of input vector.
  !   DS -- Dimension of internal vectors spaces.
  !   NS -- Number of sequential internal vector spaces.
  !   DO -- Dimension of output vector.
  SUBROUTINE NEW_MODEL(DI, DS, NS, DO)
    INTEGER, INTENT(IN) :: DI, DS, NS, DO
    ! Store the integer sizes related to this model.
    MDI = DI
    MDS = DS
    MNS = NS
    MDO = DO
    ! Deallocate any of the internal arrays if they have already been allocated.
    IF (ALLOCATED(INPUT_VECS))     DEALLOCATE(INPUT_VECS)
    IF (ALLOCATED(INPUT_SHIFT))    DEALLOCATE(INPUT_SHIFT)
    IF (ALLOCATED(INTERNAL_VECS))  DEALLOCATE(INTERNAL_VECS)
    IF (ALLOCATED(INTERNAL_SHIFT)) DEALLOCATE(INTERNAL_SHIFT)
    IF (ALLOCATED(OUTPUT_VECS))    DEALLOCATE(OUTPUT_VECS)
    IF (ALLOCATED(OUTPUT_SHIFT))   DEALLOCATE(OUTPUT_SHIFT)
    ! Allocate space for all model variables.
    ALLOCATE( &
         INPUT_VECS(DI,DS), &
         INPUT_SHIFT(DS), &
         INTERNAL_VECS(DS, DS, NS-1), &
         INTERNAL_SHIFT(DS, NS-1), &
         OUTPUT_VECS(DS, DO), &
         OUTPUT_SHIFT(DO) &
         )
  END SUBROUTINE NEW_MODEL


  ! Initialize the weights for a model, optionally provide a random seed.
  SUBROUTINE INIT_MODEL(SEED)
    INTEGER, INTENT(IN), OPTIONAL :: SEED
    !  Storage for seeding the random number generator (for repeatability).
    INTEGER, DIMENSION(:), ALLOCATABLE :: SEED_ARRAY
    ! Local iterator.
    INTEGER :: I
    ! Set a random seed, if one was provided (otherwise leave default).
    IF (PRESENT(SEED)) THEN
       CALL RANDOM_SEED(SIZE=I)
       ALLOCATE(SEED_ARRAY(I))
       SEED_ARRAY(:) = SEED
       CALL RANDOM_SEED(PUT=SEED_ARRAY(:))
    END IF
    ! Generate well spaced random unit-length vectors (no scaling biases)
    !  for all initial variables in the input, internal, and output.
    CALL RANDOM_UNIT_VECTORS(INPUT_VECS(:,:))
    DO I = 1, MNS-1
       CALL RANDOM_UNIT_VECTORS(INTERNAL_VECS(:,:,I))
    END DO
    CALL RANDOM_UNIT_VECTORS(OUTPUT_VECS(:,:))
    ! Make the output vectors have very small magnitude initially.
    OUTPUT_VECS(:,:) = OUTPUT_VECS(:,:) * INITIAL_OUTPUT_SCALE
    ! Generate random shifts for inputs and internal layers, zero
    !  shift for the output layer (first two will be rescaled).
    CALL RANDOM_NUMBER(INPUT_SHIFT(:))
    CALL RANDOM_NUMBER(INTERNAL_SHIFT(:,:))
    OUTPUT_SHIFT(:) = 0.0_RT
    ! Assuming unit standard deviaiton, a shift range will capture N
    !  standard deviations of data in the first layer, unclear in later
    !  later layers. TODO: research the distribution of values observed
    !  at any given component given this initialization scheme.
    INPUT_SHIFT(:) = INITIAL_SHIFT_RANGE * (INPUT_SHIFT(:) * 2.0_RT - 1.0_RT)
    INTERNAL_SHIFT(:,:) = INITIAL_SHIFT_RANGE * (INTERNAL_SHIFT(:,:) * 2.0_RT - 1.0_RT)
  END SUBROUTINE INIT_MODEL

  ! Generate randomly distributed vectors on the N-sphere.
  SUBROUTINE RANDOM_UNIT_VECTORS(COLUMN_VECTORS)
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: COLUMN_VECTORS
    REAL(KIND=RT), DIMENSION(SIZE(COLUMN_VECTORS,1), SIZE(COLUMN_VECTORS,2)) :: TEMP_VECS
    REAL(KIND=RT), PARAMETER :: PI = 3.141592653589793
    INTEGER :: I 
    ! Generate random numbers in the range [0,1].
    CALL RANDOM_NUMBER(COLUMN_VECTORS(:,:))
    CALL RANDOM_NUMBER(TEMP_VECS(:,:))
    ! Map the random uniform numbers to a normal distribution.
    COLUMN_VECTORS(:,:) = SQRT(-LOG(COLUMN_VECTORS(:,:))) * COS(PI * TEMP_VECS(:,:))
    ! Make the vectors uniformly distributed on the unit ball (for dimension > 1).
    IF (SIZE(COLUMN_VECTORS,1) .GT. 1) THEN
       ! Normalize all vectors to have unit length.
       DO I = 1, SIZE(COLUMN_VECTORS,2)
          COLUMN_VECTORS(:,I) = COLUMN_VECTORS(:,I) / NORM2(COLUMN_VECTORS(:,I))
       END DO
    END IF
    ! Orthogonalize the first components of the column
    !  vectors to ensure those are well spaced.
    I = MIN(SIZE(COLUMN_VECTORS,1), SIZE(COLUMN_VECTORS,2))
    IF (I .GT. 1) CALL ORTHOGONALIZE(COLUMN_VECTORS(:,1:I))

  CONTAINS

    ! Orthogonalize A with Householder reflections.
    SUBROUTINE ORTHOGONALIZE(A)
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:) :: A
      ! Local variables (perform orthogonalization by computing
      !  Householder reflectors that make A upper diagonal).
      REAL(KIND=RT), DIMENSION(SIZE(A,1), SIZE(A,2)) :: Q
      REAL(KIND=RT), DIMENSION(SIZE(A,1), SIZE(A,1)) :: QI
      INTEGER :: I, J, M, N, L
      REAL(KIND=RT) :: R
      M = SIZE(A,1)
      N = SIZE(A,2)
      L = MIN(M, N)
      ! Set Q to have ones on the diagonal.
      Q(:,:) = 0.0_RT
      FORALL (I=1:L) Q(I,I) = 1.0_RT
      ! Set Qi to the identity.
      QI(:,:) = 0.0_RT
      FORALL (I=1:M) QI(I,I) = 1.0_RT
      ! Compute Householder reflectors.
      DO I = 1, L-1
         ! Copy the column vector from A and compute its length.
         QI(I:M,I) = A(I:M,I)
         QI(I,I) = NORM2(QI(I:M,I))
         ! Make the norm have the opposite sign of the first component.
         IF (A(I,I) .LT. 0.0_RT) QI(I,I) = -QI(I,I)
         ! Subtract this computed value from the first component and
         !  normalize the vector to have unit length.
         QI(I,I) = A(I,I) - QI(I,I)
         QI(I:M,I) = QI(I:M,I) / NORM2(QI(I:M,I))
         ! Perform the outer product between the new vector and itself.
         FORALL (J=I+1:M)
            QI(I:M,J) = QI(J,I) * QI(I:M,I)
         END FORALL
         QI(I:M,I) = QI(I,I) * QI(I:M,I)
         ! Subtract twice this value all from the identity matrix.
         QI(I:M,I:M) = -2.0_RT * QI(I:M,I:M)
         FORALL (J=I:M) QI(J,J) = QI(J,J) + 1.0_RT
         ! Update the matrix as well as Q.
         !  TODO: Use SGEMM instead of builtin MATMUL.
         A(:,:) = MATMUL(QI(:,:), A(:,:))
         Q(:,:) = MATMUL(QI(:,:), Q(:,:))
         ! Reset the current row and column of QI to be from the identity.
         QI(I,I) = 1.0_RT
         QI(I+1:M,I) = 0.0_RT
         QI(I,I+1:M) = 0.0_RT
      END DO
      ! Overwrite A with the orthognal Q.
      A(1:M,1:L) = Q(1:M,1:L)
      ! Randomly flip the sign of some of the vectors in A.
      DO I = 1, L
         CALL RANDOM_NUMBER(R)
         IF (R .LT. 0.5_RT) THEN
            A(:,I) = -A(:,I)
         END IF
      END DO
    END SUBROUTINE ORTHOGONALIZE

  END SUBROUTINE RANDOM_UNIT_VECTORS



  ! Evaluate the piecewise linear regression model.
  SUBROUTINE EVALUATE(INPUTS, OUTPUTS)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: OUTPUTS
    ! Internal values.
    REAL(KIND=RT), DIMENSION(:,:), ALLOCATABLE :: TEMP_VALUES, INTERNAL_VALUES
    INTEGER :: I, D, L, N
    ! Make sure that this model has been initialized.
    IF (.NOT. ALLOCATED(INPUT_VECS)) THEN
       OUTPUTS(:,:) = 0.0
       RETURN
    END IF
    ! Get the number of points.
    N = SIZE(INPUTS,2)
    ! Allocate storage for the internal values.
    ALLOCATE(INTERNAL_VALUES(1:N,1:MDS), TEMP_VALUES(1:N,1:MDS))
    ! Compute the input transformation.
    CALL GEMM('T', 'N', N, MDS, MDI, 1.0_RT, &
         INPUTS(:,:), SIZE(INPUTS,1), &
         INPUT_VECS(:,:), SIZE(INPUT_VECS,1), &
         0.0_RT, INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1))
    DO D = 1, MDS
       INTERNAL_VALUES(:,D) = MAX( &
            INTERNAL_VALUES(:,D) + INPUT_SHIFT(D), &
            DISCONTINUITY )
    END DO
    ! Compute the next set of internal values with a rectified activation.
    DO L = 1, MNS-1
       ! Compute all vectors.
       CALL GEMM('N', 'N', N, MDS, MDS, 1.0_RT, &
            INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
            INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
            0.0_RT, TEMP_VALUES(:,:), SIZE(TEMP_VALUES,1))
       INTERNAL_VALUES(:,:) = TEMP_VALUES(:,:)
       DO D = 1, MDS
          INTERNAL_VALUES(:,D) = MAX( &
               INTERNAL_VALUES(:,D) + INTERNAL_SHIFT(D,L), &
               DISCONTINUITY )
       END DO
    END DO
    ! Return the final embedded layer if the outputs have the size of the embedding.
    IF ((MDS .NE. MDO) .AND. (SIZE(OUTPUTS,1) .EQ. MDS)) THEN
       ! Do the necessary transpose operation one dimension at a time, produce embedding.
       DO D = 1, MDS
          OUTPUTS(D,1:N) = INTERNAL_VALUES(1:N,D)
       END DO
    ! Otherwise, assume regular output computation.
    ELSE
       CALL GEMM('T', 'T', MDO, N, MDS, 1.0_RT, &
            OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
            INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
            0.0_RT, OUTPUTS(:,:), SIZE(OUTPUTS,1))
       DO D = 1, MDO
          OUTPUTS(D,1:N) = OUTPUT_SHIFT(D) + OUTPUTS(D,1:N)
       END DO
    END IF
    ! Deallocate temporary variables.
    DEALLOCATE(INTERNAL_VALUES, TEMP_VALUES)
  END SUBROUTINE EVALUATE


  ! Compute the gradient of the sum of squared error of this regression
  ! model with respect to its variables given input and output pairs.
  SUBROUTINE SSE_GRADIENT(INPUTS, OUTPUTS, SUM_SQUARED_ERROR, &
       INPUT_VECS_GRADIENT, INPUT_SHIFT_GRADIENT, &
       INTERNAL_VECS_GRADIENT, INTERNAL_SHIFT_GRADIENT, &
       OUTPUT_VECS_GRADIENT, OUTPUT_SHIFT_GRADIENT)
    ! Data values stored with contiguous points: shape = (MDI,N)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: OUTPUTS
    ! Sum (over all data) squared error (summed over dimensions).
    REAL(KIND=RT), INTENT(INOUT) :: SUM_SQUARED_ERROR
    ! Model variable gradients.
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:)   :: INPUT_VECS_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:)     :: INPUT_SHIFT_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:,:) :: INTERNAL_VECS_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:)   :: INTERNAL_SHIFT_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:)   :: OUTPUT_VECS_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:)     :: OUTPUT_SHIFT_GRADIENT     
    ! Internal values and lists of which ones are nonzero.
    REAL(KIND=RT), DIMENSION(SIZE(INPUTS,2),MDO)     :: OUTPUT_GRADIENT
    REAL(KIND=RT), DIMENSION(SIZE(INPUTS,2),MDS,MNS) :: INTERNAL_VALUES
    LOGICAL,       DIMENSION(SIZE(INPUTS,2),MDS,MNS) :: LT_DISCONTINUITY
    REAL(KIND=RT), DIMENSION(SIZE(INPUTS,2),MDS)     :: INTERNAL_TEMP
    ! Number of points.
    REAL(KIND=RT) :: N
    INTEGER :: D, L
    ! Evaluate the model at all data points, store outputs in "OUTPUT_GRADIENT"
    CALL EVALUATE_BATCH()
    ! Compute the new gradient.
    DO D = 1, MDO
       OUTPUT_GRADIENT(:,D) = OUTPUT_GRADIENT(:,D) - OUTPUTS(D,:)
    END DO
    ! Get the number of points.
    N = REAL(SIZE(INPUTS,2), RT)
    ! Compute the current sum squared error.
    SUM_SQUARED_ERROR = SUM_SQUARED_ERROR + SUM(OUTPUT_GRADIENT(:,:)**2)
    ! Compute the gradient of variables with respect to the output gradient.
    CALL GRADIENT_BATCH()

  CONTAINS

    ! Evaluate the piecewise linear regression model.
    SUBROUTINE EVALUATE_BATCH()
      ! D   - dimension index
      ! L   - layer index
      ! LP1 - layer index "plus 1" -> "P1"
      INTEGER :: D, L, LP1
      ! Compute the input transformation.
      CALL GEMM('T', 'N', SIZE(INPUTS,2), MDS, MDI, 1.0_RT, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INPUT_VECS(:,:), SIZE(INPUT_VECS,1), &
           0.0_RT, INTERNAL_VALUES(:,:,1), SIZE(INTERNAL_VALUES,1))
      DO D = 1, MDS
         INTERNAL_VALUES(:,D,1) = MAX( &
              INTERNAL_VALUES(:,D,1) + INPUT_SHIFT(D), &
              DISCONTINUITY)
      END DO
      ! Compute the next set of internal values with a rectified activation.
      DO L = 1, MNS-1
         ! Compute vector values.
         LP1 = L+1
         CALL GEMM('N', 'N', SIZE(INPUTS,2), MDS, MDS, 1.0_RT, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
              0.0_RT, INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1))
         DO D = 1, MDS
            INTERNAL_VALUES(:,D,LP1) = MAX( &
                 INTERNAL_VALUES(:,D,LP1) + INTERNAL_SHIFT(D,L), &
                 DISCONTINUITY)
         END DO
      END DO
      ! Compute the output.
      CALL GEMM('N', 'N', SIZE(INPUTS,2), MDO, MDS, 1.0_RT, &
           INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
           0.0_RT, OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1))
      DO D = 1, MDO
         OUTPUT_GRADIENT(:,D) = OUTPUT_SHIFT(D) + OUTPUT_GRADIENT(:,D)
      END DO
    END SUBROUTINE EVALUATE_BATCH

    SUBROUTINE GRADIENT_BATCH()
      ! D   - dimension index
      ! L   - layer index
      ! LP1 - layer index "plus 1" -> "P1"
      INTEGER :: D, L, LP1
      ! Compute the average output gradient.
      OUTPUT_SHIFT_GRADIENT(:) = SUM(OUTPUT_GRADIENT(:,:), 1) / N &
           + OUTPUT_SHIFT_GRADIENT(:)
      ! Compute output variable gradients.
      CALL GEMM('T', 'N', MDS, MDO, SIZE(INPUTS,2), 1.0_RT / N, &
           INTERNAL_VALUES(:,:,MNS), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           1.0_RT, OUTPUT_VECS_GRADIENT(:,:), SIZE(OUTPUT_VECS_GRADIENT,1))
      ! Propogate the gradient back to the last internal vector space.
      CALL GEMM('N', 'T', SIZE(INPUTS,2), MDS, MDO, 1.0_RT, &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
           0.0_RT, INTERNAL_TEMP(:,:), SIZE(INTERNAL_TEMP,1))
      ! Cycle over all internal layers.
      INTERNAL_REPRESENTATIONS : DO L = MNS-1, 1, -1
         LP1 = L+1
         DO D = 1, MDS
            ! Propogate the error gradient back to the preceding vectors.
            WHERE (INTERNAL_VALUES(:,D,LP1) .GT. DISCONTINUITY)
               INTERNAL_VALUES(:,D,LP1) = INTERNAL_TEMP(:,D)
            END WHERE
         END DO
         ! Compute the shift gradient.
         INTERNAL_SHIFT_GRADIENT(:,L) = SUM(INTERNAL_VALUES(:,:,LP1), 1) / N &
              + INTERNAL_SHIFT_GRADIENT(:,L)
         ! Compute the gradient with respect to each output and all inputs.
         CALL GEMM('T', 'N', MDS, MDS, SIZE(INPUTS,2), 1.0_RT / N, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              1.0_RT, INTERNAL_VECS_GRADIENT(:,:,L), SIZE(INTERNAL_VECS_GRADIENT,1))
         ! Propogate the gradient to the immediately preceding layer.
         CALL GEMM('N', 'T', SIZE(INPUTS,2), MDS, MDS, 1.0_RT, &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
              0.0_RT, INTERNAL_TEMP(:,:), SIZE(INTERNAL_TEMP,1))
      END DO INTERNAL_REPRESENTATIONS
      ! Compute the gradients going into the first layer.
      DO D = 1, MDS
         ! Propogate the error gradient back to the preceding vectors.
         WHERE (INTERNAL_VALUES(:,D,1) .GT. DISCONTINUITY)
            INTERNAL_VALUES(:,D,1) = INTERNAL_TEMP(:,D)
         END WHERE
      END DO
      ! Compute the input shift variable gradients.
      INPUT_SHIFT_GRADIENT(:) = SUM(INTERNAL_VALUES(:,:,1), 1) / N &
           + INPUT_SHIFT_GRADIENT(:)
      ! Compute the gradient of all input variables.
      !   [the INPUTS are transposed already, shape = (MDI,N)]
      CALL GEMM('N', 'N', MDI, MDS, SIZE(INPUTS,2), 1.0_RT / N, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INTERNAL_VALUES(:,:,1), SIZE(INTERNAL_VALUES,1), &
           1.0_RT, INPUT_VECS_GRADIENT(:,:), SIZE(INPUT_VECS_GRADIENT,1))
    END SUBROUTINE GRADIENT_BATCH

  END SUBROUTINE SSE_GRADIENT


  ! Fit input / output pairs by minimizing mean squared error.
  SUBROUTINE MINIMIZE_MSE(X, Y, STEPS, BATCH_SIZE, NUM_THREADS, &
       STEP_SIZE, KEEP_BEST, MEAN_SQUARED_ERROR)
    ! Define an interface for the "LOG_STATUS" function.
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: X
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: Y
    INTEGER,       INTENT(IN) :: STEPS
    INTEGER,       INTENT(IN), OPTIONAL :: BATCH_SIZE, NUM_THREADS
    REAL(KIND=RT), INTENT(IN), OPTIONAL :: STEP_SIZE
    LOGICAL,       INTENT(IN), OPTIONAL :: KEEP_BEST
    REAL(KIND=RT), INTENT(OUT) :: MEAN_SQUARED_ERROR
    !  gradient step arrays, 4 copies of internal model (total of 5).
    REAL(KIND=RT), DIMENSION(MDI,MDS)       :: V1_GRAD, V1_MEAN, V1_CURV, BEST_V1
    REAL(KIND=RT), DIMENSION(MDS)           :: S1_GRAD, S1_MEAN, S1_CURV, BEST_S1
    REAL(KIND=RT), DIMENSION(MDS,MDS,MNS-1) :: V2_GRAD, V2_MEAN, V2_CURV, BEST_V2
    REAL(KIND=RT), DIMENSION(MDS,MNS-1)     :: S2_GRAD, S2_MEAN, S2_CURV, BEST_S2
    REAL(KIND=RT), DIMENSION(MDS,MDO)       :: V3_GRAD, V3_MEAN, V3_CURV, BEST_V3
    REAL(KIND=RT), DIMENSION(MDO)           :: S3_GRAD, S3_MEAN, S3_CURV, BEST_S3
    !  local variables
    LOGICAL :: REVERT_TO_BEST
    LOGICAL, DIMENSION(MDS,MNS) :: TO_ZERO
    INTEGER :: BS, I, L, NX, NT, BATCH_START, BATCH_END
    REAL(KIND=RT) :: RNX, BATCHES, PREV_MSE, BEST_MSE, SCALAR
    REAL(KIND=RT) :: STEP_FACTOR, STEP_MEAN_CHANGE, STEP_MEAN_REMAIN, &
         STEP_CURV_CHANGE, STEP_CURV_REMAIN
    ! Function that is defined by OpenMP.
    INTERFACE
       FUNCTION OMP_GET_MAX_THREADS()
         INTEGER :: OMP_GET_MAX_THREADS
       END FUNCTION OMP_GET_MAX_THREADS
    END INTERFACE
    ! Number of points.
    NX = SIZE(X,2)
    RNX = REAL(NX, RT)
    ! Set the number of threads.
    IF (PRESENT(NUM_THREADS)) THEN
       NT = NUM_THREADS
    ELSE
       NT = OMP_GET_MAX_THREADS()
    END IF
    ! Set the batch size.
    IF (PRESENT(BATCH_SIZE)) THEN
       BS = MIN(BATCH_SIZE, NX)
    ELSE
       BS = (NX + NT - 1) / NT
    END IF
    ! Set the step factor.
    IF (PRESENT(STEP_SIZE)) THEN
       STEP_FACTOR = STEP_SIZE
    ELSE
       STEP_FACTOR = INITIAL_STEP
    END IF
    ! Set the "keep best" boolean.
    IF (PRESENT(KEEP_BEST)) THEN
       REVERT_TO_BEST = KEEP_BEST
    ELSE
       REVERT_TO_BEST = (STEPS .GE. MIN_STEPS_TO_STABILITY)
    END IF
    ! Initial rates of change of mean and variance values.
    STEP_MEAN_CHANGE = INITIAL_STEP_MEAN_CHANGE
    STEP_MEAN_REMAIN = 1.0_RT - STEP_MEAN_CHANGE
    STEP_CURV_CHANGE = INITIAL_STEP_CURV_CHANGE
    STEP_CURV_REMAIN = 1.0_RT - STEP_CURV_CHANGE
    ! Initial mean squared error is "max representable value".
    PREV_MSE = HUGE(PREV_MSE)
    BEST_MSE = HUGE(BEST_MSE)
    ! Set the average step sizes.
    V1_MEAN(:,:)   = 0.0_RT
    S1_MEAN(:)     = 0.0_RT
    V2_MEAN(:,:,:) = 0.0_RT
    S2_MEAN(:,:)   = 0.0_RT
    V3_MEAN(:,:)   = 0.0_RT
    S3_MEAN(:)     = 0.0_RT
    ! Set the estiamted curvature in steps.
    V1_CURV(:,:)   = 0.0_RT
    S1_CURV(:)     = 0.0_RT
    V2_CURV(:,:,:) = 0.0_RT
    S2_CURV(:,:)   = 0.0_RT
    V3_CURV(:,:)   = 0.0_RT
    S3_CURV(:)     = 0.0_RT
    ! Iterate, taking steps with the average gradient over all data.
    fit_loop : DO I = 1, STEPS
       ! Compute the average gradient over all points.
       MEAN_SQUARED_ERROR = 0.0_RT
       ! Set gradients to zero initially.
       V1_GRAD(:,:)   = 0.0_RT
       S1_GRAD(:)     = 0.0_RT
       V2_GRAD(:,:,:) = 0.0_RT
       S2_GRAD(:,:)   = 0.0_RT
       V3_GRAD(:,:)   = 0.0_RT
       S3_GRAD(:)     = 0.0_RT
       ! Count the number of batches.
       BATCHES = 0.0_RT
       !$OMP PARALLEL DO NUM_THREADS(NT) PRIVATE(BATCH_END) &
       !$OMP&  REDUCTION(+: BATCHES, MEAN_SQUARED_ERROR, &
       !$OMP&  V1_GRAD, S1_GRAD, V2_GRAD, S2_GRAD, V3_GRAD, S3_GRAD)
       DO BATCH_START = 1, NX, BS
          BATCHES = BATCHES + 1.0_RT
          BATCH_END = MIN(NX, BATCH_START+BS-1)
          ! Sum the gradient over all data batches.
          CALL SSE_GRADIENT(X(:,BATCH_START:BATCH_END), &
               Y(:,BATCH_START:BATCH_END), MEAN_SQUARED_ERROR, &
               V1_GRAD(:,:), S1_GRAD(:), V2_GRAD(:,:,:), S2_GRAD(:,:), &
               V3_GRAD(:,:), S3_GRAD(:))
       END DO
       !$OMP END PARALLEL DO
       ! Convert the sum of squared errors into the mean squared error.
       MEAN_SQUARED_ERROR = MEAN_SQUARED_ERROR / RNX

       ! Update the saved "best" model based on error (only when dropout is disabled).
       IF (MEAN_SQUARED_ERROR .LT. BEST_MSE) THEN
          BEST_MSE       = MEAN_SQUARED_ERROR
          BEST_V1(:,:)   = INPUT_VECS(:,:)
          BEST_S1(:)     = INPUT_SHIFT(:)
          BEST_V2(:,:,:) = INTERNAL_VECS(:,:,:)
          BEST_S2(:,:)   = INTERNAL_SHIFT(:,:)
          BEST_V3(:,:)   = OUTPUT_VECS(:,:)
          BEST_S3(:)     = OUTPUT_SHIFT(:)
       END IF

       ! Update the step factor based on model improvement.
       IF (MEAN_SQUARED_ERROR .LE. PREV_MSE) THEN
          STEP_FACTOR = STEP_FACTOR * STEP_GROWTH_RATE
          STEP_MEAN_CHANGE = STEP_MEAN_CHANGE * STEP_GROWTH_RATE
          STEP_MEAN_REMAIN = 1.0_RT - STEP_MEAN_CHANGE
          STEP_CURV_CHANGE = STEP_CURV_CHANGE * STEP_GROWTH_RATE
          STEP_CURV_REMAIN = 1.0_RT - STEP_CURV_CHANGE
       ELSE
          STEP_FACTOR = STEP_FACTOR * STEP_SHRINK_RATE
          STEP_MEAN_CHANGE = STEP_MEAN_CHANGE * STEP_SHRINK_RATE
          STEP_MEAN_REMAIN = 1.0_RT - STEP_MEAN_CHANGE
          STEP_CURV_CHANGE = STEP_CURV_CHANGE * STEP_SHRINK_RATE
          STEP_CURV_REMAIN = 1.0_RT - STEP_CURV_CHANGE
       END IF

       ! Store the previous error for tracking the best-so-far.
       PREV_MSE = MEAN_SQUARED_ERROR
       
       ! Update the steps for all of the model variables.
       CALL COMPUTE_STEP(SIZE(INPUT_VECS(:,:)),      INPUT_VECS(:,:),      V1_GRAD(:,:),   V1_MEAN(:,:),   V1_CURV(:,:))
       CALL COMPUTE_STEP(SIZE(INPUT_SHIFT(:)),       INPUT_SHIFT(:),       S1_GRAD(:),     S1_MEAN(:),     S1_CURV(:))
       CALL COMPUTE_STEP(SIZE(INTERNAL_VECS(:,:,:)), INTERNAL_VECS(:,:,:), V2_GRAD(:,:,:), V2_MEAN(:,:,:), V2_CURV(:,:,:))
       CALL COMPUTE_STEP(SIZE(INTERNAL_SHIFT(:,:)),  INTERNAL_SHIFT(:,:),  S2_GRAD(:,:),   S2_MEAN(:,:),   S2_CURV(:,:))
       CALL COMPUTE_STEP(SIZE(OUTPUT_VECS(:,:)),     OUTPUT_VECS(:,:),     V3_GRAD(:,:),   V3_MEAN(:,:),   V3_CURV(:,:))
       CALL COMPUTE_STEP(SIZE(OUTPUT_SHIFT(:)),      OUTPUT_SHIFT(:),      S3_GRAD(:),     S3_MEAN(:),     S3_CURV(:))

       ! Maintain a constant max-norm across the magnitue of input and internal vectors.
       SCALAR = SQRT(MAXVAL(SUM(INPUT_VECS(:,:)**2, 1)))
       INPUT_VECS(:,:) = INPUT_VECS(:,:) / SCALAR
       DO L = 1, MNS-1
          SCALAR = SQRT(MAXVAL(SUM(INTERNAL_VECS(:,:,L)**2, 1)))
          INTERNAL_VECS(:,:,L) = INTERNAL_VECS(:,:,L) / SCALAR
       END DO

    END DO fit_loop

    ! Restore the best model seen so far (if enough steps were taken).
    IF (REVERT_TO_BEST) THEN
       MEAN_SQUARED_ERROR   = BEST_MSE
       INPUT_VECS(:,:)      = BEST_V1(:,:)  
       INPUT_SHIFT(:)       = BEST_S1(:)    
       INTERNAL_VECS(:,:,:) = BEST_V2(:,:,:)
       INTERNAL_SHIFT(:,:)  = BEST_S2(:,:)  
       OUTPUT_VECS(:,:)     = BEST_V3(:,:)  
       OUTPUT_SHIFT(:)      = BEST_S3(:)    
    END IF

  CONTAINS

    ! Compute the current step factor as the running average step,
    ! inversely scaled by the average magnitude of the step.
    ! 
    ! Implicit arguments to this subroutine from parent scope:
    !   BATCHES
    !   STEP_MEAN_CHANGE
    !   STEP_MEAN_REMAIN
    !   STEP_CURV_CHANGE
    !   STEP_CURV_REMAIN
    !   STEP_FACTOR
    SUBROUTINE COMPUTE_STEP(N, PARAMS, CURR_STEP, STEP_MEAN, STEP_CURV)
      INTEGER, INTENT(IN) :: N
      REAL(KIND=RT), DIMENSION(1:N), INTENT(INOUT) :: PARAMS
      REAL(KIND=RT), DIMENSION(1:N), INTENT(INOUT) :: CURR_STEP
      REAL(KIND=RT), DIMENSION(1:N), INTENT(INOUT) :: STEP_MEAN
      REAL(KIND=RT), DIMENSION(1:N), INTENT(INOUT) :: STEP_CURV
      ! Convert the summed gradients to average gradients.
      CURR_STEP(:) = CURR_STEP(:) / BATCHES
      ! Update sliding window mean and variance calculations.
      STEP_MEAN(:) = STEP_MEAN_CHANGE * CURR_STEP(:) &
           + STEP_MEAN_REMAIN * STEP_MEAN(:)
      STEP_CURV(:) = STEP_CURV_CHANGE * (STEP_MEAN(:) - CURR_STEP(:))**2 &
           + STEP_CURV_REMAIN * STEP_CURV(:)
      STEP_CURV(:) = MAX(STEP_CURV(:), EPSILON(STEP_FACTOR))
      ! Take averaged gradient descent step.
      CURR_STEP(:) = STEP_FACTOR * STEP_MEAN(:)
      ! Start scaling by step magnitude by curvature once enough data is collected.
      IF (I .GT. MIN_STEPS_TO_STABILITY) THEN
         CURR_STEP(:) = CURR_STEP(:) / SQRT(STEP_CURV(:))
      END IF
      ! Take the gradient steps (based on the computed "step" above).
      PARAMS(:) = PARAMS(:) - CURR_STEP(:)
    END SUBROUTINE COMPUTE_STEP

  END SUBROUTINE MINIMIZE_MSE

END MODULE PLRM
