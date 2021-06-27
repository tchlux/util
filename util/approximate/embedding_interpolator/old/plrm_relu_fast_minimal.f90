! TODO:
! - Study convergence pattern, it's not as reliable as I hoped.
! - Compare loss per step and loss per unit time with old version.
! - Implement custom block matrix multiplication routines.
! - Training logs.
! - Save model and load model (during training).
! - Well-spaced validation data (multiply distances in input and output).
! - Validation-based training (and greedy saving).
! - Embedding support, for integer valued inputs.
! - Embedding support for vector of assigned neighbors for each point.
! - Train only the output, internal, or input layers.
! - Test larger networks on image data.
! 
! MAYBE TODO:
! - Use ALLOCATE instead of automatic, because of large N seg faults.
! - Implement similar code in C, compare speeds.


! A piecewise linear regression model.
MODULE PLRM
  USE ISO_FORTRAN_ENV, RT => REAL32
  IMPLICIT NONE

  ! MDI = Model dimension -- input
  ! MDS = Model dimension -- states (internal)
  ! MNS = Model number -- states (internal)
  ! MDO = Model dimension -- output
  INTEGER :: MDI, MDS, MNS, MDO
  ! Model parameters, for evaluation.
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INPUT_VECS
  REAL(KIND=RT), DIMENSION(:),     ALLOCATABLE :: INPUT_SHIFT
  REAL(KIND=RT), DIMENSION(:,:,:), ALLOCATABLE :: INTERNAL_VECS
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_SHIFT
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: OUTPUT_VECS
  
  PUBLIC :: NEW_MODEL, INIT_MODEL, MINIMIZE_MSE, EVALUATE_ONE
  PRIVATE :: SSE_GRADIENT

  ! Externally provided matrix-matrix multiplication routine.
  EXTERNAL :: SGEMM

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
    ! Allocate space for all model parameters.
    ALLOCATE( &
         INPUT_VECS(DI,DS), &
         INPUT_SHIFT(DS), &
         INTERNAL_VECS(DS, DS, NS-1), &
         INTERNAL_SHIFT(DS, NS-1), &
         OUTPUT_VECS(DS, DO) &
         )
  END SUBROUTINE NEW_MODEL


  ! Initialize the weights for a model, optionally provide a random seed.
  SUBROUTINE INIT_MODEL(SEED)
    INTEGER, INTENT(IN), OPTIONAL :: SEED
    !  Storage for seeding the random number generator (for repeatability).
    INTEGER, DIMENSION(:), ALLOCATABLE :: SEED_ARRAY
    INTEGER :: I
    ! Set a random seed, if one was provided (otherwise leave default).
    IF (PRESENT(SEED)) THEN
       CALL RANDOM_SEED(SIZE=I)
       ALLOCATE(SEED_ARRAY(I))
       SEED_ARRAY(:) = SEED
       CALL RANDOM_SEED(PUT=SEED_ARRAY(:))
    END IF
    ! Generate random unit-length vectors (no scaling shiftes) for
    !  all initial parameters in the input, internal, and output.
    CALL RANDOM_UNIT_VECTORS(INPUT_VECS(:,:))
    DO I = 1, MNS-1
       CALL RANDOM_UNIT_VECTORS(INTERNAL_VECS(:,:,I))
    END DO
    CALL RANDOM_UNIT_VECTORS(OUTPUT_VECS(:,:))
    ! Generate random shiftes for inputs and internal layers, zero
    !  shift for the output layer (first two will be rescaled).
    CALL RANDOM_NUMBER(INPUT_SHIFT(:))
    CALL RANDOM_NUMBER(INTERNAL_SHIFT(:,:))
    ! Values are in range [-2,2], assuming no scaling changes,
    !  this will capture 2 standard deviations in either direction.
    INPUT_SHIFT(:) = 2.0_RT * (INPUT_SHIFT(:) * 2.0_RT - 1.0_RT)
    INTERNAL_SHIFT(:,:) = 2.0_RT * (INTERNAL_SHIFT(:,:) * 2.0_RT - 1.0_RT)

  CONTAINS

    ! Generate randomly distributed vectors on the N-sphere.
    SUBROUTINE RANDOM_UNIT_VECTORS(COLUMN_VECTORS)
      REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: COLUMN_VECTORS
      INTEGER :: I 
      ! Generate random numbers in the range [0,1].
      CALL RANDOM_NUMBER(COLUMN_VECTORS(:,:))
      ! Map those random numbers to the range [-1,1].
      COLUMN_VECTORS(:,:) = 2.0_RT * COLUMN_VECTORS(:,:) - 1.0_RT
      ! Make the vectors uniformly distributed on the unit ball (for dimension > 1).
      IF (SIZE(COLUMN_VECTORS,1) .GT. 1) THEN
         ! Compute the inverse hyperbolic tangent of all values.
         COLUMN_VECTORS(:,:) = ATANH(COLUMN_VECTORS(:,:))
         ! Normalize all vectors to have unit length.
         DO I = 1, SIZE(COLUMN_VECTORS,2)
            COLUMN_VECTORS(:,I) = COLUMN_VECTORS(:,I) / NORM2(COLUMN_VECTORS(:,I))
         END DO
      END IF
      ! Orthogonalize the first components of the column
      !  vectors to ensure those are well spaced.
      I = MIN(SIZE(COLUMN_VECTORS,1), SIZE(COLUMN_VECTORS,2))
      CALL ORTHOGONALIZE(COLUMN_VECTORS(:,1:I))
    END SUBROUTINE RANDOM_UNIT_VECTORS

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

  END SUBROUTINE INIT_MODEL


  ! Evaluate the piecewise linear regression model.
  SUBROUTINE EVALUATE_ONE(INPUT, OUTPUT)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:) :: INPUT
    ! Record of values after the last transformation.
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:) :: OUTPUT
    ! Internal values.
    REAL(KIND=RT), DIMENSION(MDS) :: INTERNAL_VALUES
    ! Compute the number of nodes (to know where "shift" value is).
    INTEGER :: I, IP1
    ! Make sure that this model has been initialized.
    IF (.NOT. ALLOCATED(INPUT_VECS)) THEN
       OUTPUT(:) = 0.0
       RETURN
    END IF
    ! Compute the input layer.
    INTERNAL_VALUES(:) = MAX(0.0_RT, INPUT_SHIFT(:) + MATMUL(INPUT(:), INPUT_VECS(:,:)))
    ! Compute the next set of internal values with a rectified activation.
    DO I = 1, MNS-1
       IP1 = I+1
       INTERNAL_VALUES(:) = MAX(0.0_RT, INTERNAL_SHIFT(:,I) + &
            MATMUL(INTERNAL_VALUES(:), INTERNAL_VECS(:,:,I)))
    END DO
    ! Compute the output.
    OUTPUT(:) = MATMUL(INTERNAL_VALUES(:), OUTPUT_VECS(:,:))
  END SUBROUTINE EVALUATE_ONE

  ! Evaluate the piecewise linear regression model.
  SUBROUTINE EVALUATE(INPUTS, OUTPUTS)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: OUTPUTS
    ! Internal values.
    REAL(KIND=RT), DIMENSION(:,:), ALLOCATABLE :: TEMP_VALUES, INTERNAL_VALUES
    ! Compute the number of nodes (to know where "shift" value is).
    INTEGER :: D, L, N
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
    CALL SGEMM('T', 'N', N, MDS, MDI, 1.0_RT, &
         INPUTS(:,:), SIZE(INPUTS,1), &
         INPUT_VECS(:,:), SIZE(INPUT_VECS,1), &
         0.0_RT, INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1))
    DO D = 1, MDS
       INTERNAL_VALUES(:,D) = MAX(0.0_RT, INPUT_SHIFT(D) + INTERNAL_VALUES(:,D))
    END DO
    ! Compute the next set of internal values with a rectified activation.
    DO L = 1, MNS-1
       CALL SGEMM('N', 'N', N, MDS, MDS, 1.0_RT, &
            INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
            INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
            0.0_RT, TEMP_VALUES(:,:), SIZE(TEMP_VALUES,1))
       INTERNAL_VALUES(:,:) = TEMP_VALUES(:,:)
       DO D = 1, MDS
          INTERNAL_VALUES(:,D) =  MAX(0.0, INTERNAL_SHIFT(D,L) + INTERNAL_VALUES(:,D))
       END DO
    END DO
    ! Return the final embedded layer if the outputs have the size of the embedding.
    IF ((MDS .NE. MDO) .AND. (SIZE(OUTPUTS,1) .EQ. MDS)) THEN
       ! Do the necessary transpose operation one dimension at a time.
       DO L = 1, MDS
          OUTPUTS(L,1:N) = INTERNAL_VALUES(1:N,L)
       END DO
    ! Otherwise, assume regular output computation.
    ELSE
       CALL SGEMM('T', 'T', MDO, N, MDS, 1.0_RT, &
            OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
            INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
            0.0_RT, OUTPUTS(:,:), SIZE(OUTPUTS,1))
    END IF
    ! Deallocate temporary variables.
    DEALLOCATE(INTERNAL_VALUES, TEMP_VALUES)
  END SUBROUTINE EVALUATE


  ! Compute the gradient of the sum of squared error of this regression
  ! model with respect to its parameters given input and output pairs.
  SUBROUTINE SSE_GRADIENT(INPUTS, OUTPUTS, SUM_SQUARED_ERROR, &
       INPUT_VECS_GRADIENT, INPUT_SHIFT_GRADIENT, &
       INTERNAL_VECS_GRADIENT, INTERNAL_SHIFT_GRADIENT, &
       OUTPUT_VECS_GRADIENT)
    ! Data values stored with contiguous points: shape = (MDI,N)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: OUTPUTS
    ! Sum (over all data) squared error (summed over dimensions).
    REAL(KIND=RT), INTENT(INOUT) :: SUM_SQUARED_ERROR
    ! Model parameter gradients.
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:)   :: INPUT_VECS_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:)     :: INPUT_SHIFT_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:,:) :: INTERNAL_VECS_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:)   :: INTERNAL_SHIFT_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:)   :: OUTPUT_VECS_GRADIENT
    ! Internal values and lists of which ones are nonzero.
    REAL(KIND=RT), DIMENSION(SIZE(INPUTS,2),MDO)     :: OUTPUT_GRADIENT
    REAL(KIND=RT), DIMENSION(SIZE(INPUTS,2),MDS,MNS) :: INTERNAL_VALUES
    LOGICAL, DIMENSION(SIZE(INPUTS,2),MDS) :: LT_DISCONTINUITY
    ! Number of points.
    REAL(KIND=RT) :: N
    INTEGER :: D
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
    ! Compute the gradient of parameters with respect to the output gradient.
    CALL GRADIENT_BATCH()

  CONTAINS

    ! Evaluate the piecewise linear regression model.
    SUBROUTINE EVALUATE_BATCH()
      ! D   - dimension index
      ! L   - layer index
      ! LP1 - layer index "plus 1" -> "P1"
      INTEGER :: D, L, LP1
      ! Compute the input transformation.
      CALL SGEMM('T', 'N', SIZE(INPUTS,2), MDS, MDI, 1.0_RT, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INPUT_VECS(:,:), SIZE(INPUT_VECS,1), &
           0.0_RT, INTERNAL_VALUES(:,:,1), SIZE(INTERNAL_VALUES,1))
      DO D = 1, MDS
         INTERNAL_VALUES(:,D,1) = MAX(0.0_RT, INPUT_SHIFT(D) + INTERNAL_VALUES(:,D,1))
      END DO
      ! Compute the next set of internal values with a rectified activation.
      DO L = 1, MNS-1
         LP1 = L+1
         CALL SGEMM('N', 'N', SIZE(INPUTS,2), MDS, MDS, 1.0_RT, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
              0.0_RT, INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1))
         DO D = 1, MDS
            INTERNAL_VALUES(:,D,LP1) =  MAX(0.0_RT, INTERNAL_SHIFT(D,L) + INTERNAL_VALUES(:,D,LP1))
         END DO
      END DO
      ! Compute the output.
      CALL SGEMM('N', 'N', SIZE(INPUTS,2), MDO, MDS, 1.0_RT, &
           INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
           0.0_RT, OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1))
    END SUBROUTINE EVALUATE_BATCH

    SUBROUTINE GRADIENT_BATCH()
      ! L   - layer index
      ! LP1 - layer index "plus 1" -> "P1"
      INTEGER :: L, LP1
      ! Compute output parameter gradients.
      CALL SGEMM('T', 'N', MDS, MDO, SIZE(INPUTS,2), 1.0_RT / N, &
           INTERNAL_VALUES(:,:,MNS), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           1.0_RT, OUTPUT_VECS_GRADIENT(:,:), SIZE(OUTPUT_VECS_GRADIENT,1))
      ! Compute the gradient for the last internal vector space values.
      LT_DISCONTINUITY(:,:) = INTERNAL_VALUES(:,:,MNS) .LT. 0.0_RT
      CALL SGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDO, 1.0_RT, &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
           0.0_RT, INTERNAL_VALUES(:,:,MNS), SIZE(INTERNAL_VALUES,1))
      WHERE (LT_DISCONTINUITY(:,:))
         INTERNAL_VALUES(:,:,MNS) = 0.0_RT
      END WHERE
      ! Cycle over all internal layers.
      INTERNAL_REPRESENTATIONS : DO L = MNS-1, 1, -1
         LP1 = L+1
         ! Compute the shift gradient.
         INTERNAL_SHIFT_GRADIENT(:,L) = SUM(INTERNAL_VALUES(:,:,LP1), 1) / N &
              + INTERNAL_SHIFT_GRADIENT(:,L)
         ! Compute the gradient with respect to each output and all inputs.
         CALL SGEMM('T', 'N', MDS, MDS, SIZE(INPUTS,2), 1.0_RT / N, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              1.0_RT, INTERNAL_VECS_GRADIENT(:,:,L), SIZE(INTERNAL_VECS_GRADIENT,1))
         ! Propogate the gradient to the immediately preceding layer.
         LT_DISCONTINUITY(:,:) = INTERNAL_VALUES(:,:,L) .LT. 0.0_RT
         CALL SGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDS, 1.0_RT, &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
              0.0_RT, INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1))
         WHERE (LT_DISCONTINUITY(:,:))
            INTERNAL_VALUES(:,:,L) = 0.0_RT
         END WHERE
      END DO INTERNAL_REPRESENTATIONS
      ! Compute the input parameter shift values.
      INPUT_SHIFT_GRADIENT(:) = SUM(INTERNAL_VALUES(:,:,1), 1) / N &
           + INPUT_SHIFT_GRADIENT(:)
      ! Compute the gradient of all input parameters.
      !   [the INPUTS are transposed already, shape = (MDI,N)]
      CALL SGEMM('N', 'N', MDI, MDS, SIZE(INPUTS,2), 1.0_RT / N, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INTERNAL_VALUES(:,:,1), SIZE(INTERNAL_VALUES,1), &
           1.0_RT, INPUT_VECS_GRADIENT(:,:), SIZE(INPUT_VECS_GRADIENT,1))
    END SUBROUTINE GRADIENT_BATCH

  END SUBROUTINE SSE_GRADIENT


  ! Fit input / output pairs by minimizing mean squared error.
  SUBROUTINE MINIMIZE_MSE(X, Y, STEPS, STEP_SIZE, BATCH_SIZE, & 
       NUM_THREADS, MEAN_SQUARED_ERROR, RECORD)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: X
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: Y
    INTEGER,       INTENT(IN) :: STEPS
    REAL(KIND=RT), INTENT(IN), OPTIONAL :: STEP_SIZE
    INTEGER,       INTENT(IN), OPTIONAL :: BATCH_SIZE, NUM_THREADS
    REAL(KIND=RT), INTENT(OUT) :: MEAN_SQUARED_ERROR
    REAL(KIND=RT), INTENT(OUT), DIMENSION(STEPS), OPTIONAL :: RECORD
    !  gradient step arrays, 4 copies of internal model (total of 5).
    REAL(KIND=RT), DIMENSION(MDI,MDS)       :: V1_GRAD, V1_MEAN, V1_VAR, BEST_V1
    REAL(KIND=RT), DIMENSION(MDS)           :: S1_GRAD, S1_MEAN, S1_VAR, BEST_S1
    REAL(KIND=RT), DIMENSION(MDS,MDS,MNS-1) :: V2_GRAD, V2_MEAN, V2_VAR, BEST_V2
    REAL(KIND=RT), DIMENSION(MDS,MNS-1)     :: S2_GRAD, S2_MEAN, S2_VAR, BEST_S2
    REAL(KIND=RT), DIMENSION(MDS,MDO)       :: V3_GRAD, V3_MEAN, V3_VAR, BEST_V3
    !  local integers
    INTEGER :: BS, I, L, NX, NT, BATCH_START, BATCH_END
    REAL(KIND=RT) :: RNX, BATCHES, STEP_FACTOR, SCALAR
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
       STEP_FACTOR = 0.001
    END IF
    ! Set the average step sizes.
    V1_MEAN(:,:)   = 0.0_RT
    S1_MEAN(:)     = 0.0_RT
    V2_MEAN(:,:,:) = 0.0_RT
    S2_MEAN(:,:)   = 0.0_RT
    V3_MEAN(:,:)   = 0.0_RT
    ! Set the average step size variances.
    V1_VAR(:,:)   = 0.0_RT
    S1_VAR(:)     = 0.0_RT
    V2_VAR(:,:,:) = 0.0_RT
    S2_VAR(:,:)   = 0.0_RT
    V3_VAR(:,:)   = 0.0_RT
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
       ! Count the number of batches.
       BATCHES = 0.0_RT
       !$OMP PARALLEL DO NUM_THREADS(NT) PRIVATE(BATCH_END) &
       !$OMP&  REDUCTION(+: BATCHES, MEAN_SQUARED_ERROR, &
       !$OMP&  V1_GRAD, S1_GRAD, V2_GRAD, S2_GRAD, V3_GRAD)
       DO BATCH_START = 1, NX, BS
          BATCHES = BATCHES + 1.0_RT
          BATCH_END = MIN(NX, BATCH_START+BS-1)
          ! Compute the gradient from all data.
          CALL SSE_GRADIENT(X(:,BATCH_START:BATCH_END), &
               Y(:,BATCH_START:BATCH_END), MEAN_SQUARED_ERROR, &
               V1_GRAD, S1_GRAD, V2_GRAD, S2_GRAD, V3_GRAD)
       END DO
       !$OMP END PARALLEL DO
       ! Convert the sum of squared errors into the mean squared error.
       MEAN_SQUARED_ERROR = MEAN_SQUARED_ERROR / RNX
       ! Store this mean squared error in the history if desired.
       IF (PRESENT(RECORD)) RECORD(I) = MEAN_SQUARED_ERROR
       ! Convert the average gradients.
       V1_GRAD = V1_GRAD / BATCHES
       S1_GRAD = S1_GRAD / BATCHES
       V2_GRAD = V2_GRAD / BATCHES
       S2_GRAD = S2_GRAD / BATCHES
       V3_GRAD = V3_GRAD / BATCHES
       ! Update the steps for all different parameters.
       CALL COMPUTE_STEP(SIZE(V1_GRAD(:,:)), V1_GRAD(:,:), V1_MEAN(:,:), V1_VAR(:,:))
       CALL COMPUTE_STEP(SIZE(S1_GRAD(:)), S1_GRAD(:), S1_MEAN(:), S1_VAR(:))
       DO L = 1, MNS-1
          CALL COMPUTE_STEP(SIZE(V2_GRAD(:,:,L)), V2_GRAD(:,:,L), V2_MEAN(:,:,L), V2_VAR(:,:,L))
          CALL COMPUTE_STEP(SIZE(S2_GRAD(:,L)), S2_GRAD(:,L), S2_MEAN(:,L), S2_VAR(:,L))
       END DO
       CALL COMPUTE_STEP(SIZE(V3_GRAD(:,:)), V3_GRAD(:,:), V3_MEAN(:,:), V3_VAR(:,:))
       ! Take the gradient steps (based on the computed "step" above).
       INPUT_VECS(:,:)      = INPUT_VECS(:,:)      - V1_GRAD(:,:)
       INPUT_SHIFT(:)       = INPUT_SHIFT(:)       - S1_GRAD(:)
       INTERNAL_VECS(:,:,:) = INTERNAL_VECS(:,:,:) - V2_GRAD(:,:,:)
       INTERNAL_SHIFT(:,:)  = INTERNAL_SHIFT(:,:)  - S2_GRAD(:,:)
       OUTPUT_VECS(:,:)     = OUTPUT_VECS(:,:)     - V3_GRAD(:,:)
       ! Maintain a constant max-norm across the magnitue of internal vectors.
       DO L = 1, MNS-1
          SCALAR = SQRT(MAXVAL(SUM(INTERNAL_VECS(:,:,L)**2, 1)))
          INTERNAL_VECS(:,:,L) = INTERNAL_VECS(:,:,L) / SCALAR
       END DO
     END DO fit_loop
   CONTAINS

    ! Compute the current step factor as the running average step,
    ! inversely scaled by the average magnitude of the step.
    SUBROUTINE COMPUTE_STEP(N, CURR_STEP, STEP_MEAN, STEP_VAR)
      INTEGER, INTENT(IN) :: N
      REAL(KIND=RT), DIMENSION(1:N), INTENT(INOUT) :: CURR_STEP
      REAL(KIND=RT), DIMENSION(1:N), INTENT(INOUT) :: STEP_MEAN
      REAL(KIND=RT), DIMENSION(1:N), INTENT(INOUT) :: STEP_VAR
      ! Update sliding window mean and variance calculations.
      STEP_MEAN = 0.1_RT * CURR_STEP + 0.9_RT * STEP_MEAN
      STEP_VAR = 0.001_RT * (STEP_MEAN - CURR_STEP)**2 + 0.999_RT * STEP_VAR
      STEP_VAR = MAX(STEP_VAR, EPSILON(STEP_FACTOR))
      ! Rescale the step by the means and variances.
      CURR_STEP = STEP_FACTOR * STEP_MEAN / SQRT(STEP_VAR)
    END SUBROUTINE COMPUTE_STEP

  END SUBROUTINE MINIMIZE_MSE


END MODULE PLRM
