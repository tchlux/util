! TODO:
! - Use ALLOCATE instead of automatic, because of large N seg faults.
! - Save model and load model.
! - Train only the output, internal, or input layers.
! - Embedding support, for integer valued inputs.
! - Implement similar code in C, compare speeds.


! A piecewise linear regression model.
MODULE PLRM
  USE ISO_FORTRAN_ENV, ONLY: REAL32
  IMPLICIT NONE
  ! Value used for rectification, default is 0.
  REAL(KIND=REAL32), PARAMETER :: DISCONTINUITY = 0.0_REAL32
  REAL(KIND=REAL32), PARAMETER :: SMALL_SLOPE = 0.01_REAL32
  ! Value for bias base (in terms of standard deviation for unit
  !  variance data), default value is 2.
  REAL(KIND=REAL32), PARAMETER :: BIAS_STDEVS = 2.0_REAL32
  ! MDI = Model dimension -- input
  ! MDS = Model dimension -- states (internal)
  ! MNS = Model number -- states (internal)
  ! MDO = Model dimension -- output
  INTEGER :: MDI, MDS, MNS, MDO
  ! Model parameters, for evaluation.
  REAL(KIND=REAL32), DIMENSION(:,:),   ALLOCATABLE :: INPUT_PARAMS
  REAL(KIND=REAL32), DIMENSION(:),     ALLOCATABLE :: INPUT_BIAS
  REAL(KIND=REAL32), DIMENSION(:,:,:), ALLOCATABLE :: INTERNAL_PARAMS
  REAL(KIND=REAL32), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_BIAS
  REAL(KIND=REAL32), DIMENSION(:,:),   ALLOCATABLE :: OUTPUT_PARAMS
  REAL(KIND=REAL32), DIMENSION(:),     ALLOCATABLE :: OUTPUT_BIAS
  
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
    IF (ALLOCATED(INPUT_PARAMS))    DEALLOCATE(INPUT_PARAMS)
    IF (ALLOCATED(INPUT_BIAS))      DEALLOCATE(INPUT_BIAS)
    IF (ALLOCATED(INTERNAL_PARAMS)) DEALLOCATE(INTERNAL_PARAMS)
    IF (ALLOCATED(INTERNAL_BIAS))   DEALLOCATE(INTERNAL_BIAS)
    IF (ALLOCATED(OUTPUT_PARAMS))   DEALLOCATE(OUTPUT_PARAMS)
    IF (ALLOCATED(OUTPUT_BIAS))     DEALLOCATE(OUTPUT_BIAS)
    ! Allocate space for all model parameters.
    ALLOCATE( &
         INPUT_PARAMS(DI,DS), &
         INPUT_BIAS(DS), &
         INTERNAL_PARAMS(DS, DS, NS-1), &
         INTERNAL_BIAS(DS, NS-1), &
         OUTPUT_PARAMS(DS, DO), &
         OUTPUT_BIAS(DO) &
         )
  END SUBROUTINE NEW_MODEL


  ! Initialize the weights for a model, optionally provide a random seed.
  SUBROUTINE INIT_MODEL(INPUTS, OUTPUTS, SEED)
    REAL(KIND=REAL32), DIMENSION(:,:), INTENT(IN), OPTIONAL :: INPUTS, OUTPUTS
    INTEGER, INTENT(IN), OPTIONAL :: SEED
    ! Storage for data-related statistics for high quality initialization.
    REAL(KIND=REAL32), DIMENSION(:),   ALLOCATABLE :: VALUES
    REAL(KIND=REAL32), DIMENSION(:,:), ALLOCATABLE :: TEMP_VALUES, INTERNAL_VALUES, OUTPUT_VALUES
    !  Storage for seeding the random number generator (for repeatability).
    INTEGER, DIMENSION(:), ALLOCATABLE :: SEED_ARRAY
    ! Dimension and layer iterates, number of points, scalar used to
    !  regularly space bias values, and selected percentile.
    !  (all only relevant if INPUTS is provided).
    INTEGER :: D, L, I
    REAL(KIND=REAL32) :: N, P, MEAN, STDEV
    ! Set a random seed, if one was provided (otherwise leave default).
    IF (PRESENT(SEED)) THEN
       CALL RANDOM_SEED(SIZE=I)
       ALLOCATE(SEED_ARRAY(I))
       SEED_ARRAY(:) = SEED
       CALL RANDOM_SEED(PUT=SEED_ARRAY(:))
    END IF
    ! Generate random unit-length vectors (no scaling biases) for
    !  all initial parameters in the input, internal, and output.
    CALL RANDOM_UNIT_VECTORS(INPUT_PARAMS(:,:))
    DO I = 1, MNS-1
       CALL RANDOM_UNIT_VECTORS(INTERNAL_PARAMS(:,:,I))
    END DO
    CALL RANDOM_UNIT_VECTORS(OUTPUT_PARAMS(:,:))
    ! Generate random biases for inputs and internal layers, zero
    !  bias for the output layer (first two will be rescaled).
    CALL RANDOM_NUMBER(INPUT_BIAS(:))
    CALL RANDOM_NUMBER(INTERNAL_BIAS(:,:))
    OUTPUT_BIAS(:) = 0.0_REAL32
    ! Adjust parameters based on sample input data.
    IF (PRESENT(INPUTS)) THEN
       ! Allocate storage for intermediate arrays.
       ALLOCATE( TEMP_VALUES(1:SIZE(INPUTS,2),1:MDS), &
            INTERNAL_VALUES(1:SIZE(INPUTS,2),1:MDS), &
            OUTPUT_VALUES(1:SIZE(INPUTS,2),1:MDO), &
            VALUES(1:SIZE(INPUTS,2)) )
       ! Get the number of points.
       N = REAL(SIZE(INPUTS,2), REAL32)
       ! Rescale all input parameters by the componentwise variance.
       !  This will redistribute the random vectors to be uniform
       !  about the unit variance ball (instead of its "squashed" form).
       DO D = 1, MDI
          MEAN = SUM(INPUTS(D,:) / N)
          STDEV = SQRT(MAX(EPSILON(STDEV), SUM((INPUTS(D,:) - MEAN)**2 / N)))
          INPUT_PARAMS(D,:) = INPUT_PARAMS(D,:) * STDEV
       END DO
       ! Renormalize the length of the input parameter vectors to be unit length.
       !    (only do this when the input dimension is greater than 1).
       IF (MDI .GT. 1) THEN
          DO D = 1, MDS
             INPUT_PARAMS(:,D) = INPUT_PARAMS(:,D) / NORM2(INPUT_PARAMS(:,D))
          END DO
       END IF
       ! Compute first layer of internal values.
       CALL SGEMM('T', 'N', SIZE(INPUTS,2), MDS, MDI, 1.0_REAL32, &
            INPUTS(:,:), SIZE(INPUTS,1), &
            INPUT_PARAMS(:,:), SIZE(INPUT_PARAMS,1), &
            0.0_REAL32, INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1))
       DO L = 1, MNS
          ! Adjust magnitude of vectors to make dimension values unit variance.
          DO D = 1, MDS
             MEAN = SUM(INTERNAL_VALUES(:,D) / N)
             STDEV = SQRT(MAX(EPSILON(STDEV), SUM((INTERNAL_VALUES(:,D) - MEAN)**2 / N)))
             IF (L .EQ. 1) THEN
                ! Rescale incoming vector to make output unit variance.
                INPUT_PARAMS(:,D) = INPUT_PARAMS(:,D) / STDEV
                ! Store the random variable used for bias creation.
                P = INPUT_BIAS(D)
             ELSE
                ! Rescale incoming vector to make output unit variance.
                I = L-1
                INTERNAL_PARAMS(:,D,I) = INTERNAL_PARAMS(:,D,I) / STDEV
                ! Store the random variabble used for bias creation.
                P = INTERNAL_BIAS(D,I)
             END IF
             ! Adjust current values based on updated scaling of input vector.
             INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) / STDEV
             MEAN = MEAN / STDEV
             ! Scale bias values into the range [-2, 2] about the mean.
             P = 2.0_REAL32 * (2.0_REAL32 * P - 1.0_REAL32)  -  MEAN  +  DISCONTINUITY
             !                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ makes random value in [-1,1].
             !   ^^^^^^^^^^ scales into [-2,2] range.
             !             centers discontinuity at the mean ^^^^^^^^^^^^^^^^^^^^^^^^^^
             ! 
             ! Store the bias.
             IF (L .EQ. 1) THEN
                INPUT_BIAS(D) = P
             ELSE
                INTERNAL_BIAS(D,I) = P
             END IF
             ! Compute the values for this dimension of this layer based on the bias.
             INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) + P
             WHERE (INTERNAL_VALUES(:,D) .LT. DISCONTINUITY)
                INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) * SMALL_SLOPE
             END WHERE
          END DO
          ! Compute the internal values for the next layer.
          IF (L .LT. MNS) THEN
             CALL SGEMM('N', 'N', SIZE(INPUTS,2), MDS, MDS, 1.0_REAL32, &
                  INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
                  INTERNAL_PARAMS(:,:,L), SIZE(INTERNAL_PARAMS,1), &
                  0.0_REAL32, TEMP_VALUES(:,:), SIZE(TEMP_VALUES,1))
             INTERNAL_VALUES(:,:) = TEMP_VALUES(:,:)
          END IF
       END DO
       ! Compute the output.
       CALL SGEMM('N', 'N', SIZE(INPUTS,2), MDO, MDS, 1.0_REAL32, &
            INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
            OUTPUT_PARAMS(:,:), SIZE(OUTPUT_PARAMS,1), &
            0.0_REAL32, OUTPUT_VALUES(:,:), SIZE(OUTPUT_VALUES,1))
       ! Adjust magnitude of vectors to make dimension values unit variance.
       DO D = 1, MDO
          ! Make this component of the output have unit variance.
          MEAN = SUM(OUTPUT_VALUES(:,D) / N)
          STDEV = SQRT(MAX(EPSILON(STDEV), SUM((OUTPUT_VALUES(:,D) - MEAN)**2 / N)))
          OUTPUT_PARAMS(:,D) = OUTPUTS(:,D) / STDEV
          OUTPUT_VALUES(:,D) = OUTPUT_VALUES(:,D) / STDEV
          ! If outputs are provided.
          IF (PRESENT(OUTPUTS)) THEN
             ! Compute the variance of the output, make results match that.
             MEAN = SUM(OUTPUTS(:,D) / N)
             STDEV = SQRT(MAX(EPSILON(STDEV), SUM((OUTPUTS(:,D) - MEAN)**2 / N)))
             OUTPUT_PARAMS(:,D) = OUTPUT_PARAMS(:,D) * STDEV
             OUTPUT_VALUES(:,D) = OUTPUT_VALUES(:,D) * STDEV
             ! Set the bias value to make the output mean match.
             OUTPUT_BIAS(D) = MEAN - SUM(OUTPUT_VALUES(:,D) / N, DIM=1)
          ELSE
             ! Make this component of output have zero mean.
             OUTPUT_BIAS(D) = SUM(OUTPUT_VALUES(:,D) / N)
          END IF
       END DO
       ! Deallocate memory for temporary arrays.
       DEALLOCATE( INTERNAL_VALUES, OUTPUT_VALUES, VALUES )
    ELSE
       ! Values are in range [-2,2], assuming no scaling changes,
       !  this will capture 2 standard deviations in either direction.
       INPUT_BIAS(:) = 2.0_REAL32 * (INPUT_BIAS(:) * 2.0_REAL32 - 1.0_REAL32)
       INTERNAL_BIAS(:,:) = 2.0_REAL32 * (INTERNAL_BIAS(:,:) * 2.0_REAL32 - 1.0_REAL32)
    END IF

  CONTAINS

    ! Generate randomly distributed vectors on the N-sphere.
    SUBROUTINE RANDOM_UNIT_VECTORS(COLUMN_VECTORS)
      REAL(KIND=REAL32), INTENT(OUT), DIMENSION(:,:) :: COLUMN_VECTORS
      INTEGER :: I 
      ! Generate random numbers in the range [0,1].
      CALL RANDOM_NUMBER(COLUMN_VECTORS(:,:))
      ! Map those random numbers to the range [-1,1].
      COLUMN_VECTORS(:,:) = 2.0_REAL32 * COLUMN_VECTORS(:,:) - 1.0_REAL32
      ! Make the vectors uniformly distributed on the unit ball (for dimension > 1).
      IF (SIZE(COLUMN_VECTORS,1) .GT. 1) THEN
         ! Compute the inverse hyperbolic tangent of all values.
         COLUMN_VECTORS(:,:) = ATANH(COLUMN_VECTORS(:,:))
         ! Normalize all vectors to have unit length.
         DO I = 1, SIZE(COLUMN_VECTORS,2)
            COLUMN_VECTORS(:,I) = COLUMN_VECTORS(:,I) / NORM2(COLUMN_VECTORS(:,I))
         END DO
      END IF
    END SUBROUTINE RANDOM_UNIT_VECTORS

  END SUBROUTINE INIT_MODEL


  ! Evaluate the piecewise linear regression model, do not store any
  !  of the intermediate states (used 
  SUBROUTINE EVALUATE_ONE(INPUT, OUTPUT)
    REAL(KIND=REAL32), INTENT(IN), DIMENSION(:) :: INPUT
    ! Record of values after the last transformation.
    REAL(KIND=REAL32), INTENT(OUT), DIMENSION(:) :: OUTPUT
    ! Internal values.
    REAL(KIND=REAL32), DIMENSION(MDS) :: INTERNAL_VALUES
    ! Compute the number of nodes (to know where "bias" value is).
    INTEGER :: I, IP1
    ! Make sure that this model has been initialized.
    IF (.NOT. ALLOCATED(INPUT_PARAMS)) THEN
       OUTPUT(:) = 0.0
       RETURN
    END IF
    ! Compute the input layer.
    INTERNAL_VALUES(:) = INPUT_BIAS(:) + MATMUL(INPUT(:), INPUT_PARAMS(:,:))
    WHERE (INTERNAL_VALUES(:) .LT. DISCONTINUITY)
       INTERNAL_VALUES(:) = INTERNAL_VALUES(:) * SMALL_SLOPE
    END WHERE
    ! Compute the next set of internal values with a rectified activation.
    DO I = 1, MNS-1
       IP1 = I+1
       INTERNAL_VALUES(:) = INTERNAL_BIAS(:,I) + &
            MATMUL(INTERNAL_VALUES(:), INTERNAL_PARAMS(:,:,I))
       WHERE (INTERNAL_VALUES(:) .LT. DISCONTINUITY)
          INTERNAL_VALUES(:) = INTERNAL_VALUES(:) * SMALL_SLOPE
       END WHERE
    END DO
    ! Compute the output.
    OUTPUT(:) = OUTPUT_BIAS(:) + &
         MATMUL(INTERNAL_VALUES(:), OUTPUT_PARAMS(:,:))
  END SUBROUTINE EVALUATE_ONE

  ! Evaluate the piecewise linear regression model, do not store any
  !  of the intermediate states (used 
  SUBROUTINE EVALUATE(INPUTS, OUTPUTS)
    REAL(KIND=REAL32), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=REAL32), INTENT(OUT), DIMENSION(:,:) :: OUTPUTS
    ! Internal values.
    REAL(KIND=REAL32), DIMENSION(:,:), ALLOCATABLE :: TEMP_VALUES, INTERNAL_VALUES
    ! Compute the number of nodes (to know where "bias" value is).
    INTEGER :: D, L, N
    ! Make sure that this model has been initialized.
    IF (.NOT. ALLOCATED(INPUT_PARAMS)) THEN
       OUTPUTS(:,:) = 0.0
       RETURN
    END IF
    ! Get the number of points.
    N = SIZE(INPUTS,2)
    ! Allocate storage for the internal values.
    ALLOCATE(INTERNAL_VALUES(1:N,1:MDS), TEMP_VALUES(1:N,1:MDS))
    ! Compute the input transformation.
    CALL SGEMM('T', 'N', N, MDS, MDI, 1.0_REAL32, &
         INPUTS(:,:), SIZE(INPUTS,1), &
         INPUT_PARAMS(:,:), SIZE(INPUT_PARAMS,1), &
         0.0_REAL32, INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1))
    DO D = 1, MDS
       INTERNAL_VALUES(:,D) = INPUT_BIAS(D) + INTERNAL_VALUES(:,D)
    END DO
    WHERE (INTERNAL_VALUES(:,:) .LT. DISCONTINUITY)
       INTERNAL_VALUES(:,:) = INTERNAL_VALUES(:,:) * SMALL_SLOPE
    END WHERE
    ! Compute the next set of internal values with a rectified activation.
    DO L = 1, MNS-1
       CALL SGEMM('N', 'N', N, MDS, MDS, 1.0_REAL32, &
            INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
            INTERNAL_PARAMS(:,:,L), SIZE(INTERNAL_PARAMS,1), &
            0.0_REAL32, TEMP_VALUES(:,:), SIZE(TEMP_VALUES,1))
       INTERNAL_VALUES(:,:) = TEMP_VALUES(:,:)
       DO D = 1, MDS
          INTERNAL_VALUES(:,D) =  INTERNAL_BIAS(D,L) + INTERNAL_VALUES(:,D)
       END DO
       WHERE (INTERNAL_VALUES(:,:) .LT. DISCONTINUITY)
          INTERNAL_VALUES(:,:) = INTERNAL_VALUES(:,:) * SMALL_SLOPE
       END WHERE
    END DO
    ! Return the final embedded layer if the outputs have the size of the embedding.
    IF ((MDS .NE. MDO) .AND. (SIZE(OUTPUTS,1) .EQ. MDS)) THEN
       ! Do the necessary transpose operation one dimension at a time.
       DO L = 1, MDS
          OUTPUTS(L,1:N) = INTERNAL_VALUES(1:N,L)
       END DO
    ! Otherwise, assume regular output computation.
    ELSE
       CALL SGEMM('T', 'T', MDO, N, MDS, 1.0_REAL32, &
            OUTPUT_PARAMS(:,:), SIZE(OUTPUT_PARAMS,1), &
            INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
            0.0_REAL32, OUTPUTS(:,:), SIZE(OUTPUTS,1))
       DO D = 1, MDO
          OUTPUTS(:,D) = OUTPUT_BIAS(D) + OUTPUTS(:,D)
       END DO
    END IF
    ! Deallocate temporary variables.
    DEALLOCATE(INTERNAL_VALUES, TEMP_VALUES)
  END SUBROUTINE EVALUATE


  ! Compute the gradient of the sum of squared error of this regression
  ! model with respect to its parameters given input and output pairs.
  SUBROUTINE SSE_GRADIENT(INPUTS, OUTPUTS, SUM_SQUARED_ERROR, &
       INPUT_PARAMS_GRADIENT, INPUT_BIAS_GRADIENT, &
       INTERNAL_PARAMS_GRADIENT, INTERNAL_BIAS_GRADIENT, &
       OUTPUT_PARAMS_GRADIENT, OUTPUT_BIAS_GRADIENT)
    ! Data values stored with contiguous points: shape = (MDI,N)
    REAL(KIND=REAL32), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=REAL32), INTENT(IN), DIMENSION(:,:) :: OUTPUTS
    ! Sum (over all data) squared error (summed over dimensions).
    REAL(KIND=REAL32), INTENT(INOUT) :: SUM_SQUARED_ERROR
    ! Model parameter gradients.
    REAL(KIND=REAL32), INTENT(INOUT), DIMENSION(:,:)   :: INPUT_PARAMS_GRADIENT
    REAL(KIND=REAL32), INTENT(INOUT), DIMENSION(:)     :: INPUT_BIAS_GRADIENT
    REAL(KIND=REAL32), INTENT(INOUT), DIMENSION(:,:,:) :: INTERNAL_PARAMS_GRADIENT
    REAL(KIND=REAL32), INTENT(INOUT), DIMENSION(:,:)   :: INTERNAL_BIAS_GRADIENT
    REAL(KIND=REAL32), INTENT(INOUT), DIMENSION(:,:)   :: OUTPUT_PARAMS_GRADIENT
    REAL(KIND=REAL32), INTENT(INOUT), DIMENSION(:)     :: OUTPUT_BIAS_GRADIENT     
    ! Internal values and lists of which ones are nonzero.
    REAL(KIND=REAL32), DIMENSION(SIZE(INPUTS,2),MDO)     :: OUTPUT_GRADIENT
    REAL(KIND=REAL32), DIMENSION(SIZE(INPUTS,2),MDS,MNS) :: INTERNAL_VALUES
    LOGICAL, DIMENSION(SIZE(INPUTS,2),MDS) :: LT_DISCONTINUITY
    ! Number of points.
    REAL(KIND=REAL32) :: N
    INTEGER :: D
    ! Evaluate the model at all data points, store outputs in "OUTPUT_GRADIENT"
    CALL EVALUATE_BATCH()
    ! Compute the new gradient.
    DO D = 1, MDO
       OUTPUT_GRADIENT(:,D) = OUTPUT_GRADIENT(:,D) - OUTPUTS(D,:)
    END DO
    ! Get the number of points.
    N = REAL(SIZE(INPUTS,2), REAL32)
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
      CALL SGEMM('T', 'N', SIZE(INPUTS,2), MDS, MDI, 1.0_REAL32, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INPUT_PARAMS(:,:), SIZE(INPUT_PARAMS,1), &
           0.0_REAL32, INTERNAL_VALUES(:,:,1), SIZE(INTERNAL_VALUES,1))
      DO D = 1, MDS
         INTERNAL_VALUES(:,D,1) = INPUT_BIAS(D) + INTERNAL_VALUES(:,D,1)
      END DO
      WHERE (INTERNAL_VALUES(:,:,1) .LT. DISCONTINUITY)
         INTERNAL_VALUES(:,:,1) = INTERNAL_VALUES(:,:,1) * SMALL_SLOPE
      END WHERE
      ! Compute the next set of internal values with a rectified activation.
      DO L = 1, MNS-1
         LP1 = L+1
         CALL SGEMM('N', 'N', SIZE(INPUTS,2), MDS, MDS, 1.0_REAL32, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_PARAMS(:,:,L), SIZE(INTERNAL_PARAMS,1), &
              0.0_REAL32, INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1))
         DO D = 1, MDS
            INTERNAL_VALUES(:,D,LP1) =  INTERNAL_BIAS(D,L) + INTERNAL_VALUES(:,D,LP1)
         END DO
         WHERE (INTERNAL_VALUES(:,:,LP1) .LT. DISCONTINUITY)
            INTERNAL_VALUES(:,:,LP1) = INTERNAL_VALUES(:,:,LP1) * SMALL_SLOPE
         END WHERE
      END DO
      ! Compute the output.
      CALL SGEMM('N', 'N', SIZE(INPUTS,2), MDO, MDS, 1.0_REAL32, &
           INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_PARAMS(:,:), SIZE(OUTPUT_PARAMS,1), &
           0.0_REAL32, OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1))
      DO D = 1, MDO
         OUTPUT_GRADIENT(:,D) = OUTPUT_BIAS(D) + OUTPUT_GRADIENT(:,D)
      END DO
    END SUBROUTINE EVALUATE_BATCH

    SUBROUTINE GRADIENT_BATCH()
      ! L   - layer index
      ! LP1 - layer index "plus 1" -> "P1"
      INTEGER :: L, LP1
      ! Compute the average output gradient.
      OUTPUT_BIAS_GRADIENT(:) = SUM(OUTPUT_GRADIENT(:,:), 1) / N &
           + OUTPUT_BIAS_GRADIENT(:)
      ! Compute output parameter gradients.
      CALL SGEMM('T', 'N', MDS, MDO, SIZE(INPUTS,2), 1.0_REAL32 / N, &
           INTERNAL_VALUES(:,:,MNS), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           1.0_REAL32, OUTPUT_PARAMS_GRADIENT(:,:), SIZE(OUTPUT_PARAMS_GRADIENT,1))
      ! Compute the gradient for the last internal vector space values.
      LT_DISCONTINUITY(:,:) = INTERNAL_VALUES(:,:,MNS) .LT. DISCONTINUITY
      CALL SGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDO, 1.0_REAL32, &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           OUTPUT_PARAMS(:,:), SIZE(OUTPUT_PARAMS,1), &
           0.0_REAL32, INTERNAL_VALUES(:,:,MNS), SIZE(INTERNAL_VALUES,1))
      WHERE (LT_DISCONTINUITY(:,:))
         INTERNAL_VALUES(:,:,MNS) = INTERNAL_VALUES(:,:,MNS) * SMALL_SLOPE
      END WHERE
      ! Cycle over all internal layers.
      INTERNAL_REPRESENTATIONS : DO L = MNS-1, 1, -1
         LP1 = L+1
         ! Compute the bias gradient.
         INTERNAL_BIAS_GRADIENT(:,L) = SUM(INTERNAL_VALUES(:,:,LP1), 1) / N &
              + INTERNAL_BIAS_GRADIENT(:,L)
         ! Compute the gradient with respect to each output and all inputs.
         CALL SGEMM('T', 'N', MDS, MDS, SIZE(INPUTS,2), 1.0_REAL32 / N, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              1.0_REAL32, INTERNAL_PARAMS_GRADIENT(:,:,L), SIZE(INTERNAL_PARAMS_GRADIENT,1))
         ! Propogate the gradient to the immediately preceding layer.
         LT_DISCONTINUITY(:,:) = INTERNAL_VALUES(:,:,L) .LT. DISCONTINUITY 
         CALL SGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDS, 1.0_REAL32, &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_PARAMS(:,:,L), SIZE(INTERNAL_PARAMS,1), &
              0.0_REAL32, INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1))
         WHERE (LT_DISCONTINUITY(:,:))
            INTERNAL_VALUES(:,:,L) = INTERNAL_VALUES(:,:,L) * SMALL_SLOPE
         END WHERE
      END DO INTERNAL_REPRESENTATIONS
      ! Compute the input parameter bias values.
      INPUT_BIAS_GRADIENT(:) = SUM(INTERNAL_VALUES(:,:,1), 1) / N &
           + INPUT_BIAS_GRADIENT(:)
      ! Compute the gradient of all input parameters.
      !   [the INPUTS are transposed already, shape = (MDI,N)]
      CALL SGEMM('N', 'N', MDI, MDS, SIZE(INPUTS,2), 1.0_REAL32 / N, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INTERNAL_VALUES(:,:,1), SIZE(INTERNAL_VALUES,1), &
           1.0_REAL32, INPUT_PARAMS_GRADIENT(:,:), SIZE(INPUT_PARAMS_GRADIENT,1))
    END SUBROUTINE GRADIENT_BATCH

  END SUBROUTINE SSE_GRADIENT


  ! Fit input / output pairs by minimizing mean squared error.
  SUBROUTINE MINIMIZE_MSE(X, Y, STEPS, BATCH_SIZE, NUM_THREADS, &
       STEP_SIZE, KEEP_BEST, MEAN_SQUARED_ERROR)
    REAL(KIND=REAL32), INTENT(IN), DIMENSION(:,:) :: X
    REAL(KIND=REAL32), INTENT(IN), DIMENSION(:,:) :: Y
    INTEGER,           INTENT(IN), OPTIONAL :: STEPS, BATCH_SIZE, NUM_THREADS
    REAL(KIND=REAL32), INTENT(IN), OPTIONAL :: STEP_SIZE
    LOGICAL,           INTENT(IN), OPTIONAL :: KEEP_BEST
    REAL(KIND=REAL32), INTENT(OUT) :: MEAN_SQUARED_ERROR
    !  gradient step arrays, 4 copies of internal model (total of 5).
    REAL(KIND=REAL32), DIMENSION(MDI,MDS)       :: P1_GRAD, P1_MEAN, P1_VAR, BEST_P1
    REAL(KIND=REAL32), DIMENSION(MDS)           :: B1_GRAD, B1_MEAN, B1_VAR, BEST_B1
    REAL(KIND=REAL32), DIMENSION(MDS,MDS,MNS-1) :: P2_GRAD, P2_MEAN, P2_VAR, BEST_P2
    REAL(KIND=REAL32), DIMENSION(MDS,MNS-1)     :: B2_GRAD, B2_MEAN, B2_VAR, BEST_B2
    REAL(KIND=REAL32), DIMENSION(MDS,MDO)       :: P3_GRAD, P3_MEAN, P3_VAR, BEST_P3
    REAL(KIND=REAL32), DIMENSION(MDO)           :: B3_GRAD, B3_MEAN, B3_VAR, BEST_B3
    !  local integers
    LOGICAL :: REVERT_TO_BEST
    INTEGER :: BS, I, L, NX, NT, S, BATCH_START, BATCH_END
    REAL(KIND=REAL32) :: RNX, BATCHES, STEP_FACTOR, PREV_MSE, BEST_MSE, SCALAR
    ! Function that is defined by OpenMP.
    INTERFACE
       FUNCTION OMP_GET_MAX_THREADS()
         INTEGER :: OMP_GET_MAX_THREADS
       END FUNCTION OMP_GET_MAX_THREADS
    END INTERFACE
    ! Number of points.
    NX = SIZE(X,2)
    RNX = REAL(NX, REAL32)
    ! Set the number of steps.
    IF (PRESENT(STEPS)) THEN
       S = STEPS
    ELSE
       S = 100
    END IF
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
    ! Set the "keep best" boolean.
    IF (PRESENT(KEEP_BEST)) THEN
       REVERT_TO_BEST = KEEP_BEST
    ELSE
       REVERT_TO_BEST = (S .GE. 10)
    END IF
    ! Initial mean squared error is "max float value".
    PREV_MSE = HUGE(PREV_MSE)
    BEST_MSE = HUGE(BEST_MSE)
    ! Set the average step sizes.
    P1_MEAN(:,:)   = 0.0_REAL32
    B1_MEAN(:)     = 0.0_REAL32
    P2_MEAN(:,:,:) = 0.0_REAL32
    B2_MEAN(:,:)   = 0.0_REAL32
    P3_MEAN(:,:)   = 0.0_REAL32
    B3_MEAN(:)     = 0.0_REAL32
    ! Set the average step size variances.
    P1_VAR(:,:)   = 0.0_REAL32
    B1_VAR(:)     = 0.0_REAL32
    P2_VAR(:,:,:) = 0.0_REAL32
    B2_VAR(:,:)   = 0.0_REAL32
    P3_VAR(:,:)   = 0.0_REAL32
    B3_VAR(:)     = 0.0_REAL32
    ! Iterate, taking steps with the average gradient over all data.
    fit_loop : DO I = 1, S
       ! Compute the average gradient over all points.
       MEAN_SQUARED_ERROR = 0.0_REAL32
       ! Set gradients to zero initially.
       P1_GRAD(:,:)   = 0.0_REAL32
       B1_GRAD(:)     = 0.0_REAL32
       P2_GRAD(:,:,:) = 0.0_REAL32
       B2_GRAD(:,:)   = 0.0_REAL32
       P3_GRAD(:,:)   = 0.0_REAL32
       B3_GRAD(:)     = 0.0_REAL32
       ! Count the number of batches.
       BATCHES = 0.0_REAL32
       !$OMP PARALLEL DO NUM_THREADS(NT) PRIVATE(BATCH_END) &
       !$OMP&  REDUCTION(+: BATCHES, MEAN_SQUARED_ERROR, &
       !$OMP&  P1_GRAD, B1_GRAD, P2_GRAD, B2_GRAD, P3_GRAD, B3_GRAD)
       DO BATCH_START = 1, NX, BS
          BATCHES = BATCHES + 1.0_REAL32
          BATCH_END = MIN(NX, BATCH_START+BS-1)
          ! Compute the gradient from all data.
          CALL SSE_GRADIENT(X(:,BATCH_START:BATCH_END), &
               Y(:,BATCH_START:BATCH_END), MEAN_SQUARED_ERROR, &
               P1_GRAD, B1_GRAD, P2_GRAD, B2_GRAD, P3_GRAD, B3_GRAD)
       END DO
       !$OMP END PARALLEL DO
       ! Convert the sum of squared errors into the mean squared error.
       MEAN_SQUARED_ERROR = MEAN_SQUARED_ERROR / RNX
       ! Update the saved "best" model based on error (only when dropout is disabled).
       IF (MEAN_SQUARED_ERROR .LT. BEST_MSE) THEN
          BEST_MSE = MEAN_SQUARED_ERROR
          BEST_P1 = INPUT_PARAMS
          BEST_B1 = INPUT_BIAS
          BEST_P2 = INTERNAL_PARAMS
          BEST_B2 = INTERNAL_BIAS
          BEST_P3 = OUTPUT_PARAMS
          BEST_B3 = OUTPUT_BIAS
       END IF
       ! Update the step factor based on model improvement.
       IF (MEAN_SQUARED_ERROR .LT. PREV_MSE) THEN
          STEP_FACTOR = STEP_FACTOR * 1.01_REAL32
       ELSE
          STEP_FACTOR = STEP_FACTOR / 1.01_REAL32
       END IF
       PREV_MSE = MEAN_SQUARED_ERROR
       ! Convert the average gradients.
       P1_GRAD = P1_GRAD / BATCHES
       B1_GRAD = B1_GRAD / BATCHES
       P2_GRAD = P2_GRAD / BATCHES
       B2_GRAD = B2_GRAD / BATCHES
       P3_GRAD = P3_GRAD / BATCHES
       B3_GRAD = B3_GRAD / BATCHES
       ! Update the steps for all different parameters.
       CALL COMPUTE_STEP(SIZE(P1_GRAD(:,:)), P1_GRAD(:,:), P1_MEAN(:,:), P1_VAR(:,:))
       CALL COMPUTE_STEP(SIZE(B1_GRAD(:)), B1_GRAD(:), B1_MEAN(:), B1_VAR(:))
       DO L = 1, MNS-1
          CALL COMPUTE_STEP(SIZE(P2_GRAD(:,:,L)), P2_GRAD(:,:,L), P2_MEAN(:,:,L), P2_VAR(:,:,L))
          CALL COMPUTE_STEP(SIZE(B2_GRAD(:,L)), B2_GRAD(:,L), B2_MEAN(:,L), B2_VAR(:,L))
       END DO
       DO L = 1, MDO
          CALL COMPUTE_STEP(SIZE(P3_GRAD(:,L)), P3_GRAD(:,L), P3_MEAN(:,L), P3_VAR(:,L))
       END DO
       CALL COMPUTE_STEP(SIZE(B3_GRAD), B3_GRAD, B3_MEAN, B3_VAR)
       ! Take the gradient steps (based on the computed "step" above).
       INPUT_PARAMS(:,:)      = INPUT_PARAMS(:,:)      - P1_GRAD(:,:)
       INPUT_BIAS(:)          = INPUT_BIAS(:)          - B1_GRAD(:)
       INTERNAL_PARAMS(:,:,:) = INTERNAL_PARAMS(:,:,:) - P2_GRAD(:,:,:)
       INTERNAL_BIAS(:,:)     = INTERNAL_BIAS(:,:)     - B2_GRAD(:,:)
       OUTPUT_PARAMS(:,:)     = OUTPUT_PARAMS(:,:)     - P3_GRAD(:,:)
       OUTPUT_BIAS(:)         = OUTPUT_BIAS(:)         - B3_GRAD(:)
       ! Enforce maxnorm on internal parameters.
       DO L = 1, MNS-1
          SCALAR = MAXVAL(INTERNAL_PARAMS(:,:,L))
          IF (SCALAR .GT. REAL(MDS,REAL32)) THEN
             SCALAR = REAL(MDS,REAL32) / SCALAR
             INTERNAL_PARAMS(:,:,L) = INTERNAL_PARAMS(:,:,L) * SCALAR
          END IF
       END DO
     END DO fit_loop
     ! Restore the best model seen so far (if a reasonable number of steps were taken).
    IF (REVERT_TO_BEST) THEN
       MEAN_SQUARED_ERROR = BEST_MSE
       INPUT_PARAMS = BEST_P1
       INPUT_BIAS = BEST_B1
       INTERNAL_PARAMS = BEST_P2
       INTERNAL_BIAS = BEST_B2
       OUTPUT_PARAMS = BEST_P3
       OUTPUT_BIAS = BEST_B3
    END IF
  CONTAINS

    ! Compute the current step factor as the running average step,
    ! inversely scaled by the average magnitude of the step.
    SUBROUTINE COMPUTE_STEP(N, CURR_STEP, STEP_MEAN, STEP_VAR)
      INTEGER, INTENT(IN) :: N
      REAL(KIND=REAL32), DIMENSION(1:N), INTENT(INOUT) :: CURR_STEP
      REAL(KIND=REAL32), DIMENSION(1:N), INTENT(INOUT) :: STEP_MEAN
      REAL(KIND=REAL32), DIMENSION(1:N), INTENT(INOUT) :: STEP_VAR
      ! Update sliding window mean and variance calculations.
      STEP_MEAN = 0.1 * CURR_STEP + 0.9 * STEP_MEAN
      STEP_VAR = 0.01 * (STEP_MEAN - CURR_STEP)**2 + 0.99 * STEP_VAR
      STEP_VAR = MAX(STEP_VAR, EPSILON(STEP_FACTOR))
      ! Take averaged gradient descent step.
      CURR_STEP = STEP_FACTOR * STEP_MEAN
      ! Start scaling by step magnitude by variance once enough data is collected.
      IF (I .GT. 10) THEN
         CURR_STEP = CURR_STEP / SQRT(STEP_VAR)
      END IF
    END SUBROUTINE COMPUTE_STEP

  END SUBROUTINE MINIMIZE_MSE


END MODULE PLRM
