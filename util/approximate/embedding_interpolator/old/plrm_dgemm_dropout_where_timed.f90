! A piecewise linear regression model.
MODULE PLRM
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE
  ! Value used for rectification, default is 0.
  REAL(KIND=REAL64), PARAMETER :: DISCONTINUITY = 0.0_REAL64
  REAL(KIND=REAL64), PARAMETER :: SMALL_SLOPE = 0.01_REAL64
  ! Value for bias base (in terms of standard deviation for unit
  !  variance data), default value is 2.
  REAL(KIND=REAL64), PARAMETER :: BIAS_VALUE = 2.0_REAL64
  ! MDI = Model dimension -- input
  ! MDS = Model dimension -- states (internal)
  ! MNS = Model number -- states (internal)
  ! MDO = Model dimension -- output
  INTEGER :: MDI, MDS, MNS, MDO
  ! Model parameters, for evaluation.
  REAL(KIND=REAL64), DIMENSION(:,:),   ALLOCATABLE :: INPUT_PARAMS
  REAL(KIND=REAL64), DIMENSION(:),     ALLOCATABLE :: INPUT_BIAS
  REAL(KIND=REAL64), DIMENSION(:,:,:), ALLOCATABLE :: INTERNAL_PARAMS
  REAL(KIND=REAL64), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_BIAS
  REAL(KIND=REAL64), DIMENSION(:,:),   ALLOCATABLE :: OUTPUT_PARAMS
  REAL(KIND=REAL64), DIMENSION(:),     ALLOCATABLE :: OUTPUT_BIAS
  REAL :: DGEMM_TIMER, OTHER_TIMER, START_TIME, STOP_TIME
  
  PUBLIC :: NEW_MODEL, INIT_MODEL, MINIMIZE_MSE, EVALUATE_ONE
  PRIVATE :: DROP_INTERNAL_NODES, SWAP_PARAMS, SSE_GRADIENT


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
  SUBROUTINE INIT_MODEL(DATA, SEED)
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN), OPTIONAL :: DATA
    INTEGER, INTENT(IN), OPTIONAL :: SEED
    INTEGER :: I
    ! Storage for data-related statistics for high quality initialization.
    REAL(KIND=REAL64), DIMENSION(:),   ALLOCATABLE :: SCALE, RANGE
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: DATA_VALUES
    !  Storage for seeding the random number generator (for repeatability).
    INTEGER, DIMENSION(:), ALLOCATABLE :: SEED_ARRAY
    ! Set a random seed, if one was provided (otherwise leave default).
    IF (PRESENT(SEED)) THEN
       CALL RANDOM_SEED(SIZE=I)
       ALLOCATE(SEED_ARRAY(I))
       SEED_ARRAY(:) = SEED
       CALL RANDOM_SEED(PUT=SEED_ARRAY(:))
    END IF
    ! 
    ! TODO:
    !  - compute the scale and range of inputs, use that to 
    !    adjust the input layer
    !  - compute the scale and range of each intermediate layer and
    !    adjust all initial parameters accordingly
    !  - study output surface of model, make sure it looks like a
    !    random valued surface with * variance and zero mean.
    ! 
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
    OUTPUT_BIAS(:) = 0.0_REAL64
    ! Set input biases so that the surface produced by the rectified linear
    !  functions produce well spaced simplices (not too many near origin).
    INPUT_BIAS(:) = BIAS_VALUE * INPUT_BIAS(:) ** (1.0_REAL64 / REAL(MDS,REAL64))
    ! Internal biases will depend on vector direction, they will
    !  need to be in this range, but this will disable many.
    INTERNAL_BIAS(:,:) = BIAS_VALUE * INTERNAL_BIAS(:,:) * 2.0_REAL64 - BIAS_VALUE

  CONTAINS

    ! Generate randomly distributed vectors on the N-sphere.
    SUBROUTINE RANDOM_UNIT_VECTORS(COLUMN_VECTORS)
      REAL(KIND=REAL64), INTENT(OUT), DIMENSION(:,:) :: COLUMN_VECTORS
      INTEGER :: I 
      ! Generate random numbers in the range [0,1].
      CALL RANDOM_NUMBER(COLUMN_VECTORS(:,:))
      ! Map those random numbers to the range [-1,1].
      COLUMN_VECTORS(:,:) = 2.0_REAL64 * COLUMN_VECTORS(:,:) - 1.0_REAL64
      ! Compute the inverse hyperbolic tangent of all values.
      COLUMN_VECTORS(:,:) = ATANH(COLUMN_VECTORS(:,:))
      ! Normalize all vectors to have unit length.
      DO I = 1, SIZE(COLUMN_VECTORS,2)
         COLUMN_VECTORS(:,I) = COLUMN_VECTORS(:,I) / NORM2(COLUMN_VECTORS(:,I))
      END DO
    END SUBROUTINE RANDOM_UNIT_VECTORS

  END SUBROUTINE INIT_MODEL


  ! Swap the parameters of nodes at indices N1 and N2 in internal layer L.
  SUBROUTINE SWAP_PARAMS(L, N1, N2)
    INTEGER, INTENT(IN) :: L, N1, N2
    INTEGER :: I
    REAL(KIND=REAL64) :: TEMP_VAL
    ! Swap the input parameters.
    IF (L .EQ. 1) THEN
       DO I = 1, MDI
          TEMP_VAL = INPUT_PARAMS(I,N1)
          INPUT_PARAMS(I,N1) = INPUT_PARAMS(I,N2)
          INPUT_PARAMS(I,N2) = TEMP_VAL
       END DO
       INPUT_BIAS(N1) = INPUT_BIAS(N2)
       INPUT_BIAS(N2) = TEMP_VAL
    ELSE
       DO I = 1, MDS
          TEMP_VAL = INTERNAL_PARAMS(I,N1,L-1)
          INTERNAL_PARAMS(I,N1,L-1) = INTERNAL_PARAMS(I,N2,L-1)
          INTERNAL_PARAMS(I,N2,L-1) = TEMP_VAL
       END DO
       TEMP_VAL = INTERNAL_BIAS(N1,L-1)
       INTERNAL_BIAS(N1,L-1) = INTERNAL_BIAS(N2,L-1)
       INTERNAL_BIAS(N2,L-1) = TEMP_VAL
    END IF
    ! Swap the output parameters.
    IF (L .EQ. MNS-1) THEN
       DO I = 1, MDO
          TEMP_VAL = OUTPUT_PARAMS(N1,I)
          OUTPUT_PARAMS(N1,I) = OUTPUT_PARAMS(N2,I)
          OUTPUT_PARAMS(N2,I) = TEMP_VAL
       END DO
    ELSE
       DO I = 1, MDS
          TEMP_VAL = INTERNAL_PARAMS(N1,I,L)
          INTERNAL_PARAMS(N1,I,L) = INTERNAL_PARAMS(N2,I,L)
          INTERNAL_PARAMS(N2,I,L) = TEMP_VAL
       END DO
    END IF
  END SUBROUTINE SWAP_PARAMS


  ! Randomly disable N internal nodes at each internal layer by
  ! "swapping" their parameters to the back of the layer.
  SUBROUTINE DROP_INTERNAL_NODES(N)
    INTEGER, INTENT(IN) :: N
    INTEGER :: L, I, D, M
    REAL :: R
    ! Cycle over all internal layers.
    LAYERS : DO L = 1, MNS-1
       ! Compute the number of nodes that can be sampled from.
       M = MDS - N
       ! Cycle over all nodes that need to be dropped.
       DROPPED_NODES : DO I = MDS, MDS-N+1
          ! Pick a random node to drop.
          CALL RANDOM_NUMBER(R)
          D = 1 + INT(R * M)
          ! Swap the node that will be dropped to the end of the
          !  "droppable" nodes range.
          IF (D .NE. M) CALL SWAP_PARAMS(L, D, M)
          ! Swap the dropped node into the "dropped" range (end of layer).
          IF (M .NE. I) CALL SWAP_PARAMS(L, M, I)
          ! Exit if there are no more new nodes to swap in.
          IF (M .EQ. 1) EXIT DROPPED_NODES
          ! Reduce the number of droppable nodes 
          M = M - 1
       END DO DROPPED_NODES
    END DO LAYERS
  END SUBROUTINE DROP_INTERNAL_NODES


  ! Evaluate the piecewise linear regression model, do not store any
  !  of the intermediate states (used 
  SUBROUTINE EVALUATE_ONE(INPUT, OUTPUT)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: INPUT
    ! Record of values after the last transformation.
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(:) :: OUTPUT
    ! Internal values.
    REAL(KIND=REAL64), DIMENSION(MDS) :: INTERNAL_VALUES
    ! Compute the number of nodes (to know where "bias" value is).
    INTEGER :: I, IP1
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
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(:,:) :: OUTPUTS
    ! Internal values.
    REAL(KIND=REAL64), DIMENSION(SIZE(INPUTS,2),MDS) :: INTERNAL_VALUES
    ! Compute the number of nodes (to know where "bias" value is).
    INTEGER :: D, L, LP1
    ! Compute the input transformation.
    CALL DGEMM('T', 'N', SIZE(INPUTS,2), MDS, MDI, 1.0_REAL64, &
         INPUTS(:,:), SIZE(INPUTS,1), &
         INPUT_PARAMS(:,:), SIZE(INPUT_PARAMS,1), &
         0.0_REAL64, INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1))
    DO D = 1, MDS
       INTERNAL_VALUES(:,D) = INPUT_BIAS(D) + INTERNAL_VALUES(:,D)
    END DO
    WHERE (INTERNAL_VALUES(:,:) .LT. DISCONTINUITY)
       INTERNAL_VALUES(:,:) = INTERNAL_VALUES(:,:) * SMALL_SLOPE
    END WHERE
    ! Compute the next set of internal values with a rectified activation.
    DO L = 1, MNS-1
       LP1 = L+1
       CALL DGEMM('N', 'N', SIZE(INPUTS,2), MDS, MDS, 1.0_REAL64, &
            INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
            INTERNAL_PARAMS(:,:,L), SIZE(INTERNAL_PARAMS,1), &
            0.0_REAL64, INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1))
       DO D = 1, MDS
          INTERNAL_VALUES(:,D) =  INTERNAL_BIAS(D,L) + INTERNAL_VALUES(:,D)
       END DO
       WHERE (INTERNAL_VALUES(:,:) .LT. DISCONTINUITY)
          INTERNAL_VALUES(:,:) = INTERNAL_VALUES(:,:) * SMALL_SLOPE
       END WHERE
    END DO
    ! Compute the output.
    CALL DGEMM('T', 'T', MDO, SIZE(INPUTS,2), MDS, 1.0_REAL64, &
         OUTPUT_PARAMS(:,:), SIZE(OUTPUT_PARAMS,1), &
         INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
         0.0_REAL64, OUTPUTS(:,:), SIZE(OUTPUTS,1))
    DO D = 1, MDO
       OUTPUTS(:,D) = OUTPUT_BIAS(D) + OUTPUTS(:,D)
    END DO
  END SUBROUTINE EVALUATE


  ! Compute the gradient of the sum of squared error of this regression
  ! model with respect to its parameters given input and output pairs.
  SUBROUTINE SSE_GRADIENT(INPUTS, OUTPUTS, SUM_SQUARED_ERROR, &
       INPUT_PARAMS_GRADIENT, INPUT_BIAS_GRADIENT, &
       INTERNAL_PARAMS_GRADIENT, INTERNAL_BIAS_GRADIENT, &
       OUTPUT_PARAMS_GRADIENT, OUTPUT_BIAS_GRADIENT)
    ! Data values stored with contiguous points: shape = (MDI,N)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: OUTPUTS
    ! Sum (over all data) squared error (summed over dimensions).
    REAL(KIND=REAL64), INTENT(INOUT) :: SUM_SQUARED_ERROR
    ! Model parameter gradients.
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:,:)   :: INPUT_PARAMS_GRADIENT
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:)     :: INPUT_BIAS_GRADIENT
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:,:,:) :: INTERNAL_PARAMS_GRADIENT
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:,:)   :: INTERNAL_BIAS_GRADIENT
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:,:)   :: OUTPUT_PARAMS_GRADIENT
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:)     :: OUTPUT_BIAS_GRADIENT     
    ! Internal values and lists of which ones are nonzero.
    REAL(KIND=REAL64), DIMENSION(SIZE(INPUTS,2),MDO)     :: OUTPUT_GRADIENT
    REAL(KIND=REAL64), DIMENSION(SIZE(INPUTS,2),MDS,MNS) :: INTERNAL_VALUES
    LOGICAL, DIMENSION(SIZE(INPUTS,2),MDS) :: UNDER_SPLITS
    ! Number of points.
    REAL(KIND=REAL64) :: N
    INTEGER :: D
    ! Evaluate the model at all data points, store outputs in "OUTPUT_GRADIENT"
    CALL EVALUATE_BATCH()
    ! Compute the new gradient.
    DO D = 1, MDO
       OUTPUT_GRADIENT(:,D) = OUTPUT_GRADIENT(:,D) - OUTPUTS(D,:)
    END DO
    ! Get the number of points.
    N = REAL(SIZE(INPUTS,2), REAL64)
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
      EXTERNAL :: DGEMM
      CALL CPU_TIME(START_TIME)
      ! Compute the input transformation.
      CALL DGEMM('T', 'N', SIZE(INPUTS,2), MDS, MDI, 1.0_REAL64, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INPUT_PARAMS(:,:), SIZE(INPUT_PARAMS,1), &
           0.0_REAL64, INTERNAL_VALUES(:,:,1), SIZE(INTERNAL_VALUES,1))
      CALL CPU_TIME(STOP_TIME)
      DGEMM_TIMER = DGEMM_TIMER + (STOP_TIME - START_TIME)
      CALL CPU_TIME(START_TIME)
      DO D = 1, MDS
         INTERNAL_VALUES(:,D,1) = INPUT_BIAS(D) + INTERNAL_VALUES(:,D,1)
      END DO
      WHERE (INTERNAL_VALUES(:,:,1) .LT. DISCONTINUITY)
         INTERNAL_VALUES(:,:,1) = INTERNAL_VALUES(:,:,1) * SMALL_SLOPE
      END WHERE
      CALL CPU_TIME(STOP_TIME)
      OTHER_TIMER = OTHER_TIMER + (STOP_TIME - START_TIME)
      ! Compute the next set of internal values with a rectified activation.
      DO L = 1, MNS-1
         LP1 = L+1
         CALL CPU_TIME(START_TIME)
         CALL DGEMM('N', 'N', SIZE(INPUTS,2), MDS, MDS, 1.0_REAL64, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_PARAMS(:,:,L), SIZE(INTERNAL_PARAMS,1), &
              0.0_REAL64, INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1))
         CALL CPU_TIME(STOP_TIME)
         DGEMM_TIMER = DGEMM_TIMER + (STOP_TIME - START_TIME)
         CALL CPU_TIME(START_TIME)
         DO D = 1, MDS
            INTERNAL_VALUES(:,D,LP1) =  INTERNAL_BIAS(D,L) + INTERNAL_VALUES(:,D,LP1)
         END DO
         WHERE (INTERNAL_VALUES(:,:,LP1) .LT. DISCONTINUITY)
            INTERNAL_VALUES(:,:,LP1) = INTERNAL_VALUES(:,:,LP1) * SMALL_SLOPE
         END WHERE
         CALL CPU_TIME(STOP_TIME)
         OTHER_TIMER = OTHER_TIMER + (STOP_TIME - START_TIME)
      END DO
      CALL CPU_TIME(START_TIME)
      ! Compute the output.
      CALL DGEMM('N', 'N', SIZE(INPUTS,2), MDO, MDS, 1.0_REAL64, &
           INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_PARAMS(:,:), SIZE(OUTPUT_PARAMS,1), &
           0.0_REAL64, OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1))
      CALL CPU_TIME(STOP_TIME)
      DGEMM_TIMER = DGEMM_TIMER + (STOP_TIME - START_TIME)
      CALL CPU_TIME(START_TIME)
      DO D = 1, MDO
         OUTPUT_GRADIENT(:,D) = OUTPUT_BIAS(D) + OUTPUT_GRADIENT(:,D)
      END DO
      CALL CPU_TIME(STOP_TIME)
      OTHER_TIMER = OTHER_TIMER + (STOP_TIME - START_TIME)
    END SUBROUTINE EVALUATE_BATCH

    SUBROUTINE GRADIENT_BATCH()
      ! L   - layer index
      ! LP1 - layer index "plus 1" -> "P1"
      INTEGER :: L, LP1
      EXTERNAL :: DGEMM
      ! Compute the average output gradient.
      CALL CPU_TIME(START_TIME)
      OUTPUT_BIAS_GRADIENT(:) = SUM(OUTPUT_GRADIENT(:,:), 1) / N &
           + OUTPUT_BIAS_GRADIENT(:)
      CALL CPU_TIME(STOP_TIME)
      OTHER_TIMER = OTHER_TIMER + (STOP_TIME - START_TIME)
      ! Compute output parameter gradients.
      CALL CPU_TIME(START_TIME)
      CALL DGEMM('T', 'N', MDS, MDO, SIZE(INPUTS,2), 1.0_REAL64 / N, &
           INTERNAL_VALUES(:,:,MNS), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           1.0_REAL64, OUTPUT_PARAMS_GRADIENT(:,:), SIZE(OUTPUT_PARAMS_GRADIENT,1))
      CALL CPU_TIME(STOP_TIME)
      DGEMM_TIMER = DGEMM_TIMER + (STOP_TIME - START_TIME)
      ! Compute the gradient for the last internal vector space values.
      CALL CPU_TIME(START_TIME)
      UNDER_SPLITS(:,:) = INTERNAL_VALUES(:,:,MNS) .LT. DISCONTINUITY
      CALL CPU_TIME(STOP_TIME)
      OTHER_TIMER = OTHER_TIMER + (STOP_TIME - START_TIME)
      CALL CPU_TIME(START_TIME)
      CALL DGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDO, 1.0_REAL64, &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           OUTPUT_PARAMS(:,:), SIZE(OUTPUT_PARAMS,1), &
           0.0_REAL64, INTERNAL_VALUES(:,:,MNS), SIZE(INTERNAL_VALUES,1))
      CALL CPU_TIME(STOP_TIME)
      DGEMM_TIMER = DGEMM_TIMER + (STOP_TIME - START_TIME)
      CALL CPU_TIME(START_TIME)
      WHERE (UNDER_SPLITS(:,:))
         INTERNAL_VALUES(:,:,MNS) = INTERNAL_VALUES(:,:,MNS) * SMALL_SLOPE
      END WHERE
      CALL CPU_TIME(STOP_TIME)
      OTHER_TIMER = OTHER_TIMER + (STOP_TIME - START_TIME)
      ! Cycle over all internal layers.
      INTERNAL_REPRESENTATIONS : DO L = MNS-1, 1, -1
         LP1 = L+1
         CALL CPU_TIME(START_TIME)
         ! Compute the bias gradient.
         INTERNAL_BIAS_GRADIENT(:,L) = SUM(INTERNAL_VALUES(:,:,LP1), 1) / N &
              + INTERNAL_BIAS_GRADIENT(:,L)
         CALL CPU_TIME(STOP_TIME)
         OTHER_TIMER = OTHER_TIMER + (STOP_TIME - START_TIME)         
         CALL CPU_TIME(START_TIME)
         ! Compute the gradient with respect to each output and all inputs.
         CALL DGEMM('T', 'N', MDS, MDS, SIZE(INPUTS,2), 1.0_REAL64 / N, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              1.0_REAL64, INTERNAL_PARAMS_GRADIENT(:,:,L), SIZE(INTERNAL_PARAMS_GRADIENT,1))
         CALL CPU_TIME(STOP_TIME)
         DGEMM_TIMER = DGEMM_TIMER + (STOP_TIME - START_TIME)
         ! Propogate the gradient to the immediately preceding layer.
         CALL CPU_TIME(START_TIME)
         UNDER_SPLITS(:,:) = INTERNAL_VALUES(:,:,L) .LT. DISCONTINUITY 
         CALL CPU_TIME(STOP_TIME)
         OTHER_TIMER = OTHER_TIMER + (STOP_TIME - START_TIME)
         CALL CPU_TIME(START_TIME)
         CALL DGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDS, 1.0_REAL64, &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_PARAMS(:,:,L), SIZE(INTERNAL_PARAMS,1), &
              0.0_REAL64, INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1))
         CALL CPU_TIME(STOP_TIME)
         DGEMM_TIMER = DGEMM_TIMER + (STOP_TIME - START_TIME)
         CALL CPU_TIME(START_TIME)
         WHERE (UNDER_SPLITS(:,:))
            INTERNAL_VALUES(:,:,L) = INTERNAL_VALUES(:,:,L) * SMALL_SLOPE
         END WHERE
         CALL CPU_TIME(STOP_TIME)
         OTHER_TIMER = OTHER_TIMER + (STOP_TIME - START_TIME)
      END DO INTERNAL_REPRESENTATIONS
      ! Compute the input parameter bias values.
      CALL CPU_TIME(START_TIME)
      INPUT_BIAS_GRADIENT(:) = SUM(INTERNAL_VALUES(:,:,1), 1) / N &
           + INPUT_BIAS_GRADIENT(:)
      CALL CPU_TIME(STOP_TIME)
      OTHER_TIMER = OTHER_TIMER + (STOP_TIME - START_TIME)
      ! Compute the gradient of all input parameters.
      !   [the INPUTS are transposed already, shape = (MDI,N)]
      CALL CPU_TIME(START_TIME)
      CALL DGEMM('N', 'N', MDI, MDS, SIZE(INPUTS,2), 1.0_REAL64 / N, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INTERNAL_VALUES(:,:,1), SIZE(INTERNAL_VALUES,1), &
           1.0_REAL64, INPUT_PARAMS_GRADIENT(:,:), SIZE(INPUT_PARAMS_GRADIENT,1))
      CALL CPU_TIME(STOP_TIME)
      DGEMM_TIMER = DGEMM_TIMER + (STOP_TIME - START_TIME)
    END SUBROUTINE GRADIENT_BATCH

  END SUBROUTINE SSE_GRADIENT


  ! Fit input / output pairs by minimizing mean squared error.
  SUBROUTINE MINIMIZE_MSE(X, Y, STEPS, BATCH_SIZE, NUM_THREADS, &
       DROPOUT, MEAN_SQUARED_ERROR)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: X
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: Y
    INTEGER,           INTENT(IN), OPTIONAL :: STEPS, BATCH_SIZE, NUM_THREADS
    REAL,              INTENT(IN), OPTIONAL :: DROPOUT
    REAL(KIND=REAL64), INTENT(OUT) :: MEAN_SQUARED_ERROR
    !  gradient step arrays
    REAL(KIND=REAL64), DIMENSION(MDI,MDS)       :: INPUT_PARAMS_GRAD,    IPG_MEAN,  IPG_SQ_MEAN
    REAL(KIND=REAL64), DIMENSION(MDS)           :: INPUT_BIAS_GRAD,      IBG_MEAN,  IBG_SQ_MEAN
    REAL(KIND=REAL64), DIMENSION(MDS,MDS,MNS-1) :: INTERNAL_PARAMS_GRAD, INPG_MEAN, INPG_SQ_MEAN
    REAL(KIND=REAL64), DIMENSION(MDS,MNS-1)     :: INTERNAL_BIAS_GRAD,   INBG_MEAN, INBG_SQ_MEAN
    REAL(KIND=REAL64), DIMENSION(MDS,MDO)       :: OUTPUT_PARAMS_GRAD,   OPG_MEAN,  OPG_SQ_MEAN
    REAL(KIND=REAL64), DIMENSION(MDO)           :: OUTPUT_BIAS_GRAD,     OBG_MEAN,  OBG_SQ_MEAN
    !  local integers
    INTEGER :: BS, I, NX, NT, S, BATCH_START, BATCH_END, TO_DROP, STEPS_WITH_DROPOUT
    REAL(KIND=REAL64) :: RNX, BATCHES
    ! Function that is defined by OpenMP.
    INTERFACE
       FUNCTION OMP_GET_MAX_THREADS()
         INTEGER :: OMP_GET_MAX_THREADS
       END FUNCTION OMP_GET_MAX_THREADS
    END INTERFACE
    ! Number of points.
    NX = SIZE(X,2)
    RNX = REAL(NX, REAL64)
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
    ! Set the dropout rate and number of training steps with dropout.
    IF (PRESENT(DROPOUT)) THEN
       IF ((DROPOUT .GT. 0.0) .AND. (MDS .GT. 1) .AND. (S .GT. 1)) THEN
          TO_DROP = MAX(1,INT(DROPOUT * REAL(MDS)))
          STEPS_WITH_DROPOUT = MAX(1, 1 + S/2)
          MDS = MDS - TO_DROP
       ELSE
          TO_DROP = 0
          STEPS_WITH_DROPOUT = 0
       END IF
    ELSE
       TO_DROP = 0
       STEPS_WITH_DROPOUT = 0
    END IF
    ! Set the average step sizes.
    IPG_MEAN(:,:) = 0.0_REAL64
    IBG_MEAN(:)   = 0.0_REAL64
    INPG_MEAN(:,:,:) = 0.0_REAL64
    INBG_MEAN(:,:)   = 0.0_REAL64
    OPG_MEAN(:,:) = 0.0_REAL64
    OBG_MEAN(:)   = 0.0_REAL64
    ! Set the average step size magnitudes.
    IPG_SQ_MEAN(:,:) = 0.0_REAL64
    IBG_SQ_MEAN(:)   = 0.0_REAL64
    INPG_SQ_MEAN(:,:,:) = 0.0_REAL64
    INBG_SQ_MEAN(:,:)   = 0.0_REAL64
    OPG_SQ_MEAN(:,:) = 0.0_REAL64
    OBG_SQ_MEAN(:)   = 0.0_REAL64
    ! Iterate, taking steps with the average gradient over all data.
    DO I = 1, S
       ! Do dropout if that is requested.
       IF (I .LT. STEPS_WITH_DROPOUT) THEN
          CALL DROP_INTERNAL_NODES(TO_DROP)
       ELSE IF (I .EQ. STEPS_WITH_DROPOUT) THEN
          MDS = SIZE(INTERNAL_PARAMS,1)
       END IF
       ! Compute the average gradient over all points.
       MEAN_SQUARED_ERROR = 0.0_REAL64
       ! Set gradients to zero initially.
       INPUT_PARAMS_GRAD(:,:) = 0.0_REAL64
       INPUT_BIAS_GRAD(:) = 0.0_REAL64
       INTERNAL_PARAMS_GRAD(:,:,:) = 0.0_REAL64
       INTERNAL_BIAS_GRAD(:,:) = 0.0_REAL64
       OUTPUT_PARAMS_GRAD(:,:) = 0.0_REAL64
       OUTPUT_BIAS_GRAD(:) = 0.0_REAL64
       ! Count the number of batches.
       BATCHES = 0.0_REAL64
       !$OMP PARALLEL DO NUM_THREADS(NT) SCHEDULE(STATIC) &
       !$OMP&  PRIVATE(BATCH_END, START_TIME, STOP_TIME) &
       !$OMP&  REDUCTION(+: DGEMM_TIMER, OTHER_TIMER, &
       !$OMP&  BATCHES, MEAN_SQUARED_ERROR, &
       !$OMP&  INPUT_PARAMS_GRAD, INPUT_BIAS_GRAD, &
       !$OMP&  INTERNAL_PARAMS_GRAD, INTERNAL_BIAS_GRAD, &
       !$OMP&  OUTPUT_PARAMS_GRAD, OUTPUT_BIAS_GRAD)
       DO BATCH_START = 1, NX, BS
          BATCHES = BATCHES + 1.0_REAL64
          BATCH_END = MIN(NX, BATCH_START+BS-1)
          ! Compute the gradient from all data.
          CALL SSE_GRADIENT(X(:,BATCH_START:BATCH_END), &
               Y(:,BATCH_START:BATCH_END), &
               MEAN_SQUARED_ERROR, &
               INPUT_PARAMS_GRAD, INPUT_BIAS_GRAD, &
               INTERNAL_PARAMS_GRAD, INTERNAL_BIAS_GRAD, &
               OUTPUT_PARAMS_GRAD, OUTPUT_BIAS_GRAD)
       END DO
       !$OMP END PARALLEL DO
       ! Convert the sum of squared errors into the mean squared error.
       MEAN_SQUARED_ERROR = MEAN_SQUARED_ERROR / RNX
       ! Convert the average gradients.
       INPUT_PARAMS_GRAD = INPUT_PARAMS_GRAD / BATCHES
       INPUT_BIAS_GRAD = INPUT_BIAS_GRAD / BATCHES
       INTERNAL_PARAMS_GRAD = INTERNAL_PARAMS_GRAD / BATCHES
       INTERNAL_BIAS_GRAD = INTERNAL_BIAS_GRAD / BATCHES
       OUTPUT_PARAMS_GRAD = OUTPUT_PARAMS_GRAD / BATCHES
       OUTPUT_BIAS_GRAD = OUTPUT_BIAS_GRAD / BATCHES
       ! Update the steps for all different parameters.
       CALL COMPUTE_STEP(SIZE(INPUT_PARAMS_GRAD), INPUT_PARAMS_GRAD, IPG_MEAN, IPG_SQ_MEAN)
       CALL COMPUTE_STEP(SIZE(INPUT_BIAS_GRAD), INPUT_BIAS_GRAD, IBG_MEAN, IBG_SQ_MEAN)
       CALL COMPUTE_STEP(SIZE(INTERNAL_PARAMS_GRAD), INTERNAL_PARAMS_GRAD, INPG_MEAN, INPG_SQ_MEAN)
       CALL COMPUTE_STEP(SIZE(INTERNAL_BIAS_GRAD), INTERNAL_BIAS_GRAD, INBG_MEAN, INBG_SQ_MEAN)
       CALL COMPUTE_STEP(SIZE(OUTPUT_PARAMS_GRAD), OUTPUT_PARAMS_GRAD, OPG_MEAN, OPG_SQ_MEAN)
       CALL COMPUTE_STEP(SIZE(OUTPUT_BIAS_GRAD), OUTPUT_BIAS_GRAD, OBG_MEAN, OBG_SQ_MEAN)
       ! Shrink the gradient to take smaller steps, and take
       !  the gradient step based on the average gradient.
       INPUT_PARAMS    = INPUT_PARAMS    - INPUT_PARAMS_GRAD
       INPUT_BIAS      = INPUT_BIAS      - INPUT_BIAS_GRAD
       INTERNAL_PARAMS = INTERNAL_PARAMS - INTERNAL_PARAMS_GRAD
       INTERNAL_BIAS   = INTERNAL_BIAS   - INTERNAL_BIAS_GRAD
       OUTPUT_PARAMS   = OUTPUT_PARAMS   - OUTPUT_PARAMS_GRAD
       OUTPUT_BIAS     = OUTPUT_BIAS     - OUTPUT_BIAS_GRAD
    END DO

  CONTAINS

    ! Compute the current step factor as the running average step,
    ! inversely scaled by the average magnitude of the step.
    SUBROUTINE COMPUTE_STEP(N, CURR_STEP, STEP_MEAN, STEP_SQ_MEAN)
      INTEGER, INTENT(IN) :: N
      REAL(KIND=REAL64), DIMENSION(1:N), INTENT(INOUT) :: CURR_STEP
      REAL(KIND=REAL64), DIMENSION(1:N), INTENT(INOUT) :: STEP_MEAN
      REAL(KIND=REAL64), DIMENSION(1:N), INTENT(INOUT) :: STEP_SQ_MEAN
      STEP_MEAN = 0.1 * CURR_STEP + 0.9 * STEP_MEAN
      STEP_SQ_MEAN = 0.001 * CURR_STEP**2 + 0.999 * STEP_SQ_MEAN
      WHERE (STEP_SQ_MEAN .EQ. 0.0)
         STEP_SQ_MEAN = EPSILON(STEP_SQ_MEAN(1))
      END WHERE
      ! Start scaling by step magnitude once enough data is collected.
      IF (I .GT. 10) THEN
         CURR_STEP = 0.001 * STEP_MEAN / SQRT(STEP_SQ_MEAN)
      ! Otherwise, use a scaled-down gradient descent step.
      ELSE
         CURR_STEP = 0.01 * STEP_MEAN
      END IF
    END SUBROUTINE COMPUTE_STEP

  END SUBROUTINE MINIMIZE_MSE


END MODULE PLRM


!2021-03-19 20:45:48
!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PRINT *, ''                                       !
    ! PRINT *, 'S ', S                                  !
    ! PRINT *, 'NT', NT                                 !
    ! PRINT *, 'BS', BS                                 !
    ! PRINT *, 'TO_DROP           ', TO_DROP            !
    ! PRINT *, 'STEPS_WITH_DROPOUT', STEPS_WITH_DROPOUT !
    ! PRINT *, 'NX ', NX                                !
    ! PRINT *, 'RNX', RNX                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
