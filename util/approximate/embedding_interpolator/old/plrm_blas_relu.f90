! A piecewise linear regression model.
MODULE PLRM
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE
  ! Value used for rectification, default is 0.
  REAL(KIND=REAL64), PARAMETER :: RECTIFY_AT = 0.0_REAL64
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


  ! Evaluate the piecewise linear regression model, do not store any
  !  of the intermediate states (used 
  SUBROUTINE EVALUATE_ONE(INPUT, OUTPUT)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: INPUT
    ! Record of values after the last transformation.
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(:) :: OUTPUT
    ! Internal values and lists of which ones are nonzero.
    REAL(KIND=REAL64), DIMENSION(MDS) :: INTERNAL_VALUES
    INTEGER, DIMENSION(MDS) :: NONZEROS
    ! Compute the number of nodes (to know where "bias" value is).
    INTEGER :: I, IP1, J, K
    ! Compute the input layer.
    INTERNAL_VALUES(:) = MAX(RECTIFY_AT, INPUT_BIAS(:) + &
         MATMUL(INPUT(:), INPUT_PARAMS(:,:)))
    ! Store the number of nonzero values and their locations.
    K = 0
    DO J = 1, MDS
       IF (INTERNAL_VALUES(J) .NE. 0.0_REAL64) THEN
          K = K + 1
          NONZEROS(K) = J
       END IF
    END DO
    ! Compute the next set of internal values with a rectified activation.
    DO I = 1, MNS-1
       IP1 = I+1
       INTERNAL_VALUES(:) = MAX(RECTIFY_AT, INTERNAL_BIAS(:,I) + &
            MATMUL(INTERNAL_VALUES(NONZEROS(1:K)), &
            INTERNAL_PARAMS(NONZEROS(1:K),:,I)))
       ! Store the number of nonzero values and their locations.
       K = 0
       DO J = 1, MDS
          IF (INTERNAL_VALUES(J) .NE. 0.0_REAL64) THEN
             K = K + 1
             NONZEROS(K) = J
          END IF
       END DO
    END DO
    ! Compute the output.
    OUTPUT(:) = OUTPUT_BIAS(:) + &
         MATMUL(INTERNAL_VALUES(NONZEROS(1:K)), &
         OUTPUT_PARAMS(NONZEROS(1:K),:))
  END SUBROUTINE EVALUATE_ONE


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
    LOGICAL, DIMENSION(SIZE(INPUTS,2),MDS) :: ZEROS
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
      ! Compute the input transformation.
      ! INTERNAL_VALUES(:,:,1) = MATMUL(TRANSPOSE(INPUTS(:,:)), INPUT_PARAMS(:,:))
      CALL DGEMM('T', 'N', SIZE(INPUTS,2), MDS, MDI, 1.0_REAL64, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INPUT_PARAMS(:,:), SIZE(INPUT_PARAMS,1), &
           0.0_REAL64, INTERNAL_VALUES(:,:,1), SIZE(INTERNAL_VALUES,1))
      DO D = 1, MDS
         INTERNAL_VALUES(:,D,1) = MAX(RECTIFY_AT, &
              INPUT_BIAS(D) + INTERNAL_VALUES(:,D,1))
      END DO
      ! Compute the next set of internal values with a rectified activation.
      DO L = 1, MNS-1
         LP1 = L+1
         ! INTERNAL_VALUES(:,:,LP1) = MATMUL(INTERNAL_VALUES(:,:,L), INTERNAL_PARAMS(:,:,L))
         CALL DGEMM('N', 'N', SIZE(INPUTS,2), MDS, MDS, 1.0_REAL64, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_PARAMS(:,:,L), SIZE(INTERNAL_PARAMS,1), &
              0.0_REAL64, INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1))
         DO D = 1, MDS
            INTERNAL_VALUES(:,D,LP1) = MAX(RECTIFY_AT, &
                 INTERNAL_BIAS(D,L) + INTERNAL_VALUES(:,D,LP1))
         END DO
      END DO
      ! Compute the output.
      ! OUTPUT_GRADIENT(:,:) = MATMUL(INTERNAL_VALUES(:,:,MNS), OUTPUT_PARAMS(:,:))
      CALL DGEMM('N', 'N', SIZE(INPUTS,2), MDO, MDS, 1.0_REAL64, &
           INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_PARAMS(:,:), SIZE(OUTPUT_PARAMS,1), &
           0.0_REAL64, OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1))
      DO D = 1, MDO
         OUTPUT_GRADIENT(:,D) = OUTPUT_BIAS(D) + OUTPUT_GRADIENT(:,D)
      END DO
    END SUBROUTINE EVALUATE_BATCH

    SUBROUTINE GRADIENT_BATCH()
      ! L   - layer index
      ! LP1 - layer index "plus 1" -> "P1"
      INTEGER :: L, LP1
      EXTERNAL :: DGEMM
      ! Compute the average output gradient.
      OUTPUT_BIAS_GRADIENT(:) = SUM(OUTPUT_GRADIENT(:,:), 1) / N &
           + OUTPUT_BIAS_GRADIENT(:)
      ! Compute output parameter gradients.
      CALL DGEMM('T', 'N', MDS, MDO, SIZE(INPUTS,2), 1.0_REAL64 / N, &
           INTERNAL_VALUES(:,:,MNS), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           1.0_REAL64, OUTPUT_PARAMS_GRADIENT(:,:), SIZE(OUTPUT_PARAMS_GRADIENT,1))
      ! Compute the gradient for the last internal vector space values.
      ZEROS(:,:) = INTERNAL_VALUES(:,:,MNS) .LE. RECTIFY_AT
      CALL DGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDO, 1.0_REAL64, &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           OUTPUT_PARAMS(:,:), SIZE(OUTPUT_PARAMS,1), &
           0.0_REAL64, INTERNAL_VALUES(:,:,MNS), SIZE(INTERNAL_VALUES,1))
      WHERE (ZEROS(:,:))
         INTERNAL_VALUES(:,:,MNS) = 0.0_REAL64
      END WHERE
      ! Cycle over all internal layers.
      INTERNAL_REPRESENTATIONS : DO L = MNS-1, 1, -1
         LP1 = L+1
         ! Compute the bias gradient.
         INTERNAL_BIAS_GRADIENT(:,L) = SUM(INTERNAL_VALUES(:,:,LP1), 1) / N &
              + INTERNAL_BIAS_GRADIENT(:,L)
         ! Compute the gradient with respect to each output and all inputs.
         CALL DGEMM('T', 'N', MDS, MDS, SIZE(INPUTS,2), 1.0_REAL64 / N, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              1.0_REAL64, INTERNAL_PARAMS_GRADIENT(:,:,L), SIZE(INTERNAL_PARAMS_GRADIENT,1))
         ! Propogate the gradient to the immediately preceding layer.
         ZEROS(:,:) = INTERNAL_VALUES(:,:,L) .LE. RECTIFY_AT 
         CALL DGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDS, 1.0_REAL64, &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_PARAMS(:,:,L), SIZE(INTERNAL_PARAMS,1), &
              0.0_REAL64, INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1))
         WHERE (ZEROS(:,:))
            INTERNAL_VALUES(:,:,L) = 0.0_REAL64
         END WHERE
      END DO INTERNAL_REPRESENTATIONS
      ! Compute the input parameter bias values.
      INPUT_BIAS_GRADIENT(:) = SUM(INTERNAL_VALUES(:,:,1), 1) / N &
           + INPUT_BIAS_GRADIENT(:)
      ! Compute the gradient of all input parameters.
      !   [the INPUTS are transposed already, shape = (MDI,N)]
      CALL DGEMM('N', 'N', MDI, MDS, SIZE(INPUTS,2), 1.0_REAL64 / N, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INTERNAL_VALUES(:,:,1), SIZE(INTERNAL_VALUES,1), &
           1.0_REAL64, INPUT_PARAMS_GRADIENT(:,:), SIZE(INPUT_PARAMS_GRADIENT,1))
    END SUBROUTINE GRADIENT_BATCH

  END SUBROUTINE SSE_GRADIENT


  ! Fit input / output pairs by minimizing mean squared error.
  SUBROUTINE MINIMIZE_MSE(X, Y, STEPS, BATCH_SIZE, NUM_THREADS, MEAN_SQUARED_ERROR)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: X
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: Y
    INTEGER, INTENT(IN), OPTIONAL :: STEPS, BATCH_SIZE, NUM_THREADS
    REAL(KIND=REAL64), INTENT(OUT) :: MEAN_SQUARED_ERROR
    ! Local variables.
    REAL(KIND=REAL64), PARAMETER :: STEP_FACTOR = 0.01_REAL64
    !  gradient step arrays
    REAL(KIND=REAL64), DIMENSION(MDI,MDS)       :: INPUT_PARAMS_GRADIENT_STEP
    REAL(KIND=REAL64), DIMENSION(MDS)           :: INPUT_BIAS_GRADIENT_STEP
    REAL(KIND=REAL64), DIMENSION(MDS,MDS,MNS-1) :: INTERNAL_PARAMS_GRADIENT_STEP
    REAL(KIND=REAL64), DIMENSION(MDS,MNS-1)     :: INTERNAL_BIAS_GRADIENT_STEP
    REAL(KIND=REAL64), DIMENSION(MDS,MDO)       :: OUTPUT_PARAMS_GRADIENT_STEP
    REAL(KIND=REAL64), DIMENSION(MDO)           :: OUTPUT_BIAS_GRADIENT_STEP
    ! Internal values and lists of which ones are nonzero.
    REAL(KIND=REAL64), DIMENSION(SIZE(X,2),MDS,MNS) :: INTERNAL_VALUES
    !  local integers
    INTEGER :: BS, I, NX, NT, S, BATCH_START, BATCH_END
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
       NT = OMP_GET_MAX_THREADS() / 2
    END IF
    ! Set the batch size.
    IF (PRESENT(BATCH_SIZE)) THEN
       BS = MIN(BATCH_SIZE, NX)
    ELSE
       BS = (NX + NT - 1) / NT
    END IF
    ! Iterate, taking steps with the average gradient over all data.
    DO I = 1, S
       ! Compute the average gradient over all points.
       MEAN_SQUARED_ERROR = 0.0_REAL64
       ! Set gradients to zero initially.
       INPUT_PARAMS_GRADIENT_STEP(:,:) = 0.0_REAL64
       INPUT_BIAS_GRADIENT_STEP(:) = 0.0_REAL64
       INTERNAL_PARAMS_GRADIENT_STEP(:,:,:) = 0.0_REAL64
       INTERNAL_BIAS_GRADIENT_STEP(:,:) = 0.0_REAL64
       OUTPUT_PARAMS_GRADIENT_STEP(:,:) = 0.0_REAL64
       OUTPUT_BIAS_GRADIENT_STEP(:) = 0.0_REAL64
       ! Count the number of batches.
       BATCHES = 0.0_REAL64
       !$OMP PARALLEL DO PRIVATE(BATCH_END) NUM_THREADS(NT) SCHEDULE(STATIC) &
       !$OMP&  REDUCTION(+: BATCHES, MEAN_SQUARED_ERROR, &
       !$OMP&  INPUT_PARAMS_GRADIENT_STEP, INPUT_BIAS_GRADIENT_STEP, &
       !$OMP&  INTERNAL_PARAMS_GRADIENT_STEP, INTERNAL_BIAS_GRADIENT_STEP, &
       !$OMP&  OUTPUT_PARAMS_GRADIENT_STEP, OUTPUT_BIAS_GRADIENT_STEP)
       DO BATCH_START = 1, NX, BS
          BATCHES = BATCHES + 1.0_REAL64
          BATCH_END = MIN(NX, BATCH_START+BS-1)
          ! Compute the gradient from all data.
          CALL SSE_GRADIENT(X(:,BATCH_START:BATCH_END), &
               Y(:,BATCH_START:BATCH_END), &
               MEAN_SQUARED_ERROR, &
               INPUT_PARAMS_GRADIENT_STEP, INPUT_BIAS_GRADIENT_STEP, &
               INTERNAL_PARAMS_GRADIENT_STEP, INTERNAL_BIAS_GRADIENT_STEP, &
               OUTPUT_PARAMS_GRADIENT_STEP, OUTPUT_BIAS_GRADIENT_STEP)
       END DO
       !$OMP END PARALLEL DO

       ! Convert the sum of squared errors into the mean squared error.
       MEAN_SQUARED_ERROR = MEAN_SQUARED_ERROR / REAL(NX, REAL64)
       ! Shrink the gradient to take smaller steps, and take
       !  the gradient step based on the average gradient.
       INPUT_PARAMS    = INPUT_PARAMS    - STEP_FACTOR * INPUT_PARAMS_GRADIENT_STEP    / BATCHES
       INPUT_BIAS      = INPUT_BIAS      - STEP_FACTOR * INPUT_BIAS_GRADIENT_STEP      / BATCHES
       INTERNAL_PARAMS = INTERNAL_PARAMS - STEP_FACTOR * INTERNAL_PARAMS_GRADIENT_STEP / BATCHES
       INTERNAL_BIAS   = INTERNAL_BIAS   - STEP_FACTOR * INTERNAL_BIAS_GRADIENT_STEP   / BATCHES
       OUTPUT_PARAMS   = OUTPUT_PARAMS   - STEP_FACTOR * OUTPUT_PARAMS_GRADIENT_STEP   / BATCHES
       OUTPUT_BIAS     = OUTPUT_BIAS     - STEP_FACTOR * OUTPUT_BIAS_GRADIENT_STEP     / BATCHES
    END DO
  END SUBROUTINE MINIMIZE_MSE


END MODULE PLRM
