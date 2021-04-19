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
  ! Internal values needed for training.
  REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: INTERNAL_VALUES
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: NONZEROS
  INTEGER, DIMENSION(:), ALLOCATABLE :: NUM_NONZERO
  ! Set all arrays related to gradient computaiton as private to Open MP threads.
  !$OMP THREADPRIVATE(INTERNAL_VALUES)
  !$OMP THREADPRIVATE(NONZEROS)
  !$OMP THREADPRIVATE(NUM_NONZERO)
  
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
    ! Do all de/allocations for internal values in a "parallel"
    !  region so that each thread correctly allocates its own memory.
    !$OMP PARALLEL
    IF (ALLOCATED(INTERNAL_VALUES)) DEALLOCATE(INTERNAL_VALUES)
    IF (ALLOCATED(NONZEROS)) DEALLOCATE(NONZEROS)
    IF (ALLOCATED(NUM_NONZERO)) DEALLOCATE(NUM_NONZERO)
    ALLOCATE( &
         INTERNAL_VALUES(DS, NS), &
         NONZEROS(DS, NS), &
         NUM_NONZERO(NS) &
         )
    !$OMP END PARALLEL
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
  END SUBROUTINE INIT_MODEL


  ! Evaluate the piecewise linear regression model.
  SUBROUTINE EVALUATE(INPUT, OUTPUT)
    USE ISO_FORTRAN_ENV, ONLY: REAL64
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: INPUT
    ! Record of values after the last transformation.
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(:) :: OUTPUT
    ! Compute the number of nodes (to know where "bias" value is).
    INTEGER :: I, IP1, J, K
    ! Set all values to zero.
    INTERNAL_VALUES(:,:) = 0.0_REAL64
    ! Compute the input layer.
    INTERNAL_VALUES(:,1) = MAX(RECTIFY_AT, INPUT_BIAS(:) + &
         MATMUL(INPUT(:), INPUT_PARAMS(:,:)))
    ! Store the number of nonzero values and their locations.
    K = 0
    DO J = 1, MDS
       IF (INTERNAL_VALUES(J,1) .NE. 0.0_REAL64) THEN
          K = K + 1
          NONZEROS(K,1) = J
       END IF
    END DO
    NUM_NONZERO(1) = K
    ! Compute the next set of internal values with a rectified activation.
    DO I = 1, MNS-1
       IP1 = I+1
       INTERNAL_VALUES(:,IP1) = MAX(RECTIFY_AT, INTERNAL_BIAS(:,I) + &
            MATMUL(INTERNAL_VALUES(NONZEROS(1:K,I),I), &
            INTERNAL_PARAMS(NONZEROS(1:K,I),:,I)))
       ! Store the number of nonzero values and their locations.
       K = 0
       DO J = 1, MDS
          IF (INTERNAL_VALUES(J,IP1) .NE. 0.0_REAL64) THEN
             K = K + 1
             NONZEROS(K,IP1) = J
          END IF
       END DO
       NUM_NONZERO(IP1) = K
    END DO
    ! Compute the output.
    OUTPUT(:) = OUTPUT_BIAS(:) + &
         MATMUL(INTERNAL_VALUES(NONZEROS(1:K,MNS),MNS), &
         OUTPUT_PARAMS(NONZEROS(1:K,MNS),:))
  END SUBROUTINE EVALUATE


  ! Compute the gradient of this regression model with respect to its parameters.
  SUBROUTINE GRADIENT(INPUTS, OUTPUT_GRADIENT, &
       INPUT_PARAMS_GRADIENT, INPUT_BIAS_GRADIENT, &
       INTERNAL_PARAMS_GRADIENT, INTERNAL_BIAS_GRADIENT, &
       OUTPUT_PARAMS_GRADIENT, OUTPUT_BIAS_GRADIENT)
    USE ISO_FORTRAN_ENV, ONLY: REAL64
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: INPUTS
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: OUTPUT_GRADIENT
    ! Model parameter gradients.
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:,:)   :: INPUT_PARAMS_GRADIENT
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:)     :: INPUT_BIAS_GRADIENT
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:,:,:) :: INTERNAL_PARAMS_GRADIENT
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:,:)   :: INTERNAL_BIAS_GRADIENT
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:,:)   :: OUTPUT_PARAMS_GRADIENT
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:)     :: OUTPUT_BIAS_GRADIENT     
    INTEGER :: I, IP1, J, J1, K, K1
    ! Transfer the output bias graident.
    OUTPUT_BIAS_GRADIENT(:) = OUTPUT_BIAS_GRADIENT(:) + OUTPUT_GRADIENT(:)
    ! Compute output parameter gradients.
    K = NUM_NONZERO(MNS)
    DO I = 1, MDO
       OUTPUT_PARAMS_GRADIENT(NONZEROS(1:K,MNS),I) = &
            OUTPUT_PARAMS_GRADIENT(NONZEROS(1:K,MNS),I) + &
            INTERNAL_VALUES(NONZEROS(1:K,MNS),MNS) * OUTPUT_GRADIENT(I)
    END DO
    ! Compute the gradient for the last internal vector space values.
    INTERNAL_VALUES(NONZEROS(1:K,MNS),MNS) = MATMUL(&
         OUTPUT_PARAMS(NONZEROS(1:K,MNS),:), OUTPUT_GRADIENT(:))
    ! Cycle over all internal layers.
    INTERNAL_REPRESENTATIONS : DO I = MNS-1, 1, -1
       IP1 = I+1
       K = NUM_NONZERO(IP1)       
       ! Compute the bias gradient.
       INTERNAL_BIAS_GRADIENT(NONZEROS(1:K,IP1),I) = &
            INTERNAL_BIAS_GRADIENT(NONZEROS(1:K,IP1),I) + &
            INTERNAL_VALUES(NONZEROS(1:K,IP1),IP1)
       ! Compute the parameter gradient.
       OUTPUT_DIM : DO J1 = 1, NUM_NONZERO(IP1)
          J = NONZEROS(J1,IP1)
          INPUT_DIM : DO K1 = 1, NUM_NONZERO(I)
             K = NONZEROS(K1,I)
             INTERNAL_PARAMS_GRADIENT(K,J,I) = &
                  INTERNAL_PARAMS_GRADIENT(K,J,I) + &
                  INTERNAL_VALUES(K,I) * INTERNAL_VALUES(J,IP1)
          END DO INPUT_DIM
       END DO OUTPUT_DIM
       ! Propogate the gradient backwards to the next previous layer.
       K = NUM_NONZERO(I)
       INTERNAL_VALUES(NONZEROS(1:K,I),I) = 0.0_REAL64
       DO J1 = 1, NUM_NONZERO(IP1)
          J = NONZEROS(J1,IP1)
          DO K1 = 1, NUM_NONZERO(I)
             K = NONZEROS(K1,I)
             INTERNAL_VALUES(K,I) = INTERNAL_VALUES(K,I) + &
                  INTERNAL_PARAMS(K,J,I) * INTERNAL_VALUES(J,IP1)
          END DO
       END DO
       ! TODO: Update the list of nonzer values in this layer?
    END DO INTERNAL_REPRESENTATIONS
    K = NUM_NONZERO(1)
    ! Compute the input parameter bias values.
    INPUT_BIAS_GRADIENT(NONZEROS(1:K,1)) = &
         INPUT_BIAS_GRADIENT(NONZEROS(1:K,1)) + &
         INTERNAL_VALUES(NONZEROS(1:K,1),1)
    ! Compute the next set of internal values with a rectified activation.
    DO J = 1, K
       I = NONZEROS(J,1)
       INPUT_PARAMS_GRADIENT(:,I) = INPUT_PARAMS_GRADIENT(:,I) + &
            INPUTS(:) * INTERNAL_VALUES(I,1)
    END DO
  END SUBROUTINE GRADIENT


  ! Fit input / output pairs by minimizing mean squared error.
  SUBROUTINE FIT_DATA(X, Y, STEPS, MEAN_SQUARED_ERROR)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: X
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: Y
    INTEGER, INTENT(IN), OPTIONAL :: STEPS
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
    !  storage for the gradient of the error function
    REAL(KIND=REAL64), DIMENSION(MDO) :: OUTPUT_GRADIENT
    !  local integers
    INTEGER :: I, J, N, S
    ! Set the number of steps.
    IF (PRESENT(STEPS)) THEN ; S = STEPS
    ELSE                     ; S = 1000
    END IF
    ! Set the number of points (for convenience).
    N = SIZE(X,2)
    ! Iterate, taking steps with the average gradient over all data.
    DO I = 1, S
       ! Compute the average gradient over all points.
       INPUT_PARAMS_GRADIENT_STEP    = 0.0_REAL64
       INPUT_BIAS_GRADIENT_STEP      = 0.0_REAL64
       INTERNAL_PARAMS_GRADIENT_STEP = 0.0_REAL64
       INTERNAL_BIAS_GRADIENT_STEP   = 0.0_REAL64
       OUTPUT_PARAMS_GRADIENT_STEP   = 0.0_REAL64
       OUTPUT_BIAS_GRADIENT_STEP     = 0.0_REAL64
       MEAN_SQUARED_ERROR            = 0.0_REAL64
       ! Average the gradient steps over all points.
       !$OMP PARALLEL DO PRIVATE(OUTPUT_GRADIENT) &
       !$OMP& REDUCTION(+:MEAN_SQUARED_ERROR, &
       !$OMP& INPUT_PARAMS_GRADIENT_STEP, INPUT_BIAS_GRADIENT_STEP, &
       !$OMP& INTERNAL_PARAMS_GRADIENT_STEP, INTERNAL_BIAS_GRADIENT_STEP, &
       !$OMP& OUTPUT_PARAMS_GRADIENT_STEP, OUTPUT_BIAS_GRADIENT_STEP)
       DO J = 1, N
          ! Evaluate the model at this point.
          CALL EVALUATE(X(:,J), OUTPUT_GRADIENT(:))
          ! Compute the gradient of mean squared error with respect
          !  to the output produced by the model at this data point.
          OUTPUT_GRADIENT(:) = OUTPUT_GRADIENT(:) - Y(:,J)
          ! Update the running average of squared error.
          MEAN_SQUARED_ERROR = SUM(OUTPUT_GRADIENT(:)**2) + MEAN_SQUARED_ERROR
          ! Compute the gradient with respect to this observation.
          CALL GRADIENT(X(:,J), OUTPUT_GRADIENT(:), &
               INPUT_PARAMS_GRADIENT_STEP, INPUT_BIAS_GRADIENT_STEP, &
               INTERNAL_PARAMS_GRADIENT_STEP, INTERNAL_BIAS_GRADIENT_STEP, &
               OUTPUT_PARAMS_GRADIENT_STEP, OUTPUT_BIAS_GRADIENT_STEP)
       END DO
       !$OMP END PARALLEL DO
       MEAN_SQUARED_ERROR = MEAN_SQUARED_ERROR / REAL(N,REAL64)
       ! Shrink the gradient to take smaller steps.
       INPUT_PARAMS_GRADIENT_STEP    = STEP_FACTOR * INPUT_PARAMS_GRADIENT_STEP   / REAL(N,REAL64)
       INPUT_BIAS_GRADIENT_STEP      = STEP_FACTOR * INPUT_BIAS_GRADIENT_STEP     / REAL(N,REAL64)
       INTERNAL_PARAMS_GRADIENT_STEP = STEP_FACTOR * INTERNAL_PARAMS_GRADIENT_STEP/ REAL(N,REAL64)
       INTERNAL_BIAS_GRADIENT_STEP   = STEP_FACTOR * INTERNAL_BIAS_GRADIENT_STEP  / REAL(N,REAL64)
       OUTPUT_PARAMS_GRADIENT_STEP   = STEP_FACTOR * OUTPUT_PARAMS_GRADIENT_STEP  / REAL(N,REAL64)
       OUTPUT_BIAS_GRADIENT_STEP     = STEP_FACTOR * OUTPUT_BIAS_GRADIENT_STEP    / REAL(N,REAL64)
       ! Take the gradient step based on the average gradient.
       INPUT_PARAMS    = INPUT_PARAMS    - INPUT_PARAMS_GRADIENT_STEP   
       INPUT_BIAS      = INPUT_BIAS      - INPUT_BIAS_GRADIENT_STEP     
       INTERNAL_PARAMS = INTERNAL_PARAMS - INTERNAL_PARAMS_GRADIENT_STEP
       INTERNAL_BIAS   = INTERNAL_BIAS   - INTERNAL_BIAS_GRADIENT_STEP  
       OUTPUT_PARAMS   = OUTPUT_PARAMS   - OUTPUT_PARAMS_GRADIENT_STEP  
       OUTPUT_BIAS     = OUTPUT_BIAS     - OUTPUT_BIAS_GRADIENT_STEP    
    END DO
  END SUBROUTINE FIT_DATA


END MODULE PLRM
