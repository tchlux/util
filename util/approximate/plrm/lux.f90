! TODO:
! 
! - Integer assigned to embedding, unknown integer replaced with the
!   average embedding seen at training time.
! 
! - Zero mean and unit variance the input and output data inside the
!   fit routine, then insert those linear operators into the weights
!   of the model (to reduce work that needs to be done outside).
! 
! - Get stats on the internal values within the network during training.
!   - step size progression
!   - shift values
!   - vector magnitudes for each node
!   - output weights magnitude for each node
!   - internal node contributions to MSE
!   - data distributions at internal nodes (% less and greater than 0)
! 
! - Replace model variables with single internal container type.
! - Replace model configuration with single internal container type.
! - Replace model optimization variables with single internal container type.
! - Create one "PLRM" container type that holds these three types.
!   ^^ Simplify all internal routine signatures with these, allows for
!      storing the state of the optimizer when saving the model.
!      Also allows for the module to operate on multiple different
!      networks interchangeably (right now it only supports one).
! 
! - Implement contextual input transformation that applies a PLRM model
!   to all input components concatenated with learnable "position"
!   variables. Sum all input components across the outputs of this
!   model and pass that sum-vector as the input for each point.
!   Gradient of full model should pass back through this component.
!   For very large inputs, devise a way to sample subsets of positions
!    and only train on that subset of components.
! 
! - Training support for vector of assigned neighbors for each point.
! - categorical outputs with trainable vector representations
! - Train only the output, internal, or input layers.



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
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: A
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: B
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: C
    ! Call external single-precision matrix-matrix multiplication
    !  (should be provided by hardware manufacturer, if not use custom).
    EXTERNAL :: SGEMM 
    CALL SGEMM(OP_A, OP_B, OUT_ROWS, OUT_COLS, INNER_DIM, &
       AB_MULT, A, A_ROWS, B, B_ROWS, C_MULT, C, C_ROWS)
    ! ! Fortran intrinsic version of general matrix multiplication routine,
    ! !   first compute the initial values in the output matrix,
    ! C(:,:) = C_MULT * C(:)
    ! !   then compute the matrix multiplication.
    ! IF (OP_A .EQ. 'N') THEN
    !    IF (OP_B .EQ. 'N') THEN
    !       C(:,:) = C(:,:) + AB_MULT * MATMUL(A(:,:), B(:,:))
    !    ELSE
    !       C(:,:) = C(:,:) + AB_MULT * MATMUL(A(:,:), TRANSPOSE(B(:,:)))
    !    END IF
    ! ELSE
    !    IF (OP_B .EQ. 'N') THEN
    !       C(:,:) = C(:,:) + AB_MULT * MATMUL(TRANSPOSE(A(:,:)), B(:,:))
    !    ELSE
    !       C(:,:) = C(:,:) + AB_MULT * MATMUL(TRANSPOSE(A(:,:)), TRANSPOSE(B(:,:)))
    !    END IF
    ! END IF
  END SUBROUTINE GEMM

END MODULE MATRIX_MULTIPLICATION


! A piecewise linear regression model.
MODULE PLRM
  USE ISO_FORTRAN_ENV, ONLY: RT => REAL32
  USE MATRIX_MULTIPLICATION, ONLY: GEMM
  IMPLICIT NONE

  PUBLIC :: MODEL_CONFIG, NEW_MODEL_CONFIG, INIT_MODEL, &
       RANDOM_UNIT_VECTORS, RESET_MIN_MAX, SET_MIN_MAX, MINIMIZE_MSE
  PRIVATE :: MODEL_GRADIENT

  ! Parameters that define model properties.
  !   MDI = Model dimension -- input
  !   MDS = Model dimension -- states (internal)
  !   MNS = Model number    -- states (internal)
  !   MDO = Model dimension -- output
  TYPE, BIND(C) :: MODEL_CONFIG
     ! Size related parameters.
     INTEGER :: MDI
     INTEGER :: MDS = 32
     INTEGER :: MNS = 8
     INTEGER :: MDO
     ! Approximation settings for integer inputs.
     INTEGER :: MNE = 0
     INTEGER :: MDE = 0
     ! Summary numbers that are computed.
     INTEGER :: TOTAL_SIZE
     INTEGER :: NUM_VARS
     ! Index subsets of total size vector naming scheme:
     !   S__ -> start,   E__ -> end
     !   _I_ -> input,   _S_ -> states, _O_ -> output, _E_ -> embedding
     !   __V -> vectors, __S -> shifts,  __N -> min,   __X -> max
     INTEGER :: SIV, EIV, SIS, EIS, SIN, EIN, SIX, EIX
     INTEGER :: SSV, ESV, SSS, ESS, SSN, ESN, SSX, ESX
     INTEGER :: SOV, EOV
     INTEGER :: SEV, EEV
     ! Function parameter.
     REAL(KIND=RT) :: DISCONTINUITY = 0.0_RT
     ! Initialization related parameters.
     REAL(KIND=RT) :: INITIAL_SHIFT_RANGE = 1.0_RT
     REAL(KIND=RT) :: INITIAL_OUTPUT_SCALE = 0.1_RT
     REAL(KIND=RT) :: INITIAL_STEP = 0.001_RT
     REAL(KIND=RT) :: INITIAL_STEP_MEAN_CHANGE = 0.1_RT  
     REAL(KIND=RT) :: INITIAL_STEP_CURV_CHANGE = 0.01_RT  
     ! Optimization related parameters.
     REAL(KIND=RT) :: FASTER_RATE = 1.01_RT
     REAL(KIND=RT) :: SLOWER_RATE = 0.99_RT
     INTEGER       :: MIN_STEPS_TO_STABILITY = 1
     INTEGER       :: NUM_THREADS = 1
     LOGICAL       :: KEEP_BEST = .TRUE.
     LOGICAL       :: EARLY_STOP = .TRUE.
  END TYPE MODEL_CONFIG

  ! Function that is defined by OpenMP.
  INTERFACE
     FUNCTION OMP_GET_MAX_THREADS()
       INTEGER :: OMP_GET_MAX_THREADS
     END FUNCTION OMP_GET_MAX_THREADS
  END INTERFACE

CONTAINS

  ! Generate a model configuration given state parameters for the model.
  SUBROUTINE NEW_MODEL_CONFIG(MDI, MDO, MDS, MNS, MNE, &
       MDE, NUM_THREADS, CONFIG)
     ! Size related parameters.
     INTEGER, INTENT(IN) :: MDI
     INTEGER, INTENT(IN) :: MDO
     INTEGER, OPTIONAL, INTENT(IN) :: MDS
     INTEGER, OPTIONAL, INTENT(IN) :: MNS
     INTEGER, OPTIONAL, INTENT(IN) :: MNE
     INTEGER, OPTIONAL, INTENT(IN) :: MDE
     INTEGER, OPTIONAL, INTENT(IN) :: NUM_THREADS
     ! Output
     TYPE(MODEL_CONFIG), INTENT(OUT) :: CONFIG
     ! MDS
     IF (PRESENT(MDS)) CONFIG%MDS = MDS
     ! MNS
     IF (PRESENT(MNS)) CONFIG%MNS = MNS
     ! MNE
     IF (PRESENT(MNE)) CONFIG%MNE = MNE
     ! MDE
     IF (PRESENT(MDE)) THEN
        CONFIG%MDE = MDE
     ELSE IF (CONFIG%MNE .GT. 0) THEN
        CONFIG%MDE = MAX(1, CEILING(LOG(REAL(CONFIG%MNE,RT)) / LOG(2.0_RT)))
     END IF
     ! NUM_THREADS
     IF (PRESENT(NUM_THREADS)) THEN
        CONFIG%NUM_THREADS = NUM_THREADS
     ELSE
        CONFIG%NUM_THREADS = OMP_GET_MAX_THREADS()
     END IF
     ! Declare all required configurations.
     CONFIG%MDI = MDI + CONFIG%MDE
     CONFIG%MDO = MDO
     ! Compute indices related to the parameter vector for this model.
     CONFIG%TOTAL_SIZE = 0
     !   input vecs
     CONFIG%SIV = 1 + CONFIG%TOTAL_SIZE
     CONFIG%EIV = CONFIG%SIV-1  +  CONFIG%MDI * CONFIG%MDS
     CONFIG%TOTAL_SIZE = CONFIG%EIV
     !   input shift
     CONFIG%SIS = 1 + CONFIG%TOTAL_SIZE
     CONFIG%EIS = CONFIG%SIS-1  +  CONFIG%MDS
     CONFIG%TOTAL_SIZE = CONFIG%EIS
     !   state vecs
     CONFIG%SSV = 1 + CONFIG%TOTAL_SIZE
     CONFIG%ESV = CONFIG%SSV-1  +  CONFIG%MDS * CONFIG%MDS * (CONFIG%MNS-1)
     CONFIG%TOTAL_SIZE = CONFIG%ESV
     !   state shift
     CONFIG%SSS = 1 + CONFIG%TOTAL_SIZE
     CONFIG%ESS = CONFIG%SSS-1  +  CONFIG%MDS * (CONFIG%MNS-1)
     CONFIG%TOTAL_SIZE = CONFIG%ESS
     !   output vecs
     CONFIG%SOV = 1 + CONFIG%TOTAL_SIZE
     CONFIG%EOV = CONFIG%SOV-1  +  CONFIG%MDS * CONFIG%MDO
     CONFIG%TOTAL_SIZE = CONFIG%EOV
     !   embedding vecs
     CONFIG%SEV = 1 + CONFIG%TOTAL_SIZE
     CONFIG%EEV = CONFIG%SEV-1  +  CONFIG%MDE * CONFIG%MNE
     CONFIG%TOTAL_SIZE = CONFIG%EEV
     !   number of variables
     CONFIG%NUM_VARS = CONFIG%TOTAL_SIZE
     !   input min
     CONFIG%SIN = 1 + CONFIG%TOTAL_SIZE
     CONFIG%EIN = CONFIG%SIN-1  +  CONFIG%MDS
     CONFIG%TOTAL_SIZE = CONFIG%EIN
     !   input max
     CONFIG%SIX = 1 + CONFIG%TOTAL_SIZE
     CONFIG%EIX = CONFIG%SIX-1  +  CONFIG%MDS
     CONFIG%TOTAL_SIZE = CONFIG%EIX
     !   state min
     CONFIG%SSN = 1 + CONFIG%TOTAL_SIZE
     CONFIG%ESN = CONFIG%SSN-1  +  CONFIG%MDS * (CONFIG%MNS-1)
     CONFIG%TOTAL_SIZE = CONFIG%ESN
     !   state max
     CONFIG%SSX = 1 + CONFIG%TOTAL_SIZE
     CONFIG%ESX = CONFIG%SSX-1  +  CONFIG%MDS * (CONFIG%MNS-1)
     CONFIG%TOTAL_SIZE = CONFIG%ESX
  END SUBROUTINE NEW_MODEL_CONFIG


  ! Initialize the weights for a model, optionally provide a random seed.
  SUBROUTINE INIT_MODEL(CONFIG, MODEL, SEED)
    TYPE(MODEL_CONFIG), INTENT(IN) :: CONFIG
    REAL(KIND=RT), INTENT(IN), DIMENSION(:) :: MODEL
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
    ! Unpack the model vector into its parts.
    CALL UNPACKED_INIT_MODEL(&
         MODEL(CONFIG%SIV:CONFIG%EIV), &
         MODEL(CONFIG%SIS:CONFIG%EIS), &
         MODEL(CONFIG%SIN:CONFIG%EIN), &
         MODEL(CONFIG%SIX:CONFIG%EIX), &
         MODEL(CONFIG%SSV:CONFIG%ESV), &
         MODEL(CONFIG%SSS:CONFIG%ESS), &
         MODEL(CONFIG%SSN:CONFIG%ESN), &
         MODEL(CONFIG%SSX:CONFIG%ESX), &
         MODEL(CONFIG%SOV:CONFIG%EOV), &
         MODEL(CONFIG%SEV:CONFIG%EEV))

    CONTAINS

      ! Initialize the model after unpacking it into its constituent parts.
      SUBROUTINE UNPACKED_INIT_MODEL(&
           INPUT_VECS, INPUT_SHIFT, INPUT_MIN, INPUT_MAX, &
           STATE_VECS, STATE_SHIFT, STATE_MIN, STATE_MAX, &
           OUTPUT_VECS, EMBEDDINGS)
        REAL(KIND=RT), DIMENSION(CONFIG%MDI, CONFIG%MDS) :: INPUT_VECS
        REAL(KIND=RT), DIMENSION(CONFIG%MDS) :: INPUT_SHIFT, INPUT_MIN, INPUT_MAX
        REAL(KIND=RT), DIMENSION(CONFIG%MDS, CONFIG%MDS, CONFIG%MNS-1) :: STATE_VECS
        REAL(KIND=RT), DIMENSION(CONFIG%MDS, CONFIG%MNS-1) :: STATE_SHIFT, STATE_MIN, STATE_MAX
        REAL(KIND=RT), DIMENSION(CONFIG%MDS, CONFIG%MDO) :: OUTPUT_VECS
        REAL(KIND=RT), DIMENSION(CONFIG%MDE, CONFIG%MNE) :: EMBEDDINGS
        ! Generate well spaced random unit-length vectors (no scaling biases)
        !  for all initial variables in the input, internal, output, and embedings.
        CALL RANDOM_UNIT_VECTORS(INPUT_VECS(:,:))
        DO I = 1, CONFIG%MNS-1
           CALL RANDOM_UNIT_VECTORS(STATE_VECS(:,:,I))
        END DO
        CALL RANDOM_UNIT_VECTORS(OUTPUT_VECS(:,:))
        IF (CONFIG%MNE .GT. 0) THEN
           CALL RANDOM_UNIT_VECTORS(EMBEDDINGS(:,:))
        END IF
        ! Make the output vectors have very small magnitude initially.
        OUTPUT_VECS(:,:) = OUTPUT_VECS(:,:) * CONFIG%INITIAL_OUTPUT_SCALE
        ! Generate random shifts for inputs and internal layers, zero
        !  shift for the output layer (first two will be rescaled).
        DO I = 1, CONFIG%MDS
           INPUT_SHIFT(I) = CONFIG%INITIAL_SHIFT_RANGE * (&
                2.0_RT * REAL((I-1),RT) &
                / MAX(1.0_RT, REAL((CONFIG%MDS-1), RT)))
           STATE_SHIFT(I,:) = INPUT_SHIFT(I)
        END DO
        ! Set the initial min and max values to be unbounded.
        CALL RESET_MIN_MAX(INPUT_MIN, INPUT_MAX, STATE_MIN, STATE_MAX)
      END SUBROUTINE UNPACKED_INIT_MODEL

  END SUBROUTINE INIT_MODEL

  ! Generate randomly distributed vectors on the N-sphere.
  SUBROUTINE RANDOM_UNIT_VECTORS(COLUMN_VECTORS)
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: COLUMN_VECTORS
    REAL(KIND=RT), DIMENSION(SIZE(COLUMN_VECTORS,1), SIZE(COLUMN_VECTORS,2)) :: TEMP_VECS
    REAL(KIND=RT), PARAMETER :: PI = 3.141592653589793
    INTEGER :: I, J
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

    ! Orthogonalize and normalize A.
    SUBROUTINE ORTHOGONALIZE(A)
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:) :: A
      REAL(KIND=RT), DIMENSION(SIZE(A,2)) :: MULTIPLIERS
      INTEGER :: I, J
      A(:,1) = A(:,1) / NORM2(A(:,1))
      DO I = 2, SIZE(A,2)
         MULTIPLIERS(I:) = MATMUL(A(:,I-1), A(:,I:))
         DO J = I, SIZE(A,2)
            A(:,J) = A(:,J) - MULTIPLIERS(J) * A(:,I-1)
         END DO
         A(:,I) = A(:,I) / NORM2(A(:,I))
      END DO
    END SUBROUTINE ORTHOGONALIZE

  END SUBROUTINE RANDOM_UNIT_VECTORS


  ! Given integer inputs and embedding vectors, put embeddings in
  !  place of integer inputs inside of a real matrix.
  SUBROUTINE UNPACK_EMBEDDINGS(CONFIG, EMBEDDINGS, INT_INPUTS, EMBEDDED)
    TYPE(MODEL_CONFIG), INTENT(IN) :: CONFIG
    REAL(KIND=RT), INTENT(IN), DIMENSION(CONFIG%MDE,CONFIG%MNE) :: EMBEDDINGS
    INTEGER, INTENT(IN), DIMENSION(:,:) :: INT_INPUTS
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: EMBEDDED
    INTEGER :: N, D, E
    ! Add together appropriate embedding vectors based on integer inputs.
    EMBEDDED(:,:) = 0.0_RT
    DO N = 1, SIZE(INT_INPUTS,2)
       EMBEDDED(:,N) = 0.0_RT
       DO D = 1, SIZE(INT_INPUTS,1)
          E = INT_INPUTS(D,N)
          EMBEDDED(:,N) = EMBEDDED(:,N) + EMBEDDINGS(:,E)
       END DO
    END DO
  END SUBROUTINE UNPACK_EMBEDDINGS


  ! Reset the minimum and maximum values for internal nonlinearities.
  SUBROUTINE RESET_MIN_MAX(INPUT_MIN, INPUT_MAX, STATE_MIN, STATE_MAX)
    REAL(KIND=RT), DIMENSION(:)   :: INPUT_MIN, INPUT_MAX
    REAL(KIND=RT), DIMENSION(:,:) :: STATE_MIN, STATE_MAX
    ! Set the initial bounds on minimum and maximum values to be extreme.
    INPUT_MAX(:) =  HUGE(INPUT_MAX(1))
    INPUT_MIN(:) = -HUGE(INPUT_MIN(1))
    STATE_MAX(:,:) =  HUGE(STATE_MAX(1,1))
    STATE_MIN(:,:) = -HUGE(STATE_MIN(1,1))
  END SUBROUTINE RESET_MIN_MAX


  ! Evaluate the piecewise linear regression model and store the minimum
  !  and maximum values observed at each internal piecewise linear function.
  SUBROUTINE SET_MIN_MAX(CONFIG, MODEL, X, XI)
    TYPE(MODEL_CONFIG), INTENT(IN) :: CONFIG
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:) :: MODEL
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: X
    INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: XI
    ! Internal values.
    REAL(KIND=RT), DIMENSION(:,:), ALLOCATABLE :: XXI, TEMP_VALUES, STATE_VALUES
    INTEGER :: I, D, L, N, ES
    ! Get the number of points.
    N = SIZE(X,2)
    ! Allocate storage for the internal values.
    ALLOCATE(STATE_VALUES(1:N,1:CONFIG%MDS), TEMP_VALUES(1:N,1:CONFIG%MDS))
    IF (PRESENT(XI)) THEN
       ! Convert the integer embeddings into their real counterparts.
       ALLOCATE(XXI(CONFIG%MDI, SIZE(X,2)))
       ES = SIZE(X,1)+1
       ! Add the real inputs to the front of each vector.
       XXI(:ES-1,:) = XI(:,:)
       ! Put the embeddings into back of each vector.
       CALL UNPACK_EMBEDDINGS(CONFIG, MODEL(CONFIG%SEV:CONFIG%EEV), &
            XI, XXI(ES:,:))
       ! Unpack the model.
       CALL UNPACKED_SET_MIN_MAX(XXI(:,:), &
            MODEL(CONFIG%SIV:CONFIG%EIV), &
            MODEL(CONFIG%SIS:CONFIG%EIS), &
            MODEL(CONFIG%SIN:CONFIG%EIN), &
            MODEL(CONFIG%SIX:CONFIG%EIX), &
            MODEL(CONFIG%SSV:CONFIG%ESV), &
            MODEL(CONFIG%SSS:CONFIG%ESS), &
            MODEL(CONFIG%SSN:CONFIG%ESN), &
            MODEL(CONFIG%SSX:CONFIG%ESX), &
            MODEL(CONFIG%SOV:CONFIG%EOV))
       ! Deallocate the full inputs that included embeddings.
       DEALLOCATE(XXI)
    ELSE
       ! Unpack the model.
       CALL UNPACKED_SET_MIN_MAX(X(:,:), &
            MODEL(CONFIG%SIV:CONFIG%EIV), &
            MODEL(CONFIG%SIS:CONFIG%EIS), &
            MODEL(CONFIG%SIN:CONFIG%EIN), &
            MODEL(CONFIG%SIX:CONFIG%EIX), &
            MODEL(CONFIG%SSV:CONFIG%ESV), &
            MODEL(CONFIG%SSS:CONFIG%ESS), &
            MODEL(CONFIG%SSN:CONFIG%ESN), &
            MODEL(CONFIG%SSX:CONFIG%ESX), &
            MODEL(CONFIG%SOV:CONFIG%EOV))
    END IF
    ! Deallocate temporary variables.
    DEALLOCATE(STATE_VALUES, TEMP_VALUES)

  CONTAINS

    ! Initialize the model after unpacking it into its constituent parts.
    SUBROUTINE UNPACKED_SET_MIN_MAX(INPUTS, &
         INPUT_VECS, INPUT_SHIFT, INPUT_MIN, INPUT_MAX, &
         STATE_VECS, STATE_SHIFT, STATE_MIN, STATE_MAX, &
         OUTPUT_VECS)
      REAL(KIND=RT), DIMENSION(:,:) :: INPUTS
      ! Model.
      REAL(KIND=RT), DIMENSION(CONFIG%MDI, CONFIG%MDS) :: INPUT_VECS
      REAL(KIND=RT), DIMENSION(CONFIG%MDS) :: INPUT_SHIFT, INPUT_MIN, INPUT_MAX
      REAL(KIND=RT), DIMENSION(CONFIG%MDS, CONFIG%MDS, CONFIG%MNS-1) :: STATE_VECS
      REAL(KIND=RT), DIMENSION(CONFIG%MDS, CONFIG%MNS-1) :: STATE_SHIFT, STATE_MIN, STATE_MAX
      REAL(KIND=RT), DIMENSION(CONFIG%MDS, CONFIG%MDO) :: OUTPUT_VECS
      ! Compute the input transformation.
      CALL GEMM('T', 'N', N, CONFIG%MDS, CONFIG%MDI, 1.0_RT, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INPUT_VECS(:,:), SIZE(INPUT_VECS,1), &
           0.0_RT, STATE_VALUES(:,:), SIZE(STATE_VALUES,1))
      DO D = 1, CONFIG%MDS
         STATE_VALUES(:,D) = STATE_VALUES(:,D) + INPUT_SHIFT(D)
         INPUT_MAX(D) = MAXVAL(STATE_VALUES(:,D))
         INPUT_MIN(D) = MINVAL(STATE_VALUES(:,D))
         STATE_VALUES(:,D) = MAX( &
              STATE_VALUES(:,D), &
              CONFIG%DISCONTINUITY )
      END DO
      ! Compute the next set of internal values with a rectified activation.
      DO L = 1, CONFIG%MNS-1
         ! Compute all vectors.
         CALL GEMM('N', 'N', N, CONFIG%MDS, CONFIG%MDS, 1.0_RT, &
              STATE_VALUES(:,:), SIZE(STATE_VALUES,1), &
              STATE_VECS(:,:,L), SIZE(STATE_VECS,1), &
              0.0_RT, TEMP_VALUES(:,:), SIZE(TEMP_VALUES,1))
         STATE_VALUES(:,:) = TEMP_VALUES(:,:)
         DO D = 1, CONFIG%MDS
            STATE_VALUES(:,D) = STATE_VALUES(:,D) + STATE_SHIFT(D,L)
            STATE_MAX(D,L) = MAXVAL(STATE_VALUES(:,D))
            STATE_MIN(D,L) = MINVAL(STATE_VALUES(:,D))
            STATE_VALUES(:,D) = MAX( &
                 STATE_VALUES(:,D), &
                 CONFIG%DISCONTINUITY )
         END DO
      END DO
    END SUBROUTINE UNPACKED_SET_MIN_MAX

  END SUBROUTINE SET_MIN_MAX

  
  ! Evaluate the piecewise linear regression model.
  SUBROUTINE EVALUATE(CONFIG, MODEL, X, XI, Y, POSITIONS, EMBEDDINGS)
    TYPE(MODEL_CONFIG), INTENT(IN) :: CONFIG
    REAL(KIND=RT), INTENT(IN), DIMENSION(:) :: MODEL
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: X
    INTEGER, OPTIONAL, INTENT(IN), DIMENSION(:,:) :: XI
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: Y
    REAL(KIND=RT), OPTIONAL, INTENT(OUT), DIMENSION(:,:,:) :: EMBEDDINGS ! N, MDS, MNS
    INTEGER,       OPTIONAL, INTENT(OUT), DIMENSION(:,:,:) :: POSITIONS  ! N, MDS, MNS
    ! Internal values.
    REAL(KIND=RT), DIMENSION(:,:), ALLOCATABLE :: XXI, TEMP_VALUES, STATE_VALUES
    INTEGER :: I, N, BATCH, NB, BN, BS, BE, BT, ES
    ! Get the number of points.
    N = SIZE(X,2)
    ! Allocate storage for the internal values.
    ALLOCATE(STATE_VALUES(1:N,1:CONFIG%MDS), TEMP_VALUES(1:N,1:CONFIG%MDS))
    ! Set up batching for parallelization.
    NB = MIN(N, CONFIG%NUM_THREADS)
    ! Compute the number of points in each batch.
    BN = CEILING(REAL(N) / REAL(NB))
    ! One pluse the number of input components (in case embeddings are used).
    ES = SIZE(X,1)+1
    !$OMP PARALLEL DO PRIVATE(BS, BE, BT, XXI) IF(NB > 1)
    batch_evaluation : DO BATCH = 1, NB
       BS = BN*(BATCH-1) + 1
       BE = MIN(N, BN*BATCH)
       BT = BE-BS+1
       IF (BT .LE. 0) CYCLE batch_evaluation
       ! Generate copies needed to pass embeddings into model.
       IF (PRESENT(XI)) THEN
          ! Convert the integer embeddings into their real counterparts.
          ALLOCATE(XXI(CONFIG%MDI, BT))
          ! Add the real inputs to the front of each vector.
          XXI(:ES-1,:) = X(:,BS:BE)
          ! Put the embeddings into back of each vector.
          CALL UNPACK_EMBEDDINGS(CONFIG, MODEL(CONFIG%SEV:CONFIG%EEV), &
               XI(:,BS:BE), XXI(ES:,:))
          ! Unpack the model vector into its parts.
          CALL UNPACKED_EVALUATE(BS, BE, BT, &
               MODEL(CONFIG%SIV:CONFIG%EIV), &
               MODEL(CONFIG%SIS:CONFIG%EIS), &
               MODEL(CONFIG%SIN:CONFIG%EIN), &
               MODEL(CONFIG%SIX:CONFIG%EIX), &
               MODEL(CONFIG%SSV:CONFIG%ESV), &
               MODEL(CONFIG%SSS:CONFIG%ESS), &
               MODEL(CONFIG%SSN:CONFIG%ESN), &
               MODEL(CONFIG%SSX:CONFIG%ESX), &
               MODEL(CONFIG%SOV:CONFIG%EOV), &
               XXI(:,:), Y(:,BS:BE))
          ! Deallocate memory for embedding copy.
          DEALLOCATE(XXI)
       ELSE
          ! Unpack the model vector into its parts.
          CALL UNPACKED_EVALUATE(BS, BE, BT, &
               MODEL(CONFIG%SIV:CONFIG%EIV), &
               MODEL(CONFIG%SIS:CONFIG%EIS), &
               MODEL(CONFIG%SIN:CONFIG%EIN), &
               MODEL(CONFIG%SIX:CONFIG%EIX), &
               MODEL(CONFIG%SSV:CONFIG%ESV), &
               MODEL(CONFIG%SSS:CONFIG%ESS), &
               MODEL(CONFIG%SSN:CONFIG%ESN), &
               MODEL(CONFIG%SSX:CONFIG%ESX), &
               MODEL(CONFIG%SOV:CONFIG%EOV), &
               X(:,BS:BE), Y(:,BS:BE))
       END IF
    END DO batch_evaluation
    ! Deallocate temporary variables.
    DEALLOCATE(STATE_VALUES, TEMP_VALUES)


  CONTAINS

    SUBROUTINE UNPACKED_EVALUATE(BS, BE, BT, &
         INPUT_VECS, INPUT_SHIFT, INPUT_MIN, INPUT_MAX, &
         STATE_VECS, STATE_SHIFT, STATE_MIN, STATE_MAX, &
         OUTPUT_VECS, INPUTS, OUTPUTS)
      INTEGER, INTENT(IN) :: BS, BE, BT
      REAL(KIND=RT), DIMENSION(CONFIG%MDI, CONFIG%MDS) :: INPUT_VECS
      REAL(KIND=RT), DIMENSION(CONFIG%MDS) :: INPUT_SHIFT, INPUT_MIN, INPUT_MAX
      REAL(KIND=RT), DIMENSION(CONFIG%MDS, CONFIG%MDS, CONFIG%MNS-1) :: STATE_VECS
      REAL(KIND=RT), DIMENSION(CONFIG%MDS, CONFIG%MNS-1) :: STATE_SHIFT, STATE_MIN, STATE_MAX
      REAL(KIND=RT), DIMENSION(CONFIG%MDS, CONFIG%MDO) :: OUTPUT_VECS
      REAL(KIND=RT), INTENT(IN),  DIMENSION(:,:) :: INPUTS
      REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: OUTPUTS
      ! Local variables to evaluating a single batch.
      REAL(KIND=RT), DIMENSION(BT,CONFIG%MDS) :: STATE_VALUES, TEMP_VALUES
      INTEGER :: D, L
      ! Compute the input transformation.
      CALL GEMM('T', 'N', BT, CONFIG%MDS, CONFIG%MDI, 1.0_RT, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INPUT_VECS(:,:), SIZE(INPUT_VECS,1), &
           0.0_RT, STATE_VALUES(:,:), BT)
      DO D = 1, CONFIG%MDS
         STATE_VALUES(:,D) = STATE_VALUES(:,D) + INPUT_SHIFT(D)
         ! Store the unique positions each evaluation point land at
         !   relative to the internal basis functions.
         IF (PRESENT(POSITIONS)) THEN
            WHERE (STATE_VALUES(:,D) .LT. INPUT_MIN(D))
               POSITIONS(BS:BE,D,1) = 0
            ELSEWHERE (STATE_VALUES(:,D) .LT. CONFIG%DISCONTINUITY)
               POSITIONS(BS:BE,D,1) = 1
            ELSEWHERE (STATE_VALUES(:,D) .LE. INPUT_MAX(D))
               POSITIONS(BS:BE,D,1) = 2
            ELSEWHERE
               POSITIONS(BS:BE,D,1) = 3
            END WHERE
         END IF
         ! Clip the values based on the observerd minimums and maximums.
         STATE_VALUES(:,D) = MIN(STATE_VALUES(:,D), INPUT_MAX(D))
         STATE_VALUES(:,D) = MAX(STATE_VALUES(:,D), INPUT_MIN(D))
         ! Apply the rectification.
         STATE_VALUES(:,D) = MAX( &
              STATE_VALUES(:,D), &
              CONFIG%DISCONTINUITY )
      END DO
      ! Store the embeddings for the first layer.
      IF (PRESENT(EMBEDDINGS)) EMBEDDINGS(BS:BE,:,1) = STATE_VALUES(:,:)
      ! Compute the next set of internal values with a rectified activation.
      DO L = 1, CONFIG%MNS-1
         ! Compute all vectors.
         CALL GEMM('N', 'N', BT, CONFIG%MDS, CONFIG%MDS, 1.0_RT, &
              STATE_VALUES(:,:), BT, &
              STATE_VECS(:,:,L), SIZE(STATE_VECS,1), &
              0.0_RT, TEMP_VALUES(:,:), BT)
         STATE_VALUES(:,:) = TEMP_VALUES(:,:)
         ! Compute all piecewise lienar functions.
         DO D = 1, CONFIG%MDS
            STATE_VALUES(:,D) = STATE_VALUES(:,D) + STATE_SHIFT(D,L)
            ! Store the unique positions each evaluation point land at
            !   relative to the internal basis functions.
            IF (PRESENT(POSITIONS)) THEN
               WHERE (STATE_VALUES(:,D) .LT. STATE_MIN(D,L))
                  POSITIONS(BS:BE,D,L+1) = 0
               ELSEWHERE (STATE_VALUES(:,D) .LT. CONFIG%DISCONTINUITY)
                  POSITIONS(BS:BE,D,L+1) = 1
               ELSEWHERE (STATE_VALUES(:,D) .LE. STATE_MAX(D,L))
                  POSITIONS(BS:BE,D,L+1) = 2
               ELSEWHERE
                  POSITIONS(BS:BE,D,L+1) = 3
               END WHERE
            END IF
            ! Clip the values based on the observerd minimums and maximums.
            STATE_VALUES(:,D) = MIN(STATE_VALUES(:,D), STATE_MAX(D,L))
            STATE_VALUES(:,D) = MAX(STATE_VALUES(:,D), STATE_MIN(D,L))
            ! Apply the rectification.
            STATE_VALUES(:,D) = MAX( &
                 STATE_VALUES(:,D), &
                 CONFIG%DISCONTINUITY )
         END DO
         ! Store the embeddings for this layer.
         IF (PRESENT(EMBEDDINGS)) EMBEDDINGS(BS:BE,:,L+1) = STATE_VALUES(:,:)
      END DO
      ! Return the final embedded layer if the outputs have the size of the embedding.
      CALL GEMM('T', 'T', CONFIG%MDO, BT, CONFIG%MDS, 1.0_RT, &
           OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
           STATE_VALUES(:,:), BT, &
           0.0_RT, OUTPUTS(:,:), SIZE(OUTPUTS,1))
    END SUBROUTINE UNPACKED_EVALUATE
  END SUBROUTINE EVALUATE

  ! Compute the gradient of the sum of squared error of this regression
  ! model with respect to its variables given input and output pairs.
  SUBROUTINE MODEL_GRADIENT(CONFIG, MODEL, X, Y, SUM_SQUARED_ERROR, &
       MODEL_GRAD, XI, INPUT_GRAD)
    TYPE(MODEL_CONFIG), INTENT(IN) :: CONFIG
    REAL(KIND=RT), INTENT(IN), DIMENSION(:) :: MODEL
    ! Data values stored with contiguous points: shape = (MDI,N)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: X
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: Y
    ! Sum (over all data) squared error (summed over dimensions).
    REAL(KIND=RT), INTENT(INOUT) :: SUM_SQUARED_ERROR
    INTEGER, OPTIONAL, INTENT(IN), DIMENSION(:,:) :: XI
    ! Gradient of the model parameters.
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:) :: MODEL_GRAD
    REAL(KIND=RT), OPTIONAL, INTENT(OUT), DIMENSION(:,:) :: INPUT_GRAD
    ! Local allocations for computing gradient.
    REAL(KIND=RT), DIMENSION(:,:), ALLOCATABLE :: XXI
    REAL(KIND=RT), DIMENSION(SIZE(X,2),CONFIG%MDO)            :: Y_GRADIENT
    REAL(KIND=RT), DIMENSION(SIZE(X,2),CONFIG%MDS,CONFIG%MNS) :: STATE_VALUES
    REAL(KIND=RT), DIMENSION(SIZE(X,2),CONFIG%MDS)            :: STATE_TEMP
    INTEGER :: ES

    ! TODO:
    ! - Unpack integer embeddings into associated real vector values, sum
    !   them, and pass them into the model while computing the input gradient.
    ! - Average the input gradient over common batches of values to compute the
    !   gradient with respect to the individual embeddings assigned to inputs.
    ! - Add the zero mean and unit length code into the `MINIMIZE_MSE` routine
    !   then embed those shifts into the first layer of the model at training end.
    ! - Add optional rank normalization to model input and output, with a
    !   fixed number of rank positions (100). Include index sorting code. Use
    !   binary search mixed with walking sorted to transform new points.
    ! - Add well spaced sphere design code, replace random sphere usage
    !   with the well spaced sphere. Especially for embedding initialization.

    IF (.NOT. PRESENT(XI)) THEN
       ! Call an internal operator for this routine that takes the unpacked model.
       CALL UNPACKED_MODEL_GRADIENT(X(:,:), &
            MODEL(CONFIG%SIV:CONFIG%EIV), &
            MODEL(CONFIG%SIS:CONFIG%EIS), &
            MODEL(CONFIG%SSV:CONFIG%ESV), &
            MODEL(CONFIG%SSS:CONFIG%ESS), &
            MODEL(CONFIG%SOV:CONFIG%EOV), &
            MODEL_GRAD(CONFIG%SIV:CONFIG%EIV), &
            MODEL_GRAD(CONFIG%SIS:CONFIG%EIS), &
            MODEL_GRAD(CONFIG%SSV:CONFIG%ESV), &
            MODEL_GRAD(CONFIG%SSS:CONFIG%ESS), &
            MODEL_GRAD(CONFIG%SOV:CONFIG%EOV), &
            INPUT_GRAD)
    ELSE
       ! Convert the integer embeddings into their real counterparts.
       ALLOCATE(XXI(CONFIG%MDI, SIZE(X,2)))
       ES = SIZE(X,1)+1
       ! Add the real inputs to the front of each vector.
       XXI(:ES-1,:) = X(:,:)
       ! Put the embeddings into back of each vector.
       CALL UNPACK_EMBEDDINGS(CONFIG, MODEL(CONFIG%SEV:CONFIG%EEV), &
            XI(:,:), XXI(ES:,:))
       ! Call the regular gradient computation step, ensuring to get input gradients.
       CALL UNPACKED_MODEL_GRADIENT(XXI(:,:), &
            MODEL(CONFIG%SIV:CONFIG%EIV), &
            MODEL(CONFIG%SIS:CONFIG%EIS), &
            MODEL(CONFIG%SSV:CONFIG%ESV), &
            MODEL(CONFIG%SSS:CONFIG%ESS), &
            MODEL(CONFIG%SOV:CONFIG%EOV), &
            MODEL_GRAD(CONFIG%SIV:CONFIG%EIV), &
            MODEL_GRAD(CONFIG%SIS:CONFIG%EIS), &
            MODEL_GRAD(CONFIG%SSV:CONFIG%ESV), &
            MODEL_GRAD(CONFIG%SSS:CONFIG%ESS), &
            MODEL_GRAD(CONFIG%SOV:CONFIG%EOV), &
            INPUT_GRAD=XXI) ! <- overwrite "full inputs" with grad info
       ! Convert the computed input gradients into average gradients for each embedding.
       CALL UNPACKED_EMBEDDING_GRAD(&
            XI, XXI(ES:,:), MODEL_GRAD(CONFIG%SEV:CONFIG%EEV))
       ! Deallocate memory for holding embedded inputs.
       DEALLOCATE(XXI)
    END IF

  CONTAINS
    
    ! Average the gradients over all points to create the embedding gradient.
    SUBROUTINE UNPACKED_EMBEDDING_GRAD(INT_INPUTS, GRAD, EMBEDDING_GRAD)
      INTEGER, INTENT(IN), DIMENSION(:,:) :: INT_INPUTS
      REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: GRAD
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(CONFIG%MDE,CONFIG%MNE) :: EMBEDDING_GRAD
      REAL(KIND=RT), DIMENSION(CONFIG%MNE) :: COUNTS
      INTEGER :: N, D, E
      ! Accumulate the gradients for all embedding vectors.
      COUNTS(:) = 0.0_RT
      EMBEDDING_GRAD(:,:) = 0.0_RT
      DO N = 1, SIZE(INT_INPUTS,2)
         DO D = 1, SIZE(INT_INPUTS,1)
            E = INT_INPUTS(D,N)
            COUNTS(E) = COUNTS(E) + 1.0_RT
            EMBEDDING_GRAD(:,E) = EMBEDDING_GRAD(:,E) + GRAD(:,N)
         END DO
      END DO
      ! Average the embedding gradient by dividing by the sum of occurrences.
      DO E = 1, CONFIG%MNE
         IF (COUNTS(E) .GT. 0.0_RT) THEN
            EMBEDDING_GRAD(:,E) = EMBEDDING_GRAD(:,E) / COUNTS(E)
         END IF
      END DO
    END SUBROUTINE UNPACKED_EMBEDDING_GRAD

    ! Compute the model gradient.
    SUBROUTINE UNPACKED_MODEL_GRADIENT( INPUTS, &
         INPUT_VECS, INPUT_SHIFT, &
         STATE_VECS, STATE_SHIFT, OUTPUT_VECS, &
         INPUT_VECS_GRADIENT, INPUT_SHIFT_GRADIENT, &
         STATE_VECS_GRADIENT, STATE_SHIFT_GRADIENT, &
         OUTPUT_VECS_GRADIENT, INPUT_GRAD )
      REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: INPUTS
      ! Model variables.
      REAL(KIND=RT), INTENT(IN), DIMENSION(CONFIG%MDI,CONFIG%MDS) :: INPUT_VECS
      REAL(KIND=RT), INTENT(IN), DIMENSION(CONFIG%MDS) :: INPUT_SHIFT
      REAL(KIND=RT), INTENT(IN), DIMENSION(CONFIG%MDS,CONFIG%MDS,CONFIG%MNS-1) :: STATE_VECS
      REAL(KIND=RT), INTENT(IN), DIMENSION(CONFIG%MDS,CONFIG%MNS-1) :: STATE_SHIFT
      REAL(KIND=RT), INTENT(IN), DIMENSION(CONFIG%MDS,CONFIG%MDO) :: OUTPUT_VECS
      ! Model variable gradients.
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(CONFIG%MDI,CONFIG%MDS) :: INPUT_VECS_GRADIENT
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(CONFIG%MDS) :: INPUT_SHIFT_GRADIENT
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(CONFIG%MDS,CONFIG%MDS,CONFIG%MNS-1) :: STATE_VECS_GRADIENT
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(CONFIG%MDS,CONFIG%MNS-1) :: STATE_SHIFT_GRADIENT
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(CONFIG%MDS,CONFIG%MDO) :: OUTPUT_VECS_GRADIENT
      ! Input gradient (optional).
      REAL(KIND=RT), INTENT(INOUT), OPTIONAL, DIMENSION(:,:) :: INPUT_GRAD
      ! D   - dimension index
      ! L   - layer index
      ! LP1 - layer index "plus 1" -> "P1"
      INTEGER :: D, L, LP1
      ! Number of points.
      REAL(KIND=RT) :: N
      ! Get the number of points.
      N = REAL(SIZE(INPUTS,2), RT)
      ! -------------------------------------------------------------------------
      ! Evaluate the model at all data points, store outputs in "Y_GRADIENT"
      ! 
      ! Compute the input transformation.
      CALL GEMM('T', 'N', SIZE(INPUTS,2), CONFIG%MDS, CONFIG%MDI, 1.0_RT, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           INPUT_VECS(:,:), SIZE(INPUT_VECS,1), &
           0.0_RT, STATE_VALUES(:,:,1), SIZE(STATE_VALUES,1))
      DO D = 1, CONFIG%MDS
         STATE_VALUES(:,D,1) = MAX( STATE_VALUES(:,D,1) &
              + INPUT_SHIFT(D), CONFIG%DISCONTINUITY)
      END DO
      ! Compute the next set of internal values with a rectified activation.
      DO L = 1, CONFIG%MNS-1
         ! Compute vector values.
         LP1 = L+1
         CALL GEMM('N', 'N', SIZE(INPUTS,2), CONFIG%MDS, CONFIG%MDS, 1.0_RT, &
              STATE_VALUES(:,:,L), SIZE(STATE_VALUES,1), &
              STATE_VECS(:,:,L), SIZE(STATE_VECS,1), &
              0.0_RT, STATE_VALUES(:,:,LP1), SIZE(STATE_VALUES,1))
         DO D = 1, CONFIG%MDS
            STATE_VALUES(:,D,LP1) = MAX( STATE_VALUES(:,D,LP1) &
                 + STATE_SHIFT(D,L), CONFIG%DISCONTINUITY)
         END DO
      END DO
      ! Compute the output.
      CALL GEMM('N', 'N', SIZE(INPUTS,2), CONFIG%MDO, CONFIG%MDS, 1.0_RT, &
           STATE_VALUES(:,:,L), SIZE(STATE_VALUES,1), &
           OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
           0.0_RT, Y_GRADIENT(:,:), SIZE(Y_GRADIENT,1))
      ! ------------------------------------------------------------------------
      ! Compute the gradient of the model outputs, overwriting "Y_GRADIENT"
      CALL SQUARED_ERROR_GRADIENT(Y, Y_GRADIENT)
      ! Compute the current sum squared error.
      SUM_SQUARED_ERROR = SUM_SQUARED_ERROR + SUM(Y_GRADIENT(:,:)**2)
      ! ------------------------------------------------------------------------
      ! Compute the gradient of variables with respect to the "output gradient"
      CALL GEMM('T', 'N', CONFIG%MDS, CONFIG%MDO, SIZE(INPUTS,2), 1.0_RT / N, &
           STATE_VALUES(:,:,CONFIG%MNS), SIZE(STATE_VALUES,1), &
           Y_GRADIENT(:,:), SIZE(Y_GRADIENT,1), &
           1.0_RT, OUTPUT_VECS_GRADIENT(:,:), SIZE(OUTPUT_VECS_GRADIENT,1))
      ! Propogate the gradient back to the last internal vector space.
      CALL GEMM('N', 'T', SIZE(INPUTS,2), CONFIG%MDS, CONFIG%MDO, 1.0_RT, &
           Y_GRADIENT(:,:), SIZE(Y_GRADIENT,1), &
           OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
           0.0_RT, STATE_TEMP(:,:), SIZE(STATE_TEMP,1))
      ! Cycle over all internal layers.
      STATE_REPRESENTATIONS : DO L = CONFIG%MNS-1, 1, -1
         LP1 = L+1
         DO D = 1, CONFIG%MDS
            ! Propogate the error gradient back to the preceding vectors.
            WHERE (STATE_VALUES(:,D,LP1) .GT. CONFIG%DISCONTINUITY)
               STATE_VALUES(:,D,LP1) = STATE_TEMP(:,D)
            END WHERE
         END DO
         ! Compute the shift gradient.
         STATE_SHIFT_GRADIENT(:,L) = SUM(STATE_VALUES(:,:,LP1), 1) / N &
              + STATE_SHIFT_GRADIENT(:,L)
         ! Compute the gradient with respect to each output and all inputs.
         CALL GEMM('T', 'N', CONFIG%MDS, CONFIG%MDS, SIZE(INPUTS,2), 1.0_RT / N, &
              STATE_VALUES(:,:,L), SIZE(STATE_VALUES,1), &
              STATE_VALUES(:,:,LP1), SIZE(STATE_VALUES,1), &
              1.0_RT, STATE_VECS_GRADIENT(:,:,L), SIZE(STATE_VECS_GRADIENT,1))
         ! Propogate the gradient to the immediately preceding layer.
         CALL GEMM('N', 'T', SIZE(INPUTS,2), CONFIG%MDS, CONFIG%MDS, 1.0_RT, &
              STATE_VALUES(:,:,LP1), SIZE(STATE_VALUES,1), &
              STATE_VECS(:,:,L), SIZE(STATE_VECS,1), &
              0.0_RT, STATE_TEMP(:,:), SIZE(STATE_TEMP,1))
      END DO STATE_REPRESENTATIONS
      ! Compute the gradients going into the first layer.
      DO D = 1, CONFIG%MDS
         ! Propogate the error gradient back to the preceding vectors.
         WHERE (STATE_VALUES(:,D,1) .GT. CONFIG%DISCONTINUITY)
            STATE_VALUES(:,D,1) = STATE_TEMP(:,D)
         END WHERE
      END DO
      ! Compute the input shift variable gradients.
      INPUT_SHIFT_GRADIENT(:) = SUM(STATE_VALUES(:,:,1), 1) / N &
           + INPUT_SHIFT_GRADIENT(:)
      ! Compute the gradient of all input variables.
      !   [the INPUTS are transposed already, shape = (CONFIG%MDI,N)]
      CALL GEMM('N', 'N', CONFIG%MDI, CONFIG%MDS, SIZE(INPUTS,2), 1.0_RT / N, &
           INPUTS(:,:), SIZE(INPUTS,1), &
           STATE_VALUES(:,:,1), SIZE(STATE_VALUES,1), &
           1.0_RT, INPUT_VECS_GRADIENT(:,:), SIZE(INPUT_VECS_GRADIENT,1))
      ! Compute the gradient at the input if that was requested.
      IF (PRESENT(INPUT_GRAD)) THEN
         CALL GEMM('N', 'T', CONFIG%MDI, SIZE(INPUT_GRAD,2), CONFIG%MDS, 1.0_RT, &
              INPUT_VECS(:,:), SIZE(INPUT_VECS,1), &
              STATE_VALUES(:,:,1), SIZE(STATE_VALUES,1), &
              0.0_RT, INPUT_GRAD(:,:), SIZE(INPUT_GRAD,1))
      END IF
      ! --------------------------------------------------------------
    END SUBROUTINE UNPACKED_MODEL_GRADIENT

    ! Compute the sum of squared error, store the gradient in the OUTPUTS.
    !   TARGETS - row vectors containing target values
    !   OUTPUTS - column vectors containing model predictions
    SUBROUTINE SQUARED_ERROR_GRADIENT(TARGETS, OUTPUTS)
      REAL(KIND=RT), INTENT(IN),    DIMENSION(:,:) :: TARGETS
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:) :: OUTPUTS
      INTEGER :: D
      DO D = 1, SIZE(OUTPUTS,2)
         OUTPUTS(:,D) = OUTPUTS(:,D) - TARGETS(D,:)
      END DO
    END SUBROUTINE SQUARED_ERROR_GRADIENT

  END SUBROUTINE MODEL_GRADIENT


  ! Fit input / output pairs by minimizing mean squared error.
  SUBROUTINE MINIMIZE_MSE(CONFIG, MODEL, X, Y, STEPS, XI, &
       SUM_SQUARED_ERROR, RECORD)
    TYPE(MODEL_CONFIG), INTENT(IN) :: CONFIG
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:) :: MODEL
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: X
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: Y
    INTEGER,       INTENT(IN) :: STEPS
    INTEGER,       INTENT(IN), DIMENSION(:,:), OPTIONAL :: XI
    REAL(KIND=RT), INTENT(OUT) :: SUM_SQUARED_ERROR
    REAL(KIND=RT), INTENT(OUT), DIMENSION(2,STEPS), OPTIONAL :: RECORD
    !  Local variables.
    !    gradient step arrays, 4 copies of model + (num threads - 1)
    REAL(KIND=RT), DIMENSION(CONFIG%NUM_VARS) :: &
         MODEL_GRAD, MODEL_GRAD_MEAN, MODEL_GRAD_CURV, BEST_MODEL
    !    singletons
    LOGICAL :: REVERT_TO_BEST, DID_PRINT
    INTEGER :: BS, I, NX, NS, BATCH_START, BATCH_END
    REAL(KIND=RT) :: RNX, BATCHES, PREV_MSE, MSE, BEST_MSE
    REAL(KIND=RT) :: STEP_FACTOR, STEP_MEAN_CHANGE, STEP_MEAN_REMAIN, &
         STEP_CURV_CHANGE, STEP_CURV_REMAIN
    REAL :: LAST_PRINT_TIME, CURRENT_TIME, WAIT_TIME
    CHARACTER(LEN=*), PARAMETER :: RESET_LINE = REPEAT(CHAR(8),25)
    ! Number of points.
    NX = SIZE(X,2)
    RNX = REAL(NX, RT)
    ! Set the step factor.
    STEP_FACTOR = CONFIG%INITIAL_STEP
    ! Set the initial "number of steps taken since best" counter.
    NS = 0
    ! Set the batch size.
    BS = (NX + CONFIG%NUM_THREADS - 1) / CONFIG%NUM_THREADS
    ! Only "revert" to the best model seen if some steps are taken.
    REVERT_TO_BEST = REVERT_TO_BEST .AND. (STEPS .GT. 0)
    ! Store the start time of this routine (to make sure updates can
    !    be shown to the user at a reasonable frequency).
    CALL CPU_TIME(LAST_PRINT_TIME)
    DID_PRINT = .FALSE.
    WAIT_TIME = 2.0 * CONFIG%NUM_THREADS ! 2 seconds (times number of threads)
    ! Initial rates of change of mean and variance values.
    STEP_MEAN_CHANGE = CONFIG%INITIAL_STEP_MEAN_CHANGE
    STEP_MEAN_REMAIN = 1.0_RT - STEP_MEAN_CHANGE
    STEP_CURV_CHANGE = CONFIG%INITIAL_STEP_CURV_CHANGE
    STEP_CURV_REMAIN = 1.0_RT - STEP_CURV_CHANGE
    ! Initial mean squared error is "max representable value".
    PREV_MSE = HUGE(PREV_MSE)
    BEST_MSE = HUGE(BEST_MSE)
    ! Set the average step sizes.
    MODEL_GRAD_MEAN(:) = 0.0_RT
    ! Set the estiamted curvature in steps.
    MODEL_GRAD_CURV(:) = 0.0_RT
    ! Iterate, taking steps with the average gradient over all data.
    fit_loop : DO I = 1, STEPS
       ! Compute the average gradient over all points.
       SUM_SQUARED_ERROR = 0.0_RT
       ! Set gradients to zero initially.
       MODEL_GRAD(:) = 0.0_RT
       ! Count the number of batches.
       BATCHES = 0.0_RT
       !$OMP PARALLEL DO NUM_THREADS(CONFIG%NUM_THREADS) PRIVATE(BATCH_END) &
       !$OMP&  REDUCTION(+: BATCHES, SUM_SQUARED_ERROR, MODEL_GRAD)
       DO BATCH_START = 1, NX, BS
          BATCHES = BATCHES + 1.0_RT
          BATCH_END = MIN(NX, BATCH_START+BS-1)
          ! Sum the gradient over all data batches.
          IF (PRESENT(XI)) THEN
             CALL MODEL_GRADIENT(CONFIG, MODEL(:), &
                  X(:,BATCH_START:BATCH_END), &
                  Y(:,BATCH_START:BATCH_END), &
                  SUM_SQUARED_ERROR, MODEL_GRAD(:), &
                  XI=XI(:,BATCH_START:BATCH_END))
          ELSE
             CALL MODEL_GRADIENT(CONFIG, MODEL(:), &
                  X(:,BATCH_START:BATCH_END), &
                  Y(:,BATCH_START:BATCH_END), &
                  SUM_SQUARED_ERROR, MODEL_GRAD(:))
          END IF
       END DO
       !$OMP END PARALLEL DO
       ! Convert the sum of squared errors into the mean squared error.
       MSE = SUM_SQUARED_ERROR / RNX
       ! Update the step factor based on model improvement.
       IF (MSE .LE. PREV_MSE) THEN
          STEP_FACTOR = STEP_FACTOR * CONFIG%FASTER_RATE
          STEP_MEAN_CHANGE = STEP_MEAN_CHANGE * CONFIG%SLOWER_RATE
          STEP_MEAN_REMAIN = 1.0_RT - STEP_MEAN_CHANGE
          STEP_CURV_CHANGE = STEP_CURV_CHANGE * CONFIG%SLOWER_RATE
          STEP_CURV_REMAIN = 1.0_RT - STEP_CURV_CHANGE
       ELSE
          STEP_FACTOR = STEP_FACTOR * CONFIG%SLOWER_RATE
          STEP_MEAN_CHANGE = STEP_MEAN_CHANGE * CONFIG%FASTER_RATE
          STEP_MEAN_REMAIN = 1.0_RT - STEP_MEAN_CHANGE
          STEP_CURV_CHANGE = STEP_CURV_CHANGE * CONFIG%FASTER_RATE
          STEP_CURV_REMAIN = 1.0_RT - STEP_CURV_CHANGE
       END IF
       ! Store the previous error for tracking the best-so-far.
       PREV_MSE = MSE

       ! Record that a step was taken.
       NS = NS + 1
       ! Update the saved "best" model based on error.
       IF (MSE .LT. BEST_MSE) THEN
          NS = 0
          BEST_MSE = MSE
          IF (REVERT_TO_BEST) THEN
             BEST_MODEL(:) = MODEL(1:CONFIG%NUM_VARS)
          END IF
       ! Early stop if we don't expect to see a better solution
       !  by the time the fit operation is complete.
       ELSE IF (CONFIG%EARLY_STOP .AND. (NS .GT. STEPS - I)) THEN
          EXIT fit_loop
       END IF

       ! Convert the summed gradients to average gradients.
       MODEL_GRAD(:) = MODEL_GRAD(:) / BATCHES
       ! Update sliding window mean and variance calculations.
       MODEL_GRAD_MEAN(:) = STEP_MEAN_REMAIN * MODEL_GRAD_MEAN(:) &
            + STEP_MEAN_CHANGE * MODEL_GRAD(:)
       MODEL_GRAD_CURV(:) = STEP_CURV_REMAIN * MODEL_GRAD_CURV(:) &
            + STEP_CURV_CHANGE * (MODEL_GRAD_MEAN(:) - MODEL_GRAD(:))**2
       MODEL_GRAD_CURV(:) = MAX(MODEL_GRAD_CURV(:), EPSILON(STEP_FACTOR))
       ! Set the step as the mean direction (over the past few steps).
       MODEL_GRAD(:) = MODEL_GRAD_MEAN(:)
       ! Start scaling by step magnitude by curvature once enough data is collected.
       IF (I .GE. CONFIG%MIN_STEPS_TO_STABILITY) THEN
          MODEL_GRAD(:) = MODEL_GRAD(:) / SQRT(MODEL_GRAD_CURV(:))
       END IF
       ! Take the gradient steps (based on the computed "step" above).
       MODEL(1:CONFIG%NUM_VARS) = MODEL(1:CONFIG%NUM_VARS) - MODEL_GRAD(:) * STEP_FACTOR

       ! Record the 2-norm of the step that was taken (the GRAD variables were updated).
       IF (PRESENT(RECORD)) THEN
          ! Store the mean squared error at this iteration.
          RECORD(1,I) = MSE
          ! Store the norm of the step that was taken.
          RECORD(2,I) = SUM(MODEL_GRAD(:)**2)
       END IF

       ! Maintain a constant max-norm across the magnitue of input and internal vectors.
       CALL UNIT_MAX_NORM(MODEL(CONFIG%SIV:CONFIG%EIV), &
            MODEL(CONFIG%SSV:CONFIG%ESV))

       ! Write an update about step and convergence to the command line.
       CALL CPU_TIME(CURRENT_TIME)
       IF (CURRENT_TIME - LAST_PRINT_TIME .GT. WAIT_TIME) THEN
          IF (DID_PRINT) WRITE (*,'(A)',ADVANCE='NO') RESET_LINE
          WRITE (*,'(I6,"  (",F6.3,") [",F6.3,"]")', ADVANCE='NO') I, MSE, BEST_MSE
          LAST_PRINT_TIME = CURRENT_TIME
          DID_PRINT = .TRUE.
       END IF

    END DO fit_loop

    ! Restore the best model seen so far (if enough steps were taken).
    IF (REVERT_TO_BEST) THEN
       MSE                      = BEST_MSE
       MODEL(1:CONFIG%NUM_VARS) = BEST_MODEL(:)
    END IF

    ! Set the minimum and maximum values observved at fit time.
    CALL SET_MIN_MAX(CONFIG, MODEL(:), X(:,:), XI(:,:))

    ! Erase the printed message if one was produced.
    IF (DID_PRINT) WRITE (*,'(A)',ADVANCE='NO') RESET_LINE


  CONTAINS

    ! Set the input vectors and the state vectors to 
    SUBROUTINE UNIT_MAX_NORM(INPUT_VECS, STATE_VECS)
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(CONFIG%MDI,CONFIG%MDS)              :: INPUT_VECS
      REAL(KIND=RT), INTENT(INOUT), DIMENSION(CONFIG%MDS,CONFIG%MDS,CONFIG%MNS-1) :: STATE_VECS
      REAL(KIND=RT) :: SCALAR
      INTEGER :: L
      SCALAR = SQRT(MAXVAL(SUM(INPUT_VECS(:,:)**2, 1)))
      INPUT_VECS(:,:) = INPUT_VECS(:,:) / SCALAR
      DO L = 1, SIZE(STATE_VECS,3)
         SCALAR = SQRT(MAXVAL(SUM(STATE_VECS(:,:,L)**2, 1)))
         STATE_VECS(:,:,L) = STATE_VECS(:,:,L) / SCALAR
      END DO
    END SUBROUTINE UNIT_MAX_NORM

  END SUBROUTINE MINIMIZE_MSE

END MODULE PLRM

