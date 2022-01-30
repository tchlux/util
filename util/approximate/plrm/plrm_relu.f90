! TODO:
! 
! - Validation-based training and greedy saving / model selection.
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
! - Well-spaced validation data (multiply distances in input and output?).
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
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: A, B
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: C
    ! Call external single-precision matrix-matrix multiplication
    !  (should be provided by hardware manufacturer, if not use custom).
    EXTERNAL :: SGEMM 
    CALL SGEMM(OP_A, OP_B, OUT_ROWS, OUT_COLS, INNER_DIM, &
       AB_MULT, A, A_ROWS, B, B_ROWS, C_MULT, C, C_ROWS)
    ! TODO: debug the following, appears to not be working correctly.
    ! 
    ! ! Fortran intrinsic version of general matrix multiplication routine.
    ! IF (OP_A .EQ. 'N') THEN
    !    IF (OP_B .EQ. 'N') THEN
    !       C(:,:) = MATMUL(A(:,:), B(:,:)) * AB_MULT + C_MULT * C(:,:)
    !    ELSE
    !       C(:,:) = MATMUL(A(:,:), TRANSPOSE(B(:,:))) * AB_MULT + C_MULT * C(:,:)
    !    END IF
    ! ELSE
    !    IF (OP_B .EQ. 'N') THEN
    !       C(:,:) = MATMUL(TRANSPOSE(A(:,:)), B(:,:)) * AB_MULT + C_MULT * C(:,:)
    !    ELSE
    !       C(:,:) = MATMUL(TRANSPOSE(A(:,:)), TRANSPOSE(B(:,:))) * AB_MULT + C_MULT * C(:,:)
    !    END IF
    ! END IF
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
  REAL(KIND=RT), PARAMETER :: INITIAL_OUTPUT_SCALE = 0.1_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_STEP = 0.001_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_STEP_MEAN_CHANGE = 0.1_RT  
  REAL(KIND=RT), PARAMETER :: INITIAL_STEP_CURV_CHANGE = 0.01_RT  
  REAL(KIND=RT), PARAMETER :: FASTER_RATE = 1.01_RT
  REAL(KIND=RT), PARAMETER :: SLOWER_RATE = 0.99_RT
  INTEGER,       PARAMETER :: MIN_STEPS_TO_STABILITY = 1

  ! MDI = Model dimension -- input
  ! MDS = Model dimension -- states (internal)
  ! MNS = Model number    -- states (internal)
  ! MDO = Model dimension -- output
  INTEGER :: MDI, MDS, MNS, MDO
  ! Model variables, for evaluation.
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INPUT_VECS     ! MDI, MDS
  REAL(KIND=RT), DIMENSION(:),     ALLOCATABLE :: INPUT_SHIFT    ! MDS
  REAL(KIND=RT), DIMENSION(:),     ALLOCATABLE :: INPUT_MAX      ! MDS
  REAL(KIND=RT), DIMENSION(:),     ALLOCATABLE :: INPUT_MIN      ! MDS
  REAL(KIND=RT), DIMENSION(:,:,:), ALLOCATABLE :: INTERNAL_VECS  ! MDS, MDS, MNS-1
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_SHIFT ! MDS, MNS-1
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_MAX   ! MDS, MNS-1
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_MIN   ! MDS, MNS-1
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: OUTPUT_VECS    ! MDS, MDO
  REAL(KIND=RT), DIMENSION(:),     ALLOCATABLE :: OUTPUT_SHIFT   ! MDO

  PUBLIC :: NEW_MODEL, INIT_MODEL, RANDOM_UNIT_VECTORS, MINIMIZE_MSE, &
       RESET_MIN_MAX, SET_MIN_MAX
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
    IF (ALLOCATED(INPUT_MAX))      DEALLOCATE(INPUT_MAX)
    IF (ALLOCATED(INPUT_MIN))      DEALLOCATE(INPUT_MIN)
    IF (ALLOCATED(INTERNAL_VECS))  DEALLOCATE(INTERNAL_VECS)
    IF (ALLOCATED(INTERNAL_SHIFT)) DEALLOCATE(INTERNAL_SHIFT)
    IF (ALLOCATED(INTERNAL_MAX))   DEALLOCATE(INTERNAL_MAX)
    IF (ALLOCATED(INTERNAL_MIN))   DEALLOCATE(INTERNAL_MIN)
    IF (ALLOCATED(OUTPUT_VECS))    DEALLOCATE(OUTPUT_VECS)
    IF (ALLOCATED(OUTPUT_SHIFT))   DEALLOCATE(OUTPUT_SHIFT)
    ! Allocate space for all model variables.
    ALLOCATE( &
         INPUT_VECS(DI,DS), &
         INPUT_SHIFT(DS), &
         INPUT_MAX(DS), &
         INPUT_MIN(DS), &
         INTERNAL_VECS(DS, DS, NS-1), &
         INTERNAL_SHIFT(DS, NS-1), &
         INTERNAL_MAX(DS, NS-1), &
         INTERNAL_MIN(DS, NS-1), &
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
    DO I = 1, MDS
       INPUT_SHIFT(I) = INITIAL_SHIFT_RANGE * (&
            2.0_RT * REAL((I-1),RT) / REAL((MDS-1), RT) - 1.0_RT)
       INTERNAL_SHIFT(I,:) = INPUT_SHIFT(I)
    END DO
    OUTPUT_SHIFT(:) = 0.0_RT
    ! Set the initial min and max values to be unbounded.
    CALL RESET_MIN_MAX()
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


  ! Reset the minimum and maximum values for internal nonlinearities.
  SUBROUTINE RESET_MIN_MAX()
    ! Set the initial bounds on minimum and maximum values to be extreme.
    INPUT_MAX(:) =  HUGE(INPUT_MAX(1))
    INPUT_MIN(:) = -HUGE(INPUT_MIN(1))
    INTERNAL_MAX(:,:) =  HUGE(INTERNAL_MAX(1,1))
    INTERNAL_MIN(:,:) = -HUGE(INTERNAL_MIN(1,1))
  END SUBROUTINE RESET_MIN_MAX


  ! Evaluate the piecewise linear regression model and store the minimum
  !  and maximum values observed at each internal piecewise linear function.
  SUBROUTINE SET_MIN_MAX(INPUTS)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: INPUTS
    ! Internal values.
    REAL(KIND=RT), DIMENSION(:,:), ALLOCATABLE :: TEMP_VALUES, INTERNAL_VALUES
    INTEGER :: I, D, L, N
    ! Make sure that this model has been initialized.
    IF (.NOT. ALLOCATED(INPUT_VECS)) THEN
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
       INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) + INPUT_SHIFT(D)
       INPUT_MAX(D) = MAXVAL(INTERNAL_VALUES(:,D))
       INPUT_MIN(D) = MINVAL(INTERNAL_VALUES(:,D))
       INTERNAL_VALUES(:,D) = MAX( &
            INTERNAL_VALUES(:,D), &
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
          INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) + INTERNAL_SHIFT(D,L)
          INTERNAL_MAX(D,L) = MAXVAL(INTERNAL_VALUES(:,D))
          INTERNAL_MIN(D,L) = MINVAL(INTERNAL_VALUES(:,D))
          INTERNAL_VALUES(:,D) = MAX( &
               INTERNAL_VALUES(:,D), &
               DISCONTINUITY )
       END DO
    END DO
    ! Deallocate temporary variables.
    DEALLOCATE(INTERNAL_VALUES, TEMP_VALUES)
  END SUBROUTINE SET_MIN_MAX


  ! Evaluate the piecewise linear regression model.
  SUBROUTINE EVALUATE(INPUTS, OUTPUTS, NUM_THREADS, POSITIONS, EMBEDDINGS)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: OUTPUTS
    INTEGER, OPTIONAL, INTENT(IN) :: NUM_THREADS
    REAL(KIND=RT), OPTIONAL, INTENT(OUT), DIMENSION(:,:,:) :: EMBEDDINGS ! N, MDS, MNS
    INTEGER,       OPTIONAL, INTENT(OUT), DIMENSION(:,:,:) :: POSITIONS ! N, MDS, MNS
    ! Internal values.
    REAL(KIND=RT), DIMENSION(:,:), ALLOCATABLE :: TEMP_VALUES, INTERNAL_VALUES
    INTEGER :: I, D, L, N, BATCH, NB, BN, BS, BE, BT
    ! OMP interface for determining number of batches.
    INTERFACE
       FUNCTION OMP_GET_MAX_THREADS()
         INTEGER :: OMP_GET_MAX_THREADS
       END FUNCTION OMP_GET_MAX_THREADS
    END INTERFACE
    ! Make sure that this model has been initialized.
    IF (.NOT. ALLOCATED(INPUT_VECS)) THEN
       OUTPUTS(:,:) = 0.0
       IF (PRESENT(POSITIONS)) POSITIONS(:,:,:) = 0
       IF (PRESENT(EMBEDDINGS)) EMBEDDINGS(:,:,:) = 0.0_RT
       RETURN
    END IF
    ! Get the number of points.
    N = SIZE(INPUTS,2)
    ! Allocate storage for the internal values.
    ALLOCATE(INTERNAL_VALUES(1:N,1:MDS), TEMP_VALUES(1:N,1:MDS))
    ! Set up batching for parallelization.
    IF (PRESENT(NUM_THREADS)) THEN
       NB = NUM_THREADS
    ELSE
       NB = OMP_GET_MAX_THREADS()
    END IF
    NB = MIN(N, NB)
    BN = CEILING(REAL(N) / REAL(NB))
    !$OMP PARALLEL DO PRIVATE(D, L, BS, BE, BT) IF(NB > 1)
    batch_evaluation : DO BATCH = 1, NB
       BS = BN*(BATCH-1) + 1
       BE = MIN(N, BN*BATCH)
       BT = BE-BS+1
       IF (BT .LE. 0) CYCLE batch_evaluation
       ! Compute the input transformation.
       CALL GEMM('T', 'N', BT, MDS, MDI, 1.0_RT, &
            INPUTS(:,BS:BE), SIZE(INPUTS,1), &
            INPUT_VECS(:,:), SIZE(INPUT_VECS,1), &
            0.0_RT, INTERNAL_VALUES(BS:BE,:), BT)
       DO D = 1, MDS
          INTERNAL_VALUES(BS:BE,D) = INTERNAL_VALUES(BS:BE,D) + INPUT_SHIFT(D)
          ! Store the unique positions each evaluation point land at
          !   relative to the internal basis functions.
          IF (PRESENT(POSITIONS)) THEN
             WHERE (INTERNAL_VALUES(BS:BE,D) .LT. INPUT_MIN(D))
                POSITIONS(BS:BE,D,1) = 0
             ELSEWHERE (INTERNAL_VALUES(BS:BE,D) .LT. DISCONTINUITY)
                POSITIONS(BS:BE,D,1) = 1
             ELSEWHERE (INTERNAL_VALUES(BS:BE,D) .LE. INPUT_MAX(D))
                POSITIONS(BS:BE,D,1) = 2
             ELSEWHERE
                POSITIONS(BS:BE,D,1) = 3
             END WHERE
          END IF
          ! Clip the values based on the observerd minimums and maximums.
          INTERNAL_VALUES(BS:BE,D) = MIN(INTERNAL_VALUES(BS:BE,D), INPUT_MAX(D))
          INTERNAL_VALUES(BS:BE,D) = MAX(INTERNAL_VALUES(BS:BE,D), INPUT_MIN(D))
          ! Apply the rectification.
          INTERNAL_VALUES(BS:BE,D) = MAX( &
               INTERNAL_VALUES(BS:BE,D), &
               DISCONTINUITY )
       END DO
       ! Store the embeddings for the first layer.
       IF (PRESENT(EMBEDDINGS)) EMBEDDINGS(BS:BE,:,1) = INTERNAL_VALUES(BS:BE,:)
       ! Compute the next set of internal values with a rectified activation.
       DO L = 1, MNS-1
          ! Compute all vectors.
          CALL GEMM('N', 'N', BT, MDS, MDS, 1.0_RT, &
               INTERNAL_VALUES(BS:BE,:), BT, &
               INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
               0.0_RT, TEMP_VALUES(BS:BE,:), BT)
          INTERNAL_VALUES(BS:BE,:) = TEMP_VALUES(BS:BE,:)
          ! Compute all piecewise lienar functions.
          DO D = 1, MDS
             INTERNAL_VALUES(BS:BE,D) = INTERNAL_VALUES(BS:BE,D) + INTERNAL_SHIFT(D,L)
             ! Store the unique positions each evaluation point land at
             !   relative to the internal basis functions.
             IF (PRESENT(POSITIONS)) THEN
                WHERE (INTERNAL_VALUES(BS:BE,D) .LT. INTERNAL_MIN(D,L))
                   POSITIONS(BS:BE,D,L+1) = 0
                ELSEWHERE (INTERNAL_VALUES(BS:BE,D) .LT. DISCONTINUITY)
                   POSITIONS(BS:BE,D,L+1) = 1
                ELSEWHERE (INTERNAL_VALUES(BS:BE,D) .LE. INTERNAL_MAX(D,L))
                   POSITIONS(BS:BE,D,L+1) = 2
                ELSEWHERE
                   POSITIONS(BS:BE,D,L+1) = 3
                END WHERE
             END IF
             ! Clip the values based on the observerd minimums and maximums.
             INTERNAL_VALUES(BS:BE,D) = MIN(INTERNAL_VALUES(BS:BE,D), INTERNAL_MAX(D,L))
             INTERNAL_VALUES(BS:BE,D) = MAX(INTERNAL_VALUES(BS:BE,D), INTERNAL_MIN(D,L))
             ! Apply the rectification.
             INTERNAL_VALUES(BS:BE,D) = MAX( &
                  INTERNAL_VALUES(BS:BE,D), &
                  DISCONTINUITY )
          END DO
          ! Store the embeddings for this layer.
          IF (PRESENT(EMBEDDINGS)) EMBEDDINGS(BS:BE,:,L+1) = INTERNAL_VALUES(BS:BE,:)
       END DO
       ! Return the final embedded layer if the outputs have the size of the embedding.
       IF ((MDS .NE. MDO) .AND. (SIZE(OUTPUTS,1) .EQ. MDS)) THEN
          ! Do the necessary transpose operation one dimension at a time, produce embedding.
          DO D = 1, MDS
             OUTPUTS(D,BS:BE) = INTERNAL_VALUES(BS:BE,D)
          END DO
       ! Otherwise, assume regular output computation.
       ELSE
          CALL GEMM('T', 'T', MDO, BT, MDS, 1.0_RT, &
               OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
               INTERNAL_VALUES(BS:BE,:), BT, &
               0.0_RT, OUTPUTS(:,BS:BE), SIZE(OUTPUTS,1))
          DO D = 1, MDO
             OUTPUTS(D,BS:BE) = OUTPUT_SHIFT(D) + OUTPUTS(D,BS:BE)
          END DO
       END IF
    END DO batch_evaluation
    ! Deallocate temporary variables.
    DEALLOCATE(INTERNAL_VALUES, TEMP_VALUES)
  END SUBROUTINE EVALUATE


  ! Compute the gradient of the sum of squared error of this regression
  ! model with respect to its variables given input and output pairs.
  SUBROUTINE SSE_GRADIENT(INPUTS, OUTPUTS, SUM_SQUARED_ERROR, &
       INPUT_VECS_GRADIENT, INPUT_SHIFT_GRADIENT, &
       INTERNAL_VECS_GRADIENT, INTERNAL_SHIFT_GRADIENT, &
       OUTPUT_VECS_GRADIENT, OUTPUT_SHIFT_GRADIENT )
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
    REAL(KIND=RT), DIMENSION(SIZE(INPUTS,2),MDS)     :: INTERNAL_TEMP
    ! Number of points.
    REAL(KIND=RT) :: N
    INTEGER :: D, L
    ! Get the number of points.
    N = REAL(SIZE(INPUTS,2), RT)
    ! Evaluate the model at all data points, store outputs in "OUTPUT_GRADIENT"
    CALL EVALUATE_BATCH()
    ! Compute the new gradient.
    DO D = 1, MDO
       OUTPUT_GRADIENT(:,D) = OUTPUT_GRADIENT(:,D) - OUTPUTS(D,:)
    END DO
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
         INTERNAL_VALUES(:,D,1) = MAX( INTERNAL_VALUES(:,D,1) &
              + INPUT_SHIFT(D), DISCONTINUITY)
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
            INTERNAL_VALUES(:,D,LP1) = MAX( INTERNAL_VALUES(:,D,LP1) &
                 + INTERNAL_SHIFT(D,L), DISCONTINUITY)
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
      INTEGER :: D, L, LP1, J1, J2
      REAL(KIND=RT) :: DISTANCE
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
       STEP_SIZE, KEEP_BEST, VALIDATION_SIZE, EARLY_STOP, &
       SUM_SQUARED_ERROR, RECORD)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: X
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: Y
    INTEGER,       INTENT(IN) :: STEPS
    INTEGER,       INTENT(IN), OPTIONAL :: BATCH_SIZE, NUM_THREADS, VALIDATION_SIZE
    REAL(KIND=RT), INTENT(IN), OPTIONAL :: STEP_SIZE
    LOGICAL,       INTENT(IN), OPTIONAL :: KEEP_BEST, EARLY_STOP
    REAL(KIND=RT), INTENT(OUT) :: SUM_SQUARED_ERROR
    REAL(KIND=RT), INTENT(OUT), DIMENSION(3,STEPS), OPTIONAL :: RECORD
    !  gradient step arrays, 4 copies of internal model (total of 5).
    REAL(KIND=RT), DIMENSION(MDI,MDS)       :: V1_GRAD, V1_MEAN, V1_CURV, BEST_V1
    REAL(KIND=RT), DIMENSION(MDS)           :: S1_GRAD, S1_MEAN, S1_CURV, BEST_S1
    REAL(KIND=RT), DIMENSION(MDS,MDS,MNS-1) :: V2_GRAD, V2_MEAN, V2_CURV, BEST_V2
    REAL(KIND=RT), DIMENSION(MDS,MNS-1)     :: S2_GRAD, S2_MEAN, S2_CURV, BEST_S2
    REAL(KIND=RT), DIMENSION(MDS,MDO)       :: V3_GRAD, V3_MEAN, V3_CURV, BEST_V3
    REAL(KIND=RT), DIMENSION(MDO)           :: S3_GRAD, S3_MEAN, S3_CURV, BEST_S3
    !  local variables
    LOGICAL :: REVERT_TO_BEST, DID_PRINT, ES
    LOGICAL, DIMENSION(MDS,MNS) :: TO_ZERO
    INTEGER :: BS, NV, I, L, NX, NS, NT, BATCH_START, BATCH_END, J ! TODO: remove J
    REAL(KIND=RT), DIMENSION(:,:), ALLOCATABLE :: Y_VALID
    REAL(KIND=RT) :: RNX, BATCHES, PREV_MSE, MSE, BEST_MSE, SCALAR, TOTAL_GRAD_NORM
    REAL(KIND=RT) :: STEP_FACTOR, STEP_MEAN_CHANGE, STEP_MEAN_REMAIN, &
         STEP_CURV_CHANGE, STEP_CURV_REMAIN
    REAL :: LAST_PRINT_TIME, CURRENT_TIME, WAIT_TIME
    ! Function that is defined by OpenMP.
    INTERFACE
       FUNCTION OMP_GET_MAX_THREADS()
         INTEGER :: OMP_GET_MAX_THREADS
       END FUNCTION OMP_GET_MAX_THREADS
    END INTERFACE
    ! Number of points.
    NX = SIZE(X,2)
    RNX = REAL(NX, RT)
    ! Set the "validation size".
    IF (PRESENT(VALIDATION_SIZE)) THEN
       NV = MAX(0, MIN(VALIDATION_SIZE, NX))
    ELSE
       NV = INT(0.1_RT * RNX)
    END IF
    ! Allocate storage space for evaluating the validation data.
    IF (NV .GT. 0) THEN
       ALLOCATE(Y_VALID(MDO,NV))
    END IF
    ! Update NX and RNX based on the validation size.
    NX = NX - NV
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
    ! Set the earliest allowed early-stopping step.
    IF (PRESENT(EARLY_STOP)) THEN
       ES = EARLY_STOP
    ELSE
       ES = .TRUE.
    END IF
    ! Set the initial "number of steps taken since best" counter.
    NS = 0
    ! Set the "keep best" boolean.
    IF (PRESENT(KEEP_BEST)) THEN
       REVERT_TO_BEST = KEEP_BEST
    ELSE
       REVERT_TO_BEST = (STEPS .GE. MIN_STEPS_TO_STABILITY)
    END IF
    ! Only "revert" to the best model seen if some steps are taken.
    REVERT_TO_BEST = REVERT_TO_BEST .AND. (STEPS .GT. 0)
    ! Store the start time of this routine (to make sure updates can
    !    be shown to the user at a reasonable frequency).
    CALL CPU_TIME(LAST_PRINT_TIME)
    DID_PRINT = .FALSE.
    WAIT_TIME = 2.0 * NT ! 2 seconds (times number of threads)
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
       SUM_SQUARED_ERROR = 0.0_RT
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
       !$OMP&  REDUCTION(+: BATCHES, SUM_SQUARED_ERROR, &
       !$OMP&  V1_GRAD, S1_GRAD, V2_GRAD, S2_GRAD, V3_GRAD, S3_GRAD)
       DO BATCH_START = 1, NX, BS
          BATCHES = BATCHES + 1.0_RT
          BATCH_END = MIN(NX, BATCH_START+BS-1)
          ! Sum the gradient over all data batches.
          CALL SSE_GRADIENT(X(:,BATCH_START:BATCH_END), &
               Y(:,BATCH_START:BATCH_END), SUM_SQUARED_ERROR, &
               V1_GRAD(:,:), S1_GRAD(:), &
               V2_GRAD(:,:,:), S2_GRAD(:,:), &
               V3_GRAD(:,:), S3_GRAD(:))
       END DO
       !$OMP END PARALLEL DO

       ! Convert the sum of squared errors into the mean squared error.
       MSE = SUM_SQUARED_ERROR / RNX
       ! Update the step factor based on model improvement.
       IF (MSE .LE. PREV_MSE) THEN
          STEP_FACTOR = STEP_FACTOR * FASTER_RATE
          STEP_MEAN_CHANGE = STEP_MEAN_CHANGE * SLOWER_RATE
          STEP_MEAN_REMAIN = 1.0_RT - STEP_MEAN_CHANGE
          STEP_CURV_CHANGE = STEP_CURV_CHANGE * SLOWER_RATE
          STEP_CURV_REMAIN = 1.0_RT - STEP_CURV_CHANGE
       ELSE
          STEP_FACTOR = STEP_FACTOR * SLOWER_RATE
          STEP_MEAN_CHANGE = STEP_MEAN_CHANGE * FASTER_RATE
          STEP_MEAN_REMAIN = 1.0_RT - STEP_MEAN_CHANGE
          STEP_CURV_CHANGE = STEP_CURV_CHANGE * FASTER_RATE
          STEP_CURV_REMAIN = 1.0_RT - STEP_CURV_CHANGE
       END IF
       ! Store the previous error for tracking the best-so-far.
       PREV_MSE = MSE
       
       ! Compute the mean squared error on the validation points.
       IF (NV .GT. 0) THEN
          CALL EVALUATE(X(:,NX+1:NX+NV), Y_VALID(:,:))
          Y_VALID(:,:) = (Y_VALID(:,:) - Y(:,NX+1:NX+NV))**2
          MSE = SUM(Y_VALID) / NV
       END IF

       ! Record that a step was taken.
       NS = NS + 1
       ! Update the saved "best" model based on error.
       IF (MSE .LT. BEST_MSE) THEN
          NS             = 0
          BEST_MSE       = MSE
          IF (REVERT_TO_BEST) THEN
             BEST_V1(:,:)   = INPUT_VECS(:,:)
             BEST_S1(:)     = INPUT_SHIFT(:)
             BEST_V2(:,:,:) = INTERNAL_VECS(:,:,:)
             BEST_S2(:,:)   = INTERNAL_SHIFT(:,:)
             BEST_V3(:,:)   = OUTPUT_VECS(:,:)
             BEST_S3(:)     = OUTPUT_SHIFT(:)
          END IF
       END IF

       ! Update the steps for all of the model variables.
       CALL COMPUTE_STEP(SIZE(INPUT_VECS(:,:)),      INPUT_VECS(:,:),      V1_GRAD(:,:),   V1_MEAN(:,:),   V1_CURV(:,:))
       CALL COMPUTE_STEP(SIZE(INPUT_SHIFT(:)),       INPUT_SHIFT(:),       S1_GRAD(:),     S1_MEAN(:),     S1_CURV(:))
       CALL COMPUTE_STEP(SIZE(INTERNAL_VECS(:,:,:)), INTERNAL_VECS(:,:,:), V2_GRAD(:,:,:), V2_MEAN(:,:,:), V2_CURV(:,:,:))
       CALL COMPUTE_STEP(SIZE(INTERNAL_SHIFT(:,:)),  INTERNAL_SHIFT(:,:),  S2_GRAD(:,:),   S2_MEAN(:,:),   S2_CURV(:,:))
       CALL COMPUTE_STEP(SIZE(OUTPUT_VECS(:,:)),     OUTPUT_VECS(:,:),     V3_GRAD(:,:),   V3_MEAN(:,:),   V3_CURV(:,:))
       CALL COMPUTE_STEP(SIZE(OUTPUT_SHIFT(:)),      OUTPUT_SHIFT(:),      S3_GRAD(:),     S3_MEAN(:),     S3_CURV(:))

       ! Record the 2-norm of the step that was taken (the GRAD variables were updated).
       IF (PRESENT(RECORD)) THEN
          ! Store the mean squared error at this iteration.
          RECORD(1,I) = PREV_MSE
          ! Store the validation mean squared error.
          RECORD(2,I) = MSE
          ! Store the norm of the step that was taken.
          TOTAL_GRAD_NORM = 0.0_RT
          TOTAL_GRAD_NORM = TOTAL_GRAD_NORM + SUM(V1_GRAD(:,:)**2)
          TOTAL_GRAD_NORM = TOTAL_GRAD_NORM + SUM(S1_GRAD(:)**2)
          TOTAL_GRAD_NORM = TOTAL_GRAD_NORM + SUM(V2_GRAD(:,:,:)**2)
          TOTAL_GRAD_NORM = TOTAL_GRAD_NORM + SUM(S2_GRAD(:,:)**2)
          TOTAL_GRAD_NORM = TOTAL_GRAD_NORM + SUM(V3_GRAD(:,:)**2)
          TOTAL_GRAD_NORM = TOTAL_GRAD_NORM + SUM(S3_GRAD(:)**2)
          RECORD(3,I) = SQRT(TOTAL_GRAD_NORM)
       END IF

       ! Early stop if we don't expect to see a better solution
       !  by the time the fit operation is complete.
       IF (ES .AND. (NS .GT. STEPS - I)) EXIT fit_loop

       ! Maintain a constant max-norm across the magnitue of input and internal vectors.
       SCALAR = SQRT(MAXVAL(SUM(INPUT_VECS(:,:)**2, 1)))
       INPUT_VECS(:,:) = INPUT_VECS(:,:) / SCALAR
       DO L = 1, MNS-1
          SCALAR = SQRT(MAXVAL(SUM(INTERNAL_VECS(:,:,L)**2, 1)))
          INTERNAL_VECS(:,:,L) = INTERNAL_VECS(:,:,L) / SCALAR
       END DO

       ! Write an update about step and convergence to the command line.
       CALL CPU_TIME(CURRENT_TIME)
       IF (CURRENT_TIME - LAST_PRINT_TIME .GT. WAIT_TIME) THEN
          IF (DID_PRINT) THEN
             DO J = 1, 25
                WRITE (*,'(A)', ADVANCE='NO') CHAR(8) ! <- backspace over the message
             END DO
          END IF
          WRITE (*,'(I6,"  (",F6.3,") [",F6.3,"]")', ADVANCE='NO') I, MSE, BEST_MSE
          LAST_PRINT_TIME = CURRENT_TIME
          DID_PRINT = .TRUE.
       END IF

    END DO fit_loop

    ! Restore the best model seen so far (if enough steps were taken).
    IF (REVERT_TO_BEST) THEN
       MSE                  = BEST_MSE
       INPUT_VECS(:,:)      = BEST_V1(:,:)  
       INPUT_SHIFT(:)       = BEST_S1(:)    
       INTERNAL_VECS(:,:,:) = BEST_V2(:,:,:)
       INTERNAL_SHIFT(:,:)  = BEST_S2(:,:)  
       OUTPUT_VECS(:,:)     = BEST_V3(:,:)  
       OUTPUT_SHIFT(:)      = BEST_S3(:)    
    END IF

    ! Set the minimum and maximum values observved at fit time.
    CALL SET_MIN_MAX(X(:,:))

    ! Erase the printed message if one was produced.
    IF (DID_PRINT) THEN
       DO J = 1, 25
          WRITE (*,'(A)', ADVANCE='NO') CHAR(8)
       END DO
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
      ! Set the step as the mean direction (over the past few steps).
      CURR_STEP(:) = STEP_MEAN(:)
      ! Start scaling by step magnitude by curvature once enough data is collected.
      IF (I .GE. MIN_STEPS_TO_STABILITY) THEN
         CURR_STEP(:) = CURR_STEP(:) / SQRT(STEP_CURV(:))
      END IF
      ! Take the gradient steps (based on the computed "step" above).
      PARAMS(:) = PARAMS(:) - CURR_STEP(:) * STEP_FACTOR
    END SUBROUTINE COMPUTE_STEP

  END SUBROUTINE MINIMIZE_MSE

END MODULE PLRM


!2022-01-30 10:19:19
!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ! CALL RANDOM_NUMBER(INPUT_SHIFT(:))                                                  !
    ! ! CALL RANDOM_NUMBER(INTERNAL_SHIFT(:,:))                                             !
    ! ! ! Assuming unit standard deviaiton, a shift range will capture N                    !
    ! ! !  standard deviations of data in the first layer, unclear in later                 !
    ! ! !  later layers.                                                                    !
    ! ! INPUT_SHIFT(:) = INITIAL_SHIFT_RANGE * (INPUT_SHIFT(:) * 2.0_RT - 1.0_RT)           !
    ! ! INTERNAL_SHIFT(:,:) = INITIAL_SHIFT_RANGE * (INTERNAL_SHIFT(:,:) * 2.0_RT - 1.0_RT) !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
