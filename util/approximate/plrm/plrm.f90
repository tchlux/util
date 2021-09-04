! TODO:
! - Validation-based training and greedy saving / model selection.
! - Well-spaced validation data (multiply distances in input and output).
! - batch reduction for accelerated training (well spaced in different layer values and large gradients)
! - Support for vectorizing integer valued inputs (initially regular simplex, trainable).
! - Measure pairwise nearness between internal vectors based on point values.
! - Setting that produces convergence status on timed interval (3 seconds?).
! - Training support for vector of assigned neighbors for each point.
! - Support for huge input spaces (1M+ components), automatically
!   construct vectors that select a subset of input components to
!   process.
! - Support for huge number of input points (find well spaced sample, reduce).
! 
! - Get stats on the internal values within the network during training.
!   - shift values
!   - vector magnitudes for each node
!   - output weights magnitude for each node
!   - internal node contributions to MSE
!   - data distributions at internal nodes (% less and greater than 0)
! 
! - Measure pairwise distance between vectors (including shift) at
!   each data representation.
! - Measure pairwise distance between vector *values* (after applying
!   discontinuity) at each representation.
! - Add forcing function that pushes vectors away from each other based
!   on how similar their outputs at all data points are?
! - Pick a point (with probability corresponding to ratio of total error)
!   as the target for a new vector, pick the shift that captures it.
! 
! - Keep the length of every internal vector fixed at 1.
! - Orthogonalize all internal vectors at every step.
! - Substitute least value vector at each layer with one that captures
!   points with the current largest error (regress error given the
!   representation at that level).
! 
! - categorical outputs with trainable vector representations
! - Save model and load model (during training).
! - Train only the output, internal, or input layers.
! - Implement similar code in C, compare speeds.


! A piecewise linear regression model.
MODULE PLRM
  USE ISO_FORTRAN_ENV, ONLY: RT => REAL32
  IMPLICIT NONE
  ! Parameters that define model properties.
  REAL(KIND=RT), PARAMETER :: DISCONTINUITY = 0.0_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_SHIFT_RANGE = 0.9_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_FLEX = 0.001_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_OUTPUT_SCALE = 0.001_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_STEP = 0.001_RT
  REAL(KIND=RT), PARAMETER :: INITIAL_STEP_MEAN_CHANGE = 0.1_RT  
  REAL(KIND=RT), PARAMETER :: INITIAL_STEP_VAR_CHANGE = 0.01_RT  
  REAL(KIND=RT), PARAMETER :: STEP_GROWTH_RATE = 1.001_RT
  REAL(KIND=RT), PARAMETER :: STEP_SHRINK_RATE = 0.999_RT
  REAL(KIND=RT), PARAMETER :: DECAY_FACTOR = 1.0_RT
  INTEGER,       PARAMETER :: MIN_STEPS_TO_STABILITY = 1

  ! MDI = Model dimension -- input
  ! MDS = Model dimension -- states (internal)
  ! MNS = Model number -- states (internal)
  ! MDO = Model dimension -- output
  INTEGER :: MDI, MDS, MNS, MDO
  ! Model parameters, for evaluation.
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INPUT_VECS
  REAL(KIND=RT), DIMENSION(:),     ALLOCATABLE :: INPUT_SHIFT
  REAL(KIND=RT), DIMENSION(:),     ALLOCATABLE :: INPUT_FLEX
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_MEAN
  REAL(KIND=RT), DIMENSION(:,:,:), ALLOCATABLE :: INTERNAL_VECS
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_SHIFT
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_FLEX
  REAL(KIND=RT), DIMENSION(:,:),   ALLOCATABLE :: OUTPUT_VECS
  REAL(KIND=RT), DIMENSION(:),     ALLOCATABLE :: OUTPUT_SHIFT  
  
  PUBLIC :: NEW_MODEL, INIT_MODEL, MINIMIZE_MSE
  PRIVATE :: SSE_GRADIENT

  ! Externally provided matrix-matrix multiplication routine.
  EXTERNAL :: SGEMM
  ! EXTERNAL :: DGEMM

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
    IF (ALLOCATED(INPUT_FLEX))     DEALLOCATE(INPUT_FLEX)
    IF (ALLOCATED(INTERNAL_MEAN))  DEALLOCATE(INTERNAL_MEAN)
    IF (ALLOCATED(INTERNAL_VECS))  DEALLOCATE(INTERNAL_VECS)
    IF (ALLOCATED(INTERNAL_SHIFT)) DEALLOCATE(INTERNAL_SHIFT)
    IF (ALLOCATED(INTERNAL_FLEX))  DEALLOCATE(INTERNAL_FLEX)
    IF (ALLOCATED(OUTPUT_VECS))    DEALLOCATE(OUTPUT_VECS)
    IF (ALLOCATED(OUTPUT_SHIFT))   DEALLOCATE(OUTPUT_SHIFT)
    ! Allocate space for all model parameters.
    ALLOCATE( &
         INPUT_VECS(DI,DS), &
         INPUT_SHIFT(DS), &
         INPUT_FLEX(DS), &
         INTERNAL_MEAN(DS, NS), &
         INTERNAL_VECS(DS, DS, NS-1), &
         INTERNAL_SHIFT(DS, NS-1), &
         INTERNAL_FLEX(DS, NS-1), &
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
    !  for all initial parameters in the input, internal, and output.
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
    ! Values are in range [-2,2], assuming no scaling changes,
    !  this will capture 2 standard deviations in either direction.
    INPUT_SHIFT(:) = INITIAL_SHIFT_RANGE * (INPUT_SHIFT(:) * 2.0_RT - 1.0_RT)
    INTERNAL_SHIFT(:,:) = INITIAL_SHIFT_RANGE * (INTERNAL_SHIFT(:,:) * 2.0_RT - 1.0_RT)
    ! Set initial *_FLEX to be equal to default small value.
    INPUT_FLEX(:) = INITIAL_FLEX
    INTERNAL_FLEX(:,:) = INITIAL_FLEX

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
  SUBROUTINE EVALUATE(INPUTS, OUTPUTS)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:,:) :: OUTPUTS
    ! Internal values.
    REAL(KIND=RT), DIMENSION(:,:), ALLOCATABLE :: TEMP_VALUES, &
         TEMP_OUTPUT, INTERNAL_VALUES
    INTEGER :: I, D, L, N
    ! Make sure that this model has been initialized.
    IF (.NOT. ALLOCATED(INPUT_VECS)) THEN
       OUTPUTS(:,:) = 0.0
       RETURN
    END IF
    ! Get the number of points.
    N = SIZE(INPUTS,2)
    ! Allocate storage for the internal values.
    ALLOCATE(INTERNAL_VALUES(1:N,1:MDS), TEMP_VALUES(1:N,1:MDS), &
         TEMP_OUTPUT(1:N,1:MDO) )
    ! Compute the input transformation.
    CALL SGEMM('T', 'N', N, MDS, MDI, 1.0_RT, &
         INPUTS(:,:), SIZE(INPUTS,1), &
         INPUT_VECS(:,:), SIZE(INPUT_VECS,1), &
         0.0_RT, INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1))
    DO D = 1, MDS
       INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) + INPUT_SHIFT(D)
       WHERE (INTERNAL_VALUES(:,D) .LT. DISCONTINUITY)
          INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) * INPUT_FLEX(D)
       END WHERE
    END DO
    ! Compute the next set of internal values with a rectified activation.
    DO L = 1, MNS-1
       ! Make all internal components have a mean of zero.
       DO D = 1, MDS
          INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) - INTERNAL_MEAN(D,L)
       END DO
       ! Compute all vectors.
       CALL SGEMM('N', 'N', N, MDS, MDS, 1.0_RT, &
            INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
            INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
            0.0_RT, TEMP_VALUES(:,:), SIZE(TEMP_VALUES,1))
       INTERNAL_VALUES(:,:) = TEMP_VALUES(:,:)
       DO D = 1, MDS
          INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) + INTERNAL_SHIFT(D,L)
          WHERE (INTERNAL_VALUES(:,D) .LT. DISCONTINUITY)
             INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) * INTERNAL_FLEX(D,L)
          END WHERE
       END DO
    END DO
    ! Make all internal components have a mean of zero.
    L = MNS
    DO D = 1, MDS
       INTERNAL_VALUES(:,D) = INTERNAL_VALUES(:,D) - INTERNAL_MEAN(D,L)
    END DO
    ! Return the final embedded layer if the outputs have the size of the embedding.
    IF ((MDS .NE. MDO) .AND. (SIZE(OUTPUTS,1) .EQ. MDS)) THEN
       ! Do the necessary transpose operation one dimension at a time, produce embedding.
       DO D = 1, MDS
          OUTPUTS(D,1:N) = INTERNAL_VALUES(1:N,D)
       END DO
    ! Otherwise, assume regular output computation.
    ELSE
       ! CALL SGEMM('T', 'T', MDO, N, MDS, 1.0_RT, &
       !      OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
       !      INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
       !      0.0_RT, OUTPUTS(:,:), SIZE(OUTPUTS,1))
       ! DO D = 1, MDO
       !    OUTPUTS(D,1:N) = OUTPUT_SHIFT(D) + OUTPUTS(D,1:N)
       ! END DO
       CALL SGEMM('N', 'N', N, MDO, MDS, 1.0_RT, &
            INTERNAL_VALUES(:,:), SIZE(INTERNAL_VALUES,1), &
            OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
            0.0_RT, TEMP_OUTPUT(:,:), SIZE(TEMP_OUTPUT,1))
       DO D = 1, MDO
          OUTPUTS(D,:) = OUTPUT_SHIFT(D) + TEMP_OUTPUT(:,D)
       END DO
    END IF
    ! Deallocate temporary variables.
    DEALLOCATE(INTERNAL_VALUES, TEMP_VALUES, TEMP_OUTPUT)
  END SUBROUTINE EVALUATE


  ! Compute the gradient of the sum of squared error of this regression
  ! model with respect to its parameters given input and output pairs.
  SUBROUTINE SSE_GRADIENT(INPUTS, OUTPUTS, SUM_SQUARED_ERROR, &
       INPUT_VECS_GRADIENT, INPUT_SHIFT_GRADIENT, INPUT_FLEX_GRADIENT,&
       INTERNAL_VECS_GRADIENT, INTERNAL_SHIFT_GRADIENT, INTERNAL_FLEX_GRADIENT,&
       OUTPUT_VECS_GRADIENT, OUTPUT_SHIFT_GRADIENT)
    ! Data values stored with contiguous points: shape = (MDI,N)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: OUTPUTS
    ! Sum (over all data) squared error (summed over dimensions).
    REAL(KIND=RT), INTENT(INOUT) :: SUM_SQUARED_ERROR
    ! Model parameter gradients.
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:)   :: INPUT_VECS_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:)     :: INPUT_SHIFT_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:)     :: INPUT_FLEX_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:,:) :: INTERNAL_VECS_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:)   :: INTERNAL_SHIFT_GRADIENT
    REAL(KIND=RT), INTENT(INOUT), DIMENSION(:,:)   :: INTERNAL_FLEX_GRADIENT
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
         INTERNAL_VALUES(:,D,1) = INTERNAL_VALUES(:,D,1) + INPUT_SHIFT(D)
         LT_DISCONTINUITY(:,D,1) = INTERNAL_VALUES(:,D,1) .LT. DISCONTINUITY
         WHERE (LT_DISCONTINUITY(:,D,1))
            INTERNAL_VALUES(:,D,1) = INTERNAL_VALUES(:,D,1) * INPUT_FLEX(D)
         END WHERE
      END DO
      ! Compute the next set of internal values with a rectified activation.
      DO L = 1, MNS-1
         ! Update the MEAN for this layer.
         DO D = 1, MDS
            ! INTERNAL_MEAN(D,L) = SUM(INTERNAL_VALUES(:,D,L)) / N
            INTERNAL_MEAN(D,L) = (MAXVAL(INTERNAL_VALUES(:,D,L)) + MINVAL(INTERNAL_VALUES(:,D,L))) / 2
            ! INTERNAL_MEAN(D,L) = 0.0_RT
            INTERNAL_VALUES(:,D,L) = INTERNAL_VALUES(:,D,L) - INTERNAL_MEAN(D,L)
         END DO
         ! Compute vector values.
         LP1 = L+1
         CALL SGEMM('N', 'N', SIZE(INPUTS,2), MDS, MDS, 1.0_RT, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
              0.0_RT, INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1))
         DO D = 1, MDS
            INTERNAL_VALUES(:,D,LP1) = INTERNAL_VALUES(:,D,LP1) + INTERNAL_SHIFT(D,L)
            LT_DISCONTINUITY(:,D,LP1) = INTERNAL_VALUES(:,D,LP1) .LT. DISCONTINUITY
            WHERE (LT_DISCONTINUITY(:,D,LP1))
               INTERNAL_VALUES(:,D,LP1) = INTERNAL_VALUES(:,D,LP1) * INTERNAL_FLEX(D,L)
            END WHERE
         END DO
      END DO
      ! Update the MEAN for the output layer.
      DO D = 1, MDS
         ! INTERNAL_MEAN(D,L) = SUM(INTERNAL_VALUES(:,D,L)) / N
         INTERNAL_MEAN(D,L) = (MAXVAL(INTERNAL_VALUES(:,D,L)) + MINVAL(INTERNAL_VALUES(:,D,L))) / 2
         ! INTERNAL_MEAN(D,L) = 0.0_RT
         INTERNAL_VALUES(:,D,L) = INTERNAL_VALUES(:,D,L) - INTERNAL_MEAN(D,L)
      END DO
      ! Compute the output.
      CALL SGEMM('N', 'N', SIZE(INPUTS,2), MDO, MDS, 1.0_RT, &
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
      ! Compute output parameter gradients.
      CALL SGEMM('T', 'N', MDS, MDO, SIZE(INPUTS,2), 1.0_RT / N, &
           INTERNAL_VALUES(:,:,MNS), SIZE(INTERNAL_VALUES,1), &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           1.0_RT, OUTPUT_VECS_GRADIENT(:,:), SIZE(OUTPUT_VECS_GRADIENT,1))
      ! Propogate the gradient back to the last internal vector space.
      CALL SGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDO, 1.0_RT, &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
           0.0_RT, INTERNAL_TEMP(:,:), SIZE(INTERNAL_TEMP,1))
      ! Cycle over all internal layers.
      INTERNAL_REPRESENTATIONS : DO L = MNS-1, 1, -1
         LP1 = L+1
         DO D = 1, MDS
            ! Compute the flexure gradient.
            IF (INTERNAL_FLEX(D,L) .NE. 0.0_RT) THEN
               INTERNAL_FLEX_GRADIENT(D,L) = SUM(&
                    INTERNAL_TEMP(:,D) * INTERNAL_VALUES(:,D,LP1), 1, &
                    LT_DISCONTINUITY(:,D,LP1)) / INTERNAL_FLEX(D,L)
            ELSE
               INTERNAL_FLEX_GRADIENT(D,L) = 0.0_RT
            END IF
            ! Propogate the error gradient back to the preceding vectors.
            WHERE (LT_DISCONTINUITY(:,D,LP1))
               INTERNAL_VALUES(:,D,LP1) = INTERNAL_TEMP(:,D) * INTERNAL_FLEX(D,L)
            ELSEWHERE
               INTERNAL_VALUES(:,D,LP1) = INTERNAL_TEMP(:,D)
            END WHERE
         END DO
         ! Compute the shift gradient.
         INTERNAL_SHIFT_GRADIENT(:,L) = SUM(INTERNAL_VALUES(:,:,LP1), 1) / N &
              + INTERNAL_SHIFT_GRADIENT(:,L)
         ! Compute the gradient with respect to each output and all inputs.
         CALL SGEMM('T', 'N', MDS, MDS, SIZE(INPUTS,2), 1.0_RT / N, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              1.0_RT, INTERNAL_VECS_GRADIENT(:,:,L), SIZE(INTERNAL_VECS_GRADIENT,1))
         ! Propogate the gradient to the immediately preceding layer's flexure.
         CALL SGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDS, 1.0_RT, &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
              0.0_RT, INTERNAL_TEMP(:,:), SIZE(INTERNAL_TEMP,1))
      END DO INTERNAL_REPRESENTATIONS
      ! Compute the gradients going into the first layer.
      DO D = 1, MDS
         ! Compute the flexure gradient.
         IF (INPUT_FLEX(D) .NE. 0.0_RT) THEN
            INPUT_FLEX_GRADIENT(D) = SUM(&
                 INTERNAL_TEMP(:,D) * INTERNAL_VALUES(:,D,1), 1, &
                 LT_DISCONTINUITY(:,D,1)) / INPUT_FLEX(D)
         ELSE
            INPUT_FLEX_GRADIENT(D) = 0.0_RT
         END IF
         ! Propogate the error gradient back to the preceding vectors.
         WHERE (LT_DISCONTINUITY(:,D,1))
            INTERNAL_VALUES(:,D,1) = INTERNAL_TEMP(:,D) * INPUT_FLEX(D)
         ELSEWHERE
            INTERNAL_VALUES(:,D,1) = INTERNAL_TEMP(:,D)
         END WHERE
      END DO
      ! Compute the input parameter gradients.
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
  SUBROUTINE MINIMIZE_MSE(X, Y, STEPS, BATCH_SIZE, NUM_THREADS, &
       STEP_SIZE, KEEP_BEST, MEAN_SQUARED_ERROR, RECORD)
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: X
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: Y
    INTEGER,       INTENT(IN) :: STEPS
    INTEGER,       INTENT(IN), OPTIONAL :: BATCH_SIZE, NUM_THREADS
    REAL(KIND=RT), INTENT(IN), OPTIONAL :: STEP_SIZE
    LOGICAL,       INTENT(IN), OPTIONAL :: KEEP_BEST
    REAL(KIND=RT), INTENT(OUT) :: MEAN_SQUARED_ERROR
    REAL(KIND=RT), INTENT(OUT), DIMENSION(STEPS), OPTIONAL :: RECORD
    !  gradient step arrays, 4 copies of internal model (total of 5).
    REAL(KIND=RT), DIMENSION(MDS,MNS)       :: MEAN_BEST
    REAL(KIND=RT), DIMENSION(MDI,MDS)       :: V1_GRAD, V1_MEAN, V1_VAR, BEST_V1
    REAL(KIND=RT), DIMENSION(MDS)           :: S1_GRAD, S1_MEAN, S1_VAR, BEST_S1
    REAL(KIND=RT), DIMENSION(MDS)           :: F1_GRAD, F1_MEAN, F1_VAR, BEST_F1
    REAL(KIND=RT), DIMENSION(MDS,MDS,MNS-1) :: V2_GRAD, V2_MEAN, V2_VAR, BEST_V2
    REAL(KIND=RT), DIMENSION(MDS,MNS-1)     :: S2_GRAD, S2_MEAN, S2_VAR, BEST_S2
    REAL(KIND=RT), DIMENSION(MDS,MNS-1)     :: F2_GRAD, F2_MEAN, F2_VAR, BEST_F2
    REAL(KIND=RT), DIMENSION(MDS,MDO)       :: V3_GRAD, V3_MEAN, V3_VAR, BEST_V3
    REAL(KIND=RT), DIMENSION(MDO)           :: S3_GRAD, S3_MEAN, S3_VAR, BEST_S3
    REAL(KIND=RT), DIMENSION(SIZE(Y,1),SIZE(Y,2)) :: SAMPLE_OUTPUTS
    !  local variables
    LOGICAL :: REVERT_TO_BEST
    INTEGER :: BS, I, J, L, NX, NT, BATCH_START, BATCH_END
    REAL(KIND=RT) :: RNX, BATCHES, PREV_MSE, BEST_MSE, SCALAR
    REAL(KIND=RT) :: STEP_FACTOR, STEP_MEAN_CHANGE, STEP_MEAN_REMAIN, &
         STEP_VAR_CHANGE, STEP_VAR_REMAIN
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
    STEP_VAR_CHANGE = INITIAL_STEP_VAR_CHANGE
    STEP_VAR_REMAIN = 1.0_RT - STEP_VAR_CHANGE
    ! Initial mean squared error is "max representable value".
    PREV_MSE = HUGE(PREV_MSE)
    BEST_MSE = HUGE(BEST_MSE)
    ! Set the best shift and scale (for internal consistency) to initial values.
    MEAN_BEST(:,:) = 0.0_RT
    ! Set the average step sizes.
    V1_MEAN(:,:)   = 0.0_RT
    S1_MEAN(:)     = 0.0_RT
    F1_MEAN(:)     = 0.0_RT
    V2_MEAN(:,:,:) = 0.0_RT
    S2_MEAN(:,:)   = 0.0_RT
    F2_MEAN(:,:)   = 0.0_RT
    V3_MEAN(:,:)   = 0.0_RT
    S3_MEAN(:)     = 0.0_RT
    ! Set the average step size variances.
    V1_VAR(:,:)   = 0.0_RT
    S1_VAR(:)     = 0.0_RT
    F1_VAR(:)     = 0.0_RT
    V2_VAR(:,:,:) = 0.0_RT
    S2_VAR(:,:)   = 0.0_RT
    F2_VAR(:,:)   = 0.0_RT
    V3_VAR(:,:)   = 0.0_RT
    S3_VAR(:)     = 0.0_RT
    ! Iterate, taking steps with the average gradient over all data.
    fit_loop : DO I = 1, STEPS
       ! Compute the average gradient over all points.
       MEAN_SQUARED_ERROR = 0.0_RT
       ! Set gradients to zero initially.
       V1_GRAD(:,:)   = 0.0_RT
       S1_GRAD(:)     = 0.0_RT
       F1_GRAD(:)     = 0.0_RT
       V2_GRAD(:,:,:) = 0.0_RT
       S2_GRAD(:,:)   = 0.0_RT
       F2_GRAD(:,:)   = 0.0_RT
       V3_GRAD(:,:)   = 0.0_RT
       S3_GRAD(:)     = 0.0_RT
       ! Count the number of batches.
       BATCHES = 0.0_RT
       !$OMP PARALLEL DO NUM_THREADS(NT) PRIVATE(BATCH_END) &
       !$OMP&  REDUCTION(+: BATCHES, MEAN_SQUARED_ERROR, &
       !$OMP&  V1_GRAD, S1_GRAD, F1_GRAD, V2_GRAD, S2_GRAD, F2_GRAD, &
       !$OMP&  V3_GRAD, S3_GRAD)
       DO BATCH_START = 1, NX, BS
          BATCHES = BATCHES + 1.0_RT
          BATCH_END = MIN(NX, BATCH_START+BS-1)
          ! Sum the gradient over all data batches.
          CALL SSE_GRADIENT(X(:,BATCH_START:BATCH_END), &
               Y(:,BATCH_START:BATCH_END), MEAN_SQUARED_ERROR, &
               V1_GRAD(:,:), S1_GRAD(:), F1_GRAD(:), &
               V2_GRAD(:,:,:), S2_GRAD(:,:), F2_GRAD(:,:), &
               V3_GRAD(:,:), S3_GRAD(:))
       END DO
       !$OMP END PARALLEL DO
       ! Convert the sum of squared errors into the mean squared error.
       MEAN_SQUARED_ERROR = MEAN_SQUARED_ERROR / RNX

       ! Evaluate the network to see what its outputs should look like.
       CALL EVALUATE(X(:,:), SAMPLE_OUTPUTS(:,:))       
       BATCHES = SUM((SAMPLE_OUTPUTS(:,:) - Y(:,:))**2) / RNX
       IF ( (ABS(MEAN_SQUARED_ERROR - BATCHES) / &
            (1.0_RT + MEAN_SQUARED_ERROR)) .GT. &
            SQRT(SQRT(EPSILON(1.0_RT))) ) THEN
          PRINT *, 'Mean squared errors do not aggree between the SSE'
          PRINT *, '  computation and actual forward evaluation of network.'
          PRINT *, ' ', MEAN_SQUARED_ERROR
          PRINT *, ' ', SAMPLE_OUTPUTS(1,1)
          MEAN_SQUARED_ERROR = -1.0_RT
          RETURN
       END IF

       ! Update the saved "best" model based on error (only when dropout is disabled).
       IF (MEAN_SQUARED_ERROR .LT. BEST_MSE) THEN
          BEST_MSE       = MEAN_SQUARED_ERROR
          MEAN_BEST(:,:) = INTERNAL_MEAN(:,:)
          BEST_V1(:,:)   = INPUT_VECS(:,:)
          BEST_S1(:)     = INPUT_SHIFT(:)
          BEST_F1(:)     = INPUT_FLEX(:)
          BEST_V2(:,:,:) = INTERNAL_VECS(:,:,:)
          BEST_S2(:,:)   = INTERNAL_SHIFT(:,:)
          BEST_F2(:,:)   = INTERNAL_SHIFT(:,:)
          BEST_V3(:,:)   = OUTPUT_VECS(:,:)
          BEST_S3(:)     = OUTPUT_SHIFT(:)
       END IF
       ! Update the step factor based on model improvement.
       IF (MEAN_SQUARED_ERROR .LT. PREV_MSE) THEN
          STEP_FACTOR = STEP_FACTOR * STEP_GROWTH_RATE
          STEP_MEAN_CHANGE = STEP_MEAN_CHANGE * STEP_GROWTH_RATE
          STEP_MEAN_REMAIN = 1.0_RT - STEP_MEAN_CHANGE
          STEP_VAR_CHANGE = STEP_VAR_CHANGE * STEP_GROWTH_RATE
          STEP_VAR_REMAIN = 1.0_RT - STEP_VAR_CHANGE
       ELSE
          STEP_FACTOR = STEP_FACTOR * STEP_SHRINK_RATE
          STEP_MEAN_CHANGE = STEP_MEAN_CHANGE * STEP_SHRINK_RATE
          STEP_MEAN_REMAIN = 1.0_RT - STEP_MEAN_CHANGE
          STEP_VAR_CHANGE = STEP_VAR_CHANGE * STEP_SHRINK_RATE
          STEP_VAR_REMAIN = 1.0_RT - STEP_VAR_CHANGE
       END IF
       ! Store this error if a record is being kept.
       IF (PRESENT(RECORD)) RECORD(I) = MEAN_SQUARED_ERROR
       ! Store the previous error for tracking the best-so-far.
       PREV_MSE = MEAN_SQUARED_ERROR
       ! Convert the summed gradients to average gradients.
       V1_GRAD(:,:)   = V1_GRAD(:,:)   / BATCHES
       S1_GRAD(:)     = S1_GRAD(:)     / BATCHES
       F1_GRAD(:)     = F1_GRAD(:)     / BATCHES
       V2_GRAD(:,:,:) = V2_GRAD(:,:,:) / BATCHES
       S2_GRAD(:,:)   = S2_GRAD(:,:)   / BATCHES
       F2_GRAD(:,:)   = F2_GRAD(:,:)   / BATCHES
       V3_GRAD(:,:)   = V3_GRAD(:,:)   / BATCHES
       S3_GRAD(:)     = S3_GRAD(:)     / BATCHES
       ! Update the steps for all different parameters.
       CALL COMPUTE_STEP(SIZE(V1_GRAD(:,:)),   V1_GRAD(:,:),   V1_MEAN(:,:),   V1_VAR(:,:))
       CALL COMPUTE_STEP(SIZE(S1_GRAD(:)),     S1_GRAD(:),     S1_MEAN(:),     S1_VAR(:))
       CALL COMPUTE_STEP(SIZE(F1_GRAD(:)),     F1_GRAD(:),     F1_MEAN(:),     F1_VAR(:))
       CALL COMPUTE_STEP(SIZE(V2_GRAD(:,:,:)), V2_GRAD(:,:,:), V2_MEAN(:,:,:), V2_VAR(:,:,:))
       CALL COMPUTE_STEP(SIZE(S2_GRAD(:,:)),   S2_GRAD(:,:),   S2_MEAN(:,:),   S2_VAR(:,:))
       CALL COMPUTE_STEP(SIZE(F2_GRAD(:,:)),   F2_GRAD(:,:),   F2_MEAN(:,:),   F2_VAR(:,:))
       CALL COMPUTE_STEP(SIZE(V3_GRAD(:,:)),   V3_GRAD(:,:),   V3_MEAN(:,:),   V3_VAR(:,:))
       CALL COMPUTE_STEP(SIZE(S3_GRAD(:)),     S3_GRAD(:),     S3_MEAN(:),     S3_VAR(:))
       ! Take the gradient steps (based on the computed "step" above).
       INPUT_VECS(:,:)      = DECAY_FACTOR * (INPUT_VECS(:,:)      - V1_GRAD(:,:))
       INPUT_SHIFT(:)       = DECAY_FACTOR * (INPUT_SHIFT(:)       - S1_GRAD(:))
       INPUT_FLEX(:)        = DECAY_FACTOR * (INPUT_FLEX(:)        - F1_GRAD(:))
       INTERNAL_VECS(:,:,:) = DECAY_FACTOR * (INTERNAL_VECS(:,:,:) - V2_GRAD(:,:,:))
       INTERNAL_SHIFT(:,:)  = DECAY_FACTOR * (INTERNAL_SHIFT(:,:)  - S2_GRAD(:,:))
       INTERNAL_FLEX(:,:)   = DECAY_FACTOR * (INTERNAL_FLEX(:,:)   - F2_GRAD(:,:))
       OUTPUT_VECS(:,:)     = DECAY_FACTOR * (OUTPUT_VECS(:,:)     - V3_GRAD(:,:))
       OUTPUT_SHIFT(:)      = DECAY_FACTOR * (OUTPUT_SHIFT(:)      - S3_GRAD(:))
       ! Maintain a constant max-norm across the magnitue of internal vectors.
       DO L = 1, MNS-1
          SCALAR = SQRT(MAXVAL(SUM(INTERNAL_VECS(:,:,L)**2, 1)))
          INTERNAL_VECS(:,:,L) = INTERNAL_VECS(:,:,L) / SCALAR
          SCALAR = MAXVAL(ABS(INTERNAL_FLEX(:,L)))
          IF (SCALAR .GT. 1) THEN
             INTERNAL_FLEX(:,L) = INTERNAL_FLEX(:,L) / SCALAR
          END IF
       END DO

    END DO fit_loop

    ! Restore the best model seen so far (if enough steps were taken).
    IF (REVERT_TO_BEST) THEN
       MEAN_SQUARED_ERROR   = BEST_MSE
       INTERNAL_MEAN(:,:)   = MEAN_BEST(:,:)
       INPUT_VECS(:,:)      = BEST_V1(:,:)  
       INPUT_SHIFT(:)       = BEST_S1(:)    
       INPUT_FLEX(:)        = BEST_F1(:)    
       INTERNAL_VECS(:,:,:) = BEST_V2(:,:,:)
       INTERNAL_SHIFT(:,:)  = BEST_S2(:,:)  
       INTERNAL_FLEX(:,:)   = BEST_F2(:,:)  
       OUTPUT_VECS(:,:)     = BEST_V3(:,:)  
       OUTPUT_SHIFT(:)      = BEST_S3(:)    
    END IF

  CONTAINS

    ! Compute the current step factor as the running average step,
    ! inversely scaled by the average magnitude of the step.
    SUBROUTINE COMPUTE_STEP(N, CURR_STEP, STEP_MEAN, STEP_VAR)
      INTEGER, INTENT(IN) :: N
      REAL(KIND=RT), DIMENSION(1:N), INTENT(INOUT) :: CURR_STEP
      REAL(KIND=RT), DIMENSION(1:N), INTENT(INOUT) :: STEP_MEAN
      REAL(KIND=RT), DIMENSION(1:N), INTENT(INOUT) :: STEP_VAR
      ! Update sliding window mean and variance calculations.
      STEP_MEAN(:) = STEP_MEAN_CHANGE * CURR_STEP(:) &
           + STEP_MEAN_REMAIN * STEP_MEAN(:)
      STEP_VAR(:)  = STEP_VAR_CHANGE  * (STEP_MEAN(:) - CURR_STEP(:))**2 &
           + STEP_VAR_REMAIN  * STEP_VAR(:)
      STEP_VAR(:) = MAX(STEP_VAR(:), EPSILON(STEP_FACTOR))
      ! Take averaged gradient descent step.
      CURR_STEP(:) = STEP_FACTOR * STEP_MEAN(:)
      ! Start scaling by step magnitude by variance once enough data is collected.
      IF (I .GT. MIN_STEPS_TO_STABILITY) THEN
         CURR_STEP(:) = CURR_STEP(:) / SQRT(STEP_VAR(:))
      END IF
    END SUBROUTINE COMPUTE_STEP

  END SUBROUTINE MINIMIZE_MSE

END MODULE PLRM


! Use this module to visualize the values and error gradients within
!  the PLRM model for different parameter types.
MODULE VIS_PLRM

  ! Pieces of PLRM that are used here.
  USE ISO_FORTRAN_ENV, ONLY: RT => REAL32
  USE PLRM, ONLY: INPUT_VECS, INPUT_SHIFT, INPUT_FLEX, &
       INTERNAL_MEAN, INTERNAL_VECS, INTERNAL_SHIFT, INTERNAL_FLEX, &
       OUTPUT_VECS, OUTPUT_SHIFT, MDI, MDS, MDO, MNS, DISCONTINUITY, &
       EVALUATE
  IMPLICIT NONE

CONTAINS

  ! Disable the 'layer, position' node (1 indexed) in the network and
  !  compute the mean squared error "MSE". Use layer=-1 to compute
  !  MSE without disabling any of the internal representations.
  RECURSIVE SUBROUTINE DISABLE_AND_COMPUTE_MSE(LAYER, POSITION, INPUTS, OUTPUTS, MSE)
    INTEGER, INTENT(IN) :: LAYER, POSITION
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: INPUTS
    REAL(KIND=RT), INTENT(IN), DIMENSION(:,:) :: OUTPUTS
    REAL(KIND=RT), INTENT(OUT) :: MSE
    ! Internal values.
    REAL(KIND=RT), DIMENSION(:,:), ALLOCATABLE :: RESULTS
    REAL(KIND=RT) :: FULL_MSE
    REAL(KIND=RT) :: TEMP_SHIFT, TEMP_FLEX, N
    ! Allocate storage for the internal values.
    ALLOCATE(RESULTS(1:SIZE(OUTPUTS,1), 1:SIZE(OUTPUTS,2)))
    ! Store the flex and shift for that node temporarily.
    IF (LAYER .EQ. 1) THEN
       TEMP_SHIFT = INPUT_SHIFT(POSITION)
       TEMP_FLEX = INPUT_FLEX(POSITION)
       INPUT_SHIFT(POSITION) = -SQRT(HUGE(1.0_RT))
       INPUT_FLEX(POSITION) = 0.0_RT
    ELSE IF (LAYER .GT. 1) THEN
       TEMP_SHIFT = INTERNAL_SHIFT(POSITION,LAYER-1)
       TEMP_FLEX = INTERNAL_FLEX(POSITION,LAYER-1)
       INTERNAL_SHIFT(POSITION,LAYER-1) = -SQRT(HUGE(1.0_RT))
       INTERNAL_FLEX(POSITION,LAYER-1) = 0.0_RT
    END IF
    ! Evaluate the network (with the appropriate node disabled).
    CALL EVALUATE(INPUTS(:,:), RESULTS(:,:))
    ! Compute the diference between the results and the outputs.
    N = REAL(SIZE(INPUTS,2), RT)
    MSE = SUM((RESULTS(:,:) - OUTPUTS(:,:))**2) / N
    ! Deallocate temporary variables.
    DEALLOCATE(RESULTS)
    ! Restore the flex and shift to the original values.
    IF (LAYER .EQ. 1) THEN
       INPUT_SHIFT(POSITION) = TEMP_SHIFT
       INPUT_FLEX(POSITION) = TEMP_FLEX
    ELSE IF (LAYER .GT. 1) THEN
       INTERNAL_SHIFT(POSITION,LAYER-1) = TEMP_SHIFT
       INTERNAL_FLEX(POSITION,LAYER-1) = TEMP_FLEX
    END IF
  END SUBROUTINE DISABLE_AND_COMPUTE_MSE

  ! Given index, compute value and gradient at internal position in model.
  SUBROUTINE COMPUTE_VALUES(LAYER, POSITION, INPUTS, OUTPUTS, VALS, GRADS)
    ! Arguments
    INTEGER, INTENT(IN) :: LAYER, POSITION
    REAL(KIND=RT), INTENT(IN),  DIMENSION(:,:) :: INPUTS, OUTPUTS
    REAL(KIND=RT), INTENT(OUT), DIMENSION(:) :: VALS, GRADS
    ! Internal values and lists of which ones are nonzero.
    REAL(KIND=RT), DIMENSION(SIZE(INPUTS,2),MDO)     :: OUTPUT_GRADIENT
    REAL(KIND=RT), DIMENSION(SIZE(INPUTS,2),MDS,MNS) :: INTERNAL_VALUES
    LOGICAL,       DIMENSION(SIZE(INPUTS,2),MDS,MNS) :: LT_DISCONTINUITY
    REAL(KIND=RT), DIMENSION(SIZE(INPUTS,2),MDS)     :: INTERNAL_TEMP
    ! Number of points.
    REAL(KIND=RT) :: N
    INTEGER :: D
    ! Routine for matrix multiplication.
    EXTERNAL :: SGEMM
    ! Evaluate the model at all data points, store outputs in "OUTPUT_GRADIENT"
    CALL EVALUATE_BATCH()
    ! Get the internal values at all points at the specific node.
    VALS(:) = INTERNAL_VALUES(:,POSITION,LAYER)
    PRINT *, '      ', SUM(OUTPUT_GRADIENT(:,:)**2)
    ! Compute the new gradient.
    DO D = 1, MDO
       OUTPUT_GRADIENT(:,D) = OUTPUT_GRADIENT(:,D) - OUTPUTS(D,:)
    END DO
    ! Get the number of points.
    N = REAL(SIZE(INPUTS,2), RT)
    PRINT *, 'SS OG:', SUM(OUTPUT_GRADIENT(:,:)**2)
    PRINT *, '      ', SUM(OUTPUTS(:,:)**2)
    ! Compute the gradient of parameters with respect to the output gradient.
    CALL GRADIENT_BATCH()
    ! Get the internal gradients at the nodes.
    GRADS(:) = INTERNAL_VALUES(:,POSITION,LAYER)

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
         INTERNAL_VALUES(:,D,1) = INTERNAL_VALUES(:,D,1) + INPUT_SHIFT(D)
         LT_DISCONTINUITY(:,D,1) = INTERNAL_VALUES(:,D,1) .LT. DISCONTINUITY
         WHERE (LT_DISCONTINUITY(:,D,1))
            INTERNAL_VALUES(:,D,1) = INTERNAL_VALUES(:,D,1) * INPUT_FLEX(D)
         END WHERE
      END DO
      ! Compute the next set of internal values with a rectified activation.
      DO L = 1, MNS-1
         ! Update the MEAN for this layer.
         DO D = 1, MDS
            INTERNAL_VALUES(:,D,L) = INTERNAL_VALUES(:,D,L) - INTERNAL_MEAN(D,L)
         END DO

         LP1 = L+1
         CALL SGEMM('N', 'N', SIZE(INPUTS,2), MDS, MDS, 1.0_RT, &
              INTERNAL_VALUES(:,:,L), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
              0.0_RT, INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1))
         DO D = 1, MDS
            INTERNAL_VALUES(:,D,LP1) = INTERNAL_VALUES(:,D,LP1) + INTERNAL_SHIFT(D,L)
            LT_DISCONTINUITY(:,D,LP1) = INTERNAL_VALUES(:,D,LP1) .LT. DISCONTINUITY
            WHERE (LT_DISCONTINUITY(:,D,LP1))
               INTERNAL_VALUES(:,D,LP1) = INTERNAL_VALUES(:,D,LP1) * INTERNAL_FLEX(D,L)
            END WHERE
         END DO
      END DO
      ! Update the MEAN for the output layer.
      DO D = 1, MDS
         INTERNAL_VALUES(:,D,MNS) = INTERNAL_VALUES(:,D,MNS) - INTERNAL_MEAN(D,L)
      END DO
      ! Compute the output.
      CALL SGEMM('N', 'N', SIZE(INPUTS,2), MDO, MDS, 1.0_RT, &
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
      ! Propogate the gradient back to the last internal vector space.
      CALL SGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDO, 1.0_RT, &
           OUTPUT_GRADIENT(:,:), SIZE(OUTPUT_GRADIENT,1), &
           OUTPUT_VECS(:,:), SIZE(OUTPUT_VECS,1), &
           0.0_RT, INTERNAL_TEMP(:,:), SIZE(INTERNAL_TEMP,1))
      ! Cycle over all internal layers.
      INTERNAL_REPRESENTATIONS : DO L = MNS-1, 1, -1
         LP1 = L+1
         DO D = 1, MDS
            ! Propogate the error gradient back to the preceding vectors.
            WHERE (LT_DISCONTINUITY(:,D,LP1))
               INTERNAL_VALUES(:,D,LP1) = INTERNAL_TEMP(:,D) * INTERNAL_FLEX(D,L)
            ELSEWHERE
               INTERNAL_VALUES(:,D,LP1) = INTERNAL_TEMP(:,D)
            END WHERE
         END DO
         ! Propogate the gradient to the immediately preceding layer's flexure.
         CALL SGEMM('N', 'T', SIZE(INPUTS,2), MDS, MDS, 1.0_RT, &
              INTERNAL_VALUES(:,:,LP1), SIZE(INTERNAL_VALUES,1), &
              INTERNAL_VECS(:,:,L), SIZE(INTERNAL_VECS,1), &
              0.0_RT, INTERNAL_TEMP(:,:), SIZE(INTERNAL_TEMP,1))
      END DO INTERNAL_REPRESENTATIONS
      ! Compute the gradients going into the first layer.
      DO D = 1, MDS
         ! Propogate the error gradient back to the preceding vectors.
         WHERE (LT_DISCONTINUITY(:,D,1))
            INTERNAL_VALUES(:,D,1) = INTERNAL_TEMP(:,D) * INPUT_FLEX(D)
         ELSEWHERE
            INTERNAL_VALUES(:,D,1) = INTERNAL_TEMP(:,D)
         END WHERE
      END DO
    END SUBROUTINE GRADIENT_BATCH
    
  END SUBROUTINE COMPUTE_VALUES


END MODULE VIS_PLRM


! These go in the SSE_GRADIENT function.
! ----------------------------------------------------------------
! REAL(KIND=RT), DIMENSION(MDS,MDS,MNS) :: PAIRWISE_DIFFS
! ----------------------------------------------------------------
! ----------------------------------------------------------------
! ! Compute pairwise differences between the outputs of the internal vectors.
! PAIRWISE_DIFFS(:,:,:) = 0.0_RT
! DO L = 1, MNS
!    DO D = 1, MDS
!       IF (D .GT. 1) &
!            CONTINUE
!       IF (D .LT. MDS) &
!            CONTINUE
!    END DO
! END DO
! ----------------------------------------------------------------
