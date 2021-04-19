! A piecewise linear regression model.
MODULE PLRM
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE
  ! Value used for rectification, default is 0.
  REAL(KIND=REAL64) :: RECTIFY_AT
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
  ! Model parameter gradients.
  REAL(KIND=REAL64), DIMENSION(:,:),   ALLOCATABLE :: INPUT_PARAMS_GRADIENT
  REAL(KIND=REAL64), DIMENSION(:),     ALLOCATABLE :: INPUT_BIAS_GRADIENT
  REAL(KIND=REAL64), DIMENSION(:,:,:), ALLOCATABLE :: INTERNAL_PARAMS_GRADIENT
  REAL(KIND=REAL64), DIMENSION(:,:),   ALLOCATABLE :: INTERNAL_BIAS_GRADIENT
  REAL(KIND=REAL64), DIMENSION(:,:),   ALLOCATABLE :: OUTPUT_PARAMS_GRADIENT
  REAL(KIND=REAL64), DIMENSION(:),     ALLOCATABLE :: OUTPUT_BIAS_GRADIENT     
  ! Internal values needed for training.
  REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: INTERNAL_VALUES
  ! Set all arrays related to gradient computaiton as private to Open MP threads.
  !$OMP THREADPRIVATE(INPUT_PARAMS_GRADIENT    )
  !$OMP THREADPRIVATE(INPUT_BIAS_GRADIENT      )
  !$OMP THREADPRIVATE(INTERNAL_PARAMS_GRADIENT )
  !$OMP THREADPRIVATE(INTERNAL_BIAS_GRADIENT   )
  !$OMP THREADPRIVATE(OUTPUT_PARAMS_GRADIENT   )
  !$OMP THREADPRIVATE(OUTPUT_BIAS_GRADIENT     )
  !$OMP THREADPRIVATE(INTERNAL_VALUES)
  
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

  ! Initialize a new model.
  SUBROUTINE NEW_MODEL(DI, DS, NS, DO, SEED, RECTIFIER, BIAS)
    INTEGER, INTENT(IN) :: DI, DS, NS, DO
    INTEGER, INTENT(IN), OPTIONAL :: SEED
    REAL(KIND=REAL64), INTENT(IN), OPTIONAL :: RECTIFIER
    REAL(KIND=REAL64), INTENT(IN), OPTIONAL :: BIAS
    ! Internal parameters.
    INTEGER :: I
    REAL(KIND=REAL64) :: BIAS_VALUE
    !  Storage for seeding the random number generator (for repeatability).
    INTEGER, DIMENSION(:), ALLOCATABLE :: SEED_ARRAY
    ! Set a random seed, if one was provided (otherwise leave default).
    IF (PRESENT(SEED)) THEN
       CALL RANDOM_SEED(SIZE=I)
       ALLOCATE(SEED_ARRAY(I))
       SEED_ARRAY(:) = SEED
       CALL RANDOM_SEED(PUT=SEED_ARRAY(:))
    END IF
    ! Set the default rectifier value.
    IF (PRESENT(RECTIFIER)) THEN
       RECTIFY_AT = RECTIFIER
    ELSE
       RECTIFY_AT = 0.0_REAL64
    END IF
    ! Set the default bias value.
    IF (PRESENT(BIAS)) THEN
       BIAS_VALUE = BIAS
    ELSE
       BIAS_VALUE = 2.0_REAL64 + RECTIFY_AT
    END IF
    ! Store the integer sizes related to this model.
    MDI = DI
    MDS = DS
    MNS = NS
    MDO = DO
    ! Deallocate any of the internal arrays if they have already been allocated.
    !   parameters
    IF (ALLOCATED(INPUT_PARAMS))    DEALLOCATE(INPUT_PARAMS)
    IF (ALLOCATED(INPUT_BIAS))      DEALLOCATE(INPUT_BIAS)
    IF (ALLOCATED(INTERNAL_PARAMS)) DEALLOCATE(INTERNAL_PARAMS)
    IF (ALLOCATED(INTERNAL_BIAS))   DEALLOCATE(INTERNAL_BIAS)
    IF (ALLOCATED(OUTPUT_PARAMS))   DEALLOCATE(OUTPUT_PARAMS)
    IF (ALLOCATED(OUTPUT_BIAS))     DEALLOCATE(OUTPUT_BIAS)
    ALLOCATE( &
         ! Parameters.
         INPUT_PARAMS(DI,DS), &
         INPUT_BIAS(DS), &
         INTERNAL_PARAMS(DS, DS, NS-1), &
         INTERNAL_BIAS(DS, NS-1), &
         OUTPUT_PARAMS(DS, DO), &
         OUTPUT_BIAS(DO) &
         )
    ! Do all (de)allocations for gradients and internal values in
    !  a "parallel" region so that each thread correctly allocates.
    !$OMP PARALLEL
    !   gradients
    IF (ALLOCATED(INPUT_PARAMS_GRADIENT))    DEALLOCATE(INPUT_PARAMS_GRADIENT)
    IF (ALLOCATED(INPUT_BIAS_GRADIENT))      DEALLOCATE(INPUT_BIAS_GRADIENT)
    IF (ALLOCATED(INTERNAL_PARAMS_GRADIENT)) DEALLOCATE(INTERNAL_PARAMS_GRADIENT)
    IF (ALLOCATED(INTERNAL_BIAS_GRADIENT))   DEALLOCATE(INTERNAL_BIAS_GRADIENT)
    IF (ALLOCATED(OUTPUT_PARAMS_GRADIENT))   DEALLOCATE(OUTPUT_PARAMS_GRADIENT)
    IF (ALLOCATED(OUTPUT_BIAS_GRADIENT))     DEALLOCATE(OUTPUT_BIAS_GRADIENT)
    !   internal state holder
    IF (ALLOCATED(INTERNAL_VALUES)) DEALLOCATE(INTERNAL_VALUES)
    ! Allocate new storage.
    ALLOCATE( &
         ! Gradients.
         INPUT_PARAMS_GRADIENT(DI,DS), &
         INPUT_BIAS_GRADIENT(DS), &
         INTERNAL_PARAMS_GRADIENT(DS, DS, NS-1), &
         INTERNAL_BIAS_GRADIENT(DS, NS-1), &
         OUTPUT_PARAMS_GRADIENT(DS, DO), &
         OUTPUT_BIAS_GRADIENT(DO), &
         ! Internal values.
         INTERNAL_VALUES(DS, NS) &
         )
    !$OMP END PARALLEL

    ! Generate random unit vectors with bias values proportional to
    ! surface area of the sphere at different radii. This should create
    ! a roughly well-spaced set of intersecting simplices formed by
    ! the ReLU activaiton functions over the input data.
    CALL RANDOM_UNIT_VECTORS(INPUT_PARAMS(:,:))
    DO I = 1, NS-1
       CALL RANDOM_UNIT_VECTORS(INTERNAL_PARAMS(:,:,I))
    END DO
    ! Transform bias values based on the inputs provided such that
    !  there are more rectified units that capture all points, and
    !  fewer that are pushed towards the middle to capture fewer
    !  points.
    IF (BIAS_VALUE .NE. 0.0_REAL64) THEN
       ! Set input biases so that the surface produced by the rectified linear
       !  functions produce well spaced simplices (not too many near middle).
       CALL RANDOM_NUMBER(INPUT_BIAS(:))
       INPUT_BIAS(:) = BIAS_VALUE * INPUT_BIAS(:) ** (1.0_REAL64 / REAL(DS,REAL64))
       ! Internal biases will depend on vector direction, they will
       !  need to be in this range, but this will disable many.
       CALL RANDOM_NUMBER(INTERNAL_BIAS(:,:))
       INTERNAL_BIAS(:,:) = BIAS_VALUE * INTERNAL_BIAS(:,:) * 2.0_REAL64 - BIAS_VALUE
    ELSE
       INPUT_BIAS(:) = 0.0_REAL64
       INTERNAL_BIAS(:,:) = 0.0_REAL64
    END IF
    ! Output layer always default to a bias of zero.
    CALL RANDOM_UNIT_VECTORS(OUTPUT_PARAMS(:,:))
    OUTPUT_BIAS(:) = 0.0_REAL64

  END SUBROUTINE NEW_MODEL

  ! Evaluate the piecewise linear regression model.
  SUBROUTINE EVALUATE(INPUT, OUTPUT)
    USE ISO_FORTRAN_ENV, ONLY: REAL64
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: INPUT
    ! Record of values after the last transformation.
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(:) :: OUTPUT
    ! Compute the number of nodes (to know where "bias" value is).
    INTEGER :: I
    ! Compute the input layer.
    INTERNAL_VALUES(:,1) = MAX(RECTIFY_AT, INPUT_BIAS(:) + &
         MATMUL(INPUT(:), INPUT_PARAMS(:,:)))
    ! Compute the next set of internal values with a rectified activation.
    DO I = 1, MNS-1
       INTERNAL_VALUES(:,I+1) = MAX(RECTIFY_AT, INTERNAL_BIAS(:,I) + &
            MATMUL(INTERNAL_VALUES(:,I), INTERNAL_PARAMS(:,:,I)))
    END DO
    ! Compute the output.
    OUTPUT(:) = OUTPUT_BIAS(:) + &
         MATMUL(INTERNAL_VALUES(:,MNS), OUTPUT_PARAMS(:,:))
  END SUBROUTINE EVALUATE

  ! Compute the gradient of this regression model with respect to its parameters.
  SUBROUTINE GRADIENT(INPUTS, OUTPUT_GRADIENT)
    USE ISO_FORTRAN_ENV, ONLY: REAL64
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: INPUTS
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: OUTPUT_GRADIENT
    INTEGER :: I, J, VEC_START, VEC_END, BIAS_START, BIAS_END
    ! Transfer the output bias graident.
    OUTPUT_BIAS_GRADIENT(:) = OUTPUT_GRADIENT(:)
    ! Compute output parameter gradients.
    DO I = 1, MDO
       OUTPUT_PARAMS_GRADIENT(:,I) = INTERNAL_VALUES(:,MNS) * OUTPUT_GRADIENT(I)
    END DO
    ! Compute the gradient for the last internal vector space values.
    DO I = 1, MDS
       IF (INTERNAL_VALUES(I,MNS) .GT. RECTIFY_AT) THEN
          INTERNAL_VALUES(I,MNS) = DOT_PRODUCT(OUTPUT_PARAMS(I,:), OUTPUT_GRADIENT(:))
       ELSE
          INTERNAL_VALUES(I,MNS) = 0.0_REAL64
       END IF
    END DO
    ! Cycle over all internal layers.
    INTERNAL_REPRESENTATIONS : DO I = MNS-1, 1, -1
       ! Compute the bias values .
       INTERNAL_BIAS_GRADIENT(:,I) = INTERNAL_VALUES(:,I+1)
       ! Compute the indices that will be modified for this layer.
       OUTPUT_DIM : DO J = 1, MDS
          ! There are no updates for parameters with output gradient 0.
          IF (INTERNAL_VALUES(J,I+1) .NE. 0.0_REAL64) THEN
             ! Compute the gradient with respect to one output and
             !  all of its inputs (on the forward eval).
             INTERNAL_PARAMS_GRADIENT(:,J,I) = &
                  INTERNAL_VALUES(:,I) * INTERNAL_VALUES(J,I+1)
          ELSE
             INTERNAL_PARAMS_GRADIENT(:,J,I) = 0.0_REAL64
          END IF
       END DO OUTPUT_DIM
       ! Propogate the gradient backwards to the next previous layer.
       DO J = 1, MDS
          IF (INTERNAL_VALUES(J,I) .GT. RECTIFY_AT) THEN
             INTERNAL_VALUES(J,I) = DOT_PRODUCT( &
                  INTERNAL_PARAMS(J,:,I), INTERNAL_VALUES(:,I+1) )
          ELSE
             INTERNAL_VALUES(J,I) = 0.0_REAL64
          END IF
       END DO
    END DO INTERNAL_REPRESENTATIONS
    ! Compute the input parameter bias values.
    INPUT_BIAS_GRADIENT(:) = INTERNAL_VALUES(:,1)
    ! Compute the next set of internal values with a rectified activation.
    DO I = 1, MDS
       INPUT_PARAMS_GRADIENT(:,I) = INPUTS(:) * INTERNAL_VALUES(I,1)
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
    ! Interface for open mp functions
    INTERFACE
       FUNCTION omp_get_thread_num()
         INTEGER :: omp_get_thread_num
       END FUNCTION omp_get_thread_num
    END INTERFACE
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
          ! Compute the error.
          OUTPUT_GRADIENT(:) = OUTPUT_GRADIENT(:) - Y(:,J)
          ! Update the running average of squared error.
          MEAN_SQUARED_ERROR = SUM(OUTPUT_GRADIENT(:)**2) + MEAN_SQUARED_ERROR
          ! Compute the gradient with respect to this observation.
          CALL GRADIENT(X(:,J), OUTPUT_GRADIENT(:))
          ! Update the running average for the gradient "steps".
          INPUT_PARAMS_GRADIENT_STEP = INPUT_PARAMS_GRADIENT + INPUT_PARAMS_GRADIENT_STEP
          INPUT_BIAS_GRADIENT_STEP   = INPUT_BIAS_GRADIENT   + INPUT_BIAS_GRADIENT_STEP
          INTERNAL_PARAMS_GRADIENT_STEP = INTERNAL_PARAMS_GRADIENT + INTERNAL_PARAMS_GRADIENT_STEP
          INTERNAL_BIAS_GRADIENT_STEP   = INTERNAL_BIAS_GRADIENT   + INTERNAL_BIAS_GRADIENT_STEP
          OUTPUT_PARAMS_GRADIENT_STEP = OUTPUT_PARAMS_GRADIENT + OUTPUT_PARAMS_GRADIENT_STEP
          OUTPUT_BIAS_GRADIENT_STEP   = OUTPUT_BIAS_GRADIENT   + OUTPUT_BIAS_GRADIENT_STEP
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


!2021-02-14 10:27:33
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! Internal variables.                                                                             !
! REAL(KIND=REAL64), PARAMETER :: PI = 3.141592653589793_REAL64                                     !
! REAL(KIND=REAL64), PARAMETER :: PI2 = 2.0_REAL64 * PI                                             !
! REAL(KIND=REAL64), DIMENSION(SIZE(COLUMN_VECTORS,1)-1,SIZE(COLUMN_VECTORS,2)) :: SIGNS            !
! REAL(KIND=REAL64), DIMENSION(SIZE(COLUMN_VECTORS,2)) :: SIN_PRODUCT, TEMP                         !
! INTEGER :: D, I                                                                                   !
! ! Store the dimension.                                                                            !
! D = SIZE(COLUMN_VECTORS,1) - 1                                                                    !
! ! Generate random seed values.                                                                    !
! !  randomly sign the values (1 or -1)                                                             !
! CALL RANDOM_NUMBER(SIGNS(:,:))                                                                    !
! SIGNS(:,:) = SIGNS(:,:) - 0.5_REAL64                                                              !
! SIGNS(:,:) = SIGNS(:,:) / ABS(SIGNS(:,:))                                                         !
! !  randomly get magnitudes proportional to cross-sectional area                                   !
! CALL RANDOM_NUMBER(COLUMN_VECTORS(1:D,:))                                                         !
! COLUMN_VECTORS(1:D,:) = ACOS(SIGNS(:,:) * (1.0_REAL64 - COLUMN_VECTORS(1:D,:)**(1.0_REAL64 / D))) !
! !  randomly assign shifts of +PI                                                                  !
! CALL RANDOM_NUMBER(SIGNS(:,:))                                                                    !
! SIGNS(:,:) = MAX(SIGNS(:,:) - 0.5_REAL64, 0.0_REAL64)                                             !
! SIGNS(:,:) = PI * (SIGNS(:,:) / SIGNS(:,:))                                                       !
! COLUMN_VECTORS(1:D,:) = COLUMN_VECTORS(1:D,:) + SIGNS(:,:)                                        !
!                                                                                                   !
! ! ! Translate the seed values to the range [0,2*pi].                                              !
! ! COLUMN_VECTORS(1:D,:) = PI2 * COLUMN_VECTORS(1:D,:)                                             !
!                                                                                                   !
! ! Initialize the running product of all sin values.                                               !
! SIN_PRODUCT(:) = SIN(COLUMN_VECTORS(1,:))                                                         !
! ! First column is equal to the cos.                                                               !
! COLUMN_VECTORS(1,:) = COS(COLUMN_VECTORS(1,:))                                                    !
! ! Each next column is the product of previous sin values times                                    !
! ! the cos value for that column.                                                                  !
! DO I = 2, D                                                                                       !
!    TEMP(:) = SIN(COLUMN_VECTORS(I,:))                                                             !
!    COLUMN_VECTORS(I,:) = SIN_PRODUCT(:) * COS(COLUMN_VECTORS(I,:))                                !
!    SIN_PRODUCT(:) = SIN_PRODUCT(:) * TEMP(:)                                                      !
! END DO                                                                                            !
! ! The last column is equal to the product of all sin values.                                      !
! COLUMN_VECTORS(D+1,:) = SIN_PRODUCT(:)                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!2021-02-15 01:09:43
!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! WHERE (INTERNAL_VALUES(:,MNS) .GT. 0.0_REAL64)                             !
    !    INTERNAL_VALUES(:,MNS) = MATMUL(OUTPUT_PARAMS(:,:), OUTPUT_GRADIENT(:)) !
    ! ELSEWHERE                                                                  !
    !    INTERNAL_VALUES(:,MNS) = 0.0_REAL64                                     !
    ! END WHERE                                                                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!2021-02-15 01:11:07
!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! WHERE (INTERNAL_VALUES(:,I) .GT. 0.0_REAL64)                                     !
       !    INTERNAL_VALUES(:,I) = MATMUL(INTERNAL_PARAMS(:,:,I), INTERNAL_VALUES(:,I+1)) !
       ! ELSEWHERE                                                                        !
       !    INTERNAL_VALUES(:,I) = 0.0_REAL64                                             !
       ! END WHERE                                                                        !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!2021-02-15 22:53:47
!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! ! ! Store the average gradients.                           !
       ! ! INPUT_PARAMS_GRADIENT    = INPUT_PARAMS_GRADIENT_STEP    !
       ! ! INPUT_BIAS_GRADIENT      = INPUT_BIAS_GRADIENT_STEP      !
       ! ! INTERNAL_PARAMS_GRADIENT = INTERNAL_PARAMS_GRADIENT_STEP !
       ! ! INTERNAL_BIAS_GRADIENT   = INTERNAL_BIAS_GRADIENT_STEP   !
       ! ! OUTPUT_PARAMS_GRADIENT   = OUTPUT_PARAMS_GRADIENT_STEP   !
       ! ! OUTPUT_BIAS_GRADIENT     = OUTPUT_BIAS_GRADIENT_STEP     !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!2021-02-15 23:02:35
!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ! ! Glorot uniform uses:                                                        !
    ! ! !   bound = sqrt(6 / (di + do))                                               !
    ! ! ! Pytorch uniform usese:                                                      !
    ! ! !   bound = sqrt(6 / di)                                                      !
    ! ! !                                                                             !
    ! ! ! Uniform [-bound, bound]                                                     !
    ! ! !                                                                             !
    ! ! ! Initiazlie input parameters.                                                !
    ! ! CALL RANDOM_NUMBER(INPUT_PARAMS(:,:))                                         !
    ! ! INPUT_PARAMS(:,:) = (2.0_REAL64 * INPUT_PARAMS(:,:) - 1.0_REAL64) &           !
    ! !      * SQRT(6.0_REAL64 / REAL(DI+DS, REAL64))                                 !
    ! ! INPUT_BIAS(:) = 0.0_REAL64                                                    !
    ! ! ! Initialize internal parameters.                                             !
    ! ! CALL RANDOM_NUMBER(INTERNAL_PARAMS(:,:,:))                                    !
    ! ! INTERNAL_PARAMS(:,:,:) = (2.0_REAL64 * INTERNAL_PARAMS(:,:,:) - 1.0_REAL64) & !
    ! !      * SQRT(6.0_REAL64 / REAL(DS+DS, REAL64))                                 !
    ! ! INTERNAL_BIAS(:,:) = 0.0_REAL64                                               !
    ! ! ! Initiazlie output parameters.                                               !
    ! ! CALL RANDOM_NUMBER(OUTPUT_PARAMS(:,:))                                        !
    ! ! OUTPUT_PARAMS(:,:) = (2.0_REAL64 * OUTPUT_PARAMS(:,:) - 1.0_REAL64) &         !
    ! !      * SQRT(6.0_REAL64 / REAL(DS+DO, REAL64))                                 !
    ! ! OUTPUT_BIAS(:) = 0.0_REAL64                                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  ! ! Evaluate the piecewise linear regression model.
  ! SUBROUTINE EVALUATE(INPUT, OUTPUT)
  !   USE ISO_FORTRAN_ENV, ONLY: REAL64
  !   REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: INPUT
  !   ! Record of values after the last transformation.
  !   REAL(KIND=REAL64), INTENT(OUT), DIMENSION(:) :: OUTPUT
  !   ! Compute the number of nodes (to know where "bias" value is).
  !   INTEGER :: I, J, K
  !   INTEGER, DIMENSION(MDS) :: NONZERO
  !   ! Set all values to zero.
  !   INTERNAL_VALUES(:,:) = 0.0_REAL64
  !   ! Compute the input layer.
  !   INTERNAL_VALUES(:,1) = MAX(RECTIFY_AT, INPUT_BIAS(:) + &
  !        MATMUL(INPUT(:), INPUT_PARAMS(:,:)))
  !   K = 0
  !   DO J = 1, MDS
  !      IF (INTERNAL_VALUES(J,1) .NE. 0.0_REAL64) THEN
  !         K = K + 1
  !         NONZERO(K) = J
  !      END IF
  !   END DO
  !   ! Compute the next set of internal values with a rectified activation.
  !   DO I = 1, MNS-1
  !      INTERNAL_VALUES(:,I+1) = MAX(RECTIFY_AT, INTERNAL_BIAS(:,I) + &
  !           MATMUL(INTERNAL_VALUES(NONZERO(1:K),I), INTERNAL_PARAMS(NONZERO(1:K),:,I)))
  !      K = 0
  !      DO J = 1, MDS
  !         IF (INTERNAL_VALUES(J,I) .NE. 0.0_REAL64) THEN
  !            K = K + 1
  !            NONZERO(K) = J
  !         END IF
  !      END DO
  !   END DO
  !   ! Compute the output.
  !   OUTPUT(:) = OUTPUT_BIAS(:) + &
  !        MATMUL(INTERNAL_VALUES(NONZERO(1:K),MNS), OUTPUT_PARAMS(NONZERO(1:K),:))
  ! END SUBROUTINE EVALUATE
