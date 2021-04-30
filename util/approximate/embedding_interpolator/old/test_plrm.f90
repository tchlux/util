
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
PROGRAM TEST
  USE ISO_FORTRAN_ENV, ONLY: REAL32
  USE PLRM
  USE TIMER, ONLY: START_TIMER, STOP_TIMER, TOTAL, RESET_TIMER
  IMPLICIT NONE
  ! Global settings.
  INTEGER, PARAMETER :: N = 2**3, DI = 2**1, DO = 1
  REAL(KIND=REAL32), PARAMETER :: PI = 3.1415926535897932
  INTEGER, PARAMETER :: DS = 3, NS = 2, STEPS = 20
  ! Local variables.
  REAL(KIND=REAL32), DIMENSION(:,:), ALLOCATABLE :: X, Y
  REAL(KIND=REAL32), DIMENSION(:), ALLOCATABLE :: MEAN, STDEV
  REAL(KIND=REAL32) :: MSE
  INTEGER :: I
  INTEGER, DIMENSION(:), ALLOCATABLE :: SEED_ARRAY
  WRITE (*,*) ' '
  CALL RESET_TIMER()

  ! Set a random seed, if one was provided (otherwise leave default).
  CALL RANDOM_SEED(SIZE=I)
  ALLOCATE(SEED_ARRAY(I))
  SEED_ARRAY(:) = 0
  CALL RANDOM_SEED(PUT=SEED_ARRAY(:))

  WRITE (*,*) 'Allocating storage..' 
  CALL START_TIMER()
  ALLOCATE(X(1:DI,1:N), Y(1:DO,1:N), MEAN(1:DI), STDEV(1:DI))
  CALL STOP_TIMER()
  ! WRITE (*,'("  ",1F10.6,/)') TOTAL()

  WRITE (*,*) 'Creating and initializing new model..'
  CALL START_TIMER()
  CALL NEW_MODEL(DI, DS, NS, DO)
  CALL INIT_MODEL(SEED=0)
  CALL STOP_TIMER()
  ! WRITE (*,'("  ",1F10.6,/)') TOTAL()

  WRITE (*,*) 'Generating training data..'
  CALL START_TIMER()
  ! Generate random X and Y values.
  CALL RANDOM_NUMBER(X(:,:))
  X(:,:) = 2.0_REAL32 * PI * X(:,:)
  DO I = 1, DO
     Y(DO,:) = REAL(I-1,REAL32) + SIN(NORM2(X(:,:), DIM=1))
  END DO
  CALL STOP_TIMER()
  ! WRITE (*,'("  ",1F10.6,/)') TOTAL()

  WRITE (*,*) 'Normalizing input data to have zero mean and unit variance..'
  CALL START_TIMER()
  ! Make the X data unit mean.
  MEAN(:) = SUM(X(:,:), DIM=2) / REAL(N,REAL32)
  DO I = 1, DI
     X(I,:) = X(I,:) - MEAN(I)
  END DO
  ! Make the X data unit variance.
  STDEV(:) = SQRT(SUM(X(:,:)**2, DIM=2) / REAL(N,REAL32))
  DO I = 1, DI
     IF (STDEV(I) .EQ. 0) CYCLE
     X(I,:) = X(I,:) / STDEV(I)
  END DO
  CALL STOP_TIMER()
  ! WRITE (*,'("  ",1F10.6,/)') TOTAL()

  WRITE (*,*) 'Creating and initializing new model based on data..'
  CALL START_TIMER()
  CALL NEW_MODEL(DI, DS, NS, DO)
  CALL INIT_MODEL(INPUTS=X, OUTPUTS=Y, SEED=0)
  CALL STOP_TIMER()
  ! WRITE (*,'("  ",1F10.6,/)') TOTAL()

  WRITE (*,*) 'Training model..'
  CALL START_TIMER()
  CALL MINIMIZE_MSE(X(:,:), Y(:,:), STEPS=STEPS, MEAN_SQUARED_ERROR=MSE, &
       NUM_THREADS=2)
  CALL STOP_TIMER()
  ! WRITE (*,'("  ",1F10.6,/)') TOTAL()
  WRITE (*,*) MSE

  ! WRITE (*,*) 'Time spent:'
  ! WRITE (*,"(A14,F10.6)")   '   DGEMM time: ', DGEMM_TIMER
  ! WRITE (*,"(A14,F10.6,/)") '   OTHER time: ', OTHER_TIMER

END PROGRAM TEST
