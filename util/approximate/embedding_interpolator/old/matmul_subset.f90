PROGRAM MATMUL_SUBSET
  IMPLICIT NONE
  INTEGER, PARAMETER :: ROWS = 3, COLS = 2
  REAL :: A(ROWS,COLS), B(ROWS,COLS), C(ROWS,COLS)
  REAL, ALLOCATABLE, DIMENSION(:,:) :: TEMP
  LOGICAL :: LI(3,2)
  INTEGER :: NR(3*2), NC(3*2)
  INTEGER :: ROW, COL, COUNT
  REAL :: T1, T2

  ! Initialize the matrices with random numbers.
  CALL RANDOM_NUMBER(A(:,:))
  CALL RANDOM_NUMBER(B(:,:))

  WHERE (A .LE. 0.5) A(:,:) = 0.0
  WHERE (B .LE. 0.5) B(:,:) = 0.0

  ! Print the matrices.
  PRINT *, 'A:' 
  DO ROW = 1, SIZE(A,1)
     PRINT *, '  ', A(ROW,:)
  END DO
  PRINT *, 'B:'
  DO ROW = 1, SIZE(B,1)
     PRINT *, '  ', B(ROW,:)
  END DO

  ! Get the logical indices where A is greater than 0.5
  LI(:,:) = A(:,:) .GT. 0.5
  PRINT *, 'I:'
  DO ROW = 1, SIZE(LI,1)
     PRINT *, '  ', LI(ROW,:)
  END DO

  ! Get the numeric indices where A is greater than 0.5
  COUNT = 0
  DO ROW = 1, SIZE(A,1)
     DO COL = 1, SIZE(A,2)
        IF (A(ROW,COL) .GT. 0.5) THEN
           COUNT = COUNT + 1
           NR(COUNT) = ROW
           NC(COUNT) = COL
        END IF
     END DO
  END DO

  PRINT *, 'COUNT:', COUNT
  PRINT *, 'NR:   ', NR(1:COUNT)
  PRINT *, 'NC:   ', NC(1:COUNT)
  ! TEMP = A(:,NC(1:COUNT))
  ! PRINT *, 'A(NR,NC):', SHAPE(TEMP) ! (NR(1:COUNT),:)

  PRINT *, 'AS:' 
  DO ROW = 1, SIZE(A,1)
     PRINT *, '  ', A(ROW,:) .GT. 0.5
  END DO
  PRINT *, 'BS:'
  DO ROW = 1, SIZE(B,1)
     PRINT *, '  ', B(ROW,:) .GT. 0.5
  END DO

  PRINT *, ' trying ..'


  TEMP = TRANSPOSE(0.0 * A(:,:))

  WHERE ((A .GT. 0.5) .AND. (B .GT. 0.5))
     TEMP = TRANSPOSE(A)
     ! TEMP = MATMUL(TRANSPOSE(A),B)
  END WHERE

  PRINT *, ' success!'

  PRINT *, 'TEMP:', SHAPE(TEMP)
  DO ROW = 1, SIZE(TEMP,1)
     PRINT *, '  ', TEMP(ROW,:)
  END DO

END PROGRAM MATMUL_SUBSET
