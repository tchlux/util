PROGRAM TEST
  INTEGER, DIMENSION(6) :: C
  INTEGER, DIMENSION(3,2) :: A, B
  INTEGER :: I, J

  DO J = 1, SIZE(A,2)
     DO I = 1, SIZE(A,1)
        A(I,J) = I + (J-1)*SIZE(A,1)
     END DO
  END DO
  
  B(:,:) = 0

  C(:) = -1

  PRINT *, 'A(:,:)', A(:,:)
  PRINT *, 'C(:)  ', C(:)

  CALL ADD(RESHAPE(C,[3,2]), A)

  PRINT *, 'A(:,:)', A(:,:)
  PRINT *, 'C(:)  ', C(:)
  
CONTAINS

  SUBROUTINE ADD(A, B)
    INTEGER, DIMENSION(:,:) :: A, B
    A(:,:) = A(:,:) + B(:,:)
  END SUBROUTINE ADD

END PROGRAM TEST
