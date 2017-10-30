PROGRAM test_bounds
    INTEGER :: length
    REAL, DIMENSION(10) :: a
    REAL, DIMENSION(5)  :: b
    REAL :: sum
    
    ! Initialize arrays
    call RANDOM_NUMBER(a)
    b = 0

    PRINT *, 'A:', a
    length = 20

    CALL array_add(5, a, b, sum)
    PRINT *, 'TRUE SUM:', sum
    CALL array_add(length, a, b, sum)
    PRINT *, 'BAD SUM: ', sum
  

CONTAINS

  SUBROUTINE array_add(length, a, b, sum)
    INTEGER, INTENT(IN) :: length
    REAL, INTENT(IN), DIMENSION(length) :: a
    REAL, INTENT(IN), DIMENSION(SIZE(a)) :: b
    REAL, INTENT(OUT) :: sum
    INTEGER :: i
    sum = 0
    DO i = 1, length
       sum = sum + a(i) + b(i)
    END DO
  END SUBROUTINE array_add

END PROGRAM test_bounds
