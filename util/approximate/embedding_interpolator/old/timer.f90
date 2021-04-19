! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!  A module for timing the execution of blocks of code. Uses
!  OMP_GET_WTIME to measure elapsed time, if this returns zero
!  elapsed time then it uses CPU_TIME as a fallback.
MODULE TIMER
  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 1
  REAL, DIMENSION(N) :: CPU_START, CPU_END, WALL_START, WALL_END
  INTEGER, DIMENSION(N) :: HITS

  ! Private variables.
  PRIVATE :: N

  ! Interface for the wall-clock time (use OpenMP standard).
  INTERFACE
     FUNCTION OMP_GET_WTIME()
       USE ISO_FORTRAN_ENV, ONLY: REAL64
       REAL(KIND=REAL64) :: OMP_GET_WTIME
     END FUNCTION OMP_GET_WTIME
  END INTERFACE

CONTAINS  

  ! Reset either a specific timer or all timers (if no arguments given).
  SUBROUTINE RESET_TIMER(I)
    INTEGER, INTENT(IN), OPTIONAL :: I    
    IF (PRESENT(I)) THEN
       CPU_START(I) = 0.0
       CPU_END(I) = 0.0
       WALL_START(I) = 0.0
       WALL_END(I) = 0.0
    ELSE
       CPU_START(:) = 0.0
       CPU_END(:) = 0.0
       WALL_START(:) = 0.0
       WALL_END(:) = 0.0
    END IF
  END SUBROUTINE RESET_TIMER

  ! Start a specific timer or the first timer if no argument is given.
  SUBROUTINE START_TIMER(I)
    INTEGER, INTENT(IN), OPTIONAL :: I
    REAL :: CPU_DIFF, WALL_DIFF
    IF (PRESENT(I)) THEN
       CPU_DIFF = (CPU_END(I) - CPU_START(I))
       WALL_DIFF = (WALL_END(I) - WALL_START(I))
       CALL CPU_TIME(CPU_START(I))
       WALL_START(I) = OMP_GET_WTIME()
       CPU_START(I) = CPU_START(I) - CPU_DIFF
       WALL_START(I) = WALL_START(I) - WALL_DIFF
    ELSE
       CPU_DIFF = (CPU_END(1) - CPU_START(1))
       WALL_DIFF = (WALL_END(1) - WALL_START(1))
       CALL CPU_TIME(CPU_START(1))
       WALL_START(1) = OMP_GET_WTIME()
       CPU_START(1) = CPU_START(1) - CPU_DIFF
       WALL_START(1) = WALL_START(1) - WALL_DIFF
    END IF
  END SUBROUTINE START_TIMER

  ! End a specific timer or the first timer if no argument is given.
  SUBROUTINE STOP_TIMER(I)
    INTEGER, INTENT(IN), OPTIONAL :: I
    IF (PRESENT(I)) THEN
       CALL CPU_TIME(CPU_END(I))
       WALL_END(I) = OMP_GET_WTIME()
       HITS(I) = HITS(I) + 1
    ELSE
       CALL CPU_TIME(CPU_END(1))
       WALL_END(1) = OMP_GET_WTIME()
       HITS(1) = HITS(1) + 1
    END IF
  END SUBROUTINE STOP_TIMER

  FUNCTION TOTAL()
    REAL, DIMENSION(N) :: TOTAL
    TOTAL = WALL_END(:) - WALL_START(:)
    WHERE (TOTAL(:) .LE. 0.0)
       TOTAL(:) = CPU_END(:) - CPU_START(:)
    END WHERE
  END FUNCTION TOTAL

END MODULE TIMER
! --------------------------------------------------------------------
