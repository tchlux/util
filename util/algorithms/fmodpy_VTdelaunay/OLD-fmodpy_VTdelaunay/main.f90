
! Driver code for computing a delaunay simplex containing a given point.
! 
! Written to FORTRAN 2003 standard.
!
! Tested w/ gfortran 4.8.5 on x86_64-redhat-linux       
! AUTHOR: Tyler Chang                                   
! Last Update: August, 2017                               

PROGRAM test
!
! Driver code: Read a dataset from a file and interpolate along a dense
! grid as specified at command line
! USAGE: ./vtdel d n minp maxp step filename
! INTEGER : d = dimension of problem
! INTEGER : n = lines in file
! REAL : minp = minimum value for p
! REAL : maxp = maximum value for p
! REAL : step = step size for p
! CHAR(32) : filename = input file path (32 chars max)
!
USE VTdelaunay
USE VTgeometry
USE real_precision
IMPLICIT NONE

! data
INTEGER :: d, n
REAL(KIND=R8), ALLOCATABLE :: pts(:,:), weights(:), p(:)
REAL(KIND=R8) :: minp, maxp, step
INTEGER, ALLOCATABLE :: simp(:)
INTEGER :: error, i, j
CHARACTER(len=32) :: arg

! get d, n, minp, maxp, step, filename
CALL GET_COMMAND_ARGUMENT(1, arg)
READ(arg, *) d
CALL GET_COMMAND_ARGUMENT(2, arg)
READ(arg, *) n
CALL GET_COMMAND_ARGUMENT(3, arg)
READ(arg, *) minp
CALL GET_COMMAND_ARGUMENT(4, arg)
READ(arg, *) maxp
CALL GET_COMMAND_ARGUMENT(5, arg)
READ(arg, *) step
CALL GET_COMMAND_ARGUMENT(6, arg)

! allocate
ALLOCATE(pts(d,n), simp(d+1), weights(d+1), p(d), STAT=error)
IF (error /= 0) PRINT *, 'Warning: Malloc Fail!'

! read data
OPEN(1, FILE=TRIM(arg))
DO i = 1, n
        READ(1, *) pts(:, i)
END DO
CLOSE(1)

! create interpolation point
p = minp

! get results for all points
DO WHILE(ALL(p <= maxp))
        ! get results
        CALL delaunayp(d, n, pts, p, simp, weights, error)
        ! Successful interpolation
        IF(error == 0) THEN
                PRINT *, 'Interpolation point: ', p
                PRINT *, 'Simplex: ', simp
                PRINT *, 'Weights: ', weights
        ! Successful extrapolation
        ELSE IF(error == 1) THEN
                PRINT *, 'Extrapolation point: ', p
                PRINT *, 'Simplex: ', simp
                PRINT *, 'Weights: ', weights
        ! error code
        ELSE
                PRINT *, 'Error. Code = ', error
        END IF
        ! increment point
        DO i = 1, d-1
                IF(p(i) + step <= maxp) EXIT
                IF(i < d) p(i) = minp
        END DO
        p(i) = p(i) + step
END DO

! free heap
DEALLOCATE(pts, simp, weights, p, STAT=error)
IF (error /= 0) PRINT *, 'Warning: Mem Free Fail!'

! done
RETURN
END PROGRAM test
