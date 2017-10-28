
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
IMPLICIT NONE

! data
INTEGER :: d, n, m
REAL(KIND=R8), ALLOCATABLE :: pts(:,:), weights(:,:), p(:,:), work(:)
REAL(KIND=R8) :: minp, maxp, step
INTEGER, ALLOCATABLE :: simp(:,:), error(:)
INTEGER :: i, j
CHARACTER(len=80) :: arg

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

! compute problem dimensions
m = (FLOOR((maxp - minp) / step) + 1)**d

! allocate
ALLOCATE(error(m), pts(d,n), simp(d+1,m), weights(d+1,m), p(d,m), &
        work(MAX(5*d,m*d)), STAT=i)
IF (i /= 0) PRINT *, 'Warning: Malloc Fail!'

! read data
OPEN(1, FILE=TRIM(arg))
DO i = 1, n
        READ(1, *) pts(:, i)
END DO
CLOSE(1)

! create interpolation point
p(:,1) = minp
DO i = 2, m
        ! increment value by step
        p(:,i) = p(:,i-1) 
        p(d,i) = p(d,i) + step
        ! check if its in range
        j = d
        DO WHILE(p(j,i) > maxp)
                ! if not increment next dimension
                p(j,i) = minp
                j = j-1
                p(j,i) = p(j,i) + step
        END DO
END DO

! get results
CALL delaunayp(d, pts, p, work, simp, weights, error)
! print results
DO i = 1, m
IF(error(i) == 0) THEN
        PRINT *, 'Interpolation point: ', p(:,i)
        PRINT *, 'Simplex: ', simp(:,i)
        PRINT *, 'Weights: ', weights(:,i)
ELSE IF(error(i) == 1) THEN
        PRINT *, 'Extrapolation point: ', p(:,i)
        PRINT *, 'Simplex: ', simp(:,i)
        PRINT *, 'Weights: ', weights(:,i)
ELSE
        PRINT *, 'Error. Code = ', error(i)
END IF
END DO

! free heap
DEALLOCATE(error, pts, simp, weights, p, work, STAT=i)
IF (i /= 0) PRINT *, 'Warning: Mem Free Fail!'

! done
RETURN
END PROGRAM test
