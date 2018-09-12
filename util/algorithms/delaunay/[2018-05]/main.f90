
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
! CHAR(32) : filename = input file path (32 chars max)
! REAL : p = d dimensional vector, space separated components
!
USE VTdelaunay_MOD
IMPLICIT NONE

! data
INTEGER :: d, n, m
REAL(KIND=R8), ALLOCATABLE :: pts(:,:), weights(:,:), p(:,:)
REAL(KIND=R8) :: minp, maxp, step, rnorm(1)
INTEGER, ALLOCATABLE :: simp(:,:), error(:)
INTEGER :: i
CHARACTER(len=80) :: arg

! get d, n, minp, maxp, step, filename
CALL GET_COMMAND_ARGUMENT(1, arg)
READ(arg, *) d
CALL GET_COMMAND_ARGUMENT(2, arg)
READ(arg, *) n
CALL GET_COMMAND_ARGUMENT(3, arg)

! number of interpolation points
m = 1

! allocate
ALLOCATE(error(1), pts(d,n), simp(d+1,1), weights(d+1,1), p(d,1), STAT=i)
IF (i /= 0) PRINT *, 'Warning: Malloc Fail!'

! read data
OPEN(1, FILE=TRIM(arg))
DO i = 1, n
        READ(1, *) pts(:, i)
END DO
CLOSE(1)

! create interpolation point
DO i = 1, d
        CALL GET_COMMAND_ARGUMENT(3+i, arg)
        READ(arg, *) p(i,1)
END DO

! get results
CALL delaunayinterp(d, n, pts, 1, p, simp, weights, error, extrap_opt=2.0_R8,&
        rnorm_opt=rnorm)
! print results
IF(error(1) == 0) THEN
        PRINT *, 'Interpolation point: ', p(:,1)
        PRINT *, 'Simplex: ', simp(:,1)
        PRINT *, 'Weights: ', weights(:,1)
ELSE IF(error(1) == 1) THEN
        PRINT *, 'Extrapolation point: ', p(:,1)
        PRINT *, 'Simplex: ', simp(:,1)
        PRINT *, 'Weights: ', weights(:,1)
        PRINT *, 'Residual: ', rnorm
ELSE
        PRINT *, 'Error. Code = ', error(1)
END IF

! free heap
DEALLOCATE(error, pts, simp, weights, p, STAT=i)
IF (i /= 0) PRINT *, 'Warning: Mem Free Fail!'

! done
RETURN
END PROGRAM test
