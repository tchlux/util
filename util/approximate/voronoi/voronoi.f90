!==========================================
!     Full Naive Evaluation of Voronoi     
!==========================================

SUBROUTINE EVAL_MESH(int_prods, app_int_prods, mesh)
  USE ISO_FORTRAN_ENV, ONLY : REAL64
  IMPLICIT NONE
  ! Subroutine for evaluating a standard voronoi mesh at a set of
  ! points given the products have already been calculated.
  REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: int_prods
  REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: app_int_prods
  REAL(KIND=REAL64), DIMENSION(SIZE(app_int_prods,1),SIZE(int_prods,1)), &
       INTENT(OUT) :: mesh
  ! Local variables
  REAL(KIND=REAL64) :: max_ratio, ratio, numerator, denominator
  INTEGER :: y, xi, c
  ! Loop for computing the mesh
  approx_points : DO y = 1, SIZE(app_int_prods,1)
     !$omp parallel do private(xi,c,max_ratio,ratio,numerator,denominator)
     interp_points : DO xi = 1, SIZE(int_prods,1)
        max_ratio = 0
        cell_centers : DO c = 1, SIZE(int_prods,1)
           ! Cycle when we are considering a cell's center and itself
           IF (c .EQ. xi) CYCLE
           ! Calculate the ratio 
           numerator = app_int_prods(y,c) - app_int_prods(y,xi) - &
                int_prods(xi,c) + int_prods(xi,xi)
           denominator = int_prods(c,c) - int_prods(c,xi) - &
                int_prods(xi,c) + int_prods(xi,xi)
           valid_denom : IF (ABS(denominator) .GT. EPSILON(ratio)) THEN
              ratio = numerator / denominator
              track_max : IF (ratio > max_ratio) THEN
                 max_ratio = ratio
              END IF track_max
           END IF valid_denom
        END DO cell_centers
        ! Store the ratio for this cell in the mesh
        mesh(y,xi) = MAX(REAL(1,KIND=REAL64) - max_ratio, REAL(0,KIND=REAL64))
     END DO interp_points
     !$omp end parallel do
     ! Normalize the mesh across the interpolation points (to be convex combinations)
     mesh(y,:) = mesh(y,:) / SUM(mesh(y,:))
  END DO approx_points
END SUBROUTINE EVAL_MESH

!===============================================
!     Cached Dynamic Computation of Voronoi     
!===============================================

SUBROUTINE MAKE_HUGE(DOTS)
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  ! Given a REAL64 matrix, assign the maximum value to all elements.
  REAL(KIND=REAL64), INTENT(OUT), DIMENSION(:,:) :: DOTS  
  DOTS = HUGE(DOTS(1,1))
END SUBROUTINE MAKE_HUGE

SUBROUTINE PREDICT(POINTS, DOTS, EVAL_PT, WEIGHTS, ERROR)
  ! Given column-vector POINTS, pairwise dot product matrix DOTS, and
  ! a single evaluation point EVAL_PT, compute the WEIGHTS associated
  ! with each input point.
  USE ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE
  ! Inputs and outputs
  REAL(KIND=REAL64), INTENT(IN),    DIMENSION(:,:)                           :: POINTS
  REAL(KIND=REAL64), INTENT(IN),    DIMENSION(SIZE(POINTS,1))                :: EVAL_PT
  REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(SIZE(POINTS,2),SIZE(POINTS,2)) :: DOTS
  REAL(KIND=REAL64), INTENT(OUT),   DIMENSION(SIZE(POINTS,2))                :: WEIGHTS
  INTEGER,           INTENT(OUT)                                             :: ERROR
  ! Local variables
  INTEGER :: I_CENTER, I_OTHER, IC, IO, I
  INTEGER, DIMENSION(SIZE(POINTS,2))           :: INDICES
  REAL(KIND=REAL64), DIMENSION(SIZE(POINTS,2)) :: PROJECTIONS, DISTANCES
  REAL(KIND=REAL64), DIMENSION(SIZE(POINTS,1)) :: TO_EVAL, TO_OTHER
  LOGICAL, DIMENSION(SIZE(POINTS,2)) :: SKIP
  REAL(KIND=REAL64) :: DIST_C_TO_PT, DIST_C_TO_O, RATIO, MAX_RATIO
  REAL(KIND=REAL64), PARAMETER :: MAX_VAL = HUGE(MAX_VAL)
  REAL(KIND=REAL64), PARAMETER :: MIN_VAL = EPSILON(MIN_VAL)
  ! Initiallly assume no points should be skipped
  SKIP = .FALSE.
  WEIGHTS = 0.
  ERROR = 0
  ! Distances to all points from evaluation point.
  DO I = 1, SIZE(POINTS, 2)
     DISTANCES(I) = SUM((EVAL_PT - POINTS(:,I)) ** 2)
  END DO
  ! Get the sorted list of points (by distance to evaluation point)
  CALL QSORTC(DISTANCES, INDICES)
  ! Identification of weights
  centers : DO I_CENTER = 1, SIZE(INDICES)
     ! Get the actual index of the center point.
     IC = INDICES(I_CENTER)
     IF (SKIP(IC)) CYCLE
     MAX_RATIO = 0.
     others : DO I_OTHER = 1, SIZE(INDICES)
        ! Get the actual index of the boundary point.
        IO = INDICES(I_OTHER)
        IF (IC .EQ. IO) CYCLE others
        ! Check to see if the dot products have been computed (cache them).
        IF (DOTS(IC,IC) .EQ. MAX_VAL) THEN
           DOTS(IC,IC) = SUM(POINTS(:,IC)**2)
        END IF
        IF (DOTS(IO,IO) .EQ. MAX_VAL) THEN
           DOTS(IO,IO) = SUM(POINTS(:,IO)**2)
        END IF
        IF (DOTS(IC,IO) .EQ. MAX_VAL) THEN
           DOTS(IC, IO) = SUM(POINTS(:,IC) * POINTS(:,IO))
           DOTS(IO, IC) = DOTS(IC, IO)
        END IF
        ! Compute the projected distance (center -> interpolation point)
        DIST_C_TO_PT = SUM(EVAL_PT(:) * POINTS(:,IO)) + DOTS(IC,IC) - &
             (SUM(EVAL_PT(:) * POINTS(:,IC)) + DOTS(IC,IO))
        ! Compute the projected distance (center -> other)
        DIST_C_TO_O = DOTS(IO,IO) + DOTS(IC,IC) - &
             (DOTS(IO,IC) + DOTS(IC,IO))
        ! If the denominator is not zero (should never happen for
        ! unique interpolation points) then...
        IF (ABS(DIST_C_TO_O) .GT. MIN_VAL) THEN
           RATIO = DIST_C_TO_PT / DIST_C_TO_O
           IF (RATIO .GT. MAX_RATIO) THEN; MAX_RATIO = RATIO
           ELSE IF (RATIO .LT. 0.) THEN;   SKIP(IO) = .TRUE.
           ELSE IF (RATIO .GE. 1.) THEN
              SKIP(IC) = .TRUE.
              WEIGHTS(IC) = 0.
              CYCLE centers
           END IF
        ELSE
           ! Duplicated center points, error! Return!
           ERROR = 1
           RETURN
        END IF
     END DO others
     WEIGHTS(IC) = MAX(1. - MAX_RATIO, 0.)
  END DO centers
  ! Normalize the weights to be convex.
  WEIGHTS = WEIGHTS(:) / SUM(WEIGHTS)
  RETURN
  
CONTAINS

  SUBROUTINE QSORTC(A, IDX)
    ! This is a QuickSort routine adapted from Orderpack 2.0.
    !
    ! Also, this implementation incorporates ideas from "A Practical Introduction
    ! to Data Structures and Algorithm Analysis", by Clifford Shaffer.
    !
    ! It sorts real numbers into ascending numerical order
    ! and keeps an index of the value's original array position.
    !
    ! Author: Will Thacker, Winthrop University, July 2013.
    !
    ! QSORTC sorts the real array A and keeps the index of the value's
    ! original array position along with the value (integer array IDX).
    !
    ! On input:
    !
    ! A(:) is the array to be sorted.
    !
    ! On output:
    !
    ! A(:) is sorted.
    !
    ! IDX(1:SIZEOF(A)) contains the original positions of the sorted values.
    !    I.e., sorted(i) = orginal_unsorted(IDX(i)).
    !
    !
    REAL(KIND=REAL64), DIMENSION(:), INTENT(IN OUT):: A
    INTEGER, DIMENSION(SIZE(A)), INTENT(OUT):: IDX
    ! Local variables
    INTEGER:: I   ! Loop iteration variable.
    ! Initialize the array of original positions.
    FORALL (I=1:SIZE(A)) IDX(I)=I
    ! Call the sorting helper
    CALL QSORTC_HELPER(A, IDX, 1, SIZE(A))
  END SUBROUTINE QSORTC
  
  
  RECURSIVE SUBROUTINE QSORTC_HELPER(A, IDX, ISTART, ISTOP)
    ! This internal recursive subroutine performs the recursive quicksort
    ! algorithm.  It is needed because the boundaries of the part of the
    ! array being sorted change and the initial call to the sort routine
    ! does not need to specify boundaries since, generally, the user will
    ! want to sort the entire array passed.
    !
    ! On input:
    !
    ! A(:) contains a subpart to be sorted.
    !
    ! IDX(i) contains the initial position of the value A(i) before sorting.
    !
    ! ISTART is the starting position of the subarray to be sorted.
    !
    ! ISTOP is the ending position of the subarray to be sorted.
    !
    ! On output:
    !
    ! A(ISTART:ISTOP) will be sorted.
    !
    ! IDX(i) contains the original position for the value at A(i).
    !
    !
    REAL(KIND=REAL64), DIMENSION(:), INTENT (IN OUT):: A
    INTEGER, DIMENSION(SIZE(A)), INTENT(IN OUT):: IDX
    INTEGER, INTENT(IN):: ISTART, ISTOP

    !  Local variables
    INTEGER:: ILEFT ! A position on the left to be swapped with value at IRIGHT.
    INTEGER:: IMID ! The middle position used to select the pivot.
    INTEGER:: IRIGHT ! A position on the right to be swapped with value at ILEFT.
    INTEGER:: ITEMP  ! Used for swapping within IDX.
    REAL(KIND=REAL64):: ATEMP ! Used for swapping.
    REAL(KIND=REAL64):: PIVOT ! Holds the temporary pivot.

    ! INSMAX is used to stop recursively dividing the array and to instead
    ! use a sort that is more efficient for small arrays than quicksort.
    !
    ! The best cutoff point is system dependent.

    INTEGER, PARAMETER:: INSMAX=24

    ! Check to see if we have enough values to make quicksort useful.
    ! Otherwise let the insertion sort handle it.

    IF ((ISTOP - ISTART) < INSMAX) THEN
       CALL INSERTION(A, IDX, ISTART, ISTOP)
    ELSE

       ! Use the median of the first, middle and last items for the pivot
       ! and place the median (pivot) at the end of the list.
       ! Putting it at the end of the list allows for a guard value to keep
       ! the loop from falling off the right end of the array (no need to
       ! check for at the end of the subarray EACH time through the loop).

       IMID = (ISTART + ISTOP)/2

       IF (A(ISTOP) < A(ISTART)) THEN
          ATEMP = A(ISTART)
          A(ISTART) = A(ISTOP)
          A(ISTOP) = ATEMP

          ITEMP = IDX(ISTART)
          IDX(ISTART) = IDX(ISTOP)
          IDX(ISTOP) = ITEMP
       END IF

       IF (A(IMID) < A(ISTOP)) THEN
          ATEMP = A(ISTOP)
          A(ISTOP) = A(IMID)
          A(IMID) = ATEMP

          ITEMP = IDX(ISTOP)
          IDX(ISTOP) = IDX(IMID)
          IDX(IMID) = ITEMP

          IF (A(ISTOP) < A(ISTART)) THEN
             ATEMP = A(ISTOP)
             A(ISTOP) = A(ISTART)
             A(ISTART) = ATEMP

             ITEMP = IDX(ISTOP)
             IDX(ISTOP) = IDX(ISTART)
             IDX(ISTART) = ITEMP
          END IF
       END IF

       ! Now, the first position has a value that is less or equal to the
       ! partition. So, we know it belongs in the left side of the partition
       ! and we can skip it. Also, the pivot is at the end.  So, there is
       ! no need to compare the pivot with itself.

       PIVOT = A(ISTOP)
       ILEFT = ISTART + 1
       IRIGHT = ISTOP - 1

       DO WHILE (ILEFT < IRIGHT)
          ! Find a value in the left side that is bigger than the pivot value.
          ! Pivot is at the right end so ILEFT will not fall off the end
          ! of the subarray.

          DO WHILE (A(ILEFT) < PIVOT)
             ILEFT = ILEFT + 1
          END DO

          DO WHILE (IRIGHT .NE. ILEFT)
             IF (A(IRIGHT) .LT. PIVOT) EXIT
             IRIGHT = IRIGHT - 1
          END DO

          ! Now we have a value bigger than pivot value on the left side that can be
          ! swapped with a value smaller than the pivot on the right side.
          !
          ! This gives us all values less than pivot on the left side of the
          ! array and all values greater than the pivot on the right side.

          ATEMP = A(IRIGHT)
          A(IRIGHT) = A(ILEFT)
          A(ILEFT) = ATEMP

          ITEMP = IDX(IRIGHT)
          IDX(IRIGHT) = IDX(ILEFT)
          IDX(ILEFT) = ITEMP

       END DO
       !
       ! The last swap was in error (since the while condition is not checked
       ! until after the swap is done) so we swap again to fix it.

       ! This is done (once) rather than having an if (done many times) in the
       ! loop to prevent the swapping.


       ATEMP = A(IRIGHT)
       A(IRIGHT) = A(ILEFT)
       A(ILEFT) = ATEMP

       ITEMP = IDX(IRIGHT)
       IDX(IRIGHT) = IDX(ILEFT)
       IDX(ILEFT) = ITEMP

       ! Put the pivot value in its correct spot (between the 2 partitions)
       ! When the WHILE condition finishes, ILEFT is greater than IRIGHT.
       ! So, ILEFT has the position of the first value in the right side.
       ! This is where we can put the pivot (and where it will finally rest,
       ! so no need to look at it again).  Also, place the first value of
       ! the right side (being displaced by the pivot) at the end of the
       ! subarray (since it is bigger than the pivot).

       ATEMP = A(ISTOP)
       A(ISTOP) = A(ILEFT)
       A(ILEFT) = ATEMP

       ITEMP = IDX(ISTOP)
       IDX(ISTOP) = IDX(ILEFT)
       IDX(ILEFT) = ITEMP

       CALL QSORTC_HELPER(A, IDX, ISTART, ILEFT-1)
       CALL QSORTC_HELPER(A, IDX, ILEFT+1, ISTOP)
    END IF

  END SUBROUTINE QSORTC_HELPER


  SUBROUTINE INSERTION(A, IDX, ISTART, ISTOP)
    ! This subroutine performs an insertion sort used for sorting
    ! small subarrays efficiently.
    !
    ! This subroutine sorts a subarray of A (between positions ISTART
    ! and ISTOP) keeping a record of the original position (array IDX).

    ! On input:
    !
    ! A(:) contains a subpart to be sorted.
    !
    ! IDX(i) contains the initial position of the value A(i) before sorting.
    !
    ! ISTART is the starting position of the subarray to be sorted.
    !
    ! ISTOP is the ending position of the subarray to be sorted.
    !
    ! On output:
    !
    ! A(ISTART:ISTOP) will be sorted.
    !
    ! IDX(i) contains the original position for the value at A(i).
    !

    REAL(KIND=REAL64), DIMENSION(:), INTENT (IN OUT):: A
    INTEGER, DIMENSION(SIZE(A)), INTENT(IN OUT)::IDX
    INTEGER, INTENT(IN):: ISTART, ISTOP

    ! Local variables.

    REAL(KIND=REAL64):: AMIN  ! Temporary minimum.
    REAL(KIND=REAL64):: ATEMP  ! The value to be inserted.
    INTEGER:: I    ! Index variable.
    INTEGER:: IABOVE ! Index to find insertion point.
    INTEGER:: IMIN ! Temporary minimum position.
    INTEGER:: ITEMP ! Temporary for swapping.

    IF (ISTOP .EQ. ISTART) THEN
       RETURN
    END IF

    ! Find the smallest and put it at the top as a "guard" so there is
    ! no need for the DO WHILE to check if it is going past the top.

    AMIN = A(ISTART)
    IMIN = ISTART

    DO I=ISTOP,ISTART+1,-1
       IF (A(I) < AMIN) THEN
          AMIN = A(I)
          IMIN = I
       END IF
    END DO

    A(IMIN) = A(ISTART)
    A(ISTART) = AMIN

    ITEMP = IDX(ISTART)
    IDX(ISTART) = IDX(IMIN)
    IDX(IMIN) = ITEMP

    ! Insertion sort the rest of the array.
    DO I=ISTART+2,ISTOP
       ATEMP = A(I)
       ITEMP = IDX(I)
       IABOVE = I - 1
       IF (ATEMP < A(IABOVE)) THEN
          A(I) = A(IABOVE)
          IDX(I) = IDX(IABOVE)
          IABOVE = IABOVE - 1

          ! Stop moving items down when the position for "insertion" is found.
          !
          ! Do not have to check for "falling off" the beginning of the
          ! array since the smallest value is a guard value in the first position.

          DO WHILE (ATEMP < A(IABOVE))
             A(IABOVE+1) = A(IABOVE)
             IDX(IABOVE+1) = IDX(IABOVE)
             IABOVE = IABOVE - 1
          END DO
       END IF
       A(IABOVE+1) = ATEMP
       IDX(IABOVE+1) = ITEMP
    END DO
  END SUBROUTINE INSERTION

END SUBROUTINE PREDICT
