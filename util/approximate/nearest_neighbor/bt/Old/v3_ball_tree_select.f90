MODULE BALL_TREE
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64
  IMPLICIT NONE
  INTEGER :: DIST_CALCS
  ! REAL :: STARTT, STOPT
  ! REAL, DIMENSION(10) :: TIMES = 0.0

CONTAINS

  ! Re-arrange elements of POINTS into a binary ball tree.
  RECURSIVE SUBROUTINE BUILD_TREE(POINTS, SQ_SUMS, RADII, ORDER,&
       ROOT, LEAF_SIZE, COMPUTED_SQ_SUMS)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:,:) :: POINTS
    REAL(KIND=REAL64),   INTENT(OUT),   DIMENSION(:) :: SQ_SUMS, RADII
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL :: ROOT, LEAF_SIZE
    LOGICAL, INTENT(IN), OPTIONAL :: COMPUTED_SQ_SUMS
    ! Local variables
    INTEGER(KIND=INT64) :: CENTER_IDX, MAX_IDX, I, J, LS
    REAL(KIND=REAL64), DIMENSION(SIZE(POINTS,1)) :: PT
    REAL(KIND=REAL64), DIMENSION(SIZE(ORDER)) :: SQ_DISTS
    REAL(KIND=REAL64) :: MAX_SQ_DIST, SQ_DIST
    ! Set the leaf size to 1 by default (most possible work required,
    ! but guarantees successful use with any leaf size).
    IF (PRESENT(LEAF_SIZE)) THEN ; LS = LEAF_SIZE
    ELSE                         ; LS = 1
    END IF
    ! If no squared sums were provided, compute them.
    IF (.NOT. PRESENT(COMPUTED_SQ_SUMS) .OR. &
         .NOT. COMPUTED_SQ_SUMS) THEN
       DO I = 1, SIZE(POINTS,2)
          SQ_SUMS(I) = SUM(POINTS(:,I)**2)
       END DO
    END IF
    ! Set the index of the 'split' 
    IF (PRESENT(ROOT)) THEN ; CENTER_IDX = ROOT
    ELSE
       ! 1) Compute distances between first point (random) and all others.
       ! 2) Pick the furthest point (on conv hull) from first as the center node.
       J = ORDER(1)
       PT(:) = POINTS(:,J)
       CENTER_IDX = 1
       MAX_SQ_DIST = 0.0_REAL64
       DO I = 2, SIZE(ORDER)
          SQ_DIST = SQ_SUMS(J) + SQ_SUMS(ORDER(I)) - &
               2 * DOT_PRODUCT(POINTS(:,ORDER(I)), PT(:))
          ! If this is a new max distance, record this point.
          IF (SQ_DIST .GT. MAX_SQ_DIST) THEN
             MAX_SQ_DIST = SQ_DIST
             CENTER_IDX = I
          END IF
       END DO
       ! Now CENTER_IDX is the selected center for this node in tree.
    END IF

    ! Move the "center" to the first position.
    CALL SWAPI(ORDER(1), ORDER(CENTER_IDX))
    ! Measure squared distance beween "center" node and all other points.
    J = ORDER(1)
    PT(:) = POINTS(:,J)
    SQ_DISTS(1) = 0.0_REAL64
    CENTER_TO_ALL : DO I = 2, SIZE(ORDER)
       SQ_DISTS(I) = SQ_SUMS(J) + SQ_SUMS(ORDER(I)) - &
            2 * DOT_PRODUCT(POINTS(:,ORDER(I)), PT(:))
    END DO CENTER_TO_ALL
    ! Base case for recursion, once we have few enough points, exit.
    IF (SIZE(ORDER) .LE. LS) THEN
       RADII(ORDER(1)) = SQRT(MAXVAL(SQ_DISTS))
       RADII(ORDER(2:)) = 0.0_REAL64
       RETURN
    END IF
    ! Rearrange "SQ_DISTS" about the median value.
    ! Compute the last index that will belong "inside" this node.
    MAX_IDX = (SIZE(ORDER) + 1) / 2
    CALL SELECT(SQ_DISTS, ORDER, MAX_IDX)
    ! Now ORDER has been rearranged such that the median distance
    ! element of POINTS is at the median location.
    ! Identify the furthest point (must be in second half of list).
    I = MAX_IDX + MAXLOC(SQ_DISTS(MAX_IDX+1:),1)
    ! Store the "radius" of this ball, the furthest point.
    RADII(ORDER(1)) = SQRT(SQ_DISTS(I))
    ! Move the furthest point into the spot after the median (outer root).
    CALL SWAPI(ORDER(I), ORDER(MAX_IDX+1))
    ! Move the median point (furthest "interior") to the front (inner root).
    CALL SWAPI(ORDER(2), ORDER(MAX_IDX))

    ! Recurisively create this tree.
    !   build a tree with the root being the furthest from this center
    !   for the remaining "interior" points of this center node.
    CALL BUILD_TREE(POINTS, SQ_SUMS, RADII, ORDER(2:MAX_IDX), 1_INT64, LS, .TRUE.)
    !   build a tree with the root being the furthest from this center
    !   for the remaining "exterior" points of this center node.
    !   Only perform this operation if there are >0 points available.
    IF (MAX_IDX < SIZE(ORDER)) &
         CALL BUILD_TREE(POINTS, SQ_SUMS, RADII, &
         ORDER(MAX_IDX+1:), 1_INT64, LS, .TRUE.)
  END SUBROUTINE BUILD_TREE

  ! Compute the K nearest elements of TREE to each point in POINTS.
  SUBROUTINE NEAREST(POINTS, K, TREE, SQ_SUMS, RADII, ORDER, &
       LEAF_SIZE, INDICES, DISTS)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: POINTS, TREE
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: SQ_SUMS, RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)               :: K, LEAF_SIZE
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(K,SIZE(POINTS,2)) :: DISTS
    ! Local variables.
    INTEGER :: I
    INTEGER(KIND=INT64), DIMENSION(K+LEAF_SIZE+2) :: INDS_BUFFER
    REAL(KIND=REAL64),   DIMENSION(K+LEAF_SIZE+2) :: DISTS_BUFFER
    ! For each point in this set, use the recursive branching
    ! algorithm to identify the nearest elements of TREE.

    DO I = 1, SIZE(POINTS,2)
       DIST_CALCS = 0
       CALL PT_NEAREST(POINTS(:,I), K, TREE, SQ_SUMS, RADII, ORDER, &
            LEAF_SIZE, INDS_BUFFER, DISTS_BUFFER)
       ! Sort the first K elements of the temporary arry for return.
       INDICES(:,I) = INDS_BUFFER(:K)
       DISTS(:,I) = DISTS_BUFFER(:K)
       ! PRINT *, 'BALL_TREE.F90: DIST CALCS =', DIST_CALCS
    END DO


    ! PRINT *, ''
    ! PRINT '("sq sums   = ",f6.3,"s")',TIMES(1)
    ! PRINT '("select    = ",f6.3,"s")',TIMES(2)
    ! PRINT '("recursion = ",f6.3,"s")',TIMES(5)
    ! PRINT *, ''

  END SUBROUTINE NEAREST

  ! Compute the K nearest elements of TREE to each point in POINTS.
  RECURSIVE SUBROUTINE PT_NEAREST(POINT, K, TREE, SQ_SUMS, RADII, &
       ORDER, LEAF_SIZE, INDICES, DISTS, FOUND, PT_SS)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: POINT
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: TREE
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:)   :: SQ_SUMS, RADII
    INTEGER(KIND=INT64), INTENT(IN), DIMENSION(:) :: ORDER
    INTEGER(KIND=INT64), INTENT(IN)               :: K, LEAF_SIZE
    INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(:) :: INDICES
    REAL(KIND=REAL64),   INTENT(OUT), DIMENSION(:) :: DISTS
    INTEGER(KIND=INT64), INTENT(INOUT), OPTIONAL   :: FOUND
    REAL(KIND=REAL64),   INTENT(IN),    OPTIONAL   :: PT_SS
    ! Local variables
    INTEGER(KIND=INT64) :: F, I, I1, I2
    REAL(KIND=REAL64)   :: D, D1, D2, PS
    ! Initialize FOUND for first call, if FOUND is present then
    ! this must not be first and all are present.
    INITIALIZE : IF (PRESENT(FOUND)) THEN
       F = FOUND
       PS = PT_SS
    ELSE
       ! Start at index 0 (added onto current index). Compute squared sum.
       PS = SUM(POINT(:)**2)
       ! Measure distance to root.
       INDICES(1) = ORDER(1)
       DISTS(1) = SQRT(PS + SQ_SUMS(ORDER(1)) - &
            2*DOT_PRODUCT(POINT(:), TREE(:,ORDER(1))))
       DIST_CALCS = DIST_CALCS + 1
       ! Set the "points found" to be 1.
       F = 1
    END IF INITIALIZE
    ! If this is NOT a leaf node, then recurse.
    BRANCH_OR_LEAF : IF (SIZE(ORDER) .GT. LEAF_SIZE) THEN
       ! Measure distance to inner child.
       I1 = ORDER(2)
       D1 = SQRT(PS + SQ_SUMS(I1) - 2*DOT_PRODUCT(POINT(:),TREE(:,I1)))
       DIST_CALCS = DIST_CALCS + 1
       ! Store this distance calculation and index.
       F = F + 1_INT64
       INDICES(F) = I1
       DISTS(F) = D1
       ! Measure distance to outer child the same as above.
       I = (SIZE(ORDER) + 1) / 2_INT64
       I2 = ORDER(I+1_INT64)
       D2 = SQRT(PS + SQ_SUMS(I2) - 2*DOT_PRODUCT(POINT(:),TREE(:,I2)))
       DIST_CALCS = DIST_CALCS + 1
       ! Store this distance calculation and index.
       F = F + 1_INT64
       INDICES(F) = I2
       DISTS(F) = D2
       ! Re-organize the list of closest points, pushing them to first K spots.
       CALL SELECT(DISTS(:F), INDICES(:F), MIN(K,F))
       F = MIN(K,F)
       ! Store the maximum distance.
       D = MAXVAL(DISTS(:F),1)
       ! Determine which child to search (depth-first search) based
       ! on which one is closer to the approximation point.
       INNER_CHILD_CLOSER : IF (D1 .LE. D2) THEN
          ! Search the inner child if it could contain a nearer point.
          SEARCH_INNER1 : IF ((F .LT. K) .OR. (D .GT. D1 - RADII(I1))) THEN
             CALL PT_NEAREST(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(2:I), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_INNER1
          ! Search the outer child if it could contain a nearer point.
          SEARCH_OUTER1 : IF ((F .LT. K) .OR. (D .GT. D2 - RADII(I2))) THEN
             CALL PT_NEAREST(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(I+1_INT64:), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_OUTER1
       ELSE
          ! Search the outer child if it could contain a nearer point.
          SEARCH_OUTER2 : IF ((F .LT. K) .OR. (D .GT. D2 - RADII(I2))) THEN
             CALL PT_NEAREST(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(I+1_INT64:), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_OUTER2
          ! Search the inner child if it could contain a nearer point.
          SEARCH_INNER2 : IF ((F .LT. K) .OR. (D .GT. D1 - RADII(I1))) THEN
             CALL PT_NEAREST(POINT, K, TREE, SQ_SUMS, RADII, &
                  ORDER(2:I), LEAF_SIZE, INDICES, DISTS, F, PS)
          END IF SEARCH_INNER2
       END IF INNER_CHILD_CLOSER

    ! Since this is a leaf node, we measure distance to all children.
    ELSE
       DIST_TO_CHILDREN : DO I = 2, SIZE(ORDER)
          ! Measure distance to all children of this node.
          D = SQRT(PS + SQ_SUMS(ORDER(I)) - &
               2*DOT_PRODUCT(POINT(:),TREE(:,ORDER(I))))
          DIST_CALCS = DIST_CALCS + 1
          ! Store this distance.
          F = F + 1_INT64
          DISTS(F) = D
          INDICES(F) = ORDER(I)
       END DO DIST_TO_CHILDREN
       ! Reduce the kept points to only those that are closest.
       CALL SELECT(DISTS(:F), INDICES(:F), MIN(K,F))
       F = MIN(K, F)
    END IF BRANCH_OR_LEAF

    ! Handle closing operations..
    SORT_K : IF (PRESENT(FOUND)) THEN
       ! This is not the root, We need to pass the updated value of
       ! FOUND back up the recrusion stack.
       FOUND = F
    ELSE
       ! This is the root, initial caller. Sort the distances for return.
       CALL SELECT(DISTS(:F), INDICES(:F), K)
       CALL QSORTC(DISTS(:K), INDICES(:K), .FALSE.)
    END IF SORT_K

  END SUBROUTINE PT_NEAREST


  ! TODO -- Write the following function.
  ! 
  ! ! Swap all nodes at or above "LEVEL" to the front of the TREE,
  ! ! leaving the rest in the tail. The front portion of the resulting
  ! ! set of points is a valid TREE for use in the future.
  ! SUBROUTINE PRUNE(TREE, RADII, LEVEL)
  !   ! Establish the set of indices that will be kept recursively.
  !   ! Establish the set of indices that will be overwritten recursively.
  !   ! Execute the swap (placing all to-keep points into free indices).
  ! END SUBROUTINE PRUNE


  ! Short subroutines for swapping out two values.
  SUBROUTINE SWAPI(V1, V2)
    INTEGER(KIND=INT64), INTENT(INOUT) :: V1, V2
    ! Local temp
    INTEGER(KIND=INT64) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAPI

  SUBROUTINE SWAPR(V1, V2)
    REAL(KIND=REAL64), INTENT(INOUT) :: V1, V2
    ! Local temp
    REAL(KIND=REAL64) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAPR

  SUBROUTINE SWAPP(V1, V2)
    REAL(KIND=REAL64), DIMENSION(:) :: V1, V2
    ! Local temp
    INTEGER(KIND=INT64) :: I
    REAL(KIND=REAL64) :: TEMP
    DO I = 1, SIZE(V1)
       TEMP = V1(I)
       V1(I) = V2(I)
       V2(I) = TEMP
    END DO
  END SUBROUTINE SWAPP

  ! ------------------------------------------------------------------
  !                       FastSelect method
  ! 
  ! Given VALUES list of numbers, rearrange the elements of VALUES
  ! such that the element at index K has rank K (holds its same
  ! location as if all of VALUES were sorted). Symmetrically rearrange
  ! array INDICES to keep track of prior indices.
  ! 
  ! This algorithm uses the same conceptual approach as Floyd-Rivest,
  ! but instead of standard-deviation based selection of bounds for
  ! recursion, a rank-based method is used to pick the subset of
  ! values that is searched. This simplifies the code and improves
  ! interpretability, while achieving the same tunable performance.
  ! 
  ! Arguments:
  ! 
  !   VALUES   --  A 1D array of real numbers.
  !   INDICES  --  A 1D array of original indices for elements of VALUES.
  !   K        --  A positive integer for the rank index about which
  !                VALUES should be rearranged.
  ! Optional:
  ! 
  !   DIVISOR  --  A positive integer >= 2 that represents the
  !                division factor used for large VALUES arrays.
  !   MAX_SIZE --  An integer >= DIVISOR that represents the largest
  !                sized VALUES for which the worst-case pivot value
  !                selection is tolerable. A worst-case pivot causes
  !                O( SIZE(VALUES)^2 ) runtime. This value should be
  !                determined heuristically based on compute hardware.
  ! 
  ! Output:
  ! 
  !   The elements of the array VALUES are rearranged such that the
  !   element at position VALUES(K) is in the same location it would
  !   be if all of VALUES were in sorted order. Also known as,
  !   VALUES(K) has rank K.
  ! 
  RECURSIVE SUBROUTINE SELECT(VALUES, INDICES, K, DIVISOR, MAX_SIZE)
    ! Arguments
    REAL(KIND=REAL64), INTENT(INOUT), DIMENSION(:) :: VALUES
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:) :: INDICES
    INTEGER(KIND=INT64), INTENT(IN)                  :: K
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL        :: DIVISOR, MAX_SIZE
    ! Locals
    INTEGER(KIND=INT64) :: LEFT, RIGHT, L, R, MS, D
    REAL(KIND=REAL64) :: P
    ! Initialize the divisor (for making subsets).
    IF (PRESENT(DIVISOR)) THEN ; D = DIVISOR
    ELSE IF (SIZE(VALUES) .GE. 2**23) THEN ; D = 2**5
    ELSE IF (SIZE(VALUES) .GE. 2**20) THEN ; D = 2**3
    ELSE                                   ; D = 2**2
    END IF
    ! Initialize the max size (before subsets are created).
    IF (PRESENT(MAX_SIZE)) THEN ; MS = MAX_SIZE
    ELSE                        ; MS = 2**10
    END IF
    ! Initialize LEFT and RIGHT to be the entire array.
    LEFT = 1
    RIGHT = SIZE(VALUES)
    ! Loop until done finding the K-th element.
    DO WHILE (LEFT .LT. RIGHT)
       ! Use SELECT recursively to improve the quality of the
       ! selected pivot value for large arrays.
       IF (RIGHT - LEFT .GT. MS) THEN
          ! Compute how many elements should be left and right of K
          ! to maintain the same percentile in a subset.
          L = K - K / D
          R = L + (SIZE(VALUES) / D)
          ! Perform fast select on an array a fraction of the size about K.
          CALL SELECT(VALUES(L:R), INDICES(L:R), K - L + 1, DIVISOR, MAX_SIZE)
       END IF
       ! Pick a partition element at position K.
       P = VALUES(K)
       L = LEFT
       R = RIGHT
       ! Move the partition element to the front of the list.
       CALL SWAPR(VALUES(LEFT), VALUES(K))
       CALL SWAPI(INDICES(LEFT), INDICES(K))
       ! Pre-swap the left and right elements (temporarily putting a
       ! larger element on the left) before starting the partition loop.
       IF (VALUES(RIGHT) .GT. P) THEN
          CALL SWAPR(VALUES(LEFT), VALUES(RIGHT))
          CALL SWAPI(INDICES(LEFT), INDICES(RIGHT))
       END IF
       ! Now partition the elements about the pivot value "T".
       DO WHILE (L .LT. R)
          CALL SWAPR(VALUES(L), VALUES(R))
          CALL SWAPI(INDICES(L), INDICES(R))
          L = L + 1
          R = R - 1
          DO WHILE (VALUES(L) .LT. P) ; L = L + 1 ; END DO
          DO WHILE (VALUES(R) .GT. P) ; R = R - 1 ; END DO
       END DO
       ! Place the pivot element back into its appropriate place.
       IF (VALUES(LEFT) .EQ. P) THEN
          CALL SWAPR(VALUES(LEFT), VALUES(R))
          CALL SWAPI(INDICES(LEFT), INDICES(R))
       ELSE
          R = R + 1
          CALL SWAPR(VALUES(R), VALUES(RIGHT))
          CALL SWAPI(INDICES(R), INDICES(RIGHT))
       END IF
       ! adjust left and right towards the boundaries of the subset
       ! containing the (k - left + 1)th smallest element
       IF (R .LE. K) LEFT = R + 1
       IF (K .LE. R) RIGHT = R - 1
    END DO
  END SUBROUTINE SELECT
  ! ------------------------------------------------------------------


  ! ------------------------------------------------------------------
  !                QSORT method (very fast on average)
  ! 
  SUBROUTINE QSORTC(A, IDX, INIT)
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
    REAL(KIND=REAL64), DIMENSION(:), INTENT(INOUT):: A
    INTEGER(KIND=INT64), DIMENSION(SIZE(A)), INTENT(INOUT):: IDX
    LOGICAL, INTENT(IN), OPTIONAL :: INIT
    ! Local variables
    LOGICAL :: SHOULD_INIT
    INTEGER(KIND=INT64) :: I   ! Loop iteration variable.
    IF (.NOT. PRESENT(INIT) .OR. INIT) THEN
       ! Initialize the array of original positions.
       FORALL (I=1:SIZE(A)) IDX(I)=I
    END IF
    ! Call the sorting helper
    CALL QSORTC_HELPER(A, IDX, 1_INT64, INT(SIZE(A),INT64))
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
    INTEGER(KIND=INT64), DIMENSION(SIZE(A)), INTENT(IN OUT):: IDX
    INTEGER(KIND=INT64), INTENT(IN):: ISTART, ISTOP

    !  Local variables
    INTEGER(KIND=INT64):: ILEFT ! A position on the left to be swapped with value at IRIGHT.
    INTEGER(KIND=INT64):: IMID ! The middle position used to select the pivot.
    INTEGER(KIND=INT64):: IRIGHT ! A position on the right to be swapped with value at ILEFT.
    INTEGER(KIND=INT64):: ITEMP  ! Used for swapping within IDX.
    REAL(KIND=REAL64):: ATEMP ! Used for swapping.
    REAL(KIND=REAL64):: PIVOT ! Holds the temporary pivot.

    ! INSMAX is used to stop recursively dividing the array and to instead
    ! use a sort that is more efficient for small arrays than quicksort.
    !
    ! The best cutoff point is system dependent.

    INTEGER(KIND=INT64), PARAMETER:: INSMAX=24

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

    REAL(KIND=REAL64), DIMENSION(:), INTENT (IN OUT):: A
    INTEGER(KIND=INT64), DIMENSION(SIZE(A)), INTENT(IN OUT)::IDX
    INTEGER(KIND=INT64), INTENT(IN):: ISTART, ISTOP

    ! Local variables.

    REAL(KIND=REAL64):: AMIN  ! Temporary minimum.
    REAL(KIND=REAL64):: ATEMP  ! The value to be inserted.
    INTEGER(KIND=INT64):: I    ! Index variable.
    INTEGER(KIND=INT64):: IABOVE ! Index to find insertion point.
    INTEGER(KIND=INT64):: IMIN ! Temporary minimum position.
    INTEGER(KIND=INT64):: ITEMP ! Temporary for swapping.

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



END MODULE BALL_TREE
