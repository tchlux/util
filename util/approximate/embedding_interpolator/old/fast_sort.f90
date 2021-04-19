! TITLE:
!   Fast sort (and select)
!
! 
! PURPOSE:
!  This file contains the following modules and subroutines for quickly
!  sorting (either partially or fully) arrays in Fortran.
! 
!   SWAP
!    SWAP_I64
!    SWAP_R64
! 
!   FAST_SELECT
!    ARGSELECT_R64
!    ARGSELECT_I64
! 
!   FAST_SORT 
!    ARGSORT_R64
!    ARGPARTITION_R64
!    INSERTION_ARGSORT_R64
! 
! Respectively, the purposes of the three modules are:
! 
!   SWAP - Swap the value of two variables through a temporary.
! 
!   FAST_SELECT - Do a fast "select" operation on an array of numbers,
!     resulting in the value at index K being as if the entire array
!     were sorted. All values before K will be less or equal, but not
!     sorted. All values after K will be greater or equal, but not
!     sorted. This operation is O(N).
! 
!   FAST_SORT - Do a fast "sort" operation that uses methods
!     appropriate for the size of the array that will maximize speed.
! 
! 
! AUTHOR:
!   Thomas C.H. Lux (thomas.ch.lux@gmail.com)
! 
! 
! NOTES:
!   The algorithm for SELECT was inspired by the Floyd-Rivest method
!   of performing selection. The algorithm for SORT was inspired by
!   Data Structures and Algorithms Analysis by Dr. Clifford A. Shaffer
!   (Virginia Polytechnic Institute and State University).


! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MODULE SWAP
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64
  IMPLICIT NONE

CONTAINS
  SUBROUTINE SWAP_I64(V1, V2)
    INTEGER(KIND=INT64), INTENT(INOUT) :: V1, V2
    ! Local temp
    INTEGER(KIND=INT64) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAP_I64

  SUBROUTINE SWAP_R64(V1, V2)
    REAL(KIND=REAL64), INTENT(INOUT) :: V1, V2
    ! Local temp
    REAL(KIND=REAL64) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAP_R64
END MODULE SWAP
! ____________________________________________________________________



! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MODULE FAST_SELECT
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64
  USE SWAP, ONLY: SWAP_I64, SWAP_R64
  IMPLICIT NONE

CONTAINS
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
  RECURSIVE SUBROUTINE ARGSELECT_R64(VALUES, INDICES, K, DIVISOR, MAX_SIZE)
    ! Arguments
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:) :: VALUES
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
          CALL ARGSELECT_R64(VALUES(L:R), INDICES(L:R), K - L + 1, DIVISOR, MAX_SIZE)
       END IF
       ! Pick a partition element at position K.
       P = VALUES(K)
       L = LEFT
       R = RIGHT
       ! Move the partition element to the front of the list.
       CALL SWAP_R64(VALUES(LEFT), VALUES(K))
       CALL SWAP_I64(INDICES(LEFT), INDICES(K))
       ! Pre-swap the left and right elements (temporarily putting a
       ! larger element on the left) before starting the partition loop.
       IF (VALUES(RIGHT) .GT. P) THEN
          CALL SWAP_R64(VALUES(LEFT), VALUES(RIGHT))
          CALL SWAP_I64(INDICES(LEFT), INDICES(RIGHT))
       END IF
       ! Now partition the elements about the pivot value "P".
       DO WHILE (L .LT. R)
          CALL SWAP_R64(VALUES(L), VALUES(R))
          CALL SWAP_I64(INDICES(L), INDICES(R))
          L = L + 1
          R = R - 1
          DO WHILE (VALUES(L) .LT. P) ; L = L + 1 ; END DO
          DO WHILE (VALUES(R) .GT. P) ; R = R - 1 ; END DO
       END DO
       ! Place the pivot element back into its appropriate place.
       IF (VALUES(LEFT) .EQ. P) THEN
          CALL SWAP_R64(VALUES(LEFT), VALUES(R))
          CALL SWAP_I64(INDICES(LEFT), INDICES(R))
       ELSE
          R = R + 1
          CALL SWAP_R64(VALUES(R), VALUES(RIGHT))
          CALL SWAP_I64(INDICES(R), INDICES(RIGHT))
       END IF
       ! adjust left and right towards the boundaries of the subset
       ! containing the (k - left + 1)th smallest element
       IF (R .LE. K) LEFT = R + 1
       IF (K .LE. R) RIGHT = R - 1
    END DO
  END SUBROUTINE ARGSELECT_R64

  ! ------------------------------------------------------------------

  RECURSIVE SUBROUTINE ARGSELECT_I64(VALUES, INDICES, K, DIVISOR, MAX_SIZE)
    ! Arguments
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:) :: VALUES
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(:) :: INDICES
    INTEGER(KIND=INT64), INTENT(IN)                  :: K
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL        :: DIVISOR, MAX_SIZE
    ! Locals
    INTEGER(KIND=INT64) :: LEFT, RIGHT, L, R, MS, D, P
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
          CALL ARGSELECT_I64(VALUES(L:R), INDICES(L:R), K - L + 1, DIVISOR, MAX_SIZE)
       END IF
       ! Pick a partition element at position K.
       P = VALUES(K)
       L = LEFT
       R = RIGHT
       ! Move the partition element to the front of the list.
       CALL SWAP_I64(VALUES(LEFT), VALUES(K))
       CALL SWAP_I64(INDICES(LEFT), INDICES(K))
       ! Pre-swap the left and right elements (temporarily putting a
       ! larger element on the left) before starting the partition loop.
       IF (VALUES(RIGHT) .GT. P) THEN
          CALL SWAP_I64(VALUES(LEFT), VALUES(RIGHT))
          CALL SWAP_I64(INDICES(LEFT), INDICES(RIGHT))
       END IF
       ! Now partition the elements about the pivot value "T".
       DO WHILE (L .LT. R)
          CALL SWAP_I64(VALUES(L), VALUES(R))
          CALL SWAP_I64(INDICES(L), INDICES(R))
          L = L + 1
          R = R - 1
          DO WHILE (VALUES(L) .LT. P) ; L = L + 1 ; END DO
          DO WHILE (VALUES(R) .GT. P) ; R = R - 1 ; END DO
       END DO
       ! Place the pivot element back into its appropriate place.
       IF (VALUES(LEFT) .EQ. P) THEN
          CALL SWAP_I64(VALUES(LEFT), VALUES(R))
          CALL SWAP_I64(INDICES(LEFT), INDICES(R))
       ELSE
          R = R + 1
          CALL SWAP_I64(VALUES(R), VALUES(RIGHT))
          CALL SWAP_I64(INDICES(R), INDICES(RIGHT))
       END IF
       ! adjust left and right towards the boundaries of the subset
       ! containing the (k - left + 1)th smallest element
       IF (R .LE. K) LEFT = R + 1
       IF (K .LE. R) RIGHT = R - 1
    END DO
  END SUBROUTINE ARGSELECT_I64
END MODULE FAST_SELECT
! ____________________________________________________________________



! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MODULE FAST_SORT
  USE ISO_FORTRAN_ENV, ONLY: REAL64, INT64
  USE SWAP, ONLY: SWAP_I64, SWAP_R64
  IMPLICIT NONE

CONTAINS
  ! ------------------------------------------------------------------
  !                         FastSort
  ! 
  ! This routine uses a combination of QuickSort (with modestly
  ! intelligent pivot selection) and Insertion Sort (for small arrays)
  ! to achieve very fast average case sort times for both random and
  ! partially sorted data. The pivot is selected for QuickSort as the
  ! median of the first, middle, and last values in the array.
  ! 
  ! Arguments:
  ! 
  !   VALUES   --  A 1D array of real numbers.
  !   INDICES  --  A 1D array of original indices for elements of VALUES.
  ! 
  ! Optional:
  ! 
  !   MIN_SIZE --  An positive integer that represents the largest
  !                sized VALUES for which a partition about a pivot
  !                is used to reduce the size of a an unsorted array.
  !                Any size less than this will result in the use of
  !                INSERTION_ARGSORT instead of ARGPARTITION.
  ! 
  ! Output:
  ! 
  !   The elements of the array VALUES are sorted and all elements of
  !   INDICES are sorted symmetrically (given INDICES = 1, ...,
  !   SIZE(VALUES) beforehand, final INDICES will show original index
  !   of each element of VALUES before the sort operation).
  ! 
  RECURSIVE SUBROUTINE ARGSORT_R64(VALUES, INDICES, MIN_SIZE)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:)            :: VALUES
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(SIZE(VALUES)) :: INDICES
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL                   :: MIN_SIZE
    ! Local variables
    INTEGER(KIND=INT64) :: I, MS
    IF (PRESENT(MIN_SIZE)) THEN ; MS = MIN_SIZE
    ELSE                        ; MS = 2**6
    END IF
    ! Base case, return.
    IF (SIZE(VALUES) .LT. MS) THEN
       CALL INSERTION_ARGSORT_R64(VALUES, INDICES)
       ! Call this function recursively after pivoting about the median.
    ELSE
       ! ---------------------------------------------------------------
       ! If you are having slow runtime with the selection of pivot values 
       ! provided by ARGPARTITION, then consider using ARGSELECT instead.
       I = ARGPARTITION_R64(VALUES, INDICES)
       ! ---------------------------------------------------------------
       ! I = SIZE(VALUES) / 2
       ! CALL ARGSELECT_R64(VALUES, INDICES, I)
       ! ! Requires 'USE FAST_SELECT' at top of subroutine or module.
       ! ---------------------------------------------------------------
       CALL ARGSORT_R64(VALUES(:I-1), INDICES(:I-1), MS)
       CALL ARGSORT_R64(VALUES(I+1:), INDICES(I+1:), MS)
    END IF
  END SUBROUTINE ARGSORT_R64

  ! This function efficiently partitions values based on the median
  ! of the first, middle, and last elements of the VALUES array. This
  ! function returns the index of the pivot.
  FUNCTION ARGPARTITION_R64(VALUES, INDICES) RESULT(LEFT)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:)            :: VALUES
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(SIZE(VALUES)) :: INDICES
    INTEGER(KIND=INT64) :: LEFT, MID, RIGHT
    REAL(KIND=REAL64)   :: PIVOT
    ! Use the median of the first, middle, and last element as the
    ! pivot. Place the pivot at the end of the array.
    MID = (1 + SIZE(VALUES)) / 2
    ! Swap the first and last elements (if the last is smaller).
    IF (VALUES(SIZE(VALUES)) < VALUES(1)) THEN
       CALL SWAP_R64(VALUES(1),  VALUES(SIZE(VALUES)))
       CALL SWAP_I64(INDICES(1), INDICES(SIZE(VALUES)))
    END IF
    ! Swap the middle and first elements (if the middle is smaller).
    IF (VALUES(MID) < VALUES(SIZE(VALUES))) THEN
       CALL SWAP_R64(VALUES(MID),  VALUES(SIZE(VALUES)))
       CALL SWAP_I64(INDICES(MID), INDICES(SIZE(VALUES)))       
       ! Swap the last and first elements (if the last is smaller).
       IF (VALUES(SIZE(VALUES)) < VALUES(1)) THEN
          CALL SWAP_R64(VALUES(1),  VALUES(SIZE(VALUES)))
          CALL SWAP_I64(INDICES(1), INDICES(SIZE(VALUES)))
       END IF
    END IF
    ! Set the pivot, LEFT index and RIGHT index (skip the smallest,
    ! which is in location 1, and the pivot at the end).
    PIVOT = VALUES(SIZE(VALUES))
    LEFT  = 2
    RIGHT = SIZE(VALUES) - 1
    ! Partition all elements to the left and right side of the pivot
    ! (left if they are smaller, right if they are bigger).
    DO WHILE (LEFT < RIGHT)
       ! Loop left until we find a value that is greater or equal to pivot.
       DO WHILE (VALUES(LEFT) < PIVOT)
          LEFT = LEFT + 1
       END DO
       ! Loop right until we find a value that is less or equal to pivot (or LEFT).
       DO WHILE (RIGHT .NE. LEFT)
          IF (VALUES(RIGHT) .LT. PIVOT) EXIT
          RIGHT = RIGHT - 1
       END DO
       ! Now we know that [VALUES(RIGHT) < PIVOT < VALUES(LEFT)], so swap them.
       CALL SWAP_R64(VALUES(LEFT),  VALUES(RIGHT))
       CALL SWAP_I64(INDICES(LEFT), INDICES(RIGHT))
    END DO
    ! The last swap was done even though LEFT == RIGHT, we need to undo.
    CALL SWAP_R64(VALUES(LEFT),  VALUES(RIGHT))
    CALL SWAP_I64(INDICES(LEFT), INDICES(RIGHT))
    ! Finally, we put the pivot back into its proper location.
    CALL SWAP_R64(VALUES(LEFT),  VALUES(SIZE(VALUES)))
    CALL SWAP_I64(INDICES(LEFT), INDICES(SIZE(VALUES)))
  END FUNCTION ARGPARTITION_R64

  ! Insertion sort (best for small lists).
  SUBROUTINE INSERTION_ARGSORT_R64(VALUES, INDICES)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:)            :: VALUES
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(SIZE(VALUES)) :: INDICES
    ! Local variables.
    REAL(KIND=REAL64)   :: TEMP_VAL
    INTEGER(KIND=INT64) :: I, BEFORE, AFTER, TEMP_IND
    ! Return for the base case.
    IF (SIZE(VALUES) .LE. 1) RETURN
    ! Put the smallest value at the front of the list.
    I = MINLOC(VALUES,1)
    CALL SWAP_R64(VALUES(1),  VALUES(I))
    CALL SWAP_I64(INDICES(1), INDICES(I))
    ! Insertion sort the rest of the array.
    DO I = 3, SIZE(VALUES)
       TEMP_VAL = VALUES(I)
       TEMP_IND = INDICES(I)
       ! Search backwards in the list, 
       BEFORE = I - 1
       AFTER  = I
       DO WHILE (TEMP_VAL .LT. VALUES(BEFORE))
          VALUES(AFTER)  = VALUES(BEFORE)
          INDICES(AFTER) = INDICES(BEFORE)
          BEFORE = BEFORE - 1
          AFTER  = AFTER - 1
       END DO
       ! Put the value into its place (where it is greater than the
       ! element before it, but less than all values after it).
       VALUES(AFTER)  = TEMP_VAL
       INDICES(AFTER) = TEMP_IND
    END DO
  END SUBROUTINE INSERTION_ARGSORT_R64
  ! ------------------------------------------------------------------

END MODULE FAST_SORT
