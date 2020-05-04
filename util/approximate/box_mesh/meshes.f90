! NAME:   Bootstrapped box splines
! AUTHOR: Thomas C.H. Lux
! EMAIL:  tchlux@vt.edu
! 
! DESCRIPTION: This file (box_spline_basis.f95) contains the module
!              box_spline_basis that declares the subroutines
!              (compute_boxes, evaluate_boxes) for generating and
!              evaluating a box spline hypercube basis model for
!              continuous multidimensional (d<=20) data with real
!              response values.
! 
! The following LAPACK routines are used:
!    + DGELS
! 

MODULE meshes
  USE ISO_FORTRAN_ENV, ONLY : INT64, REAL64
  IMPLICIT NONE

CONTAINS

  !==============================================================
  !                         Max Box Mesh                         
  !==============================================================

  SUBROUTINE build_mbm(pts, box_widths)
    ! Given a set of "points", record the upper and lower widths of
    ! boxes in "structure". Assume that all boxes except the last
    ! were correctly sized before the addition of the last
    ! point. Reshape all boxes that need to be reshaped because of
    ! the single new addition to the end of "points".

    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: pts
    REAL(KIND=REAL64), DIMENSION(2*SIZE(pts,1), SIZE(pts,2)),&
         & INTENT(OUT) :: box_widths
    ! Holder for the upper and lower widths surrounding each box
    REAL(KIND=REAL64), DIMENSION(SIZE(pts,2)) :: distances
    ! Distances between points
    INTEGER, DIMENSION(SIZE(pts,2)) :: distance_indices
    ! Indices of distances (for use in sorting)
    LOGICAL, DIMENSION(SIZE(pts,1)) :: in_lower, in_upper
    ! Holder for checking if a point is inside of a box
    INTEGER :: step, box, dim, idx, other_pt_idx, other_pt
    ! Miscellaneous integers for stepping

    ! Assume all boxes except the last were already built correctly,
    ! add a new box and reshape old boxes that may have been impacted
    compute_box_shapes : DO box = 1,SIZE(pts,2)
       other_pt_idx = SIZE(pts, 2)
       ! Logical testing if the new point is inside a box:
       !  -- lower width undefined .OR. (box - other) < lower width
       in_lower = (box_widths(:SIZE(pts,1),box) .LT. 0.0_REAL64) .OR. &
            (pts(:,box)-pts(:,other_pt_idx) < box_widths(:SIZE(pts,1),box))
       !  -- upper width undefined .OR. (other - box) < upper width
       in_upper = (box_widths(SIZE(pts,1)+1:,box) .LT. 0.0_REAL64) .OR. &
            (box_widths(SIZE(pts,1)+1:,box) > pts(:,other_pt_idx)-pts(:,box))
       ! If the box contains the new point, then rebuild the box
       rebuild_box : IF (ALL(in_lower .AND. in_upper)) THEN         
          ! We are rebuilding, so reset the box widths
          box_widths(:,box) = -1.0_REAL64
          ! Find the infinity norm distances to each other point
          linf_distances : DO CONCURRENT (idx = 1: SIZE(pts,2))
             distances(idx) = MAXVAL(ABS(pts(:,idx) - pts(:,box)))
          END DO linf_distances
          ! Sort the Linf distances between pts
          CALL QSORTC(distances, distance_indices)
          ! Identify the max-Linf-width bounding rectangle
          calculate_box_width : DO other_pt = 1, SIZE(pts,2)
             other_pt_idx = distance_indices(other_pt)
             ! Don't try to use the current box center as a boundary
             IF (other_pt_idx .EQ. box) CYCLE calculate_box_width
             ! Logical testing if a point is inside a box:
             !  -- lower width undefined .OR. (box - other) < lower width
             in_lower = (box_widths(:SIZE(pts,1),box) .LT. 0.0_REAL64) .OR. &
                  (pts(:,box)-pts(:,other_pt_idx) < box_widths(:SIZE(pts,1),box))
             !  -- upper width undefined .OR. (other - box) < upper width
             in_upper = (box_widths(SIZE(pts,1)+1:,box) .LT. 0.0_REAL64) .OR. &
                  (box_widths(SIZE(pts,1)+1:,box) > pts(:,other_pt_idx)-pts(:,box))
             ! If all conditions are satisfied, the point is in the box
             IF (ALL(in_lower .AND. in_upper)) THEN
                ! Find the dimension of max magnitude (Linf) difference
                dim = MAXLOC(ABS(pts(:,other_pt_idx) - pts(:,box)),1)
                ! Check if we are adjusting the upper or lower width of
                ! the box, and adjust the boundary for this point
                IF (pts(dim,box) > pts(dim,other_pt_idx)) THEN
                   box_widths(dim,box) = &
                        pts(dim,box) - pts(dim,other_pt_idx)
                ELSE
                   box_widths(dim+SIZE(pts,1),box) = &
                        pts(dim,other_pt_idx) - pts(dim,box)
                END IF
             END IF
          END DO calculate_box_width
       END IF rebuild_box
    END DO compute_box_shapes
  END SUBROUTINE build_mbm

  !                    Box Mesh Computations                    
  !=============================================================

  SUBROUTINE eval_box_mesh(boxes, widths, x_points, box_vals)
    ! Subroutine for computing the values of all the boxes at an x-point
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: boxes
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(SIZE(boxes,1)*2, &
         SIZE(boxes,2)) :: widths
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: x_points
    REAL (KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(x_points,2)&
         &,SIZE(boxes,2)) :: box_vals
    ! ^^ Holder for the evaluation of each box at the given x-point
    INTEGER :: pt_num, box_num, dim
    REAL (KIND=REAL64), DIMENSION(SIZE(boxes,1)) :: shifted_box
    ! ^^ local variables for computing the box evaluations

    ! Compute the linear box spline for each cube
    box_evals : DO CONCURRENT (pt_num=1:SIZE(x_points,2), &
         box_num=1:SIZE(boxes,2))
       ! Calculate the normalized distance from center of this box
       normalize_box : DO dim = 1, SIZE(boxes,1)
          IF (x_points(dim,pt_num) .LT. boxes(dim,box_num)) THEN
             shifted_box(dim) = &
                  (boxes(dim,box_num) - x_points(dim,pt_num)) / &
                  widths(dim,box_num)
          ELSE
             shifted_box(dim) = &
                  (x_points(dim,pt_num) - boxes(dim,box_num)) / &
                  widths(SIZE(boxes,1)+dim,box_num)
          END IF
       END DO normalize_box
       ! Store the evaluation of this box for this x-data point

       ! Order 1 approximation
       box_vals(pt_num, box_num) = PRODUCT(&
            MAX(1.0_REAL64 - shifted_box, 0.0_REAL64))

       ! ! Order 2 approximation
       ! shifted_box = 1.5_REAL64 + 1.5_REAL64 * shifted_box
       ! WHERE ((1.0_REAL64 .LE. shifted_box) .AND. &
       !      (shifted_box .LT. 2.0_REAL64)) shifted_box = (&
       !      -(2*shifted_box**2 - 6*shifted_box + 3) )
       ! WHERE ((2.0_REAL64 .LE. shifted_box) .AND. &
       !      (shifted_box .LE. 3.0_REAL64)) shifted_box = (&
       !      (shifted_box - 3)**2)
       ! WHERE (shifted_box .GT. 3.0_REAL64) shifted_box = 0.0_REAL64
       ! box_vals(pt_num, box_num) = PRODUCT(shifted_box)
    END DO box_evals
  END SUBROUTINE eval_box_mesh


  !                          Utilities                          
  !=============================================================

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
    DO CONCURRENT (I=1:SIZE(A)); IDX(I)=I; END DO
    
    CALL QSORTC_HELPER(A, IDX, 1, SIZE(A))
    RETURN

  CONTAINS
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
  END SUBROUTINE QSORTC

END MODULE meshes
