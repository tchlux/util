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

  SUBROUTINE build_ibm(pts, sq_sums, box_widths)
    ! Given a set of "points", record the upper and lower widths of
    ! boxes in "structure". Assume that all boxes except the last
    ! were correctly sized before the addition of the last
    ! point. Reshape all boxes that need to be reshaped because of
    ! the single new addition to the end of "points".


    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: pts
    REAL(KIND=REAL64), DIMENSION(:),   INTENT(IN) :: sq_sums
    REAL(KIND=REAL64), DIMENSION(2*SIZE(pts,1), SIZE(pts,2)),&
         & INTENT(INOUT) :: box_widths
    ! Holder for the upper and lower widths surrounding each box
    REAL(KIND=REAL64), DIMENSION(SIZE(pts,2)) :: distances
    ! Distances between points
    INTEGER(KIND=INT64), DIMENSION(SIZE(pts,2)) :: distance_indices
    REAL(KIND=REAL64), DIMENSION(SIZE(pts,1)) :: center
    ! Indices of distances (for use in sorting)
    LOGICAL, DIMENSION(SIZE(pts,1)) :: in_lower, in_upper
    ! Holder for checking if a point is inside of a box
    INTEGER(KIND=INT64) :: dim, idx, b1, b2
    REAL(KIND=REAL64) :: ss
    ! Miscellaneous integers for stepping

    ! Start with the center-most box and build outwards with the
    ! iterative-box technique. For each new point:
    ! 
    !   - sort all existing boxes by their distance to that point
    !   - walk through existing boxes, closest to nearest
    !   - if a box contains the new point, bound both mutually on max dim
    ! 
    loop_over_new_boxes : DO b2 = 2, SIZE(pts,2)
       ! Measure the distance to all already-constructed boxes.
       ss = sq_sums(b2)

       !$omp parallel do private(b1)
       DO b1 = 1, b2-1
          distances(b1) = ss + sq_sums(b1) - 2 * DOT_PRODUCT(pts(:,b2), pts(:,b1))
          distance_indices(b1) = b1
       END DO
       !$omp end parallel do

       ! Sort the squared distances from each constructed box.
       CALL ARGSORT(distances(:b2-1), distance_indices(:b2-1))

       ! Find out if any existing boxes contain this new point.
       loop_over_finished_boxes : DO idx = 1, b2-1
          b1 = distance_indices(idx)

          ! Logical testing if the box_2 (new) is inside box_1 (old):
          !  -- lower width undefined .OR. (box - other) < lower width
          in_lower = (box_widths(:SIZE(pts,1),b1) .LT. 0.0_REAL64) .OR. &
               (pts(:,b1)-pts(:,b2) < box_widths(:SIZE(pts,1),b1))
          !  -- upper width undefined .OR. (other - box) < upper width
          in_upper = (box_widths(SIZE(pts,1)+1:,b1) .LT. 0.0_REAL64) .OR. &
               (box_widths(SIZE(pts,1)+1:,b1) > pts(:,b2)-pts(:,b1))
          ! If the already-constructed box (1) contains this new box (2),
          ! then reshape the old box not to contain the new box.
          rebuild_box_1 : IF (ALL(in_lower .AND. in_upper)) THEN
             ! Find the dimension of max magnitude (Linf) difference
             dim = MAXLOC(ABS(pts(:,b2) - pts(:,b1)),1)
             ! Check if we are adjusting the upper or lower width of
             ! the box, and adjust the boundary for this point
             IF (pts(dim,b1) > pts(dim,b2)) THEN
                ! Set the lower bound on the old (1)
                box_widths(dim,b1) = pts(dim,b1) - pts(dim,b2)
             ELSE
                ! Set the upper bound on the old (1)
                box_widths(dim+SIZE(pts,1),b1) = pts(dim,b2) - pts(dim,b1)
             END IF
          END IF rebuild_box_1

          ! Logical testing if the box_1 (old) is inside box_2 (new):
          !  -- lower width undefined .OR. (box - other) < lower width
          in_lower = (box_widths(:SIZE(pts,1),b2) .LT. 0.0_REAL64) .OR. &
               (pts(:,b2)-pts(:,b1) < box_widths(:SIZE(pts,1),b2))
          !  -- upper width undefined .OR. (other - box) < upper width
          in_upper = (box_widths(SIZE(pts,1)+1:,b2) .LT. 0.0_REAL64) .OR. &
               (box_widths(SIZE(pts,1)+1:,b2) > pts(:,b1)-pts(:,b2))
          ! If the new box (2) contains this already-constructed box (1),
          ! then reshape the new box not to include the old box.
          rebuild_box_2 : IF (ALL(in_lower .AND. in_upper)) THEN
             ! Find the dimension of max magnitude (Linf) difference
             dim = MAXLOC(ABS(pts(:,b1) - pts(:,b2)),1)
             ! Check if we are adjusting the upper or lower width of
             ! the box, and adjust the boundary for this point
             IF (pts(dim,b2) > pts(dim,b1)) THEN
                ! Set the lower bound on the new (2)
                box_widths(dim,b2) = pts(dim,b2) - pts(dim,b1)
             ELSE
                ! Set the upper bound on the new (2)
                box_widths(dim+SIZE(pts,1),b2) = pts(dim,b1) - pts(dim,b2)
             END IF
          END IF rebuild_box_2

       END DO loop_over_finished_boxes
    END DO loop_over_new_boxes
  END SUBROUTINE build_ibm


  !                    Box Mesh Computations                    
  !=============================================================

  SUBROUTINE eval_order_1_box_mesh(boxes, widths, x_points, box_vals)
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
    ! Variable for ensuring that weights don't get incorrectly dropped to zero.
    REAL (KIND=REAL64) :: small, ss
    SMALL = SQRT(EPSILON(small))
    ! ^^ local variables for computing the box evaluations

    !$omp parallel do private(pt_num,box_num,dim,shifted_box)
    pt_evals : DO pt_num=1, SIZE(x_points,2)
       ! Compute the linear box spline over all boxes.
       box_evals : DO box_num=1, SIZE(boxes,2)
          ! Calculate the normalized distance from center of this box
          normalize_box : DO dim = 1, SIZE(boxes,1)
             ! Use the lower width if the point is less than the center.
             IF (x_points(dim,pt_num) .LT. boxes(dim,box_num)) THEN
                shifted_box(dim) = (boxes(dim,box_num) - x_points(dim,pt_num)) &
                     / widths(dim,box_num)
             ! Use the upper width otherwise (equal or greater value).
             ELSE
                shifted_box(dim) = (x_points(dim,pt_num) - boxes(dim,box_num)) &
                     / widths(SIZE(boxes,1)+dim,box_num)
             END IF
          END DO normalize_box
          ! Compute the order 1 approximation and store the evaluation of
          ! this box for this x-data point 
          shifted_box(:) = MAX(1.0_REAL64 - shifted_box(:), 0.0_REAL64)
          box_vals(pt_num, box_num) = PRODUCT(shifted_box)
          ! If the product forced a nonzero weight box to zero, make it nonzero again.
          IF ((MINVAL(shifted_box) .GT. 0_REAL64) .AND. &
               (box_vals(pt_num,box_num) .LT. small)) &
               box_vals(pt_num, box_num) = small
       END DO box_evals
    END DO pt_evals
    !$omp end parallel do

    ! Make the weights convex.
    convexify: DO pt_num=1, SIZE(box_vals,1)
       ss = SUM(box_vals(pt_num,:))
       ! Make the weights sum to 1.
       IF (ss .GT. 0_REAL64) THEN
          box_vals(pt_num,:) = box_vals(pt_num,:) / ss
       ELSE
          ! If there was no coverage, average all values.
          ss = SIZE(box_vals,2)
          box_vals(pt_num,:) = 1_REAL64 / ss
       END IF
    END DO convexify

  END SUBROUTINE eval_order_1_box_mesh


  SUBROUTINE eval_order_2_box_mesh(boxes, widths, x_points, box_vals)
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
    REAL (KIND=REAL64) :: ss
    ! ^^ local variables for computing the box evaluations

    ! Compute the linear box spline for each cube
    box_evals : DO CONCURRENT (pt_num=1:SIZE(x_points,2), box_num=1:SIZE(boxes,2))
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

       ! ! Order 1 approximation
       ! box_vals(pt_num, box_num) = PRODUCT(&
       !      MAX(1.0_REAL64 - shifted_box, 0.0_REAL64))

       ! Order 2 approximation
       shifted_box = 1.5_REAL64 + 1.5_REAL64 * shifted_box
       WHERE ((1.0_REAL64 .LE. shifted_box) .AND. &
            (shifted_box .LT. 2.0_REAL64)) shifted_box = (&
            -(2*shifted_box**2 - 6*shifted_box + 3) )
       WHERE ((2.0_REAL64 .LE. shifted_box) .AND. &
            (shifted_box .LE. 3.0_REAL64)) shifted_box = (&
            (shifted_box - 3)**2)
       WHERE (shifted_box .GT. 3.0_REAL64) shifted_box = 0.0_REAL64
       box_vals(pt_num, box_num) = PRODUCT(shifted_box)
    END DO box_evals
    ! Make the weights convex.
    convexify: DO pt_num=1, SIZE(box_vals,1)
       ss = SUM(box_vals(pt_num,:))
       ! Make the weights sum to 1.
       IF (ss .GT. 0_REAL64) THEN
          box_vals(pt_num,:) = box_vals(pt_num,:) / ss
       ELSE
          ! If there was no coverate, average all values.
          box_vals(pt_num,:) = 1 / SIZE(box_vals,2)
       END IF
    END DO convexify
  END SUBROUTINE eval_order_2_box_mesh


  ! ------------------------------------------------------------------
  !                        FastSort method
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
  RECURSIVE SUBROUTINE ARGSORT(VALUES, INDICES, MIN_SIZE)
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
       CALL INSERTION_ARGSORT(VALUES, INDICES)
       ! Call this function recursively after pivoting about the median.
    ELSE
       ! ---------------------------------------------------------------
       ! If you are having slow runtime with the selection of pivot values 
       ! provided by ARGPARTITION, then consider using ARGSELECT instead.
       I = ARGPARTITION(VALUES, INDICES)
       ! ---------------------------------------------------------------
       ! I = SIZE(VALUES) / 2
       ! CALL ARGSELECT(VALUES, INDICES, I)
       ! ! Requires 'USE FAST_SELECT' at top of module.
       ! ---------------------------------------------------------------
       CALL ARGSORT(VALUES(:I-1), INDICES(:I-1), MS)
       CALL ARGSORT(VALUES(I+1:), INDICES(I+1:), MS)
    END IF
  END SUBROUTINE ARGSORT

  ! This function efficiently partitions values based on the median
  ! of the first, middle, and last elements of the VALUES array. This
  ! function returns the index of the pivot.
  FUNCTION ARGPARTITION(VALUES, INDICES) RESULT(LEFT)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:)            :: VALUES
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(SIZE(VALUES)) :: INDICES
    INTEGER(KIND=INT64) :: LEFT, MID, RIGHT
    REAL(KIND=REAL64)   :: PIVOT
    ! Use the median of the first, middle, and last element as the
    ! pivot. Place the pivot at the end of the array.
    MID = (1 + SIZE(VALUES)) / 2
    ! Swap the first and last elements (if the last is smaller).
    IF (VALUES(SIZE(VALUES)) < VALUES(1)) THEN
       CALL SWAPR64(VALUES(1),  VALUES(SIZE(VALUES)))
       CALL SWAPI64(INDICES(1), INDICES(SIZE(VALUES)))
    END IF
    ! Swap the middle and first elements (if the middle is smaller).
    IF (VALUES(MID) < VALUES(SIZE(VALUES))) THEN
       CALL SWAPR64(VALUES(MID),  VALUES(SIZE(VALUES)))
       CALL SWAPI64(INDICES(MID), INDICES(SIZE(VALUES)))       
       ! Swap the last and first elements (if the last is smaller).
       IF (VALUES(SIZE(VALUES)) < VALUES(1)) THEN
          CALL SWAPR64(VALUES(1),  VALUES(SIZE(VALUES)))
          CALL SWAPI64(INDICES(1), INDICES(SIZE(VALUES)))
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
       CALL SWAPR64(VALUES(LEFT),  VALUES(RIGHT))
       CALL SWAPI64(INDICES(LEFT), INDICES(RIGHT))
    END DO
    ! The last swap was done even though LEFT == RIGHT, we need to undo.
    CALL SWAPR64(VALUES(LEFT),  VALUES(RIGHT))
    CALL SWAPI64(INDICES(LEFT), INDICES(RIGHT))
    ! Finally, we put the pivot back into its proper location.
    CALL SWAPR64(VALUES(LEFT),  VALUES(SIZE(VALUES)))
    CALL SWAPI64(INDICES(LEFT), INDICES(SIZE(VALUES)))
  END FUNCTION ARGPARTITION

  ! Insertion sort (best for small lists).
  SUBROUTINE INSERTION_ARGSORT(VALUES, INDICES)
    REAL(KIND=REAL64),   INTENT(INOUT), DIMENSION(:)            :: VALUES
    INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(SIZE(VALUES)) :: INDICES
    ! Local variables.
    REAL(KIND=REAL64)   :: TEMP_VAL
    INTEGER(KIND=INT64) :: I, BEFORE, AFTER, TEMP_IND
    ! Return for the base case.
    IF (SIZE(VALUES) .LE. 1) RETURN
    ! Put the smallest value at the front of the list.
    I = MINLOC(VALUES,1)
    CALL SWAPR64(VALUES(1),  VALUES(I))
    CALL SWAPI64(INDICES(1), INDICES(I))
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
  END SUBROUTINE INSERTION_ARGSORT
  ! ------------------------------------------------------------------

  ! Routines for swapping values (cleans up code, will be optimized out).
  SUBROUTINE SWAPI64(V1, V2)
    INTEGER(KIND=INT64), INTENT(INOUT) :: V1, V2
    ! Local temp
    INTEGER(KIND=INT64) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAPI64

  SUBROUTINE SWAPR64(V1, V2)
    REAL(KIND=REAL64), INTENT(INOUT) :: V1, V2
    ! Local temp
    REAL(KIND=REAL64) :: TEMP
    TEMP = V1
    V1 = V2
    V2 = TEMP
  END SUBROUTINE SWAPR64


END MODULE meshes
