! NAME:   Voronoi Mesh
! AUTHOR: Thomas C.H. Lux
! EMAIL:  tchlux@vt.edu
! 

MODULE voronoi_mesh
  USE ISO_FORTRAN_ENV, ONLY : INT64, REAL64
  IMPLICIT NONE

CONTAINS
  SUBROUTINE train_vm(points, dots)
    ! Given the set of interpolation points, calculate the dot product
    ! of all pairs of points (because this will save time during evaluation)

    ! Input / Output variables
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: points
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2),SIZE(points,2)), &
         INTENT(OUT) :: dots
    ! Local variables
    INTEGER :: step1, step2

    ! Calculate all dot products between control points
    find_all_dots : DO step1 = 1, SIZE(points,2)
       DO step2 = 1, SIZE(points,2)
          IF (step1 .LE. step2) THEN
             dots(step1, step2) = SUM(points(:,step1)*points(:,step2))
          ELSE
             dots(step1, step2) = dots(step2, step1)
          END IF
       END DO
    END DO find_all_dots
  END SUBROUTINE train_vm

  SUBROUTINE predict_vm(control_points, control_dots, &
       approximation_points, approximation_weights)
    ! Given control points, pairwise dot products of control points,
    ! and approximation points, generate the coefficient matrix
    ! representing the convex combination of control points that can
    ! be used to predict each approximation point.

    ! Input / Output variables
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: control_points
    REAL(KIND=REAL64), DIMENSION(SIZE(control_points,2), &
         & SIZE(control_points,2)), INTENT(IN) :: control_dots
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: approximation_points
    REAL(KIND=REAL64), DIMENSION(SIZE(approximation_points,2),&
         SIZE(control_points,2)), INTENT(OUT) :: approximation_weights
    ! Local variables
    INTEGER :: step, step2
    REAL(KIND=REAL64), DIMENSION(SIZE(approximation_points,2), &
         SIZE(control_points,2)) :: point_dots

    ! Calculate all vectors between all points and all pairs of dot products
    find_all_dots : DO step = 1, SIZE(approximation_points,2)
       DO step2 = 1, SIZE(control_points,2)
          point_dots(step, step2) = &
               SUM(approximation_points(:,step)*control_points(:,step2))
       END DO
    END DO find_all_dots

    ! Evaluate the voronoi mesh at all approximation points
    CALL eval_vm(control_points, approximation_points, &
         control_dots, point_dots,  approximation_weights)
    ! ! Make sure the approximation weights are convex
    ! DO step = 1, SIZE(approximation_points,2)
    !    ! Make the weights convex
    !    approximation_weights(step,:) = approximation_weights(step,:) / &
    !         SUM(approximation_weights(step,:))
    ! END DO
  END SUBROUTINE predict_vm

  SUBROUTINE eval_vm(control_points, points, control_dots, &
       point_dots, mesh_values)
    ! Given the voronoi mesh control points and approximation
    ! points, calculate the basis function values for all cells and
    ! all points.

    ! Input / Output variables
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: control_points
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: points
    REAL(KIND=REAL64), DIMENSION(SIZE(control_points,2), &
         & SIZE(control_points,2)), INTENT(IN) :: control_dots
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2), &
         & SIZE(control_points,2)), INTENT(IN) :: point_dots
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2), &
         & SIZE(control_points,2)), INTENT(OUT) :: mesh_values

    ! Local variables
    INTEGER :: control_ind, point_ind, step
    LOGICAL, DIMENSION(SIZE(control_points,2)) :: ignore
    REAL(KIND=REAL64), DIMENSION(SIZE(control_points,2)) :: dists
    INTEGER, DIMENSION(SIZE(control_points,2)):: indices

    eval_points : DO point_ind = 1, SIZE(points, 2)
       ignore = .FALSE.
       ! Sort a list of all cells by increasing distance from current
       ! Sort a list of all cells by increasing distance from current
       DO control_ind = 1, SIZE(control_points,2)
          dists(control_ind) = SUM(points(:,point_ind)**2) - &
               2*point_dots(point_ind,control_ind) +&
               control_dots(control_ind,control_ind)
       END DO
       CALL QSORTC(dists, indices)
       ! Start with the closest cell, work outwards (the closest cells
       ! will rule out the most potential neighbors, saving time)
       eval_cells : DO step = 1, SIZE(control_points, 2)
          control_ind = indices(step)
          IF (ignore(control_ind)) THEN
             mesh_values(point_ind, control_ind) = 0
          ELSE
             mesh_values(point_ind, control_ind) = dist_to_boundary(&
                  control_points, points, control_ind, point_ind,&
                  ignore, control_dots, point_dots)
             mesh_values(point_ind, control_ind) = &
                  linear_basis_value(mesh_values(point_ind, control_ind))
          END IF
       END DO eval_cells
    END DO eval_points
  END SUBROUTINE eval_vm

  !                  Voronoi Mesh Computations                  
  !=============================================================

  FUNCTION dist_to_boundary(control_points, points, control_ind, &
       & point_ind, ignore, control_dots, point_dots)
    ! Measure the distance to the voronoi cell boundary along the
    ! vector (point - center) given all other control points. Ignore
    ! other control points that are 0-distance away.
    ! Input / Output variables
    INTEGER, INTENT(IN) :: control_ind, point_ind
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: control_points
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: points
    LOGICAL, DIMENSION(SIZE(control_points,2)), INTENT(INOUT) :: ignore
    REAL(KIND=REAL64), DIMENSION(SIZE(control_points,2), &
         & SIZE(control_points,2)), INTENT(IN) :: control_dots
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2), &
         & SIZE(control_points,2)), INTENT(IN) :: point_dots
    ! Local variables
    INTEGER :: step
    REAL(KIND=REAL64) :: dist_to_boundary, center_projected,&
         & bound_projected, point_projected, sign, ratio
    ! Initialization of local variables
    dist_to_boundary = 0.0_REAL64
    ! Find the closest boundary
    find_boundary : DO step = 1, SIZE(control_points,2)
       IF (step .EQ. control_ind) CYCLE find_boundary
       ! Project the center, bound, and point onto the boundary vector
       center_projected = control_dots(control_ind, step) - &
            control_dots(control_ind, control_ind)
       bound_projected = control_dots(step, step) - &
            control_dots(step, control_ind)
       point_projected = point_dots(point_ind, step) - &
            point_dots(point_ind, control_ind)
       ! Identify whether "point" is on the same side of "center" as "bound"
       sign = (bound_projected - center_projected) * &
            &(point_projected - center_projected)
       ! If "point" is on the same side of "center" as "bound"
       IF (sign .GT. 0) THEN
          ! Find the normalized distance to the boundary
          ratio = (point_projected - center_projected) / &
               (bound_projected - center_projected)
          ! If the boundary is the closest one, the store this boundary
          update_closest : IF (ratio > dist_to_boundary) THEN
             dist_to_boundary = ratio
          END IF update_closest
       ELSE IF (point_projected .EQ. bound_projected) THEN
          ! This means the point is on the boundary, which means the
          ! distance to a boundary is at least 1.0
          IF (dist_to_boundary .LT. 1) dist_to_boundary = 1.0
       ELSE IF (sign .LT. 0) THEN
          ignore(step) = .TRUE.
       END IF
    END DO find_boundary
  END FUNCTION dist_to_boundary

  !                          Utilities                          
  !=============================================================

  FUNCTION linear_basis_value(location)
    ! Compute the linear basis function value given a normalized
    ! location supporting strictly a range of (-1,1)
    REAL(KIND=REAL64), INTENT(IN) :: location
    REAL(KIND=REAL64) :: linear_basis_value
    linear_basis_value = MAX(1.0_REAL64 - ABS(location), 0.0_REAL64)
  END FUNCTION linear_basis_value

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

END MODULE voronoi_mesh

