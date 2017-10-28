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

MODULE bootstrapped_box_splines
  USE ISO_FORTRAN_ENV, ONLY : INT64, REAL64
  IMPLICIT NONE

CONTAINS
  SUBROUTINE compute_boxes(x_data, response, boxes, widths, weights, batch_size)
    ! This is a serial implementation of the hypercube box-spline
    ! basis approximation algorithm.
    ! 
    ! compute_boxes(x_data, response, boxes, widths, weights)
    !   x_data   -- matrix of n-dimensions by m-data points (nxm) stored as
    !               column vectors. x_data(:,1) should be the first point.
    !   response -- array of length m holding response values at each x-point
    !   boxes    -- matrix of n-dimensions by p-boxes (n x p), will be
    !               overwritten with the output boxes
    !   widths   -- array of length p holding the width of each box
    !   weights  -- array of length p holding the weight for each box function
    
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: x_data
    ! The matrix of data points (stored in rows as x-vectors)

    REAL (KIND=REAL64), INTENT(IN), DIMENSION(SIZE(x_data,2)) :: response
    ! The vector of response values associated with each point

    REAL (KIND=REAL64), INTENT(INOUT), DIMENSION(:,:) :: boxes
    ! Holder for the center coordinates of each box-spline support region

    REAL (KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(boxes,2)) :: widths
    ! Holder for the width of each box-spline support region

    REAL (KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(boxes,2)) :: weights
    ! Holder for the multiplier weight of each box-spline basis
    !   (solved for with LAPACK least squares)

    INTEGER (KIND=INT64), INTENT(IN), OPTIONAL :: batch_size
    ! Number of boxes to add per iteration

    !     Local workspace variables     
    !===================================

    INTEGER (KIND=INT64) :: local_batch_size
    ! Number of boxes to add per iteration

    REAL (KIND=REAL64), DIMENSION(:), ALLOCATABLE :: dgels_work_array
    ! Self-explanatory name

    REAL (KIND=REAL64), DIMENSION(SIZE(boxes,2),SIZE(x_data,2)) :: coef_matrix
    ! The matrix of coefficients in the least squares solve, holds
    ! the evaluation of all box-splines at each x data point

    REAL (KIND=REAL64), DIMENSION(SIZE(x_data,2)) :: local_response
    ! Local copy of 'response' so that the original is not tampered with

    REAL (KIND=REAL64), DIMENSION((SIZE(x_data,2)*(SIZE(x_data,2)-1))/2) :: distances
    ! Holder for the pairwise distances between data points (m*(m-1))/2

    INTEGER, DIMENSION((SIZE(x_data,2)*(SIZE(x_data,2)-1))/2) :: distance_indices
    ! Holder for the original un-sorted indices of pairwise distances

    REAL(KIND=REAL64), DIMENSION(2*SIZE(x_data,1),SIZE(boxes,2)) :: box_widths
    ! Storage for hyperbox dimensions before they are converted to
    ! hypercube, box_widths(:SIZE(x_data,1),:) holds the lower
    ! widths, box_widths(SIZE(x_data,1)+1:,:) holds the upper widths

    LOGICAL, DIMENSION(SIZE(x_data,2)) :: is_center_point
    ! Holder for which data points are box centers

    INTEGER :: num_box_centers, step, idx, first, second, dim, info
    ! num_box_centers -- Number of box centers, SUM(IF is_center_point 1)
    ! step -- Iterator used for stepping through various algorithms
    ! idx -- An index in the list of pairwise distances or boxes
    ! first -- The index of the first point in a pair
    ! second -- The index of the second point in a pair
    ! dim -- Index of which dimension is being checked for a data point
    ! info -- Used to send / receive messages from DGELS

    REAL (KIND=REAL64) :: random_num
    ! random_num -- Holder for random numbers in the range [0,1]

    INTEGER :: batch, new_boxes
    ! batch -- Current batch number
    ! new_boxes -- Counter for the number of new boxes that have been added
    INTEGER, DIMENSION(SIZE(boxes,2)) :: batch_errors
    ! Array of error achieved for each batch
    INTEGER, DIMENSION(SIZE(response,1)) :: error_indices
    ! Indices for sorting errors
    REAL(KIND=REAL64), DIMENSION(SIZE(response,1)) :: errors
    ! Container for actual errors
    INTEGER, DIMENSION(SIZE(boxes,2)) :: box_indices

    !     Initialization     
    !========================
    ! Set the batch size according to whether or not it was given
    IF (PRESENT(batch_size)) THEN
       local_batch_size = batch_size
    ELSE
       local_batch_size = 1
    END IF

    PRINT *, "SIZE(boxes,2): ", SIZE(boxes,2)
    PRINT *, "local_batch_size: ", local_batch_size
    PRINT *, "SIZE(x_data,1): ", SIZE(x_data,1)
    PRINT *, "SIZE(x_data,2): ", SIZE(x_data,2)
    PRINT *, "SIZE(response): ", SIZE(response)
    PRINT *, "bootstrapped_box_splines.f90 Line 119: "

    ! Initialization of local arrays
    errors = response
    local_response = response
    num_box_centers = SIZE(x_data,2)
    is_center_point = .TRUE.
    ! Initialize the list of indices of box-centers
    DO CONCURRENT (idx = 1: SIZE(box_indices,1))
       box_indices(idx) = idx
    END DO

    box_initialization : IF (local_batch_size > 1) THEN
       !     Step 1: Identify initial box centers     
       !=============================================
       !     n = local_batch_size ! "number of boxes desired"
       !     m = SIZE(x_data,2) ! "number of x-data points"
       !     d = SIZE(x_data,1) ! "dimension of the x-data"
       !     Computation overview:
       !         d m (m-1) / 2    ! pairwise L2 distance between points
       !       + m log_2(m)       ! sorting the pairwise distances
       !       + n                ! picking the box centers to keep
       
       ! Pairwise L2 distance between points
       step = 1
       pairwise_outer_loop : DO first = 1, SIZE(x_data,2)-1
          pairwise_inner_loop : DO second = first+1, SIZE(x_data,2)
             distances(step) = SUM((x_data(:,first) - x_data(:,second))**2)
             step = step + 1
          END DO pairwise_inner_loop
       END DO pairwise_outer_loop

       ! Sorting the pairwise distances
       CALL QSORTC(distances, distance_indices)

       PRINT *, "num_box_centers: ", num_box_centers
       ! Picking the box centers to keep (well spaced, QNSTOP algorithm)
       step = 0
       ! Starting with the closest pair and working to more spaced pairs,
       ! randomly select points from each pair to be dropped until there
       ! are only "n" points left to use as the centers of boxes.
       identify_centers : DO WHILE (num_box_centers > local_batch_size)
          step = step + 1
          ! Identify the two points that compose the next closest pair
          idx = distance_indices(step)
          ! Map (pair index) -> (first point index, second point index)
          CALL idx_to_pair(idx, SIZE(x_data,2), first, second)
          ! If both are still in the 'to keep' set, then remove one randomly
          IF (is_center_point(first) .AND. is_center_point(second)) THEN
             CALL RANDOM_NUMBER(random_num)
             IF (random_num > 0.5) THEN
                is_center_point(first) = .FALSE.
             ELSE
                is_center_point(second) = .FALSE.
             END IF
             num_box_centers = num_box_centers - 1
          END IF
       END DO identify_centers

    ELSE
       ! WARNING: Re-using memory, this is actually sorting response values
       CALL QSORTC(errors, error_indices)
       ! Set all center points to false
       is_center_point = .False.
       ! Make the median response valued x-point the initial box center.
       is_center_point(error_indices( SIZE(error_indices,1)/2 )) = .TRUE.
       num_box_centers = 1
    END IF box_initialization

    PRINT *, "bootstrapped_box_splines.f90 Line 179: BOX CENTERS"
    ! Initialize all box centers basd on the "is_center_point" array
    idx = 1
    step = 1
    initialize_boxes : DO step=1, SIZE(is_center_point,1)
       IF (is_center_point(step)) THEN
          boxes(:,idx) = x_data(:,step)
          PRINT *, "boxes(:,",idx,"): ", boxes(:,idx)
          idx = idx + 1
       END IF
    END DO initialize_boxes

    !     Step 2: Calculate box-spline hypercube widths     
    !=======================================================
    !     n = SIZE(boxes,2)  ! "number of boxes desired"
    !     d = SIZE(x_data,1) ! "number of x-data points"
    !     Computation overview:
    !       n *                ! "for each box-center"
    !      (  n d              ! "calculate Linf distance to other centers" 
    !       + n log_2(n)       ! "sort other boxes by Linf distances"
    !       + n d        )     ! "identify largest box width in each 
    !                          !  dimension allowed without overlap"
    !       + n d              ! "identify smallest cube containing
    !                             the max-Linf bounding rectangle"

    info = 0
    ! Call LAPACK DGELS to query the appropriate size for the work array
    CALL DGELS('T', SIZE(coef_matrix,1), SIZE(coef_matrix,2), 1, &
         coef_matrix(:,:), SIZE(coef_matrix,1), &
         local_response, SIZE(local_response,1), local_response, &
         -1, info)
    ! Check for errors in identifying amount of storage space necessary
    IF (info .EQ. 0) THEN
       dim = local_response(1)
    ELSE
       PRINT *, 'ERROR: unsuccessful query of work-array size for DGELS.'
       dim = 2 * num_box_centers * SIZE(x_data,2)
    END IF
    ! Allocate space for the working array
    ALLOCATE( dgels_work_array(1:dim) )

    ! Primary loop for performing the bootstrapping
    ! TODO: Update runtime analysis of bootstrapping process
    bootstrapping_main_loop : DO batch = 1, SIZE(batch_errors)
       ! Calculate the shapes of the max-boxes
       CALL max_boxes(boxes(:,:num_box_centers), &
            box_indices(:num_box_centers), &
            box_widths(:,:num_box_centers))
       ! Use the smallest bounding hyper-cube for each box size
       widths(:num_box_centers) = MAXVAL(box_widths(:,:num_box_centers),1)

       PRINT *, "bootstrapped_box_splines.f90 Line 233: "
       PRINT *, "batch: ", batch
       !     Step 3: Perform least-squares fit of data     
       !===================================================
       !     d = SIZE(x_data,1) ! "number of dimensions for x-data"
       !     m = SIZE(x_data,2) ! "number of x-data points"
       !     n = SIZE(boxes,2)  ! "number of boxes desired"
       !     A = [ m x n ] matrix ! "coefficient matrix of box spline
       !                             function values at all x-data points"
       !     y = [ m ] vector     ! "true response values at each x "
       !     Computation overview:
       !       n d m              ! "computing the box-splines at datum"
       !       + n^2 m            ! "computing A^t A"
       !       + n m              ! "computing A^t y"
       !       + n^3              ! "computing LU decomposition of A"

       ! Evaluate all box-spline values at all of the x-data points,
       ! generate the "coef_matrix" to be used for least squares
       CALL eval_box_coefs(boxes(:,:num_box_centers), widths(:num_box_centers), &
            &x_data, coef_matrix(:,:num_box_centers))

       PRINT *, "bootstrapped_box_splines.f90 Line 254: "
       ! Copy in the true respose values for the DGELS solve
       info = 0
       local_response = response
       ! Call the LAPACK routine for doing least squares minimization of Ax = B
       CALL DGELS('T', SIZE(coef_matrix,1), num_box_centers, 1, &
            coef_matrix(:,:num_box_centers), SIZE(coef_matrix,1), &
            local_response, SIZE(local_response,1), &
            dgels_work_array, dim, info)
       !   'T'                 -- Transpose of A necessary
       !   SIZE(coef_matrix,1) -- number of rows in A
       !   SIZE(coef_matrix,2) -- number of columns in A
       !   1                   -- Number of right hand sides
       !   coef_matrix         -- A
       !   SIZE(coef_matrix,1) -- Leading dimension of A (number of rows)
       !   local_response      -- B (will be overwritten to hold X after call)
       !   SIZE(response,1)    -- Leading dimension of B (number of rows)
       !   dgels_work_array    -- Workspace array for DGELS
       !   dim                 -- Size of dgels_work_array
       !   info                -- For verifying successful execution

       PRINT *, "bootstrapped_box_splines.f90 Line 275: "
       ! Extract the weights from the output array for DGELS
       weights = local_response(1:num_box_centers)

       !     Step 4: Evaluate current model error     
       !==============================================

       PRINT *, "bootstrapped_box_splines.f90 Line 282: "
       CALL eval_boxes(boxes(:,:num_box_centers), widths(:num_box_centers), &
            weights(:num_box_centers), x_data, errors)
       errors = - (errors - response)**2
       CALL QSORTC(errors, error_indices)
       ! Now we have the indices (sorted from least to largest
       ! magnitude to smallest magnitude)

       batch_errors(batch) = -SUM(errors)

       ! TODO: Stopping condition based on change in error
       ! IF (should_stop(batch_errors(:batch))) THEN
       !    EXIT bootstrapping_main_loop
       ! END IF

       !     Step 5: Add new boxes where error is large     
       !====================================================

       PRINT *, "bootstrapped_box_splines.f90 Line 300: "
       idx = 1
       box_add_loop : DO new_boxes = 1, local_batch_size
          ! If we cannot add another box, the break
          IF (num_box_centers + new_boxes .GT. SIZE(boxes,2)) THEN
             ! Compensate for the false "new_box" that was never added
             num_box_centers = num_box_centers - 1
             EXIT bootstrapping_main_loop
          END IF
          ! Find suitable new error point to add
          find_suitable_new_point : DO
             IF (idx > SIZE(errors)) THEN
                ! Compensate for the false "new_box" that was never added
                num_box_centers = num_box_centers - 1
                ! If we have tried all possible points, then done
                EXIT box_add_loop
             ELSE IF (boxes_contain_point(&
                  boxes(:,num_box_centers+1:num_box_centers+new_boxes), &
                  widths(num_box_centers+1:num_box_centers+new_boxes), &
                  x_data(:,idx))) THEN
                ! While the current idx error point is in an already
                ! added box, keep searching
                idx = idx + 1
             ELSE
                ! If we have found a new error point not inside any
                ! new boxes, 
                EXIT find_suitable_new_point
             END IF
          END DO find_suitable_new_point

          ! Add new box to the list of box centers
          num_box_centers = num_box_centers + new_boxes
          boxes(:,num_box_centers) = x_data(:,idx)

          ! Re-calculate the box width of *just the new box*
          CALL max_boxes(boxes(:,:num_box_centers), &
               box_indices(num_box_centers:num_box_centers), &
               box_widths(:,:num_box_centers))
          ! Use the smallest bounding hyper-cube for this new box size
          widths(num_box_centers:num_box_centers) = &
               MAXVAL(box_widths(:,num_box_centers:num_box_centers),1)
          ! Revert the number of box centers to be that minus new boxes
          num_box_centers = num_box_centers - new_boxes
       END DO box_add_loop

       PRINT *, "bootstrapped_box_splines.f90 Line 347: "
       ! Recalculate the shapes of the current max-boxes
       CALL max_boxes(boxes(:,:num_box_centers), &
            box_indices(:num_box_centers), &
            box_widths(:,:num_box_centers))
       ! Use the smallest bounding hyper-cube for each box size
       widths(:num_box_centers) = MAXVAL(box_widths(:,:num_box_centers),1)

       PRINT *, "bootstrapped_box_splines.f90 Line 355: "
    END DO bootstrapping_main_loop

    PRINT *, "bootstrapped_box_splines.f90 Line 358: "
    !     Final Step: Recalculate Full Box Mesh and Fit     
    !=======================================================
    
    ! Evaluate all box-spline values at all of the x-data points,
    ! generate the "coef_matrix" to be used for least squares
    CALL eval_box_coefs(boxes, widths, x_data, coef_matrix)
    ! Copy in the true respose values for the DGELS solve
    info = 0
    local_response = response
    ! Call the LAPACK routine for doing least squares minimization of Ax = B
    CALL DGELS('T', SIZE(coef_matrix,1), num_box_centers, 1, &
         coef_matrix(:,:num_box_centers), SIZE(coef_matrix,1), &
         local_response, SIZE(local_response,1), &
         dgels_work_array, dim, info)
    ! Extract the weights from the output array for DGELS
    weights = local_response(1:num_box_centers)

    PRINT *, "bootstrapped_box_splines.f90 Line 373: "

  CONTAINS

    ! Function for evaluating whether the bootstrapping algorithm
    ! should stop based on the progression of errors seen
    FUNCTION should_stop(errors)
      LOGICAL :: should_stop
      REAL(KIND=REAL64), DIMENSION(:) :: errors
      should_stop = ( ABS(errors(1) - errors(SIZE(errors,1))) &
           / errors(1) ) .GT. 0.9_REAL64
    END FUNCTION should_stop

    FUNCTION boxes_contain_point(boxes, widths, pt)
      LOGICAL :: boxes_contain_point
      REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: boxes
      REAL(KIND=REAL64), INTENT(IN), DIMENSION(SIZE(boxes,2)) :: widths
      REAL(KIND=REAL64), INTENT(IN), DIMENSION(SIZE(boxes,1)) :: pt
      INTEGER :: idx
      ! Initialize to assuming that no boxes contain the point
      boxes_contain_point = .FALSE.
      ! Search through all of the given boxes
      contains_search : DO idx = 1, SIZE(boxes,2)
         IF (ALL(ABS(boxes(:,idx) - pt) < widths(idx))) THEN
            ! If any box contains the point, set true and break
            boxes_contain_point = .TRUE.
            EXIT contains_search
         END IF
      END DO contains_search
    END FUNCTION boxes_contain_point

    ! Function for calculating the index number of a point based on
    ! the index of the pair (assuming pairs are [1,2], [1,3]... [2,3]...)
    SUBROUTINE idx_to_pair(idx, num_individuals, first, second)
      INTEGER, INTENT(IN) :: idx, num_individuals
      INTEGER, INTENT(OUT) :: first, second
      REAL (KIND=REAL64) :: pair_num, i
      ! pair_num -- The number whose integer part represents the index
      !             of the first point and fractional part the index of
      !             second point in "distance_indices"
      ! i        -- holds the converted index (for the math to work)

      ! First figure out how many pairs there will be, convert to
      ! backwards index (required to use functional form)
      i = ( num_individuals * (num_individuals - 1) ) / 2 - idx
      ! Get the pair number (using function whos whole + remainder form
      ! the sequence [1,2], [1,3], ... [1,num_individuals], [2,3], ...)
      pair_num = (1.0_REAL64 + SQRT(1.0_REAL64 + 8.0_REAL64 * i)) / 2.0_REAL64
       ! 'first' is determined by the whole part of "pair_num"
      first = FLOOR(pair_num)
      ! 'second' is determined by the (rounded) remainder part of "pair_num"
      second = (pair_num - first) * first + 0.5_REAL64
      ! Shift first and second to start at the front of the list
      second = num_individuals - second
      first = num_individuals - first
    END SUBROUTINE idx_to_pair

  END SUBROUTINE compute_boxes


  ! Subroutine for computing the max-box around given points
  !  where the max-box is the box defined by a center, lower, and
  !  upper widths such that MIN(lower,upper) is maximized
  SUBROUTINE max_boxes(pts, box_centers, box_widths)
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: pts
    ! Column vectors of points
    INTEGER, INTENT(IN), DIMENSION(:) :: box_centers
    ! Box center points (indices of points in "pts")
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(2*SIZE(box_centers,1),&
         SIZE(box_centers,1)) :: box_widths
    ! Holder for the upper and lower widths surrounding each box
    REAL(KIND=REAL64), DIMENSION(SIZE(pts,2)) :: distances
    ! Distances between points
    INTEGER, DIMENSION(SIZE(pts,2)) :: distance_indices
    ! Indices of distances (for use in sorting)
    LOGICAL, DIMENSION(SIZE(pts,1)) :: in_lower, in_upper
    ! Holder for checking if a point is inside of a box
    INTEGER :: step, box, dim, idx, other_pt_idx, other_pt
    ! Miscellaneous integers for stepping

    ! Initialize all box widths to be 'undefined'
    ! box_widths = HUGE(box_widths(1,1))
    box_widths= -1.0_REAL64
    ! Iterate over box centers (points)
    compute_box_shapes : DO box = 1,SIZE(box_centers,1)
       step = box_centers(box)
       ! Find the infinity norm distances to each other point
       linf_distances : DO CONCURRENT (idx = 1: SIZE(pts,2))
          distances(idx) = MAXVAL(ABS(pts(:,idx) - pts(:,step)))
       END DO linf_distances
       ! Sort the Linf distances between pts
       CALL QSORTC(distances, distance_indices)
       ! Identify the max-Linf-width bounding rectangle
       calculate_box_width : DO other_pt = 1, SIZE(pts,2)
          other_pt_idx = distance_indices(other_pt)
          ! Don't try to use the current box center as a boundary
          IF (other_pt_idx .EQ. step) CYCLE
          ! Logical testing if a point is inside a box:
          !  -- lower width undefined .OR. (box - other) < lower width
          in_lower = (box_widths(:SIZE(pts,1),step) .LT. 0.0_REAL64) .OR. &
               (pts(:,step)-pts(:,other_pt_idx) < box_widths(:SIZE(pts,1),step))
          !  -- upper width undefined .OR. (other - box) < upper width
          in_upper = (box_widths(SIZE(pts,1)+1:,step) .LT. 0.0_REAL64) .OR. &
               (box_widths(SIZE(pts,1)+1:,step) > pts(:,other_pt_idx)-pts(:,step))
          ! If all conditions are satisfied, the point is in the box
          IF (ALL(in_lower .AND. in_upper)) THEN
             ! Find the dimension of max magnitude (Linf) difference
             dim = MAXLOC(ABS(pts(:,other_pt_idx) &
                  - pts(:,step)),1)
             ! Check if we are adjusting the upper or lower width of
             ! the box, and adjust the boundary for this point
             IF (pts(dim,step) > pts(dim,other_pt_idx)) THEN
                box_widths(dim,step) = &
                     pts(dim,step) - pts(dim,other_pt_idx)
             ELSE
                box_widths(dim+SIZE(pts,1),step) = &
                     pts(dim,other_pt_idx) - pts(dim,step)
             END IF
          END IF
       END DO calculate_box_width
    END DO compute_box_shapes
  END SUBROUTINE max_boxes


  SUBROUTINE eval_box_coefs(boxes, widths, x_points, box_vals)
    ! Subroutine for computing the values of all the boxes at an x-point
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: boxes
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:)   :: widths
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: x_points
    REAL (KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(boxes,2)&
         &,SIZE(x_points,2)) :: box_vals
    ! ^^ Holder for the evaluation of each box at the given x-point
    INTEGER :: pt_num, box_num
    REAL (KIND=REAL64), DIMENSION(SIZE(boxes,1)) :: shifted_box
    ! ^^ local variables for computing the box evaluations

    ! Compute the linear box spline for each cube
    box_evals : DO CONCURRENT (pt_num=1:SIZE(x_points,2), &
         box_num=1:SIZE(boxes,2))
       ! Calculate the normalized distance from center of this box
       shifted_box = ABS((boxes(:,box_num) - x_points(:,pt_num)) / &
            widths(box_num))
       ! Store the evaluation of this box for this x-data point

       ! Order 1 approximation
       box_vals(box_num,pt_num) = PRODUCT(&
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
       ! box_vals(box_num,pt_num) = PRODUCT(shifted_box)

    END DO box_evals
  END SUBROUTINE eval_box_coefs


  SUBROUTINE eval_boxes(boxes, widths, weights, x_points, response)
    ! Subroutine for evaluating an approximation at a set of x-points
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: boxes
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:)   :: widths
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:)   :: weights
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: x_points
    REAL (KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(x_points,2))  :: response
    ! Local variables
    INTEGER :: x_pt
    REAL (KIND=REAL64), DIMENSION(SIZE(boxes,2),SIZE(x_points,2)) :: box_coefs

    ! Get the box-spline values for each box at each point
    CALL eval_box_coefs(boxes, widths, x_points, box_coefs)

    DO CONCURRENT (x_pt=1:SIZE(x_points,2))
       ! Calculate the value using the box-splines and weights
       response(x_pt) = SUM(box_coefs(:,x_pt) * weights)
    END DO
  END SUBROUTINE eval_boxes


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

END MODULE bootstrapped_box_splines


