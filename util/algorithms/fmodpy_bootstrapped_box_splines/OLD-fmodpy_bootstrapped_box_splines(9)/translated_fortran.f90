
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

MODULE BOOTSTRAPPED_BOX_SPLINES
USE ISO_FORTRAN_ENV , ONLY : INT64 , REAL64

SUBROUTINE COMPUTE_BOXES ( X_DATA , RESPONSE , BOXES , WIDTHS , WEIGHTS , BATCH_SIZE )
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

REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : X_DATA
! The matrix of data points (stored in rows as x-vectors)

REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( SIZE ( X_DATA , 2 ) ) : : RESPONSE
! The vector of response values associated with each point

REAL ( KIND = REAL64 ) , INTENT ( INOUT ) , DIMENSION ( : , : ) : : BOXES
! Holder for the center coordinates of each box-spline support region

REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( SIZE ( BOXES , 2 ) ) : : WIDTHS
! Holder for the width of each box-spline support region

REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( SIZE ( BOXES , 2 ) ) : : WEIGHTS
! Holder for the multiplier weight of each box-spline basis
!   (solved for with LAPACK least squares)

INTEGER ( KIND = INT64 ) , INTENT ( IN ) , OPTIONAL : : BATCH_SIZE
! Number of boxes to add per iteration

!     Local workspace variables
!===================================

INTEGER ( KIND = INT64 ) : : LOCAL_BATCH_SIZE
! Number of boxes to add per iteration

REAL ( KIND = REAL64 ) , DIMENSION ( : ) , ALLOCATABLE : : DGELS_WORK_ARRAY
! Self-explanatory name

REAL ( KIND = REAL64 ) , DIMENSION ( SIZE ( BOXES , 2 ) , SIZE ( X_DATA , 2 ) ) : : COEF_MATRIX
! The matrix of coefficients in the least squares solve, holds
! the evaluation of all box-splines at each x data point

REAL ( KIND = REAL64 ) , DIMENSION ( SIZE ( X_DATA , 2 ) ) : : LOCAL_RESPONSE
! Local copy of 'response' so that the original is not tampered with

REAL ( KIND = REAL64 ) , DIMENSION ( ( SIZE ( X_DATA , 2 ) * ( SIZE ( X_DATA , 2 ) - 1 ) ) / 2 ) : : DISTANCES
! Holder for the pairwise distances between data points (m*(m-1))/2

INTEGER , DIMENSION ( ( SIZE ( X_DATA , 2 ) * ( SIZE ( X_DATA , 2 ) - 1 ) ) / 2 ) : : DISTANCE_INDICES
! Holder for the original un-sorted indices of pairwise distances

REAL ( KIND = REAL64 ) , DIMENSION ( 2 * SIZE ( X_DATA , 1 ) , SIZE ( BOXES , 2 ) ) : : BOX_WIDTHS
! Storage for hyperbox dimensions before they are converted to
! hypercube, box_widths(:SIZE(x_data,1),:) holds the lower
! widths, box_widths(SIZE(x_data,1)+1:,:) holds the upper widths

LOGICAL , DIMENSION ( SIZE ( X_DATA , 2 ) ) : : IS_CENTER_POINT
! Holder for which data points are box centers

INTEGER : : NUM_BOX_CENTERS , STEP , IDX , FIRST , SECOND , DIM , INFO
! num_box_centers -- Number of box centers, SUM(IF is_center_point 1)
! step -- Iterator used for stepping through various algorithms
! idx -- An index in the list of pairwise distances or boxes
! first -- The index of the first point in a pair
! second -- The index of the second point in a pair
! dim -- Index of which dimension is being checked for a data point
! info -- Used to send / receive messages from DGELS

REAL ( KIND = REAL64 ) : : RANDOM_NUM
! random_num -- Holder for random numbers in the range [0,1]

INTEGER : : BATCH , NEW_BOXES
! batch -- Current batch number
! new_boxes -- Counter for the number of new boxes that have been added
INTEGER , DIMENSION ( SIZE ( BOXES , 2 ) ) : : BATCH_ERRORS
! Array of error achieved for each batch
INTEGER , DIMENSION ( SIZE ( RESPONSE , 1 ) ) : : ERROR_INDICES
! Indices for sorting errors
REAL ( KIND = REAL64 ) , DIMENSION ( SIZE ( RESPONSE , 1 ) ) : : ERRORS
! Container for actual errors
INTEGER , DIMENSION ( SIZE ( BOXES , 2 ) ) : : BOX_INDICES

!     Initialization
!========================
! Set the batch size according to whether or not it was given

! Initialization of local arrays
! Initialize the list of indices of box-centers

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

! Sorting the pairwise distances

! Picking the box centers to keep (well spaced, QNSTOP algorithm)
! Starting with the closest pair and working to more spaced pairs,
! randomly select points from each pair to be dropped until there
! are only "n" points left to use as the centers of boxes.
! Identify the two points that compose the next closest pair
! Map (pair index) -> (first point index, second point index)
! If both are still in the 'to keep' set, then remove one randomly

! WARNING: Re-using memory, this is actually sorting response values
! Set all center points to false
! Make the median response valued x-point the initial box center.

! Initialize all box centers basd on the "is_center_point" array

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

! Calculate the shapes of the max-boxes
! Use the smallest bounding hyper-cube for each box size

! Call LAPACK DGELS to query the appropriate size for the work array
! Check for errors in identifying amount of storage space necessary
! Allocate space for the working array

! Primary loop for performing the bootstrapping
! TODO: Update runtime analysis of bootstrapping process
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

! Copy in the true respose values for the DGELS solve
! Call the LAPACK routine for doing least squares minimization of Ax = B
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

! Extract the weights from the output array for DGELS

!     Step 4: Evaluate current model error
!==============================================

! Now we have the indices (sorted from least to largest
! magnitude to smallest magnitude)


! TODO: Stopping condition based on change in error
! IF (should_stop(batch_errors(:batch))) THEN
!    EXIT bootstrapping_main_loop
! END IF

!     Step 5: Add new boxes where error is large
!====================================================

! If we cannot add another box, the break
! Compensate for the false "new_box" that was never added
! Find suitable new error point to add
! Compensate for the false "new_box" that was never added
! If we have tried all possible points, then done
! While the current idx error point is in an already
! added box, keep searching
! If we have found a new error point not inside any
! new boxes,

! Add new box to the list of box centers

! Re-calculate the box width of *just the new box*
! Use the smallest bounding hyper-cube for this new box size
! Revert the number of box centers to be that minus new boxes

! Recalculate the shapes of the current max-boxes
! Use the smallest bounding hyper-cube for each box size


!     Final Step: Recalculate Full Box Mesh and Fit
!=======================================================

! Copy in the true respose values for the DGELS solve
! Call the LAPACK routine for doing least squares minimization of Ax = B
! Extract the weights from the output array for DGELS



! Function for evaluating whether the bootstrapping algorithm
! should stop based on the progression of errors seen
FUNCTION SHOULD_STOP ( ERRORS )
LOGICAL : : SHOULD_STOP
REAL ( KIND = REAL64 ) , DIMENSION ( : ) : : ERRORS
END FUNCTION SHOULD_STOP

FUNCTION BOXES_CONTAIN_POINT ( BOXES , WIDTHS , PT )
LOGICAL : : BOXES_CONTAIN_POINT
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : BOXES
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( SIZE ( BOXES , 2 ) ) : : WIDTHS
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( SIZE ( BOXES , 1 ) ) : : PT
INTEGER : : IDX
! Initialize to assuming that no boxes contain the point
! Search through all of the given boxes
! If any box contains the point, set true and break
END FUNCTION BOXES_CONTAIN_POINT

! Function for calculating the index number of a point based on
! the index of the pair (assuming pairs are [1,2], [1,3]... [2,3]...)
SUBROUTINE IDX_TO_PAIR ( IDX , NUM_INDIVIDUALS , FIRST , SECOND )
INTEGER , INTENT ( IN ) : : IDX , NUM_INDIVIDUALS
INTEGER , INTENT ( OUT ) : : FIRST , SECOND
REAL ( KIND = REAL64 ) : : PAIR_NUM , I
! pair_num -- The number whose integer part represents the index
!             of the first point and fractional part the index of
!             second point in "distance_indices"
! i        -- holds the converted index (for the math to work)

! First figure out how many pairs there will be, convert to
! backwards index (required to use functional form)
! Get the pair number (using function whos whole + remainder form
! the sequence [1,2], [1,3], ... [1,num_individuals], [2,3], ...)
! 'first' is determined by the whole part of "pair_num"
! 'second' is determined by the (rounded) remainder part of "pair_num"
! Shift first and second to start at the front of the list
END SUBROUTINE IDX_TO_PAIR

END SUBROUTINE COMPUTE_BOXES


! Subroutine for computing the max-box around given points
!  where the max-box is the box defined by a center, lower, and
!  upper widths such that MIN(lower,upper) is maximized
SUBROUTINE MAX_BOXES ( PTS , BOX_CENTERS , BOX_WIDTHS )
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : PTS
! Column vectors of points
INTEGER , INTENT ( IN ) , DIMENSION ( : ) : : BOX_CENTERS
! Box center points (indices of points in "pts")
REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( 2 * SIZE ( BOX_CENTERS , 1 ) , SIZE ( BOX_CENTERS , 1 ) ) : : BOX_WIDTHS
! Holder for the upper and lower widths surrounding each box
REAL ( KIND = REAL64 ) , DIMENSION ( SIZE ( PTS , 2 ) ) : : DISTANCES
! Distances between points
INTEGER , DIMENSION ( SIZE ( PTS , 2 ) ) : : DISTANCE_INDICES
! Indices of distances (for use in sorting)
LOGICAL , DIMENSION ( SIZE ( PTS , 1 ) ) : : IN_LOWER , IN_UPPER
! Holder for checking if a point is inside of a box
INTEGER : : STEP , BOX , DIM , IDX , OTHER_PT_IDX , OTHER_PT
! Miscellaneous integers for stepping

! Initialize all box widths to be 'undefined'
! box_widths = HUGE(box_widths(1,1))
! Iterate over box centers (points)
! Find the infinity norm distances to each other point
! Sort the Linf distances between pts
! Identify the max-Linf-width bounding rectangle
! Don't try to use the current box center as a boundary
! Logical testing if a point is inside a box:
!  -- lower width undefined .OR. (box - other) < lower width
!  -- upper width undefined .OR. (other - box) < upper width
! If all conditions are satisfied, the point is in the box
! Find the dimension of max magnitude (Linf) difference
! Check if we are adjusting the upper or lower width of
! the box, and adjust the boundary for this point
END SUBROUTINE MAX_BOXES


SUBROUTINE EVAL_BOX_COEFS ( BOXES , WIDTHS , X_POINTS , BOX_VALS )
! Subroutine for computing the values of all the boxes at an x-point
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : BOXES
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : ) : : WIDTHS
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : X_POINTS
REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( SIZE ( BOXES , 2 ) , SIZE ( X_POINTS , 2 ) ) : : BOX_VALS
! ^^ Holder for the evaluation of each box at the given x-point
INTEGER : : PT_NUM , BOX_NUM
REAL ( KIND = REAL64 ) , DIMENSION ( SIZE ( BOXES , 1 ) ) : : SHIFTED_BOX
! ^^ local variables for computing the box evaluations

! Compute the linear box spline for each cube
! Calculate the normalized distance from center of this box
! Store the evaluation of this box for this x-data point

! Order 1 approximation

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

END SUBROUTINE EVAL_BOX_COEFS


SUBROUTINE EVAL_BOXES ( BOXES , WIDTHS , WEIGHTS , X_POINTS , RESPONSE )
! Subroutine for evaluating an approximation at a set of x-points
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : BOXES
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : ) : : WIDTHS
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : ) : : WEIGHTS
REAL ( KIND = REAL64 ) , INTENT ( IN ) , DIMENSION ( : , : ) : : X_POINTS
REAL ( KIND = REAL64 ) , INTENT ( OUT ) , DIMENSION ( SIZE ( X_POINTS , 2 ) ) : : RESPONSE
! Local variables
INTEGER : : X_PT
REAL ( KIND = REAL64 ) , DIMENSION ( SIZE ( BOXES , 2 ) , SIZE ( X_POINTS , 2 ) ) : : BOX_COEFS

! Get the box-spline values for each box at each point

! Calculate the value using the box-splines and weights
END SUBROUTINE EVAL_BOXES


SUBROUTINE QSORTC ( A , IDX )
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
REAL ( KIND = REAL64 ) , DIMENSION ( : ) , INTENT ( IN OUT ) : : A
INTEGER , DIMENSION ( SIZE ( A ) ) , INTENT ( OUT ) : : IDX

! Local variables

INTEGER : : I

! Initialize the array of original positions.


RECURSIVE SUBROUTINE QSORTC_HELPER ( A , IDX , ISTART , ISTOP )
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
REAL ( KIND = REAL64 ) , DIMENSION ( : ) , INTENT ( IN OUT ) : : A
INTEGER , DIMENSION ( SIZE ( A ) ) , INTENT ( IN OUT ) : : IDX
INTEGER , INTENT ( IN ) : : ISTART , ISTOP

!  Local variables
INTEGER : : ILEFT
INTEGER : : IMID
INTEGER : : IRIGHT
INTEGER : : ITEMP
REAL ( KIND = REAL64 ) : : ATEMP
REAL ( KIND = REAL64 ) : : PIVOT

! INSMAX is used to stop recursively dividing the array and to instead
! use a sort that is more efficient for small arrays than quicksort.
!
! The best cutoff point is system dependent.

INTEGER , PARAMETER : : INSMAX = 24

! Check to see if we have enough values to make quicksort useful.
! Otherwise let the insertion sort handle it.


! Use the median of the first, middle and last items for the pivot
! and place the median (pivot) at the end of the list.
! Putting it at the end of the list allows for a guard value to keep
! the loop from falling off the right end of the array (no need to
! check for at the end of the subarray EACH time through the loop).








! Now, the first position has a value that is less or equal to the
! partition. So, we know it belongs in the left side of the partition
! and we can skip it. Also, the pivot is at the end.  So, there is
! no need to compare the pivot with itself.


! Find a value in the left side that is bigger than the pivot value.
! Pivot is at the right end so ILEFT will not fall off the end
! of the subarray.



! Now we have a value bigger than pivot value on the left side that can be
! swapped with a value smaller than the pivot on the right side.
!
! This gives us all values less than pivot on the left side of the
! array and all values greater than the pivot on the right side.



!
! The last swap was in error (since the while condition is not checked
! until after the swap is done) so we swap again to fix it.

! This is done (once) rather than having an if (done many times) in the
! loop to prevent the swapping.




! Put the pivot value in its correct spot (between the 2 partitions)
! When the WHILE condition finishes, ILEFT is greater than IRIGHT.
! So, ILEFT has the position of the first value in the right side.
! This is where we can put the pivot (and where it will finally rest,
! so no need to look at it again).  Also, place the first value of
! the right side (being displaced by the pivot) at the end of the
! subarray (since it is bigger than the pivot).




END SUBROUTINE QSORTC_HELPER


SUBROUTINE INSERTION ( A , IDX , ISTART , ISTOP )
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

REAL ( KIND = REAL64 ) , DIMENSION ( : ) , INTENT ( IN OUT ) : : A
INTEGER , DIMENSION ( SIZE ( A ) ) , INTENT ( IN OUT ) : : IDX
INTEGER , INTENT ( IN ) : : ISTART , ISTOP

! Local variables.

REAL ( KIND = REAL64 ) : : AMIN
REAL ( KIND = REAL64 ) : : ATEMP
INTEGER : : I
INTEGER : : IABOVE
INTEGER : : IMIN
INTEGER : : ITEMP


! Find the smallest and put it at the top as a "guard" so there is
! no need for the DO WHILE to check if it is going past the top.





! Insertion sort the rest of the array.

! Stop moving items down when the position for "insertion" is found.
!
! Do not have to check for "falling off" the beginning of the
! array since the smallest value is a guard value in the first position.

END SUBROUTINE INSERTION
END SUBROUTINE QSORTC

END MODULE BOOTSTRAPPED_BOX_SPLINES



