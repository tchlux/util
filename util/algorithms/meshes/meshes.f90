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
  !                         Voronoi Mesh                         
  !==============================================================

  SUBROUTINE train_vm(points, values, control_points,&
       & control_values, mesh_coefficients, control_dots, &
       & num_control_points, error_tolerance)  
    ! Given a set of approximation points along with their response
    ! values, fit a voronoi mesh with maximum error 'error_tolerance'.

    ! Input / Output variables
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: points
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2)), INTENT(IN) :: values
    REAL(KIND=REAL64), DIMENSION(SIZE(points,1),SIZE(points,2)),&
         & INTENT(OUT) :: control_points 
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2)), INTENT(OUT) :: &
         & control_values, mesh_coefficients
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2),SIZE(points,2)), &
         INTENT(OUT) :: control_dots
    INTEGER(KIND=INT64), INTENT(OUT) :: num_control_points
    REAL(KIND=REAL64), INTENT(IN), OPTIONAL :: error_tolerance
    ! Local variables
    REAL(KIND=REAL64) :: local_error_tolerance, error
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2), SIZE(points,2)) :: mesh_values
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2), SIZE(points,2)) :: point_dots
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2)) :: predictions
    INTEGER, DIMENSION(SIZE(points,2)) :: control_indices
    LOGICAL, DIMENSION(SIZE(points,2)) :: used_points
    INTEGER :: step, step2

   
    ! Declaration of optional "error_tolerance"
    IF (PRESENT(error_tolerance)) THEN
       local_error_tolerance = error_tolerance * (MAXVAL(values) - MINVAL(values))
    ELSE
       ! Default maximum error tolerance is 5% error
       local_error_tolerance = 0.05 * (MAXVAL(values) - MINVAL(values))
    END IF
    ! Initialization of 'used points' to be all false
    used_points = .FALSE.
    mesh_values = 0.0
    point_dots = 0
    control_dots = 0
    control_points = 0
    ! Initialize the first control point to be the most central point
    step = most_central_point(points)
    num_control_points = 1
    control_points(:,num_control_points) = points(:,step)
    control_values(num_control_points) = values(step)
    control_indices(num_control_points) = step
    used_points(step) = .TRUE.

    ! Begin the bootstrapping process
    bootstrap_vm : DO
       ! Calculate the next coloumn of dot products between points and controls
       find_all_dots : DO step = 1, SIZE(points,2)
          point_dots(step, num_control_points) = &
               SUM(points(:,step)*control_points(:,num_control_points))
       END DO find_all_dots
       ! Calculate the next column and row of the control dot products
       DO step = 1, num_control_points
          control_dots(num_control_points,step) = &
               point_dots(control_indices(step),num_control_points)
          control_dots(step,num_control_points) = &
               control_dots(num_control_points,step)
       END DO
       ! Evaluate the mesh at all points
       CALL eval_vm(control_points(:,:num_control_points), points,&
            control_dots(:num_control_points,:num_control_points), &
            point_dots(:,:num_control_points), mesh_values)
       ! Fit the mesh to the points
       CALL fit_mesh(mesh_values(:,:num_control_points), values, &
            mesh_coefficients(:num_control_points))
       ! Compute an approximation with the fit coefficients
       compute_approximation : DO step = 1, SIZE(points,2)
          ! Predictions made with fitting
          predictions(step) = SUM(mesh_values(step,:num_control_points) &
               * mesh_coefficients(:num_control_points))
          ! ! Predictions made with NURBS-style blending
          ! predictions(step) = SUM(mesh_values(step,:num_control_points) &
          !      * control_values(:num_control_points)) / &
          !      SUM(mesh_values(step,:num_control_points))
       END DO compute_approximation
       error = NORM2(values - predictions)
       ! error = length(values - predictions)
       stopping_condition : IF ((error < local_error_tolerance) .OR. &
            (num_control_points .GE. SIZE(points,2))) THEN
          EXIT bootstrap_vm
       END IF stopping_condition
       ! Re-use the variable to store the errors associated with predictions
       predictions = ABS(values - predictions)
       ! Find the point with the highest error that is not already a
       ! control point and add that as a control point
       find_next_addition : DO 
          step = MAXLOC(predictions, 1)
          IF (used_points(step)) THEN
             predictions(step) = -1
          ELSE
             EXIT find_next_addition
          END IF
       END DO find_next_addition
       used_points(step) = .TRUE.
       num_control_points = num_control_points + 1
       control_points(:,num_control_points) = points(:,step)
       control_values(num_control_points) = values(step)
       control_indices(num_control_points) = step
    END DO bootstrap_vm
  END SUBROUTINE train_vm

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

    eval_points : DO point_ind = 1, SIZE(points, 2)
       ignore = .FALSE.
       eval_cells : DO control_ind = 1, SIZE(control_points, 2)
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

  !TODO:  Update python call to predict_vm
  SUBROUTINE predict_vm(control_points, values, coefficients,&
       & control_dots, approximation_points, approximations)
    ! Given control points, response values for each control point,
    ! coefficients for each control point, and new approximation
    ! points, calculate approximate response values for each point
    ! using the voronoi mesh.

    ! Input / Output variables
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: control_points
    REAL(KIND=REAL64), DIMENSION(SIZE(control_points,2)), &
         INTENT(IN) :: values, coefficients
    REAL(KIND=REAL64), DIMENSION(SIZE(control_points,2), &
         & SIZE(control_points,2)), INTENT(IN) :: control_dots
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: approximation_points
    REAL(KIND=REAL64), DIMENSION(SIZE(approximation_points,2)), &
         INTENT(OUT) :: approximations
    ! Local variables
    REAL(KIND=REAL64), DIMENSION(SIZE(approximation_points,2), &
         SIZE(control_points,2)) :: mesh_values
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
         control_dots, point_dots,  mesh_values)
    ! Make predictions using the fit coefficients
    DO step = 1, SIZE(approximation_points,2)
       ! Predictions made with fitting
       approximations(step) = SUM(mesh_values(step,:) * coefficients)
       ! ! Predictions made with NURBS-style blending
       ! approximations(step) = (SUM(mesh_values(step,:) * values) / &
       !      SUM(mesh_values(step,:)))
    END DO

  END SUBROUTINE predict_vm

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

  !==============================================================
  !                         Max Box Mesh                         
  !==============================================================

  SUBROUTINE train_mbm(points, values, control_points,&
       & control_values, mesh_structure, mesh_coefficients, &
       & num_control_points, error_tolerance) 
    ! Given a set of approximation points along with their response
    ! values, fit a max box mesh with maximum error 'error_tolerance'
    ! using a batch size of 1.

    ! Input / Output variables
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: points
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2)), INTENT(IN) :: values
    REAL(KIND=REAL64), DIMENSION(SIZE(points,1),SIZE(points,2)),&
         & INTENT(OUT) :: control_points 
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2)), INTENT(OUT) :: &
         & control_values, mesh_coefficients
    REAL(KIND=REAL64), DIMENSION(2*SIZE(points,1), SIZE(points,2)), &
         INTENT(OUT) :: mesh_structure
    INTEGER(KIND=INT64), INTENT(OUT) :: num_control_points
    REAL(KIND=REAL64), INTENT(IN), OPTIONAL :: error_tolerance
    ! Local variables
    REAL(KIND=REAL64) :: local_error_tolerance, error
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2), SIZE(points,2)) :: mesh_values
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2)) :: predictions
    LOGICAL, DIMENSION(SIZE(points,2)) :: used_points
    INTEGER :: step
    ! Declaration of optional "error_tolerance"
    IF (PRESENT(error_tolerance)) THEN
       local_error_tolerance = error_tolerance * (MAXVAL(values) - MINVAL(values))
    ELSE
       ! Default maximum error tolerance is 5% error
       local_error_tolerance = 0.05 * (MAXVAL(values) - MINVAL(values))
    END IF
    ! Initialization of 'used points' to be all false
    mesh_structure = -1.0_REAL64
    used_points = .FALSE.
    mesh_values = 0.0
    ! Initialize the first control point to be the most central point
    step = most_central_point(points)
    control_points(:,1) = points(:,step)
    control_values(1) = values(step)
    used_points(step) = .TRUE.
    num_control_points = 1
    ! Begin the bootstrapping process
    bootstrap_mbm : DO
       ! BUILD THE MESH (assuming only 1 point has been added)
       CALL build_mbm(control_points(:,:num_control_points),&
            & mesh_structure(:,:num_control_points))
       ! EVALUATE MESH
       CALL eval_box_mesh(control_points(:,:num_control_points),&
            & mesh_structure(:,:num_control_points), points, mesh_values)
       ! FIT THE MESH
       CALL fit_mesh(mesh_values(:,:num_control_points), values, mesh_coefficients)
       ! EVALUATE APPROXIMATION
       compute_approximation : DO step = 1, SIZE(points,2)
          ! Predictions made with fitting
          predictions(step) = SUM(mesh_values(step,:num_control_points) &
               * mesh_coefficients(:num_control_points))
          ! ! Predictions made with NURBS-style blending
          ! predictions(step) = SUM(mesh_values(step,:num_control_points) &
          !      * control_values(:num_control_points)) / &
          !      SUM(mesh_values(step,:num_control_points))
       END DO compute_approximation
       error = length(values - predictions)
       ! CHECK STOPPING CONDITION
       stopping_condition : IF ((error < local_error_tolerance) .OR. &
            (num_control_points .GE. SIZE(points,2))) THEN
          EXIT bootstrap_mbm
       END IF stopping_condition
       ! ADD NEW CONTROL POINT TO MESH
       !   Find the point with the highest error that is not already
       !   a control point and add that as a control point
       predictions = ABS(values - predictions)
       find_next_addition : DO 
          step = MAXLOC(predictions, 1)
          IF (used_points(step)) THEN
             predictions(step) = -1
          ELSE
             EXIT find_next_addition
          END IF
       END DO find_next_addition
       used_points(step) = .TRUE.
       num_control_points = num_control_points + 1
       control_points(:,num_control_points) = points(:,step)
       control_values(num_control_points) = values(step)
    END DO bootstrap_mbm
  END SUBROUTINE train_mbm

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

  !==============================================================
  !                      Iterative Box Mesh                      
  !==============================================================

  SUBROUTINE train_ibm(points, values, control_points,&
       & control_values, mesh_structure, mesh_coefficients, &
       & num_control_points, error_tolerance) 
    ! Given a set of approximation points along with their response
    ! values, fit a max box mesh with maximum error 'error_tolerance'
    ! using a batch size of 1.

    ! Input / Output variables
    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: points
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2)), INTENT(IN) :: values
    REAL(KIND=REAL64), DIMENSION(SIZE(points,1),SIZE(points,2)),&
         & INTENT(OUT) :: control_points 
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2)), INTENT(OUT) :: &
         & control_values, mesh_coefficients
    REAL(KIND=REAL64), DIMENSION(2*SIZE(points,1), SIZE(points,2)), &
         INTENT(OUT) :: mesh_structure
    INTEGER(KIND=INT64), INTENT(OUT) :: num_control_points
    REAL(KIND=REAL64), INTENT(IN), OPTIONAL :: error_tolerance
    ! Local variables
    REAL(KIND=REAL64) :: local_error_tolerance, error
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2), SIZE(points,2)) :: mesh_values
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2)) :: predictions
    LOGICAL, DIMENSION(SIZE(points,2)) :: used_points
    INTEGER :: step
    ! Declaration of optional "error_tolerance"
    IF (PRESENT(error_tolerance)) THEN
       local_error_tolerance = error_tolerance * (MAXVAL(values) - MINVAL(values))
    ELSE
       ! Default maximum error tolerance is 5% error
       local_error_tolerance = 0.05 * (MAXVAL(values) - MINVAL(values))
    END IF
    ! Initialization of 'used points' to be all false
    mesh_structure = -1.0_REAL64
    used_points = .FALSE.
    mesh_values = 0.0
    ! Initialize the first control point to be the most central point
    step = most_central_point(points)
    control_points(:,1) = points(:,step)
    control_values(1) = values(step)
    used_points(step) = .TRUE.
    num_control_points = 1
    ! Begin the bootstrapping process
    bootstrap_ibm : DO
       ! BUILD THE MESH (assuming only 1 point has been added)
       CALL build_ibm(control_points(:,:num_control_points),&
            & mesh_structure(:,:num_control_points))
       ! EVALUATE MESH
       CALL eval_box_mesh(control_points(:,:num_control_points),&
            & mesh_structure(:,:num_control_points), points, mesh_values)
       ! FIT THE MESH
       CALL fit_mesh(mesh_values(:,:num_control_points), values, &
            mesh_coefficients(:num_control_points))
       ! EVALUATE APPROXIMATION
       compute_approximation : DO step = 1, SIZE(points,2)
          ! Predictions made with fitting
          predictions(step) = SUM(mesh_values(step,:num_control_points) &
               * mesh_coefficients(:num_control_points))
          ! ! Predictions made with NURBS-style blending
          ! predictions(step) = SUM(mesh_values(step,:num_control_points) &
          !      * control_values(:num_control_points)) / &
          !      SUM(mesh_values(step,:num_control_points))
       END DO compute_approximation
       error = length(values - predictions)
       ! CHECK STOPPING CONDITION
       stopping_condition : IF ((error < local_error_tolerance) .OR. &
            (num_control_points .GE. SIZE(points,2))) THEN
          EXIT bootstrap_ibm
       END IF stopping_condition
       ! ADD NEW CONTROL POINT TO MESH
       !   Find the point with the highest error that is not already
       !   a control point and add that as a control point
       predictions = ABS(values - predictions)
       find_next_addition : DO 
          step = MAXLOC(predictions, 1)
          IF (used_points(step)) THEN
             predictions(step) = -1
          ELSE
             EXIT find_next_addition
          END IF
       END DO find_next_addition
       used_points(step) = .TRUE.
       num_control_points = num_control_points + 1
       control_points(:,num_control_points) = points(:,step)
       control_values(num_control_points) = values(step)
    END DO bootstrap_ibm
  END SUBROUTINE train_ibm

  SUBROUTINE build_ibm(pts, box_widths)
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
    compute_box_shapes : DO box = 1,SIZE(pts,2)-1
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
          ! Find the dimension of max magnitude (Linf) difference
          dim = MAXLOC(ABS(pts(:,other_pt_idx) - pts(:,box)),1)
          ! Check if we are adjusting the upper or lower width of
          ! the box, and adjust the boundary for this point
          IF (pts(dim,box) > pts(dim,other_pt_idx)) THEN
             box_widths(dim,box) = &
                  pts(dim,box) - pts(dim,other_pt_idx)
             box_widths(dim+SIZE(pts,1),other_pt_idx) = &
                  pts(dim,box) - pts(dim,other_pt_idx)
          ELSE
             box_widths(dim+SIZE(pts,1),box) = &
                  pts(dim,other_pt_idx) - pts(dim,box)
             box_widths(dim,other_pt_idx) = &
                  pts(dim,box) - pts(dim,other_pt_idx)
          END IF
       END IF rebuild_box
    END DO compute_box_shapes
  END SUBROUTINE build_ibm

  ! SUBROUTINE train_ibm()
  !   ! Given a set of approximation points along with their response
  !   ! values, fit an iterative box mesh with maximum error
  !   ! 'error_tolerance'.

  !   ! Add a central box
  !   build_mesh : DO 
  !      ! Evaluate basis functions
  !      ! Fit coefficients
  !      ! IF error < tolerance, EXIT
  !      ! Identify point with largest error, add that box to mesh,
  !      !   reshape all boxes that contain that point.
  !   END DO build_mesh
  ! END SUBROUTINE train_ibm

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


  SUBROUTINE predict_box_mesh(boxes, widths, x_points, box_coefs, values)
    ! Given box centers, box widths, response values for each box,
    ! coefficients for each box, and new approximation points,
    ! calculate approximate response values for each point using the
    ! fit box mesh.

    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: boxes
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(SIZE(boxes,1)*2, &
         SIZE(boxes,2)) :: widths
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: x_points
    REAL (KIND=REAL64), INTENT(IN), DIMENSION(SIZE(boxes,2)) :: box_coefs
    REAL (KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(x_points,2)) :: values

    REAL (KIND=REAL64), DIMENSION(SIZE(x_points,2)&
         &,SIZE(boxes,2)) :: box_vals
    INTEGER :: step

    CALL eval_box_mesh(boxes, widths, x_points, box_vals)
    ! EVALUATE APPROXIMATION
    compute_approximation : DO step = 1, SIZE(x_points,2)
       ! Predictions made with fitting
       values(step) = SUM(box_vals(step,:) * box_coefs(:))
       ! ! Predictions made with NURBS-style blending
       ! predictions(step) = SUM(box_vals(step,:) &
       !      * control_values(:num_control_points)) / &
       !      SUM(mesh_values(step,:num_control_points))
    END DO compute_approximation
  END SUBROUTINE predict_box_mesh

  !                          Utilities                          
  !=============================================================

  SUBROUTINE fit_mesh(coef_matrix, response, weights)
    ! Given control points, approximation points and their response
    ! values, identify the coefficients for each control point via a
    ! least-squares solve.
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: coef_matrix
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: response
    REAL(KIND=REAL64), INTENT(OUT), DIMENSION(SIZE(coef_matrix,2)) :: weights
    ! Local variables
    REAL (KIND=REAL64), DIMENSION(:), ALLOCATABLE :: dgels_work_array
    REAL(KIND=REAL64), DIMENSION(SIZE(response,1)) :: local_response
    REAL(KIND=REAL64), DIMENSION(SIZE(coef_matrix,1), &
         & SIZE(coef_matrix,2)) :: local_coef_matrix
    INTEGER :: dim, info
    ! Initialize local variables
    local_coef_matrix = coef_matrix
    local_response = 0
    info = 0
    dim = 0
    ! 
    !     Query Size     
    ! 
    ! Call LAPACK DGELS to query the appropriate size for the work array
    CALL DGELS('N', SIZE(coef_matrix,1), SIZE(coef_matrix,2), 1, &
         local_coef_matrix(:,:), SIZE(coef_matrix,1), local_response, &
         SIZE(local_response,1), local_response, -1, info)
    ! Check for errors in identifying amount of storage space necessary
    IF (info .EQ. 0) THEN
       dim = local_response(1)
    ELSE
       PRINT *, 'ERROR: unsuccessful query of work-array size for DGELS.'
       dim = 2 * SIZE(coef_matrix,2) * SIZE(coef_matrix,1)
    END IF
    ! Allocate space for the working array
    ALLOCATE( dgels_work_array(1:dim) )
    ! 
    !     Solve System     
    ! 
    ! Copy in the true respose values for the DGELS solve
    info = 0
    local_response = response
    ! Call the LAPACK routine for doing least squares minimization of Ax = B
    CALL DGELS('N', SIZE(coef_matrix,1), SIZE(coef_matrix,2), 1, &
         local_coef_matrix, SIZE(coef_matrix,1), local_response, &
         SIZE(local_response,1), dgels_work_array, dim, info)
    !   'N'                 -- No transpose of A is necessary
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
    weights = local_response(1:SIZE(coef_matrix,2))
  END SUBROUTINE fit_mesh

  FUNCTION most_central_point(points)
    ! Given a matrix of points, identify the index of the point that
    ! is most central to those points.
    REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: points
    INTEGER(KIND=INT64) :: most_central_point
    ! Local variables
    INTEGER :: index
    REAL(KIND=REAL64), DIMENSION(SIZE(points,1)) :: center
    REAL(KIND=REAL64), DIMENSION(SIZE(points,2)) :: distances
    center = (MAXVAL(points,2) + MINVAL(points,2)) / 2
    DO index = 1,SIZE(points,2)
       distances(index) = NORM2(points(:,index) - center)
       ! distances(index) = length(points(:,index) - center)
    END DO
    most_central_point = MINLOC(distances,1)
  END FUNCTION most_central_point

  FUNCTION length(vector, norm)
    ! Compute the length of a vector with a given norm, defaults to 2.
    REAL(KIND=REAL64), DIMENSION(:), INTENT(IN) :: vector
    INTEGER(KIND=INT64), INTENT(IN), OPTIONAL :: norm
    REAL(KIND=REAL64) :: length
    ! Local variable for optional norm
    INTEGER(KIND=INT64) :: local_norm
    IF (PRESENT(norm)) THEN
       local_norm = norm
    ELSE 
       local_norm = 2
    END IF
    ! Compute the length
    length = SUM(vector**local_norm) ** (1.0_REAL64 /&
         & REAL(local_norm,REAL64))
  END FUNCTION length

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

END MODULE meshes
