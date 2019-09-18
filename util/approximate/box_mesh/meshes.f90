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

  SUBROUTINE build_ibm(pts, box_widths)
    ! Given a set of "points", record the upper and lower widths of
    ! boxes in "structure". Assume that all boxes except the last
    ! were correctly sized before the addition of the last
    ! point. Reshape all boxes that need to be reshaped because of
    ! the single new addition to the end of "points".

    REAL(KIND=REAL64), DIMENSION(:,:), INTENT(IN) :: pts
    REAL(KIND=REAL64), DIMENSION(2*SIZE(pts,1), SIZE(pts,2)),&
         & INTENT(INOUT) :: box_widths
    ! Holder for the upper and lower widths surrounding each box
    REAL(KIND=REAL64), DIMENSION(SIZE(pts,2)) :: distances
    ! Distances between points
    INTEGER, DIMENSION(SIZE(pts,2)) :: distance_indices
    ! Indices of distances (for use in sorting)
    LOGICAL, DIMENSION(SIZE(pts,1)) :: in_lower, in_upper
    ! Holder for checking if a point is inside of a box
    INTEGER :: step, box, dim, idx, other_pt_idx, other_pt
    ! Miscellaneous integers for stepping

    ! Cycle through all pairs of points and make sure the boxes do not
    ! contain neighbors. This search could be made more efficient,
    ! but that constitutes future work on this code.
    ! 
    ! For example, one could build the mesh one-step-at-a-time and
    ! assume all boxes except the last were already built correctly,
    ! add a new box and reshape old boxes that may have been impacted.
    compute_box_shapes : DO box = 1,SIZE(pts,2)
       compute_other_points : DO other_pt_idx = 1,SIZE(pts,2)
          IF (box .EQ. other_pt_idx) CYCLE compute_other_points
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
             ELSE
                box_widths(dim+SIZE(pts,1),box) = &
                     pts(dim,other_pt_idx) - pts(dim,box)
             END IF
          END IF rebuild_box
       END DO compute_other_points
    END DO compute_box_shapes
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
    ! Make the weights convex.
    convexify: DO pt_num=1, SIZE(box_vals,1)
       box_vals(pt_num,:) = box_vals(pt_num,:) / SUM(box_vals(pt_num,:))
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
       box_vals(pt_num,:) = box_vals(pt_num,:) / SUM(box_vals(pt_num,:))
    END DO convexify

  END SUBROUTINE eval_order_2_box_mesh


END MODULE meshes
