
MODULE VTdelaunay
  !
  ! Module containing DelaunayP() subroutine for computing a delaunay
  ! simplex containing a given point.
  !
  ! Written to FORTRAN 2003 standard. 
  ! Requies BLAS and LAPACK
  !
  ! Tested w/ gfortran 4.8.5 on x86_64-redhat-linux       
  ! AUTHOR: Tyler Chang                                   
  ! Last Update: September, 2017                               
  !
  USE real_precision

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: DelaunayP

CONTAINS

  SUBROUTINE DelaunayP(d, n, pts, p, simp, weights, err, &
       eps_opt, rnorm_opt, budget_opt, extrap_opt)
    !
    ! Compute the simplex in the delaunay triangulation containing the point
    ! p. Inspired by DeWall algorithm. Compute a delaunay simplex "close to"
    ! p by solving a series of LS problems. Proceed to connect faces of the
    ! triangulation until a simplex is found containing p.                  
    !
    ! INTEGER (IN) d : Dimension of space.                                  
    ! INTEGER (IN) n : Number of points in set.                             
    ! DOUBLE (IN) pts(d,n) : Set of n points in d-space.                    
    ! DOUBLE (IN) p(d) : A single point in d-space to interpolate using     
    !       the delaunay triangulation.                        
    ! INTEGER (OUT) simp(d+1) : The indices of the points which make up the 
    !       vertices of the simplex containing p.       
    ! DOUBLE (OUT) weights(d+1) : The weights for the vertices of simp for  
    !       interpolating at p.                       
    ! INTEGER (OUT) err : Error flag:
    !       SUCCESSFUL INTERPOLATION = 0
    !       SUCCESSFUL EXTRAPOLATION = 1
    !       See other codes in ErrorCodes.txt file.
    ! OPTIONAL, DOUBLE (IN) eps_opt : A small number defining the tolerance for
    !       accepting "nearly" containing simplices. 0 < eps <<< scale of data.
    !       Default value: eps = 0.1 * DIAMETER(PTS) * SQRT(EPSILON(DOUBLE)).
    ! OPTIONAL, DOUBLE (OUT) rnorm_opt : The 2-norm of the residual vector when
    !       projecting p onto the convex hull of pts. This value is always
    !       calculated when extrapolation is turned on, but only returned when
    !       this value is present.
    ! OPTIONAL, INTEGER (IN) budget_opt : The maximum number of simplices that
    !       the user is willing to construct when searching for the simplex
    !       containing p. In practice, this should occur in no more than 20~30
    !       iterations, with a typical case of 1~10. The default value for
    !       budget is: budget = 1000. If failure to converge in this many
    !       iterations, it is most likely that eps is too small for the scale
    !       of this problem.
    ! OPTIONAL, LOGICAL (IN) extrap_opt : If extrap_opt = .FALSE., then no
    !       extrapolation will be done when p is outside the convex hull and
    !       an error flag will be returned. By default, extrap = .TRUE.
    !
    ! inputs
    INTEGER, INTENT(IN) :: d, n
    REAL(KIND=R8), INTENT(IN) :: pts(d,n), p(d)
    ! outputs
    INTEGER, INTENT(OUT) :: simp(d+1)
    REAL(KIND=R8), INTENT(OUT) :: weights(d+1)
    INTEGER, INTENT(OUT) :: err
    ! optional arguments
    REAL(KIND=R8), INTENT(IN), OPTIONAL:: eps_opt
    REAL(KIND=R8), INTENT(OUT), OPTIONAL :: rnorm_opt
    INTEGER, INTENT(IN), OPTIONAL :: budget_opt
    LOGICAL, INTENT(IN), OPTIONAL :: extrap_opt
    ! local vars
    REAL(KIND=R8) :: currRad, diam, eps, minRad, side1, side2, rnorm
    REAL(KIND=R8), ALLOCATABLE :: A(:,:), b(:), center(:), LQ(:,:), &
         plane(:), proj(:), work(:)
    INTEGER :: budget, i, j, k, tmp
    INTEGER, ALLOCATABLE :: piv(:)
    LOGICAL :: extrap
    ! blas functions
    REAL(KIND=R8), EXTERNAL :: DDOT
    ! external subroutine
    EXTERNAL :: Project

    ! check for common errors
    IF (d < 1 .OR. n < 1) THEN
       err = 10
       RETURN
    ELSE IF (n < d+1) THEN
       err = 11
       RETURN
    END IF
    ! compute diameter
    CALL GetDiameter(d, n, pts, minRad, diam)
    ! check for optional inputs arguments
    IF (PRESENT(eps_opt)) THEN
       eps = eps_opt
       IF (eps <= 0) THEN
          err = 10
          RETURN
       END IF
    ELSE
       eps = diam*SQRT(EPSILON(eps))
    END IF
    IF (PRESENT(budget_opt)) THEN
       budget = budget_opt
       IF (budget < 1) THEN
          err = 10
          RETURN
       END IF
    ELSE
       budget = 1000
    END IF
    IF (PRESENT(extrap_opt)) THEN
       extrap = extrap_opt
    ELSE
       extrap = .TRUE.
    END IF
    ! check for degeneracies in point spacing
    IF (minRad < eps) THEN
       err = 13
       RETURN
    END IF
    ! allocate data
    ALLOCATE(A(d,d), b(d), center(d), LQ(d,d), piv(d), &
         plane(d+1), proj(d), work(5*d), STAT=err)
    IF (err /= 0) THEN
       err = 40
       RETURN
    END IF

    ! initialize the projection
    proj = p
    rnorm = 0.0_R8
    ! construct a delaunay simplex "near" p
    CALL MakeFirstSimp()
    IF (err /= 0) RETURN
    ! begin search for simplex containing p
    DO k = 1, budget
       ! check if contains p
       IF (PtInSimp()) EXIT
       IF (err /= 0) RETURN
       ! swap out least weighted vertex
       i = MINLOC(weights(1:d+1), DIM=1)
       tmp = simp(i)
       simp(i) = simp(d+1)
       ! if i /= 1, simply drop a row from linear system
       IF(i /= 1) THEN
          A(i-1,:) = A(d,:)
          b(i-1) = b(d)
          ! if i == 1, must reconstruct A & b
       ELSE
          DO j=1,d
             A(j,:) = pts(:,simp(j+1)) - pts(:,simp(1))
             b(j) = DDOT(d, A(j,:), 1, A(j,:), 1) / 2.0_R8 
          END DO
       END IF
       ! compute next simplex
       CALL MakeSimplex()
       IF (err /= 0) RETURN
       ! check for extrapolation condition
       IF (simp(d+1) == 0) THEN
          IF (.NOT. extrap) THEN
             err = 3
             RETURN
          END IF
          ! project p onto convex hull
          CALL Project(d, n, pts, p, proj, rnorm, err)
          IF (err /= 0) RETURN
          ! remake previous simplex
          simp(d+1) = tmp
          A(d,:) = pts(:,tmp) - pts(:,simp(1))
          b(d) = DDOT(d, A(d,:), 1, A(d,:), 1) / 2.0_R8
       END IF
    END DO
    ! free heap
    DEALLOCATE(A, b, center, LQ, piv, plane, proj, work, STAT=err)
    IF (err /= 0) THEN
       err = 41
       RETURN
    END IF
    ! check for budget violation
    IF (k > budget) THEN 
       err = 20
       RETURN
    END IF
    ! set extrapolation error flag
    IF (rnorm > eps) THEN
       IF (rnorm < (diam * 0.1_R8)) THEN
          err = 1
       ELSE
          err = 2
       END IF
    END IF
    ! optional outputs
    IF (PRESENT(rnorm_opt)) THEN
       rnorm_opt = rnorm
    END IF
    ! done
    RETURN
    ! end subroutine

    ! internal subroutines & functions
  CONTAINS

    SUBROUTINE MakeFirstSimp()
      !
      ! Build first face by solving sequence of LS problems ranging from 2Xd  
      ! up to [d-1]Xd. At each step the center is given by the minimum norm   
      ! solution to Ac = b. Then the radius is given by r = ||c||_2.          
      ! A^(m) = [ x2 - x1 | x3 - x1 | ... | xm - x1 ] where x1 - xm-1 are the 
      ! current points in the face. The xm which results in the minimum r is  
      ! added to A^(m+1) until an entire face (d points) has been constructed.
      !
      ! Dummy variables U and VT for storing U and V^T when computing SVD.
      ! Values are never initialized or used.
      REAL(KIND=R8), ALLOCATABLE :: U(:), VT(:)
      ! blas functions & lapack subroutines
      REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
      EXTERNAL :: DGELS, DGESVD
      ! find first point
      simp = 0
      minRad = HUGE(minRad)
      DO i = 1, n
         ! check distance to p
         currRad = DNRM2(d, pts(:,i) - proj(:), 1)
         IF (currRad < minRad) THEN
            minRad = currRad
            simp(1) = i
         END IF
      END DO
      ! check that a point was found
      IF (simp(1) == 0) THEN
         err = 81
         RETURN
      END IF
      ! find second point
      minRad = HUGE(minRad)
      DO i = 1, n
         ! avoid repeats
         IF (i == simp(1)) CYCLE
         ! check radius
         currRad = DNRM2(d, pts(:,i)-pts(:,simp(1)), 1)
         IF (currRad < minRad) THEN
            minRad = currRad
            simp(2) = i
         END IF
      END DO
      ! check that a point was found
      IF (simp(2) == 0) THEN
         err = 81
         RETURN
      END IF
      ! set up first row of LS system
      A(1,:) = pts(:,simp(2)) - pts(:,simp(1))
      b(1) = DDOT(d, A(1,:), 1, A(1,:), 1) / 2.0_R8
      ! loop over remaining d-1 points in first simplex
      DO i = 2, d
         ! re-init radius for each iteration
         minRad = HUGE(minRad)
         ! find next point to add
         DO j = 1, n
            ! check point is not already in face
            IF (ANY(simp == j)) CYCLE
            ! add curr point to LS system
            A(i,:) = pts(:,j) - pts(:,simp(1))
            b(i) = DDOT(d, A(i,:), 1, A(i,:), 1) / 2.0_R8
            ! solve LS problem for min norm solution
            LQ(1:i,:) = A(1:i,:)
            center(1:i) = b(1:i)
            ! solve LS problem using QR/LQ factorization
            CALL DGELS('N', i, d, 1, LQ, d, center, d, &
                 work, 5*d, err)
            ! check for errors
            IF (err < 0) THEN
               err = 80
               RETURN
               ! indicates rank-deficiency detected
            ELSE IF (err > 0) THEN
               CYCLE
            END IF
            ! calculate radius
            currRad = DNRM2(d, center, 1)
            ! check for new minimum
            IF (currRad < minRad) THEN
               minRad = currRad
               simp(i+1) = j
            END IF
         END DO
         ! check that a point was found
         IF (simp(i+1) == 0) THEN
            EXIT
         END IF
         ! add next point to LS system
         A(i,:) = pts(:,simp(i+1)) - pts(:,simp(1))
         b(i) = DDOT(d, A(i,:), 1, A(i,:), 1) / 2.0_R8 
      END DO
      ! If no rank deficiency was detected, double-check
      IF (i > d) THEN
         ! compute the SVD
         LQ = A
         CALL DGESVD('N','N',d,d,LQ,d,plane(1:d),U,1,VT,1,work,5*d,err)
         ! check for errors
         IF (err < 0) THEN
            err = 80
            RETURN
         ELSE IF (err > 0) THEN
            err = 99
            RETURN
         END IF
         ! Check for rank deficiency up to working precision.
         ! If found, re-compute first face using SVD.
         IF (plane(d)/plane(1) < eps) THEN
            CALL MakeFirstSimp_RankDef()
            IF (err /= 0) RETURN
         END IF
         ! Otherwise rank deficiency already detected
      ELSE
         CALL MakeFirstSimp_RankDef()
         IF (err /= 0) RETURN

      END IF
      ! done
      err = 0
      RETURN
    END SUBROUTINE MakeFirstSimp

    SUBROUTINE MakeFirstSimp_RankDef()
      !
      ! Build first face by solving sequence of LS problems ranging from 2Xd  
      ! up to [d-1]Xd. At each step the center is given by the minimum norm   
      ! solution to Ac = b. Then the radius is given by r = ||c||_2.          
      ! A^(m) = [ x2 - x1 | x3 - x1 | ... | xm - x1 ] where x1 - xm-1 are the 
      ! current points in the face. The xm which results in the minimum r is  
      ! added to A^(m+1) until an entire face (d points) has been constructed.
      ! Use the SVD of A to check for rank defficiency at each step.
      !
      ! blas functions & lapack subroutines
      REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
      EXTERNAL :: DGELSS
      ! We already found the first 2 points, their indices are in simp(1)
      ! and simp(2) and the LS system components A(1,:) and b(1) are
      ! constructed.
      simp(3:d+1) = 0
      ! Loop over remaining d-1 points in first simplex
      DO i = 2, d
         ! re-init radius for each iteration
         minRad = HUGE(minRad)
         ! find next point to add
         DO j = 1, n
            ! check point is not already in face
            IF (ANY(simp == j)) CYCLE
            ! add curr point to LS system
            A(i,:) = pts(:,j) - pts(:,simp(1))
            b(i) = DDOT(d, A(i,:), 1, A(i,:), 1) / 2.0_R8
            ! solve LS problem for min norm solution
            LQ(1:i,:) = A(1:i,:)
            center(1:i) = b(1:i)
            ! Solve using SVD. Store effective rank in tmp.
            CALL DGELSS(i, d, 1, LQ, d, center, d, plane, eps, &
                 tmp, work, 5*d, err)
            ! check for errors
            IF (err < 0) THEN
               err = 80
               RETURN
            ELSE IF (err > 0) THEN
               err = 99
               RETURN
            END IF
            ! if rank deficient, cycle
            IF (tmp < i) CYCLE
            ! calculate radius
            currRad = DNRM2(d, center, 1)
            ! check for new minimum
            IF (currRad < minRad) THEN
               minRad = currRad
               simp(i+1) = j
            END IF
         END DO
         ! check that a point was found
         IF (simp(i+1) == 0) THEN
            err = 12
            RETURN
         END IF
         ! add next point to LS system
         A(i,:) = pts(:,simp(i+1)) - pts(:,simp(1))
         b(i) = DDOT(d, A(i,:), 1, A(i,:), 1) / 2.0_R8 
      END DO
      ! done
      err = 0
      RETURN
    END SUBROUTINE MakeFirstSimp_RankDef

    SUBROUTINE MakeSimplex()
      !
      ! Given a face of a delaunay simplex, complete the simplex by adding a  
      ! point on the same side of the face as p. The new point must minimize  
      ! the radius of the circumcircle. Determine the radius by solving the   
      ! linear system Ac = b, then taking r = ||c||_2. Note that              
      ! A = [ x2 - x1 | x3 - x1 | ... | x(d+1) - x1 ] where x1 ... xd are the 
      ! points in the face and x(d+1) is a potential new point (completing the
      ! simplex). b_i = dot(A_i,A_i) / -2.                                    
      !
      ! blas functions & lapack subroutines
      REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
      EXTERNAL :: DGESV
      ! calculate halfspace
      CALL MakePlane()
      IF(err /= 0) RETURN
      ! calculate side of the viable points
      side1 =  DDOT(d,plane(1:d),1,proj(:),1) - plane(d+1)
      ! normalize magnitude of p
      side1 = SIGN(1.0_R8,side1)
      ! initialize center, radius, and simplex
      simp(d+1) = 0
      center = 0_R8
      minRad = HUGE(minRad)
      ! loop through all points
      DO i = 1, n
         ! check point(i) is inside current ball
         IF (DNRM2(d, pts(:,i) - center(:), 1) > minRad) CYCLE
         ! check point(i) on viable side of halfspace
         side2 = DDOT(d,plane(1:d),1,pts(:,i),1) - plane(d+1)
         IF (side1 * side2 < eps .OR. ANY(simp == i)) CYCLE
         ! add point i to linear system
         A(d,:) = pts(:,i) - pts(:,simp(1))
         b(d) = DDOT(d, A(d,:), 1, A(d,:), 1) / 2.0_R8
         LQ = A
         center = b
         ! solve linear system to get center
         CALL DGESV(d, 1, LQ, d, piv, center, d, err)
         IF (err < 0) THEN
            err = 80
            RETURN
         ELSE IF (err > 0) THEN
            err = 21
            RETURN
         END IF
         ! update radius, center, and simplex
         currRad = DNRM2(d, center, 1)
         minRad = currRad
         center = center + pts(:,simp(1))
         simp(d+1) = i
      END DO
      ! reset error flag
      err = 0
      ! check for extrapolation condition
      IF(simp(d+1) == 0) RETURN
      ! otherwise add point to linear system
      A(d,:) = pts(:,simp(d+1)) - pts(:,simp(1))
      b(d) = DDOT(d, A(d,:), 1, A(d,:), 1) / 2.0_R8
      ! done
      RETURN
    END SUBROUTINE MakeSimplex

    SUBROUTINE MakePlane()
      !
      ! Construct a plane intersecting the current d points. The coefficients 
      ! are determined by the normal vector: v and an intercept, stored in    
      ! plane(1:d) and plane(d+1) respectively. To find the normal vector, the
      ! current d points are transformed into d-1 vectors in the plane.       
      ! Together, these d-1 d-vectors have nullity rank-1. The normal vector  
      ! is any nontrivial vector in the nullspace.                            
      !
      ! blas functions & lapack subroutines
      REAL(KIND=R8), EXTERNAL :: DDOT, DNRM2
      EXTERNAL :: DGEQP3, DTRSM
      ! check that d-1 > 0
      IF (d > 1) THEN
         ! get QR
         piv=0
         LQ = A
         CALL DGEQP3(d-1,d,LQ,d,piv,plane,work,5*d,err)
         IF(err < 0) THEN
            err = 80
            RETURN
         ELSE IF (err > 0) THEN
            err = 21
            RETURN
         END IF
         ! do triangular solve to get plane
         ! store in work array
         work(1:d-1) = LQ(1:d-1,d)
         CALL DTRSM('L','U','N','N',d-1,1,-1.0_R8,LQ,d,work,d)
         work(d) = 1.0_R8
         ! undo pivots
         DO i = 1,d
            plane(piv(i)) = work(i)
         END DO
         ! normalize
         plane(1:d) = plane(1:d) / DNRM2(d,plane(1:d),1)
         ! calculate intercept
         plane(d+1) = DDOT(d,plane(1:d),1,pts(:,simp(1)),1)
      ELSE
         ! compute plane for d=1
         plane(1) = 1.0_R8
         plane(2) = pts(1,simp(1))
      END IF
      ! done 
      RETURN
    END SUBROUTINE MakePlane

    FUNCTION PtInSimp() RESULT(tf)
      !
      ! Determine if a point is in the current simplex. The simplex can also  
      ! be thought of as a cone of vectors s(i) - s(1) where each s(k) is a   
      ! vertex for the simplex. The d vectors defined thusly span the entire  
      ! d-dimensional space. If the linear combination which creates p - s(1) 
      ! is convex (all values are positive), then the cone contains p. Also,  
      ! the weights for s(2:d+1) for interpolating p will be given by the d   
      ! linear coefficients.                                                  
      !
      ! return value
      LOGICAL :: tf
      ! lapack subroutines
      EXTERNAL :: DGESV
      ! init result
      tf = .FALSE.
      ! solve linear system to get weights
      LQ = TRANSPOSE(A)
      weights(2:d+1) = proj(:) - pts(:,simp(1))
      CALL DGESV(d, 1, LQ, d, piv, weights(2:d+1), d, err)
      IF (err < 0) THEN
         err = 80
         RETURN
      ELSE IF (err > 0) THEN
         err = 21
         RETURN
      END IF
      ! get last interpolation weight
      weights(1) = 1.0_R8 - SUM(weights(2:d+1))
      ! check if weights define a convex combination
      IF (ALL(weights >= -eps)) tf = .TRUE.
      ! done
      RETURN
    END FUNCTION PtInSimp

  END SUBROUTINE DelaunayP

  SUBROUTINE GetDiameter(d, n, pts, minD, maxD)
    !
    ! Compute the diameter of a point set by brute force.
    ! O(n^2) complexity.
    ! INTEGER (IN) d : Dimension of space.                                  
    ! INTEGER (IN) n : Number of points in set.                             
    ! DOUBLE (IN) pts(d,n) : Set of n points in d-space.                    
    ! DOUBLE (OUT) minD : Minimum distance between 2 points
    ! DOUBLE (OUT) maxD : Maximum distance between 2 points
    ! INTEGER(OUT) error : Error flag
    !
    ! inputs
    INTEGER, INTENT(IN) :: d, n
    REAL(KIND=R8), INTENT(IN) :: pts(d,n)
    ! outputs
    REAL(KIND=R8), INTENT(OUT) :: minD, maxD
    ! local vars
    REAL(KIND=R8) :: distance
    INTEGER :: i, j
    ! blas function
    REAL(KIND=R8), EXTERNAL :: DNRM2
    ! initialize values
    minD = HUGE(minD)
    maxD = 0
    ! compute min and max distance in dataset
    DO i = 1, n-1
       DO j = i+1, n
          distance = DNRM2(d, pts(1:d,i) - pts(1:d,j), 1)
          IF (distance > maxD) THEN
             maxD = distance
          END IF
          IF (distance < minD) THEN
             minD = distance
          END IF
       END DO
    END DO
    ! done
    RETURN
  END SUBROUTINE GetDiameter

END MODULE VTdelaunay


