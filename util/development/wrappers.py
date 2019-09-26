# Shepard Method
SHEPARD_WEIGHT_FUNCTION = lambda dist: 1 / max(dist,SMALL_NUMBER)


BBS_DEFAULT_BATCH_SIZE = None
BBS_DEFAULT_NUM_BOXES = float("inf")
FBM_DEFAULT_NUM_BOXES = float("inf")

MBM_CLOSE_TO_SHOW = 10
MBM_DEFAULT_PADDING = 0.01


# ==============================================
#      Arbitrary Nearest Plane Approximator     
# ==============================================

# Class for computing an interpolation between the nearest n neighbors
class NearestPlane(Approximator):
    def __init__(self):
        self.points = None
        self.values = None

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values):
        # Process and store local information
        self.points = control_points.copy()
        self.values = values.copy()

    # Use fortran code to evaluate the computed boxes
    def predict(self, x):
        # Use the nearest point 
        response = []
        for pt in x:
            distances = np.sum((self.points - pt)**2, axis=1)
            closest = np.argsort(distances)
            plane_points = [closest[0], closest[1]]
            to_try = 2
            while ( len(plane_points) <= self.points.shape[1] ):
                if to_try >= len(closest):
                    raise(Exception(("ERROR: No set of %i linearly "+
                                     "independant points in the "+
                                     "training data.")%(self.points.shape[1]+1)))
                new_points = plane_points + [closest[to_try]]
                system_points = [self.points[new_points[0]] - self.points[p]
                                 for p in new_points[1:] ]
                U,s,V = np.linalg.svd(system_points)
                if (min(singular_value / s[0] for singular_value in s[1:])
                    > SMALL_NUMBER):
                    plane_points = new_points[:]
                else:
                    # The new point was not linearly independent
                    pass
                # Increment the attempt index
                to_try += 1
            # Solve for the weights (hyperplane interpolation)
            system = np.concatenate((self.points[plane_points],
                            np.ones((len(plane_points),1))), axis=1)
            pt = np.concatenate((pt,[1.0]))
            weights = np.linalg.solve(system, pt)
            value = np.dot(self.values[plane_points],weights)
            response.append(value)
        # Return the final response
        return response


# ================================
#      MARS (fortran) Wrapper     
# ================================

# In order to use fmodpy to generate the appropriate wrapper, the
# following modifications need to be made to 'marspack.f' before
# calling fmodpy:

# COMMENT OUT LINES 551 - 552, INSERT:
#      real, intent(out) :: fm(*)                                            1.1
#      integer, intent(out) :: im(*)                                         1.2
#      integer p,lx(p),mm(*)                                                 2
#      real x(n,p),y(n),w(n),sp(*)                                           3

# COMMENT OUT LINE 612 (+4), INSERT:
#      real, intent(out) :: f(*)                                            62.5
#      real x(n,:),fm(*),sp(:,2)                                            63

# Afterwards, these lines need to be reverted to their original form.
# Finally, call "Make" and the final correct output should be produced.

# If you want to make further calls to mars custom code, include them
# in the file 'marspack_wrapper.f03'

# Wrapper class for using the MARS fortran code
class FortMARS(Approximator):
    from .marspack import mars, fmod
    def __init__(self):
        self.fm = self.im = None

    # Use fortran code to compute the MARS model for given data
    #  See 'marspack.f' for details as to what these poorly named variables mean.
    def fit(self, control_points, values, max_bases=MARS_MAX_BASES,
            max_interaction=MARS_MAX_INTERACTION):
        # Normalize the control points to the unit hypercube
        control_points = control_points.copy()
        self.shift = np.min(control_points,axis=0)
        control_points -= self.shift
        self.scale = np.max(control_points,axis=0)
        control_points /= self.scale
        # Translated input parameters
        n = control_points.shape[0]
        p = control_points.shape[1]
        x = np.asfortranarray(control_points, dtype=float)
        y = np.asfortranarray(values, dtype=float)
        # Customizable parameters
        nk = min(n,max_bases)
        mi = min(p,max_interaction)
        # Hard coded parameters
        w = np.ones(shape=(n,), dtype=float, order="F") # Equal weighting
        lx = np.ones(shape=(p,), dtype=int, order="F")  # All ordinal
        nmcv = ntcv = 0 # No categoricals
        # Model storage variables (copied straight from 'marspack.f')
        fm_shape = 3+nk*(5*mi+nmcv+6)+2*p+ntcv
        im_shape = 21+nk*(3*mi+8)
        self.fm = np.ones(shape=(fm_shape,), dtype=float, order="F")
        self.im = np.ones(shape=(im_shape,), dtype=int,   order="F")
        # Runtime storage for MARS to execute (see 'marspack.f')
        sp_shape = n*(max(nk+1,2)+3)+max(3*n+5*nk+p,2*p,4*n)+2*p+4*nk
        dp_shape = max(n*nk,(nk+1)*(nk+1))+max((nk+2)*(nmcv+3),4*nk)
        mm_shape = n*p+2*max(mi,nmcv)
        sp = np.ones(shape=(sp_shape,), dtype=float, order="F")
        dp = np.ones(shape=(dp_shape,), dtype=float, order="F")
        mm = np.ones(shape=(mm_shape,), dtype=int,   order="F")
        # Make the call to MARS, in-place update to self.fm and self.im
        FortMARS.mars(n,p,x,y,w,nk,mi,lx,self.fm,self.im)#,sp,dp,mm)
        # NOTE: Any custom "call" parameters are hard-coded in 'marspack_wrapper.f03'

    # Use fortran code to evaluate the computed boxes
    def predict(self, x, model_flag=MARS_MODEL_FLAG):
        # Normalize the incoming x-points
        x -= self.shift
        x /= self.scale
        # Calculate all of the response values
        m = model_flag
        n = x.shape[0]
        x = np.asfortranarray(x,  dtype=float)
        f = np.ones(shape=(n,),   dtype=float, order="F")
        sp = np.ones(shape=(n,2), dtype=float, order="F")
        FortMARS.fmod(m,n,x,self.fm,self.im,f,sp)
        # Return the response values
        return f

# train = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/3-F_Size--R_Size--Freq/10.42--89.58/Train-0.csv"
# test = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/3-F_Size--R_Size--Freq/10.42--89.58/Test-0.csv"
# # ^^ FAILURE CASE

# The following can be added to "marspack.pyx" for testing

# # Extra variables for printing test fortran code
# cdef long nmcv = 0
# cdef long ntcv = 0
# print("")
# print("PROGRAM TEST_MARS")
# print("  IMPLICIT NONE")
# print("  INTEGER :: N, P, NK, MI")
# print("  REAL,             DIMENSION(%i,%i) :: X"%(n,p))
# print("  REAL,             DIMENSION(%i)    :: Y"%(n))
# print("  REAL,             DIMENSION(%i)    :: W"%(n))
# print("  INTEGER,          DIMENSION(%i)    :: LX"%(p))
# print("  REAL,             DIMENSION(%i)    :: FM"%(fm_0))
# print("  INTEGER,          DIMENSION(%i)    :: IM"%(im_0))
# print("  REAL,             DIMENSION(%i)    :: SP"%(sp_0))
# print("  DOUBLE PRECISION, DIMENSION(%i)    :: DP"%(dp_0))
# print("  INTEGER,          DIMENSION(%i)    :: MM"%(mm_0))
# print("  ")
# print("  N = %i    ! Number of points"%(n))
# print("  P = %i    ! Number of dimensions"%(p))
# print("  NK = %i   ! Maximum number of basis functions"%(nk))    
# print("  MI = %i   ! Maximum interaction between dimensions"%(mi))
# print("  W = 1.0   ! Weights for each point")
# print("  LX = 1    ! Flags for each dimension (1 = ordinal)")
# print("  FM = 1.0  ! Holder for MARS model")
# print("  IM = 1    ! Holder for MARS model")
# print("  SP = 1.0  ! Real workspace array for MARS")
# print("  DP = 1.0  ! Double precision workspace array for MARS")
# print("  MM = 1    ! Integer workspace array for MARS")
# print("  ")
# print("  ! Size of fm = 3+nk*(5*mi+nmcv+6)+2*p+ntcv")
# print("  !            = 3+%i*(5*%i+%i+6)+2*%i+%i)"%(nk,mi,nmcv,p,ntcv))
# print("  !            = %i"%(3+nk*(5*mi+nmcv+6)+2*p+ntcv))
# print("  ! Size of im = 21+nk*(3*mi+8)")
# print("  !            = 21+%i*(3*%i+8)"%(nk,mi))
# print("  !            = %i"%(21+nk*(3*mi+8)))
# print("  ! Size of sp = n*(max(nk+1,2)+3)+max(3*n+5*nk+p,2*p,4*n)+2*p+4*nk")
# print("  !            = n*(max(%i+1,2)+3)+max(3*%i+5*%i+%i,2*%i,4*%i)+2*%i+4*%i"%(nk,n,nk,p,p,n,p,nk))
# print("  !            = %i"%(n*(max(nk+1,2)+3)+max(3*n+5*nk+p,2*p,4*n)+2*p+4*nk))
# print("  ! Size of dp = max(n*nk,(nk+1)*(nk+1))+max((nk+2)*(nmcv+3),4*nk)")
# print("  !            = max(%i*%i,(%i+1)*(%i+1))+max((%i+2)*(%i+3),4*%i)"%(n,nk,nk,nk,nk,nmcv,nk))
# print("  !            = %i"%(max(n*nk,(nk+1)*(nk+1))+max((nk+2)*(nmcv+3),4*nk)))
# print("  ! Size of mm = n*p+2*max(mi,nmcv)")
# print("  !            = %i*%i+2*max(%i,%i)"%(n,p,mi,nmcv))
# print("  !            = %i"%(n*p+2*max(mi,nmcv)))
# print("  ")
# for step in range(1,n+1):
#     print("  X(%i,:) = (/ %s /)"%(step,", ".join(list(map(str,x[step-1,:])))))
# print("  ")
# for step in range(1,n+1):
#     print("  Y(%i) = %s"%(step,y[step-1]))
# print("  ")
# print("  CALL MARS(N,P,X,Y,W,NK,MI,LX,FM,IM,SP,DP,MM)")
# print("END PROGRAM TEST_MARS")
# print("")



# ==========================================
#      Bootstrapped Box Splines Wrapper     
# ==========================================

# Wrapper class for using the box_spline_basis fortran module
class BBS(Approximator):
    from .bootstrapped_box_splines import compute_boxes, eval_boxes
    def __init__(self):
        self.boxes = self.widths = self.weights = None
        self.shift = self.scale = None

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values,
            num_boxes=BBS_DEFAULT_NUM_BOXES,
            batch_size=BBS_DEFAULT_BATCH_SIZE):
        # Process and store local information
        x_data = np.asfortranarray(control_points.T, dtype=float)
        response = np.asfortranarray(values)
        num_boxes = min(num_boxes, x_data.shape[1])
        self.boxes = np.ones(shape=(x_data.shape[0],num_boxes), dtype=float, order="F")
        self.widths = np.ones(shape=(self.boxes.shape[1],),     dtype=float, order="F")
        self.weights = np.ones(shape=(self.boxes.shape[1],),    dtype=float, order="F")
        # In-place update of boxes, widths, and weights
        BBS.compute_boxes(x_data, response, self.boxes,
                          self.widths, self.weights, batch_size)

    # Use fortran code to evaluate the computed boxes
    def predict(self, x):
        x = np.array(x.T, order="F")
        response = np.ones(shape=(x.shape[1],), dtype=float, order="F")
        BBS.eval_boxes(self.boxes, self.widths, self.weights,
                       x, response)
        return response

# ============================
#      FitBoxMesh Wrapper     
# ============================

# Wrapper class for using the box_spline_basis fortran module
class FitBoxMesh(Approximator):
    from .fit_box_mod import compute_boxes, eval_boxes
    def __init__(self):
        self.boxes = self.widths = self.weights = None
        self.shift = self.scale = None

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values, num_boxes=FBM_DEFAULT_NUM_BOXES):
        # Process and store local information
        x_data = np.asfortranarray(control_points.T, dtype=float)
        response = np.asfortranarray(values)
        num_boxes = min(num_boxes, x_data.shape[1])
        self.boxes = np.ones(shape=(x_data.shape[0],num_boxes), dtype=float, order="F")
        self.widths = np.ones(shape=(self.boxes.shape[1],),     dtype=float, order="F")
        self.weights = np.ones(shape=(self.boxes.shape[1],),    dtype=float, order="F")
        # In-place update of boxes, widths, and weights
        FitBoxMesh.compute_boxes(x_data, response, self.boxes,
                                 self.widths, self.weights)

    # Use fortran code to evaluate the computed boxes
    def predict(self, x):
        x = np.array(x.T, order="F")
        response = np.ones(shape=(x.shape[1],), dtype=float, order="F")
        FitBoxMesh.eval_boxes(self.boxes, self.widths, self.weights,
                              x, response)
        return response


# ============================
#      MaxBoxMesh Wrapper     
# ============================

# Function for creating a surface over a region defined by the
# axis-aligned bounding box of all data points
class MaxBoxMesh(Approximator):
    from .max_box_mod import max_boxes, linear_eval

    def __init__(self):
        self.boxes = None
        self.values = None
        self.box_widths = None

    # Given a set of control points, this constructs the maximum box
    # mesh, that is to say the boxes around each control point have
    # the largest minimum side length possible.
    # Computational complexity: O(n^2*log(n) * d)
    def fit(self, control_points, values, padding=MBM_DEFAULT_PADDING):
        # Store the box-centers and the associated values
        self.boxes = np.asarray(control_points.T, dtype=float, order="F")
        self.values = np.asarray(values, dtype=float, order="F")
        # Initialize storage for the box-widths
        self.box_widths = np.ones(
            shape=(self.boxes.shape[0]*2,self.boxes.shape[1]),
            dtype=float, order="F")
        # Compute the max-boxes in-place
        MaxBoxMesh.max_boxes(self.boxes, self.box_widths)
        # Add a percentage of padding to all of the box widths
        self.box_widths += self.box_widths * padding

    # Function for displaying useful information
    def debug(self, x_pt):
        dists = np.sum((self.boxes.T - x_pt)**2,axis=1)**(1/2)
        orig_indices = np.argsort(dists)
        boxes = (self.boxes.T)[orig_indices]
        values = self.values[orig_indices]
        widths = (self.box_widths.T)[orig_indices]
        low_widths = widths[:,:widths.shape[1]//2]
        upp_widths = widths[:,widths.shape[1]//2:]
        print()
        print("BAD POINT:")
        print(CLEAN_ARRAY_STRING(x_pt))
        print()
        print("BOXES CONTAINING POINT: (center, lower, upper, shifted)")
        for i,(b,v,l,u) in enumerate(zip(boxes,values,low_widths,upp_widths)):
            if i >= MBM_CLOSE_TO_SHOW: break
            in_low = np.any( [(b-l < x_pt),(l < 0)], axis=0)
            in_upp = np.any( [(b+u > x_pt),(u < 0)], axis=0)
            in_box = all(np.all((in_low,in_upp),axis=0))
            print("  NUMBER %i"%(orig_indices[i]+1),v)
            print(CLEAN_ARRAY_STRING(b))
            print(CLEAN_ARRAY_STRING(l))
            print(CLEAN_ARRAY_STRING(u))
            shifted_box = np.where(b > x_pt, (b - x_pt) / l,
                                   (x_pt - b) / u)
            shifted_box = np.where(shifted_box < 0, 0.0, shifted_box)
            print(CLEAN_ARRAY_STRING(shifted_box))
            # print(" ",in_low)
            # print(" ",in_upp)
            print(" ",in_box)
            print()


    # Evaluate all of the box splines and return the final function
    # value at the requested point.
    def predict(self, x):
        # Transform x to be the right shape (and fortran ordered)
        x = np.asarray(x.T, dtype=float, order="F")
        response = np.ones(shape=(x.shape[1]), dtype=float, order="F")
        # Compute the linear boxes in-place, error is the last return value
        error = MaxBoxMesh.linear_eval(self.boxes, self.box_widths,
                                       self.values, x, response)[-1]
        if error > 0:
            x_pt = x.T[error-1]
            self.debug(x_pt)
            raise(Exception("ERROR: Could not compute boxes for these x points."))
        return response
