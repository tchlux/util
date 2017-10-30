# Makes available all custom algorithms that I have worked with
import numpy as np

CLEAN_ARRAY_STRING = lambda arr: " "+str(arr).replace(",","").replace("[","").replace("]","")
DEFAULT_RESPONSE = -1 #None
SMALL_NUMBER = 1.4901161193847656*10**(-8) 
#              ^^ SQRT(EPSILON( 0.0_REAL64 ))

NN_DEFAULT_NUM_NEIGHBORS = 1

MARS_MAX_BASES = float('inf')
MARS_MAX_INTERACTION = float('inf')
MARS_MODEL_FLAG = 1
# ^^ 1 -> Piecewise linear model
#    2 -> Piecewise cubic model

BBS_DEFAULT_BATCH_SIZE = None
BBS_DEFAULT_NUM_BOXES = float("inf")
FBM_DEFAULT_NUM_BOXES = float("inf")

MBM_CLOSE_TO_SHOW = 10
MBM_DEFAULT_PADDING = 0.01

BT_NUM_TREES = 100
# Linear complexity increase, ~2.8 seconds per 100 trees

DT_NUM_PARTS = 500
DT_MODEL = "linear"

TGP_BURNIN_TOTAL_ENERGY = (20, 120, 1)

# Get the name of a class as listed in python source file.
CLASS_NAME = lambda obj: (repr(obj)[1:].split(" ")[0]).split(".")[-1]

# Generic class definition for creating an algorithm that can perform
# regression in the multi-variate setting. In these classes 'x' refers
# to a potentially multi-dimensional set of euclydian coordinates with
# associated single-dimensional response values referred to as 'y'.
# 
# The "fit" function is given both 'x' and associated 'y' and creates
# any necessary (or available) model for predicting 'y' at new 'x'.
# 
# The "predict" function predicts 'y' at potentially novel 'x' values.
class Approximator:
    def __init__(self):
        pass

    def fit(self, x:np.ndarray, y:np.ndarray):
        pass

    def predict(self, x:np.ndarray, **kwargs):
        pass
    
    # Wrapper for 'predict' that returns a single value for a single
    # prediction, or an array of values for an array of predictions
    def __call__(self, x:np.ndarray, *args, **kwargs):
        single_response = len(x.shape) == 1
        if single_response:
            x = np.array([x])
        if len(x.shape) != 2:
            raise(Exception("ERROR: Bad input shape."))
        response = np.asarray(self.predict(x, *args, **kwargs), dtype=float)
        # Return the response values
        return response[0] if single_response else response


# ====================================================
#      Simple Linear Model Implemented with Numpy     
# ====================================================

# Class for creating a linear fit of D-dimensional data (hyperplane)
class LinearModel(Approximator):
    def __init__(self):
        self.model = None

    # Fit a linear model to the data
    def fit(self, control_points, values):
        # Process and store local information
        ones_column = np.ones( (control_points.shape[0],1) )
        coef_matrix = np.concatenate((control_points, ones_column), axis=1)
        # Returns (model, residuals, rank, singular values), we want model
        self.model, self.residuals = np.linalg.lstsq(coef_matrix, values) [:2]

    # Use fortran code to evaluate the computed boxes
    def predict(self, x):
        # Use the linear model to produce an estimate
        x = np.concatenate((x,np.ones( (x.shape[0],1) )), axis=1)
        return np.dot(x, self.model)


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


# ===========================================
#      Simple Nearest Neighbor Regressor     
# ===========================================

# Class for computing an interpolation between the nearest n neighbors
class NearestNeighbor(Approximator):
    def __init__(self):
        self.points = None
        self.values = None
        self.num_neighbors = None

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values,
            num_neighbors=NN_DEFAULT_NUM_NEIGHBORS):
        # Process and store local information
        self.points = control_points.copy()
        self.values = values.copy()
        self.num_neighbors = num_neighbors

    # Use fortran code to evaluate the computed boxes
    def predict(self, x):
        # Use the nearest point 
        response = []
        for pt in x:
            distances = np.sum((self.points - pt)**2, axis=1)
            closest = np.argsort(distances)[:self.num_neighbors]
            distances = distances[closest]
            if np.min(distances) <= 0:
                response.append( self.values[closest][np.argmin(distances)] )
            else:
                distances = 1 / distances
                weights = distances / sum(distances)
                response.append( np.sum(weights * self.values[closest]) )
        return response


# =======================
#      qHullDelaunay     
# =======================

# Wrapper class for using the Delaunay scipy code
class qHullDelaunay(Approximator):
    from scipy.spatial import Delaunay as Delaunay_via_qHull
    # For 1-D case, use regular interpolation
    from scipy import interpolate

    def __init__(self):
        self.pts = None
        self.values = None

    # Use fortran code to compute the boxes for the given data
    def fit(self, x, y):
        self.pts = x.copy()
        self.values = y.copy()
        # Handle special case for 1-D data, use regular linear
        #  interpolation, otherwise use standard Delaunay
        if (self.pts.shape[1] > 1):
            self.surf = qHullDelaunay.Delaunay_via_qHull(self.pts)
        else:
            self.surf = qHullDelaunay.interpolate.interp1d(self.pts[:,0],
                                                           self.values,
                                                           fill_value="extrapolate")

    # Use scipy code in order to evaluate delaunay triangulation
    def predict(self, x, debug=False, verbose=False):
        # Compute the response values
        response = []
        for i,x_pt in enumerate(x):
            # Regular linear interpolation for 1-D data
            if (self.pts.shape[1] == 1):
                value = self.surf(x_pt[0])
                response.append(value)
                continue
            else:
                # Solve for the weights in the old Delaunay model
                simp_ind = self.surf.find_simplex(x_pt)
                # If a point is outside the convex hull, use the
                # closest simplex to extrapolate the value
                if simp_ind < 0: 
                    # Calculate the distance between the x_pt each simplex
                    simp_dists = self.surf.plane_distance(x_pt)
                    # Find the index of the closest simplex
                    simp_ind = np.argmax(simp_dists)
                # Solve for the response value
                simp = self.surf.simplices[simp_ind]
                system = np.concatenate((self.pts[simp],
                    np.ones((simp.shape[0],1))), axis=1).T
                x_pt = np.concatenate((x_pt,[1]))
                weights = np.linalg.solve(system, x_pt)
                value = np.dot(self.values[simp],weights)
                response.append(value)
                if debug:
                    if verbose:
                        print("QHull")
                        print("Simplex:", sorted(simp))
                        print("Training Points:", self.pts.shape)
                        # print(CLEAN_ARRAY_STRING(self.pts))
                        print("Test Point:")
                        print(CLEAN_ARRAY_STRING(x_pt[:-1]))
                        print("Weights:")
                        print(CLEAN_ARRAY_STRING(weights[:-1]))
                        print("Value:  ", value)
                    return [simp]
        return response


# ============================
#      VTDelaunay Wrapper     
# ============================

# Wrapper class for using the Delaunay fortran code
class Delaunay(Approximator):
    from .fmodpy_VTdelaunay.VTdelaunay import delaunayp
    def __init__(self):
        self.pts = None
        self.values = None
        self.errs = {}

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values=None):
        self.pts = np.asarray(control_points.T, dtype=np.float64, order="F")
        # Store values if they were given
        if type(values) != type(None):
            self.values = np.asarray(values.reshape((1,values.shape[0])),
                                     dtype=np.float64, order="F")
        else:
            self.values = None

    # Return just the points and the weights associated with them for
    # creating the correct interpolation
    def points_and_weights(self, x):
        single_response = len(x.shape) == 1
        if single_response:
            x = np.array([x])
        if len(x.shape) != 2:
            raise(Exception("ERROR: Bad input shape."))
        points = []
        weights = []
        for i,x_pt in enumerate(x):
            x_pt = np.asarray(x_pt, dtype=float, order="F")
            simp, w, err, _ = Delaunay.delaunayp(self.pts.shape[0],
                                                 self.pts.shape[1],
                                                 self.pts, x_pt)
            if err != 0: print("Delaunay error:", err)
            points.append(self.pts.T[simp-1])
            weights.append(w)

        # Conver the lists of points and weights into numpy arrays
        points = np.asarray(points, dtype=float)
        weights = np.asarray(weights, dtype=float)
        # Return the appropriate shaped pair of points and weights
        return (points[0], weights[0]) if single_response else (points, weights)
        

    # Use fortran code to evaluate the computed boxes
    def predict(self, x, debug=False, verbose=False):
        if type(self.values) == type(None):
            raise(Exception("ERROR: Must provide values in order to get predictions."))
        # Calculate all of the response values
        response = []
        # Get the predictions from VTdelaunay
        p_in = np.asarray(x.T, dtype=np.float64, order="F")
        work_in = np.ones(shape=(p_in.shape[0]*p_in.shape[1],), 
                          dtype=np.float64, order="F")
        simp_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                           dtype=np.int32, order="F")
        weights_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                              dtype=np.float64, order="F")
        error_out = np.ones(shape=(p_in.shape[1],), 
                            dtype=np.int32, order="F")
        interp_in = self.values
        interp_out = np.ones(shape=(self.values.shape[0], p_in.shape[1]),
                             dtype=np.float64, order="F")
        Delaunay.delaunayp(self.pts.shape[0], self.pts, p_in, work_in,
                           simp_out, weights_out, error_out, 
                           extrap_opt=1000,
                           interp_in_opt=interp_in,
                           interp_out_opt=interp_out)
        return interp_out[0,:]

        # for i,x_pt in enumerate(x):
        #     x_pt = np.asarray(x_pt, dtype=float, order="F")
        #     simp, weights, err, _ = Delaunay.delaunayp(self.pts.shape[0],
        #                                                self.pts.shape[1],
        #                                                self.pts, x_pt)
        #     # Assume execution went well, collect values
        #     resp = np.sum( self.values[simp-1] * weights )
        #     response.append(resp)
        #     self.errs[err] = self.errs.get(err,0) + 1
        #     # Print some debug statements
        #     if debug:
        #         simp = simp - 1
        #         if verbose:
        #             print("Delaunay")
        #             print("Simplex:", sorted(simp))
        #             print("Training Points:", self.pts.shape)
        #             # print(Delaunay.CLEAN_ARRAY_STRING(self.pts.T))
        #             print("Test Point:")
        #             print(CLEAN_ARRAY_STRING(x_pt))
        #             print("Weights:")
        #             print(CLEAN_ARRAY_STRING(weights))
        #             print("Value:  ", resp)
        #             print("Err:    ", err)
        #         return [simp]
        # return response


# ==============================
#      MARS (Earth) Wrapper     
# ==============================

# Was using the fortran, but found a bug... Using Py-Earth instead

""" In order to download and install py-earth:
git clone git://github.com/scikit-learn-contrib/py-earth.git
cd py-earth
sudo python setup.py install
"""

# Wrapper class for using the open source MARS, py-earth
class MARS(Approximator):
    from pyearth import Earth
    def __init__(self):
        self.model = None

    # Use fortran code to compute the MARS model for given data
    #  See 'marspack.f' for details as to what these poorly named variables mean.
    def fit(self, control_points, values, max_bases=MARS_MAX_BASES,
            max_interaction=MARS_MAX_INTERACTION):
        max_bases = min(control_points.shape[0],max_bases)
        max_interaction = min(control_points.shape[1],max_interaction)
        self.model = MARS.Earth(max_terms=max_bases,
                                max_degree=max_interaction,
                                use_fast=True)
        self.model.fit(control_points, values)

    # Use Earth to evaluate the computed boxes
    def predict(self, x):
        return self.model.predict(x)


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
    from .fmodpy_marspack.marspack import mars, fmod
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


# =======================
#      LSHEP Wrapper     
# =======================

# Wrapper class for using the LSHEP fortran code
class LSHEP(Approximator):
    from .fmodpy_linear_shepard.linear_shepard import lshep, lshepval
    ierrors = {}
    def __init__(self):
        self.x = self.f = self.a = self.rw = None

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values):
        # Local operations
        self.x = np.asfortranarray(control_points.T)
        self.f = np.asfortranarray(values)
        self.a = np.ones(shape=self.x.shape, order="F")
        self.rw = np.ones(shape=(self.x.shape[1],), order="F")
        # In-place update of self.a and self.rw
        a,rw,ier = LSHEP.lshep(self.x.shape[0], self.x.shape[1],
                               self.x, self.f, self.a, self.rw)

    # Use fortran code to evaluate the computed boxes
    def predict(self, x):
        # Calculate all of the response values
        response = []
        for x_pt in x:
            x_pt = np.array(np.reshape(x_pt,(self.x.shape[0],)), order="F")
            ier, resp = LSHEP.lshepval(x_pt, self.x.shape[0],
                                       self.x.shape[1], self.x,
                                       self.f, self.a, self.rw)
            response.append(resp)
            LSHEP.ierrors[ier] = LSHEP.ierrors.get(ier,0) + 1
        # Return the response values
        return response

    # Print and return a summary of the errors experienced
    def errors(self):
        print("%i normal returns."%LSHEP.ierrors.get(0,0))
        print(("%i points outside the radius of influence "+
              "of all other points.")%LSHEP.ierrors.get(1,0))
        print(("%i points at which hyperplane normals are significantly"+
               " different.\n  This is only hazardous for piecewise "+
               "linear underlying functions.")%LSHEP.ierrors.get(2,0))
        return LSHEP.ierrors.get(0,0), LSHEP.ierrors.get(1,0), LSHEP.ierrors.get(2,0)


# ==========================================
#      Bootstrapped Box Splines Wrapper     
# ==========================================

# Wrapper class for using the box_spline_basis fortran module
class BBS(Approximator):
    from .fmodpy_bootstrapped_box_splines.bootstrapped_box_splines import compute_boxes, eval_boxes
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
    from .fmodpy_fit_box_mod.fit_box_mod import compute_boxes, eval_boxes
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
    from .fmodpy_max_box_mod.max_box_mod import max_boxes, linear_eval

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


# ===============================
#      BayesTree (R) Wrapper     
# ===============================

class BayesTree(Approximator):
    # For managing and executing R codes, as well as importing packages
    from rpy2.robjects.vectors import FloatVector
    from rpy2.robjects.packages import importr
    from rpy2.rinterface import RRuntimeError
    from rpy2.robjects import r
    from rpy2.robjects.numpy2ri import numpy2ri

    def __init__(self):
        self.train_x = None
        self.train_y = None
        self.BayesTree = None
        self.ntree = None

    # Simply store the data for later evaluation
    def fit(self, control_points, values, num_trees=BT_NUM_TREES):
        # Local operations
        self.train_x = BayesTree.r.matrix(BayesTree.numpy2ri(control_points),
                                          nrow=control_points.shape[0],
                                          ncol=control_points.shape[1])
        self.train_y = BayesTree.FloatVector(values)
        self.BayesTree = BayesTree.importr('BayesTree')
        self.ntree = num_trees

    # Fit the training data and get the returned values
    def predict(self, x):
        extra = False
        # BayesTree cannot handle single test points, for some reason
        # there must be at least 2 points. This adds an extra null point.
        if x.shape[0] == 1:
            extra = True
            x = np.concatenate((np.ones(shape=(1,x.shape[1])),x))
        # Convert the x input to an R data type
        x = BayesTree.r.matrix(BayesTree.numpy2ri(x),
                               nrow=x.shape[0], ncol=x.shape[1])
        # Calculate the response for the given x with R
        try:
            response = self.BayesTree.bart(self.train_x, self.train_y, x,
                                           verbose=False, ntree=self.ntree)
            # The 7th return value is the estimated response
            response = np.asarray(response[7], dtype=float)
            # except:
            #     response = np.array([DEFAULT_RESPONSE] * x.dim[0])
        except BayesTree.RRuntimeError:
            response = [DEFAULT_RESPONSE] * x.dim[0]
        # Remove the extra point if it was added
        if extra: response = response[-1:]
        return response

# WARNING: Example test data that fails
# train_x = [[0, 0],
#            [0, 1],
#            [1, 1]]
# train_y = [<any 3 vector>]
# test_x = [[0.2, 0.0],
#           [0.2, 0.2],
#           [1.0, 0.2]]

# I have not found a way around this bug, return DEFAULT_RESPONSE instead


# ==============================
#      dynaTree (R) Wrapper     
# ==============================

class dynaTree(Approximator):
    # For managing and executing R codes, as well as importing packages
    from rpy2.robjects.vectors import FloatVector
    from rpy2.robjects.packages import importr
    from rpy2.robjects import r

    def __init__(self):
        self.dT = None
        self.dynaTree = None

    # Simply store the data for later evaluation
    def fit(self, control_points, values, num_parts=DT_NUM_PARTS, model=DT_MODEL):
        # Local operations
        X = dynaTree.r.matrix(control_points,
                              nrow=control_points.shape[0],
                              ncol=control_points.shape[1])
        y = dynaTree.FloatVector(values)
        self.dT = dynaTree.importr('dynaTree')
        self.dynaTree = self.dT.dynaTree(X,y,N=num_parts, model=model,
                                         verb=0, minp=len(y)//2)

    # Fit the training data and get the returned values
    def predict(self, x):
        # Convert the x input to an R data type
        x = dynaTree.r.matrix(x, nrow=x.shape[0], ncol=x.shape[1])
        # Calculate the response for the given x with R
        response = self.dT.predict_dynaTree(self.dynaTree, x)
        # The 12th return value is the estimated response
        response = np.asarray(response[12], dtype=float)
        # Return the response value(s)
        return response


# =========================
#      tgp (R) Wrapper     
# =========================

class tgp(Approximator):
    # For managing and executing R codes, as well as importing packages
    from rpy2.robjects.vectors import FloatVector
    from rpy2.robjects.packages import importr
    from rpy2.robjects import r

    def __init__(self):
        self.train_x = None
        self.train_y = None
        self.tgp = None
        self.bte = None

    # Simply store the data for later evaluation
    def fit(self, control_points, values, bte=TGP_BURNIN_TOTAL_ENERGY):
        # Local operations
        self.train_x = tgp.r.matrix(control_points,
                                    nrow=control_points.shape[0],
                                    ncol=control_points.shape[1])
        self.train_y = tgp.FloatVector(values)
        self.tgp = importr('tgp')
        self.bte = bte
        self.step = 0

    # Fit the training data and get the returned values
    def predict(self, x):
        # Convert the x input to an R data type
        x = tgp.r.matrix(x, nrow=x.shape[0], ncol=x.shape[1])
        # Calculate the response for the given x with R
        bte = tgp.FloatVector(self.bte)
        response = self.tgp.btgp(self.train_x, self.train_y, x,
                                 BTE=bte, verb=0)
        # response = self.tgp.btlm(self.train_x, self.train_y, x,
        #                          BTE=bte, verb=0)
        # response = self.tgp.bcart(self.train_x, self.train_y, x,
        #                           BTE=bte, verb=0)

        # The 15th return value is the estimated response
        return np.asarray(response[15], dtype=float)


# =================================================
#           Optimization and Minimization          
# =================================================


DEFAULT_SEARCH_SIZE = 2.0**10
DEFAULT_MAX_TIME_SEC = 5
DEFAULT_MIN_STEPS = 10
DEFAULT_MIN_IMPROVEMENT = 0.0001

import time
from util.algorithms.DiRect import DiRect as optimize

# Class for tracking convergence of optimization algorithm
class Tracker:
    # Initialization of storage attributes
    def __init__(self, objective, max_time=DEFAULT_MAX_TIME_SEC,
                 min_steps=DEFAULT_MIN_STEPS, 
                 min_improvement=DEFAULT_MIN_IMPROVEMENT, display=False):
        self.max_time = max_time
        self.min_steps = min_steps
        self.min_change = min_improvement
        self.display = display
        self.obj = objective
        self.max_obj = -float('inf')
        self.min_obj = float('inf')
        self.best_sol = None
        self.record = []
        self.tries = 0
        self.start = None

    # Wrapper for the objective function that tracks progress
    def check(self,sol):
        # Initialize starting time if it hasn't been done
        if self.start == None: self.start = time.time()
        # Don't execute the objective function if time has been exceeded
        if time.time() - self.start > DEFAULT_MAX_TIME_SEC: return float('inf')
        # Calculate the new objective value
        new_obj = self.obj(sol)
        self.tries += 1
        if new_obj > self.max_obj: self.max_obj = new_obj
        if new_obj <= self.min_obj:
            self.record.append(new_obj)
            self.min_obj = new_obj
            self.best_sol = sol
            if self.display: print(new_obj, sol)
        return new_obj

    # Function that returns "True" when optimization is over
    def done(self):
        # Most importantly, check for time
        time_done = (time.time() - self.start > self.max_time)
        if time_done: return time_done
        # Next, check for convergence
        converged = False
        if len(self.record) > self.min_steps:
            change = [abs(v - self.min_obj) for v in self.record[-self.min_steps:]]
            divisor = max(self.record) - min(self.record)
            if divisor == 0: 
                converged = False
            else:
                converged = max(change) / divisor < self.min_change
        return converged

# Minimization function
def minimize(objective, solution, bounds=None,
             max_time=DEFAULT_MAX_TIME_SEC,
             min_steps=DEFAULT_MIN_STEPS,
             min_improvement=DEFAULT_MIN_IMPROVEMENT, plot=False):
    solution = np.array(list(map(float,solution)))

    # Generate some arbitrary bounds
    if type(bounds) == type(None):
        upper = solution + np.array([DEFAULT_SEARCH_SIZE]*len(solution))
        lower = solution - np.array([DEFAULT_SEARCH_SIZE]*len(solution))
        bounds = list(zip(lower, upper))

    # Initialize a tracker for halting the optimization
    t = Tracker(objective, max_time, min_steps, min_improvement, plot)

    # Call the optimization function and get the best solution
    optimize(t.check, t.done, bounds, solution)

    if plot:
        import pygal
        CUSTOM_CSS_FILE = 'file://'+os.path.expanduser("~")+'/Git_Analytics/Analytics-Vis/Code/ViewData/Value_Distribution/no_dot.css'
        config = pygal.Config()

        config.css.append(CUSTOM_CSS_FILE)
        plot = pygal.XY(config)
        plot.add("Objective", list(zip(range(len(t.record)),t.record)))
        plot.render_in_browser()
        

    return t.best_sol
