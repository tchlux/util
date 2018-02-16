# Makes available all custom algorithms that I have worked with
import os
import numpy as np
import fmodpy

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))

# Get the name of a class as listed in python source file.
CLASS_NAME = lambda obj: (repr(obj)[1:-2].split(".")[-1])
CLEAN_ARRAY_STRING = lambda arr: " "+str(arr).replace(",","").replace("[","").replace("]","")
DEFAULT_RESPONSE = -1 #None
SMALL_NUMBER = 1.4901161193847656*10**(-8) 
#              ^^ SQRT(EPSILON( 0.0_REAL64 ))

# Nearest neighbor
NN_DEFAULT_NUM_NEIGHBORS = 1

# Shepard Method
SHEPARD_WEIGHT_FUNCTION = lambda dist: 1 / max(dist,SMALL_NUMBER)

# MLP Regressor
MLP_R_ACTIVATION_FUNC = "relu"
MLP_R_SOLVER = "lbfgs"

# Support Vector Regressor
SVR_KERNEL = "poly"

# MARS
MARS_MAX_BASES = float('inf')
MARS_MAX_INTERACTION = float('inf')
MARS_MODEL_FLAG = 1
# ^^ 1 -> Piecewise linear model
#    2 -> Piecewise cubic model

# Bayes Tree
BT_NUM_TREES = 100
# Linear complexity increase, ~2.8 seconds per 100 trees

# DynaTree
DT_NUM_PARTS = 500
DT_MODEL = "linear"

# Treed Gaussian Process
TGP_BURNIN_TOTAL_ENERGY = (20, 120, 1)

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

    def __str__(self):
        return CLASS_NAME(self)

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


# ==========================================
#      Multi-Layer Perceptron Regressor     
# ==========================================
class MLPRegressor(Approximator):
    def __init__(self, *args, activation=MLP_R_ACTIVATION_FUNC,
                 solver=MLP_R_SOLVER, **kwargs):
        from sklearn.neural_network import MLPRegressor
        self.MLPRegressor = MLPRegressor
        kwargs.update(dict(activation=activation, solver=solver))
        self.mlp = self.MLPRegressor(*args, **kwargs)
    def fit(self, *args, **kwargs):
        return self.mlp.fit(*args, **kwargs)
    def predict(self, *args, **kwargs):
        return self.mlp.predict(*args, **kwargs)

    
# ==================================
#      Support Vector Regressor     
# ==================================
class SVR(Approximator):
    def __init__(self, *args, kernel=SVR_KERNEL, **kwargs):
        from sklearn.svm import SVR
        self.SVR = SVR
        kwargs.update(dict(kernel=kernel))
        self.svr = self.SVR(*args, **kwargs)
    def fit(self, *args, **kwargs):
        return self.svr.fit(*args, **kwargs)
    def predict(self, *args, **kwargs):
        return self.svr.predict(*args, **kwargs)

    
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
        self.model, self.residuals = np.linalg.lstsq(coef_matrix, values)[:2]

    # Use fortran code to evaluate the computed boxes
    def predict(self, x):
        # Use the linear model to produce an estimate
        x = np.concatenate((x,np.ones( (x.shape[0],1) )), axis=1)
        return np.dot(x, self.model)

    
# ===========================================
#      Simple Nearest Neighbor Regressor     
# ===========================================

# Class for computing an interpolation between the nearest n neighbors
class NearestNeighbor(Approximator):
    def __init__(self, num_neighbors=NN_DEFAULT_NUM_NEIGHBORS):
        self.points = None
        self.values = None
        self.num_neighbors = num_neighbors

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values, num_neighbors=None):
        if type(num_neighbors) != type(None):
            self.num_neighbors = num_neighbors
        # Process and store local information
        self.points = control_points.copy()
        self.values = values.copy()

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
    def __init__(self):
        from scipy.spatial import Delaunay as Delaunay_via_qHull
        # For 1-D case, use regular interpolation
        from scipy import interpolate
        # Store modules in self for later access
        self.Delaunay_via_qHull = Delaunay_via_qHull
        self.interpolate = interpolate
        
        self.pts = None
        self.values = None

    # Use fortran code to compute the boxes for the given data
    def fit(self, x, y):
        self.pts = x.copy()
        self.values = y.copy()
        # Handle special case for 1-D data, use regular linear
        #  interpolation, otherwise use standard Delaunay
        if (self.pts.shape[1] > 1):
            self.surf = self.Delaunay_via_qHull(self.pts)
        else:
            self.surf = self.interpolate.interp1d(self.pts[:,0],
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
    def __init__(self):
        self.delaunayp = fmodpy.fimport(
            os.path.join(CWD,"VTdelaunay","VTdelaunay.f90"),
            module_link_args=["-lgfortran","-lblas","-llapack"], 
            output_directory=CWD, ).delaunayp
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

        # Get the predictions from VTdelaunay
        p_in = np.asarray(x.T, dtype=np.float64, order="F")
        work_in = np.ones(shape=(max(5*p_in.shape[0],
                                     p_in.shape[0]*p_in.shape[1]),), 
                          dtype=np.float64, order="F")
        simp_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                           dtype=np.int32, order="F")
        weights_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                              dtype=np.float64, order="F")
        error_out = np.ones(shape=(p_in.shape[1],), 
                            dtype=np.int32, order="F")
        self.delaunayp(self.pts.shape[0], self.pts, p_in, work_in,
                       simp_out, weights_out, error_out, 
                       extrap_opt=1000)
        error_out = np.where(error_out != 1, error_out, 0)
        # Handle any errors that may have occurred.
        if (sum(error_out) != 0):
            unique_errors = sorted(np.unique(error_out))
            print(" [Delaunay errors:",end="")
            for e in unique_errors:
                if (e in {0,1}): continue
                print(" %3i"%e,"at",""+",".join(tuple(
                    str(i) for i in range(len(error_out))
                    if (error_out[i] == e)))+"}", end=";")
            print("] ",end="")
            # Reset the errors to simplex of 1s (to be 0) and weights of 0s.
            bad_indices = (error_out > 1)
            simp_out[:,bad_indices] = 1
            weights_out[:,bad_indices] = 0
        # Adjust the output simplices and weights to be expected shape
        points  = simp_out.T - 1
        weights = weights_out.T
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
        work_in = np.ones(shape=(max(5*p_in.shape[0],
                                     p_in.shape[0]*p_in.shape[1]),), 
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
        self.delaunayp(self.pts.shape[0], self.pts, p_in, work_in,
                       simp_out, weights_out, error_out, 
                       extrap_opt=1000,
                       interp_in_opt=interp_in,
                       interp_out_opt=interp_out)
        if debug:
            print("VTDelaunay")
            print("P Input:   ", p_in.T)
            print("Simp Out:  ", simp_out.T - 1)
            print("Error Out: ", error_out)
            print("Interp Out:", interp_out.T)
        return interp_out[0,:]


# ======================
#      Voronoi Mesh     
# ======================
class VoronoiMesh(Approximator):
    def __init__(self):
        self.voronoi_mesh = fmodpy.fimport(
            os.path.join(CWD,"voronoi_mesh","voronoi_mesh.f90"),
            module_link_args=["-lgfortran","-lblas","-llapack"], 
            output_directory=CWD)
        self.points = None
        self.values = None

    # Fit a set of points
    def fit(self, points, values=None):
        self.points = np.asarray(points.T.copy(), order="F")
        if type(values) != type(None):
            self.values = np.asarray(values.copy(), order="F")
        self.dots = np.ones((points.shape[0],)*2, order="F")
        # Compute all of the pairwise dot products between control points
        self.voronoi_mesh.train_vm(self.points, self.dots)

    # Calculate the weights associated with a set of points
    def points_and_weights(self, x, convex=True):
        single_response = len(x.shape) == 1
        if single_response:
            x = np.array([x])
        if len(x.shape) != 2:
            raise(Exception("ERROR: Bad input shape."))
        x = np.array(x.T, order="F")
        weights = np.ones((x.shape[1], self.points.shape[1]), 
                          dtype=np.float64, order="F")
        self.voronoi_mesh.predict_vm(self.points, self.dots,
                                     x, weights)
        if convex:
            # Make all the weights convex
            weights = (weights.T / np.sum(weights,axis=1)).T
        # Generate the list of point indices
        points = np.zeros(weights.shape, dtype=int) + np.arange(self.points.shape[1])
        # Return the appropriate shaped pair of points and weights
        return points, weights

    # Generate a prediction for a new point
    def predict(self, xs):
        if type(self.values) == type(None):
            raise(Exception("Must provide values in order to make predicitons."))
        weights = self.points_and_weights(xs)
        print(weights)
        # Generate the predictions
        predictions = np.matmul(weights, self.values)
        print(predictions)
        exit()
        return predictions

# =========================
#      MaxBoxMesh Mesh     
# =========================
class MaxBoxMesh(Approximator):
    def __init__(self):
        self.max_box_mesh = fmodpy.fimport(
            os.path.join(CWD,"max_box_mesh","max_box_mesh.f90"),
            module_link_args=["-lgfortran","-lblas","-llapack"], 
            output_directory=CWD)
        self.points = None
        self.values = None

    # Fit a set of points
    def fit(self, points, values=None):
        self.points = np.asarray(points.T.copy(), order="F")
        if type(values) != type(None):
            self.values = np.asarray(values.copy(), order="F")
        self.shapes = np.ones((self.points.shape[0]*2,self.points.shape[1]), order="F")
        # Compute all of the max-boxes
        self.max_box_mesh.max_boxes(self.points, self.shapes)
        # Expand all boxes (to make the mesh space-filling)
        self.shapes += SMALL_NUMBER

    # Calculate the weights associated with a set of points
    def points_and_weights(self, x):
        single_response = len(x.shape) == 1
        if single_response:
            x = np.array([x])
        if len(x.shape) != 2:
            raise(Exception("ERROR: Bad input shape."))

        x = np.array(x.T, order="F")
        weights = np.ones((x.shape[1], self.points.shape[1]), 
                          dtype=np.float64, order="F")
        self.max_box_mesh.eval_box_mesh(self.points, self.shapes,
                                        x, weights)

        # Make all weights convex
        for i,w in enumerate(np.sum(weights,axis=1)):
            weights[i,:] /= w
        points = np.zeros(weights.shape, dtype=int) + np.arange(self.points.shape[1])
        # Return the appropriate shaped pair of points and weights
        return points, weights

    # Generate a prediction for a new point
    def predict(self, xs):
        if type(self.values) == type(None):
            raise(Exception("Must provide values in order to make predicitons."))
        weights = self.points_and_weights(xs)
        predictions = np.matmul(weights, self.values)
        return predictions



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
    def __init__(self, max_bases=MARS_MAX_BASES,
                 max_interaction=MARS_MAX_INTERACTION):
        from pyearth import Earth
        self.Earth = Earth
        self.model = None
        self.max_bases = max_bases
        self.max_interaction = max_interaction

    # Use fortran code to compute the MARS model for given data
    #  See 'marspack.f' for details as to what these poorly named variables mean.
    def fit(self, control_points, values, max_bases=None, max_interaction=None):
        if type(max_bases) != type(None):
            self.max_bases = max_bases
        if type(max_interaction) != type(None):
            self.max_interaction = max_interaction
        max_bases = min(control_points.shape[0],self.max_bases)
        max_interaction = min(control_points.shape[1],self.max_interaction)
        self.model = self.Earth(max_terms=max_bases,
                                max_degree=max_interaction,
                                use_fast=True)
        self.model.fit(control_points, values)

    # Use Earth to evaluate the computed boxes
    def predict(self, x):
        return self.model.predict(x)

# =======================
#      LSHEP Wrapper     
# =======================

# Wrapper class for using the LSHEP fortran code
class LSHEP(Approximator):
    def __init__(self):
        self.linear_shepard = fmodpy.fimport(
            os.path.join(CWD,"linear_shepard","linear_shepard.f95"),
            module_link_args=["-lgfortran","-lblas","-llapack"], 
            output_directory=CWD)
        self.lshep = self.linear_shepard.lshep
        self.lshepval = self.linear_shepard.lshepval
        self.ierrors = {}
        self.x = self.f = self.a = self.rw = None

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values):
        # Local operations
        self.x = np.asfortranarray(control_points.T)
        self.f = np.asfortranarray(values)
        self.a = np.ones(shape=self.x.shape, order="F")
        self.rw = np.ones(shape=(self.x.shape[1],), order="F")
        # In-place update of self.a and self.rw
        self.lshep(self.x.shape[0], self.x.shape[1],
                   self.x, self.f, self.a, self.rw)

    # Use fortran code to evaluate the computed boxes
    def predict(self, x):
        # Calculate all of the response values
        response = []
        for x_pt in x:
            x_pt = np.array(np.reshape(x_pt,(self.x.shape[0],)), order="F")
            resp = self.lshepval(x_pt, self.x.shape[0],
                                 self.x.shape[1], self.x,
                                 self.f, self.a, self.rw)
            response.append(resp)
            # self.ierrors[ier] = self.ierrors.get(ier,0) + 1
        # Return the response values
        return response

    # Print and return a summary of the errors experienced
    def errors(self):
        print("%i normal returns."%self.ierrors.get(0,0))
        print(("%i points outside the radius of influence "+
              "of all other points.")%self.ierrors.get(1,0))
        print(("%i points at which hyperplane normals are significantly"+
               " different.\n  This is only hazardous for piecewise "+
               "linear underlying functions.")%self.ierrors.get(2,0))
        return self.ierrors.get(0,0), self.ierrors.get(1,0), self.ierrors.get(2,0)


# ===============================
#      BayesTree (R) Wrapper     
# ===============================

class BayesTree(Approximator):
    def __init__(self, num_trees=BT_NUM_TREES):
        try:    import rpy2
        except: raise(Exception("BayesTree not implemented, missing 'rpy2' package."))
        # For managing and executing R codes, as well as importing packages
        from rpy2.robjects.vectors import FloatVector
        from rpy2.robjects.packages import importr
        from rpy2.rinterface import RRuntimeError
        from rpy2.robjects import r
        from rpy2.robjects.numpy2ri import numpy2ri
        # Store modules in self
        self.FloatVector = FloatVector
        self.importr = importr
        self.RRuntimeError = RRuntimeError
        self.r = r
        self.numpy2ri = numpy2ri
        # Local variables
        self.train_x = None
        self.train_y = None
        self.BayesTree = None
        self.ntree = num_trees

    # Simply store the data for later evaluation
    def fit(self, control_points, values, num_trees=None):
        if type(num_trees) != type(None):
            self.ntree = num_trees
        # Local operations
        self.train_x = self.r.matrix(BayesTree.numpy2ri(control_points),
                                     nrow=control_points.shape[0],
                                     ncol=control_points.shape[1])
        self.train_y = self.FloatVector(values)
        self.BayesTree = self.importr('BayesTree')

    # Fit the training data and get the returned values
    def predict(self, x):
        extra = False
        # BayesTree cannot handle single test points, for some reason
        # there must be at least 2 points. This adds an extra null point.
        if x.shape[0] == 1:
            extra = True
            x = np.concatenate((np.ones(shape=(1,x.shape[1])),x))
        # Convert the x input to an R data type
        x = self.r.matrix(self.numpy2ri(x),
                          nrow=x.shape[0], ncol=x.shape[1])
        # Calculate the response for the given x with R
        try:
            response = self.BayesTree.bart(self.train_x, self.train_y, x,
                                           verbose=False, ntree=self.ntree)
            # The 7th return value is the estimated response
            response = np.asarray(response[7], dtype=float)
            # except:
            #     response = np.array([DEFAULT_RESPONSE] * x.dim[0])
        except self.RRuntimeError:
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
    def __init__(self, num_parts=DT_NUM_PARTS, model=DT_MODEL):
        try:    import rpy2
        except: raise(Exception("dynaTree not implemented, missing 'rpy2' package."))
        # For managing and executing R codes, as well as importing packages
        from rpy2.robjects.vectors import FloatVector
        from rpy2.robjects.packages import importr
        from rpy2.robjects import r
        # Store modules in self
        self.FloatVector = FloatVector
        self.importr = importr
        self.r = r
        # Local variables
        self.dT = None
        self.dynaTree = None
        self.num_parts = num_parts
        self.model = DT_MODEL

    # Simply store the data for later evaluation
    def fit(self, control_points, values, num_parts=None, model=None):
        if type(num_parts) != type(None):
            self.num_parts = num_parts
        if type(model) != type(None):
            self.model = model
        # Local operations
        X = self.r.matrix(control_points,
                          nrow=control_points.shape[0],
                          ncol=control_points.shape[1])
        y = self.FloatVector(values)
        self.dT = self.importr('dynaTree')
        self.dynaTree = self.dT.dynaTree(X,y,N=self.num_parts, model=self.model,
                                         verb=0, minp=len(y)//2)

    # Fit the training data and get the returned values
    def predict(self, x):
        # Convert the x input to an R data type
        x = self.r.matrix(x, nrow=x.shape[0], ncol=x.shape[1])
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
    def __init__(self, bte=TGP_BURNIN_TOTAL_ENERGY):
        try:    import rpy2
        except: raise(Exception("tgp not implemented, missing 'rpy2' package."))
        # For managing and executing R codes, as well as importing packages
        from rpy2.robjects.vectors import FloatVector
        from rpy2.robjects.packages import importr
        from rpy2.robjects import r
        # Store modules in self
        self.FloatVector = FloatVector
        self.importr = importr
        self.r = r
        # Local variables
        self.train_x = None
        self.train_y = None
        self.tgp = None
        self.bte = self.FloatVector(bte)

    # Simply store the data for later evaluation
    def fit(self, control_points, values, bte=None):
        if type(bte) != type(None):
            self.bte = self.FloatVector(bte)
        # Local operations
        self.train_x = self.r.matrix(control_points,
                                     nrow=control_points.shape[0],
                                     ncol=control_points.shape[1])
        self.train_y = self.FloatVector(values)
        self.tgp = self.importr('tgp')
        self.step = 0

    # Fit the training data and get the returned values
    def predict(self, x):
        # Convert the x input to an R data type
        x = self.r.matrix(x, nrow=x.shape[0], ncol=x.shape[1])
        # Calculate the response for the given x with R
        response = self.tgp.btgp(self.train_x, self.train_y, x,
                                 BTE=self.bte, verb=0)
        # response = self.tgp.btlm(self.train_x, self.train_y, x,
        #                          BTE=bte, verb=0)
        # response = self.tgp.bcart(self.train_x, self.train_y, x,
        #                           BTE=bte, verb=0)
        # The 15th return value is the estimated response
        return np.asarray(response[15], dtype=float)


if __name__ == "__main__":
    np.random.seed(0)

    # import time
    # print("Starting")
    # for n in range(10, 20011, 1000):
    #     n = 8000
    #     # for d in range(20,21,10):
    #     d = 5
    #     print("%6s"%n, "%3s"%d)
    #     x = np.random.random((n,d))
    #     y = np.random.random((n,))
    #     # model = NearestNeighbor(2)
    #     # model = LinearModel()
    #     # model = MARS()
    #     # model = LSHEP()
    #     # model = MLPRegressor()
    #     # model = VoronoiMesh()
    #     model = Delaunay()
    #     # model = MaxBoxMesh()
    #     # model = qHullDelaunay()
    #     start = time.time()
    #     model.fit(x, y)
    #     print("%5.2f"%(time.time() - start),end=" ")
    #     start = time.time()
    #     model(x[0])
    #     print("%3.2f"%((time.time() - start)*10000))
    #     print()
    #     exit()
    # exit()

    from util.plotly import Plot
    mult = 5
    fun = lambda x: np.cos(x[0]*mult) + np.sin(x[1]*mult)
    p = Plot()
    low = 0
    upp = 1
    dim = 2
    plot_points = 2000
    N = 2
    random = True
    if random:
        x = np.random.random(size=(N,dim))
    else:
        N = int(round(N ** (1/dim)))
        x = np.array([r.flatten() for r in np.meshgrid(np.linspace(low,upp,N), np.linspace(low,upp,N))]).T
    y = np.array([fun(v) for v in x])
    p.add("Training Points", *x.T, y)
    
    # surf = Delaunay()
    # surf = MaxBoxMesh()
    surf = VoronoiMesh()
    surf.fit(x,y)
    p.add_func(str(surf), surf, *([(low,upp)]*dim), plot_points=plot_points)
    p.plot(file_name="test_plot.html")

