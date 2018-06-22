DEFAULT_RESPONSE = -1 #None

# Bayes Tree
BT_NUM_TREES = 100
# Linear complexity increase, ~2.8 seconds per 100 trees

# DynaTree
DT_NUM_PARTS = 500
DT_MODEL = "linear"

# Treed Gaussian Process
TGP_BURNIN_TOTAL_ENERGY = (20, 120, 1)


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
