import numpy as np
from util.approximate import Approximator

# ====================================================
#      Simple Linear Model Implemented with Numpy     
# ====================================================

# Class for creating a linear fit of D-dimensional data (hyperplane)
class LinearModel(Approximator):
    def __init__(self):
        self.model = None

    # Fit a linear model to the data
    def _fit(self, control_points, values, weights=None):
        # Process and store local information
        ones_column = np.ones( (control_points.shape[0],1) )
        coef_matrix = np.concatenate((control_points, ones_column), axis=1)
        # Weight the points if necessary
        if type(weights) != type(None): coef_matrix[:,-1] *= weights
        # Returns (model, residuals, rank, singular values), we want model
        self.model, self.residuals = np.linalg.lstsq(coef_matrix, values, rcond=None)[:2]

    # Use fortran code to evaluate the computed boxes
    def _predict(self, x):
        if (len(x.shape) == 1): x = x[None,:]
        # Use the linear model to produce an estimate
        x = np.concatenate((x,np.ones( (x.shape[0],1) )), axis=1)
        return np.matmul(x, self.model)
