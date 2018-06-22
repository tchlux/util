import numpy as np
from util.algorithms import WeightedApproximator

# Nearest neighbor
NN_DEFAULT_NUM_NEIGHBORS = 1
    
# ===========================================
#      Simple Nearest Neighbor Regressor     
# ===========================================

# Class for computing an interpolation between the nearest n neighbors
class NearestNeighbor(WeightedApproximator):
    def __init__(self, num_neighbors=NN_DEFAULT_NUM_NEIGHBORS):
        self.points = None
        self.num_neighbors = num_neighbors

    # Use fortran code to compute the boxes for the given data
    def _fit(self, control_points, num_neighbors=None):
        if type(num_neighbors) != type(None):
            self.num_neighbors = num_neighbors
        # Process and store local information
        self.points = control_points.copy()

    # Function that returns the indices of points and the weights that
    # should be used to make associated predictions for each point in
    # "points".
    def _predict(self, points):
        # Body
        pts = []
        wts = []
        for pt in points:
            distances = np.sum((self.points - pt)**2, axis=1)
            closest = np.argsort(distances)[:self.num_neighbors]
            distances = distances[closest]
            if np.min(distances) <= 0:
                pts.append( [closest[0]] )
                wts.append( [1.0] )
            else:
                # distances = 1 / distances
                # weights = distances / sum(distances)
                weights = np.ones(self.num_neighbors) / self.num_neighbors
                points = closest
                pts.append( closest )
                wts.append( weights )
        # Convert into array form
        pts = np.array(pts)
        wts = np.array(wts)
        # Return the appropriate shaped pair of points and weights
        return (pts, wts)
