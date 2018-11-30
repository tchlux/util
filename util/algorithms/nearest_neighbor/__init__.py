import numpy as np
from util.algorithms import WeightedApproximator

# Construct an approximation algorithm that only returns the average
# of the fit points.
class Average(WeightedApproximator):
    def _fit(self, points):
        self.points = points.copy()
    def _predict(self, points):
        indices = [list(range(len(self.points)))] * len(points)
        w = 1 / len(self.points)
        weights = [[w]*len(self.points)] * len(points)
        return np.array(indices), np.array(weights)
# Import another approximation algorithm that reasonably blends points.
from util.algorithms import Voronoi

# Nearest neighbor
NN_DEFAULT_NUM_NEIGHBORS = 1
    
# ===========================================
#      Simple Nearest Neighbor Regressor     
# ===========================================

# Class for computing an interpolation between the nearest n neighbors
class NearestNeighbor(WeightedApproximator):
    def __init__(self, k=NN_DEFAULT_NUM_NEIGHBORS, method=Average):
        self.points = None
        self.method = method
        self.num_neighbors = k

    # Use fortran code to compute the boxes for the given data
    def _fit(self, control_points, k=None):
        if type(k) != type(None):
            self.num_neighbors = k
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
                # Use another weighted approximator to handle final prediction.
                model = self.method()
                model.fit(self.points[closest])
                points, weights = model.predict(pt)
                # Return the indices and weights associated with each.
                pts.append( closest[points] )
                wts.append( weights )
        # Convert into array form
        pts = np.array(pts)
        wts = np.array(wts)
        # Return the appropriate shaped pair of points and weights
        return (pts, wts)
