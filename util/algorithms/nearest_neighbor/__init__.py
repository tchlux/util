import numpy as np
from util.algorithms import WeightedApproximator
from util.math import is_none, abs_diff

# Automatically select the value for "k" for using nearest neighbor
# based on a random sample of the points provided.
def auto(points, values, metric=abs_diff, max_k=None, samples=100,
         mean=False, k_step=1, display=False):
    from util.stats import random_range
    # Make the maximum value for "k" the nearest power of 2 that
    # contains less than or equal to half of the provided data.
    if is_none(max_k):
        from math import floor, log2
        max_k = 2**floor(log2(len(points)//2))
    # Compute up to 'max_k' nearest neighbors about selected points.
    # Add "+1" to exclude the current active point as a neighbor.
    model = NearestNeighbor(k=max_k+1)
    model.fit(points)
    # Randomly pick a set of points as the "checks".
    indices = [i for i in random_range(len(points), count=samples)]
    neighbors = np.array([i[1:] for (i,w) in model(points[indices])])
    differences = np.array([[metric(values[i1], values[i2]) for i2 in neighbors[i]]
                            for i,i1 in enumerate(indices)])
    k_values = {}
    # Pick the function for identifying the best selection of "k".
    for k_pow in range(0, int(log2(max_k))+1, k_step):
        k = 2**k_pow
        if mean: k_values[k] = np.mean(differences[:,:k])
        else:    k_values[k] = np.mean(np.min(differences[:,:k], axis=1))
    if (2**k_pow != max_k):
        k = 2**int(log2(max_k))
        if mean: k_values[k] = np.mean(differences[:,:k])
        else:    k_values[k] = np.mean(np.min(differences[:,:k], axis=1))
    # Find the k with the lowest mean error.
    best_k = min(k_values.items(), key=lambda i: i[1])[0]
    if display:
        name = "mean" if mean else "minimum"
        from math import log10, ceil
        print('-'*52)
        print(" Estimated "+name+" error for various choices of 'k':")
        for k in sorted(k_values):
            extra = "  <-- chosen 'k'" if k == best_k else ""
            print(f"  k = {k:{ceil(log10(max_k))}d} ~ {k_values[k]:.4e}"+extra)
        print('-'*52)
    # Return the "k" with the minimum mean difference 
    return best_k

# Construct an approximation algorithm that only returns the average
# of the fit points.
class Average(WeightedApproximator):
    def _fit(self, points):
        self.points = points
    def _predict(self, points):
        indices = [list(range(len(self.points)))] * len(points)
        w = 1 / len(self.points)
        weights = [[w]*len(self.points)] * len(points)
        return np.array(indices), np.array(weights)
# Import another approximation algorithm that reasonably blends points.
from util.algorithms import Voronoi

# ===========================================
#      Simple Nearest Neighbor Regressor     
# ===========================================

# Class for computing an interpolation between the nearest n neighbors
class NearestNeighbor(WeightedApproximator):
    def __init__(self, k=None, method=Average, **auto_kwargs):
        self.points = None
        self.method = method
        self.num_neighbors = k
        self.auto_kwargs = auto_kwargs

    # Use fortran code to compute the boxes for the given data
    def _fit(self, control_points, k=None, **kwargs):
        if not is_none(k): self.num_neighbors = k
        # Process and store local information
        self.points = control_points.copy()
        # Automatically select the value for "k" if appropriate and
        # the response values are available for the points.
        if is_none(self.num_neighbors):
            if (not is_none(self.y)):
                self.auto_kwargs.update(kwargs)
                # If "mean" was not provided, pick based on problem type.
                if "mean" not in self.auto_kwargs:
                    self.auto_kwargs["mean"] = not self.classifier
                # Begin the estimation of best value for 'k'.
                print("Nearest neighbor, estimating best value for 'k'..", end="\r", flush=True)
                self.num_neighbors = auto(self.points, self.y, **self.auto_kwargs)
                print("                                                 ", end="\r", flush=True)
            else:
                self.num_neighbors = 1

    # Function that returns the indices of points and the weights that
    # should be used to make associated predictions for each point in
    # "points".
    def _predict(self, points, display=True):
        if display: import time
        # Body
        pts = []
        wts = []
        if display:
            update = end = ""
            start = time.time()
        for i,pt in enumerate(points):
            if display and (time.time() - start > 1):
                update = "\b"*len(end)
                end = f" nearest neighbor predicting {i+1} of {len(points)}.."
                print(update, end=end, flush=True)
                start = time.time()
            distances = np.sum((self.points - pt)**2, axis=1)
            closest = np.argsort(distances)[:self.num_neighbors]
            distances = distances[closest]
            # Use another weighted approximator to handle final prediction.
            model = self.method()
            model.fit(self.points[closest])
            guess_pts, guess_wts = model.predict(pt)
            # Return the indices and weights associated with each.
            pts.append( closest[guess_pts] )
            wts.append( guess_wts )
        if display: print("\b"*len(end), end="", flush=True)
        # Convert into array form
        pts = np.array(pts)
        wts = np.array(wts)
        # Return the appropriate shaped pair of points and weights
        return (pts, wts)


if __name__ == "__main__":
    d = 10
    n = 511
    points = np.random.random(size=(n,d))
    values = np.random.random(size=(n,d))
    k = auto(points, values, display=True)

