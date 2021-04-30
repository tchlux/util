import numpy as np
from util.approximate import WeightedApproximator
from util.math import abs_diff

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

# ===========================================
#      Simple Nearest Neighbor Regressor     
# ===========================================

# Class for computing an interpolation between the nearest n neighbors
class NearestNeighbor(WeightedApproximator):
    def __init__(self, k=1, method=Average, **auto_kwargs):
        self.points = None
        self.method = method
        self.num_neighbors = k
        self.auto_kwargs = auto_kwargs

    # Use fortran code to compute the boxes for the given data
    def _fit(self, control_points, k=None, display=True, **kwargs):
        # from balltree import BallTree
        from sklearn.neighbors import BallTree
        if (k is not None): self.num_neighbors=k
        # Process and store local information
        self.points = control_points.copy()
        self.tree = BallTree(self.points)
        # Update the associated 'y' values if they exist (since points
        # were shuffled on the construction of the tree).
        if (self.y is not None) and hasattr(self.tree, "index_mapping"):
            self.y = [self.y[i] for i in self.tree.index_mapping]
        # Automatically select the value for "k" if appropriate and
        # the response values are available for the points.
        if (self.num_neighbors is None):
            if (self.y is not None):
                self.auto_kwargs.update(kwargs)
                # If "mean" was not provided, pick based on problem type.
                if "mean" not in self.auto_kwargs:
                    self.auto_kwargs["mean"] = not self.classifier
                # Begin the estimation of best value for 'k'.
                end = ("\n" if display else "\r")
                print("Nearest neighbor, estimating best value for 'k'..", end=end, flush=True)
                self.num_neighbors = auto(self.points, self.y, **self.auto_kwargs)
                print(f"  chose k = {self.num_neighbors}", end=end, flush=True)
                print("                                                 ", end=end, flush=True)
            else:
                self.num_neighbors = 1

    # Function that returns the indices of points and the weights that
    # should be used to make associated predictions for each point in
    # "points".
    def _predict(self, points, max_size=2**20, display=True):
        if display: import time
        # Body
        pts = []
        wts = []
        if display:
            update = end = ""
            start = time.time()
        # Cycle through points and make prediction.
        for i,pt in enumerate(points):
            if display and (time.time() - start > 1):
                update = "\b"*len(end)
                end = f" nearest neighbor predicting {i+1} of {len(points)}.."
                print(update, end=end, flush=True)
                start = time.time()
            closest = self.tree.query(pt[None,:], k=self.num_neighbors,
                                      return_distance=False)[0]
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

# Automatically select the value for "k" for using nearest neighbor
# based on a random sample of the points provided.
def auto(points, values, metric=abs_diff, max_k=None, samples=100,
         mean=True, k_step=1, model=NearestNeighbor, display=False):
    from util.random import random_range
    # Make the maximum value for "k" the nearest power of 2 that
    # contains less than or equal to half of the provided data.
    if (max_k is None):
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


if __name__ == "__main__":

    # np.random.seed(0)

    # pts = np.random.random(size=(10000,5))
    # values = np.random.random(size=(10000,))
    # test = np.random.random(size=(1000,5))
    # print("Fitting model..", flush=True)
    # m = NearestNeighbor()
    # m.fit(pts, values)
    # print("Evaluating model..", flush=True)
    # m(test)
    # print("done")
    # exit()

    TEST_AUTO = True
    if TEST_AUTO:
        d = 10
        n = 511
        points = np.random.random(size=(n,d))
        values = np.random.random(size=(n,d))
        k = auto(points, values, display=True)

    # from util.approximate import condition
    # m = condition(NearestNeighbor, method="PLR")()
    # m = condition(NearestNeighbor, method="MPCA", dim=1)()
    # m = condition(NearestNeighbor, method="MPCA", dim=2)()
    # m = condition(NearestNeighbor, method="PCA", dim=2)()
    # m = NearestNeighbor()
    from util.approximate.testing import test_plot
    # p,x,y = test_plot(m, N=30, fun=lambda x: x[0]**3, plot_points=4000)
    p,x,y = test_plot(m, N=300)
    p.show()
