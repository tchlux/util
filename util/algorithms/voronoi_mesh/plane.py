import random
import numpy as np
from scipy.optimize import minimize

SMALL_NUM = .00001
BIG_NUM = 1000000

L1   = lambda x: np.sum(abs(x)**1) **(1/1)
L2   = lambda x: np.sum(abs(x)**2) **(1/2)
L3   = lambda x: np.sum(abs(x)**3) **(1/3)
L4   = lambda x: np.sum(abs(x)**4) **(1/4)
L5   = lambda x: np.sum(abs(x)**5) **(1/5)
L6   = lambda x: np.sum(abs(x)**6) **(1/6)
LINF = lambda x: max(abs(x))

LINEAR    = lambda r: max(0,1 - r)
QUADRATIC = lambda r: max(0,1 - r)**2
WEIRD     = lambda r: max(0,1 - r)**1.5

BASIS_SIZE = 2.0

class Linn:
    # Initialize the Linn regressor.
    def __init__(self, radius=LINEAR, length=L2, basis_size=BASIS_SIZE):
        self.radius = radius
        self.length = length
        self.basis_size = basis_size

    # Fit a set of points
    def fit(self, points, values):
        self.points = points.copy()
        self.values = values.copy()

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

    # Generate a prediction for a new point
    def predict(self, xs):
        predictions = []
        for x in xs:
            weights = []
            to_skip = set()
            nearest_indices = np.argsort([np.sum((x-p)**2) for p in self.points])
            # Cycle through all training points
            for nearest_ind in range(len(self.points)): #nearest_indices:
                # Skip points that we know will not support x
                if nearest_ind in to_skip: 
                    weights.append(0)
                    continue
                boundary = (float('inf'),)
                # Get the nearest point
                nearest = self.points[nearest_ind]
                vector = (x - nearest)
                # dist_to_x = self.length(vector)
                # Calculate the vectors defining separating planes for
                # all other training points, use the closest one in
                # the direction of 'vector' as the boundary.
                for other_ind in range(len(self.points)):
                    if (nearest_ind == other_ind): continue
                    other = self.points[other_ind]
                    divide = (other - nearest)
                    nearest_to_divide = np.dot(nearest,divide)
                    other_to_divide = np.dot(other,divide)
                    x_to_divide = np.dot(x,divide)
                    # Check to make sure that x is on the same side of
                    # nearest as other is.
                    sign = (other_to_divide - nearest_to_divide) * (x_to_divide - nearest_to_divide)
                    if (sign > 0):
                        # # Identify the midpoint between (nearest + other)
                        # middle = (nearest_to_divide + other_to_divide) / 2
                        # # Identify how much farther x is than the midpoint
                        # ratio = (x_to_divide-nearest_to_divide) / (middle-nearest_to_divide)
                        ratio = 2*(x_to_divide-nearest_to_divide)/(other_to_divide-nearest_to_divide)
                        # Calculate the distance to the boundary
                        # (using the known distance to the point)
                        # dist_to_boundary = ratio
                        dist_to_boundary = self.length(vector)/ratio
                        if (dist_to_boundary < boundary[0]):
                            # ratio = dist_to_x / dist_to_boundary
                            boundary = (dist_to_boundary, ratio)
                    else:
                        # We know that 'other' is blocked from
                        # reaching x by 'nearest' and should not check
                        # for weighting because it must be 0.
                        to_skip.add(other_ind)
                # If there was no boundary between x and nearest
                if (boundary[0] == float('inf')):
                    weights.append( 1 )
                else:
                    new = False
                    if new:
                        # Calculate the weight as an approximation to the
                        # natural neighbor measure (how large is the
                        # difference between the natural boundary and the
                        # new boundary that would be defined by x)
                        weights.append( max(0,boundary[0] - boundary[1]/2) )
                    else:
                        # boundary = (self.length(vector/boundary[0]), boundary[1])

                        # Calculate the weight of this training point response
                        ratio = boundary[1]
                        weights.append( self.radius(ratio / self.basis_size) )
                
            if sum(weights) > 0:
                guess = sum(np.array(weights)*self.values)/sum(weights)
            else:
                guess = 0
            predictions.append( guess )
        # print("weights: ",weights)
        # print("predictions: ",predictions)
        return predictions


if __name__ == "__main__":
    from util.plotly import Plot
    mult = 5
    fun = lambda x: np.cos(x[0]*mult) + np.sin(x[1]*mult)
    np.random.seed(0)
    p = Plot()
    low = 0
    upp = 1
    dim = 2
    plot_points = 2000
    N = 4
    random = True
    if random:
        x = np.random.random(size=(N,dim))
    else:
        N = int(round(N ** (1/dim)))
        x = np.array([r.flatten() for r in np.meshgrid(np.linspace(low,upp,N), np.linspace(low,upp,N))]).T
    y = np.array([fun(v) for v in x])
    p.add("Training Points", *x.T, y)
    
    surf = Linn()
    surf.fit(x,y)
    p.add_func("VMesh", surf, *([(low,upp)]*dim), plot_points=plot_points)
    p.plot(file_name="vmesh.html")

