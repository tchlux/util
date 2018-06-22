# Makes available all custom algorithms that I have worked with
import os
import numpy as np
import fmodpy

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))
# Get the name of a class as listed in python source file.
CLASS_NAME = lambda obj: (repr(obj)[1:-2].split(".")[-1])
# Pretty print an array
CLEAN_ARRAY_STRING = lambda arr: " "+str(arr).replace(",","").replace("[","").replace("]","")
# Check if a variable equals none
IS_NONE = lambda v: type(v) == type(None)
# Useful for checking near-equal cases when roundoff error is involved.
SMALL_NUMBER = 1.4901161193847656*10**(-8) 
#              ^^ SQRT(EPSILON( 0.0_REAL64 ))

# Exceptions for raising useful messages to users.
class ImproperUsage(Exception): pass
class UnexpectedType(Exception): pass
class UnexpectedShape(Exception): pass
class MissingOperator(Exception): pass
class UnderdefinedModel(Exception): pass

# Generic class definition for creating an algorithm that can perform
# regression in the multi-variate setting. In these classes 'x' refers
# to a potentially multi-dimensional set of euclydian coordinates with
# associated single-dimensional response values referred to as 'y'.
# 
# The "fit" function is given both 'x' and associated 'y' and creates
# any necessary (or available) model for predicting 'y' at new 'x'.
# 
# The "predict" function predicts 'y' at potentially novel 'x' values.
# 
# Subclasses are expected to implement _fit and _predict methods which
# take 2D numpy arrays of points as input.
class Approximator:
    _response_dim = 2

    def __str__(self):
        return CLASS_NAME(self)

    # Fit a 2D x and a 2D y with this model.
    def fit(self, x, y, *args, **kwargs):
        if ((type(x) != np.ndarray) or (len(x.shape) != 2)):
            raise(UnexpectedType("Provided 'x' should be a 2D numpy array."))
        if (type(y) != np.ndarray):
            raise(UnexpectedType("Provided 'y' should be a 1D or 2D numpy array."))
        elif (len(y.shape) == 1):
            y = np.reshape(y, (y.shape[0], 1))
            self._response_dim = 1
        elif (len(y.shape) == 2): 
            pass
        else:
            raise(UnexpectedShape("Provided 'y' should be a 1D or 2D numpy array."))
        return self._fit(x, y, *args, **kwargs)
        
    # Predict y at new x locations from fit model.
    def predict(self, x, *args, **kwargs):
        print("Using Approximator predictions..")
        if ((type(x) != np.ndarray) or (0 >= len(x.shape) < 2)):
            raise(UnexpectedType("Provided 'x' should be a 1D or 2D numpy array."))
        single_response = len(x.shape) == 1
        if single_response:
            x = np.reshape(x, (1,len(x)))
        response = np.asarray(self._predict(x, *args, **kwargs), dtype=float)
        # Reduce to one response value if that's what was trained.
        if self._response_dim == 1: response = response[:,0]
        # Reduce to one approximation point if that's what was provided.
        if single_response:         response = response[0]
        # Return the response
        return response

    # Wrapper for 'predict' that returns a single value for a single
    # prediction, or an array of values for an array of predictions
    def __call__(self, *args, **kwargs): return self.predict(*args, **kwargs)


# Generic class definition for an approximator that uses an algorithm
# which predicts new values as a linear combination of provided values.
# These types of approximators are convenient for their ability to
# predict any response value that has defined addition and
# multiplication operators.
# 
# Subclasses should implement "_fit" operator which takes 2D numpy
# arrays of row-points, and the "_predict" operator which takes 2D
# numpy arrays of row-points and returns the indices of source points
# and associated weights required to make a prediction.
class WeightedApproximator(Approximator):
    WeightedApproximator = True
    y = None

    # Fit a 2D x and a 2D y with this model.
    def fit(self, x, y=None, **kwargs):
        if ((type(x) != np.ndarray) or (len(x.shape) != 2)):
            raise(UnexpectedType("Provided 'x' should be a 2D numpy array."))
        if (not IS_NONE(y)):
            if (not hasattr(y, "__len__")):
                raise(MissingOperator("Provided 'y' to fit must have defined '__len__' operator."))
            elif (not hasattr(y, "__getitem__")):
                raise(MissingOperator("Provided 'y' to fit must have defined '__getitem__' operator."))
            elif (not hasattr(y[0], "__add__")):
                raise(MissingOperator("Elements of provided 'y' must have defined '__add__' operator."))
            elif (not hasattr(y[0], "__mul__")):
                raise(MissingOperator("Elements of provided 'y' must have defined '__mul__' operator."))
            self.y = y
        # Fit the provided x values.
        return self._fit(x, **kwargs)
        
    # Predict y at new x locations from fit model.
    def predict(self, x, *args, **kwargs):
        print("Using Weighted Approximator predictions..")
        if ((type(x) != np.ndarray) or (0 >= len(x.shape) < 2)):
            raise(UnexpectedType("Provided 'x' should be a 1D or 2D numpy array."))
        single_response = len(x.shape) == 1
        if single_response:
            x = np.reshape(x, (1,len(x)))
        indices, weights = self._predict(x, *args, **kwargs)
        # Return the indices and weights if no y values were provided.
        if (IS_NONE(self.y)): 
            response = [(ids,wts) for (ids,wts) in zip(indices, weights)]
        else:
            # Collect response values via weighted sums of self.y values
            response = []        
            for ids, wts in zip(indices, weights):
                response.append( sum(self.y[i]*w for (i,w) in zip(ids, wts)) )
        # Reduce to one approximation point if that's what was provided.
        if single_response: response = response[0]
        # Return the response
        return response

    # Wrapper for 'predict' that returns a single value for a single
    # prediction, or an array of values for an array of predictions
    def __call__(self, *args, **kwargs): return self.predict(*args, **kwargs)    


# A wrapper for an approximator that must fit unique points.
def uniqueify(weighted_approximator):
    # A class for wrapping an approximator.
    class UniquePoints(weighted_approximator):
        original_points = None
        unique_points = None
        unique_indices = None
        original_values = None

        # Wrap the fit method to capture only unique points.
        def fit(self, points, values=None, *args, **kwargs):
            if ((type(points) != np.ndarray) or (len(points.shape) != 2)):
                raise(UnexpectedType("Expected 2D numpy array as first argument."))
            self.original_points = points.copy()
            self.unique_points = {}
            for i,pt in enumerate(points):
                pt = tuple(pt)
                self.unique_points[pt] = self.unique_points.get(pt, []) + [i]
            # Store the indices of the first occurrence of each unique point.
            self.unique_indices = np.array(sorted(
                self.unique_points[pt][0] for pt in self.unique_points))
            # Average the response value for the points that are identical.
            if (not IS_NONE(values)):
                self.original_values = values
                to_add = set(self.unique_points)
                avg_values = []
                for pt in self.original_points:
                    pt = tuple(pt)
                    if pt in to_add:
                        indices = self.unique_points[pt]
                        wt = 1. / len(indices)
                        avg_values.append( sum(values[i]*wt for i in indices) )
                        to_add.remove(pt)
                args = args + (avg_values,)
            # Call the fit method on parent with unique points only.
            return super().fit(self.original_points[self.unique_indices,:], 
                               *args, **kwargs)
            
        # Wrap the predict method to return original points.
        def predict(self, points, *args, **kwargs):
            if ((type(points) != np.ndarray) or (0 >= len(points.shape) < 2)):
                raise(UnexpectedType("Provided 'points' should be a 1D or 2D numpy array."))
            # If values were provided, return usual prediction.
            elif not IS_NONE(self.original_values):
                return super().predict(points)
            # Otherwise we are getting the points and indices in original data.
            single_response = len(points.shape) == 1
            if single_response:
                points = np.reshape(points, (1,len(points)))
            # Otherwise, return points and weights in original indices.
            indices, weights = self._predict(points, *args, **kwargs)
            response = []
            for ids, wts in zip(indices, weights):
                orig_ids = []
                orig_wts = []
                for (i,w) in zip(ids, wts):
                    pt = tuple(self.original_points[self.unique_indices[i]])
                    orig_ids += self.unique_points[pt]
                    # w /= len(self.unique_points[pt]) # <- Equally weight unique points.
                    orig_wts += [w] * len(self.unique_points[pt])
                # Normalize sum of weights, giving repeated points higher 'weight'
                orig_wts_sum = sum(orig_wts)
                orig_wts = [w / orig_wts_sum for w in orig_wts]
                response.append( (orig_ids, orig_wts) )
            if single_response: response = response[0]
            return response
        
    # Return the newly constructed class
    return UniquePoints


# Given a model, test the time it takes to make a prediction
def test_time(model, ns=range(10,10011,1000), ds=range(2,11,2)):
    import time
    print("Starting")
    for d in ds:
        print(f"D: {d}")
        print()
        for n in ns:
            print(f" N: {n}")
            x = np.random.random((n,d))
            y = np.random.random((n,))
            start = time.time()
            model.fit(x, y)
            print(f"  fit:     {time.time() - start:5.2f}")
            start = time.time()
            model(x[0])
            print(f"  predict: {time.time() - start:5.2f}")
            print()
        print()

# Given a model, generate a test plot demonstrating the surface it procudes
def test_plot(model, low=0, upp=1, plot_points=3000, p=None,
              fun=lambda x: 3*x[0]+.5*np.cos(8*x[0])+np.sin(5*x[-1]),
              N=20, D=2, random=True, seed=0):
    np.random.seed(seed)
    # Generate x points
    if random:
        x = np.random.random(size=(N,D))
    else:
        N = int(round(N ** (1/D)))
        x = np.array([r.flatten() for r in np.meshgrid(*[np.linspace(0,1,N)]*D)]).T
    # Calculate response values
    y = np.array([fun(v) for v in x])
    # Fit the model to the points
    model.fit(x,y)
    # Generate the plot
    from util.plot import Plot
    if type(p) == type(None): p = Plot()
    p.add("Training Points", *x.T, y)
    p.add_func(str(model), model, *([(low-.1,upp+.1)]*D),
               plot_points=plot_points, vectorized=True)
    return p, x, y

# Given a model, use the method "points_and_weights" to identify the
# regions of support for each of the interpolation points.
def test_support(model, low=0, upp=1, plot_points=3000, p=None,
                 fun=lambda x: 3*x[0]+.5*np.cos(8*x[0])+np.sin(5*x[-1]),
                 N=20, D=2, random=True, seed=0):
    # Force D to be 2
    D = 2
    np.random.seed(seed)
    # Generate x points
    if random:
        x = np.random.random(size=(N,D))
    else:
        N = int(round(N ** (1/D)))
        x = np.array([r.flatten() for r in np.meshgrid(*[np.linspace(0,1,N)]*D)]).T
    # Calculate response values
    y = np.array([fun(v) for v in x])
    # Fit the model to the points
    model.fit(x,y)
    # Generate the plot
    from util.plotly import Plot
    if type(p) == type(None): p = Plot()
    p.add("Training Points", *x.T, color=p.color(len(x)))
    for i in range(len(x)):
        name = f"{i+1}" 
        p.add(name, [x[i][0]], [x[i][1]], group=name)
        def supported(pt):
            pts, wts = model.points_and_weights(pt)
            return (i in pts)
        p.add_region(name+" region", supported, *([(low-.1,upp+.1)]*D),
                     color=p.color(p.color_num), group=name,
                     plot_points=plot_points, show_in_legend=False) 
    # p.add_func(str(model), model, *([(low-.1,upp+.1)]*D),
    #            plot_points=plot_points, vectorized=True)
    return p, x, y


# Import all defined algorithms
from util.algorithms.voronoi import Voronoi
from util.algorithms.nearest_neighbor import NearestNeighbor
from util.algorithms.delaunay import Delaunay, qHullDelaunay

if __name__ == "__main__":
    print("Adding surface to plot..")
    p,_,_ = test_plot(Voronoi(), N=20, D=2, low=-.1, upp=1.1,
                      random=True, plot_points=4000) # 6, 8
    print("Generating plot HTML..")
    p.show()

    # # ================================================
    # #      Testing Regions of Support for Voronoi     
    # # ================================================

    # p,_,_ = test_support(Voronoi(), N=3, D=2, low=-.4, upp=1.4,
    #                      random=True, seed=8) # 6
    # p.show(width=700, height=700)
    # p,_,_ = test_support(VoronoiMesh(), N=3, D=2, low=-.4, upp=1.4,
    #                      random=True, seed=8) # 6
    # p.show(width=700, height=700)

    # p, x, y = test_plot(Delaunay, N=9, low=-.2, upp=1.2, random=False)
    # p.plot("VT Delaunay", file_name="bad_extrap.html", show=False)
    # p, x, y = test_plot(qHullDelaunay, N=9, low=-.2, upp=1.2, random=False)
    # p.plot("QHull Delaunay", file_name="bad_extrap.html", show=False, append=True)

    # # ==================================
    # #      TOMS_Delaunay Error Case     
    # # ==================================

    # from util.plotly import Plot
    # surf = Delaunay()
    # surf.fit(x,y)
    # surf2 = qHullDelaunay()
    # surf2.fit(x,y)

    # # bad = np.array([.2, -.1])
    # bad = np.array([-.01, .01])
    # sp, sw = surf.points_and_weights(bad)
    # s2p, s2w = surf2.points_and_weights(bad)

    # np.savetxt("bad_train.csv", x)
    # print("extrap_point = (\ %s \)"%(", ".join(map(str,bad))))
    # print("OBTAINED simplex:", sp + 1)
    # print("         weights:", sw)
    # print("EXPECTED simplex:", s2p + 1)
    # print("         weights:", s2w)

    # p = Plot()
    # p.add("Given Points", x[:,0], x[:,1], shade=False)
    # p.add("Extrap Point", [bad[0]], [bad[1]], color="rgb(0,0,0)", shade=False)
    # p.add("Obtained Simplex", (surf.pts.T)[sp][:,0],
    #       (surf.pts.T)[sp][:,1], color=p.color(3), shade=False)
    # p.add("Expected Simplex", surf2.pts[s2p][:,0],
    #       surf2.pts[s2p][:,1], color=p.color(1), shade=False)
    # p.show(file_name="bad_extrap.html", append=True)

