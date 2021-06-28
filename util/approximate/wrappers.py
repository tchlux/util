import os
import numpy as np
from util.math import abs_diff
from util.approximate.base import Approximator

MAX_DIM = None
MAX_SAMPLES = None
DEFAULT_COND_SCALE = True
DEFAULT_CONDITIONER = "PLR" # "PCA" # MPCA
PATH_TO_PLRM = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "embedding_interpolator", "plrm.f90")

class Split(Approximator):
    def __init__(self, model_type, k=5, seed=0, *args, **kwargs):
        # Initialize "k" internal models.
        self._seed = seed
        self._count = k
        self._models = []
        for i in range(k): self._models.append( model_type(*args, **kwargs) )

    # Fit multiple internal models.
    def _fit(self, x, y):
        np.random.seed(self._seed)
        splits = np.arange(len(x))
        np.random.shuffle(splits)
        for i in range(self._count):
            indices = splits[i:len(x):self._count]
            self._models[i].fit( x[indices], y[indices] )
        
    # Average over the predictions of multiple models.
    def _predict(self, z):
        outcomes = []
        for m in self._models:
            outcomes.append( m(z) )
        return np.mean(np.array(outcomes), axis=0)


# A wrapper for an approximator that must fit unique points. Use this
# wrapper on an un-initialized "weighted approximator" to inheret all
# methods while modifying the "fit" and "predict" methods.
def unique(weighted_approximator):
    # ----------------------------------------------------------------
    # A class for wrapping an approximator.
    class UniquePointsApproximator(weighted_approximator):
        original_points = None
        unique_points = None
        unique_indices = None
        original_values = None

        # Wrap the fit method to capture only unique points.
        def fit(self, points, values=None, *args, **kwargs):
            if ((not issubclass(type(points), np.ndarray)) or (len(points.shape) != 2)):
                from util.approximate.base import UnexpectedType
                raise(UnexpectedType("Expected 2D numpy array as first argument."))
            self.original_points = points
            self.unique_points = {}
            for i,pt in enumerate(points):
                pt = tuple(pt)
                self.unique_points[pt] = self.unique_points.get(pt, []) + [i]
            # Store the indices of the first occurrence of each unique point.
            self.unique_indices = np.array(sorted(
                self.unique_points[pt][0] for pt in self.unique_points))
            # Average the response value for the points that are identical.
            if (values is not None):
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
            if ((not issubclass(type(points), np.ndarray)) or (0 >= len(points.shape) < 2)):
                from util.approximate.base import UnexpectedType
                raise(UnexpectedType("Provided 'points' should be a 1D or 2D numpy array."))
            # If values were provided, return usual prediction.
            elif (self.original_values is not None):
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
    # ----------------------------------------------------------------        
    # Return the newly constructed class
    return UniquePointsApproximator


# A wrapper for an approximator that must fit unique points. Use this
# wrapper on an un-initialized "weighted approximator" to inheret all
# methods while modifying the "fit" and "predict" methods.
def condition(approximator, metric=abs_diff, method=DEFAULT_CONDITIONER,
              dim=MAX_DIM, samples=MAX_SAMPLES, scale=DEFAULT_COND_SCALE, 
              display=False, seed=None, **cond_kwargs):
    if method == "PLR":
        import fmodpy
        plrm = fmodpy.fimport(PATH_TO_PLRM, blas=True, omp=True, verbose=False,
                              output_dir=os.path.dirname(PATH_TO_PLRM)).plrm
        if (seed is not None):
            num_threads = 2
    elif method == "PCA":
        from util.stats import pca, normalize_error
    elif method == "MPCA":
        from util.stats import mpca
    # ----------------------------------------------------------------
    class ConditionedApproximator(approximator):
        # Wrapped "fit" method, this incorporates dimension reduction
        # and approximation conditioning based on the selected method.
        def _fit(self, x, y, *args, num_comps=None, layers=8, steps=1000, **kwargs):
            # Set the number of components appropriately.
            if (num_comps is None): num_comps = min(x.shape + (16,))
            if (dim is not None):   num_comps = min(dim, num_comps)
            # Compute the components and the values.
            if method == "PLR":
                original_y = y
                if (len(y.shape) == 1): y = y[:,None]
                self._num_comps = num_comps
                # Use a trained piecewise linear regressor to condition.
                self._x_mean = x.mean(axis=0)
                self._x_stdev = x.std(axis=0)
                # Normalize the X and Y data for this model.
                x = np.asarray((x - self._x_mean) / self._x_stdev, dtype="float32", order="C")
                y_mean = y.mean(axis=0)
                y_stdev = y.std(axis=0)
                y = np.asarray((y - y_mean) / y_stdev, dtype="float32")
                # Fit neural network embedder.
                di = x.shape[1]
                do = y.shape[1]
                ds = num_comps
                ns = layers
                plrm.new_model(di, ds, ns, do)
                plrm.init_model(seed=seed)
                plrm.minimize_mse(x.T, y.T, steps=steps, num_threads=num_threads)
                # Embed to get the internal x.
                x, old_x = np.zeros((x.shape[0], ds), dtype="float32"), x
                plrm.evaluate(old_x.T, x.T)
                # Convert back into 64 bit precision.
                x = np.asarray(x, dtype=float)
                return super()._fit(x, original_y, *args, **kwargs)
            elif method == "PCA":
                # Compute the principle components as the new axes.
                components, values = pca(x, num_components=num_comps, display=display, **cond_kwargs)
                # Compute the values so that the transformed points have unit metric slope.
                values = normalize_error(np.matmul(x, components.T), y, metric, display=display)
            elif method == "MPCA":
                # Use metric PCA to compute components and values.
                components, values = mpca(x, y, metric=metric,
                                          num_components=num_comps,
                                          num_vecs=samples,
                                          display=display,
                                          **cond_kwargs)
            # Reset the values if scale should not be used.
            if not scale: values[:] = 1.
            if display:
                np.set_printoptions(precision=3, sign=" ")
                print("\nComponents and values:")
                for (c,v) in zip(components, values): print(f" {v:.2f}  {c}")
                print()
                np.set_printoptions(precision=8, sign="-")
            # Generate the conditioning matrix.
            self.conditioner = np.matmul(np.diag(values), components).T
            # Return the normal fit operation.
            return super()._fit(np.matmul(x, self.conditioner), y, *args, **kwargs)

        # Return the predictions made on the similarly conditioned points.
        def _predict(self, x, *args, **kwargs):
            if (method == "PLR"):
                # Embed the incoming X values.
                x = np.asarray((x - self._x_mean) / self._x_stdev, dtype="float32", order="C")
                embedded_x = np.zeros((x.shape[0], self._num_comps), dtype="float32")
                plrm.evaluate(x.T, embedded_x.T)
                # Return the model approximation over embedded X.
                return super()._predict(np.asarray(embedded_x, dtype=float), *args, **kwargs)
            else:
                return super()._predict(np.matmul(x, self.conditioner), *args, **kwargs)
    # ----------------------------------------------------------------
    return ConditionedApproximator

