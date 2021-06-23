


from util.approximate import NearestNeighbor
_default_interpolant = NearestNeighbor()

# from util.approximate import Delaunay
# _default_interpolant = Delaunay()

# from util.approximate import Voronoi
# _default_interpolant = Voronoi()

# from util.approximate import LSHEP
# _default_interpolant = LSHEP()

# from util.approximate import ShepMod
# _default_interpolant = ShepMod()

from util.approximate.base import Approximator

class PLRM(Approximator):
    def __init__(self):
        import os
        import fmodpy
        this_dir = os.path.dirname(os.path.abspath(__file__))
        source_code = os.path.join(this_dir, "plrm.f90")
        output_dir = os.path.join(this_dir, "plrm")
        self.plrm = fmodpy.fimport(source_code, blas=True, omp=False,
                                   verbose=False, output_dir=output_dir).plrm

    def _fit(self, x, y, ds=64, ns=8, steps=1000, seed=None):
        # If a random seed is provided, then only 2 threads can be used
        #  because nondeterministic behavior is exhibited otherwise.
        if (seed is not None): num_threads = 2
        else:                  num_threads = None
        # Core numpy utilities.
        from numpy import asarray, zeros, where
        # Store the X and Y data for this model.
        self.x_mean = x.mean(axis=0)
        self.x_stdev = x.std(axis=0)
        self.y = y.copy()
        # Normalize the X and Y data for this model.
        x = asarray((x - self.x_mean) / self.x_stdev, dtype="float32")
        self.y_mean = y.mean(axis=0)
        self.y_stdev = y.std(axis=0)
        self.y_stdev = where(self.y_stdev == 0, 1.0, self.y_stdev)
        y = asarray((y - self.y_mean) / self.y_stdev, dtype="float32")
        # Fit internal neural network embedder.
        di = x.shape[1]
        do = y.shape[1]
        self.plrm.new_model(di, ds, ns, do)
        self.plrm.init_model(seed=seed)
        self.plrm.minimize_mse(x.T, y.T, steps=steps, num_threads=num_threads)


    def _predict(self, x):
        # Embed the incoming X values.
        from numpy import asarray, zeros
        x = asarray((x - self.x_mean) / self.x_stdev, dtype="float32")
        y = zeros((x.shape[0], self.y_mean.shape[0]), dtype="float32")
        self.plrm.evaluate(x.T, y.T)
        # Denormalize the output values and return them.
        y = (y * self.y_stdev) + self.y_mean
        return y




class EmbeddingInterpolator(Approximator):
    def __init__(self, model=_default_interpolant):
        import fmodpy
        self.plrm = fmodpy.fimport("plrm.f90", blas=True, omp=False, verbose=False).plrm
        self.model = model

    def _fit(self, x, y, ds=8, ns=4, seed=None, steps=500):
        # If a random seed is provided, then only 2 threads can be used
        #  because nondeterministic behavior is exhibited otherwise.
        if (seed is not None): num_threads = 2
        else:                  num_threads = None
        # Core numpy utilities.
        from numpy import asarray, zeros
        # Store the X and Y data for this model.
        self.x_mean = x.mean(axis=0)
        self.x_stdev = x.std(axis=0)
        self.y = y.copy()
        # Normalize the X and Y data for this model.
        x = asarray((x - self.x_mean) / self.x_stdev, dtype="float32")
        y_mean = y.mean(axis=0)
        y_stdev = y.std(axis=0)
        y = asarray((y - y_mean) / y_stdev, dtype="float32")
        # Fit internal neural network embedder.
        di = x.shape[1]
        do = y.shape[1]
        self.plrm.new_model(di, ds, ns, do)
        self.plrm.init_model(seed=seed)
        self.plrm.minimize_mse(x.T, y.T, steps=steps, num_threads=num_threads)
        # Embed to get the internal x.
        from numpy import zeros
        self.x = zeros((x.shape[0], ds), dtype="float32")
        self.plrm.evaluate(x.T, self.x.T)
        # Fit the internal model to the embedded data.
        self.model.fit(asarray(self.x, dtype=float), asarray(self.y, dtype=float))


    def _predict(self, x):
        # Embed the incoming X values.
        from numpy import asarray, zeros
        x = asarray((x - self.x_mean) / self.x_stdev, dtype="float32")
        embedded_x = zeros((x.shape[0], self.x.shape[1]), dtype="float32")
        self.plrm.evaluate(x.T, embedded_x.T)
        # Return the internal model approximation over embedded X.
        return self.model.predict(asarray(embedded_x, dtype=float))


if __name__ == "__main__":
    from util.approximate.testing import test_plot
    m = EmbeddingInterpolator()
    # m = PLRM()
    p, x, y = test_plot(m, random=True, N=40, plot_points=5000)
    p.show()

