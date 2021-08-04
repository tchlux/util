


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


class EmbeddingInterpolator(Approximator):
    def __init__(self, model=_default_interpolant):
        print("Loading PLRM for embedding interplator..")
        import fmodpy
        self.plrm = fmodpy.fimport("plrm.f90", blas=True, omp=False, verbose=True).plrm
        self.model = model

    def _fit(self, x, y, ds=16, ns=8, seed=None, steps=1000):
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
        self.record = zeros(steps, dtype="float32", order="F")
        self.plrm.minimize_mse(x.T, y.T, steps=steps,
                               num_threads=num_threads,
                               record=self.record)
        # Embed to get the internal x.
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
    p, x, y = test_plot(m, random=True, N=40, plot_points=5000)
    p.show(show=False)
    from util.plot import Plot
    p2 = Plot()
    p2.add("MSE", list(range(len(m.record))), m.record, color=1, mode="lines")
    p2.show(append=True)

