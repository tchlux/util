


from util.approximate import NearestNeighbor
_default_interpolant = NearestNeighbor()

# from util.approximate import Delaunay
# _default_interpolant = Delaunay()

# from util.approximate import Voronoi
# _default_interpolant = Voronoi()

from util.approximate.base import Approximator
class EmbeddingInterpolator(Approximator):
    def __init__(self, model=_default_interpolant):
        import fmodpy
        self.plrm = fmodpy.fimport("plrm.f90", blas=True, omp=True, verbose=False).plrm
        self.model = model

    def _fit(self, x, y, ds=8, ns=4, seed=0, steps=500):
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
        self.plrm.init_model(inputs=x.T, outputs=y.T, seed=seed)
        self.plrm.minimize_mse(x.T, y.T, steps=steps)
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
    p, x, y = test_plot(m, random=True, N=200)
    p.show()

