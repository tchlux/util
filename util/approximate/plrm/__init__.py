from util.approximate.base import Approximator

# SOURCE_FILE = "plrm.f90"
# SOURCE_FILE = "stable_flexure.f90"
SOURCE_FILE = "stable_relu.f90"

class PLRM(Approximator):

    def __init__(self):
        # from plrm import plrm
        # self.plrm = plrm
        import os
        import fmodpy
        this_dir = os.path.dirname(os.path.abspath(__file__))
        source_code = os.path.join(this_dir, SOURCE_FILE)
        plrm = fmodpy.fimport(source_code, blas=True, omp=True, wrap=True,
                              verbose=False, output_dir=this_dir)
        self.plrm = plrm.plrm

    def _fit(self, x, y, ds=32, ns=8, steps=1000, seed=None):
        # If a random seed is provided, then only 2 threads can be used
        #  because nondeterministic behavior is exhibited otherwise.
        if (seed is not None): num_threads = 2
        else:                  num_threads = None
        # Core numpy utilities.
        from numpy import where, zeros
        # Store the X and Y data for this model.
        self.x_mean = x.mean(axis=0)
        self.x_stdev = x.std(axis=0)
        
        # Normalize the X and Y data for this model.
        x = ((x - self.x_mean) / self.x_stdev).astype("float32")
        self.y_mean = y.mean(axis=0)
        self.y_stdev = y.std(axis=0)
        self.y_stdev = where(self.y_stdev == 0, 1.0, self.y_stdev)
        y = ((y - self.y_mean) / self.y_stdev).astype("float32")
        # Fit internal piecewise linear regression model.
        di = x.shape[1]
        do = y.shape[1]
        self.plrm.new_model(di, ds, ns, do)
        self.plrm.init_model(seed=seed)
        mse = self.plrm.minimize_mse(x.T, y.T, steps=steps,
                                     num_threads=num_threads,
                                     keep_best=True)


    def _predict(self, x, embed=False):
        from numpy import zeros, asarray
        # Evaluate the model at all data.
        x = asarray((x - self.x_mean) / self.x_stdev, dtype="float32", order='C')
        if embed:
            y = zeros((x.shape[0], self.plrm.mds), dtype="float32", order='C')
        else:
            y = zeros((x.shape[0], self.y_mean.shape[0]), dtype="float32", order='C')
        # Call the unerlying library.
        self.plrm.evaluate(x.T, y.T)
        # Denormalize the output values and return them.
        if (not embed):
            y = (y * self.y_stdev) + self.y_mean
        return y


if __name__ == "__main__":
    import numpy as np
    from util.approximate.testing import test_plot
    m = PLRM()
    p, x, y = test_plot(m, random=True, N=40, plot_points=500) #plot_points=5000)
    print("Generating surface plot..")
    p.show(show=True)
    print("", "done.", flush=True)
