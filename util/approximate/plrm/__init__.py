from util.approximate.base import Approximator

class PLRM(Approximator):

    def __init__(self):
        from plrm import plrm
        self.plrm = plrm
        # import os
        # import fmodpy
        # this_dir = os.path.dirname(os.path.abspath(__file__))
        # source_code = os.path.join(this_dir, "plrm.f90")
        # self.plrm = fmodpy.fimport(source_code, blas=True, omp=False, wrap=False,
        #                            verbose=True, output_dir=this_dir).plrm

    def _fit(self, x, y, ds=64, ns=8, steps=2000, seed=None):
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
        self.record = zeros(steps, dtype="float32", order="F")
        self.plrm.minimize_mse(x.T, y.T, steps=steps,
                               num_threads=num_threads,
                               record=self.record)


    def _predict(self, x):
        # Embed the incoming X values.
        from numpy import zeros
        x = ((x - self.x_mean) / self.x_stdev).astype("float32")
        y = zeros((x.shape[0], self.y_mean.shape[0]), dtype="float32")
        # Call the unerlying library.
        self.plrm.evaluate(x.T, y.T)
        # Denormalize the output values and return them.
        y = (y * self.y_stdev) + self.y_mean
        return y



if __name__ == "__main__":
    from util.approximate.testing import test_plot
    m = PLRM()
    p, x, y = test_plot(m, random=True, N=40, plot_points=5000)
    p.show(show=False)
    temp = m._predict(x)
    from util.plot import Plot
    p2 = Plot()
    p2.add("MSE", list(range(len(m.record))), m.record, color=1, mode="lines")
    p2.show(append=True)


# Gfortran compile, open BALS
#   gfortran -O3 -lblas -fopenmp -fPIC -shared -o plrm.so plrm.f90 plrm_c_wrapper.f90"
# 
# Gfortran compile, link Intel BLAS
#   gfortran -O3 -fopenmp -c plrm.f90 plrm_c_wrapper.f90
#   gfortran -fPIC -shared -I/opt/intel/oneapi/intelpython/latest/include -L/opt/intel/oneapi/intelpython/latest/lib -liomp5 -lmkl_rt -o plrm.so plrm.o plrm_c_wrapper.o"
