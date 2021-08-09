from util.approximate.base import Approximator

class PLRM(Approximator):

    def __init__(self):
        # from plrm import plrm
        # self.plrm = plrm
        import os
        import fmodpy
        this_dir = os.path.dirname(os.path.abspath(__file__))
        source_code = os.path.join(this_dir, "plrm.f90")
        plrm = fmodpy.fimport(source_code, blas=True, omp=True, wrap=True,
                                   verbose=True, output_dir=this_dir)
        self.plrm = plrm.plrm
        self.vis = plrm.vis_plrm

    def _fit(self, x, y, ds=16, ns=8, steps=2000, seed=None):
        # If a random seed is provided, then only 2 threads can be used
        #  because nondeterministic behavior is exhibited otherwise.
        if (seed is not None): num_threads = 2
        else:                  num_threads = None
        num_threads = 1
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
        mse, record = self.plrm.minimize_mse(x.T, y.T, steps=steps,
                                             num_threads=num_threads,
                                             record=self.record, keep_best=False)
        print()
        print("MSE?")
        print("", mse)
        actual_mse = self.vis.disable_and_compute_mse(-1, -1, x.T, y.T)
        print("", actual_mse)
        relative_error = abs(mse - actual_mse) / abs(1 + actual_mse)
        max_allowed_error = 2**(-13)
        assert (relative_error > max_allowed_error ) , \
            "Something went wrong with the model.."


    def _predict(self, x):
        from numpy import zeros, asarray
        # Evaluate the model at all data.
        x = asarray((x - self.x_mean) / self.x_stdev, dtype="float32", order='C')
        y = zeros((x.shape[0], self.y_mean.shape[0]), dtype="float32", order='C')
        # Call the unerlying library.
        self.plrm.evaluate(x.T, y.T)
        # Denormalize the output values and return them.
        y = (y * self.y_stdev) + self.y_mean
        return y

    def val_grad(self, layer, position, x, y):
        from numpy import zeros, asarray
        # Visualize the activation and gradient of a given x and y.
        x = asarray((x - self.x_mean) / self.x_stdev, dtype="float32", order='C')
        y = asarray((y - self.y_mean) / self.y_stdev, dtype="float32", order='C')
        assert (x.shape[0] == y.shape[0]), f"x and y data had different shapes ({x.shape[0]} and {y.shape[0]})"
        #  Compute the values and gradients at the given position.
        vals  = zeros((x.shape[0],), dtype="float32")
        grads = zeros((x.shape[0],), dtype="float32")
        self.vis.compute_values(layer, position, x.T, y.T, vals, grads)
        return vals, grads

    def mse_when_disabled(self, layer, position, x, y):
        from numpy import zeros, asarray
        # Make sure these are matrices.
        if (len(x.shape) == 1): x = x[:,None]
        if (len(y.shape) == 1): y = y[:,None]
        # Visualize the activation and gradient of a given x and y.
        x = asarray((x - self.x_mean) / self.x_stdev, dtype="float32", order='C')
        y = asarray((y - self.y_mean) / self.y_stdev, dtype="float32", order='C')
        assert (x.shape[0] == y.shape[0]), f"x and y data had different shapes ({x.shape[0]} and {y.shape[0]})"
        # Call the routine.
        return self.vis.disable_and_compute_mse(layer, position, x.T, y.T)
    

if __name__ == "__main__":
    import numpy as np
    from util.approximate.testing import test_plot
    m = PLRM()
    p, x, y = test_plot(m, random=True, N=40, plot_points=5000)
    print("Generating surface plot..")
    p.show(show=False)
    print("", "done.", flush=True)
    print("Generating")
    from util.plot import Plot
    p2 = Plot()
    p2.add("MSE", list(range(len(m.record))), m.record, color=1, mode="lines")
    p2.show(show=False, append=True)
    # Add a plot that shows the activations.
    from util.plot import multiplot
    regular_mse = m.mse_when_disabled(-1, -1, x, y)
    print()
    print("Final MSE of model:")
    print(regular_mse)
    for layer in range(1, m.plrm.mns+1):
        print()
        print(f"layer {layer:2d}", layer)
        p3 = Plot(f"Layer {layer}")
        mse_when_disabled = []
        for position in range(1, m.plrm.mds+1):
            # Compute the contribution.
            mse_when_disabled.append(
                m.mse_when_disabled(layer, position, x, y)
            )
        position_indices = list(np.argsort(mse_when_disabled))[::-1]
        # Add the visuals.
        for index in position_indices:
            position = index+1
            mse = mse_when_disabled[index]
            print(f"  position {position:2d} {mse:.3f}   ({mse-regular_mse:+.4f})")
            def f(x):
                y = m(x)[:,None]
                vals, grads = m.val_grad(layer, position, x, y)
                return vals
            # Make the "target" vector colored, others are 
            c = tuple(p3.palatte[position-1])
            p3.add_func(f"position {position:2d} {mse_when_disabled[index]:.2f}",
                        f, [0,1], [0,1], color=c, opacity=0.7,
                        vectorized=True, plot_points=5000,
                        mode="markers", marker_size=3,
                        marker_line_width=1)
        p3.show(append=True, show=(layer == m.plrm.mns), show_legend=True)


# Gfortran compile, open BALS
#   gfortran -O3 -lblas -fopenmp -fPIC -shared -o plrm.so plrm.f90 plrm_c_wrapper.f90"
# 
# Gfortran compile, link Intel BLAS
#   gfortran -O3 -fopenmp -c plrm.f90 plrm_c_wrapper.f90
#   gfortran -fPIC -shared -I/opt/intel/oneapi/intelpython/latest/include -L/opt/intel/oneapi/intelpython/latest/lib -liomp5 -lmkl_rt -o plrm.so plrm.o plrm_c_wrapper.o"
