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
        # TODO: Remove this vis code when done.
        self.vis = plrm.vis_plrm

    # def _fit(self, x, y, ds=64, ns=12, steps=10000, seed=None):
    def _fit(self, x, y, ds=16, ns=8, steps=2000, seed=None):
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
        # TODO: Remove most of the following code when testing is done.
        self.record = zeros(steps, dtype="float32", order="F")
        self.logs = zeros((8,ds,ns,steps), dtype="float32", order="F")
        mse, record, logs = self.plrm.minimize_mse(x.T, y.T, steps=steps,
                                                   num_threads=num_threads,
                                                   record=self.record,
                                                   logs=self.logs,
                                                   keep_best=True)
        print()
        print("MSE?")
        print("", mse)
        estimate_mse = self.vis.disable_and_compute_mse(-1, -1, x.T, y.T)
        print("", estimate_mse)
        relative_error = abs(mse - estimate_mse) / (1 + abs(mse))
        max_allowed_error =  2**(-13)
        if (relative_error > max_allowed_error):
            print(f"Something went wrong with the model..\n  {relative_error}\n  {max_allowed_error}")


    def _predict(self, x, embed=False):
        from numpy import zeros, asarray
        # Evaluate the model at all data.
        x = asarray((x - self.x_mean) / self.x_stdev, dtype="float32", order='C')
        if embed:
            y = zeros((x.shape[0], self.plrm.mds), dtype="float32", order='C')
            print(y.shape)
        else:
            y = zeros((x.shape[0], self.y_mean.shape[0]), dtype="float32", order='C')
        # Call the unerlying library.
        self.plrm.evaluate(x.T, y.T)
        # Denormalize the output values and return them.
        if (not embed):
            y = (y * self.y_stdev) + self.y_mean
        return y

    # TODO: Remove the following two methods after testing is done.
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
    p, x, y = test_plot(m, random=True, N=40, plot_points=500) #plot_points=5000)
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
    steps = np.arange(m.logs.shape[-1])
    print()
    print("Final MSE of model:")
    print(regular_mse)
    for layer in range(1, m.plrm.mns+1):
        print()
        print(f"layer {layer:2d}", layer)
        p3 = Plot(f"Layer {layer}") # values WRT to input
        p4 = Plot(f"Layer {layer}") # global contributions
        p5 = Plot(f"Layer {layer} 1D slices") # 1D slices of fits
        # Get the previous layer to this component.
        if (layer == 1):
            prev_layer = m.plrm.input_vecs
        else:
            prev_layer = m.plrm.internal_vecs[:,:,layer-2]
        # Get the next layer to this component.
        if (layer == m.plrm.mns):
            next_layer = m.plrm.output_vecs
        else:
            next_layer = m.plrm.internal_vecs[:,:,layer-1]
        # Track lists with values of interest.
        mse_when_disabled = []
        input_magnnitude = []
        for position in range(1, m.plrm.mds+1):
            # Compute the contribution.
            mse_when_disabled.append(
                m.mse_when_disabled(layer, position, x, y)
            )
            input_magnnitude.append(
                np.linalg.norm( prev_layer[:,position-1] )
            )
        position_indices = list(np.argsort(mse_when_disabled))[::-1]
        # Add the visuals.
        for index in position_indices:
            position = index+1
            mse = mse_when_disabled[index]
            inp_mag = input_magnnitude[index]
            vals, grads = m.val_grad(layer, position, x, y[:,None])
            grad_mag = np.linalg.norm(grads)
            print(f"  position {position:2d} {mse:.3f}   ({mse-regular_mse:+.4f}) [{inp_mag:+.4f}] [{grad_mag:+.4f}]")

            # Add the "value" plot for this node.
            def f(x):
                y = m(x)[:,None]
                vals, grads = m.val_grad(layer, position, x, y)
                return vals
            c = tuple(p3.palette[position-1])
            p3.add_func(f"position {position:2d} ({mse:.4f}) [{inp_mag:.4f}] [{grad_mag:.4f}]",
                        f, [0,1], [0,1], color=c, opacity=0.7,
                        vectorized=True, plot_points=5000,
                        mode="markers", marker_size=3,
                        marker_line_width=1, group=(layer,index))

            # Add all the log values.
            p4.add(f"{position:2d} MSE contribution",  steps, m.logs[0,position-1,layer-1,:], color=1, group=position, mode="lines", show_in_legend=True)
            # p4.add(f"{position:2d} Grad step MSE",     steps, m.logs[1,position-1,layer-1,:], color=2, group=position, mode="lines", show_in_legend=(index==position_indices[0]))
            p4.add(f"{position:2d} Magnitude",         steps, m.logs[2,position-1,layer-1,:], color=3, group=position, mode="lines", show_in_legend=(index==position_indices[0]))
            # p4.add(f"{position:2d} Gradient step mag", steps, m.logs[3,position-1,layer-1,:], color=4, group=position, mode="lines", show_in_legend=(index==position_indices[0]))
            cumulative_angle_change = np.cumsum( m.logs[4,position-1,layer-1,:] )
            p4.add(f"{position:2d} Angle change",      steps, cumulative_angle_change,        color=6, group=position, mode="lines", show_in_legend=(index==position_indices[0]))
            # p4.add(f"{position:2d} Flexure",           steps, m.logs[5,position-1,layer-1,:], color=6, group=position, mode="lines", show_in_legend=(index==position_indices[0]))
            p4.add(f"{position:2d} Shift",             steps, m.logs[6,position-1,layer-1,:], color=7, group=position, mode="lines", show_in_legend=(index==position_indices[0]))

            # vals, grads = m.val_grad(layer, position, x, y[:,None])
            # print(vals.shape)
            # print(grads.shape)
            # exit()
            # p5.add()

            # # Add the gradient plot for this node.
            # def df(x):
            #     y = m(x)
            #     vals, grads = m.val_grad(layer, position, x, y)
            #     return grads
            # c = tuple(p3.palette[position-1] * 1.3)
            # print()
            # print(y.shape)
            # _, df = m.val_grad(layer, position, x, y[:,None])
            # print(df.shape)
            # p3.add(f"position GRAD {position:2d}   ({mse-regular_mse:+.4f})",
            #        *x.T, df, color=c, opacity=0.7,
            #        mode="markers", marker_size=2,
            #        marker_line_width=1)
            

        p4.show(append=True, show=(layer == m.plrm.mns))
        # p3.show(append=True, show=(layer == m.plrm.mns), show_legend=True)

        # multiplot([ p4, p3 ],  append=True, show=(layer == m.plrm.mns))

# Gfortran compile, open BALS
#   gfortran -O3 -lblas -fopenmp -fPIC -shared -o plrm.so plrm.f90 plrm_c_wrapper.f90"
# 
# Gfortran compile, link Intel BLAS
#   gfortran -O3 -fopenmp -c plrm.f90 plrm_c_wrapper.f90
#   gfortran -fPIC -shared -I/opt/intel/oneapi/intelpython/latest/include -L/opt/intel/oneapi/intelpython/latest/lib -liomp5 -lmkl_rt -o plrm.so plrm.o plrm_c_wrapper.o"


# Should be able to reduce a large set of input tokens down to those
# that are "most relevant" to a given prediction task. That ranking
# process should be achievable at very large scale (and allow for the
# idea of "focus" shifting around to relevant parts of the input space
# over sequential predictions).
# 
# A "well spaced" set of tokens have varying types of overlap, but
# also have varying input spatial patterns too: some clustering
# locally, some covering broad ranges globally, ...

# Need to use real normal distribution instead of ATANH, that is not
#   radially symmetric (it is more dense on the corners still.
