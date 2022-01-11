from util.approximate.base import Approximator

SOURCE_FILE = "plrm.f90"
# SOURCE_FILE = "stable_flexure.f90"
# SOURCE_FILE = "stable_relu.f90"
# SOURCE_FILE = "stable_leaky_relu.f90"

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

    def _fit(self, x, y, ds=32, ns=8, steps=1000, seed=None,
             normalize_x=True, normalize_y=True, new_model=True,
             keep_best=True, num_threads=None, **kwargs):
        # If a random seed is provided, then only 2 threads can be used
        #  because nondeterministic behavior is exhibited otherwise.
        if (seed is not None):
            if (num_threads is None):
                num_threads = 2
            else:
                num_threads = min(num_threads, 2)
        # Core numpy utilities.
        from numpy import where, zeros, ones
        # Store the X and Y normalization factors for the model.
        if normalize_x:
            self.x_mean = x.mean(axis=0)
            self.x_stdev = x.std(axis=0)
            self.x_stdev = where(self.x_stdev != 0, self.x_stdev, 1.0)
        else:
            self.x_mean = zeros(x.shape[1], dtype="float32")
            self.x_stdev = ones(x.shape[1], dtype="float32")
        if normalize_y:
            self.y_mean = y.mean(axis=0)
            self.y_stdev = y.std(axis=0)
            self.y_stdev = where(self.y_stdev != 0, self.y_stdev, 1.0)
        else:
            self.y_mean = zeros(y.shape[1], dtype="float32")
            self.y_stdev = ones(y.shape[1], dtype="float32")
        # Normalize the X and Y data for this model.
        x = ((x - self.x_mean) / self.x_stdev).astype("float32")
        y = ((y - self.y_mean) / self.y_stdev).astype("float32")
        # Fit internal piecewise linear regression model.
        if new_model:
            di = x.shape[1]
            do = y.shape[1]
            self.plrm.new_model(di, ds, ns, do)
            self.plrm.init_model(seed=seed)
        self.record = zeros(steps, dtype="float32")
        mse = self.plrm.minimize_mse(x.T, y.T, steps=steps,
                                     num_threads=num_threads,
                                     keep_best=keep_best,
                                     record=self.record, **kwargs)

    def _predict(self, x, embed=False, **kwargs):
        from numpy import zeros, asarray
        # Evaluate the model at all data.
        x = asarray((x - self.x_mean) / self.x_stdev, dtype="float32", order='C')
        if embed:
            y = zeros((x.shape[0], self.plrm.mds), dtype="float32", order='C')
        else:
            y = zeros((x.shape[0], self.y_mean.shape[0]), dtype="float32", order='C')
        # Call the unerlying library.
        self.plrm.evaluate(x.T, y.T, **kwargs)
        # Denormalize the output values and return them.
        if (not embed):
            y = (y * self.y_stdev) + self.y_mean
        return y

    # Save this model to a path.
    def save(self, path):
        import json
        with open(path, "w") as f:
            all_attributes = set(dir(self.plrm))
            to_save = [n[4:] for n in sorted(all_attributes)
                       if (n[:4] == "set_") and (n[4:] in all_attributes)]
            f.write(json.dumps([
                # Create a dictionary of all attributes that can be "set_".
                {
                    n : getattr(self.plrm, n).tolist() if
                    hasattr(getattr(self.plrm, n), "tolist") else
                    getattr(self.plrm, n)
                    for n in to_save
                },
                # Create a dictionary of the known python attributes.
                {
                    "x_mean": self.x_mean.tolist() if hasattr(self, "x_mean") else None,
                    "x_stdev": self.x_stdev.tolist() if hasattr(self, "x_stdev") else None,
                    "y_mean": self.y_mean.tolist() if hasattr(self, "y_mean") else None,
                    "y_stdev": self.y_stdev.tolist() if hasattr(self, "y_stdev") else None,
                    "record": self.record.tolist() if hasattr(self, "record") else None
                }
            ]))

    # Load this model from a path (after having been saved).
    def load(self, path):
        # Read the file.
        from numpy import asarray
        import json
        with open(path, "r") as f:
            model, attrs = json.loads(f.read())
        # Load the attributes of the model.
        for key in model:
            # Try setting all possible model keys.
            try:
                if (model[key] is None): continue
                elif (type(model[key]) is int):
                    setattr(self.plrm, key, model[key])
                else:
                    setattr(self.plrm, key, asarray(model[key], dtype="float32", order="F"))
            # If this is a parameter, ignore it.
            except NotImplementedError: pass
        # Load the attributes of this class.        
        for key in attrs:
            value = attrs[key]
            if (type(value) is list): 
                value = asarray(value, dtype="float32", order="F")
            setattr(self, key, value)
        # Return self in case an assignment was made.
        return self


if __name__ == "__main__":
    import numpy as np
    from util.approximate.testing import test_plot
    m = PLRM()
    # Try saving an untrained model.
    m.save("testing_empty_save.json")
    m.load("testing_empty_save.json")
    # Create the test plot.
    p, x, y = test_plot(m, random=True, N=40, plot_points=500) #plot_points=5000)
    # Try saving the trained model and applying it after loading.
    m.save("testing_real_save.json")
    m.load("testing_real_save.json")
    p.add("Loaded values", *x.T, m(x)+0.05, color=1, marker_size=4)
    # Generate the visual.
    print("Generating surface plot..")
    p.show(show=False)
    print("Generating loss plot..")
    p = type(p)("Mean squared error")
    p.add("MSE", list(range(m.record.shape[0])), m.record, color=1, mode="lines")
    p.show(append=True, show=True)
    print("", "done.", flush=True)
    # Remove the save files.
    import os
    os.remove("testing_empty_save.json")
    os.remove("testing_real_save.json")
