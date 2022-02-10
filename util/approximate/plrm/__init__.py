from util.approximate.base import Approximator

SOURCE_FILE = "lux.f90"
# SOURCE_FILE = "plrm.f90"
# SOURCE_FILE = "stable_flexure.f90"
# SOURCE_FILE = "stable_relu.f90"
# SOURCE_FILE = "stable_leaky_relu.f90"

class PLRM(Approximator):
    def __str__(self): return str(self.model_unpacked())

    def __init__(self, **kwargs):
        # from plrm import plrm
        # self.plrm = plrm
        import os
        from numpy import zeros, ones
        import fmodpy
        this_dir = os.path.dirname(os.path.abspath(__file__))
        source_code = os.path.join(this_dir, SOURCE_FILE)
        plrm = fmodpy.fimport(source_code, blas=True, omp=True, wrap=True,
                              verbose=False, output_dir=this_dir)
        # Store the Fortran module as an attribute.
        self.PLRM = plrm.plrm
        # Set defaults for standard internal parameters.
        self.config = None
        self.model = zeros(0, dtype="float32")
        self.x_mean = zeros(1, dtype="float32")
        self.x_var = ones(1, dtype="float32")
        self.y_mean = zeros(1, dtype="float32")
        self.y_var = ones(1, dtype="float32")
        self.record = zeros(0, dtype="float32")
        # Initialize the attributes of the model that can be initialized.
        self._init_model(**kwargs)


    # Initialize a model, if possible.
    def _init_model(self, **kwargs):
        di = kwargs.pop("di", None)
        ds = kwargs.pop("ds", 32)
        ns = kwargs.pop("ns", 8)
        do = kwargs.pop("do", None)
        mde = kwargs.pop("mde", None)
        mne = kwargs.pop("mne", None)
        num_threads = kwargs.pop("num_threads", None)
        self.seed = kwargs.pop("seed", None)
        self.steps = kwargs.pop("steps", 1000)
        if (None not in {di, ds, ns, do}):
            from numpy import zeros, ones
            self.config = self.PLRM.new_model_config(di, do, ds, ns,
                                mne=mne, mde=mde, num_threads=num_threads)
            # Set any configuration keyword arguments given at initialization
            #   that were not ppassed to "new_model_config".
            for n in ({n for (n,t) in self.config._fields_} & set(kwargs)):
                setattr(self.config, n, kwargs[n])
            # Set all internal arrays and initialize the model.
            self.x_mean = zeros(di, dtype="float32")
            self.x_var = ones(di, dtype="float32")
            self.y_mean = zeros(do, dtype="float32")
            self.y_var = ones(do, dtype="float32")
            self.model = zeros(self.config.total_size, dtype="float32")
            self.PLRM.init_model(self.config, self.model, seed=self.seed)


    # Unpack the model (which is in one array) into it's constituent parts.
    def model_unpacked(self):
        class ModelUnpacked:
            config = self.config
            model  = self.model
            embeddings   = self.model[self.config.sev-1:self.config.eev].reshape(self.config.mde, self.config.mne, order="F")
            input_vecs   = self.model[self.config.siv-1:self.config.eiv].reshape(self.config.mdi, self.config.mds, order="F")
            input_shift  = self.model[self.config.sis-1:self.config.eis].reshape(self.config.mds, order="F")
            input_min    = self.model[self.config.sin-1:self.config.ein].reshape(self.config.mds, order="F")
            input_max    = self.model[self.config.six-1:self.config.eix].reshape(self.config.mds, order="F")
            state_vecs   = self.model[self.config.ssv-1:self.config.esv].reshape(self.config.mds, self.config.mds, self.config.mns-1, order="F")
            state_shift  = self.model[self.config.sss-1:self.config.ess].reshape(self.config.mds, self.config.mns-1, order="F")
            state_min    = self.model[self.config.ssn-1:self.config.esn].reshape(self.config.mds, self.config.mns-1, order="F")
            state_max    = self.model[self.config.ssx-1:self.config.esx].reshape(self.config.mds, self.config.mns-1, order="F")
            output_vecs  = self.model[self.config.sov-1:self.config.eov].reshape(self.config.mds, self.config.mdo, order="F")
            __getitem__  = lambda *args, **kwargs: getattr(*args, **kwargs)
            def __str__(self):
                # Calculate the byte size of this model (excluding python descriptors).
                byte_size = len(self.config._fields_)*4 + self.model.dtype.itemsize*self.model.size
                if (byte_size < 2**10):
                    byte_size = f"{byte_size} bytes"
                elif (byte_size < 2**20):
                    byte_size = f"{byte_size//2**10:.1f}KB"
                elif (byte_size < 2**30):
                    byte_size = f"{byte_size//2**20:.1f}MB"
                else:
                    byte_size = f"{byte_size//2**30:.1f}GB"
                # Create a function that prints the actual contents of the arrays.
                to_str = lambda arr: "\n    " + "\n    ".join(str(arr).split("\n")) + "\n"
                # Provide details (and some values where possible).
                return (
                    f"PLRM model ({self.config.total_size} parameters) [{byte_size}]\n"+
                    f"  input dimension  {self.config.mdi}\n"+
                    f"  output dimension {self.config.mdo}\n"+
                    f"  state dimension  {self.config.mds}\n"+
                    f"  number of states {self.config.mns}\n"+
                     "\n"+
                    f"  embeddings   {self.embeddings.shape}  "+to_str(self.embeddings)+
                    f"  input vecs   {self.input_vecs.shape}  "+to_str(self.input_vecs)+
                    f"  input shift  {self.input_shift.shape} "+to_str(self.input_shift)+
                    f"  input min    {self.input_min.shape}   "+to_str(self.input_min)+
                    f"  input max    {self.input_max.shape}   "+to_str(self.input_max)+
                    f"  state vecs   {self.state_vecs.shape}  "+to_str(self.state_vecs)+
                    f"  state shift  {self.state_shift.shape} "+to_str(self.state_shift)+
                    f"  state min    {self.state_min.shape}   "+to_str(self.state_min)+
                    f"  state max    {self.state_max.shape}   "+to_str(self.state_max)+
                    f"  output vecs  {self.output_vecs.shape} "+to_str(self.output_vecs)
                )
        return ModelUnpacked()


    # Fit this model.
    def _fit(self, x, y, xi=None, normalize_x=True, normalize_y=True,
             new_model=False, **kwargs):
        # Get the number of steps for training.
        steps = kwargs.get("steps", self.steps)
        # If a random seed is provided, then only 2 threads can be used
        #  because nondeterministic behavior is exhibited otherwise.
        seed = kwargs.get("seed", self.seed)
        if (seed is not None):
            num_threads = kwargs.get("num_threads", self.config.num_threads)
            if ((num_threads is None) or (num_threads > 2)):
                import warnings
                warnings.warn("Seeding a PLRM model will deterministically initialize weights, but num_threads > 2 will result in nondeterministic model updates.")
        # Check that the categories meet expectations.
        if (xi is not None):
            minval = xi.min()
            maxval = xi.max()
            num_categories = maxval - minval + 1
        # Verify the shape of the model matches the data.
        di = x.shape[1]
        do = y.shape[1]
        mne = 0
        mde = 0
        if (xi is not None):
            mne = max(self.config.mne, num_categories)
            mde = self.config.mde
            kwargs.update({"mne":mne, "mde":mde})
        # If the shape of the model does not match the data, reinitialize.
        if new_model or (self.config.mdi != di+mde) or (self.config.mdo != do):
            kwargs.update({"di":di, "do":do})
            self._init_model(**kwargs)

        # Set any configuration keyword arguments given at initialization.
        for n in ({n for (n,t) in self.config._fields_} & set(kwargs)):
            setattr(self.config, n, kwargs[n])
        # Verify that the necessary number of embeddings fit into the model.
        if (xi is not None):
            assert (minval == 1), f"Expected minimum category in 'xi' to be 1, received {minval}."
            assert (maxval <= self.config.mne), f"Expected largest category in 'xi' to be at most {self.config.mne}, received {maxval}."
        # Core numpy utilities.
        from numpy import where, zeros, ones
        # Store the X and Y normalization factors for the model.
        if normalize_x:
            self.x_mean = x.mean(axis=0).astype("float32")
            self.x_var = x.var(axis=0).astype("float32")
            self.x_var = where(self.x_var != 0, self.x_var, 1.0)
        else:
            self.x_mean = zeros(x.shape[1], dtype="float32")
            self.x_var = ones(x.shape[1], dtype="float32")
        if normalize_y:
            self.y_mean = y.mean(axis=0)
            self.y_var = y.var(axis=0)
            self.y_var = where(self.y_var != 0, self.y_var, 1.0)
        else:
            self.y_mean = zeros(y.shape[1], dtype="float32")
            self.y_var = ones(y.shape[1], dtype="float32")
        # Normalize the X and Y data for this model.
        x = ((x - self.x_mean) / self.x_var).astype("float32")
        y = ((y - self.y_mean) / self.y_var).astype("float32")
        # Minimize the mean squared error.
        self.record = zeros((steps,2), dtype="float32", order='C')
        self.PLRM.minimize_mse(self.config, self.model,
                               x.T, y.T, steps=steps,
                               xi=xi, record=self.record.T)

    # Make predictions for new data.
    def _predict(self, x, xi=None, embeddings=False, positions=False, **kwargs):
        from numpy import zeros, asarray
        # Evaluate the model at all data.
        x = asarray((x - self.x_mean) / self.x_var, dtype="float32", order='C')
        y = zeros((x.shape[0], self.y_mean.shape[0]), dtype="float32", order='C')
        # Initialize arrays for holding the embeddings and positions of points.
        if (embeddings):
            self.embeddings = zeros((x.shape[0], self.config.mds, self.config.mns), dtype="float32", order='F')
            kwargs["embeddings"] = self.embeddings
        if (positions):
            self.positions = zeros((x.shape[0], self.config.mds, self.config.mns), dtype="int32", order='F')
            kwargs["positions"] = self.positions
        # Check that the categories meet expectations.
        if (xi is not None):
            xi = asarray(xi, dtype="int32", order='C').T
            minval = xi.min()
            maxval = xi.max()
            assert (0 <= minval <= self.config.mne), f"Expected minimum category in 'xi' to be 0 (for unknown) or 1, received {minval}."
            assert (0 <= maxval <= self.config.mne), f"Expected largest category in 'xi' to be at most {self.config.mne}, received {maxval}."
        print("xi:", xi.shape)
        # Call the unerlying library.
        self.PLRM.evaluate(self.config, self.model, x.T, y.T, xi=xi, **kwargs)
        # Denormalize the output values and return them.
        return (y * self.y_var) + self.y_mean


    # Save this model to a path.
    def save(self, path):
        import json
        with open(path, "w") as f:
            # Get the config as a Python type.
            if (self.config is None): config = None
            else: config = {n:getattr(self.config, n) for (n,t) in self.config._fields_}
            # Write the JSON file with Python types.
            f.write(json.dumps({
                # Create a dictionary of the known python attributes.
                "config" : config,
                "model"  : self.model.tolist(),
                "x_mean" : self.x_mean.tolist(),
                "x_var": self.x_var.tolist(),
                "y_mean" : self.y_mean.tolist(),
                "y_var": self.y_var.tolist(),
                "record" : self.record.tolist(),
            }))


    # Load this model from a path (after having been saved).
    def load(self, path):
        # Read the file.
        import json
        from numpy import asarray
        with open(path, "r") as f:
            attrs = json.loads(f.read())
        # Load the attributes of the model.
        for key in attrs:
            value = attrs[key]
            if (type(value) is list): 
                value = asarray(value, dtype="float32")
            setattr(self, key, value)
        # Convert the the dictionary configuration into the correct type.
        if (type(self.config) is dict):
            self.config = self.PLRM.MODEL_CONFIG(**self.config)
        # Return self in case an assignment was made.
        return self


if __name__ == "__main__":
    import numpy as np
    from util.approximate.testing import test_plot
    layer_dim = 20
    num_layers = 4
    m = PLRM()
    # Try saving an untrained model.
    m.save("testing_empty_save.json")
    m.load("testing_empty_save.json")
    m = PLRM(di=2, ds=layer_dim, ns=num_layers, do=1, mde=4, mne=2, seed=0,
             steps=1000, initial_output_scale=1.0, num_threads=1,
             revert_to_best=False)
    from util.approximate import LSHEP

    # Create the test plot.
    p, x, y = test_plot(LSHEP(), random=True, N=100, plot_points=1000) #plot_points=5000)
    # # Try saving the trained model and applying it after loading.
    # m.save("testing_real_save.json")
    # m.load("testing_real_save.json")
    # p.add("Loaded values", *x.T, m(x)+0.05, color=1, marker_size=4)

    print("Building model..")
    # Fit a new model that has categories for two different functions.
    m = PLRM(di=2, ds=layer_dim, ns=num_layers, do=1, mde=4, mne=2, seed=0, steps=1000, num_threads=1)
    all_x = np.concatenate((x, x), axis=0)
    all_y = np.concatenate((y, np.cos(np.linalg.norm(x,axis=1))), axis=0)
    all_xi = np.concatenate((np.ones(len(x), dtype="int32"),2*np.ones(len(x), dtype="int32"))).reshape((1,-1))
    print("all_x.shape: ",all_x.shape)
    print("all_y.shape: ",all_y.shape)
    print("all_xi.shape: ",all_xi.shape)
    print(m)
    m.fit(all_x, all_y, xi=all_xi, num_threads=1)
    print("m.config.mne: ", m.config.mne)
    print("Evaluating y1..")
    xi1 = np.ones((len(x),1),dtype="int32")
    y1 = m(x, xi=xi1)
    print("y1.shape: ", y1.shape)
    print("Evaluating y2..")
    y2 = m(x, xi=2*xi1)
    print("Adding to plot..")
    p.add("t0", *all_x.T, all_y, color=2)
    p.add("y1", *x.T, y1, color=2)
    p.add("y2", *x.T, y2, color=3)

    # Generate the visual.
    print("Generating surface plot..")
    p.show(show=False)
    print("Generating loss plot..")
    p = type(p)("Mean squared error")
    # Rescale the columns of the record for visualization.
    record = m.record
    record /= record.max(axis=0)
    p.add("MSE", list(range(record.shape[0])), record[:,0], color=1, mode="lines")
    p.add("Step sizes", list(range(record.shape[0])), record[:,1], color=2, mode="lines")
    p.show(append=True, show=True)
    print("", "done.", flush=True)
    # Remove the save files.
    import os
    os.remove("testing_empty_save.json")
    os.remove("testing_real_save.json")

    VISUALIZE_EMBEDDINGS = False
    if VISUALIZE_EMBEDDINGS:
        # Create a blank-slate model, to see what it looks like.
        n = PLRM(di=2, ds=layer_dim, ns=num_layers, do=1, seed=0)
        from util.random import latin
        nx = latin(2000, 2)
        #   visialize the initial model at all layers.
        n(nx, embeddings=True, positions=True)
        embeddings = n.embeddings
        for i in range(num_layers):
            p = type(p)(f"Layer {i+1} initializations")
            p.add("Data", *x.T, y / 4, shade=True)
            for j in range(layer_dim):
                vals = embeddings[:,j,i]
                vals -= vals.min()
                if (vals.max() > 0):
                    vals /= vals.max()
                else:
                    vals -= .1
                p.add(f"Component {j+1}", *nx.T, vals, marker_size=3, use_gradient=True)
            p.show(append=(i>0), show=(i in {0, num_layers-1}))

        #   visualize the trained model at all layers
        n.fit(x, y)
        n(nx, embeddings=True, positions=True)
        embeddings = n.embeddings
        for i in range(num_layers):
            p = type(p)(f"Layer {i+1} final form")
            p.add("Data", *x.T, y / 4, shade=True)
            for j in range(layer_dim):
                vals = embeddings[:,j,i]
                vals -= vals.min()
                if (vals.max() > 0):
                    vals /= vals.max()
                else:
                    vals -= .1
                p.add(f"Component {j+1}", *nx.T, vals, marker_size=3, use_gradient=True)
            p.show(append=(i>0), show=(i in {0, num_layers-1}))


# 2022-02-09 20:25:55
# 
            ######################################################################
            # # Set any configuration keyword arguments given at initialization. #
            # for n in ({n for (n,t) in self.config._fields_} & set(kwargs)):    #
            #     setattr(self.config, n, kwargs[n])                             #
            ######################################################################
