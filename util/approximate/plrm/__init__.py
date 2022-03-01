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
                    (f"  embedding dimension  {self.config.mde}\n"+
                     f"  number of embeddings {self.config.mne}\n"
                     if self.config.mne > 0 else "")+"\n"+
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
        # If this model isn't initialized, then do that.
        if (self.config is None):
            self._init_model(di=x.shape[1], do=y.shape[1])
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
        # Verify the shape of the model matches the data.
        di = x.shape[1]
        do = y.shape[1]
        mne = 0
        mde = 0
        # Check that the categories meet expectations.
        if (xi is not None):
            xi = np.asarray(xi, dtype="int32", order="C").T
            minval = xi.min()
            maxval = xi.max()
            num_categories = int(maxval - minval + 1)
            mne = max(self.config.mne, num_categories)
            mde = self.config.mde
            if (mne > 0): kwargs.update({"mne":mne})
            if (mde > 0): kwargs.update({"mde":mde})
        # If the shape of the model does not match the data, reinitialize.
        if (new_model or (self.config.mdi != di+mde) or
            (self.config.mdo != do) or ((mne > 0) and (mde <= 0))):
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
            assert (xi.shape[1] == x.shape[0]), f"Provided 'x' {x.shape} and 'xi' {xi.T.shape} do not match."
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
    from util.stats import pca
    from tlux.plot import Plot
    from tlux.random import well_spaced_ball

    # Define a wrapper convenience function for computing principal components.
    from sklearn.decomposition import PCA        
    def pca(x, num_components=None):
        pca = PCA(n_components=num_components)
        if (num_components is None): num_components = min(*x.shape)
        else: num_components = min(num_components, *x.shape)
        pca.fit(x)
        return pca.components_, pca.singular_values_

    # TODO:
    #  - bounds on curvature based on the layer in the model
    #  - way to increase basis function diversity (different shift distribution?)
    #    it gets very low with depth of model using uniform and random ortho on sphere
    # 

    layer_dim = 100
    num_layers = 10

    TEST_SPECTRUM = True
    if TEST_SPECTRUM:
    #     layer_dim = 40
    #     num_layers = 20
        n = 10000
        d_in = 3
        d_out = 5
        # Control array print options for readability.
        np.set_printoptions(edgeitems=50, linewidth=100000, 
            formatter=dict(float=lambda v: f"{v:.2f}"))

        # Plot the values (project onto 3D if there are more).
        visualize = False
        def plot_vals(x, name="", d=3, show=False):
            if (not visualize): return
            if (x.shape[1] > d):
                vecs, mags = pca(x, num_components=d)
                x = np.matmul(x, vecs.T)
            # Create the plot.
            p = Plot(name)
            p.add("", *x.T, color=1, marker_size=4, marker_line_width=1)
            p.plot(show=show, append=True, show_legend=False)
        # Initialize a model and data.
        m = PLRM(di=d_in, ds=layer_dim, ns=num_layers, do=d_out, initial_shift_range=2.0)
        x = well_spaced_ball(n, d_in)
        y = np.vstack([
            np.cos(2*np.pi*(i+1)*np.linalg.norm(x,axis=1))
            for i in range(d_out)]).T
        # Show the input points.
        plot_vals(x, "Input")
        # Generate embeddings for all points.
        m(x, embeddings=True)
        for l in range(num_layers):
            emb = m.embeddings[:,:,l]
            print(l)
            print("  min", emb.min().round(2))#, emb.min(axis=0))
            print("  max", emb.max().round(2))#, emb.max(axis=0))
            print("  point spectrum:", pca(emb)[1])
            print("  basis spectrum:", pca(emb.T)[1])
            print()
            plot_vals(emb, f"Layer {l+1}", show=(l==num_layers-1))
        # Fit the model, then look at the new spectrum.
        print('-'*70)
        m.fit(x, y, steps=2000)
        m(x, embeddings=True)
        print('-'*70)
        #Show the spectrum.
        for l in range(num_layers):
            emb = m.embeddings[:,:,l]
            print(l)
            print("  min", emb.min().round(2))#, emb.min(axis=0))
            print("  max", emb.max().round(2))#, emb.max(axis=0))
            print("  point spectrum:", pca(emb)[1])
            print("  basis spectrum:", pca(emb.T)[1])
            print()
            plot_vals(emb, f"Fit Layer {l+1}", show=(l==num_layers-1))



    TEST_SAVE_LOAD = False
    if TEST_SAVE_LOAD:
        # Try saving an untrained model.
        m = PLRM()
        m.save("testing_empty_save.json")
        m.load("testing_empty_save.json")
        m = PLRM(di=2, ds=layer_dim, ns=num_layers, do=1, mde=4, mne=2, seed=0,
                 steps=1000, initial_output_scale=1.0, num_threads=2)
        # Create the test plot.
        p, x, y = test_plot(m, random=True, N=100, plot_points=1000) #plot_points=5000)
        # Try saving the trained model and applying it after loading.
        m.save("testing_real_save.json")
        m.load("testing_real_save.json")
        p.add("Loaded values", *x.T, m(x)+0.05, color=1, marker_size=4)
        p.plot(show=False)
        # Add plot showing the training loss.
        print("Generating loss plot..")
        p = type(p)("Mean squared error")
        # Rescale the columns of the record for visualization.
        record = m.record
        p.add("MSE", list(range(record.shape[0])), record[:,0], color=1, mode="lines")
        p.add("Step sizes", list(range(record.shape[0])), record[:,1], color=2, mode="lines")
        p.show(append=True, show=True)
        print("", "done.", flush=True)
        # Remove the save files.
        import os
        try: os.remove("testing_empty_save.json")
        except: pass
        try: os.remove("testing_real_save.json")
        except: pass


    TEST_INT_INPUT = False
    if TEST_INT_INPUT:
        print("Building model..")
        # Fit a small throw-away model to get the test plot points.
        m = PLRM(di=2, ds=1, ns=1, steps=10)
        p, x, y = test_plot(m, random=True, N=100)
        # Initialize a new model.
        m = PLRM(di=2, ds=layer_dim, ns=num_layers, do=1, mne=2, seed=0, steps=1000, num_threads=2)
        all_x = np.concatenate((x, x), axis=0)
        all_y = np.concatenate((y, np.cos(np.linalg.norm(x,axis=1))), axis=0)
        all_xi = np.concatenate((np.ones(len(x)),2*np.ones(len(x)))).reshape((-1,1)).astype("int32")
        print("fitting model that should reshape")
        m.fit(all_x, all_y, xi=all_xi)
        # Create an evaluation set that evaluates the model that was built over two differnt functions.
        xi1 = np.ones((len(x),1),dtype="int32")
        y1 = m(x, xi=xi1)
        y2 = m(x, xi=2*xi1)
        print("Adding to plot..")
        p = type(p)()
        p.add("xi=1 true", *x.T, all_y[:len(all_y)//2], color=0)
        p.add("xi=2 true", *x.T, all_y[len(all_y)//2:], color=1)
        p.add_func("xi=1", lambda x: m(x, np.ones(len(x), dtype="int32").reshape((1,-1))), [0,1], [0,1], color=3, mode="markers", shade=True)
        p.add_func("xi=2", lambda x: m(x, 2*np.ones(len(x), dtype="int32").reshape((1,-1))), [0,1], [0,1], color=2, mode="markers", shade=True)

        # Generate the visual.
        print("Generating surface plot..")
        p.show(show=False)
        print("Generating loss plot..")
        p = type(p)("Mean squared error")
        # Rescale the columns of the record for visualization.
        record = m.record
        p.add("MSE", list(range(record.shape[0])), record[:,0], color=1, mode="lines")
        p.add("Step sizes", list(range(record.shape[0])), record[:,1], color=2, mode="lines")
        p.show(append=True, show=True)
        print("", "done.", flush=True)


    VISUALIZE_EMBEDDINGS = False
    if VISUALIZE_EMBEDDINGS:
        # Create a blank-slate model, to see what it looks like.
        n = PLRM(di=2, ds=layer_dim, ns=num_layers, do=1, seed=0, num_threads=2)
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
