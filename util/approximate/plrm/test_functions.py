# Define a set of test functions for approximation.
import numpy as np
# Imports for testing.
import os, time, json, gzip
from itertools import product
# Generating well-spaced points over the unit sphere (or ball) in any dimension.
from util.random import latin_sphere


REAL64_PRECISION = 2 ** (-52)
SQRT_REAL64_PRECISION = 2 ** (-26)


# Given a list of values, resize the list of values to the desired
#  size by linearly interpolating the percentiles of the data.
def linsize(values, new_size):
    return np.percentile(values, np.linspace(0, 100, new_size))


# Given a set of y values representing a 1D function, smooth the
#  function values by convolving a normal distirbution over the data.
def smoothy(x, y, stdev=1/30):
    if (stdev <= 0): return y
    x = np.asarray(x)
    y = np.asarray(y)
    new_y = np.zeros(y.shape)
    stdev *= max(SQRT_REAL64_PRECISION, x.max() - x.min())
    for i in range(y.size):
        # Generate convex weights from a normal distribution.
        weights = np.exp(
            -((x - x[i]) / stdev)**2 / 2
        ) / (stdev * np.sqrt(2*np.pi))
        # Adjust the weights at the ends of the interval to account 
        #   for clipping using the normal cumulative distirbution.
        left_tail_weight = (1.0 + np.math.erf((x[0]-x[i]) / (stdev * np.sqrt(2.0)))) / 2.0
        right_tail_weight = 1.0 - (1.0 + np.math.erf((x[-1]-x[i]) / (stdev * np.sqrt(2.0)))) / 2.0
        # Rescale the tail weights so their contribution is 
        #   proportional to their probability relative to the body.
        body_weight = 1.0 - left_tail_weight - right_tail_weight
        weight_scale = weights.sum() / max(SQRT_REAL64_PRECISION, body_weight)
        left_tail_weight *= weight_scale
        right_tail_weight *= weight_scale
        weights[0] += left_tail_weight
        weights[-1] += right_tail_weight
        # Make the weights convex.
        weights /= max(SQRT_REAL64_PRECISION, weights.sum())
        # Create the new y value (by convolving the normal distribution).
        new_y[i] = np.dot(y, weights)
    return new_y


# Produce a CDF fit of some values.
def cdf_fit(values, smooth=False, smooth_n=200, stdev=1/10, eps=SQRT_REAL64_PRECISION):
    assert len(values) > 0, "Provided values array has length 0."
    cdf_x = np.array(values)
    cdf_x.sort()
    # Ensure that all x are at least "eps" apart from each other.
    for i in range(len(cdf_x)-1):
        cdf_x[i+1] = max(cdf_x[i+1], cdf_x[i]+eps)
    cdf_y = np.linspace(0, 1, cdf_x.size)
    # Smooth the fit, if desired.
    if (smooth and (smooth_n > 0)):
        new_x = linsize(cdf_x, smooth_n)
        new_y = np.interp(new_x, cdf_x, cdf_y)
        new_y = smoothy(new_x, new_y, stdev=stdev)
        new_y -= new_y.min()
        new_y /= new_y.max()
        cdf_x, cdf_y = new_x, new_y
    # Generate the fit and its inverse.
    min_max = [cdf_x.min(), cdf_x.max()]
    fit = lambda *x: np.interp(x[0], cdf_x, cdf_y) if (len(x) > 0) else min_max
    fit.inverse = lambda y: np.interp(y, cdf_y, cdf_x)
    # Generate the approximate deirvative of the fit.
    x_widths = (cdf_x[1:] - cdf_x[:-1])
    mid_slopes_x = (cdf_x[:-1] + cdf_x[1:]) / 2
    mid_slopes_y = (cdf_y[1:] - cdf_y[:-1]) / x_widths
    slopes_y = cdf_x.copy()
    slopes_y[1:-1] = (mid_slopes_y[1:] + mid_slopes_y[:-1]) / 2
    slopes_y[0] = mid_slopes_y[0]
    slopes_y[-1] = mid_slopes_y[-1]
    # Normalize the area under the curve to be 1.
    slope_area = x_widths * (slopes_y[1:] + slopes_y[:-1]) / 2
    slopes_y /= slope_area.sum()
    fit.derivative = lambda x: np.interp(x, cdf_x, slopes_y)
    # Return the fit.
    return fit


# Test the CDF fit function.
def _test_linsize_smoothy_cdf_fit():
    from util.plot import Plot
    np.random.seed(0)
    # Generate a piecewise linear fit of a uniform distribution.
    data = np.random.uniform(size=(10000))
    cdf = cdf_fit(data)
    p = Plot()
    p.add_histogram("data", data, color=(0,0,0,0.3))
    p.add_func("CDF fit", cdf, cdf(), color=1)
    p.add_func("CDF inverse", cdf.inverse, [0, 1], color=2)
    p.show(show=False)
    # Generate a smoothed fit of a normal distribution.
    data = np.random.normal(size=(100))
    cdf = cdf_fit(data, smooth=True)
    p = Plot()
    p.add_histogram("data", data, color=(0,0,0,0.3))
    p.add_func("CDF fit", cdf, cdf(), color=1)
    p.add_func("CDF inverse", cdf.inverse, [0, 1], color=2)
    p.add_func("PDF", cdf.derivative, cdf(), color=0, line_width=0, fill="tozeroy")
    p.show(append=True)


# Return a random CDF.
def rand_cdf(nodes=3, power=1.0, smooth=200):
    # Generate random steps in the x and y direction (that sum to 1).
    cdf_x = np.linspace(0, 1, nodes+2)
    # Randomly assign a probability to each bin (this does not produce diverse CDFs).
    cdf_y = np.random.random(size=cdf_x.size) ** power
    cdf_y[0] = 0.0
    cdf_y /= cdf_y.sum()
    cdf_y = cdf_y.cumsum()
    # Smooth the CDF y values if desired.
    if smooth > 0:
        new_x = linsize(cdf_x, smooth)
        new_y = np.interp(new_x, cdf_x, cdf_y)
        new_y = smoothy(new_x, new_y)
        new_y -= new_y.min()
        new_y /= new_y.max()
        # Assign the new CDF points.
        cdf_x, cdf_y = new_x, new_y
    # Generate the CDF fit and return it.
    cdf = lambda x: np.interp(x, cdf_x, cdf_y)
    cdf.inverse = lambda x: np.interp(x, cdf_y, cdf_x)
    return cdf


# Verify that the CDF functions generated have nice behaviors.
def _test_rand_cdf():
    from util.stats import plot_percentiles
    from util.plot import Plot
    p = Plot()
    # Generate samples from a lot of random CDFs.
    x_points = np.linspace(0, 1, 1000)
    y_points = []
    samples = 1000
    for i in range(samples):
        print(f"\r{i:5d} : {samples:5d}", end="")
        y_points.append(
            rand_cdf()(x_points)
        )
    print(" "*50)
    percentiles = list(range(0,101,10))
    y_points = np.asarray(y_points).T
    # Plot all of the percentiles of those CDFs to see the
    #  distribution of distributions that is generated.
    p.color_num += 1
    p = plot_percentiles(p, "Random CDFs", x_points, y_points,
                         percentiles=percentiles)
    # Add some examples for good measure.
    for i in range(5):
        f = rand_cdf()
        p.add_func(f"Sample {i+1}", f, [0, 1])
    # Generate the plot.
    p.show()


# Compute the cosine of the 2-norm of the x points (row vectors).
def cos_norm(x, radius=1, peaks=2, power=2):
    return np.cos(np.linalg.norm(x / radius, axis=1)**power
                  * max(0, peaks * 2 - 1) * np.pi).reshape((-1,1))


# Test the pure function that takes the cosine of the norm of x.
def _test_cos_norm():
    from util.plot import Plot
    p = Plot()
    n = 5000
    d = 2
    x = latin_sphere(n, d, inside=True)
    y = cos_norm(x)
    p.add("cos 2-norm", *x.T, y[:,0], use_gradient=True)
    # p.add_func("numpy", x1_cdf, x1_cdf())
    # p.add_func("lcg", x2_cdf, x2_cdf())
    p.show(z_range=[-3,5])


# A pure approximation problem. Should serve as a viable test for
# approximation algorithms. Returns a test approximation function that
# returns two arrays, x and y, given x inputs for approximation.
# 
#   d (integer)
#      Dimension of the function input.
#   df (optional integer)
#      Dimension of the function output. Defaults to 1.
#   df_intrinsic (optional integer)
#      Intrinsic dimension of the output. I.e., number of principal
#      components with nonzero variance. Default is max(1, df // 2).
#   shift (real)
#      Magnitude of the random shift applied to the input data for
#      each intrinsic dimension of output.
#   distort (optional boolean)
#      If True, apply a random distribution transformation that 
#      distorts the axis-aligned input point distribution. 
#      Defaults to False.
#   skew (optional boolean)
#      If True, apply a random linear transformation that rotates and
#      skews the input points. Defaults to False.
#   ortho (optional boolean)
#      This flag determines whether orthogonality is enforced over the
#      outputs of the underlying function in the intrinsic dimension,
#      necessary for orthogonality when distort=True of skew=True.
#      Default value is True.
#      WARNING: Setting this to True will result in different spatial
#               function evaluations depending on the number of inputs.
# 
def pure(d, df=1, df_intrinsic=None, shift=0.2, complexity=2,
         distort=False, skew=False, ortho=True, seed=0):
    np.random.seed(seed)
    # Set the intrinsic dimension of the output if it was not provided.
    if (df_intrinsic is None):
        df_intrinsic = max(1, df // 2)
    # If output is in higher dimension than the intrinsic, it will
    #  be lifted with a linear projection onto well spaced vectors.
    projection = None
    if (df != df_intrinsic):
        projection = latin_sphere(df, df_intrinsic).T
    # Add some centered shift to the input data of the function.
    shift_vecs = np.vstack((
        np.zeros((1, d)),
        latin_sphere(df_intrinsic-1, d) * shift
    ))
    # Apply a random distribution distortion to all input components.
    cdfs = None
    if (distort):
        cdfs = [rand_cdf() for i in range(d)]
    # Apply a random linear projection to x to "skew" it.
    if (skew):
        skew = latin_sphere(d, d).T
    else:
        skew = None
    # Create a function with orthogonal outputs (in the intrinsic dimension).
    def pure_func(x, complexity=complexity, shift=shift, shift_vecs=shift_vecs,
                  projection=projection, cdfs=cdfs, skew=skew, ortho=ortho):
        # Check quality of x input.
        assert (len(x.shape) == 2), f"Expected x with shape (n,{d}) but received x with shape {x.shape}."
        assert (x.shape[1] == d), f"Expected row vectors in x with dimension {d} but received x with shape {x.shape}."
        # Apply a distribution distortion to input components
        #   assuming ideal x component range is [-1,1].
        if (cdfs is not None):
            x = x.copy()
            for i in range(d):
                x[:,i] = cdfs[i].inverse((x[:,i] + 1) / 2) * 2 - 1
        # Apply a linear skew to the input.
        if (skew is not None):
            x = np.matmul(x, skew)
        # Compute a function over inputs that has orthogonal outputs.
        y = [cos_norm(x + shift_vecs[i], radius=1+shift, peaks=complexity) 
             for i in range(df_intrinsic)]
        y = np.asarray(y).reshape((df_intrinsic, -1)).T
        # Make the y outputs orthogonal in the intrinsic dimension.
        if (ortho and (df_intrinsic > 1)):
            old_y, (y, _) = y, np.linalg.qr(y)
            numerator = 1.0 if (np.dot(old_y[:,0], y[:,0]) >= 0) else -1.0
            y_scale = numerator / np.abs(y[:,0]).max()
            y *= y_scale
        # Project the output into a higher dimension.
        if (projection is not None):
            y = np.matmul(y, projection)
        # Return the final function output.
        return x, y
    # Reset the random seed for later processes.
    np.random.seed()
    # Return the pure function.
    return pure_func


# Test the random uniform by plotting the CDF.
def _test_pure(normalize=False):
    from util.plot import Plot
    p = Plot()
    n = 500
    d = 2
    df = 1
    f = pure(d, df=df)
    x = latin_sphere(n, d, inside=True)
    new_x, y = f(x)
    # Zero mean and unit variance the data like the neural network would.
    if normalize:
        x -= x.mean(axis=0)
        x /= x.var(axis=0)
        y -= y.mean(axis=0)
        y /= y.var(axis=0)
    # Visualize the pure function components.
    for i in range(df):
        p.add(f"pure {i} x", *x.T, y[:,i] * 0.0 + y[:,i].min() - 0.7, marker_size=2, color=(200,200,200),
              marker_line_color=(0,0,0), marker_line_width=1, group=i, show_in_legend=False)
        p.add(f"pure {i}", *x.T, y[:,i], marker_size=4, use_gradient=True, marker_line_width=1, group=i)
    p.show(z_range=[-3,5])


# Create a function that samples an underlying surface. The ground
#  truth prediction task has two parts, identify the type of
#  underlying surface and identify the center of that type.
def sampled(d, df=4, dp=2, df_intrinsic=None,
            max_shift=0.5, power=3, return_info=False, seed=0):
    # By default, make all of the "modes" that can be observed maximally different.
    if (df_intrinsic is None):
        df_intrinsic = 2 * df
    # Define the output dimension of the pure function based on the total output dimension.
    dp = min(dp, df)
    fdf = max(1, df - dp)
    # Generate a "truth function" that creates orthogonal output
    #  components that can be combined to form "category predictions".
    f = pure(dp, df=fdf, df_intrinsic=df_intrinsic, seed=seed)
    # Generate positions for all of the dimensions in some space.
    np.random.seed(seed)
    positions = latin_sphere(d, dp, inside=True)
    np.random.seed()
    # Rescale the function output so that the final dimension can be matched.
    def sampled_func(inputs, d=d, dp=dp, df=df, fdf=fdf, max_shift=max_shift,
                     power=power, f=f, positions=positions, seed=seed):
        assert len(inputs.shape) == 2, f"Expected input array to have shape (n,{d}), received {inputs.shape}."
        assert inputs.shape[1] == d,   f"Expected input array to have shape (n,{d}), received {inputs.shape}."
        np.random.seed(seed)
        # Get the number of points.
        n = inputs.shape[0]
        # Assign random center offsets to all of the points (the y values).
        offsets = latin_sphere(n, dp, inside=True) * max_shift
        # Generate a set of evaluation points (at all positions, for all offsets).
        x = np.concatenate([(positions - o) for o in offsets])
        # Evaluate the function at all positions, this will produce
        #   fdf channels of output for every (position,offset) pair.
        new_x, y = f(x)
        # Generate a set of weights, this will be used to blend the modes
        #   of the output from the function into single "observations".
        #   Ideally, the "most observable" output is the predicted one.
        if (fdf + dp == df):
            weights = (latin_sphere(n, fdf) + 1) ** power
        else:
            weights = np.ones((n, fdf))
        weights = (weights.T / np.sum(weights, axis=1)).T
        # Overwrite the inputs, collapsing the observations using the weights.
        for i in range(n):
            inputs[i,:] = np.dot(y[i*d:(i+1)*d], weights[i])
        # The outputs are the offset terms and weights concatenated,
        #  both are something that should be predicted by an 
        #  approximator. Unless there is not room in the output.
        if (fdf + dp == df):
            outputs = np.concatenate((offsets, weights), axis=1)
        else:
            outputs = offsets
        # Return the pairs of points.
        np.random.seed()
        return inputs, outputs
    # Store some configurations in accessible attributes.
    sampled_func.f = f
    sampled_func.positions = positions
    # Return the function.
    return sampled_func


# Test the function that mimics the sampling behavior.
def _test_sampled():
    from util.plot import Plot
    n = 10
    d = 3000
    f = sampled(d, return_info=True)
    positions = f.positions
    x = latin_sphere(n, d)
    _, y = f(x)
    p = Plot()
    for i in range(n):
        p.add(f"post {i+1}", [y[i,0],y[i,0]], [y[i,1], y[i,1]], [-2, 2],
              mode="lines", line_width=4, color=1, group=i, show_in_legend=False)
        p.add(f"Offset {i+1} ({y[i,0]:.3f}, {y[i,1]:.3f}) [{y[i,2:].round(3)}]",
              *positions.T, x[i], use_gradient=True,
              marker_size=4, marker_line_width=1, group=i)
    p.show(z_range=[-3,5])


# Test all functions.
def _test_all():
    _test_linsize_smoothy_cdf_fit()
    _test_rand_cdf()
    _test_cos_norm()
    _test_pure()
    _test_sampled()


# Define functions that map from a configuration to identifier strings.
def model_str(M):
    return str(M).split("'")[1].replace(".", "-").split("-")[-1]
def function_str(f):
    return str(str(f).split()[1:2])[2:-2]
def config_str(config):
    M, f, n, d, df, ds, ns, s, test_size, func_seed, model_seed = config
    return f"{func_seed:02d}-{model_seed:02d} {function_str(f):8s} {n:5d} {d:5d} {df} {model_str(M)}-{ds:02d}-{ns:02d} ({s} {test_size})"
def test_id_str(config):
    M, f, n, d, df, ds, ns, s, test_size, func_seed, model_seed = config
    return f"{model_seed:04d}-{func_seed:04d}-{function_str(f)}"

# Given a model class, instantiate and test various problems.
def test_model_class(M, suffix="", test_size=1000,
                     func_seeds=list(range(1,5+1)),
                     model_seeds=list(range(1,10+1)),
                     funcs=[pure, sampled],
                     num_points=[50, 100, 200, 500],
                     dimension_in=[2, 10, 100],
                     dimension_out=[1, 5],
                     model_internal_dimension=[32],
                     model_internal_states=[8],
                     model_fit_steps=[3000]):
    # Run a single test given a model class M. Return the results.
    def run_test(config):
        # Extract the configuration variables.
        M, f, n, d, df, ds, ns, s, test_size, func_seed, model_seed = config
        np.random.seed(func_seed)
        # Generate random test data.
        func = f(d, df=df, seed=func_seed)
        x = latin_sphere(n+test_size, d, inside=True)
        new_x, y = func(x)
        # Componentwise zero mean and unit variance all data.
        x_mean = x.mean(axis=0)
        x_stdev = x.std(axis=0)
        y_mean = y.mean(axis=0)
        y_stdev = y.std(axis=0)
        x = (x - x_mean) / np.where(x_stdev > 0, x_stdev, 1.0)
        y = (y - y_mean) / np.where(y_stdev > 0, y_stdev, 1.0)
        # Assign training and testing points.
        train_inds = np.arange(n+test_size)
        np.random.shuffle(train_inds)
        train_inds, test_inds = train_inds[:n], train_inds[n:]
        # Initialize the model.
        model = M(di=d, do=df, ds=ds, ns=ns, seed=model_seed)
        # Compute the error of the model before the fit operation.
        output = model(x[test_inds])
        prefit_error = ((output - y[test_inds])**2).sum(axis=1).mean()
        # Fit the model.
        fit_start = time.time()
        model.fit(x[train_inds], y[train_inds], new_model=False,
                  normalize_x=False, normalize_y=False, steps=s)
        fit_time = time.time() - fit_start
        # Predict with the model.
        predict_start = time.time()
        output = model(x[test_inds])
        predict_time = time.time() - predict_start
        # Measure the error (mean squared and percentiles).
        error = ((output - y[test_inds])**2).sum(axis=1)
        error_deciles = np.percentile(error, list(range(0,101,10)))
        # Return the test results.
        return {
            "model": model_str(M),
            "function": function_str(f),
            "func_seed": func_seed,
            "model_seed": model_seed,
            "n": n,
            "d": d,
            "df": df,
            "ds": ds,
            "ns": ns,
            "steps": s,
            "test_size": test_size,
            "fit_time": fit_time,
            "predict_time": predict_time,
            "initial_mean_squared_error": float(prefit_error),
            "mean_squared_error": float(error.mean()),
            "error_deciles": error_deciles.tolist(),
            "fit_record": getattr(model, "record", np.zeros(0)).tolist(),
        }
        
    # Load in any existing data that's already been collected.
    data_dir = os.path.join("data", f"testing_{model_str(M)}{suffix}")
    os.makedirs(data_dir, exist_ok=True)
    test_configs = set()
    for test_id in os.listdir(data_dir):
        with gzip.open(os.path.join(data_dir, test_id), "rb") as f:
            test_results = json.loads(str(f.read(), "utf8"))
        test_configs |= set(test_results)
    # Make all of the configurations for testing, skipping already covered tests.
    configs = []
    for (fs, ms, f, n, d, df, ds, ns, s) in product(
            func_seeds, model_seeds, funcs, num_points, dimension_in, 
            dimension_out, model_internal_dimension, 
            model_internal_states, model_fit_steps, 
    ):
        config = (M, f, n, d, df, ds, ns, s, test_size, fs, ms)
        conf_str = config_str(config)
        if (conf_str in test_configs):
            print(conf_str)
        else:
            configs.append(config)
    # Execute all of the tests.
    previous_time = time.time()
    previous_id = None
    test_results = {}
    for config in configs:
        # If the ID of this test doesn't match the previous ID, then
        #  save all previous results to a file and reset the result list..
        test_id = test_id_str(config)
        if (test_id != previous_id):
            if (previous_id is not None):
                with gzip.open(os.path.join(data_dir, f"{previous_id}.json.gz"), "wb") as f:
                    f.write(bytes(json.dumps(test_results), "utf8"))
            # Check for existing data from the new test config, load it if it exists.
            test_data_path = os.path.join(data_dir, f"{test_id}.json.gz")
            if os.path.exists(test_data_path):
                with gzip.open(test_data_path, "rb") as f:
                    test_results = json.loads(str(f.read(), "utf8"))
            else:
                test_results = {}
            # Reset the previous ID and the previous time.
            previous_id = test_id
            previous_time = time.time()
        # Run this new test.
        result = run_test(config)
        # Store test result.
        conf_str = config_str(config)
        test_results[conf_str] = result
        # Print an update.
        print(f"{conf_str}  {time.ctime()}  ({time.time() - previous_time:0.2f} seconds)")
        previous_time = time.time()
    # Save the last test iteration.
    with gzip.open(os.path.join(data_dir, f"{previous_id}.json.gz"), "wb") as f:
        f.write(bytes(json.dumps(test_results), "utf8"))



# Visualize the test results for a model class.
def view_model_class(M, suffix=""):
    from util.data import Data
    # Declare a file for the data to be read quickly.
    data_file = f"data_{model_str(M)}{suffix}.pkl"
    if (not os.path.exists(data_file)):
        print("Reading data from directory..")
        # Load in any existing data that's already been collected.
        data_dir = os.path.join("data", f"testing_{model_str(M)}{suffix}")
        names = ['n', 'd', 'df', 'function', 'model', 'ds', 'ns',
                 'steps', 'test_size', 'func_seed', 'model_seed',
                 'fit_time', 'predict_time',
                 'initial_mean_squared_error', 'mean_squared_error',
                 'error_deciles', 'fit_record']
        configs = set(names[:9])
        trials = set(names[9:11])
        # Read the data.
        d = Data(names=names)
        for test_id in sorted(os.listdir(data_dir)):
            with gzip.open(os.path.join(data_dir, test_id), "rb") as f:
                test_results = json.loads(str(f.read(), "utf8"))
            for conf_str in sorted(test_results):
                d.append([test_results[conf_str][n] for n in names])
        # Stack the data across the axes that are repeated for each trial.
        to_stack = [n for n in d.names if n not in configs]
        print("  reorganizing data, stacking across trials..")
        d.stack(to_stack)
        print("  deleting columns that only have 1 unique value..")
        import pickle
        hash_func = lambda v: hash(pickle.dumps(v))
        num_unique = {n:len(set(map(hash_func, d[n]))) for n in d.names}
        for n in list(d.names):
            if (num_unique[n] == 1):
                print("  ", n, d[n][0])
                d.pop(n)
        print("  saving data for faster loading..")
        d.save(data_file)
    else:
        print("Loading data from file..")
        d = Data.load(data_file)
    
    print()
    print(d)
    print()

    # Get the names of the columns that are for test configuration.
    conf_cols = [n for n,t in zip(d.names, d.types) if t is not list][::-1]
    print("conf_cols: ", conf_cols)
    trial_cols = [n for n,t in zip(d.names, d.types) if t is list]
    print("trial_cols: ", trial_cols)
    print()
    d.sort(key=lambda row: tuple(row[conf_cols]))

    # Create plots for each config.
    from util.plot import Plot
    from util.stats.plotting import plot_percentiles
    for i,test in enumerate(d):
        config = test[conf_cols]
        conf_str = " ".join(list(map(str, config)))
        print("", conf_str, end="   ")
        # Create a visual of the model training error and test error.
        p = Plot(conf_str, "step", "mean squared error")
        # Plot the percentiles of the training error.
        records = np.asarray(test["fit_record"]).T
        step = np.arange(1, 1+records.shape[0])
        p = plot_percentiles(p, "MSE", step, records, color=1)
        # Plot the testing error distribution.
        def sideways_pdf(name, values, color=0):
            cdf = cdf_fit(values, smooth=True)
            pdf_x = cdf.inverse(np.linspace(0, 1, 1000))
            pdf_y = cdf.derivative(pdf_x)
            #    reorient the CDF to be on the far right
            pdf_x, pdf_y = pdf_y, pdf_x
            pdf_x -= min(pdf_x)
            pdf_x /= -max(pdf_x)
            pdf_x -= min(pdf_x)
            pdf_x *= (max(step) - min(step)) / 20
            pdf_x += max(step)
            pdf_x = pdf_x.tolist()
            pdf_y = pdf_y.tolist()
            pdf_x = [max(pdf_x)] + pdf_x + [max(pdf_x)]
            pdf_y = [pdf_y[0]] + pdf_y + [pdf_y[-1]]
            # Plot the distributions.
            p.add(f"{name}", pdf_x, pdf_y, mode="lines", color=color,
                  fill="toself", line_width=0, group="test cdf")
            pdf_x = pdf_x[1:-1]
            pdf_y = pdf_y[1:-1]
            p.add(f"{name} Line", pdf_x, pdf_y, mode="lines", color=(0,0,0,0.8),
                  line_width=0.5, show_in_legend=False, group="test cdf")
        sideways_pdf("Prefit Error", test["initial_mean_squared_error"], color=0)
        sideways_pdf("Test Error", test["mean_squared_error"], color=3)
        # Show the plot (only the first and the last).
        p.show(append=True, show=(i==0 or i==len(d)-1))


from util.approximate import PLRM

suffix = "_101-step-adaptation"

# # Run tests for a model class.
# test_model_class(PLRM, suffix=suffix)

# from notifications import send_email
# send_email("", "Done testing")

# View the test results for a model class.
view_model_class(PLRM, suffix=suffix)


# 2022-01-24 14:39:23
# 
################################################################
# # # Add distributions for the fit and predict times.         #
# # fit_cdf = cdf_fit(test["fit_time"])                        #
# # predict_cdf = cdf_fit(test["predict_time"])                #
# # p.add_func("F -- " +conf_str, fit_cdf, fit_cdf(),          #
# #            color=i, group="fits", show_in_legend=i==0)     #
# # p.add_func("P -- " + conf_str, predict_cdf, predict_cdf(), #
# #            color=i, group="predicts", show_in_legend=i==0) #
################################################################
 


# 2022-01-24 14:42:19
# 
        ########################################################
        # # Add distributions for the testing errors.          #
        # error_cdf = cdf_fit(test["mean_squared_error"])      #
        # p.add_func("E -- "+conf_str, error_cdf, error_cdf(), #
        #            color=i) #, group="errors")               #
        ########################################################


# 2022-01-24 16:51:51
# 
        ##############################################################################
        # n_smooth, deviations, stdev, dist_points = smooth                          #
        # n_smooth += (n_smooth+1) % 2 # Force n_smooth to be odd.                   #
        # # Create smoothing evaluation points (always an odd number).               #
        # smooth_points = np.linspace(-deviations*stdev, deviations*stdev, n_smooth) #
        # # Create weights for the smoothing points based on a normal distribution.  #
        # smooth_weights = np.exp(                                                   #
        #     -(smooth_points / stdev)**2 / 2                                        #
        # ) / (stdev * np.sqrt(2*np.pi))                                             #
        # smooth_weights /= smooth_weights.sum()                                     #
        # # Compute new x and y points for the smoothed distribution.                #
        # new_cdf_x = np.linspace(0, 1, dist_points)                                 #
        # cdf_y = np.asarray([                                                       #
        #     np.interp(x + smooth_points, cdf_x, cdf_y)                             #
        #     for x in new_cdf_x                                                     #
        # ])                                                                         #
        # cdf_y = np.dot(cdf_y, smooth_weights)                                      #
        # cdf_x = new_cdf_x                                                          #
        # cdf_y -= cdf_y.min()                                                       #
        # cdf_y /= cdf_y.max()                                                       #
        ##############################################################################


# 2022-01-25 22:14:38
# 
        ######################################################################
        # print("  identifying which columns are configuration columns..")   #
        # # num_unique = {n:len(set(map(hash_func, d[n]))) for n in d.names} #
        ######################################################################
