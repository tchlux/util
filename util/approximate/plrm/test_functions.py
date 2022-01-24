# Define a set of test functions for approximation.
import numpy as np
from itertools import combinations
from util.random import latin_sphere


# Return a random CDF.
def rand_cdf(nodes=3, power=1.0, smooth=(100, 5, 0.1, 1000)):
    # Generate random steps in the x and y direction (that sum to 1).
    cdf_x = np.linspace(0, 1, nodes+2)
    # Randomly assign a probability to each bin (this does not produce diverse CDFs).
    cdf_y = np.random.random(size=cdf_x.size) ** power
    cdf_y[0] = 0.0
    cdf_y /= cdf_y.sum()
    cdf_y = cdf_y.cumsum()
    # Smooth the CDF y values if desired.
    if smooth is not None:
        n_smooth, deviations, stdev, dist_points = smooth
        n_smooth += (n_smooth+1) % 2 # Force n_smooth to be odd.
        # Create smoothing evaluation points (always an odd number).
        smooth_points = np.linspace(-deviations*stdev, deviations*stdev, n_smooth)
        # Create weights for the smoothing points based on a normal distribution.
        smooth_weights = np.exp(
            -(smooth_points / stdev)**2 / 2
        ) / (stdev * np.sqrt(2*np.pi))
        smooth_weights /= smooth_weights.sum()
        # Compute new x and y points for the smoothed distribution.
        new_cdf_x = np.linspace(0, 1, dist_points)
        cdf_y = np.asarray([
            np.interp(x + smooth_points, cdf_x, cdf_y)
            for x in new_cdf_x
        ])
        cdf_y = np.dot(cdf_y, smooth_weights)
        cdf_x = new_cdf_x
        cdf_y -= cdf_y.min()
        cdf_y /= cdf_y.max()
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
    n = 10000
    d = 2
    df = 3
    df_intrinsic = df
    shift = 0.2
    distort = False
    skew = False
    ortho = True
    f = pure(d, df=df, df_intrinsic=df_intrinsic,
             shift=shift, distort=distort, skew=skew, ortho=ortho)
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
def sampled(n, d, dp=2, df=1, nf=2, nf_intrinsic=None,
            max_shift=0.5, power=3, return_info=False, seed=0):
    # By default, make all of the "modes" that can be observed maximally different.
    if (nf_intrinsic is None):
        nf_intrinsic = 2 * nf
    # Generate a "truth function" that creates orthogonal output
    #  components that can be combined to form "category predictions".
    f = pure(dp, df=nf, df_intrinsic=nf_intrinsic, seed=seed)
    # Generate positions for all of the dimensions in some space.
    positions = latin_sphere(d, dp, inside=True)
    # Assign random center offsets to all of the points (the y values).
    offsets = latin_sphere(n, dp, inside=True) * max_shift
    # Generate a set of evaluation points (at all positions, for all offsets).
    x = np.concatenate([(positions - o) for o in offsets])
    # Evaluate the function at all positions, this will produce
    #   nf channels of output for every (position,offset) pair.
    new_x, y = f(x)
    # Generate a set of weights, this will be used to blend the modes
    #   of the output from the function into single "observations".
    #   Ideally, the "most observable" output is the predicted one.
    weights = (latin_sphere(n, nf) + 1) ** power
    weights = (weights.T / np.sum(weights, axis=1)).T
    # Initialize the inputs, then collapse the observations using the weights.
    inputs = np.ones((n, d))
    for i in range(n):
        inputs[i,:] = np.dot(y[i*d:(i+1)*d], weights[i])
    # The outputs are the offset terms and weights concatenated,
    #  both are something that should be predicted by an approximator.
    outputs = np.concatenate((offsets, weights), axis=1)
    # If desired, return the true positions.
    if (return_info):
        return inputs, outputs, positions
    # Return the pairs of points.
    return inputs, outputs


# Test the function that mimics the sampling behavior.
def _test_sampled():
    from util.plot import Plot
    n = 10
    d = 3000
    x, y, positions = sampled(n, d, return_info=True)
    p = Plot()
    for i in range(n):
        p.add(f"post {i+1}", [y[i,0],y[i,0]], [y[i,1], y[i,1]], [-2, 2],
              mode="lines", line_width=4, color=1, group=i, show_in_legend=False)
        p.add(f"Offset {i+1} ({y[i,0]:.3f}, {y[i,1]:.3f}) [{y[i,2:].round(3)}]",
              *positions.T, x[i], use_gradient=True,
              marker_size=4, marker_line_width=1, group=i)
    p.show(z_range=[-3,5])


_test_rand_cdf()
_test_cos_norm()
_test_pure()
_test_sampled()


# 2022-01-22 15:23:20
# 
####################################################################
# # Generate unique subsets of input components to apply the pure  #
# #   funciton to, that will create partially independent outputs. #
# slice_generators = []                                            #
# for num_dims in range(d-1, 0, -1):                               #
#     slice_generators.append( combinations(range(d), num_dims) )  #
# i = -1                                                           #
# while (len(y) < df) and (len(slice_generators) > 0):             #
#     i = (i + 1) % len(slice_generators)                          #
#     # Get the next set of dimensions.                            #
#     try:                                                         #
#         dims = next(slice_generators[i])                         #
#         y.append( f(x[:,dims]) )                                 #
#     except StopIteration:                                        #
#         slice_generators.pop(i)                                  #
#         i -= 1                                                   #
# # If y needs to be in higher dimension, simply raise the         #
# #   dimension of the existing output for the function.           #
# if (y.shape[1] < df):                                            #
#     projection = latin_sphere(df, y.shape[1]).T                  #
#     print("projection.shape: ", projection.shape)                #
#     y = np.matmul(y, projection)                                 #
####################################################################
