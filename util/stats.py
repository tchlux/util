import os
from util.math import abs_diff, SMALL
from util.random import pairs

# Import "numpy" by modifying the system path to remove conflicts.
import sys
_ = []
for i in range(len(sys.path)-1,-1,-1):
    if "util/util" in sys.path[i]: _ = [sys.path.pop(i)] + _
import numpy as np
sys.path = _ + sys.path


# =============================================
#      Metric Principle Component Analysis     
# =============================================

# Generate vector between scaled by metric difference. Give the metric
# the indices of "vectors" in the provided matrix.
def gen_random_metric_diff(matrix, index_metric, power=2, count=None):
    # Iterate over random pairs.
    for (p1, p2) in pairs(len(matrix), count):
        vec = matrix[p1] - matrix[p2]
        length = np.sqrt(np.sum(vec**2))
        if length > 0: vec /= length**power
        yield vec * index_metric(p1, p2)

# Given a set of row-vectors, compute a convex weighting that is
# proportional to the inverse total variation of metric distance
# between adjacent points along each vector (component).
def normalize_error(points, values, metric, display=True):
    # Compute the magnitudes using total variation.
    if display: print(" estimating error slope per axis.. ", end="\r", flush=True)
    avg_slope = np.zeros(points.shape[1])
    update = end = ""
    for axis in range(points.shape[1]):
        update = "\b"*len(end)
        end = f"{axis+1} of {points.shape[1]}"
        if display: print(update, end=end, flush=True)
        # Sort points according to this axis.
        ordering = np.argsort(points[:,axis])
        for i in range(len(ordering)-1):
            p1, p2 = ordering[i], ordering[i+1]
            diff = (points[p2,axis] - points[p1,axis])
            if (diff > 0): avg_slope[axis] += metric(values[p1], values[p2]) / (diff * points.shape[1])
    if display: print("\r                                   ", end="\r", flush=True)
    # If there are dimensions with no slope, then they are the only ones needed.
    if (min(avg_slope) <= 0.0): avg_slope = np.where(avg_slope == 0, 1., float('inf'))
    # We want to minimize the expected error associated with the 2-norm.
    # 
    # E[f(x) - f(y)] = sum( E[f(x)1 - f(y)1]^2 + ... + E[f(x)d - f(y)d]^2 )^(1/2)
    #               -> { E[f(x)1-f(y)1]^(-2), ..., E[f(x)d-f(y)d]^(-2) }
    # 
    # This is the normalizing vector we want, scaled to unit determinant.
    error_fix = avg_slope**(-2)
    error_fix /= np.median(error_fix) # Make product of all values stable.
    return error_fix / np.prod(error_fix)**(1/points.shape[1])
    # ^^ This line was numerically unstable (for small 

# Compute the metric PCA (pca of the between vectors scaled by 
# metric difference slope).
def mpca(points, values, metric=abs_diff, num_components=None,
         num_vecs=None, display=True):
    if display: print(" normalizing axes by expected error..", end="\r", flush=True)
    # Set default values for num_components and num_vecs.
    if type(num_components) == type(None): num_components = min(points.shape)
    if (type(num_vecs) == type(None)): num_vecs = (len(points)**2-len(points))//2
    num_components = min(num_components, *points.shape)
    num_vecs = min(num_vecs, (len(points)**2-len(points))//2)
    if display: print(" allocating memory for metric between vectors..", end="\r", flush=True)
    m_vecs = np.zeros((num_vecs, points.shape[1]))
    if display: print(" generating metric vectors..", end="\r", flush=True)
    # Function that takes two indices and returns metric difference.
    index_metric = lambda i1, i2: metric(values[i1], values[i2])
    # Generator that produces "between vectors".
    vec_gen = gen_random_metric_diff(points, index_metric, count=num_vecs)
    for i,vec in enumerate(vec_gen): m_vecs[i,:] = vec
    # Compute the principle components of the between vectors.
    if display: print(" computing principle components..", end="\r", flush=True)
    components, _ = pca(m_vecs, num_components=num_components)
    if display: print(" normalizing components by metric slope..", end="\r", flush=True)
    weights = normalize_error(np.matmul(points, components.T), values, metric, display)
    if display: print("                                               ", end="\r", flush=True)
    # Return the principle components of the metric slope vectors.
    return components, weights

# Compute the principle components using sklearn.
def pca(points, num_components=None, display=True):
    from util.math import is_none
    from sklearn.decomposition import PCA        
    pca = PCA(n_components=num_components)
    if is_none(num_components): num_components = min(*points.shape)
    else:       num_components = min(num_components, *points.shape)
    if display: print(f"Computing {num_components} principle components..",end="\r", flush=True)
    pca.fit(points)
    if display: print( "                                                          ",end="\r", flush=True)
    principle_components = pca.components_
    magnitudes = pca.singular_values_
    # Normalize the component magnitudes to have sum 1.
    magnitudes /= np.sum(magnitudes)
    return principle_components, magnitudes

# ===================================
#      CDF and PDF Fit Functions     
# ===================================

# A class for handling standard operators over distributions. Supports:
#    addition with another distribution
#    subtraction from another distribution (KS-statistic)
#    multiplication by a constant in [0,1]
class Distribution():
    def __init__(self, func):
        self._function = func
        standard_attributes = ["min", "max", "derivative", 
                               "integral", "inverse", "nodes"]
        # Copy over the attributes of the distribution function.
        for attr in standard_attributes:
            if hasattr(func, attr): setattr(self, attr, getattr(func, attr))

    # Define addition with another distribution.
    def __radd__(self, func):
        # Special case, when distributions are used in a "sum" then an
        # integer will be added as the "start" of the sum. Ignore this.
        if (type(func) == int): return self
        # Get the new minimum and maximum of this distribution and "func".
        new_min = min(self.min, func.min)
        new_max = max(self.max, func.max)
        def new_func(x=None, func=func):
            if (type(x) == type(None)): return (new_min, new_max)
            else:                       return self._function(x) + func(x)
        new_func.min = new_min
        new_func.max = new_max
        # Copy over the attributes of the distribution function.
        standard_attributes = ["derivative", "integral"]
        for attr in standard_attributes:
            if (hasattr(self, attr) and hasattr(func, attr)):
                mine = getattr(self, attr)
                other = getattr(func, attr)
                setattr(new_func, attr, lambda x: mine(x) + other(x))
        # Return a distribution object.
        return Distribution(new_func)
    __add__ = __radd__

    # Define multiplication by a number (usually a float).
    def __rmul__(self, number):
        def new_func(x=None, m=number, f=self._function):
            if (type(x) == type(None)): return f(None)
            else:                       return m * f(x)
        new_dist = Distribution(self)
        new_dist._function = new_func
        # Copy over the attributes of the distribution function.
        standard_attributes = ["derivative", "integral"]
        for attr in standard_attributes:
            if hasattr(self, attr):
                f = getattr(self, attr)
                setattr(new_dist, attr, lambda x: number * f(x))
        return new_dist
    __mul__ = __rmul__

    # Define subtraction from another distribution. (KS statistic)
    def __rsub__(self, func): return ks_diff(self, func)
    __sub__ = __rsub__
        
    # Define how to "call" a distribution.
    def __call__(self, *args, **kwargs):
        return self._function(*args, *kwargs)

# Construct a version of a function that has been smoothed with a gaussian.
def gauss_smooth(func, stdev=1., n=100):
    # We must incorporate some smoothing.
    from scipy.stats import norm
    # Construct a set of gaussian weights (to apply nearby on the function).
    eval_pts = np.linspace(-5 * stdev, 5 * stdev, n)
    weights = norm(0, stdev).pdf(eval_pts)
    weights /= sum(weights)
    # Generate a smoothed function (with the [min,max] no-arg convenience).
    def smooth_func(x=None):
        if (type(x) == type(None)): return func()
        else:
            try:
                wts = (np.ones((len(x), n)) * weights)
                pts = (x + (eval_pts * np.ones((len(x), n))).T).T
                return np.sum(wts * func(pts), axis=1)
            except:
                return sum(weights * func(eval_pts + x))
    # Construct the smooth derivative as well.
    eps = .01 * (func.max - func.min)
    def deriv(x): return (smooth_func(x + eps) - smooth_func(x - eps))
    smooth_func.derivative = deriv
    smooth_func.min = func.min
    smooth_func.max = func.max
    # Return the smoothed function.
    return Distribution(smooth_func)

# Given a list of numbers, generate two lists, one of the CDF x points
# and one of the CDF y points (for fitting).
def cdf_points(data):
    from util.math import SMALL
    # Sort the data (without changing it) and get the min and max
    data = np.array(sorted(data))
    min_pt = data[0]
    max_pt = data[-1]
    # Initialize a holder for all CDF points
    data_vals = [0]
    # Add all of the CDF points for the data set
    for i,val in enumerate(data):
        if ((i > 0) and (val == data[i-1])): data_vals[-1] = (i+1)/len(data)
        else:                                data_vals.append( (i+1)/len(data) )
    # Add the 100 percentile point if it was not added already
    if (data_vals[-1] != 1.0): data_vals[-1] = 1.0
    # Convert data into its unique values.
    data = np.array([data[0] - SMALL] + sorted(set(data)))
    # Convert it all to numpy format for ease of processing
    data_vals = np.array(data_vals)
    # Return the two lists that define the CDF points.
    return data, data_vals


# Returns the CDF function value for any given x-value
#   fit -- "linear" linearly interpolates between Emperical Dist points.
#          "cubic" constructs a monotonic piecewise cubic interpolating spline
#          <other> returns the raw Emperical Distribution Function.
def cdf_fit(data, fit="linear", smooth=None):
    fit_type = fit
    # Scipy functions for spline fits.
    from scipy.interpolate import splrep, splev
    from scipy.interpolate import PchipInterpolator
    # Check for function usage. Either the user gave (x, y) points of
    # a CDF that they want to fit or they gave a list of numbers.
    if (type(data) == tuple) and (len(data) == 2):
        data, data_vals = data
    else:
        data, data_vals = cdf_points(data)
    # Retrieve the minimum and maximum points.
    min_pt, max_pt = data[0], data[-1]
    # Generate a fit for the data points
    if (fit_type == "linear"):
        # Generate linear function
        fit_params = splrep(data, data_vals, k=1)
        fit = lambda x_val: splev(x_val, fit_params, ext=3)
        fit.derivative = lambda x_val: splev(x_val, fit_params, der=1, ext=3)
        # Generate inverse linear function
        inv_fit_params = splrep(data_vals, data, k=1)
        inv_fit = lambda y_val: splev(y_val, inv_fit_params, ext=3)
        inv_fit.derivative = lambda y_val: splev(y_val, inv_fit_params, der=1, ext=3)
        fit.inverse = inv_fit
    elif (fit_type == "cubic"):
        # Generate piecewise cubic monotone increasing spline
        spline_fit = PchipInterpolator(data, data_vals)
        # Define a fit function that extrapolates properly.
        def fit(x):
            try:
                # If "x" is an array, then use array ops to fix boundaries.
                len(x)
                out = spline_fit(x)
                out[x < min_pt] = 0
                out[x > max_pt] = 1
            except:
                # If "x" is not an array, then fix specific test point.
                if   (x < min_pt): out = 0.
                elif (x > max_pt): out = 1.
                else:              out = float(spline_fit(x))
            return out
        # Define and save an inverse function.
        from util.optimize import zero
        def inverse(x):
            # Define a function that computes the inverse for a single value.
            single_inverse = lambda x: zero(lambda v: fit(v) - x, min_pt, max_pt, max_steps=50)
            # Record vector result if vector was given.
            if hasattr(x, "__len__"):
                result = []
                for xi in x: result.append( single_inverse(xi) )
                result = np.array(result)
            else: result = single_inverse(x)
            # Return computed inverse.
            return result
        fit.inverse = inverse
        # Store a derivative (that is a PDF).
        def derivative(x=None, deriv=spline_fit.derivative(1)):
            if type(x) == type(None): return (min_pt, max_pt)
            try:
                # If "x" is an array, then use array ops to fix boundaries.
                len(x)
                out = deriv(x)
                out[x <= min_pt] = 0
                out[x >= max_pt] = 0
            except:
                # If "x" is not an array, then fix specific test point.
                if (x <= min_pt) or (x >= max_pt): out = 0.
                else:                              out = float(deriv(x))
            return out
        # Store a second derivative.
        def second_derivative(x=None, deriv=spline_fit.derivative(2)):
            return derivative(x, deriv)
        # Set details for the derivative function.
        derivative.min = min_pt
        derivative.max = min_pt
        derivative.integral = fit
        derivative.derivative = second_derivative
        # Set the derivative of the fit.
        fit.derivative = derivative
    else:
        # Construct the empirical distribution fit.
        def fit(x_val):
            try:
                # Handle the vector valued case
                len(x_val)
                vals = []
                for x in x_val:
                    index = np.searchsorted(data, x, side="right") - 1
                    index = max(index, 0)
                    vals.append(data_vals[index])
                return np.array(vals)
            except TypeError:
                index = np.searchsorted(data, x, side="right") - 1
                index = max(index, 0)
                # Handle the single value case
                return data_vals[index]
        # Construct an inverse function for the EDF.
        def inverse(perc, data=data, data_vals=data_vals):
            try:
                # Handle the vector valued case
                len(perc)
                indices = []
                for p in perc:
                    l_index = np.searchsorted(data_vals, max(0,min(p,1)), side="left")
                    r_index = np.searchsorted(data_vals, max(0,min(p,1)), side="right")
                    index = (l_index + r_index) // 2
                    indices.append( data[min(index, len(data)-1)] )
                return np.array(indices)
            except TypeError:
                # Handle the single value case
                l_index = np.searchsorted(data_vals, max(0,min(perc,1)), side="left")
                r_index = np.searchsorted(data_vals, max(0,min(perc,1)), side="right")
                index = (l_index + r_index) // 2
                return data[min(index, len(data)-1)]
        # Construct a derivative function for the EDF.
        def derivative(x_val):
            try:
                # Handle the vector valued case
                len(x_val)
                vals = []
                for x in x_val:
                    vals.append( float(x in data) )
                return np.array(vals)
            except TypeError:
                # Handle the single value case
                return float(x_val in data)
        # Store the inverse and derivative as attributes.
        fit.inverse = inverse
        fit.derivative = derivative
    
    # Generate a function that computes this CDF's points
    def cdf_func(x_val=None):
        # Make the function return a min and max if nothing is provided.
        if type(x_val) == type(None):
            return (min_pt, max_pt)
        else:
            # Treat this like a normal spline interpolation.
            return fit(x_val)

    # Set the min and max of the data.
    cdf_func.min = min_pt
    cdf_func.max = max_pt
    # Set the derivative function
    cdf_func.derivative = fit.derivative
    # Set the inverse function
    cdf_func.inverse = fit.inverse
    # Smooth the CDF if requested by user.
    if type(smooth) != type(None):
        # Smooth the pdf fit if that was requested.
        width = (cdf_func.max - cdf_func.min)
        stdev = smooth * width
        cdf_func = gauss_smooth(cdf_func, stdev)
    # Store the original nodes.
    cdf_func.nodes = list(zip(data, data_vals))
    # Return the custom function for this set of points
    return Distribution(cdf_func)

# Return a PDF fit function that is smoothed with a gaussian kernel.
# "smooth" is the percentage of the data width to use as the standard
# deviation of the gaussian kernel. "n" is the number of samples to
# make when doing the gaussian kernel smoothing.
def pdf_fit(data, smooth=.05, n=1000, **cdf_fit_kwargs):
    cdf_func = cdf_fit(data, **cdf_fit_kwargs)
    if smooth:
        # Smooth the pdf fit if that was requested.
        width = (cdf_func.max - cdf_func.min)
        stdev = smooth * width
        cdf_func = gauss_smooth(cdf_func, stdev, n)
    # Return the PDF as the derivative of the CDF.
    pdf_func = cdf_func.derivative
    # Take the first derivative of the CDF function to get the PDF
    def pdf(x=None):
        # Return the min and max when nothing is provided.
        if (type(x) == type(None)): return (cdf_func.min, cdf_func.max)
        try:
            # Handle the vector valued case.
            len(x)                
            # Return array of results.
            return np.array([pdf_func(v) for v in x])
        except TypeError:
            # Return the PDF function evaluation.
            return pdf_func(x)
        return pdf_func(x)
    # Set the "min" and "max" attributes of this function.
    pdf.min = cdf_func.min
    pdf.max = cdf_func.max
    pdf.integral = cdf_func
    return Distribution(pdf)

# ===============================================
#      Statistical Measurements Involving KS     
# ===============================================

# Calculate the maximum difference between two CDF functions (two sample).
def ks_diff(test_func, true_func, method=100):
    # Cycle through the functions to find the min and max of all ranges
    min_test, max_test = test_func()
    min_true, max_true = true_func()
    min_val, max_val = max(min_test, min_true), min(max_test, max_true)
    if method in {"scipy", "util"}:
        diff_func = lambda x: -abs(test_func(x) - true_func(x))
        if method == "scipy":
            from scipy.optimize import minimize
            # METHOD 1:
            #  Use scipy to maximize the difference function between
            #  the two cdfs in order to find the greatest difference.
            sol = minimize(diff_func, [(max_val - min_val) / 2],
                           bounds=[(min_val,max_val)], method='L-BFGS-B').x
        elif method == "util":
            from util.optimize import minimize
            # METHOD 2 (default):
            #  Use the default minimizer in "optimize" to maximize the
            #  difference function between the two cdfs in order to
            #  find the greatest difference.
            sol = minimize(diff_func, [(max_val - min_val) / 2],
                           bounds=[(min_val,max_val)])[0]
        greatest_diff = abs(test_func(sol) - true_func(sol))
    else:
        # METHOD 3:
        #  Generate a large set of x-points and find the difference
        #  between the functions at all of those points. Generate a
        #  grid of points and "zoom in" around the greatest difference
        #  points to identify the spot with largest gap.
        n = 100
        if (type(method) == int): n = method
        width = (max_val - min_val)
        x_points = np.linspace(min_val, max_val, n)
        diff = abs(test_func(x_points) - true_func(x_points))
        # Cycle zooming in about the greatest difference.
        greatest_diff = -float('inf')
        while (diff[np.argmax(diff)] > greatest_diff):
            lower = max(np.argmax(diff) - 1, 0)
            upper = min(np.argmax(diff) + 1, n-1)
            min_pt = max(x_points[lower], min_val)
            max_pt = min(x_points[upper], max_val)
            x_points = np.linspace(min_pt, max_pt, n)
            diff = abs(test_func(x_points) - true_func(x_points))
            greatest_diff = max(max(diff), greatest_diff)
    return greatest_diff

# Given a ks-statistic and the sample sizes of the two distributions
# compared, return the largest confidence with which the two
# distributions can be said to be the same.
def ks_p_value(ks_stat, n1, n2=float('inf')):
    # By definition of the KS-test: 
    # 
    #    KS > c(a) (1/n1 + 1/n2)^(1/2)
    # 
    #   where KS is the KS statistic, n1 and n2 are the respective
    #   sample sizes for each distribution and c(a) is defined as
    # 
    #    c(a) = ( -ln(a/2)/2 )^(1/2)
    # 
    #   is the standard for testing the probability with which two
    #   distributions come from different underlying distributions. If
    #   we want the distributions to be the same, we want the KS test
    #   to only pass with large values for "a" (large 'p-value'). The
    #   above check can be reversed to compute "a", which provides the
    #   largest p-value for which the KS test states the two
    #   distributions are not certainly different.
    # 
    #    c^-1(b) = 2 e^(-2 b^2)
    # 
    #    a = c^-1 ( (KS / (1/n1 + 1/n2)^(1/2)) )
    #      = 2 e^( -2 (KS / (1/n1 + 1/n2)^(1/2))^2 )
    #      = 2 e^( -2 KS^2 / |1/n1 + 1/n2| )
    #    
    #    and "a" cannot be larger than 1. Therefore, finally we have
    # 
    #    a = min(1, 2 e^( -2 KS^2 / |1/n1 + 1/n2| ))
    #                                                                 
    return min(1.0, 2 * np.exp( -2 * ( ks_stat**2 / abs((1/n1) + (1/n2)) )))

# Return the "confidence" that the provided distribution is normal.
def normal_confidence(distribution):
    from scipy.stats import kstest
    # # Make the distribution 0 mean and unit variance (unit standard deviation)
    # new_distribution = (distribution - np.mean(distribution)) / np.var(distribution)
    # Compare the distribution with a normal distribution
    new_distribution = distribution
    ks_statistic = kstest(new_distribution, "norm").statistic
    return ks_p_value(ks_statistic, len(distribution))


# =====================================================
#      Approximating the desired number of samples     
# =====================================================


# Compute the confidence of "num_samples" samples having less than or
# equal to "max_error" at exactly 1/2. This is achieved by factoring
# the expression (sum_{i=j}^k choose(n,i) / 2^n) to reduce operations.
def _half_confidence(num_samples, max_error):
    from util.math import Fraction, choose
    # Calculate the maximum number of steps away from the mid-point
    # that we can take *before* violating the error condition.
    odd = num_samples % 2
    start = (num_samples+1) // 2
    #         initial_error + steps * step_error <= max_error
    #   (max_error - initial_error) / step_error >= steps
    #  (max_error - initial_error) * num_samples >= steps
    # knowing that                    step_error  = 1 / num_samples
    min_error = odd * Fraction(1, (2*num_samples))    
    steps = int((max_error - min_error) * num_samples)
    # Put steps into allowable bounds.
    steps = max(0, min(num_samples - start, steps))
    # Handle two cases where there is no functioning confidence bound.
    if (odd) and (max_error < min_error): return Fraction()
    # Compute the fraction.
    numerator   = 1
    denominator = 2**(num_samples)
    # First compute the constant multiple outside of the sum.
    for i in range(start+steps+1, num_samples+1): numerator   *= i
    for i in range(2,             start+1):       denominator *= i
    # Compute the sum of the inner parts of the distribution.
    total = 0
    for i in range(start, start+steps+1):
        v = 1
        for j in range(i+1,             start+steps+1): v *= j
        for j in range(1+num_samples-i, start+1):       v *= j
        if (v > 1): total += v * (2 if ((i != start) or odd) else 1)
    # Perform the final numerator update.
    if (total > 0): numerator *= total
    # Return the (reduced) fraction.
    return Fraction(numerator, denominator)


# Given a list of numbers, return True if the given values provide an
# estimate to the underlying distribution with confidence bounded error.
def samples(size=None, error=None, confidence=None, at=None):
    from util.system import Timer
    t = Timer()

    # Determine what to calculate based on what was provided.
    from util.math import is_none, choose, Fraction
    if   is_none(size):       to_calculate = "samples"
    elif is_none(error):      to_calculate = "error"
    elif is_none(confidence): to_calculate = "confidence"
    else:                     to_calculate = "verify"
    # Default evaluation point is at (1/2), where the error is greatest.
    if is_none(at): at = Fraction(1, 2)
    else:           at = Fraction(at)
    # Set the default values for other things that were not provided.
    if type(error) == type(None):      error      = Fraction(10, 100)
    if type(confidence) == type(None): confidence = Fraction(95, 100)
    # Convert error and confidence to fraction types if necessary.
    if not type(error)      == Fraction: error      = Fraction(error)
    if not type(confidence) == Fraction: confidence = Fraction(confidence)
    # If the user provided something with a length, use that number.
    if hasattr(size, "__len__"): size = len(size)
    # \sum_{i=0}^n choose(n, i) * ( at^i (1-at)^(n-i) )
    if not is_none(size):
        # Compute the probability of any given observed EDF value.
        prob = lambda i: choose(size, i) * (at**i * (1-at)**(size-i))
        # If we are calculating the confidence or verifying, compute confidence.
        if to_calculate in {"confidence", "verify"}:
            if (at == 1/2): conf = _half_confidence(size, error)
            else:
                conf = Fraction()
                steps = 0
                # Sum those probabilities that are closer than "error" distance.
                for i in range(size+1):
                    p = Fraction(i, size)
                    if (abs(p - at) <= error):
                        steps += 1
                        conf += prob(i)
            # Return the total confidence.
            if to_calculate == "confidence": return float(conf)
            else:                            return conf >= confidence
        elif to_calculate == "error":
            # Store the "contained" outcomes by "allowed error".
            error = Fraction()
            contained = Fraction()
            # Sort the percentiles by their distance from "at".
            i_p = sorted(enumerate(Fraction(i,size,_normalize=False)
                                   for i in range(size+1)),
                         key=lambda ip: abs(ip[1]-at))
            # Cycle through percentiles, starting closest to "at" and moving out.
            for step in range(len(i_p)):
                # If this step has the same probability as the last, skip.
                if (i_p[step][1] == i_p[step-1][1]): continue
                i, p = i_p[step]
                # Compute the amount of data contained by this step away.
                next_contained = contained + prob(i)
                # If the distance from "at" is the same for two steps, take two.
                if (step+1 < len(i_p)) and (abs(at-i_p[step][1]) == abs(at-i_p[step+1][1])):
                    next_contained += prob( i_p[step+1][0] )
                # Only update the "allowed error" if confidence is maintained.
                if next_contained < confidence:
                    contained = next_contained
                    error = abs(i_p[step][1] - at)
                else: break
            return float(error)
    else:
        # Compute the number of samples required.
        size, step = 2**10, 2**9
        # print("Desired ----------------")
        # print("error:      ",error)
        # print("confidence: ",confidence)
        # for size in range(2, 500):
        #     conf_below = samples(size-1, error=error, at=at)
        #     conf_at = samples(size, error=error, at=at)
        #     print("", "size: ",size, float(f"{conf_below:.2e}"), float(f"{conf_at:.2e}"))
        # exit()

        under, over = 0, None
        # We have the right size when any smaller size is not passable.
        conf_below = samples(size=size-1, error=error, at=at)
        conf_at = samples(size=size, error=error, at=at)
        print("", "size: ",size, float(f"{conf_below:.2e}"), float(f"{conf_at:.2e}"))
        while not (conf_below < confidence <= conf_at):
            if conf_at < confidence:
                # Update "under". Scale up if we haven't found "over".
                under = max(under, size)
                if (over == None): step *= 2
                # Take the step.
                size += step
            else:
                # Update "over". Take step. Scale down.
                over = min(over if (over != None) else float('inf'), size)
                size = size - step
                step = step // 2
            # Recompute the confidence at and below this step size.
            conf_below = samples(size-1, error=error, at=at)
            conf_at = samples(size, error=error, at=at)
            print("", "size: ",size, float(f"{conf_below:.2e}"), float(f"{conf_at:.2e}"))
            # Correct for strange sample size error that can happen bc
            # of alignment of "at" and one of the sample values.
            if conf_at < conf_below:
                size = size-1
                conf_at = conf_below
                conf_below = samples(size-1, error=error, at=at)
                print("", "size: ",size, float(f"{conf_below:.2e}"), float(f"{conf_at:.2e}"))
        # Return the computed best sample size.
        return size



# ====================================================
#      Measuring the Difference Between Sequences     
# ====================================================

# Generator for returning pairs of EDF points as (val1, val2, density).
def edf_pair_gen(seq):
    # Compute the 'extra' slope lost at the beginning of the function.
    for i in range(len(seq)):
        if seq[i] != seq[0]: break
    extra = (i/len(seq))
    shift = lambda v: (v - extra) / (1 - extra)
    # Manually set the first slope value.
    x1, y1 = -float('inf'), 0
    for i in range(len(seq)):
        # Cycle until we hit the last occurrence of the next value.
        if ((i+1) < len(seq)) and (seq[i+1] == seq[i]): continue
        x2, y2 = seq[i], shift( (i+1)/len(seq) )
        yield (x1, x2, y2 - y1)
        # Cycle values.
        x1, y1 = x2, y2
    # Manually yield the last element in the sequence.
    yield (x1, float('inf'), 0)

# Use a linear interpolation of the EDF points to compute the integral
# of the absolute difference between the empirical PDF's of two sequences.
def epdf_diff(seq_1, seq_2):
    # Construct generators for the EDF point pairs.
    gen_1 = edf_pair_gen(seq_1)
    gen_2 = edf_pair_gen(seq_2)
    # Get the initial values from the generators.
    low_1, upp_1, density_1 = gen_1.__next__()
    low_2, upp_2, density_2 = gen_2.__next__()
    shared = 0.
    # Cycle until completing the integral difference between pEDFs.
    while (upp_1 < float('inf')) or (upp_2 < float('inf')):
        # Compute the width of the interval and the percentage of data
        # within the width for each of the distributions.
        width = min(upp_1, upp_2) - max(low_1, low_2)
        if width < float('inf'):
            # Convert the EDF points into densities over the interval.
            sub_density_1 = density_1 * (width / (upp_1-low_1))
            sub_density_2 = density_2 * (width / (upp_2-low_2))
            if abs(sub_density_1 - sub_density_2) >= 0:
                shared += min(sub_density_1, sub_density_2)
        # Cycle whichever range has the lower maximum.
        if (upp_1 > upp_2):
            low_2, upp_2, density_2 = gen_2.__next__()
        elif (upp_2 > upp_1):
            low_1, upp_1, density_1 = gen_1.__next__()
        else:
            low_1, upp_1, density_1 = gen_1.__next__()
            low_2, upp_2, density_2 = gen_2.__next__()
    return max(0, 1 - shared)

# Compute the PDF of a sequence of categories.
def categorical_pdf(seq, unique=None):
    from util.system import hash, sorted_unique
    if (type(unique) == type(None)): unique = sorted_unique(seq)
    try:              return {c:sum(v==c for v in seq) / len(seq) for c in unique}
    except TypeError: return {hash(c): sum(v==c for v in seq) / len(seq) for c in unique}

# Given the PDF (in a dictionary object) for two different sets of
# categorical values, compute the total change in the PDF sum(abs(diffs)).
def categorical_diff(pdf_1, pdf_2):
    # Check for equally lengthed sequences.
    class CannotCompare(Exception): pass
    def near_one(val): return (abs(val - 1) < 1e-5)
    if not near_one(sum(pdf_1.values())):
        raise(CannotCompare(f"Sum of PDF values for sequence one '{sum(pdf_1.values()):.2f}' is not 1."))
    if not near_one(sum(pdf_2.values())):
        raise(CannotCompare(f"Sum of PDF values for sequence one '{sum(pdf_2.values()):.2f}' is not 1."))
    # Sum the differences in categories.
    all_values = sorted(set(pdf_1).union(set(pdf_2)))
    return sum(abs(pdf_1.get(v,0.) - pdf_2.get(v,0.)) for v in all_values) / 2

# Compute the effect between two sequences. Use the following:
#    number vs. number     -- Correlation coefficient between the two sequences.
#    category vs. number   -- "method" 1-norm difference between full distribution
#                             and conditional distributions for categories.
#    category vs. category -- "method" total difference between full distribution
#                             of other sequence given value of one sequence.
def effect(seq_1, seq_2, method="mean"):
    from util.system import hash, sorted_unique
    from util.math import is_numeric
    # Check for equally lengthed sequences.
    class CannotCompare(Exception): pass
    if len(seq_1) != len(seq_2): raise(CannotCompare("Provided sequences must have same length for comparison to be made."))
    # Check for index accessibility.
    try:              seq_1[0], seq_2[0]
    except TypeError: raise(CannotCompare("Provided sequences must both be index-accessible."))
    # Get the set of attributes a numeric type should have.
    # The three different outcomes of comparison could be:
    num_num = (is_numeric(seq_1[0]) and is_numeric(seq_2[0]))
    cat_cat = (not is_numeric(seq_1[0])) and (not is_numeric(seq_2[0]))
    num_cat = (not cat_cat) and (not num_num)
    # Execute the appropriate comparison according to the types.
    if num_num:
        # Compute the correlation (highly correlated have no difference).
        return float(np.corrcoef(seq_1, seq_2)[0,1])
    elif num_cat:
        # Compute the change in numeric distribution caused by a chosen category.
        nums, cats = (seq_1, seq_2) if is_numeric(seq_1[0]) else (seq_2, seq_1)
        # Compute the sorted order for "nums" by index.
        sort = sorted(range(len(nums)), key=lambda i: nums[i])
        nums, cats = [nums[i] for i in sort], [cats[i] for i in sort]
        # Consider each unique categorical value, compute epdf diff.
        diffs = {}
        for cat in sorted_unique(cats):
            main_nums = [n for (n,c) in zip(nums,cats) if c != cat]
            sub_nums = [n for (n,c) in zip(nums,cats) if c == cat]
            dist_diff = epdf_diff(main_nums, sub_nums)
            try:              diffs[cat]       = dist_diff
            except TypeError: diffs[hash(cat)] = dist_diff
        if   (method == "max"):  return max(diffs.values())
        elif (method == "min"):  return min(diffs.values())
        elif (method == "mean"): return sum(diffs.values()) / len(diffs)
        else:                    return diffs
    elif cat_cat:
        diffs = {}
        # Generate unique lists of values (hash those that aren't).
        unique_1 = sorted_unique(seq_1)
        unique_2 = sorted_unique(seq_2)
        # Compute the percentage of each unique category for full sequences.
        full_1 = categorical_pdf(seq_1, unique_1)
        full_2 = categorical_pdf(seq_2, unique_2)
        # Identify the difference when selecting sequence 1 values.
        for cat in unique_1:
            main_seq = categorical_pdf([c2 for (c1,c2) in zip(seq_1,seq_2) if c1 != cat])
            sub_seq = categorical_pdf([c2 for (c1,c2) in zip(seq_1,seq_2) if c1 == cat])
            if (len(main_seq) > 0) and (len(sub_seq) > 0):
                dist_diff = categorical_diff( main_seq, sub_seq )
            else: dist_diff = 0.
            try:              diffs[(cat,None)]       = dist_diff
            except TypeError: diffs[(hash(cat),None)] = dist_diff
        # Identify the difference when selecting sequence 2 values.
        for cat in unique_2:
            main_seq = categorical_pdf([c1 for (c1,c2) in zip(seq_1,seq_2) if c2 != cat])
            sub_seq = categorical_pdf([c1 for (c1,c2) in zip(seq_1,seq_2) if c2 == cat])
            dist_diff = categorical_diff( main_seq, sub_seq )
            try:              diffs[(None,cat)]       = dist_diff
            except TypeError: diffs[(None,hash(cat))] = dist_diff
        # Return the desired measure of difference between the two.
        if   (method == "max"):  return max(diffs.values())
        elif (method == "min"):  return min(diffs.values())
        elif (method == "mean"): return sum(diffs.values()) / len(diffs)
        else:                    return diffs


# ===============================================
#      Plotting Percentile Functions as CDFS     
# ===============================================

# Linearly interpolate to guess y values for x between provided data
def linear_fit_func(x_points, y_points):
    from scipy.interpolate import splrep, splev

    # Generate the 'range' of the function
    min_max = (min(x_points), max(x_points))

    # Fit the data with degree 1 splines (linear interpolation)
    fit = splrep(x_points, y_points, k=1)

    # Generate the function to return (gives the 'range' of
    # application when called without x values)
    def func(x_val=None, fit=fit, min_max=min_max):
        if type(x_val) == type(None):
            return min_max
        else:
            # ext 0 -> extrapolate   ext 1 -> return 0
            # ext 2 -> raise errror  ext 3 -> return boundary value
            return splev(x_val, fit, ext=3)
    def deriv(x_val, fit=fit):
        # der 1 -> compute the 1st order derivative
        # ext 0 -> extrapolate   ext 1 -> return 0
        # ext 2 -> raise errror  ext 3 -> return boundary value
        return splev(x_val, fit, der=1, ext=3)

    # Make the attribute "derivative" return a function for the derivative
    setattr(func, "derivative", lambda: deriv)

    # Return the linear interpolation function
    return func

# A percentile function that can handle series with non-numbers (by ignoring).
def robust_percentile(values, perc):
    import numbers
    no_none = [v for v in values if isinstance(v, numbers.Number)]
    if len(no_none) > 0:
        return np.percentile(no_none,perc)
    else:
        return None

# Given a (N,) array x_points, and a (N,M) array of y_points, generate
# a list of (M,) linear functions that represent the "percentiles" of
# y_points at each point in x_points.
def percentile_funcs(x_points, y_points, percentiles):
    # Generate the points for the percentile CDF's
    perc_points = np.array([
        [robust_percentile([row[i] for row in y_points], p) 
         for i in range(len(x_points))]
        for p in percentiles ])
    # Generate the CDF functions for each percentile and
    # return the fit functions for each of the percentiles
    return [linear_fit_func(x_points, pts) for pts in perc_points]

# Given a (N,) array x_points, and a (N,M) array of y_points, generate
# a list of (M,) linear functions that represent the "percentiles" of
# y_points at each point in x_points.
def percentile_points(x_points, y_points, percentiles):
    # Generate the points for the percentile CDF's
    perc_points = [ [robust_percentile(row, p) for row in y_points]
                    for p in percentiles ]
    # Generate the CDF functions for each percentile and
    # return the fit functions for each of the percentiles
    return perc_points

# Given a plot, this will add the 'percentiles' cloud given a (N,)
# array x_points, and a (N,M) array of y_points, generate a list of
# (M,) linear functions that represent the "percentiles" of y_points
# at each point in x_points.
def plot_percentiles(plot, name, x_points, y_points, color=None, 
                     percentiles=[0,20,40,60,80,100], line_width=1,
                     # center_color = np.array([255,70,0,0.6]),
                     # outer_color = np.array([255,150,50,0.3])):
                     center_color=None, outer_color=None, **kwargs):
    from util.plot import color_string_to_array
    # If there are multiple percentiles, use shading to depict them all.
    if (len(percentiles) > 1):
        # Generate the color spectrum if necessary
        plot.color_num +=1
        color = color_string_to_array(plot.color())
        if type(center_color) == type(None):
            center_color = color.copy()
            center_color[-1] = center_color[-1]*0.6
        if type(outer_color) == type(None):
            outer_color = color.copy()
            outer_color[:-1] = outer_color[:-1]*1.5 + 50
            outer_color[-1] = outer_color[-1]*0.3
        # Generate a group ID for the percentiles and ensure the
        # percentiles are sorted (for calculating gaps and colors)
        group_id = name + "_percentiles"
        percentiles.sort() # Ensure they are sorted
        # Find the min and max values and generate functions
        perc_pts = percentile_points(x_points, y_points, percentiles)
        # Identify the width of gaps between percentiles (for coloring)
        gaps = [percentiles[i] - percentiles[i-1] for i in range(1,len(percentiles))]
        text_color = 'rgb(100,100,100)'
        textfont = dict(color=text_color)
        # Add the last function (bounding the set)
        plot.add("%ith percentile"%percentiles[0], x_points,
                 perc_pts[0], color=text_color, line_width=line_width,
                 group=group_id, mode="lines", show_in_legend=False,
                 textfont=textfont, **kwargs)
        for pts,p,g in zip(perc_pts[1:], percentiles[1:], gaps):
            ratio = abs((50 - (p - g/2))/50)
            color = center_color * abs(1.0 - ratio) + outer_color * ratio
            color = 'rgba(%i,%i,%i,%f)'%tuple(color)
            plot.add("%ith percentile"%p, x_points, pts,
                     line_width=line_width, fill='toprevy',
                     color=text_color, fill_color=color,
                     group=group_id, show_in_legend=False,
                     mode="lines", textfont=textfont, **kwargs)
        # Add a master legend entry.
        plot.add(name + " Percentiles", [None], [None],
                 color='rgba(%i,%i,%i,%f)'%tuple(center_color),
                 group=group_id, **kwargs)
    else:
        if type(color) == type(None):
            color = plot.color()
            plot.color_num +=1
        perc_pts = percentile_points(x_points, y_points, percentiles)
        show_name = name+f" {percentiles[0]}th percentile"
        plot.add(show_name, x_points, perc_pts[0], color=color,
                 mode="lines", group=show_name, **kwargs)
    return plot



# ============================================================
#      Automatic detection of modes (based on confidence)     
# ============================================================

def modes(data, confidence=.99, tol=1/1000):
    from util.optimize import zero
    num_samples = len(data)
    error = 2*samples(num_samples, confidence=confidence)
    print()
    print("Smallest allowed mode: ",error)
    print()
    # Get the CDF points (known to be true based on data).
    x, y = cdf_points(data)
    cdf = cdf_fit((x,y), fit="linear")
    x, y = x[1:], y[1:]
    # Generate the candidate break points based on linear interpolation.
    slopes = [(y[i+1] - y[i]) / (x[i+1] - x[i]) for i in range(len(x)-1)]
    candidates = [i for i in range(1,len(slopes)-1)
                  if (slopes[i] < slopes[i-1]) and (slopes[i] < slopes[i+1])]
    # Sort candidates by their 'width', the widest being the obvious divisors.
    candidates = sorted(candidates, key=lambda i: -(x[i+1] - x[i]))
    print("candidates: ",candidates)
    print("slopes:     ",[slopes[c] for c in candidates])
    # Break the data at candidates as much as can be done with confidence.
    breaks = [min(x), max(x)]
    sizes = [1.]
    chosen = []
    print()
    print("breaks: ",breaks)
    print("sizes:  ",sizes)
    print("chosen: ",chosen)
    print()
    # Loop until there are no candidate break points left.
    while len(candidates) > 0:
        new_break_idx = candidates.pop(0)
        new_break = (x[new_break_idx + 1] + x[new_break_idx]) / 2
        b_idx = np.searchsorted(breaks, new_break, side="right")
        # Compute the CDF values at the upper, break, and lower positions.
        upp = cdf(breaks[b_idx+1]) if (b_idx < len(breaks)-1) else cdf.max
        mid = cdf(new_break)
        low = cdf(breaks[b_idx-1])
        print()
        print("new_break: ", new_break)
        print("b_idx: ",b_idx)
        print("  upp: ",upp, upp - mid)
        print("  mid: ",mid)
        print("  low: ",low, mid - low)
        # Compute the size of the smallest mode resulting from the break.
        smallest_result = min(upp - mid, mid - low)
        # Skip the break if it makes a mode smaller than error.
        if smallest_result < error: continue

        print()
        print("Num_modes: ", len(sizes))
        print("breaks:    ", breaks)
        print("sizes:     ", sizes)
        print()

        # Update the "breaks" and "sizes" lists.
        breaks.insert(b_idx, new_break)
        sizes.insert(b_idx, upp - mid)
        sizes[b_idx-1] = mid - low
        chosen.append(new_break_idx)
        
    # From the "breaks" and "sizes", construct a list of "modes".
    # Consider the most dense point between breaks the "mode".
    modes = []
    for i in range(len(chosen)):
        low = 0 if (i == 0) else chosen[i-1]
        upp = len(slopes)-1 if (i == len(chosen)-1) else chosen[i+1]
        mode_idx = low + np.argmax(slopes[low:upp])
        modes.append( (x[mode_idx+1] + x[mode_idx]) / 2 )


    from util.plot import Plot
    p = Plot()
    pdf = pdf_fit(data)
    p.add_func("PDF", pdf, pdf(), color=p.color(1))
    # Generate the mode lines.
    mode_lines = [[],[]]
    for z in modes:
        mode_lines[0] += [z,z,None]
        mode_lines[1] += [0,.2,None]
    p.add("modes", *mode_lines, color=p.color(0), mode="lines", group="modes")
    # Generate the antimode lines.
    break_lines = [[],[]]
    for z in breaks:
        break_lines[0] += [z,z,None]
        break_lines[1] += [0,.2,None]
    p.add("seperator", *break_lines, color=p.color(3,alpha=.3), mode="lines", group="seperator")
    p.show()
    # Show CDF
    p = Plot()
    pdf = pdf_fit(data)
    p.add_func("CDF", cdf, cdf(), color=p.color(1))
    # Generate the mode lines.
    mode_lines = [[],[]]
    for z in modes:
        mode_lines[0] += [z,z,None]
        mode_lines[1] += [0,1,None]
    p.add("modes", *mode_lines, color=p.color(0), mode="lines", group="modes")
    # Generate the antimode lines.
    break_lines = [[],[]]
    for z in breaks:
        break_lines[0] += [z,z,None]
        break_lines[1] += [0,1,None]
    p.add("seperator", *break_lines, color=p.color(3,alpha=.3), mode="lines", group="seperator")
    p.show(append=True)

            

    # Try and break the data up at the lowest density regions,
    # continue breaking the data until we cannot break it any more
    # based on confidence.




# ../development/testing/test_stats.py 
if __name__ == "__main__":
    from util.random import cdf
    np.random.seed(1)
    data = cdf()
    


    from util.random import cdf
    np.random.seed(1)
    data = cdf(nodes=3).inverse(np.random.random(100))
    modes(data)

    # Look for "largest" regions of separation 



# Backwards compatibility with warning for deprecation.
def cdf_fit_func(*args, **kwargs):
    print("\nWARNING: 'cdf_fit_func' is a deprecated function. Use 'cdf_fit' instead.\n")
    return cdf_fit(*args, **kwargs)

def pdf_fit_func(*args, **kwargs):
    print("\nWARNING: 'pdf_fit_func' is a deprecated function. Use 'pdf_fit' instead.\n")
    return pdf_fit(*args, **kwargs)


