import numpy as np
import os

# A numeric attirbute should have all of these operations.
NUMERIC_ATTRIBUTES = {'__abs__', '__sub__', '__pos__', '__neg__',
                      '__floordiv__', '__truediv__', '__pow__', 
                      '__radd__', '__rsub__', '__rfloordiv__',
                      '__rtruediv__', '__rpow__'}

# =============================================
#      Metric Principle Component Analysis     
# =============================================

# Custom excepion for improper user-specified ranges.
class InvalidRange(Exception):
    def __init__(self, start, stop, step):
        return super().__init__("\n\nError: random_range(start=%s,stop=%s,step=%s) contains no elements."%(start,stop,step))

# Prime number used for generating random sequences.
ABS_DIFF = lambda p1,p2: abs(p1-p2)
        
# Return all unique prime factors of a number in sorted order.
def primes(num):
    factors = {}
    candidate = 2
    while candidate**2 <= num:
        while (num % candidate) == 0:
            factors[candidate] = factors.get(candidate,0) + 1
            num //= candidate
        candidate += 1
    if num > 1: factors.add(num)
    return sorted((k,factors[k]) for k in factors)

# Return a randomized "range" using a Linear Congruential Generator
# to produce the number sequence. Parameters are the same as for 
# python builtin "range".
#   Memory  -- O(1) storage for a few integers, regardless of parameters.
#   Compute -- O(n) at most 2 times the number of steps in the range, n.
# 
def random_range(start, stop=None, step=None):
    from random import randint
    from math import ceil, log2
    # Set a default values the same way "range" does.
    if (stop == None): start, stop = 0, start
    if (step == None): step = 1
    # Use a mapping to convert a standard range into the desired range.
    mapping = lambda i: (i*step) + start
    # Compute the number of numbers in this range.
    num_steps = (stop - start) // step
    # Check for a usage error.
    if (num_steps == 0): raise(InvalidRange(start, stop, step))
    # Seed range with a random integer.
    value = randint(0,num_steps)
    # 
    # Construct an offset, multiplier, and modulus for a linear
    # congruential generator. These generators are cyclic and
    # non-repeating when they maintain the properties:
    # 
    #   1) "modulus" and "offset" are relatively prime.
    #   2) ["multiplier" - 1] is divisible by all prime factors of "modulus".
    #   3) ["multiplier" - 1] is divisible by 4 if "modulus" is divisible by 4.
    # 
    offset = randint(0,num_steps) * 2 + 1                 # Pick a random odd-valued offset.
    multiplier = 4*(num_steps + randint(0,num_steps)) + 1 # Pick a multiplier 1 greater than a multiple of 4.
    modulus = 2**ceil(log2(num_steps))                    # Pick a modulus just big enough to generate all numbers (power of 2).
    # Track how many random numbers have been returned.
    found = 0
    while found < num_steps:
        # If this is a valid value, yield it in generator fashion.
        if value < num_steps:
            found += 1
            yield mapping(value)
        # Calculate the next value in the sequence.
        value = (value*multiplier + offset) % modulus
    # yield (offset, multiplier) # <-- For testing purposes.

# This function maps an index in the range [0, (count**2 - count) // 2] 
# to a tuple of integers in the range [0,count). The mapping is complete.
def index_to_pair(index):
    val = int(((1/4 + 2*index)**(1/2) + 1/2))
    remainder = index - val*(val - 1)//2
    return (val, remainder)

# Generate "count" random indices of pairs that are within "length" bounds.
def gen_random_pairs(length, count=None):
    # Compute the hard maximum in the number of pairs.
    max_pairs = length*(length - 1) // 2
    # Initialize count if it wasn't provided.
    if type(count) == type(None): count = max_pairs
    count = min(count, max_pairs)
    # Get a random set of pairs (no repeats).
    for c,i in enumerate(random_range(count)):
        if (i >= count): break
        yield index_to_pair(i)
    print(" "*40, end="\r")

# Generate vector between scaled by metric difference. Give the metric
# the indices of "vectors" in the provided matrix.
def gen_random_metric_diff(matrix, index_metric, count=None):
    # Iterate over random pairs.
    for (p1, p2) in gen_random_pairs(len(matrix), count):
        vec = matrix[p1] - matrix[p2]
        length = np.sqrt(np.sum(vec**2))
        if length > 0: vec /= length**2
        yield index_metric(p1, p2) * vec

# Compute the metric PCA (pca of the between vectors scaled by 
# metric difference slope).
def mpca(points, values, metric=ABS_DIFF, num_components=None,
         num_vecs=None, display=True):
    # Set default values for num_components and num_vecs.
    if type(num_components) == type(None): num_components = min(points.shape)
    if (type(num_vecs) == type(None)): num_vecs = (len(points)**2-len(points))//2
    if display: print(" allocating memory for metric between vectors..", end="\r", flush=True)
    m_vecs = np.zeros((num_vecs, num_components))
    if display: print(" generating metric vectors..", end="\r", flush=True)
    # Function that takes two indices and returns metric difference.
    index_metric = lambda i1, i2: metric(values[i1], values[i2])
    # Generator that produces "between vectors".
    vec_gen = gen_random_metric_diff(points, index_metric, count=num_vecs)
    for i,vec in enumerate(vec_gen): m_vecs[i,:] = vec
    # Compute the principle components of the between vectors.
    if display: print(" computing principle components..", end="\r", flush=True)
    components, _ = pca(m_vecs, num_components=num_components)
    # Compute the magnitudes using total variation.
    if display: print(" computing total variation per component..", end="\r", flush=True)
    avg_diffs = np.zeros(num_components)
    for i in range(num_components):
        ordering = np.argsort(np.matmul(points, components[i]))
        for j in range(len(ordering)-1):
            p1, p2 = ordering[j], ordering[j+1]
            avg_diffs[i] += metric(values[p1], values[p2])
    # Invert the total variation.
    if (min(avg_diffs) <= 0.0):
        avg_diffs = np.where(avg_diffs == 0, 1., 0.)
    else:
        avg_diffs = 1 / avg_diffs
    avg_diffs /= np.sum(avg_diffs)
    # Re-order the components according to inverse total variation.
    order = np.argsort(avg_diffs)[::-1]
    if display: print(" reordering the components by total variation..", end="\r", flush=True)
    # If they are not already ordered correctly, re-order the returns.
    if not all(order[i] < order[i+1] for i in range(len(order)-1)):
        components, avg_diffs = components[order], avg_diffs[order]
    if display: print("                                               ", end="\r", flush=True)
    # Return the principle components of the metric slope vectors.
    return components, avg_diffs

# Compute the principle components using sklearn.
def pca(points, num_components=None):
    from sklearn.decomposition import PCA        
    pca = PCA(n_components=num_components)
    pca.fit(points)
    principle_components = pca.components_
    magnitudes = pca.singular_values_
    # Normalize the component magnitudes to have sum 1.
    magnitudes /= sum(magnitudes)
    return principle_components, magnitudes

# ===================================
#      CDF and PDF Fit Functions     
# ===================================

# Construct a version of a function that has been smoothed with a gaussian.
def gauss_smooth(func, stdev=1., n=100):
    # We must incorporate ssome smoothing.
    from scipy.stats import norm
    # Construct a set of gaussian weights (to apply nearby on the function).
    eval_pts = np.linspace(-5 * stdev, 5 * stdev, n)
    weights = norm(0, stdev).pdf(eval_pts)
    # Generate a smoothed function.
    def smooth_func(x): return sum(weights * func(eval_pts + x)) / sum(weights)
    # Transfer attributes over to the smooth function.
    standard_attributes = set(dir(smooth_func))
    for attr in dir(func):
        if attr not in standard_attributes:
            setattr(smooth_func, attr, getattr(func, attr))
    # Return the smoothed function.
    return smooth_func

# Returns the CDF function value for any given x-value
#   fit -- "linear" linearly interpolates between Emperical Dist points.
#          "cubic" constructs a monotonic piecewise cubic interpolating spline
#          <other> returns the raw Emperical Distribution Function.
def cdf_fit(data, fit="linear"):
    fit_type = fit
    # Scipy functions for spline fits.
    from scipy.interpolate import splrep, splev
    from scipy.interpolate import PchipInterpolator

    # Sort the data (without changing it) and get the min and max
    data = np.array(sorted(data))
    min_pt = data[0]
    max_pt = data[-1]

    # Initialize a holder for all CDF points
    data_vals = []

    # Add all of the CDF points for the data set
    for i,val in enumerate(data):
        if ((i > 0) and (val == data[i-1])): continue
        data_vals.append( i/(len(data) - 1) )
    # Add the 100 percentile point if it was not added already
    if (data_vals[-1] != 1.0): data_vals[-1] = 1.0

    # Convert data into its unique values.
    data = np.array(sorted(set(data)))
    # Convert it all to numpy format for ease of processing
    data_vals = np.array(data_vals)

    # Generate a fit for the data points
    if (fit_type == "linear"):
        # Generate linear function
        fit_params = splrep(data, data_vals, k=1)
        fit = lambda x_val: splev(x_val, fit_params, ext=3)
        fit.derivative = lambda d=1: lambda x_val: splev(x_val, fit_params, der=d, ext=3)
        # Generate inverse linear function
        inv_fit_params = splrep(data_vals, data, k=1)
        inv_fit = lambda y_val: splev(y_val, inv_fit_params, ext=3)
        inv_fit.derivative = lambda d=1: lambda y_val: splev(y_val, inv_fit_params, der=d, ext=3)
        fit.inverse = inv_fit
    elif (fit_type == "cubic"):
        # Generate piecewise cubic monotone increasing spline
        fit = PchipInterpolator(data, data_vals)
        fit.inverse = PchipInterpolator(data_vals, data)
    else:
        # Construct the empirical distribution fit.
        def fit(x_val, data=data, data_vals=data_vals):
            try:
                # Handle the vector valued case
                len(x_val)
                return np.array([np.searchsorted(data, x, side="right") 
                                 / len(data) for x in x_val])
            except TypeError:
                # Handle the single value case
                return np.searchsorted(data, x_val, side="right") / len(data)
        # Construct an inverse function for the EDF.
        def inverse(perc, data=data, data_vals=data_vals):
            try:
                # Handle the vector valued case
                len(perc)
                indices = []
                for p in perc:
                    index = np.searchsorted(data_vals, max(0,min(p,1)), side="right")
                    indices.append( data[min(index, len(data)-1)] )
                return np.array(indices)
            except TypeError:
                # Handle the single value case
                index = np.searchsorted(data_vals, max(0,min(perc,1)), side="right")
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
        fit.derivative = lambda d=1: (derivative) if (d == 1) else (lambda: 0.)
    
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

    # Return the custom function for this set of points
    return cdf_func

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
        smooth_cdf_func = gauss_smooth(cdf_func, stdev, n)
        eps = .01 * width
        def pdf_func(x):
            return (smooth_cdf_func(x + eps) - smooth_cdf_func(x - eps))
    else:
        # Return the PDF as the derivative of the CDF.
        pdf_func = cdf_func.derivative(1)
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
    return pdf


# ===============================================
#      Statistical Measurements Involving KS     
# ===============================================

# Calculate the maximum difference between two CDF functions (two sample)
def ks_diff(test_func, true_func, method="util"):
    from util.optimize import minimize
    # Cycle through the functions to find the min and max of all ranges
    min_pt, max_pt = true_func()
    if method in {"scipy", "util"}:
        diff_func = lambda x: -abs(test_func(x) - true_func(x))
        if method == "scipy":
            # METHOD 1:
            #  Use scipy to maximize the difference function between
            #  the two cdfs in order to find the greatest difference.
            sol = minimize(diff_func, [(max_pt - min_pt) / 2],
                           bounds=[(min_pt,max_pt)], method='L-BFGS-B').x
        elif method == "util":
            # METHOD 2 (default):
            #  Use the default minimizer in "optimize" to maximize the
            #  difference function between the two cdfs in order to
            #  find the greatest difference.
            sol = minimize(diff_func, [(max_pt - min_pt) / 2],
                           bounds=[(min_pt,max_pt)])
        greatest_diff = abs(test_func(sol) - true_func(sol))
    else:
        # METHOD 3:
        #  Generate a large set of x-points and find the difference
        #  between the functions at all of those points. Generate a
        #  grid of points and "zoom in" around the greatest difference
        #  points to identify the spot with largest gap.
        x_points = np.linspace(min_pt, max_pt, 1000)
        diff = abs(test_func(x_points) - true_func(x_points))

        greatest_diff = -float('inf')
        while (diff[np.argmax(diff)] > greatest_diff):
            lower = np.argmax(diff) - 1
            upper = np.argmax(diff) + 1
            min_pt = x_points[max(lower, 0)] - (
                1 if lower < 0 else 0)
            max_pt = x_points[min(upper,len(x_points)-1)] + (
                1 if upper >= len(x_points) else 0)
            x_points = np.linspace(min_pt, max_pt, 1000)
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

# Compute the pairwise difference between two sequences. Use the following:
#    category vs. category -- Average total difference between full distribution
#                             and conditional distributions for any setting.
#    category vs. number   -- Average 1-norm difference between full distribution
#                             and conditional distributions for categories.
#    number vs. number     -- Correlation coefficient between the two numbers.
def diff(seq_1, seq_2):
    # Get the set of attributes a numeric type should have.
    # Return a boolean "is_numeric"
    def is_numeric(obj):
        attributes = set(dir(obj))
        return all(attr in attributes for attr in NUMERIC_ATTRIBUTES)
    # The three different outcomes of comparison could be:
    cat_cat = (is_numeric(seq_1) and is_numeric(seq_2))
    num_num = (not is_numeric(seq_1)) and (not is_numeric(seq_2))
    num_cat = (not cat_cat) and (not num_num)
    # Execute the appropriate comparison according to the types.
    if num_num:
        pass
        # number vs. number
        # - compute mean of each sequence
        # - compute covariance
        # - divide by product of covariances
    elif num_cat:
        pass
        # number vs. categorical
        # - compute total EDF for the numer sequence
        # - compute 1-norm difference between EDF
    elif cat_cat:
        pass
        # categorical vs. categorical


# diff([1,2,3,4],['a','b','c'])
# exit()

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



# Backwards compatibility with warning for deprecation.
def cdf_fit_func(*args, **kwargs):
    print("\nWARNING: 'cdf_fit_func' is a deprecated function. Use 'cdf_fit' instead.\n")
    return cdf_fit(*args, **kwargs)

def pdf_fit_func(*args, **kwargs):
    print("\nWARNING: 'pdf_fit_func' is a deprecated function. Use 'pdf_fit' instead.\n")
    return pdf_fit(*args, **kwargs)



if __name__ == "__main__":
    # # Use a linear interpolation of the EDF points to compute the integral
    # # of the absolute difference between the empirical PDF's of two sequences.
    # def epdf_diff(seq_1, seq_2):
    #     # Construct generators for the EDF point pairs.
    #     gen_1 = edf_pair_gen(seq_1)
    #     gen_2 = edf_pair_gen(seq_2)
    #     # Get the initial values from the generators.
    #     low_1, upp_1, density_1 = gen_1.__next__()
    #     low_2, upp_2, density_2 = gen_2.__next__()
    #     shared = 0.
    #     # Cycle until completing the integral difference between pEDFs.
    #     while (upp_1 < float('inf')) or (upp_2 < float('inf')):
    #         # Compute the width of the interval and the percentage of data
    #         # within the width for each of the distributions.
    #         width = min(upp_1, upp_2) - max(low_1, low_2)
    #         if width < float('inf'):
    #             # Convert the EDF points into densities over the interval.
    #             sub_density_1 = density_1 * (width / (upp_1-low_1))
    #             sub_density_2 = density_2 * (width / (upp_2-low_2))
    #             if abs(sub_density_1 - sub_density_2) >= 0:
    #                 shared += min(sub_density_1, sub_density_2)
    #         # Cycle whichever range has the lower maximum.
    #         if (upp_1 > upp_2):
    #             low_2, upp_2, density_2 = gen_2.__next__()
    #         elif (upp_2 > upp_1):
    #             low_1, upp_1, density_1 = gen_1.__next__()
    #         else:
    #             low_1, upp_1, density_1 = gen_1.__next__()
    #             low_2, upp_2, density_2 = gen_2.__next__()
    #     return max(0, 1 - shared)

    # # Generator for returning pairs of EDF points as (val1, val2, density).
    # def edf_pair_gen(seq):
    #     # Compute the 'extra' slope lost at the beginning of the function.
    #     for i in range(len(seq)):
    #         if seq[i] != seq[0]: break
    #     extra = (i/len(seq))
    #     shift = lambda v: (v - extra) / (1 - extra)
    #     # Manually set the first slope value.
    #     x1, y1 = -float('inf'), 0
    #     for i in range(len(seq)):
    #         # Cycle until we hit the last occurrence of the next value.
    #         if ((i+1) < len(seq)) and (seq[i+1] == seq[i]): continue
    #         x2, y2 = seq[i], shift( (i+1)/len(seq) )
    #         yield (x1, x2, y2 - y1)
    #         # Cycle values.
    #         x1, y1 = x2, y2
    #     # Manually yield the last element in the sequence.
    #     yield (x1, float('inf'), 0)


    # # ----------------------------------------------------------------
    # def demo(seq):
    #     print('-'*70)
    #     print(seq)
    #     print()
    #     total = 0
    #     for vals in edf_pair_gen(seq):
    #         total += vals[-1]
    #         print("[% 4s, % 3s] (%.2f) --"%vals, round(total,3))
    #     print('-'*70)
    #     print()
    # print()
    # seq = [0] + list(range(9))
    # demo(seq)
    # print()
    # seq = sorted(np.random.random(size=(10,)))
    # demo(seq)
    # # ----------------------------------------------------------------
    # # a = [1, 1, 3, 3, 5, 6]
    # # b = [0, 1, 2, 3, 4]
    # # 
    # n = 100
    # a = [v / n for v in list(range(n+1))]
    # b = [v+.5 for v in a]
    # # b = a
    # # b = [v+10 for v in a]

    # print(epdf_diff(a,a))
    # print(epdf_diff(a,b))
    # print(epdf_diff(b,a))
    # print()

    # # Generate a random sequence.
    # a = sorted((np.random.random(size=(10,))))
    # b = sorted((np.random.random(size=(10,)) + .5))
    # print("[%.2f, %.2f]"%(min(a),max(a)), "[%.2f, %.2f]"%(min(b),max(b)))
    # diff = epdf_diff(a, b)
    # print(diff)
    # exit()
    # # ----------------------------------------------------------------
    # from util.plot import Plot

    # p = Plot()
    # p1 = pdf_fit(a, smooth=0.00001)
    # p2 = pdf_fit(b, smooth=0.00001)
    # p.add_func("a", p1, p1()) #(-.5,2))
    # p.add_func("b", p2, p2()) #(-.5,2))
    # p.show(show=False, y_range=[-.5,1.5])

    # p = Plot()
    # p1 = cdf_fit(a)
    # p2 = cdf_fit(b)
    # p.add_func("a", p1, p1()) #(-.5,2))
    # p.add_func("b", p2, p2()) #(-.5,2))
    # p.show(append=True)
    # exit()
    # # ----------------------------------------------------------------


    import pickle, os
    import numpy as np
    from util.plot import Plot, multiplot
    
    TEST_FIT_FUNCS = False
    TEST_MPCA = True
    TEST_RANDOM_RANGE = False

    if TEST_RANDOM_RANGE:
        vals = {}
        mults = {}
        for i in range(10000):
            results = tuple(random_range(20))
            key, (val, mult) = results[:-1], results[-1]
            diff = set(key[i-1] - key[i] for i in range(len(key)))
            if len(diff) <= 2:
                vals[key] = vals.get(key, set())   .union({val})
                mults[key] = mults.get(key, set()) .union({mult})
        for v in sorted(vals): print(v, "% 10s"%sorted(vals[v]), sorted(mults[v]))
        print(len(vals))


    if TEST_FIT_FUNCS:
        # ==============================================
        #      Test the fit functions and smoothing     
        # ==============================================
        # Make data
        smooth = .1
        data = np.random.normal(size=(1000,))
        # data[:len(data)//2] += 2
        min_max = (min(data) - .1, max(data) + .1)
        print()
        print("(min, max) : (%.2f, %.2f)"%(min_max))
        print("Normal confidence: %.2f%%"%(100*normal_confidence(data)))
        print()
        # Make PDF plot
        pfit = pdf_fit(data)
        p = Plot()
        p.add_func("PDF", pfit, min_max)
        # Make smooth PDF
        pfit = pdf_fit(data, smooth=smooth)
        p.add_func("Smooth PDF", pfit, min_max)
        p.show(show=False)

        # Make CDF plota
        cfit = cdf_fit(data)
        p = Plot()
        p.add_func("CDF", cfit, min_max)
        # Make CDF whose derivative is the default PDF.
        stdev = .05 * (cfit.max - cfit.min)
        cfit = gauss_smooth(cfit, stdev)
        p.add_func("CDF for default PDF", cfit, min_max)
        # Make smooth cdf
        stdev = smooth * (cfit.max - cfit.min)
        cfit = gauss_smooth(cfit, stdev)
        p.add_func("Smooth CDF", cfit, min_max)
        p.show(append=True)

    # Example / Test case of Principle Response Analyais
    if TEST_MPCA:

        GENERATE_APPROXIMATIONS = False
        BIG_PLOTS = False
        SHOW_2D_POINTS = False

        import random
        # Generate some points for testing.
        np.random.seed(4)
        random.seed(0)
        rgen = np.random.RandomState(10)
        n = 100
        points = (rgen.rand(n,2) - .5) * 2
        # points *= np.array([.5, 1.])

        # Create some testing functions (for learning different behaviors)
        funcs = [
            lambda x: x[1]               , # Linear on y
            lambda x: abs(x[0] + x[1])   , # "V" function on 1:1 diagonal
            lambda x: abs(2*x[0] + x[1]) , # "V" function on 2:1 diagonal
            lambda x: x[0]**2            , # Quadratic on x
            lambda x: (x[0] + x[1])**2   , # Quadratic on 1:1 diagonal
            lambda x: (2*x[0] + x[1])**3 , # Cubic on 2:1 diagonal
            lambda x: (x[0]**3)          , # Cubic on x
            lambda x: rgen.rand()        , # Random function
        ]
        # Calculate the response values associated with each function.
        responses = np.vstack(tuple(tuple(map(f, points)) for f in funcs)).T

        # Reduce to just the first function
        choice = 1
        func = funcs[choice]
        response = responses[:,choice]

        # Run the princinple response analysis function.
        components, values = mpca(points, response)
        conditioner = np.matmul(components, np.diag(values))

        print()
        print("Components")
        print(components)
        print()
        print("Values")
        print(values)
        print()
        print("Conditioner")
        print(conditioner)
        print()

        # Generate a plot of the response surfaces.
        from util.plot import Plot, multiplot
        print("Generating plots of source function..")

        # Add function 1
        p1 = Plot()
        p1.add("Points", *(points.T), response, opacity=.8)
        p1.add_func("Surface", func, [-1,1], [-1,1], plot_points=100)

        if GENERATE_APPROXIMATIONS:
            from util.algorithms import NearestNeighbor
            model = NearestNeighbor()
            model.fit(points, response)
            p1.add_func("Unconditioned Approximation", model, [-1,1], [-1,1],
                        mode="markers", opacity=.8)
            # Generate a conditioned approximation
            model = NearestNeighbor()
            model.fit(np.matmul(conditioner, points), response)
            approx = lambda x: model(np.matmul(conditioner, x))
            p1.add_func("Best Approximation", approx, [-1,1], [-1,1],
                        mode="markers", opacity=.8)


        print("Generating metric principle components..")

        # Return the between vectors and the differences between those points.
        def between(x, y, unique=True):
            vecs = []
            diffs = []
            for i1 in range(x.shape[0]):
                start = i1+1 if unique else 0
                for i2 in range(start, x.shape[0]):
                    if (i1 == i2): continue
                    vecs.append(x[i2] - x[i1])
                    diffs.append(y[i2] - y[i1])
            return np.array(vecs), np.array(diffs)

        # Plot the between slopes to verify they are working.
        # Calculate the between slopes
        vecs, diffs = between(points, response)
        vec_lengths = np.sqrt(np.sum(vecs**2, axis=1))
        between_slopes = diffs / vec_lengths
        bs = ((vecs.T / vec_lengths) * between_slopes).T
        # Extrac a random subset for display
        size = 100
        random_subset = np.arange(len(bs))
        rgen.shuffle(random_subset)
        bs = bs[random_subset[:size],:]
        # Normalize the between slopes so they fit on the plot
        max_bs_len = np.max(np.sqrt(np.sum(bs**2, axis=1)))
        bs /= max_bs_len
        # Get a random subset of the between slopes and plot them.
        p2 = Plot("","Metric PCA on Z","")
        p2.add("Between Slopes", *(bs.T), color=p2.color(4, alpha=.4))

        if SHOW_2D_POINTS:
            # Add the points and transformed points for demonstration.
            new_pts = np.matmul(np.matmul(conditioner, points), np.linalg.inv(components))
            p2.add("Original Points", *(points.T))
            p2.add("Transformed Points", *(new_pts.T), color=p2.color(6, alpha=.7))

        # Add the principle response components 
        for i,(vec,m) in enumerate(zip(components, values)):
            vec = vec * m
            p2.add(f"PC {i+1}", [0,vec[0]], [0,vec[1]], mode="lines")
            ax, ay = (vec / sum(vec**2)**.5) * 3
            p2.add_annotation(f"{m:.2f}", vec[0], vec[1])


        p3 = Plot("", "PCA on X", "")
        p3.add("Points", *(points.T), color=p3.color(4, alpha=.4))

        # Add the normal principle components
        components, values = pca(points)
        for i,(vec,m) in enumerate(zip(components, values)):
            vec = vec * m
            p3.add(f"PC {i+1}", [0,vec[0]], [0,vec[1]], mode="lines")
            ax, ay = (vec / sum(vec**2)**.5) * 3
            p3.add_annotation(f"{m:.2f}", vec[0], vec[1])


        if BIG_PLOTS:
            p1.plot(file_name="source_func.html", show=False)
            p2.plot(append=True, x_range=[-8,8], y_range=[-5,5])
        else:
            # Make the plots (with manual ranges)
            p1 = p1.plot(html=False, show_legend=False)
            p2 = p2.plot(html=False, x_range=[-1,1], y_range=[-1,1], show_legend=False)
            p3 = p3.plot(html=False, x_range=[-1,1], y_range=[-1,1], show_legend=False)
            # Generate the multiplot of the two side-by-side figures
            multiplot([p1,p2,p3], height=126, width=600)
