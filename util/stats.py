import numpy as np
import os

# =============================================
#      Metric Principle Component Analysis     
# =============================================

# Prime number used for generating random sequences.
ABS_DIFF = lambda p1,p2: abs(p1-p2)
        
# Return all unique prime factors of a number in sorted order.
def primes(num):
    factors = set()
    candidate = 2
    while candidate**2 <= num:
        while (num % candidate) == 0:
            factors.add(candidate)
            num //= candidate
        candidate += 1
    if num > 1: factors.add(num)
    return sorted(factors)

# Return a randomized "range" using a Linear Congruential Generator
# to produce the number sequence. Parameters are the same as for 
# python builtin "range".
#   Memory  -- storage for 8 integers, regardless of parameters.
#   Compute -- at most 2*"maximum" steps required to generate sequence.
# 
def random_range(start, stop=None, step=None):
    import random, math
    # Set a default values the same way "range" does.
    if (stop == None): start, stop = 0, start
    if (step == None): step = 1
    # Use a mapping to convert a standard range into the desired range.
    mapping = lambda i: (i*step) + start
    # Compute the number of numbers in this range.
    maximum = (stop - start) // step
    # Seed range with a random integer.
    value = random.randint(0,maximum)
    # 
    # Construct an offset, multiplier, and modulus for a linear
    # congruential generator. These generators are cyclic and
    # non-repeating when they maintain the properties:
    # 
    #   1) "modulus" and "offset" are relatively prime.
    #   2) ["multiplier" - 1] is divisible by all prime factors of "modulus".
    #   3) ["multiplier" - 1] is divisible by 4 if "modulus" is divisible by 4.
    # 
    offset = random.randint(0,maximum) * 2 + 1      # Pick a random odd-valued offset.
    multiplier = 4*(maximum//4 + random.randint(0,maximum)) + 1                 # Pick a multiplier 1 greater than a multiple of 4.
    modulus = int(2**math.ceil(math.log2(maximum))) # Pick a modulus just big enough to generate all numbers (power of 2).
    # Track how many random numbers have been returned.
    found = 0
    while found < maximum:
        # If this is a valid value, yield it in generator fashion.
        if value < maximum:
            found += 1
            yield mapping(value)
        # Calculate the next value in the sequence.
        value = (value*multiplier + offset) % modulus









total = []
for i in range(100000):
    seq = tuple()
    for k,v in enumerate(random_range(10)):
        seq += (v,)
        if (k >= 1): break
    total.append(seq)
print()
for l in sorted(set(total)):
    print(l)
exit()

# This function maps an index in the range [0, (count**2 - count) // 2] 
# to a tuple of integers in the range [0,count). The mapping is complete.
def index_to_pair(index):
    val = int(((1/4 + 2*index)**(1/2) + 1/2))
    remainder = index - val*(val - 1)//2
    return (val, remainder)

# Generate "count" random indices of pairs that are within "length" bounds.
def gen_random_pairs(length, count=None, show_time=True, timeout=1):
    import time
    start = time.time()
    # Compute the hard maximum in the number of pairs.
    max_pairs = length*(length - 1) // 2
    # Initialize count if it wasn't provided.
    if type(count) == type(None): count = max_pairs
    count = min(count, max_pairs)
    # Get a random set of pairs (no repeats).
    for i in random_range(count, maximum=max_pairs):
        if (time.time() - start) >= timeout:
            print(f"{i: 7d} : {count: 7d}", end="\r")
            start = time.time()
        yield index_to_pair(i)
    print(" "*40, end="\r")

# Generate vector between scaled by metric difference. Give the metric
# the indices of "vectors" in the provided matrix.
def gen_random_metric_diff(matrix, index_metric, count=None):
    # Iterate over random pairs.
    for (p1, p2) in gen_random_pairs(len(matrix), count):
        vec = matrix[p1] - matrix[p2]
        vec /= np.sqrt(np.sum(vec**2))
        yield index_metric(p1, p2) * vec

# Compute the metric PCA (pca of the between vectors scaled by 
# metric difference slope).
def mpca(points, values, metric=ABS_DIFF, num_components=None, num_vecs=None):
    # Set default values for steps and directions.
    if type(num_vecs) == type(None):       num_vecs = 10 * points.shape[0]
    if type(num_components) == type(None): num_components = points.shape[1]
    # Function that takes two indices and returns metric difference.
    index_metric = lambda i1, i2: metric(values[i1], values[i2])
    # Generator that produces "between vectors".
    vec_gen = gen_random_metric_diff(points, index_metric, count=num_vecs)
    # Return the principle components of the metric slope vectors.
    return pca(np.array([v for v in vec_gen]), num_components=num_components)

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

# Returns the CDF function value for any given x-value
#   fit -- "linear" linearly interpolates between Emperical Dist points.
#          "cubic" constructs a monotonic piecewise cubic interpolating spline
#          <other> returns the raw Emperical Distribution Function.
def cdf_fit_func(data, fit="linear"):
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
            y_val = fit(x_val)
            if (type(x_val) == np.ndarray):
                y_val = np.where(x_val < min_pt, 0.0, y_val)
                y_val = np.where(x_val > max_pt, 1.0, y_val)
            return y_val

    # Set the min and max of the data.
    cdf_func.min = min_pt
    cdf_func.max = max_pt
    # Set the derivative function
    cdf_func.derivative = fit.derivative
    # Set the inverse function
    cdf_func.inverse = fit.inverse

    # Return the custom function for this set of points
    return cdf_func

# Construct a version of a function that has been smoothed with a gaussian.
def gauss_smooth(func, stdev=1., n=100):
    # We must incorporate ssome smoothing.
    from scipy.stats import norm
    # Construct a set of gaussian weights (to apply nearby on the function).
    eval_pts = np.linspace(-5 * stdev, 5 * stdev, n)
    weights = norm(0, stdev).pdf(eval_pts)
    # Return the smoothed function.
    return lambda x: sum(weights * func(eval_pts + x)) / sum(weights)

# Return a PDF fit function that is smoothed with a gaussian kernel.
# "smooth" is the percentage of the data width to use as the standard
# deviation of the gaussian kernel. "n" is the number of samples to
# make when doing the gaussian kernel smoothing.
def pdf_fit_func(data, smooth=.05, n=1000, **cdf_fit_kwargs):
    cdf_fit = cdf_fit_func(data, **cdf_fit_kwargs)
    pdf_func = cdf_fit.derivative(1)
    # Take the first derivative of the CDF function to get the PDF
    def pdf_fit(x=None):
        # Return the min and max when nothing is provided.
        if (type(x) == type(None)): return (cdf_fit.min, cdf_fit.max)
        try:
            # Handle the vector valued case.
            len(x)                
            vals = np.array([pdf_func(v) for v in x])
            vals[x < cdf_fit.min] = 0.
            vals[x > cdf_fit.max] = 0.
            # Return array of results.
            return vals
        except TypeError:
            # Return zeros for queries outside the range.
            if (cdf_fit.min > x) or (cdf_fit.max < x): return 0.0
            # Return the smoothed evaluation.
            return pdf_func(x)
        return pdf_func(x)
    # Smooth the pdf fit if that was requested.
    if smooth:
        width = (cdf_fit.max - cdf_fit.min)
        stdev = smooth * width
        smooth_cdf_fit = gauss_smooth(cdf_fit, stdev, n)
        eps = .01 * width
        def pdf_fit(x):
            return smooth_cdf_fit(x + eps) - smooth_cdf_fit(x - eps)
    # Set the "min" and "max" attributes of this function.
    pdf_fit.min = cdf_fit.min
    pdf_fit.max = cdf_fit.max
    pdf_fit.integral = cdf_fit
    return pdf_fit


# ===============================================
#      Statistical Measurements Involving KS     
# ===============================================

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

def normal_confidence(distribution):
    from scipy.stats import kstest
    # # Make the distribution 0 mean and unit variance (unit standard deviation)
    # new_distribution = (distribution - np.mean(distribution)) / np.var(distribution)
    # Compare the distribution with a normal distribution
    new_distribution = distribution
    ks_statistic = kstest(new_distribution, "norm").statistic
    return ks_p_value(ks_statistic, len(distribution))

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



if __name__ == "__main__":
    import pickle, os
    
    TEST_FIT_FUNCS = True
    TEST_CDF_PDF = False
    TEST_PRA = False

    if TEST_FIT_FUNCS:
        # ==============================================
        #      Test the fit functions and smoothing     
        # ==============================================
        from util.plot import Plot
        import numpy as np
        # Make data
        smooth = .1
        data = np.random.normal(size=(1000,))
        min_max = (min(data) - .1, max(data) + .1)
        print(min_max)

        # Make PDF plot
        pfit = pdf_fit_func(data)
        p = Plot()
        p.add_func("PDF", pfit, min_max)
        # Make smooth PDF
        pfit = pdf_fit_func(data, smooth=smooth)
        p.add_func("Smooth PDF", pfit, min_max)
        p.show(show=False)

        # Make CDF plota
        cfit = cdf_fit_func(data)
        p = Plot()
        p.add_func("CDF", cfit, min_max)
        # Make smooth cdf
        stdev = smooth * (cfit.max - cfit.min)
        cfit = gauss_smooth(cfit, stdev)
        p.add_func("Smooth CDF", cfit, min_max)
        p.show(append=True)


    if TEST_CDF_PDF:
        from util.plot import Plot, multiplot

        # Testing out the mathematical routines
        num_modes = 5
        top_mode_color = "rgba(20,20,20,0.6)"
        rest_mode_color = "rgba(210,210,210,0.4)"
        size = 200
        values = np.vstack((np.random.normal(0.4, 0.05, size=(size,)),
                            np.random.normal(0.0, 0.02, size=(size,)))).flatten()
        values.sort()

        raw = Plot("", "Value", "Normalized Scale")
        cdf = Plot("", "Value", "Percent Less")
        pdf = Plot("CDF and PDF of random normal 2-mode data", "Value", "Probability of Occurrence")

        min_max = (min(values), max(values))
        bin_size = (min_max[1]-min_max[0])/100
        pdf.add("PDF Histogram", x_values=values, plot_type="histogram",group="1",
                opacity=0.7, autobinx=False, histnorm='probability',
                xbins=dict(start=min_max[0],end=min_max[1],size=bin_size))
        cdf.add_func("CDF", cdf_fit_func(values, cubic=False), min_max,
                     mode="markers+lines", marker_size=4, group="1")

        sd_pts = cdf_second_diff(values)

        # Add vertical lines for each of the mode locations
        mode_pts = mode_points(sd_pts, thresh_perc=50.0)
        # Add the rest of the modes as another series
        mode_x_vals = []
        mode_y_vals = []
        for pt in mode_pts[num_modes:]:
            mode_x_vals += [pt[0], pt[0], None]
            mode_y_vals += [0,     1,     None]
        cdf.add("Remaining %i Modes"%(len(mode_pts) - num_modes),
                mode_x_vals, mode_y_vals, mode='lines', line_width=1,
                color=rest_mode_color, group='modes')
        # Add the 'num_modes' top modes as a series
        mode_x_vals = []
        mode_y_vals = []
        for pt in mode_pts[:num_modes]:
            mode_x_vals += [pt[0], pt[0], None]
            mode_y_vals += [0,     1,     None]
        cdf.add("Top %i Modes"%num_modes,mode_x_vals, mode_y_vals,
                mode='lines', line_width=1, color=top_mode_color, group='modes')


        # Generate the discrete CDF second derivative plot
        mode_pts = [list(pt) for pt in mode_pts]
        mode_pts.sort(key=lambda pt: pt[0])
        mode_pts = np.array(mode_pts)
        raw.add("Absolute Discrete CDF''", mode_pts[:,0], mode_pts[:,1],
                mode='lines', fill='tozeroy', color=top_mode_color,
                fill_color=rest_mode_color, line_width=1, group='modes')


        # Plot all three of these stacked together
        multiplot([[raw.plot(show=False, html=False)],
                   [cdf.plot(show=False, html=False)],
                   [pdf.plot(show=False, html=False)]
        ], shared_x=True)

    # Example / Test case of Principle Response Analyais
    if TEST_PRA:
        import numpy as np
        # Generate some points for testing.
        rgen = np.random.RandomState(10)
        n = 1000
        print("Number of points:", n)
        points = (rgen.rand(n,2) - .5) * 2

        # Create some testing functions (for learning different behaviors)
        funcs = [
            lambda x: x[0]               , # Linear on x
            lambda x: abs(x[0] + x[1])   , # "V" function on 1:1 diagonal
            lambda x: abs(2*x[0] + x[1]) , # "V" function on 2:1 diagonal
            lambda x: x[0]**2            , # Quadratic on x
            lambda x: x[1]**2            , # Quadratic on y
            lambda x: (x[0] + x[1])**2   , # Quadratic on 1:1 diagonal
            lambda x: (x[0] - x[1])**2   , # Quadratic on 1:-1 diagonal
            lambda x: (2*x[0] + x[1])**3 , # Cubic on 2:1 diagonal
            lambda x: rgen.rand()        , # Random function
        ]
        # Calculate the response values associated with each function.
        responses = np.vstack(tuple(tuple(map(f, points)) for f in funcs)).T

        # Reduce to just the first function
        choice = 4
        func = funcs[choice]
        response = responses[:,choice]

        print("Calculating PRA..")
        # Run the princinple response analysis function.
        components, lengths = pra(points, response, steps=100000, parallel=False)
        conditioner = components.copy()
        unconditioner = np.linalg.inv(components)
        # unconditioner = np.matmul(np.linalg.inv(components), np.diag(values))

        # Shorten the set of points and values for display purposes.
        size = 1000
        random_subset = np.arange(len(points))
        rgen.shuffle(random_subset)
        random_subset = random_subset[:size]
        points = points[random_subset]
        response = response[random_subset]

        print()
        print("Components")
        print(components)
        print()
        print("Conditioner")
        print(conditioner)
        print()

        # # ------------------------------------------------------------
        # # Generate a plot of the response surfaces.
        # from util.plot import Plot, multiplot
        # print("Generating plots of source function..")
        # # Add function 1
        # p1 = Plot()
        # p1.add("Points", *(points.T), response, opacity=.8)
        # p1.add_func("Surface", func, [-1,1], [-1,1], plot_points=100)
        # # Generate a regular approximation
        # from util.algorithms import NearestNeighbor as Approximator
        # model = Approximator()
        # model.fit(points, response)
        # p1.add_func("Unconditioned Approximation", model, [-1,1], [-1,1],
        #             mode="markers", opacity=.8)
        # # Generate a conditioned approximation
        # model = Approximator()
        # model.fit(np.matmul(points, conditioner), response)
        # approx = lambda x: model(np.matmul(x, conditioner))
        # p1.add_func("Best Approximation", approx, [-1,1], [-1,1],
        #             mode="markers", opacity=.8)
        # p1.plot(file_name="source_func.html", show=False)
        # # ------------------------------------------------------------

        # Return the set of (unique) vectors between points and the metric
        # change in response between each (unique) pair of points.
        def between(points, responses, metric=ABS_DIFF):
            import numpy as np
            vecs = []
            diffs = []
            for p1 in range(len(points)):
                for p2 in range(p1+1, len(points)):
                    vecs.append( points[p1] - points[p2] )
                    diffs.append( metric(responses[p1], responses[p2]) )
            return np.array(vecs, dtype=float), np.array(diffs, dtype=float)

        print("Generating demo plot 2")
        # Plot the between slopes to verify they are working.
        # Calculate the between slopes
        vecs, diffs = between(points, response)
        vec_lengths = np.sqrt(np.sum(vecs**2, axis=1))
        bs = (vecs.T * diffs / vec_lengths**2).T
        # Get a random subset of the between slopes and plot them.
        random_subset = np.arange(len(bs))
        rgen.shuffle(random_subset)
        random_subset = random_subset[:size]
        # Add function 1
        p2 = Plot()
        # Add the points, transformed points, and betweens for demonstration.
        new_pts = np.matmul(np.matmul(points, conditioner), unconditioner)
        p2.add("Original Points", *(points.T))
        # p2.add("Transformed Points", *(new_pts.T), color=p2.color(6, alpha=.7))
        p2.add("Between Slopes", *(bs[random_subset].T),
               color=p2.color(4, alpha=.4))
        # Add the principle response components 
        for i,vec in enumerate(components.T):
            vec *= 2
            p2.add(f"PC {i+1}", [0,vec[0]], [0,vec[1]], mode="lines")
            ax, ay = (vec / sum(vec**2)**.5) * 3
        p2.plot(file_name="source_func.html", 
                x_range=[-8,8], y_range=[-5,5])
        # ------------------------------------------------------------
        
