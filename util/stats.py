import numpy as np
from scipy.interpolate import splrep, splev
from scipy.interpolate import PchipInterpolator
from scipy.integrate import quad as integrate
from scipy.stats import kstest
import os

CWD = os.path.dirname(os.path.abspath(__file__))

# =====================================================
#      Principle Response and Principle Components     
# =====================================================

STANDARD_METRIC = lambda p1,p2: abs(p1-p2)

# Find a set of (nearly) orthogonal vectors that maximize the one norm
# of the vectors in "vec_iterator".
def ortho_basis(vec_iterator, num_bases=None, steps=float('inf'), parallel=False):
    # Determine the functions being used based on parallelism.
    if not parallel:
        from numpy import zeros
        from util.parallel import builtin_map as map
    else:
        from util.parallel import map, killall
        from util.parallel import shared_array as zeros
    # Get the first vector from the iterator and store as the first basis.
    for vec in vec_iterator: break
    # Store the "dimension"
    dim = len(vec)
    # Automatically set the number of bases if not provided.
    if type(num_bases) == type(None): num_bases = dim
    else:                             num_bases = min(dim, num_bases)
    bases_shape = (num_bases, len(vec))
    # Create global arrays that  safe for multiprocessing (without copies)
    if parallel: global _bases, _lengths, _counts
    _bases = zeros(bases_shape)
    _lengths = zeros((num_bases,))
    _counts = zeros((num_bases,))
    # Store the first basis vector (the vector that was pulled out).
    _bases[0,:] = vec
    _counts[0] = 1.0
    _lengths[0] = np.sqrt(np.sum(vec**2))
    # Given a vec, update the basis vectors iteratively according to
    # the provided vector (ignore race conditions).
    def update_bases(vec):
        print("STATS: Starting with a vector..", vec[0], len(vec), flush=True)
        for basis in range(num_bases):
            vec_length = np.sqrt(np.sum(vec**2))
            if (vec_length <= 0): return
            # Flip the vector if appropriate.
            if (np.dot(_bases[basis], vec) < 0): vec *= -1
            # Update the current basis vector (count, then vec, the length).
            _counts[basis] += 1
            _bases[basis] += (vec - _bases[basis]) / _counts[basis]
            _lengths[basis] += vec_length / _counts[basis]
            if (np.sum(_lengths[basis]**2) <= 1e-100): return
            # Remove this basis vector from "vec" (orthogonalize).
            shift = np.dot(vec, _bases[basis]) * (_bases[basis] / np.sum(_bases[basis]**2)) 
            vec -= shift
    print("STATS: Beginning the update process..", flush=True)
    # Perform rounds of updates using a (parallel) map operation.
    step = 1
    for _ in map(update_bases, vec_iterator):
        print("STATS: Step",_,"started..", flush=True)
        step += 1
        if step >= steps: break
    # Kill all hanging processes (because we may not have exausted the iterator).
    if parallel: killall()
    # Identify those bases that were actually computed, reduce to them.
    to_keep = _counts > 0
    bases = _bases[to_keep]
    lengths = _lengths[to_keep]
    # Normalize the bases and lengths (so that users can understand relative weight).
    bases = (bases.T / np.sqrt(np.sum(bases**2, axis=1))).T
    lengths /= np.sum(lengths)
    # Delete the temporary global variables (used for parallelism)
    if parallel: del _bases, _lengths, _counts
    return bases, lengths
        
# Generate "count" random indices of pairs that are within "length" bounds.
def gen_random_pairs(length, count=None, show_time=True, timeout=1):
    import time
    start = time.time()
    if type(count) == type(None):
        count = (length**2 - length) // 2
    # Iterate over random pairs.
    for _ in range(count):
        if (time.time() - start) >= timeout:
            print(f"{_: 7d} : {count: 7d}", end="\r")
            start = time.time()
        index_1 = np.random.randint(length)
        index_2 = np.random.randint(length-1)
        if (index_2 >= index_1): index_2 += 1
        yield index_1, index_2

# Generate vector between scaled by metric difference. Give the metric
# the indices of "vectors" in the provided matrix.
def gen_random_metric_diff(matrix, index_metric, count=None):
    # Iterate over random pairs.
    for (p1, p2) in gen_random_pairs(len(matrix), count):
        vec = matrix[p1] - matrix[p2]
        vec /= np.sqrt(np.sum(vec**2))
        yield index_metric(p1, p2) * vec

# Estimate the principle response directions using an iterative
# technique to save on memory usage. This is an approximation to the
# one-norm maximizers and may not converge.
def pra(points, responses, metric=STANDARD_METRIC,
        directions=100, steps=10000, parallel=False):
    if (len(points[0].shape) != 1): raise(Exception("Points must be an indexed set of vectors."))
    # Create a local "metric" based on the indices of row vectors.
    def index_metric(p1, p2):
        return metric(responses[p1], responses[p2])
    # Create a "vector generator" and pass it to the one-norm basis finder.
    vec_gen = gen_random_metric_diff(points, index_metric, steps)
    components, lengths = ortho_basis(vec_gen, num_bases=directions,
                                      steps=steps, parallel=parallel)
    # Currently not doing anyting with the "lengths", but they can be
    # returned too later if the code updates.
    return components, lengths


# Compute the principle components using sklearn.
def sklearn_pca(points):
    from sklearn.decomposition import PCA        
    pca = PCA()
    pca.fit(points)
    principle_components = pca.components_
    magnitudes = pca.singular_values_
    # Normalize the component magnitudes to have sum 1.
    magnitudes /= sum(magnitudes)
    return principle_components, magnitudes

# Compute the principle components of row points manually using NumPy.
#  (Eigenpairs of the covariance matrix)
def pca(points):
    points = points - np.mean(points, axis=0)
    # Compute the principle components of the (row) points.
    covariance_matrix = np.cov(points, rowvar=False)
    magnitudes, principle_components = np.linalg.eig(covariance_matrix)
    # Normalize the component magnitudes to have sum 1.
    magnitudes /= np.sum(magnitudes)
    # Get the ordering of the principle components (by decreasing value).
    order = np.argsort(-magnitudes)
    # Order the principle components by (decreasing) magnitude.
    magnitudes = magnitudes[order]
    principle_components = principle_components[order]
    # Return the components and their magnitudes.
    return principle_components, magnitudes



#      Functional Representations of Data     
# ============================================

# Given three points, this will solve the equation for the quadratic
# function which intercepts all 3
def solve_quadratic(x, y):
    if len(x) != len(y): raise(Exception("X and Y must be the same length."))
    if len(x) != 3: raise(Exception("Exactly 3 (x,y) coordinates must be given."))
    x1, x2, x3 = x
    y1, y2, y3 = y
    a = -((-x2 * y1 + x3 * y1 + x1 * y2 - x3 * y2 - x1 * y3 + x2 * y3)/((-x1 + x2) * (x2 - x3) * (-x1 + x3)))
    b = -(( x2**2 * y1 - x3**2 * y1 - x1**2 * y2 + x3**2 * y2 + x1**2 * y3 - x2**2 * y3)/((x1 - x2) * (x1 - x3) * (x2 - x3)))
    c = -((-x2**2 * x3 * y1 + x2 * x3**2 * y1 + x1**2 * x3 * y2 - x1 * x3**2 * y2 - x1**2 * x2 * y3 + x1 * x2**2 * y3)/((x1 - x2) * (x1 - x3) * (x2 - x3)))
    return (a,b,c)

# Returns the set of mode points given the CDF second difference list,
# mode points are sorted in order of magnitude of relative differences
def mode_points(second_diff, thresh_perc=0.0):
    # Collect together the maximum absolute values in each
    # group of signs (find local peaks)
    sign_maxes = []
    curr_max = second_diff[0]
    thresh = np.percentile(second_diff[:,1], thresh_perc)
    for val in second_diff:
        if (abs(val[1]) > thresh) and (val[1] * curr_max[1] < 0):
            sign_maxes.append(curr_max)
            curr_max = val
        elif abs(val[1]) > abs(curr_max[1]):
            curr_max = val
    sign_maxes = np.array(sign_maxes)
    # Record the magnitude of the difference between the peaks
    # as well as the location in the middle of the two peaks    
    sign_changes = []
    for i in range(1,len(sign_maxes)):
        x_val = (sign_maxes[i][0] + sign_maxes[i-1][0]) / 2
        y_val = abs(sign_maxes[i][1] - sign_maxes[i-1][1])
        sign_changes.append([x_val, y_val])
    # Sort the valley
    sign_changes.sort(key=lambda i: -i[1])
    sign_changes = np.array(sign_changes)
    shift = min(sign_changes[:,1])
    scale = max(sign_changes[:,1]) - shift
    sign_changes[:,1] = (sign_changes[:,1] - shift) / scale
    return sign_changes

# CDF Second Difference
def cdf_second_diff(data):
    # Sort the data and get the min and max
    data = list(data); data.sort()
    data = np.array(data)
    min_pt = data[0]
    max_pt = data[-1]
    # Initialize a holder for all CDF points
    data_pts = []
    # Add all of the CDF points for the data set
    for i,val in enumerate(data):
        if ((i > 0) and (val == data[i-1])): continue
        data_pts.append( [val, i/(len(data) - 1)] )
    # Add the 100 percentile point if it was not added already
    if data_pts[-1][1] != 1.0: data_pts[-1][1] = 1.0
    # Convert it all to numpy format for ease of processing
    data_pts = np.array(data_pts)
    second_diff_pts = []
    for i in range(1,len(data_pts)-1):
        a,_,_ = solve_quadratic(data_pts[i-1:i+2,0],
                                data_pts[i-1:i+2,1])
        second_deriv = 2*a
        second_diff_pts.append( [data_pts[i,0], second_deriv] )
    # # Sort the second_diff points by the magnitude of the second derivative
    # second_diff_pts.sort(key=lambda pt: -abs(pt[1]))
    return np.array(second_diff_pts)
        

# Returns the CDF function value for any given x-value
def cdf_fit_func(data, cubic=False):
    # Sort the data (without changing it) and get the min and max
    data = sorted(data)
    min_pt = data[0]
    max_pt = data[-1]

    # Initialize a holder for all CDF points
    data_pts = []

    # Add all of the CDF points for the data set
    for i,val in enumerate(data):
        if ((i > 0) and (val == data[i-1])): continue
        data_pts.append( [val, i/(len(data) - 1)] )
    # Add the 100 percentile point if it was not added already
    if data_pts[-1][1] != 1.0: data_pts[-1][1] = 1.0

    # Convert it all to numpy format for ease of processing
    data_pts = np.array(data_pts)

    # Generate a fit for the data points
    if not cubic:
        # Generate linear function
        fit_params = splrep(data_pts[:,0], data_pts[:,1], k=1)
        fit = lambda x_val: splev(x_val, fit_params, ext=3)
        fit.derivative = lambda x_val, d=1: splev(x_val, fit_params, der=d, ext=3)
        # Generate inverse linear function
        inv_fit_params = splrep(data_pts[:,1], data_pts[:,0], k=1)
        inv_fit = lambda y_val: splev(y_val, inv_fit_params, ext=3)
        inv_fit.derivative = lambda y_val, d=1: splev(y_val, inv_fit_params, der=d, ext=3)
        fit.inverse = inv_fit
    else:
        # Generate piecewise cubic monotone increasing spline
        fit = PchipInterpolator(data_pts[:,0], data_pts[:,1])
        fit.inverse = PchipInterpolator(data_pts[:,1], data_pts[:,0])
    
    # Generate a function that computes this CDF's points
    def cdf_func(x_val=None):
        if type(x_val) == type(None):
            return (min_pt, max_pt)
        else:
            y_val = fit(x_val)
            if (type(x_val) == np.ndarray):
                y_val = np.where(x_val < min_pt, 0.0, y_val)
                y_val = np.where(x_val > max_pt, 1.0, y_val)
            return y_val

    # Set the derivative function
    cdf_func.derivative = fit.derivative
    # Set the inverse function
    cdf_func.inverse = fit.inverse

    # Return the custom function for this set of points
    return cdf_func

# Returns the PDF function for the data
def pdf_fit_func(data=None):
    # Take the first derivative of the CDF function to get the PDF
    return cdf_fit_func(data).derivative(1)

# The normal pdf function at a given x
def normal_pdf(x, mean=0, stdev=1):
    return np.exp(-((x - mean)**2 / (2 * stdev**2)) ) / (np.sqrt(2*np.pi)*stdev)

# The normal cdf function at a given x
def normal_cdf(data, mean=0, stdev=1):
    try:
        len(data)
        return np.array([normal_cdf(row,mean,stdev) for row in data])
    except TypeError:
        return integrate(lambda x: normal_pdf(x,mean,stdev),
                         -np.inf, data)[0]

# Return the function of a normal cdf that has the mean and standard deviation
def normal_cdf_func(data):
    mean = np.mean(data)
    std = np.std(data)
    def cdf_func(x=None):
        if type(x) == type(None):
            return (min(data), max(data))
        return normal_cdf(x)
    return cdf_func

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

# Provides the probability that you would see a difference in CDFs as
# large as "ks" if the two underlying distributions are the same. This
# test makes no assumptions about where in the distribution the
# distance metric was found.
def same_prob(ks, n1, n2):
    # The probability of the two distributions taking on necessary values
    prob_val = lambda val,truth: (normal_pdf(val,truth,truth*(1-truth)) * 
                                  normal_pdf(val+ks,truth,truth*(1-truth)))
    prob_truth = lambda truth: integrate(lambda v: prob_val(v,truth), 0, 1-ks)[0]
    return 2 * (integrate(prob_truth, 0, 1)[0]) / (n1*n2)**(1/2)

def avg_same_prob(pts1, pts2):
    # get the EDF points for pts1 and pts2
    # get the differences at all known points
    # get the average probability that pts1 and pts2 are the same
    pass

def other(ks, n1, n2=float('inf')):
    # The probability of the two distributions taking on necessary values
    prob_val = lambda val,truth: (normal_pdf(val,truth,truth*(1-truth)) * 
                                  normal_pdf(val+ks,truth,truth*(1-truth)))
    prob_truth = lambda truth: integrate(lambda v: prob_val(v,truth), 0, 1-ks)[0]
    return 2 * prob_val(.5-(ks/2),.5) / np.sqrt(n1*n2)

def normal_confidence(distribution):
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

# Linearly interpolate to guess y values for x between provided data
def linear_fit_func(x_points, y_points):
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
    from util.plot import Plot
    import pickle, os

    TEST_CDF_PDF = False
    TEST_PRA = True

    if TEST_CDF_PDF:
        # num_distrs = 1000
        # for n in range(10,101,10):

        # p = Plot()
        # def c(a):
        #     return (-np.log(a/2)/2)**(1/2)
        # def f(n):
        #     return ks_p_value(.1, n)
        # def g(ks):
        #     return ks_p_value(ks, 1000)
        # p.add_func("c(a) changing a", c, [0.0001,.9999])
        # p.add_func("KS=.1 changing N", f, [10,1000])
        # p.add_func("N=1000 changing KS", g, [0,1])
        # p.plot(fixed=False)
        # exit()

        # p = Plot()
        # step_size = 1000
        # steps = 10
        # dist = []
        # mean = 10.0
        # std = 10.0
        # for size in range(step_size,step_size*steps+1,step_size):
        #     dist += list(np.random.normal(mean, std, size=(step_size,)))
        #     # dist += list(np.random.random(size=(step_size,)))
        #     new_dist = np.array(dist)
        #     new_dist = (new_dist - mean) / std
        #     func = cdf_fit_func(dist, cubic=False)
        #     p.add_func("%i"%(len(dist)),func, [min(dist), max(dist)], group=str(len(dist)))
        #     print( ("%i "+"%.2f "*3)%
        #            (size, normal_confidence(dist), np.mean(dist), np.std(dist)) )
        # p.plot()
        # exit()

        from util.plotly import Plot, multiplot

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
        components, lengths = pra(points, response, steps=100000, parallel=True)
        exit()
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
        def between(points, responses, metric=STANDARD_METRIC):
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
        
