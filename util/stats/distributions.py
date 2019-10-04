from util.stats import *

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
    import numpy as np
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
    import numpy as np
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
                index = np.searchsorted(data, x_val, side="right") - 1
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
