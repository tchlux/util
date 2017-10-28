import numpy as np
# from scipy.optimize import minimize
from scipy.interpolate import PchipInterpolator, splrep, splev
from scipy.integrate import quad as integrate
from scipy.stats import kstest
from util.algorithms import minimize


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
        
        
# Linearly interpolate to guess y values for x between provided data
def linear_fit_func(x_points, y_points, min_max=None):
    # Generate the 'range' of the function if it is not provided
    if type(min_max) == type(None):
        min_max = (min(x_points), max(x_points))

    # Fit the data with degree 1 splines (linear interpolation)
    fit = splrep(x_points, y_points, k=1)

    # Generate the function to return (gives the 'range' of
    # application when called without x values)
    def func(x_val=None, fit=fit, min_max=min_max):
        if type(x_val) == type(None):
            return min_max
        else:
            return splev(x_val, fit, ext=3)
    def deriv(x_val, fit=fit):
        return splev(x_val, fit, der=1, ext=3)

    # Make the attribute "derivative" return a function for the derivative
    setattr(func, "derivative", lambda: deriv)

    # Return the linear interpolation function
    return func


# Returns the CDF function value for any given x-value
def cdf_fit_func(data, cubic=False):
    # Sort the data and get the min and max
    data.sort()
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
        fit = splrep(data_pts[:,0], data_pts[:,1], k=1)
        fit = lambda x_val, fit=fit: splev(x_val, fit)
        fit.derivative = lambda d: (lambda x_val: splev(x_val, fit, der=d))
        # Generate inverse linear function
        inv_fit = splrep(data_pts[:,1], data_pts[:,0], k=1)
        inv_fit = lambda y_val, inv_fit=inv_fit: splev(y_val, inv_fit)
        inv_fit.derivative = lambda d: (lambda y_val: splev(y_val, inv_fit, der=d))
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
    # Generate the derivative of the function
    def deriv(x_val):
        val = fit.derivative(1)(x_val)
        val = np.where(x_val < min_pt, 0, val)
        val = np.where(x_val > max_pt, 0, val)
        return val

    # Set the derivative function
    setattr(cdf_func, "derivative", lambda: deriv)
    # Set the inverse function
    setattr(cdf_func, "inverse", lambda: fit.inverse)

    # Return the custom function for this set of points
    return cdf_func

# Returns the PDF function for the data
def pdf_fit_func(data=None):
    # Take the first derivative of the CDF function to get the PDF
    return cdf_fit_func(data).derivative(1)

# The normal pdf function at a given x
def normal_pdf(x, mean=0, stdev=1):
    return (1 / (np.sqrt(2*np.pi)*stdev)) * np.exp(
        - (x - mean)**2 / (2 * stdev**2) )

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
def ks_same_confidence(ks_stat, n1, n2=float('inf')):
    # By definition of the KS-test:
    # 
    #    c(a) = ( -ln(a/2)/2 )^(1/2)
    # 
    #   The standard test for if two distributions come from the same
    #   underlying distribution with "a" confidence. AKA, if you want
    #   the distributions to be the same, you want this test to pass
    #   with the largest "a" possible.
    # 
    #    KS > c(a) x (1/n1 + 1/n2)^(1/2)
    # 
    #   where n1 and n2 are the respective sample sizes for each
    #   distribution. The above check can be reversed to compute the
    #   largest "a" with which the KS test states the two
    #   distributions are certainly the same. (technically 1 above)
    # 
    #    c^-1(b) = 2 x e^(-2 x b^2)
    # 
    #    a = c^-1 ( (KS / (1/n1 + 1/n2)^(1/2)) )
    #      = e^( (KS / (1/n1 + 1/n2)^(1/2))^2 x -2 ) x 2
    #      = 2 x e^( -2 x KS^2 / |1/n1 + 1/n2| )
    # 
    return 2 * np.exp( -2 * ( ks_stat**2 / abs((1/n1) + (1/n2)) ))

def normal_confidence(distribution):
    # # Make the distribution 0 mean and unit variance (unit standard deviation)
    # new_distribution = (distribution - np.mean(distribution)) / np.var(distribution)
    # Compare the distribution with a normal distribution
    new_distribution = distribution
    ks_statistic = kstest(new_distribution, "norm").statistic
    return ks_same_confidence(ks_statistic, len(distribution))
    

# Calculate the maximum difference between two CDF functions (two sample)
def ks_diff(test_func, true_func):
    # Cycle through the functions to find the min and max of all ranges
    min_pt, max_pt = true_func()

    diff_func = lambda x: -abs(test_func(x) - true_func(x))
    # Use scipy minimize to find the greatest difference between the functions
    
    # sol = minimize(diff_func, [(max_pt - min_pt) / 2],
    #                bounds=[(min_pt,max_pt)], method='L-BFGS-B').x
    sol = minimize(diff_func, [(max_pt - min_pt) / 2],
                   bounds=[(min_pt,max_pt)])
    greatest_diff = abs(test_func(sol) - true_func(sol))
    # # Generate a large set of x-points (for a smooth plot)
    # x_points = np.linspace(min_pt, max_pt, 1000)
    # diff = abs(test_func(x_points) - true_func(x_points))

    # greatest_diff = -float('inf')
    # while (diff[np.argmax(diff)] > greatest_diff):
    #     lower = np.argmax(diff) - 1
    #     upper = np.argmax(diff) + 1
    #     min_pt = x_points[max(lower, 0)] - (
    #         1 if lower < 0 else 0)
    #     max_pt = x_points[min(upper,len(x_points)-1)] + (
    #         1 if upper >= len(x_points) else 0)
    #     x_points = np.linspace(min_pt, max_pt, 1000)
    #     diff = abs(test_func(x_points) - true_func(x_points))
    #     greatest_diff = max(max(diff), greatest_diff)

    return greatest_diff

# Subsample a set of data "subsamples" times with each sample of size
# "subsample_size" and then generate a list of fit-functions for each
# of the percentiles in "percentiles". Return list of fit functions.
def percentile_funcs(x_points, y_points, percentiles, min_max):
    # Generate the points for the percentile CDF's
    perc_points = []    
    for p in percentiles:
        perc_points.append(
            [np.percentile(y_points[:,i], p) for i in range(len(x_points))]
        )
    perc_points = np.array(perc_points)

    # Generate the CDF functions for each percentile
    funcs = []
    for pts in perc_points:
        # Generate a linear fit function over each set of percentile points
        func = linear_fit_func(x_points, pts)
        funcs.append(func)

    # Return the fit functions for each of the percentiles
    return funcs


# Given a plot, this will add the 'percentiles' cloud
def plot_percentiles(plot, name, funcs, percentiles, min_max_x,
                     center_color = np.array([255,70,0,0.6]),
                     outer_color = np.array([255,150,50,0.3])):
    group_id = str(np.random.randint(1000000))
    percentiles.sort() # Ensure they are sorted
    gaps = [percentiles[i] - percentiles[i-1] for i in range(1,len(percentiles))]
    textfont = dict(color='rgb(40,40,40)')
    for f,p,g in zip(funcs[:-1], percentiles[:-1], gaps):
        ratio = abs((50 - p - g/2)/50)
        color = center_color * abs(1.0 - ratio) + \
                outer_color  * ratio
        color = 'rgba(%i,%i,%i,%f)'%tuple(color)
        plot.add_func("%ith percentile"%p, f, min_max_x, fill='tonexty',
                      fill_color=color, mode='none', group=group_id,
                      show_in_legend=False, textfont=textfont)

    # Add the last function (bounding the set)
    plot.add_func("%ith percentile"%percentiles[-1], funcs[-1],
                  min_max_x, mode='lines', line_width=0,
                  group=group_id, show_in_legend=False, color=color,
                  textfont=textfont)

    # Add a master legend entry
    x_val = (min_max_x[1] - min_max_x[0])/2
    y_val = funcs[len(funcs)//2](x_val)
    plot.add(name + " Percentiles", [x_val], [y_val], mode='lines',
             color='rgba(%i,%i,%i,%f)'%tuple(center_color),
             group=group_id)

    return plot



if __name__ == "__main__":
    from util.plotly import Plot

    # p = Plot()
    # def c(a):
    #     return (-np.log(a/2)/2)**(1/2)
    # def f(n):
    #     return ks_same_confidence(.1, n)
    # def g(ks):
    #     return ks_same_confidence(ks, 1000)
    # p.add_func("c(a) changing a", c, [0.0001,.9999])
    # p.add_func("KS=.1 changing N", f, [10,1000])
    # p.add_func("N=1000 changing KS", g, [0,1])
    # p.plot(fixed=False)
    # exit()

    p = Plot()
    step_size = 1000
    steps = 10
    dist = []
    mean = 10.0
    std = 10.0
    for size in range(step_size,step_size*steps+1,step_size):
        dist += list(np.random.normal(mean, std, size=(step_size,)))
        # dist += list(np.random.random(size=(step_size,)))
        new_dist = np.array(dist)
        new_dist = (new_dist - mean) / std
        func = cdf_fit_func(dist, cubic=False)
        p.add_func("%i"%(len(dist)),func, [min(dist), max(dist)], group=str(len(dist)))
        print( ("%i "+"%.2f "*3)%
               (size, normal_confidence(dist), np.mean(dist), np.std(dist)) )
    p.plot()
    exit()


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
    
    exit()




