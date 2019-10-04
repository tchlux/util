from util.stats import *

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
        import numpy as np
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
