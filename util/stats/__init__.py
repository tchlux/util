# Import all of the submodules so that users can access all from "util.stats".
from util.stats.rank import \
    insert_sorted, product, count_greater, rank_probability, \
    performance_ratio, performance_profile, data_profile

from util.stats.ks import \
    ks_diff, ks_p_value, normal_confidence

from util.stats.difference import \
    edf_pair_gen, epdf_diff, categorical_pdf, categorical_diff, effect

from util.stats.distributions import \
    Distribution, gauss_smooth, cdf_points, cdf_fit, pdf_fit

from util.stats.metric_pca import \
    gen_random_metric_diff, normalize_error, mpca, pca

from util.stats.plotting import \
    linear_fit_func, robust_percentile, percentile_funcs, \
    percentile_points, plot_percentiles


# EXPERIMENTAL
from util.stats.samples import \
    samples

from util.stats.modes import \
    modes


# BACKWARDS COMPATIBILITY
#   with warning for deprecation.

def cdf_fit_func(*args, **kwargs):
    print("\nWARNING: 'cdf_fit_func' is a deprecated function. Use 'cdf_fit' instead.\n")
    return cdf_fit(*args, **kwargs)

def pdf_fit_func(*args, **kwargs):
    print("\nWARNING: 'pdf_fit_func' is a deprecated function. Use 'pdf_fit' instead.\n")
    return pdf_fit(*args, **kwargs)


# TESTING
# ../../development/testing/test_stats.py 
if __name__ == "__main__":
    import numpy as np
    from util.random import cdf
    np.random.seed(1)
    data = cdf()
    
    from util.random import cdf
    np.random.seed(1)
    data = cdf(nodes=3).inverse(np.random.random(100))
    # modes(data)

    # Look for "largest" regions of separation 

