import numpy as np
from util.math import abs_diff

# =============================================
#      Metric Principle Component Analysis     
# =============================================

# Generate vector between scaled by metric difference. Give the metric
# the indices of "vectors" in the provided matrix.
def gen_random_metric_diff(matrix, index_metric, power=2, count=None):
    from util.random import pairs
    # Iterate over random pairs, skip those with no difference.
    for (p1, p2) in pairs(len(matrix), count):
        metric_diff = index_metric(p1, p2)
        if (metric_diff <= 0): continue
        vec = matrix[p1] - matrix[p2]
        length = np.linalg.norm(vec)**power
        if length > 0: yield metric_diff * vec / length

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
    # Make sure the first component starts with a positive value (consistency).
    if (components[0,0] < 0): components *= -1
    # Return the principle components of the metric slope vectors.
    return components, weights

# Compute the principle components using sklearn.
def pca(points, num_components=None, display=True):
    from sklearn.decomposition import PCA        
    pca = PCA(n_components=num_components)
    if (num_components is None): num_components = min(*points.shape)
    else: num_components = min(num_components, *points.shape)
    if display: print(f"Computing {num_components} principle components..",end="\r", flush=True)
    pca.fit(points)
    if display: print( "                                                          ",end="\r", flush=True)
    principle_components = pca.components_
    magnitudes = pca.singular_values_
    # Normalize the component magnitudes to have sum 1.
    magnitudes /= np.sum(magnitudes)
    return principle_components, magnitudes

