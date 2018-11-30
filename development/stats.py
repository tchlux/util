

# Compute the metric component analysis (one-norm maximizing between
# vectors scaled by metric change in response).
def ma(points, values, metric=ABS_DIFF, num_components=None, num_vecs=None):
    # Set the default values for the number of between vectors to use
    # and the number of components to keep from the between vectors.
    if (type(num_vecs) == type(None)):       num_vecs = 10 * points.shape[0]
    if (type(num_components) == type(None)): num_components = points.shape[1]
    # Function that takes two indices and returns metric difference.
    index_metric = lambda i1, i2: metric(values[i1], values[i2])
    # Generator that produces "between vectors".
    vec_gen = gen_random_metric_diff(points, index_metric, count=num_vecs)
    # Return the principle components of the metric slope vectors.
    return mas_vecs(np.array([v for v in vec_gen]), num_components=num_components)

# Compute the metric component analysis (one-norm maximizing between
# vectors scaled by metric change in response).
def mca(points, values, metric=ABS_DIFF, num_components=None, num_vecs=None):
    # Set the default values for the number of between vectors to use
    # and the number of components to keep from the between vectors.
    if (type(num_vecs) == type(None)):       num_vecs = 10 * points.shape[0]
    if (type(num_components) == type(None)): num_components = points.shape[1]
    # Function that takes two indices and returns metric difference.
    index_metric = lambda i1, i2: metric(values[i1], values[i2])
    # Generator that produces "between vectors".
    vec_gen = gen_random_metric_diff(points, index_metric, count=num_vecs)
    # Return the principle components of the metric slope vectors.
    return ca(np.array([v for v in vec_gen]), num_components=num_components)

# Compute the orthogonal vectors with the largest metric value.
def max_vals(points, values, metric=ABS_DIFF, num_components=None, num_vecs=None):
    # Set the default values for the number of between vectors to use
    # and the number of components to keep from the between vectors.
    if (type(num_vecs) == type(None)):       num_vecs = 10 * points.shape[0]
    if (type(num_components) == type(None)): num_components = points.shape[1]
    # Function that takes two indices and returns metric difference.
    index_metric = lambda i1, i2: metric(values[i1], values[i2])
    # Generator that produces "between vectors".
    vec_gen = gen_random_metric_diff(points, index_metric, count=num_vecs)
    # Return the principle components of the metric slope vectors.
    return max_vecs(np.array([v for v in vec_gen]), num_components=num_components)


# Look for the maximum length vectors that form an orthogonal basis.
def max_vecs(points, num_components=None, update=2):
    import time
    start = time.time()
    # Convert the points into their zero-mean form.
    print(" subtracting mean..", end="\r")
    points = points.copy() - np.mean(points, axis=0)
    # Set the dimension to be full automatically.
    if (type(num_components) == type(None)): num_components = points.shape[1]
    # Compue the components.
    print(" computing lengths..", end="\r")
    vecs = np.diag(np.ones(points.shape[1]))[:num_components]
    values = np.zeros(vecs.shape[0])
    orig_lens = np.sum(points**2, axis=1)
    # Get the first component.
    vecs[0,:] = points[np.argmax(orig_lens)]
    values[0] = np.linalg.norm(vecs[0,:]**2)
    vecs[0,:] /= values[0]
    # Find the rest of the components.
    for c in range(1, num_components):
        # Check to see if updates should be printed.
        if (type(update) != type(None)) and ((time.time() - start) > update):
            print(" computing component "+str(c+1)+" of "+str(num_components),end="\r")
            start = time.time()
        # Subtract out the previous vector from all vectors.
        for pt in points:
            pt -= np.dot(pt, vecs[c-1,:]) * vecs[c-1,:]
        # Save the new vector and normalize.
        lens = np.sum(points**2, axis=1)
        index = np.argmax(lens)
        vecs[c,:] = points[index]
        vecs[c,:] /= np.linalg.norm(vecs[c,:])
        values[c] = orig_lens[index]
    # Return the vectors and their values.
    return vecs, values

# Compute the component analysis (one-norm maximizing orthonormal vectors).
def ca(points, num_components=None, max_steps=100, max_seconds=60):
    import time
    print(" adjusting mean..", end="\r")
    # Convert the points into their zero-mean form.
    points = points.copy() - np.mean(points, axis=0)
    # Set the dimension to be full automatically.
    if (type(num_components) == type(None)): num_components = points.shape[1]
    # Compue the components.
    vecs = np.diag(np.ones(points.shape[1]))[:num_components]

    # Given a vector, return all points flipped onto the same side of the vector.
    def flipped_signs(vec):
        # Compute the signs and return them.
        return np.sign(np.matmul(points, vec))

    # Update each of the vectors to be the average of flipped points.
    def update_vecs():
        print("updating..", end=" ")
        for i in range(vecs.shape[0]):
            # Take the average of the appropriately flipped points
            vecs[i] = np.sum(points.T * flipped_signs(vecs[i]), axis=1)
            vec_len = np.linalg.norm(vecs[i], ord=2)
            if (vec_len > 0): vecs[i] /= vec_len

    # Get the value of each of the vectors.
    def value_vecs():
        values = []
        for v in vecs:
            value = np.sum(abs(np.matmul(v, points.T * flipped_signs(v))))
            values.append( value )
        return np.array(values)

    # Reorder the vectors by magnitude of of 1-norm sum (their values).
    def reorder_vecs():
        print("reordering..", end=" ")
        values = value_vecs()
        order = np.argsort(values)[::-1]
        vecs[:,:] = vecs[order][:,:]

    # Orthogonalize the vectors.
    def orthogonalize_vecs(flip=False):
        print("orthogonalizng..", end=" ")
        vecs[:,:] = (np.linalg.qr(vecs.T)[0].T)[:,:]

    # Flip vectors (to make them all consistent / unique solution)
    def flip_vecs():
        # Flip the vectors to be in an expected direction. (positive)
        signs = np.sign(np.matmul(vecs, np.ones(vecs.shape[1])))
        vecs[:,:] = (vecs.T * signs).T[:,:]

    # Start the iterations to find the components.
    previous_vecs = vecs.copy()
    start = time.time()
    for i in range(max_steps):
        print("\r update "+str(i+1)+"..", end=" ")
        update_vecs()
        reorder_vecs()
        orthogonalize_vecs()
        # Check for stagnation (if the updates didn't make a change).
        if np.linalg.norm(previous_vecs - vecs) <= 0: break
        else: previous_vecs = vecs.copy()
        # Break out if too much time was taken.
        if (time.time() - start) > max_seconds: break
    # Flip the vectors for consistency.
    flip_vecs()
    print(" "*70)
    # Return the vectors and their associated values (convexified).
    values = value_vecs()
    return vecs, values / sum(values)



# Estimate the principle response directions using an iterative
# technique to save on memory usage. This is an approximation to the
# one-norm maximizers and may not converge.
def pra(points, responses, metric=ABS_DIFF,
        directions=None, steps=None, parallel=False):
    if (len(points[0].shape) != 1): raise(Exception("Points must be an indexed set of vectors."))
    # Create a local "metric" based on the indices of row vectors.
    def index_metric(p1, p2):
        return metric(responses[p1], responses[p2])
    # Create a "vector generator" and pass it to the one-norm basis finder.
    vec_gen = gen_random_metric_diff(points, index_metric, steps)
    # Notice that this "vec_gen" has mean 0! That is important for ortho_basis.
    components, lengths = ortho_basis(vec_gen, num_bases=directions,
                                      steps=steps, parallel=parallel)
    # Currently not doing anyting with the "lengths", but they can be
    # returned too later if the code updates.
    return components, lengths

# Find a set of (nearly) orthogonal vectors that maximize the one norm
# of the vectors in "vec_iterator". The elements of "vec_iterator" are
# assumed to have a mean of zero!
def ortho_basis(vec_iterator, num_bases=None, steps=float('inf'), parallel=False):
    import fmodpy
    path_to_src = os.path.join(CWD,"fort_stats.f90")
    basis_update = fmodpy.fimport(path_to_src, output_directory=CWD).basis_update
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
    # Create global arrays that are safe for multiprocessing (without copies)
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
        for basis in range(num_bases):
            _counts[basis] += 1
            # Compute the basis update (with fortran code that's faster)
            _,_,vec_length = basis_update(_counts[basis], _bases[basis], vec)
            _lengths[basis] += vec_length / _counts[basis]
    # Perform rounds of updates using a (parallel) map operation.
    step = 1
    for _ in map(update_bases, vec_iterator):
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




# Compute the principle components of row points manually using NumPy.
#  (Eigenpairs of the covariance matrix)
def numpy_pca(points):
    # Transform the points so that they have zero-mean.
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
        

# Provides the probability that you would see a difference in CDFs as
# large as "ks" if the two underlying distributions are the same. This
# test makes no assumptions about where in the distribution the
# distance metric was found.
def same_prob(ks, n1, n2):
    from scipy.integrate import quad as integrate
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


# The normal pdf function at a given x
def normal_pdf(x, mean=0, stdev=1):
    return np.exp(-((x - mean)**2 / (2 * stdev**2)) ) / (np.sqrt(2*np.pi)*stdev)

# The normal cdf function at a given x
def normal_cdf(data, mean=0, stdev=1):
    from scipy.integrate import quad as integrate
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
