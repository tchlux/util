
# Return a boolean "is_numeric"
def is_numeric(obj):
    try:
        abs((.3*obj + 1*obj) - .3*obj)
        return True
    except: return False

# Store the standard library hash function.
stdlib_hash = hash

# Generate hash values using a near-zero probability of conflict
# hashing mechanism that can handle non-hashable types. Does this
# by hashing the raw bytes (from pickle) with sha256.
def hash(something):
    try:
        return stdlib_hash(something)
    except:
        import hashlib, pickle
        return hashlib.sha256(pickle.dumps(something)).hexdigest()

# Generate a sorted unique set of values from an iterable, if the
# values are not hashable, use a hashing function to sort them.
def sorted_unique(iterable):
    try:
        # Try returning sorted unique values using hash to find unique.
        return sorted(set(iterable))
    except TypeError: 
        # Get the set of unique hashes and associated values.
        ids = {hash(v):v for v in iterable}
        # Try sorting the values directly, if that's not possible
        # then sort the values by their hash.
        try:    return sorted(ids.values())
        except: return [ids[h] for h in sorted(ids)]


# Given "d" categories that need to be converted into real space,
# generate a regular simplex in (d-1)-dimensional space. (all points
# are equally spaced from each other and the origin) This process
# guarantees that all points are not placed on a sub-dimensional
# manifold (as opposed to "one-hot" encoding).
def regular_simplex(num_categories):
    import numpy as np
    class InvalidNumberOfCategories(Exception): pass
    # Special cases for one and two categories
    if num_categories < 1:
        raise(InvalidNumberOfCategories(
            "Number of categories must be an integer greater than 0."))
    elif num_categories == 1:
        return np.array([[]])
    elif num_categories == 2:
        return np.array([[0.],[1.]])
    # Standard case for >2 categories
    d = num_categories
    # Initialize all points to be zeros.
    points = np.zeros((d,d-1))
    # Set the initial first point as 1
    points[0,0] = 1
    # Calculate all the intermediate points
    for i in range(1,d-1):
        # Set all points to be flipped from previous calculation while
        # maintaining the angle "arcos(-1/d)" between vectors.
        points[i:,i-1] = -1/(d-i) * points[i-1,i-1]
        # Compute the new coordinate using pythagorean theorem
        points[i,i] = (1 - sum(points[i,:i]**2))**(1/2)
    # Set the last coordinate of the last point as the negation of the previous
    points[i+1,i] = -points[i,i]
    # Return the regular simplex
    return points


# Given a point, determine the affine combination of categories in the
# defining regular simplex that produces the point.
def category_ratio(point):
    import numpy as np
    categories = regular_simplex(len(point)+1)
    # Add the "constant" column to the matrix of categories
    categories = np.hstack((categories, np.ones((len(point)+1,1))))
    # Add the "sum to 1" column to the affinely created point
    point = np.hstack((point, 1))
    # Calculate and return the affine weights
    return np.linalg.solve(categories.T, point)


# Return a randomized "range" using the appropriate technique based on
# the size of the range being generated. If the memory footprint is
# small (<= 32K int's) then a random sample is created and returned.
# If the memory footprint would be large, a Linear Congruential
# Generator is used to efficiently generate the sequence.
# 
# Parameters are similar to the builtin `range` with:
#   start -- int, default of 0.
#   stop  -- int > start for no / positive step, < start otherwise.
#   step  -- int (!= 0), default of 1.
#   count -- int > 0, number of samples, default (stop-start)//step.
# 
# Usage (and implied parameter ordering):
#   random_range(a)             --> range(0, a, 1)[:a]
#   random_range(a, b) [a < b]  --> range(a, b, 1)[:b-a]
#   random_range(a, b, c)       --> range(a, b, c)[:(b-a)//c]
#   random_range(a, b, c, d)    --> range(a, b, c)[:min(d, (b-a)//c)]
#   random_range(a, d) [d <= a] --> range(0, a, 1)[:d]
# 
# If the size of the range is large, a Linear Congruential Generator is used.
#   Memory  -- O(1) storage for a few integers, regardless of parameters.
#   Compute -- O(n) at most 2 times the number of steps in the range, n.
# 
def random_range(start, stop=None, step=None, count=float('inf')):
    from random import sample, randint
    from math import ceil, log2
    # Add special usage where the second argument is meant to be a count.
    if (stop != None) and (stop <= start) and ((step == None) or (step >= 0)):
        start, stop, count = 0, start, stop
    # Set a default values the same way "range" does.
    if (stop == None): start, stop = 0, start
    if (step == None): step = 1
    # Compute the number of numbers in this range, update count accordingly.
    num_steps = int((stop - start) // step)
    count = min(count, num_steps)
    # Check for a usage error.
    if (num_steps == 0) or (count <= 0): raise(InvalidRange(start, stop, step, count))
    # Use robust random method if it has a small enough memory footprint.
    # Do not use this method if there is a number too large for C implementation.
    if (num_steps <= 2**15) and (max(map(abs,(start,stop,step))) < 2**63):
        for value in sample(range(start,stop,step), count): yield value
        return
    # Use the LCG for the cases where the above is too memory intensive.
    # Enforce that all numbers are python integers (no worry of overflow).
    start, stop, step = map(int, (start, stop, step))
    # Use a mapping to convert a standard range into the desired range.
    mapping = lambda i: (i*step) + start
    # Seed range with a random integer to start.
    value = int(randint(0,num_steps))
    # 
    # Construct an offset, multiplier, and modulus for a linear
    # congruential generator. These generators are cyclic and
    # non-repeating when they maintain the properties:
    # 
    #   1) "modulus" and "offset" are relatively prime.
    #   2) ["multiplier" - 1] is divisible by all prime factors of "modulus".
    #   3) ["multiplier" - 1] is divisible by 4 if "modulus" is divisible by 4.
    # 
    offset = int(randint(0,num_steps)) * 2 + 1                 # Pick a random odd-valued offset.
    multiplier = 4*(num_steps + int(randint(0,num_steps))) + 1 # Pick a multiplier 1 greater than a multiple of 4.
    modulus = int(2**ceil(log2(num_steps)))                    # Pick a modulus just big enough to generate all numbers (power of 2).
    # Track how many random numbers have been returned.
    found = 0
    while found < count:
        # If this is a valid value, yield it in generator fashion.
        if value < num_steps:
            found += 1
            yield mapping(value)
        # Calculate the next value in the sequence.
        value = (value*multiplier + offset) % modulus


# Calculate the maximum difference between two CDF functions (two sample).
def ks_diff(l1, l2, n=100):
    import numpy as np
    # Sort both sets of points.
    l1 = sorted(l1)
    l2 = sorted(l2)
    # Get all unique points sorted.
    union = sorted(set(l1) | set(l2))
    # Compute the max difference between empirical distribution function values.
    max_diff = 0
    for v in union:
        p1_cdf = np.searchsorted(l1, v, side='right') / len(l1)
        p2_cdf = np.searchsorted(l2, v, side='right') / len(l2)
        max_diff = max(max_diff, abs(p1_cdf - p2_cdf))
    return max_diff


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
            low_2, upp_2, density_2 = next(gen_2)
        elif (upp_2 > upp_1):
            low_1, upp_1, density_1 = next(gen_1)
        else:
            low_1, upp_1, density_1 = next(gen_1)
            low_2, upp_2, density_2 = next(gen_2)
    return max(0, 1 - shared)


# Compute the effect between two sequences. Use the following:
#    number vs. number     -- Correlation coefficient between the two sequences.
#    category vs. number   -- "method" L1-norm difference between full distribution
#                             and conditional distributions for categories.
#    category vs. category -- "method" total difference between full distribution
#                             of other sequence given value of one sequence.
# 
# INPUTS:
# 
#   seq_1  -- An iterable of objects with a "__len__" property.
#   seq_2  -- An iterable of objects with a "__len__" property.
#   method -- Either "max", "min", "mean" or None. Determines how the
#             "effect" is computed for cat-num and cat-cat computations,
#             since each will yield multiple difference values.
#             (method==None) returns the full set of differences.
#   use_ks -- True or False, when measuring cat-num effect, if True, 
#             use the KS test p-value between conditional distributions
#             to compute effect, if Falsd use the e-pdf L1 difference.
# 
def effect(seq_1, seq_2, method="mean", use_ks=False):
    import numpy as np
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
            if (len(main_nums) == 0) or (len(sub_nums) == 0):
                dist_diff = 0.0
            elif (not use_ks):
                dist_diff = epdf_diff(main_nums, sub_nums)
            else:
                # Use the KS-statistic to test the likelihood that the two
                # sets come from different underlying distributions.
                dist_diff = 1.0 - ks_p_value(ks_diff(main_nums, sub_nums), 
                                             len(main_nums), len(sub_nums))
            try:              diffs[cat]       = dist_diff
            except TypeError: diffs[hash(cat)] = dist_diff
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
            if (len(main_seq) > 0) and (len(sub_seq) > 0):
                dist_diff = categorical_diff( main_seq, sub_seq )
            else: dist_diff = 0.
            try:              diffs[(None,cat)]       = dist_diff
            except TypeError: diffs[(None,hash(cat))] = dist_diff
    # Return the desired measure of difference between the two.
    if   (method == "max"):  return max(diffs.values())
    elif (method == "min"):  return min(diffs.values())
    elif (method == "mean"): return sum(diffs.values()) / len(diffs)
    else:                    return diffs

