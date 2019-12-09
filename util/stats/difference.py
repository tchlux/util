

# ====================================================
#      Measuring the Difference Between Sequences     
# ====================================================

# Generator for returning pairs of EDF points as (val1, val2, density).
def edf_pair_gen(seq):
    # Compute the 'extra' slope lost at the beginning of the function.
    for i in range(len(seq)):
        if seq[i] != seq[0]: break
    extra = (i/len(seq))
    shift = lambda v: (v - extra) / (1 - extra)
    # Manually set the first slope value.
    x1, y1 = -float('inf'), 0
    for i in range(len(seq)):
        # Cycle until we hit the last occurrence of the next value.
        if ((i+1) < len(seq)) and (seq[i+1] == seq[i]): continue
        x2, y2 = seq[i], shift( (i+1)/len(seq) )
        yield (x1, x2, y2 - y1)
        # Cycle values.
        x1, y1 = x2, y2
    # Manually yield the last element in the sequence.
    yield (x1, float('inf'), 0)

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

# Compute the PDF of a sequence of categories.
def categorical_pdf(seq, unique=None):
    from util.system import hash, sorted_unique
    if (type(unique) == type(None)): unique = sorted_unique(seq)
    try:              return {c:sum(v==c for v in seq) / len(seq) for c in unique}
    except TypeError: return {hash(c): sum(v==c for v in seq) / len(seq) for c in unique}

# Given the PDF (in a dictionary object) for two different sets of
# categorical values, compute the total change in the PDF sum(abs(diffs)).
def categorical_diff(pdf_1, pdf_2):
    # Check for equally lengthed sequences.
    class CannotCompare(Exception): pass
    def near_one(val): return (abs(val - 1) < 1e-5)
    if not near_one(sum(pdf_1.values())):
        raise(CannotCompare(f"Sum of PDF values for sequence one '{sum(pdf_1.values()):.2f}' is not 1."))
    if not near_one(sum(pdf_2.values())):
        raise(CannotCompare(f"Sum of PDF values for sequence one '{sum(pdf_2.values()):.2f}' is not 1."))
    # Sum the differences in categories.
    all_values = sorted(set(pdf_1).union(set(pdf_2)))
    return sum(abs(pdf_1.get(v,0.) - pdf_2.get(v,0.)) for v in all_values) / 2

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
    from util.system import hash, sorted_unique
    from util.math import is_numeric
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
                from util.stats import cdf_fit, ks_diff, ks_p_value
                dist_diff = 1.0 - ks_p_value(
                    ks_diff(cdf_fit(main_nums), cdf_fit(sub_nums)),
                    len(main_nums), len(sub_nums))
            try:              diffs[cat]       = dist_diff
            except TypeError: diffs[hash(cat)] = dist_diff
        if   (method == "max"):  return max(diffs.values())
        elif (method == "min"):  return min(diffs.values())
        elif (method == "mean"): return sum(diffs.values()) / len(diffs)
        else:                    return diffs
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

