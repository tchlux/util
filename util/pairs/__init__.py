
# This function maps an index in the range [0, count*(count - 1) // 2] 
# to a tuple of integers in the range [0,count). The mapping will cover
# all pairs if you use all indices between [0, count*(count - 1) // 2].
def pair_to_index(p1, p2):
    if (p1 < p2): p1, p2 = p2, p1
    return (p1 * (p1 - 1) // 2) + p2

# This function maps an index in the range [0, count*(count - 1) // 2] 
# to a tuple of integers in the range [0,count). The mapping will cover
# all pairs if you use all indices between [0, count*(count - 1) // 2].
def index_to_pair(index):
    val = int(((1/4 + 2*index)**(1/2) + 1/2))
    remainder = index - val*(val - 1)//2
    return (val, remainder)

# Compute the number of original elements when given the number of pairs.
# This is done by realizing the following:
# 
#      n(n-1) / 2 = f(n)
#   (n^2 - n) / 2 = f(n)
#         n^2 - n = 2 f(n)
#   n^2 - n + 1/4 = 2 f(n) + 1/4
#     (n - 1/2)^2 = 2 f(n) + 1/4
#         n - 1/2 = (2 f(n) + 1/4)^{1/2}
#               n = (2 f(n) + 1/4)^{1/2} + 1/2
#              2n = (8 f(n) + 1)^{1/2} + 1
#               n = ((8 f(n) + 1)^{1/2} + 1) / 2
# 
def num_from_pairs(num_pairs):
    return int(round( ((8*num_pairs + 1)**(1/2) + 1) / 2 ))

# Compute the distance between pairs, return in something that can be
# indexed like a matrix.
def pairwise_distance(x1, x2=None):
    import numpy as np
    # Define a subclass of a NumPy array that converts pair indices to vector indices.
    class Pairwise(np.ndarray):
        # Overwrite only the "[]" operator to handle pair indices.
        def __getitem__(self, index):
            try: return super().__getitem__(index)
            except:
                if (type(index) != tuple):
                    raise(IndexError("Only 1D and 2D indexing are supported."))
                elif (len(index) != 2):
                    raise(IndexError("Only 1D and 2D indexing are supported."))
                # Handle integer indices (retreival).
                if (type(index[0]) == int) and (type(index[1]) == int):
                    i1, i2 = index
                    if (i1 == i2): return 0.0
                    elif (i2 < i1): i1, i2 = i2, i1
                    if (i1 < 0) or (i2 < 0): raise(IndexError("Negative indexing is not supported."))
                    if (i1 > self.shape[0]) or (i2 >= self.shape[0]): raise(IndexError("Index out of range."))
                    return self[pair_to_index(i1, i2)]
                # Convert slice indices into lists.
                if (type(index[0]) == slice):
                    index = (list(range(len(x1))[index[0]]), index[1])
                if (type(index[1]) == slice):
                    index = (index[0], list(range(len(x1))[index[1]]))
                # Convert iterable indices into lists.
                is_iterable = lambda x: hasattr(x, "__iter__")
                if is_iterable(index[0]) and (type(index[0]) != list):
                    index = (list(index[0]), index[1])
                if is_iterable(index[1]) and (type(index[1]) != list):
                    index = (index[0], list(index[1]))
                # Verify the types of the indices.
                ok_types = {int, list}
                if (type(index[0]) not in ok_types): raise(IndexError(f"Unrecognized index at 0, '{type(index[0])}'."))
                if (type(index[1]) not in ok_types): raise(IndexError(f"Unrecognized index at 1, '{type(index[1])}'."))
                # Handle the extraction of a single row of distances.
                if (type(index[0]) == int):
                    return np.array([ self[pair_to_index(index[0], i)]
                                      if (index[0] != i) else 0.0
                                      for i in index[1] ])
                elif (type(index[1]) == int):
                    return np.array([ self[pair_to_index(i, index[1])]
                                      if (i != index[1]) else 0.0
                                      for i in index[0] ])
                # Handle the creation of a new pairwise object.
                interest_ids = sorted(set(index[0] + index[1]))
                if not (len(interest_ids) == len(set(index[0])) == len(set(index[1]))):
                    raise(IndexError("Only square selections can be made from this object."))
                if (len(interest_ids) < 1): raise(IndexError("Empty selection of indices."))
                if (len(interest_ids) == 1): return np.zeros(1)
                num = len(interest_ids)
                total = (num*(num - 1)) // 2
                dists = Pairwise((total,), dtype=float, order="F")
                for i in range(total):
                    i1, i2 = index_to_pair(i)
                    i1, i2 = interest_ids[i1], interest_ids[i2]
                    dists[i] = self[pair_to_index(i1, i2)]
                return dists
        # Overwrite the "shape" property so that it behaves in
        # correspondence with the indexing scheme.
        @property
        def shape(self): return (num_from_pairs(len(self)),) * 2


    # Import Fortran code for quickly computing pairwise distance.
    import os
    cwd = os.path.dirname(os.path.abspath(__file__))
    import fmodpy
    pairwise = fmodpy.fimport(os.path.join(cwd, "pairwise.f08"), omp=True, output_directory=cwd)
    # Define an error for having the wrong array shape.
    class ShapeError(Exception): pass
    if (len(x1.shape) != 2): raise(ShapeError("Only 2D NumPy arrays are allowed."))

    # Determine use case.
    if (type(x2) == type(None)):
        # Compute the pairwise distances.
        dists = Pairwise(((len(x1)*(len(x1) - 1)) // 2,), dtype=float, order="F")
        pairwise.full(np.asfortranarray(x1), dists)
        return dists
    else:
        if (len(x2.shape) != 2): raise(ShapeError("Only 2D NumPy arrays are allowed."))
        # Compute the pairwise distance between memebers of each set.
        return pairwise.half(np.asfortranarray(x1), np.asfortranarray(x2))\
                       .reshape((len(x1), len(x2)))


# ============================================================
#                          Test Cases     
# ============================================================

def _test_pairwise_distance(display=True):
    import numpy as np
    x = np.random.random((100,30))
    y = np.random.random((10,30))

    # Verify that the pairwise distance (between points in two sets) works.
    dists = pairwise_distance(x, y)
    for i in range(len(x)):
        for j in range(len(y)):
            truth = np.linalg.norm(x[i] - y[j])
            assert(abs(dists[i,j] - truth) < 2**(-32))

    # Verify that the pairwise distance (between points in one set) works.
    dists = pairwise_distance(x)
    for i in range(len(x)):
        for j in range(i, len(x)):
            ij = dists[i,j]
            ji = dists[j,i]
            truth = np.linalg.norm(x[i] - x[j])
            assert(ij == ji)
            assert(abs(ij - truth) < 2**(-32))

    # Verify row-extraction for pairwise distances.
    sub_dists = dists[0,:]
    for i in range(len(sub_dists)):
        assert(sub_dists[i] == dists[0,i])
    sub_dists = dists[:,1]
    for i in range(len(sub_dists)):
        assert(sub_dists[i] == dists[i,1])
    sub_dists = dists[:20:2,:20:2]
    for i_sub,i in enumerate(range(0,20,2)):
        for j_sub,j in enumerate(range(0,20,2)):
            assert(sub_dists[i_sub, j_sub] == dists[i,j])
    assert(sub_dists.shape[0] == len(range(0,20,2)))
    # Try a bad list index.
    try: sub_dists = dists[[1,2,3],[2,3,3]]
    except IndexError: pass
    else: assert(False)

if __name__ == "__main__":
    _test_pairwise_distance()

