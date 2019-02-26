import numpy as np

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

# Define a custom error for passing arrays of the wrong shape to pairwise.
class ShapeError(Exception): pass

# Compute the distance between pairs, return a matrix.
def pairwise_distance(x1, x2=None):
    # Define an error for having the wrong array shape.
    if (len(x1.shape) != 2): raise(ShapeError("Only 2D NumPy arrays are allowed."))
    # Determine use case.
    if (type(x2) == type(None)):
        # Compute the pairwise distances.
        x1_sq = np.sum(x1**2, axis=1, keepdims=True)
        d = (x1_sq + x1_sq[:,0]) - 2 * np.matmul(x1, x1.T)
        # Protect against the errors thta will occur along the diagonal.
        d[np.diag_indices(len(d))] = 1.0
        d[:,:] = np.sqrt(d[:,:])
        d[np.diag_indices(len(d))] = 0.0
        return d
    else:
        if (len(x2.shape) != 2): raise(ShapeError("Only 2D NumPy arrays are allowed."))
        # Compute the pairwise distance between memebers of each set.
        x1_sq = np.sum(x1**2, axis=1, keepdims=True)
        x2_sq = np.sum(x2**2, axis=1)
        return np.sqrt(x1_sq + x2_sq - 2 * np.matmul(x1, x2.T))


# ============================================================
#                          Test Cases     
# ============================================================

def _test_pairwise_distance(display=False):
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
        for j in range(i+1, len(x)):
            ij = dists[i,j]
            ji = dists[j,i]
            truth = np.linalg.norm(x[i] - x[j])
            assert(ij == ji)
            assert(abs(ij - truth) < 2**(-32))


    if display:
        from scipy.spatial.distance import pdist, cdist
        from util.system import Timer

        x = np.random.random((500,4000))
        y = np.random.random((5000,4000))


        t = Timer()
        d = pdist(x)
        t.stop()
        print("scipy  ", d.shape, t.total)

        t = Timer()
        d = pairwise_distance(x)
        t.stop()
        print("mine   ", d.shape, t.total)

        t = Timer()
        x2 = np.sum(x**2, axis=1, keepdims=True)
        d = np.sqrt((x2 + x2[:,0]) - 2 * np.matmul(x, x.T))
        t.stop()
        print("numpy  ", d.shape, t.total)

        print()

        t = Timer()
        d = cdist(x,y)
        t.stop()
        print("scipy  ", d.shape, t.total)

        t = Timer()
        d = pairwise_distance(x, y)
        t.stop()
        print("mine   ", d.shape, t.total)

        t = Timer()
        x2 = np.sum(x**2, axis=1, keepdims=True)
        y2 = np.sum(y**2, axis=1)
        p = 2 * np.matmul(x, y.T)
        d = np.sqrt(x2 + y2 - p)
        t.stop()
        print("numpy  ", d.shape, t.total)



if __name__ == "__main__":
    _test_pairwise_distance()

