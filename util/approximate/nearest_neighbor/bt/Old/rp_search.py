import numpy as np


class RPSearch:
    # Build this tree.
    def __init__(self, points, references=None):
        if (references is None): references = np.ceil(np.log2(len(points)))
        # Store the original points, reference indexes, and reference distances.
        self.references = references
        self.points = points
        self.ref_ids = np.zeros((references, len(points)), dtype=int)
        self.ref_dts = np.zeros((references, len(points)), dtype=float)
        # Pick (random) reference points.
        self.ref_ids[:,0] = np.random.choice(np.arange(points.shape[0]), references)
        self.ref_pts = points[self.ref_ids[:,0]]
        # Pre-compute the squared lengths of all points.
        self.sq_lens = np.sum(points**2, axis=1)
        self.ref_sq = self.sq_lens[self.ref_ids[:,0]]
        # Compute the pairwise distance between reference points and others.
        self.ref_dts[:,:] = (self.sq_lens[self.ref_ids[:,0],None] + self.sq_lens
                             - 2 * np.matmul(self.ref_pts, self.points.T))
        # Sort all reference distances.
        for r in range(references):
            sorted_ids = np.argsort(self.ref_dts[r,:])
            self.ref_ids[r,:] = sorted_ids
            self.ref_dts[r,:] = self.ref_dts[r,sorted_ids]
        # Set the first distances (correctly) to zero.
        self.ref_dts[:,0] = 0
        # Square-root the rest of the matrix to make distances correct.
        self.ref_dts[:,:] = np.sqrt(self.ref_dts[:,:])

    # Compute the nearest "k" neighbors to each query point.
    def query(self, pts, k=1, return_distance=True):
        d_calcs = 0
        # Make sure "pts" is 2D
        try:    len(pts[0])
        except: pts = [pts]
        k = min(k, len(self.points)-1)

        all_indices = []
        all_distances = []
        # Cycle through all provided points.
        for pt in pts:
            pt_sq = np.sum(pt**2)
            ref_dists = np.sqrt(pt_sq + self.ref_sq - 2*np.dot(self.ref_pts, pt))
            d_calcs += self.references
            closest = np.argmin(ref_dists)
            # Find the maximum possible "k" point distance for "pt".
            radius = ref_dists[closest] + self.ref_dts[closest, k+1]
            # Find the intersection of all "radius" bands.
            candidates = set(range(len(self.points)))
            for r in range(self.references):
                min_i = np.searchsorted(self.ref_dts[r], ref_dists[r]-radius, "left")
                max_i = np.searchsorted(self.ref_dts[r], ref_dists[r]+radius, "right")
                candidates.intersection_update( self.ref_ids[r,min_i:max_i] )
            print("ref_dists[closest]: ",ref_dists[closest])
            print("ref k dist:         ",self.ref_dts[closest, k+1])
            print("radius:             ",radius)
            print("num candidates:     ", len(candidates))
            candidates = np.array(list(candidates))
            d_calcs += len(candidates)
            indices = candidates
            distances = pt_sq + self.sq_lens[indices] - 2 * \
                        np.dot(self.points[candidates], pt)
            sorted_idxs = np.argsort(distances)[:k]
            all_indices.append( indices[sorted_idxs] )
            all_distances.append( distances[sorted_idxs] )
        # Convert all results into a NumPy array.
        indices = np.array(all_indices)
        print("Distances calculated:", d_calcs)
        # Return the appropriate info depending on usage.
        if return_distance:
            distances = np.array(all_distances)
            return (distances, indices)
        else:
            return indices
            


    from util.system import Timer
    t = Timer()
    dim = 50
    train = 100000
    test = 1
    k = 3
    x = np.random.random(size=(train,dim))
    z = np.random.random(size=(test,dim))
    print()
    print("x:", x.shape)
    print("z:", z.shape)
    print()
    print()

    print("Brute Force")
    t.start()
    d = np.sqrt(np.sum(x**2, axis=1) + np.sum(z[0]**2) - 2 * np.dot(x, z[0]))
    i = np.argsort(d)
    bt = t.stop()
    i = i[:k]
    d = d[i]
    print("Query time:", bt)
    print("d: ",d)
    print("i: ",i)
    print()

if __name__ == "__main__":
    from rp_search import RPSearch    
    print("RP Search")
    t.start()
    tree = RPSearch(x, references=40)
    ct = t.stop()
    print("Construction time:", ct)
    t.start()
    d, i = tree.query(z, k=k)
    qt = t.stop()
    print("Query time:       ", qt)
    print("d: ",d)
    print("i: ",i[0])
    print()

