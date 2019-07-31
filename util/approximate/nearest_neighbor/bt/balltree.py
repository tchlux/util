# Ball Tree wrapper for Fortran module.

import os 
import fmodpy
import numpy as np
PATH_TO_HERE = os.path.dirname(os.path.abspath(__file__))
PATH_TO_BT   = os.path.join(PATH_TO_HERE, "ball_tree.f90")
btmod = fmodpy.fimport(PATH_TO_BT, output_directory=PATH_TO_HERE,
                       autocompile_extra_files=True)


class BallTree:
    # Given points and a leaf size, construct a ball tree.
    def __init__(self, points, leaf_size=1):
        self.leaf_size = leaf_size
        self.tree = np.asarray(points.T, order='F')
        self.sq_sums = np.zeros(self.tree.shape[1], order='F')
        self.radii = np.zeros(self.tree.shape[1], order='F')
        self.order = np.arange(self.tree.shape[1]) + 1
        # Build tree (in-place operation).
        #    .tree     will not be modified
        #    .sq_sums  will contain the squared sums of each point
        #    .radii    will be modified to have the radius of specific
        #              node, or 0.0 if this point is a leaf.
        #    .order    will be the list of indices (1-indexed) that
        #              determine the structure of the ball tree.
        btmod.build_tree(self.tree, self.sq_sums, self.radii,
                         self.order, leaf_size=self.leaf_size)
        self.index_mapping = self.order.copy()-1
        # Restructure the ball tree so the points are in locally
        # contiguous blocks of memory (local by branch + leaf).
        btmod.fix_order(self.tree, self.sq_sums, self.radii, self.order)

    # Find the "k" nearest neighbors to all points in z. Uses the same
    # interface as the "BallTree.query" function, see help for more info.
    def nearest(self, *args, **kwargs): return self.query(*args, **kwargs)

    # Find the "k" nearest neighbors to all points in z.
    def query(self, z, k=1, leaf_size=None, return_distance=True):
        if (len(z.shape) != 2) or (z.shape[1] != self.tree.shape[0]):
            class WrongDimension(Exception): pass
            raise(WrongDimension(f"The provided query points must be in a row-vector matrix with dimension (..,{self.tree.shape[0]})."))
        if (leaf_size is None): leaf_size = self.leaf_size
        k = min(k, self.tree.shape[1])
        points = np.asarray(z.T, order='F')
        indices = np.ones((k, points.shape[1]), order='F', dtype=int)
        dists = np.ones((k, points.shape[1]), order='F', dtype=float)
        # Compute the nearest neighbors.
        btmod.nearest(points, k, self.tree, self.sq_sums, self.radii,
                      self.order, leaf_size, indices, dists)
        # Return the results.
        if return_distance:
            return dists.T, indices.T - 1
        else:
            return indices.T - 1

