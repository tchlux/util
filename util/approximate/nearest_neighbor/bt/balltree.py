# Ball Tree wrapper for Fortran module.

# TODO: Add another operation mode for a 'byte' tree.
# TODO: Mix test cases for different type trees, make one testing suite.
# TODO: If data is 'uint8' then add 128 to the values before building tree.
# TODO: Add code for incrementally adding points to the tree


import os 
import fmodpy
import numpy as np
PATH_TO_HERE   = os.path.dirname(os.path.abspath(__file__))
PATH_TO_BT     = os.path.join(PATH_TO_HERE, "ball_tree.f90")
PATH_TO_INT_BT = os.path.join(PATH_TO_HERE, "int_ball_tree.f90")
btmod = fmodpy.fimport(PATH_TO_BT, output_directory=PATH_TO_HERE,
                       autocompile_extra_files=True)
int_btmod = fmodpy.fimport(PATH_TO_INT_BT, output_directory=PATH_TO_HERE,
                           autocompile_extra_files=True)


class BallTree:
    # Given points and a leaf size, construct a ball tree.
    def __init__(self, points, leaf_size=1, transpose=True):
        if transpose: points = points.T
        # Handle different use cases.
        if ('int8' in points.dtype.name):
            self.leaf_size = np.int32(leaf_size)
            self.tree    = np.asarray(points, order='F', dtype='int8')
            self.sq_sums = np.zeros(self.tree.shape[1],  dtype='int64')
            self.order   = np.arange(self.tree.shape[1], dtype='int32') + np.int32(1)
            self.build_tree = int_btmod.build_tree
            self.fix_order  = int_btmod.fix_order
            self.bt_nearest = int_btmod.nearest
            self.dtype = np.int8
            self.itype = np.int32
        else:
            self.leaf_size = leaf_size
            self.tree    = np.asarray(points, order='F', dtype='float64')
            self.sq_sums = np.zeros(self.tree.shape[1],  dtype='float64')
            self.order   = np.arange(self.tree.shape[1], dtype='int64') + 1
            self.build_tree = btmod.build_tree
            self.fix_order  = btmod.fix_order
            self.bt_nearest = btmod.nearest
            self.dtype = np.float64
            self.itype = np.int64
        self.radii   = np.zeros(self.tree.shape[1],  dtype='float64')
        # Build tree (in-place operation).
        #    .tree     will not be modified
        #    .sq_sums  will contain the squared sums of each point
        #    .radii    will be modified to have the radius of specific
        #              node, or 0.0 if this point is a leaf.
        #    .order    will be the list of indices (1-indexed) that
        #              determine the structure of the ball tree.
        self.build_tree(self.tree, self.sq_sums, self.radii,
                        self.order, leaf_size=self.leaf_size)
        self.index_mapping = self.order.copy()-1
        # Restructure the ball tree so the points are in locally
        # contiguous blocks of memory (local by branch + leaf).
        self.fix_order(self.tree, self.sq_sums, self.radii, self.order)

    # Find the "k" nearest neighbors to all points in z. Uses the same
    # interface as the "BallTree.query" function, see help for more info.
    def nearest(self, *args, **kwargs): return self.query(*args, **kwargs)

    # Find the "k" nearest neighbors to all points in z.
    def query(self, z, k=1, leaf_size=None, return_distance=True, transpose=True):
        if (len(z.shape) != 2) or (z.shape[1] != self.tree.shape[0]):
            class WrongDimension(Exception): pass
            raise(WrongDimension(f"The provided query points must be in a row-vector matrix with dimension (..,{self.tree.shape[0]})."))
        if (leaf_size is None): leaf_size = self.leaf_size
        if transpose: z = z.T
        k = min(k, self.tree.shape[1])
        # Initialize holders for output.
        points  = np.asarray(z, order='F', dtype=self.dtype)
        indices = np.ones((k, points.shape[1]), order='F', dtype=self.itype)
        dists   = np.ones((k, points.shape[1]), order='F', dtype='float64')
        # Compute the nearest neighbors.
        self.bt_nearest(points, k, self.tree, self.sq_sums, self.radii,
                        self.order, leaf_size, indices, dists)
        # Return the results.
        if return_distance:
            if transpose: return dists.T, indices.T - 1
            else:         return dists,   indices - 1
        else:
            if transpose: return indices.T - 1
            else:         return indices - 1
