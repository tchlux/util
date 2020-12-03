# Makes available all custom algorithms that I have worked with
import os
import numpy as np
from util.approximate import WeightedApproximator
from util.math import SMALL

epsilon = 8 * SMALL # = 8 * SQRT(EPSILON(float))
ibudget = 10000 # 200000 
# IBUDGET - maximum number of steps to take through mesh, default 50K

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))


# =======================
#      qHullDelaunay     
# =======================

# Wrapper class for using the Delaunay scipy code
class qHullDelaunay(WeightedApproximator):
    def __init__(self):
        from scipy.spatial import Delaunay as Delaunay_via_qHull
        # For 1-D case, use regular interpolation
        from scipy import interpolate, optimize
        # Store modules in self for later access
        self.Delaunay_via_qHull = Delaunay_via_qHull
        self.interpolate = interpolate
        self.lstsq = optimize.lsq_linear
        self.pts = None

    # Use fortran code to compute the boxes for the given data
    def _fit(self, x):
        self.pts = x.copy()
        # Handle special case for 1-D data, use regular linear
        #  interpolation, otherwise use standard Delaunay
        if (self.pts.shape[1] > 1):
            self.surf = self.Delaunay_via_qHull(self.pts)
        else:
            self.sorted_pts = list(np.argsort(self.pts[:,0]))
            self.surf = None
            # self.surf = self.interpolate.interp1d(self.pts[:,0],
            #                                       self.y,
            #                                       fill_value="extrapolate")

    # Function that returns the indices of points and the weights that
    # should be used to make associated predictions for each point in
    # "points".
    def _predict(self, points):
        pts = []
        wts = []
        # Body
        for x in points:
            if hasattr(self.surf, "find_simplex"):
                # Solve for the weights in the old Delaunay model
                simp_ind = self.surf.find_simplex(x)
                # If a point is outside the convex hull, use the
                # closest simplex to extrapolate the value
                if simp_ind < 0: 
                    # Calculate the distance between the x each simplex
                    simp_dists = self.surf.plane_distance(x)
                    # Find the index of the closest simplex
                    simp_ind = np.argmax(simp_dists)
                # Solve for the response value
                simp = self.surf.simplices[simp_ind]
            else:
                # This is a 1-d problem
                simp = [np.argmin(np.sum((self.pts-x)**2,axis=1))]
                if (x[0] < self.pts[self.sorted_pts[0],0]):
                    simp += [self.sorted_pts[1]]
                    x = self.pts[self.sorted_pts[0],:]
                elif (x[0] > self.pts[self.sorted_pts[-1],0]):
                    simp += [self.sorted_pts[0]]
                    x = self.pts[self.sorted_pts[-1],:]
                elif (x[0] > self.pts[simp[0],0]):
                    simp += [self.sorted_pts[self.sorted_pts.index(simp[0])+1]]
                elif (x[0] <= self.pts[simp[0],0]):
                    simp += [self.sorted_pts[self.sorted_pts.index(simp[0])-1]]
                simp = np.array(simp)
            system = np.concatenate((self.pts[simp],
                np.ones((simp.shape[0],1))), axis=1).T
            x_pt = np.concatenate((x,[1]))
            try:
                weights = np.linalg.solve(system, x_pt)
            except np.linalg.linalg.LinAlgError:
                print(system)
                print(x_pt)
                print(simp)
                raise(Exception("Uh oh.."))
            # weights = np.where(weights > 0, weights, 0)
            # / sum(weights)
            pts.append( simp )
            wts.append( weights )
        # Convert into array form
        pts = np.array(pts)
        wts = np.array(wts)
        # Return the appropriate shaped pair of points and weights
        return (pts, wts)



# ==========================
#      Delaunay Wrapper     
# ==========================

# Wrapper class for using the Delaunay fortran code
class Delaunay(WeightedApproximator):
    # from multiprocessing import cpu_count
    # os.environ["OMP_NUM_THREADS"] = str(cpu_count())
    os.environ["OMP_NESTED"] = "TRUE"
    def __init__(self, parallel=False, pmode=None, chain=None):
        # from util.approximate.delaunay import delsparse
        import fmodpy
        self.delsparse = fmodpy.fimport("delsparse.f90", lapack=True, 
                                        omp=True, output_dir=CWD,
                                        f_compiler_args="-std=legacy -fPIC -shared -O3",
        )
        # Get the source fortran code module
        path_to_src = os.path.join(CWD,"delsparse.f90")
        # Set up the algorithm for parallel or serial evaluation.
        self.parallel = parallel
        self.delaunayp = self.delsparse.delaunaysparsep
        self.delaunays = self.delsparse.delaunaysparses
        self.pmode = pmode
        self.chain = chain
        # Initialize containers.
        self.pts = None
        self.errs = {}

    # Use fortran code to compute the boxes for the given data
    def _fit(self, control_points):
        self.pts = np.asarray(control_points.T, dtype=np.float64, order="F")

    # Return just the points and the weights associated with them for
    # creating the correct interpolation
    def _predict(self, x, allow_extrapolation=True, print_errors=True):
        # Get the predictions from VTdelaunay
        pts_in = np.asarray(self.pts.copy(), order="F")
        p_in = np.asarray(x.T, dtype=np.float64, order="F")
        simp_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                           dtype=np.int32, order="F")
        weights_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                              dtype=np.float64, order="F")
        error_out = np.ones(shape=(p_in.shape[1],), 
                            dtype=np.int32, order="F")
        if self.parallel:
            self.delaunayp(self.pts.shape[0], self.pts.shape[1],
                           pts_in, p_in.shape[1], p_in, simp_out,
                           weights_out, error_out, extrap=100.0,
                           pmode=self.pmode, ibudget=ibudget,
                           eps=epsilon, chain=self.chain)
        else:
            self.delaunays(self.pts.shape[0], self.pts.shape[1],
                           pts_in, p_in.shape[1], p_in, simp_out,
                           weights_out, error_out, extrap=100.0, 
                           ibudget=ibudget, eps=epsilon, chain=self.chain)
        # Remove "extrapolation" errors if the user doesn't care.
        if allow_extrapolation: error_out = np.where(error_out == 1, 0, error_out)
        else:
            if 1 in error_out:
                class Extrapolation(Exception): pass
                raise(Extrapolation("Encountered extrapolation point when making Delaunay prediction."))
        # Handle any errors that may have occurred.
        if (sum(error_out) != 0):
            if print_errors:
                unique_errors = sorted(np.unique(error_out))
                print(" [Delaunay errors:",end="")
                for e in unique_errors:
                    if (e == 0): continue
                    print(" %3i"%e,"at","{"+",".join(tuple(
                        str(i) for i in range(len(error_out))
                        if (error_out[i] == e)))+"}", end=";")
                print("] ")
            # Reset the errors to simplex of 1s (to be 0) and weights of 0s.
            bad_indices = (error_out > (1 if allow_extrapolation else 0))
            simp_out[:,bad_indices] = 1
            weights_out[:,bad_indices] = 0
        # Adjust the output simplices and weights to be expected shape
        points  = simp_out.T - 1
        weights = weights_out.T
        # Return the appropriate shaped pair of points and weights
        return (points, weights)
        

# Wrapper class for using the Delaunay fortran code
class DelaunayP1(Delaunay):
    def __init__(self, parallel=True, pmode=1):
        return super().__init__(parallel=parallel, pmode=pmode)

# Wrapper class for using the Delaunay fortran code
class DelaunayP2(Delaunay):
    def __init__(self, parallel=True, pmode=2):
        return super().__init__(parallel=parallel, pmode=pmode)

# Wrapper class for using the Delaunay fortran code
class DelaunayP3(Delaunay):
    def __init__(self, parallel=True, pmode=3):
        return super().__init__(parallel=parallel, pmode=pmode)

# Wrapper class for using the Delaunay fortran code
class DelaunayP3CT(Delaunay):
    def __init__(self, parallel=True, pmode=3,chain=True):
        return super().__init__(parallel=parallel, pmode=pmode, chain=chain)


if __name__ == "__main__":
    from util.approximate.testing import test_plot
    m = Delaunay()
    p, x, y = test_plot(m, random=True, N=200)
    p.show()
