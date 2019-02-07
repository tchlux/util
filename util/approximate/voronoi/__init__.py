# Makes available all custom algorithms that I have worked with
import os, time
import numpy as np
import fmodpy
from util.approximate import WeightedApproximator

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))
DISPLAY_WAIT_SEC = 3.

class DuplicateInterpolationPoints(Exception): pass

# =================
#      Voronoi     
# =================
# Class for using the voronoi mesh to make approximations. This relies
# on some OpenMP parallelized fortran source to identify the voronoi
# mesh coefficients (more quickly than can be done in python).
class Voronoi(WeightedApproximator):
    points = None
    dots = None

    # Get fortran function for calculating mesh
    def __init__(self):
        # Get the source fortran code module
        path_to_src = os.path.join(CWD,"voronoi.f90")
        # Compile the fortran source (with OpenMP for acceleration)
        self.voronoi = fmodpy.fimport(path_to_src, output_directory=CWD,
                                      f_compiler_options=["-fPIC", "-O3","-fopenmp"],
                                      module_link_args=["-lgfortran","-fopenmp"])

    # Given points, pre-calculate the inner products of all pairs.
    def _fit(self, points):
        self.points = points.copy().T
        # Allocate space for pairwise dot products between points,
        # intialize all values to be the maximum representable number.
        self.dots = np.empty((self.points.shape[1],
                              self.points.shape[1]), order='F')
        self.voronoi.make_huge(self.dots)

    # Return the source points and weights associated with all the provided points.
    def _predict(self, points, display_wait_sec=DISPLAY_WAIT_SEC):
        indices = []
        weights = []
        idx = np.arange(self.points.shape[1])
        start = time.time()
        for i,pt in enumerate(points):
            # Update user if time has elapsed
            if ((time.time() - start) > display_wait_sec):
                start = time.time()
                print(f" {100.*i/len(points):.2f}%", end="\r", flush=True)
            # Calculate the support at this point
            _, support, error = self.voronoi.predict(self.points, self.dots, 
                                                     np.asarray(pt,order='F'))
            if (error != 0):
                raise(DuplicateInterpolationPoints("Some fit points were duplicated."))
            supported_points = idx[support > 0]
            supported_weights = support[supported_points]
            indices.append( supported_points )
            weights.append( supported_weights )
        return indices, weights


if __name__ == "__main__":
    from util.approximate.testing import test_plot
    model = Voronoi()
    p, x, y = test_plot(model, random=True, N=1000)
    p.show()
