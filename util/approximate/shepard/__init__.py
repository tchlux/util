# Makes available all custom algorithms that I have worked with
import os, time
import numpy as np
import og_fmodpy as fmodpy
from util.approximate import WeightedApproximator

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))
DISPLAY_WAIT_SEC = 3.

# Get the source fortran code module
path_to_src = os.path.join(CWD,"weight.f90")
# Compile the fortran source (with OpenMP for acceleration)
shepmod = fmodpy.fimport(path_to_src, output_directory=CWD, omp=True)

# =================
#      Shepard     
# =================
# Class for using Shepard inverse distance weighting to make approximations. 
# This relies on some OpenMP parallelized fortran source to compute the
# weights to points (more quickly than can be done in python).
class Shepard(WeightedApproximator):
    points = None
    dots = None

    # Given points, pre-calculate the inner products of all pairs.
    def _fit(self, points):
        # Store the points provided at fit time.
        self.points = np.asarray(points.copy().T, order="F")

    # Return the source points and weights associated with all the provided points.
    def _predict(self, points, display_wait_sec=DISPLAY_WAIT_SEC):
        indices = np.zeros((points.shape[0], self.points.shape[1]), dtype=int)
        indices += np.arange(self.points.shape[1])
        weights = np.zeros(indices.shape)
        start = time.time()
        for i,pt in enumerate(points):
            # Update user if time has elapsed
            if ((time.time() - start) > display_wait_sec):
                start = time.time()
                print(f" {100.*i/len(points):.2f}%", end="\r", flush=True)
            # Calculate the support at this point
            weights[i,:] = shepmod.weight(self.points, pt.copy())
        return indices, weights

