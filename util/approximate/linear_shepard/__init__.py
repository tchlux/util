# Makes available all custom algorithms that I have worked with
import os, time
import numpy as np
import og_fmodpy as fmodpy
from util.approximate import Approximator

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))
# Import the source Fortran code.
lshep_mod = fmodpy.fimport(
    os.path.join(CWD,"linear_shepard.f95"),
    module_link_args=["-lblas","-llapack","-lgfortran"], 
    output_directory=CWD, autocompile_extra_files=True)


class TooFewDataPoints(Exception):
    def __init__(self, d, n):
        super().__init__(f"Only {n} points were provided but {d+2} are needed to define a modified shepard interpolant.")

# =======================
#      LSHEP Wrapper     
# =======================

# Wrapper class for using the LSHEP fortran code
class LSHEP(Approximator):
    #  radius -- Multiplier on the computed radius of influence. For
    #            data that is not uniformly dense, consider growing
    #            this number (growing the radius of influence of all
    #            points) until desired model smoothness is obtained.
    def __init__(self, radius=2):
        self.ierrors = {}
        self.radius = radius
        self.x = self.f = self.a = self.rw = None

    # Use fortran code to compute the boxes for the given data
    def _fit(self, control_points, values, **kwargs):
        # Local operations
        self.x = np.asfortranarray(control_points.T)
        # Raise an error if too few points are used for a fit.
        if (self.x.shape[1] <= self.x.shape[0]+1):
            raise(TooFewDataPoints(*self.x.shape))
        # Initialize containers.
        self.f = []
        self.a = []
        self.rw = []
        # Fit a separate LSHEP model to each dimension of output.
        for i in range(values.shape[1]):
            self.f.append( np.asfortranarray(values[:,i]) )
            self.a.append( np.ones(shape=self.x.shape, order="F") )
            self.rw.append( np.ones(shape=(self.x.shape[1],), order="F") )
            # In-place update of self.a and self.rw
            lshep_mod.lshep(self.x.shape[0], self.x.shape[1],
                            self.x, self.f[i], self.a[i], self.rw[i], **kwargs)
            # Multiply in the radius modifier.
            self.rw[-1] *= self.radius

    # Use fortran code to evaluate the computed boxes
    def _predict(self, x, display_wait_sec=.5):
        # Calculate all of the response values
        response = []
        start = time.time()
        for i,x_pt in enumerate(x):
            # Update user if time has elapsed
            if ((time.time() - start) > display_wait_sec):
                start = time.time()
                print(f" {100.*i/len(x):.2f}%", end="\r", flush=True)
            row = []
            for (f, a, rw) in zip(self.f, self.a, self.rw):
                x_pt = np.array(np.reshape(x_pt,(self.x.shape[0],)), order="F")
                ierr = 0
                resp = lshep_mod.lshepval(x_pt, self.x.shape[0], self.x.shape[1], 
                                          self.x, f, a, rw, ierr)
                self.ierrors[ierr] = self.ierrors.get(ierr, 0) + 1
                row.append(resp)
            response.append(row)
            # self.ierrors[ier] = self.ierrors.get(ier,0) + 1
        # Return the response values
        return np.array(response)

    # Print and return a summary of the errors experienced
    def errors(self):
        print("%i normal returns."%self.ierrors.get(0,0))
        print(("%i points outside the radius of influence "+
              "of all other points.")%self.ierrors.get(1,0))
        print(("%i points at which hyperplane normals are significantly"+
               " different.\n  This is only hazardous for piecewise "+
               "linear underlying functions.")%self.ierrors.get(2,0))
        return self.ierrors.get(0,0), self.ierrors.get(1,0), self.ierrors.get(2,0)


if __name__ == "__main__":
    from util.approximate.testing import test_plot
    m = LSHEP()
    p, x, y = test_plot(m, random=True, N=200)
    p.show()

