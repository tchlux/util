# Makes available all custom algorithms that I have worked with
import os, time
import numpy as np
import fmodpy
from util.approximate import Approximator

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))

# =======================
#      LSHEP Wrapper     
# =======================

# Wrapper class for using the LSHEP fortran code
class LSHEP(Approximator):
    def __init__(self):
        self.linear_shepard = fmodpy.fimport(
            os.path.join(CWD,"linear_shepard.f95"),
            module_link_args=["-lblas","-llapack","-lgfortran"], 
            output_directory=CWD, autocompile_extra_files=True)
        self.lshep = self.linear_shepard.lshep
        self.lshepval = self.linear_shepard.lshepval
        self.ierrors = {}
        self.x = self.f = self.a = self.rw = None

    # Use fortran code to compute the boxes for the given data
    def _fit(self, control_points, values, **kwargs):
        # Local operations
        self.x = np.asfortranarray(control_points.T)
        self.f = []
        self.a = []
        self.rw = []
        # Fit a separate LSHEP model to each dimension of output.
        for i in range(values.shape[1]):
            self.f.append( np.asfortranarray(values[:,i]) )
            self.a.append( np.ones(shape=self.x.shape, order="F") )
            self.rw.append( np.ones(shape=(self.x.shape[1],), order="F") )
            # In-place update of self.a and self.rw
            self.lshep(self.x.shape[0], self.x.shape[1],
                       self.x, self.f[i], self.a[i], self.rw[i], **kwargs)

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
                resp = self.lshepval(x_pt, self.x.shape[0], self.x.shape[1], 
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

