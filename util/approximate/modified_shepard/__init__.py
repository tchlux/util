# Makes available all custom algorithms that I have worked with
import os
import numpy as np
import fmodpy
from util.approximate import WeightedApproximator

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))

# Wrapper class for using the modified shepard fortran code
class ShepMod(WeightedApproximator):
    def __init__(self):
        self.modified_shepard = fmodpy.fimport(
            os.path.join(CWD,"modified_shepard.f95"),
            output_directory=CWD)
        self.shepmod = self.modified_shepard.shepmod
        self.shepmodval = self.modified_shepard.shepmodval
        self.ierrors = {}
        self.x = self.rw = None

    # Use fortran code to compute the boxes for the given data
    def _fit(self, control_points, **kwargs):
        # Local operations
        self.x = np.asfortranarray(control_points.T)
        self.rw = np.ones(shape=(self.x.shape[1],), order="F")
        # In-place update of self.rw (radius of influence).
        self.shepmod(self.x.shape[0], self.x.shape[1],
                     self.x, self.rw, **kwargs)

    # Use fortran code to evaluate the computed boxes
    def _predict(self, x):
        # Calculate all of the response values
        indices = []
        weights = []
        base = np.arange(self.x.shape[1])
        for x_pt in x:
            x_pt = np.array(np.reshape(x_pt,(self.x.shape[0],)), order="F")
            ierr = 0
            wts = np.zeros(self.x.shape[1], order="F")
            self.shepmodval(x_pt, *self.x.shape, self.x, self.rw, wts, ierr)
            # Get those indices where there is a nonzero weight.
            ids = base[wts > 0]
            self.ierrors[ierr] = self.ierrors.get(ierr, 0) + 1
            # Append the impactful indices and weights.
            indices.append( ids )
            weights.append( wts[ids] )
        # Return the source indices and weights associated with all
        # interpolation points.
        return indices, weights

    # Print and return a summary of the errors experienced
    def errors(self):
        print("%i normal returns."%self.ierrors.get(0,0))
        print(("%i points outside the radius of influence "+
              "of all other points.")%self.ierrors.get(1,0))
        return self.ierrors.get(0,0), self.ierrors.get(1,0)


if __name__ == "__main__":
    from util.approximate.testing import test_plot

    model = ShepMod()
    p, x, y = test_plot(model, random=True, N=20)
    p.show()
    print()
    model.errors()


    # # Test distribution prediction.
    # from util.plot import Plot
    # from util.stats import cdf_fit
    # np.random.seed(0)
    # x = np.random.random((100,4))
    # y = [cdf_fit(np.random.random((40,))) for _ in x]

    # model.fit(x,y)
    # test = np.random.random((20,4))
    # guesses = model(test)

    # p = Plot()
    # for i,f in enumerate(guesses):
    #     p.add_func(f"{i+1} CDF", f, f())
    # p.show()
