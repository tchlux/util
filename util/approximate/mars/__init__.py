from util.approximate import Approximator

# MARS
MARS_MAX_BASES = float('inf')
MARS_MAX_INTERACTION = float('inf')
MARS_MODEL_FLAG = 1
# ^^ 1 -> Piecewise linear model
#    2 -> Piecewise cubic model


# ==============================
#      MARS (Earth) Wrapper     
# ==============================

# Was using the fortran, but found a bug... Using Py-Earth instead

INSTALL_MESSAGE = """Failed to import `pyearth`. Download and install with:
  git clone git://github.com/scikit-learn-contrib/py-earth.git
  cd py-earth
  sudo python setup.py install --cythonize
"""

# Wrapper class for using the open source MARS, py-earth
class MARS(Approximator):
    def __init__(self, max_bases=MARS_MAX_BASES,
                 max_interaction=MARS_MAX_INTERACTION):
        # Disable the future warnings that MARS gives.
        import warnings
        warnings.filterwarnings("ignore", category=FutureWarning)
        # Warning about sklearn not working.
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        # Import the Earth package.
        try:    from pyearth import Earth
        except: raise(Exception(INSTALL_MESSAGE))
        self.Earth = Earth
        self.model = None
        self.max_bases = max_bases
        self.max_interaction = max_interaction

    # Use fortran code to compute the MARS model for given data
    #  See 'marspack.f' for details as to what these poorly named variables mean.
    def _fit(self, control_points, values, max_bases=None, max_interaction=None):
        if type(max_bases) != type(None):
            self.max_bases = max_bases
        if type(max_interaction) != type(None):
            self.max_interaction = max_interaction
        max_bases = min(control_points.shape[0],self.max_bases)
        max_interaction = min(control_points.shape[1],self.max_interaction)
        self.model = self.Earth(max_terms=max_bases,
                                max_degree=max_interaction,
                                use_fast=True)
        self.model.fit(control_points, values)

    # Use Earth to evaluate the computed boxes
    def _predict(self, x):
        response = self.model.predict(x)
        if len(response.shape) == 1: response = response[:,None]
        return response

if __name__ == "__main__":
    from util.approximate.testing import test_plot
    m = MARS()
    p, x, y = test_plot(m, random=True, N=20)
    p.show()
