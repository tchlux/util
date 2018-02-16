import numpy as np
from scipy.interpolate import PchipInterpolator, splrep, splev

DEFAULT_ORDER = 1

def spline_fit(points, values, order=1):
    # Generate linear function
    fit = splrep(points, values, k=order)
    fit = lambda x_val, fit=fit: splev(x_val, fit)
    fit.derivative = lambda d: (lambda x_val: splev(x_val, fit, der=d))
    return fit

# Class for computing an interpolation between the nearest n neighbors
class SplineMesh:
    def __init__(self):
        self.points = None
        self.values = None
        self.num_neighbors = None

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values, order=DEFAULT_ORDER):
        # Process and store local information
        self.points = control_points.copy()
        self.sorts = np.array([np.argsort(self.points[:,d].copy())
                               for d in range(self.points.shape[1])])
        print()
        print(self.points[:,0][self.sorts[0]])
        print(values[self.sorts[0]])
        print()

        self.fits = [spline_fit(self.points[:,d][sort], values[sort], order)
                     for d,sort in enumerate(self.sorts)]

        from plotly_interface import Plot
        p = Plot()
        p.add("Points", self.points[:,0], values)
        min_val = min(self.points[:,0])
        max_val = max(self.points[:,1])
        p.add_func("Fit", self.fits[0], [min_val, max_val])
        p.plot()
        exit()

        coefficients = np.array(
            [[self.fits[d](pt[d]) for d in range(self.points.shape[1])]
             for pt in self.points])
        self.weights, self.residuals = np.linalg.lstsq(coefficients, values)[:2]

    # Use fortran code to evaluate the computed boxes
    def __call__(self, x):
        if len(x.shape) == 1:
            x = np.array([x])
        if len(x.shape) != 2:
            raise(Exception("ERROR: Bad input shape."))
        # Use the nearest point 
        response = []
        for pt in x:
            value = self.fits[0](pt[0])
            # value = np.dot([self.fits[d](pt[d]) for d in range(len(pt))],self.weights)
            response.append(value)

        response = np.array(response)
        return response[0] if len(response) == 1 else response

