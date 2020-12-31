import numpy as np
from util.approximate.base import Approximator, WeightedApproximator

# This class represents a Vandermonde interpolant. This type of
# predictor generalizes the idea of polynomial interpolation to
# arbitrary sets of basis functions. It is designed to use an
# intelligent method of picking interpolation points.
# 
# This is theoretically interesting, but for polynomials chosen to be
# the full set of [x, y, x^2, xy, y^2, ..] this is the same as the
# output of "polynomial_fit" and *MUCH* more computationally expensive.
# 
class Vandermonde(WeightedApproximator):
    def __init__(self, degree=2, fekete=True):
        self.degree = degree
        self.fekete = fekete

    # "Fit" the Vandermonde interpolant.
    def _fit(self, points, funcs=None):
        # Set the default values.
        if (funcs is None):
            from util.math.points import polynomials
            funcs = polynomials(self.degree, points.shape[1])

        # Pick the subset of points that should be kept.
        if self.fekete:
            from util.math.points import polynomials, fekete_indices
            # Make sure that we have enough points to interpolate.
            degree = self.degree
            while (len(funcs) > len(points)):
                degree -= 1
                funcs = polynomials(degree, points.shape[1])
            # Now compute the Vandermonde matrix to determine keepable points.
            vm = np.array([[f(p) for p in points] for f in funcs], order="f")
            to_keep = fekete_indices(vm)[:len(funcs)]
            # Update the set of points (and "y" values) that are kept.
            points = points[to_keep]
            self.kept = to_keep
            if (self.y is not None): self.y = self.y[to_keep]
        else:
            # Make the two sets equal length, take the minimum
            usable_length = min(len(points), len(funcs))
            keep = lambda array: np.array([array[i] for i in np.random.choice(
                np.arange(len(array)), size=usable_length, replace=False)])
            # If we are shorting functions, randomly select some to keep.
            if (len(funcs)  > usable_length): funcs = keep(funcs)
            # If we are shortening points, then we need to update the associated values.
            if (len(points) > usable_length):
                keep_idxs = keep(np.arange(len(points)))
                points = points[keep_idxs]
                if (self.y is not None): self.y = self.y[keep_idxs]


        # Store the points and basis functions.
        self.points = points
        self.basis_funcs = funcs
        # Get the Vandermonde matrix.
        self.vand_mat = np.array([[f(p) for f in self.basis_funcs] for p in points])
        # Get the determinant of the Vandermonde matrix.
        self.vand_mat_det = np.linalg.det(self.vand_mat)
        # Check to see if the point set is degenerate.
        if (abs(self.vand_mat_det) <= 2**(-26)):
            if (abs(self.vand_mat_det) <= 0): # 2**(-26)):
                class BadPoints(Exception): pass
                raise(BadPoints(f"The provided points and basis functions do not determine a polynomial. ({self.vand_mat_det})"))
            # If the determinant is not zero, then simply produce a warning.
            else:
                import warnings
                class BadPoints(Warning): pass
                warnings.warn(BadPoints("Vandermonde determinant near 0, unstable approximation."))

    # Given points, make a prediciton using the Vandermonde method.
    def _predict(self, points):
        pts = []
        wts = []
        for i,pt in enumerate(points):
            weights = []
            for j in range(len(self.vand_mat)):
                temp = self.vand_mat[j,:].copy()
                self.vand_mat[j,:] = [g(pt) for g in self.basis_funcs]
                det = np.linalg.det(self.vand_mat)
                self.vand_mat[j,:] = temp
                weights.append( det / self.vand_mat_det )
            pts.append( list(range(len(self.vand_mat))) )
            wts.append( weights )
        # Return the points and weights.
        return (pts, wts)


# Construct a multivariate polynomial.
class Polynomial(Approximator):
    def __init__(self, degree=3, fekete=False):
        self.degree = degree
        self.fekete = fekete

    # Fit a polynomial function.
    def _fit(self, points, values):
        from util.math.points import polynomials
        # Get the polymomials of desired degree for the given dimension.
        self.basis_funcs = polynomials(degree=self.degree, dimension=points.shape[1])
        # Get the Vandermonde matrix (for fitting the points).
        vand_mat = np.array([[f(p) for f in self.basis_funcs] for p in points])
        # Pick the subset of points that should be kept (if this should use Fekete).
        if self.fekete and (len(self.basis_funcs) < len(points)):
            # Determine which points should be kept.
            to_keep = fekete_indices(np.asfortranarray(vand_mat.T.copy()))
            to_keep = to_keep[:len(self.basis_funcs)]
            # Update the set of points (and "y" values) that are kept.
            vand_mat = vand_mat[to_keep]
            values = values[to_keep]
            self.kept = to_keep
        # Get the weights via a least-squares fit (using SVD, minimum norm for underdetermined).
        self.weights = np.linalg.lstsq(vand_mat, values, rcond=None)[0]

    # Evaluate the multivariate polynomial at the provided points.
    def _predict(self, points):
        # Get the matrix of all basis function evaluations at all points.
        evaluations = np.array([[f(x) for f in self.basis_funcs] for x in points])
        # Return the output values.
        return np.matmul(evaluations, self.weights)
# ====================================================================


if __name__ == "__main__":
    SHOW = True
    from util.math.points import fekete_points

    def _test_approximator(n=30, dim=2, degree=4, vandermonde=False,
                           polynomial=True, points_func=fekete_points,
                           show=SHOW):

        # Construct a reasonable test function to approximate.
        test_function = lambda x: 3*x[0]+.5*np.cos(8*x[0])+np.sin(5*x[-1])

        # Construct a Vandermonde interpolant
        if (points_func is None):
            np.random.seed(0)
            points = np.random.random(size=(n, dim))
        else:
            points = points_func(n, dim)

        # Generate values to be interpolated.
        values = np.array([test_function(x) for x in points])

        # Find the boundaries of points (to only plot interesting range).
        min_x = min(points[:, 0])
        max_x = max(points[:, 0])
        min_y = min(points[:, 1])
        max_y = max(points[:, 1])

        # Construct the visual.
        from util.plot import Plot
        p = Plot(f"Multivariate Interpolation (v={vandermonde}, p={polynomial})")
        p.add("Points", *(points.T), values)
        print("Adding truth..")
        p.add_func("Truth", test_function, [min_x,max_x], [min_y,max_y],
                   mode='markers', opacity=.5)

        if vandermonde:
            # Construct a Vandermonde fit.
            m = Vandermonde(degree=degree)
            m.fit(points, values)
            print("Adding Vandermonde")
            p.add_func(f"Vandermonde", m, [min_x,max_x], [min_y,max_y], plot_points=256)
            p.add("Kept", *(points[m.kept].T), values[m.kept], color=p.color(1))

        if polynomial:
            # Construct a polynomial fit.
            m = Polynomial(degree=degree)
            m.fit(points, values)
            print("Adding polynomial..", flush=True)
            p.add_func(f"Degree {degree} polynomial", m, [min_x,max_x], [min_y,max_y])
            if hasattr(m, "kept"):
                p.add("Kept", *(points[m.kept].T), values[m.kept], color=(0,50,220))

        if show: p.show(append=True)

    # Run tests on these two approximators.
    _test_approximator(vandermonde=False, polynomial=True, points_func=None)
    _test_approximator(vandermonde=True, polynomial=False, points_func=None)
