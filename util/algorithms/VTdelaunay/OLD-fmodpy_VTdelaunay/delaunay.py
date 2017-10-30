# ============================
#      VTDelaunay Wrapper     
# ============================

# Wrapper class for using the Delaunay fortran code
class Delaunay(Approximator):
    from VTdelaunay import delaunayp
    def __init__(self):
        self.pts = None
        self.values = None
        self.errs = {}

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values):
        self.pts = np.asarray(control_points.T, dtype=float, order="F")
        self.values = np.asarray(values, dtype=float, order="F")

    # Use fortran code to evaluate the computed boxes
    def predict(self, x, debug=False, verbose=False):
        # Calculate all of the response values
        response = []
        for i,x_pt in enumerate(x):
            x_pt = np.asarray(x_pt, dtype=float, order="F")
            simp, weights, err, _ = Delaunay.delaunayp(self.pts.shape[0],
                                                       self.pts.shape[1],
                                                       self.pts, x_pt)
            # Assume execution went well, collect values
            resp = np.sum( self.values[simp-1] * weights )
            response.append(resp)
            self.errs[err] = self.errs.get(err,0) + 1
            # Print some debug statements
            if debug:
                simp = simp - 1
                if verbose:
                    print("Delaunay")
                    print("Simplex:", sorted(simp))
                    print("Training Points:", self.pts.shape)
                    # print(Delaunay.CLEAN_ARRAY_STRING(self.pts.T))
                    print("Test Point:")
                    print(Delaunay.CLEAN_ARRAY_STRING(x_pt))
                    print("Weights:")
                    print(Delaunay.CLEAN_ARRAY_STRING(weights))
                    print("Value:  ", resp)
                    print("Err:    ", err)
                return [simp]
        return response
