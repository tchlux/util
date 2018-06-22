# =======================
#      LSHEP Wrapper     
# =======================

# Wrapper class for using the LSHEP fortran code
class LSHEP(Approximator):
    def __init__(self):
        self.linear_shepard = fmodpy.fimport(
            os.path.join(CWD,"linear_shepard","linear_shepard.f95"),
            module_link_args=["-lgfortran","-lblas","-llapack"], 
            output_directory=CWD)
        self.lshep = self.linear_shepard.lshep
        self.lshepval = self.linear_shepard.lshepval
        self.ierrors = {}
        self.x = self.f = self.a = self.rw = None

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values):
        # Local operations
        self.x = np.asfortranarray(control_points.T)
        self.f = np.asfortranarray(values)
        self.a = np.ones(shape=self.x.shape, order="F")
        self.rw = np.ones(shape=(self.x.shape[1],), order="F")
        # In-place update of self.a and self.rw
        self.lshep(self.x.shape[0], self.x.shape[1],
                   self.x, self.f, self.a, self.rw)

    # Use fortran code to evaluate the computed boxes
    def predict(self, x):
        # Calculate all of the response values
        response = []
        for x_pt in x:
            x_pt = np.array(np.reshape(x_pt,(self.x.shape[0],)), order="F")
            resp = self.lshepval(x_pt, self.x.shape[0],
                                 self.x.shape[1], self.x,
                                 self.f, self.a, self.rw)
            response.append(resp)
            # self.ierrors[ier] = self.ierrors.get(ier,0) + 1
        # Return the response values
        return response

    # Print and return a summary of the errors experienced
    def errors(self):
        print("%i normal returns."%self.ierrors.get(0,0))
        print(("%i points outside the radius of influence "+
              "of all other points.")%self.ierrors.get(1,0))
        print(("%i points at which hyperplane normals are significantly"+
               " different.\n  This is only hazardous for piecewise "+
               "linear underlying functions.")%self.ierrors.get(2,0))
        return self.ierrors.get(0,0), self.ierrors.get(1,0), self.ierrors.get(2,0)

