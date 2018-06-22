import fmodpy
import numpy as np
from util.algorithms import Approximator
vtdelaunay = fmodpy.fimport("VTdelaunay.f90")

# ============================
#      VTDelaunay Wrapper     
# ============================

# Wrapper class for using the Delaunay fortran code
class Delaunay(Approximator):
    def __init__(self):
        self.delaunayinterp = vtdelaunay.delaunayinterp
        self.pts = None
        self.values = None
        self.errs = {}

    # Use fortran code to compute the boxes for the given data
    def fit(self, control_points, values=None):
        self.pts = np.asarray(control_points.T, dtype=np.float64, order="F")
        # Store values if they were given
        if type(values) != type(None):
            self.values = np.asarray(values.reshape((1,values.shape[0])),
                                     dtype=np.float64, order="F")
        else:
            self.values = None

    # Return just the points and the weights associated with them for
    # creating the correct interpolation
    def points_and_weights(self, x):
        single_response = len(x.shape) == 1
        if single_response:
            x = np.array([x])
        if len(x.shape) != 2:
            raise(Exception("ERROR: Bad input shape."))

        # Get the predictions from VTdelaunay
        pts_in = np.asarray(self.pts.copy(), order="F")
        p_in = np.asarray(x.T, dtype=np.float64, order="F")
        simp_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                           dtype=np.int32, order="F")
        weights_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                              dtype=np.float64, order="F")
        error_out = np.ones(shape=(p_in.shape[1],), 
                            dtype=np.int32, order="F")
        self.delaunayinterp(self.pts.shape[0], self.pts.shape[1],
                            pts_in, p_in.shape[1], p_in, simp_out,
                            weights_out, error_out, extrap_opt=1.0)
        error_out = np.where(error_out == 1, 0, error_out)
        # Handle any errors that may have occurred.
        if (sum(error_out) != 0):
            unique_errors = sorted(np.unique(error_out))
            print(" [Delaunay errors:",end="")
            for e in unique_errors:
                if (e in {0,1}): continue
                print(" %3i"%e,"at",""+",".join(tuple(
                    str(i) for i in range(len(error_out))
                    if (error_out[i] == e)))+"}", end=";")
            print("] ")
            # Reset the errors to simplex of 1s (to be 0) and weights of 0s.
            bad_indices = (error_out > 1)
            simp_out[:,bad_indices] = 1
            weights_out[:,bad_indices] = 0
        # Adjust the output simplices and weights to be expected shape
        points  = simp_out.T - 1
        weights = weights_out.T
        # Return the appropriate shaped pair of points and weights
        return (points[0], weights[0]) if single_response else (points, weights)
        

    # Use fortran code to evaluate the computed boxes
    def predict(self, x, debug=False, verbose=False):
        if type(self.values) == type(None):
            raise(Exception("ERROR: Must provide values in order to get predictions."))
        # Calculate all of the response values
        response = []
        for pt in x:
            pts, wts = self.points_and_weights(pt)
            response.append( sum(self.values[:,pts][0] * wts) )
        return np.array(response)


if __name__ == "__main__":
    from util.plotly import Plot
    mult = 5
    fun = lambda x: np.cos(x[0]*mult) + np.sin(x[1]*mult)
    p = Plot()
    low = 0
    upp = 1
    dim = 2
    plot_points = 2000
    N = 20
    random = True
    if random:
        x = np.random.random(size=(N,dim))
    else:
        N = int(round(N ** (1/dim)))
        x = np.array([r.flatten() for r in np.meshgrid(np.linspace(low,upp,N), np.linspace(low,upp,N))]).T
    y = np.array([fun(v) for v in x])
    p.add("Training Points", *x.T, y)
    
    # # help(vtdelaunay)
    # # exit()
    # print(x)
    # print(y)
    # exit()

    surf = Delaunay()
    surf.fit(x,y)
    p.add_func(str(surf), surf, *([(low,upp)]*dim), plot_points=plot_points)
    p.plot(file_name="test_plot.html")


