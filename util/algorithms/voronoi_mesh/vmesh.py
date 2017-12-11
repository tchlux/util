import fmodpy
vmesh = fmodpy.fimport("voronoi_mesh.f90")
import numpy as np

class VoronoiMesh:
    # Fit a set of points
    def fit(self, points, values):
        self.points = np.asarray(points.T.copy(), order="F")
        self.values = np.asarray(values.copy(), order="F")

    # Wrapper for 'predict' that returns a single value for a single
    # prediction, or an array of values for an array of predictions
    def __call__(self, x:np.ndarray, *args, **kwargs):
        single_response = len(x.shape) == 1
        if single_response:
            x = np.array([x])
        if len(x.shape) != 2:
            raise(Exception("ERROR: Bad input shape."))
        response = np.asarray(self.predict(x, *args, **kwargs), dtype=float)
        # Return the response values
        return response[0] if single_response else response

    # Generate a prediction for a new point
    def predict(self, xs):
        to_predict = np.asarray(xs.T, order="F")
        predictions = np.ones( (len(xs),) )
        vmesh.voronoi_mesh(self.points, self.values, to_predict,
                           predictions)
        return predictions
