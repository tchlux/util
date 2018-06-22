# Makes available all custom algorithms that I have worked with
import os
import numpy as np
import fmodpy
from util.algorithms import Approximator, IS_NONE

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))

# ======================
#      BoxMesh Mesh     
# ======================
class BoxMesh(Approximator):
    meshes = fmodpy.fimport(
        os.path.join(CWD,"meshes","meshes.f90"),
        module_link_args=["-lgfortran","-lblas","-llapack"], 
        output_directory=CWD)
    error_tolerance = 0.0

    # Fit a set of points
    def fit(self, points, values):
        self.points = np.asarray(points.T.copy(), order="F")
        self.values = np.asarray(values.copy(), order="F")
        self.box_sizes = np.ones((self.points.shape[0]*2,self.points.shape[1]),
                                  dtype=np.float64, order="F") * -1
        # Construct the box mesh
        for i in range(1,self.points.shape[1]+1):
            self.meshes.build_ibm(self.points[:,:i], self.box_sizes)

    # Generate a prediction for a new point
    def predict(self, xs):
        to_predict = np.asarray(xs.T, order="F")
        predictions = np.ones((xs.shape[0],), dtype=np.float64, order="F")        
        self.meshes.predict_convex_mesh(self.points, self.box_sizes,
                                        to_predict, self.values,
                                        predictions)
        return predictions

