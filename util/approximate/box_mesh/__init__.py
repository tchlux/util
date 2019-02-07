# Makes available all custom algorithms that I have worked with
import os
import numpy as np
import fmodpy
from util.approximate import WeightedApproximator

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))

# ======================
#      BoxMesh Mesh     
# ======================
class BoxMesh(WeightedApproximator):
    meshes = fmodpy.fimport(
        os.path.join(CWD,"meshes.f90"), output_directory=CWD, autocompile_extra_files=True)
    error_tolerance = 0.0

    # Fit a set of points
    def _fit(self, points):
        from util.math import is_none
        # Sort points by their distance from the center of the data.
        center = (np.max(points, axis=0) - np.min(points, axis=0)) / 2
        dists = np.linalg.norm(points - center, axis=1)
        indices = np.argsort(dists)
        if not is_none(self.y): self.y = [self.y[i] for i in indices]
        # Get points in a specific order.
        self.points = np.asarray(points[indices].T, order="F")
        self.box_sizes = np.ones((self.points.shape[0]*2,self.points.shape[1]),
                                  dtype=np.float64, order="F") * -1
        self.meshes.build_ibm(self.points, self.box_sizes)

    # Generate a prediction for a new point
    def _predict(self, xs):
        to_predict = np.asarray(xs.T, order="F")
        all_weights = np.zeros((to_predict.shape[1], self.points.shape[1]), 
                               dtype=np.float64, order="F")
        self.meshes.eval_box_mesh(self.points, self.box_sizes,
                                  to_predict, all_weights)
        idx = np.arange(self.points.shape[1])
        indices = []
        weights = []
        # Normalize the weights (so they are convex).
        for i in range(len(all_weights)):
            to_keep = idx[all_weights[i,:] > 0]
            indices.append( to_keep )
            weights.append( all_weights[i,to_keep] )
        return indices, weights


if __name__ == "__main__":
    from util.plot import Plot
    from util.approximate.testing import test_plot
    model = BoxMesh()

    p,x,y = test_plot(model)
    p.plot(show=False)
    x = model.points.T
    print(x)


    # ==============================================
    #      Display the boxes that were computed     
    # ==============================================
    p = Plot()
    # Get the extreme points.
    min_x = np.min(x[:,0]) - .1
    max_x = np.max(x[:,0]) + .1
    min_y = np.min(x[:,1]) - .1
    max_y = np.max(x[:,1]) + .1
    # Get the box edges (about the centers).
    boxes = model.box_sizes.T
    boxes[:,0] = np.where(boxes[:,0] != -1, boxes[:,0], min_x)
    boxes[:,1] = np.where(boxes[:,1] != -1, boxes[:,1], min_y)
    boxes[:,2] = np.where(boxes[:,2] != -1, boxes[:,2], max_x)
    boxes[:,3] = np.where(boxes[:,3] != -1, boxes[:,3], max_y)
    # Add boxes to plot.
    for i,(pt,bounds) in enumerate(zip(x, boxes)):
        shifts = np.array([[- boxes[i,0], - boxes[i, 1]],
                           [  boxes[i,2], - boxes[i, 1]],
                           [  boxes[i,2],   boxes[i, 3]],
                           [- boxes[i,0],   boxes[i, 3]],
                           [- boxes[i,0], - boxes[i, 1]] ])
        box = (pt + shifts)
        p.add(f"Box {i+1}", *(box.T), color=p.color(i), mode="lines")
    for i,(pt,bounds) in enumerate(zip(x, boxes)):
        p.add("", [pt[0]], [pt[1]], color=p.color(i), show_in_legend=False)
    p.show(append=True)
    exit()

