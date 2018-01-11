# See why and where the prediction error is so large
#  - check 4D data first, it's bad in some spots
#  - verify the convex hull?
#  - look at plots of th 1D predictions, sanity check

import numpy as np
import time
from util.plotly import *
from util.algorithms import BBS, FitBoxMesh, MaxBoxMesh

FUN = lambda x: np.cos(x[0]) * (np.sin(x[1]) if len(x) > 1 else 1.0) + 1
FBM_NUM_BOXES = float("inf")
BBS_NUM_BOXES = float("inf")
BBS_BATCH_SIZE = 4

EXTRA_ARGS = {
    "FitBoxMesh": [FBM_NUM_BOXES],
    "BBS": [BBS_NUM_BOXES, BBS_BATCH_SIZE]
}
PLOT_POINTS = 1000

# Preparing the data
normalize_points = False
random_points = False
random_perturbation = False
tight_plot = False
plot_results = True

perturbation_multiplier = 1.0
num_points = 2
dim = 2
multiplier = 10

algorithms = [BBS()]
# algorithms = [FitBoxMesh()]
# algorithms = [MaxBoxMesh()]

# BOXES:   [[ 0.   0. ]
#           [ 0.   0.1]
#           [ 0.1  0. ]
#           [ 0.1  0.1]]
# WIDTHS:  [ 0.1  0.1  0.1  0.1]
# WEIGHTS: [ 1.          1.09983342  1.          1.09933467]

np.random.seed(0)

# Generate testing points in a grid pattern
plot_range = [[-0.2,1.2]]*dim # [[0,1]]*dim # 
if random_points:
    # Generate testing points randomly spaced
    points = np.random.random(size=(num_points**dim,dim)) * (num_points-1)*num_points
else:
    points = np.meshgrid(*[np.linspace(0,num_points-1,num_points) for i in range(dim)])
    points = np.array([p.flatten() for p in points]).T * num_points

if random_perturbation:
    points += np.random.random(size=points.shape)*perturbation_multiplier

points *= 0.05
# Calculate the associated response values
values = np.array([[FUN(pt) for pt in points]]).T
points = np.concatenate((points, values), axis=1)
# Sort the points so they resemble grid order
points = points[points[:,1].argsort()]
points = points[points[:,0].argsort()]

print()
print("Points:   %s"%(str(points[:,:-1].shape)))#,points[:10,:-1])
print("Response: %s"%(str(points[:, -1].shape)))#,points[:10,-1])

if normalize_points:
    # Normalize the points themselves
    max_val = float(np.max(points[:,:-1]))
    min_val = float(np.min(points[:,:-1]))
    points[:,:-1] = (points[:,:-1] - min_val) / (max_val - min_val)
    # Normalize the response values
    max_val = float(np.max(points[:,-1]))
    min_val = float(np.min(points[:,-1]))
    points[:,-1] = (points[:,-1] - min_val) / (max_val - min_val) + 2.0
else:
    min_val = np.min(points, axis=0)
    max_val = np.max(points, axis=0)
    plot_range = [[ plot_range[i][0] * (max_val[i]-min_val[i]) + min_val[i],
                    plot_range[i][1] * (max_val[i]-min_val[i]) + min_val[i] ]
                  for i in range(dim)]

if tight_plot:
    plot_range = [[min(points[:,i]),max(points[:,i])]
                  for i in range(dim)]


print("Beginning to fit the data...")
for s in algorithms:
    name = (str(s.__class__).split("'")[1]).split(".")[-1]
    start = time.time()
    s.fit(points[:,:-1], points[:,-1], *EXTRA_ARGS.get(name,[]))
    total = time.time()-start
    print("%s: %s seconds"%(name,total))


print(s.boxes.T)
print(s.widths.T)
print(s.weights.T)

if plot_results:
    # Plotting the results
    print()
    print("Plotting..")
    print()
    p = Plot()
    p.add("Raw data",*(points.T))
    use_gradient = len(algorithms) <= 1
    surf_x = np.meshgrid(*[
        np.linspace(lower,upper,PLOT_POINTS**(1/dim))
        for (lower,upper) in plot_range])
    surf_x = np.array([x.flatten() for x in surf_x]).T
    marker_args = {"marker_size":3, "marker_line_width":1, "opacity":0.8}
    if len(algorithms) <= 1: marker_args = {}
    for s in algorithms:
        name = (str(s.__class__).split("'")[1]).split(".")[1]
        print("Adding '%s'..."%(name))
        start = time.time()
        surf_y = s(surf_x.copy())
        total = time.time() - start
        print("  %.2f seconds."%total)
        p.add(name, *(surf_x.T), surf_y, use_gradient=use_gradient,
              plot_type='surface',
              #mode="markers", marker_size=3, marker_line_width=1,
              **marker_args,
        )
    p.plot()

