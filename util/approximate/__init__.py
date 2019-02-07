# Get the name of a class as listed in python source file.
def class_name(obj): return (repr(obj)[1:-2].split(".")[-1])

# Import the base classes for approximation and the wrappers for approximation.
from util.approximate.base import Approximator, WeightedApproximator
from util.approximate.wrappers import unique, condition

# Import all defined approximation algorithms.

# Weighted approximators (define predictions as convex weighted sums).
from util.approximate.nearest_neighbor import NearestNeighbor
from util.approximate.delaunay import Delaunay
from util.approximate.voronoi import Voronoi
from util.approximate.box_mesh import BoxMesh
from util.approximate.modified_shepard import ShepMod
from util.approximate.shepard import Shepard
# Standard real approximator.
from util.approximate.linear_shepard import LSHEP
# Regression techniques.
from util.approximate.mars import MARS
from util.approximate.neural_network import MLPRegressor
from util.approximate.support_vector_machine import SVR

# Rename neural network.
MLP = MLPRegressor
NeuralNetwork = MLPRegressor
NN = MLPRegressor
# Rename nearest neighbor.
KNN = NearestNeighbor

# Make a "NearestNeighbor" called "CNN" with automatic conditioning.
CNN = condition(NearestNeighbor)
