# Get the name of a class as listed in python source file.
def class_name(obj):
    if type(obj) == type: return (repr(obj)[1:-2].split(".")[-1])
    else:                 return (repr(obj)[1:-2].split(".")[-1]).split(" ")[0]
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
from util.approximate.neural_network import NeuralNetwork, BFGS1000, SGD10K
from util.approximate.support_vector_machine import SVR

 # Create aliases to the neural network.
MLP = NeuralNetwork
NN = SGD10K

# Rename nearest neighbor.
KNN = NearestNeighbor

# Make a "NearestNeighbor" called "CNN" with automatic conditioning.
CNN = condition(NearestNeighbor)
