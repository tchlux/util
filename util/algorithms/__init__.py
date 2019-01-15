# Get the name of a class as listed in python source file.
def class_name(obj): return (repr(obj)[1:-2].split(".")[-1])

# Import the base classes for approximation and the wrappers for approximation.
from util.algorithms.base import Approximator, WeightedApproximator
from util.algorithms.wrappers import unique, condition

# Import all defined algorithms
from util.algorithms.voronoi import Voronoi
from util.algorithms.nearest_neighbor import NearestNeighbor
from util.algorithms.delaunay import Delaunay, qHullDelaunay
from util.algorithms.neural_network import MLPRegressor
from util.algorithms.linear_shepard import LSHEP
from util.algorithms.mars import MARS
from util.algorithms.shepard import Shepard
# Rename neural network.
NeuralNetwork = MLPRegressor
NN = MLPRegressor
KNN = NearestNeighbor
# Import additional Delaunay testing modules.
from util.algorithms.delaunay import DelaunayP1, DelaunayP2, DelaunayP3

# Make a "NearestNeighbor" called "CNN" with automatic conditioning.
CNN = condition(NearestNeighbor)
