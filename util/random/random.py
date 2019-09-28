# This module provides some utilties involving random generation. The
# provided functions in this module are:
# 
#   random_ranage -- Generates a random sample from a range-like object.
#   cdf           -- Generates a random cumulative distribution function
#                    (designed to be used as a random variable).
#   ball          -- Generates random points inside N-dimension ball.
#   latin         -- Generates random points with latin cube design (well-spaced).
#   pairs         -- Generates random pairs of indices within a range.

# Custom excepion for improper user-specified ranges.
class InvalidRange(Exception):
    def __init__(self, start, stop, step, count):
        return super().__init__("\n\nError: random_range(start=%s,stop=%s,step=%s)[:count=%s] contains no elements."%(start,stop,step,count))

# Return a randomized "range" using the appropriate technique based on
# the size of the range being generated. If the memory footprint is
# small (<= 32K int's) then a random sample is created and returned.
# If the memory footprint would be large, a Linear Congruential
# Generator is used to efficiently generate the sequence.
# 
# Parameters are similar to the builtin `range` with:
#   start -- int, default of 0.
#   stop  -- int > start for no / positive step, < start otherwise.
#   step  -- int (!= 0), default of 1.
#   count -- int > 0, number of samples, default (stop-start)//step.
# 
# Usage (and implied parameter ordering):
#   random_range(a)             --> range(0, a, 1)[:a]
#   random_range(a, b) [a < b]  --> range(a, b, 1)[:b-a]
#   random_range(a, b, c)       --> range(a, b, c)[:(b-a)//c]
#   random_range(a, b, c, d)    --> range(a, b, c)[:min(d, (b-a)//c)]
#   random_range(a, d) [d <= a] --> range(0, a, 1)[:d]
# 
# If the size of the range is large, a Linear Congruential Generator is used.
#   Memory  -- O(1) storage for a few integers, regardless of parameters.
#   Compute -- O(n) at most 2 times the number of steps in the range, n.
# 
def random_range(start, stop=None, step=None, count=float('inf')):
    from random import sample, randint
    from math import ceil, log2
    # Add special usage where the second argument is meant to be a count.
    if (stop != None) and (stop <= start) and ((step == None) or (step >= 0)):
        start, stop, count = 0, start, stop
    # Set a default values the same way "range" does.
    if (stop == None): start, stop = 0, start
    if (step == None): step = 1
    # Compute the number of numbers in this range, update count accordingly.
    num_steps = int((stop - start) // step)
    count = min(count, num_steps)
    # Check for a usage error.
    if (num_steps == 0) or (count <= 0): raise(InvalidRange(start, stop, step, count))
    # Use robust random method if it has a small enough memory footprint.
    # Do not use this method if there is a number too large for C implementation.
    if (num_steps <= 2**15) and (max(map(abs,(start,stop,step))) < 2**63):
        for value in sample(range(start,stop,step), count): yield value
        return
    # Use the LCG for the cases where the above is too memory intensive.
    # Enforce that all numbers are python integers (no worry of overflow).
    start, stop, step = map(int, (start, stop, step))
    # Use a mapping to convert a standard range into the desired range.
    mapping = lambda i: (i*step) + start
    # Seed range with a random integer to start.
    value = int(randint(0,num_steps))
    # 
    # Construct an offset, multiplier, and modulus for a linear
    # congruential generator. These generators are cyclic and
    # non-repeating when they maintain the properties:
    # 
    #   1) "modulus" and "offset" are relatively prime.
    #   2) ["multiplier" - 1] is divisible by all prime factors of "modulus".
    #   3) ["multiplier" - 1] is divisible by 4 if "modulus" is divisible by 4.
    # 
    offset = int(randint(0,num_steps)) * 2 + 1                 # Pick a random odd-valued offset.
    multiplier = 4*(num_steps + int(randint(0,num_steps))) + 1 # Pick a multiplier 1 greater than a multiple of 4.
    modulus = int(2**ceil(log2(num_steps)))                    # Pick a modulus just big enough to generate all numbers (power of 2).
    # Track how many random numbers have been returned.
    found = 0
    while found < count:
        # If this is a valid value, yield it in generator fashion.
        if value < num_steps:
            found += 1
            yield mapping(value)
        # Calculate the next value in the sequence.
        value = (value*multiplier + offset) % modulus


# Generate a random cumulative distribution function by picking
# "nodes" random step sizes. Return twice differentiable CDF fit.
def cdf(nodes=9, **kwargs):
    import numpy as np
    from util.stats import cdf_fit
    if ("fit" not in kwargs): kwargs["fit"] = "cubic"
    nodes += 2
    # Generate random steps in the x and y direction (that sum to 1).
    x = np.linspace(0, 1, nodes)
    y = [0.]
    for i in range(1,nodes-1):
        y.append( y[-1] + np.random.random()**(nodes-i-1) * (1. - y[-1]) )
    y.append(1.)
    y = np.array(y)
    print("Random cdf:", x.shape, y.shape)
    # Return the CDF fit over the generate CDF points.
    return cdf_fit((x,y), **kwargs)


# Generate "count" random indices of pairs that are within "length" bounds.
def pairs(length, count=None):
    from util.math import index_to_pair
    # Compute the hard maximum in the number of pairs.
    max_pairs = length*(length - 1) // 2
    # Initialize count if it wasn't provided.
    if type(count) == type(None): count = max_pairs
    count = min(count, max_pairs)
    # Get a random set of pairs (no repeats).
    for c,i in enumerate(random_range(count)):
        if (i >= count): break
        yield index_to_pair(i)
    print(" "*40, end="\r", flush=True)


# Generate "num_points" random points in "dimension" that have uniform
# probability over the unit ball scaled by "radius" (length of points
# are in range [0, "radius"]).
def ball(num_points, dimension, radius=1):
    from numpy import random, linalg
    # First generate random directions by normalizing the length of a
    # vector of random-normal values (these distribute evenly on ball).
    random_directions = random.normal(size=(dimension,num_points))
    random_directions /= linalg.norm(random_directions, axis=0)
    # Second generate a random radius with probability proportional to
    # the surface area of a ball with a given radius.
    random_radii = random.random(num_points) ** (1/dimension)
    # Return the list of random (direction & length) points.
    return radius * (random_directions * random_radii).T


# Generate random points with a Latin Hyper Cube design.
def latin(num_points, dimension):
    from numpy import random, ones, arange
    # Generate a "num_points^dimension" grid, making sure every grid
    # cell has at least one point in one of its dimensions.
    cell_width = 1 / num_points
    cells = ones((dimension, num_points)) * arange(num_points)
    # Randomly shuffle the selected grid cells for each point.
    for _ in map(random.shuffle, cells): pass
    # Convert the random selected cells into points.
    return cell_width * (random.random((dimension,num_points)) + cells).T
