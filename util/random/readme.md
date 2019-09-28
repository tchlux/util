<h1 align="center"><code>util.random</code></h1>

Functions for generating useful random sequences ([ranges](random.py#L16), [distributions](random.py#L88), [pairs](random.py#L107), [Latin designs](random.py#L138)).

#### [`def random_range`](random.py#L16)

Behaves exactly like a standard Python `range` except that it has been randomly shuffled. Uses memory-efficient shuffling technique (linear congruential generator) for large number ranges (>32K integers) to provide a `O(1)` memory footprint.

#### [`def cdf`](random.py#L88)

Generates a random (twice differentiable) cumulative distribution function over `[0,1]`.

#### [`def pairs`](random.py#L107)

Generate random pairs of numbers (good for picking random pairs of elements from a large sequence).

#### [`def ball`](random.py#L122)

Generate random points in the unit ball for given dimension.

#### [`def latin`](random.py#L138)

Generate random points in a Latin Hypercube design for given dimension (when space is partitioned into a grid, each row / col / etc. has exactly 1 point in one of its cells).
