<h1 align="center"><code>util.math</code></h1>

#### [`def is_numeric`](__init__.py#L16)

`True` if the given object supports multiplication by `int`, `float`, and addition / subtraction.

#### [`def abs_diff`](__init__.py#L23)

Function for performing absolute difference between numbers (or vectors). Falls back to equality for things that don't support difference, absolute value, or sums.

#### [`def is_prime`](__init__.py#L33)

Return `True` if `n` is prime, `False` otherwise.

#### [`def primes`](__init__.py#L39)

Return a list of tuples of all unique prime factors and multiplicities of `n`.

#### [`def primes_up_to`](__init__.py#L53)

Return all prime numbers up to `n`.

#### [`def choose`](__init__.py#L63)

Compute the value `n choose k` using python integers, no overflow.

#### [`def flatten`](__init__.py#L10)

Given a python list of lists, flatten into a single python list.

#### [`def transpose`](__init__.py#L13)

Given a python list of lists, transpose the lists.

## [`util.math.fraction`](fraction.py)

#### [`class Fraction`](fraction.py#L24)

This is a rational number representation that supports all the standard operators. It uses python integers (no overflow) under the hood and has a `float` fallback for interactions with numbers that are not `Fraction` or `int` typed.

## [`util.math.points`](points.py)

#### [`def chebyshev`](points.py#L9)

The classical optimally-spaced points in one dimension over an interval for polynomial interpolation.

#### [`def mesh`](points.py#L22)

Create a full tensor product of all provided input lists. If a `sample` number is provided, then a random sample is drawn from the full tensor product.

#### [`def polynomials`](points.py#L62)

Construct the list of all polynomial functions in a given dimension up to a specified degree. I.e. in three dimensions, up to degree two would return `[x, y, z, xx, yy, zz, xy, xz, yz]`.

#### [`def polynomial_indices`](points.py#L75)

Similar to `polynomials`, but returns the list of degrees and the list of indices from vectors used by polynomials up to the selected degree in the given dimension.

#### [`def fekete_indices`](points.py#L91)

Return the sorted list of indices that should be used to pick Fekete points from a set. This is computed by QR column pivoting.

#### [`def fekete_points`](points.py#L106)

Given a number of points and dimension, construct and return Fekete points over the unit cube. Uses a polynomial basis to construct the Vandermonde matrix.

## [`util.math.weakly_admissable_meshes`](weakly_admissable_meshes.py)

#### [`def polar_wam`](weakly_admissable_meshes.py#L3)

#### [`def box_wam`](weakly_admissable_meshes.py#L19)

