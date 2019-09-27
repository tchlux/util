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

## [`util.math.points`](points.py)

#### [`def chebyshev`](points.py#L9)

#### [`def mesh`](points.py#L22)

#### [`def polynomials`](points.py#L62)

#### [`def polynomial_indices`](points.py#L75)

#### [`def fekete_indices`](points.py#L91)

#### [`def fekete_points`](points.py#L106)

## [`util.math.weakly_admissable_meshes`](weakly_admissable_meshes.py)

#### [`def polar_wam`](weakly_admissable_meshes.py#L3)

#### [`def box_wam`](weakly_admissable_meshes.py#L19)

