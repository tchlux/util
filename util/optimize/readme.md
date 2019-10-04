<h1 align="center"><code>util.optimize</code></h1>

Numerous optimization utilities for both gradient-free ([`Adaptive Normal`](adaptive_normal.py#L3), [`DiRect`](direct.py#L77), [`random`](pure_random.py#L3)) and gradient-based ([`SGD`](gradient_based.py#L18), [`ADAM`](gradient_based.py#L63), [`L_BFGS`](gradient_based.py#L7)) minimization of arbitrary functions.

#### [`class Tracker`](__init__.py#L23)

This class is designed to wrap an objective function and automatically track the amount of time spent optimizing, number of evaluations, objective function values, and improvement rates. It's primary purpose is to provide the `Tracker.check` method (to be used as if an objective function) and the `Tracker.done` method which returns `True` or `False` based on initialization stopping criteria. Finally, the `Tracker` provides visualizations and printouts of improvement in performance.

#### [`def minimize`](__init__.py#L121)</h4>

Minimize arbitrary objective function given an initial solution and optionally given bounds, max computation time, minimum number of steps, or a minimum improvement stopping criteria. Supports optimization with any gradient-free optimization method.

#### [`def zero`](newton.py#L1)

An implementation of the Newton secant method for one-dimensional optimization.

#### [`util.optimize.gradient_free`](gradient_free.py)

Implementations of common gradient free methods of optimization [`Adaptive Normal`](adaptive_normal.py#L3) (quasi simulated annealing), [`AMPGO`](ampgo.py#L4) (adaptive memory programming for gradient-free optimization, `L-BFGS` over a surrogate with tunneling), [`DiRect`](direct.py#L77) (divided rectangles), and [`random`](pure_random.py#L3) (pure random search).

#### [`util.optimize.gradient_based`](gradient_based.py)

Implementations of the common gradient based optimization algorithms, [`SGD`](gradient_based.py#L18) (stochastic gradient descent), [`ADAGRAD`](gradient_based.py#L39) (adaptive gradient descent), [`ADAM`](gradient_based.py#L63) (adaptive moment estimation), and [`L_BFGS`](gradient_based.py#L7) (Limited memory Broyden-Fletcher-Goldfarb-Shanno, through scipy).
