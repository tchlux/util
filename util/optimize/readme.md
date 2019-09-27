# util.optimize

Function "minimize" uses a meta-heuristic optimization algorithm to solve a global optimization problem given an arbitrary function.

##### [`class Tracker`](https://github.com/tchlux/util/blob/master/util/optimize/__init__.py#L23)

This class is designed to wrap an objective function and automatically track the amount of time spent optimizing, number of evaluations, objective function values, and improvement rates. It's primary purpose is to provide the `Tracker.check` method (to be used as if an objective function) and the `Tracker.done` method which returns `True` or `False` based on initialization stopping criteria. Finally, the `Tracker` provides visualizations and printouts of improvement in performance.

##### [`def minimize`](https://github.com/tchlux/util/blob/master/util/optimize/__init__.py#L124)</h4>

Minimize arbitrary objective function given an initial solution and optionally given bounds, max computation time, minimum number of steps, or a minimum improvement stopping criteria. Supports optimization with any gradient-free optimization method.

##### [`util.optimize.gradient_free`](https://github.com/tchlux/util/blob/master/util/optimize/gradient_free.py)

Implementations of common gradient free methods of optimization [`Adaptive Normal`](https://github.com/tchlux/util/blob/master/util/optimize/adaptive_normal.py#L3) (quasi simulated annealing), [`AMPGO`](https://github.com/tchlux/util/blob/master/util/optimize/ampgo.py#L4) (adaptive memory programming for gradient-free optimization, `L-BFGS` over a surrogate with tunneling), [`DiRect`](https://github.com/tchlux/util/blob/master/util/optimize/direct.py#L77) (divided rectangles), and [`random`](https://github.com/tchlux/util/blob/master/util/optimize/random.py#L3) (pure random search).

##### [`util.optimize.gradient_based`](https://github.com/tchlux/util/blob/master/util/optimize/gradient_based.py)

Implementations of the common gradient based optimization algorithms, [`SGD`](https://github.com/tchlux/util/blob/master/util/optimize/gradient_based.py#L18) (stochastic gradient descent), [`ADAGRAD`](https://github.com/tchlux/util/blob/master/util/optimize/gradient_based.py#L39) (adaptive gradient descent), [`ADAM`](https://github.com/tchlux/util/blob/master/util/optimize/gradient_based.py#L63) (adaptive moment estimation), and [`L_BFGS`](https://github.com/tchlux/util/blob/master/util/optimize/gradient_based.py#L7) (Limited memory Broyden-Fletcher-Goldfarb-Shanno, through scipy).

##### [`util.optimize.zero`](https://github.com/tchlux/util/blob/master/util/optimize/newton.py#L1)

An implementation of the Newton secant method for one-dimensional optimization.
