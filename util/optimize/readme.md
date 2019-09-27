# util.optimize

Function "minimize" uses a meta-heuristic optimization algorithm to solve a global optimization problem given an arbitrary function.

#### [class Tracker](https://github.com/tchlux/util/blob/master/util/optimize/__init__.py#L23)

<h4><code><a href="https://github.com/tchlux/util/blob/master/util/optimize/__init__.py#L23">class Tracker</a></code></h4>

```python
<a href="https://github.com/tchlux/util/blob/master/util/optimize/__init__.py#L23">class Tracker</a>
```

This class is designed to wrap an objective function and automatically track the amount of time spent optimizing, number of evaluations, objective function values, and improvement rates. It's primary purpose is to provide the `Tracker.check` method (to be used as if an objective function) and the `Tracker.done` method which returns `True` or `False` based on initialization stopping criteria. Finally, the `Tracker` provides visualizations and printouts of improvement in performance.

<h4><code><a href="https://github.com/tchlux/util/blob/master/util/optimize/__init__.py#L124">def minimize</a></code></h4>

Minimize arbitrary objective function given an initial solution and
optionally given bounds, max computation time, minimum number of
steps, or a minimum improvement stopping criteria.
