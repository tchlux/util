<h1 align="center"><code>util.parallel</code></h1>

Provides a parallelized implementation of the builtin Python [`map`](parallel.py#L84) function for easy drop-in parallelization.

#### [`def map`](parallel.py#L84)

A parallel implementation of the builtin `map` function in Python. Given a function and an iterable, use all cores on the current system to map the function onto all elements of the iterable. Supports both in-order (could suffer load imbalancing) and out-of-order (always load balanced) execution. Standard output and errors are redirected to process-specific log files.

#### [`def shared_array`](parallel.py#L20)

Given a shape, create a `numpy` `float64` array that will be shared (not copied) during multiprocessing (with atomic access), but behaves like a normal array.

#### [`def split`](parallel.py#L178)

Given an iterable, split it into as many chunks as there are cores.

#### [`def killall`](parallel.py#L190)

Kill all currently running parallel processes.
