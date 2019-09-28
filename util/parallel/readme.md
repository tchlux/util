<h1 align="center"><code>util.parallel</code></h1>


#### [`def shared_array`](parallel.py#L20)

Given a shape, create a `numpy` `float64` array that will be shared (not copied) during multiprocessing (with atomic access), but behaves like a normal array.


