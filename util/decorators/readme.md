<h1 align="center"><code>util.decorators</code></h1>

#### [`def same_as`](decorators.py#L24)

Decorator that copies the documentation and arguemnts of another function (specified as input).

#### [`def cache`](decorators.py#L71)

Cache (input, output) pairs and use cached values for faster re-execution. Inputs and outputs must be `pickle`-able (or `dill`-able).

#### [`def stability_lock`](decorators.py#L131)

Use a cache of (input,output) pairs to check if a function maintains the same behaviors after modifications are made.

#### [`def timeout`](decorators.py#235)

Use a system call to cancel a python function (must respond to global interpreter lock) after a certain amount of time has elapsed.

#### [`def type_check`](decorators.py#L296)

Perform (unpythonic) type checking of function inputs before executing a function.

#### [`def capture`](decorators.py#L423)

Capture all outputs to standard output and standard error for decorated function and store them in attributes `.sdtout` and `.stderr`.

#### [`def background`](decorators.py#526)

Make calls to the decorated function run asynchronously (in the background) and immediately return a `Result` object. Get the actual returned value with a blocking access to `Result.result`.
