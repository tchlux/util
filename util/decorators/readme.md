# util.decorators

Decorator function "same_as" makes the signature and documentation of a function copy another.

Decorator function "cache" generates a unique file name for (input,output) pairs from a function and stores the pair in a serialized file for faster re-execution.

Decorator function "stability_lock" uses a cache of (input,output) pairs to check if a function maintains the same behaviors after modifications are made.

Decorator function "timeout" uses a system call to cancel a python function (must respond to global interpreter lock) after a certain amount of time has elapsed.

Decorator function "type_check" performs (unpythonic) type checking of function inputs before executing a function.
