<h1 align="center"><code>util.system</code></h1>

Function [`run`](system.py#L146) is a (python+OS)-safe interface to command-line execution that cleanly handles errors. Class [`AtomicOpen`](system.py#L181) provides an atomic file operation class that uses system locking mechanisms to enforce atomic operations on files. Also provides a robust [`hash`](system.py#L1) function, easy-to-use [`save`](system.py#L24) / [`load`](system.py#L36) functions, and a [`Timer`](system.py#L235) object.

#### [`class AtomicOpen`](system.py#L181)

Class for ensuring that all file operations are atomic. Treat initialization like a standard call to `open` that happens to be atomic, must be used in a `with` block.

#### [`class Timer`](system.py#L235)

Class for timing operations. Initialize to start, call to check, use the `start` and `stop` methods / attributes to set / observe values.

#### [`def run`](system.py#L146)

Executing a blocking shell command with a subprocess. On completion provided the return code, stdout as a list of strings, and stderr as a list of strings. This works cross-platform, on Python2.7, and on Python3.x.

#### [`def hash`](system.py#L1)

Generate hash values using a near-zero probability of conflict hashing mechanism that can handle non-hashable types. Does this by hashing the raw bytes (from `pickle` serialization) with `sha256`.

#### [`def save`](system.py#L24)

Save any python object to a file (with `pickle` or `dill`).

#### [`def load`](system.py#L36)

Load python object from file created with `save`.

#### [`def disassemble`](system.py#L59)

Break a file into chunks of smaller sizes (byte-wise).

#### [`def assemble`](system.py#L94)

Undo the `disassemble` operation, reassembling a file (byte-wise) in the correct order.

#### [`def pause`](system.py#L50)

Pause execution until the user presses `<Enter>`.

#### [`def sorted_unique`](system.py#L9)

Generate a sorted unique set of values from an iterable, if the values are not hashable, use `util.system.hash` to sort them.
