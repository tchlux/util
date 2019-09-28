<h1 align="center"><code>util.system</code></h1>

Function [`run`]() is a (python+OS)-safe interface to command-line execution that cleanly handles errors. Class [`AtomicOpen`]() provides an atomic file operation class that uses system locking mechanisms to enforce atomic operations on files.

#### [`class AtomicOpen`](system.py)

Class for ensuring that all file operations are atomic. Treat initialization like a standard call to `open` that happens to be atomic, must be used in a `with` block.

#### [`class Timer`](system.py)

Class for timing operations. Initialize to start, call to check, use the `start` and `stop` methods / attributes to set / observe values.

#### [`def run`](system.py)

Executing a blocking shell command with a subprocess. On completion provided the return code, stdout as a list of strings, and stderr as a list of strings. This works cross-platform, on Python2.7, and on Python3.x.

#### [`def hash`](system.py#L1)

Generate hash values using a near-zero probability of conflict hashing mechanism that can handle non-hashable types. Does this by hashing the raw bytes (from `pickle` serialization) with `sha256`.

#### [`def save`](system.py#L24)

Save any python object to a file (with `pickle` or `dill`).

#### [`def load`](system.py#L36)

Load python object from file created with `save`.

#### [`def disassemble`](system.py)

Break a file into chunks of smaller sizes (byte-wise).

#### [`def assemble`](system.py)

Undo the `disassemble` operation, reassembling a file (byte-wise) in the correct order.

#### [`def pause`](system.py)

Pause execution until the user presses `<Enter>`.

#### [`def sorted_unique`](system.py#L9)

Generate a sorted unique set of values from an iterable, if the values are not hashable, use `util.system.hash` to sort them.
