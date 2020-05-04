# Generate hash values using a near-zero probability of conflict
# hashing mechanism that can handle non-hashable types. Does this
# by hashing the raw bytes (from pickle) with sha256.
def hash(something):
    import hashlib, pickle
    return hashlib.sha256(pickle.dumps(something)).hexdigest()


# Generate a sorted unique set of values from an iterable, if the
# values are not hashable, use a hashing function to sort them.
def sorted_unique(iterable):
    try:
        # Try returning sorted unique values using hash to find unique.
        return sorted(set(iterable))
    except TypeError: 
        # Get the set of unique hashes and associated values.
        ids = {hash(v):v for v in iterable}
        # Try sorting the values directly, if that's not possible
        # then sort the values by their hash.
        try:    return sorted(ids.values())
        except: return [ids[h] for h in sorted(ids)]


# Save "data" in a file titled "file_name" using pickle.
def save(data, file_name="_save.pkl"):
    try:
        import pickle
        with open(file_name, "wb") as f:
            pickle.dump(data, f)
    except:
        import dill as pickle
        with open(file_name, "wb") as f:
            pickle.dump(data, f)


# Load data from a pickle file titled "file_name".
def load(file_name="_save.pkl"):
    try:
        import pickle
        with open(file_name, "rb") as f:
            data = pickle.load(f)
        return data
    except:
        import dill as pickle
        with open(file_name, "rb") as f:
            data = pickle.load(f)
        return data


# Convenience function for pausing a program for input (explicitly
# named). Uses 'getpass' in order to suppress new line character.
def pause(string="press enter to continue.."):
    import os, getpass
    print(string, end="\r")
    with open(os.devnull,"w") as stream:
        response = getpass.getpass(stream=stream)
    print(" "*len(string), end="\r")

# Disassemble a file larger than "max_size" into chunks of size
# "max_size". Place in a folder with the same name as the file, but
# having extension names generated by "chunk_ext".
# 
# ARGUMENTS:
#   file_path -- A string representing an absolute or relative path to a file.
# 
# OPTIONAL:
#   max_size_bytes -- An integer representing the maximum number of
#                     bytes that a file chunk can contain. [50MB]
#   chunk_ext      -- A string. [".chunks"]
#   verbose        -- Prints out status information when True. [True]
def disassemble(file_path, max_size_bytes=50*(2**20), skip_size=1024,
                chunk_ext=".chunks", verbose=True):
    class FileTooSmall(Exception): pass
    import os, math
    size = os.path.getsize(file_path)
    chunks = math.ceil(size / max_size_bytes)
    if (chunks <= 1): raise(FileTooSmall("The provided file is smaller than 'max_size'."))
    out_path = os.path.join(os.curdir, file_path + chunk_ext)
    if verbose:
        print()
        print("File path:  ", file_path)
        print("Out path:   ", out_path)
        print("File size:  ", size)
        print("File chunks:", chunks)
    os.makedirs(out_path, exist_ok=True)
    with open(file_path, "rb") as f:
        for c in range(1,chunks+1):
            # Open the output file (with appropriate name) and read
            # out a chunk to that output file.
            with open(os.path.join(out_path,str(c)), "wb") as out:
                if verbose:
                    print(" writing chunk %i of %i in '%s'..."%(c,chunks,out.name))
                # Attempt to read the chunk (safely, replacing bad bits with 0's).
                try: chunk = f.read(max_size_bytes)
                except OSError:
                    print( "  WARNING: Encountered unreadable data, attempting partial read.")
                    # If there was a failed I/O operation, try reading
                    # a smaller number of bytesfrom the file.
                    bad_bytes = []
                    read_bytes = 0
                    chunk_size = max_size_bytes // 2
                    chunk = bytes()
                    while (read_bytes < max_size_bytes):
                        chunk_size = min(chunk_size, max_size_bytes - read_bytes)
                        print(f"    read {read_bytes} of {max_size_bytes} (trying {chunk_size})"+" "*20, end="\r")
                        try:
                            chunk += f.read(chunk_size)
                            read_bytes += chunk_size
                            chunk_size *= 2
                        except OSError:
                            if chunk_size > skip_size:
                                chunk_size = max(skip_size,chunk_size // 10)
                            elif chunk_size <= skip_size:
                                # Get the index of the bad byte (block).
                                bad_bytes.append( f.seek(0,1) )
                                # Add a null byte to the stream in place of the truth.
                                chunk += bytes(chunk_size)
                                read_bytes += chunk_size
                                # Seek forward one byte relative to current position, skip.
                                f.seek(chunk_size,1)
                    # Let the user know how the bad bytes were handled.
                    print(f"           Found {len(bad_bytes)} bad bytes and replaced them with 0's.")
                    if len(bad_bytes) < 100:
                        print(f"           Bad bytes are located at indices:")
                        print(f"             {bad_bytes}")                    
                    # Now the full amount of data has been read into "chunk".
                # Write the output chunk to file.
                out.write(chunk)

                
# This function provides the counterpart to `disassemble`.
def assemble(names=(), chunk_ext=".chunks", verbose=True):
    if verbose: print(" Beginning `assemble`")
    # Get all of the paths to the pieces of this file.
    import os
    output_dir = os.getcwd()
    if (len(names) == 0):
        names = [p for p in os.listdir() if p[-len(chunk_ext):] == chunk_ext]
        if verbose:print(f"   automatically assembling all of:\n   ",'\n    '.join(names))
    # Re-assemble all of the chunked files.
    for name in names:
        if verbose: print(f"  assembling {name}")
        output_path = os.path.join(output_dir, name[:-len(chunk_ext)])
        # Skip the file if it already exists (do not overwrite).
        if (os.path.exists(output_path)):
            print(f"  file '{output_path}' already exists, skipping")
            continue
        # Get all of the chunks of this file in order.
        chunk_paths = sorted((n for n in os.listdir(
            os.path.join(output_dir,name))), key=lambda n: int(n))
        if verbose: print(f"    {chunk_paths}")
        with open(output_path, "wb") as out:
            # Concatenate all of the files into the single output file.
            for f_path in sorted(chunk_paths):
                with open(os.path.join(name,f_path), "rb") as f: out.write(f.read())

# Execute a blocking command with a subprocess, on completion provide
# the return code, stdout as string, and stderr as string. This should
# work across both Python2.7 and Python3.x as well as cross-platform.
#  INPUT:
#   command -- A list of strings or string (space separated) describing
#              a standard command as would be given to subprocess.Popen
# 
#  OUTPUT:
#   return_code -- Straight from the subprocess returncode
#   stdout      -- A list of strings that are the lines of stdout 
#                  produced by the execution of <command>
#   stderr      -- A list of strings that are the lines of stderr
#                  produced by the execution of <command>
def run(command, **popen_kwargs):
    import sys, subprocess
    # For Python 2.x the encoding is a string by default
    # For Python 3.6 and later the encoding can be given as an arguemnt
    if sys.version_info >= (3,6):
        popen_kwargs.update( dict(encoding="UTF-8") )
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, **popen_kwargs)
    stdout, stderr = proc.communicate()
    # Before python 3.6, the encoding had to be handled after collection
    if ((sys.version_info >= (3,)) and (sys.version_info[0] == 3)):
        if (type(stdout) != str): stdout = str(stdout, encoding="UTF-8")
        if (type(stderr) != str): stderr = str(stderr, encoding="UTF-8")
    # Remove windows specific characters and split by the new line
    if stdout: stdout = stdout.replace("\r","").split("\n")
    else:      stdout = ""
    if stderr: stderr = stderr.replace("\r","").split("\n")
    else:      stderr = ""
    # Return the exit code, standard out, and standard error
    return proc.returncode, stdout, stderr


# ==========================
#      AtomicOpen Class     
# ==========================

try:
    # Posix based file locking (Linux, Ubuntu, MacOS, etc.)
    import fcntl, os
    def lock_file(f):
        try: fcntl.lockf(f, fcntl.LOCK_EX)
        except OSError: pass
    def unlock_file(f):
        fcntl.lockf(f, fcntl.LOCK_UN)
except ModuleNotFoundError:
    # Windows file locking
    import msvcrt, os
    def file_size(f):
        return os.path.getsize( os.path.realpath(f.name) )
    def lock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_RLCK, file_size(f))
    def unlock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, file_size(f))


# Class for ensuring that all file operations are atomic. Treat
# initialization like a standard call to `open` that happens to be
# atomic, must be used in a `with` block.
class AtomicOpen:
    # First acquire a lock for the file, then open the file with
    # arguments provided by user. Attempts to inherit most of the file
    # properties, but use "__enter__" to get file object directly.
    def __init__(self, path, *args, **kwargs):
        # Open the file and acquire a lock on the file before operating.
        self.file = open(path, *args, **kwargs)
        # Lock the opened file.
        lock_file(self.file)

    # Return the opened file object (knowing a lock has been obtained).
    def __enter__(self, *args, **kwargs): return self.file

    # Unlock and then close the file object.
    def __exit__(self, exc_type=None, exc_value=None, traceback=None):
        # Flush to make sure all buffered contents are written to file.
        self.file.flush()
        os.fsync(self.file.fileno())
        # Release the lock on the file, then close.
        unlock_file(self.file)
        self.file.close()
        # Handle exceptions that may have come up during execution, by
        # default any exceptions are raised to the user.
        if (exc_type != None): return False
        else:                  return True



# Class for timing operations. Initialize to start, call to check, use
# the `start` and `stop` methods / attributes to set / observe values.
class Timer:
    #-----------------------------------------------------------------
    #                            Private
    import time as _time
    _a = 0
    _b = None
    # Offer a function for explicitly beginning a timer.
    def _begin(self):
        self._a = self._time.time()
        self._b = None
        return self._a
    # End function for this timer.
    def _end(self):
        self._b = self._time.time()
        return self._b
    # Return the elapsed time since start.
    def _check(self):
        if   (self._a is None): return 0
        elif (self._b is None): return self._time.time() - self._a
        else:                   return self._b - self._a
    # Overwrite the "__call__" method of the provided object.
    def _callset(obj, func):
        class Float(type(obj)):
            def __call__(self): return func()
            def __repr__(self): return str(float(self))
        return Float(obj)
    #-----------------------------------------------------------------
    #                             Public
    # 
    # If not stopped, return elapsed time, otherwise return total time.
    def __call__(self, precision=2): return float(f"{self._check():.{precision-1}e}")

    # Return "start time" as attribute, "begin timer" as function.
    @property
    def start(self): return Timer._callset(self._a, self._begin)
    # Return "stop time" as attribute, "end timer, return total" as function.
    @property
    def stop(self):
        if (self._b == None): return Timer._callset(self._end(), self.total)
        else:                 return Timer._callset(self._b,     self.total)
    # Return the "total time from start" if running, return "total time" if finished.
    @property
    def total(self):
        if (self._b == None): return Timer._callset(self._check(), self._check)
        else:                 return Timer._callset(self._b - self._a, lambda: self._b - self._a)
