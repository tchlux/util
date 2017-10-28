import time, os, time

# ==========================
#      AtomicOpen Class     
# ==========================

try:
    # Posix based file locking (Linux, Ubuntu, MacOS, etc.)
    import fcntl
    def lock_file(f):
        try:
            fcntl.lockf(f, fcntl.LOCK_EX)
        except OSError:
            pass
    def unlock_file(f): pass
except ModuleNotFoundError:
    # Windows file locking
    import msvcrt
    def file_size(f):
        return os.path.getsize( os.path.realpath(f.name) )
    def lock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_RLCK, file_size(f))
    def unlock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, file_size(f))


# Class for ensuring that all file operations are atomic, treat
# initialization like a standard call to 'open' that happens to be atomic
class AtomicOpen:
    # First acquire a lock for the file, then open the file with
    # arguments provided by user. Attempts to inherit most of the file
    # properties, but use "__enter__" to get file object directly.
    def __init__(self, path, *args, **kwargs):
        writing = False
        if (len(args) > 0) and (args[0]) == "w": 
            writing = True
            args = tuple(a if a != "w" else "a" for a in args)
        # Open the file and acquire a lock on the file before operating
        self.file = open(path,*args, **kwargs)
        # Lock the opened file
        lock_file(self.file)
        # If the user wanted to write (overwrite existing contents),
        # make sure that we are doing that exactly (after locking)
        if writing:
            self.file.seek(0)
            self.file.truncate()

    # Return the opened file object (knowing a lock has been obtained)
    def __enter__(self, *args, **kwargs): return self.file

    # Allows users to use the 'close' function if they want, in case
    # the user did not have the AtomicOpen in a "with" block.
    def close(self): self.__exit__()

    # Unlock the file and close the file object
    def __exit__(self, exc_type=None, exc_value=None, traceback=None):        
        # Release the lock on the file
        unlock_file(self.file)
        self.file.close()
        # Handle exceptions that may have come up during execution, by
        # default any exceptions are raised to the user
        if (exc_type != None): return False
        else:                  return True        



# class PathNotFoundError(Exception): pass

# import os
# import numpy as np

# class BigStruct:
#     def __init__(self, data=None):
#         if type(data) == str:
#             # The given data is a BigStruct stored at the given path
#             self.data_path = data
#             if not os.path.exists(data):
#                 raise(NotFoundError("'%s' does not exist."))
            

#         elif type(data) == np.ndarray:
#             # The given data is stored in a numpy (structured) array
#             pass
#         elif type(data) == list:
#             # The data is stored in a 2D list
#             pass
        
#     def save(self, path="", name="BigStruct"):
#         # Set the default path to the current directory
#         if len(path) == 0: path = os.curdir
#         # Add a number to the name if necessary
#         os.makedirs(os.path.join(path,name), exists_okay=False)
        
#     def get():
#         pass


# if __name__ == "__main__":
#     import pickle

#     size = 10**8
#     a = np.ones((size,), dtype=int)

#     # filename = "test.pkl"

#     # print("Writing", flush=True)
#     # with open(filename, "wb") as f:
#     #     pickle.dump(a,f)
#     # print("Done writing", flush=True)

#     # print("Reading", flush=True)
#     # with open(filename, "rb") as f:
#     #     b = pickle.load(f)
#     # print("Done reading", flush=True)

#     print("Reading..", flush=True)
#     with open("Data.pkl", "rb") as f:
#         data = pickle.load(f)
#     print("Done reading", flush=True)
#     print()
#     print(list(data.dtype.names))
#     print()
#     print("Slicing", flush=True)
#     data = data[data["Machine"] == "m1"]
#     data = data[data["Record_Size"] < 1000]
#     data = data[data["Frequency"] < 2800000]
#     print("Done slicing", flush=True)

