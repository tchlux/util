# Get functions required for establishing multiprocessing
from multiprocessing import cpu_count, get_context, current_process
# Imports that are used by the processes.
from queue import Empty
import os
# Try to import "dill" module for better serialization, otherwise use pickle.
try:                        from dill   import dumps, loads
except ModuleNotFoundError: from pickle import dumps, loads

# Make sure new processes started are minimal (not complete copies)
MAX_PROCS = cpu_count()
LOG_DIRECTORY = os.path.abspath(os.path.curdir)
JOB_GET_TIMEOUT = 5
# Global storage for processes
MAP_PROCESSES = []
# Get a reference to the builtin "map" function, so it is not lost
# when it is overwrittend by "def map" later in this file.
builtin_map = map

# Iterator that pulls jobs from a queue and stops iterating once the
# queue is empty (and remains empty for "JOB_GET_TIMEOUT" seconds). On
# construction, requires a "job_queue" that jobs will come from.
class JobIterator:
    def __init__(self, job_queue):
        self.job_queue = job_queue
    def __iter__(self): return self
    def __next__(self):
        try:          value = self.job_queue.get(timeout=JOB_GET_TIMEOUT)
        except Empty: raise(StopIteration)
        return value

# Given an iterable object and a queue, place jobs into the queue.
def producer(jobs_iter, job_queue):
    for job in enumerate(jobs_iter):
        job_queue.put(job)

# Given an iterable object, a queue to place return values, and a
# function, apply "func" to elements of "iterable" and put return
# values into "return_queue".
# 
# INPUTS:
#   func         -- A function that takes elements of "iterable" as
#                   input. It is expected that the function has been
#                   packaged with with "dill.dumps" to handle lambda's.
#   iterable     -- An iterable object (presumably a JobIterator).
#   return_queue -- A Queue object with a 'put' method.
#   redirect     -- Boolean, True if log files should be created per
#                   process, otherwise all print goes to standard out.
def consumer(func, iterable, return_queue, redirect):
    import sys
    # Retrieve the function (because it's been dilled for transfer)
    func = loads(func)
    if redirect:
        # Set the output file for this process so that all print statments
        # by default go there instead of to the terminal
        log_file_name = f"Process-{int(current_process().name.split('-')[1])-1}"
        log_file_name = os.path.join(LOG_DIRECTORY, log_file_name)
        sys.stdout = open("%s.log"%log_file_name,"w")
    # Iterate over values and apply function.
    for (i,value) in iterable:
        try:
            # Try executing the function on the value.
            return_queue.put((i,func(value)))
        except Exception as exception:
            # If there is an exception, put it into the queue for the parent.
            return_queue.put((i,exception))
            break
    # Once the iterable has been exhaused, put in a "StopIteration".
    return_queue.put((-1, StopIteration))
    # Close the file object for output redirection
    if redirect: sys.stdout.close()
    

# A parallel implementation of the builtin `map` function in
# Python. Provided a function and an iterable that is *not* a
# generator, use all the cores on the current system to map "func" to
# all elements of "iterable".
# 
# INPUTS
#   func     -- Any function that takes elements of "iterable" as an argument.
#   iterable -- An iterable object The elements of this will be passed into "func".
# 
# OPTIONALS:
#   max_waiting_jobs    -- Integer, maximum number of jobs that can be
#                          in the job queue at any given time.
#   max_waiting_returns -- Integer, maximum number of return values in
#                          the return queue at any given time.
#   order     -- Boolean, True if returns should be made in same order
#                as the originally generated sequence.
#                WARNING: Could potentially hold up to len(iterable)
#                return values from "func" in memory.
#   save_logs -- Boolean, True if log files from individual calls to
#                "func" should *not* be deleted.
#   redirect  -- Boolean, True if log files should be created per
#                process, otherwise all print goes to standard out.
#   chunk     -- Boolean, True if `split` should be called on the
#                iterable automatically. WARNING: This will convert
#                the iterable into a python list, more memory!
# 
# RETURNS:
#   An output generator, just like the builtin "map" function.
# 
# WARNINGS:
#   If iteration is not completed, then waiting idle processes will
#   still be alive. Call "parallel.killall()" to terminate lingering
#   map processes when iteration is prematurely terminated.
# 
#   The mapped function will not have access to global variables.
def map(func, iterable, max_waiting_jobs=1, max_waiting_returns=1,
        order=True, save_logs=False, redirect=True, chunk=True):
    # If chunking should be done, then `split` the iterable.
    if chunk:
        iterable = split(iterable)
        total_size = sum([len(i) for i in iterable])
        single_func = func
        # Create a new function that returns a list of the original
        # function mapped onto the elements.
        func = lambda series: [single_func(element) for element in series]
    # Create multiprocessing context (not to muddle with others)
    ctx = get_context() # "fork" on Linux, "spawn" on Windows.
    # Create a job_queue for passing jobs to processes
    job_queue = ctx.Queue(max_waiting_jobs)
    # Create return queue for holding results
    return_queue = ctx.Queue(max_waiting_returns)
    # Create the shared set of arguments going to all consumer processes
    producer_args = (iterable, job_queue)
    consumer_args = (dumps(func), JobIterator(job_queue), return_queue, redirect)
    # Start a "producer" process that will add jobs to the queue,
    # create the set of "consumer" processes that will get jobs from the queue.
    producer_process = ctx.Process(target=producer, args=producer_args)
    consumer_processes = [ctx.Process(target=consumer, args=consumer_args)
                          for i in range(MAX_PROCS)]
    # Track (in a global) all of the active processes.
    MAP_PROCESSES.append(producer_process)
    for p in consumer_processes: MAP_PROCESSES.append(p)
    # Start the producer and consumer processes
    producer_process.start()
    for p in consumer_processes: p.start()
    # Construct a generator function for reading the returns
    def return_generator():
        # Set up variables for returning items in order if desired.
        if order:
            idx = 0
            ready_values = {}
        # Keep track of the number of finished processes there are.
        stopped_consumers = 0
        # Continue picking up results until generator is depleted
        while (stopped_consumers < MAX_PROCS):
            i, value = return_queue.get()
            if (value == StopIteration):
                stopped_consumers += 1
            elif isinstance(value, Exception):
                # Kill all processes for this "map"
                for p in consumer_processes: 
                    p.terminate()
                producer_process.terminate()
                # Raise exception without deleting log files.
                raise(value)
            elif (not order): yield value                
            else: ready_values[i] = value
            # Yield all waiting values in order (if there are any).
            if order:
                while (idx in ready_values):
                    yield ready_values.pop(idx)
                    idx += 1
        # Join all processes (to close them gracefully).
        for p in consumer_processes:
            try: p.join()
            except AssertionError: pass
            MAP_PROCESSES.remove(p)
        try: producer_process.join()
        except AssertionError: pass
        MAP_PROCESSES.remove(producer_process)
        # Delete log files if the user doess not want them.
        if not save_logs: clear_logs()
    # If this data was chunked, then re-merge the data in order.
    if chunk:
        result = list(return_generator())
        output = []
        for i in range(total_size):
            output.append( result[i%len(result)].pop(0) )
        return output
    # Return all results in a generator object.
    return return_generator()

# Given an iterable, split it into as many chunks as there are
# processes. This is not particularly efficient, as it uses python
# lists, but it is convenient!
def split(iterable, count=MAX_PROCS):
    # Fill the chunks by walking through the iterable.
    chunks = [list() for i in range(count)]
    for i,v in enumerate(iterable):
        chunks[i%count].append( v )
    # Return the list of lists that are chunks.
    return chunks

# Kill all active "map" processes.
def killall():
    # Terimate all processes and empty global record.
    while len(MAP_PROCESSES) > 0:
        proc = MAP_PROCESSES.pop(-1)
        proc.terminate()

# Clear out the log files from each process.
def clear_logs():
    for f_name in os.listdir(LOG_DIRECTORY):
        if (f_name[:8] == "Process-") and (f_name[-4:] == ".log"):
            f_name = os.path.join(LOG_DIRECTORY, f_name)
            os.remove(f_name)

# Given a shape, create a numpy double array that is prepared to be
# shared during multiprocessing, but behaves like a normal array.
def shared_array(shape):
    import numpy as np
    from multiprocessing.sharedctypes import Array
    from ctypes import c_double
    # Allocate the memory in shared space.
    memory = Array(c_double, int(np.prod(shape)))
    # Create and return a structure to access the shared memory (numpy array).
    return np.frombuffer(memory.get_obj(), dtype=float).reshape(shape)
