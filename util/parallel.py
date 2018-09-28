# Get functions required for establishing multiprocessing
from multiprocessing import cpu_count, get_context, current_process
# Imports that are used by the processes.
from queue import Empty
import sys
from dill import dumps, loads

# Make sure new processes started are minimal (not complete copies)
MAX_PROCS = cpu_count()
JOB_GET_TIMEOUT = 1
# Global storage for processes
MAP_PROCESSES = []
# Get a reference to the builtin "map" function, so it is not lost
# when it is overwrittend by "def map" later in this file.
builtin_map = map

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

# Iterator that pulls jobs from a queue and stops iterating once the
# queue is empty (and remains empty for "JOB_GET_TIMEOUT" seconds). On
# construction, requires a "job_queue" that jobs will come from.
class JobIterator:
    def __init__(self, job_queue):
        self.job_queue = job_queue
    def __iter__(self): return self
    def __next__(self):
        try:
            value = self.job_queue.get(timeout=JOB_GET_TIMEOUT)
        except Empty:
            raise(StopIteration)
        return value

# Given an iterable object and a queue, place jobs into the queue.
def producer(jobs_iter, job_queue):
    print("PARALLEL: Starting producer process..", flush=True)
    for job in jobs_iter:
        print("PARALLEL: Producer placing job in queue", flush=True)
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
def consumer(func, iterable, return_queue):
    print("PARALLEL: Starting consumer process..", flush=True)
    # Retreive the function (because it's been dilled for transfer)
    func = loads(func)
    # Set the output file for this process so that all print statments
    # by default go there instead of to the terminal
    log_file_name = f"Process-{int(current_process().name.split('-')[1])-1}"
    # sys.stdout = open("%s.log"%log_file_name,"w")
    # Iterate over values and apply function.
    for value in iterable:
        try:
            print("PARALLEL: About to evaluate function..")
            # Try executing the function on the value.
            return_queue.put(func(value))
        except Exception as exception:
            # If there is an exception, put it into the queue for the parent.
            return_queue.put(exception)
            break
    # Once the iterable has been exhaused, put in a "StopIteration".
    return_queue.put(StopIteration)
    # Close the file object for output redirection
    # sys.stdout.close()


# A parallel OUT-OF-ORDER implementation of the builtin 'map' function
# in Python. Provided a function and an iterable that is *not* a
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
# 
# RETURNS:
#   An output generator, just like the builtin "map" function, but out-of-order.
# 
# WARNINGS:
#   If iteration is not completed, then waiting idle processes will
#   still be alive. Call "parallel.killall()" to terminate lingering
#   map processes when iteration is prematurely terminated.
def map(func, iterable, max_waiting_jobs=1, max_waiting_returns=1):
    # Create multiprocessing context (not to muddle with others)
    ctx = get_context() # "fork" on Linux, "spawn" on Windows.
    # Create a job_queue for passing jobs to processes
    job_queue = ctx.Queue(max_waiting_jobs)
    # Create return queue for holding results
    return_queue = ctx.Queue(max_waiting_returns)
    # Create the shared set of arguments going to all consumer processes
    producer_args = (iterable, job_queue)
    consumer_args = (dumps(func), JobIterator(job_queue), return_queue)
    # Start a "producer" process that will add jobs to the queue,
    # create the set of "consumer" processes that will get jobs from the queue.
    producer_process = ctx.Process(target=producer, args=producer_args)
    consumer_processes = [ctx.Process(target=consumer, args=consumer_args)
                          for i in range(MAX_PROCS)]
    # Track (in a global) all of the active processes.
    MAP_PROCESSES.append(producer_process)
    for p in consumer_processes:
        MAP_PROCESSES.append(p)
    # Start the producer and consumer processes
    producer_process.start()
    for p in consumer_processes: p.start()
    print("PARALLEL: Constructing a generator for the processes..", flush=True)
    # Construct a generator function for reading the returns
    def return_generator():
        stopped_consumers = 0
        # Continue picking up results until generator is depleted
        while (stopped_consumers < MAX_PROCS):
            value = return_queue.get()
            if (value == StopIteration):
                stopped_consumers += 1
            elif isinstance(value, Exception):
                # Kill all processes for this "map"
                for p in consumer_processes: 
                    p.terminate()
                producer_process.terminate()
                # Raise exception
                raise(value)
            else:
                yield value
        # Join all processes (to close them gracefully).
        for p in consumer_processes:
            p.join()
            MAP_PROCESSES.remove(p)
        producer_process.join()
        MAP_PROCESSES.remove(producer_process)
    print("PARALLEL: Returning the generator..", flush=True)
    # Return all results in a generator object
    return return_generator()

# Kill all active "map" processes.
def killall():
    # Terimate all processes and empty global record.
    while len(MAP_PROCESSES) > 0:
        proc = MAP_PROCESSES.pop(-1)
        proc.terminate()


# Development testing code.
if __name__ == "__main__":
    import time, random
    # A slow function
    def slow_func(value):
        import time, random
        time.sleep(random.random()*.1 + .5)
        return str(value)
    # Generate a normal "map" and a parallel "map"
    start = time.time()
    s_gen = builtin_map(slow_func, [i for i in range(10)])
    print("Standard map construction time:", time.time() - start)
    start = time.time()
    p_gen = map(slow_func, [i for i in range(10)])
    print("Parllel map construction time: ", time.time() - start)

    start = time.time()
    print()
    for value in s_gen:
        print(value, end=" ", flush=True)
        if value >= "5": break
    print()
    print("Standard map time:", time.time() - start)

    start = time.time()
    print()
    for value in p_gen:
        print(value, end=" ", flush=True)
        if value > "5": break
    print()
    print("Parallel map time:", time.time() - start)

    killall()
