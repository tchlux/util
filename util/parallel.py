# Get functions required for establishing multiprocessing
from multiprocessing import cpu_count, get_context, current_process
# Imports for each process
from queue import Empty
import sys

# Make sure new processes started are minimal (not complete copies)
MAX_PROCS = cpu_count()
PROC_TIMEOUT = None # Should be adjusted per use
JOB_GET_TIMEOUT = 1

# ==============================================================
#      Dynamic Distribution multi-processing with Arguments     
# ==============================================================

# Pre:  "job_queue" is a Queue full of tuples that are arguments to
#       "func", which is a standard (non-lambda) function
#       "return_queue" is the Queue that should be filled with the
#       return values from each call to "func"
# Post: "return_queue" is filled with return values from each call
#       made to "func" with a set of arguments from "job_queue"
def worker(job_queue, return_queue, func):
    # Set the output file for this process so that all print statments
    # by default go there instead of to the terminal
    log_file_name = "Process-"+current_process().name.split("-")[1]
    sys.stdout = open("%s.log"%log_file_name,"w")
    returns = []
    while not job_queue.empty():
        try:
            args = job_queue.get(timeout=JOB_GET_TIMEOUT)
            return_queue.put(func(*args))
        except Empty: # See imported exceptions
            print("Failed 'get', job_queue empty.")

# Pre:  "func" is a standard function (NOT A LAMBDA)
#       "args_list" is a list of tuples where each tuple is a valid
#       set of arguments for "func"
# Post: "func" is run for all of the argument instances in "args_list"
#       distributed among however many processors are on the current
#       computer using a dynamic shared jobs queue
def call(func, args_list, verbose=True):
    # Create multiprocessing context (not to muddle with others)
    ctx = get_context("spawn")    
    # Create queues
    job_queue = ctx.Queue()
    return_queue = ctx.Queue()
    proc_args = (job_queue, return_queue, func)
    processes = [ctx.Process( target=worker, args=proc_args )
                 for i in range(MAX_PROCS)]    
    # Fill the job queue
    for arg in args_list:
        job_queue.put(arg)
    # Start the processes
    for p in processes:
        p.start()
    # Continue picking up results until completed
    completed = 0
    returns = []
    if verbose: print("Distributing jobs among %i processors..."%MAX_PROCS)
    if verbose: print("[%0.1f] %i completed of %i total."%
                      (100.0 * completed / len(args_list), 
                       completed, len(args_list)), end="")
    while completed != len(args_list):
        returns.append( return_queue.get() )
        completed += 1
        if verbose:
            print("\r[%0.1f] %i completed of %i total."%
                  (100.0 * completed / len(args_list), 
                   completed, len(args_list)), end="")
    if verbose: print()
    for p in processes: p.join()
    for p in processes: p.terminate()
    # Return all results
    return returns
