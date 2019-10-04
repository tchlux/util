# Decorators implemented in this file:
# 
#  same_as        -- Copy the function signature of another function.
# 
#  cache          -- Cache (input, output) pairs and use cache on
#                    repeats. Inputs and outputs must be pickle-able.
# 
#  stability_lock -- Generate cached test cases and automatically run
#                    tests when source code of function changes.
# 
#  timeout        -- Add python signal that times out function with default.
# 
#  type_check     -- Enforce input types for function.
# 
#  capture        -- Capture all standard error and standard output
#                    from most recent function call and store it in
#                    function attributes ".stdout" and ".stderr".
# 
#  background     -- Provides a decorator-like interface that runs a
#                    command in the background. Blocking wait for
#                    results is not performed until an attribute is
#                    accessed.

# ==================================================================
#                         SameAs Decorator
# 
# Decorator that copies the documentation and arguemnts of another
# function (specified as input). Useful for making decorators (:P)
# Optional "mention_usage" updates documentation when True to disclose
# the name of the function being wrapped (and the reuse of signature).
# 
# USAGE: 
# 
#   @same_as(<func_to_copy>)
#   def <function_to_decorate>(...):
#      ...
# 
#   OR
# 
#   <function> = same_as(<func_to_copy>)(<function_to_decorate>)
#   
def same_as(to_copy, mention_usage=False):
    import inspect
    # Create a function that takes one argument, a function to be
    # decorated. This will be called by python when decorating.
    def decorator_handler(func):
        if hasattr(func, "__name__"): original_name = func.__name__
        else:                         original_name = str(func)
        # Set the documentation string for this new function
        documentation = inspect.getdoc(to_copy)
        if documentation == None: 
            documentation = inspect.getcomments(to_copy)
        # Store the documentation and signature into the wrapped function
        if hasattr(to_copy, "__name__"):
            func.__name__ = to_copy.__name__
            if mention_usage:
                documentation = (
                    "\nThe function '%s' has been decorated with the signature "+
                    "of '%s'. (likely for aliasing / decoration)\n\n")%(
                        original_name, to_copy.__name__) + documentation
        # Try copying the signature if possible
        try:               func.__signature__ = inspect.signature(to_copy)
        except ValueError: pass
        # Finalize by copying documentation
        func.__doc__ = documentation
        return func
    # Return the decorator handler
    return decorator_handler


# ==================================================================
#                    "Cache in File" Decorator     
# 
# This decorator (when wrapped around a function) uses a hash of the
# string represenetation of the parameters to a function call in order
# to write a pickle file that contains the inputs to the function and
# the outputs to the function.
# 
# 
# USAGE: 
# 
#   @cache(<max_files>=10, <cache_dir>=os.curdir, <file_prefix>="Cache_[func.__name__]")
#   def <function_to_decorate>(...):
#      ...
# 
#   OR
# 
#   <function> = cache(<max_files>, <cache_dir>, <file_prefix>)(<function_to_decorate>)
#   
def cache(max_files=10, cache_dir=None, file_prefix=None):
    import os, hashlib
    # Import "dill" if it is available, otherwise use pickle.
    try:    import dill as pickle
    except: import pickle
    # Check to see if a cache directory was provided
    if (type(cache_dir) == type(None)): cache_dir = os.path.abspath(os.curdir)
    if (not os.path.exists(cache_dir)): os.mkdir(cache_dir)
    # Create a function that takes one argument, a function to be
    # decorated. This will be called by python when decorating.
    def decorator_handler(func):
        cache_prefix = file_prefix
        if (type(file_prefix) == type(None)):
            cache_prefix = "Cache_[%s]"%(func.__name__)
        def new_func(*args, **kwargs):
            # Identify a cache name via sha256 hex over the serialization
            hash_value = hashlib.sha256(pickle.dumps((args, kwargs))).hexdigest()
            cache_suffix = ".pkl"
            cache_path = os.path.join(cache_dir, cache_prefix+"_"+hash_value+cache_suffix)
            # Check to see if a matching cache file exists
            if os.path.exists(cache_path):
                with open(cache_path, "rb") as f:
                    args, kwargs, output = pickle.load(f)
            else:
                # Calculate the output normally with the function
                output = func(*args, **kwargs)
                # Identify the names of existing caches in this directory
                existing_caches = [f for f in os.listdir(cache_dir)
                                   if cache_prefix in f[:len(cache_prefix)]]
                # Only save a cached file if there are fewer than "max_files"
                if len(existing_caches) < max_files:
                    with open(cache_path, "wb") as f:
                        pickle.dump((args, kwargs, output), f)
            # Return the output (however it was achieved)
            return output
        # Return the decorated version of the function with identical documntation
        return same_as(func)(new_func)
    # Return a function that is capable of decorating
    return decorator_handler


# ==================================================================
#                    "Stability Lock" Decorator     
# 
# This decorator (when wrapped around a function) records serialized
# examples of (input, output) pairs as "test cases" to validate that
# the behavior of a function doesn't change over time. When the source
# code of the function is modified, stability checks are re-run.
# Failed test cases are neatly displayed to the user to demonstrate
# how the behavior of the function has changed with modifications.
# 
# 
# USAGE: 
# 
#   @stability_lock(<max_tests>=10, <test_dir>=".python_stability_locks", <auto_test>=True)
#   def <function_to_decorate>(...):
#      ...
# 
#   OR
# 
#   <function> = stability_lock(<max_tests>, <test_dir>, <auto_test>)(<function_to_decorate>)
#   
def stability_lock(max_tests=10, test_dir=None, auto_test=True):
    import os, inspect
    # Import "dill" if it is available, otherwise use pickle.
    try:    import dill as pickle
    except: import pickle
    # Create a custom Exception for printing failed checks.
    class FailedStabilityTest(Exception): pass
    # Check to see if a cache directory was provided, make default the
    # same file as the one this code resides in.
    if (type(test_dir) == type(None)): 
        user_home = os.path.expanduser("~")
        test_dir = os.path.join(user_home, ".python_stability_locks")
    # Make the directory if it does not exist already
    os.makedirs(test_dir, exist_ok=True)
    # Create a function that takes one argument, a function to be
    # decorated. This will be called by python when decorating.
    def decorator_handler(func):
        test_case_prefix = "Stability_Test_[%s]"%(func.__name__)
        # Decorate the given function with a cache
        cached_func = cache(max_files=max_tests, cache_dir=test_dir,
                            file_prefix=test_case_prefix)(func)
        def new_func(*args, **kwargs):
            # If auto_test is disabled, return the function value
            # without checking anything (saves execution time).
            if not auto_test: return(func(*args, **kwargs))
            # Generate a name and path for the function record file
            func_record_name = "[%s]_Record.pkl"%(func.__name__)
            func_record_path = os.path.join(test_dir, func_record_name)
            # Get the function record as the source python code
            curr_func_record = pickle.dumps(inspect.getsource(func))
            needs_checking = True
            if os.path.exists(func_record_path):
                with open(func_record_path, "rb") as f:
                    # This function needs testing if the source code has changed
                    needs_checking = (pickle.load(f) != curr_func_record)
            # If the function needs checking, then execute tests
            if needs_checking:
                print()
                print("Checking stability of '%s'..."%(func.__name__))
                print("  using test cases stored at '%s'"%(test_dir))
                # Get the paths to all existing test cases
                test_file_paths = [os.path.join(test_dir, f_name)
                                   for f_name in os.listdir(test_dir)
                                   if test_case_prefix in f_name]
                for test_path in test_file_paths:
                    with open(test_path, "rb") as f:
                        (args, kwargs, expected_output) = pickle.load(f)
                    print("  Testing on '%s'..."%(os.path.basename(test_path)), end="\r")
                    actual_output = func(*args, **kwargs)
                    # If a stability check failed, generate informative output
                    if (pickle.dumps(expected_output) != pickle.dumps(actual_output)):
                        failure_message = "\n\nFailed stability test case '%s'\n\n"%(os.path.basename(test_path))
                        failure_message += "INPUTS\n"
                        for arg_num, arg in enumerate(args):
                            failure_message += "  Argument %i:\n"%(arg_num)
                            failure_message += "    "+str(arg)+"\n\n"
                        for kwarg_key in kwargs:
                            failure_message += "  Keyword argument '%s':\n"%(kwarg_key)
                            failure_message += "    "+str(kwargs[kwarg_key])+"\n\n"
                        failure_message += "OUTPUTS\n"
                        failure_message += "  Expected output:\n"
                        failure_message += "    "+str(expected_output)+"\n\n"
                        failure_message += "  Actual output:\n"
                        failure_message += "    "+str(actual_output)+"\n"
                        # Raise the failure as an exception
                        raise(FailedStabilityTest(failure_message))
                    else:
                        print("  Passed '%s'"%(os.path.basename(test_path)))
                # Now that tests are complete and have passed, update the record
                with open(func_record_path, "wb") as f:
                    pickle.dump(curr_func_record, f)
                print("No stability tests failed, updated record file '%s'."%(
                    os.path.basename(func_record_path)))
                print()
            # If none of the test cases available failed, execute the
            # function with a caching mechanism activated.
            return cached_func(*args, **kwargs)
        # Return the decorated version of the function with identical documntation
        return same_as(func)(new_func)
    # Return a function that is capable of decorating
    return decorator_handler


# ==================================================================
#                         Timeout Decorator     
# 
# Function decorator that uses the "signal" module in order to ensure
# that a function call does not take more than "allowed_seconds" to
# complete. If the function takes longer than "allowed_seconds" then
# the value "default" is returned.
# 
# USAGE: 
# 
#   @timeout(<allowed_seconds>, <default>)
#   def <function_to_decorate>(...):
#      ...
# 
#   OR
# 
#   <function> = timeout(<sec>, <default>)(<function_to_decorate>)
#   
def timeout(allowed_seconds=1, default=None):
    import signal, inspect
    # Create a custom exception, be sure to *only* catch it
    class TimeoutError(Exception): pass
    # Create a function that takes one argument, a function to be
    # decorated. This will be called by python when decorating.
    def decorator_handler(func):
        # Create a new function (for return when decorating)
        def new_func(*args, **kwargs):
            # Send our exception as a signal if this handler is called
            def handler(signum, frame):
                raise TimeoutError()
            # Set our handler function to be called by the alarm signal
            signal.signal(signal.SIGALRM, handler) 
            # Signal an alarm (our handler) if "allowed_seconds" elapses
            signal.alarm(round(allowed_seconds+0.5))
            try:
                # Try running the decorated function
                result = func(*args, **kwargs)
            except TimeoutError:
                # The decorated function did not finish
                result = default
            finally:
                # Cancel the alarm if the <try> completes successfully
                signal.alarm(0)
            return result
        # Set the documentation string for this new function
        return same_as(func)(new_func)
    # Check which type of usage is being implemented
    if type(allowed_seconds) == type(lambda:None):
        # If the user has applied the decorator without arguments,
        #  then the first argument provided is actually "func"
        func = allowed_seconds
        # Get the default value for "allowed_seconds"
        allowed_seconds = inspect.signature(timeout).parameters['allowed_seconds'].default
        decorated_func = decorator_handler(func)
    else:
        # The user must have given some arguments, normal decoration
        decorated_func = decorator_handler
    # Return the function that will return the decorated function
    return decorated_func


# ==================================================================
#                       Type Check Decorator     
# 
# TODO: The 'list' arg_type should do checking along each of the
#       dimensions of the input, and should only require that the
#       input be iterable up to D dimensions deep. Add associated
#       errros as well.
# TODO: Clean up the arg type function to only reference arguments by
#       name in the error messages, no need to duplicate code.
# 
# Function decorator that checks the types of arguments and keyword
# arguments. It's for guaranteeing proper usage (forgive how
# un-pythonic this is, but sometimes it's useful for producing better
# error messages). Given a list of arguments that is shorter than the
# actual list, it only checks the first N provided. Given keyword
# arguments to check, it only checks them if they are
# provided. Returns a function with a for loop over the arguments
# checking for minimum number of arguments and correct types. Also
# transfers the documentation to the decorated function.
# 
# USAGE:
# 
#   @type_check([<arg_type>, ...]=[], {<arg_name>:<arg_type>, ...}={})
#   def <func_name>( ... ):
# 
# OR
# 
#   <function> = type_check(<type_check_args>)(<function_to_decorate>)
# 
# "arg_type" can be one of the following:
#    type     -> Type checking is equality based         "type(arg) == arg_type"
#    set      -> Type checking is membership based       "type(arg) in arg_type"
#    list     -> Type checking is nested                 "(type(arg) == list) and all((type(arg_i) == t) for (arg_i,t) in zip(arg,arg_type))"
#    function -> Type checking is provided by function   "arg_type(arg)"
# 
def type_check(*types, **types_by_name):
    # Errors for type checking
    class WrongType(Exception): pass
    class WrongType_NotInSet(Exception): pass
    class WrongType_NotListed(Exception): pass
    class WrongType_FailedCheck(Exception): pass
    class WrongNumberOfArguments(Exception): pass
    class WrongUsageOfTypeCheck(Exception): pass
    FUNCTION_TYPE = type(lambda:None)
    # Check for poor usage of the type-checking function (potentially bad types XD)
    for t in types + tuple(types_by_name.values()):
        if type(t) not in {type, set, list, FUNCTION_TYPE}:
            raise(WrongUsageOfTypeCheck(
                "Type checking can only handle types, sets of "+
                "types, and functions that return types. Type '%s' not allowed."%(type(t))))

    # Create a function that takes one argument, a function to be
    # decorated. This will be called by python when decorating.
    def decorator_handler(func):
        # Create a new function (for return when decorating)
        def new_func(*args, **kwargs):
            # Make sure the correct number of arguments were given
            if len(args) < len(types):
                raise(WrongNumberOfArguments("Type checked function expected %i arguments."%(len(types))))
            # Check types of all of the regular arguments
            for i,(a, t) in enumerate(zip(args, types)):
                # Check if type equals the expected type
                if type(t) == type:
                    if not (type(a) == t):
                        raise(WrongType(("Expected argument %i to be type '%s',"+
                                         "received type '%s' instead.")%(i,t,type(a))))
                # Check if type exists in a set of allowed types
                elif type(t) == set:
                    if not (type(a) in t):
                        raise(WrongType_NotInSet(("Expected argument %i to be one of types %s,"+
                                                  "received type '%s' instead.")%(i,t,type(a))))
                # Check if list argument is sub-typed correctly
                elif type(t) == list:
                    # Check to make sure that the provided value is index accessible (not a generator).
                    try:
                        for i in range(len(a)): a[i]
                        for v in a: break
                    except TypeError:
                        raise(WrongType(("Expected argument %i to be iterable and be index"+
                                         "accessible, type '%s' does not meet these criteria.")%(i,type(a))))
                    # Check for subtypes of list
                    if not all(type(ai) == ti for (ai, ti) in zip(a,t)):
                        raise(WrongType_NotListed((
                            "Types contained in iterable argument %i did not match "+
                            "match expected type listing %s.")%(i,t)))
                # Check if type passes a type-checking function
                elif type(t) == FUNCTION_TYPE:
                    if not t(a):
                        raise(WrongType_FailedCheck(("Argument %i of type '%s' did not pass "+
                                                     "required type checking function.")%(i,type(a))))
            # Check types of all of the keyword arguments
            for arg_name in kwargs:
                if arg_name in types_by_name:
                    a, t = kwargs[arg_name], types_by_name[arg_name]
                    # Check if type equals the expected type
                    if type(t) == type:
                        if not (type(a) == t):
                            raise(WrongType(("Expected keyword argument '%s' to be type '%s',"+
                                             "received type '%s' instead.")%(arg_name,t,type(a))))
                    # Check if type exists in a set of allowed types
                    elif type(t) == set:
                        if not (type(a) in t):
                            raise(WrongType_NotInSet(("Expected keyword argument '%s' to be one of types %s,"+
                                                      "received type '%s' instead.")%(arg_name,t,type(a))))
                    # Check if list argument is sub-typed correctly
                    elif type(t) == list:
                        # Check for list type
                        if not (type(a) == list):
                            raise(WrongType(("Expected keyword argument '%s' to be type '%s',"+
                                             "received type '%s' instead.")%(arg_name,list,type(a))))
                        # Check for subtypes of list
                        if not all(type(ai) == ti for (ai, ti) in zip(a,t)):
                            raise(WrongType_NotListed((
                                "Types contained in list keyword argument '%s' did not match "+
                                "match expected type listing %s.")%(arg_name,t)))
                    # Check if type passes a type-checking function
                    elif type(t) == type(lambda:None):
                        if not t(a):
                            raise(WrongType_FailedCheck(("Keyword argument for '%s' of type '%s' did not pass "+
                                                         "required type checking function.")%(arg_name,type(a))))
            return func(*args, **kwargs)
        # Set the documentation string for this new function
        return same_as(func)(new_func)
    # Return the function that will return the decorated function
    return decorator_handler


# ==================================================================
#                         Capture Decorator
# 
# Decorator that captures all of the called functions output to
# standard error and standard output. They are stored in attributes of
# the function named "stderr" and "stdout" respectively.
# 
# USAGE: 
# 
#   @capture
#   def <function_to_decorate>(...):
#      ...
# 
#   OR
# 
#   <function> = capture(<function_to_decorate>)
#   
# Define the capture function.
def capture(func):
    import tempfile, sys
    import multiprocessing.connection
    try:    import dill as pickle
    except: import pickle

    # Get a "run" command and make it execute non-blocking.
    from util.system import run
    run = background(run)
    # Custom exception.
    class FailedConnection(Exception): pass
    # Get the capture worker program template from file.
    from util.decorators import CAPTURE_WORKER
    with open(CAPTURE_WORKER, "r") as f: capture_worker_program = f.read()

    # Retrieve the first available port in a range by tring to establish
    # listeners on different ports.
    def get_port(key, start_port=10000, attempts=1000):
        import socket
        import multiprocessing.connection        
        # Loop through potential ports.
        for i in range(start_port, start_port+attempts):
            try:
                port = i
                listener = multiprocessing.connection.Listener(
                    ('localhost', port), authkey=key)
                break
            # If the port is unavailable, keep trying.
            except Exception as ex:
                # 48 -- Port already in use.
                if ex.errno not in {48,}: raise(ex)
        else:
            raise(FailedConnection("Could not establish a listener."))
        # Free the memory (and the port since we know it's available).
        del(listener)
        return port

    stdout = []
    stderr = []
    # Create a function that runs in a sub-process and captures all output.
    def captured_func(*args, **kwargs):
        # Generat an authorization key
        authorization_str = 'jaaodfivjnlaeoi183hn08402jlkad'
        authorization = bytes(authorization_str, "ascii")
        # Get the ports for listening and sending.
        port_in = get_port(authorization)
        listener = multiprocessing.connection.Listener(('localhost', port_in), authkey=authorization)
        port_out = get_port(authorization)
        # Start child process (it should create a listener on other port) then
        # a client that needs this listener to accept it.
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as prog_file:
            temp_file_name = prog_file.name
            prog_file.write(capture_worker_program.format(
                authorization=authorization_str,
                port_in=port_out, port_out=port_in))
            # "run" the program, non-blocking in the background.
            command_output = run([sys.executable, temp_file_name])
        # Accept the connection that was requested by client.
        listener = listener.accept()
        # Request a connection to the other process.
        sender = multiprocessing.connection.Client(('localhost', port_out), authkey=authorization)
        # Send the "globals", the "func", and the "args, kwargs" pair for "func".
        sender.send_bytes(pickle.dumps(globals()))
        sender.send_bytes(pickle.dumps(func))
        sender.send_bytes(pickle.dumps((args, kwargs)))
        output = pickle.loads(listener.recv_bytes())
        # Clear the global variables used to store "stdout" and "stderr".
        stdout.clear()
        stderr.clear()
        # Retrieve the result from the "run" process executed in the background.
        code, stdout_lines, stderr_lines = command_output.result
        # Terminate any lingering processes.
        run.terminate()
        # Update the contents of stdout and stderr.
        for l in stdout_lines: stdout.append(l)
        for l in stderr_lines: stderr.append(l)
        # Return the output to the caller.
        return output
    # Add attributes to the function that hold stdout and stderr.
    setattr(captured_func, "stdout", stdout)
    setattr(captured_func, "stderr", stderr)
    # Return the captured function.
    return captured_func


# ==================================================================
#                      Background Decorator
# 
# This class provides a decorator-like interface that runs a command
# in the background. This uses a separate process so any global
# variable changes made by the funtion will not be propagated back to
# the parent process!
# 
# Once the output wants to be retrieved, use the ".result" attribute.
# 
# USAGE: 
# 
#   @background
#   def <function_to_decorate>(...):
#      ...
# 
#   OR
# 
#   <function> = background(<function_to_decorate>)
#   
# Results from the function when called will be returned in a
# 'Results' object that behaves (to a reasonable extent) like the
# actual underlying return value. If the user wants actual results
# then that must be retrieved with '<function output>.result'.
# 
# First define a "result" data-descriptor for returned results from
# the 'background' decorator. This will allow for many basic operators
# to be used on the result as if it is literally the result.
# 
class Result:
    # Initialize with a record of the parent and job number.
    def __init__(self, parent, job_num):
        super().__setattr__("_parent", parent)
        super().__setattr__("_job_num", job_num)
        super().__setattr__("_retrieved", False)
        super().__setattr__("_result", None)
    # Overwrite all attribute access to return the "result" attributes.
    def __getattr__(self, attr_name=""):
        return getattr(self.result, attr_name)
    # Make the retrieval of "result" blocking.
    @property
    def result(self):
        retrieved = super().__getattribute__('_retrieved')
        if (not retrieved):
            parent  = super().__getattribute__('_parent')
            job_num = super().__getattribute__('_job_num')
            result = parent.retrieve(job_num)
            super().__setattr__("_result", result)
            super().__setattr__("_retrieved", True)
        return super().__getattribute__('_result')
# Modify the "Result" class to behave like the "value" it contains.
_ = {"__init__", "__getattr__", "__getattribute__", "__class__",
     "__new__", "__dict__", "_retrieved", "_result", "_parent",
     "_job_num", "result"}
_ = [n for n in dir(Result) if (n not in _)]
# Create a method generating function.
def method_gen(name):
    def method(*args, **kwargs):
        self = args[0]
        return self.__getattr__(name)(*args[1:], **kwargs)
    return method
# Cycle through and overwrite all the standard methods.
for _ in _: setattr(Result, _, method_gen(_))
del(_)
# 
# 
# Beginning of actual 'Background' deocrator class definition.
class background:                              
    _job_count = 0
    _res_count = 0
    _completed = None
    _results = None
    _jobs = None

    # Define a function for running a job.
    def _run_job(q, job_num, func, args, kwargs):
        q.put( (job_num, background.loads(func)(*args, **kwargs)) )

    # Define custom exceptions.
    class NotCalled(Exception): pass
    class UsageError(Exception): pass

    # Import methods to dump and load functions to pass between processes.
    from dill import dumps, loads

    # Define a property (with no setter) for looking at the function.
    @property
    def func(self): return background.loads(self._func)
    @property
    def completed(self): return len(self._completed)
    @property
    def running(self): return len(self._jobs)
    @property
    def pending(self): return len(self._results)

    # Define this so that it behaves like a decorator.
    def __init__(self, func):
        from multiprocessing import get_context
        # Initialize storage for all of the details about this function.
        self._completed = []
        self._jobs = {}
        self._results = {}
        # Store the function in a serialized format.
        self._func = background.dumps(func)
        # Create multiprocessing context (not to muddle with others)
        self._ctx = get_context() # "fork" on Linux, "spawn" on Windows.
        # Create a queue for receiving function output.
        self._queue = self._ctx.Queue(1)
        # Register this background function for "cleanup" on python exit.
        import atexit
        atexit.register(self.terminate)

    # This is the method that's supposed to "call" the decorated function.
    def __call__(self, *args, **kwargs):
        # Serialize all the arguments to pass them to the subprocess.
        arguments = (self._queue, self._job_count, self._func, args, kwargs)
        self._jobs[self._job_count] = self._ctx.Process(
            target=background._run_job, args=arguments)
        # Start the job in the background.
        self._jobs[self._job_count].start()
        # Increment the job counter.
        self._job_count += 1
        # Return a Result object.
        return Result(self, self._job_count-1)

    # Print out information about this background function.
    def __str__(self):
        try:    name = "'" + self.func.__name__ + "'"
        except: name = ""
        string = f"Background function {name}\n"
        string += f"  completed jobs:  {self.completed}\n"
        string += f"  running jobs:    {self.running}\n"
        string += f"  pending results: {self.pending}\n"
        string += f" Get background function description with 'help(<this>.func)'.\n"
        string += f" All results are wrapped 'Result' instances that can be\n"
        string += f"  unpacked by accessing '<output>.result'."
        return string

    # If del is called, clean up all processes and terminate.
    def terminate(self):
        names = list(self._jobs.keys())
        for n in names: self._jobs.pop(n).terminate()

    # Define a "get" so that once the value is retrieved the process
    # is waited on for results.
    def retrieve(self, job):
        if (self._job_count == 0): raise(background.NotCalled(
                "This function has not been called yet."))
        if (job in self._completed): raise(background.UsageError(
                f"The job {job} has already been retrieved and deleted from this background instance."))
        # Continue getting jobs until the desired one is retrieved.
        while (job not in self._results):
            job_num, result = self._queue.get()
            self._results[job_num] = result
            # If an exception was encountered, terminate and raise.
            if isinstance(result, Exception):
                self.terminate()
                raise(value)
            # Join the job (end process) since it must be done.
            self._jobs.pop(job_num).join()
        # Increase the counter (the next result has been found and returned).
        self._completed.append(job)
        # Return the next result to the caller.
        return self._results.pop(job)

