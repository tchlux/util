import signal, inspect

# Function decorator that uses the "signal" module in order to ensure
# that a function call does not take more than "allowed_seconds" to
# complete. If the function takes longer than "allowed_seconds" then
# the value "default" is returned.
# 
# USAGE: 
# 
#   @timeout_after(<allowed_seconds>, <default>)
#   def <function_to_decorate>(...):
#      ...
# 
#   OR
# 
#   <function> = timeout_after(<sec>, <default>)(<function_to_decorate>)
#   
def timeout_after(allowed_seconds=1, default=None):
    # Create a function that takes one argument, a function to be
    # decorated. This will be called by python when decorating.
    def decorator_handler(func):
        # Create a new function (for return when decorating)
        def new_func(*args, **kwargs):
            # Create a custom exception, be sure to *only* catch it
            class TimeoutError(Exception): pass
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
        documentation = inspect.getdoc(func)
        if documentation == None: 
            comments = inspect.getcomments(func)
            if comments != None:
                documentation = "\n".join([line.lstrip("#") for line in
                                           comments.split("\n")])
            else:
                documentation = ""
        new_func.__doc__ = documentation
        return new_func
    # Check which type of usage is being implemented
    if type(allowed_seconds) == type(lambda:None):
        # If the user has applied the decorator without arguments,
        #  then the first argument provided is actually "func"
        func = allowed_seconds
        # Get the default value for "allowed_seconds"
        allowed_seconds = inspect.signature(timeout_after).parameters['allowed_seconds'].default
        decorated_func = decorator_handler(func)
    else:
        # The user must have given some arguments, normal decoration
        decorated_func = decorator_handler
    # Return the function that will return the decorated function
    return decorated_func

# This is a testing function for the timeout_after decorator. These
# comments will still be included in the documentation for the function.
@timeout_after
def test_timeout_after(sec):
    import time
    print("Going to sleep for '%.1f' seconds."%(sec))
    time.sleep(sec)
    print("Finished sleeping!")

