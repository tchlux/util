from util.decorators import *

def test_background():
    # Make a "result" object that gets a specific result from the
    # parent function. Make that object callable and give it good help
    # documentation.

    # Define a function that can take a little while.
    @background
    def sleep(time):
        from time import sleep
        sleep(time)
        return time

    print()
    print("Testing out 'background' function..", end=" ")
    import time
    start = time.time()
    out1 = sleep(1)
    out2 = sleep(.5)
    assert(time.time() - start < .01)
    # Initially there should be "0 completed, 2 running, 0 pending".
    assert(sleep.completed == 0)
    assert(sleep.running == 2)
    assert(sleep.pending == 0)
    # Call a method of "out1" knowing the result must be computed.
    out1 = out1.result
    # Now there should be "1 completed, 0 running, 1 pending".
    assert(sleep.completed == 1)
    assert(sleep.running == 0)
    assert(sleep.pending == 1)
    start = time.time()
    # Call a method of "out2" knowing the result must be computed.
    out2 = out2.result
    # No time should elapse because "out2" should already be done.
    assert(time.time() - start < .001)
    # Now there are "2 completed, 0 running, 0 pending".
    assert(sleep.completed == 2)
    assert(sleep.running == 0)
    assert(sleep.pending == 0)
    # Store the actual results instead of the wrapped results.
    assert(type(out1) == int)
    assert(type(out2) == float)
    print("passed.")


def test_capture():
    print()
    print("Testing out 'capture' function..", end=" ")

    @capture
    def test():
        import sys
        print("Hello standard output world!")
        print("Hello standard error world!", file=sys.stderr)
        return "done"
    assert(test() == "done")
    assert(test.stdout[0].strip() == "Hello standard output world!")
    assert(test.stderr[0].strip() == "Hello standard error world!")
    print("passed.")

def test_cache(temp_dir):
    print()
    # Documentation for "simple_add"
    @cache(max_files=10, cache_dir=temp_dir)
    def simple_add(a, b):
        return a.union(b)
    # Run the function many times (with a non-hashable object) to
    # demonstrate that it all works
    for i in range(10):
        print(simple_add({1,2},{i}))

def test_stability_lock(temp_dir):
    @stability_lock(max_tests=5, test_dir=temp_dir)
    def func(a, b):
        # c = 10
        return a + b

    max_val = 20
    for i in range(1,max_val+1):
        a = list(range(1,i))
        b = list(range(i+1,max_val+1))
        func(a,b)

def test_timeout():
    # This is a testing function for the timeout decorator. These
    # comments will still be included in the documentation for the function.
    @timeout(2,"Failed...")
    def sample_func(sec):
        import time
        print(" Going to sleep for '%.1f' seconds."%(sec), end="\r")
        time.sleep(sec)
        print(" Finished sleeping!", end="\r")
    print("@timeout short wait passed:",sample_func(1.5) == None)
    print("@timeout long wait passed: ",sample_func(2.5) == "Failed...")
    # Re-run with the alternative style of decoration.
    def sample_func(sec):
        import time
        print(" Going to sleep for '%.1f' seconds."%(sec), end="\r")
        time.sleep(sec)
        print(" Finished sleeping!", end="\r")
    sample_func = timeout(2,"Failed...")(sample_func)
    print("@timeout short wait passed:",sample_func(1.5) == None)
    print("@timeout long wait passed: ",sample_func(2.5) == "Failed...")
    print()

# This is a testing function for the type_check decorator.
def test_type_check():
    class TestCaseFailed(Exception): pass

    @type_check(int, {int,float}, [float, int], lambda arg: hasattr(arg, "__call__"))
    def func(arg_int, arg_int_or_float, arg_list_float_int, arg_callable):
        print("Passed", type(arg_int), type(arg_int_or_float), 
              type(arg_list_float_int), arg_callable)

    func(1, 1, [.0,1], lambda: None)
    func(0, 1.0, [.0,1], lambda: None)
    class test:
        def __call__(): pass
    func(1, 1, [.0,1], test)

    # Test case 0
    try:
        func()
        raise(TestCaseFailed())        
    except Exception as exc:
        if ("WrongNumberOfArguments" in str(type(exc))):
            print("Passed correct exception for wrong number of arguments.")
        else:
            raise(TestCaseFailed("An unexpected exception was raised for wrong number of arguments.\n\n"+str(type(exc))))
    # Test case 1
    try:
        func(1.0, 0, [.0, 1], 0)
        raise(TestCaseFailed())
    except Exception as exc:
        if ("WrongType" in str(type(exc))):
            print("Passed correct exception for wrong argument type.")
        else:
            raise(TestCaseFailed("An unexpected exception was raised for wrong argument type.\n\n"+str(type(exc))))
    # Test case 2
    try:
        func(0, 0, [.0, 1], 0)
        raise(TestCaseFailed())
    except Exception as exc:
        if ("WrongType_FailedCheck" in str(type(exc))):
            print("Passed correct exception for wrong keyword argument type.")
        else:
            raise(TestCaseFailed("An unexpected exception was raised for wrong keyword argument type.\n\n"+str(type(exc))))
    # Test case 3
    try:
        func(0, 0, 0, 0)
        raise(TestCaseFailed())
    except Exception as exc:
        if ("WrongType" in str(type(exc))):
            print("Passed correct exception for wrong type for list argument.")
        else:
            raise(TestCaseFailed("An unexpected exception was raised for wrong type for list argument.\n\n"+str(type(exc))))
    # Test case 4
    try:
        func(0, 0, [0,0], 0)
        raise(TestCaseFailed())
    except Exception as exc:
        if ("WrongType_NotListed" in str(type(exc))):
            print("Passed correct exception for wrong listed argument type.")
        else:
            raise(TestCaseFailed("An unexpected exception was raised for wrong listed argument type.\n\n"+str(type(exc))))
    # Test case 5
    try:
        func = type_check((int, float), float)(func)
        raise(TestCaseFailed())
    except Exception as exc:
        if ("WrongUsageOfTypeCheck" in str(type(exc))):
            print("Passed correct exception for improper developer usage.")
        else:
            raise(TestCaseFailed("An unexpected exception was raised for improper developer usage.\n\n"+str(type(exc))))



# ==================================================================
#                   Code for testing decorators
# 
if __name__ == "__main__":
    import os
    temp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testing_decorators")

    test_background()
    test_capture()
    test_cache(temp_dir)
    test_stability_lock(temp_dir)
    test_timeout()
    test_type_check()

    # Remove the testing directory.
    import shutil
    shutil.rmtree(temp_dir)

