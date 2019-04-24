from util.parallel import *


def test_map():
    import time, random
    # A slow function
    def slow_func(value, seconds=.45):
        import time, random
        random.seed(value)
        time.sleep(random.random()*.1 + seconds)
        return str(value)

    print()
    print("Max processes:", MAX_PROCS)
    print()
    # Generate a normal "map" and a parallel "map"
    values = list(range(1,10+1))
    start = time.time()
    s_gen = builtin_map(slow_func, values)
    print("Standard map construction time:", time.time() - start)
    start = time.time()
    po_gen = map(slow_func, values, order=True)
    print("Parllel ordered map construction time: ", time.time() - start)
    start = time.time()
    pu_gen = map(slow_func, values, order=False)
    print("Parllel unordered map construction time: ", time.time() - start)

    start = time.time()
    print()
    for value in s_gen:
        print(value, end=" ", flush=True)
    print()
    print("Standard map time:", time.time() - start)

    start = time.time()
    print()
    for value in po_gen:
        print(value, end=" ", flush=True)
    print()
    print("Ordered Parallel map time:", time.time() - start)

    start = time.time()
    print()
    for value in pu_gen:
        print(value, end=" ", flush=True)
    print()
    print("Unordered Parallel map time:", time.time() - start)

    killall()
    clear_logs()


# Development testing code.
if __name__ == "__main__":
    test_map()
