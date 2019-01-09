
def test_Timer():
    print("Testing Timer..", end=" ")
    import time
    from util.system import Timer
    a = Timer()
    first = a.start
    assert(a.check() < .01)
    time.sleep(1)
    assert(first - a.start == 0)
    assert(abs(a.check() - 1) < .01)
    a.start()
    assert(abs(first - a.start) > 1)
    assert(a.check() < .001)    
    print("passed.")


def test_AtomicOpen(display=False):
    print("Testing AtomicOpen..", end=" ")
    import time, os
    from util.parallel import map
    from util.system import AtomicOpen, Timer

    test_file = "_test_.txt"
    # Testing without atomic writes.
    def write_test(string="testing..", atomic=False):
        open_file = AtomicOpen if atomic else open
        with open_file(test_file, "w") as f:
            if display: print(string, "file opened")
            print(string, file=f)
            time.sleep(2)
        if display: print(string, "file closed")
    # Testing for atomic writes.
    def atomic_write_test(string): return write_test(string, atomic=True)

    p = 3
    import util.parallel
    util.parallel.MAX_PROCS = 50

    t = Timer()
    if display: print("Testing parallel regular write..")
    list(map(write_test, map(str,range(p)), redirect=False))
    if display:
        with open(test_file) as f: print(f'{"-"*70}\n{f.read()}{"-"*70}')
        print(t())
    assert( t() < 4 )

    t.start()
    if display: print("Testing parallel atomic write..")
    list(map(atomic_write_test, map(str,range(p)), redirect=False))
    if display:
        with open(test_file) as f: print(f'{"-"*70}\n{f.read()}{"-"*70}')
        print(t())
    assert( t() > 6 )

    os.remove(test_file)
    print("passed.")



if __name__ == "__main__":
    test_Timer()
    test_AtomicOpen()
