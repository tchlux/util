

def test_algorithms():
    from util.algorithms.testing import test_plot
    from util.algorithms import condition, LSHEP, Voronoi, Delaunay, \
        NearestNeighbor, KNN, Shepard, DelaunayP1, DelaunayP2, DelaunayP3

    print("Adding surface to plot..")
    # model = LSHEP()
    # model = condition(Voronoi, method="MPCA")
    # model = condition(Delaunay, method="MPCA", scale=True)
    model = DelaunayP3()
    # model = NearestNeighbor(k=4, method=Voronoi)
    # model = condition(KNN, display=True)(display=True)
    # model = condition(Shepard, display=True)()
    # model = condition(Voronoi, method="MPCA", display=True)()
    # model = LSHEP()
    # f = lambda x: (x[0] - .5)**2 + x[1]
    p,_,_ = test_plot(model, N=100, D=2, low=-.1, upp=1.1, noise=0, # fun=f,
                      random=False, plot_points=4000, classifier=False) # 6, 8
    # model.errors()
    print("Generating plot HTML..")
    p.show()

    # import random
    # model = KNN()
    # n = 100
    # x = np.random.random(size=(n,10))
    # y = [random.choice(['a', 'b', 'c', 'd', 'e']) for i in range(n)]
    # model.fit(x,y)
    # from util.math import is_numeric
    # print(model(np.random.random(10,)))

    


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
    a.stop()
    print("passed.")


def test_AtomicOpen(display=False):
    print("Testing AtomicOpen..", end=" ")
    import time, os
    from util.parallel import map
    from util.system import AtomicOpen, Timer

    # Testing without atomic writes.
    def write_test(string="testing..", atomic=False):
        open_file = AtomicOpen if atomic else open
        with open_file(test_file, "r") as f:
            if display: print(string, "file opened")
            if display: print(string, f.read().strip())
            # print(string, file=f)
            time.sleep(2)
        if display: print(string, "file closed")
    # Testing for atomic writes.
    def atomic_write_test(string): return write_test(string, atomic=True)

    p = 3
    import util.parallel
    util.parallel.MAX_PROCS = 50

    t = Timer()
    if display: print("Testing parallel regular write..")

    if display: print()
    test_file = "_test_.txt"
    with open("_test_.txt", "w") as f:
        print("'conents'", file=f)

    list(map(write_test, map(str,range(p)), redirect=False))
    if display:
        with open(test_file) as f: print(f'{"-"*70}\n{f.read()}{"-"*70}')
        print(t())
    # assert( t() < 4 )

    t.start()
    if display: print("Testing parallel atomic write..")

    if display: print()
    test_file = "_test_.txt"
    with open("_test_.txt", "w") as f:
        print("'conents'", file=f)

    list(map(atomic_write_test, map(str,range(p)), redirect=False))
    if display:
        with open(test_file) as f: print(f'{"-"*70}\n{f.read()}{"-"*70}')
        print(t())
    # assert( t() > 6 )

    os.remove(test_file)
    print("passed.")



if __name__ == "__main__":
    # test_algorithms()
    test_Timer()
    # test_AtomicOpen()
