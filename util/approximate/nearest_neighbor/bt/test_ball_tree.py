import numpy as np

LARGE_TEST = True
TEST_SORT = False
TREE_TEST = False
COMPARE_AGAINST_SKLEARN = True


if TEST_SORT:
    if LARGE_TEST: N = 1000000
    else:          N = 10
    from . import btmod
    # Import a timer.
    from util.system import Timer
    t = Timer()
    # Generate test numbers.
    print()
    print(f"Generating {N} numbers..", flush=True)
    x = np.random.random(size=N)
    i = np.arange(len(x)) + 1
    print()
    # Test the fortran code.
    pts = x.copy()
    ids = i.copy()
    t.start()
    pts, ids = btmod.argsort(pts, ids)
    t.stop()
    # Check for correctness.
    ids_match = np.all(x[ids-1] == pts)
    is_sorted = np.all(np.diff(pts)>=0)
    try: assert(ids_match and is_sorted)
    except:
        print("ERROR")
        print(" ids_match: ",ids_match)
        print(" is_sorted: ",is_sorted)
        print()
        print("pts:", pts)
        print("ids:", ids)
        print()
        print(x)
        print(x[ids-1])
        exit()
    print("argsort:", t())
    # Test the NumPy code.
    pts = x.copy()
    ids = i.copy()
    t.start()
    ids = pts.argsort()
    t.stop()
    print("numpy:  ", t())


if TREE_TEST:
    from balltree import BallTree
    # from sklearn.neighbors import BallTree

    np.random.seed(0)
    if LARGE_TEST: size = (100000,1000)
    else:          size = (7,2)
    print()
    print(f"Allocating array.. {size}", flush=True)
    x = np.random.random(size=size)
    # x = np.random.random(size=size)


    if len(x) < 20:
        print()
        print(x.T.shape)
        print()
    else:
        print("Building tree..", flush=True)
        from util.system import Timer
        t = Timer()
        t.start()

    leaf_size = 1
    tree = BallTree(x, leaf_size=leaf_size)

    if len(x) < 20:
        print('-'*70)
        print()
        print(tree.tree.T)
        print()
        print(tree.order)
        print(tree.radii)
        print(tree.sq_sums)
        print()
        print('-'*70)
    else:
        t.stop()
        print(f"done in {t()} seconds.", flush=True)
        print()

    z = np.random.random(size=(1,size[1]))
    k = 7
    d,i = tree.query(z, k=k)

    true_dists = np.linalg.norm(x - z[0,:], axis=1)
    print("np.argmin(true_dists): ",np.argmin(true_dists))
    print("np.min(true_dists):    ",np.min(true_dists))
    print()
    print("i: ",i)
    print("d: ",d)

    if len(x) < 20:
        d = np.linalg.norm(x - z[0], axis=1)
        i = np.argsort(d)
        print()
        print("i: ",i)
        print("d: ",d[i])


if COMPARE_AGAINST_SKLEARN:
    print()
    print("="*70)

    from util.system import Timer
    t = Timer()

    if LARGE_TEST: train, dim = 100000, 1000
    else:          train, dim = 7, 2
    test = 1
    k = 5
    print("Initializing data..", flush=True)
    np.random.seed(0)
    x = np.random.random(size=(train,dim))
    z = np.random.random(size=(test,dim))
    print()
    print("x:", x.shape)
    print("z:", z.shape)

    from balltree import BallTree as BT
    print()
    print("Fortran Ball Tree")
    t.start()
    tree = BT(x, leaf_size=10)
    ct = t.stop()
    print("Construction time:", ct)
    t.start()
    d, i = tree.query(z, k=k)
    qt = t.stop()
    print("Query time:       ", qt)
    print("d: ",d[0])
    print("i: ",i[0])
    d1, i1 = d[0].copy(), i[0].copy()

    from sklearn.neighbors import BallTree
    print()
    print("Sklearn Ball Tree")
    t.start()
    tree = BallTree(x, leaf_size=10)
    ct = t.stop()
    print("Construction time:", ct)
    t.start()
    d, i = tree.query(z, k=k)
    qt = t.stop()
    print("Query time:       ", qt)
    print("d: ",d[0])
    print("i: ",i[0])
    d2, i2 = d[0].copy(), i[0].copy()

    print()
    print("Brute Force")
    t.start()
    d = np.sqrt(np.sum(x**2, axis=1) + np.sum(z[0]**2) - 2 * np.dot(x, z[0]))
    i = np.argsort(d)
    qt = t.stop()
    i = i[:k]
    d = d[i]
    print("Query time:", qt)
    print("d: ",d)
    print("i: ",i)
    d3, i3 = d.copy(), i.copy()

    max_diff = max(max(abs(d1 - d2)), max(abs(d1-d3)))
    ds_match = max_diff < 2**(-26)
    is_match = np.all(i1 == i2) and np.all(i1 == i3)
    print()
    print(f"Max difference in distance calculations:\n   {max_diff:.32}")
    try: assert(ds_match and is_match)
    except:
        print()
        print("ERROR")
        print("  ds_match: ",ds_match)
        print("  is_match: ",di_match)

