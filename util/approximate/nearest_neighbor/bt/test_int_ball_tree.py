import numpy as np

LARGE_TEST = True
TEST_SORT = False
TREE_TEST = False
COMPARE_AGAINST_SKLEARN = True


if TEST_SORT:
    if LARGE_TEST: N = 1000000
    else:          N = 10
    from balltree import int_btmod
    # Import a timer.
    from util.system import Timer
    t = Timer()
    # Generate test numbers.
    print()
    print(f"Generating {N} numbers..", flush=True)
    x = np.asarray(np.random.random(size=N), dtype='float32')
    i = np.asarray(np.arange(len(x)) + 1, dtype='int32')
    print()
    # Test the fortran code.
    pts = x.copy()
    ids = i.copy()
    t.start()
    pts, ids = int_btmod.argsort(pts, ids)
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
    if LARGE_TEST: size = (1000,100000)
    else:          size = (7,2)
    print()
    print(f"Allocating array.. {size}", flush=True)
    x = np.asarray((np.random.random(size=size)-1/2)*128, order='f', dtype='int8')

    # Add 128 to convert uint8 into correctly ordered int8

    if len(x) < 20:
        print()
        print(x.shape)
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

    k = 7
    z = np.asarray((np.random.random(size=(1,size[1]))-1/2)*128,
                    order='f', dtype='int8')
    print()
    print("z: ",z)
    print()
    d,i = tree.query(z, k=k)

    true_dists = np.linalg.norm(x[tree.index_mapping] - z[0,:], axis=1)
    print("np.argmin(true_dists): ",np.argmin(true_dists))
    print("np.min(true_dists):    ",np.min(true_dists))
    print()
    print("  i: ",i)
    print("  d: ",d)

    if len(x) < 20:
        d = np.linalg.norm(x[tree.index_mapping].astype('float32') -
                           z[0].astype('float32'), axis=1)
        i = np.argsort(d)
        print()
        print("Truth:")
        print("  i: ",i)
        print("  d: ",d[i])

if COMPARE_AGAINST_SKLEARN:
    print()
    print("="*70)

    from util.system import Timer
    t = Timer()

    if LARGE_TEST: train, dim = 1000, 100000
    else:          train, dim = 7, 2
    test = 1
    k = 5
    print("Initializing data..", flush=True)
    np.random.seed(0)
    x = np.asarray((np.random.random(size=(train,dim))-1/2)*128, order='F', dtype='int8')
    z = np.asarray((np.random.random(size=(test, dim))-1/2)*128, order='F', dtype='int8')
    print()
    print("x:", x.shape)
    print("z:", z.shape)
    if (not LARGE_TEST):
        print()
        print("x: ",x)
        print("z: ",z)        
        print()

    from balltree import BallTree as BT
    print()
    print("Fortran Ball Tree", flush=True)
    # Transpose x to make watching memory consumption easier.
    t.start()
    tree = BT(x, leaf_size=10)
    ct = t.stop()
    print("Construction time:", ct, tree.dtype)
    t.start()
    d, i = tree.query(z, k=k)
    qt = t.stop()
    print("Query time:       ", qt)
    print("d: ",d[0])
    print("i: ",i[0])
    d1, i1 = d[0].copy(), i[0].copy()

    # Re-order "x" to be the same for all trees.
    _ = x[tree.index_mapping]
    del(x)
    x = _

    from sklearn.neighbors import BallTree
    print()
    print("Sklearn Ball Tree", flush=True)
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
    # Convert 'x' and 'z' to types that will be exact for brute force approach.
    x = x.astype('float64')
    z = z.astype('float64')
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
    print(f"Max difference in distance calculations:\n   {max_diff:.3e}")
    try: assert(ds_match and is_match)
    except:
        print()
        print("ERROR")
        print( "  is_match:",is_match)
        print(f"  ds_match: {ds_match} {max(abs(d1-d3)):.3e} {max(abs(d1 - d2)):.3e}")

