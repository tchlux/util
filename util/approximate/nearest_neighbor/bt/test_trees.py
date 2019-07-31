

if __name__ == "__main__":
    import numpy as np

    from util.system import Timer
    t = Timer()
    dim = 200
    # train = 100000
    train = 100000
    test = 1
    k = 10
    ls = (2**17) // dim # Leaf size
    ls = 40
    print("leaf size:", ls)
    print("Allocating array..", flush=True)
    np.random.seed(2)
    x = np.random.random(size=(train,dim))
    z = np.random.random(size=(test,dim))
    print()
    print("x:", x.shape)
    print("z:", z.shape)
    print()
    print()

    print("Fortran Ball Tree")
    from balltree import BallTree as BT
    t.start()
    tree = BT(x, leaf_size=ls)
    ct = t.stop()
    print("Construction time:", ct)
    t.start()
    d, i = tree.query(z, k=k, leaf_size=ls)
    qt = t.stop()
    print("Query time:       ", qt)
    print("d: ",d)
    print("i: ",i[0])
    print()
    i1 = i[0].copy()

    print("Brute Force")
    t.start()
    d = np.sqrt(np.sum(x**2, axis=1) + np.sum(z[0]**2) - 2 * np.dot(x, z[0]))
    # d = np.linalg.norm(x - z[0], axis=1)
    i = np.argsort(d)
    qt = t.stop()
    i = i[:k]
    d = d[i]
    print("Query time:", qt)
    print("d: ",d)
    print("i: ",i)
    print()

    # print("VP Tree")
    # from vp_tree import VPTree
    # t.start()
    # tree = VPTree(x, leaf_size=ls)
    # ct = t.stop()
    # print("Construction time:", ct)
    # t.start()
    # d, i = tree.query(z, k=k)
    # qt = t.stop()
    # print("Query time:       ", qt)
    # print("d: ",d)
    # print("i: ",i[0])
    # print()

    from sklearn.neighbors import KDTree, BallTree

    print("KD Tree")
    t.start()
    tree = KDTree(x, leaf_size=ls)
    ct = t.stop()
    print("Construction time:", ct)
    t.start()
    d, i = tree.query(z, k=k)
    qt = t.stop()
    print("Query time:       ", qt)
    print("d: ",d)
    print("i: ",i[0])
    print()
    i2 = i[0].copy()

    print("Sklearn Ball Tree")
    t.start()
    tree = BallTree(x, leaf_size=ls)
    ct = t.stop()
    print("Construction time:", ct)
    t.start()
    d, i = tree.query(z, k=k)
    qt = t.stop()
    print("Query time:       ", qt)
    print("d: ",d)
    print("i: ",i[0])
    print()
    i3 = i[0].copy()


    assert(np.all(i1 == i2))
    assert(np.all(i1 == i3))
