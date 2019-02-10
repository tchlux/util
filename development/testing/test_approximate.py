from util.approximate import *

def test_approximate():
    from util.approximate.testing import test_plot
    from util.approximate import condition, LSHEP, Voronoi, Delaunay, \
        NearestNeighbor, KNN, Shepard
    from util.approximate.delaunay import DelaunayP1, DelaunayP2

    print("Adding surface to plot..")
    # model = LSHEP()
    # model = condition(Voronoi, method="MPCA")
    # model = condition(Delaunay, method="MPCA", scale=True)
    model = DelaunayP2()
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

    
def test():
    print("Testing 'util.approximate'..", end=" ")
    test_approximate()
    print("passed.")
