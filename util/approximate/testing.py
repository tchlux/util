import numpy as np

# Given a model, test the time it takes to make a prediction
def test_time(model, ns=range(10,10011,1000), ds=range(2,11,2)):
    import time
    print("Starting")
    for d in ds:
        print(f"D: {d}")
        print()
        for n in ns:
            print(f" N: {n}")
            x = np.random.random((n,d))
            y = np.random.random((n,))
            start = time.time()
            model.fit(x, y)
            print(f"  fit:     {time.time() - start:5.2f}")
            start = time.time()
            model(x[0])
            print(f"  predict: {time.time() - start:5.2f}")
            print()
        print()

# Given a model, generate a test plot demonstrating the surface it procudes
def test_plot(model, low=0, upp=1, plot_points=3000, p=None,
              fun=lambda x: 3*x[0]+.5*np.cos(8*x[0])+np.sin(5*x[-1]),
              N=20, D=2, noise=0., random=True, seed=0, x=None, 
              y=None, classifier=False):
    np.random.seed(seed)
    provided_points = (type(x) != type(None)) and (type(y) != type(None))
    if (type(x) == type(None)):
        # Generate x points
        if random:
            x = np.random.random(size=(N,D))
        else:
            N = int(round(N ** (1/D)))
            x = np.array([r.flatten() for r in np.meshgrid(*[np.linspace(0,1,N)]*D)]).T
    if (type(y) == type(None)):
        # Calculate response values
        y = np.array([round(fun(v) + np.random.random()*noise)
                      if classifier else
                      fun(v)+np.random.random()*noise for v in x])
    # Fit the model to the points
    model.fit(x,y, classifier=classifier)
    # Generate the plot
    from util.plot import Plot
    if type(p) == type(None): p = Plot()
    if not provided_points: p.add("Training Points", *x.T, y)
    p.add_func(str(model), model, *([(low-.1,upp+.1)]*D),
               plot_points=plot_points, vectorized=True)
    # p.add_func("truth", fun, *([(low-.1,upp+.1)]*D),
    #            plot_points=plot_points)
    return p, x, y

# Given a model, use the method "points_and_weights" to identify the
# regions of support for each of the interpolation points.
def test_support(model, low=0, upp=1, plot_points=3000, p=None,
                 fun=lambda x: 3*x[0]+.5*np.cos(8*x[0])+np.sin(5*x[-1]),
                 N=20, D=2, random=True, seed=0):
    # Force D to be 2
    D = 2
    np.random.seed(seed)
    # Generate x points
    if random:
        x = np.random.random(size=(N,D))
    else:
        N = int(round(N ** (1/D)))
        x = np.array([r.flatten() for r in np.meshgrid(*[np.linspace(0,1,N)]*D)]).T
    # Calculate response values
    y = np.array([fun(v) for v in x])
    # Fit the model to the points
    model.fit(x,y)
    # Generate the plot
    from util.plotly import Plot
    if type(p) == type(None): p = Plot()
    p.add("Training Points", *x.T, color=p.color(len(x)))
    for i in range(len(x)):
        name = f"{i+1}" 
        p.add(name, [x[i][0]], [x[i][1]], group=name)
        def supported(pt):
            pts, wts = model.points_and_weights(pt)
            return (i in pts)
        p.add_region(name+" region", supported, *([(low-.1,upp+.1)]*D),
                     color=p.color(p.color_num), group=name,
                     plot_points=plot_points, show_in_legend=False) 
    # p.add_func(str(model), model, *([(low-.1,upp+.1)]*D),
    #            plot_points=plot_points, vectorized=True)
    return p, x, y

