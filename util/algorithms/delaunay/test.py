if __name__ == "__main__":
    import numpy as np
    from util.plot import Plot
    from util.algorithms import Delaunay

    # Settings for creating the example.
    N           = 30
    dim         = 2
    random      = False
    plot_points = 1000

    # Function for testing.
    mult = 5
    low  = 0
    upp  = 1
    fun = lambda x: np.cos(x[0]*mult) + np.sum(np.sin(x[1:]*mult))
    np.random.seed(0)
    # Generate the X and Y for the points.
    if random:
        x = np.random.random(size=(N,dim))
    else:
        width = int(round(N ** (1/dim)))
        x = np.array([r.flatten() for r in np.meshgrid(
            * ((np.linspace(low,upp,width),)*dim) )]).T
        x = x[:N]

    y = np.array([fun(v) for v in x])

    # Shift to add padding
    padding = .4
    low -= (upp - low)*padding
    upp += (upp - low)*padding


    # PMODE=1, once with PMODE=2, once with PMODE=3, and once with PMODE unspecified and CHUNKSIZE=10.
    # Fit the Delaunay model to the points
    surf = Delaunay()
    surf.fit(x,y)

    # Manual test at a specific point.
    # test = np.array([[.3587, .422]])
    test = np.array([[.25, .25]])

    # # Parallel surface (wrong)
    # guess_surf = Delaunay(parallel=True)
    # guess_surf.fit(x,y)
    # guess_pts, guess_wts = guess_surf._predict(test.copy())
    # # Serial surface (right)
    # true_surf = Delaunay(parallel=False)
    # true_surf.fit(x,y)
    # true_pts, true_wts = true_surf._predict(test.copy())
    # # Display results
    # print()
    # print("True points: ", true_pts)
    # print("    weights: ", true_wts)
    # print("      value: ", np.sum(y[true_pts] * true_wts))
    # print()
    # print("Guess points:", guess_pts)
    # print("     weights:", guess_wts)
    # print("      value: ", np.sum(y[guess_pts] * guess_wts))
    # print()
    # exit()

    # Create the surface in the plot (calling Delaunay).
    p = Plot()
    p.add("Training Points", *x.T, y)
    p.add_func(str(surf), surf, *([(low,upp)]*dim), plot_points=plot_points)
    p.plot(file_name="test_plot.html")


