if __name__ == "__main__":
    import numpy as np
    from util.plot import Plot
    from util.algorithms import Delaunay

    # Settings for creating the example.
    N           = 20
    dim         = 2
    random      = True
    plot_points = 5000

    # Function for testing.
    mult = 5
    low  = 0
    upp  = 1
    fun = lambda x: np.cos(x[0]*mult) + np.sum(np.sin(x[1:]*mult))

    # Generate the X and Y for the points.
    if random:
        x = np.random.random(size=(N,dim))
    else:
        N = int(round(N ** (1/dim)))
        x = np.array([r.flatten() for r in np.meshgrid(
            * ((np.linspace(low,upp,N),)*dim) )]).T
    y = np.array([fun(v) for v in x])

    # Shift to add padding
    padding = .4
    low -= (upp - low)*padding
    upp += (upp - low)*padding

    # Fit the Delaunay model to the points
    surf = Delaunay()
    surf.fit(x,y)

    # Create the surface in the plot (calling Delaunay).
    p = Plot()
    p.add("Training Points", *x.T, y)
    p.add_func(str(surf), surf, *([(low,upp)]*dim), plot_points=plot_points)
    p.plot(file_name="test_plot.html")


