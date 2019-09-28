<h1 align="center"><code>util.plot</code></h1>

Provides an extensive interface to `HTML` plotting through `plotly`. Simplifies the usage of *offline* python `plotly` plotting. Produce plots without ever interacting directly with the dictionary objects that plotly expects. This module currently supports 2D and 3D scatter plots with numerical axes, histograms, subplots (with varying numbers of plots in each row), animations, box-plots, and plot annotations.

This module was designed so that the `plot.py` file can be downloaded and used standalone, rather than making it entirely it's own Python package. The inline documentation is extensive.

#### [`class Plot`](plot.py#L163)

The object that is used to store relevant data for the creation of a single data visualization.

##### [`def Plot.color`](plot.py#L410)

Interface to the automatic palatte-based color scheme for this plot. This method produces an rgb string compatible with standard plotly `"rgba(%i,%i,%i,%f)"%(<red>,<green>,<blue>,<alpha>)`. Also accepts indices for color numbers (`[0,1,...]`).

##### [`def Plot.add`](plot.py#L717)

This single generic function can be used as a full-fledged interface for generating all plot components. Everything else is built off of it.

##### [`def Plot.add_region`](plot.py#L469)

Given a boolean function (over two dimensions), this uses a grid to identify the boundary of the region of True points. For convex shapes it draws the boundary, for nonconvex shapes it uses square markers to draw.

##### [`def Plot.add_function`](plot.py#L517)

Automatically evaluates a function in a mesh grid over provided intervals. Adds the function to a 2D or 3D plot automatically.

##### [`def Plot.add_histogram`](plot.py#L604)

Given a list of values, constructs an automatically binned histogram of those values in a 2D plot.

##### [`def Plot.add_box`](plot.py#L674)

Give a list of lists of values, generates a set of box plots.

##### [`def Plot.add_annotation`](plot.py#L939)

Add a textual annotation to a 2D or 3D plot.

##### [`def Plot.show`](plot.py#L1038)

Once all components have been added to a plot, use this function to generate the `plotly` internals used to create the final `HTML` render. If `show=True` then this method will also generate the HTML and open the visualization in your default web browser.

##### [`def Plot.add_node`](plot.py#L1288)

Add a single labeled marker, usually useful when creating a `graph` visualization.

##### [`def Plot.add_edge`](plot.py#L1321)

Add a line between labeled node(s), usually useful when creating a `graph` visualization.

##### [`def Plot.graph`](plot.py#L1268)

A wrapped version of the `show` method that sets up the visual to look better for graphs.

#### [`def iplot`](plot.py#1362)

Convenience function for generating (the usual HTML) interactive plots in a Jupyter notebook.

#### [`def create_html`](plot.py#1384)

Given a figure (as defined for `plotly`), generate the `HTML` that can be rendered nicely by a web browser. Removes excessive `plotly` logos and links.

#### [`def multiplot`](plot.py#1478)

Given multiple `Plot` objects, create a single-page visualization that contains all plots.
