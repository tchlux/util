# This module serves to provide a simplified interface to *offline*
# python plotly plotting. The user can produce plots without ever
# interacting directly with the dictionary objects that plotly
# expects. This module currently supports 2D and 3D scatter plots with
# numerical axes, histograms, subplots (with varying numbers of plots
# in each row), and plot annotations.

import random, numbers, os, webbrowser, imp, sys
import numpy as np
import colorlover as cl
from scipy.spatial import Delaunay, ConvexHull
from scipy.spatial.qhull import QhullError

# Importer function to make sure that the a local file name doesn't
# conflict with the true installed package name
def import_package(name, custom_name=None):
    # Use a custom name if provided, otherwise just name
    custom_name = custom_name or name
    # Find and open a file that is the non-local package name
    f, pathname, desc = imp.find_module(name, [p for p in sys.path if
                                               p not in os.getcwd()])
    # Load the module
    module = imp.load_module(custom_name, f, pathname, desc)
    return module

# Load the pypi package "plotly" that interfaces with plotly.js
plotly = import_package("plotly")

PLOT_POINTS = 1000
BRIGHTNESS_RANGE = 0.6

PALATTE = np.array(cl.to_numeric(cl.scales['5']['qual']['Set2']))
PALATTE = PALATTE**2
PALATTE = PALATTE / np.max(PALATTE) * 255
# Re-order the palatte so that the colors appear better
PALATTE = np.concatenate((PALATTE[1:], [PALATTE[0]]))
# Save the color palatte for plotting a gradient
DEFAULT_GRADIENT = np.array(cl.to_numeric(cl.scales['11']['div']['Spectral']))[::-1]
MIN_COLORS = 40

# Expand the palatte using random combinations of existing colors
RANDOM_SEED = 0
random.seed(RANDOM_SEED)
palatte_size = len(PALATTE)
for i in range(MIN_COLORS - palatte_size):
    # Create lots of extra colors
    c = np.array([random.choice(PALATTE[:palatte_size,0]),
                  random.choice(PALATTE[:palatte_size,1]),
                  random.choice(PALATTE[:palatte_size,2])])
    # Add this new random color to the palatte
    PALATTE = np.concatenate( (PALATTE, [c]), axis=0 )
# Re-seed the random number generator so that it is not tainted
random.seed()

DEFAULT_CAMERA_POSITION = dict(
    up=dict(x=0, y=0, z=1),
    center=dict(x=0, y=0, z=0),
    eye=dict(x=-1.0, y=-2.0, z=0.7)
)

#      Coloring Data     
# =======================

# Given some data, color the data according to a palatte with uniform
# interpolation between the colors in the palatte from the minimum
# value provided to the maximimum value provided
def color_data(values, palatte=DEFAULT_GRADIENT, opacity=1.0):
    no_none = [v for v in values if type(v) != type(None)]
    shift = min(no_none)
    scale = (max(no_none) - shift) * 1.11
    if (scale == 0): scale = 1.0
    def color(value):
        if value == None: return None
        # Generate the index as a float (for interpolating)
        index = len(palatte) * (value-shift) / scale
        # Get the exact colors on either side of this index
        lower = int(index)
        upper = lower + 1
        if (lower > len(palatte)-1): lower = len(palatte)-1
        if (upper > len(palatte)-1): upper = len(palatte)-1
        index -= lower
        # Interpolate between the lower and upper colors
        c = tuple(palatte[lower]*(1-index) + palatte[upper]*(index))
        # Return the interpolated color.
        return 'rgba(%i,%i,%i,%f)'%(c+(opacity,))
    return list(map(color, values))



#      Class for building a plotly plot     
# ==========================================

# Class that serves as an interface to the standard "data & layout"
# containers that need to be managed in order to produce Plotly plots.
# This class uses the offline modes of plotly to produce local HTML
# files rather than the standard web-based ones. This class also
# attempts to strip web-related features (such as the Plotly logo)
# from the upper-right hand corner plot interface.
# 
# All functionality is encapsulated in the "Plot.add" command, which
# allows for all standard plotly options to be controlled in the
# construction of data, along with the "Plot.plot" command, which
# allows for all standard plotly options that controll layout.
# 
# Additional methods that are effectively decorated versions of the
# "add" command include: 
#    add_histogram  --  For quickly creating vertically oriented or
#                       horizontally oriented histograms.
#    add_function   --  For passing a function and automatically
#                       sampling it across a meshgrid and plotting.
#    add_region     --  For drawing convex regions in 2D by providing 
#                       a boolean function that is True inside the
#                       region and False outside of the region.
# 
# The "add_annotation" function can be used to add text descriptions
# with arrows over points of interest in a plot.
# 
# The "plot" function is also capable of appending to existing HTML
# files by setting the keyword argument "append=True". This is nice
# for producing a single scrollable multi-page HTML report of plots.
# 
# The "multiplot" function, provided in this module and not part of
# the class, allows for the produciton of single pages that fit
# multiple plots. See documentation of "multiplot" for more detials.
class Plot:
    def __init__(self, title="", x_title="x", y_title="y",
                 z_title="z", mode="markers", palatte=PALATTE):
        self.title = title
        self.x_title = x_title
        self.y_title = y_title
        self.z_title = z_title
        self.x_min_max = [float('inf'), -float('inf')]
        self.y_min_max = [float('inf'), -float('inf')]
        self.z_min_max = [float('inf'), -float('inf')]
        # Specific booleans for tracking internal state
        self.is_3d = False
        self.to_reverse = []
        # Data for tracking default plot settings
        self.color_num = -1
        self.data = []
        self.annotations = []
        self.mode = mode
        self.palatte = palatte
        self.palatte_size = len(palatte)

    # Return an appropriate face color of a simplex given the simplex,
    # data z values, and either (color index and opaicty, or a list of
    # colors associated with each data value.
    def _simp_color(self, simp, z, color_ind=None, opacity=1.0, colors=None):
        shift = max(z)
        scale = shift - min(z)
        has_none = type(None) in (type(v) for v in z[simp])
        if (scale > 0) and (not has_none): 
            # If colors were provided, then average them to produce out color
            if type(colors) != type(None):
                # Get the color of each node in the simplex as a numpy array
                colors = [colors[i] for i in simp]
                colors = [c[c.index('(')+1:c.index(')')].split(',')
                          for c in colors]
                colors = np.array([list(map(float,c)) for c in colors])
                if colors.shape[1] != 4:
                    colors = np.concatenate((
                        colors,np.ones(shape=(colors.shape[0],1))),
                                            axis=1)
                # return the average color of points in the simplex
                return 'rgba(%f,%f,%f,%f)'%tuple(np.sum(colors,axis=0) / len(simp))
            else:
                simp_avg = sum(z[simp]) / len(simp)
                brightness = (1.0-BRIGHTNESS_RANGE/2) + ((simp_avg - shift) / scale) * BRIGHTNESS_RANGE
        else:
            brightness = 1.0
        return self.color(color_ind, brightness, opacity)


    # Prepare all annotations for the type of plot being  presented.
    def _clean_annotations(self, annotations):
        if not self.is_3d:
            for a in annotations:
                a.pop('z', '')
        else:
            for a in annotations:
                if type(a['z']) == type(None):
                    a['z'] = 0
        return annotations

    # Prepares all the data sets to be plotted in whatever dimension
    # is highest (2 or 3). Creates 3D meshes for all surfaces. Should
    # be stable if called multiple times, but this code is still in
    # development stage.
    def _clean_data(self, data):
        # Remove all references to 3D layout if this is a 2D plot
        if not self.is_3d:
            # 2D PLOT SETUP
            for d in data:
                d.pop('z','')
                d.pop('hoverinfo','')
                d.pop('text','')
                # Special case for plotting histograms
                if d['type'] == 'histogram':
                    if type(d.get('y','')) == type(None):
                        d.pop('y','')
                    if type(d.get('x','')) == type(None):
                        d.pop('x','')
                    d.pop('line','')
                    d.pop('mode','')
                    d.pop('fill','')
                    d['opacity'] = d['marker'].pop('opacity','')
                    d['marker'].pop('symbol','')
                    d['marker'].pop('size','')
                    d['marker']['color'] = d.pop('fillcolor','')
        else:
            # 3D PLOT SETUP
            for ind,d in enumerate(data):
                # Add z values to all scatters that may have been added
                if d['type'] == 'scatter':
                    d['z'] = np.zeros(len(d['x']))
                    d['type'] = 'scatter3d'                            
                    if d['marker']['size'] == None:
                        d['marker']['size'] = 5
                # Convert fill and / or lines into surfaces
                conv_2d = (not self.is_3d) and ('lines' in d['mode'])
                if (d.get('fill','') == 'toself') or conv_2d:
                    print("WARNING: Converting 2D to 3D automatically.")
                    d['type'] = 'surface'
                    # Get the opacity of the surface
                    if d.get('fill','') != None:
                        d['opacity'] = float(d['fillcolor'].split(',')[-1].strip(')'))
                    else:
                        d['opacity'] = float(d['line']['color'].split(',')[-1].strip(')'))
                # If user wants a surface, construct one! 
                # (plotly default surfaces are not very cooperative)
                if ('surface' in d['type']):
                    points_2D = np.vstack([d['x'], d['y']]).T
                    try:
                        simps = Delaunay(points_2D).simplices
                        d['type'] = 'mesh3d'
                        d['i'] = simps[:,0]
                        d['j'] = simps[:,1]
                        d['k'] = simps[:,2]
                        # Generate face colors with average simplex z-coordinate
                        d['facecolor'] = list(map(
                            lambda simp: self._simp_color(
                                simp, d['z'], ind,
                                d['marker']['opacity'],
                                d['marker']['color']), simps
                        ))
                        if 'opacity' not in d:
                            d['opacity'] = d['marker']['opacity']
                        d.pop('marker','')
                        d.pop('mode','')
                        d.pop('text','')
                    except QhullError:
                        d['type'] = 'scatter3d'
                        if 'mode' not in d:
                            d['mode'] = 'lines'
                # Pop out the unnecessary attributes for 3D plots
                d.pop('fill','')
                d.pop('fillcolor','')
                if 'line' not in d.get('mode',''):
                    d.pop('line','')

    # Manage plotly reverse order bug (only happens with "tonext[xy]")
    def _reorder_data(self, data):
        start = end = None
        # Cycle through the elements of data
        for i,tr in enumerate(self.to_reverse):
            if (tr and start==None):
                start = i
            if (not tr and start!=None):
                end = i+1
                # Reverse that group of plot series
                data[start:end] = data[start:end][::-1]
                start = end = None
        # Reverse the final group when self.to_reverse[-1] == True
        if (start!=None):
            end = len(data)
            data[start:end] = data[start:end][::-1]

        self.to_reverse = [False] * len(data)

    # ===================================
    #      User accessible functions     
    # ===================================

    # Interface to the automatic palatte-based color scheme for this
    # plot. This method produces an rgb string compatible with
    # standard plotly "rgba(%i,%i,%i,%f)"%(<red>,<green>,<blue>,<alpha>). 
    # 
    # This method takes the following arguments:
    #   Arg Name (Default) -- Description
    # 
    #   number (None)      -- Index in the palatte of the desired color.
    #   brightness (1.0)   -- Value ranging from 0.0 to 1.0 that can
    #                         be used to produces shades of the same color.
    #   alpha (1.0)        -- Opacity of the color produced. Note, 0.0
    #                         for this argument will cause the color
    #                         to be invisible.
    #   color (None)       -- String ".*([0-9]+,[0-9]+,[0-9][,0-9\.]*).*", 
    #                         list, tuple, or numpy array meant to be
    #                         converted into the standard rgba string.
    def color(self, number=None, brightness=1.0, alpha=None, color=None):
        if color == None:
            if (number == None): number = self.color_num
            if (number < len(self.palatte)):
                # If we have fewer entries than the palette size
                c = self.palatte[number]
            else:
                # Otherwise we have to create a new palette entry
                c = np.array([random.choice(self.palatte[:self.palatte_size,0]),
                              random.choice(self.palatte[:self.palatte_size,1]),
                              random.choice(self.palatte[:self.palatte_size,2])])
                # Add this new random color to the palatte
                self.palatte = np.concatenate( (self.palatte, [c]), axis=0 )
        elif type(color) == str:
            # Get the color as a list of numbers
            c = color[color.index('(')+1:color.index(')')].split(',')
            # Make sure the color only has [red, green, blue, alpha]
            c = np.array(list(map(float,c)))
            if (len(c) > 3) and (type(alpha) == type(None)):
                alpha = c[-1]
            c = c[:3]
        elif (type(color) == tuple or 
              type(color) == list or 
              type(color) == np.ndarray):
            c = np.array(color[:3])
            if (len(color) > 3) and (type(alpha) == type(None)): 
                alpha=color[-1]
        else:
            raise(Exception("ERROR: Color must either be a string, tuple, list, or numpy array."))
        # Define a default alpha if necessary
        if type(alpha) == type(None):
            alpha = 1.0
        # Apply the brightness to the color
        c = c*brightness
        c = np.where(c > 255, 255, c)
        c = np.where(c < 0, 0, c)
        # Return the color as a plotly color string
        return 'rgba(%i,%i,%i,%f)'%(tuple(c)+(alpha,))

    # Decorated "add" function that automatically attempts to
    # find the edge of a convex 2D region given a function that is
    # True inside and False outside. Uses a meshgrid of "plot_points"
    # points in order to approximate the boundary of the region.
    # 
    #  name        -- The string name of the series being added
    #  func        -- A function that, given a single (x,y) point 
    #                 returns True or False.
    #  min_max_x   -- A length-2 iterable for the x-range over which
    #                 to apply the meshgrid.
    #  min_max_y   -- A length-2 iterable for the y-range over which
    #                 to apply the meshgrid.
    #  plot_points -- The number of plot points in the
    #                 meshgrid. Higher numbers will yield more precise
    #                 boundaries for the region.
    #  ... <standard "add" arguments with adjusted defaults> ...
    def add_region(self, name, func, min_max_x=None, min_max_y=None,
                   plot_points=PLOT_POINTS, mode="lines", opacity=0.1,
                   fill="toself", line_width=0, **kwargs):
        if self.is_3d: raise(Exception("ERROR: Regions only work for 2D plots."))
        if type(min_max_x) == type(None):
            min_max_x = self.x_min_max.copy()
        if type(min_max_y) == type(None):
            min_max_y = self.y_min_max.copy()
        if max(map(abs,min_max_x+min_max_y)) == float('inf'):
            raise(Exception("ERROR: Invalid x or y range."))
        # Round up the number of plot points per axis
        plot_points = int(plot_points**(0.5) + 0.5)
        # Calculate the mesh grid of x and y values
        x_vals = (np.linspace(*min_max_x, num=plot_points),)
        y_vals = (np.linspace(*min_max_y, num=plot_points),)
        x,y = np.meshgrid(x_vals, y_vals)
        test_pts = np.vstack((x.flatten(), y.flatten())).T
        region_pts = np.array([pt for pt in test_pts if func(pt)])
        # Try reducing to the set of convex hull points for the region
        # and plotting that, if it fails simply print an error message.
        try:
            hull_pts = region_pts[ConvexHull(region_pts).vertices]
            self.add(name, hull_pts[:,0], hull_pts[:,1], mode=mode,
                     opacity=opacity, fill=fill,
                     line_width=line_width, **kwargs)
        except:
            print("ERROR: Could not add region.")

    # Decorated "add" function that automatically generates the
    # response values for a given "func" over a meshgrid using
    # "plot_points" points (works for 2D or 3D plotting depending on
    # how many "min_max..." ranges are provided).
    # 
    #  name        -- The string name of the series being added
    #  func        -- A function that, given a single (x[,y]) point 
    #                 returns a numeric type object.
    #  min_max_x   -- A length-2 iterable for the x-range over which
    #                 to apply the meshgrid.
    #  min_max_y   -- A length-2 iterable for the y-range over which
    #                 to apply the meshgrid. (only provided for 3D)
    #  grid_lines  -- Whether or not to add lines whose intersections
    #                 show where plot points were placed (only works
    #                 for 3D plotting).
    #  plot_points -- The number of plot points in the meshgrid.
    #  ... <standard "add" arguments with adjusted defaults> ...
    def add_func(self, name, func, min_max_x, min_max_y=[],
                 grid_lines=True, plot_points=PLOT_POINTS, mode=None, 
                 plot_type=None, **kwargs):
        if len(min_max_y) > 0: self.is_3d = True
        # If we have two control axes, square root the plot points
        if self.is_3d:
            plot_points = int(plot_points**(0.5) + 0.5)
            # If no y was provided, set it to default value
            if len(min_max_y) == 0: min_max_y = [0.0,0.0]
            if mode == None: plot_type = 'surface'
        else:
            if mode == None: mode = 'lines'

        # Generate the input points
        x_vals = (np.linspace(*min_max_x, num=plot_points),)
        if self.is_3d:
            x_vals += (np.linspace(*min_max_y, num=plot_points),)
        x_vals = tuple(x.flatten() for x in np.meshgrid(*x_vals))

        # Get the response values
        try:
            response = np.array([func(x) for x in np.vstack(x_vals).T]).flatten()
        except SystemExit: exit()
        except:
            # Provide a useful error message if the response values
            # could not be computed correctly
            try:
                sample = func((np.vstack(x_vals).T)[0])
                raise(Exception("Error in return value of provided function. Expected number, got '%s'"%(type(sample))))
            except:
                raise(Exception("Error computing the provided function."))

        # For analyzing the specific outputs
        # temp = np.concatenate((x_vals,[response])).T
        # np.savetxt("/Users/thomaslux/Desktop/temp.txt", temp)

        # Call the standard plot function
        self.add(name, *x_vals, response, mode=mode,
                 hoverinfo="name+x+y+z", plot_type=plot_type, **kwargs)
            
        # If this is a 3D surface plot and grid_lines=True, add grid lines
        if (self.is_3d and plot_type == 'surface') and grid_lines:
            opacity = kwargs.get("opacity",1.0)
            line_color = kwargs.get("line_color",'rgb(0,0,0)')
            for row in range(plot_points):
                x = x_vals[0][row*plot_points:(row+1)*plot_points]
                y = x_vals[1][row*plot_points:(row+1)*plot_points]
                z = response[row*plot_points:(row+1)*plot_points]
                self.add("", x,y,z, show_in_legend=False,
                         group=name+" (lines)", mode="lines",
                         line_width=opacity, opacity=opacity, 
                         color=line_color, hoverinfo="none")

                indices = np.arange(plot_points)*plot_points + row
                x = x_vals[0][indices]
                y = x_vals[1][indices]
                z = response[indices]
                self.add("", x,y,z, show_in_legend=False,
                         group=name+" (lines)", mode="lines",
                         line_width=opacity, opacity=opacity,
                         color=line_color, hoverinfo="none")

    # Decorated "add" function that automatically sets the options
    # necessary for plotting an N-bin PDF histogram of a given set of
    # values. By default the bars are separated along "bar_spacing"
    # axis, and the area of all bars together adds to 1.
    # 
    #  name        -- The string name of the series being added
    #  values      -- A list of ints or floats.
    #  bar_spacing -- "x" if the x-axis should be bins and y-axis
    #                 probabilities, "y" for transposing the setup.
    #  num_bins    -- The number of evenly spaced bins to use when
    #                 generating the histogram.
    #  padding     -- The amount of spacing on the min and max sides
    #                 of the histogram that is produced.
    #  histnorm    -- Standard plotly "histnorm" argument, can be
    #                 "probability" or "count" most commonly.
    #  barmode     -- Standard plotly "barmode" argument. When set to
    #                 "", plotly default will be used where
    #                 multi-series histograms will be non-overlapping.
    #                 When set "overlay", histogram series can overlap.
    #  opacity     -- See "add" function.
    def add_histogram(self, name, values, bar_spacing="x", num_bins=100,
                      padding=0.03, opacity=0.7, histnorm='probability',
                      barmode='overlay', **kwargs):
        # Check for errors in usage.
        if bar_spacing not in ("x", "y"):
            raise(Exception("ERROR: Invalid 'bar_spacing', only 'x' or 'y' are acceptable."))
        if num_bins <= 0:
            raise(Exception("ERROR: Invalid 'num_bins', must be a positive integer."))
        if len(values) == 0:
            raise(Exception("ERROR: Empty list passed in for 'values'."))
        values_name = bar_spacing + "_values"
        autobin = "autobin" + bar_spacing
        bins = bar_spacing + "bins"
        self.histogram_barmode = barmode
        # Calculate the range of the histogram
        hist_value_range = max(values) - min(values)
        hist_start_val = min(values) - hist_value_range * padding
        hist_end_val   = max(values) + hist_value_range * padding
        # Provide necessary keyword arguments (that the user has not already)
        if (values_name not in kwargs):
            kwargs[values_name] = values
        kwargs['histnorm'] = histnorm
        if (autobin not in kwargs):
            kwargs[autobin] = False
            if (bins not in kwargs):
                bin_settings = dict( start=hist_start_val,
                                     end=hist_end_val,
                                     size=hist_value_range/num_bins )
                kwargs[bins] = bin_settings
        # Store the correct extrema to be used for plotting
        min_max = getattr(self, bar_spacing+"_min_max").copy()
        min_max[0] = min(hist_start_val, min_max[0])
        min_max[1] = max(hist_end_val,   min_max[1])
        # Call the 'add' function with updated arguments
        self.add(name, plot_type='histogram', opacity=opacity, **kwargs)
        # Make sure min_max were not wrongly changed, use the extrema
        # of the desired bins as the range, not the extrema of values
        getattr(self, bar_spacing+"_min_max")[0] = min_max[0]
        getattr(self, bar_spacing+"_min_max")[1] = min_max[1]

    # Primary function for simplifying the interface to plotly
    # plotting. This single generic function can be used as a
    # full-fledged interface for generating standard plotly "data"
    # dictionary object. It can be used for both 2D and 3D plotting,
    # and allows for control of all aspects of plot styling.
    # 
    # STANDARD ARGUMENTS: The combination of these that is provided
    # determines whether a 2D or 3D plot is produced. "x_values" are
    # optional because histograms may only have y-values given. For
    # most standard usage, (x,y) will be given for 2D, (x,y,z) for 3D.
    # 
    #  name     -- Name of the series to be plotted
    #  x_values -- The x-values associated with the series
    #  y_values -- The y-values associated with the series
    #  z_values -- The z-values associated with the series
    # 
    # HIGH-LEVEL STYLING: 
    #  mode           -- The plotly series mode, "lines",  "markers",
    #                    "text", or combinations with a "+" between.
    #  plot_type      -- The plotly plot_type, "scatter[3d]" for plots
    #                    of lines and dots, "surface" for 3D surfaces,
    #                    "histogram" for producing histograms.
    #  group          -- The legend-series group name. This is used
    #                    for the simultaneous hide/show of multiple
    #                    series. This will cause increased legend spacing.
    #  show_in_legend -- True or False for if this series should show
    #                    in the legend. Currently plotly legends do
    #                    *not* support 3D surfaces in legends.
    #  shade          -- True or False if the given data series should
    #                    be shaded with different brightnesses based
    #                    on magnitude.
    #  use_gradient   -- True or False if a gradient coloring should
    #                    be applied to the given data series.
    #  palatte        -- The palatte to use when creating a gradient
    #                    of colors for the "use_gradient" option.
    #  text           -- A list of the text strings that should be
    #                    shown for each data point when a user hovers
    #                    with their mouse over that data point.
    # 
    # LOW-LEVEL STYLING:
    #  color        -- The series color as a tuple/list/array of 
    #                  (<red>,<green>,<blue>[,<alpha>])
    #                  rgb in [0,255], alpha in [0,1]
    #  opacity      -- Transparency constant for series color, 0 is
    #                  completely transparent, 1 is completely opaque.
    #                  this value is overwritten if "color" has 4 numbers.
    #  line_color   -- The color of the line for this series
    #  line_width   -- The width of the line for this series
    #  fill         -- Almost exactly the plotly "fill" argument,
    #                  options include "toprevy" "tozeroy" "toself"
    #                  and the same for x. If "tonext[xy]" is used,
    #                  the legened will be reversed. (plotly bug)
    #  fill_color   -- The color to use for the fill if active.
    #  fill_opacity -- The opacity of the fill color.
    #  symbol       -- The marker symbol, standard plotly. "circle",
    #                  "square", and a lot more on their website.
    #  dash         -- Standard plotly "dash" option. "solid", "dot",
    #                  "dash", or "1px,2px,5px[,[0-9]*px]*" list of lengths
    #  marker_size       -- The size (in pixels) of markers
    #  marker_colors     -- The color of markers
    #  marker_line_width -- The width of the bounding line of markers
    #  marker_line_color -- The color of the bounding line of markers
    #  hoverinfo         -- The information displayed when the user's
    #                       mouse hovers over the plot. Options include
    #                       "x" "y" "z" "text" "name", combined with "+"
    # 
    #  ... <any additional plotly data-dictionary args> ...
    def add(self, name, x_values=None, y_values=None, z_values=None,
            mode=None, plot_type=None, group=None,
            show_in_legend=True, shade=True, use_gradient=False,
            palatte=DEFAULT_GRADIENT, text=None, color=None,
            opacity=1.0, line_color=None, line_width=None, fill=None,
            fill_color=None, fill_opacity=0.6, symbol='circle',
            dash=None, marker_size=None, marker_colors=None,
            marker_line_width=0, marker_line_color='rgba(50,50,50,0.8)', 
            hoverinfo='name+text', **kwargs):

        # Convert the x, y (and z) values into numpy arrays and
        # store 'values' for creating marker colors based on magnitude
        if type(x_values) != type(None):
            x_values = np.array(x_values)
            values = x_values
            no_none = [v for v in x_values if type(v) != type(None)]
            if len(no_none) != 0:
                self.x_min_max = [min(min(no_none), self.x_min_max[0]),
                                  max(max(no_none), self.x_min_max[1])]
        if type(y_values) != type(None):
            y_values = np.array(y_values)
            values = y_values
            no_none = [v for v in y_values if type(v) != type(None)]
            if len(no_none) != 0:
                self.y_min_max = [min(min(no_none), self.y_min_max[0]),
                                  max(max(no_none), self.y_min_max[1])]
        if type(z_values) != type(None):
            self.is_3d = True
            z_values = np.array(z_values)
            values = z_values
            no_none = [v for v in z_values if type(v) != type(None)]
            if len(no_none) != 0:
                self.z_min_max = [min(min(no_none), self.z_min_max[0]),
                                  max(max(no_none), self.z_min_max[1])]

        # Make a nice pretty gradient of color
        if use_gradient:
            marker_colors = color_data(values, palatte)

        # Define z-values if none were given and we need them, and plot type
        if self.is_3d:
            if plot_type == None:
                plot_type = 'scatter3d'
            if type(z_values) == type(None):
                z_values = np.zeros(len(x_values))
            # Define text for all the data points
            if (hoverinfo != None) and ("text" in hoverinfo) and (text == None):
                # hoverinfo = None
                # text = None
                # WARNING: Sometimes this is causing problems where
                # the hoverinfo labels do not update on scroll, it
                # looks like another bug in the python plotly.
                text = ["%s: %s<br>%s: %s<br>%s: %s"%(
                    self.x_title,x, self.y_title,y, self.z_title,z)
                        for (x,y,z) in zip(x_values,y_values,z_values)]
        else:
            if plot_type == None:
                plot_type = 'scatter'

        # Process mode
        if type(mode) == type(None):
            mode = self.mode
        # Set the color if none was provided
        if color == None:
            self.color_num += 1
            color = self.color(self.color_num, alpha=opacity)
        else:
            shade = False
        if line_color == None:
            line_color = color
        if fill_color == None:
            fill_color = self.color(color=color, alpha=fill_opacity)
        else: 
            fill_color = self.color(color=fill_color, alpha=fill_opacity)
        if not marker_colors:
            if shade:
                marker_colors = []
                no_none = [v for v in values if v != None]
                if len(no_none) > 0:
                    shift = min(no_none)
                    scale = max(no_none) - shift
                    if scale == 0: scale = 1.0
                    for v in values:
                        if not isinstance(v,numbers.Number):
                            raise(Exception((
                                "ERROR: '%s' not permitted. Only "+
                                "numbers are allowed as values.")%(v)))
                        brightness = ((1.0-BRIGHTNESS_RANGE/2) +
                                      ((v - shift) / scale) *
                                      BRIGHTNESS_RANGE)
                        marker_colors.append( self.color(color=color,
                                                         brightness=brightness, 
                                                         alpha=opacity) )
            else:
                marker_colors = [color]*len(values)

        # Special plotly failure mode, need to reverse data for
        # 'tonext' to actually mean 'next' instead of 'previous'. This
        # bug has been reported, but no one in the plotly community is
        # addressing it (or even noticing it) as a problem.
        self.to_reverse.append((type(fill) == str) and ("tonext" in fill))

        # Now add the standard plotly "data" object to local storage
        self.data.append(dict(
            type = plot_type,
            name = name,
            x = x_values,
            y = y_values,
            z = z_values,
            hoverinfo = hoverinfo,
            text = text,
            marker = dict(
                # Generate colors based on point magnitude
                color = color if ("lines" in mode) else marker_colors,
                # color = marker_colors,
                size = marker_size,
                opacity = opacity,
                symbol = symbol,
                line = dict(
                    width = marker_line_width,
                    color = marker_line_color
                )),
            line = dict(
                width = line_width,
                color = line_color,
                dash = dash
            ),
            mode = mode,
            fill = fill,
            fillcolor = fill_color,
            legendgroup = group,    
            showlegend = show_in_legend
        ))

        # Update the newly created dictionary with any custom user settings
        self.data[-1].update(kwargs)


    # Add an annotation to the plot. These will be text boxes
    # stationed in the absolute foreground of the plot, disregarding
    # occlusion in 3D plots.
    # 
    # STANDARD ARGUMENTS
    #  text -- The text to display in the annotation.
    #  x    -- The x coordinate of the arrow for the annotation
    #  y    -- The y coordinate of the arrow for the annotation
    #  z    -- The z coordinate (if applicable) of the arrow for the annotation
    # 
    # ANNOTATION CONTROL
    #  ax        -- The x screen pixels offset for the anntation box (+ is right)
    #  ay        -- The y screen pixels offset for the annotaiton box (+ is down)
    #  opacity   -- The transparency of the entire annotation
    #  textangle -- The angle of the annotation (and bounding box)
    #  align     -- The alignment of text within the annotation box
    #  xanchor   -- The box-x anchor point for the extending arrow
    #  yanchor   -- The box-y anchor point for the extending arrow
    # 
    # FONT CONTROL
    #  font_family -- The family of font used in the annotation
    #  font_color  -- The color of the font used in the annotation
    #  font_size   -- The size of the font used in the annotation
    # 
    # BORDER CONTROL
    #  border_color -- The color of the border of the annotation box
    #  border_width -- The thickness of the border of the annotation box
    #  border_pad   -- The padding between the annotation text and box
    #  bg_color     -- The background color of the annotation box
    # 
    # ARROW CONTROL
    #  show_arrow  -- Whether or not to show an arrow at all
    #  arrow_color -- The color of the arrow
    #  arrow_size  -- The size of the arrow head
    #  arrow_width -- The width of the arrow line
    #  arrow_head  -- The type of arrow head. 0 -> None, 1-5 -> Arrows,
    #                 6 -> Dot, 7 -> Box, >7 -> None
    # 
    #  ... <any additional plotly annotation-dictionary args> ...
    def add_annotation(self, text, x, y, z=None, ax=None, ay=None,
                       opacity=0.8, text_angle=0, align="left",
                       x_anchor="center", y_anchor="bottom",
                       font_family="Arial", font_color="#0a0a0a",
                       font_size=12, border_color="#1a1a1a",
                       border_width=0, border_pad=4,
                       bg_color="#f0f0f0", show_arrow=True,
                       arrow_color="#666", arrow_size=1,
                       arrow_width=1, arrow_head=7, **kwargs):
        # Add computed values for the annotation x and y
        if show_arrow:
            if ax == None: ax = 10
            if ay == None: ay = -20
        else:
            if ax == None: ax = 0
            if ay == None: ay = 0
        # Add the annotation
        self.annotations.append(dict(
            text=text,
            # Target location
            x = x,
            y = y,
            z = z,
            ax = ax,
            ay = ay,
            # Annotation text control
            opacity = opacity,
            textangle = text_angle,
            align = align,
            # Anchor and shift
            xanchor = x_anchor,
            yanchor = y_anchor,
            xshift = 0,
            yshift = 0,
            # Font
            font = dict(
                family = font_family,
                color = font_color,
                size = font_size
            ),
            # Border control
            bordercolor = border_color,
            borderwidth = border_width,
            borderpad = border_pad,
            bgcolor = bg_color,
            # Arrow control
            showarrow = show_arrow,
            arrowcolor = arrow_color,
            arrowsize = arrow_size,
            arrowwidth = arrow_width,
            arrowhead = arrow_head, 
        ))
        self.annotations[-1].update(kwargs)

    # Second part to the simplified plotly interface. This creates the
    # layout-dictionary object and (optionally) produces the HTML and
    # opens a browser to view the plot.
    # 
    # COMMON ARGUMENTS:
    #  title       -- Title to display for this plot. (can include
    #                 HTML line break <br> and bold <b>text</b>)
    #  x_range     -- The range of x-values to default to displaying,
    #                 automatically determined by data if possible
    #  y_range     -- The range of y-values to default to displaying,
    #                 automatically determined by data if possible
    #  z_range     -- The range of z-values to default to displaying,
    #                 automatically determined by data if possible
    #  xaxis       -- Dictionary of xaxis settings, useful for type = 
    #                 "log", "date", "category" axes settings.
    #  yaxis       -- Same as xaxis, but for y.
    #  zaxis       -- Same as xaxis, but for z.
    #  fixed       -- False if plotly should automatically rescale the
    #                 plot when series are hidden/shown, True if
    #                 plotly should not rescale on hide/show.
    #  show_legend -- True if the legend should be included.
    # 
    # LAYOUT CONTROL:
    #  layout          -- Update to be performed to the plotly
    #                     layout-dictionary that is generated.
    #  aspect_mode     -- For 3D plotting, standard plotly.
    #  scene_settings  -- Standard plotly, for updating the "scene"
    #                     dictionary for 3D plotting.
    #  axis_settings   -- Controls for all of the axes. Include things
    #                     like showgrid, zeroline, showline,
    #                     showticklabels (all boolean) or ticks="<str>"
    #  camera_position -- A dictionary of dictionaries of x,y,z
    #                     values, "up" is relative up vector, "center"
    #                     is the point about which a 3D plot rotates, 
    #                     and "eye" is the camera coordinate.
    # 
    # OUTPUT CONTROL:
    #  html      -- True if "create_html" should be called.
    #  file_name -- See "create_html".
    #  show      -- See "create_html".
    #  append    -- See "create_html".
    # 
    #  ... <any additional plotly.offline.plot keyword arguments> ...
    def plot(self, title=None, x_range=None, y_range=None,
             z_range=None, xaxis={}, yaxis={}, zaxis={},
             fixed=True, show_legend=True, layout={}, 
             aspect_mode='cube', scene_settings={}, axis_settings={},
             camera_position=DEFAULT_CAMERA_POSITION,
             html=True, file_name="temp-plot.html", show=True,
             append=False, **kwargs):
        # Update title, and all plot axis ranges
        if title == None:
            title = self.title
        if (fixed and x_range == None and
            max(map(abs,self.x_min_max)) != float('inf')):
            x_width = self.x_min_max[1] - self.x_min_max[0]
            x_range = [self.x_min_max[0] - 0.05*x_width,
                       self.x_min_max[1] + 0.05*x_width]
        if (fixed and y_range == None and
            max(map(abs,self.y_min_max)) != float('inf')):
            y_width = self.y_min_max[1] - self.y_min_max[0]
            y_range = [self.y_min_max[0] - 0.05*y_width,
                       self.y_min_max[1] + 0.05*y_width]
        if (fixed and z_range == None and
            max(map(abs,self.z_min_max)) != float('inf')):
            z_width = self.z_min_max[1] - self.z_min_max[0]
            z_range = [self.z_min_max[0] - 0.05*z_width,
                       self.z_min_max[1] + 0.05*z_width]
        # Generate the layout (titles and legend)
        plot_layout = dict(
            title = title,
            showlegend = show_legend,
       )
        # Set the barmode for histograms if necessary
        if (hasattr(self, 'histogram_barmode') and
            len(self.histogram_barmode) > 0):
            plot_layout['barmode'] = self.histogram_barmode
        # Clean all annotations so they are ready for plotting
        annotations = [a.copy() for a in self.annotations]
        self._clean_annotations(annotations)
        # Setup for the axes of the plot
        scene = dict(
            xaxis = dict(title = self.x_title, range=x_range, **axis_settings),
            yaxis = dict(title = self.y_title, range=y_range, **axis_settings),
            zaxis = dict(title = self.z_title, range=z_range, **axis_settings),
        )
        # Setup the plot layout (different for 2D and 3D plots)
        if not self.is_3d:
            plot_layout.update(scene)
            plot_layout.pop('zaxis')
            plot_layout.update(dict(annotations=annotations))
        else:
            scene['aspectmode'] = aspect_mode
            scene['camera'] = camera_position
            scene.update(scene_settings)
            scene.update(dict(annotations=annotations))
            plot_layout['scene'] = scene

        # Update the plot layout with any specific user settings given
        plot_layout.update(layout)
        # Make sure all the data entries are prepared to be plotted
        # Make a deep copy of the locally stored data that can be
        # cleaned and prepared for plotting (without risk of deleting
        # information that may be necessary for re-plotting)
        data = [d.copy() for d in self.data]
        self._clean_data(data)
        # Manage plotly reverse order bug (only happens with "tonext_")
        self._reorder_data(data)
        # Generate the figure
        fig = plotly.graph_objs.Figure(data=data, layout=plot_layout)
        # Create the html file and show in browser if appropriate
        if html: create_html(fig, file_name, show, append, **kwargs)
        # Return the figure
        return fig

        
#      Functions for manipulation produces plots     
# ===================================================

# Generates the HTML file and fixes some javascript so that the
# plot does not have unwanted buttons and links.
# 
#  fig       -- A plotly figure-dictionary, with all necessary keys.
#  file_name -- The name of the output file to generate.
#  show      -- True if the output file should be opened in a
#               webbrowser after writing is complete.
#  append    -- True if this HTML code should be appended to
#               "file_name" if it already exists. This creates a
#               scrollable page, where each plot takes a full screen.
# 
#  ... <any additional plotly.offline.plot keyword arguments> ...
def create_html(fig, file_name="temp-plot.html", show=True,
                append=False, **kwargs):
    if (not append):
        print("Creating new", end=" ")
    else:
        print("Appending to", end=" ")
    print("file named: '%s'"%file_name)
    # Store the old file contents if we are appending
    if (append and os.path.exists(file_name)):
        with open(file_name) as f:
            old_contents = f.read()
    # Generate the plot offline 
    plotly.offline.plot(fig, auto_open=False, filename=file_name,
                        show_link=False, **kwargs)
    # Remove unnecessary modebar buttons and the plotly logo link
    with open(file_name) as f:
        file_string = f.read()
    file_string = file_string.replace(
        'displaylogo:!0', 'displaylogo:!1')
    file_string = file_string.replace(
        'modeBarButtonsToRemove:[]',
        'modeBarButtonsToRemove:["sendDataToCloud", "select2d", "lasso2d"]')
    file_string += "\n\n"
    # If appending, put the old contents back in front of the new
    if append: file_string = old_contents + file_string
    with open(file_name, "w") as f:
        f.write(file_string)
    # Open the plot in a webbrowser if the user wants that
    if show: webbrowser.open("file://"+os.path.abspath(file_name))
    return file_name


# Make multiple plots fit onto one browser window, options for sharing
# axes as well for naming purposes. Mixed plot types allowed too!
# Supports different number of columns per row, but no other mismatch
# 
#  plots     -- A 2D list of plots in the desired grid layout. Rows
#               can have varying numbers of plots, columns cannot.
#  x_domains -- A 2D list of pairs (3D list) each pair is [start,end]
#               where 0 <= start < end <= 1. This controls the width
#               of each column of plots. Same 2D shape as "plots".
#  y_domains -- A 2D list of pairs (3D list) each pair is [start,end]
#               where 0 <= start < end <= 1. This controls the width
#               of each row of plots. Same 2D shape as "plots".
#  shared_y  -- True if the y-axis is shared for plots in same row.
#  shared_x  -- True if the x-axis is shared for plots in same column.
#  gap       -- The amount of space between the plots.
#  specs     -- A 2D list (same shape as "plots") of dictionaries
#               representing plotly subplots "specs". Mostly for
#               telling plotly which plots are 3D and which are 2D.
#  html      -- True if "create_html" should be called.
#  show      -- See "create_html".
#  append    -- See "create_html".
# 
#  ... <any additional plotly.offline.plot keyword arguments> ...
def multiplot(plots, x_domains=None, y_domains=None, html=True,
              show=True, append=False, specs=None, shared_y=False,
              shared_x=False, gap=0.12, **kwargs): 
    # Make sure the plots array is 2D
    try:    plots[0][0]
    except: plots = [plots]
    # Convert given plots into figures (if figures were not given
    for r in plots:
        for c in range(len(r)):
            if type(r[c]) == Plot:
                r[c] = r[c].plot(html=False, show=False)

    # Count the number of rows and columns
    rows = len(plots)
    cols = [len(r) for r in plots]
    max_cols = max(c for c in cols)
    # Generate/Process the specs
    if type(specs) != type(None):
        try:    specs[0][0]
        except: specs = [specs]
    else:
        specs = [[None]*max_cols for r in range(rows)]
        for r,row in enumerate(plots):
            for c,plot in enumerate(row):
                if type(plot) == type(None): continue
                sample_data = plots[r][c]['data'][0]
                specs[r][c] = {"is_3d": 'z' in sample_data}
    # Generate the x and y domains if they are not provided by the user
    if x_domains == None:                
        x_domains = []
        for r in range(rows):
            width = (1 - (cols[r]-1)*gap) / cols[r]
            x_domains.append(
                [[c*(width+gap), c*(width+gap) + width]
                 for c in range(cols[r])])
    if y_domains == None:                
        height = (1 - (rows-1)*gap) / rows
        y_domains = [[r*(height+gap), r*(height+gap) + height]
                     for r in range(rows)]
    # Identify the number of dimensions provided in x an y domains, if
    # too few, then make sure it is the same shape as the plots
    try: x_domains[0][0][0]
    except TypeError:
        x_domains = [x_domains for r in range(rows)]
    try: y_domains[0][0][0]
    except TypeError:
        y_domains = [[y_domains[r]]*cols[r] for r in range(rows)]
            
    # Fix y-domains so that they are specified from bottom to top
    flipped_y = []
    gap = y_domains[1][0][0] - y_domains[0][0][1] if len(y_domains) > 1 else 0
    for r in range(rows):
        start = 0.0 if r == 0 else flipped_y[-1][1] + gap
        width = y_domains[rows-r-1][0][1] - y_domains[rows-r-1][0][0]
        flipped_y.append([start, start+width])
    y_domains = [[flipped_y[r]]*cols[len(cols)-1-r] for r in range(rows)][::-1]

    # Generate the holder for the multiplot
    fig = plotly.tools.make_subplots(rows=rows, cols=max_cols,
                                     specs=specs,
                                     shared_yaxes=shared_y,
                                     shared_xaxes=shared_x)
    # Generate the multi plot!
    counter_2d = 0
    counter_3d = 0
    for r,row in enumerate(plots):
        for c,plot in enumerate(row):
            # Allows for empty spaces
            if type(plot) == type(None): continue
            # Otherwise, continue assuming we have a figure!
            for d in plot['data']:
                fig.append_trace(d, r+1, c+1)
            # Extract the annotations for this plot
            plot_annotations = plot['layout'].pop('annotations',[])
            # Handle 3D and 3D differently
            if specs[r][c]['is_3d']:
                counter_3d += 1
                scene_name = 'scene' + str(counter_3d)
                fig['layout'][scene_name].update(plot['layout']['scene'])
                fig['layout'][scene_name]['domain']['x'] = x_domains[r][c]
                fig['layout'][scene_name]['domain']['y'] = y_domains[r][c]
            else:
                counter_2d += 1
                x_name = 'xaxis'+str(counter_2d)
                y_name = 'yaxis'+str(counter_2d)
                # For shared axes, only add the first entry of column or row
                # Update the domains as specified by the user
                if (not shared_x) or (r == 0):
                    fig['layout'][x_name].update(plot['layout'].pop('xaxis'))
                    fig['layout'][x_name]['domain'] = x_domains[r][c]
                if (not shared_y) or (c == 0):
                    fig['layout'][y_name].update(plot['layout'].pop('yaxis'))
                    fig['layout'][y_name]['domain'] = y_domains[r][c]
                for a in plot_annotations:
                    a['xref'] = "x" + str(counter_2d)
                    a['yref'] = "y" + str(counter_2d)
                    fig['layout']['annotations'] = fig['layout'].get(
                        'annotations',[]) + [a]
                # Ensure that no axis layouts make it into the plot that shouldn't
                plot['layout'].pop('xaxis','')
                plot['layout'].pop('yaxis','')
            fig['layout'].update(plot['layout'])
            # Return the annotations to the plot now that the figure
            # has been updated (and is not at risk of overwriting annotations)
            if len(plot_annotations) > 0:
                plot['layout']['annotations'] = plot_annotations
            # Remove the 'scene' if there is one left over
            if specs[r][c]['is_3d']: fig['layout'].pop('scene','')


    # Create the html plot if the user wants that (pass extra arguments)
    if html: create_html(fig, show=show, append=append, **kwargs)
    # Return the figure to be plotted
    return fig


# ================================================
#      Example Usage of This Plotly Interface     
# ================================================

if __name__ == "__main__":
    print()
    print("Creating a demonstration of most of the available (and useful) features!")
    print()

    # Testing code for the plotting interface
    fun = lambda x: np.sum(x**2) / 10
    # fun = lambda x: x[-1]*x[-2]
    x = np.linspace(-10,10,100)
    y = x**2 / 10    

    # # Decide randomly whether to make a 2D, 3D, or both
    # selection = np.random.randint(3)
    selection = 2
    
    # Plot the 2D
    if selection in [0,2]:
        plot = Plot("2D Title")
        # Adding a 2D function
        plot.add_func("Test Func 2D", fun,[-10,10], opacity=0.5, dash="dot")
        # Adding lines with dots
        plot.add("V Line", [0,0], [min(y), max(y)], mode="lines+markers")
        # Adding a filled region
        plot.add("Square", [-2,-2,2,2], [5,10,10,5], opacity=0.8,
                 mode="none", fill="toself")
        # Adding lines in arbitrary directions
        plot.add("H Line", [-5,5], [1,1], mode="lines+markers",
                 symbol='square', dash="1px,3px,1px")
        plot.add("H Line 2", [-5,5], [2,2], mode="lines")
        plot.add_annotation("2D Annotation", 10+.1, 10-.1, ax=40, ay=20,
                            arrow_head=2, y_anchor="top")
        plot1 = plot

    # Plot the 3D
    if selection in [1,2]:
        plot = Plot("3D Title","X Axis", "Y Axis", "Z Axis")
        rand_x = list(range(-5,6,2))
        rand_y = np.random.randint(-3,4,size=6)
        rand_z = np.random.randint(3,8,size=6)
        # Adding a 3D line
        plot.add("3D Line", rand_x, rand_y, rand_z, mode='lines')
        dens = 5
        x, y = np.meshgrid(np.linspace(-5,5,dens), np.linspace(-5,5,dens))
        x = x.flatten()
        y = y.flatten()
        fun = lambda x: -.3*x[1] + 1/2*x[0] + 1
        z = np.array(list(map(fun, zip(x,y))))
        # Adding a 3D function, and demonstrating different marker styles
        plot.add("3D Above", x, y, z+1.5, marker_size=3,
                 marker_line_width=1, group="Small")
        plot.add("3D Below", x, y, z-1.5, marker_size=2,
                 marker_line_width=1, group="Small")
        plot.add("3D B Next", x, y, z-1, marker_size=5, opacity=0.7,
                 marker_line_width=1, group="Big" )
        plot.add("3D A Next", x, y, z+1, marker_size=7, opacity=0.4,
                 marker_line_width=1, group="Big")
        plot.add_func("3D Surface", fun, [min(x),max(x)],
                      [min(y),max(y)], opacity=0.7, use_gradient=True)
        x_val, y_val = x[-5], y[-5]
        plot.add_annotation("3D Annotation", x_val, y_val, 
                            fun([x_val,y_val])+1.5, ax=-15)
        plot2 = plot

    # Adding a histogram, notice they don't have the same ranges and
    # that will reflect in their respective bin sizes.
    plot3 = Plot("Last Title", "x stuff", "y stuff")
    plot3.add_histogram("Histogram Series 1", np.random.normal(0,3,size=(400,)))
    plot3.add_histogram("Histogram Series 2", np.random.normal(15,1, size=(200,)))
    plot3.add_annotation("Histogram annotation", 0, 0.005)

    # Showing multiple plots on one screen, a grid layout with the
    # option for varying numbers of elements on each row.
    multiplot([[plot1, plot2],[plot3]], gap=0.1)

    # Demonstrate how to put a full-screen plot beneath the first.
    plot2.plot(title="'append=True' Plotting",append=True)
    # Demonstrate allowing plotly to auto-scale when series are
    # activated and deactivated (try turning off Histogram Series 1)
    plot3.plot(title="'fixed=False' Plotting", fixed=False, append=True)