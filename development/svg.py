'''
Stick Plot is an SVG plotting and animation package written in
python. It does *not* use any external libraries in an attempt to
remove its users from (sometimes vicious) update cycles of popular
web-based plotting libraries like D3. It provides slightly more
features than the nearest neighbor library "pygal", most notably
animations and human-readable SVG files are provided by stickplot.

The goal of this library is to produce pretty human-readable (and
editable) SVG files that allow for web-based display of a diverse
class of data. This library allows simplified access to most all SVG
plotting features.

Enjoy and please file issues with minimum working examples for bug
reporting and feature requests.
'''

# - Construct general coordinates for "in plot" objects.
# - Construct general coordinates for "in legend" objects.
# - If "no legend" then plot right = canvas right
# - If "no title" then plot top = canvas top
# - If "no xaxis" then plot bottom = canvas bottom + numbers
# - If "no yaxis" then plot left = canvas left + numbers
# - Place x-axis and y-axis titles
# - Place x-axis and y-axis number labels
# - Make one number for "plot padding"

# 3D plots have 3D coordinates
# 3D plots have animations

# Make axes grow out of their origin.
# Make line plots build in from line origin to end.
# Make data "shiver" a little when hovering over it
# Make data "build in" from nothing (flat vertical / horizontal, circle)
# Make a logo stick figure that waves
# Make a logo stick figure that runs

# https://plot.ly/python/line-and-scatter/
# https://developer.mozilla.org/en-US/docs/Web/SVG/Tutorial/Paths
# https://developer.mozilla.org/en-US/docs/Web/SVG/Tutorial/Basic_Shapes
# https://stackoverflow.com/questions/25150865/beginning-and-pausing-svg-animations-on-hover
# http://stylecampaign.com/blog/2014/02/svg-animation/
# https://stackoverflow.com/questions/20998840/svg-rotation-in-3d
# http://stylecampaign.com/blog/2014/03/3d-to-svg/



# For 3D, we need to generate the entire set of rotations and
# calculate the z values for all objects at each rotation.

# Each "type" plot should have it's own class that constructs the
# individual SVG objects (and animation key frames) from data.

# Rendering will be done as soon as the plot is actually requested.

# Each series type should be able to "render_data" and "render_legend".

# "render_data" will take a (time_start) and (time_end) that
# determines when in the [0,1] animation frames should be drawn, and a
# "project" function that takes data points as input and returns the
# "x,y" coordinates in terms of the drawable area.

# The "Layout" will decide where the legend is and where the series
# should be rendered relative to the scene.

# The "Layout" will sequentially render each object at different
# viewpoints constructing the animation frames.

# Animation labels provided as strings will be uniformly spaced 
# (with loop / bounce options), numeric labels will be used to space
# frames relative to how long they should last.

# All animation labels will be converted to the range [0,1] (closest
# label to current time position will be drawn during animation construction)

# If there is a 3D plot, the loop around the plot will live on
# animation frames in the range [0,1]

# If frame labels are given with a 3D plot, the rotation will proceed
# normally during the frame label rendering.

# If frames are given one function call at a time, they will be
# collected into existing objects. If multiple frames are provided for
# the same object, all of them will be rendered (raise warning?).

# Legend is based entirely on the unique names. A legend will always
# contain the list of all uniquely named series that will be displayed
# during an animation.

# 3D plot animation trajectory is determined by a [path] and [focus].

# At first, 3D surfaces will use quads generated directly from a
# function. Providing data and rendering a surface will use
# scipy.spatial.Delaunay.




# True if warning should be given to the user at runtime for various
# usages that may result in suboptimal outcomes.
WARNINGS = True

# The unit used that varies with the aspect ration of the figure.
DEFAULT_DIFFERENTIAL_UNIT = "%"

# The unit used uniformly dependent on "size", does not vary with ratio.
DEFAULT_UNIFORM_UNIT = "vh"

# The font used for titles, axes, and information displays.
DEFAULT_FONT = "Monaco"

# Attributes that are numbers with these names are assumed to be
# uniform and will use the "DEFAULT_UNIFORM_UNIT"
UNIFORM_ATTRS = ["stroke", "size"]

# Character replacements needed to translate python attribtes to XML attributes.
ATTR_REPLACEMENTS = {"_":"-"}

# ====================================================================

# Base class for all the XML elements that will be in the SVG.
class Element:
    id = ""
    content = ""

    # Default initialize transfers settings to attributes.
    def __init__(self, **kwargs):
        for k in kwargs: setattr(self, k, kwargs[k])

    # Return the lower-case version of the class name.
    def type(self):
        full_name = str(type(self))
        name_no_prefix = full_name.split(".")[-1]
        name_no_syntax = name_no_prefix.split("'")[0]
        return name_no_syntax.lower()

    # Render the XML string for this element and return it.
    def render(self, **kwargs):
        # Update the keyword arguments in this object before rendering.
        for k in kwargs: setattr(self, k, kwargs[k])
        # Initialize a XML formatted string with keyword placeholders.
        string = "\n<{type}{attributes}>{content}</{type}>\n"
        attributes = [""]
        for attr in dir(self):
            # Skip attributes that will not translate to XML text.
            if (type(getattr(self, attr))) not in (str, int, float): continue
            # Skip built-in / reserved attributes.
            if attr[0] == "_": continue
            # Skip the "content" attribute that is special to rendering.
            if attr == "content": continue
            # Otherwise, add the attribute to the settings.
            attribute = getattr(self, attr)
            # Convert numeric attributes into strings.
            if (type(attribute) != str):
                if any(ua in attr for ua in UNIFORM_ATTRS):
                    attribute = str(attribute) + DEFAULT_UNIFORM_UNIT
                else:
                    attribute = str(attribute) + DEFAULT_DIFFERENTIAL_UNIT
            # Replace characters to translate python attribtes to XML attributes.
            for find in ATTR_REPLACEMENTS:
                attr = attr.replace(find, ATTR_REPLACEMENTS[find])
            # Only append attributes that actually have values.
            if len(attribute) > 0:
                attributes.append( attr +'="'+ attribute + '"' )
        # Return the formatted string that produces XML output.
        return string.format(type=self.type(),
                             attributes=" ".join(attributes),
                             content=self.content)

class SVG(Element):
    version = "1.1"
    xmlns = "http://www.w3.org/2000/svg"
    height = 100.
    width = 100.

class Line(Element):
    x1 = 0.
    x2 = 0.
    y1 = 1.
    y2 = 1.
    stroke = "#000"
    stroke_width = .1

class Rect(Element):
    x = 0.
    y = 0.
    width = 100.
    height = 100.
    fill = "#000"
    
class Text(Element):
    x = 0.
    y = 0.
    fill = "black"
    content = "Text"
    font_size = 1.
    text_anchor = "middle"
    font_family = DEFAULT_FONT


# ====================================================================


# Base class for all "Series" objects that have their own legend entry.
class Series:
    name = None
    x = None
    y = None
    is_3d = None
    animated = None
    
    # Given a function "project" that converts points into 2D
    # coordinates for rendering, generate the SVG text for this series.
    # "start" and "end" are for animation, if this series is animated,
    # then the contents of the render will be determined by the
    # appropriate animated frames.
    def render_data(self, project, start=0., end=1.):
        pass

    # Given the coordinates to place the text entry for this series,
    # render a legend entry for this series at the location (x,y) and
    # make sure the entry has specified max width and height.
    # "start" and "end" are provided in case the legend entry for this
    # series transforms through time.
    def render_legend(self, x, y, start=0., end=1.,
                      max_width=100., max_height=100.):
        pass


# "Line" is a series of points connected by lines.
class Line(Series):
    def __init__(self, name, x, y, z=None, **kwargs):
        self.name = name
        self.x = x
        self.y = y
        if type(z) != type(None):
            self.z = z
            self.is_3d = True
        else:
            self.is_3d = False

        # Copy in any extra keyword arguments.
        for k in kwargs:
            setattr(self, k, kwargs[k])

# "Scatter" is a series of shapes at specified points.
class Scatter(Series):
    def __init__(self, name, x, y, z=None, **kwargs):
        self.name = name
        self.x = x
        self.y = y
        if type(z) != type(None):
            self.z = z
            self.is_3d = True
        else:
            self.is_3d = False

        # Copy in any extra keyword arguments.
        for k in kwargs:
            setattr(self, k, kwargs[k])

# ====================================================================


# This is the object that orchestrates rendering the entire SVG from
# all of the data. It will create the appropriate set of elements and
# size / place them into the scene according to the correct coordinates.
class Layout:
    font_family = DEFAULT_FONT
    title_font_size = 4.
    axis_title_font_size = 3.
    axis_tick_font_size = 2.
    axis_stroke_width = .3
    plot_padding = 2.5
    canvas_bg_color = "rgba(0,0,0,.1)"
    grid_lines = True
    num_grid_lines = [4,8]
    grid_line_stroke_width = ".1vh"
    axis_color = "rgba(0,0,0,.95)"
    # Title controls
    title = False
    x_title = False
    y_title = False
    legend = True
    legend_location = "right"

    # Update the attributes stored in this "Compiler"
    def set(self, **kwargs):
        for k in kwargs: setattr(self, k, kwargs[k])

    # Return the rendered SVG plot as a string.
    def compile(self, svg=None, canvas=None, **kwargs):
        # Initialize an SVG if one was not provided.
        if (type(svg) == type(None)): svg = SVG()
        # Initialzie a plotting canvas if one was not provided.
        if (type(canvas) == type(None)):
            canvas = Rect(fill=self.canvas_bg_color,
                          x=self.plot_padding, y=self.plot_padding,
                          width=(100.-self.plot_padding*2),
                          height=(100.-self.plot_padding*2))
        # Return the rendered SVG string.
        return svg.render(content=canvas.render())

        
# ====================================================================

class Plot():
    def __init__(self, **kwargs):
        self.width = None
        self.height = None
        self.ratio = None
        self.background_color = None
        # Title controls
        self.title_font = None
        self.title_font_size = None
        # Axis controls
        self.axis_font = None
        self.axis_font_size = None
        self.axis_color = None
        self.sub_axis_color = None
        self.sub_axis_ticks = None
        # Legend controls
        self.legend_font = None
        self.legend_font_size = None
        self.legend_position = None
        # Object control
        self.max_objects = 2000
        self.total_objects = 0
        # Container for all data
        self.title = ""
        self.xaxis = "X"
        self.yaxis = "Y"
        self.zaxis = "Z"
        self.data = []
        self.is_3d = False
        # Container for compiler options controlling figure layout.
        self.layout = Layout()
        # Copy in any extra keyword arguments.
        for k in kwargs:
            setattr(self, k, kwargs[k])


    # =================================
    #      Adding data to the plot     
    # =================================

    # Construct a wrapper for the most common "add" operation.
    def add(self, *args, **kwargs):
        return self.add_scatter(*args, **kwargs)

    def add_line(self, name, x, y, z=None):
        pass

    def add_region(self, name, func, x_range, y_range):
        pass

    def add_scatter(self, name, x, y, z=None):
        pass

    def add_surface(self, name, x, y, z):
        pass

    def add_histogram(self, name, x):
        pass

    def add_boxplot(self, name, x, values):
        pass

    def add_pie(self, name, values):
        pass

    def add_annotation(self, contents, x, y, z=None):
        pass

    def add_node(self, label, x, y, z=None):
        pass

    def add_edge(self, labels):
        pass


    # ===========================================
    #      Compiling and Displaying the Plot     
    # ===========================================

    # Render the plot string, save it to file, and show the plot in the browser.
    def plot(self, file_name="out.svg", dir_name=".", show=True, **compile_kwargs):
        # Update the file name to be the absolute path to the file.
        import os
        file_name = os.path.join(os.path.abspath(dir_name), file_name)
        print(f"Saving SVG at path '{file_name}'..")

        # Save the SVG string to file.
        string = self.layout.compile()
        with open(file_name, "w") as f:
            f.write(string)

        # Show the plot in the default web browser.
        if show:
            import webbrowser
            webbrowser.open("file://" + file_name)

# =================
#      Testing     
# =================

if __name__ == "__main__":
    import numpy as np
    x = np.linspace(-3,10,10)
    y = x**3
    # Create a plot
    p = Plot(title="Does this thing work?", xaxis="X VALUES", yaxis="y u do this?")
    p.add("First series", x, y)
    p.plot()
    
