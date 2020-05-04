# Compile and build the fekete point generation code.
import os
import og_fmodpy as fmodpy
CWD = os.path.dirname(os.path.abspath(__file__))
fp_mod = fmodpy.fimport(os.path.join(CWD,"fekete.f90"),
                        module_link_args=["-lblas","-llapack"],
                        output_directory=CWD)


# Given an "n", construct the 1D Chebyshev-Guass-Lobatto nodes
# (equally spaced on a unit semicircle).
def chebyshev(n, d=1, sample=None):
    from numpy import array, cos, pi
    x = []
    for k in range(1,n+1):
        x.append(cos( (2*k - 1) / (2*n) * pi ))
    # Reormalize the points to be in the range [0,1].
    x = (array(x) + 1) / 2
    # Create a tensor product to fill the desired dimension.
    return mesh(x, multiplier=d, sample=sample)


# Compute a grid of "n" points in "d" dimensions. Randomly subsample
# the full grid if it is larger than "n".
def grid(n, d):
    from random import sample
    per_dim = round( n**(1/d) + .5 )
    pts = mesh(list(range(per_dim)), d) / (per_dim-1)
    if n < len(pts): return pts[sample(range(len(pts)), n)]
    else:            return pts

# Create a tensor mesh of "len(nodes) * multiplier" dimension of all
# provided lists inside of "nodes". If "sample" is provided,
# then a random sample of size "sample" is drawn from that mesh.
# 
#    `mesh([0, 1], 2) = [[0, 0], [1, 0], [0, 1], [1, 1]]`
# 
def mesh(*nodes, multiplier=None, sample=None):
    # Extract the passed `sample` and `multiplier` from the end of the
    # argument list if they were not given as keyword arguments.
    if (sample is None) and (type(nodes[-1]) == bool):
        nodes, sample = nodes[:-1], nodes[-1]
    if (multiplier is None) and (type(nodes[-1]) == int):
        nodes, multiplier = nodes[:-1], nodes[-1]
    # Start of the actual function.
    from numpy import vstack, meshgrid
    if sample:
        # Compute the multiplicative product of a list of integers.
        # Make sure everything is PYTHON INTEGER to avoid overflow.
        def product(integers):
            p = int(integers[0])
            for i in range(1,len(integers)): p *= int(integers[i])
            return p
        mesh_size = product(list(map(len,nodes))) ** int(multiplier)
        # Only continue if the sample does *not* include the whole grid.
        if (sample < mesh_size):
            from util.random import random_range
            from numpy import zeros
            # Expand out the full set of nodes (along each component).
            nodes *= multiplier
            # Compute the total mesh size (and track the size of each component).
            sizes = list(map(len, nodes))
            # Randomly draw points from the mesh.
            points = zeros((sample, len(sizes)), dtype=float)
            for i,index in enumerate(random_range(mesh_size, count=sample)):
                smaller_prod = mesh_size
                # Decode the index into a set of indices for each component.
                for s in range(len(sizes)):
                    smaller_prod = int(smaller_prod // sizes[s])
                    curr_ind = int(index // smaller_prod )
                    # Handle the special case where remaining sizes are 1.
                    if (smaller_prod == 1):
                        curr_ind = min(curr_ind, sizes[s]-1)
                    # Subtract out the effect of the chosen index.
                    index -= curr_ind * smaller_prod
                    points[i,s] = nodes[s][curr_ind]
            # Return the randomly selected mesh points.
            return points
    # Default operation, return the full tensor product of "nodes".
    return vstack([_.flatten() for _ in meshgrid(*(nodes * multiplier))]).T

# Given a degree and a dimension, construct a list of polynomial functions.
def polynomials(degree=1, dimension=1):
    import numpy as np
    from itertools import combinations_with_replacement as combinations
    # Initialize with only a constant function.
    funcs = [lambda x: 1]
    # Construct polynomials of all degree.
    for deg in range(1,degree+1):
        for group in combinations(range(dimension), deg):
            f = lambda x, g=np.array(group): np.prod(x[g])
            funcs.append( f )
    return funcs

# Given a degree and a dimension, construct a list of polynomial
# functions in terms of the degree used in the polynomial and the indices.
def polynomial_indices(degree=1, dimension=1):
    import numpy as np
    from itertools import combinations_with_replacement as combinations
    # Initialize with only a constant function.
    degrees = [ 0 ]
    indices = [ [-1] * degree ]
    # Construct polynomials of all degree.
    for deg in range(1,degree+1):
        for group in combinations(range(dimension), deg):
            degrees.append( deg )
            indices.append( list(group) + [-1] * (degree - deg) )
    return ( np.asfortranarray(np.array(degrees, dtype=int)),
             np.asfortranarray(np.array(indices, dtype=int).T) )

# Construct Fekete points over the unit box. This operation will
# destroy the contents of "vmt" (Vandermonde matrix transpose).
def fekete_indices(vmt):
    import numpy as np
    # If the system is not overdetermined, then all points should be kept.
    if (vmt.shape[1] <= vmt.shape[0]): return np.arange(vmt.shape[1])
    # Otherwise, identify which indices should be kept by performing a
    # QR with column pivoting (and track which columns are first).
    to_keep, info = fp_mod.fekete_indices(vmt)
    if (info != 0):
        import warnings
        warning.warn(f"FEKETE POINT QR ERROR:  {info}")
    # Return the smaller vandermonde matrix (with fewer points).
    return to_keep - 1

# Construct Fekete points over the unit hypercube.
def fekete_points(num_points, dimension, min_per_dim=3,
                  max_func_ratio=1, max_point_ratio=2):
    from math import log, ceil
    import numpy as np
    assert(dimension >= 1)
    # Construct a set of functions at least as long as the number of points.
    degree = 1
    # Special case where we know exactly the necessary degree.
    if (dimension == 1): degree = num_points
    funcs = polynomials(degree, dimension)
    # Otherwise, search iteratively (should require few iterations).
    while (len(funcs) < num_points):
        del funcs
        degree += 1
        funcs = polynomials(degree, dimension)
    # Generate a warning if the user is requesting a potentially unstable point set.
    if (degree > 10):
        import warnings
        class UnstablePolynomial(Warning): pass
        warnings.warn(UnstablePolynomial("Polynomials of high degree (>10) are notoriously unstable. Use resulting Fekete points with caution."))
    # Get the set of polynomial degrees and the indices (for fast
    # Fortran evaluation of the Vandermonde matrix).
    degrees, indices = polynomial_indices(degree, dimension)
    # Find the minimum number of points that should be drawn per
    # dimension to satisfy "num_points", but draw at least
    # "min_per_dim" points per dimension.
    points = chebyshev(max(min_per_dim,ceil(num_points**(1/dimension))),
                       dimension, sample=max_point_ratio*num_points)
    # Only keep the necessary basis functions for making a Vandermonde matrix.
    if (len(funcs) > max_func_ratio * num_points):
        to_keep = max_func_ratio * num_points
        funcs = funcs[:to_keep]
        degrees = np.array(degrees[:to_keep], order='f')
        indices = np.array(indices[:,:to_keep], order='f')
    # Create a warning when a large memory footprint is imminent.
    if (len(points)*len(degrees) * 8) > (1 * 2**30):
        gb_size = (len(points)*len(degrees) * 8) / (2**30)
        shape = f"(matrix size {len(points)} x {len(degrees)})"
        import warnings
        class LargeMemoryFootprint(Warning): pass
        warnings.warn(LargeMemoryFootprint(f"\n  The memory footprint will be large (~{gb_size:.0f}GB) while constructing these Fekete points.\n  The following QR factorization {shape} will likely take time (minutes) as well.\n"))
    # Allocate space for Vandermonde matrix (tranpose)
    vmt = np.zeros((len(degrees), len(points)), dtype=float, order='f')
    points = np.asfortranarray(points.T) # <- this should require no copy
    # Construct the Vandermonde matrix (transpose).
    fp_mod.evaluate_vandermonde(points, degrees, indices+1, vmt)
    f_indices = fekete_indices(vmt)
    del vmt
    to_keep = f_indices[:num_points]
    del f_indices
    f_points = points[:,to_keep].T
    del points
    # Find the "Fekete" points.
    return f_points


# Given "d" categories that need to be converted into real space,
# generate a regular simplex in (d-1)-dimensional space. (all points
# are equally spaced from each other and the origin) This process
# guarantees that all points are not placed on a sub-dimensional
# manifold (as opposed to "one-hot" encoding).
def regular_simplex(d, volume=None):
    import numpy as np
    # Special cases for one and two categories
    if d < 1:
        class InvalidDimension(Exception): pass
        raise(InvalidDimension("Dimension must be an integer greater than 0."))
    elif d == 1: return np.array([[0.],[1.]])
    # Initialize all points to be zeros.
    points = np.zeros((d+1,d))
    # Set the initial first point as 1
    points[0,0] = 1
    # Calculate all the intermediate points
    for i in range(1,d):
        # Set all points to be flipped from previous calculation while
        # maintaining the angle "arcos(-1/d)" between vectors.
        points[i:,i-1] = -1/(1+d-i) * points[i-1,i-1]
        # Compute the new coordinate using pythagorean theorem
        points[i,i] = (1 - sum(points[i,:i]**2))**(1/2)
    # Set the last coordinate of the last point as the negation of the previous
    points[i+1,i] = -points[i,i]
    # Normalize the edge lengths to all be 1.
    points *= (2 + 2/d)**(-1/2)
    # If no volume was specified, return the vertices of the regular simplex.
    if volume is None: return points
    elif (volume > 0):
        # Return a regular simplex with volume specified by first converting
        # to unit edge length, then using the formula from [1]:
        # 
        #   volume  =   (length^d / d!) ([d+1] / 2^d)^(1/2)
        # 
        # where "volume" is the volume of the simplex defined by d+1
        # vertices, "d" is the dimension the simplex lives in (number
        # of vertices minus one), and "length" is the length of all
        # edges between vertices of the simplex. This gives:
        # 
        #  length  =  [ volume d! ([d+1] / 2^d)^(-1/2) ]^(1/d)
        #          =  (volume d!)^(1/d) ([d+1] / 2^d)^(-1/(2d))
        # 
        # Reference:
        #  [1] https://www.researchgate.net/publication/267127613_The_Volume_of_an_n-Simplex_with_Many_Equal_Edges
        # 
        from math import factorial
        edge_length = (volume*factorial(d))**(1/d) * ((d+1) / 2**d)**(-1/(2*d))
        # get_v = lambda l: (l**d / factorial(d)) * ((d+1) / 2**d)**(1/2)
        # get_l = lambda v: (v*factorial(d))**(1/d) * ((d+1) / 2**d)**(-1/(2*d))
        # print("Estimated volume: ",get_v(edge_length))
        return points * edge_length
        

# Construct a mesh formed from equilateral simplices grown out from
# the middle, attempt to have nearly "n" points in a "d" dimensional
# unit hypercube.
def regular_mesh(n, d):
    from math import factorial
    # Get a regular simplex with a prescribed volume.
    seed_simplex = regular_simplex(d, volume=1/n)
    # Shrink the simplex so that the growing procedure will create
    # approximately "n" points in total.
    import numpy as np
    print("Edge lengths:")
    print("",np.round(np.linalg.norm(seed_simplex[1:] - seed_simplex[0],axis=1),2))
    print()
    print("Vertices:")
    for i in range(len(seed_simplex)):
        print("", str(np.round(seed_simplex[i], 3)).replace("\n "," "))
    print()
    if d in {2,3}:
        from util.plot import Plot
        p = Plot()
        p.add("simplex", *(seed_simplex.T), mode="markers+lines", color=(0,0,0,.5))
        p.show(x_range=[-1, 1], y_range=[-1,1], z_range=[-1,1],
               width=500, height=500, show_legend=False)

    # Find a way to rotate simplex so it fits inside the unit hypercube,
    #   will be tough in high dimension, but should be possible by pointing
    #   long edges in the direction of many simultaneous components.
    # At that point, all vertices should be inside the unit hypercube.
    # Queue up all faces (and opposite vertices).
    # For each flip in queue, record the new vertex and all face
    #   vertices as new faces that need to be flipped across (if new
    #   vertex is in cube).



if __name__ == "__main__":
    regular_mesh(1, 20)
    # regular_mesh(10e10, 20)

    show = False

    def _test_mesh(show=True):
        from weakly_admissable_meshes import polar_wam, box_wam
        from util.plot import Plot
        p = Plot("Weakly admissable mesh")
        p.add("Polar", *(polar_wam(8).T), color=p.color(0))
        p.add("Box", *(box_wam(8).T), color=p.color(1))
        if show: p.show(width=1000*.7, height=900*.7, append=True)


    def _test_num_polynomials(dim=2):
        print()
        print("Num of polynomial functions: ")
        for d in range(1, 20):
            print("",f"degree {d}: {len(polynomials(dimension=dim, degree=d))}")
        print()

    def _test_fekete(max_n=2000, max_d=1000, start_n=1000, start_d=10, steps=2, show=True):
        if show:
            # Generate some points in 3 dimensions for showing.
            pts = fekete_points(47, 3)
            from util.plot import Plot
            p = Plot("Fekete")
            p.add("fekete points", *(pts.T))
            p.show(append=True)

        from util.system import Timer
        t = Timer()
        steps -= 1
        for n in range(start_n,max_n+1, (max_n-start_n-1) // steps):
            for d in range(start_d, max_d+1, (max_d-start_d-1) // steps):
                t.start()
                print()
                print('-'*70)
                print(f"  {n}, {d}", flush=True)
                pts = fekete_points(n, d)
                print(f"     {pts.shape}")
                print(f"     {t()} seconds", flush=True)
                t.stop()
        print()
        print()



    # _test_mesh(show=show)
    # _test_num_polynomials()
    # _test_fekete(show=show)
