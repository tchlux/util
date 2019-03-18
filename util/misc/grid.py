
# Given a "min" a "max" (values or vectors of length "dim"), generate
# a full grid of points with as near to "num" positions as possible.
# Return results in a numpy matrix.
def grid(min=0, max=1, dim=1, num=None):
    import numpy as np
    # Default value for "num" is 10 points per dimension. Convert the
    # "num" into the number of points per dimension.
    if (num == None): num = 11
    else:             num = int(round(num**(1/dim)))
    vals = []
    # Set the minimum and maximum to be arrays of length "dim".
    try:    len(min)
    except: min = np.ones(dim)*min
    try:    len(max)
    except: max = np.ones(dim)*max
    # Cycle through and generate grid points.
    for low,upp in zip(min,max):
        vals.append( np.linspace(low, upp, num=num) )
    # Combine points into a meshgrid.
    vals = np.meshgrid(*vals)
    # Return the grid of point values.
    return np.vstack([v.flatten() for v in vals]).T


if __name__ == "__main__":
    pts = grid(dim=6, num=300)
    print(pts.shape)
    print()
    print(pts)
