import numpy as np

# Given "n" points and "d" dimensions, generate a weakly admissible
# mesh. WARNING: Forced to dimension = 2.
def polar_wam(n, dimension=2):
    # Initialize points only containing the center.
    points = [[0,0]]
    for length in range(1,n):
        l = 1/2 * (1 + np.cos((length * np.pi) / n))
        for angle in range(0,2*n+1):
            a = (2 * np.pi * angle) / (2 * n + 1)
            # Convert the length and angle into a 2D point.
            x = l * np.cos(a)
            y = l * np.sin(a)
            points.append( [x,y] )
    # Convet the points into a row-vector matrix and shift into the unit box.
    return (np.array(points) + 1) / 2

# Given "n" points and "d" dimensions, generate a weakly admissible
# mesh. WARNING: Forced to dimension = 2.
def box_wam(n, dimension=2):
    # Initialize points only containing the center.
    points = []
    for length_1 in range(0,n+1):
        l1 = 1/2 * (1 + np.cos((length_1 * np.pi) / n))
        for length_2 in range(0,n+1):
            l2 = 1/2 * (1 + np.cos((length_2 * np.pi) / n))
            # Convert the length and angle into a 2D point.
            points.append( [l1,l2] )
    # Convet the points into a row-vector matrix.
    return np.array(points)
