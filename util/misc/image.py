# Import python default packages
import sys

# Define globals
PERSPECTIVE_RANGE = [-.2, .2] # Range of real-world replicable changes
ROTATION_RANGE = [-20, 20] # Steepest roads in America
SCALE_RANGE = [-3,5] # Some common and reasonable zoom ranges

# Define some errors that might be raised in this code
class CommandError(Exception): pass
class MissingPackage(Exception): pass

# Check the python version
if (sys.version_info[0] < 3):
    print("WARNING: This code was written and tested for Python3 only.")

#      Import External Packages     
# ==================================
try:
    import numpy as np
except:
    raise(MissingPackage("Missing required pacakage 'numpy'. Use 'pip install --user numpy'"))
try:
    from PIL import Image
except:
    raise(MissingPackage("Missing required pacakage 'pillow'. Use 'pip install --user Pillow'"))

# Find the coefficients of the trasnformation matrix necessary to
# perform an affine trasnformation that maps the points in
# "input_points" to the points in "output_points"
def find_coeffs(input_points, output_points):
    matrix = []
    for p1, p2 in zip(input_points, output_points):
        matrix.append([p1[0], p1[1], 1, 0, 0, 0, -p2[0]*p1[0], -p2[0]*p1[1]])
        matrix.append([0, 0, 0, p1[0], p1[1], 1, -p2[1]*p1[0], -p2[1]*p1[1]])
    A = np.matrix(matrix, dtype=np.float)
    B = np.array(output_points).reshape(8)
    res = np.dot(np.linalg.inv(A.T * A) * A.T, B)
    return np.array(res).reshape(8)

# Perspective change matrix generator
def perspective(img, magnitude):
    left_change = right_change = 0
    if magnitude < 0:
        left_change = abs(img.size[1] * magnitude)
    elif magnitude > 0:
        right_change = abs(img.size[1] * magnitude)    
    input_points = [(0,0+left_change),
                    (img.size[0],0+right_change),
                    (img.size[0],img.size[1]-right_change),
                    (0,img.size[0]-left_change),]
    output_points = [(0,0),
                     (img.size[0],0),
                     (img.size[0], img.size[1]),
                     (0,img.size[0]),]
    return find_coeffs(input_points, output_points)

def random_transformation(img, p_range=PERSPECTIVE_RANGE, 
                          r_range=ROTATION_RANGE, s_range=SCALE_RANGE):
    p_change = np.random.random()*(np.diff(p_range)[0])+p_range[0]
    r_change = np.random.random()*(np.diff(r_range)[0])+r_range[0]
    s_change = np.random.random()*(np.diff(s_range)[0])+s_range[0]
    # Transform perspective
    original_size = img.size
    # Add an alpha layer in case it doesn't exist (so margins can be hidden)
    img.putalpha(255)
    # Change perspective and rotate the image randomly
    img = img.transform(img.size, Image.PERSPECTIVE,
                        perspective(img,p_change))
    img = img.rotate(r_change, expand=True)
    img = img.resize(original_size)
    # Scale the image randomly
    if s_change < 1: s_change = 1 / (2-s_change)
    new_size = tuple(int(v+0.5) for v in s_change * np.array(img.size))
    new_img = img.resize(new_size)
    # Generate a maximally sized output image (for consistency)
    out_img = Image.new("RGBA", tuple(int(i*s_range[-1]+.5) for i in img.size))
    # Generate a random location in which to put the new random image
    shift_x = np.random.randint(0,max(1,out_img.size[0]-new_size[0]))
    shift_y = np.random.randint(0,max(1,out_img.size[1]-new_size[1]))
    out_img.paste(new_img, (shift_x, shift_y))
    # Return the output image
    return out_img
