# Import python default packages
import os, sys, time
import numpy as np
from PIL import Image

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


# Function for displaying an image directly from a vector.
# 
#    vec -- A numpy ndarray containing image data of any numeric type.
#  shape -- A tuple representing the shape of the image. Assumed square.
#   gray -- True if color image should be converted to grayscale.
#   clip -- True if values out of range should be clipped, otherwise rescaled.
# fliplr -- True if the image vector should be flipped left-right.
# flipud -- True if the image vector should be flipped up-down.
#   name -- The file name to save the temporary file as.
# 
def show_img(vec, shape=None, gray=False, clip=True, rotate=False,
             fliplr=False, flipud=False, name=TEMP_IMG):
    from PIL import Image
    vec = vec.copy()
    # Try and intelligently infer the shape. Assume square.
    if type(shape) == type(None):
        if (vec.size % 3) == 0:
            shape = (-1, int(round((vec.size//3)**(1/2))), 3)
        else: shape = (-1, int(round((vec.size//3)**(1/2))))
    # Make sure the image is ranged [0,255] everywhere.
    if clip:
        vec[vec < 0] = 0
        vec[vec > 255] = 255
    else:
        vec -= np.min(vec)
        vec /= np.max(vec)
        vec *= 255
    # Make sure the shape of the image is correct.
    vec = vec.reshape(shape)
    # Flip the image over the x and y axis if desired.
    if flipud: vec = np.flipud(vec)
    if fliplr: vec = np.fliplr(vec)
    # Make the image grayscale if that's desired.
    if gray: vec = np.mean(vec, axis=-1)
    # Make the type of the image compatible with PIL.
    vec = np.array(vec, dtype='uint8')
    # Produce an image.
    img = Image.fromarray(vec)
    if rotate: img = img.rotate(rotate)
    img.save(name)
    os.system(f"open {name}")


BLUE = np.array((100,100,255), dtype=np.uint8)
RED = np.array((255,100,100), dtype=np.uint8)
WHITE = np.array((255,255,255), dtype=np.uint8)

# ===============================================================
#     Convert a matrix into an image for display
# 
def matrix_to_img(matrix, file_name=None, directory=".", rescale=False,
                  pos_color=BLUE, mid_color=WHITE, neg_color=RED):
    # If no file name is given, use a temporary file and turn on "show".
    if (type(file_name) == type(None)):
        from tempfile import NamedTemporaryFile
        f = NamedTemporaryFile(delete=False, suffix='.png')
        file_name = f.name
        f.close()
        show = True
    # Handle a 2-d matrix.
    if len(matrix.shape) == 2:
        matrix = np.asarray(matrix, dtype="float64")
        img = np.zeros(matrix.shape + (3,), dtype="uint8")
        max_val = np.max(matrix)
        min_val = abs(np.min(matrix))
        if rescale:
            matrix = matrix - max_val
            matrix /= max_val
            matrix *= 2
            matrix -= 1
            min_val = -1
            max_val = 1
        start = time.time()
        # Rescale the matrix of numbers to range [-1,1]
        for row in range(matrix.shape[0]):
            if (time.time() - start > 1) and (not row%10):
                print(row,":",matrix.shape[0],end="\r")
                start = time.time()
            for col in range(matrix.shape[1]):
                # Rescale positive values to range from white to blue
                if matrix[row,col] >= 0:
                    val = matrix[row,col] / max_val
                    img[row,col,:] = val * pos_color + (1-val) * mid_color
                # Rescale negative values to range from white to red
                else:
                    val = - matrix[row,col] / min_val
                    img[row,col,:] = val * neg_color + (1-val) * mid_color
    elif len(matrix.shape) == 3:
        img = np.array(matrix, dtype="uint8")
    else:
        class UnrecognizedFormat(Exception): pass
        raise(UnrecognizedFormat("Unrecognized matrix shape..",matrix.shape))
    # Convert the scaled values into an actual image and save to desktop.
    img = Image.fromarray(img)
    file_name = os.path.join(directory, file_name)
    print("saving image at",file_name+"..")
    img.save(file_name)
    if show: os.system(f"open {file_name}")
# ===============================================================
