import os 
import numpy as np

# ===============================================================
#     Convert a matrix into an image for display
# 
# Hard coded colors
BLUE = np.array((100,100,255), dtype=np.uint8)
RED = np.array((255,100,100), dtype=np.uint8)
WHITE = np.array((255,255,255), dtype=np.uint8)
BLACK = np.array((0,0,0), dtype=np.uint8)
# File interactions
DIRECTORY = os.path.abspath(".")
TEMP_IMG = "temp.png"
def matrix_to_img(matrix, file_name=TEMP_IMG, directory=DIRECTORY,
                  pos_color=WHITE, mid_color=BLACK, neg_color=RED, show=False):
    from PIL import Image
    if len(matrix.shape) == 2:
        matrix = np.asarray(matrix, dtype="float64")
        img = np.zeros(matrix.shape + (3,), dtype=np.uint8)
        max_val = np.max(matrix)
        min_val = abs(np.min(matrix))
        # Rescale the matrix of numbers to range [-1,1]
        for row in range(matrix.shape[0]):
            for col in range(matrix.shape[1]):
                # Rescale positive values to range from white to blue
                if matrix[row,col] >= 0:
                    val = matrix[row,col] / max_val
                    img[row,col,:] = val * pos_color + (1-val) * mid_color
                # Rescale negative values to range from white to red
                else:
                    val = -matrix[row,col] / min_val
                    img[row,col,:] = val * neg_color + (1-val) * mid_color
    elif len(matrix.shape) == 3: img = matrix.copy()
    else: raise(Exception(f"Unrecognized matrix shape, {matrix.shape}."))
    # Convert the scaled values into an actual image and save to desktop.
    img = Image.fromarray(img)
    file_name = os.path.join(directory, file_name)
    img.save(file_name)
    if show: os.system(f"open {file_name}")
# ===============================================================


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
