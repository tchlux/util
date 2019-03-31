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
