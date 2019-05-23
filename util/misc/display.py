BLUE = np.array((100,100,255), dtype=np.uint8)
RED = np.array((255,100,100), dtype=np.uint8)
WHITE = np.array((255,255,255), dtype=np.uint8)

# ===============================================================
#     Convert a matrix into an image for display
# 
def matrix_to_img(matrix, file_name, directory=".", pos_color=BLUE,
                  mid_color=WHITE, neg_color=RED, show=False):
    import time
    from PIL import Image
    import numpy as np
    if len(matrix.shape) == 2:
        matrix = np.asarray(matrix, dtype="float64")
        img = np.zeros(matrix.shape + (3,), dtype="uint8")
        max_val = np.max(matrix)
        min_val = abs(np.min(matrix))
        start = time.time()
        # Rescale the matrix of numbers to range [-1,1]
        for row in range(matrix.shape[0]):
            if (time.time() - start > 1):
                print(row,":",matrix.shape[0],end="\r")
                start = time.time()
            for col in range(matrix.shape[1]):
                # Rescale positive values to range from white to blue
                if matrix[row,col] >= 0:
                    val = matrix[row,col] / max_val
                    img[row,col,:] = val * pos_color + (1-val) * mid_color
                # Rescale negative values to range from white to red
                else:
                    val = -matrix[row,col] / min_val
                    img[row,col,:] = val * neg_color + (1-val) * mid_color
    elif len(matrix.shape) == 3:
        img = np.array(matrix, dtype="uint8")
    else:
        class UnrecognizedFormat(Exception): pass
        raise(UnrecognizedFormat("Unrecognized matrix shape.."))
    # Convert the scaled values into an actual image and save to desktop.
    img = Image.fromarray(img)
    file_name = os.path.join(directory, file_name)
    print("saving image at",file_name+"..")
    img.save(file_name)
    if show: os.system(f"open {file_name}")
# ===============================================================
