import numpy as np

# Reduce a number to its minimal display form (no unnecessary 0's on right)
CLEAN_NUMBER_STRING = lambda number: str(number).rstrip("0").rstrip(".") \
                      if "." in str(number) else str(number)

# Print a numpy array space seperated so that it can be taken by fortran
def CLEAN_STRING(arr, string=""):
    for num in arr: string = string + " " + CLEAN_NUMBER_STRING(num)
    return string

CLEAN_ARRAY_STRING = lambda arr: (
    "\n".join([CLEAN_ARRAY_STRING(row) for row in arr]) if
    (len(arr.shape) > 1) else CLEAN_STRING(arr) )


