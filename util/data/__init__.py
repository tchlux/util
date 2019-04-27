import numpy as np

QUOTES = {'"'}
DEFAULT_DISPLAY_WAIT = 3.
COMMON_SEPARATORS = [",", " ", "	", ";"]
UPDATE_FREQUENCY = .5  # How much time (in seconds) between updates when reading data.
MAX_ERROR_PRINTOUT = 10 # Only print out this many errors when processing data

# Definition of separators by file extension when saving.
SEPARATORS = {"csv":",", "tsv":"\t", "ssv":" "}

# Declare some constants.
NP_TYPES = {str:(np.str_, 16),    # Which numpy types correspond to 
            int:(np.int64,),      # the python types of variables
            float:(np.float64,)}
PY_TYPES = {value[0]:key for (key,value) in NP_TYPES.items()}
# 
GENERATOR_TYPE = type(_ for _ in ())

FILE_SAMPLE_SIZE = 2**20
MISSING_SAMPLE_SIZE = 10000

is_none = lambda v: type(v) == type(None)

# coding: future_fstrings
# from future import print_function

from util.data.data import Data
# from util.data.image import *
# from util.data.video import *
# from util.data.categories import *
# from util.data.read import *
# from util.data.test import *
