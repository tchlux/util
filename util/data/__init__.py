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

FILE_SAMPLE_SIZE = 2**20 # <- Default amount of files read by Data (in bytes).
MISSING_SAMPLE_SIZE = 100000 # <- The max size for which "missing" data
#                                 calculations become an estimate instead
#                                 of exact numbers of missing values.
MAX_DISPLAY_COLS = 20 # <- maximum number of displayed columnns from data on read

is_none = lambda v: type(v) == type(None)

# coding: future_fstrings
# from future import print_function

from util.data.data import Data
from util.data.video import IndexedVideo
from util.data.categories import regular_simplex, category_ratio
# from util.data.test import *
