import numpy as np

DEFAULT_MAX_STR_LEN = 30 # max column width for printouts
DEFAULT_MAX_DISPLAY = 10 # max number of rows displayed for various printouts
DEFAULT_WAIT = 3. # seconds to wait before displaying progress
UPDATE_FREQUENCY = .5  # seconds between updates when reading data
MAX_ERROR_PRINTOUT = 10 # max errors printed when processing data
FILE_SAMPLE_SIZE = 2**20 # <- Default amount of files read by Data (in bytes).
MISSING_SAMPLE_SIZE = 100000 # <- The max size for which "missing" data
#                                 calculations become an estimate instead
#                                 of exact numbers of missing values.
MAX_DISPLAY_COLS = 20 # <- maximum number of displayed columnns from data on read
EXTRA_CATEGORY_KEY = "Data - other" # <- key used for values that are "other" encoded

# Definition of separators by file extension when saving.
COMMON_SEPARATORS = [",", " ", "\t", ";"]
SEPARATORS = {
    "csv" : ",",
    "txt" : " ",
    "tsv" : "\t",
}

# Allowable "quote" characters in files.
QUOTES = {'"'}

# Declare some constants.
NP_TYPES = {str:(np.str_, 16),    # Which numpy types correspond to 
            int:(np.int64,),      # the python types of variables
            float:(np.float64,)}
PY_TYPES = {value[0]:key for (key,value) in NP_TYPES.items()}
# 
GENERATOR_TYPE = type(_ for _ in ())

from .data import Data
from .video import IndexedVideo

# Automatically execute testing code?
# from util.data.test import *
