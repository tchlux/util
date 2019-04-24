# Define the path to the file that holds the worker program for "capture".
import os
CAPTURE_WORKER = os.path.join(os.path.dirname(os.path.abspath(__file__)),"capture_worker.py")

# Import all other decorators.
from util.decorators.decorators import *
