'''This Python code is an automatically generated wrapper
for Fortran code made by 'fmodpy'. The original documentation
for the Fortran source code follows.


'''

import os
import ctypes
import numpy

# --------------------------------------------------------------------
#               CONFIGURATION
# 
verbose = False
fort_compiler = "gfortran"
shared_object_name = "weight.so"
this_directory = os.path.dirname(os.path.abspath(__file__))
path_to_lib = os.path.join(this_directory, shared_object_name)
compile_options = "-fPIC -shared -O3 -fopenmp"
ordered_dependencies = ['weight.f90', 'weight_bind_c.f90']
# 
# --------------------------------------------------------------------
#               AUTO-COMPILING
#
# Try to import the existing object. If that fails, recompile and then try.
try:
    clib = ctypes.CDLL(path_to_lib)
except:
    # Remove the shared object if it exists, because it is faulty.
    if os.path.exists(shared_object_name):
        os.remove(shared_object_name)
    # Compile a new shared object.
    command = " ".join([fort_compiler, compile_options]
                        + ordered_dependencies + ["-o", path_to_lib])
    if verbose:
        print("Running command")
        print("  ", command)
    # Run the compilation command.
    import subprocess
    subprocess.run(command, shell=True, cwd=this_directory)
    # Remove all ".mod" files that were created to reduce clutter.
    all_mods = [f for f in os.listdir(os.curdir) if f[-4:] == ".mod"]
    for m in all_mods: os.remove(m)
    # Import the shared object file as a C library with ctypes.
    clib = ctypes.CDLL(path_to_lib)

# --------------------------------------------------------------------


# ----------------------------------------------
# Wrapper for the Fortran subroutine WEIGHT

def weight(pts, pt, wts=None):
    ''''''
    
    # Setting up "pts"
    if ((not issubclass(type(pts), numpy.ndarray)) or
        (not numpy.asarray(pts).flags.f_contiguous) or
        (not (pts.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'pts' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        pts = numpy.asarray(pts, dtype=ctypes.c_double, order='F')
    pts_dim_1 = ctypes.c_int(pts.shape[0])
    pts_dim_2 = ctypes.c_int(pts.shape[1])
    
    # Setting up "pt"
    if ((not issubclass(type(pt), numpy.ndarray)) or
        (not numpy.asarray(pt).flags.f_contiguous) or
        (not (pt.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'pt' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        pt = numpy.asarray(pt, dtype=ctypes.c_double, order='F')
    pt_dim_1 = ctypes.c_int(pt.shape[0])
    
    # Setting up "wts"
    if (wts is None):
        wts = numpy.zeros(shape=(pts.shape[1]), dtype=ctypes.c_double, order='F')
    wts_dim_1 = ctypes.c_int(wts.shape[0])

    # Call C-accessible Fortran wrapper.
    clib.c_weight(ctypes.byref(pts_dim_1), ctypes.byref(pts_dim_2), ctypes.c_void_p(pts.ctypes.data), ctypes.byref(pt_dim_1), ctypes.c_void_p(pt.ctypes.data), ctypes.byref(wts_dim_1), ctypes.c_void_p(wts.ctypes.data))

    # Return final results, 'INTENT(OUT)' arguments only.
    return wts

