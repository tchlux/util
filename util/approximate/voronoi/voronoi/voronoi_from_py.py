'''This Python code is an automatically generated wrapper
for Fortran code made by 'fmodpy'. The original documentation
for the Fortran source code follows.

!==========================================
!     Full Naive Evaluation of Voronoi
!==========================================
'''

import os
import ctypes
import numpy

# --------------------------------------------------------------------
#               CONFIGURATION
# 
verbose = False
fort_compiler = "gfortran"
shared_object_name = "voronoi.so"
this_directory = os.path.dirname(os.path.abspath(__file__))
path_to_lib = os.path.join(this_directory, shared_object_name)
compile_options = "-fPIC -shared -O3 -fopenmp"
ordered_dependencies = ['voronoi.f90', 'voronoi_bind_c.f90']
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
# Wrapper for the Fortran subroutine EVAL_MESH

def eval_mesh(int_prods, app_int_prods, mesh=None):
    ''''''
    
    # Setting up "int_prods"
    if ((not issubclass(type(int_prods), numpy.ndarray)) or
        (not numpy.asarray(int_prods).flags.f_contiguous) or
        (not (int_prods.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'int_prods' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        int_prods = numpy.asarray(int_prods, dtype=ctypes.c_double, order='F')
    int_prods_dim_1 = ctypes.c_int(int_prods.shape[0])
    int_prods_dim_2 = ctypes.c_int(int_prods.shape[1])
    
    # Setting up "app_int_prods"
    if ((not issubclass(type(app_int_prods), numpy.ndarray)) or
        (not numpy.asarray(app_int_prods).flags.f_contiguous) or
        (not (app_int_prods.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'app_int_prods' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        app_int_prods = numpy.asarray(app_int_prods, dtype=ctypes.c_double, order='F')
    app_int_prods_dim_1 = ctypes.c_int(app_int_prods.shape[0])
    app_int_prods_dim_2 = ctypes.c_int(app_int_prods.shape[1])
    
    # Setting up "mesh"
    if (mesh is None):
        mesh = numpy.zeros(shape=(app_int_prods.shape[0], int_prods.shape[0]), dtype=ctypes.c_double, order='F')
    mesh_dim_1 = ctypes.c_int(mesh.shape[0])
    mesh_dim_2 = ctypes.c_int(mesh.shape[1])

    # Call C-accessible Fortran wrapper.
    clib.c_eval_mesh(ctypes.byref(int_prods_dim_1), ctypes.byref(int_prods_dim_2), ctypes.c_void_p(int_prods.ctypes.data), ctypes.byref(app_int_prods_dim_1), ctypes.byref(app_int_prods_dim_2), ctypes.c_void_p(app_int_prods.ctypes.data), ctypes.byref(mesh_dim_1), ctypes.byref(mesh_dim_2), ctypes.c_void_p(mesh.ctypes.data))

    # Return final results, 'INTENT(OUT)' arguments only.
    return mesh


# ----------------------------------------------
# Wrapper for the Fortran subroutine MAKE_HUGE

def make_huge(dots):
    ''''''
    
    # Setting up "dots"
    dots_dim_1 = ctypes.c_int(dots.shape[0])
    dots_dim_2 = ctypes.c_int(dots.shape[1])

    # Call C-accessible Fortran wrapper.
    clib.c_make_huge(ctypes.byref(dots_dim_1), ctypes.byref(dots_dim_2), ctypes.c_void_p(dots.ctypes.data))

    # Return final results, 'INTENT(OUT)' arguments only.
    return dots


# ----------------------------------------------
# Wrapper for the Fortran subroutine PREDICT

def predict(points, dots, eval_pt, weights=None):
    '''! Given column-vector POINTS, pairwise dot product matrix DOTS, and
! a single evaluation point EVAL_PT, compute the WEIGHTS associated
! with each input point.'''
    
    # Setting up "points"
    if ((not issubclass(type(points), numpy.ndarray)) or
        (not numpy.asarray(points).flags.f_contiguous) or
        (not (points.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'points' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        points = numpy.asarray(points, dtype=ctypes.c_double, order='F')
    points_dim_1 = ctypes.c_int(points.shape[0])
    points_dim_2 = ctypes.c_int(points.shape[1])
    
    # Setting up "dots"
    if ((not issubclass(type(dots), numpy.ndarray)) or
        (not numpy.asarray(dots).flags.f_contiguous) or
        (not (dots.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'dots' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        dots = numpy.asarray(dots, dtype=ctypes.c_double, order='F')
    dots_dim_1 = ctypes.c_int(dots.shape[0])
    dots_dim_2 = ctypes.c_int(dots.shape[1])
    
    # Setting up "eval_pt"
    if ((not issubclass(type(eval_pt), numpy.ndarray)) or
        (not numpy.asarray(eval_pt).flags.f_contiguous) or
        (not (eval_pt.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'eval_pt' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        eval_pt = numpy.asarray(eval_pt, dtype=ctypes.c_double, order='F')
    eval_pt_dim_1 = ctypes.c_int(eval_pt.shape[0])
    
    # Setting up "weights"
    if (weights is None):
        weights = numpy.zeros(shape=(points.shape[1]), dtype=ctypes.c_double, order='F')
    weights_dim_1 = ctypes.c_int(weights.shape[0])
    
    # Setting up "error"
    error = ctypes.c_int()

    # Call C-accessible Fortran wrapper.
    clib.c_predict(ctypes.byref(points_dim_1), ctypes.byref(points_dim_2), ctypes.c_void_p(points.ctypes.data), ctypes.byref(dots_dim_1), ctypes.byref(dots_dim_2), ctypes.c_void_p(dots.ctypes.data), ctypes.byref(eval_pt_dim_1), ctypes.c_void_p(eval_pt.ctypes.data), ctypes.byref(weights_dim_1), ctypes.c_void_p(weights.ctypes.data), ctypes.byref(error))

    # Return final results, 'INTENT(OUT)' arguments only.
    return dots, weights, error.value

