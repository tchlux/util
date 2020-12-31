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
shared_object_name = "modified_shepard.so"
this_directory = os.path.dirname(os.path.abspath(__file__))
path_to_lib = os.path.join(this_directory, shared_object_name)
compile_options = "-fPIC -shared -O3"
ordered_dependencies = ['modified_shepard.f95', 'modified_shepard_bind_c.f90']
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


class modified_shepard:
    ''''''

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine SHEPMOD
    
    def shepmod(self, m, n, x, rw=None):
        '''! This subroutine computes a set of parameters defining a function that
! interpolates N interpolation nodes at scattered X(i) in M dimensions. The
! interpolant may be evaluated at an arbitrary point by the function SHEPMODVAL.
! The interpolation scheme is a modified Shepard method, and can
! use either the Shepard weighted fit or standarad shepard method.
!
! Input parameters:
!   M is the dimension of the data.
!   N is the number of nodes.
!   X(M, N) contains the coordinates of the nodes.
! Output parameters:
!   RW(N) is an array containing the radius of influence for each point.
!   IER
!   = 0, if no errors were encountered.
!   = 1, if N is too small relative to M.'''
        
        # Setting up "m"
        if (type(m) is not ctypes.c_int): m = ctypes.c_int(m)
        
        # Setting up "n"
        if (type(n) is not ctypes.c_int): n = ctypes.c_int(n)
        
        # Setting up "x"
        if ((not issubclass(type(x), numpy.ndarray)) or
            (not numpy.asarray(x).flags.f_contiguous) or
            (not (x.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'x' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            x = numpy.asarray(x, dtype=ctypes.c_double, order='F')
        x_dim_1 = ctypes.c_int(x.shape[0])
        x_dim_2 = ctypes.c_int(x.shape[1])
        
        # Setting up "rw"
        if (rw is None):
            rw = numpy.zeros(shape=(n), dtype=ctypes.c_double, order='F')
        rw_dim_1 = ctypes.c_int(rw.shape[0])
        
        # Setting up "ier"
        ier = ctypes.c_int()
    
        # Call C-accessible Fortran wrapper.
        clib.c_shepmod(ctypes.byref(m), ctypes.byref(n), ctypes.byref(x_dim_1), ctypes.byref(x_dim_2), ctypes.c_void_p(x.ctypes.data), ctypes.byref(rw_dim_1), ctypes.c_void_p(rw.ctypes.data), ctypes.byref(ier))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return rw, ier.value

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine SHEPMODVAL
    
    def shepmodval(self, xp, m, n, x, rw, wts=None):
        '''! SHEPMODVAL returns the modified Shepard approximation at the
! point XP, using the points computed by SHEPMOD.
!
! Input parameters:
!  XP is the point at which the modified Shepard interpolant
!     function approximation function is to be evaluated.
!  M is the dimension of the data.
!  N is the number of interpolation points.
!  X contains the interpolation nodes, by column.
!  RW contains the radius of influence about each interpolation
!    node returned by SHEPMOD.
! Output parameter:
!  WTS contains the weight associated with each interpolation node.
!  IER
!   = 0, normal returns.
!   = 1, if the point XP is outside the radius of influence RW(i) for all
!        nodes, in which case SHEPMODVAL is computed using the original Shepard
!        algorithm with the M+1 closest points.'''
        
        # Setting up "xp"
        if ((not issubclass(type(xp), numpy.ndarray)) or
            (not numpy.asarray(xp).flags.f_contiguous) or
            (not (xp.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'xp' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            xp = numpy.asarray(xp, dtype=ctypes.c_double, order='F')
        xp_dim_1 = ctypes.c_int(xp.shape[0])
        
        # Setting up "m"
        if (type(m) is not ctypes.c_int): m = ctypes.c_int(m)
        
        # Setting up "n"
        if (type(n) is not ctypes.c_int): n = ctypes.c_int(n)
        
        # Setting up "x"
        if ((not issubclass(type(x), numpy.ndarray)) or
            (not numpy.asarray(x).flags.f_contiguous) or
            (not (x.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'x' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            x = numpy.asarray(x, dtype=ctypes.c_double, order='F')
        x_dim_1 = ctypes.c_int(x.shape[0])
        x_dim_2 = ctypes.c_int(x.shape[1])
        
        # Setting up "rw"
        if ((not issubclass(type(rw), numpy.ndarray)) or
            (not numpy.asarray(rw).flags.f_contiguous) or
            (not (rw.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'rw' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            rw = numpy.asarray(rw, dtype=ctypes.c_double, order='F')
        rw_dim_1 = ctypes.c_int(rw.shape[0])
        
        # Setting up "wts"
        if (wts is None):
            wts = numpy.zeros(shape=(n), dtype=ctypes.c_double, order='F')
        wts_dim_1 = ctypes.c_int(wts.shape[0])
        
        # Setting up "ier"
        ier = ctypes.c_int()
    
        # Call C-accessible Fortran wrapper.
        clib.c_shepmodval(ctypes.byref(xp_dim_1), ctypes.c_void_p(xp.ctypes.data), ctypes.byref(m), ctypes.byref(n), ctypes.byref(x_dim_1), ctypes.byref(x_dim_2), ctypes.c_void_p(x.ctypes.data), ctypes.byref(rw_dim_1), ctypes.c_void_p(rw.ctypes.data), ctypes.byref(wts_dim_1), ctypes.c_void_p(wts.ctypes.data), ctypes.byref(ier))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return wts, ier.value

modified_shepard = modified_shepard()

