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
shared_object_name = "fekete.so"
this_directory = os.path.dirname(os.path.abspath(__file__))
path_to_lib = os.path.join(this_directory, shared_object_name)
compile_options = "-fPIC -shared -O3 -lblas -llapack"
ordered_dependencies = ['fekete.f90', 'fekete_bind_c.f90']
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
# Wrapper for the Fortran subroutine FEKETE_INDICES

def fekete_indices(vm_t, inds=None):
    '''! This function identifies the set of row vectors in a Vandermonde
! matrix that should be kept. It outputs the row indices, however
! it expects the Vandermonde matrix to be provided in its transpose.
! Inspiration for greedy selection technique from:
!
!    Bos, Len, et al. "Computing multivariate Fekete and Leja points
!    by numerical linear algebra." SIAM Journal on Numerical Analysis
!    48.5 (2010): 1984-1999. DOI 10.1137/090779024
!
! USES LAPACK:
!    DGEQP3  --  (double precision QR decomposition with column pivoting)
!
! INPUT:
!    VM_T  --  The transposed Vandermonde matrix. Given functions
!              {f_1, ..., f_k} and points {p_1, ..., p_n}, we have
!                   VM_T_{i,j} = f_i(p_j)
!
! OUTPUT:
!    INDS  --  The sequence of indices (of the points {p_i}) that
!              are most useful in determining an interpolant. E.g.
!              to make a Vandermonde matrix square that was
!              previously long (more points than functions), use the
!              points at indices JPVT(1:K) as the Fekete points.
!'''
    
    # Setting up "vm_t"
    if ((not issubclass(type(vm_t), numpy.ndarray)) or
        (not numpy.asarray(vm_t).flags.f_contiguous) or
        (not (vm_t.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'vm_t' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        vm_t = numpy.asarray(vm_t, dtype=ctypes.c_double, order='F')
    vm_t_dim_1 = ctypes.c_int(vm_t.shape[0])
    vm_t_dim_2 = ctypes.c_int(vm_t.shape[1])
    
    # Setting up "inds"
    if (inds is None):
        inds = numpy.zeros(shape=(vm_t.shape[1]), dtype=ctypes.c_int, order='F')
    inds_dim_1 = ctypes.c_int(inds.shape[0])
    
    # Setting up "info"
    info = ctypes.c_int()

    # Call C-accessible Fortran wrapper.
    clib.c_fekete_indices(ctypes.byref(vm_t_dim_1), ctypes.byref(vm_t_dim_2), ctypes.c_void_p(vm_t.ctypes.data), ctypes.byref(inds_dim_1), ctypes.c_void_p(inds.ctypes.data), ctypes.byref(info))

    # Return final results, 'INTENT(OUT)' arguments only.
    return inds, info.value


# ----------------------------------------------
# Wrapper for the Fortran subroutine EVALUATE_VANDERMONDE

def evaluate_vandermonde(points, degrees, indices, vandermonde=None):
    '''! Evaluate polynomials at points to construct a Vandermonde matrix,
! the lengths are the number of terms in each polynomial, the indices
! are the active index terms in each polynomial.
!
! This should be super fast, aka even a matrix with 100M elements should
! take less than one second. For some reason it appears to slow down when
! (SIZE(POINTS, 1) == 1). That is the reason for including schedule(dynamic),
! which does not have a huge negative effect otherwise.'''
    
    # Setting up "points"
    if ((not issubclass(type(points), numpy.ndarray)) or
        (not numpy.asarray(points).flags.f_contiguous) or
        (not (points.dtype == numpy.dtype(ctypes.c_double)))):
        import warnings
        warnings.warn("The provided argument 'points' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
        points = numpy.asarray(points, dtype=ctypes.c_double, order='F')
    points_dim_1 = ctypes.c_int(points.shape[0])
    points_dim_2 = ctypes.c_int(points.shape[1])
    
    # Setting up "degrees"
    if ((not issubclass(type(degrees), numpy.ndarray)) or
        (not numpy.asarray(degrees).flags.f_contiguous) or
        (not (degrees.dtype == numpy.dtype(ctypes.c_long)))):
        import warnings
        warnings.warn("The provided argument 'degrees' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
        degrees = numpy.asarray(degrees, dtype=ctypes.c_long, order='F')
    degrees_dim_1 = ctypes.c_int(degrees.shape[0])
    
    # Setting up "indices"
    if ((not issubclass(type(indices), numpy.ndarray)) or
        (not numpy.asarray(indices).flags.f_contiguous) or
        (not (indices.dtype == numpy.dtype(ctypes.c_long)))):
        import warnings
        warnings.warn("The provided argument 'indices' was not an f_contiguous NumPy array of type 'ctypes.c_long' (or equivalent). Automatically converting (probably creating a full copy).")
        indices = numpy.asarray(indices, dtype=ctypes.c_long, order='F')
    indices_dim_1 = ctypes.c_int(indices.shape[0])
    indices_dim_2 = ctypes.c_int(indices.shape[1])
    
    # Setting up "vandermonde"
    if (vandermonde is None):
        vandermonde = numpy.zeros(shape=(degrees.size, points.shape[1]), dtype=ctypes.c_double, order='F')
    vandermonde_dim_1 = ctypes.c_int(vandermonde.shape[0])
    vandermonde_dim_2 = ctypes.c_int(vandermonde.shape[1])

    # Call C-accessible Fortran wrapper.
    clib.c_evaluate_vandermonde(ctypes.byref(points_dim_1), ctypes.byref(points_dim_2), ctypes.c_void_p(points.ctypes.data), ctypes.byref(degrees_dim_1), ctypes.c_void_p(degrees.ctypes.data), ctypes.byref(indices_dim_1), ctypes.byref(indices_dim_2), ctypes.c_void_p(indices.ctypes.data), ctypes.byref(vandermonde_dim_1), ctypes.byref(vandermonde_dim_2), ctypes.c_void_p(vandermonde.ctypes.data))

    # Return final results, 'INTENT(OUT)' arguments only.
    return vandermonde

