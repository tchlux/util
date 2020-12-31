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
shared_object_name = "linear_shepard.so"
this_directory = os.path.dirname(os.path.abspath(__file__))
path_to_lib = os.path.join(this_directory, shared_object_name)
compile_options = "-fPIC -shared -O3 -lblas -llapack"
ordered_dependencies = ['linear_shepard.f95', 'linear_shepard_bind_c.f90']
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


class real_precision:
    ''''''

    # Declare 'r8'
    def get_r8(self):
        r8 = ctypes.c_int()
        clib.real_precision_get_r8(ctypes.byref(r8))
        return r8.value
    def set_r8(self, r8):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    r8 = property(get_r8, set_r8)

real_precision = real_precision()


class linear_shepard_mod:
    '''! INTERFACE
!    SUBROUTINE DGELSS(M, N, NRHS, A, LDA, B, LDB, SIGMA, RCOND, RANK, &
!         WORK, LWORK, INFO)
!      USE REAL_PRECISION
!      INTEGER, INTENT(IN) :: M, N, NRHS, LDA, LDB, LWORK
!      INTEGER, INTENT(OUT) :: RANK, INFO
!      REAL(KIND = R8), INTENT(IN) :: RCOND
!      REAL(KIND = R8), DIMENSION(LDA, *), INTENT(INOUT) :: A
!      REAL(KIND = R8), DIMENSION(LDB, *), INTENT(INOUT) :: B
!      REAL(KIND = R8), DIMENSION(N), INTENT(OUT) :: SIGMA
!      REAL(KIND = R8), DIMENSION(LWORK), INTENT(OUT) :: WORK
!    END SUBROUTINE DGELSS
! END INTERFACE'''

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine LSHEP
    
    def lshep(self, m, n, x, f, a=None, rw=None, rlshep=None):
        '''! This subroutine computes a set of parameters defining a function that
! interpolates N data values F(i) at scattered nodes X(i) in M dimensions. The
! interpolant may be evaluated at an arbitrary point by the function LSHEPVAL.
! The interpolation scheme is a modified linear Shepard method, and
! can use either the Shepard weighted least squares fit or robust
! M-estimation for each local approximation.
!
! Input parameters:
!   M is the dimension of the data.
!   N is the number of nodes.
!   X(M, N) contains the coordinates of the nodes.
!   F(N) contains the function values of the nodes.
! Output parameters:
!   A(M, N) contains the coefficients of the local fits.
!   RW(N) is an array containing the radius of influence for each point.
!   IER
!   = 0, if no errors were encountered.
!   = 1, if N is too small relative to M.
!   = 2, if any least squares problem is rank deficient, in which
!        case an SVD based minimum norm solution is returned; the
!        local fit should still be reasonable.
!   = 3, if the IRLS subroutine returns an error, meaning the IRLS
!        iteration has failed, and the returned local fit may
!        be bad.
! Optional input:
!   RLSHEP specifies by its presence that robust M-estimation is to be used to
!          compute the local approximation.'''
        
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
        
        # Setting up "f"
        if ((not issubclass(type(f), numpy.ndarray)) or
            (not numpy.asarray(f).flags.f_contiguous) or
            (not (f.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'f' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            f = numpy.asarray(f, dtype=ctypes.c_double, order='F')
        f_dim_1 = ctypes.c_int(f.shape[0])
        
        # Setting up "a"
        if (a is None):
            a = numpy.zeros(shape=(m, n), dtype=ctypes.c_double, order='F')
        a_dim_1 = ctypes.c_int(a.shape[0])
        a_dim_2 = ctypes.c_int(a.shape[1])
        
        # Setting up "rw"
        if (rw is None):
            rw = numpy.zeros(shape=(n), dtype=ctypes.c_double, order='F')
        rw_dim_1 = ctypes.c_int(rw.shape[0])
        
        # Setting up "ier"
        ier = ctypes.c_int()
        
        # Setting up "rlshep"
        rlshep_present = ctypes.c_bool(True)
        if (rlshep is None):
            rlshep_present = ctypes.c_bool(False)
            rlshep = ctypes.c_bool()
        if (type(rlshep) is not ctypes.c_bool): rlshep = ctypes.c_bool(rlshep)
    
        # Call C-accessible Fortran wrapper.
        clib.c_lshep(ctypes.byref(m), ctypes.byref(n), ctypes.byref(x_dim_1), ctypes.byref(x_dim_2), ctypes.c_void_p(x.ctypes.data), ctypes.byref(f_dim_1), ctypes.c_void_p(f.ctypes.data), ctypes.byref(a_dim_1), ctypes.byref(a_dim_2), ctypes.c_void_p(a.ctypes.data), ctypes.byref(rw_dim_1), ctypes.c_void_p(rw.ctypes.data), ctypes.byref(ier), ctypes.byref(rlshep_present), ctypes.byref(rlshep))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return a, rw, ier.value, (rlshep.value if rlshep_present else None)

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine RIPPLE
    
    def ripple(self, m, n, x, f, a=None, rw=None):
        '''! The subroutine RIPPLE, residual initiated polynomial-time piecewise linear
! estimation, computes a set of parameters defining a function that
! interpolates N data values F(i) at scattered nodes X(i) in M dimensions.
! The interpolant may be evaluated at an arbitrary point by the function
! LSHEPVAL.
! The interpolation scheme is a modified linear Shepard method.
! The time complexity of the algorithm RIPPLE is between the robust
! M-estimation and LMS estimation.
!
! Input parameters:
!   M is the dimension of the data.
!   N is the number of nodes.
!   X(M, N) contains the coordinates of the nodes.
!   F(N) contains the function values of the nodes.
! Output parameters:
!   A(M, N) contains the coefficients of the local fits.
!   RW(N) is an array containing the radius of influence for each point.
!   IER
!   = 0, if no errors were encountered.
!   = 1, if N is too small relative to M.
!   = 2, if the IRLS subroutine returns an error.'''
        
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
        
        # Setting up "f"
        if ((not issubclass(type(f), numpy.ndarray)) or
            (not numpy.asarray(f).flags.f_contiguous) or
            (not (f.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'f' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            f = numpy.asarray(f, dtype=ctypes.c_double, order='F')
        f_dim_1 = ctypes.c_int(f.shape[0])
        
        # Setting up "a"
        if (a is None):
            a = numpy.zeros(shape=(m, n), dtype=ctypes.c_double, order='F')
        a_dim_1 = ctypes.c_int(a.shape[0])
        a_dim_2 = ctypes.c_int(a.shape[1])
        
        # Setting up "rw"
        if (rw is None):
            rw = numpy.zeros(shape=(n), dtype=ctypes.c_double, order='F')
        rw_dim_1 = ctypes.c_int(rw.shape[0])
        
        # Setting up "ier"
        ier = ctypes.c_int()
    
        # Call C-accessible Fortran wrapper.
        clib.c_ripple(ctypes.byref(m), ctypes.byref(n), ctypes.byref(x_dim_1), ctypes.byref(x_dim_2), ctypes.c_void_p(x.ctypes.data), ctypes.byref(f_dim_1), ctypes.c_void_p(f.ctypes.data), ctypes.byref(a_dim_1), ctypes.byref(a_dim_2), ctypes.c_void_p(a.ctypes.data), ctypes.byref(rw_dim_1), ctypes.c_void_p(rw.ctypes.data), ctypes.byref(ier))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return a, rw, ier.value

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine LSHEPVAL
    
    def lshepval(self, xp, m, n, x, f, a, rw):
        '''! LSHEPVAL returns the linear Shepard approximation at the point XP, using the
! local linear approximations computed by LSHEP.
!
! Input parameters:
!  XP is the point at which the linear Shepard interpolant function
!     approximation function is to be evaluated.
!  M is the dimension of the data.
!  N is the number of interpolation points.
!  X contains the interpolation nodes, by column.
!  F contains the function values of the interpolation nodes.
!  A is an M by N matrix containing the coefficients for linear nodal functions
!    returned by LSHEP.
!  RW contains the radius of influence about each interpolation node returned
!    by LSHEP.
! Output parameter:
!  IER
!   = 0, normal returns.
!   = 1, if the point XP is outside the radius of influence RW(i) for all
!        nodes, in which case LSHEPVAL is computed using the original Shepard
!        algorithm with the M+1 closest points.
!   = 2, if the hyperplane normals of the local approximations with positive
!        weights are significantly different. For a nonlinear underlying
!        function f(x), e.g., quadratic f(x), very different normals are
!        typical. For a piecewise linear underlying function f(x), IER = 2
!        signals a potentially large error in LSHEPVAL, since local
!        approximations from different facets of f(x) have been used.'''
        
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
        
        # Setting up "f"
        if ((not issubclass(type(f), numpy.ndarray)) or
            (not numpy.asarray(f).flags.f_contiguous) or
            (not (f.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'f' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            f = numpy.asarray(f, dtype=ctypes.c_double, order='F')
        f_dim_1 = ctypes.c_int(f.shape[0])
        
        # Setting up "a"
        if ((not issubclass(type(a), numpy.ndarray)) or
            (not numpy.asarray(a).flags.f_contiguous) or
            (not (a.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'a' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            a = numpy.asarray(a, dtype=ctypes.c_double, order='F')
        a_dim_1 = ctypes.c_int(a.shape[0])
        a_dim_2 = ctypes.c_int(a.shape[1])
        
        # Setting up "rw"
        if ((not issubclass(type(rw), numpy.ndarray)) or
            (not numpy.asarray(rw).flags.f_contiguous) or
            (not (rw.dtype == numpy.dtype(ctypes.c_double)))):
            import warnings
            warnings.warn("The provided argument 'rw' was not an f_contiguous NumPy array of type 'ctypes.c_double' (or equivalent). Automatically converting (probably creating a full copy).")
            rw = numpy.asarray(rw, dtype=ctypes.c_double, order='F')
        rw_dim_1 = ctypes.c_int(rw.shape[0])
        
        # Setting up "ier"
        ier = ctypes.c_int()
        
        # Setting up "lshepval_result"
        lshepval_result = ctypes.c_double()
    
        # Call C-accessible Fortran wrapper.
        clib.c_lshepval(ctypes.byref(xp_dim_1), ctypes.c_void_p(xp.ctypes.data), ctypes.byref(m), ctypes.byref(n), ctypes.byref(x_dim_1), ctypes.byref(x_dim_2), ctypes.c_void_p(x.ctypes.data), ctypes.byref(f_dim_1), ctypes.c_void_p(f.ctypes.data), ctypes.byref(a_dim_1), ctypes.byref(a_dim_2), ctypes.c_void_p(a.ctypes.data), ctypes.byref(rw_dim_1), ctypes.c_void_p(rw.ctypes.data), ctypes.byref(ier), ctypes.byref(lshepval_result))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return ier.value, lshepval_result.value

linear_shepard_mod = linear_shepard_mod()

