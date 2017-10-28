
import cython
import numpy


cdef extern:
    void c_lshep( long *m, long *n, double *x, double *f, double *a, double *rw, long *ier, double *rlshep )

@cython.boundscheck(False)
@cython.wraparound(False)
def lshep( long m, long n, double[:,:] x, double[:] f, o_a=None, o_rw=None, o_ier=None, o_rlshep=None ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (type(x) != type(None)) and (not numpy.asarray(x).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(f) != type(None)) and (not numpy.asarray(f).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_a) != type(None)) and (not numpy.asarray(o_a).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_rw) != type(None)) and (not numpy.asarray(o_rw).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if type(o_a) == type(None):
        o_a = (1) * numpy.ones(shape=(m,n),dtype=float,order='F')
    cdef double [:,:] a = o_a
    if type(o_rw) == type(None):
        o_rw = (1) * numpy.ones(shape=(n),dtype=float,order='F')
    cdef double [:] rw = o_rw
    if type(o_ier) == type(None):
        o_ier = (1)
    cdef long ier = o_ier
    if type(o_rlshep) == type(None):
        o_rlshep = (1)
    cdef double rlshep = o_rlshep
    
    c_lshep(&m, &n, &x[0][0], &f[0], &a[0][0], &rw[0], &ier, &rlshep)
    return numpy.asarray(a, order='F'), numpy.asarray(rw, order='F'), ier


cdef extern:
    void c_lshepval( double *xp, long *m, long *n, double *x, double *f, double *a, double *rw, long *ier, double *lshepval_output )

@cython.boundscheck(False)
@cython.wraparound(False)
def lshepval( double[:] xp, long m, long n, double[:,:] x, double[:] f, double[:,:] a, double[:] rw, o_ier=None, o_lshepval_output=None ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (type(xp) != type(None)) and (not numpy.asarray(xp).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(x) != type(None)) and (not numpy.asarray(x).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(f) != type(None)) and (not numpy.asarray(f).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(a) != type(None)) and (not numpy.asarray(a).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(rw) != type(None)) and (not numpy.asarray(rw).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if type(o_ier) == type(None):
        o_ier = (1)
    cdef long ier = o_ier
    if type(o_lshepval_output) == type(None):
        o_lshepval_output = (1)
    cdef double lshepval_output = o_lshepval_output
    
    c_lshepval(&xp[0], &m, &n, &x[0][0], &f[0], &a[0][0], &rw[0], &ier, &lshepval_output)
    return ier, lshepval_output

