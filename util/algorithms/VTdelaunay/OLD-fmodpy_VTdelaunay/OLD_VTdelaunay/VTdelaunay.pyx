
import cython
import numpy


cdef extern:
    void c_delaunayp( long *d, long *n, double *pts, double *p, double *eps, long *simp, double *weights, long *err )

@cython.boundscheck(False)
@cython.wraparound(False)
def delaunayp( long d, long n, double[:,:] pts, double[:] p, double eps, o_simp=None, o_weights=None, o_err=None ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (type(pts) != type(None)) and (not numpy.asarray(pts).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(p) != type(None)) and (not numpy.asarray(p).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_simp) != type(None)) and (not numpy.asarray(o_simp).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_weights) != type(None)) and (not numpy.asarray(o_weights).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if type(o_simp) == type(None):
        o_simp = (1) * numpy.ones(shape=(d+1),dtype=int,order='F')
    cdef long [:] simp = o_simp
    if type(o_weights) == type(None):
        o_weights = (1) * numpy.ones(shape=(d+1),dtype=float,order='F')
    cdef double [:] weights = o_weights
    if type(o_err) == type(None):
        o_err = (1)
    cdef long err = o_err
    
    c_delaunayp(&d, &n, &pts[0][0], &p[0], &eps, &simp[0], &weights[0], &err)
    return numpy.asarray(simp, order='F'), numpy.asarray(weights, order='F'), err

