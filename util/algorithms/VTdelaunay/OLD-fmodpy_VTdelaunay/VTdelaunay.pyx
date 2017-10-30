
import cython
import numpy

cdef extern:
    void c_delaunayp( long* d, long* n, double* pts, double* p, long* simp, double* weights, long* err, double* eps_opt, bint* eps_opt_present, double* rnorm_opt, bint* rnorm_opt_present, long* budget_opt, bint* budget_opt_present, bint* extrap_opt, bint* extrap_opt_present )

@cython.boundscheck(False)
@cython.wraparound(False)
def delaunayp( long d, long n, double[:,:] pts, double[:] p, o_simp=None, o_weights=None, o_err=None, o_eps_opt=None, o_rnorm_opt=None, o_budget_opt=None, o_extrap_opt=None ):
    '''Created automatically by fmodpy. Documentation coming soon.'''

    if (not numpy.asarray(pts).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    
    if (not numpy.asarray(p).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    
    if (type(o_simp) == type(None)):
        o_simp = numpy.ones(shape=(d+1),dtype=int,order='F')
    cdef long[:] simp = o_simp
    if (not numpy.asarray(o_simp).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    
    if (type(o_weights) == type(None)):
        o_weights = numpy.ones(shape=(d+1),dtype=float,order='F')
    cdef double[:] weights = o_weights
    if (not numpy.asarray(o_weights).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    
    if (type(o_err) == type(None)):
        o_err = 1
    cdef long err = o_err
    
    cdef bint eps_opt_present = True
    if (type(o_eps_opt) == type(None)):
        eps_opt_present = False
        o_eps_opt = 1
    cdef double eps_opt = o_eps_opt
    
    cdef bint rnorm_opt_present = True
    if (type(o_rnorm_opt) == type(None)):
        rnorm_opt_present = False
        o_rnorm_opt = 1
    cdef double rnorm_opt = o_rnorm_opt
    
    cdef bint budget_opt_present = True
    if (type(o_budget_opt) == type(None)):
        budget_opt_present = False
        o_budget_opt = 1
    cdef long budget_opt = o_budget_opt
    
    cdef bint extrap_opt_present = True
    if (type(o_extrap_opt) == type(None)):
        extrap_opt_present = False
        o_extrap_opt = 1
    cdef bint extrap_opt = o_extrap_opt
    
    
    c_delaunayp(&d, &n, &pts[0][0], &p[0], &simp[0], &weights[0], &err, &eps_opt, &eps_opt_present, &rnorm_opt, &rnorm_opt_present, &budget_opt, &budget_opt_present, &extrap_opt, &extrap_opt_present)
    
    return numpy.asarray(simp, order='F'), numpy.asarray(weights, order='F'), err, rnorm_opt


