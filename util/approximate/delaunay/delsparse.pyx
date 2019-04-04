'''! This module contains the REAL_PRECISION R8 data type for 64-bit arithmetic
! and interface blocks for the DELAUNAYSPARSES and DELAUNAYSPARSEP
! subroutines for computing the Delaunay simplices containing interpolation
! points Q in R^D given data points PTS.
! Interface for serial subroutine DELAUNAYSPARSES.'''


import cython
import numpy

class NotFortranCompatible(Exception): pass

#      Wrapper for fortran function delaunaysparses     
# =================================================

cdef extern:
    void c_delaunaysparses( int* d, int* n, int* pts_0, int* pts_1, double* pts, int* m, int* q_0, int* q_1, double* q, int* simps_0, int* simps_1, int* simps, int* weights_0, int* weights_1, double* weights, int* ierr_0, int* ierr, int* interp_in_0, int* interp_in_1, double* interp_in, bint* interp_in_present, int* interp_out_0, int* interp_out_1, double* interp_out, bint* interp_out_present, double* eps, bint* eps_present, double* extrap, bint* extrap_present, int* rnorm_0, double* rnorm, bint* rnorm_present, int* ibudget, bint* ibudget_present, bint* chain, bint* chain_present )

@cython.boundscheck(False)
@cython.wraparound(False)
def delaunaysparses( int d, int n, double[:,:] pts, int m, double[:,:] q, int[:,:] simps, double[:,:] weights, int[:] ierr, interp_in=None, interp_out=None, eps=None, extrap=None, rnorm=None, ibudget=None, chain=None ):
    ''''''
    # Prepare for fortran function call (initialize optionals)
    if (not numpy.asarray(pts).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int pts_0 = pts.shape[0]
    cdef int pts_1 = pts.shape[1]
    
    if (not numpy.asarray(q).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int q_0 = q.shape[0]
    cdef int q_1 = q.shape[1]
    
    if (not numpy.asarray(simps).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int simps_0 = simps.shape[0]
    cdef int simps_1 = simps.shape[1]
    
    if (not numpy.asarray(weights).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int weights_0 = weights.shape[0]
    cdef int weights_1 = weights.shape[1]
    
    cdef int ierr_0 = ierr.shape[0]
    
    cdef bint interp_in_present = True
    if (type(interp_in) == type(None)):
        interp_in_present = False
        interp_in = numpy.zeros(shape=(1,1),dtype=numpy.float64,order='F')
    cdef double[:,:] local_interp_in = interp_in
    if (not numpy.asarray(local_interp_in).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int interp_in_0 = interp_in.shape[0]
    cdef int interp_in_1 = interp_in.shape[1]
    
    cdef bint interp_out_present = True
    if (type(interp_out) == type(None)):
        interp_out_present = False
        interp_out = numpy.zeros(shape=(1,1),dtype=numpy.float64,order='F')
    cdef double[:,:] local_interp_out = interp_out
    if (not numpy.asarray(local_interp_out).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int interp_out_0 = interp_out.shape[0]
    cdef int interp_out_1 = interp_out.shape[1]
    
    cdef bint eps_present = True
    if (type(eps) == type(None)):
        eps_present = False
        eps = 1
    cdef double local_eps = eps
    
    cdef bint extrap_present = True
    if (type(extrap) == type(None)):
        extrap_present = False
        extrap = 1
    cdef double local_extrap = extrap
    
    cdef bint rnorm_present = True
    if (type(rnorm) == type(None)):
        rnorm_present = False
        rnorm = numpy.zeros(shape=(1),dtype=numpy.float64,order='F')
    cdef double[:] local_rnorm = rnorm
    cdef int rnorm_0 = rnorm.shape[0]
    
    cdef bint ibudget_present = True
    if (type(ibudget) == type(None)):
        ibudget_present = False
        ibudget = 1
    cdef int local_ibudget = ibudget
    
    cdef bint chain_present = True
    if (type(chain) == type(None)):
        chain_present = False
        chain = 1
    cdef bint local_chain = chain
    
    
    # Make fortran function call
    c_delaunaysparses(&d, &n, &pts_0, &pts_1, &pts[0][0], &m, &q_0, &q_1, &q[0][0], &simps_0, &simps_1, &simps[0][0], &weights_0, &weights_1, &weights[0][0], &ierr_0, &ierr[0], &interp_in_0, &interp_in_1, &local_interp_in[0][0], &interp_in_present, &interp_out_0, &interp_out_1, &local_interp_out[0][0], &interp_out_present, &local_eps, &eps_present, &local_extrap, &extrap_present, &rnorm_0, &local_rnorm[0], &rnorm_present, &local_ibudget, &ibudget_present, &local_chain, &chain_present)
    # Return appropriate values based on "INTENT" from fortran code
    
    return numpy.asarray(pts, order='F'), numpy.asarray(q, order='F'), numpy.asarray(simps, order='F'), numpy.asarray(weights, order='F'), numpy.asarray(ierr, order='F'), numpy.asarray(local_interp_out, order='F'), numpy.asarray(local_rnorm, order='F')



#      Wrapper for fortran function delaunaysparsep     
# =================================================

cdef extern:
    void c_delaunaysparsep( int* d, int* n, int* pts_0, int* pts_1, double* pts, int* m, int* q_0, int* q_1, double* q, int* simps_0, int* simps_1, int* simps, int* weights_0, int* weights_1, double* weights, int* ierr_0, int* ierr, int* interp_in_0, int* interp_in_1, double* interp_in, bint* interp_in_present, int* interp_out_0, int* interp_out_1, double* interp_out, bint* interp_out_present, double* eps, bint* eps_present, double* extrap, bint* extrap_present, int* rnorm_0, double* rnorm, bint* rnorm_present, int* ibudget, bint* ibudget_present, bint* chain, bint* chain_present, int* pmode, bint* pmode_present )

@cython.boundscheck(False)
@cython.wraparound(False)
def delaunaysparsep( int d, int n, double[:,:] pts, int m, double[:,:] q, int[:,:] simps, double[:,:] weights, int[:] ierr, interp_in=None, interp_out=None, eps=None, extrap=None, rnorm=None, ibudget=None, chain=None, pmode=None ):
    ''''''
    # Prepare for fortran function call (initialize optionals)
    if (not numpy.asarray(pts).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int pts_0 = pts.shape[0]
    cdef int pts_1 = pts.shape[1]
    
    if (not numpy.asarray(q).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int q_0 = q.shape[0]
    cdef int q_1 = q.shape[1]
    
    if (not numpy.asarray(simps).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int simps_0 = simps.shape[0]
    cdef int simps_1 = simps.shape[1]
    
    if (not numpy.asarray(weights).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int weights_0 = weights.shape[0]
    cdef int weights_1 = weights.shape[1]
    
    cdef int ierr_0 = ierr.shape[0]
    
    cdef bint interp_in_present = True
    if (type(interp_in) == type(None)):
        interp_in_present = False
        interp_in = numpy.zeros(shape=(1,1),dtype=numpy.float64,order='F')
    cdef double[:,:] local_interp_in = interp_in
    if (not numpy.asarray(local_interp_in).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int interp_in_0 = interp_in.shape[0]
    cdef int interp_in_1 = interp_in.shape[1]
    
    cdef bint interp_out_present = True
    if (type(interp_out) == type(None)):
        interp_out_present = False
        interp_out = numpy.zeros(shape=(1,1),dtype=numpy.float64,order='F')
    cdef double[:,:] local_interp_out = interp_out
    if (not numpy.asarray(local_interp_out).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int interp_out_0 = interp_out.shape[0]
    cdef int interp_out_1 = interp_out.shape[1]
    
    cdef bint eps_present = True
    if (type(eps) == type(None)):
        eps_present = False
        eps = 1
    cdef double local_eps = eps
    
    cdef bint extrap_present = True
    if (type(extrap) == type(None)):
        extrap_present = False
        extrap = 1
    cdef double local_extrap = extrap
    
    cdef bint rnorm_present = True
    if (type(rnorm) == type(None)):
        rnorm_present = False
        rnorm = numpy.zeros(shape=(1),dtype=numpy.float64,order='F')
    cdef double[:] local_rnorm = rnorm
    cdef int rnorm_0 = rnorm.shape[0]
    
    cdef bint ibudget_present = True
    if (type(ibudget) == type(None)):
        ibudget_present = False
        ibudget = 1
    cdef int local_ibudget = ibudget
    
    cdef bint chain_present = True
    if (type(chain) == type(None)):
        chain_present = False
        chain = 1
    cdef bint local_chain = chain
    
    cdef bint pmode_present = True
    if (type(pmode) == type(None)):
        pmode_present = False
        pmode = 1
    cdef int local_pmode = pmode
    
    
    # Make fortran function call
    c_delaunaysparsep(&d, &n, &pts_0, &pts_1, &pts[0][0], &m, &q_0, &q_1, &q[0][0], &simps_0, &simps_1, &simps[0][0], &weights_0, &weights_1, &weights[0][0], &ierr_0, &ierr[0], &interp_in_0, &interp_in_1, &local_interp_in[0][0], &interp_in_present, &interp_out_0, &interp_out_1, &local_interp_out[0][0], &interp_out_present, &local_eps, &eps_present, &local_extrap, &extrap_present, &rnorm_0, &local_rnorm[0], &rnorm_present, &local_ibudget, &ibudget_present, &local_chain, &chain_present, &local_pmode, &pmode_present)
    # Return appropriate values based on "INTENT" from fortran code
    
    return numpy.asarray(pts, order='F'), numpy.asarray(q, order='F'), numpy.asarray(simps, order='F'), numpy.asarray(weights, order='F'), numpy.asarray(ierr, order='F'), numpy.asarray(local_interp_out, order='F'), numpy.asarray(local_rnorm, order='F')



