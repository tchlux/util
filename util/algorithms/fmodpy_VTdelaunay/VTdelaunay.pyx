'''!
! Module containing DelaunayP() subroutine for computing a delaunay
! simplex containing a given point and R8 data type for ~64 bit
! arithmetic on most known machines.
!
! Written to FORTRAN 2003 standard.
! Requires BLAS and LAPACK
!
! Tested w/ gfortran 4.8.5 on x86_64-redhat-linux
! AUTHOR: Tyler Chang
! Last Update: October, 2017
!'''


import cython
import numpy

class NotFortranCompatible(Exception): pass

#      Wrapper for fortran function delaunayp     
# =================================================

cdef extern:
    void c_delaunayp( int* d, int* pts_0, int* pts_1, double* pts, int* p_0, int* p_1, double* p, int* work_0, double* work, int* simp_0, int* simp_1, int* simp, int* weights_0, int* weights_1, double* weights, int* err_0, int* err, int* interp_in_opt_0, int* interp_in_opt_1, double* interp_in_opt, bint* interp_in_opt_present, int* interp_out_opt_0, int* interp_out_opt_1, double* interp_out_opt, bint* interp_out_opt_present, double* eps_opt, bint* eps_opt_present, int* rnorm_opt_0, double* rnorm_opt, bint* rnorm_opt_present, int* budget_opt, bint* budget_opt_present, double* extrap_opt, bint* extrap_opt_present )

@cython.boundscheck(False)
@cython.wraparound(False)
def delaunayp( int d, double[:,:] pts, double[:,:] p, double[:] work, int[:,:] simp, double[:,:] weights, int[:] err, interp_in_opt=None, interp_out_opt=None, eps_opt=None, rnorm_opt=None, budget_opt=None, extrap_opt=None ):
    '''!
    ! Compute a simplex in the delaunay triangulation containing the point(s)
    ! in p. Compute a delaunay simplex "close to" each point in p by solving
    ! a series of LS problems. Proceed to connect faces of the triangulation
    ! until a simplex is found containing each point in p.
    !
    ! INTEGER (IN) d : Dimension of space.
    !
    ! DOUBLE (IN) pts(d,n) : Set of n points in R^d to triangulate.
    !       Note: n is assumed based on SIZE(pts).
    !
    ! DOUBLE (IN) p(d,m) : Set of m points in R^d to interpolate using the
    !       delaunay triangulation.
    !       Note: m is assumed based on SIZE(p)
    !
    ! DOUBLE (INOUT) work(lwork) : Work array. The work array must have size
    !       of AT LEAST 5*d with optimal performance at size m*d.
    !
    ! INTEGER (OUT) simp(d+1) : The indices of the points which make up the
    !       vertices of the simplex containing p.
    !
    ! DOUBLE (OUT) weights(d+1) : The weights for the vertices of simp for
    !       interpolating at p.
    !
    ! INTEGER (OUT) err(m) : Vector of error flags corresponding to each
    !       point in p:
    !       SUCCESSFUL INTERPOLATION = 0
    !       SUCCESSFUL EXTRAPOLATION = 1
    !       OVER EXTRAPOLATION, NO SOLN COMPUTED = 2
    !       See other codes in ErrorCodes.txt file.
    !       Note: m must match SIZE(p,2)
    !
    ! -------------------------OPTIONAL INPUTS/OUTPUTS-------------------------
    !
    ! OPTIONAL, DOUBLE (IN) interp_in_opt(z,n) : A vector of (z) response values
    !       for each of the n triangulation points in pts. Note: z is assumed
    !       based on SIZE(interp_in_opt).
    !
    ! OPTIONAL, DOUBLE (out) interp_out_opt(z,m) : A vector of (z) interpolated
    !       values for each of the m interpolation points in p. Note: z and m
    !       must match the dimensions of interp_in_opt and p respectively.
    !
    ! OPTIONAL, DOUBLE (IN) eps_opt : A small number defining the tolerance for
    !       accepting "nearly" containing simplices. 0 < eps <<< scale of data.
    !       Default value: eps = DIAMETER(PTS) * SQRT(EPSILON(KIND=R8)).
    ! OPTIONAL, DOUBLE (IN) extrap_opt : If extrap_opt = 0 (to working
    !       precision), then no extrapolation is performed. Otherwise,
    !       extrapolation will be done when p is outside the convex hull but
    !       within extrap_opt * diameter(pts) distance of the convex hull.
    !       By default, extrap_opt = 0.1 (10% of diameter extrapolation).
    ! OPTIONAL, DOUBLE (OUT) rnorm_opt(m) : The 2-norm of the residual vector
    !       when projecting p onto the convex hull of pts. This value is always
    !       calculated when extrapolation is turned on, but only returned when
    !       this parameter is present.
    ! OPTIONAL, INTEGER (IN) budget_opt : The maximum number of simplices that
    !       the user is willing to construct when searching for the simplex
    !       containing p. In practice, solutions should be found in no more than
    !       (d)log(n) iterations, with a typical case of 1~10. The default value
    !       for budget is: budget = 1000. If failure to converge in this many
    !       iterations, it is most likely the case that eps is too small for the
    !       scale of this problem.
    !
    ! inputs'''
    # Prepare for fortran function call (initialize optionals)
    if (not numpy.asarray(pts).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int pts_0 = pts.shape[0]
    cdef int pts_1 = pts.shape[1]
    
    if (not numpy.asarray(p).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int p_0 = p.shape[0]
    cdef int p_1 = p.shape[1]
    
    if (not numpy.asarray(work).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int work_0 = work.shape[0]
    
    if (not numpy.asarray(simp).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int simp_0 = simp.shape[0]
    cdef int simp_1 = simp.shape[1]
    
    if (not numpy.asarray(weights).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int weights_0 = weights.shape[0]
    cdef int weights_1 = weights.shape[1]
    
    if (not numpy.asarray(err).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int err_0 = err.shape[0]
    
    cdef bint interp_in_opt_present = True
    if (type(interp_in_opt) == type(None)):
        interp_in_opt_present = False
        interp_in_opt = numpy.ones(shape=(1,1),dtype=numpy.float64,order='F')
    cdef double[:,:] local_interp_in_opt = interp_in_opt
    if (not numpy.asarray(local_interp_in_opt).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int interp_in_opt_0 = interp_in_opt.shape[0]
    cdef int interp_in_opt_1 = interp_in_opt.shape[1]
    
    cdef bint interp_out_opt_present = True
    if (type(interp_out_opt) == type(None)):
        interp_out_opt_present = False
        interp_out_opt = numpy.ones(shape=(1,1),dtype=numpy.float64,order='F')
    cdef double[:,:] local_interp_out_opt = interp_out_opt
    if (not numpy.asarray(local_interp_out_opt).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int interp_out_opt_0 = interp_out_opt.shape[0]
    cdef int interp_out_opt_1 = interp_out_opt.shape[1]
    
    cdef bint eps_opt_present = True
    if (type(eps_opt) == type(None)):
        eps_opt_present = False
        eps_opt = 1
    cdef double local_eps_opt = eps_opt
    
    cdef bint rnorm_opt_present = True
    if (type(rnorm_opt) == type(None)):
        rnorm_opt_present = False
        rnorm_opt = numpy.ones(shape=(1),dtype=numpy.float64,order='F')
    cdef double[:] local_rnorm_opt = rnorm_opt
    if (not numpy.asarray(local_rnorm_opt).flags.f_contiguous):
        raise(NotFortranCompatible('Only use numpy arrays that are f_contiguous.'))
    cdef int rnorm_opt_0 = rnorm_opt.shape[0]
    
    cdef bint budget_opt_present = True
    if (type(budget_opt) == type(None)):
        budget_opt_present = False
        budget_opt = 1
    cdef int local_budget_opt = budget_opt
    
    cdef bint extrap_opt_present = True
    if (type(extrap_opt) == type(None)):
        extrap_opt_present = False
        extrap_opt = 1
    cdef double local_extrap_opt = extrap_opt
    
    
    # Make fortran function call
    c_delaunayp(&d, &pts_0, &pts_1, &pts[0][0], &p_0, &p_1, &p[0][0], &work_0, &work[0], &simp_0, &simp_1, &simp[0][0], &weights_0, &weights_1, &weights[0][0], &err_0, &err[0], &interp_in_opt_0, &interp_in_opt_1, &local_interp_in_opt[0][0], &interp_in_opt_present, &interp_out_opt_0, &interp_out_opt_1, &local_interp_out_opt[0][0], &interp_out_opt_present, &local_eps_opt, &eps_opt_present, &rnorm_opt_0, &local_rnorm_opt[0], &rnorm_opt_present, &local_budget_opt, &budget_opt_present, &local_extrap_opt, &extrap_opt_present)
    # Return appropriate values based on "INTENT" from fortran code
    
    return numpy.asarray(work, order='F'), numpy.asarray(simp, order='F'), numpy.asarray(weights, order='F'), numpy.asarray(err, order='F'), numpy.asarray(local_interp_out_opt, order='F'), numpy.asarray(local_rnorm_opt, order='F')



