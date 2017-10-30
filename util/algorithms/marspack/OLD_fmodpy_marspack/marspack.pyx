
import cython
import numpy


cdef extern:
    void c_mars( long *n, long *p, double *x, double *y, double *w, long *nk, long *mi, long *lx, long *fm_0, double *fm, long *im_0, long *im )
#                 , long *sp_0, double *sp, long *dp_0, double *dp, long *mm_0, long *mm )

@cython.boundscheck(False)
@cython.wraparound(False)
def mars( long n, long p, double[:,:] x, double[:] y, double[:] w, long nk, long mi, long[:] lx, double[:] fm, long[:] im ):
    #          , double[:] sp, double[:] dp, long[:] mm ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (((type(x) != type(None)) and (not numpy.asarray(x).flags.f_contiguous))
        or ((type(y) != type(None)) and (not numpy.asarray(y).flags.f_contiguous))
        or ((type(w) != type(None)) and (not numpy.asarray(w).flags.f_contiguous))
        or ((type(lx) != type(None)) and (not numpy.asarray(lx).flags.f_contiguous))
        or ((type(fm) != type(None)) and (not numpy.asarray(fm).flags.f_contiguous))
        or ((type(im) != type(None)) and (not numpy.asarray(im).flags.f_contiguous))):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    # if (type(sp) != type(None)) and (not numpy.asarray(sp).flags.f_contiguous):
    #     raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    # if (type(dp) != type(None)) and (not numpy.asarray(dp).flags.f_contiguous):
    #     raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    # if (type(mm) != type(None)) and (not numpy.asarray(mm).flags.f_contiguous):
    #     raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    cdef long fm_0 = fm.shape[0]
    cdef long im_0 = im.shape[0]
    # cdef long sp_0 = sp.shape[0]
    # cdef long dp_0 = dp.shape[0]
    # cdef long mm_0 = mm.shape[0]
    
    # CODE FOR GENERATING FORTRANT TEST-PROGRAMS

    # # Extra variables for printing test fortran code
    # cdef long nmcv = 0
    # cdef long ntcv = 0
    # print("")
    # print("PROGRAM TEST_MARS")
    # print("  IMPLICIT NONE")
    # print("  INTEGER :: N, P, NK, MI")
    # print("  REAL,             DIMENSION(%i,%i) :: X"%(n,p))
    # print("  REAL,             DIMENSION(%i)    :: Y"%(n))
    # print("  REAL,             DIMENSION(%i)    :: W"%(n))
    # print("  INTEGER,          DIMENSION(%i)    :: LX"%(p))
    # print("  REAL,             DIMENSION(%i)    :: FM"%(fm_0))
    # print("  INTEGER,          DIMENSION(%i)    :: IM"%(im_0))
    # print("  REAL,             DIMENSION(%i)    :: SP"%(sp_0))
    # print("  DOUBLE PRECISION, DIMENSION(%i)    :: DP"%(dp_0))
    # print("  INTEGER,          DIMENSION(%i)    :: MM"%(mm_0))
    # print("  ")
    # print("  N = %i    ! Number of points"%(n))
    # print("  P = %i    ! Number of dimensions"%(p))
    # print("  NK = %i   ! Maximum number of basis functions"%(nk))    
    # print("  MI = %i   ! Maximum interaction between dimensions"%(mi))
    # print("  W = 1.0   ! Weights for each point")
    # print("  LX = 1    ! Flags for each dimension (1 = ordinal)")
    # print("  FM = 1.0  ! Holder for MARS model")
    # print("  IM = 1    ! Holder for MARS model")
    # print("  SP = 1.0  ! Real workspace array for MARS")
    # print("  DP = 1.0  ! Double precision workspace array for MARS")
    # print("  MM = 1    ! Integer workspace array for MARS")
    # print("  ")
    # print("  ! Size of fm = 3+nk*(5*mi+nmcv+6)+2*p+ntcv")
    # print("  !            = 3+%i*(5*%i+%i+6)+2*%i+%i)"%(nk,mi,nmcv,p,ntcv))
    # print("  !            = %i"%(3+nk*(5*mi+nmcv+6)+2*p+ntcv))
    # print("  ! Size of im = 21+nk*(3*mi+8)")
    # print("  !            = 21+%i*(3*%i+8)"%(nk,mi))
    # print("  !            = %i"%(21+nk*(3*mi+8)))
    # print("  ! Size of sp = n*(max(nk+1,2)+3)+max(3*n+5*nk+p,2*p,4*n)+2*p+4*nk")
    # print("  !            = n*(max(%i+1,2)+3)+max(3*%i+5*%i+%i,2*%i,4*%i)+2*%i+4*%i"%(nk,n,nk,p,p,n,p,nk))
    # print("  !            = %i"%(n*(max(nk+1,2)+3)+max(3*n+5*nk+p,2*p,4*n)+2*p+4*nk))
    # print("  ! Size of dp = max(n*nk,(nk+1)*(nk+1))+max((nk+2)*(nmcv+3),4*nk)")
    # print("  !            = max(%i*%i,(%i+1)*(%i+1))+max((%i+2)*(%i+3),4*%i)"%(n,nk,nk,nk,nk,nmcv,nk))
    # print("  !            = %i"%(max(n*nk,(nk+1)*(nk+1))+max((nk+2)*(nmcv+3),4*nk)))
    # print("  ! Size of mm = n*p+2*max(mi,nmcv)")
    # print("  !            = %i*%i+2*max(%i,%i)"%(n,p,mi,nmcv))
    # print("  !            = %i"%(n*p+2*max(mi,nmcv)))
    # print("  ")
    # for step in range(1,n+1):
    #     print("  X(%i,:) = (/ %s /)"%(step,", ".join(list(map(str,x[step-1,:])))))
    # print("  ")
    # for step in range(1,n+1):
    #     print("  Y(%i) = %s"%(step,y[step-1]))
    # print("  ")
    # print("  CALL MARS(N,P,X,Y,W,NK,MI,LX,FM,IM,SP,DP,MM)")
    # print("END PROGRAM TEST_MARS")
    # print("")

    c_mars(&n, &p, &x[0][0], &y[0], &w[0], &nk, &mi, &lx[0], &fm_0, &fm[0], &im_0, &im[0])
    #       , &sp_0, &sp[0], &dp_0, &dp[0], &mm_0, &mm[0])
    return numpy.asarray(fm, order='F'), numpy.asarray(im, order='F')


cdef extern:
    void c_fmod( long *m, long *n, long *x_1, double *x, long *fm_0, double *fm, long *im_0, long *im, long *f_0, double *f, long *sp_0, double *sp )

@cython.boundscheck(False)
@cython.wraparound(False)
def fmod( long m, long n, double[:,:] x, double[:] fm, long[:] im, double[:] f, double[:,:] sp ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (type(x) != type(None)) and (not numpy.asarray(x).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(fm) != type(None)) and (not numpy.asarray(fm).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(im) != type(None)) and (not numpy.asarray(im).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(f) != type(None)) and (not numpy.asarray(f).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(sp) != type(None)) and (not numpy.asarray(sp).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    cdef long x_1 = x.shape[1]
    cdef long fm_0 = fm.shape[0]
    cdef long im_0 = im.shape[0]
    cdef long f_0 = f.shape[0]
    cdef long sp_0 = sp.shape[0]
    
    c_fmod(&m, &n, &x_1, &x[0][0], &fm_0, &fm[0], &im_0, &im[0], &f_0, &f[0], &sp_0, &sp[0][0])
    return numpy.asarray(f, order='F')

