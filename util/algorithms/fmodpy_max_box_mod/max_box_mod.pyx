
import cython
import numpy


cdef extern:
    void c_max_boxes( long *pts_0, long *pts_1, double *pts, double *box_widths )

@cython.boundscheck(False)
@cython.wraparound(False)
def max_boxes( double[:,:] pts, o_box_widths=None ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (type(pts) != type(None)) and (not numpy.asarray(pts).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_box_widths) != type(None)) and (not numpy.asarray(o_box_widths).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    cdef long pts_0 = pts.shape[0]
    cdef long pts_1 = pts.shape[1]
    if type(o_box_widths) == type(None):
        o_box_widths = (1) * numpy.ones(shape=(2*pts.shape[1],pts.shape[2]),dtype=float,order='F')
    cdef double [:,:] box_widths = o_box_widths
    
    c_max_boxes(&pts_0, &pts_1, &pts[0][0], &box_widths[0][0])
    return numpy.asarray(box_widths, order='F')


cdef extern:
    void c_linear_eval( long *boxes_0, long *boxes_1, double *boxes, double *box_widths, double *box_values, long *pts_0, long *pts_1, double *pts, double *values, long *error )

@cython.boundscheck(False)
@cython.wraparound(False)
def linear_eval( double[:,:] boxes, double[:,:] box_widths, double[:] box_values, double[:,:] pts, o_values=None, o_error=None ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (type(boxes) != type(None)) and (not numpy.asarray(boxes).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(box_widths) != type(None)) and (not numpy.asarray(box_widths).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(box_values) != type(None)) and (not numpy.asarray(box_values).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(pts) != type(None)) and (not numpy.asarray(pts).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_values) != type(None)) and (not numpy.asarray(o_values).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    cdef long boxes_0 = boxes.shape[0]
    cdef long boxes_1 = boxes.shape[1]
    cdef long pts_0 = pts.shape[0]
    cdef long pts_1 = pts.shape[1]
    if type(o_values) == type(None):
        o_values = (1) * numpy.ones(shape=(pts.shape[2]),dtype=float,order='F')
    cdef double [:] values = o_values
    if type(o_error) == type(None):
        o_error = (1)
    cdef long error = o_error
    
    c_linear_eval(&boxes_0, &boxes_1, &boxes[0][0], &box_widths[0][0], &box_values[0], &pts_0, &pts_1, &pts[0][0], &values[0], &error)
    return numpy.asarray(values, order='F'), error


cdef extern:
    void c_qsortc( long *a_0, double *a, long *idx )

@cython.boundscheck(False)
@cython.wraparound(False)
def qsortc( double[:] a, o_idx=None ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (type(a) != type(None)) and (not numpy.asarray(a).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_idx) != type(None)) and (not numpy.asarray(o_idx).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    cdef long a_0 = a.shape[0]
    if type(o_idx) == type(None):
        o_idx = (1) * numpy.ones(shape=(a.size),dtype=int,order='F')
    cdef long [:] idx = o_idx
    
    c_qsortc(&a_0, &a[0], &idx[0])
    return numpy.asarray(a, order='F'), numpy.asarray(idx, order='F')

