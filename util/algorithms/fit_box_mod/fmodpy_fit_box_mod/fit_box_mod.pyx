
import cython
import numpy


cdef extern:
    void c_compute_boxes( long *x_data_0, long *x_data_1, double *x_data, double *response, long *boxes_0, long *boxes_1, double *boxes, double *widths, double *weights )

@cython.boundscheck(False)
@cython.wraparound(False)
def compute_boxes( double[:,:] x_data, double[:] response, double[:,:] boxes, o_widths=None, o_weights=None ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (type(x_data) != type(None)) and (not numpy.asarray(x_data).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(response) != type(None)) and (not numpy.asarray(response).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(boxes) != type(None)) and (not numpy.asarray(boxes).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_widths) != type(None)) and (not numpy.asarray(o_widths).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_weights) != type(None)) and (not numpy.asarray(o_weights).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    cdef long x_data_0 = x_data.shape[0]
    cdef long x_data_1 = x_data.shape[1]
    cdef long boxes_0 = boxes.shape[0]
    cdef long boxes_1 = boxes.shape[1]

    if type(o_widths) == type(None):
        o_widths = (1) * numpy.ones(shape=(boxes.shape[2]),dtype=float,order='F')
    cdef double [:] widths = o_widths

    if type(o_weights) == type(None):
        o_weights = (1) * numpy.ones(shape=(boxes.shape[2]),dtype=float,order='F')
    cdef double [:] weights = o_weights
    
    c_compute_boxes(&x_data_0, &x_data_1, &x_data[0][0], &response[0], &boxes_0, &boxes_1, &boxes[0][0], &widths[0], &weights[0])
    return numpy.asarray(boxes, order='F'), numpy.asarray(widths, order='F'), numpy.asarray(weights, order='F')


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
    void c_eval_box_coefs( long *boxes_0, long *boxes_1, double *boxes, long *widths_0, double *widths, long *x_points_0, long *x_points_1, double *x_points, double *box_vals )

@cython.boundscheck(False)
@cython.wraparound(False)
def eval_box_coefs( double[:,:] boxes, double[:] widths, double[:,:] x_points, o_box_vals=None ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (type(boxes) != type(None)) and (not numpy.asarray(boxes).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(widths) != type(None)) and (not numpy.asarray(widths).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(x_points) != type(None)) and (not numpy.asarray(x_points).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_box_vals) != type(None)) and (not numpy.asarray(o_box_vals).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    cdef long boxes_0 = boxes.shape[0]
    cdef long boxes_1 = boxes.shape[1]
    cdef long widths_0 = widths.shape[0]
    cdef long x_points_0 = x_points.shape[0]
    cdef long x_points_1 = x_points.shape[1]
    if type(o_box_vals) == type(None):
        o_box_vals = (1) * numpy.ones(shape=(boxes.shape[2],x_points.shape[2]),dtype=float,order='F')
    cdef double [:,:] box_vals = o_box_vals
    
    c_eval_box_coefs(&boxes_0, &boxes_1, &boxes[0][0], &widths_0, &widths[0], &x_points_0, &x_points_1, &x_points[0][0], &box_vals[0][0])
    return numpy.asarray(box_vals, order='F')


cdef extern:
    void c_eval_boxes( long *boxes_0, long *boxes_1, double *boxes, long *widths_0, double *widths, long *weights_0, double *weights, long *x_points_0, long *x_points_1, double *x_points, double *response )

@cython.boundscheck(False)
@cython.wraparound(False)
def eval_boxes( double[:,:] boxes, double[:] widths, double[:] weights, double[:,:] x_points, o_response=None ):
    '''Created automatically by fmodpy. Documentation coming soon.'''
    if (type(boxes) != type(None)) and (not numpy.asarray(boxes).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(widths) != type(None)) and (not numpy.asarray(widths).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(weights) != type(None)) and (not numpy.asarray(weights).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(x_points) != type(None)) and (not numpy.asarray(x_points).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    if (type(o_response) != type(None)) and (not numpy.asarray(o_response).flags.f_contiguous):
        raise(Exception('ERROR: Only use numpy arrays that are f_contiguous.'))
    cdef long boxes_0 = boxes.shape[0]
    cdef long boxes_1 = boxes.shape[1]
    cdef long widths_0 = widths.shape[0]
    cdef long weights_0 = weights.shape[0]
    cdef long x_points_0 = x_points.shape[0]
    cdef long x_points_1 = x_points.shape[1]
    if type(o_response) == type(None):
        o_response = (1) * numpy.ones(shape=(x_points.shape[2]),dtype=float,order='F')
    cdef double [:] response = o_response
    
    c_eval_boxes(&boxes_0, &boxes_1, &boxes[0][0], &widths_0, &widths[0], &weights_0, &weights[0], &x_points_0, &x_points_1, &x_points[0][0], &response[0])
    return numpy.asarray(response, order='F')


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

