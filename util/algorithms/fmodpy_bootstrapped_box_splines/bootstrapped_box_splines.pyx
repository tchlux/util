'''! NAME:   Bootstrapped box splines
! AUTHOR: Thomas C.H. Lux
! EMAIL:  tchlux@vt.edu
!
! DESCRIPTION: This file (box_spline_basis.f95) contains the module
!              box_spline_basis that declares the subroutines
!              (compute_boxes, evaluate_boxes) for generating and
!              evaluating a box spline hypercube basis model for
!              continuous multidimensional (d<=20) data with real
!              response values.
!
! The following LAPACK routines are used:
!    + DGELS
!'''


import cython
import numpy

#      Wrapper for fortran function compute_boxes     
# =================================================

cdef extern:
    void c_compute_boxes( int* x_data_0, int* x_data_1, double* x_data, double* response, int* boxes_0, int* boxes_1, double* boxes, double* widths, double* weights, long* batch_size, bint* batch_size_present )

@cython.boundscheck(False)
@cython.wraparound(False)
def compute_boxes( double[:,:] x_data, double[:] response, double[:,:] boxes, widths=None, weights=None, batch_size=None ):
    '''! This is a serial implementation of the hypercube box-spline
    ! basis approximation algorithm.
    !
    ! compute_boxes(x_data, response, boxes, widths, weights)
    !   x_data   -- matrix of n-dimensions by m-data points (nxm) stored as
    !               column vectors. x_data(:,1) should be the first point.
    !   response -- array of length m holding response values at each x-point
    !   boxes    -- matrix of n-dimensions by p-boxes (n x p), will be
    !               overwritten with the output boxes
    !   widths   -- array of length p holding the width of each box
    !   weights  -- array of length p holding the weight for each box function'''
    # Prepare for fortran function call (initialize optionals)
    if (not numpy.asarray(x_data).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int x_data_0 = x_data.shape[0]
    cdef int x_data_1 = x_data.shape[1]
    
    if (not numpy.asarray(boxes).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int boxes_0 = boxes.shape[0]
    cdef int boxes_1 = boxes.shape[1]
    
    if (type(widths) == type(None)):
        widths = numpy.ones(shape=(boxes.shape[2-1]),dtype=numpy.float64,order='F')
    cdef double[:] local_widths = widths
    
    if (type(weights) == type(None)):
        weights = numpy.ones(shape=(boxes.shape[2-1]),dtype=numpy.float64,order='F')
    cdef double[:] local_weights = weights
    
    cdef bint batch_size_present = True
    if (type(batch_size) == type(None)):
        batch_size_present = False
        batch_size = 1
    cdef long local_batch_size = batch_size
    
    print("python batch_size: ", batch_size)
    
    # Make fortran function call
    c_compute_boxes(&x_data_0, &x_data_1, &x_data[0][0], &response[0], &boxes_0, &boxes_1, &boxes[0][0], &local_widths[0], &local_weights[0], &local_batch_size, &batch_size_present)
    # Return appropriate values based on "INTENT" from fortran code
    
    return numpy.asarray(boxes, order='F'), numpy.asarray(local_widths, order='F'), numpy.asarray(local_weights, order='F')



#      Wrapper for fortran function max_boxes     
# =================================================

cdef extern:
    void c_max_boxes( int* pts_0, int* pts_1, double* pts, int* box_centers_0, int* box_centers, double* box_widths )

@cython.boundscheck(False)
@cython.wraparound(False)
def max_boxes( double[:,:] pts, int[:] box_centers, box_widths=None ):
    ''''''
    # Prepare for fortran function call (initialize optionals)
    if (not numpy.asarray(pts).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int pts_0 = pts.shape[0]
    cdef int pts_1 = pts.shape[1]
    
    if (not numpy.asarray(box_centers).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int box_centers_0 = box_centers.shape[0]
    
    if (type(box_widths) == type(None)):
        box_widths = numpy.ones(shape=(2*box_centers.shape[1-1],box_centers.shape[1-1]),dtype=numpy.float64,order='F')
    cdef double[:,:] local_box_widths = box_widths
    
    
    # Make fortran function call
    c_max_boxes(&pts_0, &pts_1, &pts[0][0], &box_centers_0, &box_centers[0], &local_box_widths[0][0])
    # Return appropriate values based on "INTENT" from fortran code
    
    return numpy.asarray(local_box_widths, order='F')



#      Wrapper for fortran function eval_box_coefs     
# =================================================

cdef extern:
    void c_eval_box_coefs( int* boxes_0, int* boxes_1, double* boxes, int* widths_0, double* widths, int* x_points_0, int* x_points_1, double* x_points, double* box_vals )

@cython.boundscheck(False)
@cython.wraparound(False)
def eval_box_coefs( double[:,:] boxes, double[:] widths, double[:,:] x_points, box_vals=None ):
    '''! Subroutine for computing the values of all the boxes at an x-point'''
    # Prepare for fortran function call (initialize optionals)
    if (not numpy.asarray(boxes).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int boxes_0 = boxes.shape[0]
    cdef int boxes_1 = boxes.shape[1]
    
    if (not numpy.asarray(widths).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int widths_0 = widths.shape[0]
    
    if (not numpy.asarray(x_points).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int x_points_0 = x_points.shape[0]
    cdef int x_points_1 = x_points.shape[1]
    
    if (type(box_vals) == type(None)):
        box_vals = numpy.ones(shape=(boxes.shape[2-1],x_points.shape[2-1]),dtype=numpy.float64,order='F')
    cdef double[:,:] local_box_vals = box_vals
    
    
    # Make fortran function call
    c_eval_box_coefs(&boxes_0, &boxes_1, &boxes[0][0], &widths_0, &widths[0], &x_points_0, &x_points_1, &x_points[0][0], &local_box_vals[0][0])
    # Return appropriate values based on "INTENT" from fortran code
    
    return numpy.asarray(local_box_vals, order='F')



#      Wrapper for fortran function eval_boxes     
# =================================================

cdef extern:
    void c_eval_boxes( int* boxes_0, int* boxes_1, double* boxes, int* widths_0, double* widths, int* weights_0, double* weights, int* x_points_0, int* x_points_1, double* x_points, double* response )

@cython.boundscheck(False)
@cython.wraparound(False)
def eval_boxes( double[:,:] boxes, double[:] widths, double[:] weights, double[:,:] x_points, response=None ):
    '''! Subroutine for evaluating an approximation at a set of x-points'''
    # Prepare for fortran function call (initialize optionals)
    if (not numpy.asarray(boxes).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int boxes_0 = boxes.shape[0]
    cdef int boxes_1 = boxes.shape[1]
    
    if (not numpy.asarray(widths).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int widths_0 = widths.shape[0]
    
    if (not numpy.asarray(weights).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int weights_0 = weights.shape[0]
    
    if (not numpy.asarray(x_points).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int x_points_0 = x_points.shape[0]
    cdef int x_points_1 = x_points.shape[1]
    
    if (type(response) == type(None)):
        response = numpy.ones(shape=(x_points.shape[2-1]),dtype=numpy.float64,order='F')
    cdef double[:] local_response = response
    
    
    # Make fortran function call
    c_eval_boxes(&boxes_0, &boxes_1, &boxes[0][0], &widths_0, &widths[0], &weights_0, &weights[0], &x_points_0, &x_points_1, &x_points[0][0], &local_response[0])
    # Return appropriate values based on "INTENT" from fortran code
    
    return numpy.asarray(local_response, order='F')



#      Wrapper for fortran function qsortc     
# =================================================

cdef extern:
    void c_qsortc( int* a_0, double* a, int* idx )

@cython.boundscheck(False)
@cython.wraparound(False)
def qsortc( double[:] a, idx=None ):
    '''! This is a QuickSort routine adapted from Orderpack 2.0.
    !
    ! Also, this implementation incorporates ideas from "A Practical Introduction
    ! to Data Structures and Algorithm Analysis", by Clifford Shaffer.
    !
    ! It sorts real numbers into ascending numerical order
    ! and keeps an index of the value's original array position.
    !
    ! Author: Will Thacker, Winthrop University, July 2013.
    !
    ! QSORTC sorts the real array A and keeps the index of the value's
    ! original array position along with the value (integer array IDX).
    !
    ! On input:
    !
    ! A(:) is the array to be sorted.
    !
    ! On output:
    !
    ! A(:) is sorted.
    !
    ! IDX(1:SIZEOF(A)) contains the original positions of the sorted values.
    !    I.e., sorted(i) = orginal_unsorted(IDX(i)).
    !
    !'''
    # Prepare for fortran function call (initialize optionals)
    if (not numpy.asarray(a).flags.f_contiguous):
        raise(Exception('Only use numpy arrays that are f_contiguous.'))
    cdef int a_0 = a.shape[0]
    
    if (type(idx) == type(None)):
        idx = numpy.ones(shape=(a.size),dtype=numpy.int32,order='F')
    cdef int[:] local_idx = idx
    
    
    # Make fortran function call
    c_qsortc(&a_0, &a[0], &local_idx[0])
    # Return appropriate values based on "INTENT" from fortran code
    
    return numpy.asarray(a, order='F'), numpy.asarray(local_idx, order='F')



