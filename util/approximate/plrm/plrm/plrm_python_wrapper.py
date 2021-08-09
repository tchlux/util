'''This Python code is an automatically generated wrapper
for Fortran code made by 'fmodpy'. The original documentation
for the Fortran source code follows.

! TODO:
! - Compute value and error gradient with respect to a set of points
!   for a single parameter in the model (vec, shift, or flex).
! - Compute the "importance" of different vectors by zeroing them
!   out and determining the magnitude change in mean squared error.
! - Well-spaced validation data (multiply distances in input and output).
! - Validation-based training and greedy saving / model selection.
!
! - Setting that produces convergence status on timed interval (3 seconds?).
! - Support for vectorizing integer valued inputs (initially regular simplex, trainable).
! - Training support for vector of assigned neighbors for each point.
! - Support for huge input spaces (1M+ components), automatically
!   construct vectors that select a subset of input components to
!   process.
!
! - Save model and load model (during training).
! - Support for vectorizing integer valued outputs. (use regular simplex)
! - Train only the output, internal, or input layers.
! - Consider convergence when adding new vectors.
! - Implement similar code in C, compare speeds.
'''

import os
import ctypes
import numpy

# --------------------------------------------------------------------
#               CONFIGURATION
# 
_verbose = True
_fort_compiler = "gfortran"
_shared_object_name = "plrm.so"
_this_directory = os.path.dirname(os.path.abspath(__file__))
_path_to_lib = os.path.join(_this_directory, _shared_object_name)
_compile_options = ['-fPIC', '-shared', '-O3', '-lblas', '-fopenmp']
_ordered_dependencies = ['plrm.f90', 'plrm_c_wrapper.f90']
# 
# --------------------------------------------------------------------
#               AUTO-COMPILING
#
# Try to import the existing object. If that fails, recompile and then try.
try:
    clib = ctypes.CDLL(_path_to_lib)
except:
    # Remove the shared object if it exists, because it is faulty.
    if os.path.exists(_shared_object_name):
        os.remove(_shared_object_name)
    # Compile a new shared object.
    _command = " ".join([_fort_compiler] + _compile_options + ["-o", _shared_object_name] + _ordered_dependencies)
    if _verbose:
        print("Running system command with arguments")
        print("  ", _command)
    # Run the compilation command.
    import subprocess
    subprocess.run(_command, shell=True, cwd=_this_directory)
    # Import the shared object file as a C library with ctypes.
    clib = ctypes.CDLL(_path_to_lib)
# --------------------------------------------------------------------


class plrm:
    '''! A piecewise linear regression model.'''

    # Declare 'discontinuity'
    def get_discontinuity(self):
        discontinuity = ctypes.c_float()
        clib.plrm_get_discontinuity(ctypes.byref(discontinuity))
        return discontinuity.value
    def set_discontinuity(self, discontinuity):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    discontinuity = property(get_discontinuity, set_discontinuity)

    # Declare 'initial_shift_range'
    def get_initial_shift_range(self):
        initial_shift_range = ctypes.c_float()
        clib.plrm_get_initial_shift_range(ctypes.byref(initial_shift_range))
        return initial_shift_range.value
    def set_initial_shift_range(self, initial_shift_range):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    initial_shift_range = property(get_initial_shift_range, set_initial_shift_range)

    # Declare 'initial_flex'
    def get_initial_flex(self):
        initial_flex = ctypes.c_float()
        clib.plrm_get_initial_flex(ctypes.byref(initial_flex))
        return initial_flex.value
    def set_initial_flex(self, initial_flex):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    initial_flex = property(get_initial_flex, set_initial_flex)

    # Declare 'initial_output_scale'
    def get_initial_output_scale(self):
        initial_output_scale = ctypes.c_float()
        clib.plrm_get_initial_output_scale(ctypes.byref(initial_output_scale))
        return initial_output_scale.value
    def set_initial_output_scale(self, initial_output_scale):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    initial_output_scale = property(get_initial_output_scale, set_initial_output_scale)

    # Declare 'initial_step'
    def get_initial_step(self):
        initial_step = ctypes.c_float()
        clib.plrm_get_initial_step(ctypes.byref(initial_step))
        return initial_step.value
    def set_initial_step(self, initial_step):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    initial_step = property(get_initial_step, set_initial_step)

    # Declare 'initial_step_mean_change'
    def get_initial_step_mean_change(self):
        initial_step_mean_change = ctypes.c_float()
        clib.plrm_get_initial_step_mean_change(ctypes.byref(initial_step_mean_change))
        return initial_step_mean_change.value
    def set_initial_step_mean_change(self, initial_step_mean_change):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    initial_step_mean_change = property(get_initial_step_mean_change, set_initial_step_mean_change)

    # Declare 'initial_step_var_change'
    def get_initial_step_var_change(self):
        initial_step_var_change = ctypes.c_float()
        clib.plrm_get_initial_step_var_change(ctypes.byref(initial_step_var_change))
        return initial_step_var_change.value
    def set_initial_step_var_change(self, initial_step_var_change):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    initial_step_var_change = property(get_initial_step_var_change, set_initial_step_var_change)

    # Declare 'step_growth_rate'
    def get_step_growth_rate(self):
        step_growth_rate = ctypes.c_float()
        clib.plrm_get_step_growth_rate(ctypes.byref(step_growth_rate))
        return step_growth_rate.value
    def set_step_growth_rate(self, step_growth_rate):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    step_growth_rate = property(get_step_growth_rate, set_step_growth_rate)

    # Declare 'step_shrink_rate'
    def get_step_shrink_rate(self):
        step_shrink_rate = ctypes.c_float()
        clib.plrm_get_step_shrink_rate(ctypes.byref(step_shrink_rate))
        return step_shrink_rate.value
    def set_step_shrink_rate(self, step_shrink_rate):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    step_shrink_rate = property(get_step_shrink_rate, set_step_shrink_rate)

    # Declare 'min_steps_to_stability'
    def get_min_steps_to_stability(self):
        min_steps_to_stability = ctypes.c_int()
        clib.plrm_get_min_steps_to_stability(ctypes.byref(min_steps_to_stability))
        return min_steps_to_stability.value
    def set_min_steps_to_stability(self, min_steps_to_stability):
        raise(NotImplementedError('Module attributes with PARAMETER status cannot be set.'))
    min_steps_to_stability = property(get_min_steps_to_stability, set_min_steps_to_stability)

    # Declare 'mdi'
    def get_mdi(self):
        mdi = ctypes.c_int()
        clib.plrm_get_mdi(ctypes.byref(mdi))
        return mdi.value
    def set_mdi(self, mdi):
        mdi = ctypes.c_int(mdi)
        clib.plrm_set_mdi(ctypes.byref(mdi))
    mdi = property(get_mdi, set_mdi)

    # Declare 'mds'
    def get_mds(self):
        mds = ctypes.c_int()
        clib.plrm_get_mds(ctypes.byref(mds))
        return mds.value
    def set_mds(self, mds):
        mds = ctypes.c_int(mds)
        clib.plrm_set_mds(ctypes.byref(mds))
    mds = property(get_mds, set_mds)

    # Declare 'mns'
    def get_mns(self):
        mns = ctypes.c_int()
        clib.plrm_get_mns(ctypes.byref(mns))
        return mns.value
    def set_mns(self, mns):
        mns = ctypes.c_int(mns)
        clib.plrm_set_mns(ctypes.byref(mns))
    mns = property(get_mns, set_mns)

    # Declare 'mdo'
    def get_mdo(self):
        mdo = ctypes.c_int()
        clib.plrm_get_mdo(ctypes.byref(mdo))
        return mdo.value
    def set_mdo(self, mdo):
        mdo = ctypes.c_int(mdo)
        clib.plrm_set_mdo(ctypes.byref(mdo))
    mdo = property(get_mdo, set_mdo)

    # Declare 'input_vecs'
    def get_input_vecs(self):
        input_vecs_allocated = ctypes.c_bool(False)
        input_vecs_dim_1 = ctypes.c_int()
        input_vecs_dim_2 = ctypes.c_int()
        input_vecs = ctypes.c_void_p()
        clib.plrm_get_input_vecs(ctypes.byref(input_vecs_allocated), ctypes.byref(input_vecs_dim_1), ctypes.byref(input_vecs_dim_2), ctypes.byref(input_vecs))
        if (not input_vecs_allocated.value): return None
        input_vecs_size = (input_vecs_dim_1.value) * (input_vecs_dim_2.value)
        if (input_vecs_size > 0):
            input_vecs = numpy.array(ctypes.cast(input_vecs, ctypes.POINTER(ctypes.c_float*input_vecs_size)).contents, copy=False)
        else:
            input_vecs = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        input_vecs = input_vecs.reshape(input_vecs_dim_2.value,input_vecs_dim_1.value).T
        return input_vecs
    def set_input_vecs(self, input_vecs):
        if ((not issubclass(type(input_vecs), numpy.ndarray)) or
            (not numpy.asarray(input_vecs).flags.f_contiguous) or
            (not (input_vecs.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'input_vecs' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            input_vecs = numpy.asarray(input_vecs, dtype=ctypes.c_float, order='F')
        input_vecs_dim_1 = ctypes.c_int(input_vecs.shape[0])
        input_vecs_dim_2 = ctypes.c_int(input_vecs.shape[1])
        clib.plrm_set_input_vecs(ctypes.byref(input_vecs_dim_1), ctypes.byref(input_vecs_dim_2), ctypes.c_void_p(input_vecs.ctypes.data))
    input_vecs = property(get_input_vecs, set_input_vecs)

    # Declare 'input_shift'
    def get_input_shift(self):
        input_shift_allocated = ctypes.c_bool(False)
        input_shift_dim_1 = ctypes.c_int()
        input_shift = ctypes.c_void_p()
        clib.plrm_get_input_shift(ctypes.byref(input_shift_allocated), ctypes.byref(input_shift_dim_1), ctypes.byref(input_shift))
        if (not input_shift_allocated.value): return None
        input_shift_size = (input_shift_dim_1.value)
        if (input_shift_size > 0):
            input_shift = numpy.array(ctypes.cast(input_shift, ctypes.POINTER(ctypes.c_float*input_shift_size)).contents, copy=False)
        else:
            input_shift = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        return input_shift
    def set_input_shift(self, input_shift):
        if ((not issubclass(type(input_shift), numpy.ndarray)) or
            (not numpy.asarray(input_shift).flags.f_contiguous) or
            (not (input_shift.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'input_shift' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            input_shift = numpy.asarray(input_shift, dtype=ctypes.c_float, order='F')
        input_shift_dim_1 = ctypes.c_int(input_shift.shape[0])
        clib.plrm_set_input_shift(ctypes.byref(input_shift_dim_1), ctypes.c_void_p(input_shift.ctypes.data))
    input_shift = property(get_input_shift, set_input_shift)

    # Declare 'input_flex'
    def get_input_flex(self):
        input_flex_allocated = ctypes.c_bool(False)
        input_flex_dim_1 = ctypes.c_int()
        input_flex = ctypes.c_void_p()
        clib.plrm_get_input_flex(ctypes.byref(input_flex_allocated), ctypes.byref(input_flex_dim_1), ctypes.byref(input_flex))
        if (not input_flex_allocated.value): return None
        input_flex_size = (input_flex_dim_1.value)
        if (input_flex_size > 0):
            input_flex = numpy.array(ctypes.cast(input_flex, ctypes.POINTER(ctypes.c_float*input_flex_size)).contents, copy=False)
        else:
            input_flex = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        return input_flex
    def set_input_flex(self, input_flex):
        if ((not issubclass(type(input_flex), numpy.ndarray)) or
            (not numpy.asarray(input_flex).flags.f_contiguous) or
            (not (input_flex.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'input_flex' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            input_flex = numpy.asarray(input_flex, dtype=ctypes.c_float, order='F')
        input_flex_dim_1 = ctypes.c_int(input_flex.shape[0])
        clib.plrm_set_input_flex(ctypes.byref(input_flex_dim_1), ctypes.c_void_p(input_flex.ctypes.data))
    input_flex = property(get_input_flex, set_input_flex)

    # Declare 'internal_vecs'
    def get_internal_vecs(self):
        internal_vecs_allocated = ctypes.c_bool(False)
        internal_vecs_dim_1 = ctypes.c_int()
        internal_vecs_dim_2 = ctypes.c_int()
        internal_vecs_dim_3 = ctypes.c_int()
        internal_vecs = ctypes.c_void_p()
        clib.plrm_get_internal_vecs(ctypes.byref(internal_vecs_allocated), ctypes.byref(internal_vecs_dim_1), ctypes.byref(internal_vecs_dim_2), ctypes.byref(internal_vecs_dim_3), ctypes.byref(internal_vecs))
        if (not internal_vecs_allocated.value): return None
        internal_vecs_size = (internal_vecs_dim_1.value) * (internal_vecs_dim_2.value) * (internal_vecs_dim_3.value)
        if (internal_vecs_size > 0):
            internal_vecs = numpy.array(ctypes.cast(internal_vecs, ctypes.POINTER(ctypes.c_float*internal_vecs_size)).contents, copy=False)
        else:
            internal_vecs = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        internal_vecs = internal_vecs.reshape(internal_vecs_dim_3.value,internal_vecs_dim_2.value,internal_vecs_dim_1.value).T
        return internal_vecs
    def set_internal_vecs(self, internal_vecs):
        if ((not issubclass(type(internal_vecs), numpy.ndarray)) or
            (not numpy.asarray(internal_vecs).flags.f_contiguous) or
            (not (internal_vecs.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'internal_vecs' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            internal_vecs = numpy.asarray(internal_vecs, dtype=ctypes.c_float, order='F')
        internal_vecs_dim_1 = ctypes.c_int(internal_vecs.shape[0])
        internal_vecs_dim_2 = ctypes.c_int(internal_vecs.shape[1])
        internal_vecs_dim_3 = ctypes.c_int(internal_vecs.shape[2])
        clib.plrm_set_internal_vecs(ctypes.byref(internal_vecs_dim_1), ctypes.byref(internal_vecs_dim_2), ctypes.byref(internal_vecs_dim_3), ctypes.c_void_p(internal_vecs.ctypes.data))
    internal_vecs = property(get_internal_vecs, set_internal_vecs)

    # Declare 'internal_shift'
    def get_internal_shift(self):
        internal_shift_allocated = ctypes.c_bool(False)
        internal_shift_dim_1 = ctypes.c_int()
        internal_shift_dim_2 = ctypes.c_int()
        internal_shift = ctypes.c_void_p()
        clib.plrm_get_internal_shift(ctypes.byref(internal_shift_allocated), ctypes.byref(internal_shift_dim_1), ctypes.byref(internal_shift_dim_2), ctypes.byref(internal_shift))
        if (not internal_shift_allocated.value): return None
        internal_shift_size = (internal_shift_dim_1.value) * (internal_shift_dim_2.value)
        if (internal_shift_size > 0):
            internal_shift = numpy.array(ctypes.cast(internal_shift, ctypes.POINTER(ctypes.c_float*internal_shift_size)).contents, copy=False)
        else:
            internal_shift = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        internal_shift = internal_shift.reshape(internal_shift_dim_2.value,internal_shift_dim_1.value).T
        return internal_shift
    def set_internal_shift(self, internal_shift):
        if ((not issubclass(type(internal_shift), numpy.ndarray)) or
            (not numpy.asarray(internal_shift).flags.f_contiguous) or
            (not (internal_shift.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'internal_shift' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            internal_shift = numpy.asarray(internal_shift, dtype=ctypes.c_float, order='F')
        internal_shift_dim_1 = ctypes.c_int(internal_shift.shape[0])
        internal_shift_dim_2 = ctypes.c_int(internal_shift.shape[1])
        clib.plrm_set_internal_shift(ctypes.byref(internal_shift_dim_1), ctypes.byref(internal_shift_dim_2), ctypes.c_void_p(internal_shift.ctypes.data))
    internal_shift = property(get_internal_shift, set_internal_shift)

    # Declare 'internal_flex'
    def get_internal_flex(self):
        internal_flex_allocated = ctypes.c_bool(False)
        internal_flex_dim_1 = ctypes.c_int()
        internal_flex_dim_2 = ctypes.c_int()
        internal_flex = ctypes.c_void_p()
        clib.plrm_get_internal_flex(ctypes.byref(internal_flex_allocated), ctypes.byref(internal_flex_dim_1), ctypes.byref(internal_flex_dim_2), ctypes.byref(internal_flex))
        if (not internal_flex_allocated.value): return None
        internal_flex_size = (internal_flex_dim_1.value) * (internal_flex_dim_2.value)
        if (internal_flex_size > 0):
            internal_flex = numpy.array(ctypes.cast(internal_flex, ctypes.POINTER(ctypes.c_float*internal_flex_size)).contents, copy=False)
        else:
            internal_flex = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        internal_flex = internal_flex.reshape(internal_flex_dim_2.value,internal_flex_dim_1.value).T
        return internal_flex
    def set_internal_flex(self, internal_flex):
        if ((not issubclass(type(internal_flex), numpy.ndarray)) or
            (not numpy.asarray(internal_flex).flags.f_contiguous) or
            (not (internal_flex.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'internal_flex' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            internal_flex = numpy.asarray(internal_flex, dtype=ctypes.c_float, order='F')
        internal_flex_dim_1 = ctypes.c_int(internal_flex.shape[0])
        internal_flex_dim_2 = ctypes.c_int(internal_flex.shape[1])
        clib.plrm_set_internal_flex(ctypes.byref(internal_flex_dim_1), ctypes.byref(internal_flex_dim_2), ctypes.c_void_p(internal_flex.ctypes.data))
    internal_flex = property(get_internal_flex, set_internal_flex)

    # Declare 'output_vecs'
    def get_output_vecs(self):
        output_vecs_allocated = ctypes.c_bool(False)
        output_vecs_dim_1 = ctypes.c_int()
        output_vecs_dim_2 = ctypes.c_int()
        output_vecs = ctypes.c_void_p()
        clib.plrm_get_output_vecs(ctypes.byref(output_vecs_allocated), ctypes.byref(output_vecs_dim_1), ctypes.byref(output_vecs_dim_2), ctypes.byref(output_vecs))
        if (not output_vecs_allocated.value): return None
        output_vecs_size = (output_vecs_dim_1.value) * (output_vecs_dim_2.value)
        if (output_vecs_size > 0):
            output_vecs = numpy.array(ctypes.cast(output_vecs, ctypes.POINTER(ctypes.c_float*output_vecs_size)).contents, copy=False)
        else:
            output_vecs = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        output_vecs = output_vecs.reshape(output_vecs_dim_2.value,output_vecs_dim_1.value).T
        return output_vecs
    def set_output_vecs(self, output_vecs):
        if ((not issubclass(type(output_vecs), numpy.ndarray)) or
            (not numpy.asarray(output_vecs).flags.f_contiguous) or
            (not (output_vecs.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'output_vecs' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            output_vecs = numpy.asarray(output_vecs, dtype=ctypes.c_float, order='F')
        output_vecs_dim_1 = ctypes.c_int(output_vecs.shape[0])
        output_vecs_dim_2 = ctypes.c_int(output_vecs.shape[1])
        clib.plrm_set_output_vecs(ctypes.byref(output_vecs_dim_1), ctypes.byref(output_vecs_dim_2), ctypes.c_void_p(output_vecs.ctypes.data))
    output_vecs = property(get_output_vecs, set_output_vecs)

    # Declare 'output_shift'
    def get_output_shift(self):
        output_shift_allocated = ctypes.c_bool(False)
        output_shift_dim_1 = ctypes.c_int()
        output_shift = ctypes.c_void_p()
        clib.plrm_get_output_shift(ctypes.byref(output_shift_allocated), ctypes.byref(output_shift_dim_1), ctypes.byref(output_shift))
        if (not output_shift_allocated.value): return None
        output_shift_size = (output_shift_dim_1.value)
        if (output_shift_size > 0):
            output_shift = numpy.array(ctypes.cast(output_shift, ctypes.POINTER(ctypes.c_float*output_shift_size)).contents, copy=False)
        else:
            output_shift = numpy.zeros((0,), dtype=ctypes.c_float, order='F')
        return output_shift
    def set_output_shift(self, output_shift):
        if ((not issubclass(type(output_shift), numpy.ndarray)) or
            (not numpy.asarray(output_shift).flags.f_contiguous) or
            (not (output_shift.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'output_shift' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            output_shift = numpy.asarray(output_shift, dtype=ctypes.c_float, order='F')
        output_shift_dim_1 = ctypes.c_int(output_shift.shape[0])
        clib.plrm_set_output_shift(ctypes.byref(output_shift_dim_1), ctypes.c_void_p(output_shift.ctypes.data))
    output_shift = property(get_output_shift, set_output_shift)

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine NEW_MODEL
    
    def new_model(self, di, ds, ns, do):
        '''! Allocate storage for a new model given paramters:
!   DI -- Dimension of input vector.
!   DS -- Dimension of internal vectors spaces.
!   NS -- Number of sequential internal vector spaces.
!   DO -- Dimension of output vector.'''
        
        # Setting up "di"
        if (type(di) is not ctypes.c_int): di = ctypes.c_int(di)
        
        # Setting up "ds"
        if (type(ds) is not ctypes.c_int): ds = ctypes.c_int(ds)
        
        # Setting up "ns"
        if (type(ns) is not ctypes.c_int): ns = ctypes.c_int(ns)
        
        # Setting up "do"
        if (type(do) is not ctypes.c_int): do = ctypes.c_int(do)
    
        # Call C-accessible Fortran wrapper.
        clib.c_new_model(ctypes.byref(di), ctypes.byref(ds), ctypes.byref(ns), ctypes.byref(do))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return 

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine INIT_MODEL
    
    def init_model(self, seed=None):
        '''! Initialize the weights for a model, optionally provide a random seed.'''
        
        # Setting up "seed"
        seed_present = ctypes.c_bool(True)
        if (seed is None):
            seed_present = ctypes.c_bool(False)
            seed = ctypes.c_int()
        if (type(seed) is not ctypes.c_int): seed = ctypes.c_int(seed)
    
        # Call C-accessible Fortran wrapper.
        clib.c_init_model(ctypes.byref(seed_present), ctypes.byref(seed))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return 

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine EVALUATE
    
    def evaluate(self, inputs, outputs):
        '''! Evaluate the piecewise linear regression model.'''
        
        # Setting up "inputs"
        if ((not issubclass(type(inputs), numpy.ndarray)) or
            (not numpy.asarray(inputs).flags.f_contiguous) or
            (not (inputs.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'inputs' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            inputs = numpy.asarray(inputs, dtype=ctypes.c_float, order='F')
        inputs_dim_1 = ctypes.c_int(inputs.shape[0])
        inputs_dim_2 = ctypes.c_int(inputs.shape[1])
        
        # Setting up "outputs"
        if ((not issubclass(type(outputs), numpy.ndarray)) or
            (not numpy.asarray(outputs).flags.f_contiguous) or
            (not (outputs.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'outputs' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            outputs = numpy.asarray(outputs, dtype=ctypes.c_float, order='F')
        outputs_dim_1 = ctypes.c_int(outputs.shape[0])
        outputs_dim_2 = ctypes.c_int(outputs.shape[1])
    
        # Call C-accessible Fortran wrapper.
        clib.c_evaluate(ctypes.byref(inputs_dim_1), ctypes.byref(inputs_dim_2), ctypes.c_void_p(inputs.ctypes.data), ctypes.byref(outputs_dim_1), ctypes.byref(outputs_dim_2), ctypes.c_void_p(outputs.ctypes.data))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return outputs

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine MINIMIZE_MSE
    
    def minimize_mse(self, x, y, steps, batch_size=None, num_threads=None, step_size=None, keep_best=None, record=None):
        '''! Fit input / output pairs by minimizing mean squared error.'''
        
        # Setting up "x"
        if ((not issubclass(type(x), numpy.ndarray)) or
            (not numpy.asarray(x).flags.f_contiguous) or
            (not (x.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'x' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            x = numpy.asarray(x, dtype=ctypes.c_float, order='F')
        x_dim_1 = ctypes.c_int(x.shape[0])
        x_dim_2 = ctypes.c_int(x.shape[1])
        
        # Setting up "y"
        if ((not issubclass(type(y), numpy.ndarray)) or
            (not numpy.asarray(y).flags.f_contiguous) or
            (not (y.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'y' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            y = numpy.asarray(y, dtype=ctypes.c_float, order='F')
        y_dim_1 = ctypes.c_int(y.shape[0])
        y_dim_2 = ctypes.c_int(y.shape[1])
        
        # Setting up "steps"
        if (type(steps) is not ctypes.c_int): steps = ctypes.c_int(steps)
        
        # Setting up "batch_size"
        batch_size_present = ctypes.c_bool(True)
        if (batch_size is None):
            batch_size_present = ctypes.c_bool(False)
            batch_size = ctypes.c_int()
        if (type(batch_size) is not ctypes.c_int): batch_size = ctypes.c_int(batch_size)
        
        # Setting up "num_threads"
        num_threads_present = ctypes.c_bool(True)
        if (num_threads is None):
            num_threads_present = ctypes.c_bool(False)
            num_threads = ctypes.c_int()
        if (type(num_threads) is not ctypes.c_int): num_threads = ctypes.c_int(num_threads)
        
        # Setting up "step_size"
        step_size_present = ctypes.c_bool(True)
        if (step_size is None):
            step_size_present = ctypes.c_bool(False)
            step_size = ctypes.c_float()
        if (type(step_size) is not ctypes.c_float): step_size = ctypes.c_float(step_size)
        
        # Setting up "keep_best"
        keep_best_present = ctypes.c_bool(True)
        if (keep_best is None):
            keep_best_present = ctypes.c_bool(False)
            keep_best = ctypes.c_bool()
        if (type(keep_best) is not ctypes.c_bool): keep_best = ctypes.c_bool(keep_best)
        
        # Setting up "mean_squared_error"
        mean_squared_error = ctypes.c_float()
        
        # Setting up "record"
        record_present = ctypes.c_bool(True)
        if (record is None):
            record_present = ctypes.c_bool(False)
            record = numpy.zeros(shape=(1), dtype=ctypes.c_float, order='F')
        elif (type(record) == bool) and (record):
            record = numpy.zeros(shape=(steps), dtype=ctypes.c_float, order='F')
        elif ((not issubclass(type(record), numpy.ndarray)) or
              (not numpy.asarray(record).flags.f_contiguous) or
              (not (record.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'record' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            record = numpy.asarray(record, dtype=ctypes.c_float, order='F')
        if (record_present):
            record_dim_1 = ctypes.c_int(record.shape[0])
        else:
            record_dim_1 = ctypes.c_int()
    
        # Call C-accessible Fortran wrapper.
        clib.c_minimize_mse(ctypes.byref(x_dim_1), ctypes.byref(x_dim_2), ctypes.c_void_p(x.ctypes.data), ctypes.byref(y_dim_1), ctypes.byref(y_dim_2), ctypes.c_void_p(y.ctypes.data), ctypes.byref(steps), ctypes.byref(batch_size_present), ctypes.byref(batch_size), ctypes.byref(num_threads_present), ctypes.byref(num_threads), ctypes.byref(step_size_present), ctypes.byref(step_size), ctypes.byref(keep_best_present), ctypes.byref(keep_best), ctypes.byref(mean_squared_error), ctypes.byref(record_present), ctypes.byref(record_dim_1), ctypes.c_void_p(record.ctypes.data))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return mean_squared_error.value, (record if record_present else None)

plrm = plrm()


class vis_plrm:
    '''! Use this module to visualize the values and error gradients within
!  the PLRM model for different parameter types.'''

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine DISABLE_AND_COMPUTE_MSE
    
    def disable_and_compute_mse(self, layer, position, inputs, outputs):
        '''! Disable the 'layer, position' node (1 indexed) in the network and
!  compute the mean squared error "MSE". Use layer=-1 to compute
!  MSE without disabling any of the internal representations.'''
        
        # Setting up "layer"
        if (type(layer) is not ctypes.c_int): layer = ctypes.c_int(layer)
        
        # Setting up "position"
        if (type(position) is not ctypes.c_int): position = ctypes.c_int(position)
        
        # Setting up "inputs"
        if ((not issubclass(type(inputs), numpy.ndarray)) or
            (not numpy.asarray(inputs).flags.f_contiguous) or
            (not (inputs.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'inputs' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            inputs = numpy.asarray(inputs, dtype=ctypes.c_float, order='F')
        inputs_dim_1 = ctypes.c_int(inputs.shape[0])
        inputs_dim_2 = ctypes.c_int(inputs.shape[1])
        
        # Setting up "outputs"
        if ((not issubclass(type(outputs), numpy.ndarray)) or
            (not numpy.asarray(outputs).flags.f_contiguous) or
            (not (outputs.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'outputs' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            outputs = numpy.asarray(outputs, dtype=ctypes.c_float, order='F')
        outputs_dim_1 = ctypes.c_int(outputs.shape[0])
        outputs_dim_2 = ctypes.c_int(outputs.shape[1])
        
        # Setting up "mse"
        mse = ctypes.c_float()
    
        # Call C-accessible Fortran wrapper.
        clib.c_disable_and_compute_mse(ctypes.byref(layer), ctypes.byref(position), ctypes.byref(inputs_dim_1), ctypes.byref(inputs_dim_2), ctypes.c_void_p(inputs.ctypes.data), ctypes.byref(outputs_dim_1), ctypes.byref(outputs_dim_2), ctypes.c_void_p(outputs.ctypes.data), ctypes.byref(mse))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return mse.value

    
    # ----------------------------------------------
    # Wrapper for the Fortran subroutine COMPUTE_VALUES
    
    def compute_values(self, layer, position, inputs, outputs, vals, grads):
        '''! Given index, compute value and gradient at internal position in model.
! Arguments'''
        
        # Setting up "layer"
        if (type(layer) is not ctypes.c_int): layer = ctypes.c_int(layer)
        
        # Setting up "position"
        if (type(position) is not ctypes.c_int): position = ctypes.c_int(position)
        
        # Setting up "inputs"
        if ((not issubclass(type(inputs), numpy.ndarray)) or
            (not numpy.asarray(inputs).flags.f_contiguous) or
            (not (inputs.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'inputs' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            inputs = numpy.asarray(inputs, dtype=ctypes.c_float, order='F')
        inputs_dim_1 = ctypes.c_int(inputs.shape[0])
        inputs_dim_2 = ctypes.c_int(inputs.shape[1])
        
        # Setting up "outputs"
        if ((not issubclass(type(outputs), numpy.ndarray)) or
            (not numpy.asarray(outputs).flags.f_contiguous) or
            (not (outputs.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'outputs' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            outputs = numpy.asarray(outputs, dtype=ctypes.c_float, order='F')
        outputs_dim_1 = ctypes.c_int(outputs.shape[0])
        outputs_dim_2 = ctypes.c_int(outputs.shape[1])
        
        # Setting up "vals"
        if ((not issubclass(type(vals), numpy.ndarray)) or
            (not numpy.asarray(vals).flags.f_contiguous) or
            (not (vals.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'vals' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            vals = numpy.asarray(vals, dtype=ctypes.c_float, order='F')
        vals_dim_1 = ctypes.c_int(vals.shape[0])
        
        # Setting up "grads"
        if ((not issubclass(type(grads), numpy.ndarray)) or
            (not numpy.asarray(grads).flags.f_contiguous) or
            (not (grads.dtype == numpy.dtype(ctypes.c_float)))):
            import warnings
            warnings.warn("The provided argument 'grads' was not an f_contiguous NumPy array of type 'ctypes.c_float' (or equivalent). Automatically converting (probably creating a full copy).")
            grads = numpy.asarray(grads, dtype=ctypes.c_float, order='F')
        grads_dim_1 = ctypes.c_int(grads.shape[0])
    
        # Call C-accessible Fortran wrapper.
        clib.c_compute_values(ctypes.byref(layer), ctypes.byref(position), ctypes.byref(inputs_dim_1), ctypes.byref(inputs_dim_2), ctypes.c_void_p(inputs.ctypes.data), ctypes.byref(outputs_dim_1), ctypes.byref(outputs_dim_2), ctypes.c_void_p(outputs.ctypes.data), ctypes.byref(vals_dim_1), ctypes.c_void_p(vals.ctypes.data), ctypes.byref(grads_dim_1), ctypes.c_void_p(grads.ctypes.data))
    
        # Return final results, 'INTENT(OUT)' arguments only.
        return vals, grads

vis_plrm = vis_plrm()

