import os
import fmodpy

# 1) Wrap a file that has one module and within that module
#    contains the headers for the delaunay subroutines.
# 2) Delete temporary Fortran module file, modify "_c_to_f.f90" file
#    to include the *true* module that is being wrapped.
# 3) Change the names of the generate itermediate files to share the
#    name of the file with the *true* module that is being wrapped.
# 4) Build the module using the "build_mod" function directly.

# ------------------------------ 1 -----------------------------------
# fmodpy.wrap("delaunay.f90", working_directory="fmodpy_delaunay", verbose=True)

# ------------------------------ 2 -----------------------------------
# ------------------------------ 3 -----------------------------------

# ------------------------------ 4 -----------------------------------

fmodpy.module_link_args += ["-fopenmp"]
working_directory = os.path.abspath(".")

try:
    fmodpy.build_mod("delsparse.f90", working_directory,
                     "delsparse", verbose=True)
except SystemError:
    # Because this is a sub-module, the system gets angry when a
    # relative import is attempted. Ignore this error.
    pass
