import numpy, os, sysconfig
from distutils.extension import Extension
from distutils.core import setup
from Cython.Build import build_ext

# LINK_FILES = ['linear_shepard.o','linear_shepard_wrapper.o']
# LINK_FILES = [f for f in os.listdir("/s/tchlux/TestingData/fmodpy_linear_shepard") if (f[-2:] == ".o")]
LINK_FILES = [f for f in os.listdir() if (f[-2:] == ".o")]
MODULE_NAME = 'linear_shepard'
CYTHON_SOURCE = ['linear_shepard.pyx']
COMPILE_ARGS = ['-fPIC', '-O3']
LINK_ARGS = ['-llapack', '-lblas']

# This might not work on all systems, need to investigate further
os.environ["CC"] = "gfortran"
# Get the linker options used to build python
linker_options = sysconfig.get_config_vars().get("BLDSHARED","").split(" ")[1:]
# Set the linker, with appropriate options
os.environ["LDSHARED"] = "gfortran "+" ".join(linker_options)

print("SETUP_PY: Compiling extension module with '%s'."%os.environ["CC"])

ext_modules = [Extension(MODULE_NAME,
                         CYTHON_SOURCE,
                         extra_compile_args=COMPILE_ARGS,
                         extra_link_args=LINK_ARGS+LINK_FILES)]

print("SETUP_PY: Extension module read, calling setup...")

setup(name = MODULE_NAME,
      cmdclass = {'build_ext': build_ext},
      include_dirs = [numpy.get_include()],
      # ^^ Includes the NumPy headers when compiling.
      ext_modules = ext_modules)

print("SETUP_PY: Done setting up.")
