RM=rm -f
FORT=gfortran
FORT_FLAGS=-O3 -fPIC
PY=python3
CPY=cython

all:
	$(FORT) -c $(FORT_FLAGS) real_precision.f95 linear_shepard.f95 linear_shepard_wrapper.f03
	$(CPY) linear_shepard.pyx
	$(PY) linear_shepard_setup.py build_ext --inplace

clean:
	$(RM) *.o *.mod *.c *.pyc translated_fortran.*
	$(RM) -r build
