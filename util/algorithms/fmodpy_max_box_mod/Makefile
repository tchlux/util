RM=rm -f
FORT=gfortran
FORT_FLAGS=-O3 -fPIC
PY=python3
CPY=cython

all:
	$(FORT) -c $(FORT_FLAGS) real_precision.f95 max_box_mod.f03 max_box_mod_wrapper.f03
	$(CPY) max_box_mod.pyx
	$(PY) max_box_mod_setup.py build_ext --inplace

clean:
	$(RM) *.o *.mod *.c *.pyc translated_fortran.*
	$(RM) -r build
