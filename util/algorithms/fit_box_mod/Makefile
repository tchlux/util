RM=rm -f
FORT=gfortran
FORT_FLAGS=-O3 -fPIC
PY=python3
CPY=cython

all:
	which gcc
	$(FORT) -c $(FORT_FLAGS) real_precision.f95 fit_box_mod.f03 fit_box_mod_wrapper.f03
	$(CPY) fit_box_mod.pyx
	$(PY) fit_box_mod_setup.py build_ext --inplace


clean:
	$(RM) *.o *.mod *.c *.pyc translated_fortran.*
	$(RM) -r build
