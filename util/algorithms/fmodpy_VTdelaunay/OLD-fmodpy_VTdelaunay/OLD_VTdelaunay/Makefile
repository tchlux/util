RM=rm -f
FORT=gfortran
FORT_FLAGS=-O3 -fPIC
LEGACY= -std=legacy -w
PY=python3
CPY=cython

all:
	$(FORT) -c $(FORT_FLAGS) $(LEGACY) TOMS587.f 
	$(FORT) -c $(FORT_FLAGS) real_precision.f90 Project.f90 VTdelaunay.f90 VTdelaunay_wrapper.f03
	$(CPY) VTdelaunay.pyx
	$(PY) VTdelaunay_setup.py build_ext --inplace

clean:
	$(RM) *.o *.mod *.c *.pyc translated_fortran.*
	$(RM) -r build
