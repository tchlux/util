RM=rm -f
FORT=gfortran
FORT_FLAGS=-O3 -fPIC
LEGACY= -std=legacy -w
PY=python3
CPY=cython

all:
	$(FORT) -c $(FORT_FLAGS) $(LEGACY) marspack.f
	$(FORT) -c $(FORT_FLAGS) marspack_wrapper.f03
	$(CPY) marspack.pyx
	$(PY) marspack_setup.py build_ext --inplace

marspack.o: marspack.f
	$(FORT) -c $(FORT_FLAGS) -std=legacy marspack.f

test: marspack.o test_mars.f03
	$(FORT) $(FORT_FLAGS) -o test_mars marspack.o test_mars.f03
	./test_mars

clean:
	$(RM) *.o *.mod *.c *.pyc translated_fortran.*
	$(RM) -r build
