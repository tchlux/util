FORT=gfortran
CFLAGS=-c -std=f2008 -fopenmp -fPIC
PYTHON=python3

all: delsparse.cpython-36m-darwin.so

delsparse.cpython-36m-darwin.so: make_module.py delsparse.o slatec.o lapack.o blas.o real_precision.o
	rm -f *.so
	$(PYTHON) make_module.py

delsparse.o: delsparse.f90 
	$(FORT) $(CFLAGS) delsparse.f90 -o delsparse.o

slatec.o : slatec.f
	$(FORT) $(CFLAGS) slatec.f -o slatec.o

lapack.o : lapack.f
	$(FORT) $(CFLAGS) lapack.f -o lapack.o

blas.o : blas.f
	$(FORT) $(CFLAGS) blas.f -o blas.o

real_precision.o: real_precision.f90
	$(FORT) $(CFLAGS) real_precision.f90 -o real_precision.o

clean:
	rm -f *.o *.mod delsparse delsparse.c
	rm -rf build
