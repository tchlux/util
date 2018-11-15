FORT = gfortran
CFLAGS = -c -std=f2008
LIBS = -fopenmp 

all: delsparse

delsparse: main.f90 delsparse.o slatec.o lapack.o blas.o
	$(FORT) $(LIBS) main.f90 delsparse.o slatec.o lapack.o blas.o -o delsparse

delsparse.o: delsparse.f90 
	$(FORT) $(CFLAGS) $(LIBS) delsparse.f90 -o delsparse.o

slatec.o : slatec.f
	$(FORT) $(CFLAGS) slatec.f -o slatec.o

lapack.o : lapack.f
	$(FORT) $(CFLAGS) lapack.f -o lapack.o

blas.o : blas.f
	$(FORT) $(CFLAGS) blas.f -o blas.o

clean:
	rm *.o *.mod delsparse
