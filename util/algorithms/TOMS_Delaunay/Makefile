FORT = gfortran
CFLAGS = -c

all: vtdel

vtdel: main.f90 vtdelaunay.o slatec.o lapack.o blas.o real_precision.mod
	$(FORT) main.f90 vtdelaunay.o slatec.o lapack.o blas.o -o vtdel

vtdelaunay.o: VTdelaunay.f90 real_precision.mod 
	$(FORT) $(CFLAGS) VTdelaunay.f90 -o vtdelaunay.o

slatec.o : slatec.f
	$(FORT) $(CFLAGS) slatec.f -o slatec.o

lapack.o : lapack.f
	$(FORT) $(CFLAGS) lapack.f -o lapack.o

blas.o : blas.f
	$(FORT) $(CFLAGS) blas.f -o blas.o

real_precision.mod : real_precision.f90
	$(FORT) $(CFLAGS) real_precision.f90

clean:
	rm *.o *.mod vtdel
