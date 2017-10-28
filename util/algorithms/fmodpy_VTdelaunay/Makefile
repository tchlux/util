FORT = gfortran
CFLAGS = -c -fPIC
LEGACY = -std=legacy -w
INC = -llapack -lblas

delaunay: vtdelaunay.o vtgeometry.o project.o toms587.o

all: main.f90 delaunay
	$(FORT) main.f90 vtdelaunay.o vtgeometry.o project.o toms587.o $(INC) -o vtdel

vtdelaunay.o: VTdelaunay.f90 real_precision.mod 
	$(FORT) $(CFLAGS) VTdelaunay.f90 -o vtdelaunay.o

vtgeometry.o: VTgeometry.f90 real_precision.mod
	$(FORT) $(CFLAGS) VTgeometry.f90 -o vtgeometry.o

project.o : Project.f90 real_precision.mod
	$(FORT) $(CFLAGS) Project.f90 -o project.o

toms587.o : TOMS587.f
	$(FORT) $(CFLAGS) $(LEGACY) TOMS587.f -o toms587.o

real_precision.mod : real_precision.f90
	$(FORT) $(CFLAGS) real_precision.f90

clean:
	rm -f *.o *.mod vtdel
