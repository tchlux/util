RM=rm -f
FORT=gfortran
FORT_FLAGS=-O3 -fPIC
PY=python3
CPY=cython
CFLAGS = -c -fPIC
LEGACY = -std=legacy -w
INC = -llapack -lblas

all: vtdel py

py:
	$(FORT) -c $(FORT_FLAGS) $(LEGACY) TOMS587.f
	$(FORT) -c $(FORT_FLAGS) real_precision.f90 Project.f90 VTgeometry.f90
	$(FORT) -c -w $(FORT_FLAGS) VTdelaunay.f90 VTdelaunay_wrapper.f90
	$(CPY) VTdelaunay.pyx
	$(PY) VTdelaunay_setup.py build_ext --inplace

vtdel: main.f90 vtdelaunay.o vtgeometry.o project.o toms587.o
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
	-$(RM) *.o *.mod vtdel
	-$(RM) *.o *.mod *.c *.pyc translated_fortran.*
	-$(RM) -r build


