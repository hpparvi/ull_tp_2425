FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particle.o ex2.o
PARTICULAS = geometry.o particle.o particulas.o

all: ex2 particulas

particulas: geometry.o particle.o particulas.o
	$(FC) -o $@ $(PARTICULAS)

ex2 : geometry.o particle.o ex2.o
	$(FC) -o $@ $(OBJECTS)

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm *.o *.mod ex2 particulas

test : ex2
	./ex2
	
init: particulas
	./particulas
