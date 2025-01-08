FC = mpifort
FCFLAGS = -c
OBJECTS = geometry.o particle.o tree.o create_mpi_types.o

all: ex3

ex3 : $(OBJECTS) ex3.o
	$(FC) -o $@ $(OBJECTS) ex3.o

particulas : $(OBJECTS) particulas.o
	$(FC) -o $@ $(OBJECTS) particulas.o
	./particulas

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm -f *.o *.mod
