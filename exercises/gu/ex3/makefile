FC = mpifort
FCFLAGS = -c
OBJECTS = geometry.o particle.o tree.o create_mpi_types.o

all: test

test : $(OBJECTS) test.o
	$(FC) -o $@ $(OBJECTS) test.o

particulas : $(OBJECTS) particulas.o
	$(FC) -o $@ $(OBJECTS) particulas.o
	./particulas

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm -f *.o *.mod
