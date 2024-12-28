FC = mpifort
FCFLAGS = -c
OBJECTS = geometry.o particle.o barnes_hut.o ex3.o

all: ex3

ex3 : $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm -f *.o *.mod

test : ex3
	mpirun -n 8 ex3
