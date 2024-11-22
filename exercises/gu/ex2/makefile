FC = gfortran -fopenmp
FCFLAGS = -c
OBJECTS = geometry.o particle.o tree.o

all: ex2

ex2 : $(OBJECTS) ex2.o
	$(FC) -o $@ $(OBJECTS) ex2.o

try : $(OBJECTS) try.o
	$(FC) -o $@ $(OBJECTS) try.o

particulas : $(OBJECTS) particulas.o
	$(FC) -o $@ $(OBJECTS) particulas.o
	./particulas

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm -f *.o *.mod

test : ex2
	./ex2
