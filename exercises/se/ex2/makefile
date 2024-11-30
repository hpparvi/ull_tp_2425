FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particle.o barnes_hut.o ex2.o

all: ex2

ex2 : $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm -f *.o *.mod

test : ex2
	./ex2
