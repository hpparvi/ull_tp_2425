FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particle.o ex1.o

all: ex1

ex1 : geometry.o particle.o ex1.o
	$(FC) -o $@ $(OBJECTS)

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm *.o *.mod ex1

test : ex1
	./ex1
