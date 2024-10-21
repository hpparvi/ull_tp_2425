FC = gfortran
FCFLAGS = -c -fdefault-real-8
OBJECTS = geometry.o particle.o ex1.o

all: ex1

ex1 : $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

%.o : %.f90
	$(FC) $(FCFLAGS) $<

clean :
	rm -f *.o *.mod

test : ex1
	./ex1
