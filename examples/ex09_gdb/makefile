FC = gfortran
FCFLAGS = -c
OBJECTS = m1.o ex1.o

all: ex1

test: ex1
	./ex1

ex1: $(OBJECTS)
	$(FC) -o ex1 $(OBJECTS)

%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean:
	rm -f *.o *.mod

allclean:
	rm -f ex1 *.o *.mod *~
