FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particle.o calcs.o ex1.o

all: ex1

ex1: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean: 
	del /Q *.o *.mod

test: ex1
	ex1


