FC = gfortran
FCFLAGS = -c

OBJECTS = geometry.o particle.o ex1.o


all: ex1

ex1: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean:
	del /Q *.o *.mod #it would be rm -f on linux but on windows it's del /Q

test: ex1
	ex1