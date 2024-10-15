FC = gfortran
FCFLAGS = -c
OBJECTS = geometry.o particle.o ex1.o

all: ex1

ex1: $(OBJECTS)
	$(FC) -O $@ $(OBJECTS)

%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean: 
	del /Q *.o *.mod

test: ex1
	ex1


